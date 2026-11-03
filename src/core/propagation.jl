"""
    PreparedPropagation

A fixed-sampling, envelope-only PTP, WTS or STW complex-field transform.
Construct with `prepare_propagation`. Derived phase factors are run-immutable;
FFT scratch and plans are retained execution storage. Inputs and outputs are
caller-owned. This object does not advance a `WaveFront`'s geometric state.
"""
struct PreparedPropagation{T,A,Q,R,F,B}
    before::Q
    transfer::R
    after::Q
    field::A
    forward::F
    backward::B
    direction::Int8 # 0: identity, 1: forward, -1: backward, 2: PTP pair
    scale::T
end

@inline _require_envelope_propagation(::EnvelopeOnly) = nothing
@inline _require_envelope_propagation(::CarrierPhaseStyle) =
    throw(ArgumentError("prepared field propagation requires EnvelopeOnly"))

"""
    prepare_propagation(wf, dz, ctx; kind=:ptp)

Prepare the field arithmetic of `prop_ptp`, `prop_wts` or `prop_stw` at `wf`'s
wavelength, sampling, precision and square grid. `kind` selects `:ptp`, `:wts`
or `:stw`; STW assumes spherical input and resolves zero distance to the waist.
Distance is converted to field precision.
The wavefront is not changed. Cached phase factors use the ordinary backend
arithmetic outside capture. A qualified direct spatial-phase path instead retains
its coefficients and evaluates the same expression during execution. The context's FFT
scratch/plans are shared: callers must serialize all users of that storage.
Use a new preparation when geometry or backend changes.
"""
function prepare_propagation(wf::WaveFront{T}, dz::Real, ctx::RunContext; kind::Symbol=:ptp) where {T}
    kind in (:ptp, :wts, :stw) || throw(ArgumentError("unknown propagation kind: $kind"))
    _require_envelope_propagation(carrier_phase_style(ctx))
    workspace_float_type(ctx.workspace) === T && typeof(ctx.backend) === typeof(backend_style(typeof(wf.field))) ||
        throw(ArgumentError("prepared propagation context precision/backend differs from wavefront"))
    n = _require_square_propagation_grid(wf, :prepare_propagation)
    d = T(dz)
    kind === :stw && iszero(d) && (d = T(wf.z_w0_m - wf.z_m))
    isfinite(d) && isfinite(wf.sampling_m) && wf.sampling_m > 0 &&
        isfinite(wf.wavelength_m) && wf.wavelength_m > 0 ||
        throw(ArgumentError("prepared propagation requires finite geometry and positive sampling/wavelength"))
    ws = fft_workspace(ctx)
    field = ensure_fft_scratch!(ws, n, n)
    # FFTW.MEASURE can overwrite its planning storage. Preserve a live field
    # when the caller deliberately invalidated plans on an existing scratch.
    saved = field === wf.field && !ws.plans_valid ? copy(wf.field) : nothing
    forward, backward = try
        ensure_fft_plans!(ws, n, n, fft_planning_style(ctx))
    finally
        saved === nothing || copyto!(wf.field, saved)
    end
    empty = similar(field, 0, 0)
    phase_style=Val(_direct_spatial_phase_enabled(field))
    before = after = _empty_spatial_phase(phase_style,field)
    transfer = empty
    scale = inv(T(n))
    direction = Int8(d >= 0 ? 1 : -1)
    if kind === :ptp
        direction = abs(d) < 1e-12 ? Int8(0) : Int8(2)
        if direction != 0
            transfer = similar(field)
            fill!(transfer, one(Complex{T}))
            kphase = -pi * wf.wavelength_m * float(d)
            _prepare_frequency_phase!(backend_style(typeof(transfer)), transfer, kphase, wf.sampling_m, ctx, scale)
        end
    elseif !(kind === :wts && iszero(d))
        dx = kind === :stw ? wf.wavelength_m * abs(d) / (n * wf.sampling_m) : wf.sampling_m
        phase = _prepare_spatial_phase(phase_style,field,wf,d,ctx,dx)
        if kind === :wts
            before = phase
        else
            after = phase
        end
    else
        direction = Int8(0)
    end
    direction == 0 && (scale = one(T))
    return PreparedPropagation(_compact_propagation_phase(before),
        _compact_propagation_phase(transfer), _compact_propagation_phase(after),
        field, forward, backward, direction, scale)
end

function _prepare_frequency_phase!(::CPUBackend, phase::AbstractMatrix{Complex{T}}, k, dx, ctx, scale) where {T}
    return _apply_fft_quadratic_phase_strided!(phase, T(k), T(dx), qphase_workspace(ctx), scale)
end
function _prepare_frequency_phase!(::BackendStyle, phase, k, dx, ctx, scale)
    return _apply_fft_quadratic_phase_ka!(phase, k, dx, qphase_workspace(ctx), scale)
end

# Direct Float32 spatial chirps are cheaper to recompute on the qualified
# CUDA path than to read from a large cache. Keep frequency-domain factors
# cached: their separable Float64 angle evaluation has a different cost.
struct DirectPhaseFactors{T}
    k::T
    dx::T
    active::Bool
end
Base.isempty(phase::DirectPhaseFactors) = !phase.active
_direct_spatial_phase_enabled(::AbstractArray) = false
_empty_spatial_phase(::Val{false},field) = similar(field,0,0)
_empty_spatial_phase(::Val{true},::AbstractMatrix{Complex{T}}) where {T} =
    DirectPhaseFactors(zero(T),one(T),false)
function _prepare_spatial_phase(::Val{false},field,wf,d,ctx,dx)
    phase=similar(field)
    fill!(phase,one(eltype(field)))
    phase_wf=WaveFront(phase,wf.wavelength_m,dx,zero(dx),wf.beam_diameter_m)
    prop_qphase(phase_wf,d,ctx)
    return phase
end
function _prepare_spatial_phase(::Val{true},::AbstractMatrix{Complex{T}},wf,d,ctx,dx) where {T}
    k=iszero(d) ? zero(T) : T(pi/(wf.wavelength_m*d))
    return DirectPhaseFactors(k,T(dx),!iszero(d))
end
@inline function _propagation_phase(phase::DirectPhaseFactors{T},y,x,ny,nx) where {T}
    xx=T(_ka_shifted_index_0based(x-1,nx))*phase.dx
    yy=T(_ka_shifted_index_0based(y-1,ny))*phase.dx
    return one(Complex{T})*cis(phase.k*(xx*xx+yy*yy))
end
_compact_propagation_phase(phase::DirectPhaseFactors) = phase
_propagation_phase_storage(::DirectPhaseFactors) = nothing

# Quadratic phase is even in each raw FFT coordinate separately. Retain the
# nonnegative magnitudes, including the unique even-grid Nyquist coordinate.
# This copies already evaluated factors; no trig, precision or grouping changes.
_propagation_phase_storage(phase::AbstractMatrix) = phase

function _compact_propagation_phase(phase)
    isempty(phase) && return phase
    ny,nx=size(phase)
    return copy(view(phase,1:fld(ny,2)+1,1:fld(nx,2)+1))
end

# Dense factors are retained as the numerical/performance reference. Compact
# quadratic factors use raw-coordinate magnitude; neither changes values.
@inline function _propagation_phase(phase::AbstractMatrix,y,x,ny,nx)
    iy = size(phase,1)==ny ? y : min(y,ny-y+2)
    ix = size(phase,2)==nx ? x : min(x,nx-x+2)
    @inbounds return phase[iy,ix]
end

@kernel function _propagation_transfer_kernel!(field,phase)
    y,x=@index(Global,NTuple)
    ny,nx=size(field)
    if y<=ny && x<=nx
        @inbounds field[y,x] *= _propagation_phase(phase,y,x,ny,nx)
    end
end

@inline _propagation_input(value, ::Nothing, y, x) = value
@inline _propagation_input(value, multiplier, y, x) = value * multiplier[y, x]

# Checked equal one-based shapes and a single wrap bound every access. Phase
# arrays are empty or the reflected half-axis shape constructed above.
@kernel function _propagation_pack_kernel!(field, input, multiplier, phase, shift)
    y, x = @index(Global, NTuple)
    ny, nx = size(field)
    if y <= ny && x <= nx
        sy, sx = y + shift, x + shift
        sy > ny && (sy -= ny)
        sx > nx && (sx -= nx)
        @inbounds begin
            value = _propagation_input(input[sy, sx], multiplier, sy, sx)
            field[y, x] = isempty(phase) ? value : value * _propagation_phase(phase,y,x,ny,nx)
        end
    end
end

@kernel function _propagation_unpack_kernel!(out, field, phase, scale, multiplier, shift)
    y, x = @index(Global, NTuple)
    ny, nx = size(field)
    if y <= ny && x <= nx
        sy, sx = y + shift, x + shift
        sy > ny && (sy -= ny)
        sx > nx && (sx -= nx)
        @inbounds begin
            value = field[sy, sx] * scale
            value = isempty(phase) ? value : value * _propagation_phase(phase,sy,sx,ny,nx)
            out[y, x] = value * multiplier
        end
    end
end

"""
    propagate!(out, prepared, input; centered=false, input_multiplier=nothing,
               output_multiplier=1)

Enqueue one prepared field transform. `centered=true` accepts and returns
centered `(y,x)` fields; otherwise both use raw FFT order. Optional input and
output multipliers run before and after the transform, respectively. The sequence of edge,
quadratic phase, FFT scaling and output multiplication is retained; backend
fusion and cached phase evaluation can differ from ordinary execution by roundoff. Input/output may alias each other, but must not alias FFT scratch.
All arrays must have the prepared shape, precision, backend and one-based axes.
Synchronize at the caller's completion boundary; do not overlap shared scratch.
"""
function propagate!(out, p::PreparedPropagation{T}, input;
    centered::Bool=false, input_multiplier=nothing, output_multiplier=one(Complex{T})) where {T}
    _check_propagation_fields(out, input, p)
    _check_propagation_multiplier(input_multiplier, input, p.field)
    backend = get_backend(p.field)
    n = size(p.field, 1)
    _propagation_pack_kernel!(backend, (16,16))(p.field, input, input_multiplier,
        p.before, centered ? fld(n,2) : 0; ndrange=size(input))
    _propagation_fft!(p)
    _propagation_unpack_kernel!(backend, (16,16))(out, p.field, p.after, p.scale,
        Complex{T}(output_multiplier), centered ? cld(n,2) : 0; ndrange=size(out))
    return out
end

@inline _check_propagation_multiplier(::Nothing, input, scratch) = nothing
@inline function _check_propagation_multiplier(multiplier, input, scratch)
    Base.mightalias(multiplier,scratch) && throw(ArgumentError("prepared propagation multiplier aliases FFT scratch"))
    axes(multiplier) == axes(input) || throw(DimensionMismatch("prepared propagation multiplier axes differ"))
    eltype(multiplier) === eltype(input) && typeof(get_backend(multiplier)) === typeof(get_backend(input)) ||
        throw(ArgumentError("prepared propagation multiplier precision/backend differs"))
    return nothing
end

# Join adjacent transforms while their scratch remains in raw FFT order.
# Each arithmetic step matches an unpack followed by the next pack. The shifted
# multiplier is read at the corresponding centered coordinate, including odd n.
@kernel function _propagation_join_kernel!(field, after, scale, output_multiplier,
    input_multiplier, before, shift)
    y, x = @index(Global, NTuple)
    ny, nx = size(field)
    if y <= ny && x <= nx
        sy, sx = y + shift, x + shift
        sy > ny && (sy -= ny)
        sx > nx && (sx -= nx)
        @inbounds begin
            value = field[y,x] * scale
            value = isempty(after) ? value : value * _propagation_phase(after,y,x,ny,nx)
            value *= output_multiplier
            value = _propagation_input(value, input_multiplier, sy, sx)
            field[y,x] = isempty(before) ? value : value * _propagation_phase(before,y,x,ny,nx)
        end
    end
end

"""
    propagate!(out, (first, second), input; centered=false,
        input_multipliers=(nothing,nothing), output_multipliers=(1,1))

Execute two adjacent prepared transforms using the same FFT scratch. This
preserves each transform's multiplications and scaling, while combining the
first output mapping and second input mapping into one raw-order field pass.
Input/output and scratch contracts are the same as the single-transform method.
"""
function propagate!(out, pair::Tuple{P,P}, input; centered::Bool=false,
    input_multipliers=(nothing,nothing), output_multipliers=(1,1)) where {T,P<:PreparedPropagation{T}}
    first, second = pair
    first.field === second.field || throw(ArgumentError("adjacent propagation requires shared FFT scratch"))
    _check_propagation_fields(out, input, first)
    _check_propagation_multiplier(input_multipliers[1], input, first.field)
    _check_propagation_multiplier(input_multipliers[2], input, first.field)
    backend = get_backend(first.field)
    shift = centered ? fld(size(first.field,1),2) : 0
    _propagation_pack_kernel!(backend, (16,16))(first.field, input, input_multipliers[1],
        first.before, shift; ndrange=size(input))
    _propagation_fft!(first)
    _propagation_join_kernel!(backend, (16,16))(first.field, first.after, first.scale,
        Complex{T}(output_multipliers[1]), input_multipliers[2], second.before, shift; ndrange=size(input))
    _propagation_fft!(second)
    _propagation_unpack_kernel!(backend, (16,16))(out, second.field, second.after, second.scale,
        Complex{T}(output_multipliers[2]), centered ? cld(size(first.field,1),2) : 0; ndrange=size(out))
    return out
end

function _check_propagation_fields(out, input, p)
    size(out) == size(input) == size(p.field) || throw(DimensionMismatch("prepared propagation field shapes differ"))
    Base.require_one_based_indexing(out, input)
    eltype(out) === eltype(input) === eltype(p.field) ||
        throw(ArgumentError("prepared propagation field precisions differ"))
    typeof(get_backend(out)) === typeof(get_backend(input)) === typeof(get_backend(p.field)) ||
        throw(ArgumentError("prepared propagation field backends differ"))
    (Base.mightalias(out, p.field) || Base.mightalias(input, p.field)) &&
        throw(ArgumentError("prepared propagation requires separate FFT scratch"))
    return nothing
end

@inline function _propagation_fft!(p)
    if p.direction == 2
        p.forward * p.field
        _propagation_transfer_kernel!(get_backend(p.field),(16,16))(
            p.field,p.transfer;ndrange=size(p.field))
        p.backward * p.field
    elseif p.direction == 1
        p.forward * p.field
    elseif p.direction == -1
        p.backward * p.field
    end
    return nothing
end
