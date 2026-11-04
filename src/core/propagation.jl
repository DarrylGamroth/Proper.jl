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
@inline _propagation_input(value, ::Tuple{}, y, x) = value
@inline _propagation_input(value, multipliers::Tuple, y, x) =
    _propagation_input(_propagation_input(value, first(multipliers), y, x),
        Base.tail(multipliers), y, x)

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

# Independent field batches retain the scalar kernels' multiplication order.
# The third axis labels fields, not a Fourier coordinate. Callers validate
# exactly two square planes and tuple entries, precision/backend and disjoint outputs;
# FFT plans must transform only axes (1,2). Common phase factors are evaluated
# once per coordinate and reused across the independent planes.
@kernel function _propagation_batch_pack_kernel!(field,inputs,masks,edge,before,shift)
 y,x=@index(Global,NTuple);n=size(field,1)
 if y<=n && x<=n
  sy,sx=y+shift,x+shift;sy>n && (sy-=n);sx>n && (sx-=n)
  phase=isempty(before) ? one(eltype(field)) : _propagation_phase(before,y,x,n,n)
  @inbounds for z in 1:2
   value=_propagation_input(inputs[z][sy,sx],masks[z],sy,sx)
   value=_propagation_input(value,edge,sy,sx)
   field[y,x,z]=isempty(before) ? value : value*phase
  end
 end
end
@kernel function _propagation_batch_join_kernel!(field,after,scale,gain,masks,edge,before,shift)
 y,x=@index(Global,NTuple);n=size(field,1)
 if y<=n && x<=n
  sy,sx=y+shift,x+shift;sy>n && (sy-=n);sx>n && (sx-=n)
  post=isempty(after) ? one(eltype(field)) : _propagation_phase(after,y,x,n,n)
  pre=isempty(before) ? one(eltype(field)) : _propagation_phase(before,y,x,n,n)
  @inbounds for z in 1:2
   value=field[y,x,z]*scale
   value=isempty(after) ? value : value*post
   value*=gain
   value=_propagation_input(value,masks[z],sy,sx)
   value=_propagation_input(value,edge,sy,sx)
   field[y,x,z]=isempty(before) ? value : value*pre
  end
 end
end
@kernel function _propagation_batch_transfer_kernel!(field,transfer)
 y,x=@index(Global,NTuple);n=size(field,1)
 if y<=n && x<=n
  phase=_propagation_phase(transfer,y,x,n,n)
  @inbounds for z in 1:2;field[y,x,z]*=phase;end
 end
end
@kernel function _propagation_batch_unpack_kernel!(outputs,field,after,scale,gain,shift)
 y,x=@index(Global,NTuple);n=size(field,1)
 if y<=n && x<=n
  sy,sx=y+shift,x+shift;sy>n && (sy-=n);sx>n && (sx-=n)
  phase=isempty(after) ? one(eltype(field)) : _propagation_phase(after,sy,sx,n,n)
  @inbounds for z in 1:2
   value=field[sy,sx,z]*scale
   value=isempty(after) ? value : value*phase
   outputs[z][y,x]=value*gain
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
All arrays must have the prepared shape, backend and one-based axes. Input
multipliers may be real or complex at the field's component precision.
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
@inline function _check_propagation_multiplier(multipliers::Tuple, input, scratch)
    foreach(multiplier -> _check_propagation_multiplier(multiplier,input,scratch),multipliers)
    return nothing
end
@inline function _check_propagation_multiplier(multiplier, input, scratch)
    Base.mightalias(multiplier,scratch) && throw(ArgumentError("prepared propagation multiplier aliases FFT scratch"))
    axes(multiplier) == axes(input) || throw(DimensionMismatch("prepared propagation multiplier axes differ"))
    eltype(multiplier) in (eltype(input),real(eltype(input))) &&
        typeof(get_backend(multiplier)) === typeof(get_backend(input)) ||
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
    propagate!(out, transforms::Tuple, input; centered=false,
        input_multipliers=map(_ -> nothing, transforms),
        output_multipliers=map(_ -> 1, transforms))

Execute adjacent prepared transforms using the same FFT scratch. This
preserves each transform's multiplications and scaling, while combining each
output mapping and following input mapping into one raw-order field pass.
An input multiplier may be a tuple of sampled multipliers, applied in order.
The outer input/output multiplier collections must be tuples.
Input/output and scratch contracts are the same as the single-transform method.
"""
function propagate!(out, transforms::Tuple{P,Vararg{P}}, input; centered::Bool=false,
    input_multipliers::Tuple=map(_ -> nothing,transforms),
    output_multipliers::Tuple=map(_ -> 1,transforms)) where {T,P<:PreparedPropagation{T}}
    length(input_multipliers)==length(output_multipliers)==length(transforms) ||
        throw(DimensionMismatch("prepared propagation multiplier counts differ"))
    initial=first(transforms)
    all(p -> p.field===initial.field,transforms) ||
        throw(ArgumentError("adjacent propagation requires shared FFT scratch"))
    _check_propagation_fields(out,input,initial)
    foreach(multiplier -> _check_propagation_multiplier(multiplier,input,initial.field),input_multipliers)
    backend=get_backend(initial.field)
    shift=centered ? fld(size(initial.field,1),2) : 0
    _propagation_pack_kernel!(backend,(16,16))(initial.field,input,first(input_multipliers),
        initial.before,shift;ndrange=size(input))
    _propagation_fft!(initial)
    _propagation_continue!(out,initial,Base.tail(transforms),Base.tail(input_multipliers),
        output_multipliers,centered,shift)
    return out
end

@inline function _propagation_continue!(out,previous,::Tuple{},::Tuple{},gains,centered,shift)
    _propagation_unpack_kernel!(get_backend(previous.field),(16,16))(
        out,previous.field,previous.after,previous.scale,eltype(out)(first(gains)),
        centered ? cld(size(previous.field,1),2) : 0;ndrange=size(out))
    return nothing
end
@inline function _propagation_continue!(out,previous,remaining::Tuple,inputs,gains,centered,shift)
    next=first(remaining)
    _propagation_join_kernel!(get_backend(previous.field),(16,16))(
        previous.field,previous.after,previous.scale,eltype(out)(first(gains)),
        first(inputs),next.before,shift;ndrange=size(out))
    _propagation_fft!(next)
    return _propagation_continue!(out,next,Base.tail(remaining),Base.tail(inputs),
        Base.tail(gains),centered,shift)
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

# Private raw-order train execution: the live wavefront is exactly the shared
# FFT scratch. Zero shifts make every elementwise in-place access independent.
# The public centered/input-output API deliberately keeps its stricter alias rule.
function _propagate_inplace!(prepared::Tuple{P,Vararg{P}},inputs::Tuple,tail::Tuple,empty) where {P<:PreparedPropagation}
    length(prepared)==length(inputs) || throw(ArgumentError("prepared train multiplier count differs"))
    initial=first(prepared);field=initial.field
    isempty(empty) || throw(ArgumentError("prepared train terminal before factors must be empty"))
    all(p->p.field===field,prepared) || throw(ArgumentError("prepared train requires identical FFT storage"))
    foreach(input->_check_propagation_multiplier(input,field,field),inputs)
    _check_propagation_multiplier(tail,field,field)
    if !isempty(first(inputs)) || !isempty(initial.before)
        _propagation_pack_kernel!(get_backend(field),(16,16))(
            field,field,first(inputs),initial.before,0;ndrange=size(field))
    end
    _propagation_fft!(initial)
    _propagate_inplace_continue!(initial,Base.tail(prepared),Base.tail(inputs),tail,empty)
    return field
end
@inline function _propagate_inplace_continue!(previous,remaining::Tuple,inputs,tail,empty)
    next=first(remaining);field=previous.field
    _propagation_join_kernel!(get_backend(field),(16,16))(
        field,previous.after,previous.scale,one(eltype(field)),first(inputs),next.before,0;ndrange=size(field))
    _propagation_fft!(next)
    return _propagate_inplace_continue!(next,Base.tail(remaining),Base.tail(inputs),tail,empty)
end
@inline function _propagate_inplace_continue!(previous,::Tuple{},::Tuple{},tail,empty)
    field=previous.field
    # No next before phase follows the final ordered tail multipliers.
    _propagation_join_kernel!(get_backend(field),(16,16))(
        field,previous.after,previous.scale,one(eltype(field)),tail,empty,0;
        ndrange=size(field))
    return nothing
end
