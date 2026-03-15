@inline workspace_vector(::Type{A}, ::Type{T}, n::Integer=0) where {A<:AbstractArray,T<:AbstractFloat} =
    Vector{T}(undef, n)

@inline workspace_backend_array_type(::Type{A}) where {A<:AbstractArray} = A

@inline workspace_matrix(::Type{A}, ::Type{T}, ny::Integer=0, nx::Integer=0) where {A<:AbstractArray,T<:AbstractFloat} =
    Matrix{T}(undef, ny, nx)

@inline workspace_complex_matrix(::Type{A}, ::Type{T}, ny::Integer=0, nx::Integer=0) where {A<:AbstractArray,T<:AbstractFloat} =
    Matrix{Complex{T}}(undef, ny, nx)

@inline function _ensure_workspace_vector(a::Vector{T}, n::Integer) where {T}
    length(a) == n || resize!(a, n)
    return a
end

@inline function _ensure_workspace_vector(a::AbstractVector{T}, n::Integer) where {T}
    size(a, 1) == n && return a
    return similar(a, T, n)
end

@inline function _ensure_workspace_matrix(a::Matrix{T}, ny::Integer, nx::Integer) where {T}
    size(a) == (ny, nx) && return a
    return Matrix{T}(undef, ny, nx)
end

@inline function _ensure_workspace_matrix(a::AbstractMatrix{T}, ny::Integer, nx::Integer) where {T}
    size(a) == (ny, nx) && return a
    return similar(a, T, ny, nx)
end

@inline function _fft_shifted_index_0based(p::Int, n::Int)
    c = n ÷ 2
    cut = n - c - 1
    return p <= cut ? p : (p - n)
end

mutable struct InterpWorkspace{T<:AbstractFloat,VX<:AbstractVector{T},VY<:AbstractVector{T}}
    xcoords::VX
    ycoords::VY
end

InterpWorkspace(::Type{T}=Float64) where {T<:AbstractFloat} = InterpWorkspace(Matrix, T)

function InterpWorkspace(::Type{A}, ::Type{T}=Float64) where {A<:AbstractArray,T<:AbstractFloat}
    xcoords = workspace_vector(A, T, 0)
    ycoords = workspace_vector(A, T, 0)
    return InterpWorkspace{T,typeof(xcoords),typeof(ycoords)}(xcoords, ycoords)
end

@inline function ensure_interp_axes!(ws::InterpWorkspace, nx::Integer, ny::Integer)
    ws.xcoords = _ensure_workspace_vector(ws.xcoords, nx)
    ws.ycoords = _ensure_workspace_vector(ws.ycoords, ny)
    return ws.xcoords, ws.ycoords
end

mutable struct MaskWorkspace{T<:AbstractFloat,M<:AbstractMatrix{T}}
    mask::M
    valid::Bool
    xverts::Vector{T}
    yverts::Vector{T}
end

MaskWorkspace(::Type{T}=Float64) where {T<:AbstractFloat} = MaskWorkspace(Matrix, T)

function MaskWorkspace(::Type{A}, ::Type{T}=Float64) where {A<:AbstractArray,T<:AbstractFloat}
    mask = workspace_matrix(A, T, 0, 0)
    return MaskWorkspace{T,typeof(mask)}(mask, false, Vector{T}(undef, 0), Vector{T}(undef, 0))
end

@inline function ensure_mask_buffer!(ws::MaskWorkspace{T}, ny::Integer, nx::Integer) where {T}
    if !(ws.valid && size(ws.mask) == (ny, nx))
        ws.mask = _ensure_workspace_matrix(ws.mask, ny, nx)
        ws.valid = true
    end
    return ws.mask
end

@inline function ensure_mask_vertices!(ws::MaskWorkspace{T}, nverts::Integer) where {T}
    resize!(ws.xverts, nverts)
    resize!(ws.yverts, nverts)
    return ws.xverts, ws.yverts
end

const FFTWFwdPlan2D{T} = FFTW.cFFTWPlan{Complex{T},-1,true,2,Tuple{Int,Int}}
const FFTWBwdPlan2D{T} = FFTW.cFFTWPlan{Complex{T},1,true,2,Tuple{Int,Int}}
const MaybeFFTWFwdPlan2D{T} = Union{Nothing,FFTWFwdPlan2D{T}}
const MaybeFFTWBwdPlan2D{T} = Union{Nothing,FFTWBwdPlan2D{T}}

mutable struct FFTWorkspace{
    T<:AbstractFloat,
    SCR<:AbstractMatrix{Complex{T}},
    RS<:AbstractMatrix{T},
    FP,
    BP,
}
    rho2::Matrix{T}
    scratch::SCR
    real_scratch::RS
    forward_plan::FP
    backward_plan::BP
    nx::Int
    ny::Int
    dx::T
    plan_flags::UInt32
    valid::Bool
    scratch_valid::Bool
    real_scratch_valid::Bool
    plans_valid::Bool
end

function FFTWorkspace(::Type{T}=Float64) where {T<:AbstractFloat}
    scratch = Matrix{Complex{T}}(undef, 0, 0)
    real_scratch = Matrix{T}(undef, 0, 0)
    return FFTWorkspace{T,typeof(scratch),typeof(real_scratch),MaybeFFTWFwdPlan2D{T},MaybeFFTWBwdPlan2D{T}}(
        Matrix{T}(undef, 0, 0),
        scratch,
        real_scratch,
        nothing,
        nothing,
        0,
        0,
        zero(T),
        zero(UInt32),
        false,
        false,
        false,
        false,
    )
end

function FFTWorkspace(::Type{A}, ::Type{T}=Float64) where {A<:AbstractArray,T<:AbstractFloat}
    scratch = workspace_complex_matrix(A, T, 0, 0)
    real_scratch = workspace_matrix(A, T, 0, 0)
    return FFTWorkspace{T,typeof(scratch),typeof(real_scratch),MaybeFFTWFwdPlan2D{T},MaybeFFTWBwdPlan2D{T}}(
        Matrix{T}(undef, 0, 0),
        scratch,
        real_scratch,
        nothing,
        nothing,
        0,
        0,
        zero(T),
        zero(UInt32),
        false,
        false,
        false,
        false,
    )
end

@inline function ensure_rho2_map!(ws::FFTWorkspace{T}, ny::Integer, nx::Integer, dx::Real) where {T}
    dxT = T(dx)
    if !(ws.valid && ws.nx == nx && ws.ny == ny && ws.dx == dxT)
        ws.rho2 = _ensure_workspace_matrix(ws.rho2, ny, nx)
        _fill_fft_order_rho2!(CPUBackend(), ws.rho2, dxT)
        ws.nx = nx
        ws.ny = ny
        ws.dx = dxT
        ws.valid = true
    end
    return ws.rho2
end

@inline function ensure_fft_scratch!(ws::FFTWorkspace{T}, ny::Integer, nx::Integer) where {T}
    if !(ws.scratch_valid && size(ws.scratch) == (ny, nx))
        ws.scratch = _ensure_workspace_matrix(ws.scratch, ny, nx)
        ws.scratch_valid = true
        ws.plans_valid = false
    end
    return ws.scratch
end

@inline function ensure_fft_real_scratch!(ws::FFTWorkspace{T}, ny::Integer, nx::Integer) where {T}
    if !(ws.real_scratch_valid && size(ws.real_scratch) == (ny, nx))
        ws.real_scratch = _ensure_workspace_matrix(ws.real_scratch, ny, nx)
        ws.real_scratch_valid = true
    end
    return ws.real_scratch
end

@inline fftw_flags(::FFTEstimateStyle) = FFTW.ESTIMATE
@inline fftw_flags(::FFTMeasureStyle) = FFTW.MEASURE
@inline fftw_flags(::FFTPlanningStyle) = FFTW.ESTIMATE

abstract type FFTPlanExecStyle end
struct FFTPlanAvailableStyle <: FFTPlanExecStyle end
struct FFTPlanUnavailableStyle <: FFTPlanExecStyle end

@inline fft_plan_exec_style(::Type{<:StridedMatrix}) = FFTPlanAvailableStyle()
@inline fft_plan_exec_style(::Type{<:AbstractMatrix}) = FFTPlanUnavailableStyle()

@inline function _ensure_fft_plans!(
    ::FFTPlanAvailableStyle,
    ws::FFTWorkspace{T},
    ny::Integer,
    nx::Integer,
    planning::FFTPlanningStyle,
) where {T}
    ensure_fft_scratch!(ws, ny, nx)
    flags = fftw_flags(planning)
    if !(ws.plans_valid && ws.forward_plan !== nothing && ws.backward_plan !== nothing && ws.plan_flags == flags)
        ws.forward_plan = FFTW.plan_fft!(ws.scratch; flags=flags)
        ws.backward_plan = FFTW.plan_bfft!(ws.scratch; flags=flags)
        ws.plan_flags = flags
        ws.plans_valid = true
    end
    return ws.forward_plan::FFTWFwdPlan2D{T}, ws.backward_plan::FFTWBwdPlan2D{T}
end

@inline function ensure_fft_plans!(ws::FFTWorkspace, ny::Integer, nx::Integer, planning::FFTPlanningStyle=FFTEstimateStyle())
    return _ensure_fft_plans!(fft_plan_exec_style(typeof(ws.scratch)), ws, ny, nx, planning)
end

@inline function _ensure_fft_plans!(
    ::FFTPlanUnavailableStyle,
    ws::FFTWorkspace,
    ny::Integer,
    nx::Integer,
    planning::FFTPlanningStyle,
)
    _ = ws
    _ = ny
    _ = nx
    _ = planning
    throw(ArgumentError("FFT plans are only available for strided CPU workspaces"))
end

@inline function _fill_fft_order_rho2!(::CPUBackend, rho2::AbstractMatrix{T}, dx::T) where {T<:AbstractFloat}
    ny, nx = size(rho2)
    inv_dy = inv(T(ny) * dx)
    inv_dx = inv(T(nx) * dx)
    @inbounds for j in 1:nx
        fx = T(_fft_shifted_index_0based(j - 1, nx)) * inv_dx
        fx2 = fx * fx
        for i in 1:ny
            fy = T(_fft_shifted_index_0based(i - 1, ny)) * inv_dy
            rho2[i, j] = fx2 + fy * fy
        end
    end
    return rho2
end

@inline function fill_affine_axis!(
    ::AxisFillLoopExecStyle,
    out::AbstractVector{T},
    origin::T,
    scale::T,
    offset::T,
) where {T<:AbstractFloat}
    @inbounds for i in eachindex(out)
        out[i] = (T(i - 1) - origin) * scale + offset
    end
    return out
end

@inline function fill_affine_axis!(
    ::AxisFillKAExecStyle,
    out::AbstractVector,
    origin,
    scale,
    offset,
)
    return ka_fill_affine_axis!(out, origin, scale, offset)
end

@inline function reset_workspace!(ws::InterpWorkspace)
    return ws
end

@inline function reset_workspace!(ws::MaskWorkspace)
    ws.valid = false
    resize!(ws.xverts, 0)
    resize!(ws.yverts, 0)
    return ws
end

@inline function reset_workspace!(ws::FFTWorkspace)
    ws.valid = false
    ws.scratch_valid = false
    ws.real_scratch_valid = false
    ws.plans_valid = false
    ws.forward_plan = nothing
    ws.backward_plan = nothing
    ws.nx = 0
    ws.ny = 0
    ws.dx = zero(eltype(ws.rho2))
    ws.plan_flags = zero(UInt32)
    return ws
end

mutable struct ProperWorkspace{
    T<:AbstractFloat,
    IW<:InterpWorkspace{T},
    MW<:MaskWorkspace{T},
    FW<:FFTWorkspace{T},
}
    interp::IW
    mask::MW
    fft::FW
end

ProperWorkspace(::Type{T}=Float64) where {T<:AbstractFloat} = ProperWorkspace(Matrix, T)

function ProperWorkspace(::Type{A}, ::Type{T}=Float64) where {A<:AbstractArray,T<:AbstractFloat}
    interp = InterpWorkspace(A, T)
    mask = MaskWorkspace(A, T)
    fft = FFTWorkspace(A, T)
    return ProperWorkspace{T,typeof(interp),typeof(mask),typeof(fft)}(interp, mask, fft)
end

@inline function reset_workspace!(ws::ProperWorkspace)
    reset_workspace!(ws.interp)
    reset_workspace!(ws.mask)
    reset_workspace!(ws.fft)
    return ws
end

@inline workspace_backend_array_type(ws::ProperWorkspace) = workspace_backend_array_type(typeof(field_backend_template(ws)))

@inline function fresh_workspace(ws::ProperWorkspace{T}) where {T<:AbstractFloat}
    return ProperWorkspace(workspace_backend_array_type(ws), T)
end
