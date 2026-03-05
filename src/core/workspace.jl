mutable struct InterpWorkspace{T<:AbstractFloat}
    xcoords::Vector{T}
    ycoords::Vector{T}
end

InterpWorkspace(::Type{T}=Float64) where {T<:AbstractFloat} =
    InterpWorkspace{T}(Vector{T}(undef, 0), Vector{T}(undef, 0))

@inline function ensure_interp_axes!(ws::InterpWorkspace, nx::Integer, ny::Integer)
    resize!(ws.xcoords, nx)
    resize!(ws.ycoords, ny)
    return ws.xcoords, ws.ycoords
end

mutable struct MaskWorkspace{T<:AbstractFloat}
    mask::Matrix{T}
    valid::Bool
    xverts::Vector{T}
    yverts::Vector{T}
end

MaskWorkspace(::Type{T}=Float64) where {T<:AbstractFloat} =
    MaskWorkspace{T}(Matrix{T}(undef, 0, 0), false, Vector{T}(undef, 0), Vector{T}(undef, 0))

@inline function ensure_mask_buffer!(ws::MaskWorkspace{T}, ny::Integer, nx::Integer) where {T}
    if !(ws.valid && size(ws.mask) == (ny, nx))
        ws.mask = Matrix{T}(undef, ny, nx)
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

mutable struct FFTWorkspace{T<:AbstractFloat}
    rho2::Matrix{T}
    scratch::Matrix{Complex{T}}
    real_scratch::Matrix{T}
    forward_plan::Union{Nothing,FFTWFwdPlan2D{T}}
    backward_plan::Union{Nothing,FFTWBwdPlan2D{T}}
    nx::Int
    ny::Int
    dx::T
    valid::Bool
    scratch_valid::Bool
    real_scratch_valid::Bool
    plans_valid::Bool
end

FFTWorkspace(::Type{T}=Float64) where {T<:AbstractFloat} =
    FFTWorkspace{T}(
        Matrix{T}(undef, 0, 0),
        Matrix{Complex{T}}(undef, 0, 0),
        Matrix{T}(undef, 0, 0),
        nothing,
        nothing,
        0,
        0,
        zero(T),
        false,
        false,
        false,
        false,
    )

@inline function ensure_rho2_map!(ws::FFTWorkspace{T}, ny::Integer, nx::Integer, dx::Real) where {T}
    dxT = T(dx)
    if !(ws.valid && ws.nx == nx && ws.ny == ny && ws.dx == dxT)
        ws.rho2 = fft_order_rho2_map(ny, nx, dxT)
        ws.nx = nx
        ws.ny = ny
        ws.dx = dxT
        ws.valid = true
    end
    return ws.rho2
end

@inline function ensure_fft_scratch!(ws::FFTWorkspace{T}, ny::Integer, nx::Integer) where {T}
    if !(ws.scratch_valid && size(ws.scratch) == (ny, nx))
        ws.scratch = Matrix{Complex{T}}(undef, ny, nx)
        ws.scratch_valid = true
        ws.plans_valid = false
    end
    return ws.scratch
end

@inline function ensure_fft_real_scratch!(ws::FFTWorkspace{T}, ny::Integer, nx::Integer) where {T}
    if !(ws.real_scratch_valid && size(ws.real_scratch) == (ny, nx))
        ws.real_scratch = Matrix{T}(undef, ny, nx)
        ws.real_scratch_valid = true
    end
    return ws.real_scratch
end

@inline function ensure_fft_plans!(ws::FFTWorkspace{T}, ny::Integer, nx::Integer) where {T}
    ensure_fft_scratch!(ws, ny, nx)
    if !(ws.plans_valid && ws.forward_plan !== nothing && ws.backward_plan !== nothing)
        ws.forward_plan = FFTW.plan_fft!(ws.scratch; flags=FFTW.ESTIMATE)
        ws.backward_plan = FFTW.plan_bfft!(ws.scratch; flags=FFTW.ESTIMATE)
        ws.plans_valid = true
    end
    return ws.forward_plan::FFTWFwdPlan2D{T}, ws.backward_plan::FFTWBwdPlan2D{T}
end

@inline function reset_workspace!(ws::InterpWorkspace)
    resize!(ws.xcoords, 0)
    resize!(ws.ycoords, 0)
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
    return ws
end

mutable struct ProperWorkspace{T<:AbstractFloat}
    interp::InterpWorkspace{T}
    mask::MaskWorkspace{T}
    fft::FFTWorkspace{T}
end

ProperWorkspace(::Type{T}=Float64) where {T<:AbstractFloat} =
    ProperWorkspace{T}(InterpWorkspace(T), MaskWorkspace(T), FFTWorkspace(T))

@inline function reset_workspace!(ws::ProperWorkspace)
    reset_workspace!(ws.interp)
    reset_workspace!(ws.mask)
    reset_workspace!(ws.fft)
    return ws
end
