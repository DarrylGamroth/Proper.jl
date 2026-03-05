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

mutable struct FFTWorkspace{T<:AbstractFloat}
    rho2::Matrix{T}
    scratch::Matrix{Complex{T}}
    nx::Int
    ny::Int
    dx::T
    valid::Bool
    scratch_valid::Bool
end

FFTWorkspace(::Type{T}=Float64) where {T<:AbstractFloat} =
    FFTWorkspace{T}(Matrix{T}(undef, 0, 0), Matrix{Complex{T}}(undef, 0, 0), 0, 0, zero(T), false, false)

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
    end
    return ws.scratch
end

mutable struct ProperWorkspace{T<:AbstractFloat}
    interp::InterpWorkspace{T}
    fft::FFTWorkspace{T}
end

ProperWorkspace(::Type{T}=Float64) where {T<:AbstractFloat} =
    ProperWorkspace{T}(InterpWorkspace(T), FFTWorkspace(T))
