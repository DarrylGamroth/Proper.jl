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
    nx::Int
    ny::Int
    dx::T
    valid::Bool
end

FFTWorkspace(::Type{T}=Float64) where {T<:AbstractFloat} =
    FFTWorkspace{T}(Matrix{T}(undef, 0, 0), 0, 0, zero(T), false)

@inline function ensure_rho2_map!(ws::FFTWorkspace{T}, nx::Integer, ny::Integer, dx::Real) where {T}
    dxT = T(dx)
    if !(ws.valid && ws.nx == nx && ws.ny == ny && ws.dx == dxT)
        ws.rho2 = fft_order_rho2_map(nx, ny, dxT)
        ws.nx = nx
        ws.ny = ny
        ws.dx = dxT
        ws.valid = true
    end
    return ws.rho2
end

mutable struct ProperWorkspace{T<:AbstractFloat}
    interp::InterpWorkspace{T}
    fft::FFTWorkspace{T}
end

ProperWorkspace(::Type{T}=Float64) where {T<:AbstractFloat} =
    ProperWorkspace{T}(InterpWorkspace(T), FFTWorkspace(T))
