"""Thread-capable cubic convolution entrypoint (currently delegates to scalar kernel)."""
libcconvthread(a::AbstractMatrix, y::Real, x::Real) = libcconv(a, y, x)
