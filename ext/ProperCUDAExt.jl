module ProperCUDAExt

using Proper
using CUDA
using CUDA.CUFFT

Proper.backend_style(::Type{<:CUDA.CuArray}) = Proper.CUDABackend()
Proper.fft_style(::Type{<:CUDA.CuArray}) = Proper.CUFFTStyle()
Proper.interp_style(::Type{<:CUDA.CuArray}) = Proper.CubicInterpStyle()
Proper.shift_kernel_style(::Type{<:CUDA.CuArray}) = Proper.ShiftKAStyle()
Proper.geometry_kernel_style(::Type{<:CUDA.CuArray}) = Proper.GeometryKAStyle()
Proper.sampling_kernel_style(::Type{<:CUDA.CuArray}) = Proper.SamplingKAStyle()
Proper.interp_kernel_style(::Type{<:CUDA.CuArray}) = Proper.InterpKAStyle()
Proper.workspace_vector(::Type{<:CUDA.CuArray}, ::Type{T}, n::Integer=0) where {T<:AbstractFloat} = CUDA.CuVector{T}(undef, n)
Proper.workspace_matrix(::Type{<:CUDA.CuArray}, ::Type{T}, ny::Integer=0, nx::Integer=0) where {T<:AbstractFloat} = CUDA.CuMatrix{T}(undef, ny, nx)
Proper.workspace_complex_matrix(::Type{<:CUDA.CuArray}, ::Type{T}, ny::Integer=0, nx::Integer=0) where {T<:AbstractFloat} = CUDA.CuMatrix{Complex{T}}(undef, ny, nx)

@inline Proper.ka_mask_enabled(::Type{A}, ny::Integer, nx::Integer) where {A<:CUDA.CuArray} = ny > 0 && nx > 0
@inline Proper.ka_end_enabled(::Type{A}, ny::Integer, nx::Integer) where {A<:CUDA.CuArray} = ny > 0 && nx > 0
@inline Proper.ka_geometry_enabled(::Type{A}, ny::Integer, nx::Integer) where {A<:CUDA.CuArray} = ny > 0 && nx > 0
@inline Proper.ka_sampling_enabled(::Type{A}, ny::Integer, nx::Integer) where {A<:CUDA.CuArray} = ny > 0 && nx > 0
@inline Proper.ka_cubic_grid_enabled(::Type{A}, ny::Integer, nx::Integer) where {A<:CUDA.CuArray} = ny > 0 && nx > 0
@inline Proper.ka_rotate_enabled(::Type{A}, ny::Integer, nx::Integer) where {A<:CUDA.CuArray} = ny > 0 && nx > 0

end
