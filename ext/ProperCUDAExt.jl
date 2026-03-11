module ProperCUDAExt

using Proper
using CUDA

Proper.backend_style(::Type{<:CUDA.CuArray}) = Proper.CUDABackend()
Proper.fft_style(::Type{<:CUDA.CuArray}) = Proper.CUFFTStyle()
Proper.interp_style(::Type{<:CUDA.CuArray}) = Proper.CubicInterpStyle()
Proper.shift_kernel_style(::Type{<:CUDA.CuArray}) = Proper.ShiftKAStyle()
Proper.interp_kernel_style(::Type{<:CUDA.CuArray}) = Proper.InterpKAStyle()

@inline Proper.ka_mask_enabled(::Type{A}, ny::Integer, nx::Integer) where {A<:CUDA.CuArray} = ny > 0 && nx > 0
@inline Proper.ka_end_enabled(::Type{A}, ny::Integer, nx::Integer) where {A<:CUDA.CuArray} = ny > 0 && nx > 0
@inline Proper.ka_cubic_grid_enabled(::Type{A}, ny::Integer, nx::Integer) where {A<:CUDA.CuArray} = ny > 0 && nx > 0
@inline Proper.ka_rotate_enabled(::Type{A}, ny::Integer, nx::Integer) where {A<:CUDA.CuArray} = ny > 0 && nx > 0

end
