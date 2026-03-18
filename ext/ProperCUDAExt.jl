module ProperCUDAExt

using Proper
using CUDA
using CUDA.CUFFT

const CUFFTFwdPlan2D{T} = CUDA.CUFFT.CuFFTPlan{Complex{T},Complex{T},CUDA.CUFFT.CUFFT_FORWARD,true,2,2,Nothing}
const CUFFTBwdPlan2D{T} = CUDA.CUFFT.CuFFTPlan{Complex{T},Complex{T},CUDA.CUFFT.CUFFT_INVERSE,true,2,2,Nothing}
const MaybeCUFFTFwdPlan2D{T} = Union{Nothing,CUFFTFwdPlan2D{T}}
const MaybeCUFFTBwdPlan2D{T} = Union{Nothing,CUFFTBwdPlan2D{T}}

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
Proper.fft_plan_exec_style(::Type{<:CUDA.CuArray}) = Proper.FFTPlanAvailableStyle()

function Proper.FFTWorkspace(::Type{<:CUDA.CuArray}, ::Type{T}=Float64) where {T<:AbstractFloat}
    rho2 = CUDA.CuMatrix{T}(undef, 0, 0)
    scratch = CUDA.CuMatrix{Complex{T}}(undef, 0, 0)
    real_scratch = CUDA.CuMatrix{T}(undef, 0, 0)
    return Proper.FFTWorkspace{T,typeof(rho2),typeof(scratch),typeof(real_scratch),MaybeCUFFTFwdPlan2D{T},MaybeCUFFTBwdPlan2D{T}}(
        rho2,
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

function Proper._ensure_fft_plans!(
    ::Proper.FFTPlanAvailableStyle,
    ws::Proper.FFTWorkspace{T,<:CUDA.CuMatrix{T},<:CUDA.CuMatrix{Complex{T}},<:CUDA.CuMatrix{T},MaybeCUFFTFwdPlan2D{T},MaybeCUFFTBwdPlan2D{T}},
    ny::Integer,
    nx::Integer,
    planning::Proper.FFTPlanningStyle,
) where {T<:AbstractFloat}
    _ = planning
    Proper.ensure_fft_scratch!(ws, ny, nx)
    if !(ws.plans_valid && ws.forward_plan !== nothing && ws.backward_plan !== nothing)
        ws.forward_plan = CUDA.CUFFT.plan_fft!(ws.scratch, (1, 2))
        ws.backward_plan = CUDA.CUFFT.plan_bfft!(ws.scratch, (1, 2))
        ws.plans_valid = true
    end
    return ws.forward_plan::CUFFTFwdPlan2D{T}, ws.backward_plan::CUFFTBwdPlan2D{T}
end

@inline Proper.ka_mask_enabled(::Type{A}, ny::Integer, nx::Integer) where {A<:CUDA.CuArray} = ny > 0 && nx > 0
@inline Proper.ka_end_enabled(::Type{A}, ny::Integer, nx::Integer) where {A<:CUDA.CuArray} = ny > 0 && nx > 0
@inline Proper.ka_shift_enabled(::Type{A}, ny::Integer, nx::Integer) where {A<:CUDA.CuArray} = ny > 0 && nx > 0
@inline Proper.ka_geometry_enabled(::Type{A}, ny::Integer, nx::Integer) where {A<:CUDA.CuArray} = ny > 0 && nx > 0
@inline Proper.ka_sampling_enabled(::Type{A}, ny::Integer, nx::Integer) where {A<:CUDA.CuArray} = ny > 0 && nx > 0
@inline Proper.ka_cubic_grid_enabled(::Type{A}, ny::Integer, nx::Integer) where {A<:CUDA.CuArray} = ny > 0 && nx > 0
@inline Proper.ka_rotate_enabled(::Type{A}, ny::Integer, nx::Integer) where {A<:CUDA.CuArray} = ny > 0 && nx > 0

end
