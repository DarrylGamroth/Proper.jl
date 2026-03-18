module ProperAMDGPUExt

using Proper
using AMDGPU
using FFTW

const ROCFFTFwdPlan2D{T} = AMDGPU.rocFFT.cROCFFTPlan{Complex{T},true,true,2,2,Nothing}
const ROCFFTBwdPlan2D{T} = AMDGPU.rocFFT.cROCFFTPlan{Complex{T},false,true,2,2,Nothing}
const MaybeROCFFTFwdPlan2D{T} = Union{Nothing,ROCFFTFwdPlan2D{T}}
const MaybeROCFFTBwdPlan2D{T} = Union{Nothing,ROCFFTBwdPlan2D{T}}

Proper.backend_style(::Type{<:AMDGPU.ROCArray}) = Proper.AMDGPUBackend()
Proper.fft_style(::Type{<:AMDGPU.ROCArray}) = Proper.ROCFFTStyle()
Proper.interp_style(::Type{<:AMDGPU.ROCArray}) = Proper.CubicInterpStyle()
Proper.shift_kernel_style(::Type{<:AMDGPU.ROCArray}) = Proper.ShiftKAStyle()
Proper.geometry_kernel_style(::Type{<:AMDGPU.ROCArray}) = Proper.GeometryKAStyle()
Proper.sampling_kernel_style(::Type{<:AMDGPU.ROCArray}) = Proper.SamplingKAStyle()
Proper.interp_kernel_style(::Type{<:AMDGPU.ROCArray}) = Proper.InterpKAStyle()
Proper.workspace_vector(::Type{<:AMDGPU.ROCArray}, ::Type{T}, n::Integer=0) where {T<:AbstractFloat} = AMDGPU.ROCVector{T}(undef, n)
Proper.workspace_matrix(::Type{<:AMDGPU.ROCArray}, ::Type{T}, ny::Integer=0, nx::Integer=0) where {T<:AbstractFloat} = AMDGPU.ROCMatrix{T}(undef, ny, nx)
Proper.workspace_complex_matrix(::Type{<:AMDGPU.ROCArray}, ::Type{T}, ny::Integer=0, nx::Integer=0) where {T<:AbstractFloat} = AMDGPU.ROCMatrix{Complex{T}}(undef, ny, nx)
Proper.fft_plan_exec_style(::Type{<:AMDGPU.ROCArray}) = Proper.FFTPlanAvailableStyle()

function Proper.FFTWorkspace(::Type{<:AMDGPU.ROCArray}, ::Type{T}=Float64) where {T<:AbstractFloat}
    rho2 = AMDGPU.ROCMatrix{T}(undef, 0, 0)
    scratch = AMDGPU.ROCMatrix{Complex{T}}(undef, 0, 0)
    real_scratch = AMDGPU.ROCMatrix{T}(undef, 0, 0)
    return Proper.FFTWorkspace{T,typeof(rho2),typeof(scratch),typeof(real_scratch),MaybeROCFFTFwdPlan2D{T},MaybeROCFFTBwdPlan2D{T}}(
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
    ws::Proper.FFTWorkspace{T,<:AMDGPU.ROCMatrix{T},<:AMDGPU.ROCMatrix{Complex{T}},<:AMDGPU.ROCMatrix{T},MaybeROCFFTFwdPlan2D{T},MaybeROCFFTBwdPlan2D{T}},
    ny::Integer,
    nx::Integer,
    planning::Proper.FFTPlanningStyle,
) where {T<:AbstractFloat}
    _ = planning
    Proper.ensure_fft_scratch!(ws, ny, nx)
    if !(ws.plans_valid && ws.forward_plan !== nothing && ws.backward_plan !== nothing)
        ws.forward_plan = plan_fft!(ws.scratch, (1, 2))
        ws.backward_plan = plan_bfft!(ws.scratch, (1, 2))
        ws.plans_valid = true
    end
    return ws.forward_plan::ROCFFTFwdPlan2D{T}, ws.backward_plan::ROCFFTBwdPlan2D{T}
end

@inline Proper.ka_mask_enabled(::Type{A}, ny::Integer, nx::Integer) where {A<:AMDGPU.ROCArray} = ny > 0 && nx > 0
@inline Proper.ka_end_enabled(::Type{A}, ny::Integer, nx::Integer) where {A<:AMDGPU.ROCArray} = ny > 0 && nx > 0
@inline Proper.ka_geometry_enabled(::Type{A}, ny::Integer, nx::Integer) where {A<:AMDGPU.ROCArray} = ny > 0 && nx > 0
@inline Proper.ka_sampling_enabled(::Type{A}, ny::Integer, nx::Integer) where {A<:AMDGPU.ROCArray} = ny > 0 && nx > 0
@inline Proper.ka_cubic_grid_enabled(::Type{A}, ny::Integer, nx::Integer) where {A<:AMDGPU.ROCArray} = ny > 0 && nx > 0
@inline Proper.ka_rotate_enabled(::Type{A}, ny::Integer, nx::Integer) where {A<:AMDGPU.ROCArray} = ny > 0 && nx > 0

end
