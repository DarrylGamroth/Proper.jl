abstract type BackendStyle end
struct CPUBackend <: BackendStyle end
struct CUDABackend <: BackendStyle end
struct UnknownBackend <: BackendStyle end

backend_style(::Type{<:AbstractArray}) = CPUBackend()

abstract type ArrayLayoutStyle end
struct GenericLayout <: ArrayLayoutStyle end
struct StridedLayout <: ArrayLayoutStyle end

array_layout_style(::Type{<:AbstractArray}) = GenericLayout()
array_layout_style(::Type{<:StridedArray}) = StridedLayout()

abstract type FFTStyle end
struct FFTWStyle <: FFTStyle end
struct CUFFTStyle <: FFTStyle end
struct GenericFFTStyle <: FFTStyle end

fft_style(::Type{<:AbstractArray}) = FFTWStyle()

abstract type InterpStyle end
struct GenericInterpStyle <: InterpStyle end
struct CubicInterpStyle <: InterpStyle end

interp_style(::Type{<:AbstractArray}) = GenericInterpStyle()
interp_style(::Type{<:StridedMatrix}) = CubicInterpStyle()

abstract type RNGStyle end
struct GenericRNGStyle <: RNGStyle end

rng_style(::Type) = GenericRNGStyle()

abstract type ShiftKernelStyle end
struct ShiftLoopStyle <: ShiftKernelStyle end
struct ShiftKAStyle <: ShiftKernelStyle end

shift_kernel_style(::Type{<:AbstractArray}) = ShiftLoopStyle()
shift_kernel_style(::Type{<:StridedMatrix}) = ShiftKAStyle()

const KA_MASK_MIN_ELEMS = typemax(Int)
const KA_END_MIN_ELEMS = 200_000

abstract type InterpKernelStyle end
struct InterpLoopStyle <: InterpKernelStyle end
struct InterpKAStyle <: InterpKernelStyle end

interp_kernel_style(::Type{<:AbstractArray}) = InterpLoopStyle()
interp_kernel_style(::Type{<:StridedMatrix}) = InterpKAStyle()

const KA_CUBIC_GRID_MIN_ELEMS = typemax(Int)
const KA_ROTATE_MIN_ELEMS = typemax(Int)

@inline ka_enabled(::ShiftKAStyle, ny::Integer, nx::Integer, min_elems::Int) = (ny * nx) >= min_elems
@inline ka_enabled(::InterpKAStyle, ny::Integer, nx::Integer, min_elems::Int) = (ny * nx) >= min_elems
@inline ka_enabled(::ShiftKernelStyle, ny::Integer, nx::Integer, min_elems::Int) = false
@inline ka_enabled(::InterpKernelStyle, ny::Integer, nx::Integer, min_elems::Int) = false

@inline function ka_mask_enabled(::Type{A}, ny::Integer, nx::Integer) where {A<:AbstractArray}
    return ka_enabled(shift_kernel_style(A), ny, nx, KA_MASK_MIN_ELEMS)
end

@inline function ka_end_enabled(::Type{A}, ny::Integer, nx::Integer) where {A<:AbstractArray}
    return ka_enabled(shift_kernel_style(A), ny, nx, KA_END_MIN_ELEMS)
end

abstract type GeometryKernelStyle end
struct GeometryLoopStyle <: GeometryKernelStyle end
struct GeometryKAStyle <: GeometryKernelStyle end

geometry_kernel_style(::Type{<:AbstractArray}) = GeometryLoopStyle()

const KA_GEOMETRY_MIN_ELEMS = typemax(Int)

@inline ka_enabled(::GeometryKAStyle, ny::Integer, nx::Integer, min_elems::Int) = (ny * nx) >= min_elems
@inline ka_enabled(::GeometryKernelStyle, ny::Integer, nx::Integer, min_elems::Int) = false

@inline function ka_geometry_enabled(::Type{A}, ny::Integer, nx::Integer) where {A<:AbstractArray}
    return ka_enabled(geometry_kernel_style(A), ny, nx, KA_GEOMETRY_MIN_ELEMS)
end

abstract type GeometryExecStyle end
struct GeometryLoopExecStyle <: GeometryExecStyle end
struct GeometryKAExecStyle <: GeometryExecStyle end

@inline geometry_exec_style(::Val{true}) = GeometryKAExecStyle()
@inline geometry_exec_style(::Val{false}) = GeometryLoopExecStyle()

@inline function geometry_exec_style(::Type{A}, ny::Integer, nx::Integer) where {A<:AbstractArray}
    return geometry_exec_style(Val(ka_geometry_enabled(A, ny, nx)))
end

abstract type SamplingKernelStyle end
struct SamplingLoopStyle <: SamplingKernelStyle end
struct SamplingKAStyle <: SamplingKernelStyle end

sampling_kernel_style(::Type{<:AbstractArray}) = SamplingLoopStyle()

const KA_SAMPLING_MIN_ELEMS = typemax(Int)

@inline ka_enabled(::SamplingKAStyle, ny::Integer, nx::Integer, min_elems::Int) = (ny * nx) >= min_elems
@inline ka_enabled(::SamplingKernelStyle, ny::Integer, nx::Integer, min_elems::Int) = false

@inline function ka_sampling_enabled(::Type{A}, ny::Integer, nx::Integer) where {A<:AbstractArray}
    return ka_enabled(sampling_kernel_style(A), ny, nx, KA_SAMPLING_MIN_ELEMS)
end

abstract type SamplingExecStyle end
struct SamplingLoopExecStyle <: SamplingExecStyle end
struct SamplingKAExecStyle <: SamplingExecStyle end

@inline sampling_exec_style(::Val{true}) = SamplingKAExecStyle()
@inline sampling_exec_style(::Val{false}) = SamplingLoopExecStyle()

@inline function sampling_exec_style(::Type{A}, ny::Integer, nx::Integer) where {A<:AbstractArray}
    return sampling_exec_style(Val(ka_sampling_enabled(A, ny, nx)))
end

@inline function ka_cubic_grid_enabled(::Type{A}, ny::Integer, nx::Integer) where {A<:AbstractArray}
    return ka_enabled(interp_kernel_style(A), ny, nx, KA_CUBIC_GRID_MIN_ELEMS)
end

@inline function ka_rotate_enabled(::Type{A}, ny::Integer, nx::Integer) where {A<:AbstractArray}
    return ka_enabled(interp_kernel_style(A), ny, nx, KA_ROTATE_MIN_ELEMS)
end

abstract type AxisFillExecStyle end
struct AxisFillLoopExecStyle <: AxisFillExecStyle end
struct AxisFillKAExecStyle <: AxisFillExecStyle end

@inline axis_fill_exec_style(::CPUBackend) = AxisFillLoopExecStyle()
@inline axis_fill_exec_style(::CUDABackend) = AxisFillKAExecStyle()
@inline axis_fill_exec_style(::BackendStyle) = AxisFillLoopExecStyle()

@inline same_backend_style(::B, ::B) where {B<:BackendStyle} = true
@inline same_backend_style(::BackendStyle, ::BackendStyle) = false
@inline same_backend_style(::Type{A}, ::Type{B}) where {A<:AbstractArray,B<:AbstractArray} =
    same_backend_style(backend_style(A), backend_style(B))

@inline function backend_adapt(template::AbstractArray, src::AbstractArray)
    if same_backend_style(typeof(template), typeof(src))
        return src
    end
    out = similar(template, eltype(src), size(src)...)
    copyto!(out, src)
    return out
end
