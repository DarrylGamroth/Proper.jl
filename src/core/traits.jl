abstract type BackendStyle end
struct CPUBackend <: BackendStyle end
struct CUDABackend <: BackendStyle end
struct UnknownBackend <: BackendStyle end

backend_style(::Type{<:AbstractArray}) = CPUBackend()

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

@inline function ka_mask_enabled(::Type{A}, ny::Integer, nx::Integer) where {A<:AbstractArray}
    return (shift_kernel_style(A) isa ShiftKAStyle) && (ny * nx >= KA_MASK_MIN_ELEMS)
end

@inline function ka_end_enabled(::Type{A}, ny::Integer, nx::Integer) where {A<:AbstractArray}
    return (shift_kernel_style(A) isa ShiftKAStyle) && (ny * nx >= KA_END_MIN_ELEMS)
end

@inline function ka_cubic_grid_enabled(::Type{A}, ny::Integer, nx::Integer) where {A<:AbstractArray}
    return (interp_kernel_style(A) isa InterpKAStyle) && (ny * nx >= KA_CUBIC_GRID_MIN_ELEMS)
end

@inline function ka_rotate_enabled(::Type{A}, ny::Integer, nx::Integer) where {A<:AbstractArray}
    return (interp_kernel_style(A) isa InterpKAStyle) && (ny * nx >= KA_ROTATE_MIN_ELEMS)
end

@inline function same_backend_style(::Type{A}, ::Type{B}) where {A<:AbstractArray,B<:AbstractArray}
    return typeof(backend_style(A)) === typeof(backend_style(B))
end

@inline function backend_adapt(template::AbstractArray, src::AbstractArray)
    if same_backend_style(typeof(template), typeof(src))
        return src
    end
    out = similar(template, eltype(src), size(src)...)
    copyto!(out, src)
    return out
end
