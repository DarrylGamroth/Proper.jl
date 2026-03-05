abstract type BackendStyle end
struct CPUBackend <: BackendStyle end
struct UnknownBackend <: BackendStyle end

backend_style(::Type{<:AbstractArray}) = CPUBackend()

abstract type FFTStyle end
struct FFTWStyle <: FFTStyle end
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

@inline function ka_mask_enabled(::Type{A}, ny::Integer, nx::Integer) where {A<:AbstractArray}
    return (shift_kernel_style(A) isa ShiftKAStyle) && (ny * nx >= KA_MASK_MIN_ELEMS)
end

@inline function ka_end_enabled(::Type{A}, ny::Integer, nx::Integer) where {A<:AbstractArray}
    return (shift_kernel_style(A) isa ShiftKAStyle) && (ny * nx >= KA_END_MIN_ELEMS)
end
