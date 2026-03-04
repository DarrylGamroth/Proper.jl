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

interp_style(::Type{<:AbstractArray}) = GenericInterpStyle()

abstract type RNGStyle end
struct GenericRNGStyle <: RNGStyle end

rng_style(::Type) = GenericRNGStyle()
