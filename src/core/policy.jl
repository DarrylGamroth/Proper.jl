abstract type CompatPolicy end

struct Python334Policy <: CompatPolicy end
struct CorrectedPolicy <: CompatPolicy end

@inline resolve_compat_policy(::Val{:python334}) = Python334Policy()
@inline resolve_compat_policy(::Val{:corrected}) = CorrectedPolicy()

function resolve_compat_policy(mode::Symbol)
    if mode === :python334
        return Python334Policy()
    elseif mode === :corrected
        return CorrectedPolicy()
    end
    throw(ArgumentError("Unsupported compat_mode=$(mode). Supported: :python334, :corrected"))
end
