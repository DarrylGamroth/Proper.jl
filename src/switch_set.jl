"""Compatibility helper for uppercase/lowercase boolean keyword switches."""
function switch_set(key::Union{String,Symbol}; kwargs...)
    s = Symbol(key)
    l = Symbol(lowercase(String(s)))
    if haskey(kwargs, s)
        return compat_bool(kwargs[s])
    elseif haskey(kwargs, l)
        return compat_bool(kwargs[l])
    end
    return false
end
