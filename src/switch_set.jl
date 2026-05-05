"""Compatibility helper for uppercase/lowercase boolean keyword switches."""
function switch_set(key::Union{String,Symbol}; kwargs...)
    isempty(kwargs) && return false
    s = Symbol(key)
    if haskey(kwargs, s)
        return compat_bool(kwargs[s])
    end
    l = Symbol(lowercase(String(s)))
    if haskey(kwargs, l)
        return compat_bool(kwargs[l])
    end
    return false
end
