const _BOOL_COMPAT_KEYS = Set((:NOABS, :NORM, :verbose, :noabs, :norm))

function normalize_kwargs(kwargs::Base.Iterators.Pairs)
    d = Dict{Symbol,Any}()
    for (k, v) in kwargs
        sk = Symbol(k)
        lk = Symbol(lowercase(String(sk)))
        if haskey(d, lk) && d[lk] != v
            throw(ArgumentError("Conflicting values for keyword $(sk)"))
        end
        d[lk] = v
    end
    return NamedTuple{Tuple(keys(d))}(Tuple(values(d)))
end

@inline compat_bool(v::Bool) = v
@inline compat_bool(v::Integer) = v != 0
compat_bool(v) = throw(ArgumentError("Expected Bool or 0/1 integer, got $(typeof(v))"))

@inline kw_resolve(primary, secondary, default) =
    primary === nothing ? (secondary === nothing ? default : secondary) : primary

@inline kw_resolve_bool(primary, secondary, default::Bool=false)::Bool =
    compat_bool(kw_resolve(primary, secondary, default))

@inline kw_resolve_float(primary, secondary, default::Real)::Float64 =
    float(kw_resolve(primary, secondary, default))

@inline function kw_resolve_float(primary, secondary, default::Nothing=nothing)
    v = kw_resolve(primary, secondary, default)
    return v === nothing ? nothing : float(v)
end

@inline function kw_resolve_string(primary, secondary, default::Union{Nothing,AbstractString}=nothing)
    v = kw_resolve(primary, secondary, default)
    return v === nothing ? nothing : String(v)
end

@inline kw_resolve_symbol(primary, secondary, default::Symbol)::Symbol =
    Symbol(lowercase(String(kw_resolve(primary, secondary, default))))

@inline function kw_lookup(kwargs::Base.Iterators.Pairs, key::Symbol, default=nothing)
    lk = Symbol(lowercase(String(key)))
    if haskey(kwargs, key)
        return kwargs[key]
    elseif haskey(kwargs, lk)
        return kwargs[lk]
    end
    return default
end

@inline kw_lookup_bool(kwargs::Base.Iterators.Pairs, key::Symbol, default::Bool=false)::Bool =
    compat_bool(kw_lookup(kwargs, key, default))

@inline function kw_lookup_float(kwargs::Base.Iterators.Pairs, key::Symbol, default::Real)
    return float(kw_lookup(kwargs, key, default))
end

@inline function kw_lookup_float(kwargs::Base.Iterators.Pairs, key::Symbol, default::Nothing=nothing)
    v = kw_lookup(kwargs, key, default)
    return v === nothing ? nothing : float(v)
end

@inline function kw_lookup_string(kwargs::Base.Iterators.Pairs, key::Symbol, default::Union{Nothing,AbstractString}=nothing)
    v = kw_lookup(kwargs, key, default)
    return v === nothing ? nothing : String(v)
end

@inline kw_lookup_present(kwargs::Base.Iterators.Pairs, key::Symbol)::Bool =
    haskey(kwargs, key) || haskey(kwargs, Symbol(lowercase(String(key))))
