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
