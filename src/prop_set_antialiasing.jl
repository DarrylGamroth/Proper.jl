"""Set polygon edge antialiasing subsampling factor (odd integer)."""
const _ANTIALIAS_SUBSAMPLING = Ref{Int}(11)

@inline antialias_subsampling() = _ANTIALIAS_SUBSAMPLING[]

function prop_set_antialiasing(level::Integer=11)
    if level < 1 || iseven(level)
        throw(ArgumentError("PROP_SET_ANTIALIASING: subsampling factor must be odd-valued integer."))
    end
    _ANTIALIAS_SUBSAMPLING[] = Int(level)
    return _ANTIALIAS_SUBSAMPLING[]
end

prop_set_antialiasing(level::Real) = prop_set_antialiasing(round(Int, level))
