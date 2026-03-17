"""Set polygon edge antialiasing subsampling factor (odd integer)."""
const _ANTIALIAS_SUBSAMPLING = Ref{Int}(11)

@inline antialias_subsampling() = _ANTIALIAS_SUBSAMPLING[]

"""
    prop_set_antialiasing(level=11)

Set the polygon-edge antialiasing subsampling factor.

# Arguments
- `level`: odd integer subsampling factor used by geometric mask builders

# Returns
- The stored antialiasing factor.
"""
function prop_set_antialiasing(level::Integer=11)
    if level < 1 || iseven(level)
        throw(ArgumentError("PROP_SET_ANTIALIASING: subsampling factor must be odd-valued integer."))
    end
    _ANTIALIAS_SUBSAMPLING[] = Int(level)
    return _ANTIALIAS_SUBSAMPLING[]
end

prop_set_antialiasing(level::Real) = prop_set_antialiasing(round(Int, level))
