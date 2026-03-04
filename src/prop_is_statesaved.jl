"""Check whether a serialized wavefront state file exists."""
prop_is_statesaved(path::AbstractString) = isfile(path)
