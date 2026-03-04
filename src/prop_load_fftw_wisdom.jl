"""Load FFTW wisdom placeholder; currently no-op."""
function prop_load_fftw_wisdom(path::AbstractString)
    return isfile(path)
end
