"""Load FFTW wisdom by grid/thread naming convention or explicit file path."""
function prop_load_fftw_wisdom(gridsize::Integer, nthreads::Integer=Threads.nthreads())
    path = joinpath(homedir(), ".proper_$(Int(gridsize))pix$(Int(nthreads))threads_wisdomfile")
    return prop_load_fftw_wisdom(path)
end

function prop_load_fftw_wisdom(path::AbstractString)
    if isfile(path)
        FFTW.forget_wisdom()
        FFTW.import_wisdom(path)
        return true
    end
    return false
end
