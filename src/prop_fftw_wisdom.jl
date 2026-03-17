function _wisdom_path(gridsize::Integer, nthreads::Integer)
    return joinpath(homedir(), ".proper_$(Int(gridsize))pix$(Int(nthreads))threads_wisdomfile")
end

"""
    prop_fftw_wisdom(gridsize; nthreads=Threads.nthreads())
    prop_fftw_wisdom(path)

Generate or export FFTW wisdom for PROPER FFT workloads.

# Arguments
- `gridsize`: FFT grid size used to create measured wisdom
- `nthreads`: thread count to bake into the generated wisdom
- `path`: explicit output path for the currently loaded FFTW wisdom

# Returns
- The path written to disk.
"""
function prop_fftw_wisdom(gridsize::Integer; nthreads::Integer=Threads.nthreads())
    n = Int(gridsize)
    t = max(1, Int(nthreads))

    old_threads = FFTW.get_num_threads()
    FFTW.set_num_threads(t)
    try
        data = ones(ComplexF64, n, n)
        p_f = plan_fft!(data; flags=FFTW.MEASURE)
        p_b = plan_ifft!(data; flags=FFTW.MEASURE)
        p_f * data
        p_b * data

        path = _wisdom_path(n, t)
        FFTW.export_wisdom(path)
        return path
    finally
        FFTW.set_num_threads(old_threads)
    end
end

# Backward-compatible method for the older compatibility API that accepted a path.
function prop_fftw_wisdom(path::AbstractString)
    FFTW.export_wisdom(path)
    return path
end
