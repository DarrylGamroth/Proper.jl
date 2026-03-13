try
    using CUDA
    CUDA.functional() || error("CUDA.functional() returned false")
    CUDA.allowscalar(false)
catch err
    println(stderr, sprint(showerror, err))
    exit(1)
end
