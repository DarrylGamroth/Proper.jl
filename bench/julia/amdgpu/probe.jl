try
    using AMDGPU
    AMDGPU.functional() || error("AMDGPU.functional() returned false")
    AMDGPU.functional(:rocfft) || error("AMDGPU.functional(:rocfft) returned false")
    AMDGPU.allowscalar(false)
catch err
    println(stderr, sprint(showerror, err))
    exit(1)
end
