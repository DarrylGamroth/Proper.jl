using AMDGPU
using Proper

struct PreparedAMDGPUOpticalPath{W,A,P,O,C,T}
    wavefront::W
    entrance_field::A
    amplitude::P
    opd_m::P
    output::O
    context::C
    dz_m::T
end

function propagate_amdgpu_optical_path!(path::PreparedAMDGPUOpticalPath)
    # Device-kernel assignment is HIP-stream-capture compatible.
    path.wavefront.field .= path.entrance_field
    path.wavefront.reference_surface = Proper.PLANAR
    prop_ptp(path.wavefront, path.dz_m, path.context)
    prop_multiply(path.wavefront, path.amplitude)
    prop_add_phase(path.wavefront, path.opd_m)
    prop_end!(path.output, path.wavefront; noabs=true)
    return path.output
end

function amdgpu_graph_prepared_demo(; gridsize::Integer=256, warmup::Integer=10)
    AMDGPU.functional() || error("AMDGPU is not functional on this host")
    AMDGPU.functional(:rocfft) || error("AMDGPU ROCFFT is not functional on this host")
    AMDGPU.allowscalar(false)

    T = Float32
    wavelength_m = T(1.55e-6)
    entrance_field = AMDGPU.fill(Complex{T}(1), gridsize, gridsize)
    wavefront = Proper.WaveFront(
        similar(entrance_field),
        wavelength_m,
        T(15e-6),
        zero(T),
        one(T),
    )
    path = PreparedAMDGPUOpticalPath(
        wavefront,
        entrance_field,
        AMDGPU.ones(T, gridsize, gridsize),
        AMDGPU.zeros(T, gridsize, gridsize),
        AMDGPU.zeros(Complex{T}, gridsize, gridsize),
        RunContext(typeof(entrance_field)),
        T(0.01),
    )

    graph = prepare_hip_graph(propagate_amdgpu_optical_path!, path; warmup)

    # A new frame changes existing device-buffer contents, not their identities.
    fill!(path.opd_m, T(10e-9))
    launch_hip_graph!(graph)
    AMDGPU.synchronize()
    return Array(path.output), graph
end

if abspath(PROGRAM_FILE) == @__FILE__
    output, _ = amdgpu_graph_prepared_demo()
    println("output checksum: ", sum(abs2, output))
end
