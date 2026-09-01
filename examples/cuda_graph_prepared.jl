using CUDA
using Proper

struct PreparedCUDAOpticalPath{W,A,P,O}
    wavefront::W
    entrance_field::A
    amplitude::P
    opd_m::P
    output::O
end

function propagate_optical_path!(path::PreparedCUDAOpticalPath)
    copyto!(path.wavefront.field, path.entrance_field)
    prop_multiply(path.wavefront, path.amplitude)
    prop_add_phase(path.wavefront, path.opd_m)
    prop_end!(path.output, path.wavefront; noabs=true)
    return path.output
end

function cuda_graph_prepared_demo(; gridsize::Integer=256, warmup::Integer=10)
    CUDA.functional() || error("CUDA is not functional on this host")
    CUDA.allowscalar(false)

    T = Float32
    wavelength_m = T(1.55e-6)
    entrance_field = CUDA.fill(Complex{T}(1), gridsize, gridsize)
    wavefront = Proper.WaveFront(
        similar(entrance_field),
        wavelength_m,
        T(15e-6),
        zero(T),
        one(T),
    )
    path = PreparedCUDAOpticalPath(
        wavefront,
        entrance_field,
        CUDA.ones(T, gridsize, gridsize),
        CUDA.zeros(T, gridsize, gridsize),
        CUDA.zeros(Complex{T}, gridsize, gridsize),
    )

    graph = prepare_cuda_graph(propagate_optical_path!, path; warmup)

    # A new frame changes existing device-buffer contents, not their identities.
    fill!(path.opd_m, T(10e-9))
    launch_cuda_graph!(graph)
    CUDA.synchronize()
    return Array(path.output), graph
end

if abspath(PROGRAM_FILE) == @__FILE__
    output, _ = cuda_graph_prepared_demo()
    println("output checksum: ", sum(abs2, output))
end
