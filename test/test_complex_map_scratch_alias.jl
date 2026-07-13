using Test
using FFTW

function _asymmetric_complex_map(::Type{T}, n::Integer) where {T<:AbstractFloat}
    return [
        complex(T(1) + T(0.031) * T(i) + T(0.047) * T(j), T(0.019) * T(i) - T(0.023) * T(j))
        for i in 1:n, j in 1:n
    ]
end

function _asymmetric_complex_field(::Type{T}, n::Integer) where {T<:AbstractFloat}
    return [
        complex(T(0.13) * T(i) - T(0.071) * T(j), T(0.029) * T(i * j) + T(0.11))
        for i in 1:n, j in 1:n
    ]
end

function _propagated_complex_map_alias_regression!(
    device_array,
    sync!,
    ::Type{T};
    n::Int=9,
    atol::T=zero(T),
    rtol::T=zero(T),
) where {T<:AbstractFloat}
    initial_cpu = _asymmetric_complex_field(T, n)
    map_cpu = _asymmetric_complex_map(T, n)
    map_backend = device_array(map_cpu)

    wf_mul = Proper.WaveFront(device_array(initial_cpu), T(500e-9), T(1e-4), zero(T), one(T))
    ctx_mul = RunContext(wf_mul)
    prop_ptp(wf_mul, T(0.01), ctx_mul)
    sync!()
    @test wf_mul.field === wf_mul.workspace.fft.scratch

    propagated_mul = copy(wf_mul.field)
    sync!()
    shifted = Proper.shift_center_for_wavefront!(wf_mul, map_backend; inverse=true)
    sync!()
    @test shifted !== wf_mul.field
    @test isapprox(Array(wf_mul.field), Array(propagated_mul); atol=atol, rtol=rtol)
    @test isapprox(Array(shifted), FFTW.ifftshift(map_cpu); atol=atol, rtol=rtol)

    expected_mul = Array(propagated_mul) .* FFTW.ifftshift(map_cpu)
    prop_multiply(wf_mul, map_backend)
    sync!()
    @test isapprox(Array(wf_mul.field), expected_mul; atol=atol, rtol=rtol)

    wf_div = Proper.WaveFront(device_array(initial_cpu), T(500e-9), T(1e-4), zero(T), one(T))
    ctx_div = RunContext(wf_div)
    prop_ptp(wf_div, T(0.01), ctx_div)
    sync!()
    @test wf_div.field === wf_div.workspace.fft.scratch

    propagated_div = copy(wf_div.field)
    sync!()
    expected_div = Array(propagated_div) ./ FFTW.ifftshift(map_cpu)
    prop_divide(wf_div, map_backend)
    sync!()
    @test isapprox(Array(wf_div.field), expected_div; atol=atol, rtol=rtol)
    return nothing
end

@testset "Complex centered maps preserve propagated FFT-scratch fields" begin
    _propagated_complex_map_alias_regression!(Array, () -> nothing, Float64)

    @testset "Optional CUDA" begin
        cuda_ready = false
        try
            @eval using CUDA
            cuda_ready = CUDA.functional()
        catch
            cuda_ready = false
        end

        if cuda_ready
            CUDA.allowscalar(false)
            _propagated_complex_map_alias_regression!(
                CUDA.CuArray,
                CUDA.synchronize,
                Float32;
                n=16,
                atol=2f-5,
                rtol=2f-5,
            )
        else
            @test true
        end
    end

    @testset "Optional AMDGPU" begin
        amdgpu_ready = false
        try
            @eval using AMDGPU
            amdgpu_ready = AMDGPU.functional() && AMDGPU.functional(:rocfft)
        catch
            amdgpu_ready = false
        end

        if amdgpu_ready
            AMDGPU.allowscalar(false)
            _propagated_complex_map_alias_regression!(
                AMDGPU.ROCArray,
                AMDGPU.synchronize,
                Float32;
                n=16,
                atol=4f-4,
                rtol=1f-3,
            )
        else
            @test true
        end
    end
end
