using Test, Proper, Random, FFTW, KernelAbstractions

_phase_original(value,opd,coefficient) =
    ComplexF32(ComplexF64(value)*cis(coefficient*Float64(opd)))
_phase_component_scale(value) = max(abs(Float64(real(value))),abs(Float64(imag(value))))
_phase_component_error(actual,expected) =
    max(abs(Float64(real(actual))-Float64(real(expected))),
        abs(Float64(imag(actual))-Float64(imag(expected))))

@testset "Prepared Float32 phase scalar accuracy and fallback" begin
    rng=MersenneTwister(48191)
    lower=floatmin(Float32);upper=floatmax(Float32)/4
    fields=ComplexF32[complex(0f0,0f0),complex(-0f0,-0f0),1+0.3im,-0.7+0.9im,
        complex(lower,lower),complex(lower,nextfloat(0f0)),
        complex(prevfloat(lower),-prevfloat(lower)),complex(nextfloat(0f0),0f0),
        complex(upper,-upper/2),complex(nextfloat(upper),1f0),
        complex(floatmax(Float32),floatmax(Float32)),complex(Inf32,1f0),complex(NaN32,1f0)]
    for wavelength in Float32[0.4e-6,1.55e-6,3e-6]
        coefficient=2pi/Float64(wavelength)
        angles=Float64[-0.0,0.0,-1e-8,1e-8,-1e10,1e10]
        for center in Float32[-2pi,-pi,-pi/2,pi/2,pi,2pi,-1024,1024]
            append!(angles,Float64[prevfloat(center),center,nextfloat(center)])
        end
        append!(angles,2048 .* (rand(rng,32).-0.5))
        opds=Float32.(angles./coefficient)
        threshold=Float32(1024/coefficient)
        append!(opds,Float32[prevfloat(threshold),threshold,nextfloat(threshold),
            -prevfloat(threshold),-threshold,-nextfloat(threshold),floatmax(Float32),-floatmax(Float32)])
        high=Float32.(coefficient.*Float64.(opds))
        @test any(abs.(high).<=1024f0) && any(abs.(high).>1024f0)
        for value in fields
            actual=Proper._phase_value_f32.(Ref(value),opds,coefficient)
            expected=_phase_original.(Ref(value),opds,coefficient)
            magnitude=_phase_component_scale(value)
            field_fallback=!isfinite(magnitude) || (magnitude!=0 && !(Float64(lower)<=magnitude<=Float64(upper)))
            fallback=map(h -> field_fallback || !isfinite(h) || abs(h)>1024f0,high)
            @test all(isequal(actual[i],expected[i]) for i in eachindex(opds) if fallback[i])
            @test all(_phase_component_error(actual[i],expected[i])<=8eps(Float32)*magnitude
                for i in eachindex(opds) if !fallback[i])
        end
        # CPU cis(Inf) can throw; preserve the original result or exception.
        for opd in (NaN32,Inf32,-Inf32)
            result=try
                _phase_original(1f0+0.3f0im,opd,coefficient)
            catch error
                error
            end
            if result isa Exception
                @test_throws typeof(result) Proper._phase_value_f32(1f0+0.3f0im,opd,coefficient)
            else
                @test isequal(Proper._phase_value_f32(1f0+0.3f0im,opd,coefficient),result)
            end
        end
    end
end

_phase_test_capture(::Proper.CUDABackend,operation) = Proper.prepare_cuda_graph(operation;warmup=1)
_phase_test_capture(::Proper.AMDGPUBackend,operation) = Proper.prepare_hip_graph(operation;warmup=1)
_phase_test_replay(::Proper.CUDABackend,graph) = Proper.launch_cuda_graph!(graph)
_phase_test_replay(::Proper.AMDGPUBackend,graph) = Proper.launch_hip_graph!(graph)

function check_prepared_phase_boundaries(array;captured=false)
    lower=floatmin(Float32);upper=floatmax(Float32)/4
    fields=ComplexF32[complex(0f0,0f0),complex(-0f0,-0f0),1+0.3im,
        complex(prevfloat(lower),-prevfloat(lower)),complex(lower,lower),
        complex(nextfloat(lower),lower),complex(lower,nextfloat(0f0)),
        complex(nextfloat(0f0),-nextfloat(0f0)),complex(prevfloat(upper),1f0),
        complex(upper,-upper/2),complex(nextfloat(upper),1f0),
        complex(floatmax(Float32),floatmax(Float32)),complex(NaN32,1f0),complex(Inf32,-Inf32)]
    wavelength=1.55f-6;coefficient=2pi/Float64(wavelength)
    angles=Float64[-0.0,0.0,-1e-8,1e-8,-1200,1200,-1e10,1e10]
    for center in Float32[-pi/2,pi/2]
        append!(angles,Float64[prevfloat(center),center,nextfloat(center)])
    end
    threshold=Float32(1024/coefficient)
    opds=vcat(Float32.(angles./coefficient),
        Float32[prevfloat(threshold),threshold,nextfloat(threshold),
            -prevfloat(threshold),-threshold,-nextfloat(threshold)])
    # GPU trig returns nonfinite values for infinite OPD; compare with the
    # same backend's original graph. CPU exception parity is checked above.
    array===Array || append!(opds,Float32[NaN32,Inf32,-Inf32])
    n=max(length(fields),length(opds))
    input=ComplexF32[fields[mod1(y,length(fields))] for y in 1:n,x in 1:n]
    host_phase=Float32[opds[mod1(x,length(opds))] for y in 1:n,x in 1:n]
    phase=array(host_phase)
    wf=WaveFront(array(input),wavelength,1f-3,0f0,1f0)
    original=WaveFront(array(input),wavelength,1f-3,0f0,1f0)
    prepared=Proper._PreparedPhase(wavelength,coefficient)
    operation()=Proper._apply_phase!(wf,phase,prepared)
    reference()=prop_add_phase(original,phase)
    operation();reference()
    KernelAbstractions.synchronize(get_backend(wf.field))
    backend=Proper.backend_style(typeof(wf.field))
    graph=captured ? _phase_test_capture(backend,operation) : nothing
    reference_graph=captured ? _phase_test_capture(backend,reference) : nothing
    copyto!(wf.field,array(input));copyto!(original.field,array(input))
    graph===nothing ? operation() : _phase_test_replay(backend,graph)
    reference_graph===nothing ? reference() : _phase_test_replay(backend,reference_graph)
    KernelAbstractions.synchronize(get_backend(wf.field))
    actual=Array(wf.field);expected=Array(original.field);shifted=ifftshift(host_phase)
    @test wf.field!==original.field && wf.workspace!==original.workspace
    for index in eachindex(input)
        magnitude=_phase_component_scale(input[index])
        high=Float32(coefficient*Float64(shifted[index]))
        fallback=!isfinite(high) || abs(high)>1024f0 || !isfinite(magnitude) ||
            (magnitude!=0 && !(Float64(lower)<=magnitude<=Float64(upper)))
        if fallback || !isfinite(expected[index])
            @test isequal(actual[index],expected[index])
        else
            @test _phase_component_error(actual[index],expected[index])<=8eps(Float32)*magnitude
        end
    end
    return nothing
end

function check_prepared_phase(array;captured=false)
    @testset "Prepared phase maps $array" begin
        for wavelength in Float32[0.4e-6,1.55e-6,3e-6], n in (5,6)
            coefficient=2pi/Float64(wavelength)
            field=array(zeros(ComplexF32,n,n))
            phase=array(zeros(Float32,n,n))
            wf=WaveFront(field,wavelength,1f-3,0f0,1f0)
            array===Array && @test Proper._prepare_phase(wf)===nothing
            prepared=Proper._PreparedPhase(wavelength,coefficient)
            operation()=Proper._apply_phase!(wf,phase,prepared)
            operation()
            KernelAbstractions.synchronize(get_backend(field))
            retained=wf.workspace.fft.real_scratch
            backend=Proper.backend_style(typeof(field))
            graph=captured ? _phase_test_capture(backend,operation) : nothing
            for state in (1,2,1)
                input=ComplexF32[complex(sin(0.3y+state),cos(0.2x-state)) for y in 1:n,x in 1:n]
                opd=Float32[(state==1 ? 0.17y-0.13x : 1200sin(0.2y+0.3x))/coefficient
                    for y in 1:n,x in 1:n]
                copyto!(field,array(input));copyto!(phase,array(opd))
                graph===nothing ? operation() : _phase_test_replay(backend,graph)
                KernelAbstractions.synchronize(get_backend(field))
                expected=_phase_original.(input,ifftshift(opd),coefficient)
                actual=Array(field)
                @test all(_phase_component_error(actual[i],expected[i])<=
                    8eps(Float32)*_phase_component_scale(input[i]) for i in eachindex(input))
                @test wf.field===field && wf.workspace.fft.real_scratch===retained
                @test Array(phase)==opd
            end
            before=Array(field);scratch_before=Array(retained)
            bad=array(zeros(Float32,n,n+1))
            @test_throws ArgumentError Proper._apply_phase!(wf,bad,prepared)
            @test Array(field)==before && Array(retained)==scratch_before
            wf.wavelength_m=nextfloat(wavelength)
            @test_throws ArgumentError operation()
            @test Array(field)==before && Array(retained)==scratch_before
            wf.wavelength_m=wavelength
            # The public original remains the exact Nothing-dispatch oracle.
            reference=WaveFront(copy(field),wavelength,1f-3,0f0,1f0)
            prop_add_phase(reference,phase)
            Proper._apply_phase!(wf,phase,nothing)
            KernelAbstractions.synchronize(get_backend(field))
            @test isequal(Array(field),Array(reference.field))
        end
        check_prepared_phase_boundaries(array;captured)
        wf64=WaveFront(array(ones(ComplexF64,2,2)),1.55e-6,1e-3,0.0,1.0)
        @test Proper._prepare_phase(wf64)===nothing
    end
end

if "--cuda" in ARGS
    using CUDA
    CUDA.allowscalar(false)
    check_prepared_phase(CuArray;captured=true)
elseif "--amdgpu" in ARGS
    using AMDGPU
    AMDGPU.allowscalar(false)
    check_prepared_phase(ROCArray;captured=true)
else
    check_prepared_phase(Array)
end
