using Test, Proper, FFTW, LinearAlgebra, Random

function check_prepared_propagation(array, synchronize; capture=nothing, replay=nothing)
    @testset "Prepared field propagation $array" begin
        for T in (Float32,Float64), n in (7,8,512), kind in (:ptp,:wts,:stw), sign in (-1,0,1)
            rng = MersenneTwister(19)
            d, dx, lambda = T(sign*.013), T(4e-5), T(1.6e-6)
            input = array(randn(rng,Complex{T},n,n))
            edge = array(complex.(T(.7) .+ T(.2).*rand(rng,T,n,n)))
            leading = array(complex.(T(.8) .+ T(.1).*rand(rng,T,n,n)))
            wf = WaveFront(copy(input),lambda,dx,T(.01),one(T))
            wf.z_w0_m = T(.003)
            ctx = RunContext(wf)
            p = Proper.prepare_propagation(wf,d,ctx;kind)
            @test Array(wf.field) == Array(input)
            @test wf.sampling_m == dx && wf.z_m == T(.01)
            out = similar(input)
            for centered in (false,true)
                wf.field .= input
                wf.field .*= edge
                centered && (wf.field .= ifftshift(wf.field))
                wf.sampling_m = dx
                wf.z_m = T(.01)
                wf.reference_surface = kind===:stw ? Proper.SPHERICAL : Proper.PLANAR
                getproperty(Proper,Symbol(:prop_,kind))(wf,d,ctx)
                expected = Array((centered ? fftshift(wf.field) : wf.field) .* Complex{T}(-im))
                operation() = Proper.propagate!(out,p,input;centered,
                    input_multiplier=edge,output_multiplier=-im)
                operation();synchronize()
                @test Array(out) ≈ expected rtol=128eps(T) atol=128eps(T)
                if capture!==nothing && n==7
                    graph=capture(operation)
                    replay(graph);synchronize()
                    @test Array(out) ≈ expected rtol=128eps(T) atol=128eps(T)
                    fill!(input,0)
                    replay(graph);synchronize()
                    @test all(iszero,Array(out))
                    copyto!(input,array(randn(MersenneTwister(19),Complex{T},n,n)))
                end
            end
            # Two different propagations share exactly the prepared FFT storage.
            second_wf=WaveFront(similar(input),lambda,T(6e-5),zero(T),one(T))
            second=Proper.prepare_propagation(second_wf,T(-.02),ctx;kind=:wts)
            temporary=similar(input)
            Proper.propagate!(temporary,p,input;centered=true,input_multiplier=leading,output_multiplier=-im)
            Proper.propagate!(out,second,temporary;centered=true,input_multiplier=edge,output_multiplier=im)
            synchronize();expected=Array(out)
            operation_pair()=Proper.propagate!(out,(p,second),input;centered=true,
                input_multipliers=(leading,edge),output_multipliers=(-im,im))
            operation_pair();synchronize()
            @test Array(out) ≈ expected rtol=128eps(T) atol=128eps(T)
            if capture!==nothing && n==8
                graph=capture(operation_pair)
                replay(graph);synchronize()
                @test Array(out) ≈ expected rtol=128eps(T) atol=128eps(T)
            end
            @test_throws ArgumentError Proper.propagate!(out,p,input;input_multiplier=p.field)
            @test_throws ArgumentError Proper.propagate!(out,(p,second),input;input_multipliers=(edge,p.field))
            @test_throws ArgumentError Proper.propagate!(p.field,p,input)
            @test_throws DimensionMismatch Proper.propagate!(similar(input,n-1,n),p,input)
            @test_throws DimensionMismatch Proper.propagate!(out,p,input;input_multiplier=similar(edge,n-1,n))
        end
        wf=WaveFront(array(ones(ComplexF32,8,8)),1.6f-6,4f-5,0f0,1f0)
        @test_throws ArgumentError Proper.prepare_propagation(wf,.01f0,RunContext(Matrix{ComplexF64}))
        @test_throws ArgumentError Proper.prepare_propagation(wf,0f0,RunContext(wf);kind=:unknown)
        @test_throws ArgumentError Proper.prepare_propagation(wf,Inf32,RunContext(wf))
        @test_throws ArgumentError Proper.prepare_propagation(wf,.01f0,RunContext(wf;carrier_phase=TrackCarrierPhase()))
        if array === Array
            context=RunContext(wf;fft_planning=FFTMeasureStyle())
            wf.field=Proper.ensure_fft_scratch!(Proper.fft_workspace(context),8,8)
            wf.field .= ComplexF32[sin(i)+im*cos(j) for i in 1:8,j in 1:8]
            retained=copy(wf.field)
            Proper.prepare_propagation(wf,.01f0,context)
            @test wf.field == retained
        end
    end
end

if "--cuda" in ARGS
    using CUDA
    CUDA.allowscalar(false)
    check_prepared_propagation(CuArray,CUDA.synchronize;
        capture=operation->Proper.prepare_cuda_graph(operation;warmup=1),replay=Proper.launch_cuda_graph!)
else
    check_prepared_propagation(Array,()->nothing)
end
