using Test, Proper, FFTW, LinearAlgebra, Random, KernelAbstractions

@kernel function expand_prepared_phase!(out,phase)
    y,x=@index(Global,NTuple)
    ny,nx=size(out)
    if y<=ny && x<=nx
        @inbounds out[y,x]=Proper._propagation_phase(phase,y,x,ny,nx)
    end
end

function check_phase_factors(p,wf,d,ctx,kind,synchronize)
    n=size(wf.field,1)
    factor=kind===:ptp ? p.transfer : kind===:wts ? p.before : p.after
    isempty(factor) && return nothing
    expected=similar(wf.field)
    fill!(expected,one(eltype(expected)))
    if kind===:ptp
        Proper._prepare_frequency_phase!(Proper.backend_style(typeof(expected)),expected,
            -pi*wf.wavelength_m*d,wf.sampling_m,ctx,p.scale)
    else
        effective=kind===:stw && iszero(d) ? wf.z_w0_m-wf.z_m : d
        dx=kind===:stw ? wf.wavelength_m*abs(effective)/(n*wf.sampling_m) : wf.sampling_m
        phase_wf=WaveFront(expected,wf.wavelength_m,dx,zero(dx),wf.beam_diameter_m)
        prop_qphase(phase_wf,effective,ctx)
    end
    expanded=similar(expected)
    expand_prepared_phase!(KernelAbstractions.get_backend(expanded),(16,16))(
        expanded,factor;ndrange=size(expanded))
    synchronize()
    @test isequal(Array(expanded),Array(expected))
    stored=Proper._propagation_phase_storage(factor)
    stored===nothing || @test size(stored)==(fld(n,2)+1,fld(n,2)+1)
    return nothing
end

function check_prepared_propagation(array, synchronize; capture=nothing, replay=nothing)
    @testset "Prepared field propagation $array" begin
        for T in (Float32,Float64), n in (1,2,7,8,512), kind in (:ptp,:wts,:stw), sign in (-1,0,1)
            rng = MersenneTwister(19)
            d, dx, lambda = T(sign*.013), T(4e-5), T(1.6e-6)
            input = array(randn(rng,Complex{T},n,n))
            edge = array(complex.(T(.7) .+ T(.2).*rand(rng,T,n,n)))
            leading = array(complex.(T(.8) .+ T(.1).*rand(rng,T,n,n)))
            wf = WaveFront(copy(input),lambda,dx,T(.01),one(T))
            wf.z_w0_m = T(.003)
            ctx = RunContext(wf)
            p = Proper.prepare_propagation(wf,d,ctx;kind)
            check_phase_factors(p,wf,d,ctx,kind,synchronize)
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
        # STW at the waist resolves its zero sentinel to zero: it still
        # transforms/scales the field, and its quadratic factor is identity.
        for n in (8,512)
            input=array(randn(MersenneTwister(23),ComplexF32,n,n))
            wf=WaveFront(copy(input),1.6f-6,4f-5,.01f0,1f0)
            wf.z_w0_m=wf.z_m
            wf.reference_surface=Proper.SPHERICAL
            context=RunContext(wf)
            prepared=Proper.prepare_propagation(wf,0f0,context;kind=:stw)
            prop_stw(wf,0f0,context)
            synchronize();expected=Array(wf.field)
            out=similar(input)
            operation()=Proper.propagate!(out,prepared,input)
            operation();synchronize()
            @test Array(out)≈expected rtol=128eps(Float32) atol=128eps(Float32)
            if capture!==nothing
                graph=capture(operation);replay(graph);synchronize()
                @test Array(out)≈expected rtol=128eps(Float32) atol=128eps(Float32)
            end
        end
        # Longer raw-order groups preserve every scalar gain and sampled
        # multiplier, including a mask followed by the next edge multiplier.
        for T in (Float32,Float64), n in (7,8,512), count in (3,4,5,8), centered in (false,true)
            rng=MersenneTwister(41)
            input=array(randn(rng,Complex{T},n,n))
            wf=WaveFront(copy(input),T(1.6e-6),T(4e-5),zero(T),one(T))
            context=RunContext(wf)
            kinds=(:stw,:wts,:ptp,:stw,:wts)
            transforms=ntuple(i->Proper.prepare_propagation(wf,T(isodd(i) ? .013 : -.025),context;kind=kinds[mod1(i,5)]),count)
            edge=array(Complex{T}.(0.9 .+0.02im .*randn(rng,n,n)))
            mask=array(Complex{T}.(0.8 .+0.03im .*randn(rng,n,n)))
            inputs=ntuple(i->i==3 ? (mask,edge) : (isodd(i) ? nothing : edge),count)
            gains=ntuple(i->Complex{T}(cis(T(i*.17))),count)
            out,expected=similar(input),similar(input)
            operation()=Proper.propagate!(out,transforms,input;centered,
                input_multipliers=inputs,output_multipliers=gains)
            captured=capture===nothing ? nothing : capture(operation)
            for state in 1:2
                input .*=Complex{T}(0.91+0.07im)
                copyto!(expected,input)
                for i in 1:count
                    multiplier=inputs[i]
                    if multiplier isa Tuple
                        for factor in multiplier
                            expected .*=factor
                        end
                    elseif multiplier!==nothing
                        expected .*=multiplier
                    end
                    Proper.propagate!(expected,transforms[i],expected;centered,output_multiplier=gains[i])
                end
                captured===nothing ? operation() : replay(captured)
                synchronize()
                @test Array(out)≈Array(expected) rtol=128eps(T) atol=128eps(T)
            end
            @test_throws DimensionMismatch Proper.propagate!(out,transforms,input;input_multipliers=(nothing,))
            @test_throws ArgumentError Proper.propagate!(out,transforms,input;
                input_multipliers=ntuple(_->(edge,first(transforms).field),count))
            if T===Float32 && n==8 && count==4 && centered
                saved=Array(out)
                @test_throws TypeError Proper.propagate!(out,transforms,input;input_multipliers=collect(inputs))
                @test_throws TypeError Proper.propagate!(out,transforms,input;output_multipliers=collect(gains))
                synchronize()
                @test Array(out)==saved
            end
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
