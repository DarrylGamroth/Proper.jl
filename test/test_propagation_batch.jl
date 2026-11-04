# The original paired transfer store is part of the numerical contract,
# including conversion before scalar scale, after-phase and gain.
function check_propagation_transfer_unpack(array,synchronize;capture=nothing,replay=nothing)
    @testset "Transfer and terminal unpack $array" begin
        for T in (Float32,Float64),n in (1,2,7,8,512)
            host=ntuple(2) do plane
                values=Complex{T}[Complex{T}(sin(T(.17y+.09x+plane)),cos(T(.11y-.13x-plane)))
                    for y in 1:n,x in 1:n]
                plane==1 && (values[1,1]=Complex{T}(-zero(T),zero(T)))
                values
            end
            input=map(array,host)
            work=similar(first(input),n,n,2)
            planes=ntuple(i->view(work,:,:,i),2)
            reference=map(similar,input);output=map(similar,input)
            backend=KernelAbstractions.get_backend(first(input))
            half=fld(n,2)+1
            factors=Complex{T}[Complex{T}(cos(T(.03y+.07x)),sin(T(.11y-.02x)))
                for y in 1:half,x in 1:half]
            dense=Complex{T}[factors[min(y,n-y+2),min(x,n-x+2)] for y in 1:n,x in 1:n]
            compact=array(factors)
            transfers=(compact,array(dense))
            # One promoted transfer case makes the intermediate conversion
            # observable rather than relying on same-precision temporaries.
            T===Float32 && n==7 && (transfers=(transfers...,array(ComplexF64.(dense))))
            empty=array(zeros(Complex{T},0,0))
            cases=((empty,T(.37),Complex{T}(.31,-.27)),
                (compact,T(-.73),Complex{T}(-.19,.43)))
            n<=2 && (cases=(cases...,(empty,-zero(T),Complex{T}(one(T),-zero(T)))))
            for transfer in transfers,(after,scale,gain) in cases
                oracle()=begin
                    for plane in 1:2
                        copyto!(planes[plane],input[plane])
                    end
                    Proper._propagation_batch_transfer_kernel!(backend,(16,16))(
                        work,transfer;ndrange=(n,n))
                    for plane in 1:2
                        Proper._propagation_unpack_kernel!(backend,(16,16))(
                            reference[plane],planes[plane],after,scale,gain,cld(n,2);ndrange=(n,n))
                    end
                end
                operation()=begin
                    for plane in 1:2
                        Proper._propagation_transfer_unpack_kernel!(backend,(16,16))(
                            output[plane],input[plane],transfer,after,scale,gain,cld(n,2);ndrange=(n,n))
                    end
                end
                oracle();operation();synchronize()
                @test all(isequal(Array(output[i]),Array(reference[i])) for i in 1:2)
                @test all(isequal(Array(input[i]),host[i]) for i in 1:2)
                if capture!==nothing && n==7
                    recorded=capture(operation)
                    saved=Array(first(output))
                    fill!(last(input),zero(Complex{T}))
                    oracle();replay(recorded);synchronize()
                    @test all(isequal(Array(output[i]),Array(reference[i])) for i in 1:2)
                    @test isequal(Array(first(output)),saved)
                    @test all(iszero,Array(last(output)))
                    @test isequal(Array(first(input)),first(host)) && all(iszero,Array(last(input)))
                    copyto!(last(input),array(last(host)))
                end
            end
        end
    end
end

# Private batch kernels must retain the scalar kernels' arithmetic and plane
# isolation. FFT batching itself is exercised by the shared AOS branch gate.
function check_propagation_batch_kernels(array,synchronize;capture=nothing,replay=nothing)
    check_propagation_transfer_unpack(array,synchronize;capture,replay)
    @testset "Independent field batch kernels $array" begin
        for T in (Float32,Float64),n in (1,2,7,8,512)
            input=ntuple(2) do plane
                array(Complex{T}[Complex{T}(sin(T(.17i+.09j+plane)),cos(T(.11i-.13j))) for i in 1:n,j in 1:n])
            end
            masks=ntuple(plane->array(Complex{T}[Complex{T}(.7+.01plane*i,.03j) for i in 1:n,j in 1:n]),2)
            edge=array(fill(T(.83),n,n))
            wf=WaveFront(copy(first(input)),T(1.6e-6),T(4e-5),zero(T),one(T));context=RunContext(wf)
            before=Proper.prepare_propagation(wf,T(.013),context;kind=:wts).before
            after=Proper.prepare_propagation(wf,T(-.009),context;kind=:stw).after
            transfer=Proper.prepare_propagation(wf,T(.004),context;kind=:ptp).transfer
            field=similar(first(input),n,n,2)
            expected=ntuple(_->similar(first(input)),2)
            outputs=ntuple(_->similar(first(input)),2)
            reference=ntuple(_->similar(first(input)),2)
            backend=KernelAbstractions.get_backend(field)
            scale=T(.37);gain=Complex{T}(.31,-.27)
            for plane in 1:2
                Proper._propagation_pack_kernel!(backend,(16,16))(expected[plane],input[plane],
                    (masks[plane],edge),before,fld(n,2);ndrange=(n,n))
                Proper._propagation_join_kernel!(backend,(16,16))(expected[plane],after,scale,gain,
                    (masks[plane],edge),before,fld(n,2);ndrange=(n,n))
                Proper._propagation_transfer_kernel!(backend,(16,16))(expected[plane],transfer;ndrange=(n,n))
                Proper._propagation_unpack_kernel!(backend,(16,16))(reference[plane],expected[plane],
                    after,scale,gain,cld(n,2);ndrange=(n,n))
            end
            operation()=begin
                Proper._propagation_batch_pack_kernel!(backend,(16,16))(
                    field,input,masks,edge,before,fld(n,2);ndrange=(n,n))
                Proper._propagation_batch_join_kernel!(backend,(16,16))(
                    field,after,scale,gain,masks,edge,before,fld(n,2);ndrange=(n,n))
                Proper._propagation_batch_transfer_kernel!(backend,(16,16))(
                    field,transfer;ndrange=(n,n))
                Proper._propagation_batch_unpack_kernel!(backend,(16,16))(
                    outputs,field,after,scale,gain,cld(n,2);ndrange=(n,n))
            end
            operation();synchronize()
            for plane in 1:2
                @test Array(outputs[plane])≈Array(reference[plane]) rtol=128eps(T) atol=128eps(T)
            end
            if capture!==nothing && n==7
                recorded=capture(operation);replay(recorded);synchronize()
                saved=Array(first(outputs))
                fill!(last(input),0)
                replay(recorded);synchronize()
                @test Array(first(outputs))==saved
                @test all(iszero,Array(last(outputs)))
            end
        end
    end
end
