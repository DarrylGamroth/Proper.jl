# Private batch kernels must retain the scalar kernels' arithmetic and plane
# isolation. FFT batching itself is exercised by the shared AOS branch gate.
function check_propagation_batch_kernels(array,synchronize;capture=nothing,replay=nothing)
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
