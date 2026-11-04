function check_inplace_propagation(array,synchronize;capture=nothing,replay=nothing)
    @testset "Raw in-place prepared propagation $array" begin
        for T in (Float32,Float64),n in (2,7,8),sign in (-1,1)
            input=array(Complex{T}[complex(sin(T(.2i+.1j)),cos(T(.1i-.3j))) for i in 1:n,j in 1:n])
            wf=WaveFront(copy(input),T(1.6e-6),T(4e-5),zero(T),one(T));ctx=RunContext(wf)
            plans=(Proper.prepare_propagation(wf,T(sign*.004),ctx;kind=:ptp),
                Proper.prepare_propagation(wf,T(sign*.013),ctx;kind=:wts),
                Proper.prepare_propagation(wf,T(-sign*.009),ctx;kind=:stw))
            mask=array(Complex{T}[complex(T(.8+.01i),T(.02j)) for i in 1:n,j in 1:n])
            edge=array(fill(T(.83),n,n))
            multipliers=((mask,edge),(),(edge,mask));tail=(mask,edge)
            expected=similar(input);empty=similar(input,0,0)
            Proper.propagate!(expected,plans,input;input_multipliers=multipliers)
            for factor in tail;expected .*=factor;end
            operation()=begin
                first(plans).field .=input
                Proper._propagate_inplace!(plans,multipliers,tail,empty)
            end
            operation();synchronize()
            recorded=capture===nothing ? nothing : capture(operation)
            recorded===nothing || replay(recorded)
            synchronize()
            @test Array(first(plans).field)≈Array(expected) rtol=128eps(T) atol=128eps(T)
            retained=copy(first(plans).field)
            @test_throws ArgumentError Proper._propagate_inplace!(plans,(),tail,empty)
            @test_throws ArgumentError Proper._propagate_inplace!(plans,multipliers,(first(plans).field,),empty)
            @test_throws ArgumentError Proper._propagate_inplace!(plans,multipliers,tail,input)
            synchronize()
            @test Array(first(plans).field)==Array(retained)
        end
    end
end

function check_propagation_observer()
    @testset "Ordinary propagation observer" begin
        for T in (Float32,Float64),to_plane in (false,true)
            a=prop_begin(T(.008),T(1.55e-6),8;beam_diam_fraction=T(.5))
            b=prop_begin(T(.008),T(1.55e-6),8;beam_diam_fraction=T(.5))
            ca,cb=RunContext(a),RunContext(b)
            names=Symbol[];distances=T[]
            observer=function(operation,wf,dz,ctx)
                push!(names,nameof(operation));push!(distances,T(dz))
                operation(wf,dz,ctx)
            end
            for (focal,distance) in ((.08,.08),(-.06,-.02),(.04,0.0))
                prop_lens(a,T(focal),ca);prop_lens(b,T(focal),cb)
                prop_propagate(a,T(distance),ca;TO_PLANE=to_plane)
                Proper._prop_propagate!(b,T(distance),Proper.PropagateOptions(to_plane),cb,observer)
                @test a.field==b.field
                for field in fieldnames(typeof(a))[2:end-1]
                    @test getfield(a,field)==getfield(b,field)
                end
            end
            @test !isempty(names) && all(name->name in (:prop_ptp,:prop_wts,:prop_stw),names)
            @test all(isfinite,distances)
        end
    end
end
