@testset "prop_run_multi backend scheduling" begin
    @test Proper._multi_run_exec_style(Proper.CPUBackend()) isa Proper.MultiRunThreadedExecStyle
    @test Proper._multi_run_exec_style(nothing) isa Proper.MultiRunThreadedExecStyle
    @test Proper._multi_run_exec_style(Proper.CUDABackend()) isa Proper.MultiRunSerialExecStyle
    @test Proper._multi_run_exec_style(Proper.AMDGPUBackend()) isa Proper.MultiRunSerialExecStyle
    @test Proper._multi_run_exec_style(Proper.UnknownBackend()) isa Proper.MultiRunSerialExecStyle

    ctx = RunContext(Matrix{Float32})
    independent = Proper.prepared_contexts(
        prepare_prescription((λ, n; kwargs...) -> prop_begin(1.0f0, λ, n), 0.55f0, 8; context=ctx),
        2,
    )
    @test Proper._multi_run_exec_style(independent) isa Proper.MultiRunThreadedExecStyle
    @test Proper._multi_run_exec_style([ctx, ctx]) isa Proper.MultiRunSerialExecStyle

    calls = Int[]
    runner = function (item, pass, slot)
        push!(calls, slot)
        return fill(Float32(item + pass), 2, 3), Float32(slot)
    end
    stack, samplings = Proper._prop_run_multi_items(
        Proper.MultiRunSerialExecStyle(),
        collect(1:4),
        collect(10:13),
        runner,
    )
    @test calls == collect(1:4)
    @test size(stack) == (2, 3, 4)
    @test vec(stack[1, 1, :]) == Float32[11, 13, 15, 17]
    @test samplings == Float32[1, 2, 3, 4]
end
