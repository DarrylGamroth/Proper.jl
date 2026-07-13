@testset "prop_run_multi backend scheduling" begin
    @test Proper._multi_run_exec_style(Proper.CPUBackend()) isa Proper.MultiRunThreadedExecStyle
    @test Proper._multi_run_exec_style(nothing) isa Proper.MultiRunThreadedExecStyle
    @test Proper._multi_run_exec_style(Proper.CUDABackend()) isa Proper.MultiRunSerialExecStyle
    @test Proper._multi_run_exec_style(Proper.AMDGPUBackend()) isa Proper.MultiRunSerialExecStyle
    @test Proper._multi_run_exec_style(Proper.UnknownBackend()) isa Proper.MultiRunSerialExecStyle

    @test Proper.cpu_inner_kernel_parallelism_enabled()
    @test Proper.cpu_inner_kernel_min_elems() == 0
    @test !Proper.with_cpu_inner_kernel_parallelism(false) do
        Proper.cpu_inner_kernel_parallelism_enabled()
    end
    @test Proper.cpu_inner_kernel_parallelism_enabled()
    @test_throws ErrorException Proper.with_cpu_inner_kernel_parallelism(false) do
        @test !Proper.cpu_inner_kernel_parallelism_enabled()
        error("restore inner-kernel policy")
    end
    @test Proper.cpu_inner_kernel_parallelism_enabled()

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

    inner_enabled = fill(true, 4)
    inner_min_elems = fill(0, 4)
    nested_runner = function (item, pass, slot)
        inner_enabled[slot] = Proper.cpu_inner_kernel_parallelism_enabled()
        inner_min_elems[slot] = Proper.cpu_inner_kernel_min_elems()
        return fill(Float32(item + pass), 2, 3), Float32(slot)
    end
    Proper._prop_run_multi_items(
        Proper.MultiRunThreadedExecStyle(),
        collect(1:4),
        collect(10:13),
        nested_runner,
    )
    @test all(inner_enabled)
    expected_min_elems = Threads.nthreads() > 1 && length(inner_enabled) >= Threads.nthreads() ? 262_144 : 0
    @test all(==(expected_min_elems), inner_min_elems)
end
