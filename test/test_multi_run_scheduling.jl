function _preallocated_multi_prescription(λm, n; wavefront, output)
    prop_begin!(wavefront, one(λm), λm; beam_diam_fraction=one(λm) / 2)
    return prop_end(wavefront, output)
end

function _multi_value_prescription(λm, n, value)
    T = typeof(λm)
    return fill(T(value), n, n), λm
end

function _preallocated_multi_fixture(::Type{T}, nruns::Integer; n::Integer=4) where {T<:AbstractFloat}
    stack = zeros(T, n, n, nruns)
    samplings = zeros(T, nruns)
    runs = map(1:nruns) do i
        context = RunContext(Matrix{T})
        field = Matrix{Complex{T}}(undef, n, n)
        wavelength_microns = T(0.5) + T(i - 1) * T(0.01)
        wavefront = prop_begin!(
            field,
            one(T),
            wavelength_microns * T(1e-6);
            context=context,
        )
        prepared = prepare_prescription(
            _preallocated_multi_prescription,
            wavelength_microns,
            n;
            context=context,
        )
        return prepare_run(
            prepared;
            activate_context=false,
            wavefront=wavefront,
            output=@view(stack[:, :, i]),
        )
    end
    return stack, samplings, runs
end

function _preallocated_multi_allocated(stack, samplings, runs)
    prop_run_multi!(stack, samplings, runs)
    return @allocated prop_run_multi!(stack, samplings, runs)
end

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

    # BitArray stores adjacent logical slices in shared packed words. Threaded
    # writes directly into a BitArray stack can therefore race even when each
    # task owns a distinct third-axis slice. Typed per-run staging followed by
    # serial stack assembly must preserve every bit.
    bit_items = collect(1:64)
    bit_passes = zeros(Int, length(bit_items))
    bit_expected = falses(3, 3, length(bit_items))
    for i in eachindex(bit_items)
        @views bit_expected[:, :, i] .= isodd(i)
    end
    bit_runner = function (_, _, slot)
        yield()
        return isodd(slot) ? trues(3, 3) : falses(3, 3), Float64(slot)
    end
    for _ in 1:32
        bit_stack, bit_samplings = Proper._prop_run_multi_items(
            Proper.MultiRunThreadedExecStyle(),
            bit_items,
            bit_passes,
            bit_runner,
        )
        @test bit_stack isa BitArray{3}
        @test bit_stack == bit_expected
        @test bit_samplings == Float64.(bit_items)
    end
end

@testset "prop_run_multi! caller-owned storage" begin
    single_stack, single_samplings, single_run = _preallocated_multi_fixture(Float64, 1)
    @test isconcretetype(eltype(single_run))
    @test Proper._is_exact_output_slice(single_run[1].kwargs.output, single_stack, 1)
    returned_stack, returned_samplings = @inferred prop_run_multi!(
        single_stack,
        single_samplings,
        single_run,
    )
    @test returned_stack === single_stack
    @test returned_samplings === single_samplings
    @test single_stack == ones(4, 4, 1)
    @test single_samplings == [0.5]
    single_allocated = _preallocated_multi_allocated(single_stack, single_samplings, single_run)
    @test single_allocated <= (VERSION < v"1.11.0-DEV" ? 128 : 0)

    stack, samplings, runs = _preallocated_multi_fixture(Float32, 4)
    @test Proper._multi_run_exec_style(runs) isa Proper.MultiRunThreadedExecStyle
    prop_run_multi!(stack, samplings, runs)
    @test stack == ones(Float32, 4, 4, 4)
    @test samplings == fill(0.5f0, 4)
    @test _preallocated_multi_allocated(stack, samplings, runs) <= 32_768

    allocated_stack, allocated_samplings = prop_run_multi(runs)
    @test allocated_stack == stack
    @test allocated_samplings == samplings

    wrapper_stack = zeros(Float32, 4, 4, 2)
    wrapper_samplings = zeros(Float32, 2)
    prop_run_multi!(
        wrapper_stack,
        wrapper_samplings,
        _multi_value_prescription,
        0.5f0,
        4;
        PASSVALUE=Float32[2, 3],
    )
    @test wrapper_stack[:, :, 1] == fill(2f0, 4, 4)
    @test wrapper_stack[:, :, 2] == fill(3f0, 4, 4)
    @test wrapper_samplings == fill(0.5f-6, 2)

    model = prepare_model(_multi_value_prescription, 0.5f0, 4; pool_size=2)
    fill!(wrapper_stack, 0)
    fill!(wrapper_samplings, 0)
    prop_run_multi!(
        wrapper_stack,
        wrapper_samplings,
        model;
        PASSVALUE=Float32[4, 5],
    )
    @test wrapper_stack[:, :, 1] == fill(4f0, 4, 4)
    @test wrapper_stack[:, :, 2] == fill(5f0, 4, 4)
    @test wrapper_samplings == fill(0.5f-6, 2)

    @test_throws ArgumentError prop_run_multi!(zeros(Float32, 4, 4, 3), zeros(Float32, 3), runs)
    @test_throws ArgumentError prop_run_multi!(zeros(Float32, 4, 4, 4), zeros(Float32, 3), runs)
    @test_throws ArgumentError prop_run_multi!(zeros(Float64, 4, 4, 1), zeros(Float64, 1), single_run; PASSVALUE=1)
    @test_throws ArgumentError prop_run_multi!(zeros(Float64, 4, 4, 1), zeros(Float64, 1), single_run; gain=2)
    @test_throws ArgumentError prop_run_multi!(zeros(Float32, 3, 3, 1), zeros(Float32, 1), single_run)
    @test_throws ArgumentError prop_run_multi!(zeros(Float32, 4, 4, 1), zeros(Float32, 1), single_run)

    bit_stack = falses(3, 3, 8)
    bit_samplings = zeros(Float64, 8)
    bit_style = Proper._multi_run_mutating_exec_style(
        Proper.MultiRunThreadedExecStyle(),
        bit_stack,
        bit_samplings,
    )
    @test bit_style isa Proper.MultiRunSerialExecStyle
    bit_runner = (_, _, slot) -> (isodd(slot) ? trues(3, 3) : falses(3, 3), Float64(slot))
    Proper._prop_run_multi_items!(
        bit_style,
        bit_stack,
        bit_samplings,
        collect(1:8),
        zeros(Int, 8),
        bit_runner,
    )
    @test vec(bit_stack[1, 1, :]) == isodd.(1:8)
    @test bit_samplings == Float64.(1:8)
end
