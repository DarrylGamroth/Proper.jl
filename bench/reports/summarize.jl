using JSON3

function loadjson(path)
    isfile(path) || return nothing
    return JSON3.read(read(path, String))
end

root = @__DIR__
py = loadjson(joinpath(root, "python_steady_state.json"))
jl = loadjson(joinpath(root, "julia_steady_state.json"))
ttfx = loadjson(joinpath(root, "julia_cold_start.json"))
phase2 = loadjson(joinpath(root, "phase2_kernels.json"))
refactor = loadjson(joinpath(root, "refactor_kernels.json"))
ka_interp = loadjson(joinpath(root, "ka_interp_kernels.json"))
ka_geom = loadjson(joinpath(root, "ka_geometry_sampling_kernels.json"))
examples = loadjson(joinpath(root, "example_workflows.json"))
cuda_jl = loadjson(joinpath(root, "julia_cuda_steady_state.json"))
cuda_kernels = loadjson(joinpath(root, "cuda_supported_kernels.json"))

if py !== nothing && jl !== nothing
    py_med = Float64(py["stats"]["median_ns"])
    jl_med = Float64(jl["stats"]["median_ns"])
    speedup = py_med / jl_med
    println("# Steady-State Comparison")
    println("Python median ns: ", py_med)
    println("Julia median ns: ", jl_med)
    println("Speedup (Python/Julia): ", speedup)
    println("Note: Julia TTFx excluded from this comparison.")
end

if ttfx !== nothing
    println("\n# Cold Start / TTFx")
    println("Julia first call ns: ", ttfx["first_call_ns"])
end

if phase2 !== nothing
    println("\n# Phase-2 Kernel Matrix")
    for (name, stats) in pairs(phase2["kernels"])
        println(name, " median ns: ", stats["median_ns"], " median bytes: ", stats["median_bytes"])
    end
end

if refactor !== nothing
    println("\n# Refactor Kernel Deltas")
    for (name, payload) in pairs(refactor["pairs"])
        println(
            name,
            " speedup(wrapper/mutating): ",
            payload["speedup_wrapper_over_mutating"],
            " median byte reduction: ",
            payload["median_byte_reduction"],
        )
    end
    hotspots = refactor["hotspots"]
    if haskey(hotspots, "psd_errormap_no_apply")
        hs = hotspots["psd_errormap_no_apply"]
        println(
            "psd_errormap_no_apply median ns: ",
            hs["median_ns"],
            " median bytes: ",
            hs["median_bytes"],
        )
    else
        if haskey(hotspots, "psd_errormap_no_apply_wrapper")
            hs_wrap = hotspots["psd_errormap_no_apply_wrapper"]
            println(
                "psd_errormap_no_apply_wrapper median ns: ",
                hs_wrap["median_ns"],
                " median bytes: ",
                hs_wrap["median_bytes"],
            )
        end
        if haskey(hotspots, "psd_errormap_no_apply_mutating")
            hs_mut = hotspots["psd_errormap_no_apply_mutating"]
            println(
                "psd_errormap_no_apply_mutating median ns: ",
                hs_mut["median_ns"],
                " median bytes: ",
                hs_mut["median_bytes"],
            )
        end
    end
end

if ka_interp !== nothing
    println("\n# KA Interp Pilot")
    for (name, payload) in pairs(ka_interp["pairs"])
        println(
            name,
            " speedup(loop/ka): ",
            payload["speedup_loop_over_ka"],
        )
    end
    for (name, stats) in pairs(ka_interp["public"])
        println(name, " median ns: ", stats["median_ns"], " median bytes: ", stats["median_bytes"])
    end
end

if ka_geom !== nothing
    println("\n# KA Geometry/Sampling Pilot")
    for (name, payload) in pairs(ka_geom["pairs"])
        println(
            name,
            " speedup(loop/ka): ",
            payload["speedup_loop_over_ka"],
        )
    end
end

if examples !== nothing
    println("\n# Example Workflow Matrix")
    for (name, stats) in pairs(examples["examples"])
        println(name, " median ns: ", stats["median_ns"], " median bytes: ", stats["median_bytes"])
    end
end

if cuda_jl !== nothing
    println("\n# CUDA Steady-State")
    if Bool(cuda_jl["meta"]["available"])
        println("CUDA Julia median ns: ", cuda_jl["stats"]["median_ns"])
        if haskey(cuda_jl["meta"], "device")
            println("CUDA device: ", cuda_jl["meta"]["device"])
        end
    else
        println("Skipped: ", cuda_jl["reason"])
    end
end

if cuda_kernels !== nothing
    println("\n# CUDA Supported Kernels")
    if Bool(cuda_kernels["meta"]["available"])
        for (name, stats) in pairs(cuda_kernels["kernels"])
            println(name, " median ns: ", stats["median_ns"], " median bytes: ", stats["median_bytes"])
        end
    else
        println("Skipped: ", cuda_kernels["reason"])
    end
end
