using JSON3

function loadjson(path)
    isfile(path) || return nothing
    return JSON3.read(read(path, String))
end

root = @__DIR__
py = loadjson(joinpath(root, "python_steady_state.json"))
jl = loadjson(joinpath(root, "julia_steady_state.json"))
ttfx = loadjson(joinpath(root, "julia_cold_start.json"))

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
