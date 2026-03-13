using JSON3
using Printf

function loadjson(path::AbstractString)
    isfile(path) || return nothing
    return JSON3.read(read(path, String))
end

function getkey(obj, key)
    obj === nothing && return nothing
    return haskey(obj, key) ? obj[key] : nothing
end

function getpath(obj, keys...)
    cur = obj
    for key in keys
        cur = getkey(cur, key)
        cur === nothing && return nothing
    end
    return cur
end

maybe_float(x) = x === nothing ? nothing : Float64(x)
maybe_int(x) = x === nothing ? nothing : Int(round(Float64(x)))
maybe_bool(x) = x === nothing ? false : Bool(x)

function fmt_ns(x)
    x === nothing && return "-"
    v = Float64(x)
    if v >= 1e9
        return @sprintf("%.2f s", v / 1e9)
    elseif v >= 1e6
        return @sprintf("%.2f ms", v / 1e6)
    elseif v >= 1e3
        return @sprintf("%.2f us", v / 1e3)
    end
    return @sprintf("%.0f ns", v)
end

function fmt_bytes(x)
    x === nothing && return "-"
    v = Float64(x)
    if v >= 1024^2
        return @sprintf("%.2f MiB", v / 1024^2)
    elseif v >= 1024
        return @sprintf("%.2f KiB", v / 1024)
    end
    return @sprintf("%.0f B", v)
end

function fmt_int(x)
    x === nothing && return "-"
    return string(maybe_int(x))
end

function fmt_ratio(x)
    x === nothing && return "-"
    return @sprintf("%.2fx", Float64(x))
end

function padcell(text::AbstractString, width::Int, align::Symbol)
    return align === :r ? lpad(text, width) : rpad(text, width)
end

function print_table(headers::Vector{String}, rows::Vector{Vector{String}}; aligns=fill(:l, length(headers)))
    isempty(rows) && return
    widths = [length(h) for h in headers]
    for row in rows
        for i in eachindex(headers)
            widths[i] = max(widths[i], length(row[i]))
        end
    end
    println(join((padcell(headers[i], widths[i], aligns[i]) for i in eachindex(headers)), "  "))
    println(join((repeat("-", widths[i]) for i in eachindex(headers)), "  "))
    for row in rows
        println(join((padcell(row[i], widths[i], aligns[i]) for i in eachindex(headers)), "  "))
    end
end

function print_section(title::AbstractString)
    println()
    println(title)
    println(repeat("=", length(title)))
end

function ordered_names(preferred::Vector{String}, dicts...)
    names = String[]
    seen = Set{String}()
    for name in preferred
        if any(d -> d !== nothing && haskey(d, name), dicts)
            push!(names, name)
            push!(seen, name)
        end
    end
    extra = String[]
    for d in dicts
        d === nothing && continue
        for name in keys(d)
            s = String(name)
            s in seen && continue
            push!(extra, s)
            push!(seen, s)
        end
    end
    append!(names, sort(extra))
    return names
end

root = @__DIR__
py = loadjson(joinpath(root, "python_steady_state.json"))
jl = loadjson(joinpath(root, "julia_steady_state.json"))
ttfx = loadjson(joinpath(root, "julia_cold_start.json"))
cpu_supported = loadjson(joinpath(root, "julia_supported_kernels.json"))
phase2 = loadjson(joinpath(root, "phase2_kernels.json"))
refactor = loadjson(joinpath(root, "refactor_kernels.json"))
ka_interp = loadjson(joinpath(root, "ka_interp_kernels.json"))
ka_geom = loadjson(joinpath(root, "ka_geometry_sampling_kernels.json"))
examples = loadjson(joinpath(root, "example_workflows.json"))
cuda_jl = loadjson(joinpath(root, "julia_cuda_steady_state.json"))
cuda_kernels = loadjson(joinpath(root, "cuda_supported_kernels.json"))

print_section("Steady-State Workload")
steady_rows = Vector{Vector{String}}()
py_med = maybe_float(getpath(py, "stats", "median_ns"))
jl_med = maybe_float(getpath(jl, "stats", "median_ns"))
cuda_available = maybe_bool(getpath(cuda_jl, "meta", "available"))
cuda_med = cuda_available ? maybe_float(getpath(cuda_jl, "stats", "median_ns")) : nothing

if py !== nothing
    push!(steady_rows, [
        "Python CPU",
        fmt_ns(py_med),
        fmt_int(getpath(py, "stats", "samples")),
        "1.00x",
        jl_med === nothing ? "-" : fmt_ratio(jl_med / py_med),
    ])
end
if jl !== nothing
    push!(steady_rows, [
        "Julia CPU",
        fmt_ns(jl_med),
        fmt_int(getpath(jl, "stats", "samples")),
        py_med === nothing ? "-" : fmt_ratio(py_med / jl_med),
        "1.00x",
    ])
end
if cuda_available
    push!(steady_rows, [
        "Julia CUDA",
        fmt_ns(cuda_med),
        fmt_int(getpath(cuda_jl, "stats", "samples")),
        py_med === nothing ? "-" : fmt_ratio(py_med / cuda_med),
        jl_med === nothing ? "-" : fmt_ratio(jl_med / cuda_med),
    ])
end
print_table(
    ["Backend", "Median", "Samples", "Vs Python", "Vs Julia CPU"],
    steady_rows;
    aligns=[:l, :r, :r, :r, :r],
)
println()
println("Note: steady-state rows exclude Julia TTFx by construction.")
println("Ratio columns are reference/row, so values greater than 1.00x mean the row is faster.")
if ttfx !== nothing
    println("Julia cold start / TTFx: ", fmt_ns(getpath(ttfx, "first_call_ns")))
end
if cuda_available
    println("CUDA device: ", getpath(cuda_jl, "meta", "device"))
elseif cuda_jl !== nothing
    println("CUDA status: skipped")
    println("Reason: ", replace(String(getpath(cuda_jl, "reason")), '\n' => ' '))
end

cpu_kernel_data = getpath(cpu_supported, "kernels")
cuda_kernel_data = cuda_available ? getpath(cuda_kernels, "kernels") : nothing
if cpu_kernel_data !== nothing
    print_section("Supported Kernels: CPU vs CUDA")
    preferred = [
        "prop_qphase",
        "prop_ptp",
        "prop_circular_aperture",
        "prop_end_mutating",
        "prop_rotate_mutating",
        "prop_magnify_quick_mutating",
        "prop_szoom_mutating",
        "prop_pixellate_mutating",
        "prop_resamplemap_mutating",
        "prop_rectangle_mutating",
        "prop_rounded_rectangle_mutating",
    ]
    rows = Vector{Vector{String}}()
    for name in ordered_names(preferred, cpu_kernel_data, cuda_kernel_data)
        cpu_stats = getkey(cpu_kernel_data, name)
        cuda_stats = getkey(cuda_kernel_data, name)
        cpu_ns = maybe_float(getpath(cpu_stats, "median_ns"))
        cuda_ns = maybe_float(getpath(cuda_stats, "median_ns"))
        push!(rows, [
            name,
            fmt_ns(cpu_ns),
            fmt_ns(cuda_ns),
            cpu_ns === nothing || cuda_ns === nothing ? "-" : fmt_ratio(cpu_ns / cuda_ns),
            fmt_int(getpath(cpu_stats, "median_allocs")),
            fmt_int(getpath(cuda_stats, "median_allocs")),
            fmt_bytes(getpath(cpu_stats, "median_bytes")),
            fmt_bytes(getpath(cuda_stats, "median_bytes")),
        ])
    end
    print_table(
        ["Kernel", "CPU", "CUDA", "CPU/CUDA", "CPU allocs", "CUDA allocs", "CPU bytes", "CUDA bytes"],
        rows;
        aligns=[:l, :r, :r, :r, :r, :r, :r, :r],
    )
end

if phase2 !== nothing
    print_section("Phase-2 CPU Kernels")
    kernel_data = getpath(phase2, "kernels")
    rows = Vector{Vector{String}}()
    for name in ordered_names(["prop_lens", "prop_qphase", "prop_ptp"], kernel_data)
        stats = getkey(kernel_data, name)
        push!(rows, [
            name,
            fmt_ns(getpath(stats, "median_ns")),
            fmt_int(getpath(stats, "median_allocs")),
            fmt_bytes(getpath(stats, "median_bytes")),
            fmt_int(getpath(stats, "samples")),
        ])
    end
    print_table(["Kernel", "Median", "Allocs", "Bytes", "Samples"], rows; aligns=[:l, :r, :r, :r, :r])
end

if refactor !== nothing
    print_section("Wrapper vs Mutating")
    rows = Vector{Vector{String}}()
    for name in ordered_names(
        ["rectangle", "ellipse", "polygon", "irregular_polygon", "rotate", "resamplemap", "magnify_quick", "psd_errormap"],
        getpath(refactor, "pairs"),
    )
        payload = getpath(refactor, "pairs", name)
        push!(rows, [
            name,
            fmt_ns(getpath(payload, "wrapper", "median_ns")),
            fmt_ns(getpath(payload, "mutating", "median_ns")),
            fmt_ratio(getpath(payload, "speedup_wrapper_over_mutating")),
            fmt_bytes(getpath(payload, "median_byte_reduction")),
        ])
    end
    print_table(
        ["Kernel", "Wrapper", "Mutating", "Wrapper/Mut", "Bytes saved"],
        rows;
        aligns=[:l, :r, :r, :r, :r],
    )
end

if ka_interp !== nothing
    print_section("KA Interpolation Pilot")
    rows = Vector{Vector{String}}()
    for name in ordered_names(["cubic_conv_grid", "rotate_cubic", "rotate_linear"], getpath(ka_interp, "pairs"))
        payload = getpath(ka_interp, "pairs", name)
        push!(rows, [
            name,
            fmt_ns(getpath(payload, "loop", "median_ns")),
            fmt_ns(getpath(payload, "ka", "median_ns")),
            fmt_ratio(getpath(payload, "speedup_loop_over_ka")),
        ])
    end
    print_table(["Kernel", "Loop", "KA", "Loop/KA"], rows; aligns=[:l, :r, :r, :r])
end

if ka_geom !== nothing
    print_section("KA Geometry/Sampling Pilot")
    rows = Vector{Vector{String}}()
    for name in ordered_names(["rectangle", "ellipse", "irregular_polygon", "rounded_rectangle", "szoom", "pixellate"], getpath(ka_geom, "pairs"))
        payload = getpath(ka_geom, "pairs", name)
        push!(rows, [
            name,
            fmt_ns(getpath(payload, "loop", "median_ns")),
            fmt_ns(getpath(payload, "ka", "median_ns")),
            fmt_ratio(getpath(payload, "speedup_loop_over_ka")),
        ])
    end
    print_table(["Kernel", "Loop", "KA", "Loop/KA"], rows; aligns=[:l, :r, :r, :r])
end

if examples !== nothing
    print_section("Example Workflows")
    rows = Vector{Vector{String}}()
    for name in ordered_names(["simple_prescription_256", "simple_telescope_256", "psdtest_128"], getpath(examples, "examples"))
        stats = getpath(examples, "examples", name)
        push!(rows, [
            name,
            fmt_ns(getpath(stats, "median_ns")),
            fmt_int(getpath(stats, "median_allocs")),
            fmt_bytes(getpath(stats, "median_bytes")),
            fmt_int(getpath(stats, "samples")),
        ])
    end
    print_table(["Workflow", "Median", "Allocs", "Bytes", "Samples"], rows; aligns=[:l, :r, :r, :r, :r])
end
