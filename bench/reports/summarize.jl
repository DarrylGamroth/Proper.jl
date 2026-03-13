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

function render_ascii_table(headers::Vector{String}, rows::Vector{Vector{String}}; aligns=fill(:l, length(headers)))
    isempty(rows) && return ""
    widths = [length(h) for h in headers]
    for row in rows
        for i in eachindex(headers)
            widths[i] = max(widths[i], length(row[i]))
        end
    end

    io = IOBuffer()
    println(io, join((padcell(headers[i], widths[i], aligns[i]) for i in eachindex(headers)), "  "))
    println(io, join((repeat("-", widths[i]) for i in eachindex(headers)), "  "))
    for row in rows
        println(io, join((padcell(row[i], widths[i], aligns[i]) for i in eachindex(headers)), "  "))
    end
    return String(take!(io))
end

function md_escape(text::AbstractString)
    return replace(text, "|" => "\\|")
end

function render_markdown_table(headers::Vector{String}, rows::Vector{Vector{String}})
    isempty(rows) && return ""
    io = IOBuffer()
    println(io, "| ", join(md_escape.(headers), " | "), " |")
    println(io, "| ", join(fill("---", length(headers)), " | "), " |")
    for row in rows
        println(io, "| ", join(md_escape.(row), " | "), " |")
    end
    return String(take!(io))
end

function csv_escape(text::AbstractString)
    return "\"" * replace(text, "\"" => "\"\"") * "\""
end

function write_csv(path::AbstractString, headers::Vector{String}, rows::Vector{Vector{String}})
    open(path, "w") do io
        println(io, join(csv_escape.(headers), ","))
        for row in rows
            println(io, join(csv_escape.(row), ","))
        end
    end
    return path
end

function append_ascii_section!(io::IO, title::AbstractString, headers::Vector{String}, rows::Vector{Vector{String}}; aligns=fill(:l, length(headers)), notes=String[])
    isempty(rows) && return
    println(io)
    println(io, title)
    println(io, repeat("=", length(title)))
    print(io, render_ascii_table(headers, rows; aligns=aligns))
    for note in notes
        println(io)
        println(io, note)
    end
end

function append_markdown_section!(io::IO, title::AbstractString, headers::Vector{String}, rows::Vector{Vector{String}}; notes=String[])
    isempty(rows) && return
    println(io, "## ", title)
    println(io)
    print(io, render_markdown_table(headers, rows))
    println(io)
    for note in notes
        println(io, note)
    end
    println(io)
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
cuda_precision = loadjson(joinpath(root, "cuda_precision_split.json"))

summary_md_path = joinpath(root, "benchmark_summary.md")
generated_paths = String[]

term = IOBuffer()
md = IOBuffer()
println(md, "# Benchmark Summary")
println(md)
println(md, "Generated from JSON reports in `bench/reports/`.")
println(md)

steady_headers = ["Backend", "Median", "Samples", "Vs Python", "Vs Julia CPU"]
steady_rows = Vector{Vector{String}}()
steady_notes = String[]
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
push!(steady_notes, "Note: steady-state rows exclude Julia TTFx by construction.")
push!(steady_notes, "Ratio columns are reference/row, so values greater than 1.00x mean the row is faster.")
if ttfx !== nothing
    push!(steady_notes, "Julia cold start / TTFx: $(fmt_ns(getpath(ttfx, "first_call_ns")))")
end
if cuda_available
    push!(steady_notes, "CUDA device: $(getpath(cuda_jl, "meta", "device"))")
elseif cuda_jl !== nothing
    reason = replace(String(getpath(cuda_jl, "reason")), '\n' => ' ')
    push!(steady_notes, "CUDA status: skipped")
    push!(steady_notes, "Reason: $reason")
end
append_ascii_section!(term, "Steady-State Workload", steady_headers, steady_rows; aligns=[:l, :r, :r, :r, :r], notes=steady_notes)
append_markdown_section!(md, "Steady-State Workload", steady_headers, steady_rows; notes=steady_notes)
push!(generated_paths, write_csv(joinpath(root, "steady_state_comparison.csv"), steady_headers, steady_rows))

cpu_kernel_data = getpath(cpu_supported, "kernels")
cuda_kernel_data = cuda_available ? getpath(cuda_kernels, "kernels") : nothing
if cpu_kernel_data !== nothing
    kernel_headers = ["Kernel", "CPU", "CUDA", "CPU/CUDA", "CPU allocs", "CUDA allocs", "CPU bytes", "CUDA bytes"]
    kernel_rows = Vector{Vector{String}}()
    preferred = [
        "prop_qphase",
        "prop_ptp",
        "prop_wts",
        "prop_stw",
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
    for name in ordered_names(preferred, cpu_kernel_data, cuda_kernel_data)
        cpu_stats = getkey(cpu_kernel_data, name)
        cuda_stats = getkey(cuda_kernel_data, name)
        cpu_ns = maybe_float(getpath(cpu_stats, "median_ns"))
        cuda_ns = maybe_float(getpath(cuda_stats, "median_ns"))
        push!(kernel_rows, [
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
    append_ascii_section!(term, "Supported Kernels: CPU vs CUDA", kernel_headers, kernel_rows; aligns=[:l, :r, :r, :r, :r, :r, :r, :r])
    append_markdown_section!(md, "Supported Kernels: CPU vs CUDA", kernel_headers, kernel_rows)
    push!(generated_paths, write_csv(joinpath(root, "supported_kernels_comparison.csv"), kernel_headers, kernel_rows))
end

cuda_precision_available = maybe_bool(getpath(cuda_precision, "meta", "available"))
if cuda_precision_available
    precision_headers = ["Metric", "FP64", "FP32", "FP64/FP32", "FP64 allocs", "FP32 allocs", "FP64 bytes", "FP32 bytes"]
    precision_rows = Vector{Vector{String}}()

    workload_fp64 = getpath(cuda_precision, "workloads", "steady_state_fp64")
    workload_fp32 = getpath(cuda_precision, "workloads", "steady_state_fp32")
    if workload_fp64 !== nothing || workload_fp32 !== nothing
        fp64_ns = maybe_float(getpath(workload_fp64, "median_ns"))
        fp32_ns = maybe_float(getpath(workload_fp32, "median_ns"))
        push!(precision_rows, [
            "steady_state_workload",
            fmt_ns(fp64_ns),
            fmt_ns(fp32_ns),
            fp64_ns === nothing || fp32_ns === nothing ? "-" : fmt_ratio(fp64_ns / fp32_ns),
            fmt_int(getpath(workload_fp64, "median_allocs")),
            fmt_int(getpath(workload_fp32, "median_allocs")),
            fmt_bytes(getpath(workload_fp64, "median_bytes")),
            fmt_bytes(getpath(workload_fp32, "median_bytes")),
        ])
    end

    for name in ordered_names(
        ["prop_qphase", "prop_ptp", "prop_wts", "prop_stw", "prop_circular_aperture", "prop_end_mutating"],
        getpath(cuda_precision, "kernels"),
    )
        payload = getpath(cuda_precision, "kernels", name)
        fp64 = getpath(payload, "fp64")
        fp32 = getpath(payload, "fp32")
        local fp64_ns = maybe_float(getpath(fp64, "median_ns"))
        local fp32_ns = maybe_float(getpath(fp32, "median_ns"))
        push!(precision_rows, [
            name,
            fmt_ns(fp64_ns),
            fmt_ns(fp32_ns),
            fp64_ns === nothing || fp32_ns === nothing ? "-" : fmt_ratio(fp64_ns / fp32_ns),
            fmt_int(getpath(fp64, "median_allocs")),
            fmt_int(getpath(fp32, "median_allocs")),
            fmt_bytes(getpath(fp64, "median_bytes")),
            fmt_bytes(getpath(fp32, "median_bytes")),
        ])
    end

    precision_notes = ["Ratio column is FP64/FP32, so values greater than 1.00x mean the FP32 row is faster."]
    append_ascii_section!(term, "CUDA Precision Split", precision_headers, precision_rows; aligns=[:l, :r, :r, :r, :r, :r, :r, :r], notes=precision_notes)
    append_markdown_section!(md, "CUDA Precision Split", precision_headers, precision_rows; notes=precision_notes)
    push!(generated_paths, write_csv(joinpath(root, "cuda_precision_split.csv"), precision_headers, precision_rows))
elseif cuda_precision !== nothing
    reason = replace(String(getpath(cuda_precision, "reason")), '\n' => ' ')
    precision_headers = ["Metric", "Status"]
    precision_rows = [["cuda_precision_split", "skipped"]]
    precision_notes = ["Reason: $reason"]
    append_ascii_section!(term, "CUDA Precision Split", precision_headers, precision_rows; aligns=[:l, :l], notes=precision_notes)
    append_markdown_section!(md, "CUDA Precision Split", precision_headers, precision_rows; notes=precision_notes)
end

if phase2 !== nothing
    phase2_headers = ["Kernel", "Median", "Allocs", "Bytes", "Samples"]
    phase2_rows = Vector{Vector{String}}()
    kernel_data = getpath(phase2, "kernels")
    for name in ordered_names(["prop_lens", "prop_qphase", "prop_ptp"], kernel_data)
        stats = getkey(kernel_data, name)
        push!(phase2_rows, [
            name,
            fmt_ns(getpath(stats, "median_ns")),
            fmt_int(getpath(stats, "median_allocs")),
            fmt_bytes(getpath(stats, "median_bytes")),
            fmt_int(getpath(stats, "samples")),
        ])
    end
    append_ascii_section!(term, "Phase-2 CPU Kernels", phase2_headers, phase2_rows; aligns=[:l, :r, :r, :r, :r])
    append_markdown_section!(md, "Phase-2 CPU Kernels", phase2_headers, phase2_rows)
    push!(generated_paths, write_csv(joinpath(root, "phase2_cpu_kernels.csv"), phase2_headers, phase2_rows))
end

if refactor !== nothing
    refactor_headers = ["Kernel", "Wrapper", "Mutating", "Wrapper/Mut", "Bytes saved"]
    refactor_rows = Vector{Vector{String}}()
    for name in ordered_names(
        ["rectangle", "ellipse", "polygon", "irregular_polygon", "rotate", "resamplemap", "magnify_quick", "psd_errormap"],
        getpath(refactor, "pairs"),
    )
        payload = getpath(refactor, "pairs", name)
        push!(refactor_rows, [
            name,
            fmt_ns(getpath(payload, "wrapper", "median_ns")),
            fmt_ns(getpath(payload, "mutating", "median_ns")),
            fmt_ratio(getpath(payload, "speedup_wrapper_over_mutating")),
            fmt_bytes(getpath(payload, "median_byte_reduction")),
        ])
    end
    append_ascii_section!(term, "Wrapper vs Mutating", refactor_headers, refactor_rows; aligns=[:l, :r, :r, :r, :r])
    append_markdown_section!(md, "Wrapper vs Mutating", refactor_headers, refactor_rows)
    push!(generated_paths, write_csv(joinpath(root, "wrapper_vs_mutating.csv"), refactor_headers, refactor_rows))
end

if ka_interp !== nothing
    interp_headers = ["Kernel", "Loop", "KA", "Loop/KA"]
    interp_rows = Vector{Vector{String}}()
    for name in ordered_names(["cubic_conv_grid", "rotate_cubic", "rotate_linear"], getpath(ka_interp, "pairs"))
        payload = getpath(ka_interp, "pairs", name)
        push!(interp_rows, [
            name,
            fmt_ns(getpath(payload, "loop", "median_ns")),
            fmt_ns(getpath(payload, "ka", "median_ns")),
            fmt_ratio(getpath(payload, "speedup_loop_over_ka")),
        ])
    end
    append_ascii_section!(term, "KA Interpolation Pilot", interp_headers, interp_rows; aligns=[:l, :r, :r, :r])
    append_markdown_section!(md, "KA Interpolation Pilot", interp_headers, interp_rows)
    push!(generated_paths, write_csv(joinpath(root, "ka_interpolation_pilot.csv"), interp_headers, interp_rows))
end

if ka_geom !== nothing
    geom_headers = ["Kernel", "Loop", "KA", "Loop/KA"]
    geom_rows = Vector{Vector{String}}()
    for name in ordered_names(["rectangle", "ellipse", "irregular_polygon", "rounded_rectangle", "szoom", "pixellate"], getpath(ka_geom, "pairs"))
        payload = getpath(ka_geom, "pairs", name)
        push!(geom_rows, [
            name,
            fmt_ns(getpath(payload, "loop", "median_ns")),
            fmt_ns(getpath(payload, "ka", "median_ns")),
            fmt_ratio(getpath(payload, "speedup_loop_over_ka")),
        ])
    end
    append_ascii_section!(term, "KA Geometry/Sampling Pilot", geom_headers, geom_rows; aligns=[:l, :r, :r, :r])
    append_markdown_section!(md, "KA Geometry/Sampling Pilot", geom_headers, geom_rows)
    push!(generated_paths, write_csv(joinpath(root, "ka_geometry_sampling_pilot.csv"), geom_headers, geom_rows))
end

if examples !== nothing
    example_headers = ["Workflow", "Median", "Allocs", "Bytes", "Samples"]
    example_rows = Vector{Vector{String}}()
    for name in ordered_names(["simple_prescription_256", "simple_telescope_256", "psdtest_128"], getpath(examples, "examples"))
        stats = getpath(examples, "examples", name)
        push!(example_rows, [
            name,
            fmt_ns(getpath(stats, "median_ns")),
            fmt_int(getpath(stats, "median_allocs")),
            fmt_bytes(getpath(stats, "median_bytes")),
            fmt_int(getpath(stats, "samples")),
        ])
    end
    append_ascii_section!(term, "Example Workflows", example_headers, example_rows; aligns=[:l, :r, :r, :r, :r])
    append_markdown_section!(md, "Example Workflows", example_headers, example_rows)
    push!(generated_paths, write_csv(joinpath(root, "example_workflows.csv"), example_headers, example_rows))
end

open(summary_md_path, "w") do io
    write(io, String(take!(md)))
end
pushfirst!(generated_paths, summary_md_path)

print(String(take!(term)))
println()
println("Generated files")
println("===============")
for path in generated_paths
    println(relpath(path, pwd()))
end
