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

function render_markdown_table(headers::Vector{String}, rows::Vector{Vector{String}})
    isempty(rows) && return ""
    io = IOBuffer()
    println(io, "| ", join(headers, " | "), " |")
    println(io, "| ", join(fill("---", length(headers)), " | "), " |")
    for row in rows
        println(io, "| ", join(row, " | "), " |")
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
jl = loadjson(joinpath(root, "julia_steady_state.json"))
cpu_supported = loadjson(joinpath(root, "julia_supported_kernels.json"))
cuda_jl = loadjson(joinpath(root, "julia_cuda_steady_state.json"))
cuda_jl_fp64 = loadjson(joinpath(root, "julia_cuda_steady_state_fp64.json"))
cuda_jl_fp32 = loadjson(joinpath(root, "julia_cuda_steady_state_fp32.json"))
cuda_kernels = loadjson(joinpath(root, "cuda_supported_kernels.json"))
cuda_precision = loadjson(joinpath(root, "cuda_precision_split.json"))
cuda_isolated = loadjson(joinpath(root, "cuda_isolated_wavefront_kernels.json"))
amdgpu_jl = loadjson(joinpath(root, "julia_amdgpu_steady_state.json"))
amdgpu_kernels = loadjson(joinpath(root, "amdgpu_supported_kernels.json"))

summary_md_path = joinpath(root, "julia_cpu_gpu_summary.md")
generated_paths = String[]

term = IOBuffer()
md = IOBuffer()
println(md, "# Julia CPU vs GPU Benchmark Summary")
println(md)
println(md, "Generated from Julia CPU and CUDA reports in `bench/reports/`.")
println(md)

steady_headers = ["Backend", "Median", "Samples", "Vs Julia CPU"]
steady_rows = Vector{Vector{String}}()
steady_notes = String[]
jl_med = maybe_float(getpath(jl, "stats", "median_ns"))
cuda_available = maybe_bool(getpath(cuda_jl, "meta", "available"))
cuda_med = cuda_available ? maybe_float(getpath(cuda_jl, "stats", "median_ns")) : nothing
amdgpu_available = maybe_bool(getpath(amdgpu_jl, "meta", "available"))
amdgpu_med = amdgpu_available ? maybe_float(getpath(amdgpu_jl, "stats", "median_ns")) : nothing

if jl !== nothing
    push!(steady_rows, [
        "Julia CPU",
        fmt_ns(jl_med),
        fmt_int(getpath(jl, "stats", "samples")),
        "1.00x",
    ])
end
if cuda_available
    push!(steady_rows, [
        "Julia CUDA",
        fmt_ns(cuda_med),
        fmt_int(getpath(cuda_jl, "stats", "samples")),
        jl_med === nothing ? "-" : fmt_ratio(jl_med / cuda_med),
    ])
end
if amdgpu_available
    push!(steady_rows, [
        "Julia AMDGPU",
        fmt_ns(amdgpu_med),
        fmt_int(getpath(amdgpu_jl, "stats", "samples")),
        jl_med === nothing ? "-" : fmt_ratio(jl_med / amdgpu_med),
    ])
end
push!(steady_notes, "Steady-state rows use BenchmarkTools and exclude Julia cold-start / TTFx by construction.")
if cuda_available
    push!(steady_notes, "CUDA device: $(getpath(cuda_jl, "meta", "device"))")
    push!(steady_notes, "Standard CUDA steady-state row is sourced from the standalone FP64 workload report.")
elseif cuda_jl !== nothing
    reason = replace(String(getpath(cuda_jl, "reason")), '\n' => ' ')
    push!(steady_notes, "CUDA status: skipped")
    push!(steady_notes, "Reason: $reason")
end
if amdgpu_available
    push!(steady_notes, "AMDGPU device: $(getpath(amdgpu_jl, "meta", "device"))")
elseif amdgpu_jl !== nothing
    reason = replace(String(getpath(amdgpu_jl, "reason")), '\n' => ' ')
    push!(steady_notes, "AMDGPU status: skipped")
    push!(steady_notes, "Reason: $reason")
end
append_ascii_section!(term, "Julia Steady-State CPU vs GPU", steady_headers, steady_rows; aligns=[:l, :r, :r, :r], notes=steady_notes)
append_markdown_section!(md, "Julia Steady-State CPU vs GPU", steady_headers, steady_rows; notes=steady_notes)
push!(generated_paths, write_csv(joinpath(root, "julia_cpu_gpu_steady_state.csv"), steady_headers, steady_rows))

cpu_kernel_data = getpath(cpu_supported, "kernels")
cuda_kernel_data = cuda_available ? getpath(cuda_kernels, "kernels") : nothing
if cpu_kernel_data !== nothing && cuda_kernels !== nothing
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
    append_ascii_section!(term, "Julia Supported Kernels: CPU vs GPU", kernel_headers, kernel_rows; aligns=[:l, :r, :r, :r, :r, :r, :r, :r])
    append_markdown_section!(md, "Julia Supported Kernels: CPU vs GPU", kernel_headers, kernel_rows)
    push!(generated_paths, write_csv(joinpath(root, "julia_cpu_gpu_supported_kernels.csv"), kernel_headers, kernel_rows))
end

amdgpu_kernel_data = amdgpu_available ? getpath(amdgpu_kernels, "kernels") : nothing
if cpu_kernel_data !== nothing && amdgpu_kernels !== nothing
    kernel_headers = ["Kernel", "CPU", "AMDGPU", "CPU/AMDGPU", "CPU allocs", "AMDGPU allocs", "CPU bytes", "AMDGPU bytes"]
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
    for name in ordered_names(preferred, cpu_kernel_data, amdgpu_kernel_data)
        cpu_stats = getkey(cpu_kernel_data, name)
        gpu_stats = getkey(amdgpu_kernel_data, name)
        cpu_ns = maybe_float(getpath(cpu_stats, "median_ns"))
        gpu_ns = maybe_float(getpath(gpu_stats, "median_ns"))
        push!(kernel_rows, [
            name,
            fmt_ns(cpu_ns),
            fmt_ns(gpu_ns),
            cpu_ns === nothing || gpu_ns === nothing ? "-" : fmt_ratio(cpu_ns / gpu_ns),
            fmt_int(getpath(cpu_stats, "median_allocs")),
            fmt_int(getpath(gpu_stats, "median_allocs")),
            fmt_bytes(getpath(cpu_stats, "median_bytes")),
            fmt_bytes(getpath(gpu_stats, "median_bytes")),
        ])
    end
    append_ascii_section!(term, "Julia Supported Kernels: CPU vs AMDGPU", kernel_headers, kernel_rows; aligns=[:l, :r, :r, :r, :r, :r, :r, :r])
    append_markdown_section!(md, "Julia Supported Kernels: CPU vs AMDGPU", kernel_headers, kernel_rows)
    push!(generated_paths, write_csv(joinpath(root, "julia_cpu_amdgpu_supported_kernels.csv"), kernel_headers, kernel_rows))
end

cuda_precision_available = maybe_bool(getpath(cuda_precision, "meta", "available"))
if cuda_precision_available
    precision_headers = ["Metric", "FP64", "FP32", "FP64/FP32", "FP64 allocs", "FP32 allocs", "FP64 bytes", "FP32 bytes"]
    precision_rows = Vector{Vector{String}}()

    workload_fp64 = maybe_bool(getpath(cuda_jl_fp64, "meta", "available")) ? getpath(cuda_jl_fp64, "stats") : nothing
    workload_fp32_stats = maybe_bool(getpath(cuda_jl_fp32, "meta", "available")) ? getpath(cuda_jl_fp32, "stats") : nothing
    if workload_fp64 !== nothing || workload_fp32_stats !== nothing
        fp64_ns = maybe_float(getpath(workload_fp64, "median_ns"))
        fp32_ns = maybe_float(getpath(workload_fp32_stats, "median_ns"))
        push!(precision_rows, [
            "steady_state_workload",
            fmt_ns(fp64_ns),
            fmt_ns(fp32_ns),
            fp64_ns === nothing || fp32_ns === nothing ? "-" : fmt_ratio(fp64_ns / fp32_ns),
            fmt_int(getpath(workload_fp64, "median_allocs")),
            fmt_int(getpath(workload_fp32_stats, "median_allocs")),
            fmt_bytes(getpath(workload_fp64, "median_bytes")),
            fmt_bytes(getpath(workload_fp32_stats, "median_bytes")),
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
    push!(generated_paths, write_csv(joinpath(root, "julia_cuda_precision_split.csv"), precision_headers, precision_rows))
end

cuda_isolated_available = maybe_bool(getpath(cuda_isolated, "meta", "available"))
if cuda_isolated_available
    isolated_headers = ["Kernel", "FP64 host", "FP64 device", "FP32 host", "FP32 device", "Host64/Dev64", "Host32/Dev32"]
    isolated_rows = Vector{Vector{String}}()
    isolated_kernel_data = getpath(cuda_isolated, "kernels")
    for name in ordered_names(
        ["prop_qphase", "prop_ptp", "prop_wts", "prop_stw", "prop_circular_aperture", "prop_end_mutating"],
        isolated_kernel_data,
    )
        payload = getpath(cuda_isolated, "kernels", name)
        fp64 = getpath(payload, "fp64")
        fp32 = getpath(payload, "fp32")
        fp64_host = maybe_float(getpath(fp64, "host", "median_ns"))
        fp64_device = maybe_float(getpath(fp64, "timing", "device", "median_ns"))
        fp32_host = maybe_float(getpath(fp32, "host", "median_ns"))
        fp32_device = maybe_float(getpath(fp32, "timing", "device", "median_ns"))
        push!(isolated_rows, [
            name,
            fmt_ns(fp64_host),
            fmt_ns(fp64_device),
            fmt_ns(fp32_host),
            fmt_ns(fp32_device),
            fp64_host === nothing || fp64_device === nothing ? "-" : fmt_ratio(fp64_host / fp64_device),
            fp32_host === nothing || fp32_device === nothing ? "-" : fmt_ratio(fp32_host / fp32_device),
        ])
    end
    isolated_notes = [
        "Host columns use BenchmarkTools wall time and device columns use CUDA.@elapsed.",
        "Host/device ratios greater than 1.00x mean launch or synchronization overhead exceeds raw device execution time.",
    ]
    append_ascii_section!(term, "CUDA Isolated Wavefront Kernels", isolated_headers, isolated_rows; aligns=[:l, :r, :r, :r, :r, :r, :r], notes=isolated_notes)
    append_markdown_section!(md, "CUDA Isolated Wavefront Kernels", isolated_headers, isolated_rows; notes=isolated_notes)
    push!(generated_paths, write_csv(joinpath(root, "julia_cuda_isolated_wavefront_kernels.csv"), isolated_headers, isolated_rows))
end

println(term, "\nGenerated artifacts:")
for path in generated_paths
    println(term, "- ", relpath(path, root))
end

open(summary_md_path, "w") do io
    write(io, String(take!(md)))
end
println(String(take!(term)))
println("Wrote ", relpath(summary_md_path, root))
