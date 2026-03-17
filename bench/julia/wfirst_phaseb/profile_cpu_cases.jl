using JSON3
using Printf

function arg_value(flag::String, default=nothing)
    idx = findfirst(==(flag), ARGS)
    idx === nothing && return default
    idx < length(ARGS) || error("missing value for $(flag)")
    return ARGS[idx + 1]
end

function loadjson(path::AbstractString)
    isfile(path) || error("missing report artifact: $(path)")
    return JSON3.read(read(path, String))
end

fmt_ms(x) = x === nothing ? "n/a" : @sprintf("%.2f ms", Float64(x) / 1e6)
function fmt_bytes(x)
    x === nothing && return "n/a"
    v = Float64(x)
    if v >= 1024^2
        return @sprintf("%.2f MiB", v / 1024^2)
    elseif v >= 1024
        return @sprintf("%.2f KiB", v / 1024)
    end
    return @sprintf("%.0f B", v)
end

function profile_row(case_name::AbstractString)
    root = joinpath(@__DIR__, "..", "..", "reports")
    report = loadjson(joinpath(root, "julia_wfirst_phaseb_$(case_name).json"))
    stats = report.stats
    return Dict(
        "case" => case_name,
        "description" => String(report.meta.description),
        "timed" => Bool(stats.timed),
        "samples" => Int(stats.samples),
        "threads" => Int(report.meta.threads),
        "median_ns" => stats.median_ns === nothing ? nothing : Int(stats.median_ns),
        "min_ns" => stats.min_ns === nothing ? nothing : Int(stats.min_ns),
        "max_ns" => stats.max_ns === nothing ? nothing : Int(stats.max_ns),
        "median_alloc_bytes" => stats.median_alloc_bytes === nothing ? nothing : Int(stats.median_alloc_bytes),
        "min_alloc_bytes" => stats.min_alloc_bytes === nothing ? nothing : Int(stats.min_alloc_bytes),
        "max_alloc_bytes" => stats.max_alloc_bytes === nothing ? nothing : Int(stats.max_alloc_bytes),
    )
end

function main()
    cases_arg = arg_value("--cases", "compact_hlc,full_hlc,compact_spc_spec_long,full_spc_spec_long")
    cases = filter!(!isempty, split(String(cases_arg), ','))
    rows = [profile_row(case) for case in cases]

    out = Dict("cases" => rows)
    outpath = joinpath(@__DIR__, "..", "..", "reports", "wfirst_phaseb_cpu_profile.json")
    open(outpath, "w") do io
        JSON3.write(io, out)
    end

    println("WFIRST Phase B CPU Profile")
    println("=========================")
    println(rpad("Case", 24), lpad("Median", 12), lpad("Min", 12), lpad("Max", 12), lpad("Median alloc", 14), lpad("Min alloc", 14), lpad("Max alloc", 14), lpad("Samples", 10), lpad("Threads", 10))
    for row in rows
        println(
            rpad(String(row["case"]), 24),
            lpad(fmt_ms(row["median_ns"]), 12),
            lpad(fmt_ms(row["min_ns"]), 12),
            lpad(fmt_ms(row["max_ns"]), 12),
            lpad(fmt_bytes(row["median_alloc_bytes"]), 14),
            lpad(fmt_bytes(row["min_alloc_bytes"]), 14),
            lpad(fmt_bytes(row["max_alloc_bytes"]), 14),
            lpad(string(row["samples"]), 10),
            lpad(string(row["threads"]), 10),
        )
    end
end

main()
