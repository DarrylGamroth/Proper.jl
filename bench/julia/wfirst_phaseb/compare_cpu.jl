using FITSIO
using JSON3
using LinearAlgebra

function arg_value(flag::String, default=nothing)
    idx = findfirst(==(flag), ARGS)
    idx === nothing && return default
    idx < length(ARGS) || error("missing value for $(flag)")
    return ARGS[idx + 1]
end

function loadjson(path::AbstractString)
    isfile(path) || return nothing
    return JSON3.read(read(path, String))
end

function load_case_stack(prefix::AbstractString)
    idx = 1
    reals = Matrix{Float64}[]
    imags = Matrix{Float64}[]
    while isfile(prefix * "_$(idx)_real.fits")
        push!(reals, FITS(prefix * "_$(idx)_real.fits") do f
            read(f[1])
        end)
        push!(imags, FITS(prefix * "_$(idx)_imag.fits") do f
            read(f[1])
        end)
        idx += 1
    end
    isempty(reals) && error("no FITS outputs found for prefix $(prefix)")
    n = length(reals)
    sy, sx = size(reals[1])
    stack = Array{ComplexF64}(undef, n, sy, sx)
    for i in 1:n
        stack[i, :, :] = ComplexF64.(reals[i] .+ im .* imags[i])
    end
    return stack
end

timed_median_ms(report) = let timed = get(report.stats, :timed, true)
    if !timed || isnothing(report.stats.median_ns)
        nothing
    else
        Float64(report.stats.median_ns) / 1e6
    end
end

ratio_or_nothing(a, b) = (a === nothing || b === nothing) ? nothing : a / b
fmt_sci(x) = x === nothing ? "n/a" : string(round(Float64(x); sigdigits=4))
fmt_threads(x) = x === nothing ? "n/a" : string(Int(x))

function case_comparison(case_name::AbstractString)
    root = joinpath(@__DIR__, "..", "..", "reports")
    py_report = loadjson(joinpath(root, "python_wfirst_phaseb_$(case_name).json"))
    jl_report = loadjson(joinpath(root, "julia_wfirst_phaseb_$(case_name).json"))
    if py_report === nothing || jl_report === nothing
        return Dict(
            "case" => case_name,
            "available" => false,
            "python_median_ms" => nothing,
            "julia_median_ms" => nothing,
            "python_over_julia" => nothing,
            "relative_l2" => nothing,
            "max_abs_diff" => nothing,
            "max_sampling_abs_diff" => nothing,
            "threads" => nothing,
            "python_timed" => false,
            "julia_timed" => false,
            "reason" => "missing report artifact",
        )
    end
    py_prefix = joinpath(root, "python_wfirst_phaseb_$(case_name)")
    jl_prefix = joinpath(root, "julia_wfirst_phaseb_$(case_name)")
    if !isfile(py_prefix * "_1_real.fits") || !isfile(jl_prefix * "_1_real.fits")
        return Dict(
            "case" => case_name,
            "available" => false,
            "python_median_ms" => timed_median_ms(py_report),
            "julia_median_ms" => timed_median_ms(jl_report),
            "python_over_julia" => ratio_or_nothing(timed_median_ms(py_report), timed_median_ms(jl_report)),
            "relative_l2" => nothing,
            "max_abs_diff" => nothing,
            "max_sampling_abs_diff" => nothing,
            "threads" => hasproperty(jl_report.meta, :threads) ? Int(jl_report.meta.threads) : nothing,
            "python_timed" => Bool(get(py_report.stats, :timed, true)),
            "julia_timed" => Bool(get(jl_report.stats, :timed, true)),
            "reason" => "missing FITS output artifact",
        )
    end
    py_stack = load_case_stack(py_prefix)
    jl_stack = load_case_stack(jl_prefix)

    diff = jl_stack .- py_stack
    py_norm = norm(vec(py_stack))
    rel_l2 = py_norm == 0 ? norm(vec(diff)) : norm(vec(diff)) / py_norm
    max_abs = maximum(abs.(diff))
    py_sampling = Float64.(py_report.output.sampling)
    jl_sampling = Float64.(jl_report.output.sampling)
    sampling_abs = maximum(abs.(jl_sampling .- py_sampling))
    python_median_ms = timed_median_ms(py_report)
    julia_median_ms = timed_median_ms(jl_report)

    return Dict(
        "case" => case_name,
        "available" => true,
        "python_median_ms" => python_median_ms,
        "julia_median_ms" => julia_median_ms,
        "python_over_julia" => ratio_or_nothing(python_median_ms, julia_median_ms),
        "relative_l2" => rel_l2,
        "max_abs_diff" => max_abs,
        "max_sampling_abs_diff" => sampling_abs,
        "threads" => Int(jl_report.meta.threads),
        "python_timed" => Bool(get(py_report.stats, :timed, true)),
        "julia_timed" => Bool(get(jl_report.stats, :timed, true)),
    )
end

fmt_ms(x) = x === nothing ? "n/a" : string(round(Float64(x); digits=2), " ms")
fmt_ratio(x) = x === nothing ? "n/a" : string(round(Float64(x); digits=2), "x")

function main()
    cases_arg = arg_value("--cases", nothing)
    cases = cases_arg === nothing ? [
        "compact_hlc",
        "full_hlc",
        "compact_spc_spec_long",
        "full_spc_spec_long",
        "compact_spc_wide",
        "full_spc_wide",
    ] : filter!(!isempty, split(String(cases_arg), ','))
    rows = [case_comparison(case) for case in cases]
    out = Dict("cases" => rows)
    outpath = joinpath(@__DIR__, "..", "..", "reports", "wfirst_phaseb_cpu_comparison.json")
    open(outpath, "w") do io
        JSON3.write(io, out)
    end

    println("WFIRST Phase B CPU Comparison")
    println("============================")
    println(rpad("Case", 16), lpad("Python", 12), lpad("Julia", 12), lpad("Py/Jl", 10), lpad("RelL2", 14), lpad("MaxAbs", 14), lpad("dSampling", 14), lpad("Threads", 10))
    for row in rows
        print(
            rpad(String(row["case"]), 16),
            lpad(fmt_ms(row["python_median_ms"]), 12),
            lpad(fmt_ms(row["julia_median_ms"]), 12),
            lpad(fmt_ratio(row["python_over_julia"]), 10),
            lpad(fmt_sci(row["relative_l2"]), 14),
            lpad(fmt_sci(row["max_abs_diff"]), 14),
            lpad(fmt_sci(row["max_sampling_abs_diff"]), 14),
            lpad(fmt_threads(row["threads"]), 10),
        )
        if haskey(row, "reason") && row["reason"] !== nothing
            print("  # ", row["reason"])
        end
        println()
    end
end

main()
