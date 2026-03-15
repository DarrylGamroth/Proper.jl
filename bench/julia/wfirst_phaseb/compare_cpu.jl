using FITSIO
using JSON3
using LinearAlgebra

function loadjson(path::AbstractString)
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

function case_comparison(case_name::AbstractString)
    root = joinpath(@__DIR__, "..", "..", "reports")
    py_report = loadjson(joinpath(root, "python_wfirst_phaseb_$(case_name).json"))
    jl_report = loadjson(joinpath(root, "julia_wfirst_phaseb_$(case_name).json"))
    py_stack = load_case_stack(joinpath(root, "python_wfirst_phaseb_$(case_name)"))
    jl_stack = load_case_stack(joinpath(root, "julia_wfirst_phaseb_$(case_name)"))

    diff = jl_stack .- py_stack
    py_norm = norm(vec(py_stack))
    rel_l2 = py_norm == 0 ? norm(vec(diff)) : norm(vec(diff)) / py_norm
    max_abs = maximum(abs.(diff))
    py_sampling = Float64.(py_report.output.sampling)
    jl_sampling = Float64.(jl_report.output.sampling)
    sampling_abs = maximum(abs.(jl_sampling .- py_sampling))

    return Dict(
        "case" => case_name,
        "python_median_ms" => Float64(py_report.stats.median_ns) / 1e6,
        "julia_median_ms" => Float64(jl_report.stats.median_ns) / 1e6,
        "python_over_julia" => Float64(py_report.stats.median_ns) / Float64(jl_report.stats.median_ns),
        "relative_l2" => rel_l2,
        "max_abs_diff" => max_abs,
        "max_sampling_abs_diff" => sampling_abs,
        "threads" => Int(jl_report.meta.threads),
    )
end

function main()
    cases = ["compact_hlc", "full_hlc"]
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
        println(
            rpad(String(row["case"]), 16),
            lpad(string(round(row["python_median_ms"]; digits=2), " ms"), 12),
            lpad(string(round(row["julia_median_ms"]; digits=2), " ms"), 12),
            lpad(string(round(row["python_over_julia"]; digits=2), "x"), 10),
            lpad(string(round(row["relative_l2"]; sigdigits=4)), 14),
            lpad(string(round(row["max_abs_diff"]; sigdigits=4)), 14),
            lpad(string(round(row["max_sampling_abs_diff"]; sigdigits=4)), 14),
            lpad(string(row["threads"]), 10),
        )
    end
end

main()
