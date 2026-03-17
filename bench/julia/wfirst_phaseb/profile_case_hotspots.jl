using Profile
using Proper
using Proper.WFIRSTPhaseBProper

function arg_value(flag::String, default=nothing)
    idx = findfirst(==(flag), ARGS)
    idx === nothing && return default
    idx < length(ARGS) || error("missing value for $(flag)")
    return ARGS[idx + 1]
end

function main()
    case_name = String(arg_value("--case", "full_spc_spec_long"))
    data_root = String(arg_value("--data-root", phaseb_default_data_root()))
    mincount = parse(Int, String(arg_value("--mincount", "25")))
    maxdepth = parse(Int, String(arg_value("--maxdepth", "32")))

    cases = phaseb_case_definitions()
    haskey(cases, case_name) || error("unsupported case $(case_name)")
    case = cases[case_name]
    models, _assets = prepare_phaseb_models(case; data_root=data_root)
    workload() = run_phaseb_case(models; threaded=false)

    workload()
    Profile.init(; n=10_000_000, delay=0.001, limitwarn=false)
    Profile.clear()
    Profile.@profile workload()

    outdir = joinpath(@__DIR__, "..", "..", "reports")
    mkpath(outdir)
    outpath = joinpath(outdir, "wfirst_phaseb_hotspots_$(case_name).txt")
    open(outpath, "w") do io
        println(io, "WFIRST Phase B Hotspots: $(case_name)")
        println(io, "data_root=$(abspath(data_root))")
        println(io, "threaded=false")
        println(io)
        Profile.print(io; format=:flat, sortedby=:count, mincount=mincount, maxdepth=maxdepth)
    end

    println(outpath)
    Profile.print(stdout; format=:flat, sortedby=:count, mincount=mincount, maxdepth=maxdepth)
end

main()
