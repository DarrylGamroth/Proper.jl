using JSON3
using Proper
using Proper.WFIRSTPhaseBProper

function arg_value(flag::String, default=nothing)
    idx = findfirst(==(flag), ARGS)
    idx === nothing && return default
    idx < length(ARGS) || error("missing value for $(flag)")
    return ARGS[idx + 1]
end

function stage_row(case_name::AbstractString, data_root::AbstractString)
    cases = phaseb_case_definitions()
    haskey(cases, case_name) || error("unsupported case $(case_name)")
    base_case = cases[case_name]
    timer = Proper.WFIRSTPhaseBProper.PhaseBStageTimer()
    prof_case = merge(base_case, (passvalue=merge(base_case.passvalue, Dict("stage_timer" => timer)),))
    models, _assets = prepare_phaseb_models(prof_case; data_root=data_root)
    run_phaseb_case(models; threaded=false)
    return Dict(
        "case" => case_name,
        "description" => String(base_case.description),
        "wavelength_count" => length(base_case.wavelengths_m),
        "stage_report" => Proper.WFIRSTPhaseBProper.phaseb_stage_report(timer),
    )
end

function main()
    cases_arg = String(arg_value("--cases", "full_spc_spec_long"))
    data_root = String(arg_value("--data-root", phaseb_default_data_root()))
    cases = filter!(!isempty, split(cases_arg, ','))
    rows = [stage_row(case_name, data_root) for case_name in cases]

    out = Dict(
        "meta" => Dict(
            "julia_version" => string(VERSION),
            "backend" => "cpu",
            "threaded" => false,
            "data_root" => abspath(data_root),
        ),
        "cases" => rows,
    )

    outpath = joinpath(@__DIR__, "..", "..", "reports", "wfirst_phaseb_stage_profile.json")
    open(outpath, "w") do io
        JSON3.write(io, out)
    end

    println(JSON3.pretty(JSON3.read(JSON3.write(out))))
end

main()
