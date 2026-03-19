using Proper

@isdefined(WFIRSTPhaseBProper) || include(joinpath(@__DIR__, "..", "reference_models", "wfirst_phaseb_proper", "__init__.jl"))
using .WFIRSTPhaseBProper

function wfirst_phaseb_reference_demo(; case_name::AbstractString="compact_hlc", threaded::Bool=true, data_root::AbstractString=phaseb_default_data_root())
    cases = phaseb_case_definitions()
    haskey(cases, case_name) || throw(ArgumentError("unsupported WFIRST Phase B case $(case_name)"))
    case = cases[case_name]
    models, _ = prepare_phaseb_models(case; data_root=data_root)
    stack, samplings = run_phaseb_case(models; threaded=threaded)
    println("WFIRST Phase B reference example")
    println("  case = ", case_name)
    println("  wavelengths = ", case.wavelengths_um, " um")
    println("  output size = ", size(stack))
    println("  samplings = ", samplings)
    println("  peak abs = ", maximum(abs, stack))
    return stack, samplings
end

if abspath(PROGRAM_FILE) == @__FILE__
    wfirst_phaseb_reference_demo()
end
