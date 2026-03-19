module WFIRSTPhaseBCPUModels

using Proper

@isdefined(WFIRSTPhaseBProper) || include(joinpath(@__DIR__, "..", "..", "reference_models", "wfirst_phaseb_proper", "__init__.jl"))
using .WFIRSTPhaseBProper: phaseb_case_definitions, phaseb_default_data_root, load_phaseb_hlc_assets,
    prepare_phaseb_models, run_phaseb_case, wfirst_phaseb_compact, wfirst_phaseb

export phaseb_case_definitions, phaseb_default_data_root, load_phaseb_hlc_assets
export prepare_phaseb_models, run_phaseb_case
export wfirst_phaseb_compact, wfirst_phaseb

end
