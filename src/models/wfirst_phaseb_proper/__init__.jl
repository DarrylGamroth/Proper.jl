module WFIRSTPhaseBProper

using FFTW
using LinearAlgebra
import ..Proper
import ..Proper: PreparedModel, PreparedAssetPool, prepare_asset_pool, prepare_model, prop_add_phase,
    prop_begin, prop_circular_aperture, prop_define_entrance, prop_end!, prop_fits_read,
    prop_get_beamradius, prop_get_fratio, prop_get_sampling, prop_get_wavelength, prop_lens,
    prop_dm,
    prop_magnify!, prop_multiply, prop_propagate, prop_run, active_run_context

const map_dir = "/maps/"
const polfile = "/pol/new_toma"
const _data_dir_ref = Ref{String}(abspath(get(ENV, "WFIRST_PHASEB_DATA_ROOT", joinpath(pwd(), ".cache", "wfirst", "data_phaseb_from_roman_preflight"))))

phaseb_default_data_root() = _data_dir_ref[]
data_dir() = _data_dir_ref[]
set_data_dir!(path::AbstractString) = (_data_dir_ref[] = abspath(path))

@inline function _phaseb_python_fits(path::AbstractString)
    data = prop_fits_read(path)
    nd = ndims(data)
    nd <= 1 && return data
    return permutedims(data, Tuple(nd:-1:1))
end

include("trim.jl")
include("ffts.jl")
include("mft2.jl")
include("polmap.jl")
include("set_data_dir.jl")
include("copy_here.jl")
include("copy_examples_here.jl")
include("common.jl")
include("wfirst_phaseb_compact.jl")
include("wfirst_phaseb.jl")

export map_dir, polfile, data_dir, set_data_dir!, set_data_dir, phaseb_default_data_root
export copy_here, copy_examples_here
export trim, ffts, ffts!, mft2
export polmap, polab
export PhaseBHLCAssets, PhaseBHLCPreparedSharedAssets, PhaseBModelWorkspace, PhaseBPreparedAssets
export phaseb_case_definitions, load_phaseb_hlc_assets, prepare_phaseb_shared_assets
export prepare_phaseb_models, run_phaseb_case
export wfirst_phaseb_compact, wfirst_phaseb

end
