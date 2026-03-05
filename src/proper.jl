module proper

using LinearAlgebra
using FFTW
using Statistics

include("core/policy.jl")
include("core/traits.jl")
include("core/types.jl")
include("core/context.jl")
include("core/keywords.jl")
include("core/grid.jl")
include("core/resample.jl")
include("internal/unimplemented.jl")

include("libcconv.jl")
include("libcconvthread.jl")
include("libszoom.jl")
include("prop_8th_order_mask.jl")
include("prop_add_phase.jl")
include("prop_add_wavefront.jl")
include("prop_begin.jl")
include("prop_circular_aperture.jl")
include("prop_circular_obscuration.jl")
include("prop_compile_c.jl")
include("prop_cubic_conv.jl")
include("prop_define_entrance.jl")
include("prop_dftidefs.jl")
include("prop_divide.jl")
include("prop_dm.jl")
include("prop_ellipse.jl")
include("prop_elliptical_aperture.jl")
include("prop_elliptical_obscuration.jl")
include("prop_end.jl")
include("prop_end_savestate.jl")
include("prop_errormap.jl")
include("prop_execute_multi.jl")
include("prop_ffti.jl")
include("prop_fftw.jl")
include("prop_fftw_wisdom.jl")
include("prop_fit_dm.jl")
include("prop_fit_zernikes.jl")
include("prop_fits_read.jl")
include("prop_fits_write.jl")
include("prop_get_amplitude.jl")
include("prop_get_beamradius.jl")
include("prop_get_distancetofocus.jl")
include("prop_get_fratio.jl")
include("prop_get_gridsize.jl")
include("prop_get_nyquistsampling.jl")
include("prop_get_phase.jl")
include("prop_get_refradius.jl")
include("prop_get_sampling.jl")
include("prop_get_sampling_arcsec.jl")
include("prop_get_sampling_radians.jl")
include("prop_get_wavefront.jl")
include("prop_get_wavelength.jl")
include("prop_get_z.jl")
include("prop_hex_wavefront.jl")
include("prop_hex_zernikes.jl")
include("prop_init_savestate.jl")
include("prop_irregular_polygon.jl")
include("prop_is_statesaved.jl")
include("prop_lens.jl")
include("prop_load_fftw_wisdom.jl")
include("prop_magnify.jl")
include("prop_multiply.jl")
include("prop_noll_zernikes.jl")
include("prop_pixellate.jl")
include("prop_polygon.jl")
include("prop_print_zernikes.jl")
include("prop_propagate.jl")
include("prop_psd_errormap.jl")
include("prop_ptp.jl")
include("prop_qphase.jl")
include("prop_radius.jl")
include("prop_readmap.jl")
include("prop_rectangle.jl")
include("prop_rectangular_aperture.jl")
include("prop_rectangular_obscuration.jl")
include("prop_resamplemap.jl")
include("prop_rotate.jl")
include("prop_rounded_rectangle.jl")
include("prop_run.jl")
include("prop_run_multi.jl")
include("prop_savestate.jl")
include("prop_select_propagator.jl")
include("prop_set_antialiasing.jl")
include("prop_shift_center.jl")
include("prop_sinc.jl")
include("prop_state.jl")
include("prop_stw.jl")
include("prop_szoom.jl")
include("prop_table.jl")
include("prop_use_ffti.jl")
include("prop_use_fftw.jl")
include("prop_wavefront.jl")
include("prop_writemap.jl")
include("prop_wts.jl")
include("prop_zernikes.jl")
include("switch_set.jl")

export CompatPolicy, Python334Policy, CorrectedPolicy, resolve_compat_policy
export BackendStyle, CPUBackend, UnknownBackend, FFTStyle, FFTWStyle, GenericFFTStyle
export InterpStyle, GenericInterpStyle, CubicInterpStyle, RNGStyle, GenericRNGStyle
export RunContext, ProperConfig, ProperRuntime, WaveFront
export prop_wavefront, prop_begin, prop_end, prop_add_phase, prop_add_wavefront
export prop_multiply, prop_divide, prop_shift_center
export prop_radius, prop_ellipse, prop_rectangle
export prop_circular_aperture, prop_circular_obscuration
export prop_elliptical_aperture, prop_elliptical_obscuration
export prop_rectangular_aperture, prop_rectangular_obscuration
export prop_qphase, prop_lens, prop_ptp, prop_wts, prop_stw, prop_propagate
export prop_magnify, prop_rotate, prop_resamplemap, prop_select_propagator
export prop_fits_read, prop_fits_write, prop_readmap, prop_writemap
export prop_errormap, prop_psd_errormap
export prop_savestate, prop_state, prop_is_statesaved, prop_init_savestate, prop_end_savestate
export prop_execute_multi, prop_use_fftw, prop_use_ffti, prop_fftw_wisdom, prop_load_fftw_wisdom
export prop_set_antialiasing, prop_table
export prop_dm, prop_sinc, prop_noll_zernikes, prop_zernikes, prop_fit_zernikes
export prop_fit_dm, prop_pixellate, prop_polygon, prop_irregular_polygon, prop_rounded_rectangle
export prop_compile_c, prop_dftidefs, prop_ffti, prop_fftw, prop_szoom
export DftiErrorMessage
export libcconv, libcconvthread, libszoom, prop_cubic_conv, prop_define_entrance
export prop_hex_wavefront, prop_hex_zernikes, prop_8th_order_mask
export prop_get_wavefront, prop_get_amplitude, prop_get_phase
export prop_get_sampling, prop_get_sampling_radians, prop_get_sampling_arcsec
export prop_get_wavelength, prop_get_z, prop_get_gridsize
export prop_get_beamradius, prop_get_distancetofocus, prop_get_fratio
export prop_get_nyquistsampling, prop_get_refradius
export prop_run, prop_run_multi
export switch_set

end # module proper
