struct LensOptions
    surface_name::String
end

@inline LensOptions(surface_name::AbstractString="") = LensOptions(String(surface_name))

@inline _lens_phase(::Val{INSIDE_TO_INSIDE}, lens_fl::Real, r_beam_old::Real, r_beam::Real) = inv(lens_fl)
@inline _lens_phase(::Val{INSIDE_TO_OUTSIDE}, lens_fl::Real, r_beam_old::Real, r_beam::Real) = inv(lens_fl) + inv(r_beam)
@inline _lens_phase(::Val{OUTSIDE_TO_INSIDE}, lens_fl::Real, r_beam_old::Real, r_beam::Real) = inv(lens_fl) - inv(r_beam_old)

@inline function _lens_phase(::Val{OUTSIDE_TO_OUTSIDE}, lens_fl::Real, r_beam_old::Real, r_beam::Real)
    if r_beam_old == 0
        return inv(lens_fl) + inv(r_beam)
    elseif r_beam == 0
        return inv(lens_fl) - inv(r_beam_old)
    end
    return inv(lens_fl) - inv(r_beam_old) + inv(r_beam)
end

@inline function _lens_phase(pt::PropagatorType, lens_fl::Real, r_beam_old::Real, r_beam::Real)
    return _lens_phase(Val(pt), lens_fl, r_beam_old, r_beam)
end

@inline function _prop_lens!(wf::WaveFront, lens_fl::Real, opts::LensOptions)
    _ = opts
    iszero(lens_fl) && throw(ArgumentError("lens_fl must be non-zero"))

    λ = wf.wavelength_m
    z = wf.z_m
    z_w0 = wf.z_w0_m

    wf.z_rayleigh_m = pi * wf.w0_m^2 / λ
    w_at_surface = wf.w0_m * sqrt(1 + ((z - z_w0) / wf.z_rayleigh_m)^2)

    gR_beam_inf = false
    gR_beam = zero(λ)

    if (z - z_w0) != 0
        gR_beam_old = (z - z_w0) + wf.z_rayleigh_m^2 / (z - z_w0)
        if gR_beam_old != lens_fl
            gR_beam = inv(inv(gR_beam_old) - inv(float(lens_fl)))
        else
            gR_beam_inf = true
        end
    else
        gR_beam = -float(lens_fl)
    end

    r_beam_old = (wf.beam_type_old === INSIDE || wf.reference_surface === PLANAR) ? zero(λ) : (z - z_w0)

    if !gR_beam_inf
        wf.z_w0_m = -gR_beam / (1 + (λ * gR_beam / (pi * w_at_surface^2))^2) + z
        wf.w0_m = w_at_surface / sqrt(1 + (pi * w_at_surface^2 / (λ * gR_beam))^2)
    else
        wf.z_w0_m = z
        wf.w0_m = w_at_surface
    end

    wf.z_rayleigh_m = pi * wf.w0_m^2 / λ
    beam_type_new = abs(wf.z_w0_m - z) < wf.rayleigh_factor * wf.z_rayleigh_m ? INSIDE : OUTSIDE
    r_beam = beam_type_new === INSIDE ? zero(λ) : (z - wf.z_w0_m)

    wf.propagator_type = propagator_transition(wf.beam_type_old, beam_type_new)
    lens_phase = _lens_phase(wf.propagator_type, float(lens_fl), r_beam_old, r_beam)

    rho = prop_radius(wf)
    prop_add_phase(wf, @. -(rho^2) * (lens_phase / 2))

    wf.reference_surface = beam_type_new === INSIDE ? PLANAR : SPHERICAL
    wf.beam_type_old = beam_type_new
    wf.current_fratio = abs(wf.z_w0_m - z) / (2 * w_at_surface)

    return wf
end

"""Alter the current wavefront as a perfect thin lens would."""
function prop_lens(wf::WaveFront, lens_fl::Real, surface_name::AbstractString="")
    return _prop_lens!(wf, lens_fl, LensOptions(surface_name))
end
