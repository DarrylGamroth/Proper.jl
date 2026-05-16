module Card00BenchmarkHelpers

using Proper

export CARD00_DM_ACTIVE_COUNTS, CARD00_ZERNIKE_NTERMS, CARD00_ZERNIKE_GRIDS
export CARD00_ZERNIKE_FIT_GRIDS, CARD00_WFIRST_CASES
export active_actuator_positions, active_dm_command, count_active_actuators
export zernike_coefficients, circular_fit_pupil, zernike_fit_wavefront

const CARD00_DM_ACTIVE_COUNTS = (277, 468, 1024, 4096)
const CARD00_ZERNIKE_NTERMS = (8, 22, 64)
const CARD00_ZERNIKE_GRIDS = (128, 512, 1024)
const CARD00_ZERNIKE_FIT_GRIDS = (256,)
const CARD00_WFIRST_CASES = (
    "compact_hlc",
    "full_hlc",
    "compact_spc_spec_long",
    "full_spc_spec_long",
)

"""
    active_actuator_positions(active_count)

Return `(side, positions)` for a deterministic centered actuator mask. The
storage grid side is `ceil(Int, sqrt(active_count))`, and positions are the
`active_count` grid coordinates nearest the geometric center.
"""
function active_actuator_positions(active_count::Integer)
    nactive = Int(active_count)
    nactive > 0 || throw(ArgumentError("active_count must be positive"))
    side = ceil(Int, sqrt(nactive))
    center = (side + 1) / 2
    ranked = Vector{Tuple{Float64,Int,Int}}(undef, side * side)
    k = 0
    @inbounds for iy in 1:side
        dy = iy - center
        for ix in 1:side
            dx = ix - center
            k += 1
            ranked[k] = (dx * dx + dy * dy, iy, ix)
        end
    end
    sort!(ranked; by=x -> x)
    positions = Vector{CartesianIndex{2}}(undef, nactive)
    @inbounds for i in 1:nactive
        _, iy, ix = ranked[i]
        positions[i] = CartesianIndex(iy, ix)
    end
    return side, positions
end

"""
    active_dm_command(active_count; T=Float64, amplitude_m=25e-9)

Build the Card 00 deterministic DM command matrix. Inactive actuators are zero;
active actuator positions are selected by `active_actuator_positions`.
"""
function active_dm_command(active_count::Integer; T::Type{<:AbstractFloat}=Float64, amplitude_m::Real=25e-9)
    side, positions = active_actuator_positions(active_count)
    cmd = zeros(T, side, side)
    center = T((side + 1) / 2)
    denom = max(T(side - 1) / T(2), one(T))
    scale = T(amplitude_m)
    nactive = T(active_count)
    @inbounds for (k, pos) in enumerate(positions)
        y = (T(pos[1]) - center) / denom
        x = (T(pos[2]) - center) / denom
        cmd[pos] = scale * (one(T) + T(0.08) * x - T(0.05) * y + T(0.03) * sinpi(T(k) / nactive))
    end
    return (
        values=cmd,
        active_positions=positions,
        active_count=Int(active_count),
        grid_side=side,
        label="$(Int(active_count)) on $(side)x$(side)",
    )
end

count_active_actuators(command::AbstractMatrix) = count(!iszero, command)

function zernike_coefficients(nterms::Integer; T::Type{<:AbstractFloat}=Float64, scale_m::Real=10e-9)
    n = Int(nterms)
    n > 0 || throw(ArgumentError("nterms must be positive"))
    coeffs = Vector{T}(undef, n)
    base = T(scale_m)
    @inbounds for j in 1:n
        sgn = isodd(j) ? one(T) : -one(T)
        coeffs[j] = sgn * base * (one(T) + T(j) / T(n))
    end
    return coeffs
end

function circular_fit_pupil(grid_n::Integer; T::Type{<:AbstractFloat}=Float64)
    n = Int(grid_n)
    n > 2 || throw(ArgumentError("grid_n must be greater than 2"))
    radius = T(n ÷ 2 - 1)
    x0 = n ÷ 2
    y0 = n ÷ 2
    pupil = zeros(T, n, n)
    @inbounds for j in 1:n
        x = T(j - 1 - x0) / radius
        for i in 1:n
            y = T(i - 1 - y0) / radius
            pupil[i, j] = hypot(x, y) < one(T) ? one(T) : zero(T)
        end
    end
    return (mask=pupil, radius=radius, xc=x0, yc=y0)
end

function zernike_fit_wavefront(grid_n::Integer, nterms::Integer; T::Type{<:AbstractFloat}=Float64)
    wf = prop_begin(T(1.0), T(0.55e-6), Int(grid_n))
    nums = collect(1:Int(nterms))
    coeffs = zernike_coefficients(nterms; T=T, scale_m=7e-9)
    return prop_zernikes(wf, nums, coeffs; no_apply=true)
end

end
