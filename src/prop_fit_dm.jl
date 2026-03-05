"""Placeholder DM fit returning input map and zero residual."""
function prop_fit_dm(dm_map::AbstractMatrix; kwargs...)
    return copy(dm_map), zero(eltype(dm_map))
end

@inline function _reflect_index(i::Int, n::Int)
    @inbounds while i < 1 || i > n
        if i < 1
            i = 1 - i
        else
            i = 2n - i + 1
        end
    end
    return i
end

@inline function _convolve_reflect(a::AbstractMatrix{T}, k::AbstractMatrix{T}) where {T<:AbstractFloat}
    ny, nx = size(a)
    ky, kx = size(k)
    cy = ky ÷ 2 + 1
    cx = kx ÷ 2 + 1
    out = similar(a)
    @inbounds for y in 1:ny
        for x in 1:nx
            acc = zero(T)
            for j in 1:ky
                yy = _reflect_index(y + (j - cy), ny)
                for i in 1:kx
                    xx = _reflect_index(x + (i - cx), nx)
                    acc += a[yy, xx] * k[j, i]
                end
            end
            out[y, x] = acc
        end
    end
    return out
end

"""Determine DM actuator pistons to fit desired surface with influence kernel."""
function prop_fit_dm(dm_map::AbstractMatrix, inf_kernel::AbstractMatrix; kwargs...)
    T = float(promote_type(eltype(dm_map), eltype(inf_kernel)))
    dm_target = T.(dm_map)
    ker = T.(inf_kernel)

    e0 = T(1e6)
    dm_cmd = copy(dm_target)
    last_good = copy(dm_cmd)
    dm_surface = _convolve_reflect(dm_cmd, ker)
    diff = dm_target .- dm_surface
    e = sqrt(sum(abs2, diff))

    # Follow upstream convergence rule.
    while e < e0 && (e > zero(T)) && ((e0 - e) / e > T(0.01))
        last_good .= dm_cmd
        e0 = e
        dm_cmd .+= diff
        dm_surface = _convolve_reflect(dm_cmd, ker)
        diff = dm_target .- dm_surface
        maximum(abs, diff) < T(1e-15) && break
        e = sqrt(sum(abs2, diff))
    end

    if e > e0
        dm_cmd .= last_good
        dm_surface = _convolve_reflect(dm_cmd, ker)
    end

    return dm_cmd, dm_surface
end
