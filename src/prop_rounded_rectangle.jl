"""Return rectangular mask with rounded corners."""
function prop_rounded_rectangle(wf::WaveFront, corner_radius::Real, width::Real, height::Real, xc::Real=0.0, yc::Real=0.0)
    ny, nx = size(wf.field)
    RT = float(promote_type(real(eltype(wf.field)), typeof(corner_radius), typeof(width), typeof(height), typeof(xc), typeof(yc), typeof(wf.sampling_m)))
    dx = RT(wf.sampling_m)
    x = coordinate_axis(nx, dx) .- RT(xc)
    y = coordinate_axis(ny, dx) .- RT(yc)

    r = max(RT(corner_radius), zero(RT))
    hw = RT(width) / RT(2)
    hh = RT(height) / RT(2)
    one_rt = one(RT)
    zero_rt = zero(RT)

    m = zeros(RT, ny, nx)
    @inbounds for j in 1:nx
        qx = abs(x[j]) - (hw - r)
        for i in 1:ny
            qy = abs(y[i]) - (hh - r)
            ax = max(qx, zero_rt)
            ay = max(qy, zero_rt)
            inside = (qx <= 0 && abs(y[i]) <= hh) || (qy <= 0 && abs(x[j]) <= hw) || (ax * ax + ay * ay <= r * r)
            m[i, j] = inside ? one_rt : zero_rt
        end
    end

    return m
end
