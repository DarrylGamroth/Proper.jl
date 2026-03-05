"""Return rectangular mask with rounded corners."""
function prop_rounded_rectangle(wf::WaveFront, corner_radius::Real, width::Real, height::Real, xc::Real=0.0, yc::Real=0.0)
    ny, nx = size(wf.field)
    dx = wf.sampling_m
    x = coordinate_axis(nx, dx) .- float(xc)
    y = coordinate_axis(ny, dx) .- float(yc)

    r = max(float(corner_radius), 0.0)
    hw = float(width) / 2
    hh = float(height) / 2

    m = zeros(Float64, ny, nx)
    @inbounds for j in 1:nx
        qx = abs(x[j]) - (hw - r)
        for i in 1:ny
            qy = abs(y[i]) - (hh - r)
            ax = max(qx, 0.0)
            ay = max(qy, 0.0)
            inside = (qx <= 0 && abs(y[i]) <= hh) || (qy <= 0 && abs(x[j]) <= hw) || (ax * ax + ay * ay <= r * r)
            m[i, j] = inside ? 1.0 : 0.0
        end
    end

    return m
end
