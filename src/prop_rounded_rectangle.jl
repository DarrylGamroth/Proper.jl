"""Return rectangular mask with rounded corners."""
function _prop_rounded_rectangle!(
    ::GeometryLoopExecStyle,
    image::AbstractMatrix,
    wf::WaveFront,
    corner_radius::Real,
    width::Real,
    height::Real,
    xc::Real,
    yc::Real,
)
    ny, nx = size(wf.field)
    size(image) == (ny, nx) || throw(ArgumentError("output size must match wavefront"))
    RT = eltype(image)
    dx = RT(wf.sampling_m)
    x = coordinate_axis(nx, dx) .- RT(xc)
    y = coordinate_axis(ny, dx) .- RT(yc)

    r = max(RT(corner_radius), zero(RT))
    hw = RT(width) / RT(2)
    hh = RT(height) / RT(2)
    one_rt = one(RT)
    zero_rt = zero(RT)

    @inbounds for j in 1:nx
        qx = abs(x[j]) - (hw - r)
        for i in 1:ny
            qy = abs(y[i]) - (hh - r)
            ax = max(qx, zero_rt)
            ay = max(qy, zero_rt)
            inside = (qx <= 0 && abs(y[i]) <= hh) || (qy <= 0 && abs(x[j]) <= hw) || (ax * ax + ay * ay <= r * r)
            image[i, j] = inside ? one_rt : zero_rt
        end
    end

    return image
end

function _prop_rounded_rectangle!(
    ::GeometryKAExecStyle,
    image::AbstractMatrix,
    wf::WaveFront,
    corner_radius::Real,
    width::Real,
    height::Real,
    xc::Real,
    yc::Real,
)
    size(image) == size(wf.field) || throw(ArgumentError("output size must match wavefront"))
    RT = eltype(image)
    dx = RT(wf.sampling_m)
    r = max(RT(corner_radius), zero(RT))
    hw = RT(width) / RT(2)
    hh = RT(height) / RT(2)
    return ka_rounded_rectangle_mask!(image, dx, RT(xc), RT(yc), r, hw, hh)
end

function _prop_rounded_rectangle!(
    image::AbstractMatrix,
    wf::WaveFront,
    corner_radius::Real,
    width::Real,
    height::Real,
    xc::Real,
    yc::Real,
)
    return _prop_rounded_rectangle!(geometry_exec_style(typeof(image), size(image, 1), size(image, 2)), image, wf, corner_radius, width, height, xc, yc)
end

function prop_rounded_rectangle!(
    image::AbstractMatrix,
    wf::WaveFront,
    corner_radius::Real,
    width::Real,
    height::Real,
    xc::Real=0.0,
    yc::Real=0.0,
)
    return _prop_rounded_rectangle!(image, wf, corner_radius, width, height, xc, yc)
end

function prop_rounded_rectangle(
    wf::WaveFront,
    corner_radius::Real,
    width::Real,
    height::Real,
    xc::Real=0.0,
    yc::Real=0.0,
)
    RT = real(eltype(wf.field))
    image = similar(wf.field, RT, size(wf.field)...)
    return _prop_rounded_rectangle!(image, wf, corner_radius, width, height, xc, yc)
end
