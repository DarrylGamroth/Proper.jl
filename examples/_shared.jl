function center_crop(a::AbstractMatrix, n::Integer)
    ny, nx = size(a)
    cy = ny ÷ 2 + 1
    cx = nx ÷ 2 + 1
    ry = (cy - n ÷ 2):(cy + (n - 1) ÷ 2)
    rx = (cx - n ÷ 2):(cx + (n - 1) ÷ 2)
    return @view a[ry, rx]
end
