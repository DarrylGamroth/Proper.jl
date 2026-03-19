@inline function trim(input_image::AbstractMatrix, output_dim::Integer)
    input_dim = size(input_image, 2)
    if input_dim == output_dim
        return input_image
    elseif output_dim < input_dim
        x1 = input_dim ÷ 2 - output_dim ÷ 2 + 1
        x2 = x1 + output_dim - 1
        return copy(@view input_image[x1:x2, x1:x2])
    end

    output_image = zeros(eltype(input_image), output_dim, output_dim)
    x1 = output_dim ÷ 2 - input_dim ÷ 2 + 1
    x2 = x1 + input_dim - 1
    @views output_image[x1:x2, x1:x2] .= input_image
    return output_image
end
