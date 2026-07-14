module ProperExamples

"""Return the source path for a ported upstream PROPER example."""
function example_path(name::Union{AbstractString,Symbol})
    stem = splitext(String(name))[1]
    return joinpath(@__DIR__, string(stem, ".jl"))
end

end
