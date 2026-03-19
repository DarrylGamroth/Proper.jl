function copy_examples_here(dest::AbstractString=pwd())
    mkpath(dest)
    for file in ("wfirst_phaseb_reference.jl",)
        cp(joinpath(dirname(dirname(@__DIR__)), "examples", file), joinpath(dest, file); force=true)
    end
    return dest
end
