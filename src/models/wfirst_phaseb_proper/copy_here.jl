function copy_here(dest::AbstractString=pwd())
    mkpath(dest)
    for file in ("wfirst_phaseb.jl", "wfirst_phaseb_compact.jl")
        cp(joinpath(@__DIR__, file), joinpath(dest, file); force=true)
    end
    return dest
end
