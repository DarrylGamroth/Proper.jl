@noinline function _not_implemented(fname::Symbol)
    throw(ErrorException("$(fname) is not implemented yet in Proper.jl"))
end
