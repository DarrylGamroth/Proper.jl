function load_example_module(path::AbstractString; name::Symbol=gensym(Symbol(splitext(basename(path))[1])))
    mod = Module(name)
    Core.eval(mod, :(include(path::AbstractString) = Base.include($mod, path)))
    Base.include(mod, path)
    return mod
end
