function passvalue_kwargs(::Nothing)
    return NamedTuple()
end

function passvalue_kwargs(passvalue::NamedTuple{K}) where {K}
    keys_tuple = Tuple(Symbol(lowercase(String(k))) for k in K)
    return NamedTuple{keys_tuple}(values(passvalue))
end

function passvalue_kwargs(passvalue::AbstractDict)
    keys_tuple = Tuple(Symbol(lowercase(String(k))) for k in keys(passvalue))
    values_tuple = Tuple(values(passvalue))
    return NamedTuple{keys_tuple}(values_tuple)
end

function passvalue_kwargs(passvalue)
    throw(ArgumentError("PASSVALUE must be a NamedTuple, dictionary, or nothing for this example; got $(typeof(passvalue))"))
end
