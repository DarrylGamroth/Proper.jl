"""Return MKL DFTI-compatible constant dictionary for compatibility shims."""
const _DFTIDEFS = Dict{Symbol,Int}(
    :DFTI_FORWARD_DOMAIN => 0,
    :DFTI_DIMENSION => 1,
    :DFTI_LENGTHS => 2,
    :DFTI_PRECISION => 3,
    :DFTI_FORWARD_SCALE => 4,
    :DFTI_BACKWARD_SCALE => 5,
    :DFTI_FORWARD_SIGN => 6,
    :DFTI_NUMBER_OF_TRANSFORMS => 7,
    :DFTI_COMPLEX_STORAGE => 8,
    :DFTI_REAL_STORAGE => 9,
    :DFTI_PLACEMENT => 11,
    :DFTI_THREAD_LIMIT => 27,
    :DFTI_NO_ERROR => 0,
    :DFTI_MEMORY_ERROR => 1,
    :DFTI_INVALID_CONFIGURATION => 2,
    :DFTI_INCONSISTENT_CONFIGURATION => 3,
)

prop_dftidefs() = copy(_DFTIDEFS)

function DftiErrorMessage(code::Integer)
    return get(
        Dict(
            0 => "No error",
            1 => "Memory error",
            2 => "Invalid configuration",
            3 => "Inconsistent configuration",
            4 => "Multithreaded error",
            5 => "Bad descriptor",
            6 => "Unimplemented",
            7 => "MKL internal error",
            8 => "Number of threads error",
            9 => "1D length exceeds int32",
        ),
        Int(code),
        "Unknown error code",
    )
end
