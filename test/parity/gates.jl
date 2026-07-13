function parity_threshold_failures(
    metrics::AbstractDict{<:AbstractString,<:Real},
    thresholds::AbstractDict{<:AbstractString,<:Real},
)
    failures = String[]
    for (metric, limit) in thresholds
        if !isfinite(limit) || limit < 0
            push!(failures, "$metric has invalid threshold $limit")
        elseif !haskey(metrics, metric)
            push!(failures, "missing metric $metric")
        else
            value = metrics[metric]
            if !isfinite(value)
                push!(failures, "$metric is not finite: $value")
            elseif value > limit
                push!(failures, "$metric=$value > $limit")
            end
        end
    end
    return failures
end

function parity_shape_failure(actual, expected)
    actual == expected && return nothing
    return "shape mismatch: got $actual, expected $expected"
end
