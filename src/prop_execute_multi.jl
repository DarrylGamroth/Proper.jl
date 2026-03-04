"""Compatibility alias for `prop_run_multi`."""
function prop_execute_multi(args...; kwargs...)
    return prop_run_multi(args...; kwargs...)
end
