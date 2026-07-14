using Documenter
using Proper

DocMeta.setdocmeta!(Proper, :DocTestSetup, :(using Proper); recursive=true)

makedocs(
    modules=[Proper],
    sitename="Proper.jl",
    root=normpath(joinpath(@__DIR__, "..")),
    source="docs",
    build=joinpath("build", "documentation"),
    clean=true,
    doctest=true,
    checkdocs=:none,
    pagesonly=true,
    warnonly=false,
    format=Documenter.HTML(prettyurls=false, size_threshold_warn=150 * 2^10),
    pages=[
        "Home" => "index.md",
        "Guides" => [
            "Documentation Map" => "README.md",
            "API Examples" => "API_EXAMPLES.md",
            "Migration" => "MIGRATION_GUIDE.md",
            "Prescription Authoring" => "PRESCRIPTION_AUTHORING_GUIDE.md",
            "Prepared Execution" => "PREPARED_EXECUTION_GUIDE.md",
            "Latency Benchmarking" => "LATENCY_BENCHMARKING.md",
            "Adaptive Optics Integration" => "ADAPTIVE_OPTICS_INTEGRATION.md",
        ],
        "Contracts" => [
            "API" => "api_contract.md",
            "Numerics" => "numerics_contract.md",
            "Backend Traits" => "backend_traits.md",
            "Parity Harness" => "parity_harness_contract.md",
            "Parity Thresholds" => "parity_thresholds.md",
            "Compatibility Decisions" => "compat_decisions.md",
        ],
        "Maintainer Reference" => [
            "Porting Checklist" => "PORTING_CHECKLIST.md",
            "GPU Implementation" => "GPU_IMPLEMENTATION_PLAN.md",
            "WFIRST Phase B Matrix" => "WFIRST_PHASEB_CONFIG_MATRIX.md",
            "Historical Parity Closure" => "PARITY_CLOSURE.md",
            "Historical Semantic Reconciliation" => "SEMANTIC_RECONCILIATION.md",
        ],
    ],
)
