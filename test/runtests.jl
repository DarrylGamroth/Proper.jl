using Test
using Proper

include("example_loader.jl")
include("test_phase1_foundation.jl")
include("test_file_mapping.jl")
include("test_phase2_core.jl")
include("test_inference_allocations.jl")
include("test_r2_trait_routing.jl")
include("test_r3_mutating_workspace.jl")
include("test_r4_state_typing.jl")
include("test_r5_performance_gates.jl")
include("test_examples_smoke.jl")
include("test_phase9_semantic_reconciliation.jl")
