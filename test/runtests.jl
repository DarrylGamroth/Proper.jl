using Test
using Proper
using Random

const TEST_RNG = MersenneTwister(20260428)
Random.seed!(TEST_RNG, 20260428)

include("example_loader.jl")
include("test_phase1_foundation.jl")
include("test_file_mapping.jl")
include("test_ka_cpu_kernels.jl")
include("test_public_helper_coverage.jl")
include("test_phase2_core.jl")
include("test_inference_allocations.jl")
include("test_r2_trait_routing.jl")
include("test_r3_mutating_workspace.jl")
include("test_r4_state_typing.jl")
include("test_r5_performance_gates.jl")
include("test_examples_smoke.jl")
include("test_doc_examples.jl")
include("test_api_contract.jl")
include("test_wfirst_phaseb_reference.jl")
include("test_phase9_semantic_reconciliation.jl")
include("test_aqua.jl")
