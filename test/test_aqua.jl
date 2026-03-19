using Aqua
using Test
using Proper

@testset "Aqua" begin
    Aqua.test_all(
        Proper;
        ambiguities=false,
        persistent_tasks=false,
        project_extras=false,
        stale_deps=false,
        deps_compat=false,
    )
end
