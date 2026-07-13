using Aqua
using Test
using Proper

@testset "Aqua" begin
    Aqua.test_all(
        Proper;
        project_extras=true,
        deps_compat=true,
    )
end
