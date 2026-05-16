include(joinpath(@__DIR__, "..", "bench", "cards", "card00", "benchmark_helpers.jl"))
using .Card00BenchmarkHelpers

@testset "Card 00 benchmark helpers" begin
    expected_sides = Dict(277 => 17, 468 => 22, 1024 => 32, 4096 => 64)
    for active_count in CARD00_DM_ACTIVE_COUNTS
        cmd = active_dm_command(active_count)
        @test cmd.grid_side == expected_sides[active_count]
        @test size(cmd.values) == (cmd.grid_side, cmd.grid_side)
        @test length(cmd.active_positions) == active_count
        @test length(unique(cmd.active_positions)) == active_count
        @test count_active_actuators(cmd.values) == active_count
        @test all(!iszero(cmd.values[pos]) for pos in cmd.active_positions)
        @test count(iszero, cmd.values) == cmd.grid_side^2 - active_count
        @test cmd.label == "$(active_count) on $(cmd.grid_side)x$(cmd.grid_side)"

        center = (cmd.grid_side + 1) / 2
        active_dist2 = Set((pos[1] - center)^2 + (pos[2] - center)^2 for pos in cmd.active_positions)
        inactive_dist2 = [
            (i - center)^2 + (j - center)^2
            for i in 1:cmd.grid_side for j in 1:cmd.grid_side
            if cmd.values[i, j] == 0
        ]
        isempty(inactive_dist2) || @test maximum(active_dist2) <= minimum(inactive_dist2)
    end

    side, positions = active_actuator_positions(5)
    @test side == 3
    @test positions[1] == CartesianIndex(2, 2)
    @test positions == active_actuator_positions(5)[2]

    pupil = circular_fit_pupil(32)
    @test size(pupil.mask) == (32, 32)
    @test pupil.radius == 15.0
    @test pupil.mask[pupil.yc + 1, pupil.xc + 1] == 1.0
end
