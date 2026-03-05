@testset "Phase 1 foundation" begin
    ctx_default = RunContext()
    @test ctx_default.policy isa Python334Policy

    ctx_corrected = RunContext(; compat_mode=:corrected)
    @test ctx_corrected.policy isa CorrectedPolicy

    @test_throws ArgumentError RunContext(; compat_mode=:unknown)

    @test proper.interp_style(Matrix{Float64}) isa CubicInterpStyle
    @test proper.interp_style(AbstractMatrix{Float64}) isa GenericInterpStyle
    ctx_f32 = RunContext(Matrix{Float32}; compat_mode=:python334)
    @test ctx_f32.interp isa CubicInterpStyle

    wf = prop_begin(2.0, 550e-9, 16)
    @test wf isa WaveFront
    @test size(wf.field) == (16, 16)

    prop_add_phase(wf, zeros(16, 16))
    out_i, s_i = prop_end(wf)
    @test eltype(out_i) <: Real
    @test s_i == wf.sampling_m

    out_c, s_c = prop_end(wf; noabs=true)
    @test eltype(out_c) <: Complex
    @test s_c == wf.sampling_m

    @test prop_shift_center(reshape(1:16, 4, 4)) == [11 15 3 7; 12 16 4 8; 9 13 1 5; 10 14 2 6]
end
