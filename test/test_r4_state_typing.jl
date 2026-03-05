using Test

@testset "R4 state typing and transition selectors" begin
    wf = prop_begin(1.0, 500e-9, 32)

    @test wf.reference_surface isa Proper.ReferenceSurface
    @test wf.beam_type_old isa Proper.BeamType
    @test wf.propagator_type isa Proper.PropagatorType

    @test wf.reference_surface === Proper.PLANAR
    wf.reference_surface = Proper.SPHERICAL
    wf.beam_type_old = Proper.OUTSIDE
    wf.propagator_type = Proper.OUTSIDE_TO_OUTSIDE
    @test wf.reference_surface === Proper.SPHERICAL
    @test wf.beam_type_old === Proper.OUTSIDE
    @test wf.propagator_type === Proper.OUTSIDE_TO_OUTSIDE

    @test @inferred(Proper.propagator_transition(Proper.INSIDE, Proper.INSIDE)) === Proper.INSIDE_TO_INSIDE
    @test @inferred(Proper.propagator_transition(Proper.INSIDE, Proper.OUTSIDE)) === Proper.INSIDE_TO_OUTSIDE
    @test @inferred(Proper.propagator_transition(Proper.OUTSIDE, Proper.INSIDE)) === Proper.OUTSIDE_TO_INSIDE
    @test @inferred(Proper.propagator_transition(Proper.OUTSIDE, Proper.OUTSIDE)) === Proper.OUTSIDE_TO_OUTSIDE

    @test @inferred(Proper.to_plane_propagator(Proper.INSIDE_TO_OUTSIDE)) === Proper.INSIDE_TO_INSIDE
    @test @inferred(Proper.to_plane_propagator(Proper.OUTSIDE_TO_OUTSIDE)) === Proper.OUTSIDE_TO_INSIDE

    wf2 = prop_begin(1.0, 500e-9, 32)
    dzw = @inferred prop_select_propagator(wf2, 0.1)
    @test dzw isa Float64
    @test wf2.propagator_type isa Proper.PropagatorType
end
