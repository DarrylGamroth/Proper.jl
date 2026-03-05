using Test

@testset "R4 state typing and transition selectors" begin
    wf = prop_begin(1.0, 500e-9, 32)

    @test wf.reference_surface isa Proper.ReferenceSurface
    @test wf.beam_type_old isa Proper.BeamType
    @test wf.propagator_type isa Proper.PropagatorType

    # Compatibility-facing symbol comparisons/assignments remain usable.
    @test wf.reference_surface == :PLANAR
    wf.reference_surface = :SPHERI
    wf.beam_type_old = :OUTSIDE
    wf.propagator_type = :OUTSIDE_to_OUTSIDE
    @test wf.reference_surface === Proper.SPHERI
    @test wf.beam_type_old === Proper.OUTSIDE
    @test wf.propagator_type === Proper.OUTSIDE_to_OUTSIDE

    @test @inferred(Proper.propagator_transition(Proper.INSIDE_, Proper.INSIDE_)) === Proper.INSIDE__to_INSIDE_
    @test @inferred(Proper.propagator_transition(Proper.INSIDE_, Proper.OUTSIDE)) === Proper.INSIDE__to_OUTSIDE
    @test @inferred(Proper.propagator_transition(Proper.OUTSIDE, Proper.INSIDE_)) === Proper.OUTSIDE_to_INSIDE_
    @test @inferred(Proper.propagator_transition(Proper.OUTSIDE, Proper.OUTSIDE)) === Proper.OUTSIDE_to_OUTSIDE

    @test @inferred(Proper.to_plane_propagator(Proper.INSIDE__to_OUTSIDE)) === Proper.INSIDE__to_INSIDE_
    @test @inferred(Proper.to_plane_propagator(Proper.OUTSIDE_to_OUTSIDE)) === Proper.OUTSIDE_to_INSIDE_

    wf2 = prop_begin(1.0, 500e-9, 32)
    dzw = @inferred prop_select_propagator(wf2, 0.1)
    @test dzw isa Float64
    @test wf2.propagator_type isa Proper.PropagatorType
end
