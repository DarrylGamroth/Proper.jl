@testset "Entrance pupil normalization" begin
    for T in (Float32, Float64)
        input = reshape(
            Complex{T}.(1:16) .+ im .* reverse(Complex{T}.(1:16)),
            4,
            4,
        )
        wf = Proper.WaveFront(copy(input), T(500e-9), T(1e-3), zero(T), one(T))
        expected = input ./ sqrt(sum(abs2, input))

        @test (@inferred prop_define_entrance(wf)) === wf
        @test eltype(wf.field) === Complex{T}
        @test isapprox(wf.field, expected; atol=T(8) * eps(T), rtol=T(8) * eps(T))
        @test isapprox(sum(abs2, wf.field), one(T); atol=T(16) * eps(T), rtol=T(16) * eps(T))
    end

    for invalid in (
        zeros(ComplexF64, 4, 4),
        fill(ComplexF64(NaN, 0), 4, 4),
        fill(ComplexF64(Inf, 0), 4, 4),
    )
        wf = Proper.WaveFront(copy(invalid), 500e-9, 1e-3, 0.0, 1.0)
        field_before = copy(wf.field)
        @test_throws DomainError prop_define_entrance(wf)
        @test isequal(wf.field, field_before)
    end
end
