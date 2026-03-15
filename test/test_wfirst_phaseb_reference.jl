using Test
using FFTW
using Proper.WFIRSTPhaseBProper

@testset "WFIRST Phase B reference helpers" begin
    old_root = data_dir()
    try
        set_data_dir!(joinpath(pwd(), ".cache", "wfirst", "dummy"))
        @test phaseb_default_data_root() == data_dir()
    finally
        set_data_dir!(old_root)
    end

    a = reshape(collect(1.0:16.0), 4, 4)
    @test trim(a, 2) == [6.0 10.0; 7.0 11.0]
    @test size(trim(a, 6)) == (6, 6)

    c = ComplexF64.(reshape(1:16, 4, 4))
    @test ffts(copy(c), -1) ≈ circshift(fft(circshift(c, (-2, -2))) ./ length(c), (2, 2))

    m = mft2(c, 0.1, 2.0, 4, -1)
    @test size(m) == (4, 4)

    cases = phaseb_case_definitions()
    @test Set(keys(cases)) == Set(["compact_hlc", "full_hlc"])
    @test cases["compact_hlc"].func === wfirst_phaseb_compact
    @test cases["full_hlc"].func === wfirst_phaseb
end
