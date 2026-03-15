using Test
using FFTW
using Statistics
using Proper
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

    mktempdir() do d
        old_root = data_dir()
        try
            set_data_dir(d)
            @test data_dir() == abspath(d)

            copied = joinpath(d, "copied")
            copy_here(copied)
            @test isfile(joinpath(copied, "wfirst_phaseb.jl"))
            @test isfile(joinpath(copied, "wfirst_phaseb_compact.jl"))

            copied_examples = joinpath(d, "examples")
            copy_examples_here(copied_examples)
            @test isfile(joinpath(copied_examples, "wfirst_phaseb_reference.jl"))

            polroot = joinpath(d, "coeffs")
            zamp = zeros(Float64, 2, 2, 6, 22)
            zpha = zeros(Float64, 2, 2, 6, 22)
            zamp[1, 1, :, 1] .= 2.0
            prop_fits_write(polroot * "_amp.fits", zamp)
            prop_fits_write(polroot * "_pha.fits", zpha)

            amp, pha = polab(polroot, 550e-9, 6.0, -1)
            @test size(amp) == size(pha)
            @test all(iszero, pha)
            @test isapprox(mean(amp), 2.0; atol=1e-12)

            wf = prop_begin(1.0, 550e-9, 8)
            field_before = copy(wf.field)
            polmap(wf, polroot, 6.0, -1)
            @test wf.field ≈ field_before .* 2
        finally
            set_data_dir!(old_root)
        end
    end
end
