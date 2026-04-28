using Test

@testset "One-to-one mapping from Python" begin
    pyproper_root = get(ENV, "PYPROPER_ROOT", normpath(joinpath(@__DIR__, "..", "..", "proper_v3.3.4_python")))
    py_root = normpath(joinpath(pyproper_root, "proper"))
    if !isdir(py_root)
        @test_skip "Python PROPER baseline not available; set PYPROPER_ROOT to enable file mapping coverage"
    else
        py_files = sort(filter(f -> endswith(f, ".py"), readdir(py_root)))
        jl_files = Set(filter(f -> endswith(f, ".jl"), readdir(joinpath(@__DIR__, "..", "src"))))
        intentionally_unmapped = Set([
            "libcconvthread.jl",
            "libszoom.jl",
            "prop_compile_c.jl",
            "prop_execute_multi.jl",
            "prop_ffti.jl",
            "prop_fftw.jl",
            "prop_table.jl",
            "prop_use_ffti.jl",
            "prop_use_fftw.jl",
        ])

        for py in py_files
            py == "__init__.py" && continue
            expected = replace(py, ".py" => ".jl")
            expected in intentionally_unmapped && continue
            @test expected in jl_files
        end
    end
end

@testset "Example mapping from Python" begin
    pyproper_root = get(ENV, "PYPROPER_ROOT", normpath(joinpath(@__DIR__, "..", "..", "proper_v3.3.4_python")))
    py_examples = normpath(joinpath(pyproper_root, "proper", "examples"))
    if !isdir(py_examples)
        @test_skip "Python PROPER example baseline not available; set PYPROPER_ROOT to enable example mapping coverage"
    else
        py_files = sort(filter(f -> endswith(f, ".py"), readdir(py_examples)))
        jl_files = Set(filter(f -> endswith(f, ".jl"), readdir(joinpath(@__DIR__, "..", "examples"))))

        for py in py_files
            py == "__init__.py" && continue
            expected = replace(py, ".py" => ".jl")
            @test expected in jl_files
        end
    end
end
