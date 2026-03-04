using Test

@testset "One-to-one mapping from Python" begin
    py_root = normpath(joinpath(@__DIR__, "..", "..", "proper_v3.3.4_python", "proper"))
    @test isdir(py_root)

    py_files = sort(filter(f -> endswith(f, ".py"), readdir(py_root)))
    jl_files = Set(filter(f -> endswith(f, ".jl"), readdir(joinpath(@__DIR__, "..", "src"))))

    for py in py_files
        py == "__init__.py" && continue
        expected = replace(py, ".py" => ".jl")
        @test expected in jl_files
    end
end

@testset "Example mapping from Python" begin
    py_examples = normpath(joinpath(@__DIR__, "..", "..", "proper_v3.3.4_python", "proper", "examples"))
    @test isdir(py_examples)

    py_files = sort(filter(f -> endswith(f, ".py"), readdir(py_examples)))
    jl_files = Set(filter(f -> endswith(f, ".jl"), readdir(joinpath(@__DIR__, "..", "examples"))))

    for py in py_files
        py == "__init__.py" && continue
        expected = replace(py, ".py" => ".jl")
        @test expected in jl_files
    end
end
