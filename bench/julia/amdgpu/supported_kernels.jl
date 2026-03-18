using BenchmarkTools
using JSON3
using Proper
using AMDGPU
include(joinpath(@__DIR__, "..", "..", "common", "metadata.jl"))
include(joinpath(@__DIR__, "..", "..", "common", "amdgpu_support.jl"))
include(joinpath(@__DIR__, "..", "..", "common", "wavefront_state.jl"))
include(joinpath(@__DIR__, "..", "..", "common", "amdgpu_wavefront_kernel_cases.jl"))
using .BenchMetadata

const RUN_TAG = "steady_state_amdgpu_supported_kernels"
const REPORT_PATH = joinpath(@__DIR__, "..", "..", "reports", "amdgpu_supported_kernels.json")

AMDGPU.allowscalar(false)

function bench_amdgpu_supported_kernels()
    nprop = 512
    nmap = 256
    samples = 20
    wavefront_stats = benchmark_amdgpu_wavefront_kernel_stats(Float64; grid_n=nprop, samples=samples)

    img = amdgpu_rand(Float32, nmap, nmap)
    ctx_img = RunContext(typeof(img))
    rot_out = similar(img)
    mag_out = similar(img)
    szoom_out = similar(img)
    pix_out = similar(img, nmap ÷ 2, nmap ÷ 2)

    wf_map = amdgpu_wavefront_begin(1.0, 0.55e-6, nmap)
    dmap = amdgpu_rand(Float32, nmap, nmap)
    ctx_map = RunContext(wf_map)
    res_out = amdgpu_zeros(Float64, nmap, nmap)
    res_opts = Proper.ResampleMapOptions(wf_map, wf_map.sampling_m, nmap / 2, nmap / 2)
    rect_out = amdgpu_zeros(Float64, nmap, nmap)
    round_out = similar(rect_out)

    prop_rotate!(rot_out, img, 12.0, ctx_img)
    prop_magnify!(mag_out, img, 1.1, ctx_img; QUICK=true)
    prop_szoom!(szoom_out, img, 1.1, ctx_img)
    Proper._prop_pixellate_factor!(pix_out, img, 2)
    prop_resamplemap!(res_out, wf_map, dmap, res_opts, ctx_map)
    prop_rectangle!(rect_out, wf_map, 0.4, 0.2, 0.03, -0.05; ROTATION=22.0, NORM=true)
    prop_rounded_rectangle!(round_out, wf_map, 0.05, 0.3, 0.2, 0.01, -0.02)
    prop_rotate!(rot_out, img, 12.0, ctx_img)
    prop_magnify!(mag_out, img, 1.1, ctx_img; QUICK=true)
    prop_szoom!(szoom_out, img, 1.1, ctx_img)
    Proper._prop_pixellate_factor!(pix_out, img, 2)
    prop_resamplemap!(res_out, wf_map, dmap, res_opts, ctx_map)
    prop_rectangle!(rect_out, wf_map, 0.4, 0.2, 0.03, -0.05; ROTATION=22.0, NORM=true)
    prop_rounded_rectangle!(round_out, wf_map, 0.05, 0.3, 0.2, 0.01, -0.02)
    amdgpu_sync()

    r = run(@benchmarkable begin
        prop_rotate!($rot_out, $img, 12.0, $ctx_img)
        amdgpu_sync()
    end evals=1 samples=samples)

    m = run(@benchmarkable begin
        prop_magnify!($mag_out, $img, 1.1, $ctx_img; QUICK=true)
        amdgpu_sync()
    end evals=1 samples=samples)

    sz = run(@benchmarkable begin
        prop_szoom!($szoom_out, $img, 1.1, $ctx_img)
        amdgpu_sync()
    end evals=1 samples=samples)

    px = run(@benchmarkable begin
        Proper._prop_pixellate_factor!($pix_out, $img, 2)
        amdgpu_sync()
    end evals=1 samples=samples)

    rs = run(@benchmarkable begin
        prop_resamplemap!($res_out, $wf_map, $dmap, $res_opts, $ctx_map)
        amdgpu_sync()
    end evals=1 samples=samples)

    rc = run(@benchmarkable begin
        prop_rectangle!($rect_out, $wf_map, 0.4, 0.2, 0.03, -0.05; ROTATION=22.0, NORM=true)
        amdgpu_sync()
    end evals=1 samples=samples)

    rr = run(@benchmarkable begin
        prop_rounded_rectangle!($round_out, $wf_map, 0.05, 0.3, 0.2, 0.01, -0.02)
        amdgpu_sync()
    end evals=1 samples=samples)

    return Dict(
        "meta" => merge(amdgpu_report_meta(RUN_TAG; device=amdgpu_device_label()), Dict("prop_grid_n" => nprop, "map_grid_n" => nmap)),
        "policy" => "steady-state supported AMDGPU kernel timings with per-sample wavefront state restore; TTFx excluded; per-sample synchronization included",
        "kernels" => Dict(
            "prop_qphase" => wavefront_stats["prop_qphase"],
            "prop_ptp" => wavefront_stats["prop_ptp"],
            "prop_wts" => wavefront_stats["prop_wts"],
            "prop_stw" => wavefront_stats["prop_stw"],
            "prop_circular_aperture" => wavefront_stats["prop_circular_aperture"],
            "prop_end_mutating" => wavefront_stats["prop_end_mutating"],
            "prop_rotate_mutating" => trial_stats(r),
            "prop_magnify_quick_mutating" => trial_stats(m),
            "prop_szoom_mutating" => trial_stats(sz),
            "prop_pixellate_mutating" => trial_stats(px),
            "prop_resamplemap_mutating" => trial_stats(rs),
            "prop_rectangle_mutating" => trial_stats(rc),
            "prop_rounded_rectangle_mutating" => trial_stats(rr),
        ),
    )
end

write_benchmark_report(REPORT_PATH, bench_amdgpu_supported_kernels())
