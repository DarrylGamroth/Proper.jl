# API Examples

These examples are written in `jldoctest` style and mirrored by
`test/test_doc_examples.jl` so they remain runnable in this repository even
without a full Documenter pipeline.

## `prop_run`

```jldoctest
julia> using Proper

julia> function demo_prescription(λm, n)
           wf = prop_begin(1.0, λm, n)
           prop_circular_aperture(wf, 0.5)
           prop_define_entrance(wf)
           return prop_end(wf)
       end
demo_prescription (generic function with 1 method)

julia> psf, sampling = prop_run(demo_prescription, 0.55, 32);

julia> size(psf), sampling > 0
((32, 32), true)
```

## `prepare_prescription`

```jldoctest
julia> using Proper

julia> demo_prescription(λm, n) = prop_end(prop_begin(1.0, λm, n));

julia> prepared = prepare_prescription(demo_prescription, 0.55, 16);

julia> psf, sampling = prop_run(prepared);

julia> size(psf), sampling > 0
((16, 16), true)
```

## `prepare_prescription_batch` and `prop_run_multi`

```jldoctest
julia> using Proper

julia> function batch_demo(λm, n, pass)
           wf = prop_begin(1.0 + pass, λm, n)
           return prop_end(wf)
       end
batch_demo (generic function with 1 method)

julia> prepared = prepare_prescription(batch_demo, 0.55, 16);

julia> batch = prepare_prescription_batch(prepared; pool_size=2);

julia> stack, samplings = prop_run_multi(batch; PASSVALUE=[0.0, 0.5]);

julia> size(stack), length(samplings)
((16, 16, 2), 2)
```

## `prepare_model`

```jldoctest
julia> using Proper

julia> function model_demo(λm, n, pass; gain=1.0)
           psf = fill(Float64(gain + pass), n, n)
           return psf, 1.0e-3
       end
model_demo (generic function with 1 method)

julia> model = prepare_model(:model_demo, model_demo, 0.55, 8; pool_size=2, assets=(gain=2.0,));

julia> psf, sampling = prop_run(model; slot=1, PASSVALUE=3.0);

julia> psf[1, 1], sampling
(5.0, 0.001)
```
