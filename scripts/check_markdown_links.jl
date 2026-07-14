#!/usr/bin/env julia

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))

function live_markdown_files()
    files = String[
        joinpath(REPO_ROOT, "README.md"),
        joinpath(REPO_ROOT, "PORTING_PLAN.md"),
        joinpath(REPO_ROOT, "bench", "suites", "README.md"),
    ]
    append!(
        files,
        sort!(filter(path -> endswith(path, ".md"), readdir(joinpath(REPO_ROOT, "docs"); join=true))),
    )
    return files
end

function local_target(raw_target::AbstractString)
    raw = strip(raw_target)
    isempty(raw) && return nothing
    any(prefix -> startswith(raw, prefix), ("http://", "https://", "mailto:", "data:", "#")) &&
        return nothing

    target = if startswith(raw, "<") && occursin('>', raw)
        first(split(raw[2:end], '>'; limit=2))
    else
        first(split(raw))
    end
    target = first(split(first(split(target, '#'; limit=2)), '?'; limit=2))
    isempty(target) && return nothing
    return replace(target, "%20" => " ")
end

function check_markdown_links(files)
    failures = String[]
    link_pattern = r"!?\[[^\]]*\]\(([^)\n]+)\)"
    for file in files
        source = read(file, String)
        for match in eachmatch(link_pattern, source)
            target = local_target(match.captures[1])
            isnothing(target) && continue
            resolved = startswith(target, '/') ?
                normpath(joinpath(REPO_ROOT, target[2:end])) :
                normpath(joinpath(dirname(file), target))
            ispath(resolved) && continue
            push!(
                failures,
                "$(relpath(file, REPO_ROOT)): $(match.captures[1]) -> $(relpath(resolved, REPO_ROOT))",
            )
        end
    end
    return failures
end

files = live_markdown_files()
failures = check_markdown_links(files)
if !isempty(failures)
    foreach(failure -> println(stderr, "broken local Markdown target: ", failure), failures)
    exit(1)
end
println("Validated local Markdown targets in $(length(files)) live documents.")
