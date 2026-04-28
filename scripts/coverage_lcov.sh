#!/usr/bin/env bash
set -euo pipefail

output="${1:-lcov.info}"

find src test examples reference_models -name '*.cov' -delete

julia --project=. --code-coverage=user -e 'using Pkg; Pkg.test(coverage=true)'

julia --project=coverage -e '
using Pkg
Pkg.instantiate()

using Coverage

coverage = process_folder("src")
LCOV.writefile(ARGS[1], coverage)

covered, total = get_summary(coverage)
percent = total == 0 ? 0.0 : 100 * covered / total
println("Coverage: $(round(percent; digits=2))% ($(covered)/$(total))")
' "${output}"

if [[ "${KEEP_JULIA_COV:-0}" != "1" ]]; then
  find src test examples reference_models -name '*.cov' -delete
fi
