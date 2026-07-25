#!/usr/bin/env sh
set -eu

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
repo_root=$(CDPATH= cd -- "$script_dir/.." && pwd)
fake_stata="$repo_root/tests/fake_stata.sh"

run_fake() {
    mode=$1
    shift
    NRH_FAKE_STATA_MODE=$mode \
        STATA_BIN="$fake_stata" \
        STATA_BATCH_FLAG=-e \
        "$repo_root/scripts/run_verify.sh" "$@"
}

run_fake success

if run_fake stata_failure; then
    echo "Verify launcher accepted a nonzero Stata return code." >&2
    exit 1
fi

if run_fake invalid_status; then
    echo "Verify launcher accepted a malformed status sidecar." >&2
    exit 1
fi

if run_fake missing_status; then
    echo "Verify launcher accepted a missing status sidecar." >&2
    exit 1
fi

if run_fake process_failure; then
    echo "Verify launcher ignored an unexpected process failure." >&2
    exit 1
fi

if run_fake success 'invalid"argument'; then
    echo "Verify launcher accepted an unsafe quoted argument." >&2
    exit 1
fi

if "$repo_root/scripts/run_smoke.sh"; then
    echo "Retired smoke launcher returned success." >&2
    exit 1
fi

echo "Verify-launcher status tests passed."
