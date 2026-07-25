#!/usr/bin/env sh
set -eu

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
repo_root=$(CDPATH= cd -- "$script_dir/.." && pwd)

stata_bin=${STATA_BIN:-stata-mp}
if [ -n "${STATA_BATCH_FLAG:-}" ]; then
    batch_flag=$STATA_BATCH_FLAG
elif [ "$(uname -s)" = "Darwin" ]; then
    batch_flag=-e
else
    batch_flag=-b
fi

case "$stata_bin" in
    */*)
        if [ ! -x "$stata_bin" ]; then
            echo "Could not execute Stata: $stata_bin" >&2
            exit 127
        fi
        ;;
    *)
        if ! command -v "$stata_bin" >/dev/null 2>&1; then
            echo "Could not find Stata executable: $stata_bin" >&2
            echo "Set STATA_BIN to the Stata executable path." >&2
            exit 127
        fi
        ;;
esac

temp_dir=$(mktemp -d "${TMPDIR:-/tmp}/nrh-verify.XXXXXX")
status_file="$temp_dir/stata-status.txt"

cleanup() {
    rm -f -- "$status_file" \
        "$temp_dir/run_verify.log" \
        "$temp_dir/run_verify.smcl"
    rmdir -- "$temp_dir" 2>/dev/null || true
}
trap cleanup EXIT HUP INT TERM

cd "$temp_dir"

stata_command="do \"$repo_root/scripts/run_verify.do\" \"$repo_root\" \"$status_file\""
for stata_arg in "$@"; do
    case "$stata_arg" in
        *\"*)
            echo "Verify-launcher arguments may not contain a double quote." >&2
            exit 198
            ;;
    esac
    stata_command="$stata_command \"$stata_arg\""
done

set +e
"$stata_bin" "$batch_flag" "$stata_command"
stata_process_rc=$?
set -e

if [ ! -f "$status_file" ]; then
    echo "Stata did not produce the verify-run status sidecar." >&2
    if [ "$stata_process_rc" -ne 0 ]; then
        exit "$stata_process_rc"
    fi
    exit 1
fi

stata_run_rc=$(sed -n '1p' "$status_file")
case "$stata_run_rc" in
    ''|*[!0-9]*)
        echo "Stata produced an invalid verify-run status: $stata_run_rc" >&2
        exit 1
        ;;
esac

if [ "$stata_run_rc" -ne 0 ]; then
    echo "NRH no-data verification failed with Stata return code $stata_run_rc." >&2
    exit 1
fi

if [ "$stata_process_rc" -ne 0 ]; then
    echo "Stata exited unexpectedly with process status $stata_process_rc." >&2
    exit "$stata_process_rc"
fi

echo "NRH no-data verification completed successfully."
