#!/usr/bin/env sh
set -eu

command_text=${2:-}
status_file=$(printf '%s\n' "$command_text" | awk -F'"' '{print $6}')

if [ -z "$status_file" ]; then
    echo "Fake Stata could not resolve the launcher-owned status path." >&2
    exit 2
fi

case "${NRH_FAKE_STATA_MODE:-success}" in
    success)
        printf '0\n' >"$status_file"
        ;;
    stata_failure)
        printf '459\n' >"$status_file"
        ;;
    invalid_status)
        printf 'not-a-return-code\n' >"$status_file"
        ;;
    missing_status)
        ;;
    process_failure)
        printf '0\n' >"$status_file"
        exit 7
        ;;
    *)
        echo "Unknown fake-Stata mode." >&2
        exit 2
        ;;
esac
