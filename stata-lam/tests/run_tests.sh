#!/bin/sh
set -eu

# Invoke from any directory. Override STATA_BIN for a different installation.
package_dir=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
workspace_dir=$(dirname -- "$package_dir")
stata_bin=${STATA_BIN:-}
if [ -z "$stata_bin" ]; then
    for candidate in stata-mp stata-se stata; do
        if command -v "$candidate" >/dev/null 2>&1; then
            stata_bin=$(command -v "$candidate")
            break
        fi
    done
fi
if [ -z "$stata_bin" ] && [ -x /Applications/StataNow/StataMP.app/Contents/MacOS/stata-mp ]; then
    stata_bin=/Applications/StataNow/StataMP.app/Contents/MacOS/stata-mp
fi
if [ -z "$stata_bin" ] || [ ! -x "$stata_bin" ]; then
    printf '%s\n' 'Set STATA_BIN to a licensed Stata command-line executable.' >&2
    exit 1
fi

# Batch startup logs stay temporary because they contain local license details.
scratch_dir=$(mktemp -d "${TMPDIR:-/tmp}/lam-tests.XXXXXX")
trap 'rm -rf "$scratch_dir"' EXIT HUP INT TERM
run_check() {
    test_name=$1
    pass_marker=$2
    rm -f "$package_dir/tests/$test_name.log"
    (cd "$scratch_dir" && "$stata_bin" -bq do "$package_dir/tests/$test_name.do" "$workspace_dir")
    if ! grep -q "$pass_marker" "$package_dir/tests/$test_name.log"; then
        printf '%s\n' "FAILED $test_name" >&2
        tail -n 25 "$package_dir/tests/$test_name.log" >&2
        exit 1
    fi
    printf '%s\n' "PASS $test_name"
}
run_check test_lam 'ALL CORE TESTS PASSED'
run_check test_derivatives 'DERIVATIVE AND TAIL TESTS PASSED'
run_check test_separation 'SEPARATION AND PERIOD RECODING TESTS PASSED'
run_check test_postestimation 'POSTESTIMATION TESTS PASSED'
if [ -f "$workspace_dir/lam-main/lam-china.dta" ] || [ -f "$workspace_dir/lam-china.dta" ]; then
    run_check test_china 'CHINA REPLICATION PASSED'
else
    printf '%s\n' 'SKIP China replication because the original dataset is absent.'
fi
printf '%s\n' 'All available LAM checks passed.'
