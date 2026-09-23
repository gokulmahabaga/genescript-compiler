#!/usr/bin/env bash
# Test runner for gsc (called by `make test`).
#
# 1. Every program in examples/ and tests/cases/ is compiled and run with
#    BOTH front ends (Bison and recursive-descent). Each run's output must
#    match the saved <name>.expected file next to it, so the two parsers
#    are checked against each other as well as against the right answer.
# 2. Every program in tests/errors/ must be REJECTED (non-zero exit)
#    by both front ends, with an error message.
#
# A test can pass extra gsc options by putting them in <name>.flags
# (e.g. "--show-opt --dump-tac"), so the dump/optimizer output is tested too.
#
# To (re)generate an expected file after deliberately changing output:
#   tests/run_tests.sh --update

set -u
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
GSC="$ROOT/gsc"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

UPDATE=0
[ "${1:-}" = "--update" ] && UPDATE=1

pass=0; fail=0

# Compile+run one file in the scratch dir, output without gsc's own [gsc] log lines.
run_gsc() {
    local flags=""
    [ -f "${1%.gs}.flags" ] && flags="$(cat "${1%.gs}.flags")"
    # shellcheck disable=SC2086
    (cd "$WORK" && "$GSC" "$1" --frontend="$2" $flags 2>&1) | grep -v '^\[gsc\]'
    return "${PIPESTATUS[0]}"
}

for f in "$ROOT"/examples/*.gs "$ROOT"/tests/cases/*.gs; do
    [ -e "$f" ] || continue
    expected="${f%.gs}.expected"
    name="${f#$ROOT/}"

    if [ "$UPDATE" = 1 ]; then
        run_gsc "$f" bison > "$expected"
        echo "updated $expected"
        continue
    fi
    if [ ! -f "$expected" ]; then
        echo "MISSING  $name (no .expected file -- run tests/run_tests.sh --update)"
        fail=$((fail + 1)); continue
    fi

    for fe in bison rd; do
        actual="$(run_gsc "$f" "$fe")"
        if [ "$actual" = "$(cat "$expected")" ]; then
            echo "PASS     $name [$fe]"; pass=$((pass + 1))
        else
            echo "FAIL     $name [$fe]"; fail=$((fail + 1))
            diff <(cat "$expected") <(echo "$actual") | sed 's/^/           /'
        fi
    done
done

[ "$UPDATE" = 1 ] && exit 0

for f in "$ROOT"/tests/errors/*.gs; do
    [ -e "$f" ] || continue
    name="${f#$ROOT/}"
    for fe in bison rd; do
        out="$(run_gsc "$f" "$fe")"; code=$?
        if [ "$code" -ne 0 ] && echo "$out" | grep -q "Error"; then
            echo "PASS     $name [$fe] rejected: $(echo "$out" | grep Error | head -1)"
            pass=$((pass + 1))
        else
            echo "FAIL     $name [$fe] should have been rejected (exit $code)"
            fail=$((fail + 1))
        fi
    done
done

echo
echo "$pass passed, $fail failed"
[ "$fail" -eq 0 ]
