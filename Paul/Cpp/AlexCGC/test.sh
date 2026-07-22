#!/usr/bin/env bash
#
# test.sh -- for SU(4), list every irrep mu in the tensor product
#            [S,S,0,0] (x) [S,S,0,0] and report how long it takes to compute
#            the highest-weight-state CGCs of each mu (no coefficient dumps).
#
# Usage:
#     bash test.sh S [reps]
#
#   S      the rectangular-rep parameter: nu = [S,S,0,0]
#   reps   timed repetitions per mu (default 51); reported time is the minimum
#
# Example:  bash test.sh 2
#
set -euo pipefail
cd "$(dirname "$0")"

if [ "$#" -lt 1 ]; then
    echo "Usage: bash test.sh S [reps]" >&2
    exit 1
fi
S="$1"
REPS="${2:-51}"
NU="$S $S 0 0"

# Build the benchmark if it isn't there yet.
if [ ! -x ./cgc_hw_bench ]; then
    echo "building cgc_hw_bench ..." >&2
    bash build.sh >/dev/null
fi

# Format a raw number of seconds (always in seconds; LC_ALL=C so the decimal
# point parses on a non-English locale).
fmt() { LC_ALL=C awk -v s="$1" 'BEGIN{ printf "%.9f s", s }'; }

echo "SU(4):  irreps mu in  [${S},${S},0,0] (x) [${S},${S},0,0]"
echo "        highest-weight-state CGC timing  (reps=${REPS}, headline = minimum)"
echo
printf '%-16s %8s %6s %14s %14s\n' "mu" "dim" "mult" "hw-solve" "total(+LR)"
printf '%-16s %8s %6s %14s %14s\n' "----------------" "--------" "------" "--------------" "--------------"

n_mu=0
tot_solve=0
# Enumerate the mu's (fast: decomposition only), then time each one.
while read -r _tag m1 m2 m3 m4 _kw_mult mult _kw_dim dim; do
    line=$(./cgc_hw_bench 4 $NU $NU "$m1" "$m2" "$m3" "$m4" "$REPS" machine)
    hw=$(printf '%s\n'  "$line" | sed -n 's/.*hw_solve_s=\([^ ]*\).*/\1/p')
    tot=$(printf '%s\n' "$line" | sed -n 's/.*total_s=\([^ ]*\).*/\1/p')
    printf '%-16s %8s %6s %14s %14s\n' \
        "($m1 $m2 $m3 $m4)" "$dim" "$mult" "$(fmt "$hw")" "$(fmt "$tot")"
    n_mu=$((n_mu + 1))
    tot_solve=$(LC_ALL=C awk -v a="$tot_solve" -v b="$hw" 'BEGIN{printf "%.9g", a+b}')
done < <(./cgc_hw_bench decompose 4 $NU $NU)

echo
echo "irreps mu: ${n_mu}    total highest-weight-solve time: $(fmt "$tot_solve")"
