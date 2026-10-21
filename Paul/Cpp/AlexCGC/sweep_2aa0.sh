#!/usr/bin/env bash
#
# sweep_2aa0.sh -- SU(3): highest-weight-state CGC timing for the target
#                  mu = [2a,a,0] inside [2a,a,0] (x) [2a,a,0], for a = 1..amax,
#                  run sequentially (a = 1 is the adjoint 8 in 8 (x) 8).
#
# Usage:
#     bash sweep_2aa0.sh [amax] [reps] [limit_s]
#
#   amax     largest a to run (default 3)
#   reps     timed repetitions per a (default 5); headline = minimum.
#            Note cgc_hw_bench also does one untimed warm-up solve, so the
#            wall-clock per a is roughly (reps + 1) x hw-solve.
#   limit_s  stop the sweep once a single hw-solve exceeds this many seconds
#            (default 60), so a steep blow-up doesn't run for hours.
#
# Output: one table row per a, plus the growth factor hw-solve(a)/hw-solve(a-1).
# Also appended as CSV to results/sweep_2aa0.csv (a,dim,mult,hw_solve_s,total_s).
#
set -euo pipefail
cd "$(dirname "$0")"

AMAX="${1:-3}"
REPS="${2:-5}"
LIMIT="${3:-60}"

[ -x ./cgc_hw_bench ] || bash build.sh >/dev/null

CSV="results/sweep_2aa0.csv"
[ -f "$CSV" ] || echo "a,dim,mult,hw_solve_s,total_s" > "$CSV"

echo "SU(3):  mu = [2a,a,0]  in  [2a,a,0] (x) [2a,a,0]    (reps=${REPS}, headline = minimum)"
echo "        $(LC_ALL=C grep -m1 'model name' /proc/cpuinfo 2>/dev/null | sed 's/.*: //')"
echo
printf '%-4s %-12s %8s %6s %16s %16s %10s %10s\n' "a" "mu" "dim" "mult" "hw-solve" "total(+LR)" "growth" "wall"
printf '%-4s %-12s %8s %6s %16s %16s %10s %10s\n' "----" "------------" "--------" "------" "----------------" "----------------" "----------" "----------"

prev=""
for ((a = 1; a <= AMAX; ++a)); do
    R="$((2*a)) $a 0"
    w0=$(date +%s.%N)
    line=$(./cgc_hw_bench 3 $R $R $R "$REPS" machine)
    w1=$(date +%s.%N)
    dim=$(printf '%s\n'  "$line" | sed -n 's/.*dim=\([^ ]*\).*/\1/p')
    mult=$(printf '%s\n' "$line" | sed -n 's/.*mult=\([^ ]*\).*/\1/p')
    hw=$(printf '%s\n'   "$line" | sed -n 's/.*hw_solve_s=\([^ ]*\).*/\1/p')
    tot=$(printf '%s\n'  "$line" | sed -n 's/.*total_s=\([^ ]*\).*/\1/p')
    LC_ALL=C awk -v a="$a" -v mu="($R)" -v d="$dim" -v m="$mult" -v hw="$hw" -v t="$tot" \
                 -v p="$prev" -v w0="$w0" -v w1="$w1" 'BEGIN{
        g = (p == "" || p == 0) ? "-" : sprintf("x%.1f", hw / p);
        printf "%-4s %-12s %8s %6s %14.9f s %14.9f s %10s %8.1f s\n", a, mu, d, m, hw, t, g, w1 - w0 }'
    echo "$a,$dim,$mult,$hw,$tot" >> "$CSV"
    prev="$hw"
    if LC_ALL=C awk -v hw="$hw" -v L="$LIMIT" 'BEGIN{ exit !(hw > L) }'; then
        echo
        echo "hw-solve at a=$a exceeded ${LIMIT} s -- stopping the sweep here."
        break
    fi
done
