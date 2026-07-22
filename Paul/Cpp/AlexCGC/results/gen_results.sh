#!/usr/bin/env bash
#
# Regenerate the highest-weight CGC timing results for SU(4)
#   mu in [S,S,0,0] (x) [S,S,0,0],  for S = 1 .. Smax  (reps = 1).
#
# Usage:
#     bash results/gen_results.sh [Smax]      # default Smax = 8
#
# Writes results/hw_timings.txt (full tables) and results/summary.csv
# (one row per S: singlet time, #irreps, total highest-weight-solve time).
#
# WARNING: cost grows steeply. Rough single-run wall-clock on a modern laptop:
#   S<=5 : seconds     S=6 : ~5 s     S=7 : ~30 s     S=8 : ~3 min
# Going beyond S=8 can take many minutes (the singlet row dominates).
#
set -euo pipefail
cd "$(dirname "$0")/.."          # -> the AlexCGC directory
[ -x ./cgc_hw_bench ] || bash build.sh >/dev/null

SMAX="${1:-8}"
OUT="results/hw_timings.txt"
CSV="results/summary.csv"

{
  echo "Highest-weight-state Clebsch-Gordan coefficient timings"
  echo "SU(4):   mu appearing in   [S,S,0,0] (x) [S,S,0,0]"
  echo
  echo "Generated : $(date -u '+%Y-%m-%d %H:%M UTC')"
  echo "CPU       : $(LC_ALL=C grep -m1 'model name' /proc/cpuinfo | sed 's/.*: //')  ($(nproc) threads)"
  echo "Compiler  : $(g++ --version | head -1)"
  echo "LAPACK    : liblapack-dev $(dpkg -l liblapack-dev 2>/dev/null | awk '/^ii/{print $3}')"
  echo "Tool      : cgc_hw_bench   (built from ClebschGordan.cpp; see ../README.md)"
  echo "Setting   : reps = 1  (single run per mu; times are seconds)"
  echo
  echo "Columns   : hw-solve   = compute_highest_weight_coeffs() only (the CGC solve)"
  echo "            total(+LR) = hw-solve + Littlewood-Richardson decomposition"
  echo "Note      : cost tracks the size of the highest-weight linear system, which"
  echo "            is LARGEST for the singlet (0 0 0 0) -- so the singlet is the"
  echo "            slowest row despite being 1-dimensional."
  echo "======================================================================================"
  for ((S = 1; S <= SMAX; ++S)); do
    echo
    bash test.sh "$S" 1
    echo "--------------------------------------------------------------------------------------"
  done
} > "$OUT"

# Build a compact per-S summary from the full results (portable awk / mawk-safe).
echo "S,n_irreps,singlet_hw_solve_s,total_hw_solve_s" > "$CSV"
LC_ALL=C awk '
  /\(x\)/ && /\[/ { s = $0; sub(/.*\[/, "", s); sub(/,.*/, "", s); S = s }
  /^\(0 0 0 0\)/  { singlet = $7 }              # hw-solve value of the singlet row
  /^irreps mu:/   { printf "%s,%s,%s,%s\n", S, $3, singlet, $7 }
' "$OUT" >> "$CSV"

echo "wrote $OUT and $CSV"
