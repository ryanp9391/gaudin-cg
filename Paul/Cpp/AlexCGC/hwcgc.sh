#!/usr/bin/env bash
#
# hwcgc.sh — print ONLY the highest-weight-state Clebsch-Gordan coefficients
# of a chosen target irrep S'' in the product S (x) S'.
#
# The highest-weight state of S'' is the unique state annihilated by every
# raising operator; in this program it is the state with the LARGEST index
# Q(M'') = dim(S''). This script runs the calculator (menu option 5) and keeps
# only those rows.
#
# Usage:
#     bash hwcgc.sh N "S" "S'" "S''" [alpha]
#
# where N is the N of SU(N); S, S', S'' are irreps given as N space-separated
# integers (quote them); alpha is the outer-multiplicity index (default 1,
# only relevant when S'' occurs more than once in the decomposition).
#
# Examples:
#     bash hwcgc.sh 3 "2 1 0" "2 1 0" "3 0 0"
#     bash hwcgc.sh 4 "1 0 0 0" "1 0 0 0" "2 0 0 0"
#
set -e
cd "$(dirname "$0")"

if [ "$#" -lt 4 ]; then
    echo "Usage: bash hwcgc.sh N \"S\" \"S'\" \"S''\" [alpha]" >&2
    exit 1
fi

N="$1"; S="$2"; SP="$3"; SPP="$4"; ALPHA="${5:-1}"

# Drive menu option 5: choice, N, S, S', S'', alpha, blank filename (screen), quit.
OUT=$(printf '5\n%s\n%s\n%s\n%s\n%s\n\n0\n' "$N" "$S" "$SP" "$SPP" "$ALPHA" | ./clebsch 2>&1)

# Keep only the numeric coefficient rows (4 columns: Q(M) Q(M') Q(M'') CGC),
# find the maximum Q(M'') (= dim(S''), the highest-weight state), print those.
# LC_ALL=C so the decimal points parse correctly regardless of system locale.
echo "Highest-weight-state CGCs of S''=($SPP) in ($S) (x) ($SP),  alpha=$ALPHA"
echo "Q(M)"$'\t'"Q(M')"$'\t'"Q(M'')"$'\t'"CGC"
printf '%s\n' "$OUT" | LC_ALL=C awk '
    $1 ~ /^[0-9]+$/ && $2 ~ /^[0-9]+$/ && $3 ~ /^[0-9]+$/ && NF==4 {
        rows[NR]=$1"\t"$2"\t"$3"\t"$4; q[NR]=$3+0; if ($3+0>mx) mx=$3+0
    }
    END {
        for (i=1;i<=NR;i++) if (i in rows && q[i]==mx) print rows[i]
    }'
