#!/usr/bin/env bash
# G7 batch A: Munich as distributed (cgc_hw_bench via test.sh), reps=1 for S=1..8, reps=3 for S=4,6,8.  Single thread, LC_ALL=C.
export LC_ALL=C OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$(dirname "$0")/../.."
R=results/su4-same-machine
ulimit -v 12000000
run() { # tag S reps script [env]
  local tag=$1 S=$2 reps=$3 scr=$4
  echo "[$(date +%T)] load before $(cat /proc/loadavg)" >> $R/load_log.txt
  /usr/bin/time -v timeout 1800 bash $scr $S $reps > $R/${tag}.txt 2> $R/${tag}.time
  echo "[$(date +%T)] $tag done rc=$? load after $(cat /proc/loadavg) maxrss_kB=$(awk '/Maximum resident/{print $6}' $R/${tag}.time) wall=$(awk '/Elapsed/{print $8}' $R/${tag}.time)" >> $R/load_log.txt
}
for S in 1 2 3 4 5 6 7 8; do run orig_S$S $S 1 test.sh; done
for S in 4 6 8; do run orig3_S$S $S 3 test.sh; done
echo BATCH_A_DONE >> $R/load_log.txt
