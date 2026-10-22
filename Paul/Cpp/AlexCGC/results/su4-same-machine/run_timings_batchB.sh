#!/usr/bin/env bash
# G7 batch B: repaired Munich (cgc_hw_prof_ob via test_prof.sh).  Single thread, LC_ALL=C.
export LC_ALL=C OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$(dirname "$0")/../.."
R=results/su4-same-machine
ulimit -v 12000000
run() { # tag mode S reps
  local tag=$1 mode=$2 S=$3 reps=$4
  echo "[$(date +%T)] load before $(cat /proc/loadavg)" >> $R/load_log.txt
  HW_MODE=$mode /usr/bin/time -v timeout 1800 bash test_prof.sh $S $reps > $R/${tag}.txt 2> $R/${tag}.time
  echo "[$(date +%T)] $tag done rc=$? load after $(cat /proc/loadavg) maxrss_kB=$(awk '/Maximum resident/{print $6}' $R/${tag}.time) wall=$(awk '/Elapsed/{print $8}' $R/${tag}.time)" >> $R/load_log.txt
}
for S in 1 2 3 4 5 6 7 8 9 10 11 12; do run eig_S$S eig $S 1; done
for S in 1 2 3 4 5 6 7 8; do run compact_S$S compact $S 1; done
for S in 1 2 3 4 5 6 7; do run proforig_S$S orig $S 1; done
for S in 4 6 8 10; do run eig3_S$S eig $S 3; done
for S in 4 6 8; do run compact3_S$S compact $S 3; done
echo BATCH_B_DONE >> $R/load_log.txt
