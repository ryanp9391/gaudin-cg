#!/usr/bin/env bash
export LC_ALL=C OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
cd /home/paul/git_projects/gaudin-CG/Paul/Python/gaudin_gl4; R=../../Cpp/AlexCGC/results/su4-same-machine; PY=../.venv/bin/python
echo "load before $(cat /proc/loadavg)" > $R/ours_load.txt
$PY -X utf8 $R/ours_products.py > $R/ours_products.txt 2>&1
echo "after products $(cat /proc/loadavg)" >> $R/ours_load.txt
$PY -X utf8 slice_engine_gl4.py bench 12 > $R/ours_bench12_single.txt 2>&1
echo "after bench12 $(cat /proc/loadavg)" >> $R/ours_load.txt
$PY -X utf8 slice_ops_gl4.py bench > $R/ours_slice_ops_bench_single.txt 2>&1
echo "after slice_ops bench $(cat /proc/loadavg)" >> $R/ours_load.txt
unset OPENBLAS_NUM_THREADS OMP_NUM_THREADS MKL_NUM_THREADS
$PY -X utf8 slice_engine_gl4.py bench 12 > $R/ours_bench12_multithread.txt 2>&1
echo "after bench12 multithread $(cat /proc/loadavg)" >> $R/ours_load.txt
echo DONE >> $R/ours_load.txt
