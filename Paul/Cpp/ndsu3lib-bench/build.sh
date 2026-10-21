#!/bin/bash
# Reproducible build of ndsu3lib (GSL for SU(2) coefficients, double precision, no OpenMP) + driver.
set -e
cd "$(dirname "$0")"
export MAMBA_ROOT_PREFIX=$PWD/tools/mamba
if [ ! -x tools/bin/micromamba ]; then
  mkdir -p tools && curl -sL https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xj -C tools bin/micromamba
fi
[ -d env ] || tools/bin/micromamba create -y -q -p $PWD/env -c conda-forge gfortran gxx gcc cmake make gsl libopenblas liblapack
export PATH=$PWD/env/bin:$PATH
[ -d src_dl/ndsu3lib ] || { mkdir -p src_dl; GIT_TERMINAL_PROMPT=0 git clone https://github.com/nd-nuclear-theory/ndsu3lib src_dl/ndsu3lib; git -C src_dl/ndsu3lib checkout dbd38b30c3342d3deeb1e9c2f1885590926ecc0e; }
cmake -S src_dl/ndsu3lib -B build -DCMAKE_BUILD_TYPE=Release -DSU2COEF_LIBRARY=gsl -DPRECISION=double \
  -DCMAKE_PREFIX_PATH=$PWD/env -DCMAKE_Fortran_FLAGS="-O3" -DCMAKE_C_FLAGS="-O3" -DCMAKE_CXX_FLAGS="-O3"
cmake --build build -j8
g++ -O3 -std=c++17 -I src_dl/ndsu3lib driver.cpp -o bench_driver build/libndsu3lib.a -L$PWD/env/lib -lgfortran -lgsl -lgslcblas -llapack -lblas -lm -Wl,-rpath,$PWD/env/lib
