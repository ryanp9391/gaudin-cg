#!/usr/bin/env bash
#
# Rebuilds the `clebsch` program from ClebschGordan.cpp.
# Run it with:   bash build.sh
#
# Requires the LAPACK and BLAS libraries. If the build complains that it
# cannot find -llapack or -lblas, install them once with:
#     sudo apt-get install -y liblapack-dev libblas-dev
#
set -e
cd "$(dirname "$0")"

echo "Compiling ClebschGordan.cpp ..."
g++ ClebschGordan.cpp -O2 -o clebsch -llapack -lblas
echo "  built: ./clebsch          (interactive calculator)"

echo "Compiling cgc_hw_bench.cpp ..."
g++ cgc_hw_bench.cpp -O2 -o cgc_hw_bench -llapack -lblas
echo "  built: ./cgc_hw_bench     (highest-weight-state timing benchmark)"

echo "Done."
