# AlexCGC — SU(N) Clebsch–Gordan coefficient calculator

C++ program from Alex, Kalus, Huckleberry & von Delft, *"A numerical algorithm
for the explicit calculation of SU(N) and SL(N,C) Clebsch–Gordan coefficients"*
(arXiv:1009.0437). It computes, for any N, the decomposition of a tensor product
S ⊗ S′ of two su(N) irreps and the associated Clebsch–Gordan coefficients, using
the Gelfand–Tsetlin pattern calculus.

Kept here as an independent cross-check for the su(4) CG work in this repo
(companion to the Python port under `Paul/Python/`). This is upstream code, not
ours — the algorithm lives in the paper.

## Files

- `ClebschGordan.cpp` — the authors' original source, downloaded from
  http://homepages.physik.uni-muenchen.de/~vondelft/Papers/ClebschGordan/ClebschGordan.cpp
  (unmodified).
- `build.sh` — one-command rebuild.
- `clebsch` — the compiled program (rebuild it if you're on a different machine;
  a compiled binary is not portable across systems).

## One-time setup

The program relies on the standard LAPACK and BLAS numerical libraries. Install
them once (asks for your password):

    sudo apt-get install -y liblapack-dev libblas-dev

## Build

    bash build.sh

(Equivalently, by hand: `g++ ClebschGordan.cpp -O2 -o clebsch -llapack -lblas`.)

## Run

    ./clebsch

It shows a menu. The two most useful options:

- **6** — all Clebsch–Gordan coefficients for a product S ⊗ S′.
- **5** — coefficients for one specific target irrep S″ in S ⊗ S′.

It then asks for:

- **N** — the N in SU(N), e.g. `4`.
- **Irrep S** and **Irrep S′** — each entered as N space-separated integers, a
  highest-weight / Young-diagram row-length label. Examples: `2 1 0` for SU(3),
  `1 0 0 0` (the defining **4**) for SU(4).
- **File name** — press Enter to print to the screen, or type a filename to save.

To **quit**, choose `0` at the menu (or press Ctrl+D).

### Reading the output

For a product it first prints the Littlewood–Richardson decomposition, one line
per irrep S″:

    [irrep index] S'' (outer multiplicity) {dimension}

Then, for each S″, a table of nonzero coefficients:

    Q(M)  Q(M')  Q(M'')  CGC

where Q(M), Q(M′), Q(M″) are integer indices labelling the basis states of S, S′
and S″ respectively (menu options 3 and 4 convert between these indices and the
explicit Gelfand–Tsetlin patterns), and CGC is the coefficient.

## Highest-weight-state CGCs only

The highest-weight state of a target irrep S″ is the unique state killed by all
raising operators. Because of the index convention (App. C.2), it is always the
state with the **largest** index `Q(M″) = dim(S″)`, and its coefficients are the
ones the algorithm computes directly (Eq. 35) — every other state is derived from
them by lowering. So to get just those, run the calculator and keep the rows with
the maximum `Q(M″)`.

The helper `hwcgc.sh` does exactly this:

    bash hwcgc.sh N "S" "S'" "S''" [alpha]

e.g. `bash hwcgc.sh 4 "1 0 0 0" "1 0 0 0" "1 1 0 0"` prints the highest-weight
CGCs of the **6** in 4 ⊗ 4. `alpha` (default 1) selects the copy when S″ has outer
multiplicity > 1.

Note: post-processing the numbers with `awk` etc. needs `LC_ALL=C` on a
non-English locale, or decimal points are misread (the program's own output is
unaffected).

## Benchmarking just the highest-weight state (`cgc_hw_bench`)

`cgc_hw_bench` is a modified copy of the program (see `cgc_hw_bench.cpp`) that
computes **only** the highest-weight-state CGCs of one target irrep S″ and
reports how long that takes — *without* building the descendant states. This
isolates the cost of the highest-weight computation for comparison against other
methods. The original algorithm already computes the highest-weight coefficients
as a distinct first step (`compute_highest_weight_coeffs`, the raising-operator
null-space / SVD); this build simply times that step and skips the rest.

    bash build.sh                       # builds ./clebsch and ./cgc_hw_bench
    ./cgc_hw_bench N  s1..sN  p1..pN  q1..qN  [reps] [verify]

- `N` — the N of SU(N).
- `s1..sN`, `p1..pN`, `q1..qN` — the irreps S, S′, and target S″, each as N
  integers.
- `reps` — timed repetitions (default 11); the headline figure is the **minimum**
  (least affected by OS jitter). Increase it for sub-millisecond cases.
- `verify` — also run the full (all-descendants) computation once, time it,
  report the speedup, and check the highest-weight coefficients agree exactly.

It reports three timings so the comparison is fair regardless of what your own
method measures:

1. **highest-weight solve only** — the pure CGC computation for the h.w. state.
2. **decomposition (LR) only** — the prerequisite step that confirms S″ occurs.
3. **total** — labels → highest-weight CGCs (1 + 2).

It then prints the highest-weight coefficients themselves (per outer-multiplicity
copy α), with a `sum of CGC^2 = 1` normalization check and, when α > 1, the
overlaps between copies (≈ 0).

Example:

    ./cgc_hw_bench 4  3 2 1 0  3 2 1 0  6 4 2 0  51 verify

For that SU(4) case (target dim 729) the highest-weight solve is ~0.2 ms while
the full all-states computation is ~0.65 s — a ~3000× difference, confirming the
descendants dominate the original program's cost.

### Rectangular-rep family sweep (`test.sh`)

`test.sh` is a convenience wrapper for the specific SU(4) family of this project:
it lists every irrep mu appearing in the tensor product `[S,S,0,0] (x) [S,S,0,0]`
and reports the highest-weight-state CGC timing for each (no coefficient dumps).

    bash test.sh S [reps]

- `S` — the rectangular-rep parameter; nu = `[S,S,0,0]`.
- `reps` — timed repetitions per mu (default 51); headline is the minimum.

It prints one row per mu with its dimension, outer multiplicity, the
highest-weight-solve time, and the total (solve + Littlewood-Richardson), plus a
summed highest-weight-solve time across all mu. Example (`bash test.sh 2`,
20 x 20 = 400):

    mu                    dim   mult       hw-solve     total(+LR)
    (0 0 0 0)               1      1  0.000141885 s  0.000144160 s
    (2 1 1 0)              15      1  0.000131038 s  0.000135079 s
    (2 2 0 0)              20      1  0.000050382 s  0.000052579 s
    (4 2 2 0)              84      1  0.000071416 s  0.000074655 s
    (4 3 1 0)             175      1  0.000070702 s  0.000074115 s
    (4 4 0 0)             105      1  0.000041731 s  0.000044133 s

Note the timing is NOT monotonic in dim(mu): it is governed by the size of the
highest-weight linear system (the number of product states `|M> (x) |M'>` whose
weights add up to the highest weight of mu), which is largest for the singlet
`(0 0 0 0)` -- hence the singlet is consistently the slowest row.

### Saved results (`results/`)

- `results/hw_timings.txt` — full per-mu timing tables for `[S,S,0,0]²`, S = 1..8,
  with a metadata header (machine, compiler, LAPACK version, date).
- `results/summary.csv` — one row per S: `S, n_irreps, singlet_hw_solve_s,
  total_hw_solve_s` (handy for plotting the scaling).
- `results/gen_results.sh` — regenerates both from scratch:
  `bash results/gen_results.sh [Smax]` (default Smax=8; note S=8 takes ~3 min,
  and cost grows steeply beyond that — the singlet row dominates).

Snapshot (13th Gen i5-1345U, reps=1): the singlet's highest-weight solve grows
0.0004 s → 0.099 s → 18.7 s across S = 2, 5, 8 (~×5 per step in S), and is
consistently the slowest row of each family.

### What was changed vs. the original

`cgc_hw_bench.cpp` is `ClebschGordan.cpp` with three additions, all clearly
commented with `// added`:
1. a `highest_weight_only` flag on the `coefficients` constructor that returns
   right after `compute_highest_weight_coeffs()` instead of building descendants;
2. a scope timer that records the duration of `compute_highest_weight_coeffs()`;
3. a new `main()` (the benchmark driver) in place of the interactive menu.
The algorithm/numerics are untouched.

## Verified examples

- SU(3): `(2 1 0) ⊗ (2 1 0)` reproduces Eq. (32) of the paper —
  `(0,0,0) ⊕ (2,1,0)×2 ⊕ (3,0,0) ⊕ (3,3,0) ⊕ (4,2,0)`, dimensions
  1+8+8+10+10+27 = 64 = 8×8, with (2,1,0) appearing at outer multiplicity 2.
- SU(4): `(1 0 0 0) ⊗ (1 0 0 0)` gives `(1,1,0,0) ⊕ (2,0,0,0)` — the **6** and
  **10**, with 6+10 = 16 = 4×4.
