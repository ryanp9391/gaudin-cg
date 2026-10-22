# Spec G7: the Munich comparison for gl(4) on the Ubuntu laptop — correctness of the HW Clebsch–Gordan vectors (multiplicity-free AND multiplicity sectors) and a same-machine timing table (written 2026-10-22 on Windows; to be run on the Ubuntu laptop)

For a general-purpose subagent (SONNET) in a session on the **Ubuntu laptop** (i5-1345U, 12 threads, 15 GB, Ubuntu, g++ 11.4, reference LAPACK/BLAS 3.10 installed in July 2026;
Python venv `Paul/Python/.venv` with Python 3.10.12, numpy 2.2.6 (bundled OpenBLAS, multithreaded), scipy 1.15.3, sympy 1.14.0, mpmath 1.3.0 — the venv used for the gl(3)
timings in `Paul/Cpp/ndsu3lib-bench/python-timings.md`). **No Wolfram, no notebooks, no git commits unless Paul asks** (read-only `git status`/`git log` are fine).
Never edit the Windows-verified modules (`glN_reps.py`, `slice_engine_gl4.py`, `gaudin_gl4_ops.py`, `slice_ops_gl4.py`, `oper_gl4.py`, `structure_gl4.py`, the tests, any RESULTS/SPEC)
except for a Python-version shim if one is unavoidable (then say exactly what and why in RESULTS). Never overwrite `Paul/Cpp/AlexCGC/results/hw_timings.txt` or `summary.csv`
(the published reference numbers). Files you own: `Paul/Python/gaudin_gl4/munich_bridge.py`, `RESULTS-G7-munich-comparison.md`, `munich_bridge_output.txt`,
`Paul/Cpp/AlexCGC/test_prof.sh`, `Paul/Cpp/AlexCGC/results/su4-same-machine/` (new directory: every raw output file), and one new block at the top of `HANDOFF.md` when done.

Read first: this file; `RESULTS-G1-slice-engine-gl4.md` §0, §6, §9 (the Munich bridge: Q(M) = our GT index + 1, orthonormal GT basis = ours rescaled by 1/√g, trivial sign gauge,
the conversion CGC = σ·Ψ[ν,q]·√(g_ν g_q/GNORM2), the comparer `munichcmp`, the 19 commands); `RESULTS-G3-slice-ops-gl4.md` §8, §9, §13 (canonical copies, G⊗G-orthogonality,
`copies`, `copy_float_slices`, the bench stage); `benchmark_gl4.md` (the table this package completes with same-machine columns); `Paul/Cpp/AlexCGC/README.md`, `hwcgc.sh`,
`test.sh`, `results/gen_results.sh`, `cgc_hw_prof.cpp` (header comments: `HW_MODE` = orig | compact | eig); `Paul/Cpp/ndsu3lib-bench/python-timings.md` (how the gl(3)
same-machine timings were taken on this laptop: venv, load recording, min of 5, single process) and the HANDOFF block "DONE (2026-10-21 night, Ubuntu): ndsu3lib benchmark".

## 0. What this settles

The Windows session built and verified the gl(4) slice engine (G1), the Gaudin missing-label operators (G2), their slice reduction and the canonical copies (G3) and the oper
side (G5), but the only general-N reference code, Munich (`ClebschGordan.cpp`, arXiv:1009.0437), could not be built on Windows. What was settled without it: the two GT bases
coincide (pattern order, orthonormal lowering coefficients, trivial sign gauge, G1 §9), and a comparer exists and is self-tested on synthetic Munich output. **What this package
delivers:** (D1) the literal comparison of Munich's highest-weight CG coefficients with ours — multiplicity-free irreps vector by vector up to one global sign, multiplicity
sectors subspace by subspace with the rotation between Munich's Gram–Schmidt copies and our Gaudin copies; (D2) a same-machine timing table for SU(4) [S,S,0,0]² (Munich as
distributed, Munich with its linear algebra repaired, our engine); (D3) Munich timings on multiplicity products against our G3 numbers on the same laptop; (D4) the report.

## 1. Setup (do these in order; stop and report if a step fails)

1. `cd` to the repo root (`gaudin-cg`), `git log -1`, `git status --short | head`: the gl(4) material must be present — `Paul/Python/gaudin_gl4/` (modules + RESULTS-G1..G6),
   `Paul/Mathematica/Gaudin/ref/gl4-*.txt` (17 files), `Paul/Mathematica/Gaudin/gl4-slice-engine-notes.{tex,md}`. If missing, Paul has to push from Windows and pull here; stop.
2. Python: use the venv explicitly everywhere: `PY=../.venv/bin/python` from `Paul/Python/gaudin_gl4` (`$PY --version` must print 3.10.x). Run the five test files
   (`$PY -X utf8 test_glN_reps.py`, `test_gaudin_gl4_ops.py`, `test_slice_ops_gl4.py`, `test_oper_gl4.py`, `test_structure_gl4.py`; ≈ 2 + 5 + 2 + 5 + 6 s on Windows): all must
   pass. The code was checked to parse as Python 3.10 but was only ever run under 3.13 on Windows; a runtime difference would show here — report it before anything else.
   Then `$PY -X utf8 slice_engine_gl4.py munich` (prints the 19 commands of step 3 and runs the comparer self-test: 19 accepted, 39 corrupted rejected).
3. Munich: `cd Paul/Cpp/AlexCGC && bash build.sh` (builds `clebsch` and `cgc_hw_bench` with `-llapack -lblas`). Then the repaired profiler, two variants:
   `g++ cgc_hw_prof.cpp -O2 -o cgc_hw_prof -llapack -lblas` (reference BLAS) and `g++ cgc_hw_prof.cpp -O2 -o cgc_hw_prof_ob -lopenblas` (if `libopenblas-dev` is not
   installed and you cannot install it, use the OpenBLAS of the ndsu3lib micromamba env: `-L../ndsu3lib-bench/env/lib -lopenblas -Wl,-rpath,$PWD/../ndsu3lib-bench/env/lib`;
   `ldd cgc_hw_prof_ob | grep -i blas` must show it). `HW_MODE=eig ./cgc_hw_prof_ob 4 1 1 0 0 1 1 0 0 2 2 0 0 5 machine` must print one RESULT line. The three binaries are
   gitignored except `cgc_hw_prof*` — add them to `.gitignore` only if Paul asks; do not commit binaries.
4. Discipline for every timing: `export LC_ALL=C` (French locale breaks awk/sort), `export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1` (single-thread
   comparison: Munich's published table used single-threaded reference BLAS), nothing else running, `cat /proc/loadavg` before and after each batch and quote it, `/usr/bin/time -v`
   for maximum RSS on every Munich run, `timeout 1800` on every single Munich invocation, and before any Munich original run compute its memory need
   n_W · dim V_λ · dim V_λ′ · 8 bytes (n_W from our Racah–Speiser / `glN_reps.nW_of`; the original allocates that much: 9.5 GB for the S = 8 singlet, 33 GB at S = 9) and skip
   runs above 10 GB; `ulimit -v 12000000` as a guard. Record the environment once (`lscpu | head -20`, `g++ --version`, `dpkg -l | grep -i 'lapack\|blas'`, `$PY -c "import numpy; numpy.show_config()"`).

## 2. D1 — correctness: Munich vs ours

**Conventions (from G1 §9).** Munich's basis index Q(M) = our GT index + 1 in every irrep (both sites); Munich's orthonormal basis = our plain GT basis rescaled by 1/√g with a
trivial sign gauge; so for a multiplicity-free irrep the rows `Q(M) Q(M') Q(M'')=dim(μ) CGC` printed by `hwcgc.sh` equal σ·Ψ[ν,q]·√(g_ν g_q/GNORM2) for ONE global sign σ.
Munich's irrep label is ours minus μ₄·(1,1,1,1) (`glN_reps.munich_label_to_ours` converts); `hwcgc.sh N "S" "S'" "S''" [alpha]` drives menu option 5 of `clebsch` (it computes the
FULL table of S'' including descendants and keeps the HW rows, so it is slow for big targets — correctness only) and `alpha` selects Munich's outer-multiplicity copy.

**2a. Multiplicity-free (all m = 1): vector by vector.** From `Paul/Cpp/AlexCGC`:
    bash hwcgc.sh 4 "S S 0 0" "S S 0 0" "m1 m2 m3 m4" > results/su4-same-machine/hw_S_m1_m2_m3_m4.txt
    cd ../../Python/gaudin_gl4 && $PY -X utf8 slice_engine_gl4.py munichcmp S m1 m2 m3 m4 ../../Cpp/AlexCGC/results/su4-same-machine/hw_S_m1_m2_m3_m4.txt
for EVERY irrep of [S,S,0,0]², S = 1, 2, 3, 4 (3 + 6 + 10 + 15 = 34 irreps; the S ≤ 3 commands are printed by `slice_engine_gl4.py munich`; for S = 4 call
`slice_engine_gl4.munich_commands(Slist=(4,))` from Python, it prints the pairs), and for the following irreps of S = 6 and S = 8 (Munich labels): S = 6: `0 0 0 0`, `2 1 1 0`,
`6 4 2 0`, `12 12 0 0`; S = 8: `0 0 0 0` (Munich 18.7 s for the HW solve; menu 5 then builds the descendants of a 1-dimensional irrep, cheap), `2 1 1 0`, `8 6 2 0`, `16 16 0 0`
(`timeout 1800`; if a run exceeds it, say so and drop it). **Acceptance: `MATCH` (same support, max deviation ≤ 1e-9 after one global sign) on all 34 + 8.** Record σ per irrep.

**2b. Multiplicity sectors: subspace by subspace.** Write `munich_bridge.py` (new file; imports `slice_engine_gl4` (`parse_munich_rows`, `to_munich`, `hw_from_slice`,
`prv_space`), `slice_ops_gl4` (`copies`, `copy_float_slices`), `glN_reps` (`racah_speiser`, `munich_label_to_ours`, `GLNRep`) read-only) with
`compare_sector(lam, lamp, mu_ours, dir)`: (i) runs `bash hwcgc.sh 4 "lam" "lamp" "mu_munich" alpha` for alpha = 1..m (m = Racah–Speiser), saves each output, parses the HW rows
(those with Q(M'') = dim μ), assembles Munich's m vectors M_α over the index pairs (Q(M), Q(M')); checks they are orthonormal (Munich's α-copies are Gram–Schmidt orthonormal);
(ii) builds our m canonical copies (`so.copies(lam, lamp, mu, dps=30)` → `so.copy_float_slices` → `se.hw_from_slice` → `se.to_munich(lam, lamp, hw, k)` for k = 0..m−1), which
are G⊗G-orthogonal (G3 §8), normalises them to orthonormal O_k; (iii) computes the m×m overlap Q = Oᵀ M, the principal angles between span(M) and span(O) (SVD of Q: all
singular values = 1 ⇔ same subspace), and reports: max |1 − σ_i(Q)| (acceptance ≤ 1e-8), the matrix Q itself (3 decimals; it is the rotation from Munich's Gram–Schmidt copies in
Munich's index order to the Gaudin eigenbasis — a result Paul wants to see), and whether Q is a signed permutation (it will not be in general). For m = 1 the function must
reduce to 2a (|Q| = 1). Products and sectors: **all** irreps (m ≥ 1) of [2,1,0,0]², [2,1,1,0]², [2,1,0,0]⊗[1,1,0,0], [1,1,0,0]⊗[2,1,0,0], [2,2,1,0]⊗[2,1,1,0], [3,2,1,0]²
(16 irreps, m ≤ 4, incl. (5,4,2,1) where V₃ is degenerate and the joint diagonalisation matters), [4,2,0,0]² (24 irreps, m ≤ 3), [4,2,1,0]² (27 irreps, m ≤ 5; Munich needs
n_W·dim²·8 B ≤ 1 GB here — check), and, if `hwcgc.sh` runs within the timeout, the five largest-m sectors of [6,4,2,0]² with `cgc_hw_bench`/`cgc_hw_prof` instead of `clebsch`
(see §3c: the bench binaries print the HW coefficients per α when asked with `verify`? — check their output format; if they do not print coefficients, skip [6,4,2,0]² and say so).
**Acceptance: max |1 − σ_i(Q)| ≤ 1e-8 and Munich orthonormality ≤ 1e-9 on every sector; the Q matrices tabulated.** CLI: `$PY -X utf8 munich_bridge.py compare` (all of 2a + 2b from
one command, writing every raw Munich output to `results/su4-same-machine/`), `munich_bridge.py selftest` (synthetic Munich outputs from our own vectors, incl. a random rotation
of the copies for m ≥ 2: accepted; corrupted: rejected).

## 3. D2/D3 — timings on this machine (single thread, load ≈ 1, min over runs)

**3a. Munich as distributed (`cgc_hw_bench`, the binary of the published table):** `bash test.sh S 1 > results/su4-same-machine/orig_S<S>.txt` for S = 1..8 (reps = 1 as in
`gen_results.sh`; S = 8 ≈ 3 min; S = 9 needs 33 GB: do not run). Then re-run S = 4, 6, 8 with reps = 3 (`bash test.sh 8 3`, ≈ 10 min) and quote the minimum. Compare with
`results/hw_timings.txt` (same machine, July 2026): the ratios must be ≈ 1 — if not, the laptop's state (load, thermal) is different; say so.
**3b. Munich repaired:** `test_prof.sh` = a copy of `test.sh` that calls `./cgc_hw_prof_ob` with `HW_MODE` taken from the environment (write it; keep the same output format), run
`HW_MODE=eig bash test_prof.sh S 1` and `HW_MODE=compact bash test_prof.sh S 1` for S = 1..8, then S = 9, 10 with `eig` if `/usr/bin/time -v` shows it fits (report max RSS per S);
also `HW_MODE=orig` for one S as a consistency check against 3a. For gl(3) the repair gave ×19 at [14,7,0]²; the SU(4) factor is new information — report it per S.
**3c. Munich on multiplicity products:** for each product P ∈ {[2,1,1,0]², [3,2,1,0]², [4,2,0,0]², [4,2,1,0]², [4,3,1,0]²} and each irrep from
`./cgc_hw_bench decompose 4 l1 l2 l3 l4 l1 l2 l3 l4`: `./cgc_hw_bench 4 <λ> <λ> <μ_munich> 5 machine` (original) and `HW_MODE=eig ./cgc_hw_prof_ob ... 5 machine`; sum the
hw-solve times per product; for [6,4,2,0]² only the `eig`/`compact` variants and only if the memory estimate allows (the original needs 20 GB on its largest sectors) — a
per-irrep table with the sectors that could not be run marked.
**3d. Ours, same laptop:** from `Paul/Python/gaudin_gl4` with the venv and the single-thread environment of §1.4: `$PY -X utf8 slice_engine_gl4.py bench 12` (G1's protocol: min
of 5 cold runs, t_a/t_b per S, writes its table; ≈ 3 min) and `$PY -X utf8 slice_ops_gl4.py bench` (G3's families: [2a,a,0,0]² a ≤ 6 and [3a,2a,a,0]² a ≤ 2(3); check its docstring
for the exact stage name and the `ext` option; ≈ 5 min), plus the per-product t_a/t_b for the products of 3c (call `slice_ops_gl4`'s timing function on those products, or time
`so.triple_mats` + `so.copies` + `se.hw_from_slice` over all irreps yourself with `perf_counter`, cold caches, min of 3). Then one multithreaded repeat of `bench 12` with the thread
variables unset, reported separately.
**3e. The table.** `RESULTS-G7` §3: for S = 1..12: Munich orig (published July / re-run now), Munich eig, Munich compact, ours t_a, ours t_b (single-thread; multithread in a
second column), ratios orig/t_b and eig/t_b, max RSS of each Munich variant; for the multiplicity products: Munich orig / eig (sum over irreps) vs ours t_a / t_b, with #irreps,
Σm, max m, max n_W. Every number with its load average. Then the honest sentence for the notes: the ×860 of `benchmark_gl4.md` (two machines) becomes ×(…) on one machine, and the
repaired Munich is ×(…) — and whether the repaired code also hits the memory wall at S = 9 or 10.

## 4. D4 — report and hand-off

`RESULTS-G7-munich-comparison.md` (style of RESULTS-G1): environment; the test-suite run under Python 3.10 (any difference); D1 tables (per irrep: σ, support sizes, max deviation;
per multiplicity sector: m, singular values of Q, the Q matrix); D2/D3 tables; surprises; what could not be run and why (memory, timeout); the paragraph for the notes (§8 of
`gl4-slice-engine-notes` lists the literal Munich comparison as "not run" — give the replacement text; do not edit the notes). Raw outputs in `results/su4-same-machine/`
(gitignored? no — they are small text files; keep them). Printed output of every run in `munich_bridge_output.txt`. Add a block at the top of `HANDOFF.md`: "DONE (date,
Ubuntu): gl(4) Munich comparison — …" with the headline numbers and the file list; do not commit unless Paul asks (list the files to commit).

## 5. Acceptance for the whole package

All 42 multiplicity-free irreps `MATCH` at 1e-9; every multiplicity sector of the eight products in 2b has max |1 − σ_i(Q)| ≤ 1e-8 (≥ 100 sectors, Σm ≥ 170); the same-machine
table is complete for S ≤ 8 (orig, eig, compact) and S ≤ 12 (ours), with RSS and load recorded; the multiplicity-product table for the five products; the report. If a comparison
fails, do not stop: record the failing irrep/sector with the numbers (support difference, deviation), check the obvious causes (Munich label conversion, factor order site 1/site 2,
`alpha` count = m, Gram–Schmidt vs our copies), and report what you found — a genuine discrepancy would be the most important result of the package.

## 6. Pitfalls

* `hwcgc.sh` selects the HW rows by the LARGEST Q(M'') = dim(S''); for λ ≠ λ′ pass the factors in OUR order (site 1 = first argument = V_λ, site 2 = V_λ′): Munich's Q(M) is the
  first factor's index, Q(M') the second's; `to_munich` uses the same order.
* Munich's copies for m ≥ 2 are an orthonormal basis chosen by its own procedure (Gram–Schmidt in index order), not canonical; only the span is comparable. Our copies are
  G⊗G-orthogonal with distinct joint eigenvalues, so Q must be orthogonal; the ordering of our copies is by (ev₃, ev₄ᵃ, ev₄ᵇ) ascending (G3 §8).
* Munich prints 15 decimals: tolerance 1e-9, not tighter. One global sign per vector is arbitrary in both codes.
* Timing: the published `hw_timings.txt` was reps = 1 (noisy); min of 3–5 is the better statistic; say which you quote. Keep the two Munich binaries' outputs in separate files.
* Locale: `LC_ALL=C` for every awk/sort/printf; the `test.sh` table parsing silently mis-sums otherwise.
* Do not let a Munich run swap the machine to death: the memory formula in §1.4 and `ulimit -v`.
