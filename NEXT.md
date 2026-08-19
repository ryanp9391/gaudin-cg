# Next — gaudin-CG

*Updated 2026-08-15 · derived from `HANDOFF.md` (last session 2026-07-22)*

**Focus:** CG coefficients directly in the XXX spin chain (`Paul/Mathematica/XXX/`). The Gaudin
subproject proper is parked; the h→0 scaling limit is the live connection point.

1. **Extract explicit CG coefficients from the generating function** — the one piece needed to
   make the head-to-head against the C++ benchmark identical (C++ outputs explicit
   coefficients; GF timing is currently build-only). Also the natural next step now that τ2 is
   closed. `Paul/Mathematica/Gaudin/Clean/CGGenFun.wb`.
2. **Drop GF per-μ timings into a CSV** beside the C++ `Paul/Cpp/AlexCGC/results/` for a direct
   side-by-side. (Benchmark result so far: GF ~470× faster at S=8; C++ cannot run S=10 at all.)
3. **Promote `baxter_SU4_L2_XXX.wb`** (Casimir-form τ1, τ2, τ3 + gauge-fixed Baxter solver)
   from `Experiments/` to `Clean/`.
4. **Write up in `paper/xxxCG.tex`** — still single-site only; no L=2 factorization, scalar
   product, or CG-overlap dictionary written up yet.
5. **Push the factorization to s2≠0 / general L.**

## Housekeeping

- 12 uncommitted paths, including all of `Paul/Cpp/AlexCGC/` — deliberately uncommitted
  (Paul's choice: place, don't commit; the binary is machine-specific). Decide whether the
  sources should now be tracked.
- Kernel-state loss recurred all last session — re-run cells top-to-bottom, check `DownValues`
  with `?symbolName` rather than trusting a previous turn's definitions.
