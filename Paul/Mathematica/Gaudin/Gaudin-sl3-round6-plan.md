# sl(3) Gaudin, round 6 plan: beat the Munich CGC code with Q-function technology for general reps — [14,7,0]² as the benchmark (2026-10-21)

Status: DISPATCHED 2026-10-21 (five SONNET agents, round 6 = §3). Paul: "the goal of the project is to beat the Munich code using
Q-function technology for general reps, e.g. [14,7,0]² is a good example. Set up a complete plan of attack, then get cheaper sonnet
agents to implement." Rounds 1–5 are in `Gaudin-sl3-round{2,3,4,5}-plan.md`, HANDOFF.md and `Gaudin-sl3-CG-notes.md`.

## 0. The target, and what "beat" has to mean

**The problem.** For Λ = [λ₁,λ₂,0] (benchmark Λ = [14,7,0], dim 512, |λ|=21) compute the highest-weight (HW) Clebsch–Gordan vectors
Ψ₀ ∈ V_Λ⊗V_Λ for every irrep μ of the decomposition (113 irreps, 344 copies with multiplicity at a = 7 where Λ = [2a,a,0]), in the
GT basis, up to normalisation and sign, with the multiplicity copies in the SoV eigenbasis (eigenvectors of the leading-order Gaudin
Hamiltonian — the notebook's vertex operator V, equivalently the Laurent datum a₁ of the oper).

**The competitors and their complexity** (n_W = dimension of the weight-μ subspace of V_Λ⊗V_Λ; ~2 000 near the singlet at a = 7,
grows like |λ|⁴; m = multiplicity; w = dim V_Λ[μ−Λ] = a weight multiplicity of V_Λ, ≤ 8 at a = 7, ≤ 15 at a = 14):

| method | what it does | cost | a = 7, all 113 irreps | a = 8 | a = 14 |
|---|---|---|---|---|---|
| Munich (Alex–Kalus–Huckleberry–von Delft 1009.0437, `Paul/Cpp/AlexCGC`) | dense n_coeff × dim² array + full `dgesvd` null space | O(n_W³), memory n_W·dim² | 206 s (i5), OOM-free | OOM (> 12 GB) | impossible |
| Munich with the linear algebra fixed (`cgc_hw_prof`, HW_MODE=eig + OpenBLAS) | compact array, lowest eigenvectors of AᵀA | O(n_W²·…) | 10.9 s | feasible | no |
| our sparse null space + V (`Minimal_complete_working-fast.wl`, `fastEig`) | sparse Arnoldi on AᵀA, V on the HW space | O(n_W·nnz) | 12.7 s (i7) | ? | ? |
| **slice engine (new, §2.1; chain side, no integrability)** | PRV/Zhelobenko: HW vectors ↔ a w-dimensional slice; m×m eigenproblem of the reduced V | O(w³ + m·n_W) (n_W only to write the vector out) | expected ≪ 1 s | seconds | seconds |
| **Q-route (§2.2; the goal of the project)** | Q-system → Laurent data → jet formula → the slice | poly(|λ|), **independent of n_W** | to be measured | | |

Two consequences that shape the plan:
1. **"Beating Munich" on wall-clock at a = 7 is trivially achieved without integrability** (the slice engine; even fastEig ties it).
   So the Q-route is not measured against Munich's 206 s but against the real complexity floor: *the data that distinguishes the
   copies is w-dimensional*, and a method is only interesting if it produces it with cost independent of n_W. The Q-route must
   deliver the **slice** s_B ∈ V_Λ[μ−Λ] (w numbers, cf. §2.1) per branch B from Q-data alone; the full vector is representation theory
   shared by every method (and output-size bound).
2. **Where Munich dies, poly(|λ|) wins automatically:** at a ≥ 8 the Munich code is out of memory on a 15 GB laptop, at a = 14 the
   weight spaces have n_W ~ 10⁵ and no dense method exists. The benchmark that matters is therefore the *scaling* a = 7 → 14 of the
   Q-route and the slice engine, with the Munich column ending at a = 7.

**Deliverables of round 6** (all Sonnet-sized, parallel): the slice engine as the reference generator and complexity floor (R); the
Q-system spectrum at scale (Q0); the kernel at λ₁ = 14 (K); the ε-series solver at scale (S); the Gauss dictionary for general Λ (D).
Round 7 assembles them into `cg_from_Q` for general Λ, validates at [4,2,0], [6,3,0], [8,4,0], [14,7,0] against R, and benchmarks.

## 1. Where every ingredient of the Q-route stands (end of round 5 + today), and what breaks at λ₁ = 14

The formula (notes §5, §7, §9; constant exactly 1 on all 26 branches of three reps, N1/W5/S1):

```
Σ_ij G_ij(εω̂₁, εω̂₁′; θ₁₂) ν_i(θ₁; ε) ν_j(θ₂; ε) = ε^{N_B} · ĉ(φ, ω̂) · ⟨φ_u(ω̂)⊗e₁|Ψ₀⟩ / X_B + O(ε^{N_B+1}),
k_B = 2λ₁ − (μ₁−μ₃),  N_B = 4λ₁ − 2λ₂ − k_B,  G = (ω₁ω₁′)^{2λ₁−λ₂} Θ(1/ω₁, 1/ω₁′, 1/θ₁₂),
φ_u(ω̂) = e₁ᵀρ_Λ(n′),  n′ = exp(a′E₁₂) exp(g′E₁₃) exp(c′E₂₃)  (a′, g′, c′ explicit in φ, ω̂₁, ω̂₁′; Λ-independent),
ĉ = (φ₁−φ₂)^{2λ₂} Π_{ω∈{ω̂₁,ω̂₁′}} (φ₂−ω)^{λ₁−λ₂}(φ₃−ω)^{λ₁}.
```
The right-hand side is the generating function of the slice: ⟨φ_u(ω̂)⊗e₁|Ψ₀⟩ = Σ_ν [e₁ᵀρ(n′(ω̂))]_ν C_ν with C_ν = ⟨e_ν⊗e₁|Ψ₀⟩, ν running
over the weight space V_Λ[μ−Λ] (w components). X_B is one constant per branch (drops out of the direction; the Q-side norm Q2 fixes it).

| ingredient | status (round 5) | what scales / what breaks at [14,7,0] | round-6 WP |
|---|---|---|---|
| (a) the ε = 0 Q-system: branches, V-spectrum, untwisted Laurent data (H⁰, a⁰, b⁰) | **done today**: `pencilSolve` (linear pencil, exact GCD of two consistency polynomials); exact multiplicity m for every irrep tested; 0.05 s (a = 7, μ = Λ), 0.2 s (a = 14) | poly(|λ|); untested over ALL irreps and for the near-singlet sectors at a = 14 (largest Q-degrees) | Q0 |
| (b) ν-formulas (Plücker ratios at θ_s as polynomials in the Laurent data) | general-Λ algorithm `nu_general.frobenius/nu_symbolic` (N1), symbolic sympy | the Frobenius recursion to order λ₁+2 is poly(λ₁); the sympy *symbolic* route (expand in 10 variables) will not scale — needs a numeric/series evaluation of the same recursion | S (numeric ν on series) |
| (c) ε-series of (H_s, a_s, b_s) | Bethe-side degenerate Fredholm solver `bethe_series.py` (S1), exact in Q(√D), 0.1–2.4 s/branch for λ₁ ≤ 3; order 0 by Gröbner (did not finish for [3,2,0]) | order 0 must come from (a) (Gröbner is hopeless); exact Q(√D) arithmetic impossible for copies of multiplicity m (degree-m algebraic numbers) → high-precision floats; **depth**: N1 measured max ε-order 9 for [3,1,0] (N_B-ish, not the adjoint's 5) — at [14,7,0] N_B ≤ 42, so expect depth ~40: per order a 6×6 solve + series products, cheap *if* the equations are evaluated as truncated series, not as expanded polynomials | S |
| (d) the kernel Θ at λ₁ = 14 | Θ^{[λ₁,λ₂,0]} = block of the one-row master kernel (block theorem); master kernels exact for λ₁ ≤ 6 from ONE universal matrix N = (Y(1)Y₀⁻¹)⁻¹ of the one-row jet basis (M2/M3); closed form only through t⁵ | [14,7,0] needs the one-row λ₁ = 14 kernel on the 92 non-forced pairs per site (l ≥ 8): 92² entries, polynomials in (a, a′, t) with t-degree ≤ 4·13 − 11 = 41, a-degree ≤ 21, ~10⁷ coefficients; N on labels ≤ 15 is 14400×14400 — exact impossible, machine precision feasible (32 GB, 20 threads here) | K |
| (e) the dictionary: ĉ, φ_u = e₁ᵀρ_Λ(n′), the slice basis | adjoint only (`untwist_limit.gf_polys`: sympy exponentials of the 8×8 generators); N1/U2 solved the polynomial identity coefficientwise (overdetermined) | n′ is the same 3×3 unipotent for every Λ: ρ_Λ(n′) in the GT basis is cheap (nilpotent exponentials of GL3Rep generators); the overdetermined polynomial fit must become an m-point evaluation (the slice has w unknowns; we need only m+1 points ω̂_k) | D |
| (f) multiplicity structure: PRV slice space, reduced V | new today (§2.1) | makes (e)'s slice basis canonical and gives the n_W-free reference | R |
| (g) absolute normalisation / sign | Q2's general-Λ norm formula; sign is a phase convention | not needed to beat Munich (normalise the reconstructed vector); keep for the paper | round 7 |

## 2. Strategy

### 2.1 The slice engine (reference, complexity floor, and the canonical slice basis) — WP-R

Write Ψ₀ = Σ_q ψ_q ⊗ q over the GT basis q of site 2; the HW condition (E_i⊗1 + 1⊗E_i)Ψ₀ = 0 (E_i = E₁₂, E₂₃) read off at the site-2
component r gives E_i ψ_r + Σ_q (E_i)_{rq} ψ_q = 0, which determines the ψ_q of lower site-2 weight from those of higher weight, level by
level from the top ψ_top := ψ_{e₁} (site 2 = HW). The only constraints on ψ_top come from the chains f_i^k e₁ = 0 at k = ⟨Λ,α_i^∨⟩+1:
**E₁₂^{λ₁−λ₂+1} ψ_top = 0 and E₂₃^{λ₂+1} ψ_top = 0** — the Zhelobenko/PRV theorem: the multiplicity space of μ in V_Λ⊗V_Λ is
M_μ = {v ∈ V_Λ[μ−Λ] : E₁₂^{λ₁−λ₂+1}v = E₂₃^{λ₂+1}v = 0}, dim M_μ = m, inside the w-dimensional weight space, and v ↦ Ψ₀(v) is the
inverse of the slice map. Hence the hard part of the CG problem is a null space of size w (≤ 8 at a = 7). The vertex operator
V = ½Σ_{σ,τ∈S₃} sgnσ sgnτ e¹_{σ₁τ₁}e¹_{σ₂τ₂}e²_{σ₃τ₃} has ONE site-2 generator per term, so the slice of VΨ₀(v) needs ψ_q only for q
at most one lowering step below e₁ (which the first level of the recursion gives from v): **V_red : M_μ → M_μ is an m×m matrix
computable without ever touching n_W**; its eigenvectors are the SoV copies, its eigenvalues the pencil roots of (a). Writing the
full vector out is one triangular recursion (output-size cost). All of this is classical representation theory and is the honest
baseline; it also provides, for the Q-route, (i) reference slices for Λ where no other reference exists ([6,3,0], [8,4,0], [14,7,0])
and (ii) the canonical basis of M_μ in which the Q-route's m+1-point dictionary is solved.

### 2.2 The Q-route, reorganised for scale

Per branch B = (μ, copy): (a) ε = 0 Q-system → a⁰ (= V eigenvalue), H⁰, b⁰ (Q0); (c) ε-series of (H_s, a_s, b_s) to the needed depth,
degenerate solver with the order-0 point from (a) (S); (b) ν_i^{(n)}(θ_s) as series from the Frobenius recursion evaluated on series
(S); (d) the kernel pieces G^{[d]}_ij(ω̂_k) at m+1 numeric points ω̂_k (K); contraction L_B(ω̂_k) = Σ_{d+n+n′=N_B} G^{[d]}ν^{(n)}ν^{(n′)};
(e) solve L_B(ω̂_k)/ĉ(ω̂_k) = Σ_j x_j ⟨φ_u(ω̂_k)|b_j⟩ for the m coordinates x of the slice in the PRV basis b_j (D, R). Cost per branch
= O(depth² · (#pairs)² + depth · λ₁² · 6) + O(m · 92² · depth): poly(|λ|), no n_W anywhere. Open scaling question (round 7, §4): the
depth. Everything else is bookkeeping that already exists for small Λ.

## 3. Round-6 work packages (dispatched; one SONNET agent each; specs are the binding briefs)

| WP | agent / brief | language / kernel | deliverable | acceptance |
|---|---|---|---|---|
| **R** slice engine | `Paul/Python/gaudin_sl3/SPEC-slice-engine.md` | Python (gl3_reps + scipy.sparse) | `slice_engine.py`, `RESULTS-slice-engine.md`, `ref/cg-slices-{420,630,840,1470}.txt` (PRV bases, V_red, eigenvalues, slices, HW residuals) | reproduces W4's 45 exact branches (`ref/untwisted-exact-*.txt`) exactly; dim M_μ = multiplicity for every irrep of [14,7,0]² (Racah–Speiser); HW residual ≤ 1e-12 after reconstruction; timing table a = 1..14 |
| **Q0** Q-system at scale | `Paul/Mathematica/Gaudin/Experiments/Gaudin-Qsystem-scale-spec.md` | Mathematica, K1 (the ONLY agent on the kernel); `gl3-GT-tools.wl` + `Minimal_complete_working-fast.wl` | `Gaudin-Qsystem-scale.wl`, `Gaudin-Qsystem-scale.md`, `ref/qsystem-spectrum-{630,840,1470}.txt` (per irrep: pencil roots, untwisted Laurent data, Q polynomials) | #roots = multiplicity for all 113 irreps at a = 7 and all irreps at a = 10, 14; roots = V-spectrum of `fastEig` (a ≤ 7) and H⁰ = W4's eta/θ₁₂ closed form; timing vs a and vs μ |
| **K** kernel at λ₁ = 14 | `SPEC-kernel-L14.md` | Python numeric (numpy/scipy; mpmath for checks) on M2/M3's operator route | `kernel_L14.py`, `RESULTS-kernel-L14.md`, `kernel_N_L{8,14}.npz` (universal N, labels ≤ L+1) and `kernel_Theta_{840,1470}.npz` (non-forced block) | reproduces the exact kernels rep_300/400/500/600 to ≥ 1e-10 relative; the universality N^{(L)}|_{labels ≤ 7} = N^{(6)}; top-state sum rule and exchange rule on the new kernels; memory/time report |
| **S** ε-series at scale | `SPEC-bethe-series-scale.md` | Python (mpmath 60–100 digits), generalising `bethe_series.py` | `bethe_series_scale.py`, `RESULTS-bethe-series-scale.md` | exact agreement with `ref/laurent-series-{210,310,320}.txt` to order 10 (all 26 branches); [4,2,0] (multiplicity 3) all branches; [6,3,0] all branches; depth and time per branch vs λ |
| **D** Gauss dictionary for general Λ | `SPEC-gauss-dictionary.md` | Python (gl3_reps, sympy/mpmath) | `gauss_dictionary.py`, `RESULTS-gauss-dictionary.md` | reproduces `untwist_limit.gf_polys` for the adjoint exactly; the m+1-point dictionary recovers N1's 20 slices of [3,1,0]/[3,2,0] from N1's L_B; ĉ and φ_u for [14,7,0] in < 1 s |

Dependencies: none within the round (D uses the PRV basis from R only as an option; it can define M_μ itself). CPU: K is the heavy
one (14400² dense); S moderate; R, D light; Q0 on the kernel. ≤ 2 heavy processes per agent. No agent commits to git.

## 4. Round 7 (after the five return): assembly, benchmark, and the two research questions

- **E2E**: `cg_from_Q_general.py` = (a)+(c)+(b)+(d)+(e) for arbitrary Λ; validate every slice of [4,2,0]², [6,3,0]², [8,4,0]² against
  R; then [14,7,0]² (all 113 irreps); the benchmark table of §0 with measured numbers for all four methods at a = 7, 8, 10, 14.
- **Research 1 — the depth.** Measure D(Λ, B) = max ε-order of ν needed (N1: 9 for [3,1,0]); find the law (conjecture: bounded by the
  kernel's t-degree + N_B − …). If D ~ 4λ₁, the Q-route stays poly but with a large constant; look for a reformulation with O(1)
  ratios (normalise by Φ_∞ at a reference point instead of x₀^G — removes X_B but not the cancellations; test whether the dressed
  Wronskian variables a^{e_i}m_i of the one-site kernel reduce the depth; the exact algebraic structure of the series — W5's "H₁
  linear in ζ", the 2×2 and pencil char-poly closed forms — may give the needed Taylor coefficients without order-by-order solving).
- **Research 2 — the kernel in closed form** (round-5 WP-A continued): M3's normal-ordered exponential exp(ΣtⁿD_n) is closed through
  t⁵ with a tree-level Hamilton–Jacobi recursion; the loop terms from n = 4 are open. With K's numeric N at L = 14 as data, fit/derive
  D_6.. and the one-loop generating rule. Not needed for the benchmark (K's numeric kernel suffices) but needed for the paper's "formula".
- Then the paper (`paper/gaudinCG.tex` body, decision still Paul's).

## 5. Risks and fallbacks

- K: Y(1) at L = 14 may be too large/ill-conditioned numerically → do L = 8 (ample for [8,4,0]) and report the scaling; fallback for
  [14,7,0]: the state-independent K-tensor solve on the chain at the m+1 needed points (one n_W-sized solve per Λ, allowed as a stopgap,
  clearly labelled). The exact-rational universality check must be done at L ≤ 6 where exact N exists.
- S: the degenerate solver's branch selection (K-rational roots) has no numeric analogue → follow all candidates and prune by
  consistency at the next order; if several survive, use (a)'s data (a⁰ from the pencil, and H⁰ = eta/θ₁₂) to pick. Depth may exceed
  40 for near-singlet sectors at a = 7: report the time per order; do not try to be clever.
- R: the PRV identification is a theorem but the GT-basis conventions (gl3_reps vs the notebook vs W4) are a trap — the acceptance
  test against W4's exact HW vectors (all 45 branches) is mandatory before any large-Λ run.
- Q0: the near-singlet irreps at a = 14 have the largest Q-degrees (M ≈ 42); the exact GCD of two degree-~40 consistency polynomials
  over Q may be slow — fall back to numeric roots with a tolerance and report.
- D: ĉ may fail to clear denominators for some components (U2's Surprise 1) — those components vanish on HW kets; assert, don't fix.

## 6. Conventions for the agents (binding)

- Python: run from `Paul/Python/gaudin_sl3` with `python -X utf8`; import pattern of `gaudin_chain.py` (`sys.path` gets the script dir
  and `..` for `gl3_reps`). Existing scripts are read-only (import them; never edit). Caches `*.pkl`/`*.npz` are gitignored and
  rebuildable. Reference data for other agents go to `Paul/Mathematica/Gaudin/ref/*.txt` (TAB-separated, Wolfram InputForm numbers,
  a `#` header documenting every record — W4's files are the model). Each script gets `RESULTS-<name>.md` beside it (files, conventions,
  headline results with numbers, tables, surprises, not done) — the same style as the round-5 RESULTS files.
- Mathematica (Q0 only): the kernel K1 via the wolfbook MCP tools; load `Paul/Mathematica/Gaudin/Experiments/gl3-GT-tools.wl` then
  `Minimal_complete_working-fast.wl` (repo root) with `Get`; write a `.wl` script and a `.md` report in `Experiments/`; nobody else
  touches the kernel while Q0 runs.
- All: no git commits/pushes; ≤ 2 heavy processes; report failures with numbers, not adjectives; the GT basis is non-orthonormal
  (compare vectors up to scale, slices componentwise after fixing the scale on the first nonzero entry).

## 7. Decision points for Paul

1. Is the slice (w numbers per branch) an acceptable Q-route deliverable, with the full vector left to the universal recursion?
   (Recommended: yes — it is the only formulation in which "cost independent of dim V" is even possible.)
2. Should the Munich comparison in the paper be against the original (206 s) or the fixed (11 s) code? (Recommended: both, with the
   diagnosis; the scaling table a = 7 → 14 is the real argument.)
3. Whether to spend round 7 on the depth question before or after the end-to-end [14,7,0] run (recommended: run first, then measure).
