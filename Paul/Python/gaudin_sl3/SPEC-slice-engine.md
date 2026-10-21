# Spec R: the slice engine — PRV/Zhelobenko slices, the reduced vertex operator V_red, and the n_W-free reference generator for the HW Clebsch–Gordan vectors of gl(3), V_Λ ⊗ V_Λ (round 6, WP-R; 2026-10-21)

Delegated to a general-purpose subagent (SONNET), Python only (no Wolfram kernel, no wolfbook tools). Read first: `Paul/Mathematica/Gaudin/Gaudin-sl3-round6-plan.md`
§0, §1 (row (f)), §2.1, §3 (row R), §5, §6 (binding conventions of the round). Environment: Python 3.13 (`python -X utf8`), numpy, scipy, sympy 1.14, mpmath 1.3 (no gmpy2, no `uv`);
run from `Paul/Python/gaudin_sl3/`. ≤ 2 heavy processes; use `time.perf_counter()` (not `time.time()`, 15 ms resolution on Windows). **Never edit existing files; import them read-only.
No git, no notebooks.** Files you own: `slice_engine.py`, `RESULTS-slice-engine.md`, `slice_engine_output.txt`, `../../Mathematica/Gaudin/ref/cg-slices-{420,630,840,1470}.txt`.

## 0. Why this exists

The project's goal is to beat the Munich CGC code (`Paul/Cpp/AlexCGC`) with Q-function technology for general reps, benchmark [14,7,0]². The honest complexity floor is NOT the
n_W-sized null space that Munich (and our `fastEig`) solve, where n_W = dim of the weight-μ subspace of V_Λ⊗V_Λ (≈ 2000 at a = 7, ≈ 25 000 at a = 14), but the w-dimensional slice
(w = a weight multiplicity of V_Λ ≤ a+1, so ≤ 8 at [14,7,0]). This work package implements that classical-representation-theory engine: (i) the exact slice space M_μ per irrep,
(ii) the m×m matrix V_red of the Gaudin vertex operator on it (m = multiplicity) — its eigenvectors are the SoV copies — and (iii) the reconstruction of the full HW vectors by a
level recursion. It serves as the reference generator (slices of [6,3,0], [8,4,0], [14,7,0] exist nowhere else), as the canonical slice basis for the Q-route's dictionary (WP-D, round 7),
and as the benchmark baseline. Λ = [λ₁,λ₂,λ₃] with λ₃ = 0 in all runs (the code may allow λ₃ ≥ 0).

## 1. Conventions (binding)

* **One site**: GT states in the order of `gl3_reps.GL3Rep(Λ).patterns` (patterns are tuples of rows (singleton row, row of 2, Λ); the order is ascending lexicographic in (m₁,m₂,k) with the
  middle row (m₁,m₂) and the singleton k; **the highest-weight (HW) state is the LAST index**). Matrices act on kets as columns, E e_col = Σ_row M[row,col] e_row. E₁₂,E₂₃,E₂₁,E₃₂ are the
  `GL3Rep.generator` matrices (rational, non-unitary GT normalisation); **E₁₃ = [E₁₂,E₂₃], E₃₁ = [E₃₂,E₂₁]** (as `gaudin_chain.Chain.__init__` and the notebook). Weights are the
  (E₁₁,E₂₂,E₃₃) eigenvalues, α₁ = (1,−1,0), α₂ = (0,1,−1). Reuse `gl3_reps._raise/_lower/_weight/gt_patterns` (read-only) but note that `GL3Rep` builds dense sympy matrices:
  for dim up to 3375 (a = 14) build sparse exact matrices (Fractions) from `_raise`/`_lower`, dropping targets outside the interlacing range exactly as `GL3Rep._adjacent_up/_down` do,
  and test them against `GL3Rep.generator` on small reps.
* **Invariant metric**: the diagonal G with E_ijᵀ G = G E_ji, G[HW] = 1 (W4's `GMET1`; `nu_general.gram_diag`). The *orthonormal GT basis* is ẽ_q = e_q/√g_q (E-matrices there: Ê[r,q] = E[r,q]√(g_r/g_q),
  Êᵀ_ij = Ê_ji); do the floating-point recursion in it (conditioning) and convert back.
* **Two sites**: Ψ = Σ Ψ[ν,q] e_ν ⊗ e_q, **site 1 = first factor (ν, "slow"), site 2 = second factor (q, "fast")**; flat index (0-based) = ν·dim + q. Total generator E^{tot} = E⊗1 + 1⊗E.
  Site-2 HW state q_top = the HW GT pattern (last index). **Slice** of Ψ = ψ_top := Ψ[·,q_top] ∈ V_Λ (site 1), a vector of weight μ − Λ if Ψ has weight μ. W4's `cg` record is exactly this slice.
* **W4/notebook order** (files `ref/untwisted-exact-<l1><l2>0.txt`, parsed by `nu_general.load_w4`): descending lexicographic in (k,m₁,m₂) of the patterns {{k},{m₁,m₂},{l₁,l₂,l₃}} (index 1 = HW state);
  W4's two-site flat index is (i−1)·DIM + j (1-based, i = site 1, j = site 2). The map to the GL3Rep order is by pattern: `nu_general.w4_to_mine` (or match `GTBASIS`). Every comparison with W4 goes through it.
* **The vertex operator** (notebook `V`, `gl3-GT-tools.wl`, `Minimal_complete_working-fast.wl::vTriples`):
  **V = ½ Σ_{σ,τ∈S₃} sgn σ sgn τ  e¹_{σ₁τ₁} e¹_{σ₂τ₂} e²_{σ₃τ₃}**, e^s_{ij} = E_ij acting on site s, the product ordered as written (e¹_{σ₁τ₁} acts last), 36 terms. V commutes with the diagonal gl(3)
  action; its eigenvalues on the HW vectors of one irrep are the "V-spectrum" (= pencil roots of `pencilSolve`).
* Normalisation of copies: the slice is scaled so that **its first nonzero component (ascending GT index) is exactly 1**; the full vector carries the same scale. Eigenvectors are ordered by increasing eigenvalue.
  Signs/phases are conventions; compare vectors only up to scale (the GT basis is non-orthonormal).
* Racah–Speiser: multiplicity of V_μ in V_Λ⊗V_Λ = Σ_{w ∈ wt(V_Λ)} m_Λ(w) · sgn(s) over those w for which Λ + w + ρ is regular with s(Λ+w+ρ) = μ + ρ (ρ = (2,1,0)); m_Λ(w) = number of GT patterns of weight w.

## 2. The mathematics (plan §2.1, made precise)

**Theorem (Zhelobenko / PRV).** For a dominant μ, the map Ψ ↦ ψ_top is an isomorphism from the space of HW vectors of weight μ in V_Λ⊗V_Λ onto
**M_μ = { v ∈ V_Λ[μ−Λ] : E₁₂^{λ₁−λ₂+1} v = 0 , E₂₃^{λ₂−λ₃+1} v = 0 }** ⊂ V_Λ[μ−Λ] (weight space of dimension w ≤ min(λ₁−λ₂,λ₂)+1, = a+1 at most for Λ = [2a,a,0]); in particular dim M_μ = m = the Racah–Speiser multiplicity.

**Level recursion (the inverse map).** The HW conditions (E_i⊗1 + 1⊗E_i)Ψ = 0, i = 1,2 (E_1 = E₁₂, E_2 = E₂₃), read at the site-2 component r:  E_i ψ_r + Σ_q (E_i)_{rq} ψ_q = 0.
Label the weights of site 2 by (a,b) (here a,b ≥ 0 are weight labels, not the parameter of [2a,a,0]): wt(q) = Λ − aα₁ − bα₂, Q(a,b) = the GT states of that weight (|Q| ≤ min(λ₁−λ₂,λ₂)+1), X_{ab} = (ψ_q)_{q∈Q(a,b)} (each ψ_q ∈ V_Λ[μ − Λ + aα₁ + bα₂]).
For (a,b) ≠ (0,0):   **A_ab X_ab = −[ E₁₂^{site 1} X_{a−1,b} ; E₂₃^{site 1} X_{a,b−1} ]**,  A_ab = [ (E₁₂)_{Q(a−1,b),Q(a,b)} ; (E₂₃)_{Q(a,b−1),Q(a,b)} ] (stacked site-2 raising blocks).
A_ab has full column rank (a vector of V_Λ killed by E₁₂ and E₂₃ is the HW vector, and (a,b) ≠ (0,0)) and **does not depend on μ or on the slice**: one left inverse per (Λ,(a,b)).
So ψ_top = v determines everything, block by block in lexicographic (a,b); the overdetermined rows are consistent iff v ∈ M_μ (the theorem) — in exact arithmetic the leftover rows give exact consistency conditions
(T·RHS = 0) that must vanish; they also include the "dangling" conditions at the end of the α₁/α₂ strings (E₁₂ψ_r = 0 when Λ−aα₁ leaves the weights), i.e. the PRV constraints.
Total work O(n_W) per vector (n_W = Σ_{(a,b)} |Q| · dim V_Λ[…]); **the engine never solves an n_W-sized system**.

**V_red.** Raising generators annihilate Ψ, so projecting E_ab^{tot}Ψ = 0 (a<b, including E₁₃) on site-2 = q_top gives Σ_q (E_ab)_{top,q}ψ_q = −E_ab ψ_top; for a = b it is λ_a ψ_top; for a > b it is 0
(a lowering operator never produces the top state). Hence, for v ∈ M_μ, the slice of VΨ(v) is
**V_red v = ½ Σ_{σ,τ} sgn σ sgn τ E_{σ₁τ₁} E_{σ₂τ₂} Φ_{σ₃τ₃} v,   Φ_{cc} = λ_c,  Φ_{ab} = −E_ab (a<b),  Φ_{ab} = 0 (a>b)**  (only the slice and 24 of the 36 terms; **no level-1 ψ_q and no n_W**;
the plan's "level-1 ψ_q needed" is not necessary). V_red maps M_μ into M_μ (V preserves HW vectors); its matrix in the PRV basis is exact rational; its eigenvectors are the SoV copies;
its eigenvalues are real. **Cross-check requirement:** V (all 36 terms, site-2 operators as X·Eᵀ) applied to the reconstructed full vectors must equal Σ_j (V_red)_{ji}Ψ_j on ALL entries.

## 3. Tasks and acceptance (every number goes into RESULTS-slice-engine.md)

**R1 `prv_space(lam, mu)`** — exact basis of M_μ: nullspace over ℚ of [E₁₂^{λ₁−λ₂+1}; E₂₃^{λ₂−λ₃+1}] on the weight space V_Λ[μ−Λ] (w ≤ 15: Fractions are fine; skip a constraint whose image weight does not exist).
Canonical basis: RREF with columns in ascending GT index, one basis vector per free column (1 there, 0 at the other free columns) — then the coordinates of any v ∈ M_μ are its free entries.
Write the tiny Racah–Speiser routine (`racah_speiser(lam)`; check Σ m dim V_μ = (dim V_Λ)²). **Acceptance:** dim M_μ = RS multiplicity for EVERY irrep of V_Λ⊗V_Λ (and for all dominant μ = Λ + w) for
Λ = [2,1,0], [3,1,0], [3,2,0], [4,2,0], [6,3,0], [8,4,0], [14,7,0]: table (#irreps, Σm, max m, max w, mismatches = 0); expected: [14,7,0] has 113 irreps, 344 copies, max m = 8, max n_W = 2080.

**R2 `hw_from_slice(lam, mu, v)`** — the full HW vector(s) by the level recursion (floats in the orthonormal basis, left inverses of A_ab by SVD, numpy/scipy; an exact Fraction twin for small reps is wanted:
exact HW vectors and the exact consistency conditions). **Acceptance:** relative HW residual ‖(E_i^{tot})Ψ‖/‖Ψ‖ ≤ 1e-12 (i = 12, 23, 13) and total weight exactly μ for all irreps of all small reps and for ≥ 10
irreps of [14,7,0] including the singlet (14,14,14), (15,14,13), (28,14,0), (21,14,7); the exact twin gives exactly 0 residual; a slice NOT in M_μ must be detected (exact consistency condition fails).

**R3 `V_red(lam, mu, basis)`** — the exact m×m matrix (column j = V_red b_j in PRV coordinates; it must stay inside M_μ — assert via the coordinates), computed from the slice only. Cross-check against the
full 36-term V applied to the reconstructed vectors (all small reps, ≥ 5 irreps of [14,7,0]); eigen-decomposition with mpmath (≥ 40 digits): `copies()` returns eigenvalues (ascending) and the copies' slices.

**R4 Validation against W4's exact files** `ref/untwisted-exact-{110,210,220,310,320,410}.txt` — **all 45 branches**. (i) multiplicity; (ii) W4's exact `HWBASIS` kets: their slices lie in M_μ and the exact reconstruction from the
slice reproduces the whole ket EXACTLY; (iii) V_red eigenvalues vs W4's `alpha` — **find and state the exact relation** (hint: LAURMAT's A1 is the matrix of the double-pole datum α₁ on the HWBASIS kets; compare
A1 with V_red in that basis; a shift by a μ-dependent scalar is expected; check it against `eta`, `cas2` in the SECTOR record) and the copy ordering convention (copy 1 = smaller α₁); (iv) copy slices vs `cg` up to scale;
(v) full vectors vs `vec` up to scale (float recursion at 1e-12, and at ≥ 30 digits through the exact HW vectors). **Acceptance: 0 mismatches at 1e-12.**

**R5 Benchmark** (i7-12700KF, 32 GB; `perf_counter`; min of 5 cold runs, caches cleared): Λ = [2a,a,0], a = 1..14: (a) the n_W-free part = exact GT generators + RS + PRV + V_red (exact) + float eigendecomposition
for all irreps; (b) additionally the metric/float matrices/left inverses + the full reconstruction of all copies' vectors (batched over the copies of an irrep). Table: a, dim V_Λ, #irreps, Σ copies, max w, max n_W,
t_a, t_b (+ breakdown, entries written, worst HW residual). Context (HANDOFF.md, a = 7): Munich original 206 s; Munich with the linear algebra fixed 10.9 s; our sparse null space `fastEig` 12.7 s.

**R6 Reference files** `Paul/Mathematica/Gaudin/ref/cg-slices-<l1><l2>0.txt` for [4,2,0], [6,3,0], [8,4,0], [14,7,0]: TAB-separated, `#` header documenting every record and convention (the GL3Rep order and the map to the
notebook order, the definition and normalisation of V, the relation of its eigenvalue to W4's α₁, script name, date). Records: `REP`, `DIM`, `GTBASIS` (patterns in the GL3Rep order, W4 format), `WEIGHTS1`, `GMET1`, `NSECTORS`;
`SECTOR mu m w nW dimVmu cas2mu eta1`; `WSPACE mu {nu,...}` (1-based indices of the weight space V_Λ[μ−Λ]); `PRVBASIS mu {{nu,val},...}` (one line per basis vector, k = 1..m in line order, exact rationals, nonzero entries); `VRED mu {{...}}` (m×m exact, rows; column j = image of basis vector j); `CHARPOLY mu {c0,...,cm}` (exact characteristic polynomial of V_red, ascending powers, c_m = 1 — its roots are the V-spectrum that the Q-system's pencil must reproduce);
`COPY mu k eigenvalue {{nu,val},...}` (30 digits, slice with first nonzero component 1, k by increasing eigenvalue); `HWRESID mu k residual` (relative HW residual of the reconstructed copy, Python float syntax);
for the small reps also `VEC mu k {{idx,val},...}` (full copy vector, 30 digits, idx = 1-based flat index (ν−1)·DIM + q, both 1-based in the GL3Rep order, same scale as the COPY slice).

**R7 Independent re-verification (added by the first implementation; recommended):** (i) re-parse the written files from the TEXT and check them with code paths independent of the engine (`GL3Rep` float matrices for the chain conditions, the 36-term V and the HW conditions in the plain basis; `gaudin_chain.Chain` two-site operators for the VEC records);
(ii) a Munich-like brute-force null space (eigen-decomposition of AᵀA of the total E₁₂,E₂₃ on the weight-μ subspace, n_W columns) for every irrep of [4,2,0], [6,3,0], [8,4,0], [14,7,0]: null dimension = RS multiplicity, our copies lie in it and span it, V restricted to it has our eigenvalues;
(iii) the exact twin on ALL PRV basis vectors of [2a,a,0] (a = 7, 14): exact HW vectors, exact consistency conditions, exact HW check (E₁₂,E₂₃,E₁₃).

**Expected facts the implementation should confirm (so a re-implementation can check itself):** the metric equals W4's `GMET1` exactly; dim M_μ = RS multiplicity everywhere; the full 36-term V on the reconstructed vectors equals V_red to 1e-15;
α₁(W4) = eigenvalue(V_red) + η₁/2 with η₁ = |Λ|² − Ω = |Λ|² + cas2(Λ) − cas2(μ)/2 (|Λ| = λ₁+λ₂+λ₃, Ω = Σ_ij E¹_ij E²_ji = (cas2(μ) − 2cas2(Λ))/2, cas2(μ) = Σ_i μ_i(μ_i + 4 − 2i) = W4's `cas2`), i.e. `LAURMAT` A1 = V_red + η₁/2 exactly in W4's HWBASIS slice basis, for all 45 branches. Further structure found by the first implementation (cheap checks for a re-implementation): V_red, as the formula-defined operator on the whole weight space V_Λ[μ−Λ], is tridiagonal in the GT basis (every irrep tested, up to [28,14,0]), with strictly positive off-diagonal products in the w = m sectors (⇒ simple real spectrum there); V(singlet of [2a,a,0]²) = 2a²(2a+1) (a = 1..5, 7), V(top irrep (2λ₁,2λ₂,0)) = −λ₁λ₂; the first non-null eigenvalue of AᵀA of the total E₁₂,E₂₃ on a weight-μ subspace (orthonormal basis) is exactly 1. Precision note: the float64 recursion's HW residual grows like 2.2^a (1.6e-14 at a = 7, 5.3e-12 at a = 14); the exact recursion rounded to double stays at ~1e-15 (`hw_from_slice_rounded`).

## 4. Reporting

`RESULTS-slice-engine.md` (style of the round-5 RESULTS files): files, conventions, headline results with numbers, the tables of R1 and R5, the W4 validation numbers per rep, what surprised you, what is not done. Save the printed output of the runs
to `slice_engine_output.txt`. End with a concise summary (what was verified against W4 with numbers, the benchmark table, surprises). Report failures with numbers; if a test fails, do not paper over it.
