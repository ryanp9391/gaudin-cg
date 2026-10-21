# Spec: Stage B of `Two-site-rotation-brackets-plan.md` — the two-site K tensor, (c,c′)-blocks, det₆ regression, and the bracket count of ⟨Ω_w⊗Ω_{w′}| (2026-10-11)

Delegated to a wolfbook-builder subagent (sonnet); planned and reviewed by the main session. **New notebook:**
`Paul/Mathematica/new-conventions/L2-two-site-rotation-brackets.wb` (self-contained, rep-generic with `Λ` as the only
rep switch; first code cell `Quiet[ClearAll["Global`*"]]`; run [2,1,0] first, then [3,1,0]). Create it with
`wolfbook_newNotebook`; if cells of the new file will not execute natively, run the sources in K1 and store printed
results as markdown cells under each code cell (say so).

Read first: the plan §1, §3, §6, §7; `L2-site2-rotation-brackets-spec.md` and the Stage A result cell of
`L2-site2-rotation-brackets.wb` (§0 below summarises what Stage A found); `L2-210-site1-rotation-brackets.md`;
`Two-rotations-det-overlap.md` §1 (det₆, Ωc, brackets) and §4 (Part B: Ω_w⊗Ω_{w′} is not a det₆, full-rank 8×8
coefficient matrices — the statement this stage replaces by a block-rank statement). Cell-source dumps (scratchpad,
UTF-8): `ref-general-all.txt` (`L2-site1-rotation-brackets-general.wb`, 63 cells — the machinery: cells 3–6 setup,
9–10 SoV basis/labels/monomials, 13 K map, 16 tensors + Plücker tools + `decompForm` + `Bk2`, 19–21 loci, 24–25
reduced forms, 28 reach), `ref-rotation-sov-all.txt` (`L2-rotation-sov.wb`: cell 8 det₆ with `rowG`/`ins`/`wgt`
on all six poles and `Fvac` = L=2 w-vacuum direct Q's; cell 25 `Ωc`, `PppW`, the one-site vacuum direct Q's `QDv1`;
cell 82 the w′-vacuum Q's and the per-ladder bra data; cell 84 the 8×8 coefficient matrices `Mmat`), and the Stage A
notebook's own cells (dump `ref-site2-all.txt` if provided, else read it with `wolfbook_getNotebookContext`).

## 0. What is known (inputs to this stage)

- Conventions as in the general notebook: companion frame, L = 2, h = 1/2, θ = (2/7, −3/5) (t = 62/35), z = (3/2, 5/4,
  2/3), z′ = (3/2, 5/4, 3/7); `prec` = 60; SoV labels (a|bc) per site with roots θ_s + {−1−a, −b, −c}h; monomials
  M^{(s)}_{(a|bc)} = minor_s(−1−a, −b)·v₃(θ_s−ch), minor_s(n,n′) = v₂(n′)v₃(n) − v₃(n′)v₂(n) on ladder s (rows {2,3},
  k = 3), normalised by the vacuum (λ₁|λ₁λ₂); ⟨x|Ψ⟩/⟨x₀|Ψ⟩ = N_x M^{(1)}M^{(2)}/(vac). Rotation g =
  exp(αE₁₂)exp(βE₁₃)exp(γE₂₃), first row (1, α, β+αγ). Ladder brackets ⟨Q^jF⟩_s, Bk₂^{(s)}[F,G].
- One-site vacuum (companion frame): Ωc[{w₁,w₂}] = e₁ᵀ exp(E₁₂/w₁) exp(−E₁₃/(w₁w₂)) exp(((w₁+w₂)/(w₁w₂))E₂₃) =
  gaussCov[{1/w₁, −1/(w₁w₂), (w₁+w₂)/(w₁w₂)}], on the surface β+αγ = α² (plan §0; verified exactly in Stage A).
  The [3,1,0] best site-1 locus IS this surface (5 brackets; generic 7). Ωc[w]⊗Ωc[w] ∝ det₆[F_vac,w]/Ψ_Q (10⁻⁵⁰).
- Site-1 rotation ⊗ e₁ᵀ: support = site-2 labels (λ₁|λ₁ c₂), factors λ₂^{(c₂)} = v₃(θ₂−c₂h)/v₃(θ₂−λ₂h); counts
  5/7 generic, 4/5 reduced ([2,1,0]/[3,1,0]).
- **Stage A findings (`L2-site2-rotation-brackets.wb`, 36 cells, result cell 36; dump `ref-site2-all.txt`):** the SoV
  construction is strongly site-asymmetric. e₁ᵀ⊗e_νᵀ (site-2 rotation) involves ALL dim V site-1 labels on a triangular
  set of label pairs (35 of 64 states for [2,1,0], 103 of 225 for [3,1,0]; the vacuum label on either site pairs with
  every label on the other); only 15 / 29 of those pairs have x₁ = (λ₁|λ₁c₁); no additive weight rule (SoV covectors
  are not weight eigen-covectors). Per-label automatic count 17 / 40 (vs 5 / 7 at site 1); in the block form of §1 below
  (blocks A^{(c₁,c′)}, matrix-rank decomposition) the ν₁ = 1 slice has block ranks 3,3,3,2 (sum 11, 15 products) for
  [2,1,0] and 6,6,6,4 (sum 22, 35 products) for [3,1,0]. The vacuum family Ωc[w] does NOT reduce the site-2 count; the
  only cheap site-2 locus is b₃ = β+αγ = 0 (13 / 28, not CG-complete); the CG-complete reducing loci reduce by 1–2 only.
  The adjoint site-2 conic constant is cc₂(t) = (t+2)/(2(t+1)) = cc₁(1−t) with cc₁(t) = (t−3)/(2(t−2)) the site-1
  constant (exact over 15 (h,θ) points). Everything z-independent (K identical at z′, z″). Also verified exactly: Ωc[w] =
  gaussCov[{1/w₁, −1/(w₁w₂), (w₁+w₂)/(w₁w₂)}]; the adjoint quadric = {Ωc[w]ρ(D) : D torus, d₁d₃/d₂² = cc₁}.
- **Reading of the asymmetry (main session, to test in Task 2b):** the reversed-order chain t̃(u) = tr_a[ℒ^{(2)}(u−θ₂)
  ℒ^{(1)}(u−θ₁)G] is conjugate to t(u) by the Yangian R-matrix on V⊗V: ℒ^{(2)}ℒ^{(1)} = R T R⁻¹ entrywise in the auxiliary
  indices, R = R^{ΛΛ}(θ₁−θ₂) ∈ End(V⊗V), gl(3)-invariant, so t̃ = R t R⁻¹, Ψ̃ = RΨ (same Q's), B̃ = RBR⁻¹, ⟨x̃| = ⟨x|R⁻¹
  and the SoV wavefunctions ⟨x̃|Ψ̃⟩ = ⟨x|Ψ⟩ are identical. In the reversed chain the monodromy-first site is tensor
  factor 2, so by the mirror symmetry of the construction e₁ᵀ⊗g₂ is the cheap one there: ⟨e₁ᵀ⊗g₂|Ψ̃⟩ = ⟨(e₁ᵀ⊗g₂)R|Ψ⟩
  should have support only on site-1 labels (λ₁|λ₁c₁) with the site-1 counts 5 / 7 and the adjoint conic constant
  cc₁(−t) (θ₁₂ → −θ₁₂). The t → 1−t of the undressed covector would then be the R-matrix shift. (e₁ᵀ⊗e₁ᵀ)R ∝ e₁ᵀ⊗e₁ᵀ
  (HW vector of the top irrep), so x₀ is unchanged. If this holds, the cheap site-2 family is the R-dressed one,
  ⟨(e₁ᵀ⊗Ωc[w′])R(θ₁₂)|, and the natural two-site candidate is ⟨(e₁ᵀ⊗g₂)R(ρ(g₁)⊗1)| (items (g), (h) of Task 3).

## 1. Definitions (fix these conventions and keep them for Stages C/D)

- **Two-site K tensor.** For the dim V² GT covectors e_{ν₁}ᵀ⊗e_{ν₂}ᵀ (index (ν₁−1)d₁+ν₂, i.e. `kron[uv[ν₁],uv[ν₂]]`),
  e_{ν₁}ᵀ⊗e_{ν₂}ᵀ = Σ_m K_{ν₁ν₂;m} ⟨x_m|/N_m in normalised-wavefunction units, so that
  ⟨Φ|Ψ_A⟩/⟨x₀|Ψ_A⟩ = Σ_m φ_m · M^{(1)}_{x₁(m)}M^{(2)}_{x₂(m)}/(M^{(1)}_vac M^{(2)}_vac), φ_m = Σ c₁_{ν₁}c₂_{ν₂}K_{ν₁ν₂;m}
  for Φ = c₁⊗c₂. One `LinearSolve[Transpose[vecB], ·]` with all d₂ unit covectors as right-hand sides (the general
  notebook's cell 13 with d₂ instead of d₁ RHS), times N_m. Rationalise (`ratz`, tolerance 10⁻⁴⁰), report the error
  and the max denominator digits; keep the rationalised K as an exact `SparseArray`.
- **Blocks.** With x₁ = (n₁,n₁′; −c) (minor (n₁,n₁′) on ladder 1, lone point −c) and x₂ = (n₂,n₂′; −c′): for Φ,
  A^{(c,c′)}(Φ) := the nMin × nMin matrix with entry [(n₁n₁′),(n₂n₂′)] = φ_m at the state m with those labels (0 if
  no such label), minor order `mnL` of cell 16 on both sites. Then
  ⟨Φ|Ψ⟩/⟨x₀|Ψ⟩ · M^{(1)}_vac M^{(2)}_vac = Σ_{c,c′} v₃(θ₁−ch) v₃(θ₂−c′h) Σ A^{(c,c′)}_{(n₁n₁′),(n₂n₂′)} minor₁ minor₂.
  Verify this identity on all states for every Φ you use (10⁻⁴⁰ or better).
- **Counts.** For a block A (an element of Λ²(ladder-1 points) ⊗ Λ²(ladder-2 points), as a matrix): r = MatrixRank
  (the "separable" count: number of terms (2-form on ladder 1) × (2-form on ladder 2)). A product Bk₂^{(1)}·Bk₂^{(2)}
  is a wedge⊗wedge, matrix rank 1, so the number P of bracket products needed satisfies **P ≥ r** (exact lower bound)
  and **P ≤ Σ_k (rank ω_k/2)(rank ω′_k/2)** for any rank decomposition A = Σ_k ω_k⊗ω′_k (upper bound; wedge-decompose
  each ω_k, ω′_k with `decompForm`). Compute the lower bound r and the upper bound for the obvious decomposition (row
  space basis by column pivots; also try the transpose-side pivots and take the smaller) and for the greedy variant
  (repeatedly strip the rank-1-in-Λ²⊗Λ² term with the largest … — only if cheap; a tensor-rank optimisation is NOT
  asked). Report both bounds per block and their totals; "the count" = the upper bound, flagged exact when it equals r.

## 2. Tasks (each its own code cell(s) + a markdown result cell with a `Grid` of checks; save after each task)

**Task 0 — setup and SoV basis.** General cells 3–6, 9–10 (regression: label rule, monomials, `Nm`); `Ωc`, `gaussCov`;
the det₆ machinery, rep-generic as in `ref-310-two-rotations-all.txt` cell 17 (copy it: poles {s,n}, n = λ₃..λ₁ per
ladder (3 per ladder for the adjoint, 4 for [3,1,0]), `wgt`, insertions 1 and e^{2πix/h}, `rowG`, `det6`, Ψ_Q with the
vacuum triple (λ₁+1, λ₁, λ₂) i.e. (3,2,1) / (4,3,1)) with the L=2 w-vacuum direct Q's (that cell's pattern: Ω⊗Ω is the
exact eigen-covector of t_a twisted by the companion matrix G(w), `solveQ` on the direct Baxter equation with w in
place of z; degrees (0,2,4) adjoint, (0,4,6) [3,1,0]); w = (5/2, 7/3, 3/5), w′ = (7/4, −2/5, 9/7). VERIFY
Ωc[w]⊗Ωc[w] ∝ det₆[F_vac,w]/Ψ_Q on all states for the rep (10⁻⁴⁵ or better) before relying on it. Time the full setup.

**Task 1 — the two-site K tensor.** Build K (d₂ × d₂, rows ν = (ν₁,ν₂), columns m) as in §1; checks: expansion error;
rationalisation error; denominators; the slice ν₂ = 1 reproduces the general notebook's `Kg[c₂]` (rebuild them here
from the slice and compare with the general notebook's values — recompute, don't copy numbers); the slice ν₁ = 1
reproduces Stage A's site-2 K maps; z-independence (rebuild at z′: identical rationals). Rank of the full K (expect
d₂ = full). Time it.

**Task 2 — supports.** For Φ = c₁⊗c₂ with c₁ = gaussCov[g₁], c₂ = gaussCov[g₂] generic rational: the set of (x₁,x₂)
label pairs with φ_m ≠ 0; the site-1 label set as a function of x₂ and vice versa; is the support the product of the
two one-site supports {(x₁, (λ₁|λ₁c₂))} ∪ {((λ₁|λ₁c₁), x₂)} or larger (expect larger: all 225 for generic g₁, g₂;
report the count and the structure — e.g. which (c,c′) blocks are populated and the number of populated slots per
block). Same for e₁ᵀ⊗c₂ and c₁⊗e₁ᵀ (must reproduce the one-site supports) and for Ωc[w]⊗Ωc[w].

**Task 2b — the R-matrix dressing (tests the reading of the asymmetry in §0; exact; do it before Task 3 so items (g),
(h) are available).** Build the reversed-order monodromy T̃(u) = ℒ^{(2)}(u−θ₂)ℒ^{(1)}(u−θ₁) (same `lax` code, factors
swapped) and solve for R ∈ End(V⊗V) with R·T_{ab}(u) = T̃_{ab}(u)·R for all nine auxiliary entries (a,b) at two rational
u (exact sparse linear system). Reduce the unknowns first: R commutes with the diagonal gl(3) action, so it is
block-diagonal in the total-weight decomposition of V⊗V (unknowns = Σ_weights (dim)²; adjoint ≈ a few hundred, [3,1,0]
≈ a few thousand — fine for `NullSpace` on a sparse rational matrix; if too slow, impose commutation with the Cartan
and with E₁₂+E₂₃ type generators first, or restrict to one weight block at a time). Expect a 1-dim null space (report
the dimension; if > 1 add a third u). Checks: R commutes with ρ(g)⊗ρ(g) for a rational g (exact); t̃_a(u) = R t_a(u) R⁻¹
at one u; (e₁ᵀ⊗e₁ᵀ)R ∝ e₁ᵀ⊗e₁ᵀ; R(θ₁₂) at θ₁₂ → ∞ ∝ permutation-free (normalise R so that its action on x₀ is 1).
Then the test: for Φ = (e₁ᵀ⊗gaussCov[g₂])·R and ·R⁻¹ (try both; also the two-sided orderings if neither works — report
which one is right), the support in the SoV basis (expect only site-1 labels (λ₁|λ₁c₁), with the factors
λ₁^{(c₁)} = v₃(θ₁−c₁h)/v₃(θ₁−λ₂h)), the per-label automatic count (expect 5 / 7), and for the adjoint the conic locus in
g₂ with its constant (expect cc₁(−t) = (t+3)/(2(t+2)) = 229/342 at t = 62/35; also at (h,θ) = (1/3,(1/5,−4/9))) and for
[3,1,0] the best locus (expect β+αγ = α², count 5). If the expectation fails, record what the dressed support is and
move on — Task 3 does not depend on it.

**Task 3 — blocks and counts at numeric points (the core).** For each Φ in the list below: all blocks A^{(c,c′)},
c, c′ = 0..λ₂; per block: populated slots, matrix rank r, upper-bound count, the 2-form ranks of the decomposition
factors; totals (lower bound Σr, upper bound ΣP); the overlap identity check. List:
 (a) generic pair (g₁, g₂) = ((2/3, −3/7, 5/4), (−1/2, 4/9, 3/5)); (b) a random pair of covectors (control);
 (c) Ωc[w]⊗Ωc[w] (the det₆ vacuum) and Ωc[w′]⊗Ωc[w′];
 (d) Ωc[w]⊗Ωc[w′] and Ωc[w′]⊗Ωc[w] (the headline family at one point);
 (e) Ωc[w]⊗e₁ᵀ and e₁ᵀ⊗Ωc[w′] (regression: site-1 count 4/5, site-2 count from Stage A — the blocks with c′ = λ₂
     resp. c = λ₂ only, and the block ranks must reproduce the one-site column 2-form ranks: say exactly how);
 (f) gaussCov[g₁]⊗e₁ᵀ with g₁ generic (regression: 5/7);
 (g) the R-dressed site-2 vacuum (e₁ᵀ⊗Ωc[w′])R (with the R and the side found in Task 2b) — expect the site-1-type
     block structure (blocks c = λ₂ only? report) and count 4/5 if Ωc lies on the dressed best locus;
 (h) the R-dressed two-site candidates (e₁ᵀ⊗Ωc[w′])R(ρ(Ωc-rotation g₁)⊗1) = (e₁ᵀ⊗gaussCov[g₂])R(ρ(gaussCov-type
     g₁)⊗1) with g₁ = the Gauss parameters of Ωc[w], and the other order (gaussCov[g₁]⊗e₁ᵀ)R(1⊗ρ(g₂)): block ranks and
     counts, compared with (d) Ωc[w]⊗Ωc[w′] — the question is whether dressing by R makes a two-rotation covector cheap.
     ρ(g) on V is `nexp`/`MatrixExp` of the Gauss factors as in `gaussCov`; (ρ(g₁)⊗1) is `kron[ρ(g₁), Id]`.
**det₆ regression (essential, plan §7):** for (c), derive the expected structure first: det₆ = (e₂−e₁)^{…} × the
per-ladder 6×6 with rows (s,(i,j)), whose Laplace expansion along the three ladder-1 rows is Σ over the 20 3-subsets of
columns of (3×3 ladder-1 minor) × (complementary 3×3 ladder-2 minor); each 3×3 per-ladder minor is a polynomial in
the per-ladder brackets, hence in our monomials a sum over the lone point c of (2-form)_c × v₃(θ_s−ch). So the vacuum's
blocks satisfy A^{(c,c′)} = Σ_{S} ω^{(c)}_S ⊗ ω′^{(c′)}_{S^c}: r_{cc′} ≤ 20 and the total P is bounded by the Laplace
count. **Compute the actual r_{cc′} of the vacuum and compare with the generic pair (a):** the vacuum must have visibly
lower block ranks (it is a single determinant); if it does not, the block bookkeeping (minor signs, lone-point
conventions on site 2, the K normalisation) is wrong — fix before Task 4. Also check: is Σ_{c,c′} r_{cc′} for the
vacuum ≤ the number of nonzero Laplace terms, and does the per-block structure match the (2-form)_c ⊗ (2-form)_{c′}
decomposition you get by expanding det₆ directly (optional, only if quick: build the 20 Laplace terms numerically from
the bracket rows of cell 8 and compare the resulting Λ²⊗Λ² tensor with A^{(c,c′)} up to the overall constant).
**First decisive number (plan §6):** at point (d), is the total count the product of the one-site counts (4×4 = 16 for
the adjoint, 5×5 = 25 for [3,1,0], with the site-2 count from Stage A if different), lower, or higher? And how does it
compare with the vacuum (c) and the generic pair (a)?

**Task 4 — the family Ωc[w]⊗Ωc[w′] over its 4 parameters (w₁,w₂,w₁′,w₂′).** Blocks symbolic: φ_m is bilinear in
(c₁,c₂), c_s = gaussCov[1/w_s1, −1/(w_s1 w_s2), (w_s1+w_s2)/(w_s1 w_s2)] rational in w; so A^{(c,c′)}(w,w′) is a matrix
of rational functions — clear denominators (multiply by the obvious monomials in w₁w₂, w₁′w₂′) to get polynomial
blocks. (i) Generic ranks r_{cc′}(w,w′) (rank over the function field: evaluate at 3 random rational points, take the
max); (ii) rank drops: for each block, the ideal of (r×r minors) — `Factor` a few minors first and look for the obvious
factors (w₁−w₁′, w₂−w₂′, w₁−w₂′, w₂−w₁′, w₁−w₂, …); then the sub-families: w′ = w; w₁′ = w₁ (common w₁); w₂′ = w₂;
w₁′ = w₂; w₂′ = w₁; w₂ → ∞ (site-1 rotation limit, c₁ → gaussCov[{1/w₁, 0, 1/w₁}]? — compute the limit of Ωc, don't
guess) and w₂′ → ∞: rank tables and counts on each; (iii) the minimal count over the family and where it is attained;
whether the count is a product (site-1 count) × (site-2 count) on each sub-family. Gröbner only if the minors are small;
otherwise rank tables on the sub-families (say which).

**Task 5 — one-site-fixed slices (plan §3.5).** g₂ = Ωc[w′] fixed (rational), g₁ = gaussCov[α,β,γ] free: the loci
in (α,β,γ) where some block rank drops below its generic value (minors, `Factor`, dimension, genuine content) — in
particular whether the site-1 vacuum family β+αγ = α² is the only locus (then the two sites decouple at the bracket
level and the count is a product), or whether cross terms appear (loci depending on w′); and the best site-1 count
at fixed g₂. Then vice versa (g₁ = Ωc[w] fixed, g₂ free). For the adjoint, do the quadric 2(t−2)(β+αγ) = (t−3)α²
and the vacuum family both appear as site-1 loci when site 2 is Ωc[w′]?

**Task 6 — markdown summary** `### two-site rotation brackets (Stage B) — result (2026-10-11)`: per rep the K tensor
facts (size, errors, slices reproduce one-site maps, z-independence), the support structure, the block-rank tables
(lower/upper bounds) for (a)–(f), the det₆ regression outcome with the expected vs actual ranks, the decisive
comparison "Ωc[w]⊗Ωc[w′] count vs product of one-site counts vs vacuum vs generic", the 4-parameter family's rank
drops and sub-family counts, the one-site-fixed loci, timings, surprises, and a recommendation for Stage C (which
family and which blocks to write in closed form). Save. Do not write HANDOFF.

## 3. Rules and gotchas

Exact wherever possible; eigenstates/Q's at `prec`. Time one block/one LinearSolve before loops; nothing over ~3 min
without a timed smaller run; `TimeConstrained` 120 s on symbolic solves, fall back to rank tables at rational points
(say so). Markdown result cells under code cells; save after each task. Do not edit other notebooks, `paper/*.tex`,
HANDOFF.md, git. If the det₆ regression fails after two bookkeeping fixes, stop at Task 3, report the block tables
and the failure precisely (the main session decides).
- Shared kernel K1: the first cell wipes other notebooks' state (Stage A's notebook is saved; do not reuse its kernel
  state). Never write `*)` inside a comment. No `D`, `E`, `N`, `I`, `C`, `K` as symbols (use `Kt` for the tensor).
  Memoized definitions: purge with literal symbols. Semicolons between statements. Fresh names in `Table` iterators.
- `LinearSolve` with 225 right-hand sides at 60 digits is one dense solve — fine. Store K as rationals once.
- `MatrixRank` on rational matrices is exact; on symbolic matrices use `MatrixRank` at random rational points for the
  generic rank, and minors for the drop loci.
- The vacuum det₆ identity holds up to a state-independent constant; the block comparison is up to overall scale.
