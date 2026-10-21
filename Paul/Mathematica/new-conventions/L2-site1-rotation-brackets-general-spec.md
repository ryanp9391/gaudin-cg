# Spec: the site-1 rotation on the SoV vacuum as ladder-1 Q-brackets — generalisation to [3,1,0] and [3,2,0] (2026-10-10)

Delegated to a wolfbook-builder subagent (sonnet); planned and reviewed by the main session. **New notebook:**
`Paul/Mathematica/new-conventions/L2-site1-rotation-brackets-general.wb` (self-contained, rep-generic code with
`Λ` as the only rep switch; first code cell `Quiet[ClearAll["Global`*"]]`). Reference cell sources (scratchpad
files named in the task message):
- `ref-L2-210-site1-rotation-brackets-all.txt` — the adjoint notebook of this construction (all 36 cells): copy
  its structure; the generic parts (SoV basis of B, K map, tensors, Plücker tools, loci solve, reach table) are
  what you generalise. Its write-up is `L2-210-site1-rotation-brackets.md` — read it first.
- `ref-L2-310-two-rotations-cells1-12.txt`, `ref-L2-320-two-rotations-cells1-8.txt` — exact sparse transfer
  matrices for dim-15 reps, all 225 eigenstates at 200 digits, dual-Q solver with degrees μ_j, Ψ_Q with the
  rep's vacuum triple, companion/diagonal frames (rep-independent Gauss data).
- `ref-L2-rotation-sov-cells1-22-59-80.txt` — adjoint originals (B operator `Bop2`, `TMf`, `qminor`, SoV labels).

## 0. What was found for the adjoint and what to generalise

Adjoint result (`L2-210-site1-rotation-brackets.md`): with x₀ = e₁ᵀ⊗e₁ᵀ (companion frame) and a site-1
rotation g = exp(αE₁₂)exp(βE₁₃)exp(γE₂₃), the overlap ⟨e₁ᵀρ(g)⊗e₁ᵀ|Ψ⟩/⟨x₀|Ψ⟩ · Ψ₁ = Σ_x φ_x M_x over the eight
site-1 SoV labels x = (a|bc) (triple of roots θ₁ + {−1−a, −b, −c}h) with monomials M_x =
det₂[v(θ₁−(1+a)h), v(θ₁−bh)]_{rows}·v_k(θ₁−ch), i.e. a coefficient tensor C_{(nn′),c} (2-form in the minor index
(n,n′) per column c). On the quadric 2(t−2)(β+αγ) = (t−3)α² every column is a decomposable 2-form, so
N = Σ_c det₂[⟨Q^jF_c⟩,⟨Q^jG_c⟩]·Q³(θ₁−ch) — "sum over c of one bracket determinant times a point value" (R2δ);
the λ₂¹ group (site 2 at the label (λ₁|λ₁ 0)) has this form for every rotation. Generic rotations need one
more bracket determinant (column 0 has rank 4 on four points).

**Label rule to verify for general Λ = (λ₁, λ₂, 0):** SoV labels (a|bc) with λ₂ ≤ a ≤ λ₁, 0 ≤ c ≤ λ₂, c ≤ b ≤ a
(count = dim V: 15 for both [3,1,0] and [3,2,0]); roots θ_s + {−1−a, −b, −c}h; vacuum (λ₁|λ₁λ₂) with triple
(−1−λ₁, −λ₁, −λ₂) (matches the known vacuum triples (4,3,1) for [3,1,0], (4,3,2) for [3,2,0] and the minor
W(λ₁, λ₁+1)); wavefunction monomial M_{(a|bc)} = det₂[v(θ−(1+a)h), v(θ−bh)]_{rows}·v_k(θ−ch), normalised by
the vacuum. If the monomial rule fails for some labels (larger reps may need shifted points), find the correct
one from the B-eigenvalue roots and the 225-state data, and record it.

**Expected tensor structure:** columns c ∈ {0..λ₂} (λ₂+1 columns), each a 2-form on the ladder points
{−1−λ₁, …, 0} (λ₁+2 points); column c is supported on pairs (−1−a, −b) with b ≥ c, so column c lives on the
points {−1−λ₁, …, −c} ∪ … — report the support per column. A 2-form on p points has rank ≤ p (even), so a
column needs ≤ ⌈p/2⌉ bracket determinants for ANY covector: the automatic bracket count is Σ_c ⌈p_c/2⌉ per
site-2 group. That count is polynomial in (λ₁, λ₂) — Paul's "scales with boxes" criterion is met automatically;
the loci are where fewer suffice. [3,1,0]: columns c = 1 (4 points, rank ≤ 4; decomposable ⇔ 1 Plücker
quadric) and c = 0 (5 points; decomposable ⇔ 3 independent conditions, rank ≤ 4 automatic). [3,2,0]: c = 2
(points {−4,−3,−2}: always decomposable), c = 1 (4 points), c = 0 (5 points). So the full (R2δ) analogue is
expected to need 4 conditions on the 3-parameter rotation (plus scale) — likely empty or isolated for [3,1,0];
the partial loci (which columns are decomposable) are the result. Compute, don't assume.

## 1. Conventions

Companion frame, L = 2, new conventions, h = 1/2, θ = (2/7, −3/5) (t = 62/35), ket twist z = (3/2, 5/4, 2/3)
(+ z′ = (3/2, 5/4, 3/7) once for z-independence). `prec` = 60 for eigenvectors; dual Q's as in the 310 notebook
(200 digits are not needed; 60 suffice for 10⁻⁴⁰ checks). All loci exact. Rotation parametrisation: Gauss chart
(α,β,γ) as in the adjoint notebook cell 21 (`gaussAB`), and the (a,b) cone form if you need the charts at
infinity. "Genuine" = αβγ ≠ 0. Check the orbit of e₁ᵀ is 3-dimensional for both reps (Jacobian rank 3).

## 2. Tasks (rep-generic code; run the adjoint FIRST as regression — it must reproduce the quadric with
43/16 = (t−3)/(2(t−2)) — then [3,1,0], then [3,2,0]; each its own cells with a `Grid` of checks)

**Task 0 — setup per rep.** Generators `JJ`, metric, exact sparse transfer matrices (310 pattern), the 225 (64)
eigenstates of the twisted chain at z (companion-frame kets `ketC` = P_phys⁻¹⊗P_phys⁻¹ · diagonal-frame kets),
x₀ = e₁ᵀ⊗e₁ᵀ, dual Q's with degrees μ_j and the Baxter residual, `colV2[B,s,m]` = v(θ_s+mh) for m ∈
{−1−λ₁, …, 0} (and m = +1 if needed). Timing per step.

**Task 1 — SoV basis and labels.** Untwisted B(u) = T₁₁(u)·qmin_{12,12}(u+h) + T₂₁(u)·qmin_{12,13}(u+h) (rep-generic
code of the adjoint cell 7/16); left eigenvectors (simple spectrum? check; if degenerate, use two u-values
jointly); roots of each eigenvalue → per-site triples; verify the label rule of §0 (count, vacuum triple, every
triple of the form {−1−a, −b, −c}); record the bijection label ↔ triple. Wavefunction: ⟨x|Ψ_B⟩/⟨x₀|Ψ_B⟩ =
N_x Π_s M^{(s)}_x / (vacuum) for all 225 states and all labels, for the row/k choices ({1,2},1), ({2,3},3),
({1,3},2) — report which hold and the max relative deviation. If the monomial rule fails, search the
alternatives (which two of the three roots form the minor; shifts by ±h of the single point) and report.

**Task 2 — site-1 covectors ⊗ e₁ᵀ.** For a generic rotation and for the dim-V basis covectors x_ν ⊗ e₁ᵀ (x_ν the
L=1 SoV covectors — rebuild them as the left eigenvectors of the one-site untwisted B, normalised as in the
adjoint `xL1`, or simply use the GT basis covectors e_νᵀ⊗e₁ᵀ): support in the L=2 SoV basis — which site-2
labels occur (expect the vacuum (λ₁|λ₁λ₂) and the labels (λ₁|λ₁ c), c < λ₂, with factors
λ₂^{(c)} = v_k(θ₂−ch)/v_k(θ₂−λ₂h) — verify, incl. k-independence), the K map per site-2 group (dim V × dim V,
rational: rationalise and check), ranks. The overlap then is Σ_groups λ₂^{(c₂)} Σ_x φ^{(c₂)}_x ψ_x.

**Task 3 — tensors and automatic bracket count.** For each site-2 group: C^{(c₂)}_{(nn′),c} (rows = minors in a
fixed order over the λ₁+2 points, columns c = 0..λ₂), populated slots, support per column; for a generic
rotation (rational point) and a random covector: the rank of each column's 2-form (= 2 × number of bracket
determinants needed) — give the automatic count per group and in total, and compare with dim V. Verify the
bracket identity Bk₂[F,G] = Σ(F∧G)·minor on the larger point set, and the Pfaffian/Plücker tools in your minor
order (all 4×4 Pfaffians of a 2-form on p points: decomposable ⇔ all vanish; rank ≤ 4 ⇔ all 6×6 Pfaffians vanish).

**Task 4 — loci (exact, Gauss chart; check charts at infinity only if a component seems to escape).** With
C(α,β,γ) bilinear-in-(a,b) hence polynomial in (α,β,γ): for each site-2 group and each column c, the ideal
"column c decomposable" (4×4 Pfaffians): `Factor` the generators, `Reduce`/Gröbner for the locus, dimension and
components, genuine rotations on it or not; the joint locus "all columns of the λ₂^{(λ₂)} (vacuum) group
decomposable" ((R2δ) analogue) and the same for every group; the weaker "rank ≤ 4" loci where the column has
> 4 points. For the adjoint this must give back {β = 0} ∪ {2(t−2)(β+αγ) = (t−3)α²}. For [3,1,0]/[3,2,0] report
per column: the polynomial(s), their factorisation, the locus dimension in the 3-dim Gauss chart, and whether a
genuine rotation lies on it. Also test the natural guess that the adjoint quadric generalises as
"2(t−2)·b₁b₃ = (t−3)·b₂²" on the first row b of g, or with t → t − (λ₁−2)-type shifts: substitute and see if
it lies in some column's ideal (report, don't force).

**Task 5 — explicit bracket forms and verification.** (i) For a generic rotation: the automatic form (per column
a rank-decomposition of the 2-form into ⌈p/2⌉ wedges — any exact decomposition, e.g. via a symplectic-type
Gram–Schmidt over the rationals), explicit F, G vectors, 225-state check of the full overlap identity
(relative deviation; expect 10⁻⁴⁰ or better). (ii) At a rational point of the best locus found (fewest
brackets): the reduced form, F, G, verification. (iii) Count: number of bracket determinants in (i) and (ii),
vs boxes |Λ| and vs dim V.

**Task 6 — reach table (exact, finite twist, no untwisting).** As the adjoint Task 7: P_phys⁻¹ = L₀D₀U₀ (verify
e₁ᵀρ(P⁻¹) ∝ e₁ᵀρ(U₀) exactly for the rep), undone covector e₁ᵀρ(gL₀D₀); for a generic rotation and for three
rational points of each locus of Task 4 that contains genuine rotations: weights of V reached and the rank on
every multi-dimensional weight space (list the weights of V with their multiplicities first), compared with the
slice weights needed for Λ⊗Λ (the HW weights ν of Λ⊗Λ minus Λ; take them from `L2-310-untwist-CG.md` §2 for
[3,1,0]: 8 irreps/10 HW states; for [3,2,0] compute the HW sectors as in the untwist notebooks' cell 6 pattern).

**Task 7 — z-independence and t-dependence.** Loci at z′ identical? For the main locus polynomials of each rep:
coefficients at one more (h,θ) (h = 1/3, θ = (1/5,−4/9), t = 29/15) and a guess of the t-dependence if a single
coefficient is involved (as c(t) = (t−3)/(2(t−2)) for the adjoint).

**Task 8 — markdown summary** `### site-1 rotation brackets, general reps — result (2026-10-10)`: per rep the
label rule and wavefunction monomial (verified), site-2 groups and factors, the tensors' supports, automatic
bracket counts (generic rotation / random covector) vs dim V and |Λ|, the loci per column with dimensions and
genuine-rotation content, the best reduced form and its count, the reach table, z/t-dependence, timings, surprises.

## 3. Rules

Exact wherever possible; eigenstates/Q's numeric at `prec`. Time one case before loops (225-state loops: time
one state); nothing over ~3 minutes without a timed smaller run; the Gröbner solves: start with one column and
`TimeConstrained` at 120 s, fall back to `Reduce` on a chart or to numeric-point rank tables if a symbolic
solve stalls (say so). Markdown result cells under code cells; save after each task; if VS Code does not
execute cells of the new file, evaluate the cell sources in K1 and store the printed results as markdown
(say so). Do not edit other notebooks, `paper/*.tex`, git. If the label rule or the monomial rule fails for
a rep and you cannot repair it within two attempts, stop for that rep, keep what works, report.

## 4. Gotchas

- Shared kernel K1: the first cell wipes other notebooks' state. Never write `*)` inside a comment. No `D`, `E`,
  `N`, `I`, `C`, `K` as symbols. Memoized definitions: purge with literal symbols. Semicolons between statements.
- The adjoint notebook's `wedge2` had a symbol clash with the Table variable `l` — use fresh names.
- Degenerate B spectrum: for dim-15 reps B(u) at one u might have coincident eigenvalues; diagonalise a
  combination B(u₁) + πB(u₂) or use the roots jointly.
- Companion-frame kets: `ketC` = (P_phys⁻¹⊗P_phys⁻¹)·Ψ (check the direction against the adjoint cell 9/18).
- `Position`-with-level-spec and `FirstPosition` pitfalls (the earlier agent hit `Position[.., {1}]` matching the
  List head) — select by exact equality of the label.
