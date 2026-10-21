# Spec: the 2603.26617 two-rotation family for gl(4) rectangular reps [S,S,0,0], L = 2 — untwisting, and why no contamination (2026-10-10)

Delegated to a wolfbook-builder subagent (sonnet); planned and reviewed by the main session. **New notebook:**
`Paul/Mathematica/new-conventions/L2-gl4-rect-untwist.wb` (create it; self-contained; first code cell clears the
shared kernel as every new-conventions notebook does). Reference cell sources (gl(3) patterns to adapt) are in the
scratchpad file named in the task message. gl(4) analogue of `L2-210-untwist-CG-spec.md` / `L2-210-untwist-CG.wb`.

## 0. The question (Paul)

In arXiv:2603.26617 (Bargheer–Bercini–Lefundes–Ryan, App. G) the gl(4) chain in the 6-dim rep [1,1,0,0] has a
product covector ⟨W_{ω,κ} ⊗ W_ω| of two *distinct* rotated vacua whose overlap with every eigenstate is a
product of Q-determinants (their (G12)–(G18)); expanding in the free parameters after untwisting
(z, κ, ω → 1 together, their eq. (4)) generates every CG coefficient at **leading order** in the twist.
In our gl(3) work (adjoint [2,1,0] and [3,1,0], `Untwist-two-rotation-notes.md`) the analogous det₆ × det₂ family
had to be expanded to **two or three orders** in ε because several CG slices of the undone covector sit at O(ε)
or O(ε²), at the same order as the first contamination K = ⟨covector|(1−P_HW)Ψ⟩. Paul wants the gl(4)
rectangular case set up in our conventions and wants to know **why it works there and not here**.

**Main-session hypothesis (to be confirmed or refuted by the notebook):** for [S,S,0,0] the HW covector e₁ᵀ is
stabilised by the parabolic with Levi gl(2)⊕gl(2), whose nilradical span{E₁₃, E₁₄, E₂₃, E₂₄} is **abelian**.
So every rotated vacuum is e₁ᵀ exp(X), X = Σ_{i≤2<j} x_ij E_ij, and the undone covector of the family is
e₁ᵀ exp(X′) with X′_ij = x_ij (δ_j/δ_i − 1) — entrywise, no commutator cross-terms. In the joint untwisting
regime every x_ij and every (δ_j/δ_i − 1) is O(1), so **all** weight components of the undone covector (incl. the
singlet slice ∝ det X′) are O(1), every CG slice is at ε⁰, and the contamination ratio K/S is O(ε) for every
HW state. In gl(3) the nilradical is the Heisenberg algebra and the E₁₃ coefficient of the undone
n′ = Δ⁻¹nΔn⁻¹ is the commutator combination b′ − a′c′ = a′c′ (z₂−w₂)/(y−z₂), which vanishes with the twist
difference z₂−w₂ = O(ε) because the ladder-insertion point y of the family is pinned at y₀ = (t−1)/(t−2) ≠ 1
(it does not untwist). That is the degeneration that pushed slices to O(ε), O(ε²).

The notebook must establish the gl(4) side with numbers; the gl(3) contrast is one symbolic cell (Task 5).

## 1. Conventions (new conventions, as in every new-conventions notebook)

- Lax ℒ⁽ˢ⁾_ij(u) = (u−θ_s)δ_ij + h E_ij, monodromy T(u) = ℒ⁽¹⁾(u−θ₁)ℒ⁽²⁾(u−θ₂), twist on the right T^M = T·M,
  t_a(u) = Σ_{|I|=a} qmin_{I,I}(u) (quantum minor: rows antisymmetrised, column k at u−(k−1)h), a = 1..4.
  t_a(u) has degree 2a in u at L = 2 (so 9 interpolation points; keep the name `uS` only for those points).
  For a diagonal twist the (I,I) minor is (Π_{i∈I} z_i) × the untwisted minor, so compute the untwisted basis
  minors once (`tBasisC[a, I]` pattern) and any twist is a linear combination.
- Representation: GT basis with the generic `ps`/`J0`/`Jp`/`Jm`/`JJ` code (it is rank-generic; it is the same
  code as in `old-conventions/old/Gaudin/Clean/Hamiltonians_SU4.wb` for gl(4)). Λ = {SS, SS, 0, 0} with
  `SS = 1` first (dim 6; dim 20 for SS = 2). `Ee[i,j] := JJ[i,j,Λ]`, i,j = 1..4; HW ket e₁ first (check:
  `Ee[i,j].e1 = 0` for all i<j, weight of e₁ = Λ). Covectors are rows acted on from the right; e₁ᵀEe[i,j] ≠ 0
  only for i ≤ 2 < j (check this: it is the rectangular/abelian property).
- Invariant inner product: `JJ` is not Hermitian; find the diagonal metric G with Ee[i,j]ᵀ G = G Ee[j,i]
  (pattern of `L2-210-untwist-CG.wb` cell 3, `gMet`); P_HW is the G⊗G-orthogonal projector onto the HW kets.
- Numeric point P1: h = 1/2, θ = (2/7, −3/5) (t := θ₁₂/h = 62/35). Working precision pr = 200 for numerics;
  everything that can be exact stays exact (rationals in ε).
- Diagonal frame: ket twist diag(z). **Companion matrix** of a twist tw = (w₁..w₄): G_c(tw) = {{c₁, −c₂, c₃, −c₄},
  {1,0,0,0}, {0,1,0,0}, {0,0,1,0}} with c_k the elementary symmetric polynomials (same sign pattern as the gl(3)
  `Gc`). **P_aux(z)**: the 4×4 matrix whose rows are left eigenvectors of G_c(z) ordered z₁..z₄
  (so P_aux G_c(z) P_aux⁻¹ = diag(z) — verify; the main session verified this is the gl(3) convention with
  sympy: Paux.Gc[z].Paux⁻¹ = diag(z)). Normalise each row so its first entry is 1. The **w-twist in the
  z-diagonal frame** is G′(z,w) = P_aux(z) G_c(w) P_aux(z)⁻¹.
- **One-site w-vacuum Ω_zw**: the left eigenvector of the one-site transfer matrix t₁^{(1)}(u; G′) =
  Σ_{ij} ℒ_ij(u−θ₁) G′_ji with the vacuum eigenvalue (u−θ₁)Σ_k w_k + h·SS·(w₁+w₂) — compute it as the
  left null vector of t₁^{(1)}(u₀; G′) − eigenvalue at one rational u₀ (exact if z, w rational), check it is
  u₀-independent and also a left eigenvector of t₂^{(1)}(u; G′) (and t₃, t₄ if cheap), normalise (Ω)₁ = 1.
  Then **verify the abelian form**: Ω_zw = e₁ᵀ exp(x₁₃E₁₃ + x₁₄E₁₄ + x₂₃E₂₃ + x₂₄E₂₄) exactly, with x_ij read
  off from the four components of weight Λ − e_i + e_j (the matrix exponential is a finite sum since the E's are
  nilpotent; note the E_ij commute, so the ordering of the four exponentials is irrelevant — check that too).
  Get **closed forms** x_ij(z, w) (exact symbolic `Together`/`Factor` with symbolic z, w; if the symbolic
  4×4 eigenvector solve is slow, rationalise from exact rational points instead and verify at fresh points).
  Record which of w₁..w₄ enter, which are torus directions (rescalings x_ij → x_ij d_j/d_i) and what the torus
  invariant is: the cross-ratio x₁₃x₂₄/(x₁₄x₂₃) as a function of z, w. Check the gl(3)-like expectation that the
  vacuum is representation-independent (same x_ij for SS = 1 and SS = 2 — at least at one numeric point;
  for the first pass with SS = 1 just record x_ij).
- ρ(D) for a diagonal D = diag(d₁..d₄): the torus element `ρdiag[d]` = diag over GT patterns of Π d_k^{J0[k][p]}.
  δ^Λ := (δ₁δ₂)^SS so that e₁ᵀρ(diag δ) = δ^Λ e₁ᵀ.

## 2. The family (paper (G13)–(G14) in our conventions)

  v₂ = Ω_zw,   v₁ = Ω_zw · ρ(Δ),   Δ = diag(w₂′ − z_k)  (paper: diag(1 − κ z_k) with κ = 1/w₂′; same torus up to a scalar).

NB this is the **polynomial** torus diag(w₂′ − z_k), not the resolvent diag(1/(w₂′−z_k)) of the gl(3) family —
the paper's (G14) C_ij = (1−κz_j)/(1−κz_i) R_ij is reproduced by exactly this choice (check: v₁ = δ^Λ e₁ᵀ exp(X_C)
with (X_C)_ij = x_ij δ_j/δ_i, δ_k = w₂′ − z_k). w₂′ is a free parameter, independent of w (Paul's "two rotated
vacuum states"). Keep w₂′ as its own symbol (call it `wk` in code; never `κ`-related names that collide).

**Undo (exact for HW kets).** n := exp(X) ∈ the abelian unipotent group; ρ(n)⊗ρ(n) fixes every HW ket of V⊗V
(it is exp of a sum of total raising operators). Hence ⟨v₁⊗v₂|h_ν⟩ = ⟨φ_u ⊗ e₁ᵀ|h_ν⟩ with
φ_u := v₁ ρ(n)⁻¹ = δ^Λ e₁ᵀ exp(X′), X′_ij = x_ij (δ_j/δ_i − 1). Verify the exactness of the undo on all HW kets
and the formula for X′ (symbolically in δ, x if cheap; else exact rationals).

## 3. The untwisting regime

  z = 1 + εφ,  w = 1 + εω (ω = (ω₁..ω₄)),  w₂′ = 1 + εk,

fixed rational directions, ε → 0 (rational paths). Directions: φ = (1/3, −2/7, 5/11, −4/9) (primary),
φ′ = (−2/5, 3/8, 7/11, −1/6) (control); ω ∈ a small grid, e.g. ω₁ ∈ {3/4, −5/3}, ω₂ ∈ {−1/5, 9/7, 5/3},
ω₃ = 2/9, ω₄ = −7/4 fixed (vary them once to see whether they matter); k ∈ {−12/5, 7/2, 1/8}. Never let k or any
ω_i coincide with a φ_j (poles), and avoid k = ω_i. ε ∈ {10⁻⁶, 2·10⁻⁶, 4·10⁻⁶} for Richardson. Expectation (to
confirm): every x_ij(ε) is O(1) with an O(1) limit (ratios of O(ε) differences), every δ_j/δ_i − 1 =
(z_i − z_j)/(w₂′ − z_j) → (φ_i − φ_j)/(k − φ_j) is O(1), hence X′(0) is a generic 2×2 block and
φ_u(0)/δ^Λ = e₁ᵀ exp(X′(0)) has O(1) components at every weight of the orbit, with full rank over the grid.

## 4. Tasks (each its own cell(s) with a `Grid` of checks; keep `checks` as an Association; timing rule: time one
state / one parameter point before any loop; nothing over ~3 minutes without a timed smaller run)

**Task 0 — setup checks.** gl(4) relations of `JJ` on Λ; G metric; HW kets of V⊗V (common null space of the six
total raising operators E_ij⁽¹⁾+E_ij⁽²⁾, i<j, per total-weight sector); the HW sectors and multiplicities:
SS = 1: [2,2,0,0], [2,1,1,0], [1,1,1,1] (3 irreps, all multiplicity 1); SS = 2: [4,4,0,0], [4,3,1,0], [4,2,2,0],
[3,3,1,1], [3,2,2,1], [2,2,2,2] (6 irreps, multiplicity 1 each) — report what you find. Quantum determinant check
t₄ = (Πz) ν₁(u)ν₂(u−h)ν₃(u−2h)ν₄(u−3h), ν_n(x) = Π_s (x−θ_s+hΛ_n), at two u's (SS = 1; for SS = 2 only if cheap);
[t₁(u), t₂(v)] = 0. One-site vacuum: eigen-bra checks, abelian form, closed forms x_ij, torus invariant
(§1). Also record: the weights of V reached by e₁ᵀexp(X) for generic X (should be all weights of V, with the
multiplicity-2 weight (1,1,1,1) of [2,2,0,0] reached in both directions — relevant for SS = 2).

**Task 1 — localization and the product/Q-number structure (finite twist, exact rationals, z = (3/2,5/4,2/3,4/7),
w = (5/2,7/3,3/5,9/8), one or two values of w₂′).** (a) Rank scan: for the shapes T_λ built by CBR fusion from
t₁..t₄ (λ ∈ {(1,0),(1,1),(2,0),(2,1),(2,2)}; T_{(1,0)} = t₁, T_{(1,1)} = t₂, T_{(2,0)}(v) = t₁(v)t₁(v+h) − t₂(v+h),
T_{(2,1)}(v) = t₂(v)t₁(v+h) − t₃(v+h), T_{(2,2)}(v) = t₂(v)t₂(v+h) − t₃(v+h)t₁(v) + t₄(v+h)·(check the CBR
formula you use against the gl(3) convention T_λ(x) = det_{i,j} t_{λ′_j+i−j}(x+(i−1)h))) and the points v =
θ₁ + m h, m ∈ {−SS−1, …, 2}: rank of ⟨Ω_zw⊗Ω_zw|T_λ(v) as a d₁×d₁ matrix (rank 1 ⇔ localised on site 1 as
(site-1 covector)⊗Ω_zw). Expectation from the paper's App. M (their (M4)–(M8), translated): t₁ and t₂ localise at
**v = θ₁** (the point where ℒ₃₃ = ℒ₄₄ = 0 on the HW covector; the paper's t_{[1,0]}(θ₁), t_{[1,1]}(θ₁)); report
all rank-1 (λ, m). (b) The family is in the localising span: for three values of w₂′, solve exactly
⟨v₁⊗v₂| = ⟨Ω⊗Ω| (c₀ + Σ_λ c_λ T_λ(v_λ)) over the localising shapes found in (a); report the c's and check
they are polynomial in w₂′ (for SS = 1: the paper's 1 − κ t_{[1,0]} + κ² t_{[1,1]} up to our normalisations,
i.e. c₀ : c₁ : c₂ ∝ w₂′² : −w₂′·n₁ : n₂ or similar with w₂′-independent n's — find the exact form). Consequence
(state it in a markdown cell): ⟨v₁⊗v₂|Ψ_B⟩ = (c₀ + Σ c_λ τ_λ^B)·⟨Ω⊗Ω|Ψ_B⟩, a Q-number, for every w₂′ — the
Q-determinant form (paper (G7)) is NOT needed for this notebook. (c) Verify on the actual eigenstates at this
finite twist (all states of the HW sectors): ⟨v₁⊗v₂|Ψ_B⟩/⟨Ω⊗Ω|Ψ_B⟩ equals the eigenvalue combination
(relative deviation, with Precision).

**Task 2 — exact ε-structure of the undone covector (symbolic ε, rational directions).** As in
`L2-210-untwist-CG.wb` cells 11–12: build φ_u(ε)/δ^Λ in the regime of §3 with `epsS` symbolic, `Series` each
component to O(ε²); per weight of V: dimension, leading order, leading coefficient; verify (i) each x_ij(ε) is
ε-free or has a finite nonzero limit (give the limits in terms of φ, ω — factor them), (ii) δ_j/δ_i − 1 →
(φ_i−φ_j)/(k−φ_j), (iii) φ_u(0)/δ^Λ = e₁ᵀ exp(X′(0)) exactly, (iv) det X′(0) ≠ 0 and the singlet slice
component (weight (0,0,SS,SS): (0,0,1,1) for SS = 1, (0,0,2,2) for SS = 2, whose coefficient is ∝ (det X′)^SS)
is O(1). Then over the (ω, k) grid: per-weight ranks of the ε⁰ coefficient vectors (exact);
expectation: rank = dim for every weight at ε⁰ already (for SS = 2 including the dim-2 weight (1,1,1,1)).
Repeat for φ′. Slice weights ν − Λ: SS = 1: (1,1,0,0), (1,0,1,0), (0,0,1,1); SS = 2: (2,2,0,0), (2,1,1,0),
(2,0,2,0), (1,1,1,1), (1,0,2,1), (0,0,2,2).

**Task 3 — chain side: signal vs contamination in the regime (the decisive test).** For ε ∈ {1,2,4}·10⁻⁶, φ
primary, the (ω,k) grid: eigenstates of the HW sectors at z = 1+εφ (per-sector `Eigensystem` of a generic
commuting member t₁(u₀) + π t₂(u₁) restricted to the sector, u₀ = 3/7, u₁ = −2/9, at pr digits; HW-limit
branches = the dim(HW_μ) states with the smallest raising residual; normalise each Ψ_B by its largest
component; store ‖Ψ‖, ‖P_HWΨ‖/‖Ψ‖, raising residual, eigenvalue). For each branch and grid point:
R_C := ⟨v₁⊗v₂|Ψ_B⟩, S := ⟨v₁⊗v₂|P_HW Ψ_B⟩, K := R_C − S, κ := R_C/S; check S = ⟨φ_u⊗e₁ᵀ|P_HWΨ_B⟩ (exact
undo); ε-exponents of |S|, |K| (ratio between consecutive ε's, log₂), and max over the grid of |K/S| at each ε.
**Expectation to confirm or refute:** K/S = O(ε) (κ → 1 linearly) for EVERY branch, including the singlet
[1,1,1,1] (SS = 1) and all six irreps (SS = 2). If any branch has κ ↛ 1, do the `showContam` breakdown
(which site-1 weights carry K, at which ε-powers) for it and report; do not fit anything.

**Task 4 — leading-order CG extraction from the overlap alone.** For each non-top branch with a 1-dim slice
weight: C(ω,k,ε) := R_C/(δ^Λ · [φ_u(0)/δ^Λ]_{slice component}) must be (ω,k)-independent up to O(ε) and equal
to the chain-side slice ⟨m_slice⊗e₁|P_HWΨ_B⟩·(normalisation) — report the spread over the grid and the relative
error, at two ε's (both ∝ ε), and the Richardson limit vs the chain-side limit. For the dim-2 slice (1,1,1,1)
at SS = 2: the two ε⁰ functionals (the x′₁₃x′₂₄ and x′₁₄x′₂₃ monomials of e₁ᵀexp(X′) at that weight) — show
the 2×2 system over two grid points determines both slice components, with O(ε) error. Conclusion cell: at
leading order in ε, with no higher-order expansion, every CG datum of [S,S,0,0]⊗[S,S,0,0] is read off from
Q-numbers (the overlap R_C = Q-number by Task 1).

**Task 5 — the gl(3) contrast, symbolic (one cell, no chain).** With symbols: (a) gl(3), n = {{1,a,b},{0,1,c},
{0,0,1}}, Δ = diag(δ₁,δ₂,δ₃), n′ = Δ⁻¹ n Δ n⁻¹: show n′₁₃ − n′₁₂n′₂₃ = (δ₃/δ₁ − 1)(b − ac z₂₁/z₃₁) when
δ_k = 1/(w₂−z_k), and with the gl(3) vacuum (b = r₂ + r₁r₃, r's of `oneSiteVac` with (w₁, y)) that this equals
a′c′ (z₂−w₂)/(y−z₂) (a′ = n′₁₂, c′ = n′₂₃), i.e. O(ε) in the gl(3) regime because y → y₀ ≠ 1; (b) gl(4) abelian:
n = exp(X) with a symbolic 2×2 block X, n′ = Δ⁻¹nΔn⁻¹ = exp(X′) with X′_ij = x_ij(δ_j/δ_i − 1) exactly (4×4
matrices), so no entry can vanish unless a twist difference does. Put the comparison in a markdown cell: the
paper's limit (4) untwists κ (the torus point) together with z and ω; in the gl(3) family the ladder insertion
point y = w₂ y₀ cannot untwist (y₀ = (t−1)/(t−2) ≠ 1 is forced by the product condition), and the commutator
direction E₁₃ carries the factor (z₂ − w₂).

**Task 6 — SS = 2.** After everything passes for SS = 1 and the notebook is saved, duplicate the file as
`L2-gl4-rect-untwist-S2.wb`, set `SS = 2`, rerun Tasks 0–4 (Task 1 rank scan may be slow: time one shape/point
first; reduce the grid to 2×2 if the Task 3 loop exceeds ~10 minutes; report timings). Six HW branches; the
dim-2 slice (1,1,1,1) of [3,3,1,1] is the interesting one.

## 5. Reporting

End each notebook with a markdown cell `### gl(4) [S,S,0,0] untwist result (2026-10-10)` listing, per task, the
key numbers (residuals with their Precision, the per-weight order/rank tables, the κ tables with ε-exponents,
the CG-extraction errors, the localisation table, the closed forms x_ij and the torus invariant) and every
structural surprise. Save (`wolfbook_saveNotebook`). If VS Code does not execute cells of the new file (this
happened for every new notebook on 2026-10-10), run each cell's source in the shared kernel K1 via
`wolfbook_evaluateExpression` and store the printed results in a markdown cell directly under the code cell, as
the earlier agents did; say so in the report. Then report back concisely in the same structure, numbers not
adjectives. Do not edit any other notebook; do not touch `paper/*.tex`; do not commit. If a task's expectation is
refuted or a check fails and you cannot fix it within two attempts, stop at that point, keep what works, and
report what you saw.

## 6. Gotchas (all have bitten before)

- The shared kernel K1 is the only kernel on this machine; your first cell's `Quiet[ClearAll["Global`*"]]` wipes
  the live state of the other new-conventions notebooks (files are intact). Fine, but do not evaluate cells of
  other notebooks afterwards.
- Never write the two characters `*)` inside a comment (e.g. "y*)" closes it).
- Do not redefine `uS` or use `D`, `E`, `N`, `I`, `C`, `K` as symbols (K is Protected? — use `KK`/`kk` names).
- Memoized functions cache bad results: after fixing a definition, re-run the defining cell and everything
  downstream. `Do[DownValues[f] = …, {f, {…}}]` purges nothing (HoldAll) — use literal symbols.
- Semicolons between statements in multi-line cells; keep each expression on one line where possible.
- `Max[{}]` → −Infinity; floor at 0.
- Precision: compute at pr = 200, report `Precision` of residuals; `SetPrecision` before `Norm` of near-zero
  vectors; exact zero tests on rationals with `===`, on numerics never with `PossibleZeroQ`.
- Wolfbook replies truncate at `{`; read long results via single-line
  `StringReplace[ToString[.., InputForm], "\n" -> " "]`.
- `Ω_zw` degenerates when a w_i equals a z_j or when z's coincide; the (ω,k) directions above avoid this at
  small ε, but if you change a direction, check.
- Use SparseArray for the 400-dim objects at SS = 2 (E₁, E₂, quantum minors); exact rational dense 400×400
  products are slow.
- Sector `Eigensystem` at pr digits: feed `N[matrix, pr]`; for the exact HW kets use exact `NullSpace`.

## Addendum A (main session, after a sympy cross-check of §1–§3; supersedes §2–§3 where they conflict)

**Fact (sympy, 6-rep, exact at five random rational (z,w) points and at P1):** the companion-frame vacuum
Ω_zw = e₁ᵀ exp(X) has x₁₃x₂₄/(x₁₄x₂₃) = (z₁−z₄)(z₂−z₃)/((z₁−z₃)(z₂−z₄)), the cross-ratio of the **z's**,
for every w; w₃, w₄ drop out and w₁, w₂ only rescale the x_ij (torus directions). Consequently the undone
block of the κ-torus family of §2, X′_ij = x_ij(δ_j/δ_i − 1) = x_ij (z_i − z_j)/(w₂′ − z_i), has
**det X′ ≡ 0 identically** (finite twist, any w₂′): e₁ᵀexp(X′) has no (0,0,1,1) component, so the family
v₁ = Ω_zw ρ(diag(w₂′−z)), v₂ = Ω_zw never sees the singlet [1,1,1,1] (and at SS = 2 misses every slice that
needs det X′ ≠ 0: (1,0,2,1), (0,0,2,2), and one direction of (1,1,1,1)). **This is not the paper's family.**

**The paper's ⟨W_ω| (their (F4)–(F5)):** ⟨W_ω| = ⟨Z|exp(R) with R_ij = V_ij (1 − ω z_a z_b)/(1 − ω z₁z₂), where
(a,b) is the weight of the component e₁ᵀE_ij (e₁₃ ↔ (1,3): R₂₃/V₂₃ = (1−ωz₁z₃)/(1−ωz₁z₂), e₁₄: R₂₄/V₂₄ =
(1−ωz₁z₄)/(…), e₂₃: R₁₃/V₁₃ = (1−ωz₂z₃)/(…), e₂₄: R₁₄/V₁₄ = (1−ωz₂z₄)/(…)), and V = our companion-frame
vacuum block x (its z-cross-ratio is exactly what makes this consistent at the (0,0,1,1) component, where
det R·(1−ωz₁z₂) must equal det V·(1−ωz₃z₄); the main session checked this algebraically). In our diagonal
frame this is simply

    ⟨W_ω| = Ω_zw · (1 − ω ρdiag[z]) / (1 − ω z₁z₂),

ρdiag[z] the 6×6 torus matrix (weights z_az_b) — a linear combination of two group elements, not a group
element, which nevertheless maps the special rotation Ω_zw to another rotation (check this exactly: W_ω =
e₁ᵀexp(R) with the (F5) ratios, and the (0,0,1,1) component equals det R). Its cross-ratio is ω-dependent:
CR(R) = CR(z) · (1−ωz₂z₃)(1−ωz₁z₄)/((1−ωz₂z₄)(1−ωz₁z₃)).

**The paper's family, corrected §2:**

    v₂ = W_ω,   v₁ = W_ω · ρ(diag(w₂′ − z_k))   (κ = 1/w₂′; paper (G13)–(G14) with R the ω-deformed block).

Undone block: X′_ij = R_ij(δ_j/δ_i − 1), det X′ ∝ CR(R) − CR(z) ∝ ω (z₂−z₁)(z₃−z₄)/((1−ωz₂z₄)(1−ωz₁z₃)).

**Corrected regime §3:** z = 1+εφ, w = 1+εω_w (ω_w a fixed torus direction, e.g. (3/4, −1/5, 2/9, −7/4); vary once),
**ω = 1 + ε o** with o ∈ {−12/5, 7/2, 1/8, 5/3} (never o = −(φ_a+φ_b) for any pair, poles of (1−ωz_az_b)),
w₂′ = 1 + εk with k ∈ {−1/5, 9/7, 5/3} (never k = φ_j). Expectation: CR(R)/CR(z) → (o+φ₂+φ₃)(o+φ₁+φ₄)/
((o+φ₂+φ₄)(o+φ₁+φ₃)) = O(1) ≠ 1, so det X′(0) ≠ 0 and the singlet slice is O(1) — but ONLY because ω untwists
together with z (if ω stayed O(1) away from 1, CR(R)/CR(z) − 1 = O(ε²) and the singlet slice would sit two
orders down: the exact analogue of the gl(3) y₀ ≠ 1 problem). Record both: the joint regime (main runs) and
one control run with ω fixed at, say, 3/7 (expect the singlet slice at O(ε²) and K/S ↛ 0 for the singlet).

**Task 1 (localisation) stays as written but with ⟨W_ω ⊗ W_ω| in place of ⟨Ω⊗Ω| for part (b)/(c)**; also test
whether ⟨W_ω⊗W_ω| itself lies in the localising span of ⟨Ω_zw⊗Ω_zw| (exact solve over the rank-1 shapes at θ₁,
plus the global torus character) — for SS = 1 the paper's (F3) says ⟨W_ω|Ψ⟩ is a Q-determinant, so the
overlap of the full family is a Q-number for all (ω, κ); report what the exact solve gives.

**SS = 2:** define W_ω by the same rep-independent group element e₁ᵀexp(R) with the (F5) ratios applied to the
SS = 2 vacuum block (check first that the SS = 2 companion vacuum has the same x_ij as SS = 1 at a numeric
point), and test (i) whether ⟨W_ω⊗W_ω| and the κ-family are in the localising span of ⟨Ω_zw⊗Ω_zw| (if not,
say so: then the Q-number property is open for SS = 2 and the untwist analysis is chain-side only), (ii) the
untwist structure (Tasks 2–4) — now all six slices, including the dim-2 (1,1,1,1), are expected O(1).

**Task 5 stays**, with this sharpened comparison for the markdown cell: gl(4)-rect has TWO independent
localising deformations at the single point θ₁ — the 4-torus Λ²(1−κG) ↔ 1 − κt_{[1,0]} + κ²t_{[1,1]} (a group
element) and the Λ²-polynomial 1 − ω ρ(G) (not a group element, but rotation-preserving because of the special
cross-ratio of the companion vacuum) — and the product condition does not tie κ to ω, so both untwist freely and
the undone block at ε = 0 is a generic 2×2 matrix (abelian nilradical: entrywise, no commutator direction). In
gl(3)-adjoint the localising span at θ₁−2h meets the product locus only at the isolated point c*, which ties the
ladder insertion y to the resolvent point w₂ by the θ-dependent factor y₀ = (t−1)/(t−2); the commutator
direction E₁₃ (which exists only because N₊ is non-abelian) carries b′ − a′c′ = a′c′(z₂−w₂)/(y−z₂), O(ε) because
y cannot untwist. Also note the finite-twist blind spot of the pure κ-torus family (det X′ ≡ 0) as the gl(4)
shadow of the gl(3) "b′ = 0 ⟺ ω = y" blind spot.
