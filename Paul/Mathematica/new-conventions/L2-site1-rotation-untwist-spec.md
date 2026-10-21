# Spec: untwisting the site-1-rotation bracket families — adjoint quadric and [3,1,0] (2026-10-10)

Delegated to a wolfbook-builder subagent (sonnet); planned and reviewed by the main session. Two jobs, in order,
on the shared kernel:

**Job A (short) — finish `L2-site1-rotation-brackets-general.wb`.** Tasks 0–6 are complete and saved (38 cells;
cells 37–38 are placeholders; cells 31–36 are the Task 7 driver). Paul restricted the scope to [3,1,0]: run Task 7
ONLY for [2,1,0] (regression) and [3,1,0], and only at two extra parameter points (z′ = (3/2,5/4,3/7) at the base
(h,θ), and (h,θ) = (1/3, (1/5,−4/9)) at the base z), i.e. 3 points per rep instead of 8; skip the [3,2,0] cell 34
and the [3,2,0] addendum cell 36 (leave them in place, mark them "not run, scope reduced to [3,1,0]" in the
markdown). Report: are the loci polynomials z-independent; which conic constants cc change with t (give the
values at t = 29/15 next to those at t = 62/35) — the [3,1,0] best locus has cc = 1: is it t-independent?
Then write the summary cell (replace the placeholders) in the structure of the general spec's Task 8, [3,2,0]
included for what was computed. Save. **Then re-run cells 3–6 of that notebook with Λ = {2,1,0}?** No — the
untwist notebook below is self-contained; just leave the kernel as it is.

**Job B — new notebook `Paul/Mathematica/new-conventions/L2-site1-rotation-untwist.wb`** (self-contained, first
code cell `Quiet[ClearAll["Global`*"]]`, rep-generic with Λ as the switch; run [2,1,0] first, then [3,1,0]).
References (scratchpad, named in the task message): `ref-L2-210-untwist-CG-cells1-12.txt` (the adjoint untwist
pipeline: regime `reg`, HW kets and G-orthogonal projectors `PHWc`, `branchesAt` with the HW-limit branches,
exact-undo check, symbolic-ε structure of the undone covector, signal/contamination S, K, κ = R_C/S, the
leading-order CG extraction and Richardson); `ref-L2-310-two-rotations-cells1-12.txt` (sparse transfer
matrices and states for dim 15); the general notebook above (its cells 3–6, 9–10, 13, 16, 19, 24–25, 28 for the
SoV basis, K maps, tensors, loci, bracket decomposition, Gauss factors L₀D₀U₀). Write-ups to read first:
`L2-210-site1-rotation-brackets.md`, `L2-210-untwist-CG.md` (§1–§3, §5 for the reference CG numbers),
`L2-310-untwist-CG.md` (§2–§4: slice weights, reference numbers: [6,1,1] 46.6745…, [5,3,0] −36.8838…,
[5,2,1] −1325.40…, −334.33… at P1 — normalisation-dependent, compare only after matching the normalisation).

## 0. The family and what "untwist" means here

Companion-frame covector Φ(g) = e₁ᵀρ(g)⊗e₁ᵀ, g = exp(αE₁₂)exp(βE₁₃)exp(γE₂₃), on the locus
- [2,1,0]: 2(t−2)(β+αγ) = (t−3)α²  (t = θ₁₂/h; 43/16 at P1);
- [3,1,0]: β+αγ = α² (first row b = (1,α,α²)); keep the curve {243β = 280γ², 351α = 140γ} as a second option.
In the diagonal frame (where the untwist pipeline, the HW kets and the eigenstates live) the covector is
Φ_d(g;z) = (e₁ᵀρ(g)⊗e₁ᵀ)·(P_phys(z)⁻¹ ⊗ P_phys(z)⁻¹). With P_phys⁻¹ = L₀D₀U₀ (3×3 Gauss, exact, z-dependent;
`sVec`/`rsSol` of the reference cells), e₁ᵀρ(P⁻¹) = e₁ᵀρ(U₀) (D₀,₁₁ = 1) and the site-2 factor is the rotation
e₁ᵀρ(U₀); since ρ(U)⊗ρ(U) fixes every HW ket, the exact undo is

    ⟨Φ_d|h_ν⟩ = ⟨ e₁ᵀρ(g L₀ D₀) ⊗ e₁ᵀ | h_ν ⟩,   and   e₁ᵀρ(M) = (D_M)₁₁^{…}·e₁ᵀρ(U_M) for M = L_M D_M U_M (3×3 Gauss),

so the signal functional is the rotation by U_M, the upper unipotent factor of M(ε) := g(ε)·L₀(z(ε))·D₀(z(ε)) —
a 3×3 computation. Regime: z = 1 + εφ (φ = (1/3, −2/7, 5/11), control φ′ = (−2/5, 3/8, 7/11)); there is no w.
The Gauss data of P_phys⁻¹ have poles in ε (s₁ ~ ε⁻², s₂ ~ ε⁻¹), so g must be ε-dependent: g(ε) with
(α, β, γ) = (α₀ε^{p}, β₀ε^{q}, γ₀ε^{r}) on the locus. **Goal:** find (p, q, r) (and the free O(1) parameters
(α₀, γ₀), β₀ fixed by the locus) such that U_M(ε) → a generic unipotent at ε = 0: u₁₂, u₂₃, u₁₃ − u₁₂u₂₃ all O(1)
and varying independently over a grid of (α₀, γ₀) — then every weight of V is reached at ε⁰ (rank 2 on the
double weights over the grid) and K/S = O(ε) for every HW branch, and all CG data of Λ⊗Λ come from the
bracket Q-number at leading order.

## 1. Tasks (per rep; adjoint first; each its own cells with a `Grid` of checks)

**Task 0 — setup.** Generators, metric G, HW kets of V⊗V and the G-orthogonal projectors per HW sector
(reference cell 6 pattern; for [3,1,0] 10 HW states / 8 sectors), the exact sparse transfer matrices, `x0Of`
(factorised inverse, never `Inverse` at small ε), `branchesAt[z]` (HW-limit branches, normalised ⟨x₀|Ψ⟩ = 1),
the 3×3 Gauss factors L₀(z), D₀(z), U₀(z) symbolic in z (exact; verify e₁ᵀρ(P⁻¹) = e₁ᵀρ(U₀) at a rational z and
ρ(L₀)ρ(D₀)ρ(U₀) = ρ(P⁻¹) on V), the locus parametrisation g(α,γ) (β from the locus), the covector Φ_d(g;z),
and the exact undo check ⟨Φ_d|h⟩ = ⟨e₁ᵀρ(gL₀D₀)⊗e₁ᵀ|h⟩ on all HW kets at a rational z and rational g.

**Task 1 — the ε-scaling (3×3, exact symbolic ε).** With z = 1+εφ exact, compute L₀(ε), D₀(ε) and their
leading orders (print the ε-exponents of every entry). Then for M(ε) = g(α₀ε^p, β(·), γ₀ε^r)·L₀·D₀ compute the
Gauss factor U_M(ε) exactly (3×3 rational functions) and the ε-exponents of u₁₂, u₂₃, u₁₃ − u₁₂u₂₃ for all
(p, r) ∈ {−3,…,3}² (β fixed by the locus; also try β as the free one with α fixed by the locus if the locus is
not solvable for β). Tabulate; select the (p, r) with all three exponents equal (any common value: it is an
overall power) and leading coefficients that depend on (α₀, γ₀) non-trivially (rank of the map (α₀,γ₀) ↦
(u₁₂⁰, u₂₃⁰, (u₁₃−u₁₂u₂₃)⁰) = 2 at a generic point) — report ALL admissible (p, r). If none has all three at a
common order, report the best (which direction lags, by how many orders) — the analogue of the old E₁₃
degeneration — and continue with it (the deficient slices then need the next order; do not fit anything).
Also do the control: generic rotation (no locus) — what (p,r,q) work there (the locus is a constraint, the
control shows whether it costs a direction).

**Task 2 — ε-structure of the undone covector (exact, symbolic ε).** With the admissible scaling: φ_u(ε) =
e₁ᵀρ(g(ε)L₀D₀)/(overall power); per weight of V: leading order and the rank of the ε⁰ coefficient vectors over
a grid of 6 (α₀, γ₀) points (and φ′). Expectation: every weight at ε⁰, rank = multiplicity (2 on (1,1,1) for
the adjoint; on (2,1,1), (1,2,1), (1,1,2) for [3,1,0]). Control: g fixed (ε-independent) → expect the depth
table (ε^depth per weight).

**Task 3 — chain side: signal vs contamination.** ε ∈ {1,2,4}·10⁻⁶, grid of 6 (α₀,γ₀) × φ (then φ′ at one ε):
branches at z = 1+εφ; R_C = ⟨Φ_d|Ψ_B⟩, S = ⟨Φ_d|P_HWΨ_B⟩ = ⟨φ_u⊗e₁ᵀ|P_HWΨ_B⟩·(prefactor) (check the exact undo
numerically, 10⁻¹⁵⁰), K = R_C − S; ε-exponents of |S|, |K|, max|K/S| over the grid per branch. **Expectation to
confirm or refute: K/S = O(ε) for EVERY HW branch** (adjoint: 6 incl. singlet and both octets; [3,1,0]: 10).
If a branch fails: `showContam`-type breakdown (which site-1 weights carry K at which order); no fitting.

**Task 4 — CG at leading order from the overlap alone.** C := R_C/((prefactor)·[φ_u(0)]_slice) for 1-dim slices:
grid spread and error vs the chain slice at two ε's (∝ ε), Richardson limits vs chain (sanity against the
reference numbers after matching normalisation: adjoint [4,1,1] 46.67453891, [3,3,0] 10.92853892, singlet
−5118.538314 in the `L2-210-untwist-CG` normalisation — compute the chain slice in THIS notebook and compare to
it, quote the references as a cross-check only). Dim-2 slices: 2×2 systems from two grid points (adjoint octets:
mixing R vs −1.694, 1.096; [3,1,0]: [5,2,1]×2, [4,3,1]×2). Conclusion cell per rep: how many of the HW CG data are
read at leading order (expect all).

**Task 5 — the Q-number at small ε (one grid point, ε = 10⁻⁶, HW-limit branches).** The bracket expression:
R_C/⟨x₀|Ψ_B⟩ · M⁽¹⁾_vac = Σ_{c₂} λ₂^{(c₂)} Σ_c Σ_{wedges} Bk₂[F,G]·Q³(θ₁−ch) with the measures F, G decomposed
exactly from the tensors at g(ε) (general notebook cells 16/25 machinery; the K maps are z-independent, the
tensors depend on g only — so F, G at g(ε) are exact rationals in ε), the dual Q's of the branch from the
Baxter equation at this z (reference cell 9's `auditMonicQ`, degrees μ_j). Verify on the branches (relative
deviation at 200 digits), report the ε-exponents of F, G entries and of the bracket determinants (where the
cancellation sits). State the final cost: brackets count (4 adjoint, 5 [3,1,0]), Q-values needed.

**Task 6 — markdown summary** `### site-1 rotation families: untwist result (2026-10-10)` per rep: the
admissible scalings (p,q,r) and the leading Gauss parameters, the per-weight order/rank table (+ control), the
K/S table with ε-exponents, the CG errors and Richardson limits, the Q-number check and cost, surprises. Save.

## 2. Rules and gotchas

As before: exact where possible, 200 digits on the chain; time one state/one point before loops ([3,1,0]: 56
HW-sector states per ε ≈ 3 s in the earlier notebooks; symbolic-ε Gauss of 3×3: trivial); nothing over ~3
minutes without a timed smaller run; markdown result cells under code cells; save after each task; if VS Code
will not execute cells of the new file, run cell sources in K1 and store printed results as markdown (say so).
Do not edit other notebooks except Job A's notebook as specified, no `paper/*.tex`, no git. If Task 1 finds no
admissible scaling for the locus family, say so clearly, run the control, and stop after Task 2 with the
per-weight table of the best scaling (the main session decides).
- Shared kernel K1; first cell wipes other state. Never `*)` in comments. No `D`,`E`,`N`,`I`,`C`,`K` symbols.
- `x0Of` with the factorised inverse; `Inverse[Pph1]` at ε = 10⁻⁶ loses digits.
- `Do[DownValues[f] = …]` purges nothing; purge memoized definitions with literal symbols.
- `Max[{}]` = −∞; `SetPrecision` before `Norm` of near-zero vectors; exact zero tests with `===` only on rationals.
- The locus is t-dependent for the adjoint (43/16 at P1); keep h, θ as symbols where cheap so the scaling
  analysis is point-independent, else at P1 and once at P3 = (1/2, (3/7, −1/4)).
