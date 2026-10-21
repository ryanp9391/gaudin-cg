# Spec: untwisting the two-rotation det₆ × det₂ formula and reading off CG coefficients (adjoint, 2026-10-10)

Delegated to a wolfbook-builder subagent (sonnet); planned and reviewed by the main session. **New notebook:** `Paul/Mathematica/new-conventions/L2-210-untwist-CG.wb` (create it; self-contained; first code cell clears the shared kernel as every new-conventions notebook does). Companion reference cells (verbatim sources to copy from) are in the scratchpad file named in the task message.

## 0. What is known (read `Undo-site2-rotation-notes.md` §1–§4 and `Two-rotations-det-overlap.md` §5.6b, §5.8 if anything below is unclear)

gl(3), L = 2, Λ = (2,1,0) on both sites, new conventions, diagonal frame (ket twist diag(z)), GT basis with plain generators `JJ` (exact rationals), HW ket e₁ first; covectors are rows acted on from the right; `Ee[i,j]` with i<j raise.

**The family.** For a ket twist z, a bra-twist parameter w₁ and a free parameter y, with w₂(y) := y(θ₁₂−2h)/(θ₁₂−h):
- v₂ = Ω_{(w₁,y)} = e₁ᵀ ρ(n), n = [[1,r₁,r₂+r₁r₃],[0,1,r₃],[0,0,1]] ∈ N₊, r_k = `rOm`-formulas with w = (w₁, y, w₃) (w₃ drops out),
- v₁ = v₂ ρ(Δ), Δ = diag(1/(w₂(y) − z_k)) (`ρdiag[Dy[w2Of[y]]]`),
- **overlap identity (finite twist, state-independent to 10⁻⁵⁰):**
  ⟨v₁⊗v₂|Ψ_B⟩/⟨x₀|Ψ_B⟩ = c · det₆[Q_B, Q_vac,(w₁,w₂(y))]/Ψ_Q^B(x₀) · det₂^{(1)}[Q_B; y]/minor₁(B),
  det₆ = the 6×6 all-pole bracket determinant with the vacuum direct Q's of the twist (w₁, w₂(y)) (`det6[B, Fvac]`), det₂^{(1)} = ladder-1 su(2)-like 2×2 with insertion y^{−(x−θ₁)/h} (`su2det1[B, y]`), minor₁ = Q²(θ₁−2h)Q³(θ₁−3h) − Q³(θ₁−2h)Q²(θ₁−3h).
- **Undo (exact for HW kets):** ⟨v₁⊗v₂|HW⟩ = s⟨φ⊗e₁ᵀ|HW⟩, s = (v₂)₁ = 1, φ = v₁ρ(n)⁻¹ = δ^Λ e₁ᵀρ(n′), δ_k = 1/(w₂(y)−z_k), δ^Λ = δ₁²δ₂, n′ = Δ⁻¹nΔn⁻¹ = [[1,a′,b′],[0,1,c′],[0,0,1]] with
  a′ = r₁(δ₂−δ₁)/δ₁, c′ = r₃(δ₃−δ₂)/δ₂, b′ = r₁r₃(δ₁−δ₂)/δ₁ − (r₂+r₁r₃)(δ₁−δ₃)/δ₁, and b′/(a′c′) = (y−w₂(y))/(y−z₂).
  Weight components of φ/δ^Λ (cell 118 of `L2-rotation-sov.wb`): (2,1,0): 1; (1,2,0): a′; (2,0,1): c′; (1,1,1) [dim 2]: b′·e₁ᵀE₁₃-type + a′c′·e₁ᵀE₂₃E₁₂-type; (0,2,1): a′²c′p + a′b′q; (1,0,2): b′c′; (0,1,2): 3b′(a′c′−b′).
- **The CG target.** For each HW ket h_ν of V⊗V (6: [4,2,0], [4,1,1], [3,3,0], [3,2,1]×2, [2,2,2]) the family measures the slice ⟨m⊗top|h_ν⟩, m ∈ V_{ν−Λ}. The twisted eigenstate Ψ_B of the HW-limit branch of ν is **not** a HW ket; its overlap splits exactly as ⟨v₁⊗v₂|Ψ_B⟩ = S + K with S = ⟨v₁⊗v₂|P_HWΨ_B⟩ (signal) and K = ⟨v₁⊗v₂|(1−P_HW)Ψ_B⟩ (contamination), as in `Clean-CG-combinations.md` §2–3.

**The untwisting regime (decided by the main session; the reason is in §4 below).** Everything untwists together:

  z = 1 + εφ,  w₁ = 1 + εω₁,  w₂ := 1 + εω₂,  y := w₂ (θ₁₂−h)/(θ₁₂−2h)  (so that w₂(y) = w₂ exactly),

with fixed rational directions φ, ω₁, ω₂ and ε → 0 (rational paths, not exponentials; same first-order direction). Then (checked by the main session with sympy at P1): a′ = (ω₁−φ₁)(φ₂−φ₁)/((ω₁−φ₂)(ω₂−φ₂)) **exactly ε-independent**, c′ = r₃(ε)(φ₃−φ₂)/(ω₂−φ₃) = O(1), b′ = O(1), and

  b′/(a′c′) = 1 + ε (θ₁₂/h − 2)(φ₂ − ω₂) + O(ε²),

so at ε = 0 the undone covector is φ₀ = e₁ᵀ e^{a′E₁₂} e^{c′E₂₃} (b′ = a′c′): components of weights (2,1,0), (1,2,0), (2,0,1), (1,0,2) and ONE direction of (1,1,1) are O(1); the second (1,1,1) direction, (0,2,1) and the singlet direction (0,1,2) are O(ε). (At fixed y instead, the components are O(ε⁰,ε¹,ε¹,ε²,ε³,ε³,ε⁴) — the old table — and the contamination dominates; that regime is NOT used here.)

**Numeric point.** P1: h = 1/2, θ = (2/7, −3/5), so t := θ₁₂/h = 62/35 and y₀ := (θ₁₂−h)/(θ₁₂−2h) = −27/8. Directions: φ = (1/3, −2/7, 5/11) (primary), φ′ = (−2/5, 3/8, 7/11) (control); ω₁ ∈ {3/4, −5/3, 7/2}; ω₂ ∈ {−1/5, 9/7, −12/5, 5/3, 2/9, −7/4}. (Never let ω₂ or ω₁ equal a φ_k: Δ and the r's have poles there; never let y equal a z_i.) ε ∈ {10⁻⁶, 2·10⁻⁶, 4·10⁻⁶} for Richardson extrapolation; working precision pr = 200 (check `Precision` of every reported number; a "0" with zero accuracy is a trap, see `Clean-CG-combinations.md` §7).

## 1. Notebook setup (Section 0 of the notebook; copy, do not re-derive)

Copy verbatim from the reference file, adapting only what is said here:
1. GT generators, Λ, `Ee`, `E1`, `E2`, `wtBasis`, `nuP` (rotation-sov cell 2). Set `pr = 200` and keep `hN = 1/2`, `thN = {2/7, -3/5}`; do NOT hard-wire `zN`, `wN` into anything below — every object below is a function of the twists.
2. Monodromy/transfer polynomials `Tm`, `TMf`, `qminor`, `tAt`, `tPoly[M, a]` (cell 3; interpolation points `uS` — keep that name only for these points and never reuse it). Define `tZof[zz_] := Table[tPoly[DiagonalMatrix[N[zz, pr]], a], {a, 3}]` (memoize on zz). Check the quantum determinant for one numeric zz.
3. Auxiliary/physical Gauss factors `Paux`, `rsSol`, `Gc`, `χ` (cell 6), `oneSiteVac[zz, ww]` (cell 98), `ρdiag`, `Dy` (cell 10; make `Dy[y1, zz]` take zz), the physical-frame x₀ as a function of zz **using the factorized inverse** (pattern `auditX0` in diag-EH cell 96: products of MatrixExp of the negated exponents and the inverse diagonal; never `Inverse[Pph1]` at small ε), the u12 machinery `u12rs`, `u12nMat`, `u12rN`, `u12rNinv`, `u12v2`, `u12v1`, `u12phi` with zz, w₁, y as arguments and `w2Of` replaced by the regime of §0 (cell 114), and the exact HW kets `u12HW` with their weights (cell 115) plus orthonormal per-sector bases/projectors `hwB[μ]`, `PHW[μ]` (Minimal cell 8 pattern, but built from the exact `u12HW`).
4. Q-machinery: poles/weights `poles`, `xP`, `wgt`, brackets `rowG`/`det6[B, F]` (cell 8) and ladder-1 `wgt1`, `br1`, `su2det1`, `minor1`, `psiQ` (cells 7, 68), rewritten so that the state's dual Q's and the twist are passed in (an association per state: `"qu"` = three monic dual polynomials in X, `"z"`), not read from globals `qu`/`zN`. Q-functions: Q^j(x) = z_j^{x/h} q^j(x) (dual, deg q^j = μ_j), vacuum direct Q_i^vac(x) = w_i^{−x/h} q_i^vac(x) with degrees (0, 2, 4).
5. **Q solvers (fixed degree, overdetermined least squares, no Chop):** copy `auditNu`, `auditVac` (vacuum eigenvalues of t_a for a twist, explicit), `auditQdet`, `auditOp`, `auditMonicQ` (diag-EH cell 95). The state's t₁, t₂ eigenvalue polynomials come from the eigenvector (largest component, as in rotation-sov cell 4); t₃ = `auditQdet`. Degrees: dual q^j: μ_j; vacuum direct q_i: 2(i−1). Report the Baxter residual of every Q you use (pattern: `qres` in diag-EH cell 99).
6. **Eigenstates and HW-limit branches** at a numeric twist zz (Minimal cell 8 + diag-EH cell 99 pattern): per HW sector μ ∈ {(4,2,0),(4,1,1),(3,3,0),(3,2,1),(2,2,2)}, `Eigensystem` of (t₁ + π t₂)(3/7) restricted to the sector, normalize ⟨x₀|Ψ⟩ = 1, raising residual rr = ‖(E1+E2)[1,2]Ψ‖+‖(E1+E2)[2,3]Ψ‖ normalized; the HW-limit branches are the dim(HW_μ) states with the smallest rr, copies ordered by Re of the eigenvalue. Store μ, copy, Ψ, eigenvalue, rr, t₁(u), t₂(u), q^j, and ‖P_HWΨ‖/‖Ψ‖. Check: 6 branches; the top branch is e₁⊗e₁ up to scale.

Timing rule: time one state / one parameter point before any loop; nothing over ~3 minutes without a timed smaller run.

## 2. Tasks (each its own cell(s) with a `Grid` of checks; keep `checks` as an Association)

**Task 0 — exact ε-structure of the undone covector in the regime (symbolic ε, rational everything else).** With `epsS` symbolic and (φ, ω₁, ω₂) = (φ, 3/4, −1/5): build φ(ε)/δ^Λ via `u12phi` with the §0 regime (δ^Λ = δ₁²δ₂, divide it out), `Series` each of the 8 components to O(ε³). Report per weight: leading order and the leading coefficient; verify (i) a′ is exactly ε-free and equals (ω₁−φ₁)(φ₂−φ₁)/((ω₁−φ₂)(ω₂−φ₂)), (ii) b′/(a′c′) = 1 + ε(t−2)(φ₂−ω₂) + O(ε²) with t = θ₁₂/h (at these numbers the O(ε) coefficient is 24/1225), (iii) φ₀ := φ(0)/δ^Λ equals e₁ᵀ·MatrixExp[a′₀Ee[1,2]]·MatrixExp[c′₀Ee[2,3]] exactly, (iv) the expected orders (2,1,0):0, (1,2,0):0, (2,0,1):0, (1,1,1): 0 with rank 1 of 2 at ε⁰ and rank 2 including the ε¹ coefficient, (0,2,1): 1, (1,0,2): 0, (0,1,2): 1. Then, over the 3×6 grid of (ω₁, ω₂), the per-weight ranks of the ε⁰ coefficient vectors and of the ε¹ coefficient vectors (exact). Also record the O(ε) coefficient of the (0,1,2) component in factored form (it should be ∝ (a′c′₀)²(φ₂−ω₂)).

**Task 1 — chain side: signal vs contamination in the regime.** For ε ∈ {10⁻⁶, 2·10⁻⁶, 4·10⁻⁶}, φ primary, all 18 (ω₁, ω₂): build v₁, v₂ (exact rationals in ε — fine, or N[.., pr]), the 6 HW-limit branches at z = 1+εφ, and for each branch
- R_C := ⟨v₁⊗v₂|Ψ_B⟩/⟨x₀|Ψ_B⟩ (with ⟨x₀|Ψ⟩ = 1 this is just the overlap),
- S := ⟨v₁⊗v₂|P_HW Ψ_B⟩, K := R_C − S; check S = s⟨φ⊗e₁ᵀ|P_HWΨ_B⟩ with φ = `u12phi` (exact undo, relative 10⁻^{pr−20}),
- ε-exponents (from the three ε's, as `εpow` in Minimal cell 15) of |S|, |K|, ‖Ψ_B‖, and the ratio κ := R_C/S = 1 + K/S.
Report a Grid (one row per branch, κ at the three ε's for two representative (ω₁,ω₂) points, plus max over the grid of |K/S|). **Expectation to confirm or refute:** K/S = O(ε) (i.e. κ → 1 linearly in ε) for [4,1,1], [3,3,0] and both octets (whose S at ε⁰ is the single functional ∝ a′c′), and κ ≠ 1 at ε → 0 for the singlet (S and K both O(ε) relative to the ε⁰ scale of the other weights). If an expected-clean state shows κ ↛ 1, stop and report (do not proceed to Task 2 for that state).
Also a `showContam`-type breakdown (Minimal cell 15 pattern, with R := ρ(n)⊗ρ(n) of the site-2 vacuum rotation, so that K = Σ_ν φ^{(ν)}·[(R)Ψ⊥]_ν after undoing) for the singlet and for octet copy 1: which site-1 weights ν carry the contamination and at which ε-powers.

**Task 2 — Q side at small ε, and the CG numbers from the formula.** For ε ∈ {10⁻⁶, 4·10⁻⁶} and all 18 (ω₁, ω₂):
- vacuum direct Q's of the twist (w₁, w₂) (w₃ irrelevant; use w₃ = 1+ε·(2/9) or anything) from `auditVac` + `auditMonicQ` (degrees 0,2,4); the branches' dual Q's (degrees μ_j); report max Baxter residuals.
- R_Q(B) := det6[B, Fvac]/psiQ[B] × su2det1[B, y]/minor1[B] with y = w₂·y₀; check R_C(B)/R_Q(B) is state-independent over the 6 branches (and over ALL states of the five HW sectors, 21 states, at one parameter point) — report max |ratio/ratio_top − 1| and its Precision. This is the finite-twist identity at small ε; expect ~10⁻¹⁵⁰ or better.
- Define the pure Q-side number Qnum(B) := R_Q(B)/R_Q(top) × ⟨v₁⊗v₂|e₁⊗e₁⟩/⟨x₀|e₁⊗e₁⟩ (top state = e₁⊗e₁ exactly; its dual Q's are those of the top branch). Check Qnum = R_C to the same precision.
- **CG extraction from Qnum alone (leading order).** For [4,1,1] (site-1 weight (2,0,1), 1-dim) and [3,3,0] ((1,2,0), 1-dim): with the exact ε = 0 covector φ₀(ω) of Task 0, the prediction is Qnum(ω) → φ₀,μ(ω)·C_B(ε) + O(ε·…) where C_B(ε) := ⟨h_μ ⊗ top|P_HWΨ_B⟩-slice, i.e. C_B^{Q}(ω, ε) := Qnum(ω, ε)/φ₀,μ(ω) must be ω-independent up to O(ε) and equal to the chain-side C_B(ε) = (P_HWΨ_B) component to O(ε) relative. Report the spread over the 18 ω's and the relative error vs the chain side, at both ε's (expect both ∝ ε). For the two octet copies: the ε⁰ functional is a′c′₀·(u+u′)·P_HWΨ_B; check the same ω-independence of Qnum/(a′c′₀) and compare with the chain-side (u+u′)-component. For the 27: trivial check Qnum = δ^Λ-type constant.
- Richardson: from ε = 10⁻⁶, 2·10⁻⁶, 4·10⁻⁶ (run the third ε only for this step, at one (ω₁,ω₂)), extrapolate C_B^{Q}·ε^{k_B} to ε → 0 (k_B = the ε-exponent of ‖Ψ_B‖ found in Task 1) and compare with the chain-side extrapolation; report the agreement.

**Task 3 — diagnostic for the O(ε) layer (octet second direction, singlet).** Do not try to solve it; measure it. For the singlet and both octet copies, at ε = 10⁻⁶ and 4·10⁻⁶ and all 18 ω's:
- top-killed quantities: R_C^{tk} := R_C − δ^Λ s ⟨v₂⊗v₂|Ψ_B⟩/⟨x₀|Ψ_B⟩ (the second term is the plain vacuum overlap of the twist (w₁, y), itself a Q-number: det6 with the vacuum Q's of (w₁, y) over psiQ — verify this equality to 10⁻⁴⁰ at one point), and K^{tk} := K − (same subtraction);
- the ratio ρ_B(ω) := K^{tk}/S at both ε's; test (a) ε-independence (structural, as in Clean-CG), (b) φ-independence (repeat at φ′ for 4 of the 18 ω's), (c) print ρ_B(ω) as rationalized numbers (`Rationalize[.., 10^-20]`) over the grid, and the ratio (ρ_B at ω)/(ρ_B at ω′) for pairs differing only in ω₂ and only in ω₁.
Report the table; no fitting beyond that.

## 3. Reporting

End the notebook with a markdown cell `### Untwist result (2026-10-10)` listing, per task, the key numbers (residuals with their Precision, the ε-exponent tables, the κ tables, the CG-extraction errors, the ρ_B tables) and every structural surprise. Save the notebook (`wolfbook_saveNotebook`). Then report back concisely in the same structure. Do not edit any other notebook; do not touch `paper/*.tex`; do not commit. Everything must be visible in notebook cells (no scratch files). If a task's expectation is refuted or a check fails and you cannot fix it within two attempts, stop at that point, keep what works, and report what you saw (numbers, not adjectives).

## 4. Why this regime (for the record; the main session's reasoning)

With Ψ_B = ε^{−k}(h + εψ₁ + …)/norm, the contamination starts at ε¹ relative to the covector's own components. At fixed y the covector's components below the top are themselves O(ε^d) (d = depth of the weight), so contamination is at the same order as the signal for every state (the Clean-CG situation). In the regime of §0 the resolvent torus Δ = ε⁻¹diag(1/(ω₂−φ_k)) is a genuine non-scalar torus element at ε = 0, so the two rotations v₁, v₂ stay distinct and the undone covector has O(1) components; only the E₁₃-direction b′−a′c′ degenerates (the ε = 0 shadow of the old b′ = 0 blind spot, because the vacuum's r₂/(r₁r₃) and the torus Δ are built from the same φ). Hence: 27, [4,1,1], [3,3,0] and one octet functional per copy are clean at leading order; the octet mixing and the singlet live at O(ε) together with the first contamination and need a cleaning step (Task 3 measures it).

## 5. Gotchas (all have bitten before)

- The shared kernel: your first cell's `ClearAll` wipes the live state of the other new-conventions notebooks (files are intact). Fine, but do not then evaluate cells of other notebooks.
- Never write the two characters `*)` inside a comment (e.g. "y*)" closes it; write y-star).
- Do not redefine `uS` (interpolation points) or use `D`, `E`, `N`, `I`, `C` as symbols.
- Memoized functions cache bad results: after fixing a definition, re-run the defining cell and everything downstream.
- Semicolons between statements in multi-line cells.
- `Max[{}]` → −Infinity; floor at 0 (`maxc`).
- Precision: compute at pr = 200, report `Precision` of residuals; use `SetPrecision` before `Norm` of near-zero vectors (`auditNorm`).
- Wolfbook replies truncate at `{`; read long results via single-line `StringReplace[ToString[.., InputForm], "\n" -> " "]`.
- Ω_{(w₁,y)} degenerates at y = z_i (0/0); the r's have poles at w₁ = z_k; Δ at w₂ = z_k. With the directions above none of these occur, but if you change a direction, check.
