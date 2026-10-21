# Spec: the [3,1,0] k = 2 point (point A) through the untwist pipeline (2026-10-10)

Delegated to a wolfbook-builder subagent (sonnet); planned and reviewed by the main session. **New notebook:** `Paul/Mathematica/new-conventions/L2-310-untwist-CG.wb` (self-contained; first code cell clears the shared kernel). This is the [3,1,0] version of `L2-210-untwist-CG.wb` / `L2-210-untwist-CG.md` — read that write-up first (§0–§5), then `Undo-site2-rotation-notes.md` §5.2–5.4 and §6 (point A, its dressing, its det₂ identification, the Q-dictionary). Verbatim sources of the cells to copy are in the scratchpad file named in the task message.

## 0. What is known

gl(3), L = 2, Λ = (3,1,0) on both sites (dim 15, 225 states), new conventions, diagonal frame, plain GT generators `JJ` (exact sparse rationals: `Ee`, `tL2`, `TL2` of `L2-310-two-rotations.wb` cells 2–3). 10 HW states in 8 irreps: [6,2,0] (top), [6,1,1], [5,3,0], 2×[5,2,1], [4,4,0], 2×[4,3,1], [4,2,2], [3,3,2]; the site-1 weight carrying the CG slice of each is ν−Λ: (3,1,0), (3,0,1), (2,2,0), (2,1,1) [dim 2], (1,3,0), (1,2,1) [dim 2], (1,1,2) [dim 2], (0,2,2).

**Point A (k = 2)** of `L2-310-two-rotations.wb` (cells 37–39, 44–63): with y_A = w₂(θ₁₂−h)/(θ₁₂−2h) (the adjoint's y-star),
- site 2: v₂ = e₁ᵀρ(n_A), n_A the N₊ Gauss element of `rOm[z,{w₁,y_A,w₃}]`;
- site 1: v₁ ∝ e₁ᵀ(1 + βE₁₂ + γ(E₁₃ + ½E₁₂E₂₃)) ρ(n_A) ρ(Δ), Δ = diag(1/(w₂−z_k)) (`Dy[w₂]`), with the exact closed forms (cell 63)
  β = (z₁−z₂)(w₁−z₁)w₂ / (2z₂(w₁−z₂)(w₂−z₁)),
  γ = (θ₁₂−h)(z₁−z₃)(z₂−z₃)(w₁−z₁)w₂² / (4(θ₁₂−2h)z₂z₃(w₁−z₃)(w₂−z₁)(y_A−z₃));
- overlap (225 states, 5·10⁻⁴⁵): ⟨v₁⊗v₂|Ψ⟩/⟨x₀|Ψ⟩ ∝ det₆[Q_B, Q_vac,(w₁,w₂)]/Ψ_Q × det₂^{(0,1,2)}[Q_B; y_A]/W(3,4), where det₆ is the 8-pole (n = 0..3 on both ladders) bracket determinant with the vacuum direct Q's of twist (w₁,w₂) (degrees 0, 4, 6), Ψ_Q = Π_s[Q²(θ_s−3h)Q³(θ_s−4h) − Q³(θ_s−3h)Q²(θ_s−4h)]Q³(θ_s−h), det₂^{(P)}[y] = Σ_{n,m∈P} w_n w_m y^{n+m}[Q²(θ₁−nh)Q³(θ₁−(m+1)h) − Q³(θ₁−nh)Q²(θ₁−(m+1)h)] with one-site weights w_n = 1/Π_{m∈P∖n}(m−n)h on P = {0,1,2}, and W(a,b) = Q²(θ₁−ah)Q³(θ₁−bh) − Q³(θ₁−ah)Q²(θ₁−bh) (cells 18, 50, 56, 67);
- undo (exact): ⟨v₁⊗v₂|HW⟩ = s⟨φ_A⊗e₁ᵀ|HW⟩, φ_A = v₁ρ(n_A)⁻¹; at finite twist φ_A spans all 15 site-1 directions and determines all 10 HW states (cell 39).

**The regime** (as for the adjoint): z = 1+εφ, w₁ = 1+εω₁, w₂ = 1+εω₂, y_A = w₂·y₀ with y₀ = (θ₁₂−h)/(θ₁₂−2h), rational paths. In this regime β = (φ₁−φ₂)(ω₁−φ₁)/(2(ω₁−φ₂)(ω₂−φ₁)) + O(ε) is O(1) while **γ = O(ε)** (its numerator carries three ε's, its denominator two), so at ε = 0 the dressing is e₁ᵀ(1+β₀E₁₂) only, and, as for the adjoint, n′ = Δ⁻¹n_AΔn_A⁻¹ has b′/(a′c′) = 1 + ε(t−2)(φ₂−ω₂) (t = θ₁₂/h). Expect: every E₁₃-type direction is O(ε); one direction of each dim-2 weight and the weight (0,2,2) are O(ε). Verify, don't assume.

Numeric point P1: h = 1/2, θ = (2/7, −3/5), t = 62/35, y₀ = −27/8. Directions φ = (1/3, −2/7, 5/11), φ′ = (−2/5, 3/8, 7/11), ω₁ ∈ {3/4, −5/3, 7/2}, ω₂ ∈ {−1/5, 9/7, −12/5, 5/3, 2/9, −7/4}, w₃ = 1+ε·5/13; ε ∈ {10⁻⁶, 2·10⁻⁶, 4·10⁻⁶}; pr = 200 (check `Precision` of every reported number). Never let ω₁, ω₂ equal a φ_k.

## 1. Setup (copy; parametrize everything by the twists)

1. Generators, `Ee`, `E1`, `E2`, `wtBasis`, `secOf`, `nuP`, the sparse exact `tL2[a, x, M]` and polynomial-in-u transfer matrices by exact Lagrange interpolation on 7 points (cell 16 of the 310 notebook uses `uPts`; keep that name), as functions of a numeric z (memoize on z; **purge with literal symbols**, never `Do[DownValues[f] = …]`, see `L2-210-untwist-CG.wb` cell 36).
2. Invariant metric G for Λ = (3,1,0) (the plain `JJ` are not Hermitian; solve M⁻¹E_ijᵀM = E_ji with M diagonal, M₁ = 1, as `MsolEH` in `Minimal-ingredients-EH.wb` cell 2) and the G-orthogonal projector P_HW per HW sector from the exact HW kets (common null space of the total raising generators per sector).
3. x₀ with the factorized inverse (adjoint untwist cell 5), `oneSiteVac`, `ρdiag`, `Dy[y, zz]`, the u12 machinery with `w2Of` from the regime, the dressing g(β, γ) = e₁ᵀ(1 + βE₁₂ + γ(E₁₃ + ½E₁₂E₂₃)) with the closed forms above, `vA1[zz,w1,w2]` = g·ρ(n_A)·ρ(Dy[w₂]), `vA2` = e₁ᵀρ(n_A), `phiA` = vA1·ρ(n_A)⁻¹, and `reg[eps, φ, ω1, ω2]`.
4. Finite-twist regression test before anything small-ε: at P1's old z = (3/2,5/4,2/3), w = (5/2,7/3,3/5) (y_A = −63/8), rebuild v₁, v₂ from the closed forms and check (i) ⟨v₁⊗v₂|HW⟩ = s⟨φ_A⊗e₁ᵀ|HW⟩ on the 10 HW kets exactly, (ii) the det identity ⟨v₁⊗v₂|Ψ⟩/⟨x₀|Ψ⟩ ∝ det₆·det₂^{(0,1,2)}[y_A]/(Ψ_Q W(3,4)) state-independent over all 225 states at 60 digits (expect ~10⁻⁴⁵ as in cell 56). If (ii) fails, stop and report (the closed-form dressing or a convention is off).
5. Q machinery with the state's dual Q's and the twist passed in: fixed-degree monic solver (`auditMonicQ`-type, degrees: dual q^j: μ_j; vacuum direct q_i: 0, 4, 6), 8-pole det₆, det₂^{(0,1,2)}, W(3,4), Ψ_Q, exact-arithmetic Baxter residuals.
6. `branchesAt[zz]`: all states of the 8 HW sectors at numeric z (sector dimensions up to the (3,3,2) one; time one sector first), ⟨x₀|Ψ⟩ = 1, raising residual rr, HW-limit branches = the dim(HW_μ) smallest-rr states per sector, copies ordered by Re eigenvalue; store t₁, t₂ polynomials, dual Q's, Baxter residual, ‖P_HWΨ‖/‖Ψ‖ (G-norm). Expect 10 branches; top = e₁⊗e₁.

## 2. Tasks

**Task 0 — ε-structure of φ_A in the regime (exact, ε symbolic, rational φ, ω).** Per site-1 weight (12 weights): leading ε-order and the leading coefficient; ranks of the ε⁰ and ε¹ coefficient vectors over the 3×6 grid; verify γ = O(ε), b′/(a′c′) = 1 + ε(t−2)(φ₂−ω₂), and record φ_A(0)/δ^Λ in closed form (is it e₁ᵀ(1+β₀E₁₂)e^{a′₀E₁₂}e^{c′₀E₂₃}?). Expectation to check: weights (3,1,0), (2,2,0), (1,3,0), (3,0,1), (2,0,2), (1,0,3)?, and one direction each of (2,1,1), (1,2,1), (1,1,2) at ε⁰; the second directions, (0,2,2), (0,3,1)?, (0,1,3)? at ε¹ (or ε²: report what you find; anything at ε² needs flagging because it would require a third order).

**Task 1 — chain side.** For the 10 branches, at the three ε and 18 ω: R_C, S = ⟨v₁⊗v₂|P_HWΨ⟩, K; the exact-undo check; ε-exponents of S, K, ‖Ψ‖; κ = R_C/S; max |K/S| per branch. Classify: clean at leading order (κ−1 ∝ ε) vs contaminated (κ−1 = O(1)). Contamination breakdown by site-1 weight for each contaminated branch.

**Task 2 — Q side.** Vacuum Q's and branch Q's at small ε (Baxter residuals); state-independence of R_C/R_Q over the 10 branches and over all states of the HW sectors at one ω (report the number of states and the residual); Qnum = R_Q/R_Q(top)·⟨v₁⊗v₂|e₁e₁⟩/⟨x₀|e₁e₁⟩. Leading-order extraction for every clean branch: C^Q = Qnum/(δ^Λ·fac(ω)) with fac the ε⁰ covector coefficient of the branch's site-1 weight (for dim-2 weights the one direction seen at ε⁰): spread over ω and error vs the chain-side slice at both ε; Richardson limits Q vs chain.

**Task 3 — the O(ε) layer, diagnostics (no closed-form fitting yet).** Top-kill with the (w₁,y_A)-vacuum det₆ (verify it equals δ^Λ s⟨v₂⊗v₂|Ψ⟩/⟨x₀|Ψ⟩). For every branch that is not fully determined at leading order (expected: the second directions of [5,2,1]×2, [4,3,1]×2, [4,2,2]; and [3,3,2] at (0,2,2)): F = Qnum^tk/(δ^Λ·fac₀) = F₀ + εF₁ + …, three-ε Richardson; is F₁/F₀ a straight line in ω₂? ω₁-independent? slope φ-independent (φ′ at 4 ω's)? For (0,2,2)-type branches with F₀ ≡ 0: is K/S a constant over ε, ω, φ (like the adjoint singlet)? Tabulate the structural numbers (slopes, constants) at P1 with the chain-side mixing R (defined per weight as the ratio of the second-direction to the first-direction component of the CG slice; state the direction convention you use, u and u′, explicitly) so that Δ := slope + (t−2)R can be read off. Report which quantities look structural (ε-, ω-, φ-independent).

Stop rules as before: if the finite-twist regression (setup 4) fails, or a Task-0 expectation is off by more than one ε-order, stop and report with numbers. Time one state/one point before each loop; nothing over ~5 minutes without a timed smaller run (225-dim sparse exact objects: the 310 notebook builds its T's in seconds, its 225 eigenstates at 60 digits in ~10 s; at 200 digits expect a few times more).

## 3. Reporting

Markdown cell `### [3,1,0] point A untwist result (2026-10-10)` with every table; save the notebook; report back per task with numbers (residuals with Precision, ε-exponent table, κ table, extraction errors, the O(ε) structure table). Do not edit other notebooks, `paper/*.tex`, or git. Note: VS Code may not execute cells of a new file (it did not for the adjoint notebook); if so, run the cell sources in the kernel as the adjoint agent did, and say so.
