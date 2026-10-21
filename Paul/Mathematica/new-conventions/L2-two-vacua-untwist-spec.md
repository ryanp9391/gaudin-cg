# Spec: Stage D (chain side) of `Two-site-rotation-brackets-plan.md` — untwisting ⟨Ω_w ⊗ Ω_{w′}| with w, w′ → 1 independently; reach, K/S, CG at leading order (2026-10-11)

Delegated to a wolfbook-builder subagent (sonnet); planned and reviewed by the main session. **New notebook:**
`Paul/Mathematica/new-conventions/L2-two-vacua-untwist.wb` (self-contained, rep-generic with `Λ` as the switch; first
code cell `Quiet[ClearAll["Global`*"]]`; [2,1,0] first, then [3,1,0]). Create/run as the previous Stage notebooks (JSON
build, sources evaluated in K1, printed results as markdown cells; say so). This stage is the chain side of plan §5; it
does NOT depend on the bracket formula (Stages B/C) — the Q-side number is a separate later task.

Read first: plan §0 (last paragraph) and §5; `L2-210-untwist-CG.md` §1–§3 (regime, exact undo, per-weight table,
K/S, CG extraction, reference numbers); `L2-210-localising-span-spec.md` §6 (the exact commutator formula and why the
old family degenerated); `L2-site1-rotation-untwist-spec.md` §0 and Tasks 0–4 (the pipeline this spec adapts: HW kets
and G-orthogonal projectors, `x0Of`, `branchesAt`, signal/contamination, leading-order CG, Richardson). Cell-source dumps
(scratchpad, UTF-8): `ref-210-untwist-CG-all.txt` (`L2-210-untwist-CG.wb`, 45 cells: cells 1–12 setup + regime + exact
undo + symbolic-ε covector; 13–21 K/S and CG), `ref-rotation-sov-all.txt` (cell 98 `oneSiteVac[z,w]`: the diagonal-frame
Gauss parameters r₁ = (w₁−z₁)/(w₁−z₂), r₂ = (w₁−z₁)(z₁−w₂)(z₂−z₃)/((z₁−z₂)(w₁−z₃)(w₂−z₃)), r₃ = (w₁−z₂)(z₂−w₂)(z₁−z₃)/
((z₁−z₂)(w₁−z₃)(z₃−w₂)); cell 25 `Ωc`), `ref-310-two-rotations-all.txt` (dim-15 transfer matrices, states; the HW
sectors of [3,1,0]⊗[3,1,0] as in `L2-310-untwist-CG.md` §2), `ref-general-all.txt` cell 28 (Gauss factors L₀D₀U₀ of
P_phys⁻¹, `rhoLow`, `rhoUp`).

## 0. The family and the regime

Diagonal frame (where the eigenstates, HW kets and the CG extraction live): Ω_w = e₁ᵀρ(n(z,w)), n(z,w) =
exp(r₁E₁₂)exp(r₂E₁₃)exp(r₃E₂₃) with the r_k above (depends on w₁, w₂ only); Ω_w⊗Ω_w is the exact eigen-covector of
t_a with bra twist w. Companion frame: Ωc[w] = Ω_w·P_phys (check: Ωc[{w₁,w₂}] = e₁ᵀρ(n(z,w))·P_phys exactly at rational
z, w — this ties Stage A's companion-frame vacuum to the diagonal-frame one). The family: Φ = Ω_w ⊗ Ω_{w′} (four
parameters w₁, w₂, w₁′, w₂′; w₃, w₃′ irrelevant). HW kets h_ν of Λ⊗Λ are fixed by ρ(n)⊗ρ(n) for n ∈ N₊, so the EXACT
undo is ⟨Ω_w⊗Ω_{w′}|h⟩ = ⟨e₁ᵀρ(n(w)n(w′)⁻¹) ⊗ e₁ᵀ|h⟩ (no twist-matrix Gauss factors needed in this frame: the vacuum
Ω_{w′} on site 2 is itself e₁ᵀρ(n(w′)), unlike the companion-frame x₀^{(1)} = e₁ᵀρ(U₀)). The undone unipotent
n′ = n(w)n(w′)⁻¹ has Gauss parameters (a′, b′, c′) rational in (z, w, w′), and the commutator direction b′ − a′c′ is
the one that degenerated in the old det₆×det₂ family (y pinned). Here all three are expected O(1) in the joint regime.

**Regime:** z = 1+εφ, φ = (1/3, −2/7, 5/11) (control φ′ = (−2/5, 3/8, 7/11)); w = 1+εω, ω = (ω₁, ω₂, 5/13) (w₃
irrelevant but keep it ≠ others); w′ = 1+εω′. Grid: ω₁ ∈ {3/4, −5/3}, ω₂ ∈ {−1/5, 9/7, 5/3}, ω₁′ ∈ {−12/5, 7/2},
ω₂′ ∈ {1/8, −3/4} (24 points; all distinct from the φ_k and from each other); ε ∈ {1, 2, 4}·10⁻⁶; 200 digits on the
chain, `x0Of`-type factorised inverses (never `Inverse` at small ε).

## 1. Tasks (each its own cell(s) + markdown result cell with a `Grid`; save after each)

**Task 0 — setup.** Generators, metric G, HW kets of Λ⊗Λ and the G-orthogonal projectors per HW sector (adjoint: 6 HW
states / 5 irreps; [3,1,0]: 10 / 8, from `L2-310-untwist-CG.md` §2), exact sparse transfer matrices, `branchesAt[z]`
(HW-limit branches normalised ⟨x₀|Ψ⟩ = 1 in the diagonal frame as in the reference), `oneSiteVac`, `Ωd[z,w]` :=
e₁ᵀρ(n(z,w)) on V (ρ = `MatrixExp` of the Gauss factors), the companion-frame tie Ωc[w] = Ωd[z,w]·P_phys(z) (exact,
rational point), and the exact-undo check ⟨Ωd[w]⊗Ωd[w′]|h_ν⟩ = ⟨e₁ᵀρ(n(w)n(w′)⁻¹)⊗e₁ᵀ|h_ν⟩ on all HW kets at rational
(z, w, w′) (exact zero). Also verify Ωd[w]⊗Ωd[w] is the exact left eigen-covector of t_a(u) with twist diag(w) (one u).

**Task 1 — the undone unipotent, symbolic in ε (3×3, exact).** n′(ε) = n(z(ε),w(ε))·n(z(ε),w′(ε))⁻¹; print a′, c′,
b′ − a′c′ as rational functions of ε and their leading orders; the ε = 0 values in closed form in (φ, ω, ω′) (`Factor`;
expectation: all three O(1), and the map (ω₁, ω₂, ω₁′, ω₂′) ↦ (a′₀, c′₀, b′₀ − a′₀c′₀) has rank 3 at a generic point).
Controls: (i) w′ = w exactly (n′ = 1); (ii) the old regime's degeneration: w′ with w₂′ pinned at an O(1) value ≠ 1
(w₂′ = y* = w₂(t−1)/(t−2)-type) — expect b′ − a′c′ to drop an order; (iii) common w₁′ = w₁: which direction is lost
(the localising-span result says vacua of common w₁ differ only in E₂₃: b′ = a′c′ exactly? check).

**Task 2 — ε-structure of the undone covector.** φ_u(ε) = e₁ᵀρ(n′(ε)) on V: per weight of V the leading order and the
rank of the ε⁰ coefficient vectors over the grid (and φ′). Expectation: EVERY weight at ε⁰ with rank = multiplicity
(adjoint: 7 weights, rank 2 on (1,1,1); [3,1,0]: 12 weights, rank 2 on (2,1,1), (1,2,1), (1,1,2)). Control: w′ fixed
(ε-independent, e.g. w′ = (7/4, −2/5, 9/7)) → the depth table (which weights drop to which order).

**Task 3 — chain side: signal vs contamination.** At the grid × ε: R_C = ⟨Φ|Ψ_B⟩, S = ⟨Φ|P_HW Ψ_B⟩ (= ⟨φ_u⊗e₁ᵀ|P_HWΨ_B⟩
by the exact undo — check numerically to 10⁻¹⁵⁰), K = R_C − S; ε-exponents of |S|, |K|, max |K/S| per branch over the
grid. **Expectation to confirm or refute: K/S = O(ε) for EVERY HW branch** (adjoint 6 incl. singlet and both octets;
[3,1,0] all 10 incl. [4,4,0], [4,3,1]×2, [4,2,2], [3,3,2] which had ε-constant κ in the det₆×det₂ family). If a branch
fails: `showContam`-type breakdown (which site-1 weights carry K at which order), no fitting.

**Task 4 — CG at leading order from the overlap alone.** 1-dim slices: C := R_C/[φ_u(0)]_slice, grid spread and error vs
the chain slice at two ε's (∝ ε), Richardson limits vs chain (adjoint sanity against 46.67453891 ([4,1,1]),
10.92853892 ([3,3,0]), −5118.538314 (singlet) in the `L2-210-untwist-CG` normalisation — compute the chain slice here
and compare to it; quote references only as a cross-check). 2-dim slices: 2×2 systems from two grid points (adjoint
octets: mixing R vs −1.694, 1.096; [3,1,0]: [5,2,1]×2 and [4,3,1]×2). Conclusion per rep: how many of the HW CG data
are read at leading order (expect all).

**Task 5 — the det₆-type data of the family at small ε (preparation for the Q side; cheap).** For one grid point and
ε = 10⁻⁶: the vacuum direct Q's of twist w and of twist w′ (degrees as in the reference), the dual Q's of the branches,
and the four numbers ⟨Ω_w⊗Ω_w|Ψ⟩, ⟨Ω_{w′}⊗Ω_{w′}|Ψ⟩ (= det₆/Ψ_Q, check), ⟨Ω_w⊗Ω_{w′}|Ψ⟩, ⟨Ω_{w′}⊗Ω_w|Ψ⟩ on the HW-limit
branches, with their ε-exponents per branch. Store them in a markdown table for the Stage C/Q-side task.

**Task 6 — markdown summary** `### two vacua Ω_w⊗Ω_{w′}: untwist result (2026-10-11)` per rep: the closed-form ε = 0
Gauss parameters and the rank-3 statement, the per-weight order/rank table (+ controls), the K/S table with
ε-exponents, the CG errors and Richardson limits, the Task 5 table, surprises. Save. Do not write HANDOFF.

## 2. Rules and gotchas

As in the previous stage specs: exact where possible, 200 digits on the chain; time one state/one point before loops
([3,1,0]: 56 HW-sector states per ε ≈ 3 s earlier); nothing over ~3 min without a timed smaller run; markdown result
cells under code cells; save after each task; no other notebooks, no `paper/*.tex`, no git, no HANDOFF.
- Shared kernel K1; first cell wipes other state. Never `*)` in comments. No `D`,`E`,`N`,`I`,`C`,`K` symbols.
- `x0Of` with the factorised inverse; `Do[DownValues[f] = …]` purges nothing; `Max[{}]` = −∞; `SetPrecision` before
  `Norm` of near-zero vectors; exact zero tests with `===` only on rationals.
- The reference untwist notebooks used the diagonal frame with x₀ = (e₁ᵀP_phys⁻¹)^{⊗2}; here x₀ enters only through the
  branch normalisation ⟨x₀|Ψ⟩ = 1 — keep it identical to the reference so CG numbers are comparable.
