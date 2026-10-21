# Spec W3: sl(3) Gaudin chain, L = 2, Λ = [3,2,0] — the full chain-side check of every rep-generic rule, both two-vacua families, the cleaning question (2026-10-16)

Delegated to a wolfbook-builder subagent (sonnet); planned by the main session. **New notebook:**
`Paul/Mathematica/Gaudin/Experiments/Gaudin-L2-chain-320.wb` (self-contained; `Λ = {3, 2, 0}` as the switch; first code cell
`Quiet[ClearAll["Global`*"]]`). Read first: `../Gaudin-sl3-round3-plan.md` §0 (the prediction table is what you check), §2, §6; then the
round-2 plan `../Gaudin-sl3-round2-plan.md` §1 (the verified facts, incl. the "Added by W2/P3/D2/D1" blocks) and the round-1 plan §1.
**Templates:** `Gaudin-L2-chain.wb` (round 1, 48 cells) and `Gaudin-L2-chain-2.wb` (round 2, 80 cells), same directory — read their JSON
with a `python -X utf8` script and copy the code cells you need (setup, Lax, cdet, Laurent coefficients, vacuum/frames/κ^G, x0G, hwKets/PHWc,
eigenstates/branches, undo/K/S/CG, localising scan incl. the 𝓛 computation and sparse supports, closed-form cells, universality cell). Their
last cells are the summaries of what is verified. [3,2,0] has dim 15 and 225 states like [3,1,0], which both notebooks already ran, so the
timings are known (SparseArray, dense per sector at 200 digits; time one product first).

## 0. Points and conventions

Λ = {3,2,0}; θ = {2/7, −3/5}; ζ: P = {3/2, 5/4, 2/3}, P′ = {2, 5/3, −1/2}; ω = {5/2, 7/3, 3/5}, ω′ = {7/4, −2/5, 9/7}; (ω₁, ω₁′) = (5/2, 7/4) for
the (∞,∞) family, second pair (−2/5, 3); regime ζ = εφ, φ = {1/3, −2/7, 5/11}, ω̂ ∈ {3/4, −5/3}×{−1/5, 9/7, 5/3} and ω̂′ ∈ {−12/5, 7/2}×{1/8, −3/4}
for the full family (reduce to the 4-point sub-grid if the loop exceeds 10 min), ω̂₁ ∈ {3/4, −5/3}, ω̂₁′ ∈ {−12/5, 7/2} for the (∞,∞) family;
ε ∈ {1,2,4}·10⁻⁶; pr = 200. Everything exact where possible.

## 1. Tasks (code cell + markdown result cell with a `Grid` each; `checks` association; save after each)

**A1 — Hamiltonians.** cdet as in round 1 (copy); checks: [S_a(u), S_b(v)] = 0 (≥ 4 pairs), S₁ = e₁(ζ) + 5/(u−θ₁) + 5/(u−θ₂), weight conservation;
Laurent fits: the S₂ double-pole scalar (predicted 8 = λ₂(λ₁+1)), the S₃ triple pole (predicted absent), the six non-scalar data; H₁'s explicit
form with the round-1 sign; the D₂ relations: test S3c[s,1] = c_Λ·S3c[s,2](H_s − S3c[s,2]) for c_Λ ∈ {1/4, 1/3, 1/5, 1/6, 1/8, 1/(λ₁+1), 1/\|Λ\|}
and, if none, solve for the exact quadratic relation(s) in D₂ (rank of D₂; predicted ≤ 26).

**A2 — vacuum and frames.** ΩD[ζ,ω] (same r_k; ρ via the exponential series with E^{λ₁−λ₃+1} = E⁴ = 0), eigencovector of S₂^{(g′)}, S₃^{(g′)} with
g′ = PauxG[ζ].gc[ω].PauxG[ζ]⁻¹ (and the control that the transposes fail); the vacuum Laurent data vs the prediction H_s^{vac} = 2ω₁ + 3ω₂ + 5ω₃ ±
k/θ₁₂ (find k; predicted 12), a_s^{vac} = 5ω₃, b_s^{vac} = ω₃(H_s − 5ω₃); companion vacuum Ωc^G[ω] and its Gauss parameters (predicted the same
(1/ω₁, −1/(ω₁ω₂), (ω₁+ω₂)/(ω₁ω₂)) — rep-independent); frame tie and κ^G = Δ₁^{λ₁−λ₂}Δ₂^{λ₂−λ₃}Δ₃^{λ₃} = Δ₁Δ₂² (verify); x0G nonvanishing on all states.

**A3 — eigenstates, HW sectors, branches, universality.** All 225 states at P (eigenvalue table, 30 digits, for P4); the HW sectors of V⊗V with
multiplicities (predicted the dual images of [3,1,0]'s: 10 HW states, multiplicity 2 in (5,4,1) and (5,3,2)); exact HW kets and G-projectors;
HW-limit branches at the regime point (raising residuals and gaps); for each multiplicity-2 sector the matrix of S3c[1,2] on the exact HW
space at θ₁₂ = 31/35 and at θ = (1/5, −4/9) and (3/7, −1/4): is it (1/θ₁₂) × a θ-independent matrix with no scalar part (W2's finding for
[2,1,0] and [3,1,0])? Its eigenvectors = the untwisted copies; the coordinate-ratio quadratic for each sector; compare with the [3,1,0]
sectors (5,2,1) ↔ (5,4,1) and (4,3,1) ↔ (5,3,2) under duality (do the quadratics map into each other under the dual GT relabelling?).

**A4 — the two-vacua untwist, BOTH families (the cleaning question).** (a) Full family Φ = ΩD[ζ,ω]⊗ΩD[ζ,ω′] on the grid: exact undo, closed forms
of (a′, c′, g′), per-weight ranks (CG-complete?), K/S for ALL 10 branches at three ε (expect O(ε), doubling), CG at leading order vs the exact
G-norm-1 HW slices after Richardson, incl. the 2-dim slices of (5,4,1) and (5,3,2) and the mixing vs the exact copies of A3. **Report the
second singular value of every 2×2 extraction system** (the XXX "weakly conditioned" diagnostic) and whether any branch would have needed a
cleaning step (i.e. whether K/S ≠ O(ε) anywhere, or whether the extracted CG converge to a wrong limit). (b) The (∞,∞) family Φ∞ =
ΩDinf[ζ,ω₁]⊗ΩDinf[ζ,ω₁′] (W2's r_k^∞): the same list; controls ω₁′ = ω₁ and ω₁′ fixed (W2 found the fixed control FAILS at ω₂ = ∞ for the
singlet-type branches — does it fail for [3,2,0] too, and on which branches?).

**A5 — localisers.** D₁ rank (7), D₂ rank and its exact relations (A1), the exact localising subspace 𝓛 ⊂ D₂ (dim; predicted 4), 𝓛∩D₁ (predicted
{id, H₁, S3c[1,2]}), the fourth localiser: does X₄ = H₁² − ½H₁S3c[1,2] − (2/θ₁₂)W − S3c[1,1] still localise (mod D₁∩𝓛)? If not, fit X₄ =
H₁² + αH₁S3c[1,2] + βS3c[1,2]² + (γ/θ₁₂)W + δS3c[1,1] with exact α, β, γ, δ (and report their Λ-dependence against [3,1,0]'s (−½, 0, −2, −1));
sparsest supports; the fourth covector v₄ decomposed in Ω₁·(one-site operators and pair products); spectators f(H₁) (predicted 5Σζ − ζ·Λ +
12/θ₁₂), f(S3c[1,2]) (predicted 5ζ₃), f₄; site 2 mirror. Weights reached by the 3-dim and 4-dim L₁.

**A6 — closed forms.** c_H in ⟨Ω⊗Ω|H₁ = [Ω₁(c_H − C_ζ)]⊗Ω₁ (predicted 5Σζ + 12/θ₁₂, already verified by W2 for [3,2,0] — a one-line regression);
the spectator and vacuum closed forms of A2/A5 confirmed at P, P′ and two more (ζ, ω, θ) points each.

**A7 — reference tables for P4.** At P with (ω₁, ω₁′) = (5/2, 7/4): for all 225 states (G-normalised, round-1 sign convention) the eigenvalue
pair (S₂(u₀), S₃(u₁)), ⟨x0G|Ψ⟩, ⟨Φ∞|Ψ⟩, Qnum = ⟨Φ∞|Ψ⟩/⟨x0G|Ψ⟩, 30 digits; the same for the 10 HW branches at the regime point (ω̂₁, ω̂₁′) =
(3/4, −12/5), ε = 10⁻⁶.

**A8 — summary cell** `### sl(3) Gaudin L = 2, Λ = [3,2,0] (2026-10-16)`: the prediction table of the plan with a verdict per row (✔ / ✘ + the
actual value), the cleaning verdict (K/S table, 2×2 conditioning), the fourth localiser, surprises, not done.

## 2. Rules and gotchas

New-notebook workaround (JSON built by a Python script, `code <path>` once, every code cell's exact saved source evaluated in K1 via
`wolfbook_evaluateExpression`, printed `Grid` stored as a markdown cell under the code cell, re-read before every rewrite, no native outputs —
say so). K1 is yours alone (kernel `k-d8a6a147`, client `VSCode[gaudin-cg]`); first cell wipes Global`. Semicolons; never `*)` in a comment;
no symbols D, E, N, I, C, K; memoised definitions cache bad results; `===` on rationals, `Precision` on numerics; Laurent fits must not sample
at u = θ_s or at u₀ = 3/7 if θ₁ = 3/7 (round-2 pitfall); time one 225×225 product and one 200-digit sector before loops; nothing over ~3 min
without a timed smaller run. Only this notebook; no other notebooks, `paper/*.tex`, HANDOFF.md, CLAUDE.md, `Paul/Python/`, or git. Numbers
with precision, not adjectives; if an expectation is refuted and two attempts do not fix it, stop that item, keep what works, report.
