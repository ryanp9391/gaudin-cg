# Spec (Job C): what is κ(s,t)? — the adjoint (∞,∞) family as a global rotation of the vacuum times a site-2 principal torus; localisation and vacuum-Q tests (2026-10-11)

Delegated to a wolfbook-builder subagent (sonnet); planned and reviewed by the main session. **New notebook:**
`Paul/Mathematica/new-conventions/L2-210-principal-torus.wb` (adjoint; build/run as the previous stages). Stage C found
that the (∞,∞) two-vacuum family Φ(a,a′) = gaussCov[{a,0,a}]⊗gaussCov[{a′,0,a′}] has block tensors
A^{(cc′)} = a^{m+m′+e}a′^{n+n′+e′}κ^{(cc′)}(s), s = a′/a, with κ polynomial in s (degree ≤ 4) and rational in t (poles
t ∈ {0,±1}); the 10-product formula's measures are F_m = a^m f_m(s), G_m = a^{m+e}g_m(s). Paul's question: is κ(s,t) a
known object?

**Main-session observation to verify first.** With D_a := diag(1,a,a²) and g₀ := exp(E₁₂)exp(E₂₃):
D_a⁻¹ g₀ D_a has Gauss parameters (α,β,γ) → (a, 0, a) (torus conjugation rescales (α,β,γ) by (d₂/d₁, d₃/d₁, d₃/d₂)), so
e₁ᵀρ(exp(aE₁₂)exp(aE₂₃)) = e₁ᵀρ(D_a⁻¹)ρ(g₀)ρ(D_a) ∝ e₁ᵀρ(g₀)ρ(D_a). Hence
  Φ(a,a′) ∝ (e₁ᵀρ(g₀) ⊗ e₁ᵀρ(g₀)) · (ρ(D_a) ⊗ ρ(D_{a′})) = Ω₀ · (ρ(D_a)⊗ρ(D_a)) · (1 ⊗ ρ(D_s)),  Ω₀ := x₀(ρ(g₀)⊗ρ(g₀)), s = a′/a,
i.e. a FIXED global rotation Ω₀ of the SoV vacuum, a global principal torus (which should be responsible for the
a^{m+m′+…} grading of the SoV coefficients), and a site-2 principal torus ρ(D_s) = diag(1,s,s²) on V (eigenvalues
s^{−2..2} on the adjoint) carrying all the s-dependence. So κ(s,t) = the s-expansion of the site-2 torus insertion on
the globally rotated vacuum. Read first: HANDOFF top section, Stage C notebook cells 26–31, 59–64 (kernels and their
t-laws), `L2-210-localising-span.md` and the localising-span notebook (`L2-210-localising-span.wb`: the fused operators
T_{ab}(v) at the site points, the spans D₁ (7), D₂ (24), D₃ (48) on the base ⟨Ω_w⊗Ω_w|, the exact membership tests),
`Two-rotations-det-overlap.md` §5.2 (localisation on the twisted vacuum, the ω₂ direction).

## 1. Tasks (adjoint; each its own cells + markdown result; save after each)

**Task 0 — the factorisation (exact).** Stage B/C setup (SoV basis, K tensor, blocks). Verify on V: e₁ᵀρ(gaussCov-type g
with (a,0,a)) = λ(a)·e₁ᵀρ(g₀)ρ(D_a) exactly (rational a; find λ(a)); Φ(a,a′) = λλ′·Ω₀(ρ(D_a)⊗ρ(D_a))(1⊗ρ(D_s)) on V⊗V.
Then the grading: is every SoV covector ⟨x_m| an eigen-covector of the global principal torus ρ(D_a)⊗ρ(D_a) (exact test on
the 64 covectors; eigenvalue a^{γ(m)}; if yes, γ(m) vs the labels — expect γ = −(m+m′+n+n′) + const, matching Stage C's
exponents (e,e′) = (7,5),(7,4),(6,5),(6,4) per block; if not, say what the torus does to the SoV basis). Consequence: the
whole 4-parameter family's a-, a′-dependence is a grading; κ(s) is the only content.

**Task 1 — κ(s,t) as the site-2 torus on Ω₀.** Expand ρ(D_s) = Σ_{k=−2}^{2} s^k Π_k (Π_k the projectors onto the principal
grades of the adjoint); the covectors Φ_k := Ω₀(1⊗Π_k) (5 of them); their SoV coefficient vectors φ^{(k)} and blocks;
check Σ_k s^k A(Φ_k) reproduces κ^{(cc′)}(s) up to the grading factor (exact in s at t = 62/35). Ranks of the blocks of
each Φ_k and of partial sums (which k's carry the rank-3 parts). Which Φ_k read which irreps (⟨Φ_k|h_ν⟩ on the 6 HW
kets, exact: Ω₀ itself reads only the top; the k ≠ 0 pieces are what makes the family CG-complete — tabulate).

**Task 2 — localisation test (the "known object" test).** Base covector Ω₀ (companion frame; convert to the frame of the
localising-span notebook if its operators are built there — or rebuild the fused T_{ab}(v) in the companion frame, they
are frame-covariant: T_c = (P⁻¹⊗P⁻¹)… check one eigenvalue). Spans: D₁ = {Ω₀ T_{ab}(v₀)} over the six non-scalar
coefficient operators (t₁ at 2 points, t₂ at 4; see the localising-span write-up for the exact list and the site point
v₀ = θ₁ − 2h; ALSO the site-2 point θ₂ − 2h since the insertion is on site 2 — add the ladder-2 versions), D₂ = products
of two, D₃ = products of three (dimensions as in the write-up: 7, 24, 48 on the vacuum base; recompute on Ω₀). Membership
of each Φ_k, k = −2..2 (and of Φ(s) at three rational s) in D₁, D₂, D₃: exact linear algebra (rank of the augmented
matrix). If some Φ_k ∈ D_ℓ, give the expansion coefficients (rational functions of t? repeat at the second (h,θ) point
(1/3,(1/5,−4/9))). If no Φ_k is in D₃: the rational version — is Φ(s)·Q(T) ∈ D₂ for some polynomial Q(T) of the same
operators (the (2,2)-rational test of the localising-span notebook §3) — try the (2,2) and (3,1) levels; report
solution dimensions. State the conclusion: κ(s,t) is / is not an eigenvalue-polynomial (transfer-matrix) insertion on
the rotated vacuum through level 3 / rational level (2,2).

**Task 3 — vacuum-Q test.** One-site vacuum direct Q's of twist (y₁, y₂): q_i^{vac}(x) from the one-site Baxter equation
(rotation-sov cell 25 pattern), symbolic in (y₁, y₂) if feasible (degrees 0,1,2 — `Solve` with symbolic twist, h = 1/2,
θ₁), then the limit y₂ → ∞ (rescale) and y = 1/a, 1/a′. Test whether the kernel functions f_m(s), g_m(s) (Stage C cell 62,
all 10 products) are ratios/products of such vacuum-Q values at the ladder points, i.e. fit each f_m(s) and g_m(s) as
Π_i Q_i^{vac,1/a′}(θ+m_ih)^{±1} × const over small exponent patterns — a finite search; also the simpler test: are the
roots in s of the κ^{(cc′)}(s) polynomials (over Q(t)) equal to ratios of roots of the one-site vacuum q's (compute both
root sets at t = 62/35 and 29/15 exactly and compare). Report matches or none.

**Task 4 — markdown summary** `### the (∞,∞) family as rotated vacuum × site-2 principal torus — result (2026-10-11)`:
the factorisation, the grading statement, the Φ_k table (ranks, irreps read), the localisation membership results with
dimensions, the vacuum-Q test, and the answer to "is κ(s,t) a known object" with the evidence. Save; no HANDOFF.

## 2. Rules and gotchas

As before (K1 shared; first cell wipes; no `*)` in comments; no `D`,`E`,`N`,`I`,`C`,`K` symbols — call the torus `Dg`;
time before loops; `TimeConstrained`; save after each task; no other notebooks, no `paper/*.tex`, no git, no HANDOFF).
The localising-span notebook works in the diagonal frame with x₀ = (e₁ᵀP⁻¹)^{⊗2}; converting a companion-frame covector
Φ_c to the diagonal frame is Φ_d = Φ_c·(P⁻¹⊗P⁻¹) (check against Stage D cell 11's frame tie). Exact membership: use
rationalised SoV coefficients or exact covectors (Ω₀ is exact rational; the T's are exact), so ranks are exact.
