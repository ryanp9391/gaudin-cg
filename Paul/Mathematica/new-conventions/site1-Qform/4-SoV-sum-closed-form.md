# Direction 4: accept the 12-state SoV sum and put it in closed form

## What exists (exact, `../L2-rotation-sov.wb` §3–4)

- **SoV basis.** The untwisted L=2 B has a simple spectrum, and its 64 SoV states form an 8×8 grid.
- **Wavefunction.** ⟨x|Ψ⟩/⟨x₀|Ψ⟩ = N_x ∏_s det₂[v(x_{s21}), v(x_{s1})] v_k(x_{s22}) / vac, exact.
- **The site-1 rotation:**

  ⟨(e₁ᵀρ(g))⊗e₁ᵀ|Ψ⟩/⟨x₀|Ψ⟩ = Σ_ν c^{L=1}_ν(g) W_ν(B),  W_ν = Σ_x K_{νx} ψ⁽¹⁾_{x₁}(B) λ₂(B)^{[x₂≠vac]}

  - c^{L=1}(g) = c(a,b) is the L=1 bilinear classification formula.
  - K is an 8 × 12 rational matrix (table in `../Rotations-SoV-Q.pdf` §8). Its entries contain 31s and 35s from t = θ₁₂/h = 62/35.

## Goal

A complete, explicit Q-formula with no numerics:
- K_{νx} as rational functions of t (and possibly h, but h should scale out);
- the N_x normalizations in closed form;
- then a resummed form: Σ_ν c_ν(w, y₁) W_ν as an explicit expression in the L=2 Q's at ladder-1 points, plus the single site-2 quantity λ₂(B).

## Steps

1. **Second and third numeric points.**
   - Recompute K at P2 = (h, θ) = (1/3, (1/5, −4/9)), where t = 29/15, and at a third point P3.
   - K should not depend on z (check: recompute at a second z).
   - Fit each entry as a rational function of t of low degree. Entries like −35/124 suggest factors such as t/(t+…) (at P1, t = 62/35, and 124 = 2·62).
2. **Symbolic route (preferable).**
   - The 4 non-eigen covectors x_ν⊗e₁ᵀ (ν = {2|11}, {2|10}, {2|00}, {1|00}) need the L=2 B-eigenvectors with site-2 pattern {2|20}.
   - B is a polynomial matrix in u with explicit θ's, so compute these eigenvectors symbolically in θ₁₂, h: solve (c·B(u) − β(u)c) = 0 for c in the small subspace spanned by x_μ⊗e₁, x_μ⊗e₂ and their neighbours.
   - This gives K exactly, together with the N_x.
3. **Interpretation of λ₂.**
   - λ₂(B) = v(θ₂)/v(θ₂−h) is the proportionality constant of the degenerate dual columns at site 2.
   - Via the Baxter equation it should equal a T_λ(θ₂ + …) eigenvalue ratio.
   - Identify it, so the formula reads "site-1 SoV sum × (1 + transfer-matrix correction at site 2)".
4. **Resummation.**
   - At L=1, the sum Σ c_μ(w, y) ψ_μ resums to det₃ × det₂ / Ψ_Q.
   - At L=2, try to resum the 8 vacuum-site-2 terms into a ladder-1-only bracket determinant built from the L=2 Q's, and the 4 λ₂ terms into a second one. Directions 2 and 3 provide candidate determinant shapes.

## Success criteria

- K(t) and N_x verified exactly at P1, P2 and P3.
- The full formula reproduces all 64 overlaps at a fresh numeric point with no fitting.

## Cost and risks

- **Cost.** Step 1 is about 2 × (setup ~1 min + the K computation, seconds). Step 2 is symbolic eigenvectors of 64×64 polynomial matrices restricted to small subspaces, which should be fine.
- **Risk.** The rational t-dependence may need more points than expected, if the degrees are high.
