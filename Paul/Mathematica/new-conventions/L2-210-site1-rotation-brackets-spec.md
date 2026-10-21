# Spec: the generic site-1 rotation on the SoV vacuum as ladder-1 Q-brackets — exact loci in rotation space (2026-10-10)

Delegated to a wolfbook-builder subagent (sonnet); planned and reviewed by the main session. **New notebook:**
`Paul/Mathematica/new-conventions/L2-210-site1-rotation-brackets.wb` (self-contained; first code cell
`Quiet[ClearAll["Global`*"]]`). Reference cell sources: scratchpad file
`ref-L2-rotation-sov-cells1-22-59-80.txt` (named in the task message) = `L2-rotation-sov.wb` cells 1–22 (setup,
companion frame, SoV basis of the untwisted B, 12-state sum, K map) and 59–80 (2-plane variables, the exact
polynomial, the equivariant map X, Kvac/K20, the symbolic rotation c(a,b)). Read it before writing a cell; copy
what you need. Background write-up: `Two-rotations-det-overlap.md` §1–§3 (Part A).

## 0. Paul's idea

Forget the w-vacuum Ω_w. Start from the SoV vacuum itself, x₀ = e₁ᵀ⊗e₁ᵀ in the companion frame, put the most
generic rotation on site 1, Φ(g) = e₁ᵀρ(g) ⊗ e₁ᵀ (three parameters), and ask for which g the overlap
⟨Φ(g)|Ψ_B⟩/⟨x₀|Ψ_B⟩ can be written purely in terms of ladder-1 brackets of Q with some measures — "integrals of Q
with measures, like the entries of the ⟨Ω|Ψ⟩ determinant" — without insisting on a determinant. First pull out
the θ₂-ladder dependence (it enters only through λ₂ = Q^k(θ₂)/Q^k(θ₂−h)) and treat the θ₁ part alone.

**What is already exact** (`L2-rotation-sov.wb` §8, Part A): for any site-1 covector c = Σ c_ν x_ν (x_ν the 8
constant L=1 SoV covectors),

    ⟨c⊗e₁ᵀ|Ψ_B⟩/⟨x₀|Ψ_B⟩ = Σ_{x∈trip} φ⁰_x ψ_x + λ₂ Σ_x φ¹_x ψ_x,   φ⁰ = c·Kvac, φ¹ = c·K20,

with the 8 site-1 SoV wavefunctions ψ = {1, β, −α, λ₁, λ₁β, λ₁²β, −λ₁α, −λ₁²α} (trip order
(−3,−2,−1), (−3,−1,−1), (−2,−1,−1), (−3,−2,0), (−3,−1,0), (−3,0,0), (−2,−1,0), (−2,0,0)), where
v(θ₁−h) = αv(θ₁−3h) + βv(θ₁−2h), λ₁ = v_k(θ₁)/v_k(θ₁−h), v = (Q¹,Q²,Q³). Kvac has rank 8, K20 rank 4.
For a rotation, c = c(a,b) is bilinear in X = a bᵀ (b·a = 0), cell 80.

**Main-session observation (to be verified in Task 1).** With u := v(θ₁−3h), w := v(θ₁−2h) and the ladder-1
normalisation Ψ₁ := det₂[w,u]_{rows 2,3}·v₃(θ₁−h) (the ladder-1 factor of Ψ_Q(x₀)), every ψ_x·Ψ₁ is a
monomial M_x = det₂[v(n),v(n′)]_{rows 2,3} · v₃(r) in the ladder-1 Q-values:

| trip x | (n, n′; r) | ψ_x |
|---|---|---|
| (−3,−2,−1) | (−3,−2; −1) | 1 |
| (−3,−1,−1) | (−3,−1; −1) | β |
| (−2,−1,−1) | (−2,−1; −1) | −α |
| (−3,−2,0) | (−3,−2; 0) | λ₁ |
| (−3,−1,0) | (−3,−1; 0) | λ₁β |
| (−3,0,0) | (−3,0; 0) | λ₁²β |
| (−2,−1,0) | (−2,−1; 0) | −λ₁α |
| (−2,0,0) | (−2,0; 0) | −λ₁²α |

(up to the sign convention of det₂ — fix it so the table holds). So N⁰ := Ψ₁·Σφ⁰_xψ_x = Σ_x φ⁰_x M_x and
N¹ := Ψ₁·Σφ¹_xψ_x = Σ_x φ¹_x M_x are explicit polynomials of multidegree (0,1,2) in (Q¹,Q²,Q³) on the four
ladder-1 points, and the question "brackets with measures" becomes a question about the tensor
C_{(nn′),r} (6 minors × 4 points, 8 populated slots for N⁰, 4 for N¹): a ladder-1 bracket is
⟨Q^jF⟩₁ = Σ_n F(n)Q^j(θ₁−nh) (F = measure × bra data, a free function on the 4 points), a 2×2 bracket
determinant det[⟨Q²F⟩,⟨Q³F⟩;⟨Q²G⟩,⟨Q³G⟩] = Σ_{n<n′}(F∧G)_{nn′}det₂[v(n),v(n′)], and products with a third
bracket ⟨Q³H⟩ give rank-1 tensors (F∧G)⊗H. The structured forms, from most to least restrictive:

  (R1) single product: N = det₂[F,G]·⟨Q³H⟩ ⇔ C = A⊗H with A decomposable;
  (R3) one 3×3 bracket determinant of the L=1 type, det[⟨Q²F_i⟩, ⟨Q³F_i⟩, ⟨Q³H_i⟩]_{i=1..3} (shared bra
       functions, Laplace-expanded along the third column: Σ_i ±(F_j∧F_k)⊗H_i);
  (R2δ) sum of two products with point measures for the Q³ factor: N = det₂[F₋₁,G₋₁]·Q³(θ₁−h) +
       det₂[F₀,G₀]·Q³(θ₁) ⇔ both columns C_{·,−1}, C_{·,0} decomposable (C_{·,−3} = C_{·,−2} = 0 already);
  (R2) sum of two products with free H's — vacuous (every pencil of 2-forms on C⁴ contains two decomposable
       ones), do NOT report it as a result.
Counting on N⁰: C_{·,−1} lives on the points {−3,−2,−1} (always decomposable); C_{·,0} on four points
(decomposable iff one Plücker quadric, φ⁰₆φ⁰₇ − φ⁰₅φ⁰₈-type with your sign conventions). So (R2δ) should cut a
quadric hypersurface in rotation space (a 2-parameter family of rotations), (R1) needs C_{·,0} ∝ C_{·,−1}
(φ⁰₆ = φ⁰₈ = 0 and proportionality: expected only degenerate rotations), (R3) needs in addition that the two
decomposable columns share a factor (C_{·,−1}∧C_{·,0} = 0 in Λ⁴C⁴): a second quadric. On N¹ both columns live
on {−3,−2,−1} so (R2δ) is automatic and (R3) should be automatic too. The interesting composite is the full
λ₂-dependent object N(λ₂) = N⁰ + λ₂N¹: column r is the pencil A_r(λ₂) = C⁰_{·,r} + λ₂C¹_{·,r}; it is
decomposable for ALL λ₂ (so the θ₂-dependence can be absorbed into the measures F_r(λ₂), G_r(λ₂)) iff
C⁰_{·,r}, C¹_{·,r} are decomposable and C⁰_{·,r}∧C¹_{·,r} = 0 — automatic for r = −1, two conditions for
r = 0. **These exact loci in (a,b) are the deliverable.** Everything is exact rational/symbolic; no numerics
except the 64-state verifications.

## 1. Conventions

Companion frame, gl(3) adjoint, L = 2, new conventions, h = 1/2, θ = (2/7, −3/5), ket twist z = (3/2, 5/4, 2/3)
(and once z′ = (3/2, 5/4, 3/7) to confirm the loci are z-independent — Kvac, K20 are rational constants
depending on (h, θ) only; verify at z′). Precision `prec` = 60 for the eigenstates (as the reference cells),
everything else exact. Reuse: `ketC` (64 companion-frame eigenkets), `x0C`, `colV2[B,s,m]` (dual Q values
v(θ_s+mh)), the SoV basis `vecB`, `psi1`, `λ2`, `trip`, `xL1`, `Kvac`, `K20`, `Xc`, `rank1dev`, `rotFromAB`,
`cAB`, `Φvac`, `Φ20` (cells 16–21, 59–64, 80). If `rotFromAB` is defined in a cell not in the reference file,
rebuild it: c(a,b) is the unique (linear in a bᵀ) covector with X(c) = a bᵀ; alternatively parametrise the
rotation directly as c = (e₁ᵀexp(aE₁₂)exp(bE₁₃)exp(cE₂₃))·xL1⁻¹ with symbols (a, b, c) — use whichever makes
the exact solves shortest, and state which (the (a,b)-bilinear form with b·a = 0 is projective: 4 parameters
including the scale; the Gauss form is affine 3 + scale).

## 2. Tasks (each its own cell(s) with a `Grid` of checks; `checks` an Association; time one case before loops)

**Task 0 — setup and regression.** Rebuild cells 1–22 and 59–64, 80 (skip the unrelated experiments); checks:
the 12-state sum reproduces the exact overlap on all 64 states for a generic rotation and a random c (10⁻⁴⁷);
Kvac rank 8, K20 rank 4; `Φvac`, `Φ20` bilinear in (a,b); X(c(a,b)) = a bᵀ. Timing.

**Task 1 — the monomial table.** Verify on all 64 states: ψ_x·Ψ₁ = M_x for the eight (n,n′;r) of §0 with
Ψ₁ = det₂[v(−2),v(−3)]_{rows 2,3}·v₃(−1) (adjust the det₂ sign convention and the row pair {2,3} vs {1,2} to
the one your `psi1` uses — `psiSite[B,1,t,{1,2},1]` suggests rows {1,2} and component 1; take whatever makes
ψ_x·Ψ₁ = M_x hold exactly and record it). Hence N⁰ = Σφ⁰_xM_x, N¹ = Σφ¹_xM_x: build the tensors C⁰, C¹ as
6×4 matrices (rows = minors (n,n′) in the order (−3,−2),(−3,−1),(−2,−1),(−3,0),(−2,0),(−1,0); columns = r ∈
{−3,−2,−1,0}) with entries linear in φ, i.e. bilinear in (a,b). Verify on 64 states that
⟨Φ(g)|Ψ⟩/⟨x₀|Ψ⟩·Ψ₁ = Σ_{nn′,r}(C⁰+λ₂C¹)_{nn′,r}det₂[v(n),v(n′)]v_k(r) for two rotations and a random c.

**Task 2 — bracket dictionary and the forms (symbolic).** Define brackets ⟨Q^jF⟩ = Σ_n F_nQ^j(θ₁−nh) on the
four points with F a free 4-vector; show symbolically that det₂[⟨Q²F⟩,⟨Q³F⟩;⟨Q²G⟩,⟨Q³G⟩]·⟨Q³H⟩ equals
Σ(F∧G)_{nn′}H_r det₂[v(n),v(n′)]v₃(r) (fixes the orientation of F∧G relative to your minor convention), and
that the L=1-type 3×3 det[⟨Q²F_i⟩,⟨Q³F_i⟩,⟨Q³H_i⟩] expands to Σ_i ±(F_j∧F_k)⊗H_i. Write the Plücker quadric
Pf(A) = A₁₂A₃₄ − A₁₃A₂₄ + A₁₄A₂₃ in your minor order and the pairing A∧B for two 2-forms on C⁴.

**Task 3 — exact loci for N⁰ (rotations).** With C⁰(a,b) from Task 1 (impose b·a = 0 by b₃ = −(a₁b₁+a₂b₂)/a₃,
and treat a₃ ≠ 0; also check the chart a₃ = 0 separately if cheap):
  (i) (R2δ): Pf(C⁰_{·,0}) = 0 — a polynomial in (a,b); `Factor` it; describe the locus (dimension of the
      solution set in the 4-dim (a,b) cone; irreducible components; does it contain the whole rotation cone?
      — it should not, since generic rotations are not local determinants). Give a rational parametrisation
      of the main component if one exists (e.g. solve for one variable).
  (ii) (R1): C⁰_{·,0} ∝ C⁰_{·,−1} — exact solve; expect only degenerate rotations (identify them: which
      Gauss parameters vanish).
  (iii) (R3): Pf(C⁰_{·,0}) = 0 and C⁰_{·,−1}∧C⁰_{·,0} = 0 — exact solve; dimension and components; if a
      1-parameter family of genuine rotations survives, that is the "single 3×3 bracket determinant" family:
      extract F₁,F₂,F₃ (three 4-point functions) and H (supported on {−1,0}) for it symbolically, as
      functions of the family parameter. Cross-check with Part A's statement that the local 3×3 family is a
      4-dim cone and rotations are generically outside it (the intersection should be exactly this locus).
  For each locus: substitute a rational point of it, build the explicit bracket expression (F, G, H numbers on
  the 4 points), and verify the identity N⁰ = bracket form on all 64 states (relative deviation).

**Task 4 — N¹ and the λ₂-pencil.** (i) Confirm (R2δ) and (R3) are automatic for N¹ (both columns on three
points); give F, G per column. (ii) The composite: for r = 0, the conditions Pf(C⁰_{·,0}) = 0 and
C⁰_{·,0}∧C¹_{·,0} = 0; for r = −1 automatic (verify). Exact solve for the locus where N(λ₂) = Σ_{r∈{−1,0}}
det₂[F_r(λ₂),G_r(λ₂)]·Q³(θ₁−rh) for all λ₂ with λ₂-affine bra functions F_r, G_r. Dimension, components,
parametrisation; are there genuine (non-degenerate: all three Gauss parameters nonzero) rotations on it? At
a rational point of the locus, extract F_r(λ₂), G_r(λ₂) explicitly (affine in λ₂) and verify the full overlap
identity on all 64 states. (iii) Same with the (R3) shape for the composite: is N(λ₂)·(anything) a single
3×3 bracket determinant with λ₂-dependent bra data? (conditions: the four 2-forms C⁰_{·,−1}, C¹_{·,−1},
C⁰_{·,0}, C¹_{·,0} span a space of decomposable forms sharing a common vector — i.e. all pairwise ∧ vanish and
all Pf vanish; exact solve.)

**Task 5 — natural measures (report, do not force).** On each nontrivial locus found, test whether the bra
functions can be chosen "natural": F_n = wgt₁(n)·w^{−x_n/h}·(p₀ + p₁x_n) with wgt₁(n) = 1/Π_{m≠n}(m−n)h the
one-site SoV weights, x_n = θ₁−nh, (w, p₀, p₁) free (the det₆-entry type: measure × twist-exponential × linear
polynomial), for F, G (and H if present) — a 4-point function has 3 projective degrees of freedom and this
family has 2 (w and p₁/p₀), so each is one condition; solve exactly for the family parameter(s) where F AND G
are both natural (possibly with different w's — then say so). Also test the pure delta/shift dictionary
(F ∈ span of δ_n, i.e. F supported on ≤ 2 points). Report what is possible; negative results are results.

**Task 6 — z-independence and a second θ.** Re-derive Kvac, K20 at z′ (the companion-frame covectors are
z-free, K may depend on (h,θ)): check the loci of Tasks 3–4 are identical. If cheap, repeat Task 3(i) at
h = 1/3, θ = (1/5, −4/9) (P2) and say how the quadric's coefficients depend on t = θ₁₂/h.

**Task 7 — what the loci give for CG (markdown only, no chain run).** For the main (R2δ)/(R3)/pencil locus:
the undone site-1 covector relative to site 2 = x₀^{(1)}: in the diagonal frame x₀^{(1)} = e₁ᵀρ(U₀(z)) with U₀
the upper factor of the Gauss decomposition of P_phys⁻¹ (verify at z that it exists: e₁ᵀP_phys⁻¹ is of the form
e₁ᵀρ(upper unipotent) up to scale — one exact check); then Φ(g) in the diagonal frame is
e₁ᵀρ(gP⁻¹)⊗e₁ᵀρ(U₀) and the undone covector is e₁ᵀρ(gP⁻¹U₀⁻¹). Report, for a rational point of each locus,
which site-1 weights the undone covector reaches and the rank on (1,1,1) over 3 points of the locus (exact,
finite twist only). No untwisting here.

## 3. Reporting

End with a markdown cell `### site-1 rotation as ladder-1 brackets — result (2026-10-10)`: the monomial
table and sign conventions; the tensors' populated slots; for each form (R1), (R2δ), (R3), pencil, composite
3×3: the exact locus (polynomial, factorisation, dimension, components, parametrisation), the explicit bra
functions at a rational point, the 64-state verification; the natural-measure findings; z-independence;
the CG-reach table of Task 7. Save after each task. If VS Code does not execute cells of the new file, run
each cell's source in K1 via `wolfbook_evaluateExpression` and store printed results in a markdown cell under
the code cell; say so. Report back concisely, numbers not adjectives. Do not edit other notebooks, `paper/*.tex`,
git. If something is refuted (e.g. the monomial table fails), stop at that point and report what you see.

## 4. Gotchas

- The shared kernel K1 is the only kernel; the first cell wipes other notebooks' live state (files intact).
- Never write `*)` inside a comment. No `D`, `E`, `N`, `I`, `C`, `K` as symbols — use `CC0`, `CC1` for the
  tensors, `KK` etc.
- Memoized definitions cache bad results; purge with literal symbols.
- Exact Gröbner/`Solve` over (a,b) with the b·a = 0 substitution: eliminate the scale first (set a₃ = 1 and
  one b-component to 1 on a chart; check the other charts if a component of the locus could be at infinity).
- The reference cells use `prec`-digit eigenstates; keep every locus computation exact (rational K entries —
  cell 60 rationalises them and checks the error; reuse that).
- Wolfbook replies truncate at `{`; read long results with `StringReplace[ToString[.., InputForm], "\n" -> " "]`.
