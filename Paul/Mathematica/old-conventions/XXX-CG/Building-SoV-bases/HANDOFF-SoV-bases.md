# Hand-off: Building SoV bases / one-site localization (sessions of 2026-09-25)

For the next agent picking up this thread. Read this, then `main-notes.tex` (the rolling record of
results with formulas; **keep appending to it**, see memory `feedback_rolling_main_notes`). All live
computation is in `L2-localization-tests.wb` (this folder), in order, with a markdown "Findings" cell
after each block.

## Latest (2026-09-29, later session): SoV-basis operators close all the gaps

Formulas, tables and points checked are all in `main-notes.tex`, in the sections "L=1 left SoV basis for general representations" and "L=2: the L=1 SoV-basis operators on the twist-w vacuum".

**1. L=1 left B-eigenbasis for any [λ1,λ2,0]** (`Companion-B-general-rep.wb`).
- Fused transfer matrices come from the CBR formula: T_λ(u) = det T_{λ'_j+i−j,1}(u+h(i−1)), with T_{a,1} = t_a.
  This reproduces all of the hand-written adjoint t's.
- The basis states are ⟨x_μ| = ⟨lw| T_{μ11,0}(θ1) T_{λ2,0}(θ1)⁻¹ T_{μ21,μ12}(θ1), with GT labels {Λ,{μ11,μ12},{μ21}}.
- The denominator must be **T_{λ2,0}, not T_{1,0}**. The two agree when λ2 = 1.
- Checked: [2,1,0], [3,1,0], [3,2,0], [4,2,0] and [4,1,0] all give a complete eigenbasis.

**2. L=2 on the twist-w vacuum** (`L2-twistw-vacuum-SoV-ops.wb`).
- Apply the two-site O_μ to ⟨v_g| = (lw·g)⊗(lw·g), with r = r(w).
- The one-body operators are exactly 1, t11, and T_{n,0} T_{λ2,0}⁻¹ T_{λ2,λ2} for n = λ2 … λ1.
- The bare ratio T_{n,0}/T_{λ2,0} stops localizing, except in the degenerate case where w1 is a chain-twist eigenvalue.
- **The top operator O_max = T_{λ1,0} T_{λ2,0}⁻¹ T_{λ2,λ2}(θ1) reaches every irrep of Λ⊗Λ** for [3,1,0], [3,2,0], [4,2,0] and [4,1,0].
  This includes the old gaps [3,3,2] and [4,4,0] for [3,1,0], and [4,4,4] and [6,6,0] for [4,2,0], using full transfer matrices only.
- Robustness: 5 configurations of (point, w, k), including a new point 3.
  [3,1,0] and [3,2,0] are exact. [4,2,0] and [4,1,0] were computed mod two large primes, because the exact 729-dim LU is too slow.

**3. Single-site interpretation.**
- Every one-body operator acts as the same operator built for the site-1 chain alone, tensored with the identity on site 2:
  ⟨v_g|O(θ1) = c_O (⟨s|O^{L=1}(θ1)) ⊗ ⟨s|.
- c_O = ∏_{i<j}(x+ih−λ2h)(x+ih−h−λ1h) · ∏_{k=λ2}^{n−1}(x+kh−λ1h), with x = θ12 and j the column height of T_{j,j}.
  This is the site-2 lowest-weight value of L22·L33 for each column pair, and of L33 for each extra row box.
- c_O does not depend on w. The formula is fitted, not derived, but it holds in all 5 reps.

**Traps from this session.**
- Kernel globals clash between notebooks: `rOfW` already had a value from `transfer-algebra-dim.wb`. The new notebook uses `hh`, `th` and `rFromW`.
- An aborted long run left the kernel stuck "busy", and I had to restart it.
- At point 1, the single-site t1(θ1) is singular for [4,1,0]. The weight (3,2,0) gives z·m = 0.
- Mod-p code: `reachTestP` and `siteCompareP`. Both were validated against the exact [3,1,0] and [3,2,0] runs.

**4. The ratios are polynomial.** As an exact operator identity at θ1, T_{n,0}T_{λ2,λ2}(θ1) = T_{λ2,0}T_{n,λ2}(θ1), with coefficient 1.
- So O_max = T_{λ1,λ2}(θ1), and R t11 = T_{2,1}(θ1).
- It holds only at θ1: it fails at any shift s ≠ 0 and at generic u. At θ1, t3(θ1) = t3(θ1+2h) = 0.
- Checked for [2,1,0], [3,1,0] and [3,2,0] at points 1 and 2. Details are in the last section of `L2-twistw-vacuum-SoV-ops.wb` and in `main-notes.tex`.
- Not yet checked at operator level: [4,2,0] and [4,1,0] (use the mod-p code), and a symbolic proof.
- Open question: the generic-u difference. It is not of the form c·t3(u+ah)·T_{n−1,λ2−1}(u+bh). Ask Paul for the relevant Plücker/fusion relation.

**5. What "reach" means (don't overstate it).**
- The test takes one covector ⟨X| = ⟨v_g|O(θ1). An irrep ρ counts as reached if ⟨X|Π(k)|hw_ρ⟩ ≠ 0 for one generic k.
- Only the zero cases were saved. So "reaches every irrep" means there is a nonzero overlap with at least one copy of each irrep.
- It does **not** mean O generates all states, and multiplicity coverage is unknown.
- This is now explained in `main-notes.tex`, in a paragraph before the reach results. The full per-operator missed lists for [4,2,0] and [4,1,0] are in `main-notes.tex` too.

**Next steps.**
- The Q-function form of ⟨v_g|O_max|Ψ⟩. It is now just the eigenvalue of T_{λ1,λ2}(θ1), which is a Q-bilinear or CBR determinant. With the single-site result, this should reduce to an L=1 SoV functional times c_O.
- Derive c_O.
- Record multiplicity coverage: the rank of ⟨X|Π(k)|hw⟩ over each irrep's multiplicity space, not just the zero cases. **[2,1,0] (2026-09-30):** the generic-k test says full multiplicity, but that test is misleading because Π(k) mixes weights. The correct test is on the lowest-weight space, with no k. On that test {1, t11, T_{2,1}} over 8 values of w give 5 of the 6 directions (octet ×2 is OK), but they **miss the 10 [4,1,1]**: v = a g⁻¹ is identically zero at weights (2,0,1) and (1,0,2). Next: find an operator that feeds weight (1,0,2). See `main-notes.tex`.
- Check reps with λ1−λ2 ≥ 4, such as [5,1,0], if the pattern needs more confirmation.

## Earlier (2026-09-29): best result so far (Paul)

In `transfer-algebra-dim.wb` (this folder, last sections), companion frame, exact at points 1 and 2:
- Take the vacuum rotated globally by a generic g on **both** sites, (e8ᵀg)⊗(e8ᵀg), g = e^{r1E21}e^{r2E31}e^{r3E32}.
- Then ⟨v_g| R(θ1) t11(θ1) = ⟨a(r)| ⊗ ⟨e8ᵀg|. This is one-body for all r, and site 2 stays in the rotated vacuum.
- Its singlet overlap is ∝ P1(r)P2(r) ≠ 0. P1, P2 are polynomials in r whose coefficients come from the companion first row χ; the formulas are in `main-notes.tex`.
- So this is the first **one-body, singlet-reaching** covector built from a single site-1 operator.
- **Refinement:** choosing r = (−1/w1, 1/(w1w2), −(w1+w2)/(w1w2)) makes ⟨v_g| the reference state (left eigenvector) of the chain with companion twist G_c(w).
  - Then the singlet overlap is ∝ det(w1 − G_c(z))² / (w1⁶w2). It vanishes only if w1 is one of the chain's twist eigenvalues z.
  - Mechanism: Π(g) t_G Π(g)⁻¹ = t_{gGg⁻¹}, and the vacuum is an eigenvector for upper-triangular twists.
- **Why earlier attempts failed:** R(θ1) localizes on the rotated vacuum iff P1 = 0, which is exactly where R t11 misses the singlet. The old vacua all required R to localize.
- **v components:** after undoing g globally, v = a g⁻¹ has all 8 components nonzero, in both the companion and the diagonal frame. A further w0 rotation puts site 2 at e1.
- **Other reps** (full t's only, one generic r, point 1):
  - [3,1,0]²: t21(θ1) ≃ t11(θ1)t10(θ1+h) reaches the missed [4,2,2], but not [3,3,2] or [4,4,0]. Those need companion-piece triples P2 P1 P1.
  - [3,2,0]²: **t22(θ1)·t10(θ1+2h) reaches every irrep**, including the missed [4,3,3].
  - [4,2,0]²: the same operator reaches the missed [6,3,3] and [5,4,3]. [4,4,4] (the singlet analog) and [6,6,0] stay missed, even with depth-3 P triples.
  - The summary table is in `main-notes.tex`, section "Summary: rotated vacuum + full transfer matrices".
- **Saved data:** `transfer-algebra-dim-results.m` (plain Wolfram; `Get` it). It has the SoV coefficients, n_x, P1/P2, twist-w, the v components and the per-rep operator tables.
- **Kernel note:** the notebook's later cells (`buildRepV`, `scanRepV`) are self-contained. The rep scans take about 35 s for dim 225, and for [4,2,0] about 10 min at depth 2 and about 15 min for the depth-3 extension.
- Next steps:
  - Q-form of the ⟨v_g|R t11|Ψ⟩ functional (r-deformed bracket?) and of t22(θ1)t10(θ1+2h).
  - The SoV / single-site content of a(r).
  - The λ1−λ2 = 2 gaps: try larger fused diagrams (three columns) and other ratios. Ask Paul for the fusion conventions first.
  - Confirm the rep scans at point 2 and at a second r.

Same session, also in `transfer-algebra-dim.wb` + `main-notes.tex`:
- Transfer-matrix algebra dimensions at θ1+rh.
- The θ1+rh polynomial for the singlet is dense (it equals Y8).
- e1⊗e8 = two SoV terms.
- Site-1 rotation in the SoV basis (12 terms, closed form in T).
- A negative result: no contour form from the ν-measure.

## The goal (Paul)

Get Clebsch–Gordan coefficients for gl(3), two sites, adjoint ⊗ adjoint, by reducing the two-body
problem to a **one-body** problem. Mechanism: transfer-matrix actions that *localize*, i.e. map a
product covector ⟨a|⊗⟨b| to another product covector changing only one site. Then the overlap with a
Bethe state Ψ is a Q-function functional, and moving everything to one site on highest-weight states
gives CG data. Everything must stay one-body: **Paul explicitly rejects anything that needs operators
on the coupled tensor-product space**, i.e. entangled covectors.

The hard target is the **singlet** in V⊗V. All standard families miss it.

## Files

| file | content |
|---|---|
| `main-notes.tex/.pdf` | rolling record, every result with formula, status and points checked |
| `L2-localization-tests.wb` | all L=2 computations (≈145 cells) |
| `Companion-B-C-diagonalize.wb`, `diagonalise-B-C-diagonal.wb` | L=1 notebooks cleaned from Paul's `.nb` files (the `.nb` originals are deleted) |
| `new-sov-basis-tests.wb` | L=1 swapped-basis tests (xbar, ybar) |
| `../L2-mathcalA.*`, `../L2-det-mathcalA.*`, `../L2-localizing-transfer.*`, `../xxx-210-cg-test.*` | earlier site-1 results this thread builds on |
| `Paul/Papers/SoV_2.pdf` | Ryan–Volin arXiv:2002.12341, eq. (5.22) is the SoV wave function used below |

## Conventions (details in `main-notes.tex` §1)

- Lax operator L_ij(u) = u δ_ij − h E_ji, and T = L1(u−θ1)·L2(u−θ2)·G. Site 1 is the left tensor factor.
  The GT basis has the highest weight e1 first and the lowest weight e8 last.
- Two twist frames. **Diagonal**: G = diag(z). **Companion**: G_c = P·G_d·P⁻¹, with P the
  Vandermonde matrix of the (z_a², z_a, 1). They are related by t_c = Π(P)·t_d·Π(P)⁻¹ (checked).
- `t` is the fused tower at θ1. `tS` is the shifted tower at θ2+2h, which is the one that localizes
  on site 2. R = t20·t10⁻¹, and T := θ12/h.
- Q's come from the Baxter operators of `L2-mathcalA.wb`: lower Q_i = z_i^{u/h} q_i, upper
  Q^j = z_j^{−u/h} q^j, monic.
- Numerical points: point 1 is h=1/3, θ=(2/7, −5/11), z=(2, −3, 5). Point 2 is h=1/2,
  θ=(3/13, −4/17), z=(3, −7/2, 11/4). Results are checked at both unless stated.

## What worked (details and formulas in `main-notes.tex`)

1. **Site-2 localization.** tS at θ2+2h localizes on site 2, and t at θ1 localizes on site 1.
   R only localizes with the right vacuum. In the companion frame that is the unrotated lowest-weight
   state e8⊗e8; in the diagonal frame it is the **rotparams** vacuum, i.e. ⟨lw|Π(P).
2. **Site-2 Q functionals.** eig Rˢ = Q³(θ2)/Q³(θ2+h), and eig tS11 = Q²³(θ2+h)/Q²³(θ2+2h).
   The brackets are ∮Q³ k^{−u/h}, with poles θ2, θ2+h, normalized at θ2+h; and ∮Q²³ k^{−u/h}, with
   poles θ2+h, θ2+2h, normalized at θ2+2h. Hodge duality lets both be rewritten in terms of Q12 and Q1.
   ⟨Ω|tS11 = c12·(single-site tS11 on site 2), with c12 = (θ21+h)/(θ21+2h).
3. **Transfer-matrix eigenvalues as Q-bilinears.** τ1 = K/ν2 · Σ κ_a Q_a(x+h) Q^a(x−h), and
   τ2 = K′ · Σ κ′_a Q^a(x+h) Q_a(x−2h).
4. **True SoV basis, companion frame.** B = T11·QM12,12(u+h) + T21·QM12,13(u+h) with the
   **untwisted** T. The recipe states ⟨0|O_g1(θ1)O_g2(θ2), with ⟨0| = e8⊗e8, are exactly its
   eigenvectors. The wave functions match RV (5.22):
   ∏_s Q1(θs+r3 h) · D_s(r1, r2), where r1 ≤ r2 ≤ r3 are the B-roots per site.
5. **The site-1-localized space is 8-dimensional**, with explicit t-forms; table in `main-notes.tex`.
   **Row 8 reaches the singlet and is one-body:**
   ⟨0| R(θ1) [ t21(θ1) − 1/(T+1) · t11(θ1) R(θ2) ] = ⟨v| ⊗ e8ᵀ.
   Its Q-form is C · Q1(θ1+2h)/Q1(θ1+h) · [D(0,2) + D(0,1)·ρ2/(T−1)] / D(−1,0), with
   ρ2 = Q1(θ2+2h)/Q1(θ2+h). Checked on all 64 states at both points.
6. **Replacing t11(θ1) in the ⟨k1Q1⟩⟨k2Q12⟩ family.**
   - No *linear* combination of single-site quantum minors reaches the singlet. The weight argument:
     a minor has distinct rows, hence at most one E_j3.
   - *Quadratic* insertions do reach it: t11(θ1)², fused t22(θ1±h), and others. These must be
     conjugated by g0, because they act before g0 is undone.

## What did NOT work (and why; don't redo blindly)

- **Diagonal frame with generic rotation angles.** R does not localize there. You need rotparams:
  p3 fixed by p1 and p2.
- **Lowest-weight vacua with the original x, y.** This is the same pair as rotparams, giving the same
  overlap matrix.
- **Swapped bases (xbar, ybar) at L=1.** The vacuum-orthogonality condition is solvable, but the
  overlap matrix is dense.
- **Products of site-1 and site-2 families.** They reach the singlet, but every singlet-reaching
  product is **entangled** (rank 2–4). Only R·Rˢ stays a product, and it has zero singlet overlap.
  Paul rejected this route.
- **Iterating localization.** V_loc(O; w) is never invariant under O; it is a one-step property.
  Nontrivial localization only happens for w a rotated extremal covector.
- **Compact Q-forms in the recipe/diagonal basis for the dressed t10 combinations.** Single ratios,
  products of ≤3 ratios, reweighted bilinears, site-1-only data and the site-1 SoV span all fail.
  The fix was to switch to the true SoV basis (item 4 above).
- **A B built from the twisted monodromy T·G_c.** It commutes, but its wave functions do not separate.
- **A double-contour integral for the singlet functional.** Factorized rational insertions cannot give
  c = 1/(T−1). They are degree 0 in T, whereas c → 0 as T → ∞. This points to single-variable
  brackets with poles on **both** sites (L2-det-mathcalA style). That is the next ansatz, but it has
  **not been tried**.

## Traps hit this session (save yourself the time)

- **Wolfbook cell parsing.** A multi-line statement whose line ends in ` .` (Dot) is silently cut
  there. Put Dot-chains on one line or wrap them in parentheses.
- **MCP replies are truncated to the first line.** Return `ToString[..., InputForm]` on a single line,
  or slice the output.
- **`MatrixPower[nilpotent, 0]` is not the identity.** Use an explicit `mpow`.
- **Separability or proportionality tests of B-eigenvectors must allow a per-label constant.** Each
  eigenvector ⟨x| has an arbitrary scale. Test that ratios are constant across states, not that they
  equal 1. This caused a false "does not separate" once.
- **Each Ψ_A has an arbitrary normalization.** Always divide by a reference overlap before comparing.
- **Read the shift-scan direction carefully.** I once recorded the tS11 ratio upside down
  (corrected in `main-notes.tex`). The convention is C = ρ·Q(θ+bh)/Q(θ+ah).
- **Kernel state.** Everything runs in one shared kernel. `setPoint[p]` resets the memoized operators
  for a new point, but `computeSpectrum[]` must be rerun for evecs, tauPoly and the Q's. Gc, Pv and
  psiC are recomputed inside `sovPipeline` and `singletFromQ`. The kernel is currently at
  **point 2**, apart from the last few cells.
- **Diagonal-twist transfer matrices conserve weight.** Single-site t's only change weight after
  conjugation by g0.
- **Paul's working style.** Stop and ask before open-ended ansatz searches. He often knows the
  structure: he supplied the RV wave-function paper when the factorization fit stalled. Keep runs
  short.

## Open / suggested next steps

1. A contour-bracket (k-deformed) form of the singlet functional, using single-variable brackets
   with poles on both sites. Ask Paul which measure: the L2-det-mathcalA one, or FSoV from Gromov–Primi–Ryan.
2. A Q-functional that produces a quadratic one-body insertion (for example t11(θ1)²) when λ2 = 1,
   or reuse the ⟨k2 Q12⟩ bracket twice.
3. The site-2 mirror of the SoV-basis localization table.
4. The explicit one-body covector v for row 8, and the RV Gamma-factor normalization of C.
5. Symbolic θ, h, z, and other reps. Everything here is at two numerical points.
