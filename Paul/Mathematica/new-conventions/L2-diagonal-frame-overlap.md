# gl(3) XXX, L=2, adjoint⊗adjoint: the overlap formula in the diagonal frame

These notes cover `L2-diagonal-frame-overlap.wb` and its Hermitian-generator copy `L2-diagonal-frame-overlap-EH.wb`, from work on 2026-10-01. The design spec is `docs/superpowers/specs/2026-10-01-L2-diagonal-frame-overlap-design.md`. The L=1 counterpart is row 19 of `Companion-B-C-diagonalize-clean.md`.

## Files

| File | Role |
|---|---|
| `L2-diagonal-frame-overlap.wb` | E-basis generators. It clears the kernel and runs in about a minute. Its summary table has 37 checks, of which only 4 FAIL, and those 4 are deliberate (see "What does not work"). |
| `L2-diagonal-frame-overlap-EH.wb` | The same notebook with Hermitian generators EH (see "EH variant"). 39 checks, the same 4 deliberate FAILs. |

**Numeric point.** All values are exact rationals, with 60 digits wherever floats enter:
- h = 1/2;
- θ = (2/7, −3/5);
- z = (3/2, 5/4, 2/3);
- w = (5/2, 7/3, 3/5).

Only this one point has been checked.

## Setup

- **Lax and monodromy.** ℒ_ij(u) = u δ_ij + h E_ij, and T(u) = ℒ⁽¹⁾(u−θ₁) ℒ⁽²⁾(u−θ₂), a product in auxiliary space.
- **Twist.** Diagonal, T^z = T·diag(z). The companion frame is not used anywhere.
- **Transfer matrices.** t_a has degree 2a in u. Each is rebuilt exactly by Lagrange interpolation from 7 rational sample points, each an exact 64×64 product, which takes about 4 s.
- **States.** The 64 eigenstates are found one weight sector at a time; there are 19 sectors. The diagonal twist conserves the weight, so μ is simply the label of the sector.
- **The bra.** Ω = Ω_zw ⊗ Ω_zw, where Ω_zw is the L=1 vacuum of G′ = P_aux G(w) P_aux⁻¹, i.e. the w-twist written in the frame where z is diagonal. It is an exact left eigenvector of every t_a(G′), with t₁ = Σ_j w_j ν_j(u).

## Results

| # | Statement | Status |
|---|---|---|
| 1 | **Quantum determinant.** t₃ = χ₃ ν₁(u)ν₂(u−h)ν₃(u−2h), with ν_n(u) = ∏_s (u−θ_s+hλ_n). The t_a commute and conserve the weight | exact |
| 2 | **Baxter equations.** The L=1 direct and dual equations with each ν_n replaced by the product over sites. deg q_i = 2λ₁−μ_i and deg q^j = μ_j−2λ₃; the solver finds the minimal degree, so this is tested, not assumed | all 64 states |
| 3 | **QQ relation.** Q_aQ_b^{[2]} − Q_bQ_a^{[2]} = (1/z_b−1/z_a) F Q^c, with F = ∏_s (u−θ_s+2h) | all 64 states |
| 4 | **Orthogonality.** With the brackets below, the 18×6 matrix (9 pairs × 2 brackets) has rank 5 for A≠B, with σ_min/σ_max ≲ 10⁻⁶³ | all pairs |
| 5 | **Overlap formula**, given below | 50–51 digits, all 64 |

**The overlap formula:**

Ω·Ψ_B = [ (Ω·y₀) / (c₀ ⟨x₀|y₀⟩² Φ_Q^vac(y₀)) ] · det N[B, vac] / Ψ_Q^B(x₀) · ⟨x₀|Ψ_B⟩

- **x₀ and y₀** are the L=1 reference vectors, taken as products over the two sites:
  - x₀ = (e₁ᵀ P_phys⁻¹)^{⊗2};
  - y₀ = y₀(θ₁) ⊗ y₀(θ₂), where each factor is the L=1 y_{2|21} built at that site in the companion frame and rotated by P_phys.
- **Ψ_Q(x₀) = ∏_s [Q²(θ_s−2h)Q³(θ_s−3h) − Q³(θ_s−2h)Q²(θ_s−3h)] · Q³(θ_s−h)**, using the ket's dual Q's.
- **Φ_Q(y₀) = ∏_s [Q₁(θ_s−2h)Q₂(θ_s) − Q₂(θ_s−2h)Q₁(θ_s)] · Q₁(θ_s−2h)**, using the vacuum's direct Q's.
- **⟨x₀|y₀⟩ per site** is 9h⁴ with E generators and 3h⁴ with EH. The plain product over sites works; the inverse x/y overlap matrix was never needed.
- **det N** is a 6×6 matrix:
  - Rows are (i,j; k), with (i,j) ∈ {(1,3),(1,2),(2,3)} (the L=1 pairs) and k = 1, 2 labelling the two brackets.
  - Columns are the three L=1 functionals −⟨PQ⟩, ⟨PQ^{[2]}⟩ and ⟨P^{[−2]}Q⟩, each multiplied by x^{β−1} for β = 1, 2.
  - P = Q^j is the ket's dual Q and Q = Q_i is the vacuum's direct Q.
- **The brackets** (Paul's preferred choice, Section 9 of the notebook):
  - Both are residue sums over **all** poles of ρ = ∏_s ∏_{n=0}^{2} (x−θ_s+nh)⁻¹.
  - ⟨f⟩₁ carries the insertion 1, and ⟨f⟩₂ carries e^{2πix/h}.
  - On ladder s the insertion is the constant ε_s = e^{2πiθ_s/h}.
- **c₀** is the common coefficient of the 16 monomials of Ψ_Q(x₀)Φ_Q(y₀) in det N, computed symbolically:

  c₀ = (ε₁−ε₂)³ / [16 h¹² θ₁₂³ (θ₁₂²−h²)⁴ (θ₁₂²−4h²)²],  θ₁₂ = θ₁−θ₂

  - With per-ladder brackets instead (Section 7), the factor (ε₁−ε₂)³ is replaced by −1, so c₀ is (1/4h⁶)² times factors that come from the residue weights coupling the two sites.
  - The two versions differ by det[[1,1],[ε₁,ε₂]]³, as expected from the row mixing.
  - Fitted value at the numeric point: R = (Ω·y₀)/(c₀⟨x₀|y₀⟩²Φ_Q^vac(y₀)) ≈ 9.1365×10⁻⁴, the same in both bases.

## What does not work (kept in Section 6 on purpose)

- **One bracket summed over both ladders.** With insertion 1 only, plus x-power columns, the 9×6 matrix still has rank 5 for A≠B. But det N then depends on which 6 of the 9 rows are kept, and the dependence varies from state to state, even at A=B. No reference factor fits; the spread is ~10⁸.
- **The row choice {(1,3),(2,1),(3,2)}.** It is not equivalent to the L=1 pairs. At L=1, 7 of the 84 row choices failed as well.

## EH variant

- **Generators.** EH_ij = R E_ij R⁻¹, with R = √M = diag(1,1,1,√2,√6,3,2,3) in `ps` order. M is the unique diagonal solution of M⁻¹E_ijᵀM = E_ji with M₁ = 1.
- **Hermitian structure.** At L=2, t_a(u)ᵀ is the transfer matrix of the site-reversed monodromy ℒ⁽²⁾(u−θ₂)ℒ⁽¹⁾(u−θ₁), for all a. It is not t_a itself. See `Hermitian-structure-abstract.md`.
- **What changes.**
  - ⟨x₀|y₀⟩ becomes 3h⁴. The seed vectors e₁ and e₈ are the same GT states, so the overlap picks up R₁/R₈ = 1/3.
  - Everything else is basis-independent and unchanged.
- **Trap.** Exact matrices now contain surds, and Mathematica does not combine them canonically (√6/2 → √(3/2) sits next to √6). So `Union[...] === {0}` reports false failures. Use `zeroQ`, which applies `PossibleZeroQ` to the distinct entries.

## Implementation notes

- **Computing c₀.** Set every Q-symbol that is not in the target monomial to 0 *before* taking the 6×6 determinant. This leaves the target coefficients unchanged and is about 100× faster than expanding first.
- **Insertion order.** After inserting k cells after cell N, insert the next one after N+k.

## Section 10 (EH notebook only): localized action at θ₁−2h and the highest-weight CG data

All of this is in `L2-diagonal-frame-overlap-EH.wb`, Sections 10–10f (cells 28–49). It uses the kernel state from cells 3–14 (setup, the t's with twist diag(z), eigenstates, Ω). Sections 5–9 are not needed. Notation: R = e^{r₁E₁₂}e^{r₂E₁₃}e^{r₃E₂₃}, so Ω_zw = e₁ᵀR, with r_i(z,w) as in Section 4. θ₁₂ = θ₁−θ₂.

**10. Localization on site 1.**
- The transfer matrices are built with the **ket's** twist diag(z), so ⟨Ω|T_λ|Ψ_B⟩ = T_λ^B ⟨Ω|Ψ_B⟩.
- They are fused by CBR: T_(1,1)(u) = t₂(u) and T_(2,1)(u) = t₂(u)t₁(u+h) − t₃(u+h).
- **Result:** ⟨Ω|T_λ(θ₁−2h) = f_λ · ⟨v_λ| ⊗ ⟨Ω_zw|, where ⟨v_λ| = ⟨Ω_zw|T_λ^{L=1,θ₁}(θ₁−2h) (one-site chain, same twist).
  - **Exact** at the numeric point.
  - **Nontrivial:** v_λ is not proportional to Ω_zw.
  - **Not local on site 2.**
- **Site-2 factors:** f_(1,1) = (θ₁₂−h)(θ₁₂−3h) and f_(2,1) = (θ₁₂−h)²(θ₁₂−3h). These are **fitted at one point only** (θ₁₂ = 31/35, h = 1/2).

**10b. Undoing R on site 1.**
- **Identity:** ⟨v_λ|R⁻¹ = e₁ᵀ T_λ^{L=1}[g⁻ᵀ diag(z) gᵀ](θ₁−2h), with g = R in the fundamental. This follows from the invariance of ℒ = u + hΣe_ij⊗E_ij under g⁻ᵀ⊗π(g). The other three conjugations fail.
- **Support of ⟨v_(1,1)|R⁻¹:** only GT states 1, 3 and 5, i.e. weights (2,1,0), (1,2,0), (1,1,1).
- **Support of ⟨v_(2,1)|R⁻¹:** all 8 states.

**10c. Highest-weight projections of c_λ = ⟨v_λ|R⁻¹ ⊗ e₁ᵀ.**
- **Decomposition:** [2,1,0]² = [4,2,0] ⊕ [4,1,1] ⊕ [3,3,0] ⊕ 2[3,2,1] ⊕ [2,2,2], so the HW space is 6-dimensional.
- **What (1,1) reaches:** [4,2,0], [3,3,0] and [3,2,1], but not [4,1,1] or the singlet. This is forced by weights.
- **What (2,1) reaches:** everything.
- **[3,2,1]:** the (1,1) and (2,1) projections are independent.

**10d. Symbolic z, w, h.**
- **Where θ enters:** at u = θ₁−2h only u−θ₁ = −2h appears, so ⟨v_λ|R⁻¹ is independent of θ. θ enters only through f_λ.
- **Checks:** the symbolic result agrees with the numeric point, and the 10b identity holds symbolically.
- **Overlaps.** Write z_ij = z_i−z_j and X = w₂(w₁−z₁)z₁₃z₂₃/[(w₁−z₃)(w₂−z₃)]:

| μ | (1,1) | (2,1) |
|---|---|---|
| [4,2,0] | 3h²z₂z₃ | −3h³z₂z₃² |
| [4,1,1] | 0 | −3h³z₂z₃(z₂−w₁)(z₂−w₂)z₁₃z₂₃ / [√2 z₁₂(w₁−z₃)(w₂−z₃)] |
| [3,3,0] | 3h²(w₁−z₁)z₁₂z₃ / [√2(w₁−z₂)] | −3h³(w₁−z₁)z₁₂z₃² / [√2(w₁−z₂)] |
| [3,2,1]₁ | 2√(3/7) h²X | −3√(3/7) h³z₃X |
| [3,2,1]₂ | −√(3/35) h²X | −√(3/35) h³(2w₂−7z₂)z₃X/w₂ |
| [2,2,2] | 0 | 3h³w₂(w₁−z₁)²z₁₃²z₂₃² / [2√2(w₁−z₃)²(w₂−z₃)²] |

- **Basis dependence:** the [3,2,1] coordinates are in a Gram–Schmidt HW basis. Only zero versus nonzero and the rank are basis-independent.

**10e–10f. Small twists, z_k = e^{εφ_k}, w_k = e^{εω_k}.**
- **Leading orders** ε^#, for [4,2,0], [4,1,1], [3,3,0], [3,2,1], [2,2,2]:
  - (1,1): 0, —, 1, 1, —
  - (2,1): 0, 1, 1, 1, 2
- **Pattern (fit only):** # = max(4−μ₁, μ₃).
- **Directions:** at leading order, each irrep's direction does not depend on ω. Varying ω only rescales.
- **Weight sectors (10f):** e_jᵀ⊗e₁ᵀ has 7 weight sectors. All are 1-dimensional except total weight (3,2,1), GT states 4 and 5.
  - With both λ and ω sampled (φ = (1/3,−2/7,5/11), h = 1, 8 random ω), the leading-ε vectors span **every** sector.
  - The (3,2,1) sector needs both λ: each λ alone gives rank 1 there.
- **Cell 49:** HW_μ → (e_j⊗e₁ components) is injective for every μ.
  - So the leading-ε data determine all highest-weight CG data of [2,1,0]².
  - Paul's conclusion: these vectors generate all the CG coefficients of highest-weight states.

**Open from Section 10:**
- **Second numeric point.** Re-check localization and the f_λ fit at a second numeric point, or symbolically.
- **Overlaps → CG.** Assemble the full chain, ⟨Ω|Ψ_B⟩ formula → T_λ^B(θ₁−2h) → conjugation by R → ε → 0 → HW CG. Paul is doing this with another agent.
- **Other reps.**

## Section 11 (EH notebook only): small-twist limit and the CG chain

Code: `L2-diagonal-frame-overlap-EH.wb`, cells 50–78. These cells need only cells 3–14 run first; Sections 5–10 are not needed. Everything is numeric at 160–300 digits. Two chain points are used: P1 = (h, θ) = (1/2, (2/7, −3/5)), which is the main one, and P2 = (1/3, (1/5, −4/9)). Twists are z_k = e^{εφ_k} and w_k = e^{εω_k}, with t = θ₁₂/h. The states kept are those whose ε → 0 limit is a highest-weight (HW) state. There are six: [4,2,0], [4,1,1], [3,3,0], [3,2,1]± and [2,2,2]. The two [3,2,1] copies are ordered by the eigenvalue of t₁+πt₂ at u = 3/7.

**Helper functions** (cells 51–77):
- **Overlap ratio and twist dependence:** `smallTw` gives ρ and its pieces at one ε; `leadC`, `facs` and `Xpred` handle the leading coefficient, the twist monomial and the constant X_B.
- **The CG chain:** `cgData` gives the chain quantities, `cgParts` splits them by weight, and `wData` gives the projected covector W_μ.
- **Reconstruction:** `wCov` gives the localized covectors for one w, `x0Cov` gives x₀, `psiStates` gives the states with ⟨x₀|Ψ⟩ = 1, and `recon` reconstructs from many w.
- **Second chain point:** `clearTm[]` drops the memoized monodromy; the cells using it swap (h, θ) temporarily and then restore them.

### 11. The mixed-twist ratio ρ_B (cells 50–55)

The ratio is ρ_B = ⟨Ω|Ψ_B⟩/(⟨Ω|y₀⟩⟨x₀|Ψ_B⟩). By Sections 7–9 it equals det N[B,vac]/(c₀⟨x₀|y₀⟩²Φ_Q^vac Ψ_Q^B), and it does not depend on how Ψ_B is normalized. As ε → 0:

  ⟨x₀|y₀⟩² ρ_B → ε^# · X_B · ∏_{i,j}(ω_i−φ_j)^{n_ij},  with n_ij = max(d_j − D_i, 0).

- **The exponents.** d_j = 2λ₁−μ_j are the direct Q-degrees of Ψ_B, and D = (0,2,4) are those of the w-vacuum.
- **The leading power.** # = Σ n_ij = 8, 8, 8, 7, 7, 6.
- **The constant.** X_B = 1, 1, 1 for [4,2,0], [4,1,1], [3,3,0]. For the octet copies, X± = (2±√(4+5t²))/t². For the singlet, X = −8/(t²−1).
- **Evidence.** The monomial was fitted over 26 random (φ,ω) directions, with log-residuals ≲ 10⁻¹³. X_B was checked at P1 and P2 to about 10⁻¹³.
- **Where the ε⁸ comes from.** For Ψ = e₁⊗e₁ two things hold exactly: ⟨Ω|Ψ⟩ = 1 and ⟨x₀|Ψ⟩ = (z₂−z₁)². So the ε⁸ is pure reference-vector normalization: x₀ and y₀ are the companion-frame vectors carried by the singular P_phys, and ⟨Ω|y₀⟩ ~ ε⁻¹⁰.

### 11b. The y₀-free ratio R_B = ρ_B/ρ_[4,2,0] = (z₂−z₁)²⟨Ω|Ψ_B⟩/⟨x₀|Ψ_B⟩ (cells 56–59)

| B | leading R_B | order |
|---|---|---|
| [4,2,0] | 1 | ε⁰ |
| [4,1,1] | (w₁−z₂)(w₂−z₂)/[(w₁−z₃)(w₂−z₃)], **exact at finite twist** | ε⁰ |
| [3,3,0] | (w₁−z₁)/(w₁−z₂) = r₁, **exact at finite twist** | ε⁰ |
| [3,2,1]± | X± (w₁−z₁)/[(w₁−z₃)(w₂−z₃)] | ε⁻¹ |
| [2,2,2] | −8/(t²−1) · (w₁−z₁)²/[(w₁−z₃)²(w₂−z₃)²] | ε⁻² |

- **Evidence.** Checked from the vectors at P1 and P2 to about 10⁻¹³. The two exact rows hold to 54 digits at ε = 1/3.
- **Why those two rows are exact.** In the 1- and 2-dimensional weight sectors, Ω and x₀ are exactly parallel. So there R depends only on z and w, and is the same for every state of that weight.

### 11c. The CG chain at small twist (cells 60–62)

The chain is c_λ·(R⊗R)Ψ_B = T_λ^B(θ₁−2h)⟨Ω|Ψ_B⟩/f_λ, with c_λ = ⟨v_λ|R⁻¹⊗e₁ᵀ. Section 10e replaced (R⊗R)Ψ_B by its HW limit.

- **The identity is exact at small twist** (agreement to 160 digits). This also confirms f_λ at new twists.
- **x₀ and Ω only see the 27.** Both are coherent states (e₁ᵀπ(g))^{⊗2} in the dual of [4,2,0], so they see only the [4,2,0] admixture of Ψ_B. In particular ⟨x₀|h_μ⟩ = 0 exactly for μ ≠ [4,2,0].
- **Replacing (R⊗R)Ψ by its HW limit is wrong at leading order** for every state except [4,2,0]:
  - **Where the CG overlap is zero.** For (1,1) at [4,1,1] and at [2,2,2], c·h_μ is exactly 0, yet the chain gives ε¹ and ε² respectively.
  - **Everywhere else.** The ratio chain/CG tends to κ_{B,λ} ≠ 1, and κ does not depend on φ or ω.
  - **[3,3,0] (both λ) and [4,1,1] (2,1):** κ = (t+1)/2, which is 97/70 at P1 and 22/15 at P2.
  - **[2,2,2] (2,1):** κ = (t+1)(t+3)/12, **pair-fitted**.
  - **Octets:** κ is about (1.529, 1.072) for one copy and (0.642, 0.772) for the other, at P1.

### 11d. Weight decomposition with ⟨x₀|Ψ⟩ = 1 (cells 63–65)

With this normalization ‖Ψ_B‖ ~ ε^{−k}, k = 2, 3, 3, 4, 4, 6. The chain splits by weight as c·(R⊗R)Ψ = Σ_{ν≥μ} c^{(ν)}·[(R⊗R)Ψ]_ν.

- **The ν = μ term** is exactly c^{(μ)}·Ψ, and its leading part is ε^{−k} times the CG number.
- **The [4,2,0] term ties with it.** c^{(420)} is O(1) while c^{(μ)} is O(ε) or O(ε²), and the higher-weight pieces of (R⊗R)Ψ are only one order down. So the ν = [4,2,0] term always has the same order as the CG term. For example, at [3,3,0] both are ε¹·ε⁻³ and ε⁰·ε⁻², i.e. ε⁻² each.

### 11e. The projected covector W_μ = ⟨v_λ⊗Ω_zw|P_μ (cells 66–68)

- **Sizes.** ‖W_μ‖ is O(1). W_μ·h_μ equals c_λ·h_μ, which is the CG number and is of order ε¹ or ε².
- **Consequence.** The leading direction w₀ of W_μ is orthogonal to h_μ. So W_μ·Ψ = ε^{−k}(ε w₁·h_μ + ε w₀·Ψ⁽¹⁾ + …), where Ψ⁽¹⁾ is the O(ε) correction to the eigenstate. The second term is the contamination.

### 11f. Coincident φ and ω (cells 69–71)

Ω_zw does not depend on w₃. The allowed single coincidences are ω₁ = φ₁, ω₂ = φ₁ and ω₂ = φ₂. In each case either an irrep vanishes exactly on both sides, or κ is unchanged. **None of them helps.**

### 11g. Killing the [4,2,0] component (cells 72–74)

From 10d, c₂₁·e₁e₁ = −h z₃ c₁₁·e₁e₁ exactly. So C = c₂₁ + h z₃ c₁₁ has no [4,2,0] component. On the Q side it is (T₂₁/f₂₁ + h z₃ T₁₁/f₁₁)⟨Ω|Ψ⟩.

- **[4,1,1] and both [3,2,1] copies:** κ_C = 1 for 3 twist directions and at P1 and P2. These states are clean.
- **[4,2,0] and [3,3,0]:** C vanishes on both sides. [4,2,0] is already clean with either λ.
- **[2,2,2]:** κ_C = (t+1)/3, **pair-fitted** (97/105 at P1, 44/45 at P2).

### 11h. Reconstruction from all w (cells 75–78)

Here only the allowed inputs are used: W_λ(w)·Ψ for both λ and 20 random w, plus ⟨x₀|Ψ⟩ = 1. With these, Ψ_B is determined on its weight space up to the exact null space of the covector family. Rank, with x₀ included:

| sector | rank / dim | outcome |
|---|---|---|
| [4,2,0] | 1/1 | exact reconstruction |
| [4,1,1] | 2/2 | exact reconstruction, no κ |
| [3,3,0] | 2/2 | exact reconstruction, no κ |
| (3,2,1) | 5/6 | one exact null vector n(ε) |
| (2,2,2) | 7/10 | three null vectors |

- **The (3,2,1) null vector.** Its singular value is 10⁻³⁰⁸ at 300 digits, and it is present at every ε, including ε = 1/3. At finite ε it is not HW, but as ε → 0 it tends to a vector n₀ in HW_[3,2,1].
  - **So the octet mixing is undetermined by these data**, at any order in ε and for any w.
  - Fixing the free coefficient by minimizing the non-HW part fails: the direction error stays at 0.024 and 0.049.
  - The single clean functional on the octet HW space is C, and C·n₀ = 0.
- **The (2,2,2) null vectors.** They overlap the singlet, so the singlet's direction is trivial but its scale is not fixed (about 7.6% off).
- **Conclusion.** With T₍₁,₁₎ and T₍₂,₁₎ at θ₁−2h and Ω_zw for all w, every HW CG datum is recoverable except the octet mixing and the singlet scale. One more independent ingredient is needed: something site-2 localized, the left-vector/norm formula, or a T localizing at another point.

### Traps (Section 11)
- **Surds in y₀.** The exact EH surds in the one-site Lax that builds y₀ made `N` fail with `meprec` at h = 1/3. The fix is `EeN = N[Ee, 400]`.
- **Ordering the octet copies.** t₁(0) is the same for both [3,2,1] copies, so it can't order them. Use the eigenvalue of t₁+πt₂.
- **Comparing octet copies.** Don't re-sort them by value after multiplying by the twist monomial; that can flip the order.
- **Summary-table FAILs.** The ε-limit checks are O(ε) residuals of about 10⁻¹³, so they show FAIL against the summary table's 10⁻²⁵ threshold.

### Open (Section 11)
- **The octet mixing.** It needs one input beyond the site-1-localized T's. **Closed in Section 12** by site-2-localized t's.
- **Closed forms.** The singlet κ's are pair-fitted, and the octet κ's are not identified. Deriving them from first-order perturbation theory of Ψ_B is open.
- **X_B in Q's.** Writing X_B as an explicit untwisted-Q expression is open.
- **Other reps.** The n_ij rule and the rank count have not been tested for [3,1,0]².

## Section 12 (EH notebook only): site-2 localization and the full reconstruction

Code: `L2-diagonal-frame-overlap-EH.wb`, cells 79–87.
- **What it needs.** Cells 3–14, plus the Section 11 helpers used in 12c (`wCov`, `x0Cov`, `psiStates`, `hwB`, `nexpS`, `EeN`, `rOmF`, `rankTol`, `φT`, `ωMany`).
- **Conventions.** All t's carry the ket twist diag(z), as in Section 10. A covector is site-2 local when ⟨Ω|𝒪 = ⟨Ω_zw| ⊗ ⟨v|, with v ∦ Ω_zw.

### 12. Which t's localize on site 2 (cells 79–82; exact at P1)

**t\* (Paul's guess).**
- **Definition.** The transfer matrices of the inverse monodromy, T⁻¹ = ℒ⁽²⁾⁻¹ℒ⁽¹⁾⁻¹, with the twist diag(z)⁻¹ on the left. Site 2 now sits on the left in auxiliary space.
- **Reduction to ordinary t's.** t\*₁(v) ∝ t₂(v+2h) and t\*₂(v) ∝ t₁(v+2h), up to c-numbers. So searching CBR T_λ and products of t's at θ₂ + ℤh covers the t\*'s.

**CBR scan** (u = θ₂ + ch, c = −5..5, all λ ⊆ 3×3 with λ₃ = 0):
- **Site 2.** The only site-2-local one is **T₍₁,₁₎(θ₂+h) = t₂(θ₂+h) ∝ t\*₁(θ₂−h)**. No T₍₂,₁₎ localizes on site 2 at any integer shift.
- **Site 1 (control).** The same scan at θ₁ + ch returns exactly the Section 10 pair, T₍₁,₁₎ and T₍₂,₁₎ at θ₁−2h.

**Products of two t's at θ₂ + ℤh** (shifts −4..4):
- **The one new local product.** **t₂(θ₂+h)·t₁(θ₂−h) ∝ t\*₁(θ₂−h)·t\*₂(θ₂−3h)**.
- **Relation to site 1.** It plays the role of T₍₂,₁₎(θ₁−2h) = t₂(θ₁−2h)t₁(θ₁−h) − t₃. But the gap between the arguments is −2h rather than +h, so it is not a CBR-fused T\*₍₂,₁₎.
- **Adding t₃.** t₃ is a c-number and Ω_zw⊗Ω_zw is itself site-2 local, so adding multiples of t₃ never changes locality.

**Factorization against the one-site chain at θ₂.** These factors are now **proved symbolically** in Section 13d below. The original one-point fit at h = 1/2 missed a power of h in the scalar term:
- ⟨Ω|t₂(θ₂+h) = (θ₁₂−h)(θ₁₂−3h) · ⟨Ω_zw| ⊗ ⟨Ω_zw|t₂^{L=1,θ₂}(θ₂+h). The factor is the same f₍₁,₁₎ as on site 1.
- ⟨Ω|t₂(θ₂+h)t₁(θ₂−h) = ⟨Ω_zw| ⊗ [ −(θ₁₂−h)²(θ₁₂−3h) · ⟨Ω_zw|t₂^{L=1}(θ₂+h)t₁^{L=1}(θ₂−h) + **6h⁴χ₃(θ₁₂−h)(θ₁₂−3h)** · ⟨Ω_zw| ]. The old coefficient 3h³χ₃f equals the correct 6h⁴χ₃f only at h = 1/2 (apart from zeros).
- Ω_zw, v_A and v_B are linearly independent.

### 12b. Undoing R on site 2 (cells 83–84; exact at P1)
- **The 10b identity holds at θ₂.** ⟨w|R⁻¹ = e₁ᵀ T^{L=1,θ₂}[g⁻ᵀ diag(z) gᵀ] for w_A = Ω_zw t₂^{L=1}(θ₂+h) and for w_B = w_A t₁^{L=1}(θ₂−h).
- **HW reach of e₁ᵀ ⊗ wR⁻¹:**
  - **A** reaches [4,2,0], [4,1,1] and [3,2,1], but **not** [3,3,0] or the singlet. This is the mirror image of site-1 (1,1), which misses [4,1,1] and the singlet.
  - **B** reaches every irrep.
  - The four covectors (site-1 (1,1), site-1 (2,1), A, B) span the 2-dim [3,2,1] HW space.

### 12c. Reconstruction (cells 85–87; redo of 11h)

The site-2 rows are the exact two-site covectors ⟨Ω|t₂(θ₂+h)[t₁(θ₂−h)] at each ε, so the one-point fits above are not used. At ε = 1/3 and 12 random w, these covectors are site-2 local to 10⁻⁷⁴. Rank per weight space, including the x₀ row:

| sector | dim | site1 (= 11h) | site1+Ω | site2 | site2+Ω | site1+site2 |
|---|---|---|---|---|---|---|
| [4,2,0] | 1 | 1 | 1 | 1 | 1 | 1 |
| [4,1,1], [3,3,0] | 2 | 2 | 2 | 2 | 2 | 2 |
| (3,2,1) | 6 | 5 | 6 | 5 | 6 | **6** |
| (2,2,2) | 10 | 7 | 9 | 7 | 9 | **10** |

- **The octet and singlet are recovered.** With site1 + site2 and ⟨x₀|Ψ⟩ = 1, **all six HW-limit states are reconstructed exactly**: both octets and the singlet scale, with no κ.
  - The relative error is 0 except for the singlet: 8×10⁻¹²³ at ε = 10⁻⁶ and 2×10⁻⁷⁷ at ε = 10⁻⁹.
  - The 11h null vector n(ε) is lifted.
- **Side remark.** 11h never used the bare overlap ⟨Ω(w)|Ψ⟩ (the overlap formula itself) as a row. Adding it to site 1 alone already fixes the octet (rank 6), but not the singlet (9 of 10). Site 2 is needed for the singlet scale.
- **What the "data" mean here.** Section 12 uses M·Ψ_B as a stand-in for the Q side. **Section 13 now verifies the full chain using determinant overlaps and transfer eigenvalues recovered from Q functions.**

### Traps (Section 12)
- **meprec on surd covectors.** Exact EH covectors contain surd sums that cancel to 0, and N[...] on them hits meprec. Use `RootReduce` entrywise first (`hwProj`).
- **Missing Section 10 helpers.** The live kernel did not have `loc2Q` or `ΩT`. Section 12 defines its own `site1Q` and `site2Q`.

### Open (Section 12)
- **Closed in Section 13:** symbolic localization factors, exact site-2 HW overlap formulas, the Q-side reconstruction test, and the minimum number of w for the four-row family plus x₀.
- **Other reps.**

## Section 13: independent audit and Q-driven reconstruction

Code: `L2-diagonal-frame-overlap-EH.wb`, cells **90–108**. Rerun **3–14, then 90–108**; no dependency on Sections 5–12. The final audit ledger has **15 passing checks**. This audit was performed after returning to [2,1,0]² at Paul's request.

### Symbolic localization, including a corrected coefficient

Put Δ = θ₁−θ₂ and f = (Δ−h)(Δ−3h). All four identities are polynomial identities, not fits:

| two-site operator | local covector after factoring the spectator |
|---|---|
| T₍₁,₁₎(θ₁−2h) | site 1: f Ω₁ T₍₁,₁₎¹(θ₁−2h) |
| T₍₂,₁₎(θ₁−2h) | site 1: (Δ−h)f Ω₁ T₍₂,₁₎¹(θ₁−2h) |
| A = t₂(θ₂+h) | site 2: f Ω₁ t₂¹(θ₂+h) |
| B = t₂(θ₂+h)t₁(θ₂−h) | site 2: −(Δ−h)f Ω₁ t₂¹(θ₂+h)t₁¹(θ₂−h) + **6h⁴χ₃f Ω₁** |

The site-1 spectator is on the right; the site-2 spectator is on the left. Superscript 1 denotes the one-site chain at the active site's inhomogeneity.

**Proof in the notebook.** Undo R on both sites. The bra becomes e₁ᵀ⊗e₁ᵀ and the twist becomes the lower-triangular matrix g⁻ᵀ diag(z) gᵀ. With all six triangular entries independent and θ=(d,0), h=1, every entry of all four residual matrices simplifies identically to zero. Translation and homogeneity restore arbitrary θ and h. Exact checks at three independent rational chain/twist points also pass. The old β=3h³χ₃f fails at the two points with h≠1/2.

### A simpler determinant normalization

Let D_B(w)=det N[B,vac_w], P_B=Ψ_Q^B(x₀), and B₀=e₁⊗e₁. Then

\[
\boxed{\frac{\langle\Omega(w)|\Psi_B\rangle}{\langle x_0|\Psi_B\rangle}
=\frac{D_B(w)P_{B_0}}{D_{B_0}(w)P_B(z_2-z_1)^2}.}
\]

This follows by dividing the full overlap formula by its B₀ instance, using Ω·B₀=1 and x₀·B₀=(z₂−z₁)². It removes **y₀, Φ_Q and c₀** entirely. The per-ladder and exponential-insertion determinants acquire the same state-independent factor, which cancels in this ratio. Assume generic nonzero denominators; singular parameter choices need limits.

`auditQOverlap` evaluates this ratio. `auditFromQ` recovers t₁,t₂ from the dual Baxter coefficient equations, with the known scalar t₃. The reconstruction RHS is the ratio times these Q-derived operator eigenvalues, followed by the normalization row 1. **It never uses M·Ψ.** The LHS rows use only one-site 8×8 matrices and the proved factors.

### Validation and the untwisted CG limit

At the original finite twist, all 64 states reconstruct from three bra twists at 60-digit working precision, with about 48 digits retained in the resulting vectors. Baseline residuals reported as tracked numerical zero are agreement at finite precision, not exact symbolic equality.

For z=exp(εφ), w=exp(εω), the full Q-driven reconstruction recovers all six HW-limit branches, including both octets and the singlet scale:

| chain | ε | working digits | largest absolute vector error | largest relative vector error |
|---|---|---|---|---|
| P1 | 10⁻⁶ | 240 | 1.67×10⁻¹²¹ | 3.84×10⁻¹⁶² |
| P2 | 10⁻⁶ | 240 | 1.12×10⁻¹²⁰ | 3.48×10⁻¹⁶¹ |
| P1 | 10⁻⁹ | 300 | 5.44×10⁻¹³⁴ | 1.25×10⁻¹⁹² |
| P2, different φ | 10⁻⁹ | 300 | 5.81×10⁻¹³⁴ | 2.57×10⁻¹⁹¹ |
| P1, precision control | 10⁻⁶ | 160 | 3.35×10⁻³⁷ | 7.72×10⁻⁷⁸ |

The 160- and 240-digit reconstructions agree to relative 7.73×10⁻⁷⁸. P1 and P2 are the Section 11 points. The primary φ is (1/3,−2/7,5/11); the alternate is (−2/5,3/8,7/11). The three ω directions are (2/5,1/7,−3/4), (−3/5,2/3,5/7), and (4/7,−5/6,2/9).

Separately, normalized reconstructed vectors converge to the untwisted transfer eigenbasis in the exact HW spaces. At P1 the largest direction error drops from 1.87×10⁻⁷ to 1.87×10⁻¹⁰ as ε drops from 10⁻⁶ to 10⁻⁹. Each non-top branch has error ratio 1000 within 2×10⁻⁷ relatively. The P2 alternate-φ error is below 2.64×10⁻¹⁰. This is **O(ε) convergence**, distinct from the much smaller finite-twist reconstruction error.

**Scope at the end of Section 13.** These were independently validated reconstruction tests, but the Q-polynomials were obtained from Baxter equations using transfer eigenvalues from the test eigensystem. Section 14 below supersedes that input limitation for the tested adjoint-square examples with a numerical polynomial Baxter solver. The eigensystem remains a verification reference for finite-twist eigenvector checks and the untwisted limit comparison. This does not establish general completeness. Leading HW substitution remains invalid in general; reconstruct first, then take ε→0.

### Three bra twists are minimal for this data family

With the four operator rows at each w and one x₀ row, two w values give at most nine equations in the ten-dimensional (2,2,2) weight space. Hence at least three are necessary for unrestricted finite-twist reconstruction. At an exact rational point, three give ranks **1,2,2,6,10**, equal to all five relevant weight-space dimensions. Since the entries are rational functions with fixed algebraic coefficients, a nonzero maximal minor establishes generic sufficiency. The baseline w values are (5/2,7/3,3/5), (11/5,5/3,3/7), and (9/4,7/4,2/5).

This minimality statement is for the specified four rows plus x₀. It does not exclude a smaller family after adding bare Ω, derivatives or other information.

### Exact site-2 HW overlap formulas

Write the conjugated twist as
\[
Z'=\begin{pmatrix}z_a&0&0\\a&z_b&0\\b&c&z_c\end{pmatrix},\qquad K=ac-bz_b.
\]
The exact orthonormal HW basis is `auditHW[mu]` (raising-operator null space, then exact orthogonalization). The table gives the one-site rows after undoing R; spectator factors are not included:

| μ | e₁⊗(Ω₁t₂¹)R⁻¹ | e₁⊗(Ω₁t₂¹t₁¹)R⁻¹ |
|---|---|---|
| [4,2,0] | 3h²zₐzᵦ | 3h³zₐzᵦ(zₐ−z꜀) |
| [4,1,1] | −3h²zₐc/√2 | −3h³zₐ²c/√2 |
| [3,3,0] | 0 | −3h³zₐzᵦa/√2 |
| [3,2,1] | h²K(2√(3/7), −√(3/35)) | h³zₐ(3√(3/7)K, √(3/35)(2ac+5bzᵦ)) |
| [2,2,2] | 0 | 3h³bK/(2√2) |

These are exact finite-twist covector identities. They do not remove the eigenvector-correction terms discussed in Section 11.

## Open

- Other representations, starting with [3,1,0]².
- Whether the e^{2πiu/h} bracket choice generalizes to L sites as L brackets with insertions e^{2πiku/h}, k = 0..L−1.

## Section 14: polynomial Baxter spectral solver and spectral-only reconstruction

Notebook: `L2-diagonal-frame-overlap-EH.wb`, cells 109–114; the notebook was explicitly saved at revision 675 with 114 cells. The final cells 110, 112 and 114 were executed individually after fixes. No new final audit ledger was recorded for Section 14; Section 13's 15/15 checks remain its own ledger.

### Ansatz and solve

`specAnsatz` constructs monic dual polynomials qʲ of degrees μⱼ, along with t₁ and t₂ whose leading coefficients are fixed by the weight and twist, and the known quantum determinant t₃. There are ten unknowns: six Q coefficients and four transfer coefficients. The polynomial dual Baxter identities determine them without physical matrices or eigenvectors as input. In particular,

\[
t_1(u)=\chi_1u^2+[h\sum_i z_i\mu_i-(\theta_1+\theta_2)\chi_1]u+st_1,
\]

and t₂ has leading term χ₂u⁴ and u³ coefficient \(\sum_{i<j}z_iz_j[h(\mu_i+\mu_j-2)-2(\theta_1+\theta_2)]\), with three lower coefficients unknown.

`specSolve` obtains seeds with exact-rational `NSolve`. `specNewton` then performs four scaled, overdetermined Newton/least-squares corrections at working precision `pr+80`, seeded only by those Baxter roots. It rounds to the requested precision and checks the stored numerical coefficients by rationalizing them in the equations. Raising `NSolve` directly to `2pr` was too slow and was aborted. The refinement equations supply the accuracy; precision padding of seeds is not the source of it. No kernel operation remains in flight.

### Reconstruction and validation

`specReconstruct` consumes only the spectral solutions, local 8×8 covectors and determinant overlaps, three bra twists, and x₀ normalization. `specValidate` constructs the physical transfer matrices afterward and checks eigenvector equations coefficientwise with guard precision. Thus physical matrices are a downstream validation tool; a physical `Eigensystem` does not supply the spectral solutions in this finite-twist test.

At h=1, θ=(5/2,0), z=(2,3,5), all 19 weight sectors produce 64 solutions, with counts matching physical weight dimensions, and all 64 reconstruct. A run before the final Newton refinement reported maximum Baxter, physical eigenvector and normalization residuals 3.34×10⁻⁹⁴, 4.02×10⁻⁹² and 9.28×10⁻⁹³ at 100 digits. The final code rerun succeeded, but its exact residuals were not read back; characterize it as roughly 90-digit validation rather than assigning new residual values. The three bra twists are (7/3,11/4,4/7), (9/4,13/5,3/7), and (12/5,17/6,5/8).

### Near-untwisted checks

`specNear` uses exact rational paths z=1+εφ and w=1+εω. These have the same first-order directions as exponential twists, but are not exponential twists. It solves all 21 states in the five HW-weight sectors of dimensions 1,2,2,6,10; reconstruction passes for all, and raising-operator residual selects the six HW branches. Newton refinement fixes a difficult descendant that initially lost about 100 digits through cancellation. The untwisted HW eigensystem in `auditLimitError` is a verification reference only.

| Point and parameters | Max Baxter coefficient residual | Max physical eigenvector residual / norm | Max CG-limit direction error |
|---|---:|---:|---:|
| P1, ε=10⁻⁶, pr=200 | 1.873621426×10⁻¹⁷³ | 1.809978519×10⁻¹³² | 1.864923519×10⁻⁷ |
| P1, ε=10⁻⁹, pr=280 | 4.779807849×10⁻²⁴² | 3.397139631×10⁻¹⁸¹ | 1.864923680×10⁻¹⁰ |
| P2, ε=10⁻⁹, pr=280, alternate φ | 8.105060246×10⁻²⁴⁶ | 2.811004503×10⁻¹⁸¹ | 2.638405385×10⁻¹⁰ |

Here P1=(h,θ)=((1/2),(2/7,−3/5)); P2=((1/3),(1/5,−4/9)). The primary φ is (1/3,−2/7,5/11), alternate φ is (−2/5,3/8,7/11), and the three ω directions are listed in Section 13. The near-zero finite-twist residuals measure Baxter/eigenvector validation; the CG direction errors are O(ε) and should not be conflated with them.

**Scope.** Section 14 demonstrates a numerical polynomial Baxter solver on these finite examples. It is not an analytic closed-form or general completeness proof and does not implement arbitrary representations.
