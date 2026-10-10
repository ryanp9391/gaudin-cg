# `SU3_L2_Untwisting.wb` — companion reference

gl(3) XXX chain, **L = 2**, both sites in the adjoint `[2,1,0]`. Takes the exact twisted identity of
`SU3_L2_Baxter_Twisted_Final.wb` (normalized physical scalar product = normalized FSoV determinant)
to zero twist. Along `z_j = e^{εφ_j}`, `Σφ_j = 0`, the notebook establishes:

1. **Divergence.** How `⟨0|R⟩` and `⟨L|0⟩` vanish as ε → 0, and why.
2. **Twist factors.** Every twist factor in the identity is explicit.
3. **Untwisted limit.** The limit of the norm is written purely through **untwisted Q-functions**
   (Appendix N of arXiv:2603.26617, generalized).
4. **Bare eigenvectors.** How that object compares with bare untwisted eigenvectors
   (normalization as a function of `t = θ1 − θ2`).
5. **Separate limits.** The normalized right and left eigenvectors each have an untwisted limit.

- **Write-up:** `SU3_L2_Untwisting.tex` / `.pdf` (same directory). It has the fully explicit formula.
- **Twisted parent:** `SU3_L2_Baxter_Twisted_Final.wb` / `.md` (all-64 match, 36-triple invariance).
- **Older untwisted/twisted L=2 notebooks:** `SU3_Adjoint_SoV_ScalarProduct_L2.wb` (3×3 bare-row norm
  `N̂`), `SU3_Adjoint_Twisted_SoV_L2.wb` (column basis H/O, block dominance).
- **Status.** Every number below comes from a run of the cells in order, in a live kernel. There has been
  **no single cold top-to-bottom run** yet (see §0.3). Numerics: 200 digits.

---

## 0. Read this first

1. **Unimodular twist, `z_j = e^{εφ_j}`.** Three directions are used, `φ ∈ {(1,1/3,−4/3), (2,−1/2,−3/2),
   (1/2,1,−3/2)}` (`phiSetU`), and ε ∈ {1e-3,…,1e-8}. Limits use quadratic extrapolation through
   ε = 1e-6, 1e-7, 1e-8 (error ~1e-18).
2. **Chain point.** P1 is `θ = (1/3, −8/57)`, `h = 1` (hard-coded in shifts). Here `t = θ1−θ2 = 9/19`,
   where the octets are rational. The independent check point is P2: `θ = (1/5, 46/55)`.
3. **Cells 13, 16, 17 and 19 overwrite global chain state.**
   - Cell 13 rebuilds everything at P2.
   - `buildAtU[t]` (cells 16–19) rebuilds `qcoU` and `untwU` at `θ = (t, 0)`. It does **not** rebuild
     `polesU`/`weightsU`, so the functions `matUntwU`/`untwHU` are stale after it. `BexplU` has its own
     bracket and is safe.
   - Twisted data are memoized per `(ε, φ)` (`twU`, `qSideU`, `matRegU`). Cell 13 therefore uses a
     fresh `φ` so that nothing is reused across chain points.
   - Run the cells in order. A cold run needs ~10–20 min, dominated by the twisted solves in cells 5–7
     and 13.
4. **`qpolyU` is Horner** (`Fold`). The earlier `cl.y^Range[…]` gave `0^0 = Indeterminate` when a pole
   sits at u = 0 (θ2 = 0). The cell-3 source is fixed.
5. **Wolfbook replies truncate at `{`.** Read results by converting `ToString[..., InputForm]` with
   braces replaced.
6. **Kernel.** K2 is shared with other notebooks in the directory; restart before a clean run.

## 1. Conventions

| symbol | meaning |
|---|---|
| `Q_i = z_i^u q_i`, `Q^j = z_j^{−u} q^j` | twisted direct/dual; raw `q_i` monic of degree `4−μ_i`, `q^j` of degree `μ_j` |
| untwisted `q_i`, `q^j` | staggered: `deg q_i = (4−λ_i)+(i−1)`, `deg q^j = λ_j+(3−j)` (`untwU[[m,"q"/"qd"]]`) |
| `S(A,B)` | `⟨L_A|R_B⟩/(⟨L_A|0⟩⟨0|R_B⟩)`, references = product of one-site branch-I states, `a13=3/5`, `a23=4/7` |
| `⟨f⟩_α` | `Σ_{6 poles θ_s+n, n=0,1,2} Res M(u) e^{2πi(α−1)u} f(u)`, `M = 1/(ν1ν2ν3)` |
| rows | `(α,i,j)`, α∈{1,2}, `(i,j) ∈ tripleU = {(1,3),(1,2),(2,3)}` |
| columns `c1..c6` | `(−Q^jQ_i, Q^jQ_i^{[−2]}, Q^{j[2]}Q_i)` and the same times `u` |
| `H` | `(c1+c3, c5−c6+5c3, c2−c3)`: integral-of-motion directions surviving at zero twist |
| `O` | `(c4, c6+c3, c3)`: the three `O(φ)` directions |
| `V_λ(z)` | `∏_{i<j}(z_i−z_j)^{λ_i−λ_j}`; `V_λ ≈ ε^{Σ} ∏ φ_ij^{λ_i−λ_j}`, `φ_ij = φ_i−φ_j` |
| states `m = 1..6` | 27 `(4,2,0)`, 10 `(3,3,0)`, 10̄ `(4,1,1)`, 8₁, 8₂ `(3,2,1)`, 1 `(2,2,2)`; octets ordered by `[u⁰]T₂` (8₁ lower) |
| `S⊗S` | `S = diag(1,1,1,2,6,9,4,9)`: `E_ji = S^{−1}E_ij^T S` in the GT basis |

**Why this triple.** `Q_1` and `Q^3` stay finite as ε → 0. `Q_2 → q_2`, `Q^2 → q^2` need the moves
`Q_2 + xQ_1` and `Q^2 + yQ^3`, which only add an existing row (exact invariance). The row-choice
invariance of `SU3_L2_Baxter_Twisted_Final.wb` guarantees the result does not depend on the triple.

## 2. Cell-by-cell symbol reference

### Cell 2 — generators and chain
`ps, dim, J0, Jp, Jm, JJ` (GT); `repU, dsU=8, dcU=64, EEU, E1U, E2U, EtotU`; `thetaU`;
`TmonU`, `qminorU`, `qcoU[I,k]` (exact coefficients of principal quantum minors);
`TcoU[a,k,z]` (twisted transfer-matrix coefficients for any z); `nu1U, nu2U, nu3U, qdetU`.
Checks: commutativity, `T₃ = qdet` (unimodular), gl(3) invariance at z = 1.

### Cell 3 — untwisted highest-weight states
`BaxDirU[T1,T2,z]`, `BaxDualU[T1,T2,y]` (reduced Baxter operators; `z=1` untwisted), `qpolyU` (Horner),
`polyKernelU` (exact polynomial kernel, staggering gauge), `wtU`, `hwMusU`, `hwBasisU`,
`untwU` (6 records `mu, vec, T1, T2, q, qd`), `psiRefLocU`, `cRefLocU` (normalization factors for the
triple content `{a,a,b}`, `{c,c,d}`).

### Cell 4 — twisted solver
`precU = 200`, `phiU`, `zOf[ε,φ]`, `solveMonicU`, `refBra1U[z]`, `refKet1U[z]` (branch I),
`multipletOf` (nearest untwisted `(T1,T2)`), `twU[ε,φ]` (memoized; all 64 records with `vec, left,
q, qd, mult, hwLimit, braR = ⟨0|R⟩, Lket = ⟨L|0⟩, eigRes, bRes`).

### Cell 5 — divergence (step 1)
`epsListU`, `keyOf`, `slopeU`, `divTab`, `limitCheck`.

### Cell 6 — Q side of the identity
`qRU, qLU`, `polesU, weightsU`, `rowU[A,B,z,{α,{i,j}}]`, `tripleU`, `detNU`, `psiRefU`, `cRefU`,
`hwRecU[ε,(φ),m]`, `qSideU[ε,φ]` (`Cvac` and per-state `detN, cRef, psiRef, Svec`), `powU`, `qTabU`.

### Cell 7 — ε → 0 limits
`epsLimU`, `phiSetU`, `kU = {0,1,1,2,2,4}`, `limU`, `limDataU[φ]` (`X = lim ε^{2k}S_vec`,
`D = lim ε^{−(12−2k)}det N`, `P`, `K`, `Cv`), `untwFactorsU` (`Ψ^ref`, `c^ref` untwisted).

### Cell 8 — twist structure
`lamU`, `nU = {3,1,4,2,2,1}` (= λ₁−λ₂+1), `ndU = {3,4,1,2,2,1}` (= λ₂−λ₃+1), `vandPhiU`, `DhatU`, `D0U`,
`kappaU`, `psi0NormU`, `rowUntwU`, `matUntwU` (untwisted 6×6), `structTab`, `rankTab`.

### Cell 11 — block powers in the H/O basis
`colChangeU` (det = −1), `regQU` (2-point matching regularization), `rowRegU`, `matRegU[ε,φ,m]`
(rows × (H|O)), `entryPowU`, `blockPowTab`.

### Cell 12 — precise relation at P1
`splitResU`, `untwHU[m]` (untwisted periodic minor), `EformulaU[λ]`, `EnumU`, `relationTab`.

### Cell 13 — independent check at P2
Rebuilds the chain at `θ = (1/5,46/55)`; `phiP2 = (3/2,1/2,−2)`, `epsP2`, `checkP2`.

### Cell 16 — explicit Q-formula at `θ = (t,0)`
`Sdiag`, `S2U`, `buildAtU[t]`, `bracketExplU`, `BexplU[s]` (= `E_λ det𝒩/(c^ref Ψ^ref)`, the tex formula),
`tValsU`, `dataAtU`.

### Cell 17 — t-dependence and required normalization
`tMoreU` (16 rational t with rational octets), `dataMoreU`, `ratFitU`, `fitEvenU`, `octDataU`,
`octFitU`, `octFitResU`, `normTabU`.

### Cell 19 — separate right/left limits
`betaRLU`, `betaTabU`, `aRU`, `aLU`, `monoU`, `compU`, `bConstU`, `splitTabU`.

## 3. Results

### 3.1 Divergence (cells 5–6)
- `⟨0|R_A⟩ ~ ε^{k}` and `⟨L_A|0⟩ ~ ε^{k}`, with `k = ht(λ_vac − λ_A) = 0,1,1,2,2,4`.
- Twisted eigenvectors converge at `O(ε)` to the untwisted hw vectors.
- **Mechanism.** The reference covector stays finite as ε → 0 and lies in the 27*, so it annihilates
  every other multiplet at zero twist.
- **Exact identity** at every ε: `|S_Q/S_vec − 1| ≲ 1e-90`.
- **ε-powers:** `C_vac ~ ε^{−8}`; `c^ref, Ψ^ref ~ ε²`; `det N ~ ε^{12−2k}`; `S_vec ~ ε^{−2k}`.

### 3.2 Twist structure of every factor (cells 7–8; residuals ≤ 1e-18)
- `Ψ^ref → (z1−z2)² Ψ^ref_untw/n²` and `c^ref → (z2−z3)² c^ref_untw/n'²`, with `n = λ₁−λ₂+1`,
  `n' = λ₂−λ₃+1`.
- `det N_reg → V_{λ_A}(z) · D̂⁰_A`, with `D̂⁰_A` φ-independent.
- `C_vac → κ/V_{λ_vac}(z)`, with `κ = (cΨ)^untw_vac/D̂⁰_vac`.
- Hence `lim V_{λvac−λA}(z) S_vec(A,A) = [D̂⁰_A/(cΨ)^untw_A]/[D̂⁰_vac/(cΨ)^untw_vac]`.
- The untwisted 6×6 (all 18 rows too) has rank 3,4,4,4,4,6. The deficiency sits in the **columns**,
  which have state-dependent rational null vectors. Only for the singlet is `D̂⁰` the untwisted 6×6
  determinant itself.

### 3.3 The precise relation (cells 11–13) — Appendix N mechanism
- **Bare rows are at infinity.** The α=1 rows vanish. The power of row `(i,j)` in column O1 is
  `4 − deg q_i − deg q^j` (untwisted staggered degrees), exactly as residue at infinity predicts. In
  H1..3, O1, O2, O3 the profile is `(n+2,n+2,n+2,n,n,n+1)`.
- **Laplace split.** `det N_reg = −det[A1]_O · det[A2]_H · (1+O(ε))` for all 6 states.
  - Residual 1e-6.5 at ε=1e-7 and 1e-7.5 at ε=1e-8 (O(ε)).
  - `O1+O2` cancels at leading order, so no other Laplace term competes.
- **Periodic block.** `det[A2]_H → det[A2]^untw_H`, the untwisted 3×3 minor (residual 1e-7 at ε=1e-8).
- **Bare block.**
  `det[A1]_O → V_λ(z) E_λ`, with `E_λ = −(λ1−λ2+1)(λ2−λ3+1)/[(λ1−λ3+1)!(λ1−λ2)!(λ2−λ3)!]`.
  - Values: −3/160, −1/36, −1/36, −2/3, −2/3, −1.
  - The factorials are `(5−d_ij)!` over the three rows. The numerator `n n'` is the product of the leading
    coefficients of the untwisted Casoratians `q1∧q2` and `q^3∧q^2`.
- **Evidence levels.**
  - **P1:** `E` agrees with its numerical value to 1e-18; the boxed ratio matches the physical limit to 1e-19.
  - **P2** (independent chain point and φ): `E_λ` formula to 1e-13–1e-14 for all six states; split residual
    ~1e-7.
  - **Status:** `E_λ` is identified at P1 and independently validated at P2, but **not derived**.
- **Boxed result.**
  ```
  lim_{ε→0} V_{λvac}(z)/V_{λA}(z) · S(A,A) = B_{λ_A}/B_{λ_vac},    B_λ = E_λ det𝒩/(c^ref Ψ^ref)
  ```
  `𝒩` is the 3×3 matrix of **untwisted** periodic brackets, rows (1,3),(1,2),(2,3), with columns
  - `⟨q^j(u+1)q_i − q^j q_i⟩`
  - `⟨u q^j q_i(u−1) − (u−5) q^j(u+1) q_i⟩`
  - `⟨q^j q_i(u−1) − q^j(u+1) q_i⟩`

  This explicit form was checked against `untwHU` for all 6 states (exact to working precision).
- `det[A2]^untw_H` is **not** proportional to the older bare-row 3×3 `N̂`.

### 3.4 Q-expression versus bare untwisted eigenvectors (cells 16–17; exact, 16 rational t)
- **Left vectors.** `l_0 = (S⊗S) r_0`; `⟨l_0|r_0⟩ = ‖r_0‖²_S`.
- **`B_A/B_vac` as a function of `t = θ1−θ2`** (h=1), even in t:
  - 27: `1`
  - 10, 10̄: `8/(1−t²)`
  - 1: `73728/((1−t²)²(9−t²)²)`
  - 8±: `[3(4+5t²)(5t²−10t+9)(5t²+10t+9) ± 6(125t⁴−350t²+81)√(4+5t²)] / [t²(1−t²)²(9−t²)²]`,
    with 8₁ = −, 8₂ = +. This matches exact values, e.g. t=9/19 gives `17511884375/163924992` and
    `5603803/43200`.
- **At P1** (t=9/19): `B/B_vac = 1, 361/35, 361/35, 17511884375/163924992, 5603803/43200,
  16983563041/10672200` (exact).
- **Normalization.** Fix `r_0` by the component `|hw⟩⊗|extremal⟩` = 1 (27: |1,1⟩; 10: |1,3⟩; 10̄: |1,2⟩;
  1: |1,8⟩). Then `B/B_vac = ‖r_0‖²_S/(N_R N_L)`, with `N_R N_L = 1, (1−t²)/4, (1−t²)/4,
  [(1−t²)(9−t²)/32]²`. These are products of `(k±t)`, k odd, of degree 2k.
- **Octets:** no canonical component (the site-2 weight is zero, a two-dimensional space). With |2,3⟩
  the required normalization does not factorize. **Open.**

### 3.5 Separate untwisted limits of `|R_A⟩` and `⟨L_A|` (cell 19; P1, three φ)
- `⟨0|R⟩/ε^k → β_R(φ) = b_R ∏φ_ij^{aR}` and `⟨L|0⟩/ε^k → β_L(φ) = b_L ∏φ_ij^{aL}`. Both are monomials;
  `b_R`, `b_L` are φ-independent to ~1e-20.
- `aR + aL = ((λvac−λA)_i − (λvac−λA)_j)` for all states. Hence
  `|R_untw⟩ = lim ∏(z_i−z_j)^{aR}|R_A⟩ = r_0/b_R` and `⟨L_untw| = lim ∏(z_i−z_j)^{aL}⟨L_A| = (S r_0)^T/b_L`,
  and **the boxed object equals `⟨L_untw|R_untw⟩`**.
- Exponents on (φ12, φ13, φ23):

  | state | aR | aL | b_R, b_L (component of r₀ = 1) |
  |---|---|---|---|
  | 27 | (0,0,0) | (0,0,0) | 1, 1 |
  | 10 | (1,1,−1) | (1,0,0) | 21/76, 10/57 (|1,3⟩) |
  | 10̄ | (0,0,1) | (−1,1,1) | −20/133, −49/152 (|1,2⟩) |
  | 8, 8 | (1,1,0) | (0,1,1) | numerical (|2,3⟩); rationals unreliable |
  | 1 | (2,2,0) | (0,2,2) | numerical (|1,8⟩); rationals unreliable |

- **Caveats.**
  - The left/right split and `b_R`, `b_L` very likely depend on the reference parameters (a13, a23, branch);
    only the product is fixed. The negative exponents are probably an artifact of that choice.
  - Checked at one chain point only.

## 4. Pitfalls hit (cost real time)
- **Output truncation.** Wolfbook output stops at the first `{`; stringify with braces replaced.
- **`0^0` in `qpolyU`** when a pole is at u = 0 (fixed with Horner).
- **`KroneckerProduct` of vectors** is an outer matrix (flatten); inherited from the parent notebook.
- **Misread rank.** I first misread the octet untwisted rank as 5; it is 4 (singular values 1e-2 vs 1e-78).
- **Relative cutoffs vs working precision.** Always use working precision, not a fixed relative cutoff.
- **Toy model.** Untwisted polynomial Q's × exact exponentials give the right vanishing order, but a
  **wrong, φ-dependent** coefficient. The ε-corrections of the polynomial parts contribute at leading order.
  Do not use the toy.
- **Octet fits** need the √(4+5t²) branch explicitly; a plain rational fit in t fails.
- **Large-denominator `Rationalize` output** (octets/singlet at 1e-14) is not trustworthy; use exact
  evaluations (`BexplU`) instead.
- **`buildAtU` does not rebuild `polesU`/`weightsU`** (see §0.3).

## 5. Open
1. **Derive `E_λ` analytically** from the residue at infinity. `O1+O2` cancels at leading order, so the
   next order of the large-u expansion is needed, including the two-exponential structure of the
   regularized `Q_2`, `Q^2`.
2. **Natural octet normalization** of bare eigenvectors. Test: is `N_R N_L` an eigenvalue of the
   adjoint⊗adjoint R-matrix `Ř(θ1−θ2)`, or `⟨l_0|Ř|r_0⟩`? The octet block of Ř carries the same
   √(4+5t²).
3. **Closed form of `b_R`, `b_L`** and their dependence on `a13`, `a23` and branch. Does a symmetric
   split exist?
4. **Cold top-to-bottom run** of the notebook; promotion to `Clean/` after that.
5. **Relation of `det[A2]^untw_H`** (periodic rows) to the older bare-row `N̂`, and whether other
   triples give the same untwisted minor up to explicit factors.
