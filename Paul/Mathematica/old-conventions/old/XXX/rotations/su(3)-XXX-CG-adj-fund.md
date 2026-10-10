# `su(3)-XXX-CG-adj-fund.wb` — companion reference

gl(3) XXX chain, **L = 2**, site 1 in the adjoint `[2,1,0]`, site 2 in the fundamental `[1,0,0]`
(24 states). Builds the companion-twisted transfer matrices, the SoV basis, the FSoV scalar
product, and then studies the central question of the project in this setting:

> **Which covectors of the form `<0| g`, `g ∈ GL(3)`, can be written as a determinant of
> integrals of Q-functions?**

All objects considered turn out to be instances of **one** construction — a minor of the
scalar-product matrix, with a polynomial inserted in the integrand to cancel some poles of the
measure. This file records what is established, the logic behind it, and what is open.

- **Notebook:** `su(3)-XXX-CG-adj-fund.wb` (sections 0–15, cells 1–47).
- **Relatives:** `SU3_L2_Baxter_Twisted_Final.wb` (adjoint², diagonal twist, 64×64 match),
  `SU3_L1_GeneralRep_Normalization.wb` (general-rep `Ψ^ref`, `c^ref`),
  `SU3_310x210_L2_ZeroTwist.wb`, `SU3_100x100_L2_ZeroTwist.wb` (mixed / fundamental measures).
- **Status:** cells 1–16 cold-run verified; later sections run in a warm kernel after them.
  Exact rational arithmetic where stated, otherwise 80 digits.

---

## 0. Conventions

| symbol | meaning |
|---|---|
| `θ = (1/3, -8/57)`, `h = 1` | chain point (`prmX`); `t := θ₂ - θ₁ = -9/19` |
| `λ1 = {2,1,0}`, `λ2 = {1,0,0}` | site reps; `d1 = 8`, `d2 = 3`, `dtot = 24` |
| `L_ij(u) = u δ_ij - h E_ji` | Lax; `T_ij(u) = Σ_k L^{(1)}_ik(u-θ₁) L^{(2)}_kj(u-θ₂)` |
| `G = {{χ1,-χ2,χ3},{1,0,0},{0,1,0}}` | companion twist, eigenvalues `z`, `det G = χ3` |
| `τ_a(u) = Σ_{|I|=a} qmin^G_{I,I}(u)` | twisted transfer matrices (Cauchy–Binet from `qmin`) |
| `ν_i(u) = ∏_s (u - θ_s - h λ^{(s)}_i)` | Yangian weights: `ν1 = (u-θ₁-2h)(u-θ₂-h)`, `ν2 = (u-θ₁-h)(u-θ₂)`, `ν3 = (u-θ₁)(u-θ₂)` |
| `z = (2,3,1/6)`, `z' = (3/2,5/4,8/15)` | two unimodular chain twists (the second only as a test) |
| `Q_i = z_i^{u/h} q_i`, `Q^j = z_j^{-u/h} q^j` | direct / dual Baxter solutions, minimal degree |
| `B(u)` | SoV operator, built from the **untwisted** `T` |
| `x^α_{kj} = θ_α + h(μ^α_{kj} - j + 1)` | SoV roots; per site `(x₁₁,x₂₁,x₂₂) = θ + h(λ₂₁, λ₁₁, λ₂₂-1)` |
| `⟨0| = e₂₄` | lowest-weight covector = minimal-`x` eigencovector of `B` |
| `Ψ_A(x) = ∏_s Q₁(x₁₁)·[Q₁(x₂₁)Q₂(x₂₂) - Q₂(x₂₁)Q₁(x₂₂)]` | SoV wave function |
| `⟨x|Ψ_A⟩ = Ψ_A(x)/Ψ_A(x₀)` | normalisation of the SoV covectors (`x₀` = minimal label) |

**Rotation parametrisation.** `⟨0|E_ij ≠ 0` exactly for `i ≥ j`, so the open Bruhat cell is

```
<V_abc| = <0| e^{a E21} e^{b E31} e^{c E32}
```

with `(a,b,c)` readable off any covector: `V₂₄ = 1` and three components are `a`, `b/2`, `c/2`.
Substituting back and comparing all 24 components is the **rotation test**; an equivalent
invariant test is that the line stabiliser `{X ∈ gl(3) : ⟨V|ρ(X) ∈ C⟨V|}` has **dimension 6**
(controls: 6 for `e₂₄`, 1 for a random covector; 3 for a random adjoint-site covector).

---

## 1. Established results

### 1.1 Chain, B operator, SoV basis (§1–5, exact)

- Companion `τ_1, τ_2` commute; `τ_3 = χ3 ν1^{[-4]}ν2^{[-2]}ν3`; spectrum equals the diagonal-twist
  spectrum and is simple (24 distinct eigenvalues).
- `B(u) = T_{11}(u)·qmin_{12,12}(u+h) + T_{21}(u)·qmin_{12,13}(u+h)` from the **untwisted** `T` has
  degree 6, leading coefficient 1, commuting coefficients, and spectrum `∏(u - x)` with the root
  formula above; non-degenerate. `e₂₄·B = b₀(u) e₂₄` holds symbolically, so `⟨0| = e₂₄`.
- `|Ω⟩` (companion `τ₁` eigenvalue `Σ z_i ν_i`, normalised `⟨0|Ω⟩ = 1`) has `q₁ = 1`, `deg q₂ = 2`,
  `deg q₃ = 3`.
- Normalising `⟨x|` on `Ω` alone, `⟨x|Ψ_A⟩ = Ψ_A(x)/Ψ_A(x₀)` holds for **all 24 eigenstates**
  (max deviation at the 80-digit floor, ~`10^-73`). **The twist exponential is required**: dropping
  `z_i^{x/h}` gives O(1) deviations.

### 1.2 FSoV scalar product (§6, 80 digits)

Rows `(α,i,j)`, `α ∈ {1,2}` (second carries `e^{2πiu/h}`), `(i,j) ∈ {(1,3),(1,2),(2,3)}`;
columns `u^β(-Q^jQ_i, Q^jQ_i^{[-2]}, Q^{j[2]}Q_i)`, `β = 0,1`; **ladder measure**, poles
`θ₁+{0,h,2h}`, `θ₂+{0,h}` (5 simple poles, weights `1/∏_{p'≠p}(p-p')`).

- Orthogonal over all 24×24 pairs (`9×10^-72`).
- `Ψ^ref`, `c^ref` are the general-rep ones with each site's own `λ₂` (adjoint 1, fundamental 0);
  `Ψ^ref_B = Ψ_B(x₀)`, consistent with `⟨0| = e₂₄`.
- The **implied ket reference is also `e₂₄`** (solved for, not assumed).
- `C_vac det N^{AB}/(c^ref_A Ψ^ref_B) = ⟨L_A|R_B⟩/(⟨L_A|e₂₄⟩⟨e₂₄|R_B⟩)` to `2×10^-78`, with
  `1/C_vac = S_vec(Ω,Ω) = 22627/9`.

### 1.3 What "expands into wave functions" has to mean (§7)

- A **formal** expansion in abstract `Q_A` values does **not** exist: the abstract determinant has
  242 monomials, only 96 of them occur in the 24 `Ψ_A(x)`, and those carry ~3% of the weight. The
  rest live at non-SoV points and presumably reduce only through `A`'s Baxter equation.
- An **on-shell** expansion always exists (24 states, 24 labels), so by itself it tests nothing.
- **Working criterion:** solve `c_x` on-shell at two independent chain twists `z`, `z'` and demand
  they agree (up to one overall scale — a convention artefact of the exponentials). Controls: single
  `Q₁²Q₂`-monomials fail this (rank 2). The normalised `⟨x|` are identical at `z` and `z'`.

### 1.4 The rotations reachable from the scalar product (§7–11)

Dual Q's of the **vacuum at a second twist `w`** (unimodular), abstract `A`:

```
(a, b, c) = ( -w₁w₂ , w₁ , -w₁(w₂+w₃) ),      independent of θ
```

- exact on 6 twists; the covector is `⟨Ω^L_w|` (left companion vacuum at `w`);
- eliminating `w` leaves **one relation**

```
                     b = a (c - a)            (equivalently c = a + b/a)
```

so this is a **codimension-1 surface** in the 3-dimensional rotation group;
- inverse map `w = (b, -a/b, -1/a)`, automatically unimodular (exact on 4 targets);
- **nothing deforms off it:** all 225 admissible 6×6 minors give the *same* covector (span
  dimension 1); the `e^{ku}` dressing acts as `(a,b,c) → (a e^k, b e^{2k}, c e^k)`, which preserves
  the relation (Jacobian rank 2, determinant tracking the finite-difference step); the six vacua
  `τ₁ = Σ_i w_{σ(i)}ν_i` give the permuted formula, hence the same surface; dual Q's of non-vacuum
  states are not rotations at all (stabiliser 3 or 5).
- Pure `E₃₂` rotations `(0,0,c)` satisfy the relation; pure `E₂₁` rotations `(a,0,0)` do **not** —
  that is the direction to aim for.

### 1.5 Small objects: one site (§8) and two sites (§13–15)

**One site (degree 1 in `q`).** `⟨e^{ku}q₁(u+h)⟩/q₁(θ₁+h)` with measure `1/((u-θ₁)(u-θ₁-h))`
expands **formally** (unlike the 6×6) into two labels and gives a covector that is
`s_k ⊗ e₃` with `s_k ∝ e₈ exp(a(k) ρ(E₃₂))`, `a(k) = -½ e^{kh}`: a rotation **on site 1 only**.
Site-1 stabiliser 6, global stabiliser 5.

**Two sites (degree 2 in `q`).** With brackets at each site's own `x₁₁` points,

```
<f>_1 = a₁ f(θ₁+h) + b₁ f(θ₁+2h),      <f>_2 = a₂ f(θ₂) + b₂ f(θ₂+h)
D     = det [ <e^{ku} q₁>_s , <e^{ku} u q₁>_s ] / ( q₁(θ₁+h) q₁(θ₂) )
```

is a **global** `E₃₂` rotation iff

```
b₁/a₁ = c e^{-kh} · t/(t-2h),        b₂/a₂ = c e^{-kh} · (t-2h)/t,        t = θ₂-θ₁
```

(verified for `c = 1, -3/7, 2` and `k = 0, 1/2`: stabiliser 6, `(a,b,c) = (0,0,c)`, residual
`~10^-94`; equal weights give stabiliser 5). Only `c e^{-kh}` enters, so `k` is redundant.

**These weights are residues.** `ν₁ν₂` vanishes exactly at the four `x₁₁` points, and

```
<f>_s = ∮_{C_s} du/(2πi) · e^{ku} f(u) / (ν₁(u) ν₂(u))        ⇒     c = -e^{kh}
```

(`C_s` around that site's two poles). Verified at `k = 0, 1/2, -1, iπ/h`; the last gives `c = +1`.

**And this is a minor of the scalar product.** `1/(ν₁ν₂) = (u-θ₁) × (ladder measure)`, and the
rows `α = 1,2` (full contours, weights `(1,1)` and `(e^{2πiθ₁/h}, e^{2πiθ₂/h})`) span the same
2-plane as the two per-site contours. For the twist-`w` vacuum `deg q³ = 0`, so `Q³ = w₃^{-u/h}`
*is* the `e^{ku}` insertion. Hence the object is the minor on rows `{(1,1,3),(2,1,3)}` and columns
`{(r=1,β=0),(r=1,β=1)}`, with the integrand dressed by `(u-θ₁)`. Checked: same `(0,0,-e^{kh})`, and

```
det(minor, α rows) / det(site rows) = e^{2πiθ₂/h} - e^{2πiθ₁/h}
```

state-independent to `10^-73`.

**The one-site object is the same scheme.** Relative to the ladder it cancels `θ₁` **and both**
`θ₂` poles — factor `(u-θ₁)(u-θ₂)(u-θ₂-h) = ν₃(u)·(u-θ₂-h)` — leaving only `θ₁+h, θ₁+2h`. No
surviving pole at site 2, so site 2 is untouched: the rotation is local.

---

## 2. The logic

Read every object as

```
        minor of the scalar-product matrix   ×   polynomial insertion P(u) in the integrand
```

and three things follow.

**(i) `P` selects the poles, the poles select the coordinates.** The ladder has 5 poles. Each factor
of `P` cancels one. Which poles survive at a site decides what that site contributes:

| surviving poles at a site | that site's contribution |
|---|---|
| none | untouched (identity) |
| its two `x₁₁` points | `x₁₁` moves ⇒ `E₃₂`-type |
| points inside the Wronskian block | `x₂₁`, `x₂₂` move ⇒ `E₂₁`, `E₃₁` directions |

**(ii) The degree in `q` per site is the other half.** `Ψ_A(x)` factorises per site into
`Q₁(x₁₁)` and the Wronskian `Q₁(x₂₁)Q₂(x₂₂) - Q₂(x₂₁)Q₁(x₂₂)`. An object built from `q₁` alone can
only move `x₁₁`, whatever the measure — which is why every `q₁`-only object we built landed on
`E₃₂`. Moving `x₂₁`/`x₂₂` **requires `q₂` rows**, i.e. a Wronskian-shaped minor.

**(iii) Being a rotation is a pattern condition on the coefficients.** With `V = Σ_x c_x ⟨x|`, the
group element is visible in the exact coefficient patterns, e.g. on the four `x₁₁` labels

```
global e^{cE32}:        [[1, (47/28)c], [(9/28)c, c²]]
per-site e^{c1E32}⊗e^{c2E32}: [[1, (47/28)c₂], [c₁-(19/28)c₂, c₁c₂]]
```

and the constants are **pole ratios**: `47/28 = (t-2h)/(t-h)`, `9/28 = t/(t-h)`. The
weight-independent invariant is a cross-ratio,

```
c₁₁c₂₂/(c₁₂c₂₁) = cross-ratio of the four evaluation points = t(t-2h)/(t-h)² = 784/423
```

which the `E₃₂` pattern requires. This is **automatic** when the insertion and `q₁` are evaluated
at the same four points (§13), and **fails** when they are not (§12). The remaining conditions then
fix the free weights/insertions and the rotation parameter.

**Corollary (why `k` keeps reappearing as a scale).** Inserting `e^{ku}` multiplies the residue at
`p` by `e^{kp}`; ratios between points one step apart pick up `e^{kh}`, which is why every parameter
found so far is `∝ e^{kh}` (`a(k) = -½e^{kh}`, `c = -e^{kh}`, and the surface's Cartan rescaling
`(a,b,c) → (ae^k, be^{2k}, ce^k)`).

---

## 3. What is open

**The goal.** Generate more rotations — ideally *local* ones, independent on the two sites — and in
particular reach `b ≠ a(c-a)`, i.e. the `E₂₁`/`E₃₁` directions, which no object built so far
reaches.

1. **Wronskian minors — partly settled, see §6 below.** Minors of the FSoV matrix with `i = 2` rows
   are *never* SoV expansions (all 9315 of sizes 2–5 fail the two-twist test), and any object built
   from one site's `Q`-values can only be a local `E₃₂`. What is left is the **free-weight** `4×4`
   object of shape (3 rows `i=1`, 1 row `i=2`) with poles surviving at both sites, matched against
   the explicit degree-4 `E₂₁` target in §6.
2. **Local `E₃₂` rotations, cheaply — settled, see §5 below.** The product of two one-site brackets
   gives `e^{c₁E₃₂} ⊗ e^{c₂E₃₂}` only on the one-parameter locus `c₁ = ((t-2h)/(t-h))c₂`; independent
   parameters are impossible. The one-site bracket on the *fundamental* site is not a rotation at
   all.
3. **A classification.** For each subset of surviving poles and each minor shape, which covector
   patterns are attainable? The cross-ratio condition suggests the attainable set is cut out by
   pole cross-ratios, i.e. by `θ` — worth making precise, since it would say exactly which
   rotations are Q-determinants and which are not.
4. **Closed form for `(a,b,c)(w)`'s derivation.** The formula `(-w₁w₂, w₁, -w₁(w₂+w₃))` is verified
   but not derived; deriving it would probably also explain the surface relation `b = a(c-a)` as a
   statement about the companion twist's flag.
5. **Then the actual target.** Rotations are the tool, not the goal: the project wants the
   localisation/CG construction for general (non-rectangular) gl(3) reps. The relevant lesson so far
   is that the reachable rotations are constrained by the pole structure of the measure, which is
   itself rep-dependent (`Ψ^ref` uses each site's `λ₂`).

---

## 4. Pitfalls hit (do not repeat)

- **Per-site shifts matter.** `q₁(u+h)` at both sites is wrong; the fundamental site needs
  `q₁(u+h λ₂) = q₁(u)`. With the wrong shift the object samples a non-SoV point and is not an SoV
  expansion at all (twist-dependent, stabiliser 1).
- **The insertion must sit at the same points as `q`.** Otherwise the cross-ratio invariant is wrong
  and no weights or `k` can repair it (§12 needed a tuned, complex second column; §13 with the right
  points needs only real weights).
- **`⟨x|` are not product states.** Do not test "is it a product / is it a per-site rotation" by
  slicing the 24-vector as an 8×3 matrix and reading off a row and a column — that was wrong once in
  this session. Test with the stabiliser dimension and the exact coefficient patterns.
- **A formal expansion of the 6×6 does not exist** — see §1.3; use the two-twist criterion.
- **Rank-2 coefficient matrices are not evidence against a rotation**, because the `⟨x|` basis is not
  a product basis; a genuine `E₃₂` rotation also has a rank-2 matrix on those four labels.
- **`GroebnerBasis`, not `Solve`.** The "which rotations are supported inside this label set" systems
  (§6) are small, but `Solve` runs for tens of minutes on them and does not respond to `Abort`;
  `GroebnerBasis` returns instantly. Read the non-radical answer `⟨c², b, 19a+9c⟩` as `⟨c,b,a⟩`.
- **Precision 0 destroys survey loops.** An identically vanishing minor is pure cancellation, so at
  80 digits its determinant carries precision 0 and every later comparison triggers
  `$MaxExtraPrecision` retries. The §6 survey took 45+ min at 80 digits and 2 s at machine
  precision. Filter at machine precision with a deliberately loose tolerance, then re-test the
  survivors at 80 digits (a tight machine tolerance would risk false negatives, since these
  determinants are badly conditioned).
- **`AppendTo` inside a scan loop is quadratic** — use `Reap`/`Sow`.
- **Stale notebook path.** The notebook was moved into `rotations/` while VS Code had it open from
  `Experiments/`; `wolfbook_saveNotebook` then fails with `ENOENT` unless the full `rotations/` path
  is passed explicitly as the `notebook` argument.
- **Never trust an inherited warm kernel (this one bit hard).** The kernel inherited on 2026-09-18
  held caches from two epochs — the earlier session had restarted mid-way — so
  `vacAtTwist[w]["qd"]` no longer matched the cached results of §7–§11. Symptom: re-running the §10
  object gave a covector 0.1% away from `⟨Ω^L_w|` **at 71-digit precision** (so not noise), the
  rotation fit residual read `10^-2`, and recomputing `vacAtTwist` from scratch failed outright
  (`NullSpace` empty). A cold restart and a full re-run of cells 1–57 reproduced every stored result
  exactly — nothing in §1–§15 was wrong, and §16/§17's exact results were bit-identical — but two
  minor counts in §6.4 changed. The failure mode is silent and looks like physics. Cold-run before
  building on a base point.
- Wolfbook replies truncate at `{`; read results as strings with braces replaced. One cell hung the
  kernel (a contrast block, not the physics) — restart and rebuild rather than debugging in place.

---

## 5. Section 16 — products of one-site brackets (cells 48–50)

**The exact target.** The per-site rotation covector `(e₈e^{c₁E₃₂}) ⊗ (e₃e^{c₂E₃₂})` (built with
`JJ[3,2,λ_s]` per site, so it is a genuine product covector) has SoV support on **exactly** the four
`x₁₁` labels `{20,21,23,24}` = `fourLabels`, and there

```
c = [[ 1              , (47/28) c₂ ],
     [ c₁ - (19/28)c₂ , c₁ c₂      ]],        det c = -19 c₂ (28c₁ - 47c₂)/784
```

exact in `θ`; `47/28 = (t-2h)/(t-h)`, `9/28 = t/(t-h)`, `t = θ₂-θ₁`.

**The object.** Each one-site bracket is `B_s = Q₁(p_{s,1}) + r_s Q₁(p_{s,2})`, normalised by
`Q₁(p_{s,1})`, on that site's two `x₁₁` points (`ptsSite`). The ratio `r_s` is free: a two-pole
measure `1/((u-p_{s,1})(u-p_{s,2}))` dressed by `e^{k_su}` gives `r_s = -e^{k_sh}`, and complex `k_s`
covers every nonzero value. The product `B₁B₂` has the **rank-one** coefficient matrix
`c_{ij} = r_1^{i-1} r_2^{j-1}`, i.e. `[[1, r₂], [r₁, r₁r₂]]`.

**Result (superseded in part by §8 below: free weights plus a free insertion remove the constraint).**
Rank one forces `det c = 0`, whose two branches are exactly

| branch | weights | rotation |
|---|---|---|
| `c₂ = 0` | `r₁ = c₁`, `r₂ = 0` | `e^{c₁E₃₂} ⊗ 1`, `c₁` free — this is §8 |
| `c₁ = (47/28)c₂` | `r₁ = c₂`, `r₂ = c₁` | `e^{c₁E₃₂} ⊗ e^{c₂E₃₂}` with `c₁/c₂ = (t-2h)/(t-h)` fixed |

so the old conjecture (independent `c₁`, `c₂`) is **false**; what exists is a one-parameter family
whose ratio is a pole ratio. In terms of dressings, the locus is
`k₂ - k₁ = log[(t-2h)/(t-h)]/h = log(47/28) ≈ 0.5178` and then `(c₁,c₂) = (-e^{k₂h}, -e^{k₁h})`.

**Evidence.** `rank{c(z),c(z')} = 1` (twist-independent), `σ₂/σ₁` of the 8×3 matrix `~10^-96`
(the covector really is a product), per-site `fitRep` residuals `~10^-94`. Checked at
`(r₁,r₂) = (1/2, 47/56) → (c₁,c₂) = (47/56, 1/2)`, `(-3/5, -141/140) → (-141/140, -3/5)`, and via
dressings at `k₁ = 0` and `k₁ = 1/2`.

**Two negatives.**

1. Equal dressings `k₁ = k₂`, and in particular the natural residue ratios of `1/(ν₁ν₂)`
   (`{r₁,r₂} = {-9/47, -47/9}`), give an object that is neither a product nor a rotation
   (`σ₂/σ₁ = 4×10^-1`, global fit residual `O(1)`). With those same weights the 2×2 **determinant**
   of §13–§14 is the global `e^{cE₃₂}`: product and determinant of the same two brackets are
   different covectors, and only the determinant is a global rotation.
2. The one-site bracket on the **fundamental** site (`r₁ = 0`, `r₂ ≠ 0`) is not a rotation of any
   kind and not even a product (`σ₂/σ₁ = 1.9×10^-2`, global residual `4.8×10^-2`): its pattern
   `[[1,r₂],[0,0]]` would need `c₁ - (19/28)c₂ = 0` and `c₁c₂ = 0`, hence `c₂ = 0`. The adjoint
   site's `x₁₁` shift alone is a rotation; the fundamental site's is not.

*Trap.* For the non-product rows the per-site slice fits still return residuals `~10^-95`. They are
meaningless until `σ₂/σ₁ ≈ 0`; the verdict column in cell 49 applies the product test first.

---

## 6. Section 17 — Wronskian minors and two no-go results (cells 51–54)

### 6.1 Structure of what a bracket object can produce

A `2×2` minor with one `Q₁` row and one `Q₂` row, columns supplying the shifts `σ = 0, -h`, equals
`Σ_{p,p'} v_p v'_{p'} W(p+σ₁, p'+σ₂)` — the diagonal terms cancel by antisymmetry of `W`, so it is
automatically a sum of Wronskians. With the site-2 poles cancelled and only `p ∈ {θ₁, θ₁+h}`
surviving it reaches `x₂₁ ∈ {θ₁, θ₁+h}`, `x₂₂ ∈ {θ₁-h, θ₁}`, i.e. the three labels with `x₁₁`
unmoved (`threeLab` = `{24, 18, 15}`). Keeping `θ₁+2h` would produce `λ₁₁ = 2`, illegal unless `x₁₁`
also moves — which needs a `Q₁` row supplying `Q₁(x₁₁)`, i.e. a `3×3` minor of shape `(1,1,2)`.

### 6.2 No-go: an object built from one site's `Q`-values is at best a local `E₃₂` (exact)

Such an object has `c_x` supported on the **8 labels with site 2 at base** (`site1Block`), because the
site-2 factor of `Ψ_A(x)` cancels against the normalisation. Which rotations are supported inside a
label set is decided by the ideal generated by the coefficients *outside* it:

| rotations supported inside | Gröbner basis | conclusion |
|---|---|---|
| site-1-local, the 8 site-2-base labels | `⟨b, a⟩` | `a = b = 0`: **only `E₃₂`** |
| site-1-local, the 3 Wronskian labels | `⟨c, b, a⟩` | nothing |
| global `⟨0|g`, the 8 site-2-base labels | `⟨c², b, 19a+9c⟩` (radical `⟨c,b,a⟩`) | nothing |
| global `⟨0|g`, the 3 Wronskian labels | `⟨c, b, a⟩` | nothing |

Consequences: **no local `E₂₁`/`E₃₁` from a site-1-only object**, whatever the measure or minor shape;
and §8's object has global stabiliser 5 rather than 6 because no global rotation lives in that block
either. The general site-1-*local* rotation needs **12** labels, four of which (`8,11,20,23`) have
site 2 moved — a rotation acting as the identity on site 2 still needs site-2-shifted labels, because
`⟨x|` are not product states.

The site-1-local generator coefficients (cell 52, exact) coincide with the global ones at linear
order for `aa`, since `⟨e₃|E₂₁ = 0` on the fundamental site.

### 6.3 The `E₂₁` target, explicitly

`⟨0|e^{aE₂₁}` is **exactly linear** in `a` (`⟨0|E₂₁² = 0`):

```
c_{x₀} = 1,   c_{(0,1,0|0,0,0)} = -20a/9,   c_{(1,0,0|0,0,0)} = -9a/28,   c_{(0,0,0|1,0,0)} = 893a/252
```

and **no doubly-moved label**. The four labels share the factor `W₂(θ₂,θ₂-h)`; dividing it out, the
object must be the homogeneous degree-4 expression

```
Q₁(θ₁+h) W₁(θ₁,θ₁-h) Q₁(θ₂)
  + a [ -(9/28) Q₁(θ₁+2h) W₁(θ₁,θ₁-h) Q₁(θ₂)
        -(20/9)  Q₁(θ₁+h)  W₁(θ₁+h,θ₁-h) Q₁(θ₂)
        +(893/252) Q₁(θ₁+h) W₁(θ₁,θ₁-h) Q₁(θ₂+h) ]
```

— three `Q₁`'s and one `Q₂`, i.e. a `4×4` minor of shape (3 rows `i=1`, 1 row `i=2`) with site 2's
Wronskian frozen and surviving poles at **both** sites. The absence of cross terms rules out a
*product* of independent brackets (that would carry `a²`, `a³`), so the cancellation has to come from
the antisymmetry of a genuine determinant.

### 6.4 No-go: no FSoV minor with a `Q₂` row is an SoV expansion (numerical, complete)

All minors of the `6×6` scalar-product matrix containing at least one `i=2` row, sizes 2–5, against
the twist-`w` vacuum dual `Q`'s, two-twist criterion:

| shape (#`i=1`, #`i=2` rows) | minors | not identically zero | best twist mismatch |
|---|---|---|---|
| (1,1) | 540 | 532 | `2.04×10^-1` |
| (2,1) | 1800 | 880 | `4.02×10^-1` |
| (3,1) | 1800 | 9 | `7.94×10^-1` |
| (2,2) | 3375 | 64 | `4.03×10^-1` |
| (3,2) | 1800 | 0 | — |

(Counts and mismatches from the cold-rebuilt kernel of 2026-09-18. An earlier run in the inherited
warm kernel gave 884 and 73 nonzero for (2,1) and (2,2) — see the last pitfall below; the conclusion
was unaffected.)

The *best* mismatch over 9315 minors is `O(10^-1)`, so this is not a tolerance artefact: the machine
filter accepts anything below `10^-3` and re-tests it at 80 digits, and nothing survived. **Only the
full (4,2) minor — the whole `6×6` — is an SoV expansion.** This matches the degree count: `Ψ_A(x)` is
degree 6 (four `Q₁`, two `Q₂`), and a sub-determinant expands into wave functions only when the
missing factors are common to every label. Freezing a coordinate produces such a factor; cutting a
Wronskian in half does not. Two structural zeros worth keeping: shape (3,2) vanishes identically for
all 1800 minors, and (3,1) is nonzero in only 9 of 1800.

### 6.5 What is left after §7 below

Not a minor of the FSoV matrix. The remaining freedom is the one §8 and §13 already needed for
`E₃₂`: free weights. §7 below tests that freedom at first order and finds the surface rigid, so what
survives is listed at the end of §7.

---

## 7. Section 18 — free row weights: the exact linear response of the surface (cells 55–57)

### 7.1 Set-up

Write the `6×6` with **free row weights**, `M_{mc} = Σ_p v^{(m)}_p × [column c at pole p]`, the base
point being the FSoV choice `v^{(m)}_p = w_p e^{2πi(α_m-1)p/h}` which gives the rotation `⟨Ω^L_w|` on
the surface. The first-order response is **exact**, not a finite difference: a determinant is
multilinear in its rows, so replacing row `m`'s weights by a single-pole row *is* `∂D/∂v^{(m)}_p`
(`dDof` in cell 56). Both properties to be preserved are linear in the deformation `δ ∈ C^30`:

1. **Still an SoV expansion.** With `c^{(z)} = c₀ + εA` and `c^{(z')} = c₀/κ + εB`, proportionality at
   first order is `A - κB ∈ C c₀` — 23 conditions.
2. **Still a rotation.** `δV ∈ W := span{V₀ ρ(X) : X ∈ gl(3)}`, which is 4-dimensional (three orbit
   directions plus `V₀`) — 20 conditions.

Seven deformations are known to be allowed and harmless: rescaling each of the 6 rows, and the
`e^{ku}` dressing `v_p → v_p e^{kp}` on every row, which acts as `(a,b,c) → (ae^k, be^{2k}, ce^k)` and
multiplies `b - a(c-a)` by `e^{2kh}`.

### 7.2 Result

Base point (cold kernel): twist residual `0` to 69 digits, rotation fit residual `0` to 71 digits,
`(a,b,c) = (-35/6, 5/2, -263/42)` exactly `= (-w₁w₂, w₁, -w₁(w₂+w₃))`, surface value `0` to 70 digits.
The 7 known-harmless directions satisfy all 44 conditions to `10^-68`, which validates the
construction.

| quantity | value |
|---|---|
| conditions × directions | 44 × 30 |
| `dim W` | 4 (as expected) |
| **dim of the allowed deformation space** | **14** |
| rank {allowed space, 7 known directions} | 14 — the known ones lie inside, so 7 directions are new |
| `max |δ[b-a(c-a)]| / |δ(a,b,c)|` over the 14 | `0` to 57–65 digits |
| rank of the induced `(δa,δb,δc)` | **2** |
| rank `{δV} ∪ {V₀}` | 3, i.e. a 2-parameter projective family |

So free weights are a genuine extra freedom — 14 allowed directions, not 7 — but the extra freedom
moves the rotation **along** the surface. The rank-2 image is exactly the tangent of `b = a(c-a)`:
the weights reach the same two-parameter family that varying the dual twist `w` reaches, and nothing
more. This extends §1.4's "nothing deforms off it" (225 minors, `e^{ku}` dressing, six permuted
vacua) to the whole 30-parameter weight family at first order, and shows the surface is not an
artefact of the FSoV measure.

### 7.3 What survives for `E₂₁`/`E₃₁`

1. **Finite weight deformations** of the `4×4` shape (3 rows `i=1`, 1 row `i=2`) against the explicit
   degree-4 target of §6.3. First-order rigidity does not forbid a large-weight solution but makes one
   unlikely: 16 parameters against 46 conditions, so the structure must conspire. Since the
   determinant is multilinear, the *linear* span of such objects is all of `C^24` and proves nothing;
   the content is the factorisation (rank-one) condition. **Agree the ansatz with Paul before
   grinding on it.**
2. **Vary the dual state and the weights together.** §10 killed non-vacuum dual `Q`'s with FSoV
   weights; §18 killed infinitesimal weight deformations with a vacuum dual. Never both at once.
3. **Change the reference covector**, i.e. a different `B` and hence a different SoV basis. The Weyl
   conjugates of §11 stay on the same surface, but they keep the same `B`.
4. **Derive** `(a,b,c)(w)` and `b = a(c-a)`. §18 makes this more pressing: the relation now looks like
   a property of the construction, not of the measure.

---

## 8. Section 19 — independent rotations on the two sites (cells 58–60)

§5's constraint `c₁/c₂ = (t-2h)/(t-h)` came from the coefficient matrix of a *product* being rank
one. A genuine `2×2` **determinant** with free weights and a free insertion is not rank-one limited.
With site `s` carrying weights `v^{(s)}` on its two `x₁₁` points and the second column carrying `g(u)`,

```
D = det[[<Q1>_1, <g Q1>_1], [<Q1>_2, <g Q1>_2]] = Σ_{i,j} v^{(1)}_i v^{(2)}_j (g(p_{2,j}) - g(p_{1,i})) Q1(p_{1,i}) Q1(p_{2,j})
```

so with `c₁₁ = 1` and the affine normalisation `g(p_{1,1}) = 0`, `g(p_{2,1}) = 1`,

```
c₁₂ = S δ,     c₂₁ = R (1 - γ),     c₂₂ = R S (δ - γ),
R = v^{(1)}₂/v^{(1)}₁,  S = v^{(2)}₂/v^{(2)}₁,  γ = g(p_{1,2}),  δ = g(p_{2,2}).
```

Four unknowns against the three ratios of the target pattern, so `γ` stays free and each target has a
one-parameter family of realisations. Since the per-site rotation's support is *exactly* those four
labels (§5), matching the three ratios matches the whole 24-component covector.

**Result: every target is hit exactly.**

| target `(c₁,c₂)` | `γ` | `{R,S}` | rank`{c(z),c(z')}` | `σ₂/σ₁` | achieved | per-site residual |
|---|---|---|---|---|---|---|
| `(1,2)` | 2 | `{5/14, -157/140}` | 1 | `7e-94` | `(1,2)` | `2e-94` |
| `(2,1/2)` | 2 | — | 1 | `9e-95` | `(2,0.5)` | `1e-94` |
| `(-3/7,5/2)` | `-1` | — | 1 | `3e-95` | `(-0.4285714,2.5)` | `3e-94` |
| `(47/28,1)` (on §5's locus) | 2 | `{-1, 47/28}` | 1 | `8e-94` | `(1.6785714,1)` | — |
| `(1,0)` (site 1 only) | 2 | `{-1, 0}` | 1 | `1e-94` | `(1,0)` | — |
| `(0,1)` (site 2 only) | 2 | `{19/28, 47/56}` | 1 | `7e-94` | `(0,1)` | `1e-94` |

All six are twist-independent, genuine product covectors, and fit a per-site rotation in both reps.
All six have global stabiliser 5 and global rotation residual `O(1)` — correctly: a per-site rotation
with `c₁ ≠ c₂` is not of the form `⟨0|g`, so it is neither on nor off the surface `b = a(c-a)`, which
is a statement about global rotations only.

Two corrections to §5:

1. The constraint `c₁/c₂ = (t-2h)/(t-h)` belongs to the product construction, not to bracket objects.
2. `1 ⊗ e^{c₂E₃₂}` — which no product can give — *is* a Q-determinant, with `R = 19/28`,
   `S = 47/56` and the cubic insertion
   `g(u) = 23675/11421 + (31375/5076)u - (59147/5922)u² + (37183/11844)u³`. The site asymmetry of §5
   was an artefact of the product too.

**Conventions.** As in §1.5's two-site object, each site's bracket is its own two-pole contour
integral, and `g(u)` sits in the integrand at the same four points as `Q₁`, so §4's "same points"
requirement is respected.

**Open for per-site rotations:** the `E₂₁`/`E₃₁` directions at a single site, which need `x₂₁`, `x₂₂`
to move and hence `Q₂` rows. §6.2's no-go covers only objects built from one site's `Q`-values, and
§6.4 covers only minors of the FSoV matrix; a two-site free-weight Wronskian object is untested, and
free weights are exactly what unlocked this section.

---

## 9. Section 20 — (global 2x2 determinant) x (one-site bracket): nearly, but not quite (cells 61–64)

Paul's construction: multiply the §1.5/§6 global object `detCont(k)` (which gives `e^{cE32}` with
`c = -e^{kh}`) by the §8 one-site bracket at a *different* parameter `k'`, normalising by
`Q1(θ1+h)^2 Q1(θ2)`. Claim: the result is `⟨0|g1 g2|Ψ⟩` up to a state-independent constant.

Write `X1 = Q1(θ1+2h)/Q1(θ1+h)` and `X2 = Q1(θ2+h)/Q1(θ2)`; on the four `x11` labels
`Ψ_A(x)/Ψ_A(x0) = X1^{i-1} X2^{j-1}`, so any covector supported there is of degree <= 1 in each.

**What the object is.** A least-squares fit over the 24 states against the six monomials `X1^a X2^b`,
`a <= 2`, `b <= 1`, leaves residual `0` to **72 digits** at every `(k,k')` tried, so exactly

```
D(k,k') = G_c(X1,X2) (1 + c' X1),    G_c = 1 + (9c/28) X1 + (47c/28) X2 + c^2 X1 X2,
c = -e^{kh},   c' = -e^{k'h}/z1
```

(the `z1` in `c'` because §8's bracket is written in `q1`, not `Q1`; measured `c' = -1/2` at `k' = 0`
with `z1 = 2`, `-e^{1/2}/2` at `k' = 1/2`, `-e^{-1}/2` at `k' = -1`).

**What the target is.** `E32` at the two sites commutes, so `g1 g2 = e^{(c+c')E32} ⊗ e^{c E32}` and by
§5's exact pattern

```
P = 1 + (9c/28 + c') X1 + (47c/28) X2 + c (c + c') X1 X2.
```

**The difference, exactly:**

```
D - P = (c c' / 28) X1 (9 X1 + 19 X2 + 28 c X1 X2).
```

The constant, `X1` and `X2` coefficients agree to 72 digits — the construction gets those right. The
`X1X2` coefficient overshoots by `(19/28) c c'`, and two illegitimate monomials appear: `X1^2` with
`(9/28) c c'` and `X1^2 X2` with `c^2 c'`.

**Why.** The product carries `Q1` at a site-1 `x11` point **twice** (once from the determinant, once
from the bracket), whereas `Ψ_A` has exactly one `Q1(x11)` per site. `1/(Q1(θ1+h)^2 Q1(θ2))` is the
correct normalisation for the product, but no normalisation can remove a surplus power of `Q1` at one
site. Consequences, all measured on the 24 states: `rank{c(z),c(z')} = 2` (not an SoV expansion at
all), 14 of 24 labels used, line-stabiliser dimension 1, and an `O(1)` (2.4–13) deviation from
`⟨0|g1g2|Ψ_A⟩` after fitting the best state-independent constant.

**The minimal correction** (difference verified identically zero):

```
⟨0|g1 g2   <->   G_c + c' X1 (1 + c X2).
```

The one-site dressing must multiply only the **site-2 leg** `(1 + c X2)`, not the whole global object;
multiplying the whole object is what generates `X1^2`. Since `1 + c X2` is itself a one-site bracket on
site 2 at the same `c`, the correct object is `global + c' x (site-1 x11 shift) x (site-2 bracket)` — a
**sum** of two products, i.e. an antisymmetric combination, which is what a determinant gives for free.
That is exactly §8's object: there the tuned insertion `g(u)` generates the mixed term without ever
putting two `Q1`'s on the same site. So independent per-site rotations *are* Q-determinants; this
particular product is not one.
