# `rotations/` — READ THIS FIRST

This folder holds one self-contained thread of the project:

> **Which covectors `⟨0| g`, `g ∈ GL(3)`, can be written as a determinant of integrals of
> Q-functions?**

Everything here is about the gl(3) XXX chain at **L = 2**, site 1 = adjoint `[2,1,0]`, site 2 =
fundamental `[1,0,0]` (24 states).

## Files

| file | what it is |
|---|---|
| `su(3)-XXX-CG-adj-fund.wb` | the notebook. Sections 0–20, cells 1–64. Everything below was computed here. |
| `su(3)-XXX-CG-adj-fund.md` | the detailed reference: conventions, every result with its residual and evidence level, the logic, the open problems, the pitfalls. **The real documentation.** |
| `su(3)-XXX-CG-adj-fund.tex` / `.pdf` | the same content as a short paper, with the equations typeset properly. |
| `img/` | wolfbook eval logs (tooling artefacts, not content). |

**Read order:** this file, then the `.tex`/`.pdf` for the mathematical story, then the `.md`
when you need exact numbers, symbol names, or the list of things that already failed.

## The one-paragraph summary

`⟨0| = e₂₄` is the SoV reference covector; the SoV basis `⟨x|` is normalised so that
`⟨x|Ψ⟩ = Ψ(x)/Ψ(x₀)`. A determinant of Q-integrals, normalised by `Ψ_A(x₀)`, expands as
`Σ_x c_x ⟨x|Ψ_A⟩`, and the covector `⟨V| = Σ_x c_x ⟨x|` is sometimes a rotation `⟨0| g`.
**Which rotation you get is decided by which poles of the measure survive and by the degree in
`q` per site.** Two families are established: the FSoV scalar product against the vacuum at a
second twist `w` gives a 2-parameter family constrained by `b = a(c-a)`, and small
(1×1, 2×2) determinants with a pole-cancelling polynomial give the `E₃₂` rotations `(0,0,c)`.
Nothing built so far leaves the surface `b = a(c-a)`, and §17 now explains part of why: any object
built from one site's `Q`-values can *only* be a local `E₃₂`, and no minor of the FSoV matrix that
contains a `Q₂` row is an SoV expansion at all.

## What is established (all verified in the notebook)

1. **SoV set-up** (exact): companion transfer matrices; `B` built from the **untwisted** `T` has
   spectrum `∏(u - x)` with `x = θ + h(λ₂₁, λ₁₁, λ₂₂-1)` per site, non-degenerate, and
   `⟨0| = e₂₄`. `⟨x|Ψ_A⟩ = Ψ_A(x)/Ψ_A(x₀)` holds for all 24 eigenstates (`10^-73`); the twist
   exponential `z_i^{u/h}` in `Q_i` is required.
2. **FSoV scalar product** (`10^-78`): the 6×6 determinant equals the physical scalar product with
   `e₂₄` on *both* sides; the ket reference was solved for, not assumed.
3. **Rotations from the scalar product**: with the dual Q's of the vacuum at twist `w`,
   `(a,b,c) = (-w₁w₂, w₁, -w₁(w₂+w₃))`, independent of `θ`, hence the surface **`b = a(c-a)`**,
   inverse `w = (b, -a/b, -1/a)`. Nothing deforms off it: not the 225 admissible minors, not the
   `e^{ku}` dressing, not the six permuted vacua. Non-vacuum dual states stop being rotations.
4. **Small determinants give `E₃₂`**: the one-site bracket is a rotation of site 1 only; the
   two-site 2×2 with the right weights is the global `e^{cE₃₂}` with `c = -e^{kh}`. Those weights
   are the residues of `1/(ν₁ν₂)`, which is the ladder measure times `(u-θ₁)`, so the object is
   literally a **minor of the scalar-product matrix with one pole cancelled**.
5. **§16/§19 — per-site `E₃₂` pairs: the full two-parameter family is reachable.** The per-site
   rotation `(e₈e^{c₁E₃₂}) ⊗ (e₃e^{c₂E₃₂})` has SoV support on exactly the four `x₁₁` labels, with
   pattern `[[1, (47/28)c₂], [c₁-(19/28)c₂, c₁c₂]]` and determinant `-19c₂(28c₁-47c₂)/784`
   (`47/28 = (t-2h)/(t-h)`). A *product* of two one-site brackets is rank one there, so it only reaches
   the two branches `c₂ = 0` (the §8 object) and `c₁ = ((t-2h)/(t-h))c₂` — a one-parameter family
   (§16). **But rank one is a property of the product, not of bracket objects**: the `2×2`
   *determinant* with free per-site weights `R, S` and a free insertion `g(u)` has four knobs against
   three conditions, and §19 hits **every** target `(c₁,c₂)` exactly — `(1,2)`, `(2,½)`, `(-3/7,5/2)`,
   and both `(1,0)` and `(0,1)` — twist-independent, `σ₂/σ₁ ~ 10^-94` (genuine products), per-site fit
   residuals `~10^-94`. In particular `1 ⊗ e^{c₂E₃₂}`, which no product can give, *is* a
   Q-determinant (`R = 19/28`, `S = 47/56`, cubic insertion). The apparent site asymmetry of §16 was
   an artefact of the product construction. Per-site rotations are **not** of the form `⟨0|g`, so they
   are neither on nor off the surface — that surface constrains global rotations only.
6. **§17 — a no-go for one-site objects (exact, Gröbner).** An object built from site-1 `Q`-values
   has `c_x` supported on the 8 labels with site 2 at base. The rotations supported there are cut out
   by an ideal whose Gröbner basis is `⟨b, a⟩`: **only `E₃₂`**. For global rotations the basis is
   `⟨c², b, 19a+9c⟩`, i.e. nothing. (Hence §8's global stabiliser 5.) The general site-1-*local*
   rotation needs 12 labels, four of them with site 2 moved — `⟨x|` are not product states.
   A 2×2 Wronskian minor reaches only 3 labels, and no rotation lives in them (`⟨c,b,a⟩`).
7. **§17 — no FSoV minor with a `Q₂` row is an SoV expansion.** All 9315 minors of sizes 2–5 with at
   least one `i=2` row, two-twist criterion: the *best* mismatch is `O(10^-1)`, not `O(10^-3)`, so it
   is not a tolerance artefact. Only the full 6×6 works. This matches the degree count: `Ψ_A(x)` is
   degree 6 (four `Q₁`, two `Q₂`) and a sub-determinant expands into wave functions only if the
   missing factors are common to every label — half a Wronskian is not a common factor. Counts, from
   the cold-rebuilt kernel: nonzero minors 532/540, 880/1800, 9/1800, 64/3375, 0/1800 for the shapes
   (1,1), (2,1), (3,1), (2,2), (3,2); best mismatches `0.204, 0.402, 0.794, 0.403, —`. Two structural
   zeros worth keeping: the size-5 shape (3 `i=1`, 2 `i=2`) vanishes identically for all 1800 minors,
   and (3,1) is nonzero in only 9 of 1800.
8. **§18 — the surface is rigid against free row weights, which reach exactly its tangent.** Write the
   6×6 with free row weights `v^{(m)}_p` (30 parameters, base point = the FSoV measure). Because a
   determinant is multilinear in its rows, the first-order response is *exact*: replacing row `m` by a
   single-pole row **is** `∂D/∂v^{(m)}_p`. Preserving twist-independence (23 conditions) and
   rotation-ness (`δV ∈ span{V₀ρ(X)}`, 20 conditions) leaves a **14-dimensional** allowed space — twice
   the 7 obvious directions (6 row rescalings + the `e^{ku}` dressing), which sit inside it. On all 14,
   `δ[b - a(c-a)] = 0` to 57–65 digits, and the induced `(δa,δb,δc)` has rank exactly **2** = the
   dimension of the surface. So free weights move the rotation *along* `b = a(c-a)`, reaching the same
   2-parameter family as varying the dual twist `w`, and never leave it. The surface is not an artefact
   of the FSoV measure.

9. **§20 — (global $2\times2$ det at `k`) × (one-site bracket at `k'`) is *nearly* `⟨0|g₁g₂`.** Tested on
   all 24 states: the object equals, to a 72-digit fit residual,
   `D = G_c(X₁,X₂)·(1 + c'X₁)` with `G_c = 1 + (9c/28)X₁ + (47c/28)X₂ + c²X₁X₂`, `c = -e^{kh}`,
   `c' = -e^{k'h}/z₁`, where `X₁ = Q₁(θ₁+2h)/Q₁(θ₁+h)`, `X₂ = Q₁(θ₂+h)/Q₁(θ₂)`. The target is
   `P = 1 + (9c/28 + c')X₁ + (47c/28)X₂ + c(c+c')X₁X₂`, and

   ```
   D - P = (c c'/28) · X₁ · (9X₁ + 19X₂ + 28c X₁X₂).
   ```

   The constant, `X₁` and `X₂` coefficients are **exactly** right; `X₁X₂` overshoots by `(19/28)cc'`,
   and `X₁²`, `X₁²X₂` appear with coefficients `(9/28)cc'`, `c²c'`. `X₁²` is not a wave-function ratio
   (site 1 has only `x₁₁ ∈ {θ₁+h, θ₁+2h}`, so `Ψ_A(x)/Ψ_A(x₀)` is degree ≤ 1 in `X₁`): the product puts
   `Q₁` at a site-1 `x₁₁` point **twice**, once from the determinant and once from the bracket, and
   `1/(Q₁(θ₁+h)²Q₁(θ₂))` is the right normalisation for the product but cannot cure the extra power.
   Measured: rank`{c(z),c(z')} = 2` (not even an SoV expansion), 14 labels used, stabiliser 1,
   `O(1)` deviation after the best constant. **The minimal correction is exact:**
   `⟨0|g₁g₂ ↔ G_c + c'X₁(1 + cX₂)` — the one-site dressing must multiply only the *site-2 leg*, not the
   whole global object. Since `1 + cX₂` is itself a site-2 bracket at the same `c`, the right object is
   a **sum** (global + `c'`×site-1 shift×site-2 bracket), i.e. an antisymmetric combination — which is
   what §19's determinant with a tuned insertion produces.

## The mechanism (the thing worth remembering)

- A polynomial insertion `P(u)` cancels poles of the ladder measure. A site with **no** surviving
  pole is untouched (local rotation); a site keeping its two `x₁₁` points contributes `E₃₂`.
- The degree in `q` per site decides the rest: `q₁` alone can only move `x₁₁`. Moving `x₂₁`/`x₂₂`
  needs `q₂` too, because they sit inside the Wronskian `Q₁(x₂₁)Q₂(x₂₂) - Q₂(x₂₁)Q₁(x₂₂)`.
- **Degree counting is a hard constraint** (§17): a determinant of size `n` is a degree-`n`
  polynomial in Q-values, and the labels it can reach are those whose `Ψ_A(x)` share a common factor
  of degree `6-n`. Freezing a coordinate is what produces such a factor; cutting a Wronskian is not.
- Being a rotation is a **pattern condition** on the `c_x`, whose constants are pole ratios; the
  weight-independent part of it is a **cross-ratio** of the evaluation points.
- This is the rational-measure version of Paul's periodic `ρ_α` in the thesis
  (`Paul/Papers/Ryan_Thesis_2201.12057.pdf`, eq. (12.31)), which kills the poles of the other
  sites so that `⟨·⟩_α` sees only `θ_α`. Bracket and measure: (12.27)–(12.28).

## What to do next

The `E₂₁`/`E₃₁` directions — leaving `b = a(c-a)` — remain the open problem, and §17–§18 closed the
obvious routes. What survives:

1. **Per-site `E₂₁`/`E₃₁` (the natural next target, and the most promising).** §19 got *independent*
   per-site `E₃₂` rotations out of a `2×2` determinant the moment the weights and the insertion were
   both left free. The same freedom has never been tried on an object with `Q₂` rows, which is what
   moving `x₂₁`, `x₂₂` requires. §17's no-gos do **not** cover this: they rule out objects built from a
   single site's `Q`-values, and all minors of the FSoV matrix, but not a two-site free-weight
   Wronskian object.
2. **Finite (not infinitesimal) weight deformations of the `4×4` shape (3 `i=1`, 1 `i=2`).**
   `⟨0|e^{aE₂₁}` is exactly *linear* in `a` (higher powers annihilate `⟨0|`), with `c_{x₀}=1`,
   `c_{(0,1,0|0,0,0)}=-20a/9`, `c_{(1,0,0|0,0,0)}=-9a/28`, `c_{(0,0,0|1,0,0)}=893a/252` and **no
   doubly-moved label**. Dividing out the frozen site-2 Wronskian, the object must be the degree-4
   expression written out in the §17 result cell. §18 shows first-order rigidity, which does not forbid
   a large-weight solution but makes one unlikely: 16 parameters against 46 conditions, so the
   structure has to conspire — a search alone will not do it. **Agree the ansatz with Paul first.**
3. **Vary the dual state and the weights together.** §10 killed non-vacuum dual `Q`'s and §18 killed
   infinitesimal weight deformations, but only one knob at a time.
4. **Change the reference covector.** Everything is tied to `⟨0| = e₂₄` being the lowest-weight
   covector of this particular `B`. The Weyl conjugates of §11 gave the same surface; a different `B`
   (a different SoV basis) would move the question itself.
5. **Derive** `(a,b,c)(w)` and `b = a(c-a)` (verified, not derived). §18 makes this more pressing: the
   relation now looks like a property of the construction rather than of the choice of measure.
6. **Check the `θ`-dependence.** Every constant in §16–§17 is a pole ratio (`47/28 = (t-2h)/(t-h)`
   etc.). Whether the reachable *locus* is genuinely `θ`-dependent, or only its parametrisation,
   decides whether these constraints survive the general-rep limit.

## How to re-run

The notebook is self-contained but **not** a single-shot run: code cells are
`2,3,5,7,8,10,12,14,16,19,23,26,29,30,33,36,37,40,43,46,49,52,53,56,59,62,63`; run them in order (markdown cells in
between). **All of cells 1–64 were cold-run verified on 2026-09-18** (restart, then 1–16, 17–33, 34–46, 47–64 in
order, no kernel messages); the whole rebuild takes well under a minute of kernel time, cell 26 (the
225-minor survey) being the slowest at ~7 s and cell 53 (9315 minors) ~2 s.

Parameters: `θ = (1/3, -8/57)`, `h = 1`, chain twist `z = (2,3,1/6)`, control twist
`z' = (3/2,5/4,8/15)`, second twist `w = (5/2,7/3,6/35)`; 80-digit numerics where not exact.

## Traps (each of these cost time already)

- **The twist-independence criterion is the right notion of "expands into wave functions".** A
  formal expansion in abstract `Q` does **not** exist for the 6×6; an on-shell expansion always
  exists and so tests nothing. Solve `c_x` at two twists `z`, `z'` and demand agreement.
- **Per-site shifts**: `q₁(u + hλ₂^{(s)})`, i.e. `q₁(u+h)` on the adjoint site and `q₁(u)` on the
  fundamental. Using `q₁(u+h)` at both samples a non-SoV point and the object is not an SoV
  expansion at all.
- **The insertion must be evaluated at the same points as `q`**, or the cross-ratio invariant is
  wrong and no weights or `k` can repair it.
- **`⟨x|` are not product states.** Do not test "is it a per-site rotation" by slicing the
  24-vector into an 8×3 matrix and reading off a row and a column — that is wrong and produced a
  false claim once. Use the line-stabiliser dimension (6 ⟺ rotation of `e₂₄`), the rank of the 8×3
  matrix (rank 1 ⟺ the covector really is a product) and the exact coefficient patterns. §16's
  non-product rows still return per-site fit residuals of `10^-95`; they mean nothing until
  `σ₂/σ₁ ≈ 0`.
- **`GroebnerBasis`, not `Solve`.** The "which rotations are supported in this label set" systems
  are small but `Solve` runs for tens of minutes and ignores `Abort`; `GroebnerBasis` answers them
  instantly. Read `⟨c², b, 19a+9c⟩` as the radical `⟨c,b,a⟩`.
- **Precision 0 kills survey loops.** An identically vanishing minor is pure cancellation, so at 80
  digits its determinant has precision 0 and *every later comparison* triggers `$MaxExtraPrecision`
  retries: the same survey took 45+ min at 80 digits and 2 s at machine precision. Filter at machine
  precision with a loose tolerance, then re-test the survivors at 80 digits.
- **`AppendTo` in a scan loop is quadratic** — use `Reap`/`Sow`.
- **Do not trust an inherited warm kernel.** The kernel handed to the 2026-09-18 session held caches
  from two epochs (the earlier session had restarted mid-way, per the note at the end of this list), so
  `vacAtTwist[w]["qd"]` was inconsistent with the cached results of §7–§11: re-running the §10 object
  gave a covector 0.1% off `⟨Ω^L_w|` at 71-digit precision, and re-deriving `vacAtTwist` outright
  failed (`NullSpace` empty). A cold restart plus a full re-run reproduced every stored result exactly,
  so nothing in §1–§15 was wrong — but the symptom is silent and looks like physics. Cold-run before
  trusting a base point, and prefer recomputation to a cache whenever a result will be built on.
- Wolfbook replies truncate at `{`; read results as strings with braces replaced.
- The notebook was moved into `rotations/` after VS Code had it open from `Experiments/`, so the
  editor can hold a stale path: `wolfbook_saveNotebook` fails with `ENOENT` unless you pass the
  `notebook` argument with the full `rotations/` path.
