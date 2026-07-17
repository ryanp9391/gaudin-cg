# su(3) general gl(3) reps per site — design

**Date:** 2026-07-17
**Subproject:** XXX (gl(3))
**Notebook:** `Paul/Mathematica/XXX/Experiments/su3_V2_general.wb` (currently a verbatim copy of
`su3_V1.wb`; generalize in place)
**Status:** design approved 2026-07-17

## Goal

Generalize the two-site gl(3) XXX construction of `su3_V1.wb` (both Section A — reps/Lax/twist/
transfer matrices — and Section B — eigensystem + Baxter Q-functions) from **symmetric reps
`{S,0,0}` per site** to **arbitrary gl(3) reps labelled by Young diagrams** `λ = {λ1,λ2,λ3}` with
`λ1 ≥ λ2 ≥ λ3 ≥ 0`, all integers, one independent diagram per site.

Everything V1 established stays: diagonal-twist frame `tg`, `(H1,H2,n)` state labelling, numeric
spectrum with symbolic-in-`u` eigenvalues, the three Baxter Q-functions `Q_k = z_k^{u/h} q_k(u)`,
and the same hard-assert / zero-residual verification convention. The companion-twist SoV covectors
(the CG-overlap program) remain out of scope, as in V1.

## Fixed decisions (confirmed with the user)

- **Interface = λ-vectors everywhere.** Every function takes full 3-vectors `λ1,λ2` (one per site);
  the symmetric case is just `{S,0,0}`. No scalar-`S` backward-compat wrappers — every Section A/B
  cell is re-signatured `S1,S2 → λ1,λ2`.
- **Degree formula = re-derive the closed form.** Extract `M_k` numerically via the rep-agnostic
  `Mdeg` machinery across the sweep, then guess-and-verify a closed form in `(λ1,λ2,H1,H2)` — the
  general analog of V1's `M_k = N − n_k`. Headline deliverable, as in V1.
- **Validation sweep = mixed, small dims** (see the sweep list below): symmetric and non-symmetric
  reps, including reps with interior weight-multiplicity `> 1` to genuinely exercise the generalized
  `wtMult` and the harder degeneracy in `TauEigensystem3`.
- **Twist frame stays diagonal (`tg`)**, numeric `θ,h,z[k]` with symbolic-in-`u` eigenvalues —
  unchanged from V1 (`θ1=1/3, θ2=1/7, h=1`, three distinct high-precision moduli-1 twist
  eigenvalues; `χ1,χ2,χ3` their elementary symmetric functions).

## Existing state (`su3_V2_general.wb` = copy of `su3_V1.wb`)

Full V1 notebook: Section A (gl(3) GT generators `ps/dim/J0/Jp/Jm/JJ`; embeddings `Ee/id/idd/EE`;
`C1`; Lax `L[S]`; twist `G`/`z[k]`/`χ`; monodromies `T/Tg/TG`; Yangian weights `ν`/`ν1`; HWS;
quantum minors `qmin/qming/qminG`; transfer matrices `tg/tG`; `qdetT`; `Bax`) and Section B
(`tgCoeffs`, `H1op/H2op`, `wtMult`, `TauEigensystem3`, `Psi/tau1/tau2`, `BaxOp`, `alphaCoeffs`,
`Mdeg`, `QSolve3`, `qpoly`, `Qfun`, full-sweep verification). Everything parameterized by two
scalars `S1,S2`, each denoting the symmetric rep `{S,0,0}`.

**Key structural fact:** the GT engine `ps[rep]`, `dim`, `J0/Jp/Jm/JJ` already accepts an arbitrary
Young diagram `rep` — only the *embedding* layer hardcodes `{S,0,0}`. So most of the work is
re-signaturing, plus a small number of genuine math edits below.

## Change classification

### (a) Pure relabeling `S → λ`-vector (structurally unchanged)

Re-signature and swap `{S,0,0} → λ`:

- `Ee[λ_][i_,j_] := JJ[i,j,λ]`
- `id[λ_] := IdentityMatrix[dim[λ]]`
- `EE[λ1_,λ2_][i_,j_] := KroneckerProduct[Ee[λ1][i,j],id[λ2]] + KroneckerProduct[id[λ1],Ee[λ2][i,j]]`
- `C1[λ1_,λ2_]`
- `L[λ_][i_,j_][u_] := IdentityMatrix[dim[λ]] u KroneckerDelta[i,j] − h JJ[j,i,λ]`
- `T/Tg/TG[λ1_,λ2_]`, `qmin/qming/qminG[λ1_,λ2_]`, `tg/tG[λ1_,λ2_]`, `HWS[λ1_,λ2_]`
- Section B: `tgCoeffs`, `H1op/H2op`, `TauEigensystem3`, `Psi/tau1/tau2`, `BaxOp`, `alphaCoeffs`,
  `Mdeg`, `QSolve3`, `Qfun` — signature `S1,S2 → λ1,λ2`, internals otherwise unchanged. (`dim[{S,0,0}]`
  → `dim[λ]` inside `HWS` and tensor-dimension expressions.)

### (b) Genuine math changes

**1. Yangian weights `ν` — the one real physics edit.**
Derived from the Lax action on the highest-weight state: `L[λ][i,i][u] |hw⟩ = (u − h λ_i) |hw⟩`
(the weight of the HWS is `(λ1,λ2,λ3)`), so on the two-site HWS
`T[i,i](u) |hw⟩ = (u−θ1−h·λ1[[i]])(u−θ2−h·λ2[[i]]) |hw⟩`. Therefore:

```
ν[λ1_,λ2_][i_][u_] = (u − θ1 − h λ1[[i]]) (u − θ2 − h λ2[[i]])   for i = 1,2,3
```

V1's `ν[1]=(u−θ1−hS1)(u−θ2−hS2)` and `ν[2]=ν[3]=(u−θ1)(u−θ2)` are exactly the `λ=(S,0,0)` special
case. `ν1[λ1,λ2][u] := ν[λ1,λ2][1][u]` inherits the general form; `qdetT`, `Bax`, `BaxOp` keep their
V1 structure and pick up the general `ν` automatically:

```
qdetT[λ1,λ2][u] = χ3 · ν[λ1,λ2][1][u−2h] ν[λ1,λ2][2][u−h] ν[λ1,λ2][3][u]
```

**2. `wtMult` — rewritten as an honest tensor-product weight tally.**
V1's splitting combinatorics assumes each single-site weight has multiplicity 1 (true only for
symmetric reps). Replace with:

- `siteWeightMult[λ]` — tally the gl(3) weight `(w1,w2,w3)` of every GT pattern in `ps[{λ}]`
  (`w_i` = the `J0`-style row-sum difference), returning an Association `weight → multiplicity`.
- `wtMult[λ1_,λ2_][H1_,H2_]` — recover `(n1,n2,n3)` from `(H1,H2)` via `H1=n1−n2`, `H2=n2−n3`,
  `n1+n2+n3 = N` (`N =` total boxes `= Σ_i λ1_i + Σ_i λ2_i`), reject non-integer / negative, then
  return `Σ_{a+b=(n1,n2,n3)} siteWeightMult[λ1][a] · siteWeightMult[λ2][b]` (convolution of the two
  single-site tallies). This correctly yields multiplicity `> 1` and is the exact per-sector count
  `TauEigensystem3` must reproduce.

**3. Degree formula — re-derive the closed form.**
`Mdeg` (the large-`u` asymptotic master formula, cell ~29) is already rep-agnostic — it reads the
`α` coefficients off the actual eigenvalue asymptotics — so it stays and predicts the integer
degrees per state directly. New work: tabulate `M_k` vs `(λ1,λ2,H1,H2)` across the sweep and extract
the closed form. Candidate forms to test (must reduce to `N − n_k` when both reps are symmetric):

- `M_k = N − n_k` (V1 form, `N =` total boxes),
- `M_k = Λ_1 − n_k` with `Λ = λ1 + λ2` the pseudovacuum (component-wise sum of highest weights),
- forms built from partial sums of the `λ` components.

Confirm `deg q_k == M_k` (hard-assert) for the winning form across the whole sweep, keeping V1's
"`M_k` pairs with `z_k` in-index" ordering.

## Verification sweep (mixed, small dims)

Site pairs `(λ1, λ2)`:

| `λ1` | `λ2` | tensor dim | note |
|---|---|---|---|
| `{1,0,0}` | `{1,0,0}` | 9 | baseline, all mult 1 |
| `{2,0,0}` | `{1,1,0}` | 18 | symmetric ⊗ conjugate-fundamental |
| `{1,1,0}` | `{1,1,0}` | 9 | two conjugate-fundamentals |
| `{2,1,0}` | `{1,0,0}` | 24 | adjoint (8-dim, mult-2 zero weight) |
| `{2,1,0}` | `{1,1,0}` | 24 | adjoint ⊗ 3̄ |
| `{2,2,0}` | `{1,0,0}` | 18 | 6̄-type ⊗ fundamental |

Per state (hard-assert convention — a real failure `Abort`s):

1. **TQ residual** `≈ 0` (numeric, symbolic-in-`u`) for each `Q_k = z_k^{u/h} q_k(u)`.
2. **Casoratian / QQ** — the 3×3 discrete Wronskian of `{Q_1,Q_2,Q_3}` reproduces `qdetT` (up to the
   standard normalization); the 4×4 Casoratian reconstruction recovers `A1,A2,A3` (i.e. `τ1,τ2,qdetT`).
3. **Eigenvalue consistency** — `τ1,τ2` reconstructed from the `{Q_k}` match the `TauEigensystem3`
   eigenvalues.
4. **Sector counts** — states per `(H1,H2)` sector `== wtMult[λ1,λ2][H1,H2]`.
5. **Degree formula** — `deg q_k == M_k` from the extracted closed form.

Plus the Section A guard checks: `C1 == total boxes`; `qdetT == tg[·][3,1]` (catches any `ν`
shift-ordering error for non-symmetric reps); HWS eigenvalues of `T[i,i]` match `ν[i]`; `tgCoeffs`
leading coefficients are the expected `χ_a · 𝟙`.

## Risks

1. **`ν` shift-ordering in `qdetT`.** The `ν[1][u−2h] ν[2][u−h] ν[3][u]` shift assignment was only
   ever checked for symmetric reps (where `ν[2]=ν[3]`). For non-symmetric reps the three factors are
   distinct — verify against the directly-computed `tg[·][3,1]`. If it mismatches, the fix is a
   reassignment of shifts to `ν` factors, not a change to `ν` itself.
2. **`TauEigensystem3` on multiplicity-`> 1` sectors.** The generic-combination diagonalizer must
   fully split sectors that now genuinely have dimension `> 1`. Mitigated by the existing r-vector
   retry loop and the residual-genuineness / sector-count hard-asserts (it fails loudly, never
   silently passes). If a sector refuses to split because two states share *all* conserved charges,
   that is itself a finding to report, not a bug to paper over.
3. **Degree closed form may not be a single clean expression.** Fallback: the numeric-`Mdeg` degree
   guarantee (`deg q_k == Mdeg[…][k]`) still holds per state, so verification does not depend on
   finding the closed form; only the "figure it out" deliverable does.

## Out of scope

- Companion-twist SoV covectors / FSoV scalar product / CG-overlap dictionary (the later program).
- Clean promotion and companion `.md`/`.tex` reference docs.
- Lengths other than `L = 2`.
- `paper/xxxCG.tex` write-up.

## Open questions

None. Three choices confirmed with the user: λ-vectors everywhere (no scalar wrappers); re-derive
the degree closed form; mixed-small-dims sweep.
