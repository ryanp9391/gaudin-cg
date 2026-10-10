## Deleted Cell 21 — 2026-09-16 08:36:18
**Kind:** markdown
```markdown
## 11. The determinant does not expand into the 8 wave functions

An expansion `det = sum_ns c_ns Psi_B[ns]` (equivalently, after dividing by `Psi({1},{0,0})`, an expansion of the normalised determinant in the normalised `<x|Psi>`) exists **iff** every direct-side coefficient `U_m` of

$$\det N \;=\; \sum_m (\text{dual monomial } m)\; U_m(\text{directs})$$

lies in the span of the eight `Psi_B[ns]`. Dividing by `Psi({1},{0,0})` does not change this: the factor cancels from both sides.

| | rank of `det`'s direct side | rank of the 8 `Psi_B` | rank together |
|---|---|---|---|
| abstract `Q` (free sample values) | 16 | 8 | 16 |
| `Q_i = z_i^{u/h}P_i`, `deg P_i <= 2` | 16 | 8 | 16 |

`rank together = 16 > 8` in both rows, so the eight `Psi_B` span a strict 8-dimensional subspace of the 16 directions the determinant actually uses. The expansion does not exist — and imposing the L=1 degree structure `deg P_i = lambda_1 - mu_i <= 2` does **not** help.

**Diagnosis.** The premise that everything must work with abstract `Q` uses the fact that `Q^j_A` and `Q^B_i` belong to different states and are therefore unrelated. That is correct, but it is not the constraint that is missing. The missing constraint is *within* one state: `Q_1^B` and `Q_2^B` are two solutions of the **same** Baxter equation, with the same `T_{a,1}`. In the parametrisation used above their coefficients `a0,a1,a2` and `b0,b1,b2` were varied independently, which is strictly more freedom than any actual state has.

On shell the direct data of a state is rigid: for the vacuum, `Q_1 = z_1^{u/h}` and `Q_2 = z_2^{u/h}P_2` with `P_2` fixed by cell 14, so only the two overall normalisations are free, and `<x|Psi> = Psi(x)/Psi({1},{0,0})` is completely determined (cell 17). The 16 directions are an artefact of letting `Q_1` and `Q_2` roam independently.

**Consequence for the test.** With a single state the comparison is vacuous (one equation, eight unknown weights). The decisive check needs the whole spectrum: build all 8 states, their direct `Q_1, Q_2` and dual `Q^1, Q^2`, and test

$$\det N^{AB} \;\stackrel{?}{=}\; \sum_{ns} c(ns)\,\langle \Psi_A | x_{ns}\rangle\,\langle x_{ns} | \Psi_B\rangle$$

over all 64 pairs `(A,B)` with only 8 unknown weights `c(ns)` — heavily overdetermined, and therefore conclusive either way.
```

## Deleted Cell 20 — 2026-09-16 08:36:18
**Kind:** code
```wolfram
(* ---- 14. Expand the 3x3 determinant into wave functions, normalised by Psi({1},{0,0}). ---- *)
ClearAll[PsiBa, PsiAa, gtT, dirCoeffs, subDir, rndPts, rankOf, ctrl, onsh];
(* wave functions on the SoV lattice, in the convention verified in cell 17 *)
PsiBa[{a_, b_, c_}] := au[1, a] (au[1, b] au[2, c - 1] - au[2, b] au[1, c - 1]);
PsiAa[{a_, b_, c_}] := ad[1, a] (ad[1, b] ad[2, c + 1] - ad[2, b] ad[1, c + 1]);
gtT = {{1, 0, 0}, {1, 1, 0}, {1, 1, 1}, {2, 0, 0}, {2, 1, 0}, {2, 1, 1}, {2, 2, 0}, {2, 2, 1}};
(* det = sum_m (dual monomial m) * U_m(directs); an expansion det = sum_ns c_ns(duals) Psi_B[ns]
   exists iff every U_m lies in the SPAN of the Psi_B[ns] with constant coefficients *)
dirCoeffs[e_] := Module[{vd, vu, cr},
   vd = Union@Cases[e, ad[_, _], Infinity]; vu = Union@Cases[e, au[_, _], Infinity];
   cr = CoefficientRules[Expand[e], Join[vd, vu]];
   Values@Merge[Table[(Times @@ (vd^r[[1, 1 ;; Length[vd]]])) -> r[[2]] (Times @@ (vu^r[[1, Length[vd] + 1 ;;]])), {r, cr}], Total]];
(* at L=1 every Q_i = z_i^{u/h} P_i with deg P_i = lambda_1 - mu_i <= 2, so the sampled values are constrained *)
subDir = {au[1, n_] :> z1^n (a0 + a1 n + a2 n^2), au[2, n_] :> z2^n (b0 + b1 n + b2 n^2)};
rndPts = Table[{a0 -> RandomInteger[{-9, 9}], a1 -> RandomInteger[{-9, 9}], a2 -> RandomInteger[{-9, 9}], b0 -> RandomInteger[{-9, 9}], b1 -> RandomInteger[{-9, 9}], b2 -> RandomInteger[{-9, 9}], z1 -> RandomInteger[{2, 9}], z2 -> RandomInteger[{2, 9}]}, {40}];
rankOf[exprs_] := MatrixRank[Table[Together[e /. p], {e, exprs}, {p, rndPts}]];
(* control: fully abstract Q, using the polynomial span *)
ctrl = Module[{Um = dirCoeffs[detA[trG[[1]]]], Pb = PsiBa /@ gtT, cm},
   cm[polys_] := Module[{vs = Union@Cases[polys, au[_, _], Infinity], mons},
     mons = Union@Flatten[First /@ CoefficientRules[#, vs] & /@ polys, 1];
     Table[Lookup[Association[CoefficientRules[p, vs]], Key /@ mons, 0], {p, polys}]];
   {MatrixRank[cm[Um]], MatrixRank[cm[Pb]], MatrixRank[cm[Join[Pb, Um]]]}];
(* on shell: impose the L=1 degree structure and redo the same rank comparison numerically *)
onsh = Module[{Um = (dirCoeffs[detA[trG[[1]]]] /. subDir), Pb = (PsiBa /@ gtT /. subDir)},
   {rankOf[Um], rankOf[Pb], rankOf[Join[Pb, Um]]}];
Grid[{
  {"", "rank of det's direct side", "rank of the 8 Psi_B", "rank of both together"},
  {"ABSTRACT Q (free sample values)", ctrl[[1]], ctrl[[2]], ctrl[[3]]},
  {"ON SHELL (Q_i = z_i^{u/h} P_i, deg P_i <= 2)", onsh[[1]], onsh[[2]], onsh[[3]]},
  {"expansion in the 8 Psi_B exists?", "abstract: " <> ToString[ctrl[[3]] == ctrl[[2]]], "on shell: " <> ToString[onsh[[3]] == onsh[[2]]], SpanFromLeft}
  }, Frame -> All, Alignment -> Left]
```

