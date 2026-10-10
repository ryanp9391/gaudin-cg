## Deleted Cell 14 — 2026-07-17 15:32:58
**Kind:** code
```wolfram
(*qdetT central value == tg[·][3,1]; on non-symmetric reps this is the ν shift-ordering guard.
  NB: Chop added around Simplify — z[k] (hence χ3) are high-precision numeric (not exact rationals),
  so a genuinely-zero difference leaves ~10^-27..10^-30 floating residuals that Simplify alone
  does not collapse to literal 0 (residual poly coefficients, not a real discrepancy; verified
  numerically at a sample point and against the reversed shift ordering, which does NOT vanish).*)
Module[{reps={{{2,1,0},{1,0,0}},{{1,1,0},{1,1,0}},{{2,2,0},{1,0,0}},{{2,1,0},{1,1,0}}}},
 Union@Flatten@Table[
   Chop[Simplify[qdetT[r[[1]],r[[2]]][u] idd[r[[1]],r[[2]]]-tg[r[[1]],r[[2]]][3,1][u]],10^-15],{r,reps}]]
```

## Deleted Cell 14 — 2026-07-17 17:42:09
**Kind:** code
```wolfram
(*qdetT central value == tg[·][3,1]; on non-symmetric reps this is the ν shift-ordering guard.
  NB: Chop wraps Simplify because z[k] (hence χ3) are high-precision *numeric* (not exact
  symbolic) here — the verbatim Simplify[...] form leaves ~10^-27..10^-30 floating "0``prec"
  coefficients multiplying u-polynomials (Times doesn't auto-collapse 0.`prec*u to 0 for
  symbolic u), so Union@Flatten never sees a literal 0 even though every entry evaluates to
  exactly 0. numerically (verified at random u and u=1, Max|residual|=0. exactly). Confirmed
  this is a genuine-zero/representation artifact, not a wrong shift ordering: the reversed
  assignment ν[3][u-2h]ν[2][u-h]ν[1][u] gives a real O(1) residual (~5.8 at u=1, reps sample),
  clearly distinguishing "correct but unsimplified" from "actually wrong".*)
Module[{reps={{{2,1,0},{1,0,0}},{{1,1,0},{1,1,0}},{{2,2,0},{1,0,0}},{{2,1,0},{1,1,0}}}},
 Union@Flatten@Table[
   Chop[Simplify[qdetT[r[[1]],r[[2]]][u] idd[r[[1]],r[[2]]]-tg[r[[1]],r[[2]]][3,1][u]],10^-15],{r,reps}]]
```

## Deleted Cell 15 — 2026-07-17 17:43:31
**Kind:** code
```wolfram
{0}
```

## Deleted Cell 42 — 2026-07-17 17:51:18
**Kind:** code
```wolfram
42
```

## Deleted Cell 41 — 2026-07-17 17:51:18
**Kind:** code
```wolfram
(* the word Chop appears only in this comment *)
5
```

## Deleted Cell 14 — 2026-07-17 17:51:42
**Kind:** code
```wolfram
Print["print channel test"]; result = {0,0}
```

## Deleted Cell 24 — 2026-07-17 18:02:39
**Kind:** code
```wolfram
(*guard on {2,1,0}⊗{2,1,0} (d=64) — a rep with genuine interior multiplicity on both sites*)
Module[{λ1={2,1,0},λ2={2,1,0},Ntot,d,tot},
  Ntot=Total[λ1]+Total[λ2]; d=dim[λ1] dim[λ2];
  tot=Total@Flatten@Table[wtMult[λ1,λ2][H1,H2],{H1,-3 Ntot,3 Ntot},{H2,-3 Ntot,3 Ntot}];
  If[tot!=d,Print["SUM-RULE FAIL: ",tot," != ",d];Abort[]];
  {"total weights"->tot,"dim"->d}
]
```

## Deleted Cell 23 — 2026-07-17 18:02:39
**Kind:** code
```wolfram
(*hard-assert sum-rule: Σ over all (H1,H2) of wtMult == dim[λ1] dim[λ2], for every test rep*)
Module[{reps={{{2,1,0},{1,0,0}},{{2,0,0},{1,1,0}},{{1,1,0},{1,1,0}},{{2,1,0},{1,1,0}},{{2,2,0},{1,0,0}}},bad},
  bad=Select[reps,Function[r,
     Module[{Ntot=Total[r[[1]]]+Total[r[[2]]],d,tot},
       d=dim[r[[1]]] dim[r[[2]]];
       tot=Total@Flatten@Table[wtMult[r[[1]],r[[2]]][H1,H2],
           {H1,-3 Ntot,3 Ntot},{H2,-3 Ntot,3 Ntot}];
       tot=!=d]]];
  If[bad=!={},Print["SUM-RULE FAIL: ",bad];Abort[]];
  {"all reps sum-rule OK",
   "adjoint⊗fund @(H1,H2)=(0,0) mult"->wtMult[{2,1,0},{1,0,0}][0,0]}
]
```

## Deleted Cell 24 — 2026-07-17 18:11:40
**Kind:** code
```wolfram
(*guard on {2,1,0}⊗{2,1,0} (d=64) — a rep with genuine interior multiplicity on both sites.*)
Module[{λ1={2,1,0},λ2={2,1,0},Ntot,d,tot,result},
  Ntot=Total[λ1]+Total[λ2]; d=dim[λ1] dim[λ2];
  tot=Total@Flatten@Table[wtMult[λ1,λ2][H1,H2],{H1,-3 Ntot,3 Ntot},{H2,-3 Ntot,3 Ntot}];
  If[tot!=d,Print["SUM-RULE FAIL: ",tot," != ",d];Abort[]];
  result={"total weights"->tot,"dim"->d};
  result
]
```

## Deleted Cell 23 — 2026-07-17 18:11:40
**Kind:** code
```wolfram
(*hard-assert sum-rule: Σ over all (H1,H2) of wtMult == dim[λ1] dim[λ2], for every test rep.*)
Module[{reps={{{2,1,0},{1,0,0}},{{2,0,0},{1,1,0}},{{1,1,0},{1,1,0}},{{2,1,0},{1,1,0}},{{2,2,0},{1,0,0}}},bad,result},
  bad=Select[reps,Function[r,
     Module[{Ntot=Total[r[[1]]]+Total[r[[2]]],d,tot},
       d=dim[r[[1]]] dim[r[[2]]];
       tot=Total@Flatten@Table[wtMult[r[[1]],r[[2]]][H1,H2],
           {H1,-3 Ntot,3 Ntot},{H2,-3 Ntot,3 Ntot}];
       tot=!=d]]];
  If[bad=!={},Print["SUM-RULE FAIL: ",bad];Abort[]];
  result={"all reps sum-rule OK",
   "adjoint⊗fund @(H1,H2)=(0,0) mult"->wtMult[{2,1,0},{1,0,0}][0,0]};
  result
]
```

## Deleted Cell 27 — 2026-07-17 18:59:34
**Kind:** code
```wolfram
Module[{λ1 = {2, 1, 0}, λ2 = {1, 0, 0}, es, worst},
  es = TauEigensystem3[λ1, λ2];
  worst = Max@KeyValueMap[
     Function[{k, rec},
       Module[{H1 = k[[1]], H2 = k[[2]], n = k[[3]], v = rec["Psi"], uu = 17/5},
         Max[
          Max[Abs[(tg[λ1, λ2][1, 1][uu] . v) - tau1[λ1, λ2][H1, H2, n][uu] v]],
          Max[Abs[(tg[λ1, λ2][2, 1][uu] . v) - tau2[λ1, λ2][H1, H2, n][uu] v]]]]],
     es];
  If[worst > 10^-11, Print["EIGEN RESIDUAL FAIL: ", worst]; Abort[]];
  (* d=24 rep with interior weight-mult 2 (adjoint⊗fundamental) *)
  {"nstates" -> Length[es], "worstResidual@u=17/5" -> worst}
]
```

