# AI Evaluation Log — GaudinSolver

## 2026-07-22 07:28:52
**Input:**
```wolfram
ClearAll[mm, res1, res0, doRes0, doRes1, det1];
mm[u_, n_, S1_, S2_, w_] := Exp[w u] u^n / (u^(2 S1) (u - 1)^(2 S2));
res1[n_, S1_, S2_] := Exp[w] LaguerreL[2 S2 - 1, n - 2 S1 - 2 S2 + 1, -w];
res0[n_, S1_, S2_] := If[2 S1 > n, (-1)^(2 S1 + 2 S2 - n + 1) LaguerreL[2 S1 - n - 1, 1 - 2 S1 - 2 S2 + n, w], 0];
doRes0[S1_, S2_, f_] := Module[{nmax = Exponent[f, u], coeffs = CoefficientList[f, u]}, Plus @@ (coeffs Table[res0[n, S1, S2], {n, 0, nmax}])];
doRes1[S1_, S2_, f_] := Module[{nmax = Exponent[f, u], coeffs = CoefficientList[f, u]}, Plus @@ (coeffs Table[res1[n, S1, S2], {n, 0, nmax}])];
det1[S1_, S2_, mu_] := Module[{q = Q12[S1, S2, mu] /. z -> u}, Det[{{doRes0[S1, S2, q], doRes0[S1, S2, u q]}, {doRes1[S1, S2, q], doRes1[S1, S2, u q]}}]/(q /. u -> 0)/(q /. u -> 1)/Exp[w]];
"residue defs + det1 loaded"
```
**Output:**
```
Out= residue defs + det1 loaded
```

## 2026-07-22 07:29:01
**Input:**
```wolfram
(* cross-check doRes vs direct Residue for a Q12 example *)
q11 = Q12[1, 1, {2, 1, 1, 0}] /. z -> u;
brDirect0 = Residue[Exp[w u] q11/(u^(2*1) (u - 1)^(2*1)), {u, 0}];
brDirect1 = Residue[Exp[w u] (u q11)/(u^(2*1) (u - 1)^(2*1)), {u, 1}];
check0 = Simplify[doRes0[1, 1, q11] - brDirect0];
check1 = Simplify[doRes1[1, 1, u q11] - brDirect1];
(* det1 on examples, and test w-independence *)
d1 = det1[1, 1, {2, 1, 1, 0}];
d2 = det1[1, 2, {3, 2, 1, 0}];
d3 = det1[2, 1, {2, 2, 1, 1}];
{"doResVsResidue" -> {check0, check1}, "det1[1,1,{2,1,1,0}]" -> Simplify[d1], "wIndep1" -> (Simplify[D[d1, w]] === 0), "det1[1,2,{3,2,1,0}]" -> Simplify[d2], "wIndep2" -> (Simplify[D[d2, w]] === 0), "det1[2,1,{2,2,1,1}]" -> Simplify[d3], "wIndep3" -> (Simplify[D[d3, w]] === 0)}
```
**Output:**
```
Out= {"doResVsResidue" -> {0, 0}, "det1[1,1,{2,1,1,0}]" -> w^2, "wIndep1" -> False, "det1[1,2,{3,2,1,0}]" -> (w^3*(3 + 2*w))/12, "wIndep2" -> False, "det1[2,1,{2,2,1,1}]" -> (w^2*(-40 - 4*w + w^2))/6, "wIndep3" -> False}
```
