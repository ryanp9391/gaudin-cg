# AI Evaluation Log — CGGenFun

## 2026-07-22 07:17:07
**Input:**
```wolfram
ClearAll[schur, tensorDecomp, dimGL4, irreps, muhat, xw, Pw, Cas2, Cas4, degQ, gaudinODE, wsub, GaudinSolve, BuildQs, Q12, cq];
schur[lam_, x_List] := Cancel[Det[Table[x[[i]]^(lam[[j]] + Length[x] - j), {i, Length[x]}, {j, Length[x]}]]/Det[Table[x[[i]]^(Length[x] - j), {i, Length[x]}, {j, Length[x]}]]];
dimGL4[lam_] := Product[(lam[[i]] - lam[[j]] + j - i)/(j - i), {i, 4}, {j, i + 1, 4}];
tensorDecomp[lam_, mu_] := Module[{x = Array[xx, 4], prod, rules, nu, mult, out = {}}, prod = Expand[schur[lam, x] schur[mu, x]]; While[prod =!= 0, rules = Association[CoefficientRules[prod, x]]; nu = Last[Sort[Keys[rules]]]; mult = rules[nu]; AppendTo[out, nu -> mult]; prod = Expand[prod - mult schur[nu, x]]]; out];
irreps[S1_, S2_] := tensorDecomp[{S1, S1, 0, 0}, {S2, S2, 0, 0}][[All, 1]];
muhat[mu_] := Table[mu[[i]] + 4 - i, {i, 4}];
xw[mu_] := Table[mu[[i]] + 5/2 - i, {i, 4}];
Pw[mu_, k_] := Total[xw[mu]^k];
Cas2[mu_] := Pw[mu, 2] - 5;
Cas4[mu_] := Pw[mu, 4] + 4 Pw[mu, 3] - Pw[mu, 1] Pw[mu, 2] + Pw[mu, 2]/2 - 2 Pw[mu, 1]^2 - 10 Pw[mu, 1] - 51/4;
degQ[S1_, S2_, mu_] := Table[(S1 + S2 + 3) - muhat[mu][[k]], {k, 4}];
gaudinODE = (2*C2 + C2^2 - 2*C4 + 24*P + 2*C2*P + 20*P^2 + 4*P^3 + 24*S + 2*C2*S + 40*P*S + 12*P^2*S + 20*S^2 + 12*P*S^2 + 4*S^3)*Q[z] + 4*(C2 - 4*P - 2*P^2 - 4*S - 4*P*S - 2*S^2)*(2*z + P*z + S*z - w1 - P*w1 - w2 - S*w2)*Derivative[1][Q][z] - 4*(C2*z^2 - 6*P*z^2 - 4*P^2*z^2 - 6*S*z^2 - 8*P*S*z^2 - 4*S^2*z^2 - C2*z*w1 + 8*P*z*w1 + 6*P^2*z*w1 + 4*S*z*w1 + 8*P*S*z*w1 + 2*S^2*z*w1 - 2*P*w1^2 - 2*P^2*w1^2 - C2*z*w2 + 4*P*z*w2 + 2*P^2*z*w2 + 8*S*z*w2 + 8*P*S*z*w2 + 6*S^2*z*w2 + C2*w1*w2 - 4*P*w1*w2 - 2*P^2*w1*w2 - 4*S*w1*w2 - 8*P*S*w1*w2 - 2*S^2*w1*w2 - 2*S*w2^2 - 2*S^2*w2^2)*Derivative[2][Q][z] - 16*(z - w1)*(z - w2)*(P*z + S*z - P*w1 - S*w2)*Derivative[3][Q][z] + 8*(z - w1)^2*(z - w2)^2*Derivative[4][Q][z];
wsub = {w1 -> 0, w2 -> 1};
GaudinSolve[S1_, S2_, mu_] := Module[{ds, ode}, ds = degQ[S1, S2, mu]; ode = gaudinODE /. {C2 -> Cas2[mu], C4 -> Cas4[mu], S -> S1, P -> S2} /. wsub; Table[Module[{m = ds[[kk]], gauge = ds[[1 ;; kk - 1]], q, cl, sol}, q = Sum[cq[j] z^j, {j, 0, m}]; cl = CoefficientList[Expand[ode /. Q -> Function[x, q /. z -> x]], z]; sol = Solve[Join[Thread[cl == 0], {cq[m] == 1}, Thread[(cq /@ gauge) == 0]], Table[cq[j], {j, 0, m}]]; Table[cq[j], {j, 0, m}] /. sol[[1]]], {kk, 4}]];
BuildQs[S1_, S2_, mu_] := Module[{ds = degQ[S1, S2, mu], coeffs = GaudinSolve[S1, S2, mu]}, Table[Sum[coeffs[[kk, j + 1]] z^j, {j, 0, ds[[kk]]}], {kk, 4}]];
Q12[S1_, S2_, mu_] := Expand[Wronskian[BuildQs[S1, S2, mu][[{1, 2}]], z]];
"defs loaded"
```
**Output:**
```
Out= defs loaded
```

## 2026-07-22 07:17:16
**Input:**
```wolfram
odeCheck[S1_, S2_, mu_] := Module[{qs = BuildQs[S1, S2, mu], ode}, ode = gaudinODE /. {C2 -> Cas2[mu], C4 -> Cas4[mu], S -> S1, P -> S2} /. wsub; Table[Expand[ode /. Q -> Function[x, qs[[k]] /. z -> x]], {k, 4}]];
worst = Max[Abs[Flatten[Table[odeCheck[S1, S2, mu], {S1, 1, 2}, {S2, 1, 2}, {mu, irreps[S1, S2]}]] /. z -> 37]];
ds = degQ[1, 1, {2, 1, 1, 0}];
ex = Q12[1, 1, {2, 1, 1, 0}];
qs = BuildQs[1, 1, {2, 1, 1, 0}];
manual = Expand[qs[[1]] D[qs[[2]], z] - qs[[2]] D[qs[[1]], z]];
{"worstODEresidual" -> worst, "Q12ex" -> ex, "matchesManualWronskian" -> (Expand[ex - manual] === 0), "deg1deg2" -> ds[[{1, 2}]], "degQ12" -> Exponent[ex, z], "expectedDeg" -> ds[[1]] + ds[[2]] - 1, "Q12_1_2_example_deg" -> Exponent[Q12[1, 2, {3, 2, 1, 0}], z]}
```
**Output:**
```
Out= {"worstODEresidual" -> 0, "Q12ex" -> -1 + 2*z, "matchesManualWronskian" -> True, "deg1deg2" -> {0, 2}, "degQ12" -> 1, "expectedDeg" -> 1, "Q12_1_2_example_deg" -> 1}
```
## 2026-07-22 07:21:08 — ▶️ RUN CELLS 3–11
Cell 3: ✓ — (no output)
Cell 4: ✓ — (no output)
Cell 5: ✓ — (no output)
Cell 6: ✓ — (no output)
Cell 7: ✓ — (no output)
Cell 9: ✓ — (no output)
Cell 11: ✓ — (no output)


## 2026-07-22 07:21:11
**Input:**
```wolfram
{worstODE, Q12ok, Q12[1, 1, {2, 1, 1, 0}], irreps[1, 1]}
```
**Output:**
```
Out= {0, True, -1 + 2*z, {{2, 2, 0, 0}, {2, 1, 1, 0}, {1, 1, 1, 1}}}
```
## 2026-07-22 07:30:55 — 📥 BULK INSERT 7 CELL(S) at positions 10–16
1. [markdown] ## Residue primitives (copied from `FastQ.wb`)↵↵Measure `1/(u^(2 S1) (u-1)^(2 S2))`. For a polynomia
2. [code] ClearAll[mm, res1, res0];↵mm[u_, n_, S1_, S2_, w_] := Exp[w u] u^n / (u^(2 S1) (u - 1)^(2 S2));↵res1
3. [code] (* br[0][e^{w u} f] = doRes0[S1,S2,f];  br[1][e^{w u} f] = doRes1[S1,S2,f].  f a polynomial in u. *)
4. [markdown] ## `det1` -- 2x2 residue determinant↵↵$$\det{}_1[S_1,S_2,\mu]=\frac{1}{Q_{12}(0)\,Q_{12}(1)\,e^{w}}\
5. [code] ClearAll[det1];↵det1[S1_, S2_, mu_] := Module[{q = Q12[S1, S2, mu] /. z -> u},↵   Det[{{doRes0[S1, S
6. [markdown] ## Validation↵↵Solver (`worstODE`: every gauge-fixed `Q_k` in the Gaudin ODE over `(S1,S2)∈{1,2}²`, 
7. [code] worstODE = Max[Abs[Flatten[Table[↵      Module[{qs = BuildQs[S1, S2, mu], ode = gaudinODE /. {C2 ->


## 2026-07-22 07:30:56
**Input:**
```wolfram
ClearAll[mm, res1, res0];
mm[u_, n_, S1_, S2_, w_] := Exp[w u] u^n / (u^(2 S1) (u - 1)^(2 S2));
res1[n_, S1_, S2_] := Exp[w] LaguerreL[2 S2 - 1, n - 2 S1 - 2 S2 + 1, -w];
res0[n_, S1_, S2_] := If[2 S1 > n, (-1)^(2 S1 + 2 S2 - n + 1) LaguerreL[2 S1 - n - 1, 1 - 2 S1 - 2 S2 + n, w], 0];
```
**Output:**
```
(no output)
```

## 2026-07-22 07:30:56
**Input:**
```wolfram
(* br[0][e^{w u} f] = doRes0[S1,S2,f];  br[1][e^{w u} f] = doRes1[S1,S2,f].  f a polynomial in u. *)
doRes0[S1_, S2_, f_] := Module[{nmax, coeffs, res},
  nmax = Exponent[f, u];
  coeffs = CoefficientList[f, u];
  res = Plus @@ (coeffs Table[res0[n, S1, S2], {n, 0, nmax}]);
  res
  ]

doRes1[S1_, S2_, f_] := Module[{nmax, coeffs, res},
  nmax = Exponent[f, u];
  coeffs = CoefficientList[f, u];
  res = Plus @@ (coeffs Table[res1[n, S1, S2], {n, 0, nmax}]);
  res
  ]
```
**Output:**
```
(no output)
```

## 2026-07-22 07:30:56
**Input:**
```wolfram
ClearAll[det1];
det1[S1_, S2_, mu_] := Module[{q = Q12[S1, S2, mu] /. z -> u},
   Det[{{doRes0[S1, S2, q], doRes0[S1, S2, u q]},
        {doRes1[S1, S2, q], doRes1[S1, S2, u q]}}]/(q /. u -> 0)/(q /. u -> 1)/Exp[w]
   ];
```
**Output:**
```
(no output)
```

## 2026-07-22 07:30:57
**Input:**
```wolfram
worstODE = Max[Abs[Flatten[Table[
      Module[{qs = BuildQs[S1, S2, mu], ode = gaudinODE /. {C2 -> Cas2[mu], C4 -> Cas4[mu], S -> S1, P -> S2} /. wsub},
        Table[Expand[ode /. Q -> Function[x, qs[[k]] /. z -> x]], {k, 4}]],
      {S1, 1, 2}, {S2, 1, 2}, {mu, irreps[S1, S2]}]] /. z -> 37]];
qchk = Q12[1, 1, {2, 1, 1, 0}] /. z -> u;
brChk = Simplify[{doRes0[1, 1, qchk] - Residue[Exp[w u] qchk/(u^2 (u - 1)^2), {u, 0}],
    doRes1[1, 1, u qchk] - Residue[Exp[w u] (u qchk)/(u^2 (u - 1)^2), {u, 1}]}];
Print["worstODE = ", worstODE, "  (expect 0)"];
Print["br vs Residue = ", brChk, "  (expect {0,0})"];
Print["det1[1,1,{2,1,1,0}] = ", Simplify[det1[1, 1, {2, 1, 1, 0}]]];
Print["det1[1,2,{3,2,1,0}] = ", Simplify[det1[1, 2, {3, 2, 1, 0}]]];
Print["det1[2,1,{2,2,1,1}] = ", Simplify[det1[2, 1, {2, 2, 1, 1}]]];
```
**Output:**
```
(no output)
```

## 2026-07-22 07:31:11
**Input:**
```wolfram
{"worstODE" -> worstODE, "brChk" -> brChk, "det1_1_1" -> Simplify[det1[1, 1, {2, 1, 1, 0}]], "det1_1_2" -> Simplify[det1[1, 2, {3, 2, 1, 0}]], "det1_2_1" -> Simplify[det1[2, 1, {2, 2, 1, 1}]]}
```
**Output:**
```
Out= {"worstODE" -> 0, "brChk" -> {0, 0}, "det1_1_1" -> w^2, "det1_1_2" -> (w^3*(3 + 2*w))/12, "det1_2_1" -> (w^2*(-40 - 4*w + w^2))/6}
```
## 2026-07-22 07:33:51 — ✏️ EDIT CELL 11 [batch]
ClearAll[mm, res1, res0];
mm[z_, n_, S1_, S2_, w_] := Exp[w z] z^n / (z^(2 S1) (z - 1)^(2 S2));
res1[n_, S1_, S2_] := Exp[w] LaguerreL[2 S2 - 1, n - 2 S1 - 2 S2 + 1, -w];
res0[n_, S1_, S2_] := If[2 S1…


## 2026-07-22 07:33:51
**Input:**
```wolfram
ClearAll[mm, res1, res0];
mm[z_, n_, S1_, S2_, w_] := Exp[w z] z^n / (z^(2 S1) (z - 1)^(2 S2));
res1[n_, S1_, S2_] := Exp[w] LaguerreL[2 S2 - 1, n - 2 S1 - 2 S2 + 1, -w];
res0[n_, S1_, S2_] := If[2 S1 > n, (-1)^(2 S1 + 2 S2 - n + 1) LaguerreL[2 S1 - n - 1, 1 - 2 S1 - 2 S2 + n, w], 0];
```
**Output:**
```
Out= (no output)
```
## 2026-07-22 07:33:51 — ✏️ EDIT CELL 12 [batch]
(* br[0][e^{w z} f] = doRes0[S1,S2,f];  br[1][e^{w z} f] = doRes1[S1,S2,f].  f a polynomial in z. *)
doRes0[S1_, S2_, f_] := Module[{nmax, coeffs, res},
  nmax = Exponent[f, z];
  coeffs = Coefficient…


## 2026-07-22 07:33:51
**Input:**
```wolfram
(* br[0][e^{w z} f] = doRes0[S1,S2,f];  br[1][e^{w z} f] = doRes1[S1,S2,f].  f a polynomial in z. *)
doRes0[S1_, S2_, f_] := Module[{nmax, coeffs, res},
  nmax = Exponent[f, z];
  coeffs = CoefficientList[f, z];
  res = Plus @@ (coeffs Table[res0[n, S1, S2], {n, 0, nmax}]);
  res
  ]

doRes1[S1_, S2_, f_] := Module[{nmax, coeffs, res},
  nmax = Exponent[f, z];
  coeffs = CoefficientList[f, z];
  res = Plus @@ (coeffs Table[res1[n, S1, S2], {n, 0, nmax}]);
  res
  ]
```
**Output:**
```
Out= (no output)
```
## 2026-07-22 07:33:51 — ✏️ EDIT CELL 14 [batch]
ClearAll[det1];
det1[S1_, S2_, mu_] := Module[{q = Q12[S1, S2, mu]},
   Det[{{doRes0[S1, S2, q], doRes0[S1, S2, z q]},
        {doRes1[S1, S2, q], doRes1[S1, S2, z q]}}]/(q /. z -> 0)/(q /. z -> 1)/Ex…


## 2026-07-22 07:33:51
**Input:**
```wolfram
ClearAll[det1];
det1[S1_, S2_, mu_] := Module[{q = Q12[S1, S2, mu]},
   Det[{{doRes0[S1, S2, q], doRes0[S1, S2, z q]},
        {doRes1[S1, S2, q], doRes1[S1, S2, z q]}}]/(q /. z -> 0)/(q /. z -> 1)/Exp[w]
   ];
```
**Output:**
```
Out= (no output)
```
## 2026-07-22 07:33:51 — ✏️ EDIT CELL 16 [batch]
worstODE = Max[Abs[Flatten[Table[
      Module[{qs = BuildQs[S1, S2, mu], ode = gaudinODE /. {C2 -> Cas2[mu], C4 -> Cas4[mu], S -> S1, P -> S2} /. wsub},
        Table[Expand[ode /. Q -> Function[x, q…


## 2026-07-22 07:33:51
**Input:**
```wolfram
worstODE = Max[Abs[Flatten[Table[
      Module[{qs = BuildQs[S1, S2, mu], ode = gaudinODE /. {C2 -> Cas2[mu], C4 -> Cas4[mu], S -> S1, P -> S2} /. wsub},
        Table[Expand[ode /. Q -> Function[x, qs[[k]] /. z -> x]], {k, 4}]],
      {S1, 1, 2}, {S2, 1, 2}, {mu, irreps[S1, S2]}]] /. z -> 37]];
qchk = Q12[1, 1, {2, 1, 1, 0}];
brChk = Simplify[{doRes0[1, 1, qchk] - Residue[Exp[w z] qchk/(z^2 (z - 1)^2), {z, 0}],
    doRes1[1, 1, z qchk] - Residue[Exp[w z] (z qchk)/(z^2 (z - 1)^2), {z, 1}]}];
Print["worstODE = ", worstODE, "  (expect 0)"];
Print["br vs Residue = ", brChk, "  (expect {0,0})"];
Print["det1[1,1,{2,1,1,0}] = ", Simplify[det1[1, 1, {2, 1, 1, 0}]]];
Print["det1[1,2,{3,2,1,0}] = ", Simplify[det1[1, 2, {3, 2, 1, 0}]]];
Print["det1[2,1,{2,2,1,1}] = ", Simplify[det1[2, 1, {2, 2, 1, 1}]]];
```
**Output:**
```
Out= (no output)
```
## 2026-07-22 07:33:58 — ✏️ EDIT CELL 10
## Residue primitives (adapted from `FastQ.wb`, in the spectral variable `z`)

Measure `1/(z^(2 S1) (z-1)^(2 S2))`. For a polynomial `f(z)`, the residue of `e^{w z} f(z)` against this measure is
`br[a…

## 2026-07-22 07:34:04 — ✏️ EDIT CELL 13
## `det1` -- 2x2 residue determinant

$$\det{}_1[S_1,S_2,\mu]=\frac{1}{Q_{12}(0)\,Q_{12}(1)\,e^{w}}\;\det\begin{pmatrix} \mathrm{br}_0[e^{wz}Q_{12}(z)] & \mathrm{br}_0[e^{wz}\,z\,Q_{12}(z)]\\ \mathrm{…

## 2026-07-22 07:34:07 — ✏️ EDIT CELL 15
## Validation

Solver (`worstODE`: every gauge-fixed `Q_k` in the Gaudin ODE over `(S1,S2) ∈ {1,2}²`, all irreps), `br` entries vs direct `Residue` (in `z`), and a few `det1` values.


## 2026-07-22 07:34:12
**Input:**
```wolfram
{"worstODE" -> worstODE, "brChk" -> brChk, "det1_1_1" -> Simplify[det1[1, 1, {2, 1, 1, 0}]], "det1_1_2" -> Simplify[det1[1, 2, {3, 2, 1, 0}]], "det1_2_1" -> Simplify[det1[2, 1, {2, 2, 1, 1}]]}
```
**Output:**
```
Out= {"worstODE" -> 0, "brChk" -> {0, 0}, "det1_1_1" -> w^2, "det1_1_2" -> (w^3*(3 + 2*w))/12, "det1_2_1" -> (w^2*(-40 - 4*w + w^2))/6}
```

---

## 2026-07-22 07:34:47 — 🔄 KERNEL RESTART


---

## 2026-07-22 07:52:05 — 🔄 KERNEL RESTART


## 2026-07-22 07:53:02
**Input:**
```wolfram
ir = irreps[5, 5]; mu0 = ir[[Ceiling[Length[ir]/2]]];
{"numIrreps" -> Length[ir], "mu0" -> mu0, "degQ" -> degQ[5, 5, mu0],
 "t_GaudinSolve" -> First@AbsoluteTiming[GaudinSolve[5, 5, mu0];],
 "t_Q12" -> First@AbsoluteTiming[Q12[5, 5, mu0];],
 "t_det1" -> First@AbsoluteTiming[det1[5, 5, mu0];]}
```
**Output:**
```
Out= {"numIrreps" -> 2, "mu0" -> 5, "degQ" -> degQ[5, 5, 5], "t_GaudinSolve" -> 1.*^-6, "t_Q12" -> 1.*^-6, "t_det1" -> 1.*^-6}
```

## 2026-07-22 07:53:19
**Input:**
```wolfram
{"irreps11" -> irreps[1, 1], "irreps55_head" -> Head[irreps[5, 5]], "irreps55_len" -> Length[irreps[5, 5]], "irreps55_first3" -> Take[irreps[5, 5], 3], "Q12def?" -> (Head[Q12] === Symbol && DownValues[Q12] =!= {})}
```
**Output:**
```
[message] Take::take: Cannot take positions 1 through 3 in irreps[5, 5].
Out= {"irreps11" -> irreps[1, 1], "irreps55_head" -> irreps, "irreps55_len" -> 2, "irreps55_first3" -> Take[irreps[5, 5], 3], "Q12def?" -> False}
```
## 2026-07-22 07:53:35 — ▶️ RUN CELLS 1–16
Cell 3: ✓ — (no output)
Cell 4: ✓ — (no output)
Cell 5: ✓ — (no output)
Cell 6: ✓ — (no output)
Cell 7: ✓ — (no output)
Cell 9: ✓ — (no output)
Cell 11: ✓ — (no output)
Cell 12: ✓ — (no output)
Cell 14: ✓ — (no output)
Cell 15: ✓ — Out[30]= \left\{0.025927,2\,z-1\right\}
Cell 16: ✓ — (no output)


## 2026-07-22 07:53:44
**Input:**
```wolfram
ir = irreps[5, 5]; mu0 = ir[[Ceiling[Length[ir]/2]]];
q0 = Q12[5, 5, mu0];
{"numIrreps" -> Length[ir], "mu0" -> mu0, "degQ" -> degQ[5, 5, mu0], "degQ12" -> Exponent[q0, z],
 "t_GaudinSolve" -> First@AbsoluteTiming[GaudinSolve[5, 5, mu0];],
 "t_BuildQs" -> First@AbsoluteTiming[BuildQs[5, 5, mu0];],
 "t_Q12" -> First@AbsoluteTiming[Q12[5, 5, mu0];],
 "t_doRes_both" -> First@AbsoluteTiming[{doRes0[5, 5, q0], doRes0[5, 5, z q0], doRes1[5, 5, q0], doRes1[5, 5, z q0]};],
 "t_det1" -> First@AbsoluteTiming[det1[5, 5, mu0];]}
```
**Output:**
```
Out= {"numIrreps" -> 21, "mu0" -> {9, 5, 5, 1}, "degQ" -> {1, 6, 7, 12}, "degQ12" -> 6, "t_GaudinSolve" -> 0.002803, "t_BuildQs" -> 0.002599, "t_Q12" -> 0.002635, "t_doRes_both" -> 0.00299, "t_det1" -> 0.008335}
```

## 2026-07-22 07:54:39
**Input:**
```wolfram
(* baseline full sweeps *)
t55 = First@AbsoluteTiming[det1[5, 5, #] & /@ irreps[5, 5];];
t88 = First@AbsoluteTiming[det1[8, 8, #] & /@ irreps[8, 8];];
(* how much of GaudinSolve is the high-degree Q3,Q4? time a 2-Q solver *)
GaudinSolve2[S1_, S2_, mu_] := Module[{ds = degQ[S1, S2, mu], ode},
   ode = gaudinODE /. {C2 -> Cas2[mu], C4 -> Cas4[mu], S -> S1, P -> S2} /. wsub;
   Table[Module[{m = ds[[kk]], gauge = ds[[1 ;; kk - 1]], q, cl, sol},
       q = Sum[cq[j] z^j, {j, 0, m}];
       cl = CoefficientList[Expand[ode /. Q -> Function[x, q /. z -> x]], z];
       sol = Solve[Join[Thread[cl == 0], {cq[m] == 1}, Thread[(cq /@ gauge) == 0]], Table[cq[j], {j, 0, m}]];
       Table[cq[j], {j, 0, m}] /. sol[[1]]], {kk, 2}]];
mu8 = irreps[8, 8][[Ceiling[Length[irreps[8, 8]]/2]]];
{"sweep_det1_5,5(21)" -> t55, "sweep_det1_8,8(N)" -> t88, "Nirreps88" -> Length[irreps[8, 8]],
 "mu8" -> mu8, "degQ_8,8" -> degQ[8, 8, mu8],
 "t_GaudinSolve4" -> First@AbsoluteTiming[GaudinSolve[8, 8, mu8];],
 "t_GaudinSolve2only" -> First@AbsoluteTiming[GaudinSolve2[8, 8, mu8];]}
```
**Output:**
```
Out= {"sweep_det1_5,5(21)" -> 0.629618, "sweep_det1_8,8(N)" -> 4.396827, "Nirreps88" -> 45, "mu8" -> {14, 9, 7, 2}, "degQ_8,8" -> {2, 8, 11, 17}, "t_GaudinSolve4" -> 0.004248, "t_GaudinSolve2only" -> 0.001519}
```

## 2026-07-22 07:55:41
**Input:**
```wolfram
(* (A) memoized residues; (B) 2-Q Wronskian; combined fast det1 with residue tables built once *)
ClearAll[res0m, res1m, Q12fast, det1fast];
res0m[n_, S1_, S2_] := res0m[n, S1, S2] = If[2 S1 > n, (-1)^(2 S1 + 2 S2 - n + 1) LaguerreL[2 S1 - n - 1, 1 - 2 S1 - 2 S2 + n, w], 0];
res1m[n_, S1_, S2_] := res1m[n, S1, S2] = Exp[w] LaguerreL[2 S2 - 1, n - 2 S1 - 2 S2 + 1, -w];
Q12fast[S1_, S2_, mu_] := Module[{ds = degQ[S1, S2, mu], ode, qs},
   ode = gaudinODE /. {C2 -> Cas2[mu], C4 -> Cas4[mu], S -> S1, P -> S2} /. wsub;
   qs = Table[Module[{m = ds[[kk]], gauge = ds[[1 ;; kk - 1]], q, cl, sol},
       q = Sum[cq[j] z^j, {j, 0, m}];
       cl = CoefficientList[Expand[ode /. Q -> Function[x, q /. z -> x]], z];
       sol = Solve[Join[Thread[cl == 0], {cq[m] == 1}, Thread[(cq /@ gauge) == 0]], Table[cq[j], {j, 0, m}]];
       Sum[(cq[j] /. sol[[1]]) z^j, {j, 0, m}]], {kk, 2}];
   Expand[qs[[1]] D[qs[[2]], z] - qs[[2]] D[qs[[1]], z]]];
det1fast[S1_, S2_, mu_] := Module[{q = Q12fast[S1, S2, mu], c, L, R0, R1},
   c = CoefficientList[q, z]; L = Length[c];
   R0 = Table[res0m[n, S1, S2], {n, 0, L}];
   R1 = Table[res1m[n, S1, S2], {n, 0, L}];
   Det[{{c.R0[[1 ;; L]], c.R0[[2 ;; L + 1]]}, {c.R1[[1 ;; L]], c.R1[[2 ;; L + 1]]}}]/(q /. z -> 0)/(q /. z -> 1)/Exp[w]];
(* correctness: identical to det1 across full sweeps *)
chk55 = Simplify[(det1[5, 5, #] - det1fast[5, 5, #]) & /@ irreps[5, 5]] // Flatten // Union;
chk88 = Simplify[(det1[8, 8, #] - det1fast[8, 8, #]) & /@ irreps[8, 8]] // Flatten // Union;
{"match55" -> chk55, "match88" -> chk88}
```
**Output:**
```
Out= {"match55" -> {0}, "match88" -> {0}}
```

## 2026-07-22 07:56:08
**Input:**
```wolfram
clearMemo := (ClearAll[res0m, res1m];
   res0m[n_, S1_, S2_] := res0m[n, S1, S2] = If[2 S1 > n, (-1)^(2 S1 + 2 S2 - n + 1) LaguerreL[2 S1 - n - 1, 1 - 2 S1 - 2 S2 + n, w], 0];
   res1m[n_, S1_, S2_] := res1m[n, S1, S2] = Exp[w] LaguerreL[2 S2 - 1, n - 2 S1 - 2 S2 + 1, -w];);
bench[S1_, S2_] := Module[{ir = irreps[S1, S2], tb, tf},
   tb = First@AbsoluteTiming[det1[S1, S2, #] & /@ ir;];
   clearMemo;
   tf = First@AbsoluteTiming[det1fast[S1, S2, #] & /@ ir;];
   {"N" -> Length[ir], "baseline" -> tb, "fast" -> tf, "speedup" -> tb/tf}];
{"S=5" -> bench[5, 5], "S=8" -> bench[8, 8], "S=10" -> bench[10, 10]}
```
**Output:**
```
Out= {"S=5" -> {"N" -> 21, "baseline" -> 0.168544, "fast" -> 0.081866, "speedup" -> 2.0587789802848557}, "S=8" -> {"N" -> 45, "baseline" -> 1.499916, "fast" -> 0.103521, "speedup" -> 14.489002231431304}, "S=10" -> {"N" -> 66, "baseline" -> 6.725356, "fast" -> 0.143924, "speedup" -> 46.72852338734332}}
```

## 2026-07-22 07:57:02
**Input:**
```wolfram
(* memoize the REAL res0/res1 (used by doRes) + make a Q12 that solves only Q1,Q2; keep det1 calling doRes *)
memoizeRes := (ClearAll[res0, res1];
   res0[n_, S1_, S2_] := res0[n, S1, S2] = If[2 S1 > n, (-1)^(2 S1 + 2 S2 - n + 1) LaguerreL[2 S1 - n - 1, 1 - 2 S1 - 2 S2 + n, w], 0];
   res1[n_, S1_, S2_] := res1[n, S1, S2] = Exp[w] LaguerreL[2 S2 - 1, n - 2 S1 - 2 S2 + 1, -w];);
(* det1orig still defined as current notebook det1 (calls doRes0/doRes1, uses Q12) *)
det1orig[S1_, S2_, mu_] := Module[{q = Q12[S1, S2, mu]}, Det[{{doRes0[S1, S2, q], doRes0[S1, S2, z q]}, {doRes1[S1, S2, q], doRes1[S1, S2, z q]}}]/(q /. z -> 0)/(q /. z -> 1)/Exp[w]];
benchMemo[S1_, S2_] := Module[{ir = irreps[S1, S2], tMemoOrig},
   memoizeRes;  (* cold cache *)
   tMemoOrig = First@AbsoluteTiming[det1orig[S1, S2, #] & /@ ir;];
   tMemoOrig];
{"S=8_memoOnly_origDet1" -> benchMemo[8, 8], "S=10_memoOnly_origDet1" -> benchMemo[10, 10]}
```
**Output:**
```
Out= {"S=8_memoOnly_origDet1" -> 0.313715, "S=10_memoOnly_origDet1" -> 0.731402}
```

## 2026-07-22 07:57:45
**Input:**
```wolfram
(* Q12 solving only Q1,Q2 (same result), memoized res, det1 still readable via doRes *)
Q12two[S1_, S2_, mu_] := Module[{ds = degQ[S1, S2, mu], ode, qs},
   ode = gaudinODE /. {C2 -> Cas2[mu], C4 -> Cas4[mu], S -> S1, P -> S2} /. wsub;
   qs = Table[Module[{m = ds[[kk]], gauge = ds[[1 ;; kk - 1]], qq, cl, sol},
       qq = Sum[cq[j] z^j, {j, 0, m}];
       cl = CoefficientList[Expand[ode /. Q -> Function[x, qq /. z -> x]], z];
       sol = Solve[Join[Thread[cl == 0], {cq[m] == 1}, Thread[(cq /@ gauge) == 0]], Table[cq[j], {j, 0, m}]];
       Sum[(cq[j] /. sol[[1]]) z^j, {j, 0, m}]], {kk, 2}];
   Expand[qs[[1]] D[qs[[2]], z] - qs[[2]] D[qs[[1]], z]]];
det1two[S1_, S2_, mu_] := Module[{q = Q12two[S1, S2, mu]}, Det[{{doRes0[S1, S2, q], doRes0[S1, S2, z q]}, {doRes1[S1, S2, q], doRes1[S1, S2, z q]}}]/(q /. z -> 0)/(q /. z -> 1)/Exp[w]];
benchTwo[S1_, S2_] := Module[{ir = irreps[S1, S2], t}, memoizeRes; t = First@AbsoluteTiming[det1two[S1, S2, #] & /@ ir;]; t];
corr = Simplify[(det1[10, 10, #] - det1two[10, 10, #]) & /@ irreps[10, 10]] // Flatten // Union;
{"S=8_readable" -> benchTwo[8, 8], "S=10_readable" -> benchTwo[10, 10], "correctness10" -> corr}
```
**Output:**
```
Out= {"S=8_readable" -> 0.108973, "S=10_readable" -> 0.153538, "correctness10" -> {0}}
```
## 2026-07-22 07:58:51 — ✏️ EDIT CELL 6 [batch]
(* --- Gaudin ODE solver ---
   Arguments always ordered (S1,S2,mu): S1,S2 label the two input reps [S1,S1,0,0],[S2,S2,0,0]; mu
   is the output irrep. Solves the 4th-order Gaudin ODE for the four pol…


## 2026-07-22 07:58:52
**Input:**
```wolfram
(* --- Gaudin ODE solver ---
   Arguments always ordered (S1,S2,mu): S1,S2 label the two input reps [S1,S1,0,0],[S2,S2,0,0]; mu
   is the output irrep. Solves the 4th-order Gaudin ODE for the four polynomial Q's of degrees
   degQ[S1,S2,mu], returned as monic coefficient lists (coeffs of z^0,z^1,...,z^m). Optional kmax
   (default 4) returns only the first kmax Q's -- Q12/det1 need just Q1,Q2, and the higher Q's carry
   the largest degrees, so kmax=2 is much cheaper. GAUGE FIXING: the four Q-degrees are strictly
   increasing, so any lower-degree solution can be added to a higher one (q_k -> q_k + c q_j, j<k);
   removed by forcing the coefficient of z^(deg q_j) in q_k to vanish for every j<k. *)
ClearAll[GaudinSolve];
GaudinSolve[S1_, S2_, mu_] := GaudinSolve[S1, S2, mu, 4];
GaudinSolve[S1_, S2_, mu_, kmax_] := Module[{ds, ode},
   ds = degQ[S1, S2, mu];
   ode = gaudinODE /. {C2 -> Cas2[mu], C4 -> Cas4[mu], S -> S1, P -> S2} /. wsub;
   Table[Module[{m = ds[[kk]], gauge = ds[[1 ;; kk - 1]], q, cl, sol},
       q = Sum[cq[j] z^j, {j, 0, m}];
       cl = CoefficientList[Expand[ode /. Q -> Function[x, q /. z -> x]], z];
       sol = Solve[Join[Thread[cl == 0], {cq[m] == 1}, Thread[(cq /@ gauge) == 0]], Table[cq[j], {j, 0, m}]];
       Table[cq[j], {j, 0, m}] /. sol[[1]]], {kk, kmax}]];
```
**Output:**
```
Out= (no output)
```
## 2026-07-22 07:58:52 — ✏️ EDIT CELL 7 [batch]
(* --- BuildQs: same argument convention as GaudinSolve, (S1,S2,mu[,kmax]), but returns the Q's
   as explicit polynomials in z (rather than coefficient lists). kmax (default 4) limits to the
   first…


## 2026-07-22 07:58:52
**Input:**
```wolfram
(* --- BuildQs: same argument convention as GaudinSolve, (S1,S2,mu[,kmax]), but returns the Q's
   as explicit polynomials in z (rather than coefficient lists). kmax (default 4) limits to the
   first kmax Q's. *)
ClearAll[BuildQs];
BuildQs[S1_, S2_, mu_] := BuildQs[S1, S2, mu, 4];
BuildQs[S1_, S2_, mu_, kmax_] := Module[{ds = degQ[S1, S2, mu], coeffs = GaudinSolve[S1, S2, mu, kmax]},
   Table[Sum[coeffs[[kk, j + 1]] z^j, {j, 0, ds[[kk]]}], {kk, kmax}]];
```
**Output:**
```
Out= (no output)
```
## 2026-07-22 07:58:52 — ✏️ EDIT CELL 9 [batch]
(* Q12 = Wronskian of Q1,Q2 only, so build just those (kmax=2) -- skips the higher-degree Q3,Q4 solves. *)
ClearAll[Q12];
Q12[S1_, S2_, mu_] := Expand[Wronskian[BuildQs[S1, S2, mu, 2], z]];


## 2026-07-22 07:58:52
**Input:**
```wolfram
(* Q12 = Wronskian of Q1,Q2 only, so build just those (kmax=2) -- skips the higher-degree Q3,Q4 solves. *)
ClearAll[Q12];
Q12[S1_, S2_, mu_] := Expand[Wronskian[BuildQs[S1, S2, mu, 2], z]];
```
**Output:**
```
Out= (no output)
```
## 2026-07-22 07:58:52 — ✏️ EDIT CELL 11 [batch]
ClearAll[mm, res1, res0];
mm[z_, n_, S1_, S2_, w_] := Exp[w z] z^n / (z^(2 S1) (z - 1)^(2 S2));
(* res0/res1 memoized (f[x_]:=f[x]=...): for a fixed (S1,S2) sweep they are identical for every mu,
   a…


## 2026-07-22 07:58:52
**Input:**
```wolfram
ClearAll[mm, res1, res0];
mm[z_, n_, S1_, S2_, w_] := Exp[w z] z^n / (z^(2 S1) (z - 1)^(2 S2));
(* res0/res1 memoized (f[x_]:=f[x]=...): for a fixed (S1,S2) sweep they are identical for every mu,
   and each is a symbolic LaguerreL in w, so caching avoids recomputing them per mu. The ClearAll
   above clears the cache on re-run -- do re-run this cell if you ever change the measure/formulas. *)
res1[n_, S1_, S2_] := res1[n, S1, S2] = Exp[w] LaguerreL[2 S2 - 1, n - 2 S1 - 2 S2 + 1, -w];
res0[n_, S1_, S2_] := res0[n, S1, S2] = If[2 S1 > n, (-1)^(2 S1 + 2 S2 - n + 1) LaguerreL[2 S1 - n - 1, 1 - 2 S1 - 2 S2 + n, w], 0];
```
**Output:**
```
Out= (no output)
```

## 2026-07-22 07:59:48
**Input:**
```wolfram
(* cold cache: re-eval the memoized res defs to clear, then confirm values + speed *)
ClearAll[res0, res1];
res1[n_, S1_, S2_] := res1[n, S1, S2] = Exp[w] LaguerreL[2 S2 - 1, n - 2 S1 - 2 S2 + 1, -w];
res0[n_, S1_, S2_] := res0[n, S1, S2] = If[2 S1 > n, (-1)^(2 S1 + 2 S2 - n + 1) LaguerreL[2 S1 - n - 1, 1 - 2 S1 - 2 S2 + n, w], 0];
knownvals = Simplify[{det1[1, 1, {2, 1, 1, 0}], det1[1, 2, {3, 2, 1, 0}], det1[2, 1, {2, 2, 1, 1}]}];
t55 = First@AbsoluteTiming[det1[5, 5, #] & /@ irreps[5, 5];];
ClearAll[res0, res1];
res1[n_, S1_, S2_] := res1[n, S1, S2] = Exp[w] LaguerreL[2 S2 - 1, n - 2 S1 - 2 S2 + 1, -w];
res0[n_, S1_, S2_] := res0[n, S1, S2] = If[2 S1 > n, (-1)^(2 S1 + 2 S2 - n + 1) LaguerreL[2 S1 - n - 1, 1 - 2 S1 - 2 S2 + n, w], 0];
t1010 = First@AbsoluteTiming[det1[10, 10, #] & /@ irreps[10, 10];];
{"knownvals(expect {w^2, w^3(3+2w)/12, w^2(-40-4w+w^2)/6})" -> knownvals, "sweep55" -> t55, "sweep1010" -> t1010, "worstODE_still0" -> Max[Abs[Flatten[Table[Module[{qs = BuildQs[S1, S2, mu], ode = gaudinODE /. {C2 -> Cas2[mu], C4 -> Cas4[mu], S -> S1, P -> S2} /. wsub}, Table[Expand[ode /. Q -> Function[x, qs[[k]] /. z -> x]], {k, 4}]], {S1, 1, 2}, {S2, 1, 2}, {mu, irreps[S1, S2]}]] /. z -> 37]]}
```
**Output:**
```
Out= {"knownvals(expect {w^2, w^3(3+2w)/12, w^2(-40-4w+w^2)/6})" -> {w^2, (w^3*(3 + 2*w))/12, (w^2*(-40 - 4*w + w^2))/6}, "sweep55" -> 0.396498, "sweep1010" -> 6.985291, "worstODE_still0" -> 0}
```

## 2026-07-22 08:00:31
**Input:**
```wolfram
mu10 = irreps[10, 10][[Ceiling[Length[irreps[10, 10]]/2]]];
q10 = Q12[10, 10, mu10];
{"mu10" -> mu10, "degQ" -> degQ[10, 10, mu10], "degQ12" -> Exponent[q10, z],
 "t_Q12" -> First@AbsoluteTiming[Q12[10, 10, mu10];],
 "t_GaudinSolve_kmax2" -> First@AbsoluteTiming[GaudinSolve[10, 10, mu10, 2];],
 "t_GaudinSolve_default" -> First@AbsoluteTiming[GaudinSolve[10, 10, mu10];],
 "t_doRes0_cached" -> First@AbsoluteTiming[doRes0[10, 10, q10];],
 "t_doRes0_again" -> First@AbsoluteTiming[doRes0[10, 10, q10];],
 "res0_memoized?" -> (Length[DownValues[res0]] > 1),
 "t_det1_single" -> First@AbsoluteTiming[det1[10, 10, mu10];]}
```
**Output:**
```
Out= {"mu10" -> {17, 15, 5, 3}, "degQ" -> {3, 6, 17, 20}, "degQ12" -> 8, "t_Q12" -> 0.001642, "t_GaudinSolve_kmax2" -> 0.001213, "t_GaudinSolve_default" -> 0.003475, "t_doRes0_cached" -> 0.000086, "t_doRes0_again" -> 0.000039, "res0_memoized?" -> True, "t_det1_single" -> 0.03709}
```

## 2026-07-22 08:01:34
**Input:**
```wolfram
ir = irreps[10, 10];
(* warm the res cache first *)
det1[10, 10, #] & /@ ir;
(* now time: full det1 sweep, vs just the Q12 part, vs Q12+doRes (no Det/division) *)
tFull = First@AbsoluteTiming[det1[10, 10, #] & /@ ir;];
tQ12 = First@AbsoluteTiming[(Q12[10, 10, #];) & /@ ir;];
tQ12doRes = First@AbsoluteTiming[Module[{q}, q = Q12[10, 10, #]; {doRes0[10, 10, q], doRes0[10, 10, z q], doRes1[10, 10, q], doRes1[10, 10, z q]}] & /@ ir;];
(* the Det + division step alone, on precomputed entries *)
entries = Module[{q = Q12[10, 10, #]}, {{doRes0[10, 10, q], doRes0[10, 10, z q]}, {doRes1[10, 10, q], doRes1[10, 10, z q]}, q}] & /@ ir;
tDet = First@AbsoluteTiming[(Det[#[[1 ;; 2]]]/(#[[3]] /. z -> 0)/(#[[3]] /. z -> 1)/Exp[w]) & /@ entries;];
{"tFull" -> tFull, "tQ12" -> tQ12, "tQ12doRes" -> tQ12doRes, "tDet" -> tDet, "maxDegQ12" -> Max[Exponent[Q12[10, 10, #], z] & /@ ir]}
```
**Output:**
```
Out= {"tFull" -> 2.816065, "tQ12" -> 0.129974, "tQ12doRes" -> 0.20214, "tDet" -> 0.027263, "maxDegQ12" -> 20}
```

## 2026-07-22 08:02:17
**Input:**
```wolfram
ir = irreps[10, 10];
degs = Exponent[Q12[10, 10, #], z] & /@ ir;
muBig = ir[[First@Ordering[degs, -1]]];
q = Q12[10, 10, muBig];
a = First@AbsoluteTiming[d0 = doRes0[10, 10, q];];
b = First@AbsoluteTiming[d0z = doRes0[10, 10, z q];];
c = First@AbsoluteTiming[d1 = doRes1[10, 10, q];];
d = First@AbsoluteTiming[d1z = doRes1[10, 10, z q];];
e = First@AbsoluteTiming[dt = Det[{{d0, d0z}, {d1, d1z}}];];
f = First@AbsoluteTiming[res = dt/(q /. z -> 0)/(q /. z -> 1)/Exp[w];];
g = First@AbsoluteTiming[full = det1[10, 10, muBig];];
{"muBig" -> muBig, "degQ12" -> Exponent[q, z], "t_doRes0" -> a, "t_doRes0z" -> b, "t_doRes1" -> c, "t_doRes1z" -> d, "t_Det" -> e, "t_divide" -> f, "t_det1_full" -> g, "same?" -> Simplify[full - res] === 0}
```
**Output:**
```
Out= {"muBig" -> {10, 10, 10, 10}, "degQ12" -> 20, "t_doRes0" -> 0.000111, "t_doRes0z" -> 0.000086, "t_doRes1" -> 0.000134, "t_doRes1z" -> 0.00011, "t_Det" -> 0.000218, "t_divide" -> 0.000048, "t_det1_full" -> 0.060942, "same?" -> True}
```

## 2026-07-22 08:03:12
**Input:**
```wolfram
mub = {10, 10, 10, 10};
{"degQ" -> degQ[10, 10, mub],
 "len_BuildQs_2" -> Length[BuildQs[10, 10, mub, 2]],
 "len_GaudinSolve_2" -> Length[GaudinSolve[10, 10, mub, 2]],
 "t_Q12" -> Mean[Table[First@AbsoluteTiming[Q12[10, 10, mub];], {5}]],
 "t_Q12two" -> Mean[Table[First@AbsoluteTiming[Q12two[10, 10, mub];], {5}]],
 "t_GaudinSolve2" -> Mean[Table[First@AbsoluteTiming[GaudinSolve[10, 10, mub, 2];], {5}]],
 "t_GaudinSolve4" -> Mean[Table[First@AbsoluteTiming[GaudinSolve[10, 10, mub];], {5}]],
 "Q12_eq_Q12two" -> (Expand[Q12[10, 10, mub] - Q12two[10, 10, mub]] === 0)}
```
**Output:**
```
Out= {"degQ" -> {10, 11, 12, 13}, "len_BuildQs_2" -> 2, "len_GaudinSolve_2" -> 2, "t_Q12" -> 0.0055236, "t_Q12two" -> 0.0028986000000000003, "t_GaudinSolve2" -> 0.003104, "t_GaudinSolve4" -> 0.0037208, "Q12_eq_Q12two" -> True}
```

## 2026-07-22 08:03:39
**Input:**
```wolfram
ir = irreps[10, 10];
det1[10, 10, #] & /@ ir;  (* warm *)
{"det1_sweep_run1" -> First@AbsoluteTiming[det1[10, 10, #] & /@ ir;],
 "det1_sweep_run2" -> First@AbsoluteTiming[det1[10, 10, #] & /@ ir;],
 "det1two_sweep" -> First@AbsoluteTiming[det1two[10, 10, #] & /@ ir;],
 "N" -> Length[ir]}
```
**Output:**
```
Out= {"det1_sweep_run1" -> 2.863392, "det1_sweep_run2" -> 2.845619, "det1two_sweep" -> 0.137977, "N" -> 66}
```

## 2026-07-22 08:04:30
**Input:**
```wolfram
mub = {10, 10, 10, 10};
q1 = Q12[10, 10, mub];
q2 = Q12two[10, 10, mub];
{"identical" -> (q1 === q2),
 "leafCount" -> {LeafCount[q1], LeafCount[q2]},
 "head" -> {Head[q1], Head[q2]},
 "sameExpanded" -> (Expand[q1] === Expand[q2]),
 "coeffsSame" -> (CoefficientList[q1, z] === CoefficientList[q2, z]),
 "t_doRes0_on_q1" -> Mean[Table[First@AbsoluteTiming[doRes0[10, 10, q1];], {5}]],
 "t_doRes0_on_q2" -> Mean[Table[First@AbsoluteTiming[doRes0[10, 10, q2];], {5}]],
 "t_full_inline_Q12" -> First@AbsoluteTiming[Module[{q = Q12[10, 10, mub]}, Det[{{doRes0[10, 10, q], doRes0[10, 10, z q]}, {doRes1[10, 10, q], doRes1[10, 10, z q]}}]/(q /. z -> 0)/(q /. z -> 1)/Exp[w]];],
 "t_full_precomp_q1" -> First@AbsoluteTiming[Module[{q = q1}, Det[{{doRes0[10, 10, q], doRes0[10, 10, z q]}, {doRes1[10, 10, q], doRes1[10, 10, z q]}}]/(q /. z -> 0)/(q /. z -> 1)/Exp[w]];]}
```
**Output:**
```
Out= {"identical" -> True, "leafCount" -> {128, 128}, "head" -> {Plus, Plus}, "sameExpanded" -> True, "coeffsSame" -> True, "t_doRes0_on_q1" -> 0.0001788, "t_doRes0_on_q2" -> 0.000154, "t_full_inline_Q12" -> 0.007576, "t_full_precomp_q1" -> 0.001127}
```

## 2026-07-22 08:05:04
**Input:**
```wolfram
ir = irreps[10, 10];
tim = Table[{mu, First@AbsoluteTiming[det1[10, 10, mu];], First@AbsoluteTiming[Q12[10, 10, mu];], degQ[10, 10, mu][[1 ;; 2]]}, {mu, ir}];
slow = TakeLargest[tim, 6, 2] // (Function[t, {#[[1]], "t_det1" -> #[[2]], "t_Q12" -> #[[3]], "degQ12" -> #[[4]]} & /@ t]);
{"total_det1" -> Total[tim[[All, 2]]], "total_Q12" -> Total[tim[[All, 3]]], "slowest6" -> slow}
```
**Output:**
```
[message] TakeLargest::nonopt: Options expected (instead of 2) beyond position 2 in TakeLargest[{{{20, 20, 0, 0}, 0.004294, 0.000782, {0, 1}}, {{20, 19, 1, 0}, 0.005275, 0.001029, {0, 2}}, {{20, 18, 2, 0}, 0.00… [3193 chars]
[message] Part::partd: Part specification 6[[1]] is longer than depth of object.
[message] [Part::partd ×2 more suppressed — first instance shown above]
[message] General::stop: Further output of Part::partd will be suppressed during this calculation.
[message] [TakeLargest::nonopt ×1 more suppressed — first instance shown above]
Out= {"total_det1" -> 2.7965240000000002, "total_Q12" -> 0.16771999999999998, "slowest6" -> TakeLargest[{{{20, 20, 0, 0}, 0.004294, 0.000782, {0, 1}}, "t_det1" -> {{20, 19, 1, 0}, 0.005275, 0.001029, {0, 2}}, "t_Q12" -> {{20, 18, 2, 0}, 0.006743, 0.001186, {0, 3}}, "degQ12" -> {{20, 17, 3, 0}, 0.020495, 0.001332, {0, 4}}}, {6[[1]], "t_det1" -> 6[[2]], "t_Q12" -> 6[[3]], "degQ12" -> 6[[4]]}, {2[[1]], "t_det1" -> 2[[2]], "t_Q12" -> 2[[3]], "degQ12" -> 2[[4]]}]}
```

## 2026-07-22 08:06:19
**Input:**
```wolfram
(* redefine Q12 to the inline 2-Q solve (Q12two body); re-time det1 sweep *)
ClearAll[Q12];
Q12[S1_, S2_, mu_] := Module[{ds = degQ[S1, S2, mu], ode, qs},
   ode = gaudinODE /. {C2 -> Cas2[mu], C4 -> Cas4[mu], S -> S1, P -> S2} /. wsub;
   qs = Table[Module[{m = ds[[kk]], gauge = ds[[1 ;; kk - 1]], qq, cl, sol},
       qq = Sum[cq[j] z^j, {j, 0, m}];
       cl = CoefficientList[Expand[ode /. Q -> Function[x, qq /. z -> x]], z];
       sol = Solve[Join[Thread[cl == 0], {cq[m] == 1}, Thread[(cq /@ gauge) == 0]], Table[cq[j], {j, 0, m}]];
       Sum[(cq[j] /. sol[[1]]) z^j, {j, 0, m}]], {kk, 2}];
   Expand[qs[[1]] D[qs[[2]], z] - qs[[2]] D[qs[[1]], z]]];
ir = irreps[10, 10];
det1[10, 10, #] & /@ ir;
{"det1_sweep_now" -> First@AbsoluteTiming[det1[10, 10, #] & /@ ir;],
 "det1_5,5" -> First@AbsoluteTiming[det1[5, 5, #] & /@ irreps[5, 5];],
 "det1_8,8" -> First@AbsoluteTiming[det1[8, 8, #] & /@ irreps[8, 8];]}
```
**Output:**
```
Out= {"det1_sweep_now" -> 2.778754, "det1_5,5" -> 0.311915, "det1_8,8" -> 2.754547}
```

## 2026-07-22 08:06:54
**Input:**
```wolfram
Off[General::spell];
{"DV_det1" -> (DownValues[det1] // InputForm // ToString),
 "DV_det1two" -> (DownValues[det1two] // InputForm // ToString),
 "res0_memoized_entries" -> Length[DownValues[res0]],
 "res1_memoized_entries" -> Length[DownValues[res1]],
 "t_det1_88" -> First@AbsoluteTiming[det1[8, 8, #] & /@ irreps[8, 8];],
 "t_det1two_88" -> First@AbsoluteTiming[det1two[8, 8, #] & /@ irreps[8, 8];]}
```
**Output:**
```
Out= {"DV_det1" -> "{HoldPattern[det1[S1_, S2_, mu_]] :> Module[{q = Q12[S1, S2, mu]}, ExpandAll[Det[{{doRes0[S1, S2, q], doRes0[S1, S2, z*q]}, {doRes1[S1, S2, q], doRes1[S1, S2, z*q]}}]/(q /. z -> 0)/(q /. z -> 1)/Exp[w]]]}", "DV_det1two" -> "{HoldPattern[det1two[S1_, S2_, mu_]] :> Module[{q = Q12two[S1, S2, mu]}, Det[{{doRes0[S1, S2, q], doRes0[S1, S2, z*q]}, {doRes1[S1, S2, q], doRes1[S1, S2, z*q]}}]/(q /. z -> 0)/(q /. z -> 1)/Exp[w]]}", "res0_memoized_entries" -> 53, "res1_memoized_entries" -> 53, "t_det1_88" -> 2.760399, "t_det1two_88" -> 1.671012}
```

## 2026-07-22 08:08:10
**Input:**
```wolfram
(* res1 without the Exp[w] factor -> br1 "bar"; det1 with no Exp[w] at all *)
ClearAll[res1bar, doRes1bar];
res1bar[n_, S1_, S2_] := res1bar[n, S1, S2] = LaguerreL[2 S2 - 1, n - 2 S1 - 2 S2 + 1, -w];
doRes1bar[S1_, S2_, f_] := Module[{nmax = Exponent[f, z], c = CoefficientList[f, z]}, Plus @@ (c Table[res1bar[n, S1, S2], {n, 0, nmax}])];
det1A[S1_, S2_, mu_] := Module[{q = Q12[S1, S2, mu]}, ExpandAll[Det[{{doRes0[S1, S2, q], doRes0[S1, S2, z q]}, {doRes1[S1, S2, q], doRes1[S1, S2, z q]}}]/(q /. z -> 0)/(q /. z -> 1)/Exp[w]]];
det1B[S1_, S2_, mu_] := Module[{q = Q12[S1, S2, mu]}, Expand[Det[{{doRes0[S1, S2, q], doRes0[S1, S2, z q]}, {doRes1[S1, S2, q], doRes1[S1, S2, z q]}}]/(q /. z -> 0)/(q /. z -> 1)]/Exp[w] // Expand];
det1C[S1_, S2_, mu_] := Module[{q = Q12[S1, S2, mu], a0, a0z, b1, b1z}, a0 = doRes0[S1, S2, q]; a0z = doRes0[S1, S2, z q]; b1 = doRes1bar[S1, S2, q]; b1z = doRes1bar[S1, S2, z q]; Expand[(a0 b1z - a0z b1)/((q /. z -> 0) (q /. z -> 1))]];
ir = irreps[10, 10];
{"consistent_AC" -> (Simplify[(det1A[10, 10, #] - det1C[10, 10, #]) & /@ ir] // Flatten // Union),
 "t_A_ExpandAll" -> First@AbsoluteTiming[det1A[10, 10, #] & /@ ir;],
 "t_C_ExpFree" -> First@AbsoluteTiming[det1C[10, 10, #] & /@ ir;],
 "t_C_88" -> First@AbsoluteTiming[det1C[8, 8, #] & /@ irreps[8, 8];],
 "t_C_55" -> First@AbsoluteTiming[det1C[5, 5, #] & /@ irreps[5, 5];]}
```
**Output:**
```
Out= {"consistent_AC" -> {0}, "t_A_ExpandAll" -> 2.776189, "t_C_ExpFree" -> 0.197659, "t_C_88" -> 1.739919, "t_C_55" -> 0.276142}
```

## 2026-07-22 08:10:00
**Input:**
```wolfram
det1D[S1_, S2_, mu_] := Module[{q = Q12[S1, S2, mu]}, Expand[Det[{{doRes0[S1, S2, q], doRes0[S1, S2, z q]}, {doRes1[S1, S2, q], doRes1[S1, S2, z q]}}]/(q /. z -> 0)/(q /. z -> 1)/Exp[w]]];
warm := (det1D[5, 5, #] & /@ irreps[5, 5]; det1D[8, 8, #] & /@ irreps[8, 8]; det1C[8, 8, #] & /@ irreps[8, 8]; det1D[10, 10, #] & /@ irreps[10, 10]; det1C[10, 10, #] & /@ irreps[10, 10];);
warm;
{"consistent_CD_10" -> (Simplify[(det1C[10, 10, #] - det1D[10, 10, #]) & /@ irreps[10, 10]] // Flatten // Union),
 "D_Expand_55" -> Min[Table[First@AbsoluteTiming[det1D[5, 5, #] & /@ irreps[5, 5];], {3}]],
 "D_Expand_88" -> Min[Table[First@AbsoluteTiming[det1D[8, 8, #] & /@ irreps[8, 8];], {3}]],
 "D_Expand_1010" -> Min[Table[First@AbsoluteTiming[det1D[10, 10, #] & /@ irreps[10, 10];], {3}]],
 "C_ExpFree_88" -> Min[Table[First@AbsoluteTiming[det1C[8, 8, #] & /@ irreps[8, 8];], {3}]],
 "C_ExpFree_1010" -> Min[Table[First@AbsoluteTiming[det1C[10, 10, #] & /@ irreps[10, 10];], {3}]]}
```
**Output:**
```
Out= {"consistent_CD_10" -> {0}, "D_Expand_55" -> 0.297899, "D_Expand_88" -> 1.85353, "D_Expand_1010" -> 4.893411, "C_ExpFree_88" -> 1.723334, "C_ExpFree_1010" -> 4.428173}
```

## 2026-07-22 08:10:43
**Input:**
```wolfram
(* fixed hard input: highest-degree Q12 in irreps[10,10] *)
ir = irreps[10, 10]; degs = Exponent[Q12[10, 10, #], z] & /@ ir; muH = ir[[First@Ordering[degs, -1]]];
qH = Q12[10, 10, muH]; n0 = qH /. z -> 0; n1 = qH /. z -> 1;
d0 = doRes0[10, 10, qH]; d0z = doRes0[10, 10, z qH]; d1 = doRes1[10, 10, qH]; d1z = doRes1[10, 10, z qH]; d1b = doRes1bar[10, 10, qH]; d1bz = doRes1bar[10, 10, z qH];
{"muH" -> muH, "degQ12" -> Exponent[qH, z],
 "ExpandAll" -> First@RepeatedTiming[ExpandAll[Det[{{d0, d0z}, {d1, d1z}}]/n0/n1/Exp[w]]],
 "Expand_keepExp" -> First@RepeatedTiming[Expand[Det[{{d0, d0z}, {d1, d1z}}]/n0/n1/Exp[w]]],
 "Expand_ExpFree" -> First@RepeatedTiming[Expand[(d0 d1bz - d0z d1b)/(n0 n1)]],
 "sameResult" -> (Expand[ExpandAll[Det[{{d0, d0z}, {d1, d1z}}]/n0/n1/Exp[w]] - (d0 d1bz - d0z d1b)/(n0 n1)] === 0)}
```
**Output:**
```
Out= {"muH" -> {10, 10, 10, 10}, "degQ12" -> 20, "ExpandAll" -> 0.0601916875, "Expand_keepExp" -> 0.007782, "Expand_ExpFree" -> 0.001687529296875, "sameResult" -> True}
```

## 2026-07-22 08:13:09
**Input:**
```wolfram
ClearAll[det1F, res1bar, doRes1bar];
res1bar[n_, S1_, S2_] := res1bar[n, S1, S2] = LaguerreL[2 S2 - 1, n - 2 S1 - 2 S2 + 1, -w];
doRes1bar[S1_, S2_, f_] := Module[{c = CoefficientList[f, z]}, c . Table[res1bar[n, S1, S2], {n, 0, Length[c] - 1}]];
det1F[S1_, S2_, mu_] := Module[{q = Q12[S1, S2, mu], n0, n1}, n0 = q /. z -> 0; n1 = q /. z -> 1; Expand[(doRes0[S1, S2, q] doRes1bar[S1, S2, z q] - doRes0[S1, S2, z q] doRes1bar[S1, S2, q])/(n0 n1)]];
(* current notebook version for comparison *)
det1Now[S1_, S2_, mu_] := Module[{q = Q12[S1, S2, mu]}, ExpandAll[Det[{{doRes0[S1, S2, q], doRes0[S1, S2, z q]}, {doRes1[S1, S2, q], doRes1[S1, S2, z q]}}]/(q /. z -> 0)/(q /. z -> 1)/Exp[w]]];
warm := (det1F[#[[1]], #[[2]], #3] & ; Do[det1F[s, s, mu]; det1Now[s, s, mu], {s, {5, 8, 10}}, {mu, irreps[s, s]}];);
warm;
bench[s_] := {"S" -> s, "N" -> Length[irreps[s, s]],
   "now_ExpandAll" -> Min[Table[First@AbsoluteTiming[det1Now[s, s, #] & /@ irreps[s, s];], {3}]],
   "fast_ExpFree" -> Min[Table[First@AbsoluteTiming[det1F[s, s, #] & /@ irreps[s, s];], {3}]],
   "correct" -> (Simplify[(det1F[s, s, #] - det1Now[s, s, #]) & /@ irreps[s, s]] // Flatten // Union)};
{bench[5], bench[8], bench[10]}
```
**Output:**
```
Out= {{"S" -> 5, "N" -> 21, "now_ExpandAll" -> 0.312507, "fast_ExpFree" -> 0.271882, "correct" -> {0}}, {"S" -> 8, "N" -> 45, "now_ExpandAll" -> 2.705339, "fast_ExpFree" -> 1.666051, "correct" -> {0}}, {"S" -> 10, "N" -> 66, "now_ExpandAll" -> 6.996706, "fast_ExpFree" -> 4.447583, "correct" -> {0}}}
```

## 2026-07-22 08:13:41
**Input:**
```wolfram
ir = irreps[10, 10];
prof = Table[Module[{tq, trest, q},
    tq = First@AbsoluteTiming[q = Q12[10, 10, mu]];
    trest = First@AbsoluteTiming[Expand[(doRes0[10, 10, q] doRes1bar[10, 10, z q] - doRes0[10, 10, z q] doRes1bar[10, 10, q])/((q /. z -> 0) (q /. z -> 1))]];
    {mu, Exponent[q, z], tq, trest}], {mu, ir}];
top = Take[SortBy[prof, -(#[[4]]) &], 5];
{"sum_Q12" -> Total[prof[[All, 3]]], "sum_rest" -> Total[prof[[All, 4]]],
 "top5_by_rest" -> ({"mu" -> #[[1]], "degQ12" -> #[[2]], "tQ12" -> #[[3]], "tRest" -> #[[4]]} & /@ top)}
```
**Output:**
```
Out= {"sum_Q12" -> 0.092407, "sum_rest" -> 0.09892800000000002, "top5_by_rest" -> {{"mu" -> {13, 10, 10, 7}, "degQ12" -> 17, "tQ12" -> 0.00333, "tRest" -> 0.003245}, {"mu" -> {13, 11, 9, 7}, "degQ12" -> 16, "tQ12" -> 0.001818, "tRest" -> 0.002984}, {"mu" -> {16, 10, 10, 4}, "degQ12" -> 14, "tQ12" -> 0.001726, "tRest" -> 0.002353}, {"mu" -> {13, 12, 8, 7}, "degQ12" -> 15, "tQ12" -> 0.001767, "tRest" -> 0.002267}, {"mu" -> {12, 11, 9, 8}, "degQ12" -> 17, "tQ12" -> 0.00177, "tRest" -> 0.002176}}}
```

## 2026-07-22 08:14:16
**Input:**
```wolfram
ir = irreps[10, 10];
Do[det1F[10, 10, mu], {mu, ir}];  (* warm *)
{"map_suppressed" -> First@AbsoluteTiming[det1F[10, 10, #] & /@ ir;],
 "map_toNull" -> First@AbsoluteTiming[Do[det1F[10, 10, mu], {mu, ir}]],
 "sum_individual" -> Total[Table[First@AbsoluteTiming[det1F[10, 10, mu]], {mu, ir}]],
 "scan" -> First@AbsoluteTiming[Scan[det1F[10, 10, #] &, ir]]}
```
**Output:**
```
Out= {"map_suppressed" -> 0.188634, "map_toNull" -> 0.199617, "sum_individual" -> 0.194058, "scan" -> 0.1903}
```

## 2026-07-22 08:14:44
**Input:**
```wolfram
bench2[s_] := Module[{ir = irreps[s, s]},
   Do[det1Now[s, s, mu]; det1F[s, s, mu], {mu, ir}];  (* warm both *)
   {"S" -> s, "N" -> Length[ir],
    "now_ExpandAll" -> Total[Table[First@AbsoluteTiming[det1Now[s, s, mu]], {mu, ir}]],
    "fast_ExpFree" -> Total[Table[First@AbsoluteTiming[det1F[s, s, mu]], {mu, ir}]]}];
Prepend[Association /@ (bench2 /@ {5, 8, 10}), "speedups"]
```
**Output:**
```
Out= {"speedups", <|"S" -> 5, "N" -> 21, "now_ExpandAll" -> 0.078458, "fast_ExpFree" -> 0.035533999999999996|>, <|"S" -> 8, "N" -> 45, "now_ExpandAll" -> 1.146477, "fast_ExpFree" -> 0.111306|>, <|"S" -> 10, "N" -> 66, "now_ExpandAll" -> 2.797767, "fast_ExpFree" -> 0.19264799999999999|>}
```

## 2026-07-22 08:15:46
**Input:**
```wolfram
det1Dv[S1_, S2_, mu_] := Module[{q = Q12[S1, S2, mu]}, Expand[Det[{{doRes0[S1, S2, q], doRes0[S1, S2, z q]}, {doRes1[S1, S2, q], doRes1[S1, S2, z q]}}]/(q /. z -> 0)/(q /. z -> 1)/Exp[w]]];
cmp[s_] := Module[{ir = irreps[s, s]}, Do[det1Dv[s, s, mu]; det1F[s, s, mu], {mu, ir}];
   <|"S" -> s, "Expand_keepExp" -> Total[Table[First@AbsoluteTiming[det1Dv[s, s, mu]], {mu, ir}]],
     "ExpFree" -> Total[Table[First@AbsoluteTiming[det1F[s, s, mu]], {mu, ir}]],
     "same" -> (Simplify[(det1Dv[s, s, #] - det1F[s, s, #]) & /@ ir] // Flatten // Union)|>];
cmp /@ {8, 10}
```
**Output:**
```
Out= {<|"S" -> 8, "Expand_keepExp" -> 0.285939, "ExpFree" -> 0.11144300000000001, "same" -> {0}|>, <|"S" -> 10, "Expand_keepExp" -> 0.659085, "ExpFree" -> 0.20218099999999997, "same" -> {0}|>}
```
## 2026-07-22 08:16:53 — ✏️ EDIT CELL 14
(* det1: the br[1] entries carry a common factor Exp[w] (from res1) that factors out of the 2x2
   determinant and cancels the overall 1/Exp[w], leaving a plain polynomial in w. Stripping that
   Exp[…

## 2026-07-22 08:17:02 — ✏️ EDIT CELL 13
## `det1` -- 2x2 residue determinant

$$\det{}_1[S_1,S_2,\mu]=\frac{1}{Q_{12}(0)\,Q_{12}(1)\,e^{w}}\;\det\begin{pmatrix} \mathrm{br}_0[e^{wz}Q_{12}(z)] & \mathrm{br}_0[e^{wz}\,z\,Q_{12}(z)]\\ \mathrm{…

## 2026-07-22 08:17:12 — ▶️ RUN CELLS 1–14
Cell 3: ✓ — (no output)
Cell 4: ✓ — (no output)
Cell 5: ✓ — (no output)
Cell 6: ✓ — (no output)
Cell 7: ✓ — (no output)
Cell 9: ✓ — (no output)
Cell 11: ✓ — (no output)
Cell 12: ✓ — (no output)
Cell 14: ✓ — (no output)


## 2026-07-22 08:17:33
**Input:**
```wolfram
(* verify notebook det1 (uses notebook Q12 = BuildQs/Wronskian): known values + warm sweep timing *)
known = Simplify[{det1[1, 1, {2, 1, 1, 0}], det1[1, 2, {3, 2, 1, 0}], det1[2, 1, {2, 2, 1, 1}]}];
Do[det1[s, s, mu], {s, {5, 8, 10}}, {mu, irreps[s, s]}];  (* warm *)
t = Association @@ Table[s -> Total[Table[First@AbsoluteTiming[det1[s, s, mu]], {mu, irreps[s, s]}]], {s, {5, 8, 10}}];
{"known (expect {w^2, w^3(3+2w)/12, w^2(-40-4w+w^2)/6})" -> known, "sweep_times(s)" -> t, "user_cell19" -> (2540160000 det1[5, 5, {10, 8, 2, 0}] // Expand)}
```
**Output:**
```
Out= {"known (expect {w^2, w^3(3+2w)/12, w^2(-40-4w+w^2)/6})" -> {w^2, (w^3*(3 + 2*w))/12, (w^2*(-40 - 4*w + w^2))/6}, "sweep_times(s)" -> <|5 -> 0.060633, 8 -> 0.189363, 10 -> 0.369168|>, "user_cell19" -> w^16 + (25*w^18)/1296}
```

## 2026-07-22 08:18:10
**Input:**
```wolfram
Q12wr[S1_, S2_, mu_] := Expand[Wronskian[BuildQs[S1, S2, mu, 2], z]];
Q12man[S1_, S2_, mu_] := Module[{qs = BuildQs[S1, S2, mu, 2]}, Expand[qs[[1]] D[qs[[2]], z] - qs[[2]] D[qs[[1]], z]]];
ir = irreps[10, 10];
Do[Q12wr[10, 10, mu]; Q12man[10, 10, mu], {mu, ir}];
{"Q12_Wronskian_builtin" -> Total[Table[First@AbsoluteTiming[Q12wr[10, 10, mu]], {mu, ir}]],
 "Q12_manual_D" -> Total[Table[First@AbsoluteTiming[Q12man[10, 10, mu]], {mu, ir}]],
 "same" -> (Simplify[(Q12wr[10, 10, #] - Q12man[10, 10, #]) & /@ ir] // Flatten // Union)}
```
**Output:**
```
Out= {"Q12_Wronskian_builtin" -> 0.10876799999999999, "Q12_manual_D" -> 0.089402, "same" -> {0}}
```

---

## 2026-07-22 08:28:18 — 🔄 KERNEL RESTART


## 2026-07-22 08:37:02
**Input:**
```wolfram
clearAll3 := (ClearAll[res0, res1, res1bar];
   res1[n_, S1_, S2_] := res1[n, S1, S2] = Exp[w] LaguerreL[2 S2 - 1, n - 2 S1 - 2 S2 + 1, -w];
   res0[n_, S1_, S2_] := res0[n, S1, S2] = If[2 S1 > n, (-1)^(2 S1 + 2 S2 - n + 1) LaguerreL[2 S1 - n - 1, 1 - 2 S1 - 2 S2 + n, w], 0];
   res1bar[n_, S1_, S2_] := res1bar[n, S1, S2] = LaguerreL[2 S2 - 1, n - 2 S1 - 2 S2 + 1, -w];);
coldwarm[s_] := Module[{ir = irreps[s, s], tc, tw},
   clearAll3;
   tc = First@AbsoluteTiming[det1[s, s, #] & /@ ir;];
   tw = First@AbsoluteTiming[det1[s, s, #] & /@ ir;];
   {"S" -> s, "cold" -> tc, "warm" -> tw, "penalty" -> tc - tw}];
{coldwarm[8], coldwarm[10], coldwarm[12]}
```
**Output:**
```
Out= {{"S" -> 8, "cold" -> 0.190418, "warm" -> 0.184384, "penalty" -> 0.006034000000000012}, {"S" -> 10, "cold" -> 3.592713, "warm" -> 0.852881, "penalty" -> 2.739832}, {"S" -> 12, "cold" -> 6.865884, "warm" -> 2.600185, "penalty" -> 4.265699}}
```

## 2026-07-22 08:37:48
**Input:**
```wolfram
clearAll3;
(* distinct residues needed for the S=10 sweep, and cost to build the tables *)
ir = irreps[10, 10]; nmax = Max[Exponent[z Q12[10, 10, #], z] & /@ ir];
{"nmax" -> nmax,
 "t_res1bar_table" -> First@AbsoluteTiming[Table[res1bar[n, 10, 10], {n, 0, nmax}];],
 "t_res0_table" -> First@AbsoluteTiming[Table[res0[n, 10, 10], {n, 0, nmax}];],
 "t_single_LagD19_am39" -> First@AbsoluteTiming[LaguerreL[19, -39, -w]],
 "t_single_LagD19_am20" -> First@AbsoluteTiming[LaguerreL[19, -20, -w]],
 "t_single_LagD10_am30" -> First@AbsoluteTiming[LaguerreL[10, -30, w]],
 "sumFormula_vs_builtin_same" -> (Expand[LaguerreL[19, -25, -w] - Sum[(-1)^k Binomial[19 - 25, 19 - k] (-w)^k/k!, {k, 0, 19}]] === 0),
 "t_sumFormula_D19" -> First@AbsoluteTiming[Expand[Sum[(-1)^k Binomial[19 - 25, 19 - k] (-w)^k/k!, {k, 0, 19}]]]}
```
**Output:**
```
Out= {"nmax" -> 21, "t_res1bar_table" -> 0.004473, "t_res0_table" -> 0.002406, "t_single_LagD19_am39" -> 0.000223, "t_single_LagD19_am20" -> 0.000193, "t_single_LagD10_am30" -> 0.000111, "sumFormula_vs_builtin_same" -> True, "t_sumFormula_D19" -> 0.000074}
```

## 2026-07-22 08:38:29
**Input:**
```wolfram
ir = irreps[10, 10];
(* 1) cold: clear, immediately sweep *)
clearAll3; tCold = First@AbsoluteTiming[det1[10, 10, #] & /@ ir;];
(* 2) precompute: clear, fill res tables for n=0..nmax, then sweep *)
clearAll3;
nmax = 25;
tPre = First@AbsoluteTiming[Table[res0[n, 10, 10]; res1bar[n, 10, 10]; res1[n, 10, 10], {n, 0, nmax}];];
tSweepAfterPre = First@AbsoluteTiming[det1[10, 10, #] & /@ ir;];
(* 3) warm: sweep again *)
tWarm = First@AbsoluteTiming[det1[10, 10, #] & /@ ir;];
{"tCold" -> tCold, "tPrecompute" -> tPre, "tSweepAfterPrecompute" -> tSweepAfterPre, "tWarm" -> tWarm,
 "res0_cached_count" -> Length[DownValues[res0]], "res1bar_cached_count" -> Length[DownValues[res1bar]]}
```
**Output:**
```
Out= {"tCold" -> 0.546509, "tPrecompute" -> 0.012391, "tSweepAfterPrecompute" -> 0.411902, "tWarm" -> 0.399906, "res0_cached_count" -> 27, "res1bar_cached_count" -> 27}
```

## 2026-07-22 08:40:07
**Input:**
```wolfram
trial[s_] := Module[{ir = irreps[s, s], lazy, pre, sweepPre, nmax = 30},
   (* lazy first run *)
   clearAll3; lazy = First@AbsoluteTiming[det1[s, s, #] & /@ ir;];
   (* precompute then run *)
   clearAll3;
   pre = First@AbsoluteTiming[Do[res0[n, s, s]; res1bar[n, s, s], {n, 0, nmax}]];
   sweepPre = First@AbsoluteTiming[det1[s, s, #] & /@ ir;];
   {"lazy_firstrun" -> lazy, "precompute" -> pre, "sweep_after_pre" -> sweepPre, "precompute_total" -> pre + sweepPre}];
(* min over 4 trials to cut noise *)
summarize[s_] := Module[{r = Table[trial[s], {4}]},
   <|"S" -> s,
     "lazy_firstrun" -> Min[#["lazy_firstrun"] & /@ (Association /@ r)],
     "precompute_total" -> Min[#["precompute_total"] & /@ (Association /@ r)],
     "precompute_step_only" -> Min[#["precompute"] & /@ (Association /@ r)]|>];
{summarize[10], summarize[12]}
```
**Output:**
```
Out= {<|"S" -> 10, "lazy_firstrun" -> 0.39495, "precompute_total" -> 0.389238, "precompute_step_only" -> 0.008264|>, <|"S" -> 12, "lazy_firstrun" -> 0.678363, "precompute_total" -> 0.6784519999999999, "precompute_step_only" -> 0.01078|>}
```

## 2026-07-22 08:43:19
**Input:**
```wolfram
(* profile GaudinSolve internals for the hardest Q2 at S=12 *)
s = 12; ir = irreps[s, s];
mu = ir[[First@Ordering[degQ[s, s, #][[2]] & /@ ir, -1]]];  (* max deg Q2 *)
ode = gaudinODE /. {C2 -> Cas2[mu], C4 -> Cas4[mu], S -> s, P -> s} /. wsub;
ds = degQ[s, s, mu]; m = ds[[2]]; gauge = {ds[[1]]};
cvars = Table[cq[j], {j, 0, m}];
tExpand = First@AbsoluteTiming[cl = CoefficientList[Expand[ode /. Q -> Function[x, (cvars . Table[z^j, {j, 0, m}]) /. z -> x]], z];];
tSolve = First@AbsoluteTiming[sol = Solve[Join[Thread[cl == 0], {cq[m] == 1}, Thread[(cq /@ gauge) == 0]], cvars];];
tCoeffArr = First@AbsoluteTiming[mat = CoefficientArrays[cl, cvars][[2]];];
{"mu" -> mu, "degQ2" -> m, "numEqs" -> Length[cl], "t_Expand+CoeffList" -> tExpand, "t_Solve" -> tSolve, "t_CoefficientArrays" -> tCoeffArr, "matDims" -> Dimensions[mat]}
```
**Output:**
```
Out= {"mu" -> {12, 12, 12, 12}, "degQ2" -> 13, "numEqs" -> 12, "t_Expand+CoeffList" -> 0.000563, "t_Solve" -> 0.001365, "t_CoefficientArrays" -> 0.000031, "matDims" -> {12, 14}}
```

## 2026-07-22 08:44:06
**Input:**
```wolfram
ClearAll[GaudinSolveLS];
GaudinSolveLS[S1_, S2_, mu_, kmax_] := Module[{ds, ode},
   ds = degQ[S1, S2, mu];
   ode = gaudinODE /. {C2 -> Cas2[mu], C4 -> Cas4[mu], S -> S1, P -> S2} /. wsub;
   Table[Module[{m = ds[[kk]], gauge = ds[[1 ;; kk - 1]], cvars, cl, mat, freeidx, cfree, cvals},
       cvars = Table[cq[j], {j, 0, m}];
       cl = CoefficientList[Expand[ode /. Q -> Function[x, (cvars . Table[z^j, {j, 0, m}]) /. z -> x]], z];
       mat = CoefficientArrays[cl, cvars][[2]];
       freeidx = Complement[Range[0, m], Append[gauge, m]];
       cfree = LinearSolve[mat[[All, freeidx + 1]], -mat[[All, m + 1]]];
       cvals = ConstantArray[0, m + 1]; cvals[[m + 1]] = 1;
       MapThread[(cvals[[#1 + 1]] = #2) &, {freeidx, cfree}];
       cvals], {kk, kmax}]];
(* correctness vs Solve-based GaudinSolve, over full sweeps *)
chk[s_] := (GaudinSolve[s, s, #, 2] - GaudinSolveLS[s, s, #, 2]) & /@ irreps[s, s] // Flatten // Union;
{"match_S8" -> chk[8], "match_S10" -> chk[10], "match_S12" -> chk[12]}
```
**Output:**
```
[message] Part::partw: Part 2 of {{}} does not exist.
[message] Part::pkspec1: The expression Integer[] cannot be used as a part specification.
[message] Part::partd: Part specification {{}}[[2]][[All,1]] is longer than depth of object.
[message] MapThread::mptd: Object LinearSolve[{}[[Integer[]]], -{{}}[[2]][[All,1]]] at position {2, 2} in MapThread[(cvals$218540[[#1 + 1]] = #2) & , {{}, LinearSolve[{}[[Integer[]]], -{{}}[[2]][[All,1]]]}] has… [233 chars]
[message] [Part::partw ×2 more suppressed — first instance shown above]
[message] [Part::pkspec1 ×2 more suppressed — first instance shown above]
[message] General::stop: Further output of Part::partw will be suppressed during this calculation.
[message] [MapThread::mptd ×2 more suppressed — first instance shown above]
[message] [General::stop ×3 more suppressed — first instance shown above]
[message] [Part::partd ×2 more suppressed — first instance shown above]
Out= {"match_S8" -> {0}, "match_S10" -> {0}, "match_S12" -> {0}}
```

## 2026-07-22 08:45:22
**Input:**
```wolfram
ClearAll[GaudinSolveLS, Q12LS, det1LS];
GaudinSolveLS[S1_, S2_, mu_, kmax_] := Module[{ds, ode},
   ds = degQ[S1, S2, mu];
   ode = gaudinODE /. {C2 -> Cas2[mu], C4 -> Cas4[mu], S -> S1, P -> S2} /. wsub;
   Table[Module[{m = ds[[kk]], gauge = ds[[1 ;; kk - 1]], cvars, cl, mat, freeidx, cvals},
       cvars = Table[cq[j], {j, 0, m}]; cvals = ConstantArray[0, m + 1]; cvals[[m + 1]] = 1;
       freeidx = Complement[Range[0, m], Append[gauge, m]];
       If[freeidx =!= {},
         cl = CoefficientList[Expand[ode /. Q -> Function[x, (cvars . Table[z^j, {j, 0, m}]) /. z -> x]], z];
         mat = CoefficientArrays[cl, cvars][[2]];
         MapThread[(cvals[[#1 + 1]] = #2) &, {freeidx, LinearSolve[mat[[All, freeidx + 1]], -mat[[All, m + 1]]]}]];
       cvals], {kk, kmax}]];
Q12LS[S1_, S2_, mu_] := Module[{c = GaudinSolveLS[S1, S2, mu, 2], ds = degQ[S1, S2, mu], qs},
   qs = Table[Sum[c[[kk, j + 1]] z^j, {j, 0, ds[[kk]]}], {kk, 2}]; Expand[qs[[1]] D[qs[[2]], z] - qs[[2]] D[qs[[1]], z]]];
det1LS[S1_, S2_, mu_] := Module[{q = Q12LS[S1, S2, mu]}, Expand[Det[{{doRes0[S1, S2, q], doRes0[S1, S2, z q]}, {doRes1bar[S1, S2, q], doRes1bar[S1, S2, z q]}}]/(q /. z -> 0)/(q /. z -> 1)]];
bench[s_] := Module[{ir = irreps[s, s]},
   Do[det1[s, s, mu]; det1LS[s, s, mu], {mu, ir}];
   <|"S" -> s, "correct" -> (Simplify[(det1[s, s, #] - det1LS[s, s, #]) & /@ ir] // Flatten // Union),
     "Solve" -> Min[Table[Total[Table[First@AbsoluteTiming[det1[s, s, mu]], {mu, ir}]], {3}]],
     "LinearSolve" -> Min[Table[Total[Table[First@AbsoluteTiming[det1LS[s, s, mu]], {mu, ir}]], {3}]]|>];
bench /@ {8, 10, 12}
```
**Output:**
```
Out= {<|"S" -> 8, "correct" -> {0}, "Solve" -> 0.196913, "LinearSolve" -> 0.171509|>, <|"S" -> 10, "correct" -> {0}, "Solve" -> 0.37141199999999996, "LinearSolve" -> 0.337232|>, <|"S" -> 12, "correct" -> {0}, "Solve" -> 0.686871, "LinearSolve" -> 0.6000409999999999|>}
```

## 2026-07-22 08:46:11
**Input:**
```wolfram
(* verify LinearSolve solver matches Solve for full kmax=4 (higher gauge constraints) across sweeps *)
chk4[s_] := (GaudinSolve[s, s, #, 4] - GaudinSolveLS[s, s, #, 4]) & /@ irreps[s, s] // Flatten // Union;
{"match4_S5" -> chk4[5], "match4_S8" -> chk4[8], "match4_S10" -> chk4[10],
 "worstODE_LS" -> Max[Abs[Flatten[Table[
      Module[{ds = degQ[S1, S2, mu], c = GaudinSolveLS[S1, S2, mu, 4], qs, ode},
        qs = Table[Sum[c[[kk, j + 1]] z^j, {j, 0, ds[[kk]]}], {kk, 4}];
        ode = gaudinODE /. {C2 -> Cas2[mu], C4 -> Cas4[mu], S -> S1, P -> S2} /. wsub;
        Table[Expand[ode /. Q -> Function[x, qs[[k]] /. z -> x]], {k, 4}]],
      {S1, 1, 2}, {S2, 1, 2}, {mu, irreps[S1, S2]}]] /. z -> 37]]}
```
**Output:**
```
Out= {"match4_S5" -> {0}, "match4_S8" -> {0}, "match4_S10" -> {0}, "worstODE_LS" -> 0}
```
## 2026-07-22 08:46:42 — ✏️ EDIT CELL 6
(* --- Gaudin ODE solver ---
   Arguments always ordered (S1,S2,mu): S1,S2 label the two input reps [S1,S1,0,0],[S2,S2,0,0]; mu
   is the output irrep. Solves the 4th-order Gaudin ODE for the four pol…


## 2026-07-22 08:46:42
**Input:**
```wolfram
(* --- Gaudin ODE solver ---
   Arguments always ordered (S1,S2,mu): S1,S2 label the two input reps [S1,S1,0,0],[S2,S2,0,0]; mu
   is the output irrep. Solves the 4th-order Gaudin ODE for the four polynomial Q's of degrees
   degQ[S1,S2,mu], returned as monic coefficient lists (coeffs of z^0,z^1,...,z^m). Optional kmax
   (default 4) returns only the first kmax Q's -- Q12/det1 need just Q1,Q2, and the higher Q's carry
   the largest degrees, so kmax=2 is much cheaper. GAUGE FIXING: the four Q-degrees are strictly
   increasing, so any lower-degree solution can be added to a higher one (q_k -> q_k + c q_j, j<k);
   removed by forcing the coefficient of z^(deg q_j) in q_k to vanish for every j<k.
   The unknown coefficients enter the ODE linearly, so we build the (numeric, rational) coefficient
   matrix and LinearSolve it -- faster than a symbolic Solve. cq[m]=1 (monic) and the gauge coeffs=0
   are imposed by dropping those columns (the monic one is moved to the right-hand side). *)
ClearAll[GaudinSolve];
GaudinSolve[S1_, S2_, mu_] := GaudinSolve[S1, S2, mu, 4];
GaudinSolve[S1_, S2_, mu_, kmax_] := Module[{ds, ode},
   ds = degQ[S1, S2, mu];
   ode = gaudinODE /. {C2 -> Cas2[mu], C4 -> Cas4[mu], S -> S1, P -> S2} /. wsub;
   Table[Module[{m = ds[[kk]], gauge = ds[[1 ;; kk - 1]], cvars, cl, mat, freeidx, cvals},
       cvars = Table[cq[j], {j, 0, m}];
       cvals = ConstantArray[0, m + 1]; cvals[[m + 1]] = 1;   (* cq[m]=1 monic; gauge coeffs stay 0 *)
       freeidx = Complement[Range[0, m], Append[gauge, m]];   (* 0-based indices of the free coeffs *)
       If[freeidx =!= {},
         cl = CoefficientList[Expand[ode /. Q -> Function[x, (cvars . Table[z^j, {j, 0, m}]) /. z -> x]], z];
         mat = CoefficientArrays[cl, cvars][[2]];              (* cl == mat . cvars, homogeneous *)
         MapThread[(cvals[[#1 + 1]] = #2) &, {freeidx, LinearSolve[mat[[All, freeidx + 1]], -mat[[All, m + 1]]]}]];
       cvals], {kk, kmax}]];
```
**Output:**
```
Out= Null
```

## 2026-07-22 08:48:02
**Input:**
```wolfram
known = Simplify[{det1[1, 1, {2, 1, 1, 0}], det1[1, 2, {3, 2, 1, 0}], det1[2, 1, {2, 2, 1, 1}]}];
worst = Max[Abs[Flatten[Table[Module[{qs = BuildQs[S1, S2, mu], ode = gaudinODE /. {C2 -> Cas2[mu], C4 -> Cas4[mu], S -> S1, P -> S2} /. wsub}, Table[Expand[ode /. Q -> Function[x, qs[[k]] /. z -> x]], {k, 4}]], {S1, 1, 2}, {S2, 1, 2}, {mu, irreps[S1, S2]}]] /. z -> 37]];
Do[det1[s, s, mu], {s, {8, 10, 12}}, {mu, irreps[s, s]}];
t = Association @@ Table[s -> Min[Table[Total[Table[First@AbsoluteTiming[det1[s, s, mu]], {mu, irreps[s, s]}]], {3}]], {s, {8, 10, 12}}];
{"known" -> known, "worstODE" -> worst, "sweep_times" -> t}
```
**Output:**
```
Out= {"known" -> {w^2, (w^3*(3 + 2*w))/12, (w^2*(-40 - 4*w + w^2))/6}, "worstODE" -> 0, "sweep_times" -> <|8 -> 0.180288, 10 -> 0.343754, 12 -> 0.601361|>}
```

---

## 2026-07-22 08:49:48 — 🔄 KERNEL RESTART


## 2026-07-22 08:57:38
**Input:**
```wolfram
test[S1_, S2_] := Module[{ir = irreps[S1, S2], dexp = 2 (S1 + S2) - 2, fac = (2 S1 - 1)! (2 S2 - 1)!},
   {"S1,S2" -> {S1, S2}, "predDeg" -> dexp, "factor" -> fac,
    "degOK" -> (Union[(Exponent[det1[S1, S2, #], w]) & /@ ir] === {dexp}),
    "leadCoeffs" -> Union[Coefficient[det1[S1, S2, #], w, dexp] & /@ ir],
    "monicAfterScale" -> Union[Coefficient[fac det1[S1, S2, #], w, dexp] & /@ ir]}];
{test[1, 1], test[1, 2], test[2, 1], test[2, 2], test[3, 2], test[2, 4], test[3, 3]}
```
**Output:**
```
Out= {{"S1,S2" -> {1, 1}, "predDeg" -> 2, "factor" -> 1, "degOK" -> True, "leadCoeffs" -> {1}, "monicAfterScale" -> {1}}, {"S1,S2" -> {1, 2}, "predDeg" -> 4, "factor" -> 6, "degOK" -> True, "leadCoeffs" -> {1/6}, "monicAfterScale" -> {1}}, {"S1,S2" -> {2, 1}, "predDeg" -> 4, "factor" -> 6, "degOK" -> True, "leadCoeffs" -> {1/6}, "monicAfterScale" -> {1}}, {"S1,S2" -> {2, 2}, "predDeg" -> 6, "factor" -> 36, "degOK" -> True, "leadCoeffs" -> {1/36}, "monicAfterScale" -> {1}}, {"S1,S2" -> {3, 2}, "predDeg" -> 8, "factor" -> 720, "degOK" -> True, "leadCoeffs" -> {1/720}, "monicAfterScale" -> {1}}, {"S1,S2" -> {2, 4}, "predDeg" -> 10, "factor" -> 30240, "degOK" -> True, "leadCoeffs" -> {1/30240}, "monicAfterScale" -> {1}}, {"S1,S2" -> {3, 3}, "predDeg" -> 10, "factor" -> 14400, "degOK" -> True, "leadCoeffs" -> {1/14400}, "monicAfterScale" -> {1}}}
```

## 2026-07-22 08:57:53
**Input:**
```wolfram
fac = (2*5 - 1)! (2*5 - 1)!;
raw = det1[5, 5, {10, 8, 2, 0}];
monic = Expand[fac raw];
{"factor_9!^2" -> fac,
 "raw_coeffs" -> CoefficientList[raw, w],
 "monic_coeffs" -> CoefficientList[monic, w],
 "monic_is_monic" -> (Coefficient[monic, w, Exponent[monic, w]] == 1),
 "monic_denominators" -> Union[Denominator /@ CoefficientList[monic, w]],
 "userFactor_2540160000_gives" -> CoefficientList[Expand[2540160000 raw], w]}
```
**Output:**
```
Out= {"factor_9!^2" -> 131681894400, "raw_coeffs" -> {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1/2540160000, 0, 1/131681894400}, "monic_coeffs" -> {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1296/25, 0, 1}, "monic_is_monic" -> True, "monic_denominators" -> {1, 25}, "userFactor_2540160000_gives" -> {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 25/1296}}
```
## 2026-07-22 08:59:02 — ✏️ EDIT CELL 14
(* det1: the br[1] entries carry a common factor Exp[w] (from res1) that factors out of the 2x2
   determinant and cancels the overall 1/Exp[w], leaving a plain polynomial in w. Stripping that
   Exp[…


## 2026-07-22 08:59:02
**Input:**
```wolfram
(* det1: the br[1] entries carry a common factor Exp[w] (from res1) that factors out of the 2x2
   determinant and cancels the overall 1/Exp[w], leaving a plain polynomial in w. Stripping that
   Exp[w] up front -- res1bar = res1/Exp[w] -- lets us Expand a bare polynomial instead of
   ExpandAll-ing against the exponential: identical result, ~15x faster at S1=S2=10. res1bar is
   memoized like res0/res1.
   MONIC NORMALIZATION: the bare determinant has leading-in-w coefficient 1/((2S1-1)!(2S2-1)!) for
   every mu (each residue's top-w term is a pure factorial: res0[0] ~ 1/(2S1-1)!, res1bar[n] ~
   1/(2S2-1)!, and the Q12(0)Q12(1) cancels), so we multiply by (2S1-1)!(2S2-1)! to make the
   generating function MONIC in w. This also clears the large factorial denominators. *)
ClearAll[det1, res1bar, doRes1bar];
res1bar[n_, S1_, S2_] := res1bar[n, S1, S2] = LaguerreL[2 S2 - 1, n - 2 S1 - 2 S2 + 1, -w];
doRes1bar[S1_, S2_, f_] := Module[{c = CoefficientList[f, z]}, c . Table[res1bar[n, S1, S2], {n, 0, Length[c] - 1}]];
det1[S1_, S2_, mu_] := Module[{q = Q12[S1, S2, mu]},
   Expand[(2 S1 - 1)! (2 S2 - 1)! Det[{{doRes0[S1, S2, q],    doRes0[S1, S2, z q]},
               {doRes1bar[S1, S2, q], doRes1bar[S1, S2, z q]}}]/(q /. z -> 0)/(q /. z -> 1)]
   ];
```
**Output:**
```
Out= Null
```
## 2026-07-22 08:59:12 — ✏️ EDIT CELL 13
## `det1` -- 2x2 residue determinant (monic in `w`)

$$\det{}_1[S_1,S_2,\mu]=\frac{(2S_1-1)!\,(2S_2-1)!}{Q_{12}(0)\,Q_{12}(1)\,e^{w}}\;\det\begin{pmatrix} \mathrm{br}_0[e^{wz}Q_{12}(z)] & \mathrm{br}_…


## 2026-07-22 09:00:26
**Input:**
```wolfram
leadOK[S1_, S2_] := Union[Coefficient[det1[S1, S2, #], w, Exponent[det1[S1, S2, #], w]] & /@ irreps[S1, S2]];
{"monic_leadcoeffs_perS" -> Association[{S1, S2} -> leadOK[S1, S2] & @@@ {{1, 1}, {2, 3}, {3, 3}, {5, 5}, {4, 6}}],
 "examples" -> {det1[1, 1, {2, 1, 1, 0}], det1[1, 2, {3, 2, 1, 0}], det1[2, 1, {2, 2, 1, 1}]},
 "S5_example_coeffs" -> CoefficientList[det1[5, 5, {10, 8, 2, 0}], w]}
```
**Output:**
```
[message] Table::iterb: Iterator {j, 0, S1 + S2} does not have appropriate bounds.
[message] [Table::iterb ×2 more suppressed — first instance shown above]
[message] ConstantArray::ilsmn: Single or list of non-negative machine-sized integers expected at position 2 of ConstantArray[0, 1 + S1 + S2].
[message] Set::pkspec1: The expression 1 + S1 + S2 cannot be used as a part specification.
[message] Range::range: Range specification in Range[0, S1 + S2] does not have appropriate bounds.
[message] Complement::heads: Heads List and Range at positions 2 and 1 are expected to be the same.
[message] General::stop: Further output of Table::iterb will be suppressed during this calculation.
[message] CoefficientArrays::poly: 16 S1 Table[0, {j, 0, S1 + S2}] . Table[j z      , {j, 0, S1 + S2}] + 16 S1  Table[0, {j, 0, S1 + S2}] . Table[j z      , {j, 0, S1 + S2}] + 24 S1 Table[0, {j, 0, S1 + S2}] . … [2095 chars]
[message] Part::pkspec1: The expression 1 + Complement[Range[0, S1 + S2], {S1 + S2}] cannot be used as a part specification.
[message] [Part::pkspec1 ×2 more suppressed — first instance shown above]
[message] MapThread::mptd: Object Complement[Range[0, S1 + S2], {S1 + S2}] at position {2, 1} in MapThread[(cvals$31804[[#1 + 1]] = #2) & , {Complement[Range[0, S1 + S2], {S1 + S2}], LinearSolve[Table[cq[j], {j… [357 chars]
[message] [ConstantArray::ilsmn ×2 more suppressed — first instance shown above]
[message] [Set::pkspec1 ×2 more suppressed — first instance shown above]
[message] [Range::range ×2 more suppressed — first instance shown above]
[message] [Complement::heads ×2 more suppressed — first instance shown above]
[message] [CoefficientArrays::poly ×2 more suppressed — first instance shown above]
[message] [General::stop ×9 more suppressed — first instance shown above]
[message] [MapThread::mptd ×2 more suppressed — first instance shown above]
[message] Power::indet: Indeterminate expression 0  encountered.
[message] Power::infy: Infinite expression - encountered.                                  0
[message] [Power::infy ×2 more suppressed — first instance shown above]
[message] [Power::indet ×2 more suppressed — first instance shown above]
Out= {"monic_leadcoeffs_perS" -> <|{S1, S2} -> {-(((-1 + 2*S1)!*(-1 + 2*S2)!*If[2*S1 > 1, (-1)^(2*S1 + 2*S2 - 1 + 1)*LaguerreL[2*S1 - 1 - 1, 1 - 2*S1 - 2*S2 + 1, w], 0]*LaguerreL[-1 + 2*S2, 1 - 2*S1 - 2*S2, -w]*Sum[z^j*{ConstantArray[0, 1 + S1 + S2], ConstantArray[0, 2 + S1 + S2]}[[1,1 + j]], {j, 0, S1 + S2}]^2*Sum[j*z^(-1 + j)*{ConstantArray[0, 1 + S1 + S2], ConstantArray[0, 2 + S1 + S2]}[[2,1 + j]], {j, 0, 1 + S1 + S2}]^2)/((-(Sum[j*{ConstantArray[0, 1 + S1 + S2], ConstantArray[0, 2 + S1 + S2]}[[1,1 + j]], {j, 0, S1 + S2}]*Sum[{ConstantArray[0, 1 + S1 + S2], ConstantArray[0, 2 + S1 + S2]}[[2,1 + j]], {j, 0, 1 + S1 + S2}]) + Sum[{ConstantArray[0, 1 + S1 + S2], ConstantArray[0, 2 + S1 + S2]}[[1,1 + j]], {j, 0, S1 + S2}]*Sum[j*{ConstantArray[0, 1 + S1 + S2], ConstantArray[0, 2 + S1 + S2]}[[2,1 + j]], {j, 0, 1 + S1 + S2}])*(-(Sum[0^(-1 + j)*j*{ConstantArray[0, 1 + S1 + S2], ConstantArray[0, 2 + S1 + S2]}[[1,1 + j]], {j, 0, S1 + S2}]*Sum[0^j*{ConstantArray[0, 1 + S1 + S2], ConstantArray[0, 2 + S1 + S2]}[[2,1 + j]], {j, 0, 1 + S1 + S2}]) + Sum[0^j*{ConstantArray[0, 1 + S1 + S2], ConstantArray[0, 2 + S1 + S2]}[[1,1 + j]], {j, 0, S1 + S2}]*Sum[0^(-1 + j)*j*{ConstantArray[0, 1 + S1 + S2], ConstantArray[0, 2 + S1 + S2]}[[2,1 + j]], {j, 0, 1 + S1 + S2}]))) + ((-1 + 2*S1)!*(-1 + 2*S2)!*If[2*S1 > 0, (-1)^(2*S1 + 2*S2 - 0 + 1)*LaguerreL[2*S1 - 0 - 1, 1 - 2*S1 - 2*S2 + 0, w], 0]*LaguerreL[-1 + 2*S2, 2 - 2*S1 - 2*S2, -w]*Sum[z^j*{ConstantArray[0, 1 + S1 + S2], ConstantArray[0, 2 + S1 + S2]}[[1,1 + j]], {j, 0, S1 + S2}]^2*Sum[j*z^(-1 + j)*{ConstantArray[0, 1 + S1 + S2], ConstantArray[0, 2 + S1 + S2]}[[2,1 + j]], {j, 0, 1 + S1 + S2}]^2)/((-(Sum[j*{ConstantArray[0, 1 + S1 + S2], ConstantArray[0, 2 + S1 + S2]}[[1,1 + j]], {j, 0, S1 + S2}]*Sum[{ConstantArray[0, 1 + S1 + S2], ConstantArray[0, 2 + S1 + S2]}[[2,1 + j]], {j, 0, 1 + S1 + S2}]) + Sum[{ConstantArray[0, 1 + S1 + S2], ConstantArray[0, 2 + S1 + S2]}[[1,1 + j]], {j, 0, S1 + S2}]*Sum[j*{ConstantArray[0, 1 + S1 + S2], ConstantArray[0, 2 + S1 + S2]}[[2,1 + j]], {j, 0, 1 + S1 + S2}])*(-(Sum[0^(-1 + j)*j*{ConstantArray[0, 1 + S1 + S2], ConstantArray[0, 2 + S1 + S2]}[[1,1 + j]], {j, 0, S1 + S2}]*Sum[0^j*{ConstantArray[0, 1 + S1 + S2], ConstantArray[0, 2 + S1 + S2]}[[2,1 + j]], {j, 0, 1 + S1 + S2}]) + Sum[0^j*{ConstantArray[0, 1 + S1 + S2], ConstantArray[0, 2 + S1 + S2]}[[1,1 + j]], {j, 0, S1 + S2}]*Sum[0^(-1 + j)*j*{ConstantArray[0, 1 + S1 + S2], ConstantArray[0, 2 + S1 + S2]}[[2,1 + j]], {j, 0, 1 + S1 + S2}])) + (2*(-1 + 2*S1)!*(-1 + 2*S2)!*If[2*S1 > 1, (-1)^(2*S1 + 2*S2 - 1 + 1)*LaguerreL[2*S1 - 1 - 1, 1 - 2*S1 - 2*S2 + 1, w], 0]*LaguerreL[-1 + 2*S2, 1 - 2*S1 - 2*S2, -w]*Sum[j*z^(-1 + j)*{ConstantArray[0, 1 + S1 + S2], ConstantArray[0, 2 + S1 + S2]}[[1,1 + j]], {j, 0, S1 + S2}]*Sum[z^j*{ConstantArray[0, 1 + S1 + S2], ConstantArray[0, 2 + S1 + S2]}[[1,1 + j]], {j, 0, S1 + S2}]*Sum[j*z^(-1 + j)*{ConstantArray[0, 1 + S1 + S2], ConstantArray[0, 2 + S1 + S2]}[[2,1 + j]], {j, 0, 1 + S1 + S2}]*Sum[z^j*{ConstantArray[0, 1 + S1 + S2], ConstantArray[0, 2 + S1 + S2]}[[2,1 + j]], {j, 0, 1 + S1 + S2}])/((-(Sum[j*{ConstantArray[0, 1 + S1 + S2], ConstantArray[0, 2 + S1 + S2]}[[1,1 + j]], {j, 0, S1 + S2}]*Sum[{ConstantArray[0, 1 + S1 + S2], ConstantArray[0, 2 + S1 + S2]}[[2,1 + j]], {j, 0, 1 + S1 + S2}]) + Sum[{ConstantArray[0, 1 + S1 + S2], ConstantArray[0, 2 + S1 + S2]}[[1,1 + j]], {j, 0, S1 + S2}]*Sum[j*{ConstantArray[0, 1 + S1 + S2], ConstantArray[0, 2 + S1 + S2]}[[2,1 + j]], {j, 0, 1 + S1 + S2}])*(-(Sum[0^(-1 + j)*j*{ConstantArray[0, 1 + S1 + S2], ConstantArray[0, 2 + S1 + S2]}[[1,1 + j]], {j, 0, S1 + S2}]*Sum[0^j*{ConstantArray[0, 1 + S1 + S2], ConstantArray[0, 2 + S1 + S2]}[[2,1 + j]], {j, 0, 1 + S1 + S2}]) + Sum[0^j*{ConstantArray[0, 1 + S1 + S2], ConstantArray[0, 2 + S1 + S2]}[[1,1 + j]], {j, 0, S1 + S2}]*Sum[0^(-1 + j)*j*{ConstantArray[0, 1 + S1 + S2], ConstantArray[0, 2 + S1 + S2]}[[2,1 + j]], {j, 0, 1 + S1 + S2}])) - (2*(-1 + 2*S1)!*(-1 + 2*S2)!*If[2*S1 > 0, (-1)^(2*S1 + 2*S2 - 0 + 1)*LaguerreL[2*S1 - 0 - 1, 1 - 2*S1 - 2*S2 + 0, w], 0]*LaguerreL[-1 + 2*S2, 2 - 2*S1 - 2*S2, -w]*Sum[j*z^(-1 + j)*{C
[output truncated — 7932 chars total]
```

## 2026-07-22 09:00:38
**Input:**
```wolfram
leadOK[S1_, S2_] := Union[Function[mu, Coefficient[det1[S1, S2, mu], w, Exponent[det1[S1, S2, mu], w]]] /@ irreps[S1, S2]];
{"leadcoeffs_1,1" -> leadOK[1, 1], "leadcoeffs_2,3" -> leadOK[2, 3], "leadcoeffs_3,3" -> leadOK[3, 3], "leadcoeffs_5,5" -> leadOK[5, 5],
 "examples" -> {det1[1, 1, {2, 1, 1, 0}], det1[1, 2, {3, 2, 1, 0}], det1[2, 1, {2, 2, 1, 1}]},
 "S5_coeffs" -> CoefficientList[det1[5, 5, {10, 8, 2, 0}], w]}
```
**Output:**
```
Out= {"leadcoeffs_1,1" -> {1}, "leadcoeffs_2,3" -> {1}, "leadcoeffs_3,3" -> {1}, "leadcoeffs_5,5" -> {1}, "examples" -> {w^2, (3*w^3)/2 + w^4, -40*w^2 - 4*w^3 + w^4}, "S5_coeffs" -> {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1296/25, 0, 1}}
```

---

## 2026-07-22 09:01:52 — 🔄 KERNEL RESTART


---

## 2026-07-22 09:06:00 — 🔄 KERNEL RESTART

