# AI Action Log — GaudinSolver

## 2026-07-21 06:13:34 — 📥 BULK INSERT 5 CELL(S) at positions 1–5
1. [markdown] # Gaudin ODE solver (gl(4), two-site, [S,S,0,0] x [P,P,0,0])↵↵Self-contained solver for the gl(4) Ga
2. [code] (* --- Tensor-product decomposition [S,S,0,0] x [P,P,0,0] for gl(4) via Schur polynomials ---↵   The
3. [code] (* --- Casimir eigenvalues on the output irrep mu, and the Q-polynomial degrees ---↵   Column-ordere
4. [code] (* --- Gaudin ODE (hardcoded closed form; derived as the h->0 limit of the XXX Baxter equation↵   in
5. [code] (* --- Gaudin ODE solver ---↵   Input: irrep mu in [S,S,0,0]x[P,P,0,0]. Solves the 4th-order Gaudin


## 2026-07-21 06:13:34
**Input:**
```wolfram
(* --- Tensor-product decomposition [S,S,0,0] x [P,P,0,0] for gl(4) via Schur polynomials ---
   The gl(4) character of irrep lam is the Schur polynomial s_lam(x1..x4) (Weyl bialternant).
   Since s_lam s_mu = Sum_nu c^nu_{lam mu} s_nu, the irreps nu in the tensor product are exactly
   those in the product's Schur expansion; c^nu = Littlewood-Richardson multiplicity. Decompose
   by peeling off the lex-highest monomial (a dominant highest weight) with its coefficient. *)
ClearAll[schur, tensorDecomp, dimGL4, irreps];
schur[lam_, x_] := Cancel[Det[Table[x[[i]]^(lam[[j]] + Length[x] - j), {i, Length[x]}, {j, Length[x]}]]/Det[Table[x[[i]]^(Length[x] - j), {i, Length[x]}, {j, Length[x]}]]];
tensorDecomp[lam_, mu_] := Module[{x = Array[xx, 4], prod, rules, nu, mult, out = {}},
   prod = Expand[schur[lam, x] schur[mu, x]];
   While[prod =!= 0, rules = Association[CoefficientRules[prod, x]]; nu = Last[Sort[Keys[rules]]];
     mult = rules[nu]; AppendTo[out, nu -> mult]; prod = Expand[prod - mult schur[nu, x]]];
   out];
dimGL4[lam_] := Product[(lam[[i]] - lam[[j]] + j - i)/(j - i), {i, 4}, {j, i + 1, 4}];
(* the irreps mu (all multiplicity 1 for this family) appearing in [S,S,0,0] x [P,P,0,0] *)
irreps[S_, P_] := tensorDecomp[{S, S, 0, 0}, {P, P, 0, 0}][[All, 1]];
tensorDecomp[{1, 1, 0, 0}, {1, 1, 0, 0}]
```
**Output:**
```
(no output)
```

## 2026-07-21 06:13:35
**Input:**
```wolfram
(* --- Casimir eigenvalues on the output irrep mu, and the Q-polynomial degrees ---
   Column-ordered generating function (Kazakov et al conventions): G4[u] = Prod_i (u - h(mu_i+4-i)),
   whose coefficients are the Casimirs. Eigenvalues are symmetric functions of the centered
   shifted weights xw_i = mu_i + 5/2 - i (= mu_i + rho_i); Pw_k = Sum_i xw_i^k. *)
ClearAll[muhat, xw, Pw, Cas2, Cas4, degQ];
muhat[mu_] := Table[mu[[i]] + 4 - i, {i, 4}];
xw[mu_] := Table[mu[[i]] + 5/2 - i, {i, 4}];
Pw[mu_, k_] := Total[xw[mu]^k];
Cas2[mu_] := Pw[mu, 2] - 5;
Cas4[mu_] := Pw[mu, 4] + 4 Pw[mu, 3] - Pw[mu, 1] Pw[mu, 2] + Pw[mu, 2]/2 - 2 Pw[mu, 1]^2 - 10 Pw[mu, 1] - 51/4;
(* the SAME shifted weights fix the Q-degrees: deg q_k = (S+P+3) - muhat_k = (S+P) - mu_k + (k-1) *)
degQ[mu_, SS_, PP_] := Table[(SS + PP + 3) - muhat[mu][[k]], {k, 4}];
{Cas2[{2, 1, 1, 0}], Cas4[{1, 0, 0, 0}], degQ[{m1, m2, m3, m4}, S, P]}
```
**Output:**
```
(no output)
```

## 2026-07-21 06:13:35
**Input:**
```wolfram
(* --- Gaudin ODE (hardcoded closed form; derived as the h->0 limit of the XXX Baxter equation
   in baxter_SU4_L2_XXX.wb, cells 17-22 -- NOT re-derived here). Q'''' coefficient is
   8(u-th1)^2(u-th2)^2 -- double singular points at the two Gaudin sites, as expected. C2,C4 are
   the abstract Casimirs of the output irrep mu; S,P label the two input reps [S,S,0,0],[P,P,0,0].
   Verified (in the source notebook) to have polynomial solutions of exactly degree degQ[mu,S,P]
   for every irrep of every (S,P) tested. *)
ClearAll[gaudinODE];
gaudinODE = (2*C2 + C2^2 - 2*C4 + 24*P + 2*C2*P + 20*P^2 + 4*P^3 + 24*S + 2*C2*S + 40*P*S + 12*P^2*S + 20*S^2 + 12*P*S^2 + 4*S^3)*Q[u]
   + 4*(C2 - 4*P - 2*P^2 - 4*S - 4*P*S - 2*S^2)*(2*u + P*u + S*u - θ1 - P*θ1 - θ2 - S*θ2)*Derivative[1][Q][u]
   - 4*(C2*u^2 - 6*P*u^2 - 4*P^2*u^2 - 6*S*u^2 - 8*P*S*u^2 - 4*S^2*u^2 - C2*u*θ1 + 8*P*u*θ1 + 6*P^2*u*θ1 + 4*S*u*θ1 + 8*P*S*u*θ1
      + 2*S^2*u*θ1 - 2*P*θ1^2 - 2*P^2*θ1^2 - C2*u*θ2 + 4*P*u*θ2 + 2*P^2*u*θ2 + 8*S*u*θ2 + 8*P*S*u*θ2 + 6*S^2*u*θ2
      + C2*θ1*θ2 - 4*P*θ1*θ2 - 2*P^2*θ1*θ2 - 4*S*θ1*θ2 - 8*P*S*θ1*θ2 - 2*S^2*θ1*θ2 - 2*S*θ2^2 - 2*S^2*θ2^2)*Derivative[2][Q][u]
   - 16*(u - θ1)*(u - θ2)*(P*u + S*u - P*θ1 - S*θ2)*Derivative[3][Q][u]
   + 8*(u - θ1)^2*(u - θ2)^2*Derivative[4][Q][u];
(* sites fixed at theta1=0, theta2=1 for this solver *)
thsub = {θ1 -> 0, θ2 -> 1};
```
**Output:**
```
(no output)
```

## 2026-07-21 06:13:35
**Input:**
```wolfram
(* --- Gaudin ODE solver ---
   Input: irrep mu in [S,S,0,0]x[P,P,0,0]. Solves the 4th-order Gaudin ODE for the four
   polynomial Q's of degrees degQ[mu,S,P]. GAUGE FIXING: the four Q-degrees are strictly
   increasing, so any lower-degree solution can be added to a higher one (q_k -> q_k + c q_j,
   j<k); removed by forcing the coefficient of u^(deg q_j) in q_k to vanish for every j<k. *)
ClearAll[GaudinSolve];
GaudinSolve[mu_, SS_, PP_] := Module[{ds, ode},
   ds = degQ[mu, SS, PP];
   ode = gaudinODE /. {C2 -> Cas2[mu], C4 -> Cas4[mu], S -> SS, P -> PP} /. thsub;
   Table[Module[{m = ds[[kk]], gauge = ds[[1 ;; kk - 1]], q, cl, sol},
       q = Sum[cq[j] u^j, {j, 0, m}];
       cl = CoefficientList[Expand[ode /. Q -> Function[x, q /. u -> x]], u];
       sol = Solve[Join[Thread[cl == 0], {cq[m] == 1}, Thread[(cq /@ gauge) == 0]], Table[cq[j], {j, 0, m}]];
       Table[cq[j], {j, 0, m}] /. sol[[1]]], {kk, 4}]];
(* example: mu={2,1,1,0} at (S,P)=(1,1) -> q1=1, q2=u(u-1), q3=u^3, q4=u^5-(5/2)u^4 *)
GaudinSolve[{2, 1, 1, 0}, 1, 1]
```
**Output:**
```
(no output)
```
## 2026-07-21 06:13:49 — ▶️ RUN CELLS 2–5
Cell 2: ✓ — (no output)
Cell 3: ✓ — (no output)
Cell 4: ✓ — (no output)
Cell 5: ✓ — (no output)


## 2026-07-21 06:13:53
**Input:**
```wolfram
GaudinSolve[{2, 1, 1, 0}, 1, 1]
```
**Output:**
```
Out= {{1}, {0, -1, 1}, {0, 0, 0, 1}, {0, 0, 0, 0, -5/2, 1}}
```

## 2026-07-21 06:13:59
**Input:**
```wolfram
{tensorDecomp[{1, 1, 0, 0}, {1, 1, 0, 0}], {Cas2[{2, 1, 1, 0}], Cas4[{1, 0, 0, 0}]}, degQ[{m1, m2, m3, m4}, S, P]}
```
**Output:**
```
Out= {{{2, 2, 0, 0} -> 1, {2, 1, 1, 0} -> 1, {1, 1, 1, 1} -> 1}, {12, 64}, {-m1 + P + S, 1 - m2 + P + S, 2 - m3 + P + S, 3 - m4 + P + S}}
```

## 2026-07-21 06:14:05
**Input:**
```wolfram
verify[mu_, SS_, PP_] := Module[{qs, ds},
   ds = degQ[mu, SS, PP]; qs = GaudinSolve[mu, SS, PP];
   AllTrue[Range[4], (Expand[(gaudinODE /. {C2 -> Cas2[mu], C4 -> Cas4[mu], S -> SS, P -> PP} /. thsub)
      /. Q -> Function[x, (qs[[#]] . x^Range[0, ds[[#]]])]] === 0) &]];
AllTrue[irreps[3, 3], verify[#, 3, 3] &]
```
**Output:**
```
Out= False
```

## 2026-07-21 06:14:09
**Input:**
```wolfram
Table[{mu, verify[mu, 3, 3]}, {mu, irreps[3, 3]}]
```
**Output:**
```
Out= {{{6, 6, 0, 0}, False}, {{6, 5, 1, 0}, False}, {{6, 4, 2, 0}, False}, {{6, 3, 3, 0}, False}, {{5, 5, 1, 1}, False}, {{5, 4, 2, 1}, False}, {{5, 3, 3, 1}, False}, {{4, 4, 2, 2}, False}, {{4, 3, 3, 2}, False}, {{3, 3, 3, 3}, False}}
```

## 2026-07-21 06:14:15
**Input:**
```wolfram
verify[{2, 1, 1, 0}, 1, 1]
```
**Output:**
```
Out= False
```

## 2026-07-21 06:14:19
**Input:**
```wolfram
mu = {2, 1, 1, 0}; SS = 1; PP = 1;
ds = degQ[mu, SS, PP];
qs = GaudinSolve[mu, SS, PP];
odeS = gaudinODE /. {C2 -> Cas2[mu], C4 -> Cas4[mu], S -> SS, P -> PP} /. thsub;
qf1 = Function[x, qs[[1]] . x^Range[0, ds[[1]]]];
Expand[odeS /. Q -> qf1]
```
**Output:**
```
Out= {{0}}
```

## 2026-07-21 06:14:41
**Input:**
```wolfram
{"ds" -> ds, "qs" -> qs, "qf1[u]" -> qf1[u], "qf1''[u]" -> qf1''[u]}
```
**Output:**
```
Out= {"ds" -> {0, 2, 3, 5}, "qs" -> {{1}, {0, -1, 1}, {0, 0, 0, 1}, {0, 0, 0, 0, -5/2, 1}}, "qf1[u]" -> 1, "qf1''[u]" -> {0}}
```

## 2026-07-21 06:15:34
**Input:**
```wolfram
Clear[verify2];
verify2[mu_, SS_, PP_] := Module[{qs, ds, ode, qf},
   ds = degQ[mu, SS, PP]; qs = GaudinSolve[mu, SS, PP];
   ode = gaudinODE /. {C2 -> Cas2[mu], C4 -> Cas4[mu], S -> SS, P -> PP} /. thsub;
   AllTrue[Range[4], Module[{m = ds[[#]], cf = qs[[#]]},
      qf[x_] := Sum[cf[[j + 1]] x^j, {j, 0, m}];
      Expand[ode /. Q -> qf] === 0] &]];
{"{2,1,1,0}@(1,1)" -> verify2[{2, 1, 1, 0}, 1, 1],
 "(3,3) all irreps" -> AllTrue[irreps[3, 3], verify2[#, 3, 3] &]}
```
**Output:**
```
Out= {"{2,1,1,0}@(1,1)" -> True, "(3,3) all irreps" -> True}
```

## 2026-07-21 06:15:39
**Input:**
```wolfram
Grid[Table[{{SS, PP}, Length[irreps[SS, PP]], AllTrue[irreps[SS, PP], verify2[#, SS, PP] &]},
   {SS, {1, 2, 3}}, {PP, {1, 2, 3}}] // Flatten[#, 1] &, Frame -> All]
```
**Output:**
```
Out= Grid[{{{1, 1}, 3, True}, {{1, 2}, 3, True}, {{1, 3}, 3, True}, {{2, 1}, 3, True}, {{2, 2}, 6, True}, {{2, 3}, 6, True}, {{3, 1}, 3, True}, {{3, 2}, 6, True}, {{3, 3}, 10, True}}, Frame -> All]
```
## 2026-07-21 06:15:47 — 📥 BULK INSERT 1 CELL(S) at positions 6–6
1. [code] (* --- validation: every gauge-fixed q_k solves the Gaudin ODE exactly, across full (S,P)↵   decompo


## 2026-07-21 06:15:47
**Input:**
```wolfram
(* --- validation: every gauge-fixed q_k solves the Gaudin ODE exactly, across full (S,P)
   decompositions. NB: rebuild q as a Sum (matching how GaudinSolve itself builds it internally)
   -- reconstructing q via Dot[coeffs, u^Range[...]] and then differentiating symbolically has a
   degenerate edge case for length-1 coefficient lists (the deg=0 case) that spuriously fails. *)
ClearAll[verifyGaudin];
verifyGaudin[mu_, SS_, PP_] := Module[{qs, ds, ode, qf},
   ds = degQ[mu, SS, PP]; qs = GaudinSolve[mu, SS, PP];
   ode = gaudinODE /. {C2 -> Cas2[mu], C4 -> Cas4[mu], S -> SS, P -> PP} /. thsub;
   AllTrue[Range[4], Module[{m = ds[[#]], cf = qs[[#]]},
      qf[x_] := Sum[cf[[j + 1]] x^j, {j, 0, m}];
      Expand[ode /. Q -> qf] === 0] &]];
Grid[Table[{{SS, PP}, Length[irreps[SS, PP]], AllTrue[irreps[SS, PP], verifyGaudin[#, SS, PP] &]},
   {SS, {1, 2, 3}}, {PP, {1, 2, 3}}] // Flatten[#, 1] &, Frame -> All]
```
**Output:**
```
(no output)
```
## 2026-07-21 06:25:17 — ✏️ EDIT CELL 6
(* --- Gaudin ODE (hardcoded closed form; derived as the h->0 limit of the XXX Baxter equation
   in baxter_SU4_L2_XXX.wb, cells 17-22 -- NOT re-derived here). Q'''' coefficient is
   8(u-w1)^2(u-w2)^…


## 2026-07-21 06:25:17
**Input:**
```wolfram
(* --- Gaudin ODE (hardcoded closed form; derived as the h->0 limit of the XXX Baxter equation
   in baxter_SU4_L2_XXX.wb, cells 17-22 -- NOT re-derived here). Q'''' coefficient is
   8(u-w1)^2(u-w2)^2 -- double singular points at the two Gaudin sites, as expected. C2,C4 are
   the abstract Casimirs of the output irrep mu; S,P label the two input reps [S,S,0,0],[P,P,0,0].
   Verified (in the source notebook) to have polynomial solutions of exactly degree degQ[mu,S,P]
   for every irrep of every (S,P) tested. *)
ClearAll[gaudinODE];
gaudinODE = (2*C2 + C2^2 - 2*C4 + 24*P + 2*C2*P + 20*P^2 + 4*P^3 + 24*S + 2*C2*S + 40*P*S + 12*P^2*S + 20*S^2 + 12*P*S^2 + 4*S^3)*Q[u]
   + 4*(C2 - 4*P - 2*P^2 - 4*S - 4*P*S - 2*S^2)*(2*u + P*u + S*u - w1 - P*w1 - w2 - S*w2)*Derivative[1][Q][u]
   - 4*(C2*u^2 - 6*P*u^2 - 4*P^2*u^2 - 6*S*u^2 - 8*P*S*u^2 - 4*S^2*u^2 - C2*u*w1 + 8*P*u*w1 + 6*P^2*u*w1 + 4*S*u*w1 + 8*P*S*u*w1
      + 2*S^2*u*w1 - 2*P*w1^2 - 2*P^2*w1^2 - C2*u*w2 + 4*P*u*w2 + 2*P^2*u*w2 + 8*S*u*w2 + 8*P*S*u*w2 + 6*S^2*u*w2
      + C2*w1*w2 - 4*P*w1*w2 - 2*P^2*w1*w2 - 4*S*w1*w2 - 8*P*S*w1*w2 - 2*S^2*w1*w2 - 2*S*w2^2 - 2*S^2*w2^2)*Derivative[2][Q][u]
   - 16*(u - w1)*(u - w2)*(P*u + S*u - P*w1 - S*w2)*Derivative[3][Q][u]
   + 8*(u - w1)^2*(u - w2)^2*Derivative[4][Q][u];
(* sites fixed at w1=0, w2=1 for this solver *)
thsub = {w1 -> 0, w2 -> 1};
```
**Output:**
```
Out= Null
```
## 2026-07-21 06:25:26 — ✏️ EDIT CELL 9
(* --- (S,P)=(8,8): Gaudin ODE solver over the full [8,8,0,0]x[8,8,0,0] decomposition (45 irreps),
   w1=0, w2=1, per-irrep timing. Each q_k verified to solve the ODE exactly. *)
irr88 = irreps[8, 8];…


## 2026-07-21 06:25:28
**Input:**
```wolfram
(* --- (S,P)=(8,8): Gaudin ODE solver over the full [8,8,0,0]x[8,8,0,0] decomposition (45 irreps),
   w1=0, w2=1, per-irrep timing. Each q_k verified to solve the ODE exactly. *)
irr88 = irreps[8, 8];
tblG88 = Table[
   Module[{mu = irr88[[i]], tim, qs, ds, ode, ok, qf},
     ds = degQ[mu, 8, 8];
     {tim, qs} = AbsoluteTiming[GaudinSolve[mu, 8, 8]];
     ode = gaudinODE /. {C2 -> Cas2[mu], C4 -> Cas4[mu], S -> 8, P -> 8} /. thsub;
     ok = AllTrue[Range[4], (qf = qs[[#]] . u^Range[0, ds[[#]]];
        Expand[ode /. Q -> Function[x, qf /. u -> x]] === 0) &];
     {Row[mu, ","], dimGL4[mu], Row[ds, ","], NumberForm[tim, {6, 4}], ok}],
   {i, Length[irr88]}];
Print["# irreps: ", Length[irr88], "   all solve: ", AllTrue[tblG88, Last],
   "   total time (s): ", Total[tblG88[[All, 4]] /. NumberForm[x_, _] :> x]];
Grid[Prepend[tblG88, {"mu", "dim", "deg q_k", "time (s)", "solves?"}],
  Frame -> All, Alignment -> Left,
  Background -> {None, {RGBColor[.88, .95, .88], {GrayLevel[1], RGBColor[.96, .99, .96]}}}]
```
**Output:**
```
Out= Grid[{{"mu", "dim", "deg q_k", "time (s)", "solves?"}, {Row[{16, 16, 0, 0}, ","], 8721, Row[{0, 1, 18, 19}, ","], NumberForm[0.003987, {6, 4}], True}, {Row[{16, 15, 1, 0}, ","], 27455, Row[{0, 2, 17, 19}, ","], NumberForm[0.003287, {6, 4}], True}, {Row[{16, 14, 2, 0}, ","], 47424, Row[{0, 3, 16, 19}, ","], NumberForm[0.003262, {6, 4}], True}, {Row[{16, 13, 3, 0}, ","], 62700, Row[{0, 4, 15, 19}, ","], NumberForm[0.003488, {6, 4}], True}, {Row[{16, 12, 4, 0}, ","], 69825, Row[{0, 5, 14, 19}, ","], NumberForm[0.003113, {6, 4}], True}, {Row[{16, 11, 5, 0}, ","], 67431, Row[{0, 6, 13, 19}, ","], NumberForm[0.00362, {6, 4}], True}, {Row[{16, 10, 6, 0}, ","], 55860, Row[{0, 7, 12, 19}, ","], NumberForm[0.003631, {6, 4}], True}, {Row[{16, 9, 7, 0}, ","], 36784, Row[{0, 8, 11, 19}, ","], NumberForm[0.003003, {6, 4}], True}, {Row[{16, 8, 8, 0}, ","], 12825, Row[{0, 9, 10, 19}, ","], NumberForm[0.002903, {6, 4}], True}, {Row[{15, 15, 1, 1}, ","], 5440, Row[{1, 2, 17, 18}, ","], NumberForm[0.002928, {6, 4}], True}, {Row[{15, 14, 2, 1}, ","], 16575, Row[{1, 3, 16, 18}, ","], NumberForm[0.002934, {6, 4}], True}, {Row[{15, 13, 3, 1}, ","], 27489, Row[{1, 4, 15, 18}, ","], NumberForm[0.002767, {6, 4}], True}, {Row[{15, 12, 4, 1}, ","], 34476, Row[{1, 5, 14, 18}, ","], NumberForm[0.003061, {6, 4}], True}, {Row[{15, 11, 5, 1}, ","], 35700, Row[{1, 6, 13, 18}, ","], NumberForm[0.003558, {6, 4}], True}, {Row[{15, 10, 6, 1}, ","], 30855, Row[{1, 7, 12, 18}, ","], NumberForm[0.003056, {6, 4}], True}, {Row[{15, 9, 7, 1}, ","], 20825, Row[{1, 8, 11, 18}, ","], NumberForm[0.00312, {6, 4}], True}, {Row[{15, 8, 8, 1}, ","], 7344, Row[{1, 9, 10, 18}, ","], NumberForm[0.003105, {6, 4}], True}, {Row[{14, 14, 2, 2}, ","], 3185, Row[{2, 3, 16, 17}, ","], NumberForm[0.003136, {6, 4}], True}, {Row[{14, 13, 3, 2}, ","], 9295, Row[{2, 4, 15, 17}, ","], NumberForm[0.002827, {6, 4}], True}, {Row[{14, 12, 4, 2}, ","], 14580, Row[{2, 5, 14, 17}, ","], NumberForm[0.003757, {6, 4}], True}, {Row[{14, 11, 5, 2}, ","], 16940, Row[{2, 6, 13, 17}, ","], NumberForm[0.003338, {6, 4}], True}, {Row[{14, 10, 6, 2}, ","], 15625, Row[{2, 7, 12, 17}, ","], NumberForm[0.003045, {6, 4}], True}, {Row[{14, 9, 7, 2}, ","], 10935, Row[{2, 8, 11, 17}, ","], NumberForm[0.003096, {6, 4}], True}, {Row[{14, 8, 8, 2}, ","], 3920, Row[{2, 9, 10, 17}, ","], NumberForm[0.003134, {6, 4}], True}, {Row[{13, 13, 3, 3}, ","], 1716, Row[{3, 4, 15, 16}, ","], NumberForm[0.00284, {6, 4}], True}, {Row[{13, 12, 4, 3}, ","], 4719, Row[{3, 5, 14, 16}, ","], NumberForm[0.002943, {6, 4}], True}, {Row[{13, 11, 5, 3}, ","], 6825, Row[{3, 6, 13, 16}, ","], NumberForm[0.002926, {6, 4}], True}, {Row[{13, 10, 6, 3}, ","], 7020, Row[{3, 7, 12, 16}, ","], NumberForm[0.002995, {6, 4}], True}, {Row[{13, 9, 7, 3}, ","], 5200, Row[{3, 8, 11, 16}, ","], NumberForm[0.003668, {6, 4}], True}, {Row[{13, 8, 8, 3}, ","], 1911, Row[{3, 9, 10, 16}, ","], NumberForm[0.003099, {6, 4}], True}, {Row[{12, 12, 4, 4}, ","], 825, Row[{4, 5, 14, 15}, ","], NumberForm[0.003631, {6, 4}], True}, {Row[{12, 11, 5, 4}, ","], 2079, Row[{4, 6, 13, 15}, ","], NumberForm[0.003336, {6, 4}], True}, {Row[{12, 10, 6, 4}, ","], 2640, Row[{4, 7, 12, 15}, ","], NumberForm[0.003016, {6, 4}], True}, {Row[{12, 9, 7, 4}, ","], 2156, Row[{4, 8, 11, 15}, ","], NumberForm[0.002908, {6, 4}], True}, {Row[{12, 8, 8, 4}, ","], 825, Row[{4, 9, 10, 15}, ","], NumberForm[0.002951, {6, 4}], True}, {Row[{11, 11, 5, 5}, ","], 336, Row[{5, 6, 13, 14}, ","], NumberForm[0.002937, {6, 4}], True}, {Row[{11, 10, 6, 5}, ","], 735, Row[{5, 7, 12, 14}, ","], NumberForm[0.003023, {6, 4}], True}, {Row[{11, 9, 7, 5}, ","], 729, Row[{5, 8, 11, 14}, ","], NumberForm[0.002846, {6, 4}], True}, {Row[{11, 8, 8, 5}, ","], 300, Row[{5, 9, 10, 14}, ","], NumberForm[0.003071, {6, 4}], True}, {Row[{10, 10, 6, 6}, ","], 105, Row[{6, 7, 12, 13}, ","], NumberForm[0.002928, {6, 4}], True}, {Row[{10, 9, 7, 6}, ","], 175, Row[{6, 8, 11, 13}, ","], NumberForm[0.002966, {6, 4}], True}, {Row[{10, 8, 8, 6}, ","], 84, Row[{6, 9, 10, 13}, ","], NumberForm[0.00306, {6, 4}], True}, {Row[{9, 9, 7, 7}, 
[output truncated — 4481 chars total]
```
## 2026-07-21 06:25:41 — ▶️ RUN CELLS 3–9
Cell 3: ✓ — (no output)
Cell 4: ✓ — Out[167]= \begin{pmatrix}5 & 5 & 0 & 0\\5 & 4 & 1 & 0\\5 & 3 & 2 & 0\\4 & 4 & 1 & 1\\4 & 3 & 2 & 1\\3 & 3 & 2 & 2\end{pmatrix}
Cell 5: ✓ — (no output)
Cell 6: ✓ — (no output)
Cell 7: ✓ — Out[180]= \left\{\left\{1\right\},\left\{0,-1,1\right\},\left\{0,0,0,1\right\},\left\{0,0,0,0,-\frac{5}{2},1\right\}\right\}
Cell 8: ✓ — Out[183]= \begin{array}{|c|c|c|}\hline\left\{1,1\right\} & 3 & \mathrm{True}\\\hline\left\{1,2\right\} & 3 & \mathrm{True}\\\hline\left\{1,3\right\} & 3 & \mathrm{True}\\\hline\left\{2,1\right\} & 3 & \mathrm{True}\\\hline\left\{2,2\right\} & 6 & \mathrm{True}\\\hline\left\{2,3\right\} & 6 & \mathrm{True}\\\hline\left\{3,1\right\} & 3 & \mathrm{True}\\\hline\left\{3,2\right\} & 6 & \mathrm{True}\\… [464 chars total]
Cell 9: ✓ — Print: # irreps:[formula]   all solve:[formula]   total time (s):[formula] | Out[187]= \begin{array}{|c|c|c|c|c|}\hline\text{mu} & \text{dim} & \text{deg q\_k} & \text{time (s)} & \text{solves?}\\\hline16,16,0,0 & 8721 & 0,1,18,19 & \text{0.0030} & \mathrm{True}\\\hline16,15,1,0 & 27455 & 0,2,17,19 & \text{0.0029} & \mathrm{True}\\\hline16,14,2,0 & 47424 & 0,3,16,19 & \text{0.0031} & \mathrm{True}… [1324 chars total]

## 2026-07-21 06:34:35 — ✏️ EDIT CELL 6
(* --- Gaudin ODE --- 
   Q'''' coefficient is
   8(z-w1)^2(z-w2)^2 -- double singular points at the two Gaudin sites, as expected. C2,C4 are
   the abstract Casimirs of the output irrep mu; S,P label…


## 2026-07-21 06:34:35
**Input:**
```wolfram
(* --- Gaudin ODE --- 
   Q'''' coefficient is
   8(z-w1)^2(z-w2)^2 -- double singular points at the two Gaudin sites, as expected. C2,C4 are
   the abstract Casimirs of the output irrep mu; S,P label the two input reps [S,S,0,0],[P,P,0,0]*)
ClearAll[gaudinODE];
gaudinODE = (2*C2 + C2^2 - 2*C4 + 24*P + 2*C2*P + 20*P^2 + 4*P^3 + 24*S + 2*C2*S + 40*P*S + 12*P^2*S + 20*S^2 + 12*P*S^2 + 4*S^3)*Q[z]
   + 4*(C2 - 4*P - 2*P^2 - 4*S - 4*P*S - 2*S^2)*(2*z + P*z + S*z - w1 - P*w1 - w2 - S*w2)*Derivative[1][Q][z]
   - 4*(C2*z^2 - 6*P*z^2 - 4*P^2*z^2 - 6*S*z^2 - 8*P*S*z^2 - 4*S^2*z^2 - C2*z*w1 + 8*P*z*w1 + 6*P^2*z*w1 + 4*S*z*w1 + 8*P*S*z*w1
      + 2*S^2*z*w1 - 2*P*w1^2 - 2*P^2*w1^2 - C2*z*w2 + 4*P*z*w2 + 2*P^2*z*w2 + 8*S*z*w2 + 8*P*S*z*w2 + 6*S^2*z*w2
      + C2*w1*w2 - 4*P*w1*w2 - 2*P^2*w1*w2 - 4*S*w1*w2 - 8*P*S*w1*w2 - 2*S^2*w1*w2 - 2*S*w2^2 - 2*S^2*w2^2)*Derivative[2][Q][z]
   - 16*(z - w1)*(z - w2)*(P*z + S*z - P*w1 - S*w2)*Derivative[3][Q][z]
   + 8*(z - w1)^2*(z - w2)^2*Derivative[4][Q][z];
(* punctures fixed at w1=0, w2=1 for this solver *)
wsub = {w1 -> 0, w2 -> 1};
```
**Output:**
```
Out= Null
```
## 2026-07-21 06:34:44 — ✏️ EDIT CELL 7
(* --- Gaudin ODE solver ---
   Input: irrep mu in [S,S,0,0]x[P,P,0,0]. Solves the 4th-order Gaudin ODE for the four
   polynomial Q's of degrees degQ[mu,S,P], returned as explicit polynomials in z (n…


## 2026-07-21 06:34:44
**Input:**
```wolfram
(* --- Gaudin ODE solver ---
   Input: irrep mu in [S,S,0,0]x[P,P,0,0]. Solves the 4th-order Gaudin ODE for the four
   polynomial Q's of degrees degQ[mu,S,P], returned as explicit polynomials in z (not
   coefficient lists). GAUGE FIXING: the four Q-degrees are strictly increasing, so any
   lower-degree solution can be added to a higher one (q_k -> q_k + c q_j, j<k); removed
   by forcing the coefficient of z^(deg q_j) in q_k to vanish for every j<k. *)
ClearAll[GaudinSolve];
GaudinSolve[mu_, SS_, PP_] := Module[{ds, ode},
   ds = degQ[mu, SS, PP];
   ode = gaudinODE /. {C2 -> Cas2[mu], C4 -> Cas4[mu], S -> SS, P -> PP} /. wsub;
   Table[Module[{m = ds[[kk]], gauge = ds[[1 ;; kk - 1]], q, cl, sol},
       q = Sum[cq[j] z^j, {j, 0, m}];
       cl = CoefficientList[Expand[ode /. Q -> Function[x, q /. z -> x]], z];
       sol = Solve[Join[Thread[cl == 0], {cq[m] == 1}, Thread[(cq /@ gauge) == 0]], Table[cq[j], {j, 0, m}]];
       q /. sol[[1]]], {kk, 4}]];
(* example: mu={2,1,1,0} at (S,P)=(1,1) -> {1, z(z-1), z^3, z^5-(5/2)z^4} *)
GaudinSolve[{2, 1, 1, 0}, 1, 1]
```
**Output:**
```
Out= {1, -z + z^2, z^3, (-5*z^4)/2 + z^5}
```
## 2026-07-21 06:34:52 — ✏️ EDIT CELL 8
(* --- validation: every gauge-fixed Q_k(z) solves the Gaudin ODE exactly, across full (S,P)
   decompositions. GaudinSolve now returns explicit polynomials in z, so we substitute directly
   via Q ->…


## 2026-07-21 06:34:53
**Input:**
```wolfram
(* --- validation: every gauge-fixed Q_k(z) solves the Gaudin ODE exactly, across full (S,P)
   decompositions. GaudinSolve now returns explicit polynomials in z, so we substitute directly
   via Q -> Function[x, poly /. z -> x] (safe: this is the same pattern verified earlier to
   avoid the Dot-reconstruction/derivative edge case for degree-0 solutions). *)
ClearAll[verifyGaudin];
verifyGaudin[mu_, SS_, PP_] := Module[{qs, ode},
   qs = GaudinSolve[mu, SS, PP];
   ode = gaudinODE /. {C2 -> Cas2[mu], C4 -> Cas4[mu], S -> SS, P -> PP} /. wsub;
   AllTrue[qs, Expand[ode /. Q -> Function[x, # /. z -> x]] === 0 &]];
Grid[Table[{{SS, PP}, Length[irreps[SS, PP]], AllTrue[irreps[SS, PP], verifyGaudin[#, SS, PP] &]},
   {SS, {1, 2, 3}}, {PP, {1, 2, 3}}] // Flatten[#, 1] &, Frame -> All]
```
**Output:**
```
Out= Grid[{{{1, 1}, 3, True}, {{1, 2}, 3, True}, {{1, 3}, 3, True}, {{2, 1}, 3, True}, {{2, 2}, 6, True}, {{2, 3}, 6, True}, {{3, 1}, 3, True}, {{3, 2}, 6, True}, {{3, 3}, 10, True}}, Frame -> All]
```
## 2026-07-21 06:35:00 — ✏️ EDIT CELL 9
(* --- (S,P)=(8,8): Gaudin ODE solver over the full [8,8,0,0]x[8,8,0,0] decomposition (45 irreps),
   w1=0, w2=1, per-irrep timing. Each Q_k(z) verified to solve the ODE exactly. *)
irr88 = irreps[8, …


## 2026-07-21 06:35:02
**Input:**
```wolfram
(* --- (S,P)=(8,8): Gaudin ODE solver over the full [8,8,0,0]x[8,8,0,0] decomposition (45 irreps),
   w1=0, w2=1, per-irrep timing. Each Q_k(z) verified to solve the ODE exactly. *)
irr88 = irreps[8, 8];
tblG88 = Table[
   Module[{mu = irr88[[i]], tim, qs, ds, ode, ok},
     ds = degQ[mu, 8, 8];
     {tim, qs} = AbsoluteTiming[GaudinSolve[mu, 8, 8]];
     ode = gaudinODE /. {C2 -> Cas2[mu], C4 -> Cas4[mu], S -> 8, P -> 8} /. wsub;
     ok = AllTrue[qs, Expand[ode /. Q -> Function[x, # /. z -> x]] === 0 &];
     {Row[mu, ","], dimGL4[mu], Row[ds, ","], NumberForm[tim, {6, 4}], ok}],
   {i, Length[irr88]}];
Print["# irreps: ", Length[irr88], "   all solve: ", AllTrue[tblG88, Last],
   "   total time (s): ", Total[tblG88[[All, 4]] /. NumberForm[x_, _] :> x]];
Grid[Prepend[tblG88, {"mu", "dim", "deg q_k", "time (s)", "solves?"}],
  Frame -> All, Alignment -> Left,
  Background -> {None, {RGBColor[.88, .95, .88], {GrayLevel[1], RGBColor[.96, .99, .96]}}}]
```
**Output:**
```
Out= Grid[{{"mu", "dim", "deg q_k", "time (s)", "solves?"}, {Row[{16, 16, 0, 0}, ","], 8721, Row[{0, 1, 18, 19}, ","], NumberForm[0.010311, {6, 4}], True}, {Row[{16, 15, 1, 0}, ","], 27455, Row[{0, 2, 17, 19}, ","], NumberForm[0.003242, {6, 4}], True}, {Row[{16, 14, 2, 0}, ","], 47424, Row[{0, 3, 16, 19}, ","], NumberForm[0.003294, {6, 4}], True}, {Row[{16, 13, 3, 0}, ","], 62700, Row[{0, 4, 15, 19}, ","], NumberForm[0.003167, {6, 4}], True}, {Row[{16, 12, 4, 0}, ","], 69825, Row[{0, 5, 14, 19}, ","], NumberForm[0.003937, {6, 4}], True}, {Row[{16, 11, 5, 0}, ","], 67431, Row[{0, 6, 13, 19}, ","], NumberForm[0.003572, {6, 4}], True}, {Row[{16, 10, 6, 0}, ","], 55860, Row[{0, 7, 12, 19}, ","], NumberForm[0.003507, {6, 4}], True}, {Row[{16, 9, 7, 0}, ","], 36784, Row[{0, 8, 11, 19}, ","], NumberForm[0.003086, {6, 4}], True}, {Row[{16, 8, 8, 0}, ","], 12825, Row[{0, 9, 10, 19}, ","], NumberForm[0.002871, {6, 4}], True}, {Row[{15, 15, 1, 1}, ","], 5440, Row[{1, 2, 17, 18}, ","], NumberForm[0.00315, {6, 4}], True}, {Row[{15, 14, 2, 1}, ","], 16575, Row[{1, 3, 16, 18}, ","], NumberForm[0.002907, {6, 4}], True}, {Row[{15, 13, 3, 1}, ","], 27489, Row[{1, 4, 15, 18}, ","], NumberForm[0.002881, {6, 4}], True}, {Row[{15, 12, 4, 1}, ","], 34476, Row[{1, 5, 14, 18}, ","], NumberForm[0.003072, {6, 4}], True}, {Row[{15, 11, 5, 1}, ","], 35700, Row[{1, 6, 13, 18}, ","], NumberForm[0.003385, {6, 4}], True}, {Row[{15, 10, 6, 1}, ","], 30855, Row[{1, 7, 12, 18}, ","], NumberForm[0.003086, {6, 4}], True}, {Row[{15, 9, 7, 1}, ","], 20825, Row[{1, 8, 11, 18}, ","], NumberForm[0.003106, {6, 4}], True}, {Row[{15, 8, 8, 1}, ","], 7344, Row[{1, 9, 10, 18}, ","], NumberForm[0.003082, {6, 4}], True}, {Row[{14, 14, 2, 2}, ","], 3185, Row[{2, 3, 16, 17}, ","], NumberForm[0.004673, {6, 4}], True}, {Row[{14, 13, 3, 2}, ","], 9295, Row[{2, 4, 15, 17}, ","], NumberForm[0.006273, {6, 4}], True}, {Row[{14, 12, 4, 2}, ","], 14580, Row[{2, 5, 14, 17}, ","], NumberForm[0.005854, {6, 4}], True}, {Row[{14, 11, 5, 2}, ","], 16940, Row[{2, 6, 13, 17}, ","], NumberForm[0.006424, {6, 4}], True}, {Row[{14, 10, 6, 2}, ","], 15625, Row[{2, 7, 12, 17}, ","], NumberForm[0.005251, {6, 4}], True}, {Row[{14, 9, 7, 2}, ","], 10935, Row[{2, 8, 11, 17}, ","], NumberForm[0.003925, {6, 4}], True}, {Row[{14, 8, 8, 2}, ","], 3920, Row[{2, 9, 10, 17}, ","], NumberForm[0.00335, {6, 4}], True}, {Row[{13, 13, 3, 3}, ","], 1716, Row[{3, 4, 15, 16}, ","], NumberForm[0.002889, {6, 4}], True}, {Row[{13, 12, 4, 3}, ","], 4719, Row[{3, 5, 14, 16}, ","], NumberForm[0.00302, {6, 4}], True}, {Row[{13, 11, 5, 3}, ","], 6825, Row[{3, 6, 13, 16}, ","], NumberForm[0.003021, {6, 4}], True}, {Row[{13, 10, 6, 3}, ","], 7020, Row[{3, 7, 12, 16}, ","], NumberForm[0.003018, {6, 4}], True}, {Row[{13, 9, 7, 3}, ","], 5200, Row[{3, 8, 11, 16}, ","], NumberForm[0.003127, {6, 4}], True}, {Row[{13, 8, 8, 3}, ","], 1911, Row[{3, 9, 10, 16}, ","], NumberForm[0.00286, {6, 4}], True}, {Row[{12, 12, 4, 4}, ","], 825, Row[{4, 5, 14, 15}, ","], NumberForm[0.003007, {6, 4}], True}, {Row[{12, 11, 5, 4}, ","], 2079, Row[{4, 6, 13, 15}, ","], NumberForm[0.003056, {6, 4}], True}, {Row[{12, 10, 6, 4}, ","], 2640, Row[{4, 7, 12, 15}, ","], NumberForm[0.002996, {6, 4}], True}, {Row[{12, 9, 7, 4}, ","], 2156, Row[{4, 8, 11, 15}, ","], NumberForm[0.005316, {6, 4}], True}, {Row[{12, 8, 8, 4}, ","], 825, Row[{4, 9, 10, 15}, ","], NumberForm[0.005183, {6, 4}], True}, {Row[{11, 11, 5, 5}, ","], 336, Row[{5, 6, 13, 14}, ","], NumberForm[0.003386, {6, 4}], True}, {Row[{11, 10, 6, 5}, ","], 735, Row[{5, 7, 12, 14}, ","], NumberForm[0.003165, {6, 4}], True}, {Row[{11, 9, 7, 5}, ","], 729, Row[{5, 8, 11, 14}, ","], NumberForm[0.003167, {6, 4}], True}, {Row[{11, 8, 8, 5}, ","], 300, Row[{5, 9, 10, 14}, ","], NumberForm[0.002977, {6, 4}], True}, {Row[{10, 10, 6, 6}, ","], 105, Row[{6, 7, 12, 13}, ","], NumberForm[0.00299, {6, 4}], True}, {Row[{10, 9, 7, 6}, ","], 175, Row[{6, 8, 11, 13}, ","], NumberForm[0.002906, {6, 4}], True}, {Row[{10, 8, 8, 6}, ","], 84, Row[{6, 9, 10, 13}, ","], NumberForm[0.002946, {6, 4}], True}, {Row[{9, 9, 7, 7}, "
[output truncated — 4478 chars total]
```
## 2026-07-21 06:35:11 — ▶️ RUN CELLS 3–9
Cell 3: ✓ — (no output)
Cell 4: ✓ — Out[209]= \begin{pmatrix}5 & 5 & 0 & 0\\5 & 4 & 1 & 0\\5 & 3 & 2 & 0\\4 & 4 & 1 & 1\\4 & 3 & 2 & 1\\3 & 3 & 2 & 2\end{pmatrix}
Cell 5: ✓ — (no output)
Cell 6: ✓ — (no output)
Cell 7: ✓ — Out[222]= \left\{1,z^2-z,z^3,z^5-\frac{5\,z^4}{2}\right\}
Cell 8: ✓ — Out[225]= \begin{array}{|c|c|c|}\hline\left\{1,1\right\} & 3 & \mathrm{True}\\\hline\left\{1,2\right\} & 3 & \mathrm{True}\\\hline\left\{1,3\right\} & 3 & \mathrm{True}\\\hline\left\{2,1\right\} & 3 & \mathrm{True}\\\hline\left\{2,2\right\} & 6 & \mathrm{True}\\\hline\left\{2,3\right\} & 6 & \mathrm{True}\\\hline\left\{3,1\right\} & 3 & \mathrm{True}\\\hline\left\{3,2\right\} & 6 & \mathrm{True}\\… [464 chars total]
Cell 9: ✓ — Print: # irreps:[formula]   all solve:[formula]   total time (s):[formula] | Out[229]= \begin{array}{|c|c|c|c|c|}\hline\text{mu} & \text{dim} & \text{deg q\_k} & \text{time (s)} & \text{solves?}\\\hline16,16,0,0 & 8721 & 0,1,18,19 & \text{0.0033} & \mathrm{True}\\\hline16,15,1,0 & 27455 & 0,2,17,19 & \text{0.0030} & \mathrm{True}\\\hline16,14,2,0 & 47424 & 0,3,16,19 & \text{0.0030} & \mathrm{True}… [1324 chars total]

## 2026-07-21 06:38:02 — ✏️ EDIT CELL 7
(* --- Gaudin ODE solver ---
   Input: irrep mu in [S,S,0,0]x[P,P,0,0]. Solves the 4th-order Gaudin ODE for the four
   polynomial Q's of degrees degQ[mu,S,P], returned as monic coefficient lists (coe…


## 2026-07-21 06:38:02
**Input:**
```wolfram
(* --- Gaudin ODE solver ---
   Input: irrep mu in [S,S,0,0]x[P,P,0,0]. Solves the 4th-order Gaudin ODE for the four
   polynomial Q's of degrees degQ[mu,S,P], returned as monic coefficient lists (coeffs of
   z^0,z^1,...,z^m). GAUGE FIXING: the four Q-degrees are strictly increasing, so any
   lower-degree solution can be added to a higher one (q_k -> q_k + c q_j, j<k); removed
   by forcing the coefficient of z^(deg q_j) in q_k to vanish for every j<k. *)
ClearAll[GaudinSolve];
GaudinSolve[mu_, SS_, PP_] := Module[{ds, ode},
   ds = degQ[mu, SS, PP];
   ode = gaudinODE /. {C2 -> Cas2[mu], C4 -> Cas4[mu], S -> SS, P -> PP} /. wsub;
   Table[Module[{m = ds[[kk]], gauge = ds[[1 ;; kk - 1]], q, cl, sol},
       q = Sum[cq[j] z^j, {j, 0, m}];
       cl = CoefficientList[Expand[ode /. Q -> Function[x, q /. z -> x]], z];
       sol = Solve[Join[Thread[cl == 0], {cq[m] == 1}, Thread[(cq /@ gauge) == 0]], Table[cq[j], {j, 0, m}]];
       Table[cq[j], {j, 0, m}] /. sol[[1]]], {kk, 4}]];
(* example: mu={2,1,1,0} at (S,P)=(1,1) -> {{1},{0,-1,1},{0,0,0,1},{0,0,0,0,-5/2,1}} *)
GaudinSolve[{2, 1, 1, 0}, 1, 1]
```
**Output:**
```
Out= {{1}, {0, -1, 1}, {0, 0, 0, 1}, {0, 0, 0, 0, -5/2, 1}}
```
## 2026-07-21 06:38:14 — 📥 BULK INSERT 1 CELL(S) at positions 8–8
1. [code] (* --- BuildQs: same arguments as GaudinSolve, but returns the four Q's as explicit↵   polynomials i


## 2026-07-21 06:38:15
**Input:**
```wolfram
(* --- BuildQs: same arguments as GaudinSolve, but returns the four Q's as explicit
   polynomials in z (rather than coefficient lists) -- thin convenience wrapper. *)
ClearAll[BuildQs];
BuildQs[mu_, SS_, PP_] := Module[{ds = degQ[mu, SS, PP], coeffs = GaudinSolve[mu, SS, PP]},
   Table[Sum[coeffs[[kk, j + 1]] z^j, {j, 0, ds[[kk]]}], {kk, 4}]];
(* example: mu={2,1,1,0} at (S,P)=(1,1) -> {1, z^2-z, z^3, z^5-(5/2)z^4} *)
BuildQs[{2, 1, 1, 0}, 1, 1]
```
**Output:**
```
Out[232]= \left\{1,z^2-z,z^3,z^5-\frac{5\,z^4}{2}\right\}
```
## 2026-07-21 06:38:35 — ✏️ EDIT CELL 9
(* --- validation: every gauge-fixed Q_k(z) (via BuildQs) solves the Gaudin ODE exactly, across
   full (S,P) decompositions. *)
ClearAll[verifyGaudin];
verifyGaudin[mu_, SS_, PP_] := Module[{qs, ode}…


## 2026-07-21 06:38:35
**Input:**
```wolfram
(* --- validation: every gauge-fixed Q_k(z) (via BuildQs) solves the Gaudin ODE exactly, across
   full (S,P) decompositions. *)
ClearAll[verifyGaudin];
verifyGaudin[mu_, SS_, PP_] := Module[{qs, ode},
   qs = BuildQs[mu, SS, PP];
   ode = gaudinODE /. {C2 -> Cas2[mu], C4 -> Cas4[mu], S -> SS, P -> PP} /. wsub;
   AllTrue[qs, Expand[ode /. Q -> Function[x, # /. z -> x]] === 0 &]];
Grid[Table[{{SS, PP}, Length[irreps[SS, PP]], AllTrue[irreps[SS, PP], verifyGaudin[#, SS, PP] &]},
   {SS, {1, 2, 3}}, {PP, {1, 2, 3}}] // Flatten[#, 1] &, Frame -> All]
```
**Output:**
```
Out= Grid[{{{1, 1}, 3, True}, {{1, 2}, 3, True}, {{1, 3}, 3, True}, {{2, 1}, 3, True}, {{2, 2}, 6, True}, {{2, 3}, 6, True}, {{3, 1}, 3, True}, {{3, 2}, 6, True}, {{3, 3}, 10, True}}, Frame -> All]
```
## 2026-07-21 06:38:46 — ✏️ EDIT CELL 10
(* --- (S,P)=(8,8): Gaudin ODE solver over the full [8,8,0,0]x[8,8,0,0] decomposition (45 irreps),
   w1=0, w2=1, per-irrep timing (via BuildQs). Each Q_k(z) verified to solve the ODE exactly. *)
irr8…


## 2026-07-21 06:38:48
**Input:**
```wolfram
(* --- (S,P)=(8,8): Gaudin ODE solver over the full [8,8,0,0]x[8,8,0,0] decomposition (45 irreps),
   w1=0, w2=1, per-irrep timing (via BuildQs). Each Q_k(z) verified to solve the ODE exactly. *)
irr88 = irreps[8, 8];
tblG88 = Table[
   Module[{mu = irr88[[i]], tim, qs, ds, ode, ok},
     ds = degQ[mu, 8, 8];
     {tim, qs} = AbsoluteTiming[BuildQs[mu, 8, 8]];
     ode = gaudinODE /. {C2 -> Cas2[mu], C4 -> Cas4[mu], S -> 8, P -> 8} /. wsub;
     ok = AllTrue[qs, Expand[ode /. Q -> Function[x, # /. z -> x]] === 0 &];
     {Row[mu, ","], dimGL4[mu], Row[ds, ","], NumberForm[tim, {6, 4}], ok}],
   {i, Length[irr88]}];
Print["# irreps: ", Length[irr88], "   all solve: ", AllTrue[tblG88, Last],
   "   total time (s): ", Total[tblG88[[All, 4]] /. NumberForm[x_, _] :> x]];
Grid[Prepend[tblG88, {"mu", "dim", "deg q_k", "time (s)", "solves?"}],
  Frame -> All, Alignment -> Left,
  Background -> {None, {RGBColor[.88, .95, .88], {GrayLevel[1], RGBColor[.96, .99, .96]}}}]
```
**Output:**
```
Out= Grid[{{"mu", "dim", "deg q_k", "time (s)", "solves?"}, {Row[{16, 16, 0, 0}, ","], 8721, Row[{0, 1, 18, 19}, ","], NumberForm[0.002974, {6, 4}], True}, {Row[{16, 15, 1, 0}, ","], 27455, Row[{0, 2, 17, 19}, ","], NumberForm[0.003025, {6, 4}], True}, {Row[{16, 14, 2, 0}, ","], 47424, Row[{0, 3, 16, 19}, ","], NumberForm[0.003005, {6, 4}], True}, {Row[{16, 13, 3, 0}, ","], 62700, Row[{0, 4, 15, 19}, ","], NumberForm[0.003052, {6, 4}], True}, {Row[{16, 12, 4, 0}, ","], 69825, Row[{0, 5, 14, 19}, ","], NumberForm[0.003033, {6, 4}], True}, {Row[{16, 11, 5, 0}, ","], 67431, Row[{0, 6, 13, 19}, ","], NumberForm[0.00311, {6, 4}], True}, {Row[{16, 10, 6, 0}, ","], 55860, Row[{0, 7, 12, 19}, ","], NumberForm[0.003085, {6, 4}], True}, {Row[{16, 9, 7, 0}, ","], 36784, Row[{0, 8, 11, 19}, ","], NumberForm[0.003008, {6, 4}], True}, {Row[{16, 8, 8, 0}, ","], 12825, Row[{0, 9, 10, 19}, ","], NumberForm[0.00278, {6, 4}], True}, {Row[{15, 15, 1, 1}, ","], 5440, Row[{1, 2, 17, 18}, ","], NumberForm[0.002937, {6, 4}], True}, {Row[{15, 14, 2, 1}, ","], 16575, Row[{1, 3, 16, 18}, ","], NumberForm[0.002893, {6, 4}], True}, {Row[{15, 13, 3, 1}, ","], 27489, Row[{1, 4, 15, 18}, ","], NumberForm[0.002888, {6, 4}], True}, {Row[{15, 12, 4, 1}, ","], 34476, Row[{1, 5, 14, 18}, ","], NumberForm[0.00298, {6, 4}], True}, {Row[{15, 11, 5, 1}, ","], 35700, Row[{1, 6, 13, 18}, ","], NumberForm[0.002933, {6, 4}], True}, {Row[{15, 10, 6, 1}, ","], 30855, Row[{1, 7, 12, 18}, ","], NumberForm[0.002939, {6, 4}], True}, {Row[{15, 9, 7, 1}, ","], 20825, Row[{1, 8, 11, 18}, ","], NumberForm[0.003132, {6, 4}], True}, {Row[{15, 8, 8, 1}, ","], 7344, Row[{1, 9, 10, 18}, ","], NumberForm[0.003066, {6, 4}], True}, {Row[{14, 14, 2, 2}, ","], 3185, Row[{2, 3, 16, 17}, ","], NumberForm[0.003126, {6, 4}], True}, {Row[{14, 13, 3, 2}, ","], 9295, Row[{2, 4, 15, 17}, ","], NumberForm[0.00284, {6, 4}], True}, {Row[{14, 12, 4, 2}, ","], 14580, Row[{2, 5, 14, 17}, ","], NumberForm[0.002941, {6, 4}], True}, {Row[{14, 11, 5, 2}, ","], 16940, Row[{2, 6, 13, 17}, ","], NumberForm[0.00296, {6, 4}], True}, {Row[{14, 10, 6, 2}, ","], 15625, Row[{2, 7, 12, 17}, ","], NumberForm[0.002854, {6, 4}], True}, {Row[{14, 9, 7, 2}, ","], 10935, Row[{2, 8, 11, 17}, ","], NumberForm[0.003045, {6, 4}], True}, {Row[{14, 8, 8, 2}, ","], 3920, Row[{2, 9, 10, 17}, ","], NumberForm[0.003042, {6, 4}], True}, {Row[{13, 13, 3, 3}, ","], 1716, Row[{3, 4, 15, 16}, ","], NumberForm[0.002838, {6, 4}], True}, {Row[{13, 12, 4, 3}, ","], 4719, Row[{3, 5, 14, 16}, ","], NumberForm[0.002972, {6, 4}], True}, {Row[{13, 11, 5, 3}, ","], 6825, Row[{3, 6, 13, 16}, ","], NumberForm[0.003197, {6, 4}], True}, {Row[{13, 10, 6, 3}, ","], 7020, Row[{3, 7, 12, 16}, ","], NumberForm[0.00313, {6, 4}], True}, {Row[{13, 9, 7, 3}, ","], 5200, Row[{3, 8, 11, 16}, ","], NumberForm[0.002998, {6, 4}], True}, {Row[{13, 8, 8, 3}, ","], 1911, Row[{3, 9, 10, 16}, ","], NumberForm[0.003043, {6, 4}], True}, {Row[{12, 12, 4, 4}, ","], 825, Row[{4, 5, 14, 15}, ","], NumberForm[0.002953, {6, 4}], True}, {Row[{12, 11, 5, 4}, ","], 2079, Row[{4, 6, 13, 15}, ","], NumberForm[0.003015, {6, 4}], True}, {Row[{12, 10, 6, 4}, ","], 2640, Row[{4, 7, 12, 15}, ","], NumberForm[0.002933, {6, 4}], True}, {Row[{12, 9, 7, 4}, ","], 2156, Row[{4, 8, 11, 15}, ","], NumberForm[0.00283, {6, 4}], True}, {Row[{12, 8, 8, 4}, ","], 825, Row[{4, 9, 10, 15}, ","], NumberForm[0.003071, {6, 4}], True}, {Row[{11, 11, 5, 5}, ","], 336, Row[{5, 6, 13, 14}, ","], NumberForm[0.003136, {6, 4}], True}, {Row[{11, 10, 6, 5}, ","], 735, Row[{5, 7, 12, 14}, ","], NumberForm[0.003488, {6, 4}], True}, {Row[{11, 9, 7, 5}, ","], 729, Row[{5, 8, 11, 14}, ","], NumberForm[0.003086, {6, 4}], True}, {Row[{11, 8, 8, 5}, ","], 300, Row[{5, 9, 10, 14}, ","], NumberForm[0.003196, {6, 4}], True}, {Row[{10, 10, 6, 6}, ","], 105, Row[{6, 7, 12, 13}, ","], NumberForm[0.00309, {6, 4}], True}, {Row[{10, 9, 7, 6}, ","], 175, Row[{6, 8, 11, 13}, ","], NumberForm[0.003093, {6, 4}], True}, {Row[{10, 8, 8, 6}, ","], 84, Row[{6, 9, 10, 13}, ","], NumberForm[0.003019, {6, 4}], True}, {Row[{9, 9, 7, 7}, ","]
[output truncated — 4477 chars total]
```
## 2026-07-21 06:38:59 — ▶️ RUN CELLS 3–10
Cell 3: ✓ — (no output)
Cell 4: ✓ — Out[238]= \begin{pmatrix}5 & 5 & 0 & 0\\5 & 4 & 1 & 0\\5 & 3 & 2 & 0\\4 & 4 & 1 & 1\\4 & 3 & 2 & 1\\3 & 3 & 2 & 2\end{pmatrix}
Cell 5: ✓ — (no output)
Cell 6: ✓ — (no output)
Cell 7: ✓ — Out[251]= \left\{\left\{1\right\},\left\{0,-1,1\right\},\left\{0,0,0,1\right\},\left\{0,0,0,0,-\frac{5}{2},1\right\}\right\}
Cell 8: ✓ — Out[254]= \left\{1,z^2-z,z^3,z^5-\frac{5\,z^4}{2}\right\}
Cell 9: ✓ — Out[257]= \begin{array}{|c|c|c|}\hline\left\{1,1\right\} & 3 & \mathrm{True}\\\hline\left\{1,2\right\} & 3 & \mathrm{True}\\\hline\left\{1,3\right\} & 3 & \mathrm{True}\\\hline\left\{2,1\right\} & 3 & \mathrm{True}\\\hline\left\{2,2\right\} & 6 & \mathrm{True}\\\hline\left\{2,3\right\} & 6 & \mathrm{True}\\\hline\left\{3,1\right\} & 3 & \mathrm{True}\\\hline\left\{3,2\right\} & 6 & \mathrm{True}\\… [464 chars total]
Cell 10: ✓ — Print: # irreps:[formula]   all solve:[formula]   total time (s):[formula] | Out[261]= \begin{array}{|c|c|c|c|c|}\hline\text{mu} & \text{dim} & \text{deg q\_k} & \text{time (s)} & \text{solves?}\\\hline16,16,0,0 & 8721 & 0,1,18,19 & \text{0.0032} & \mathrm{True}\\\hline16,15,1,0 & 27455 & 0,2,17,19 & \text{0.0032} & \mathrm{True}\\\hline16,14,2,0 & 47424 & 0,3,16,19 & \text{0.0032} & \mathrm{True}… [1324 chars total]


---

## 2026-07-21 06:39:22 — 🔄 KERNEL RESTART


---

## 2026-07-21 06:39:45 — 🔄 KERNEL RESTART

## 2026-07-21 06:45:39 — ✏️ EDIT CELL 5
(* --- Casimir eigenvalues on the output irrep mu, and the Q-polynomial degrees ---
   Column-ordered generating function (Kazakov et al conventions): G4[u] = Prod_i (u - h(mu_i+4-i)),
   whose coeffi…


## 2026-07-21 06:45:39
**Input:**
```wolfram
(* --- Casimir eigenvalues on the output irrep mu, and the Q-polynomial degrees ---
   Column-ordered generating function (Kazakov et al conventions): G4[u] = Prod_i (u - h(mu_i+4-i)),
   whose coefficients are the Casimirs. Eigenvalues are symmetric functions of the centered
   shifted weights xw_i = mu_i + 5/2 - i (= mu_i + rho_i); Pw_k = Sum_i xw_i^k. *)
ClearAll[muhat, xw, Pw, Cas2, Cas4, degQ];
muhat[mu_] := Table[mu[[i]] + 4 - i, {i, 4}];
xw[mu_] := Table[mu[[i]] + 5/2 - i, {i, 4}];
Pw[mu_, k_] := Total[xw[mu]^k];
Cas2[mu_] := Pw[mu, 2] - 5;
Cas4[mu_] := Pw[mu, 4] + 4 Pw[mu, 3] - Pw[mu, 1] Pw[mu, 2] + Pw[mu, 2]/2 - 2 Pw[mu, 1]^2 - 10 Pw[mu, 1] - 51/4;
(* the SAME shifted weights fix the Q-degrees: deg q_k = (S+P+3) - muhat_k = (S+P) - mu_k + (k-1)
   -- convention: arguments always ordered (S,P,mu) *)
degQ[SS_, PP_, mu_] := Table[(SS + PP + 3) - muhat[mu][[k]], {k, 4}]
```
**Output:**
```
Out= Null
```
## 2026-07-21 06:45:50 — ✏️ EDIT CELL 7
(* --- Gaudin ODE solver ---
   Arguments always ordered (S,P,mu): S,P label the two input reps [S,S,0,0],[P,P,0,0]; mu is
   the output irrep. Solves the 4th-order Gaudin ODE for the four polynomial …


## 2026-07-21 06:45:50
**Input:**
```wolfram
(* --- Gaudin ODE solver ---
   Arguments always ordered (S,P,mu): S,P label the two input reps [S,S,0,0],[P,P,0,0]; mu is
   the output irrep. Solves the 4th-order Gaudin ODE for the four polynomial Q's of degrees
   degQ[S,P,mu], returned as monic coefficient lists (coeffs of z^0,z^1,...,z^m). GAUGE FIXING:
   the four Q-degrees are strictly increasing, so any lower-degree solution can be added to a
   higher one (q_k -> q_k + c q_j, j<k); removed by forcing the coefficient of z^(deg q_j) in
   q_k to vanish for every j<k. *)
ClearAll[GaudinSolve];
GaudinSolve[SS_, PP_, mu_] := Module[{ds, ode},
   ds = degQ[SS, PP, mu];
   ode = gaudinODE /. {C2 -> Cas2[mu], C4 -> Cas4[mu], S -> SS, P -> PP} /. wsub;
   Table[Module[{m = ds[[kk]], gauge = ds[[1 ;; kk - 1]], q, cl, sol},
       q = Sum[cq[j] z^j, {j, 0, m}];
       cl = CoefficientList[Expand[ode /. Q -> Function[x, q /. z -> x]], z];
       sol = Solve[Join[Thread[cl == 0], {cq[m] == 1}, Thread[(cq /@ gauge) == 0]], Table[cq[j], {j, 0, m}]];
       Table[cq[j], {j, 0, m}] /. sol[[1]]], {kk, 4}]];
```
**Output:**
```
Out= Null
```
## 2026-07-21 06:45:53 — ✏️ EDIT CELL 8
(* example: (S,P)=(1,1), mu={2,1,1,0} *)
GaudinSolve[1, 1, {2, 1, 1, 0}]


## 2026-07-21 06:45:53
**Input:**
```wolfram
(* example: (S,P)=(1,1), mu={2,1,1,0} *)
GaudinSolve[1, 1, {2, 1, 1, 0}]
```
**Output:**
```
Out= {{1}, {0, -1, 1}, {0, 0, 0, 1}, {0, 0, 0, 0, -5/2, 1}}
```
## 2026-07-21 06:45:59 — ✏️ EDIT CELL 9
(* --- BuildQs: same argument convention as GaudinSolve, (S,P,mu), but returns the four Q's
   as explicit polynomials in z (rather than coefficient lists) *)
ClearAll[BuildQs];
BuildQs[SS_, PP_, mu_]…


## 2026-07-21 06:45:59
**Input:**
```wolfram
(* --- BuildQs: same argument convention as GaudinSolve, (S,P,mu), but returns the four Q's
   as explicit polynomials in z (rather than coefficient lists) *)
ClearAll[BuildQs];
BuildQs[SS_, PP_, mu_] := Module[{ds = degQ[SS, PP, mu], coeffs = GaudinSolve[SS, PP, mu]},
   Table[Sum[coeffs[[kk, j + 1]] z^j, {j, 0, ds[[kk]]}], {kk, 4}]];
```
**Output:**
```
Out= Null
```
## 2026-07-21 06:46:01 — ✏️ EDIT CELL 10
(* example: (S,P)=(1,1), mu={2,1,1,0} *)
BuildQs[1, 1, {2, 1, 1, 0}]


## 2026-07-21 06:46:01
**Input:**
```wolfram
(* example: (S,P)=(1,1), mu={2,1,1,0} *)
BuildQs[1, 1, {2, 1, 1, 0}]
```
**Output:**
```
Out= {1, -z + z^2, z^3, (-5*z^4)/2 + z^5}
```
## 2026-07-21 06:46:17 — ✏️ EDIT CELL 11
(* --- validation: every gauge-fixed Q_k(z) (via BuildQs) solves the Gaudin ODE exactly, across
   full (S,P) decompositions. Argument convention: (S,P,mu). *)
ClearAll[verifyGaudin];
verifyGaudin[SS_…


## 2026-07-21 06:46:17
**Input:**
```wolfram
(* --- validation: every gauge-fixed Q_k(z) (via BuildQs) solves the Gaudin ODE exactly, across
   full (S,P) decompositions. Argument convention: (S,P,mu). *)
ClearAll[verifyGaudin];
verifyGaudin[SS_, PP_, mu_] := Module[{qs, ode},
   qs = BuildQs[SS, PP, mu];
   ode = gaudinODE /. {C2 -> Cas2[mu], C4 -> Cas4[mu], S -> SS, P -> PP} /. wsub;
   AllTrue[qs, Expand[ode /. Q -> Function[x, # /. z -> x]] === 0 &]];
Grid[Table[{{SS, PP}, Length[irreps[SS, PP]], AllTrue[irreps[SS, PP], verifyGaudin[SS, PP, #] &]},
   {SS, {1, 2, 3}}, {PP, {1, 2, 3}}] // Flatten[#, 1] &, Frame -> All]
```
**Output:**
```
Out= Grid[{{{1, 1}, 3, True}, {{1, 2}, 3, True}, {{1, 3}, 3, True}, {{2, 1}, 3, True}, {{2, 2}, 6, True}, {{2, 3}, 6, True}, {{3, 1}, 3, True}, {{3, 2}, 6, True}, {{3, 3}, 10, True}}, Frame -> All]
```
## 2026-07-21 06:46:27 — ✏️ EDIT CELL 12
(* --- (S,P)=(8,8): Gaudin ODE solver over the full [8,8,0,0]x[8,8,0,0] decomposition (45 irreps),
   per-irrep timing (via BuildQs, so the check substitutes actual Q_k(z) polynomials -- checking
   t…


## 2026-07-21 06:46:29
**Input:**
```wolfram
(* --- (S,P)=(8,8): Gaudin ODE solver over the full [8,8,0,0]x[8,8,0,0] decomposition (45 irreps),
   per-irrep timing (via BuildQs, so the check substitutes actual Q_k(z) polynomials -- checking
   the raw GaudinSolve coefficient lists directly against the ODE is meaningless and spuriously
   fails everything). Argument convention: (S,P,mu). *)
irr88 = irreps[8, 8];
tblG88 = Table[
   Module[{mu = irr88[[i]], tim, qs, ds, ode, ok},
     ds = degQ[8, 8, mu];
     {tim, qs} = AbsoluteTiming[BuildQs[8, 8, mu]];
     ode = gaudinODE /. {C2 -> Cas2[mu], C4 -> Cas4[mu], S -> 8, P -> 8} /. wsub;
     ok = AllTrue[qs, Expand[ode /. Q -> Function[x, # /. z -> x]] === 0 &];
     {Row[mu, ","], dimGL4[mu], Row[ds, ","], NumberForm[tim, {6, 4}], ok}],
   {i, Length[irr88]}];
Print["# irreps: ", Length[irr88], "   all solve: ", AllTrue[tblG88, Last],
   "   total time (s): ", Total[tblG88[[All, 4]] /. NumberForm[x_, _] :> x]];
Grid[Prepend[tblG88, {"mu", "dim", "deg Q", "time (s)", "solves?"}],
  Frame -> All, Alignment -> Left,
  Background -> {None, {RGBColor[.88, .95, .88], {GrayLevel[1], RGBColor[.96, .99, .96]}}}]
```
**Output:**
```
Out= Grid[{{"mu", "dim", "deg Q", "time (s)", "solves?"}, {Row[{16, 16, 0, 0}, ","], 8721, Row[{0, 1, 18, 19}, ","], NumberForm[0.002853, {6, 4}], True}, {Row[{16, 15, 1, 0}, ","], 27455, Row[{0, 2, 17, 19}, ","], NumberForm[0.00301, {6, 4}], True}, {Row[{16, 14, 2, 0}, ","], 47424, Row[{0, 3, 16, 19}, ","], NumberForm[0.002981, {6, 4}], True}, {Row[{16, 13, 3, 0}, ","], 62700, Row[{0, 4, 15, 19}, ","], NumberForm[0.003077, {6, 4}], True}, {Row[{16, 12, 4, 0}, ","], 69825, Row[{0, 5, 14, 19}, ","], NumberForm[0.003147, {6, 4}], True}, {Row[{16, 11, 5, 0}, ","], 67431, Row[{0, 6, 13, 19}, ","], NumberForm[0.003286, {6, 4}], True}, {Row[{16, 10, 6, 0}, ","], 55860, Row[{0, 7, 12, 19}, ","], NumberForm[0.003092, {6, 4}], True}, {Row[{16, 9, 7, 0}, ","], 36784, Row[{0, 8, 11, 19}, ","], NumberForm[0.002928, {6, 4}], True}, {Row[{16, 8, 8, 0}, ","], 12825, Row[{0, 9, 10, 19}, ","], NumberForm[0.002767, {6, 4}], True}, {Row[{15, 15, 1, 1}, ","], 5440, Row[{1, 2, 17, 18}, ","], NumberForm[0.002728, {6, 4}], True}, {Row[{15, 14, 2, 1}, ","], 16575, Row[{1, 3, 16, 18}, ","], NumberForm[0.002735, {6, 4}], True}, {Row[{15, 13, 3, 1}, ","], 27489, Row[{1, 4, 15, 18}, ","], NumberForm[0.00276, {6, 4}], True}, {Row[{15, 12, 4, 1}, ","], 34476, Row[{1, 5, 14, 18}, ","], NumberForm[0.003026, {6, 4}], True}, {Row[{15, 11, 5, 1}, ","], 35700, Row[{1, 6, 13, 18}, ","], NumberForm[0.002957, {6, 4}], True}, {Row[{15, 10, 6, 1}, ","], 30855, Row[{1, 7, 12, 18}, ","], NumberForm[0.003072, {6, 4}], True}, {Row[{15, 9, 7, 1}, ","], 20825, Row[{1, 8, 11, 18}, ","], NumberForm[0.002976, {6, 4}], True}, {Row[{15, 8, 8, 1}, ","], 7344, Row[{1, 9, 10, 18}, ","], NumberForm[0.00292, {6, 4}], True}, {Row[{14, 14, 2, 2}, ","], 3185, Row[{2, 3, 16, 17}, ","], NumberForm[0.002703, {6, 4}], True}, {Row[{14, 13, 3, 2}, ","], 9295, Row[{2, 4, 15, 17}, ","], NumberForm[0.00277, {6, 4}], True}, {Row[{14, 12, 4, 2}, ","], 14580, Row[{2, 5, 14, 17}, ","], NumberForm[0.002948, {6, 4}], True}, {Row[{14, 11, 5, 2}, ","], 16940, Row[{2, 6, 13, 17}, ","], NumberForm[0.003016, {6, 4}], True}, {Row[{14, 10, 6, 2}, ","], 15625, Row[{2, 7, 12, 17}, ","], NumberForm[0.003202, {6, 4}], True}, {Row[{14, 9, 7, 2}, ","], 10935, Row[{2, 8, 11, 17}, ","], NumberForm[0.00317, {6, 4}], True}, {Row[{14, 8, 8, 2}, ","], 3920, Row[{2, 9, 10, 17}, ","], NumberForm[0.00307, {6, 4}], True}, {Row[{13, 13, 3, 3}, ","], 1716, Row[{3, 4, 15, 16}, ","], NumberForm[0.002745, {6, 4}], True}, {Row[{13, 12, 4, 3}, ","], 4719, Row[{3, 5, 14, 16}, ","], NumberForm[0.002982, {6, 4}], True}, {Row[{13, 11, 5, 3}, ","], 6825, Row[{3, 6, 13, 16}, ","], NumberForm[0.003176, {6, 4}], True}, {Row[{13, 10, 6, 3}, ","], 7020, Row[{3, 7, 12, 16}, ","], NumberForm[0.002898, {6, 4}], True}, {Row[{13, 9, 7, 3}, ","], 5200, Row[{3, 8, 11, 16}, ","], NumberForm[0.002823, {6, 4}], True}, {Row[{13, 8, 8, 3}, ","], 1911, Row[{3, 9, 10, 16}, ","], NumberForm[0.002921, {6, 4}], True}, {Row[{12, 12, 4, 4}, ","], 825, Row[{4, 5, 14, 15}, ","], NumberForm[0.002841, {6, 4}], True}, {Row[{12, 11, 5, 4}, ","], 2079, Row[{4, 6, 13, 15}, ","], NumberForm[0.003095, {6, 4}], True}, {Row[{12, 10, 6, 4}, ","], 2640, Row[{4, 7, 12, 15}, ","], NumberForm[0.002834, {6, 4}], True}, {Row[{12, 9, 7, 4}, ","], 2156, Row[{4, 8, 11, 15}, ","], NumberForm[0.002815, {6, 4}], True}, {Row[{12, 8, 8, 4}, ","], 825, Row[{4, 9, 10, 15}, ","], NumberForm[0.00304, {6, 4}], True}, {Row[{11, 11, 5, 5}, ","], 336, Row[{5, 6, 13, 14}, ","], NumberForm[0.002987, {6, 4}], True}, {Row[{11, 10, 6, 5}, ","], 735, Row[{5, 7, 12, 14}, ","], NumberForm[0.002922, {6, 4}], True}, {Row[{11, 9, 7, 5}, ","], 729, Row[{5, 8, 11, 14}, ","], NumberForm[0.003191, {6, 4}], True}, {Row[{11, 8, 8, 5}, ","], 300, Row[{5, 9, 10, 14}, ","], NumberForm[0.003258, {6, 4}], True}, {Row[{10, 10, 6, 6}, ","], 105, Row[{6, 7, 12, 13}, ","], NumberForm[0.002895, {6, 4}], True}, {Row[{10, 9, 7, 6}, ","], 175, Row[{6, 8, 11, 13}, ","], NumberForm[0.002909, {6, 4}], True}, {Row[{10, 8, 8, 6}, ","], 84, Row[{6, 9, 10, 13}, ","], NumberForm[0.003333, {6, 4}], True}, {Row[{9, 9, 7, 7}, ","],
[output truncated — 4476 chars total]
```
## 2026-07-21 06:46:39 — ▶️ RUN CELLS 3–12
Cell 3: ✓ — (no output)
Cell 4: ✓ — Out[35]= \begin{pmatrix}5 & 5 & 0 & 0\\5 & 4 & 1 & 0\\5 & 3 & 2 & 0\\4 & 4 & 1 & 1\\4 & 3 & 2 & 1\\3 & 3 & 2 & 2\end{pmatrix}
Cell 5: ✓ — (no output)
Cell 6: ✓ — (no output)
Cell 7: ✓ — (no output)
Cell 8: ✓ — Out[48]= \left\{\left\{1\right\},\left\{0,-1,1\right\},\left\{0,0,0,1\right\},\left\{0,0,0,0,-\frac{5}{2},1\right\}\right\}
Cell 9: ✓ — (no output)
Cell 10: ✓ — Out[51]= \left\{1,z^2-z,z^3,z^5-\frac{5\,z^4}{2}\right\}
Cell 11: ✓ — Out[54]= \begin{array}{|c|c|c|}\hline\left\{1,1\right\} & 3 & \mathrm{True}\\\hline\left\{1,2\right\} & 3 & \mathrm{True}\\\hline\left\{1,3\right\} & 3 & \mathrm{True}\\\hline\left\{2,1\right\} & 3 & \mathrm{True}\\\hline\left\{2,2\right\} & 6 & \mathrm{True}\\\hline\left\{2,3\right\} & 6 & \mathrm{True}\\\hline\left\{3,1\right\} & 3 & \mathrm{True}\\\hline\left\{3,2\right\} & 6 & \mathrm{True}\\\… [463 chars total]
Cell 12: ✓ — Print: # irreps:[formula]   all solve:[formula]   total time (s):[formula] | Out[58]= \begin{array}{|c|c|c|c|c|}\hline\text{mu} & \text{dim} & \text{deg Q} & \text{time (s)} & \text{solves?}\\\hline16,16,0,0 & 8721 & 0,1,18,19 & \text{0.0028} & \mathrm{True}\\\hline16,15,1,0 & 27455 & 0,2,17,19 & \text{0.0031} & \mathrm{True}\\\hline16,14,2,0 & 47424 & 0,3,16,19 & \text{0.0030} & \mathrm{True}\\\h… [1323 chars total]

## 2026-07-21 06:48:21 — ✏️ EDIT CELL 11
(* --- validation: every gauge-fixed Q_k solves the Gaudin ODE exactly, across full (S,P)
   decompositions. Works directly off GaudinSolve's coefficient lists (BuildQs is just a display
   convenienc…


## 2026-07-21 06:48:21
**Input:**
```wolfram
(* --- validation: every gauge-fixed Q_k solves the Gaudin ODE exactly, across full (S,P)
   decompositions. Works directly off GaudinSolve's coefficient lists (BuildQs is just a display
   convenience, not used for verification). NB: q is rebuilt as a Sum (SetDelayed) rather than a
   Dot[coeffs, z^Range[...]] pure Function -- the latter has a degenerate edge case when
   differentiating symbolically for length-1 coefficient lists (the deg=0 case) that spuriously
   fails. Argument convention: (S,P,mu). *)
ClearAll[verifyGaudin];
verifyGaudin[SS_, PP_, mu_] := Module[{qs, ds, ode, qf},
   ds = degQ[SS, PP, mu]; qs = GaudinSolve[SS, PP, mu];
   ode = gaudinODE /. {C2 -> Cas2[mu], C4 -> Cas4[mu], S -> SS, P -> PP} /. wsub;
   AllTrue[Range[4], Module[{m = ds[[#]], cf = qs[[#]]},
      qf[x_] := Sum[cf[[j + 1]] x^j, {j, 0, m}];
      Expand[ode /. Q -> qf] === 0] &]];
Grid[Table[{{SS, PP}, Length[irreps[SS, PP]], AllTrue[irreps[SS, PP], verifyGaudin[SS, PP, #] &]},
   {SS, {1, 2, 3}}, {PP, {1, 2, 3}}] // Flatten[#, 1] &, Frame -> All]
```
**Output:**
```
Out= Grid[{{{1, 1}, 3, True}, {{1, 2}, 3, True}, {{1, 3}, 3, True}, {{2, 1}, 3, True}, {{2, 2}, 6, True}, {{2, 3}, 6, True}, {{3, 1}, 3, True}, {{3, 2}, 6, True}, {{3, 3}, 10, True}}, Frame -> All]
```
## 2026-07-21 06:48:31 — ✏️ EDIT CELL 12
(* --- (S,P)=(8,8): Gaudin ODE solver over the full [8,8,0,0]x[8,8,0,0] decomposition (45 irreps),
   per-irrep timing. Works directly off GaudinSolve's coefficient lists (BuildQs is just a display
  …


## 2026-07-21 06:48:33
**Input:**
```wolfram
(* --- (S,P)=(8,8): Gaudin ODE solver over the full [8,8,0,0]x[8,8,0,0] decomposition (45 irreps),
   per-irrep timing. Works directly off GaudinSolve's coefficient lists (BuildQs is just a display
   convenience, not used here) -- q rebuilt as a Sum (SetDelayed), which avoids the degenerate
   Dot/derivative edge case for deg=0 solutions. Argument convention: (S,P,mu). *)
irr88 = irreps[8, 8];
tblG88 = Table[
   Module[{mu = irr88[[i]], tim, qs, ds, ode, ok, qf},
     ds = degQ[8, 8, mu];
     {tim, qs} = AbsoluteTiming[GaudinSolve[8, 8, mu]];
     ode = gaudinODE /. {C2 -> Cas2[mu], C4 -> Cas4[mu], S -> 8, P -> 8} /. wsub;
     ok = AllTrue[Range[4], Module[{m = ds[[#]], cf = qs[[#]]},
        qf[x_] := Sum[cf[[j + 1]] x^j, {j, 0, m}];
        Expand[ode /. Q -> qf] === 0] &];
     {Row[mu, ","], dimGL4[mu], Row[ds, ","], NumberForm[tim, {6, 4}], ok}],
   {i, Length[irr88]}];
Print["# irreps: ", Length[irr88], "   all solve: ", AllTrue[tblG88, Last],
   "   total time (s): ", Total[tblG88[[All, 4]] /. NumberForm[x_, _] :> x]];
Grid[Prepend[tblG88, {"mu", "dim", "deg Q", "time (s)", "solves?"}],
  Frame -> All, Alignment -> Left,
  Background -> {None, {RGBColor[.88, .95, .88], {GrayLevel[1], RGBColor[.96, .99, .96]}}}]
```
**Output:**
```
Out= Grid[{{"mu", "dim", "deg Q", "time (s)", "solves?"}, {Row[{16, 16, 0, 0}, ","], 8721, Row[{0, 1, 18, 19}, ","], NumberForm[0.003092, {6, 4}], True}, {Row[{16, 15, 1, 0}, ","], 27455, Row[{0, 2, 17, 19}, ","], NumberForm[0.002889, {6, 4}], True}, {Row[{16, 14, 2, 0}, ","], 47424, Row[{0, 3, 16, 19}, ","], NumberForm[0.00319, {6, 4}], True}, {Row[{16, 13, 3, 0}, ","], 62700, Row[{0, 4, 15, 19}, ","], NumberForm[0.003868, {6, 4}], True}, {Row[{16, 12, 4, 0}, ","], 69825, Row[{0, 5, 14, 19}, ","], NumberForm[0.003258, {6, 4}], True}, {Row[{16, 11, 5, 0}, ","], 67431, Row[{0, 6, 13, 19}, ","], NumberForm[0.004135, {6, 4}], True}, {Row[{16, 10, 6, 0}, ","], 55860, Row[{0, 7, 12, 19}, ","], NumberForm[0.003583, {6, 4}], True}, {Row[{16, 9, 7, 0}, ","], 36784, Row[{0, 8, 11, 19}, ","], NumberForm[0.002956, {6, 4}], True}, {Row[{16, 8, 8, 0}, ","], 12825, Row[{0, 9, 10, 19}, ","], NumberForm[0.003364, {6, 4}], True}, {Row[{15, 15, 1, 1}, ","], 5440, Row[{1, 2, 17, 18}, ","], NumberForm[0.003049, {6, 4}], True}, {Row[{15, 14, 2, 1}, ","], 16575, Row[{1, 3, 16, 18}, ","], NumberForm[0.002943, {6, 4}], True}, {Row[{15, 13, 3, 1}, ","], 27489, Row[{1, 4, 15, 18}, ","], NumberForm[0.00275, {6, 4}], True}, {Row[{15, 12, 4, 1}, ","], 34476, Row[{1, 5, 14, 18}, ","], NumberForm[0.00506, {6, 4}], True}, {Row[{15, 11, 5, 1}, ","], 35700, Row[{1, 6, 13, 18}, ","], NumberForm[0.00589, {6, 4}], True}, {Row[{15, 10, 6, 1}, ","], 30855, Row[{1, 7, 12, 18}, ","], NumberForm[0.005425, {6, 4}], True}, {Row[{15, 9, 7, 1}, ","], 20825, Row[{1, 8, 11, 18}, ","], NumberForm[0.005561, {6, 4}], True}, {Row[{15, 8, 8, 1}, ","], 7344, Row[{1, 9, 10, 18}, ","], NumberForm[0.004499, {6, 4}], True}, {Row[{14, 14, 2, 2}, ","], 3185, Row[{2, 3, 16, 17}, ","], NumberForm[0.008099, {6, 4}], True}, {Row[{14, 13, 3, 2}, ","], 9295, Row[{2, 4, 15, 17}, ","], NumberForm[0.005801, {6, 4}], True}, {Row[{14, 12, 4, 2}, ","], 14580, Row[{2, 5, 14, 17}, ","], NumberForm[0.006174, {6, 4}], True}, {Row[{14, 11, 5, 2}, ","], 16940, Row[{2, 6, 13, 17}, ","], NumberForm[0.005787, {6, 4}], True}, {Row[{14, 10, 6, 2}, ","], 15625, Row[{2, 7, 12, 17}, ","], NumberForm[0.003227, {6, 4}], True}, {Row[{14, 9, 7, 2}, ","], 10935, Row[{2, 8, 11, 17}, ","], NumberForm[0.003249, {6, 4}], True}, {Row[{14, 8, 8, 2}, ","], 3920, Row[{2, 9, 10, 17}, ","], NumberForm[0.003032, {6, 4}], True}, {Row[{13, 13, 3, 3}, ","], 1716, Row[{3, 4, 15, 16}, ","], NumberForm[0.002756, {6, 4}], True}, {Row[{13, 12, 4, 3}, ","], 4719, Row[{3, 5, 14, 16}, ","], NumberForm[0.002896, {6, 4}], True}, {Row[{13, 11, 5, 3}, ","], 6825, Row[{3, 6, 13, 16}, ","], NumberForm[0.003276, {6, 4}], True}, {Row[{13, 10, 6, 3}, ","], 7020, Row[{3, 7, 12, 16}, ","], NumberForm[0.00289, {6, 4}], True}, {Row[{13, 9, 7, 3}, ","], 5200, Row[{3, 8, 11, 16}, ","], NumberForm[0.002962, {6, 4}], True}, {Row[{13, 8, 8, 3}, ","], 1911, Row[{3, 9, 10, 16}, ","], NumberForm[0.002926, {6, 4}], True}, {Row[{12, 12, 4, 4}, ","], 825, Row[{4, 5, 14, 15}, ","], NumberForm[0.002941, {6, 4}], True}, {Row[{12, 11, 5, 4}, ","], 2079, Row[{4, 6, 13, 15}, ","], NumberForm[0.00296, {6, 4}], True}, {Row[{12, 10, 6, 4}, ","], 2640, Row[{4, 7, 12, 15}, ","], NumberForm[0.002834, {6, 4}], True}, {Row[{12, 9, 7, 4}, ","], 2156, Row[{4, 8, 11, 15}, ","], NumberForm[0.002767, {6, 4}], True}, {Row[{12, 8, 8, 4}, ","], 825, Row[{4, 9, 10, 15}, ","], NumberForm[0.002919, {6, 4}], True}, {Row[{11, 11, 5, 5}, ","], 336, Row[{5, 6, 13, 14}, ","], NumberForm[0.003744, {6, 4}], True}, {Row[{11, 10, 6, 5}, ","], 735, Row[{5, 7, 12, 14}, ","], NumberForm[0.004015, {6, 4}], True}, {Row[{11, 9, 7, 5}, ","], 729, Row[{5, 8, 11, 14}, ","], NumberForm[0.003199, {6, 4}], True}, {Row[{11, 8, 8, 5}, ","], 300, Row[{5, 9, 10, 14}, ","], NumberForm[0.003656, {6, 4}], True}, {Row[{10, 10, 6, 6}, ","], 105, Row[{6, 7, 12, 13}, ","], NumberForm[0.003619, {6, 4}], True}, {Row[{10, 9, 7, 6}, ","], 175, Row[{6, 8, 11, 13}, ","], NumberForm[0.005399, {6, 4}], True}, {Row[{10, 8, 8, 6}, ","], 84, Row[{6, 9, 10, 13}, ","], NumberForm[0.006333, {6, 4}], True}, {Row[{9, 9, 7, 7}, ","]
[output truncated — 4476 chars total]
```
## 2026-07-21 06:48:43 — ▶️ RUN CELLS 3–12
Cell 3: ✓ — (no output)
Cell 4: ✓ — Out[64]= \begin{pmatrix}5 & 5 & 0 & 0\\5 & 4 & 1 & 0\\5 & 3 & 2 & 0\\4 & 4 & 1 & 1\\4 & 3 & 2 & 1\\3 & 3 & 2 & 2\end{pmatrix}
Cell 5: ✓ — (no output)
Cell 6: ✓ — (no output)
Cell 7: ✓ — (no output)
Cell 8: ✓ — Out[77]= \left\{\left\{1\right\},\left\{0,-1,1\right\},\left\{0,0,0,1\right\},\left\{0,0,0,0,-\frac{5}{2},1\right\}\right\}
Cell 9: ✓ — (no output)
Cell 10: ✓ — Out[80]= \left\{1,z^2-z,z^3,z^5-\frac{5\,z^4}{2}\right\}
Cell 11: ✓ — Out[83]= \begin{array}{|c|c|c|}\hline\left\{1,1\right\} & 3 & \mathrm{True}\\\hline\left\{1,2\right\} & 3 & \mathrm{True}\\\hline\left\{1,3\right\} & 3 & \mathrm{True}\\\hline\left\{2,1\right\} & 3 & \mathrm{True}\\\hline\left\{2,2\right\} & 6 & \mathrm{True}\\\hline\left\{2,3\right\} & 6 & \mathrm{True}\\\hline\left\{3,1\right\} & 3 & \mathrm{True}\\\hline\left\{3,2\right\} & 6 & \mathrm{True}\\\… [463 chars total]
Cell 12: ✓ — Print: # irreps:[formula]   all solve:[formula]   total time (s):[formula] | Out[87]= \begin{array}{|c|c|c|c|c|}\hline\text{mu} & \text{dim} & \text{deg Q} & \text{time (s)} & \text{solves?}\\\hline16,16,0,0 & 8721 & 0,1,18,19 & \text{0.0030} & \mathrm{True}\\\hline16,15,1,0 & 27455 & 0,2,17,19 & \text{0.0033} & \mathrm{True}\\\hline16,14,2,0 & 47424 & 0,3,16,19 & \text{0.0031} & \mathrm{True}\\\h… [1323 chars total]

## 2026-07-21 08:02:25 — ✏️ EDIT CELL 12
(* --- CheckAllIrreps[S,P]: Gaudin ODE solver over the FULL [S,S,0,0]x[P,P,0,0] decomposition,
   per-irrep timing. Works directly off GaudinSolve's coefficient lists (BuildQs is just a display
   con…


## 2026-07-21 08:02:25
**Input:**
```wolfram
(* --- CheckAllIrreps[S,P]: Gaudin ODE solver over the FULL [S,S,0,0]x[P,P,0,0] decomposition,
   per-irrep timing. Works directly off GaudinSolve's coefficient lists (BuildQs is just a display
   convenience, not used here) -- q rebuilt as a Sum (SetDelayed), which avoids the degenerate
   Dot/derivative edge case for deg=0 solutions. Argument convention: (S,P). *)
ClearAll[CheckAllIrreps];
CheckAllIrreps[SS_, PP_] := Module[{irr = irreps[SS, PP], tbl},
   tbl = Table[
      Module[{mu = irr[[i]], tim, qs, ds, ode, ok, qf},
        ds = degQ[SS, PP, mu];
        {tim, qs} = AbsoluteTiming[GaudinSolve[SS, PP, mu]];
        ode = gaudinODE /. {C2 -> Cas2[mu], C4 -> Cas4[mu], S -> SS, P -> PP} /. wsub;
        ok = AllTrue[Range[4], Module[{m = ds[[#]], cf = qs[[#]]},
           qf[x_] := Sum[cf[[j + 1]] x^j, {j, 0, m}];
           Expand[ode /. Q -> qf] === 0] &];
        {Row[mu, ","], dimGL4[mu], Row[ds, ","], NumberForm[tim, {6, 4}], ok}],
      {i, Length[irr]}];
   Print["# irreps: ", Length[irr], "   all solve: ", AllTrue[tbl, Last],
      "   total time (s): ", Total[tbl[[All, 4]] /. NumberForm[x_, _] :> x]];
   Grid[Prepend[tbl, {"mu", "dim", "deg Q", "time (s)", "solves?"}],
     Frame -> All, Alignment -> Left,
     Background -> {None, {RGBColor[.88, .95, .88], {GrayLevel[1], RGBColor[.96, .99, .96]}}}]];
```
**Output:**
```
Out= Null
```
## 2026-07-21 08:02:29 — 📥 BULK INSERT 1 CELL(S) at positions 13–13
1. [code] (* example: full sweep at (S,P)=(8,8) *)↵CheckAllIrreps[8, 8]


## 2026-07-21 08:02:32
**Input:**
```wolfram
(* example: full sweep at (S,P)=(8,8) *)
CheckAllIrreps[8, 8]
```
**Output:**
```
Print: # irreps:[formula]   all solve:[formula]   total time (s):[formula] | Out[101]= \begin{array}{|c|c|c|c|c|}\hline\text{mu} & \text{dim} & \text{deg Q} & \text{time (s)} & \text{solves?}\\\hline16,16,0,0 & 8721 & 0,1,18,19 & \text{0.0028} & \mathrm{True}\\\hline16,15,1,0 & 27455 & 0,2,17,19 & \text{0.0030} & \mathrm{True}\\\hline16,14,2,0 & 47424 & 0,3,16,19 & \text{0.0031} & \mathrm{True}\\\… [1324 chars total]
```

## 2026-07-21 08:02:37
**Input:**
```wolfram
CheckAllIrreps[3, 2]
```
**Output:**
```
Print:
*STR*# irreps: \034BoxData["6"]\034*STR*   all solve: \034BoxData["True"]\034*STR*   total time (s): \034BoxData["0.014884999999999999`"]
Out= Grid[{{"mu", "dim", "deg Q", "time (s)", "solves?"}, {Row[{5, 5, 0, 0}, ","], 196, Row[{0, 1, 7, 8}, ","], NumberForm[0.003575, {6, 4}], True}, {Row[{5, 4, 1, 0}, ","], 384, Row[{0, 2, 6, 8}, ","], NumberForm[0.002121, {6, 4}], True}, {Row[{5, 3, 2, 0}, ","], 300, Row[{0, 3, 5, 8}, ","], NumberForm[0.002158, {6, 4}], True}, {Row[{4, 4, 1, 1}, ","], 50, Row[{1, 2, 6, 7}, ","], NumberForm[0.00209, {6, 4}], True}, {Row[{4, 3, 2, 1}, ","], 64, Row[{1, 3, 5, 7}, ","], NumberForm[0.001937, {6, 4}], True}, {Row[{3, 3, 2, 2}, ","], 6, Row[{2, 3, 5, 6}, ","], NumberForm[0.003004, {6, 4}], True}}, Frame -> All, Alignment -> Left, Background -> {None, {RGBColor[0.88, 0.95, 0.88], {GrayLevel[1], RGBColor[0.96, 0.99, 0.96]}}}]
```
## 2026-07-21 08:02:48 — ▶️ RUN CELLS 3–13
Cell 3: ✓ — (no output)
Cell 4: ✓ — Out[107]= \begin{pmatrix}5 & 5 & 0 & 0\\5 & 4 & 1 & 0\\5 & 3 & 2 & 0\\4 & 4 & 1 & 1\\4 & 3 & 2 & 1\\3 & 3 & 2 & 2\end{pmatrix}
Cell 5: ✓ — (no output)
Cell 6: ✓ — (no output)
Cell 7: ✓ — (no output)
Cell 8: ✓ — Out[120]= \left\{\left\{1\right\},\left\{0,-1,1\right\},\left\{0,0,0,1\right\},\left\{0,0,0,0,-\frac{5}{2},1\right\}\right\}
Cell 9: ✓ — (no output)
Cell 10: ✓ — Out[123]= \left\{1,z^2-z,z^3,z^5-\frac{5\,z^4}{2}\right\}
Cell 11: ✓ — Out[126]= \begin{array}{|c|c|c|}\hline\left\{1,1\right\} & 3 & \mathrm{True}\\\hline\left\{1,2\right\} & 3 & \mathrm{True}\\\hline\left\{1,3\right\} & 3 & \mathrm{True}\\\hline\left\{2,1\right\} & 3 & \mathrm{True}\\\hline\left\{2,2\right\} & 6 & \mathrm{True}\\\hline\left\{2,3\right\} & 6 & \mathrm{True}\\\hline\left\{3,1\right\} & 3 & \mathrm{True}\\\hline\left\{3,2\right\} & 6 & \mathrm{True}\\… [464 chars total]
Cell 12: ✓ — Print: # irreps:[formula]   all solve:[formula]   total time (s):[formula] | Out[87]= \begin{array}{|c|c|c|c|c|}\hline\text{mu} & \text{dim} & \text{deg Q} & \text{time (s)} & \text{solves?}\\\hline16,16,0,0 & 8721 & 0,1,18,19 & \text{0.0030} & \mathrm{True}\\\hline16,15,1,0 & 27455 & 0,2,17,19 & \text{0.0033} & \mathrm{True}\\\hline16,14,2,0 & 47424 & 0,3,16,19 & \text{0.0031} & \mathrm{True}\\\h… [1323 chars total]
Cell 13: ✓ — Print: # irreps:[formula]   all solve:[formula]   total time (s):[formula] | Out[129]= \begin{array}{|c|c|c|c|c|}\hline\text{mu} & \text{dim} & \text{deg Q} & \text{time (s)} & \text{solves?}\\\hline16,16,0,0 & 8721 & 0,1,18,19 & \text{0.0029} & \mathrm{True}\\\hline16,15,1,0 & 27455 & 0,2,17,19 & \text{0.0029} & \mathrm{True}\\\hline16,14,2,0 & 47424 & 0,3,16,19 & \text{0.0030} & \mathrm{True}\\\… [1324 chars total]

## 2026-07-21 08:03:21 — ▶️ RUN CELL 12 (0.07 s)
**In [12]:** `(* --- CheckAllIrreps[S,P]: Gaudin ODE solver over the FULL [S,S,0,0]x[P,P,0,0] decomposition,↵   per-irrep timing. Works directly off GaudinSolve's coefficient lists (BuildQs is just a display↵   con`
**Out:** Print: # irreps:[formula]   all solve:[formula]   total time (s):[formula] | Out[87]= \begin{array}{|c|c|c|c|c|}\hline\text{mu} & \text{dim} & \text{deg Q} & \text{time (s)} & \text{solves?}\\\hline16,16,0,0 & 8721 & 0,1,18,19 & \text{0.0030} & \mathrm{True}\\\hline16,15,1,0 & 27455 & 0,2,17,19 & \t


## 2026-07-21 08:17:35
**Input:**
```wolfram
?dimGL4
```
**Output:**
```
Out= InformationData[<|"ObjectType" -> "Symbol", "Usage" -> "Global`dimGL4", "Documentation" -> None, "OwnValues" -> None, "UpValues" -> None, "DownValues" -> Information`InformationValueForm[DownValues, dimGL4, {dimGL4[lam_] :> Product[(lam[[i]] - lam[[j]] + j - i)/(j - i), {i, 4}, {j, i + 1, 4}]}], "SubValues" -> None, "DefaultValues" -> None, "NValues" -> None, "FormatValues" -> None, "Options" -> None, "Attributes" -> {}, "FullName" -> "Global`dimGL4"|>, False]
```

## 2026-07-21 08:17:40
**Input:**
```wolfram
{"dimGL4[{8,8,0,0}]" -> dimGL4[{8, 8, 0, 0}],
 "dimGL4[{16,16,0,0}]" -> dimGL4[{16, 16, 0, 0}],
 "Total[dimGL4/@irreps[8,8]]" -> Total[dimGL4[#] & /@ irreps[8, 8]],
 "dimGL4[{8,8,0,0}]^2" -> dimGL4[{8, 8, 0, 0}]^2}
```
**Output:**
```
Out= {"dimGL4[{8,8,0,0}]" -> 825, "dimGL4[{16,16,0,0}]" -> 8721, "Total[dimGL4/@irreps[8,8]]" -> 680625, "dimGL4[{8,8,0,0}]^2" -> 680625}
```

## 2026-07-21 08:23:23
**Input:**
```wolfram
AbsoluteTiming[irr88test = irreps[8, 8];]
```
**Output:**
```
Out= {1.56199, Null}
```
## 2026-07-21 08:23:35 — ✏️ EDIT CELL 12
(* --- CheckAllIrreps[S,P,irrList]: Gaudin ODE solver over a precomputed list of irreps
   (e.g. irreps[S,P]), per-irrep timing. Precomputing irrList once and reusing it across calls
   avoids redoing…


## 2026-07-21 08:23:35
**Input:**
```wolfram
(* --- CheckAllIrreps[S,P,irrList]: Gaudin ODE solver over a precomputed list of irreps
   (e.g. irreps[S,P]), per-irrep timing. Precomputing irrList once and reusing it across calls
   avoids redoing the expensive Schur-polynomial tensor decomposition every time -- irreps[8,8]
   alone takes ~1.6s, versus ~0.13s to actually solve all 45 irreps. Works directly off
   GaudinSolve's coefficient lists (BuildQs is just a display convenience, not used here) -- q
   rebuilt as a Sum (SetDelayed), which avoids the degenerate Dot/derivative edge case for deg=0
   solutions. Argument convention: (S,P,irrList). *)
ClearAll[CheckAllIrreps];
CheckAllIrreps[SS_, PP_, irrList_] := Module[{tbl},
   tbl = Table[
      Module[{mu = irrList[[i]], tim, qs, ds, ode, ok, qf},
        ds = degQ[SS, PP, mu];
        {tim, qs} = AbsoluteTiming[GaudinSolve[SS, PP, mu]];
        ode = gaudinODE /. {C2 -> Cas2[mu], C4 -> Cas4[mu], S -> SS, P -> PP} /. wsub;
        ok = AllTrue[Range[4], Module[{m = ds[[#]], cf = qs[[#]]},
           qf[x_] := Sum[cf[[j + 1]] x^j, {j, 0, m}];
           Expand[ode /. Q -> qf] === 0] &];
        {Row[mu, ","], dimGL4[mu], Row[ds, ","], NumberForm[tim, {6, 4}], ok}],
      {i, Length[irrList]}];
   Print["# irreps: ", Length[irrList], "   all solve: ", AllTrue[tbl, Last],
      "   total time (s): ", Total[tbl[[All, 4]] /. NumberForm[x_, _] :> x]];
   Grid[Prepend[tbl, {"mu", "dim", "deg Q", "time (s)", "solves?"}],
     Frame -> All, Alignment -> Left,
     Background -> {None, {RGBColor[.88, .95, .88], {GrayLevel[1], RGBColor[.96, .99, .96]}}}]];
```
**Output:**
```
Out= Null
```
## 2026-07-21 08:23:40 — ✏️ EDIT CELL 13
(* example: precompute the decomposition once, reuse it (e.g. for repeated checks
   or across sessions without re-running the expensive Schur decomposition) *)
irr88 = irreps[8, 8];
CheckAllIrreps[8,…


## 2026-07-21 08:23:42
**Input:**
```wolfram
(* example: precompute the decomposition once, reuse it (e.g. for repeated checks
   or across sessions without re-running the expensive Schur decomposition) *)
irr88 = irreps[8, 8];
CheckAllIrreps[8, 8, irr88]
```
**Output:**
```
Out= Grid[{{"mu", "dim", "deg Q", "time (s)", "solves?"}, {Row[{16, 16, 0, 0}, ","], 8721, Row[{0, 1, 18, 19}, ","], NumberForm[0.00288, {6, 4}], True}, {Row[{16, 15, 1, 0}, ","], 27455, Row[{0, 2, 17, 19}, ","], NumberForm[0.002865, {6, 4}], True}, {Row[{16, 14, 2, 0}, ","], 47424, Row[{0, 3, 16, 19}, ","], NumberForm[0.002954, {6, 4}], True}, {Row[{16, 13, 3, 0}, ","], 62700, Row[{0, 4, 15, 19}, ","], NumberForm[0.002921, {6, 4}], True}, {Row[{16, 12, 4, 0}, ","], 69825, Row[{0, 5, 14, 19}, ","], NumberForm[0.002921, {6, 4}], True}, {Row[{16, 11, 5, 0}, ","], 67431, Row[{0, 6, 13, 19}, ","], NumberForm[0.003001, {6, 4}], True}, {Row[{16, 10, 6, 0}, ","], 55860, Row[{0, 7, 12, 19}, ","], NumberForm[0.003095, {6, 4}], True}, {Row[{16, 9, 7, 0}, ","], 36784, Row[{0, 8, 11, 19}, ","], NumberForm[0.002762, {6, 4}], True}, {Row[{16, 8, 8, 0}, ","], 12825, Row[{0, 9, 10, 19}, ","], NumberForm[0.002712, {6, 4}], True}, {Row[{15, 15, 1, 1}, ","], 5440, Row[{1, 2, 17, 18}, ","], NumberForm[0.002684, {6, 4}], True}, {Row[{15, 14, 2, 1}, ","], 16575, Row[{1, 3, 16, 18}, ","], NumberForm[0.002796, {6, 4}], True}, {Row[{15, 13, 3, 1}, ","], 27489, Row[{1, 4, 15, 18}, ","], NumberForm[0.002782, {6, 4}], True}, {Row[{15, 12, 4, 1}, ","], 34476, Row[{1, 5, 14, 18}, ","], NumberForm[0.002796, {6, 4}], True}, {Row[{15, 11, 5, 1}, ","], 35700, Row[{1, 6, 13, 18}, ","], NumberForm[0.002895, {6, 4}], True}, {Row[{15, 10, 6, 1}, ","], 30855, Row[{1, 7, 12, 18}, ","], NumberForm[0.002994, {6, 4}], True}, {Row[{15, 9, 7, 1}, ","], 20825, Row[{1, 8, 11, 18}, ","], NumberForm[0.003055, {6, 4}], True}, {Row[{15, 8, 8, 1}, ","], 7344, Row[{1, 9, 10, 18}, ","], NumberForm[0.003057, {6, 4}], True}, {Row[{14, 14, 2, 2}, ","], 3185, Row[{2, 3, 16, 17}, ","], NumberForm[0.002833, {6, 4}], True}, {Row[{14, 13, 3, 2}, ","], 9295, Row[{2, 4, 15, 17}, ","], NumberForm[0.002781, {6, 4}], True}, {Row[{14, 12, 4, 2}, ","], 14580, Row[{2, 5, 14, 17}, ","], NumberForm[0.002965, {6, 4}], True}, {Row[{14, 11, 5, 2}, ","], 16940, Row[{2, 6, 13, 17}, ","], NumberForm[0.002897, {6, 4}], True}, {Row[{14, 10, 6, 2}, ","], 15625, Row[{2, 7, 12, 17}, ","], NumberForm[0.002956, {6, 4}], True}, {Row[{14, 9, 7, 2}, ","], 10935, Row[{2, 8, 11, 17}, ","], NumberForm[0.002949, {6, 4}], True}, {Row[{14, 8, 8, 2}, ","], 3920, Row[{2, 9, 10, 17}, ","], NumberForm[0.003115, {6, 4}], True}, {Row[{13, 13, 3, 3}, ","], 1716, Row[{3, 4, 15, 16}, ","], NumberForm[0.002765, {6, 4}], True}, {Row[{13, 12, 4, 3}, ","], 4719, Row[{3, 5, 14, 16}, ","], NumberForm[0.002898, {6, 4}], True}, {Row[{13, 11, 5, 3}, ","], 6825, Row[{3, 6, 13, 16}, ","], NumberForm[0.002972, {6, 4}], True}, {Row[{13, 10, 6, 3}, ","], 7020, Row[{3, 7, 12, 16}, ","], NumberForm[0.002879, {6, 4}], True}, {Row[{13, 9, 7, 3}, ","], 5200, Row[{3, 8, 11, 16}, ","], NumberForm[0.002985, {6, 4}], True}, {Row[{13, 8, 8, 3}, ","], 1911, Row[{3, 9, 10, 16}, ","], NumberForm[0.003008, {6, 4}], True}, {Row[{12, 12, 4, 4}, ","], 825, Row[{4, 5, 14, 15}, ","], NumberForm[0.002811, {6, 4}], True}, {Row[{12, 11, 5, 4}, ","], 2079, Row[{4, 6, 13, 15}, ","], NumberForm[0.002972, {6, 4}], True}, {Row[{12, 10, 6, 4}, ","], 2640, Row[{4, 7, 12, 15}, ","], NumberForm[0.002835, {6, 4}], True}, {Row[{12, 9, 7, 4}, ","], 2156, Row[{4, 8, 11, 15}, ","], NumberForm[0.002962, {6, 4}], True}, {Row[{12, 8, 8, 4}, ","], 825, Row[{4, 9, 10, 15}, ","], NumberForm[0.002903, {6, 4}], True}, {Row[{11, 11, 5, 5}, ","], 336, Row[{5, 6, 13, 14}, ","], NumberForm[0.002931, {6, 4}], True}, {Row[{11, 10, 6, 5}, ","], 735, Row[{5, 7, 12, 14}, ","], NumberForm[0.002974, {6, 4}], True}, {Row[{11, 9, 7, 5}, ","], 729, Row[{5, 8, 11, 14}, ","], NumberForm[0.002945, {6, 4}], True}, {Row[{11, 8, 8, 5}, ","], 300, Row[{5, 9, 10, 14}, ","], NumberForm[0.003057, {6, 4}], True}, {Row[{10, 10, 6, 6}, ","], 105, Row[{6, 7, 12, 13}, ","], NumberForm[0.002971, {6, 4}], True}, {Row[{10, 9, 7, 6}, ","], 175, Row[{6, 8, 11, 13}, ","], NumberForm[0.002814, {6, 4}], True}, {Row[{10, 8, 8, 6}, ","], 84, Row[{6, 9, 10, 13}, ","], NumberForm[0.003182, {6, 4}], True}, {Row[{9, 9, 7, 7},
[output truncated — 4481 chars total]
```
## 2026-07-21 08:23:52 — ▶️ RUN CELLS 3–13
Cell 3: ✓ — (no output)
Cell 4: ✓ — Out[137]= \begin{pmatrix}5 & 5 & 0 & 0\\5 & 4 & 1 & 0\\5 & 3 & 2 & 0\\4 & 4 & 1 & 1\\4 & 3 & 2 & 1\\3 & 3 & 2 & 2\end{pmatrix}
Cell 5: ✓ — (no output)
Cell 6: ✓ — (no output)
Cell 7: ✓ — (no output)
Cell 8: ✓ — Out[150]= \left\{\left\{1\right\},\left\{0,-1,1\right\},\left\{0,0,0,1\right\},\left\{0,0,0,0,-\frac{5}{2},1\right\}\right\}
Cell 9: ✓ — (no output)
Cell 10: ✓ — Out[153]= \left\{1,z^2-z,z^3,z^5-\frac{5\,z^4}{2}\right\}
Cell 11: ✓ — Out[156]= \begin{array}{|c|c|c|}\hline\left\{1,1\right\} & 3 & \mathrm{True}\\\hline\left\{1,2\right\} & 3 & \mathrm{True}\\\hline\left\{1,3\right\} & 3 & \mathrm{True}\\\hline\left\{2,1\right\} & 3 & \mathrm{True}\\\hline\left\{2,2\right\} & 6 & \mathrm{True}\\\hline\left\{2,3\right\} & 6 & \mathrm{True}\\\hline\left\{3,1\right\} & 3 & \mathrm{True}\\\hline\left\{3,2\right\} & 6 & \mathrm{True}\\… [464 chars total]
Cell 12: ✓ — Print: # irreps:[formula]   all solve:[formula]   total time (s):[formula] | Out[87]= \begin{array}{|c|c|c|c|c|}\hline\text{mu} & \text{dim} & \text{deg Q} & \text{time (s)} & \text{solves?}\\\hline16,16,0,0 & 8721 & 0,1,18,19 & \text{0.0030} & \mathrm{True}\\\hline16,15,1,0 & 27455 & 0,2,17,19 & \text{0.0033} & \mathrm{True}\\\hline16,14,2,0 & 47424 & 0,3,16,19 & \text{0.0031} & \mathrm{True}\\\h… [1323 chars total]
Cell 13: ✓ — Print: # irreps:[formula]   all solve:[formula]   total time (s):[formula] | Out[160]= \begin{array}{|c|c|c|c|c|}\hline\text{mu} & \text{dim} & \text{deg Q} & \text{time (s)} & \text{solves?}\\\hline16,16,0,0 & 8721 & 0,1,18,19 & \text{0.0029} & \mathrm{True}\\\hline16,15,1,0 & 27455 & 0,2,17,19 & \text{0.0047} & \mathrm{True}\\\hline16,14,2,0 & 47424 & 0,3,16,19 & \text{0.0037} & \mathrm{True}\\\… [1324 chars total]

