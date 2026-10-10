## Deleted Cell 43 — 2026-09-09 10:40:30
**Kind:** code
```wolfram
(* Order eps^0, done properly. At eps=0 the four twists collide and the four Baxter equations
   become ONE equation, so the four q_k must share the same leading-order solution. Imposing
   q_1^(0) = ... = q_4^(0) removes exactly the 6-parameter unipotent rotation found above and
   makes order 0 UNIQUE. This is not an arbitrary gauge choice: if it were wrong, the order-1
   linear system below would be inconsistent. *)
sol0fix = Join[First[sol0] /. Thread[free0 -> 0], Thread[free0 -> 0]];
lead0 = Table[pv[v, 0] /. sol0fix, {v, pert["unk"]}];
{"order-0 tc" -> Take[lead0, 5], "order-0 q coefficients" -> Drop[lead0, 5],
 "order-0 residual" -> Union[Expand[Flatten[pert["eqs"][[1]]] /. sol0fix]]}
```

## Deleted Cell 43 — 2026-09-09 10:40:42
**Kind:** code
```wolfram
(* Independent check of the order-1 verdict, using the EXACT (non-perturbative) solver.
   Solve the pseudovacuum with eps fully symbolic (multiplicity 1, so this is rational and exact),
   renormalise each monic q_k to q_k(u0)=1, and take eps -> 0 coefficient by coefficient
   (never on the whole expression -- see the Section 6 trap). The perturbative claim is that
   all four limits coincide, i.e. the order-0 rotation vanishes. *)
Module[{rec, qn, x},
  rec = First@TwistedSolve[1, 1, {2, 2, 0, 0}, eps1];
  qn = Table[Module[{c = rec["q"][[k]]}, Together[c/qpoly[c, u0]]], {k, 4}];
  Grid[Prepend[
    Table[{k, rec["M"][[k]], Simplify[Limit[#, eps -> 0]] & /@ qn[[k]]}, {k, 4}],
    {"k", "M_k", "coefficients of q_k(u0)=1 at eps -> 0"}], Frame -> All, Alignment -> Left]]
```

## Deleted Cell 52 — 2026-09-09 10:54:45
**Kind:** code
```wolfram
(* Twisted large-u asymptotics. For each branch, q_k(u) = u^M_k (1 + Sum_i r[k,i] u^-i); the
   twist z[j] is left fully symbolic (the characteristic identity holds automatically, checked
   in Section 1). Collect the large-u equations and see which tc they determine WITHOUT ever
   solving the Baxter equation itself. *)
ClearAll[asymTw, rr];
asymTw[SS_, PP_, μ_List, K_] := Module[{x, Ms, eqs = {}, qf, R, ser},
   Ms = Mdeg[SS, PP, μ];
   Do[qf[y_] := y^Ms[[k]] (1 + Sum[rr[k, i] y^-i, {i, K}]);
      R = BaxRes[SS, PP, μ][k][qf][x];
      ser = Series[Together[R]/x^(Ms[[k]] + 8), {x, Infinity, K + 1}];
      eqs = Join[eqs, DeleteCases[Simplify[CoefficientList[Normal[ser] /. x -> 1/w, w]], 0]],
      {k, 4}];
   eqs];

(* Concrete but generic state: (S,P)=(3,1), mu={3,3,1,1}, M={1,1,3,3}. theta1, theta2, h and
   all four z[j] stay symbolic. *)
asymEqs = asymTw[3, 1, {3, 3, 1, 1}, 3];
asymUnk = Join[tcvars, Flatten[Table[rr[k, i], {k, 4}, {i, 3}]]];
asymSolTw = Solve[Thread[asymEqs == 0], asymUnk];
{"#equations" -> Length[asymEqs], "#solutions" -> Length[asymSolTw],
 "tc determined" -> If[asymSolTw === {}, $Failed, Select[tcvars, MemberQ[First[asymSolTw][[All, 1]], #] &]],
 "tc still free" -> If[asymSolTw === {}, $Failed, Complement[tcvars, First[asymSolTw][[All, 1]]]]}
```

## Deleted Cell 51 — 2026-09-09 10:54:45
**Kind:** markdown
```markdown
### 8.3 How many transfer-matrix unknowns are there really?

The ansatz of Section 2 leaves five: `tc[1]` (constant of `tau1`), `tc[2],tc[3],tc[4]` (the `u^2,u^1,u^0` of `tau2`) and `tc[5]` (constant in `tau3`'s quadratic). Untwisted, the asymptotic solve of `baxter_SU4_L2_XXX.wb` fixes the `u^1` coefficient of `tau2` as well, leaving only the two Casimir-carrying coefficients. The question here is whether that extra relation survives the twist, i.e. whether the true count is **1 + 2 + 1 = 4**.

Test: push the large-`u` expansion of the twisted Baxter equation past the order that gave the subleading formula, with the subleading `Q`-data `r[k,i]` carried explicitly, and see which `tc` drop out.
```

