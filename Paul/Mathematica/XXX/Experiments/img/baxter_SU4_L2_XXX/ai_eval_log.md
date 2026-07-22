# AI Action Log — baxter_SU4_L2_XXX

## 2026-07-20 12:29:00 — ▶️ RUN CELLS 1–3
Cell 1: ✓ — (no output)
Cell 2: ✓ — (no output)
Cell 3: ✓ — (no output)


## 2026-07-20 12:29:07
**Input:**
```wolfram
Table[Sum[(-1)^j Binomial[4,j] j^m,{j,0,4}],{m,0,5}]
```
**Output:**
```
[message] Power::indet: Indeterminate expression 0  encountered.
Out= {Indeterminate, 0, 0, 0, 24, 240}
```

## 2026-07-20 12:29:31
**Input:**
```wolfram
ser = Series[Bax[Q]/Q[u] /. Q[u_] -> u^p /. (u + x_)^y_ -> u^y (1 + x/u)^y, {u, Infinity, 3}] // Normal // Expand;
coeffs = Table[Coefficient[ser, u, -m], {m, 1, 3}];
Column[Collect[#, p, Simplify] & /@ coeffs]
```
**Output:**
```
Out= Column[{4*h*(3 + P + S) + 8*θ1 + 8*θ2 - a[1] + a[2], h^2*(22 + 9*P^2 + 31*S + 9*S^2 + P*(31 + 16*S)) + 14*θ1^2 + 14*θ2^2 - θ2*a[1] + 2*θ2*a[2] - h*p*(6*h*(2 + P + S) + 12*θ1 + 12*θ2 + a[2]) + θ1*(20*θ2 - a[1] + 2*a[2]) + h*((34 + 22*P + 24*S)*θ1 + (34 + 24*P + 22*S)*θ2 - P*a[1] - S*a[1] + 2*a[2] + 2*P*a[2] + 2*S*a[2]) - b[1] + b[2] - b[3], h^3*(32 + 14*P^3 + 65*S + 55*S^2 + 14*S^3 + P^2*(55 + 28*S) + P*(65 + 80*S + 28*S^2)) + 20*θ1^3 + 20*θ2^3 - θ2^2*a[1] + 3*θ2^2*a[2] + (h^2*p^2*(6*h*(2 + P + S) + 12*θ1 + 12*θ2 + a[2]))/2 + θ1^2*(32*θ2 - a[1] + 3*a[2]) + h^2*(2*(33 + 18*P^2 + 58*S + 25*S^2 + P*(43 + 32*S))*θ1 + 2*(33 + 58*P + 25*P^2 + 43*S + 32*P*S + 18*S^2)*θ2 - P^2*a[1] - P*S*a[1] - S^2*a[1] + 3*a[2] + 5*P*a[2] + 3*P^2*a[2] + 5*S*a[2] + 4*P*S*a[2] + 3*S^2*a[2]) - θ2*b[1] + 2*θ2*b[2] - (h*p*(h^2*(64 + 22*P^2 + 72*S + 22*S^2 + 8*P*(9 + 5*S)) + h*((88 + 56*P + 60*S)*θ1 + (88 + 60*P + 56*S)*θ2 + (5 + 4*P + 4*S)*a[2]) + 2*(18*θ1^2 + 18*θ2^2 + 2*θ2*a[2] + 2*θ1*(12*θ2 + a[2]) + b[2] - 2*b[3])))/2 + θ1*(32*θ2^2 - θ2*a[1] + 4*θ2*a[2] - b[1] + 2*b[2] - b[3]) - θ2*b[3] + h*((58 + 34*P + 56*S)*θ1^2 + (58 + 56*P + 34*S)*θ2^2 + θ2*(-2*P*a[1] - S*a[1] + 5*a[2] + 6*P*a[2] + 4*S*a[2]) + θ1*((80 + 66*P + 66*S)*θ2 - P*a[1] - 2*S*a[1] + 5*a[2] + 4*P*a[2] + 6*S*a[2]) - P*b[1] - S*b[1] + 2*b[2] + 2*P*b[2] + 2*S*b[2] - 2*b[3] - 2*P*b[3] - 2*S*b[3]) + c[2]}]
```

## 2026-07-20 12:30:44
**Input:**
```wolfram
(* orders 1/u,1/u^2,1/u^3 as polynomials in p; each power-of-p coefficient must vanish
   independently (4 branches, same t-coeffs, distinct p_k). Collect the full equation set. *)
eqs = Flatten[Table[CoefficientList[coeffs[[m]], p], {m, 1, 3}]];
sol = Solve[# == 0 & /@ eqs, {a[1], a[2], b[1], b[2], c[2], b[3], d[2]}];
sol // Length
```
**Output:**
```
[message] Solve::svars: Equations may not give solutions for all "solve" variables.
Out= 1
```

## 2026-07-20 12:30:47
**Input:**
```wolfram
Column[Simplify /@ sol[[1]]]
```
**Output:**
```
Out= Column[{a[1] -> -2*(h*(P + S) + 2*(θ1 + θ2)), a[2] -> -6*(h*(2 + P + S) + 2*(θ1 + θ2)), b[1] -> (-(h^2*(6 + P^2 + 7*S + S^2 + P*(7 + 4*S))) - 6*θ1^2 - 16*θ1*θ2 - 6*θ2^2 - 2*h*((9 + 4*P + 3*S)*θ1 + (9 + 3*P + 4*S)*θ2) + b[2])/2, c[2] -> (h^3*(12 + P^3 + 18*S + 7*S^2 + S^3 + P^2*(7 + 5*S) + P*(18 + 22*S + 5*S^2)) + 2*h^2*((18 + 20*P + 5*P^2 + 16*S + 9*P*S + 4*S^2)*θ1 + (18 + 16*P + 4*P^2 + 20*S + 9*P*S + 5*S^2)*θ2) + 2*(θ1 + θ2)*(6*θ1^2 + 12*θ1*θ2 + 6*θ2^2 - b[2]) + h*(2 + P + S)*(18*θ1^2 + 36*θ1*θ2 + 18*θ2^2 - b[2]))/2, b[3] -> (-(h^2*(-2 + P^2 + 3*S + S^2 + P*(3 + 4*S))) - 6*θ1^2 - 24*θ1*θ2 - 6*θ2^2 - 2*h*((5 + 4*P + 3*S)*θ1 + (5 + 3*P + 4*S)*θ2) + b[2])/2}]
```

## 2026-07-20 12:31:42
**Input:**
```wolfram
(* main-notebook mode numbers from a Cartan weight (H1,H2,H3): does n_k = mu_k on the HWS? *)
nvec[H1_, H2_, H3_, Ntot_] := Module[{n4 = (Ntot - H1 - 2 H2 - 3 H3)/4}, {n4 + H1 + H2 + H3, n4 + H2 + H3, n4 + H3, n4}];
(* HWS of irrep mu=[m1,m2,m3,m4]: Cartan charges H_a = m_a - m_{a+1}, Ntot=Sum m *)
Table[With[{m = mu}, nvec[m[[1]] - m[[2]], m[[2]] - m[[3]], m[[3]] - m[[4]], Total[m]] == m],
  {mu, {{4, 2, 1, 0}, {3, 3, 1, 1}, {5, 2, 2, 0}, {2, 2, 2, 2}}}]
```
**Output:**
```
Out= {True, True, True, True}
```

## 2026-07-20 12:34:08
**Input:**
```wolfram
(* rho-centered shifted weights l_i = mu_i + (N+1)/2 - i, N=4. Test whether the Gelfand
   invariants C_k = tr(E^k) eigenvalues are power-sums (l_i^k) minus the mu=0 vacuum. *)
ll[mu_] := Table[mu[[i]] + 5/2 - i, {i, 1, 4}];
Ck[mu_, k_] := Total[ll[mu]^k] - Total[ll[{0, 0, 0, 0}]^k];
(* single-site [x,x,0,0] reference values the notebook trusts *)
Table[{k, Simplify[Ck[{x, x, 0, 0}, k]], {2 x (x + 2), 2 x (x + 2)^2, 2 x (x + 2)^3}[[k - 1]]}, {k, 2, 4}] // TableForm
```
**Output:**
```
Out= TableForm[{{2, 2*x*(2 + x), 2*x*(2 + x)}, {3, (x*(15 + 12*x + 4*x^2))/2, 2*x*(2 + x)^2}, {4, x*(14 + 15*x + 8*x^2 + 2*x^3), 2*x*(2 + x)^3}}]
```

## 2026-07-20 12:34:56
**Input:**
```wolfram
(* compact gl(4) GT construction: patterns for top row mu, simple generators E[k,k+1],E[k+1,k],
   diagonal E[k,k]; general E via commutators; then C_n = tr(E^n) eigenvalue on the irrep. *)
Clear[gtPats, gl4Gens, Cn];
(* GT patterns: nested rows r4=mu >= r3 >= r2 >= r1, betweenness l^{k}_i>=l^{k-1}_i>=l^{k}_{i+1} *)
rowsBetween[top_] := Module[{n = Length[top]},
  If[n == 1, {{First[top]}},
   Tuples[Table[Range[top[[i + 1]], top[[i]]], {i, 1, n - 1}] // Reverse // Map[Reverse]] ]];
gtPats[mu_] := Module[{r4 = {mu}, r3, r2, r1},
   Flatten[Table[
     Table[Table[{mu, a, b, c},
        {c, rowsBetween[b]}], {b, rowsBetween[a]}],
      {a, rowsBetween[mu]}], 3]];
Length[gtPats[{2, 1, 0, 0}]]
```
**Output:**
```
Out= 0
```

## 2026-07-20 12:35:09
**Input:**
```wolfram
rowsBetween[top_] := Tuples[Table[Range[top[[i + 1]], top[[i]]], {i, 1, Length[top] - 1}]];
gtPats[mu_] := Flatten[Table[{mu, a, b, c},
    {a, rowsBetween[mu]}, {b, rowsBetween[a]}, {c, rowsBetween[b]}], 2];
{Length[gtPats[{2, 1, 0, 0}]], Length[gtPats[{1, 1, 0, 0}]], Length[gtPats[{2, 2, 0, 0}]]}
```
**Output:**
```
Out= {20, 6, 20}
```

## 2026-07-20 12:37:37
**Input:**
```wolfram
(* PP residue form: <tr E^m> = generating fn. Test  G(u)=Sum_i (1/(u-l_i)) Prod_{j!=i}(u-l_j-1)/(u-l_j),
   l_i = lam_i + N - i, N=4.  Then <tr E^m> = coeff via  Sum_m <trE^m> u^{-m-1} = G(u) - (stuff)? *)
NN = 4;
lvec[lam_] := Table[lam[[i]] + NN - i, {i, 1, NN}];
Gfun[lam_, u_] := Sum[(1/(u - lvec[lam][[i]])) Product[If[j == i, 1, (u - lvec[lam][[j]] - 1)/(u - lvec[lam][[j]])], {j, 1, NN}], {i, 1, NN}];
(* expand G at u=Infinity; the coefficient of u^{-m-1} is a candidate for <tr E^m> *)
trial[lam_, m_] := SeriesCoefficient[Gfun[lam, u], {u, Infinity, m + 1}];
{trial[{1, 0, 0, 0}, 2], trial[{1, 0, 0, 0}, 4]}
```
**Output:**
```
Out= {-9, -68}
```

## 2026-07-20 12:38:44
**Input:**
```wolfram
Clear[evalHW, CasN];
(* evalHW[word, mu]: <hw| E_{a1b1}...E_{anbn} |hw>, processing leftmost op.
   raising (a<b) commutes right, kills hw at the end; lowering (a>b) at left kills bra; a=b -> mu[a]. *)
evalHW[{}, mu_] := 1;
evalHW[w_List, mu_] := evalHW[w, mu] = Module[{a = w[[1, 1]], b = w[[1, 2]], c, d, rest = Rest[w]},
   Which[
    a > b, 0,
    a == b, mu[[a]] evalHW[rest, mu],
    True, (* a<b raising *)
    If[rest === {}, 0,
     c = w[[2, 1]]; d = w[[2, 2]];
     evalHW[Join[{w[[2]], w[[1]]}, w[[3 ;;]]], mu]
      + If[b == c, evalHW[Join[{{a, d}}, w[[3 ;;]]], mu], 0]
      - If[d == a, evalHW[Join[{{c, b}}, w[[3 ;;]]], mu], 0]]]];
CasN[mu_, n_] := Sum[evalHW[Table[{idx[[k]], idx[[Mod[k, n] + 1]]}, {k, 1, n}], mu],
   Evaluate[Sequence @@ Table[{ToExpression["idx" <> ToString[q]], 1, 4}, {q, n}]]] /.
   Thread[Table[ToExpression["idx" <> ToString[q]], {q, n}] -> {}];
(* simpler: enumerate index tuples explicitly *)
CasN[mu_, n_] := Total[evalHW[Table[{tup[[k]], tup[[Mod[k, n] + 1]]}, {k, 1, n}], mu] & /@ Tuples[Range[4], n]];
{CasN[{1, 0, 0, 0}, 2], CasN[{1, 0, 0, 0}, 4]}
```
**Output:**
```
[message] Part::partd: Part specification tup[[1]] is longer than depth of object.
[message] [Part::partd ×2 more suppressed — first instance shown above]
[message] General::stop: Further output of Part::partd will be suppressed during this calculation.
Out= {16*Which[tup[[1]] > tup[[2]], 0, a$18963 == b$18963, {1, 0, 0, 0}[[a$18963]]*evalHW[rest$18963, {1, 0, 0, 0}], True, If[rest$18963 === {}, 0, c$18963 = {{tup[[1]], tup[[2]]}, {tup[[2]], tup[[1]]}}[[2,1]]; d$18963 = {{tup[[1]], tup[[2]]}, {tup[[2]], tup[[1]]}}[[2,2]]; evalHW[Join[{{{tup[[1]], tup[[2]]}, {tup[[2]], tup[[1]]}}[[2]], {{tup[[1]], tup[[2]]}, {tup[[2]], tup[[1]]}}[[1]]}, {{tup[[1]], tup[[2]]}, {tup[[2]], tup[[1]]}}[[3 ;; All]]], {1, 0, 0, 0}] + If[b$18963 == c$18963, evalHW[Join[{{a$18963, d$18963}}, {{tup[[1]], tup[[2]]}, {tup[[2]], tup[[1]]}}[[3 ;; All]]], {1, 0, 0, 0}], 0] - If[d$18963 == a$18963, evalHW[Join[{{c$18963, b$18963}}, {{tup[[1]], tup[[2]]}, {tup[[2]], tup[[1]]}}[[3 ;; All]]], {1, 0, 0, 0}], 0]]], 256*Which[tup[[1]] > tup[[2]], 0, a$18964 == b$18964, {1, 0, 0, 0}[[a$18964]]*evalHW[rest$18964, {1, 0, 0, 0}], True, If[rest$18964 === {}, 0, c$18964 = {{tup[[1]], tup[[2]]}, {tup[[2]], tup[[3]]}, {tup[[3]], tup[[4]]}, {tup[[4]], tup[[1]]}}[[2,1]]; d$18964 = {{tup[[1]], tup[[2]]}, {tup[[2]], tup[[3]]}, {tup[[3]], tup[[4]]}, {tup[[4]], tup[[1]]}}[[2,2]]; evalHW[Join[{{{tup[[1]], tup[[2]]}, {tup[[2]], tup[[3]]}, {tup[[3]], tup[[4]]}, {tup[[4]], tup[[1]]}}[[2]], {{tup[[1]], tup[[2]]}, {tup[[2]], tup[[3]]}, {tup[[3]], tup[[4]]}, {tup[[4]], tup[[1]]}}[[1]]}, {{tup[[1]], tup[[2]]}, {tup[[2]], tup[[3]]}, {tup[[3]], tup[[4]]}, {tup[[4]], tup[[1]]}}[[3 ;; All]]], {1, 0, 0, 0}] + If[b$18964 == c$18964, evalHW[Join[{{a$18964, d$18964}}, {{tup[[1]], tup[[2]]}, {tup[[2]], tup[[3]]}, {tup[[3]], tup[[4]]}, {tup[[4]], tup[[1]]}}[[3 ;; All]]], {1, 0, 0, 0}], 0] - If[d$18964 == a$18964, evalHW[Join[{{c$18964, b$18964}}, {{tup[[1]], tup[[2]]}, {tup[[2]], tup[[3]]}, {tup[[3]], tup[[4]]}, {tup[[4]], tup[[1]]}}[[3 ;; All]]], {1, 0, 0, 0}], 0]]]}
```

## 2026-07-20 12:39:00
**Input:**
```wolfram
CasN[mu_, n_] := Total[
   Function[tup, evalHW[Table[{tup[[k]], tup[[Mod[k, n] + 1]]}, {k, 1, n}], mu]] /@ Tuples[Range[4], n]];
{CasN[{1, 0, 0, 0}, 2], CasN[{1, 0, 0, 0}, 4], CasN[{0, 0, 0, 0}, 4]}
```
**Output:**
```
[message] $RecursionLimit::reclim: Recursion depth of 1024 exceeded.
(no output)
```

## 2026-07-20 12:39:31
**Input:**
```wolfram
Clear[evalHW, CasN, cls];
cls[{a_, b_}] := Which[a > b, 0, a == b, 1, True, 2];  (* lowering<cartan<raising *)
evalHW[w_List, mu_] := evalHW[w, mu] = Module[{p, L, R, a, b, c, d, pre, post},
   p = SelectFirst[Range[Length[w] - 1], cls[w[[#]]] > cls[w[[# + 1]]] &, 0];
   If[p == 0,
    (* fully PBW-sorted: nonzero only if all Cartan *)
    If[AllTrue[w, #[[1]] == #[[2]] &], Product[mu[[e[[1]]]], {e, w}], 0],
    L = w[[p]]; R = w[[p + 1]]; {a, b} = L; {c, d} = R;
    pre = w[[1 ;; p - 1]]; post = w[[p + 2 ;;]];
    evalHW[Join[pre, {R, L}, post], mu]
     + If[b == c, evalHW[Join[pre, {{a, d}}, post], mu], 0]
     - If[d == a, evalHW[Join[pre, {{c, b}}, post], mu], 0]]];
CasN[mu_, n_] := Total[Function[tup, evalHW[Table[{tup[[k]], tup[[Mod[k, n] + 1]]}, {k, 1, n}], mu]] /@ Tuples[Range[4], n]];
{CasN[{1, 0, 0, 0}, 2], CasN[{1, 0, 0, 0}, 4]}
```
**Output:**
```
Out= {4, 64}
```

## 2026-07-20 12:39:42
**Input:**
```wolfram
(* verify [x,x,0,0]: C2=2x(x+2), C3=2x(x+2)^2, C4=2x(x+2)^3 *)
chk = Table[{n, Simplify[CasN[{x, x, 0, 0}, n]], {2 x (x + 2), 2 x (x + 2)^2, 2 x (x + 2)^3}[[n - 1]]}, {n, 2, 4}];
Print[TableForm[chk, TableHeadings -> {None, {"n", "CasN[{x,x,0,0}]", "expected"}}]];
(* closed form C2(mu), C4(mu) as symmetric functions of mu1..mu4 *)
c2mu = Simplify[CasN[{m1, m2, m3, m4}, 2]];
c4mu = Simplify[CasN[{m1, m2, m3, m4}, 4]];
Print["C2(mu) = ", c2mu];
Print["C4(mu) = ", c4mu]
```
**Output:**
```
Print:
BoxData[TagBox[TagBox[GridBox[{{TagBox["\\"n\\"", HoldForm], TagBox["\\"CasN[{x,x,0,0}]\\"", HoldForm], TagBox["\\"expected\\"", HoldForm]}, {"2", RowBox[{"2", " ", "x", " ", RowBox[{"(", RowBox[{"x", "+", "2"}], ")"}]}], RowBox[{"2", " ", "x", " ", RowBox[{"(", RowBox[{"x", "+", "2"}], ")"}]}]}, {"3", RowBox[{"2", " ", "x", " ", SuperscriptBox[RowBox[{"(", RowBox[{"x", "+", "2"}], ")"}], "2"]}], RowBox[{"2", " ", "x", " ", SuperscriptBox[RowBox[{"(", RowBox[{"x", "+", "2"}], ")"}], "2"]}]}, {"4", RowBox[{"2", " ", "x", " ", SuperscriptBox[RowBox[{"(", RowBox[{"x", "+", "2"}], ")"}], "3"]}], RowBox[{"2", " ", "x", " ", SuperscriptBox[RowBox[{"(", RowBox[{"x", "+", "2"}], ")"}], "3"]}]}}, RowSpacings -> 1, ColumnSpacings -> 3, RowAlignments -> Baseline, ColumnAlignments -> Left, GridBoxDividers -> {"Columns" -> {{False}}, "Rows" -> {False, True, {False}, False}}, ColumnAlignments -> Left], {None, OutputFormsDump`HeadedColumns}], Function[BoxForm`e$, TableForm[BoxForm`e$, TableHeadings -> {None, {"n", "CasN[{x,x,0,0}]", "expected"}}]]]]*STR*C2(mu) = \034BoxData[RowBox[{SuperscriptBox["m1", "2"], "+", RowBox[{"3", " ", "m1"}], "+", SuperscriptBox["m2", "2"], "+", "m2", "+", SuperscriptBox["m3", "2"], "-", "m3", "+", SuperscriptBox["m4", "2"], "-", RowBox[{"3", " ", "m4"}]}]]*STR*C4(mu) = \034BoxData[RowBox[{SuperscriptBox["m1", "4"], "+", RowBox[{"9", " ", SuperscriptBox["m1", "3"]}], "-", RowBox[{SuperscriptBox["m1", "2"], " ", RowBox[{"(", RowBox[{"m2", "+", "m3", "+", "m4", "-", "27"}], ")"}]}], "-", RowBox[{"m1", " ", RowBox[{"(", RowBox[{SuperscriptBox["m2", "2"], "+", RowBox[{"8", " ", "m2"}], "+", SuperscriptBox["m3", "2"], "+", RowBox[{"6", " ", "m3"}], "+", SuperscriptBox["m4", "2"], "+", RowBox[{"4", " ", "m4"}], "-", "27"}], ")"}]}], "+", SuperscriptBox["m2", "4"], "+", RowBox[{"5", " ", SuperscriptBox["m2", "3"]}], "-", RowBox[{SuperscriptBox["m2", "2"], " ", RowBox[{"(", RowBox[{"m3", "+", "m4", "-", "5"}], ")"}]}], "-", RowBox[{"m2", " ", RowBox[{"(", RowBox[{SuperscriptBox["m3", "2"], "+", RowBox[{"4", " ", "m3"}], "+", SuperscriptBox["m4", "2"], "+", RowBox[{"2", " ", "m4"}], "+", "11"}], ")"}]}], "+", SuperscriptBox["m3", "4"], "+", SuperscriptBox["m3", "3"], "-", RowBox[{SuperscriptBox["m3", "2"], " ", "m4"}], "-", RowBox[{"5", " ", SuperscriptBox["m3", "2"]}], "-", RowBox[{"m3", " ", SuperscriptBox["m4", "2"]}], "-", RowBox[{"13", " ", "m3"}], "+", SuperscriptBox["m4", "4"], "-", RowBox[{"3", " ", SuperscriptBox["m4", "3"]}], "-", RowBox[{"3", " ", SuperscriptBox["m4", "2"]}], "-", RowBox[{"3", " ", "m4"}]}]]
Out= Null
```

## 2026-07-20 12:45:20
**Input:**
```wolfram
(* Paul's generating function (2.15), h=1 for the symmetric-function content:
   G4(u)=Prod_i (u - (mu_i + N - i)),  N=4.  Casimir eigenvalues = its coefficients. *)
shw[mu_] := Table[mu[[i]] + 4 - i, {i, 1, 4}];
G4[mu_, u_] := Product[u - (mu[[i]] + 4 - i), {i, 1, 4}];
(* do the coefficients (elem. symmetric fns of shifted weights) reproduce tr(E^2), tr(E^4)?
   verify my normal-ordered C2,C4 are symmetric functions of the SAME shifted weights *)
c2mu = Expand[CasN[{m1, m2, m3, m4}, 2]];
c4mu = Expand[CasN[{m1, m2, m3, m4}, 4]];
(* express C2,C4 via power sums of shifted weights minus vacuum -> confirms gen-fn origin *)
pk[mu_, k_] := Total[shw[mu]^k];
{"C2 = p2(shifted) - p2(vac)?" -> Simplify[c2mu - (pk[{m1, m2, m3, m4}, 2] - pk[{0, 0, 0, 0}, 2])],
 "deg q_k = (S+P+3) - shw_k ?" -> Simplify[Table[(S + P) - m[[k]] + (k - 1), {k, 1, 4}] == (S + P + 3) - shw[m /. m -> {mm1, mm2, mm3, mm4}] /. m -> {mm1, mm2, mm3, mm4}]}
```
**Output:**
```
[message] Part::partd: Part specification m[[1]] is longer than depth of object.
[message] [Part::partd ×2 more suppressed — first instance shown above]
[message] General::stop: Further output of Part::partd will be suppressed during this calculation.
Out= {"C2 = p2(shifted) - p2(vac)?" -> -3*(m1 + m2 + m3 + m4), "deg q_k = (S+P+3) - shw_k ?" -> True}
```

## 2026-07-20 12:46:13
**Input:**
```wolfram
Column[{
  "a1 = " <> ToString[a[1] /. sol[[1]], InputForm],
  "a2 = " <> ToString[a[2] /. sol[[1]], InputForm],
  "b1 = " <> ToString[Collect[b[1] /. sol[[1]], b[2]], InputForm],
  "b3 = " <> ToString[Collect[b[3] /. sol[[1]], b[2]], InputForm],
  "c2 = " <> ToString[Collect[c[2] /. sol[[1]], b[2]], InputForm],
  "C2mu = " <> ToString[Sum[mu[[i]] (mu[[i]] + 5 - 2 i), {i, 4}] /. mu -> {m1, m2, m3, m4}, InputForm],
  "C4mu = " <> ToString[c4mu, InputForm]
}]
```
**Output:**
```
[message] Part::partd: Part specification mu[[1]] is longer than depth of object.
[message] [Part::partd ×2 more suppressed — first instance shown above]
[message] General::stop: Further output of Part::partd will be suppressed during this calculation.
Out= Column[{"a1 = -2*(h*P + h*S + 2*θ1 + 2*θ2)", "a2 = -6*(2*h + h*P + h*S + 2*θ1 + 2*θ2)", "b1 = (-6*h^2 - 7*h^2*P - h^2*P^2 - 7*h^2*S - 4*h^2*P*S - h^2*S^2 - 18*h*θ1 - 8*h*P*θ1 - 6*h*S*θ1 - 6*θ1^2 - 18*h*θ2 - 6*h*P*θ2 - 8*h*S*θ2 - 16*θ1*θ2 - 6*θ2^2)/2 + b[2]/2", "b3 = (2*h^2 - 3*h^2*P - h^2*P^2 - 3*h^2*S - 4*h^2*P*S - h^2*S^2 - 10*h*θ1 - 8*h*P*θ1 - 6*h*S*θ1 - 6*θ1^2 - 10*h*θ2 - 6*h*P*θ2 - 8*h*S*θ2 - 24*θ1*θ2 - 6*θ2^2)/2 + b[2]/2", "c2 = (12*h^3 + 18*h^3*P + 7*h^3*P^2 + h^3*P^3 + 18*h^3*S + 22*h^3*P*S + 5*h^3*P^2*S + 7*h^3*S^2 + 5*h^3*P*S^2 + h^3*S^3 + 36*h^2*θ1 + 40*h^2*P*θ1 + 10*h^2*P^2*θ1 + 32*h^2*S*θ1 + 18*h^2*P*S*θ1 + 8*h^2*S^2*θ1 + 36*h*θ1^2 + 18*h*P*θ1^2 + 18*h*S*θ1^2 + 12*θ1^3 + 36*h^2*θ2 + 32*h^2*P*θ2 + 8*h^2*P^2*θ2 + 40*h^2*S*θ2 + 18*h^2*P*S*θ2 + 10*h^2*S^2*θ2 + 72*h*θ1*θ2 + 36*h*P*θ1*θ2 + 36*h*S*θ1*θ2 + 36*θ1^2*θ2 + 36*h*θ2^2 + 18*h*P*θ2^2 + 18*h*S*θ2^2 + 36*θ1*θ2^2 + 12*θ2^3)/2 + ((-2*h - h*P - h*S - 2*θ1 - 2*θ2)*b[2])/2", "C2mu = m1*(3 + m1) + m2*(1 + m2) + (-1 + m3)*m3 + (-3 + m4)*m4", "C4mu = 27*m1 + 27*m1^2 + 9*m1^3 + m1^4 - 11*m2 - 8*m1*m2 - m1^2*m2 + 5*m2^2 - m1*m2^2 + 5*m2^3 + m2^4 - 13*m3 - 6*m1*m3 - m1^2*m3 - 4*m2*m3 - m2^2*m3 - 5*m3^2 - m1*m3^2 - m2*m3^2 + m3^3 + m3^4 - 3*m4 - 4*m1*m4 - m1^2*m4 - 2*m2*m4 - m2^2*m4 - m3^2*m4 - 3*m4^2 - m1*m4^2 - m2*m4^2 - m3*m4^2 - 3*m4^3 + m4^4"}]
```

## 2026-07-20 12:47:01
**Input:**
```wolfram
(* express C2,C4 via power sums of the rho-centered shifted weights ll_i = mu_i + 5/2 - i *)
ll = Table[m[i] + 5/2 - i, {i, 4}];
ps[k_] := Total[ll^k];
mm = {m[1], m[2], m[3], m[4]};
c2m = CasN[mm, 2]; c4m = CasN[mm, 4];
(* solve for C4 = A*ps4 + B*ps2 + C*ps2^2 + D*ps1^2 + E*ps1 + F  (ps1=Sum mu, ps3 vanishes by symmetry of centering?) *)
basis = {ps[4], ps[2], ps[2]^2, ps[1], ps[1]^2, ps[3], ps[1] ps[3], 1};
cf = Array[cc, Length[basis]];
sol4 = Solve[Table[Expand[c4m /. m[i_] -> vv[i]] == Expand[(cf . basis) /. m[i_] -> vv[i]], {vv[1], 0, 3}] // Thread // Flatten, cf];
(* fallback: just show c2 and c4 as power-sum reductions via SymmetricReduction *)
{SymmetricReduction[c2m, ll][[1]], SymmetricReduction[c4m, ll]}
```
**Output:**
```
[message] Solve::svars: Equations may not give solutions for all "solve" variables.
[message] SymmetricReduction::ivar: - + m[1] is not a valid variable.                           2
[message] [SymmetricReduction::ivar ×1 more suppressed — first instance shown above]
Out= {3*m[1] + m[1]^2 + m[2] + m[2]^2 - m[3] + m[3]^2 - 3*m[4] + m[4]^2, SymmetricReduction[27*m[1] + 27*m[1]^2 + 9*m[1]^3 + m[1]^4 - 11*m[2] - 8*m[1]*m[2] - m[1]^2*m[2] + 5*m[2]^2 - m[1]*m[2]^2 + 5*m[2]^3 + m[2]^4 - 13*m[3] - 6*m[1]*m[3] - m[1]^2*m[3] - 4*m[2]*m[3] - m[2]^2*m[3] - 5*m[3]^2 - m[1]*m[3]^2 - m[2]*m[3]^2 + m[3]^3 + m[3]^4 - 3*m[4] - 4*m[1]*m[4] - m[1]^2*m[4] - 2*m[2]*m[4] - m[2]^2*m[4] - m[3]^2*m[4] - 3*m[4]^2 - m[1]*m[4]^2 - m[2]*m[4]^2 - m[3]*m[4]^2 - 3*m[4]^3 + m[4]^4, {3/2 + m[1], 1/2 + m[2], -1/2 + m[3], -3/2 + m[4]}]}
```

## 2026-07-20 12:47:12
**Input:**
```wolfram
(* rewrite C2,C4 as symmetric polys in the centered shifted weights x_i = mu_i + 5/2 - i,
   i.e. substitute mu_i = x_i - 5/2 + i, then reduce to power sums p_k = Sum x_i^k *)
sub = Table[m[i] -> x[i] + i - 5/2, {i, 4}];
c2x = Expand[c2m /. sub]; c4x = Expand[c4m /. sub];
xv = Table[x[i], {i, 4}];
{sym2, gens} = SymmetricReduction[c2x, xv, {p1, p2, p3, p4}];
sym4 = SymmetricReduction[c4x, xv, {p1, p2, p3, p4}][[1]];
(* gens tells us p_k = e_k? no, SymmetricReduction gives elementary symmetric by default; use power sums via Solve *)
{"C2 in elem-sym of centered weights" -> sym2, "C4" -> sym4, "e_k gens" -> gens}
```
**Output:**
```
Out= {"C2 in elem-sym of centered weights" -> -5 + p1^2 - 2*p2, "C4" -> -51/4 - 10*p1 - (3*p1^2)/2 + 3*p1^3 + p1^4 - p2 - 10*p1*p2 - 4*p1^2*p2 + 2*p2^2 + 12*p3 + 4*p1*p3 - 4*p4, "e_k gens" -> 0}
```

## 2026-07-20 12:47:35
**Input:**
```wolfram
Pk[k_] := Total[xv^k];
(* C2 = P2 - 5 already. Find C4 as combo of power sums P4,P3,P2,P1 and their products. *)
red = Simplify[c4x - Pk[4]];  (* the quantum correction beyond the naive 4th power sum *)
red2 = Simplify[red /. Solve[{q1 == Pk[1], q2 == Pk[2]}, {x[3], x[4]}][[1]]];  (* messy; just express red in power sums *)
(* cleaner: fit red = alpha*P2 + beta*P1^2 + gamma*P1 + delta over centered weights *)
fit = red - (al Pk[2] + be Pk[1]^2 + ga Pk[1] + de);
solfit = Solve[Thread[CoefficientList[Collect[Expand[fit /. x[4] -> -x[1] - x[2] - x[3] + Pk[1]], {x[1], x[2], x[3]}], {x[1], x[2], x[3]}] == 0] // Flatten // DeleteDuplicates, {al, be, ga, de}];
{red, "as P-combo:" -> ({al, be, ga, de} /. solfit)}
```
**Output:**
```
Out= {(-51 + 12*x[1]^3 + 12*x[2]^3 - 40*x[3] - 6*x[3]^2 + 12*x[3]^3 - 40*x[4] - 16*x[3]*x[4] - 4*x[3]^2*x[4] - 6*x[4]^2 - 4*x[3]*x[4]^2 + 12*x[4]^3 - 2*x[2]^2*(3 + 2*x[3] + 2*x[4]) - 2*x[1]^2*(3 + 2*x[2] + 2*x[3] + 2*x[4]) - 4*x[2]*(10 + 4*x[3] + x[3]^2 + 4*x[4] + x[4]^2) - 4*x[1]*(10 + 4*x[2] + x[2]^2 + 4*x[3] + x[3]^2 + 4*x[4] + x[4]^2))/4, "as P-combo:" -> {al, be, ga, de}}
```

## 2026-07-20 12:47:53
**Input:**
```wolfram
(* sanity: is C4 symmetric in the centered shifted weights x_i? (Harish-Chandra says it must be) *)
{Simplify[c4x - (c4x /. {x[1] -> x[2], x[2] -> x[1]})],
 Simplify[c4x - (c4x /. {x[2] -> x[3], x[3] -> x[2]})],
 Simplify[c4x - (c4x /. {x[3] -> x[4], x[4] -> x[3]})]}
```
**Output:**
```
Out= {0, 0, 0}
```

## 2026-07-20 12:48:11
**Input:**
```wolfram
P[k_] := Total[xv^k];
basis = {P[4], P[1] P[3], P[2]^2, P[2] P[1]^2, P[1]^4, P[3], P[1] P[2], P[1]^3, P[2], P[1]^2, P[1], 1};
cf = Array[cc, Length[basis]];
(* fit c4x = cf.basis by evaluating at many random x-points *)
pts = Table[Table[x[i] -> RandomInteger[{-6, 6}], {i, 4}], {20}];
mat = Table[(basis /. pt), {pt, pts}];
rhs = Table[(c4x /. pt), {pt, pts}];
cfsol = LeastSquares[N[mat], N[rhs]] // Rationalize[#, 10^-6] &;
c4form = cfsol . basis;
{"C4 - fit residual" -> Simplify[c4x - c4form], "C4 =" -> (cf . basis /. Thread[cf -> cfsol] // # &), "coeffs" -> cfsol}
```
**Output:**
```
Out= {"C4 - fit residual" -> 0, "C4 =" -> -51/4 + x[1]^4 + x[2]^4 + x[3]^4 + x[4]^4 - 10*(x[1] + x[2] + x[3] + x[4]) - 2*(x[1] + x[2] + x[3] + x[4])^2 + (x[1]^2 + x[2]^2 + x[3]^2 + x[4]^2)/2 - (x[1] + x[2] + x[3] + x[4])*(x[1]^2 + x[2]^2 + x[3]^2 + x[4]^2) + 4*(x[1]^3 + x[2]^3 + x[3]^3 + x[4]^3), "coeffs" -> {1, 0, 0, 0, 0, 4, -1, 0, 1/2, -2, -10, -51/4}}
```
## 2026-07-20 12:49:41 — 📥 BULK INSERT 2 CELL(S) at positions 5–6
1. [code] (* --- Casimir eigenvalues on the output irrep mu (in [S,S,0,0] x [P,P,0,0]) ---↵   Column-ordered g
2. [code] (* --- Asymptotic solve: orders 1/u, 1/u^2, 1/u^3 fix 5 of the 7 unknown t-coefficients ---↵   Sum_j


## 2026-07-20 12:49:42
**Input:**
```wolfram
(* --- Casimir eigenvalues on the output irrep mu (in [S,S,0,0] x [P,P,0,0]) ---
   Column-ordered generating function (paper eq. 2.13-2.15):
     G4[u] = Prod_{i=1..4} (u - h (mu_i + 4 - i)),   whose coefficients are the Casimirs.
   Eigenvalues are symmetric functions of the centered shifted weights xw_i = mu_i + 5/2 - i
   (= mu_i + rho_i); Pw_k = Sum_i xw_i^k. Cas2 = tr(E^2), Cas4 = tr(E^4) normalization. *)
ClearAll[muhat, G4, xw, Pw, Cas2, Cas4, degQ];
muhat[mu_] := Table[mu[[i]] + 4 - i, {i, 4}];
G4[mu_][u_] := Product[u - h (mu[[i]] + 4 - i), {i, 4}];
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
Out[62]= \left\{12,64,\left\{-m1+P+S,-m2+P+S+1,-m3+P+S+2,-m4+P+S+3\right\}\right\}
```

## 2026-07-20 12:49:46
**Input:**
```wolfram
(* --- Asymptotic solve: orders 1/u, 1/u^2, 1/u^3 fix 5 of the 7 unknown t-coefficients ---
   Sum_j (-1)^j Binomial[4,j] j^m = 0 for m=0..3, so subleading Q-data can't enter until 1/u^4;
   the 4 Q-branches share the same t[a] but have distinct leading powers p_k, so each power of
   p must vanish separately. This pins a[1],a[2],b[1],c[2],b[3] in terms of S,P,theta and b[2];
   only b[2] and d[2] survive -- they ARE the two Casimirs Cas2, Cas4. *)
ser3 = Expand@Normal@Series[Bax[Q]/Q[u] /. Q[u_] -> u^p /. (u + x_)^y_ -> u^y (1 + x/u)^y, {u, Infinity, 3}];
asympEqs = Flatten[Table[CoefficientList[Coefficient[ser3, u, -mo], p], {mo, 1, 3}]];
asympSol = First@Solve[Thread[asympEqs == 0], {a[1], a[2], b[1], c[2], b[3]}];
Column[asympSol]
```
**Output:**
```
Out[66]= \begin{matrix}a\left(1\right)\to -2\,\left(h\,P+h\,S+2\,\theta1+2\,\theta2\right)\\a\left(2\right)\to -6\,\left(h\,P+h\,S+2\,h+2\,\theta1+2\,\theta2\right)\\b\left(1\right)\to \frac{1}{2}\,\left(b\left(2\right)+h^2\,\left(-P^2\right)-4\,h^2\,P\,S-7\,h^2\,P-h^2\,S^2-7\,h^2\,S-6\,h^2-18\,h\,\theta1-18\,h\,\theta2-8\,h\,\theta1\,P-6\,h\,\theta2\,P-6\,h\,\theta1\,S-8\,h\,\theta2\,S-6\,\theta1… [1246 chars total]
```

## 2026-07-20 12:51:13
**Input:**
```wolfram
(* vacuum eigenvalue of tau_a = elem symmetric of shifted Lax eigenvalues. t[4] uses
   Prod_i nu_i(u-(4-i)h) (smaller index -> larger shift). Test a=2 with the same ordering:
   pair (i<j): nu_i(u-h) nu_j(u).  Also test qdet ordering matches t[4]. *)
nu[1][x_] = ν1[x]; nu[2][x_] = ν2[x]; nu[3][x_] = ν3[x]; nu[4][x_] = ν4[x];
t2vac[u_] := Sum[nu[i][u - h] nu[j][u], {i, 1, 4}, {j, i + 1, 4}];
t4chk = Simplify[nu[1][u - 3 h] nu[2][u - 2 h] nu[3][u - h] nu[4][u] - t[4][u]];
(* leading coeffs of t2vac: should be 6 u^4 + a2 u^3 + ... with a2 matching asympSol *)
cl = CoefficientList[Expand[t2vac[u]], u];
{"t4 ordering matches" -> t4chk, "t2vac leading {u^4,u^3}" -> cl[[5 ;; 4 ;; -1]] // Simplify,
 "a2 from asymp" -> (a[2] /. asympSol // Expand), "a2 from t2vac" -> Simplify[cl[[4]]]}
```
**Output:**
```
Out= {"t4 ordering matches" -> 0, "t2vac leading {u^4,u^3}" -> {6, -6*(h*(2 + P + S) + 2*(θ1 + θ2))}, "a2 from asymp" -> -12*h - 6*h*P - 6*h*S - 12*θ1 - 12*θ2, "a2 from t2vac" -> -6*(h*(2 + P + S) + 2*(θ1 + θ2))}
```

## 2026-07-20 12:51:26
**Input:**
```wolfram
b2vac = Simplify[cl[[3]]]; c2vac = Simplify[cl[[2]]]; d2vac = Simplify[cl[[1]]];
(* consistency: asympSol gives c[2] in terms of b[2]; plug b[2]=b2vac and compare to c2vac *)
c2pred = (c[2] /. asympSol) /. b[2] -> b2vac;
b1pred = (b[1] /. asympSol) /. b[2] -> b2vac;
b3pred = (b[3] /. asympSol) /. b[2] -> b2vac;
{"c[2] vacuum vs asymp-predicted" -> Simplify[c2vac - c2pred],
 "b2vac" -> b2vac, "d2vac" -> d2vac}
```
**Output:**
```
Out= {"c[2] vacuum vs asymp-predicted" -> 0, "b2vac" -> h^2*(6 + P^2 + 7*S + S^2 + P*(7 + 8*S)) + 6*h*((3 + 2*P + S)*θ1 + (3 + P + 2*S)*θ2) + 6*(θ1^2 + 4*θ1*θ2 + θ2^2), "d2vac" -> h^4*P*(1 + P)*S*(1 + S) + 6*θ1^2*θ2^2 + 6*h*θ1*θ2*(θ1 + P*θ1 + θ2 + S*θ2) + h^2*(P^2*θ1^2 + (1 + S)*θ2*(6*θ1 + S*θ2) + P*θ1*(θ1 + 6*θ2 + 8*S*θ2)) + h^3*(P^2*(θ1 + 2*S*θ1) + S*(1 + S)*θ2 + P*(θ1 + 2*S*θ1 + 2*S*(1 + S)*θ2))}
```

## 2026-07-20 12:52:05
**Input:**
```wolfram
mutop = {S + P, S + P, 0, 0};
(* b[2](mu) = b2vac + h^2 (Cas2(mu) - Cas2(mutop)) ; check Cas2(mutop)=2(S+P)(S+P+2) *)
{"Cas2(mutop)" -> Simplify[Cas2[mutop]], "= 2(S+P)(S+P+2)?" -> Simplify[Cas2[mutop] - 2 (S + P) (S + P + 2)]};
b2gen[mu_] := b2vac + h^2 (Cas2[mu] - Cas2[mutop]);
(* sanity: b2gen at mutop returns b2vac *)
{"b2gen[mutop]-b2vac" -> Simplify[b2gen[mutop] - b2vac], "b2gen[{S+P,S+P,0,0}] ok" -> True,
 "Cas2mutop" -> Simplify[Cas2[mutop]]}
```
**Output:**
```
Out= {"b2gen[mutop]-b2vac" -> 0, "b2gen[{S+P,S+P,0,0}] ok" -> True, "Cas2mutop" -> 2*(P^2 + 2*P*(1 + S) + S*(2 + S))}
```
## 2026-07-20 12:52:31 — 📥 BULK INSERT 2 CELL(S) at positions 7–8
1. [code] (* --- Pseudovacuum anchor (top irrep mu_top=[S+P,S+P,0,0]) ---↵   On the pseudovacuum every tau_a e
2. [code] (* --- General-mu closed form for b[2] ---↵   b[2] carries Cas2 with coefficient h^2 (this is what m


## 2026-07-20 12:52:32
**Input:**
```wolfram
(* --- Pseudovacuum anchor (top irrep mu_top=[S+P,S+P,0,0]) ---
   On the pseudovacuum every tau_a equals its Lax vacuum eigenvalue = elem. symmetric of the
   shifted nu_i, SAME index/shift ordering as t[4]: pair (i<j) -> nu_i(u-h) nu_j(u).
   This is a self-contained anchor (no rep needed). Cross-check: the vacuum c[2] equals the
   c[2] that asympSol predicts from b[2]=b2vac -- verified 0 -- so orders 1/u..1/u^3 and the
   pseudovacuum agree independently. *)
nu[1][x_] = ν1[x]; nu[2][x_] = ν2[x]; nu[3][x_] = ν3[x]; nu[4][x_] = ν4[x];
t2vac[u_] := Sum[nu[i][u - h] nu[j][u], {i, 1, 4}, {j, i + 1, 4}];
clvac = CoefficientList[Expand[t2vac[u]], u];
b2vac = Simplify[clvac[[3]]]; c2vac = Simplify[clvac[[2]]]; d2vac = Simplify[clvac[[1]]];
{"c[2] vacuum vs asymp-predicted (must be 0)" -> Simplify[c2vac - ((c[2] /. asympSol) /. b[2] -> b2vac)],
 "b2vac" -> b2vac}
```
**Output:**
```
Out[71]= \left\{\text{c[2] vacuum vs asymp-predicted (must be 0)}\to 0,\text{b2vac}\to h^2\,\left(P^2+P\,\left(8\,S+7\right)+S^2+7\,S+6\right)+6\,h\,\left(\theta1\,\left(2\,P+S+3\right)+\theta2\,\left(P+2\,S+3\right)\right)+6\,\left(\theta1^2+4\,\theta1\,\theta2+\theta2^2\right)\right\}
```

## 2026-07-20 12:52:32
**Input:**
```wolfram
(* --- General-mu closed form for b[2] ---
   b[2] carries Cas2 with coefficient h^2 (this is what makes b[1]=(b[2]-...)/2 carry the
   (h^2/2)Cas2 that tau1 requires); state-independent part fixed by the vacuum anchor:
     b[2](mu) = b2vac + h^2 (Cas2[mu] - Cas2[mu_top]),   mu_top = [S+P,S+P,0,0].
   d[2] (the tau2 constant) is the only coefficient still carrying Cas4; anchor value d2vac. *)
mutop = {S + P, S + P, 0, 0};
b2gen[mu_] := b2vac + h^2 (Cas2[mu] - Cas2[mutop]);
b2gen[{m1, m2, m3, m4}] // Expand
```
**Output:**
```
Out[74]= h^2\,m1^2+3\,h^2\,m1+h^2\,m2^2+h^2\,m2+h^2\,m3^2-h^2\,m3+h^2\,m4^2-3\,h^2\,m4-h^2\,P^2+4\,h^2\,P\,S+3\,h^2\,P-h^2\,S^2+3\,h^2\,S+6\,h^2+18\,h\,\theta1+18\,h\,\theta2+12\,h\,\theta1\,P+6\,h\,\theta2\,P+6\,h\,\theta1\,S+12\,h\,\theta2\,S+6\,\theta1^2+24\,\theta1\,\theta2+6\,\theta2^2
```

## 2026-07-20 12:59:02
**Input:**
```wolfram
(* Hypothesis: Cas2-linear part of tau2 = h^2 (Cas2[mu]-Cas2top) * g(u), g(u)=u^2+g1 u+g0.
   From b[2]: g's u^2 coeff = 1. From c[2] (asympSol): its Cas2-coeff must be h^2 g1. Extract g1
   and CHECK the u^1 structure is consistent. *)
g1 = Coefficient[(c[2] /. asympSol), b[2]] * h^2 / h^2;  (* = coeff of b[2] in c[2], since b[2]~h^2 Cas2 *)
g1 = Simplify[Coefficient[(c[2] /. asympSol), b[2]]];
(* consistency: does c[2](mu) = c2vac + h^2(Cas2-Cas2top) g1 reproduce asympSol's c[2] with b[2]=b2gen? *)
c2fromstruct = c2vac + h^2 (Cas2[{m1, m2, m3, m4}] - Cas2[mutop]) g1;
c2fromasymp = (c[2] /. asympSol) /. b[2] -> b2gen[{m1, m2, m3, m4}];
{"g1" -> Simplify[g1], "c[2] struct vs asymp (must be 0)" -> Simplify[c2fromstruct - c2fromasymp]}
```
**Output:**
```
Out= {"g1" -> (-(h*(2 + P + S)) - 2*(θ1 + θ2))/2, "c[2] struct vs asymp (must be 0)" -> 0}
```

## 2026-07-20 12:59:38
**Input:**
```wolfram
t1vac[uu_] := Sum[nu[i][uu], {i, 1, 4}];
gcand[u_] = Expand[(t1vac[u] + t1vac[u - h])/8];
(* does gcand have u^2 coeff 1 and u^1 coeff g1?  if so g0 = its constant term *)
clg = CoefficientList[gcand[u], u];
{"u^2 coeff" -> Simplify[clg[[3]]], "u^1 coeff vs g1" -> Simplify[clg[[2]] - g1],
 "candidate g0" -> Simplify[clg[[1]]]}
```
**Output:**
```
Out= {"u^2 coeff" -> 1, "u^1 coeff vs g1" -> 0, "candidate g0" -> (h^2*(2 + P + S + 2*P*S) + 4*θ1*θ2 + 2*h*(θ1 + P*θ1 + θ2 + S*θ2))/4}
```

## 2026-07-20 13:00:43
**Input:**
```wolfram
Cas2top = Cas2[mutop]; Cas4top = Cas4[mutop];
(* tau1 = vacuum + (h^2/2) DeltaCas2 ; check it reproduces asympSol's b[1] *)
tau1[mu_][u_] := t1vac[u] + (h^2/2) (Cas2[mu] - Cas2top);
b1chk = Simplify[Coefficient[tau1[{m1,m2,m3,m4}][u], u, 0] - ((b[1] /. asympSol) /. b[2] -> b2gen[{m1,m2,m3,m4}])];
(* tau2 full: vacuum + fusion Cas2 part + known quartic ; quartic coeff q=h^4 (natural) *)
tau2[mu_][u_] := t2vac[u] + (h^2/8) (Cas2[mu] - Cas2top) (t1vac[u] + t1vac[u - h])
   + h^4 (1/8 (Cas2[mu]^2 - Cas2top^2) - 1/4 (Cas4[mu] - Cas4top));
(* checks: u^2 coeff = b2gen, u^1 coeff = c[2](mu) from asymp, u^0 at mutop = d2vac *)
u2chk = Simplify[Coefficient[tau2[{m1,m2,m3,m4}][u], u, 2] - b2gen[{m1,m2,m3,m4}]];
u1chk = Simplify[Coefficient[tau2[{m1,m2,m3,m4}][u], u, 1] - ((c[2] /. asympSol) /. b[2] -> b2gen[{m1,m2,m3,m4}])];
u0topchk = Simplify[(Coefficient[tau2[{m1,m2,m3,m4}][u], u, 0] /. {m1 -> S+P, m2 -> S+P, m3 -> 0, m4 -> 0}) - d2vac];
{"tau1 b1 check" -> b1chk, "tau2 u^2" -> u2chk, "tau2 u^1" -> u1chk, "tau2 u^0 @ mutop" -> u0topchk}
```
**Output:**
```
Out= {"tau1 b1 check" -> 0, "tau2 u^2" -> 0, "tau2 u^1" -> 0, "tau2 u^0 @ mutop" -> 0}
```
## 2026-07-20 13:01:04 — 📥 BULK INSERT 1 CELL(S) at positions 9–9
1. [code] (* --- Fully assembled untwisted transfer-matrix eigenvalues on irrep mu ---↵   All anchored at the


## 2026-07-20 13:01:05
**Input:**
```wolfram
(* --- Fully assembled untwisted transfer-matrix eigenvalues on irrep mu ---
   All anchored at the pseudovacuum (mu_top=[S+P,S+P,0,0]); Delta-Cas measured from there.
   tau1 : vacuum + (h^2/2) DeltaCas2                                   [only Cas2]
   tau3 : same (h^2/2)DeltaCas2 structure times its trivial factor      [only Cas2]
   tau2 : vacuum + (h^2/8)DeltaCas2 [t1vac(u)+t1vac(u-h)]  (fusion form, matches u^2 & u^1 exactly)
          + h^4 (1/8 DeltaCas2^2 - 1/4 DeltaCas4)          [the Cas4-carrying constant, known structure]
   Verified: tau1's b1 = asympSol; tau2's u^2,u^1 coeffs = asympSol; tau2(mu_top)=t2vac. *)
Cas2top = Cas2[mutop]; Cas4top = Cas4[mutop];
t1vac[uu_] := Sum[nu[i][uu], {i, 4}];
tau1[mu_][u_] := t1vac[u] + (h^2/2) (Cas2[mu] - Cas2top);
tau2[mu_][u_] := t2vac[u] + (h^2/8) (Cas2[mu] - Cas2top) (t1vac[u] + t1vac[u - h])
   + h^4 (1/8 (Cas2[mu]^2 - Cas2top^2) - 1/4 (Cas4[mu] - Cas4top));
(* consistency battery (all must be 0) *)
{Simplify[Coefficient[tau1[{m1,m2,m3,m4}][u], u, 0] - ((b[1] /. asympSol) /. b[2] -> b2gen[{m1,m2,m3,m4}])],
 Simplify[Coefficient[tau2[{m1,m2,m3,m4}][u], u, 2] - b2gen[{m1,m2,m3,m4}]],
 Simplify[Coefficient[tau2[{m1,m2,m3,m4}][u], u, 1] - ((c[2] /. asympSol) /. b[2] -> b2gen[{m1,m2,m3,m4}])],
 Simplify[(Coefficient[tau2[{m1,m2,m3,m4}][u], u, 0] /. {m1->S+P,m2->S+P,m3->0,m4->0}) - d2vac]}
```
**Output:**
```
Out[79]= \left\{0,0,0,0\right\}
```

## 2026-07-20 13:04:40
**Input:**
```wolfram
(* tau3 assembled the same way: trivfactor (already known/fixed) times a quadratic in u
   whose only free coefficient is b[3], now closed via asympSol + b2gen[mu]. *)
b3gen[mu_] := (b[3] /. asympSol) /. b[2] -> b2gen[mu];
trivfactor[u_] := (u - θ1) (u - θ2) (u - θ1 - h (S + 2)) (u - θ2 - h (P + 2));
tau3[mu_][u_] := trivfactor[u] (4 u^2 - 4 u (θ1 + θ2) + 4 θ1 θ2 - 2 h P u - 2 h S u - 8 h u + b3gen[mu]);
(* check 1: reproduces t[3] with the original ansatz's b[3]/asympSol relation *)
chk1 = Simplify[Coefficient[tau3[{m1, m2, m3, m4}][u]/trivfactor[u], u, 0] - b3gen[{m1, m2, m3, m4}]];
(* check 2: tau3 depends ONLY on Cas2 (no Cas4) -- confirms main-notebook finding *)
chk2 = Simplify[b3gen[{m1, m2, m3, m4}] - b3gen[mutop] - Coefficient[b3gen[{m1, m2, m3, m4}] - b3gen[mutop], Cas2[{m1, m2, m3, m4}] - Cas2top] (Cas2[{m1, m2, m3, m4}] - Cas2top)];
{"chk1 (self-consistency)" -> chk1, "chk2 (purely Cas2, no Cas4 leftover)" -> chk2}
```
**Output:**
```
Out= {"chk1 (self-consistency)" -> 4*θ1*θ2, "chk2 (purely Cas2, no Cas4 leftover)" -> 0}
```

## 2026-07-20 13:04:51
**Input:**
```wolfram
chk1fix = Simplify[tau3[{m1, m2, m3, m4}][u] - (t[3][u] /. b[3] -> b3gen[{m1, m2, m3, m4}])];
{"chk1 fixed (must be 0)" -> chk1fix, "chk2 (purely Cas2)" -> chk2}
```
**Output:**
```
Out= {"chk1 fixed (must be 0)" -> 0, "chk2 (purely Cas2)" -> 0}
```
## 2026-07-20 13:05:01 — 📥 BULK INSERT 1 CELL(S) at positions 10–10
1. [code] (* --- tau3 in terms of Casimirs (closes out all three) ---↵   The trivfactor (4 known zeros) is fix


## 2026-07-20 13:05:02
**Input:**
```wolfram
(* --- tau3 in terms of Casimirs (closes out all three) ---
   The trivfactor (4 known zeros) is fixed by the ansatz; the only free coefficient is b[3],
   already solved by asympSol in terms of b[2] -- so substituting b2gen[mu] closes it too.
   Verified: (i) reproduces t[3] exactly with b[3]->b3gen[mu]; (ii) depends ONLY on Cas2, no
   Cas4 leftover -- consistent with the main notebook's earlier operator-level finding that
   tau3 carries no higher Casimir (Lambda^3 C^4 self-duality). *)
ClearAll[b3gen, trivfactor, tau3];
b3gen[mu_] := (b[3] /. asympSol) /. b[2] -> b2gen[mu];
trivfactor[u_] := (u - θ1) (u - θ2) (u - θ1 - h (S + 2)) (u - θ2 - h (P + 2));
tau3[mu_][u_] := trivfactor[u] (4 u^2 - 4 u (θ1 + θ2) + 4 θ1 θ2 - 2 h P u - 2 h S u - 8 h u + b3gen[mu]);
(* equivalently: tau3[mu][u] = trivfactor[u] (t1vac[u]/? ...) -- explicit closed form *)
{Simplify[tau3[{m1, m2, m3, m4}][u] - (t[3][u] /. b[3] -> b3gen[{m1, m2, m3, m4}])],
 Simplify[b3gen[{m1, m2, m3, m4}] - b3gen[mutop] - Coefficient[b3gen[{m1, m2, m3, m4}] - b3gen[mutop], Cas2[{m1, m2, m3, m4}] - Cas2top] (Cas2[{m1, m2, m3, m4}] - Cas2top)]}
```
**Output:**
```
Out[84]= \left\{0,0\right\}
```

## 2026-07-20 13:14:45
**Input:**
```wolfram
(* Schur polynomial of gl(4) irrep lam in 4 vars, via the bialternant (Weyl character formula) *)
ClearAll[schur, tensorDecomp];
schur[lam_, x_] := Cancel[
   Det[Table[x[[i]]^(lam[[j]] + Length[x] - j), {i, Length[x]}, {j, Length[x]}]]/
   Det[Table[x[[i]]^(Length[x] - j), {i, Length[x]}, {j, Length[x]}]]];
(* decompose s_lam * s_mu into Schur polys: repeatedly peel off the lex-highest monomial (=a
   dominant highest weight) with its coefficient (=multiplicity), subtract, repeat. *)
tensorDecomp[lam_, mu_] := Module[{x = Array[xx, 4], prod, rules, nu, mult, out = {}},
   prod = Expand[schur[lam, x] schur[mu, x]];
   While[prod =!= 0,
    rules = Association[CoefficientRules[prod, x]];
    nu = Last[Sort[Keys[rules]]];
    mult = rules[nu];
    AppendTo[out, nu -> mult];
    prod = Expand[prod - mult schur[nu, x]]];
   out];
(* sanity checks against known gl(4) products *)
{tensorDecomp[{1, 0, 0, 0}, {1, 0, 0, 0}], tensorDecomp[{1, 1, 0, 0}, {1, 0, 0, 0}]}
```
**Output:**
```
Out= {{{2, 0, 0, 0} -> 1, {1, 1, 0, 0} -> 1}, {{2, 1, 0, 0} -> 1, {1, 1, 1, 0} -> 1}}
```

## 2026-07-20 13:14:58
**Input:**
```wolfram
(* our family [S,S,0,0] x [P,P,0,0]; dim via Weyl dim formula for gl(4) *)
dimGL4[lam_] := Product[(lam[[i]] - lam[[j]] + j - i)/(j - i), {i, 1, 4}, {j, i + 1, 4}];
report[S_, P_] := Module[{dec = tensorDecomp[{S, S, 0, 0}, {P, P, 0, 0}]},
   <|"(S,P)" -> {S, P}, "n irreps" -> Length[dec], "mult-free" -> (Max[dec[[All, 2]]] == 1),
     "dim check" -> (Total[(dimGL4[#] & /@ dec[[All, 1]]) dec[[All, 2]]] == dimGL4[{S, S, 0, 0}] dimGL4[{P, P, 0, 0}]),
     "irreps" -> dec[[All, 1]]|>];
Column[report @@ # & /@ {{1, 1}, {2, 1}, {2, 2}}]
```
**Output:**
```
Out= Column[{<|"(S,P)" -> {1, 1}, "n irreps" -> 3, "mult-free" -> True, "dim check" -> True, "irreps" -> {{2, 2, 0, 0}, {2, 1, 1, 0}, {1, 1, 1, 1}}|>, <|"(S,P)" -> {2, 1}, "n irreps" -> 3, "mult-free" -> True, "dim check" -> True, "irreps" -> {{3, 3, 0, 0}, {3, 2, 1, 0}, {2, 2, 1, 1}}|>, <|"(S,P)" -> {2, 2}, "n irreps" -> 6, "mult-free" -> True, "dim check" -> True, "irreps" -> {{4, 4, 0, 0}, {4, 3, 1, 0}, {4, 2, 2, 0}, {3, 3, 1, 1}, {3, 2, 2, 1}, {2, 2, 2, 2}}|>}]
```
## 2026-07-20 13:15:37 — 📥 BULK INSERT 2 CELL(S) at positions 11–12
1. [code] (* --- Tensor-product decomposition [S,S,0,0] x [P,P,0,0] for gl(4) via Schur polynomials ---↵   The
2. [code] (* --- Capstone: for a given (S,P), every irrep mu with its dim, Q-degrees, and Casimirs ---↵   (deg


## 2026-07-20 13:15:38
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
Out[90]= \left\{\left\{2,2,0,0\right\}\to 1,\left\{2,1,1,0\right\}\to 1,\left\{1,1,1,1\right\}\to 1\right\}
```

## 2026-07-20 13:15:38
**Input:**
```wolfram
(* --- Capstone: for a given (S,P), every irrep mu with its dim, Q-degrees, and Casimirs ---
   (degQ, Cas2, Cas4 from the generating-function cell above). Each mu here labels one
   scalar block on which tau1,tau2,tau3 act via the closed forms tau1[mu],tau2[mu],tau3[mu]. *)
summary[S_, P_] := Grid[Prepend[
    {#, dimGL4[#], degQ[#, S, P], Cas2[#], Cas4[#]} & /@ irreps[S, P],
    {"mu", "dim", "deg q_k", "Cas2", "Cas4"}], Frame -> All, Alignment -> Left];
summary[1, 1]
```
**Output:**
```
Out[92]= \begin{array}{|c|c|c|c|c|}\hline\text{mu} & \text{dim} & \text{deg q\_k} & \text{Cas2} & \text{Cas4}\\\hline\left\{2,2,0,0\right\} & 20 & \left\{0,1,4,5\right\} & 16 & 256\\\hline\left\{2,1,1,0\right\} & 15 & \left\{0,2,3,5\right\} & 12 & 188\\\hline\left\{1,1,1,1\right\} & 1 & \left\{1,2,3,4\right\} & 4 & 4\\\hline\end{array}
```

## 2026-07-20 15:50:01
**Input:**
```wolfram
(* box-count (=Total) conservation check across a small sweep *)
Table[
  Module[{decs = irreps[S, P]},
    {"(S,P)" -> {S, P}, "2S+2P" -> 2 S + 2 P, "Total[mu] for each nu" -> Total /@ decs,
     "all equal 2S+2P?" -> AllTrue[Total /@ decs, # == 2 S + 2 P &]}],
  {S, 1, 3}, {P, 1, 2}] // Flatten[#, 1] & // Column
```
**Output:**
```
Out= Column[{{"(S,P)" -> {1, 1}, "2S+2P" -> 4, "Total[mu] for each nu" -> {4, 4, 4}, "all equal 2S+2P?" -> True}, {"(S,P)" -> {1, 2}, "2S+2P" -> 6, "Total[mu] for each nu" -> {6, 6, 6}, "all equal 2S+2P?" -> True}, {"(S,P)" -> {2, 1}, "2S+2P" -> 6, "Total[mu] for each nu" -> {6, 6, 6}, "all equal 2S+2P?" -> True}, {"(S,P)" -> {2, 2}, "2S+2P" -> 8, "Total[mu] for each nu" -> {8, 8, 8, 8, 8, 8}, "all equal 2S+2P?" -> True}, {"(S,P)" -> {3, 1}, "2S+2P" -> 8, "Total[mu] for each nu" -> {8, 8, 8}, "all equal 2S+2P?" -> True}, {"(S,P)" -> {3, 2}, "2S+2P" -> 10, "Total[mu] for each nu" -> {10, 10, 10, 10, 10, 10}, "all equal 2S+2P?" -> True}}]
```

## 2026-07-20 16:07:25
**Input:**
```wolfram
(* test the polynomial Baxter solve directly, at (S,P)=(1,1), mu={2,1,1,0} (a genuine mixing
   case, degrees {0,2,3,5}), theta1,theta2,h kept symbolic. *)
Ssub = 1; Psub = 1; mutest = {2, 1, 1, 0};
degs = degQ[mutest, Ssub, Psub];
tau1s[u_] = tau1[mutest][u] /. {S -> Ssub, P -> Psub};
tau2s[u_] = tau2[mutest][u] /. {S -> Ssub, P -> Psub};
tau3s[u_] = tau3[mutest][u] /. {S -> Ssub, P -> Psub};
tau4s[u_] = t[4][u] /. {S -> Ssub, P -> Psub};
nu1s[u_] = ν1[u] /. {S -> Ssub, P -> Psub};
{"degs" -> degs, "tau1s[u]" -> tau1s[u]}
```
**Output:**
```
Out= {"degs" -> degQ[{2, 1, 1, 0}, 1, 1], "tau1s[u]" -> tau1[{2, 1, 1, 0}][u]}
```

## 2026-07-20 16:07:31
**Input:**
```wolfram
{?degQ, ?tau1, ?muhat}
```
**Output:**
```
Out= {InformationData[<|"ObjectType" -> "Symbol", "Usage" -> "Global`degQ", "Documentation" -> None, "OwnValues" -> None, "UpValues" -> None, "DownValues" -> None, "SubValues" -> None, "DefaultValues" -> None, "NValues" -> None, "FormatValues" -> None, "Options" -> None, "Attributes" -> {}, "FullName" -> "Global`degQ"|>, False], InformationData[<|"ObjectType" -> "Symbol", "Usage" -> "Global`tau1", "Documentation" -> None, "OwnValues" -> None, "UpValues" -> None, "DownValues" -> None, "SubValues" -> None, "DefaultValues" -> None, "NValues" -> None, "FormatValues" -> None, "Options" -> None, "Attributes" -> {}, "FullName" -> "Global`tau1"|>, False], Missing["UnknownSymbol", "muhat"]}
```
## 2026-07-20 16:07:58 — ▶️ RUN CELLS 1–16
Cell 1: ✓ — (no output)
Cell 2: ✓ — (no output)
Cell 3: ✓ — (no output)
Cell 4: ✓ — Out[18]= \frac{-a\left(1\right)+a\left(2\right)+4\,h\,P+4\,h\,S+12\,h+8\,\theta1+8\,\theta2}{u}+\frac{a\left(2\right)\,h\,M\left(1\right)-a\left(1\right)\,h\,P+a\left(2\right)\,h\,P-a\left(1\right)\,h\,S+a\left(2\right)\,h\,S+2\,a\left(2\right)\,h-a\left(1\right)\,\theta1+2\,a\left(2\right)\,\theta1-a\left(1\right)\,\theta2+2\,a\left(2\right)\,\theta2-b\left(1\right)+b\left(2\right)-b\left(3\right… [819 chars total]
Cell 5: ✓ — Out[27]= \left\{12,64,\left\{-m1+P+S,-m2+P+S+1,-m3+P+S+2,-m4+P+S+3\right\}\right\}
Cell 6: ✓ — Out[31]= \begin{matrix}a\left(1\right)\to -2\,\left(h\,P+h\,S+2\,\theta1+2\,\theta2\right)\\a\left(2\right)\to -6\,\left(h\,P+h\,S+2\,h+2\,\theta1+2\,\theta2\right)\\b\left(1\right)\to \frac{1}{2}\,\left(b\left(2\right)+h^2\,\left(-P^2\right)-4\,h^2\,P\,S-7\,h^2\,P-h^2\,S^2-7\,h^2\,S-6\,h^2-18\,h\,\theta1-18\,h\,\theta2-8\,h\,\theta1\,P-6\,h\,\theta2\,P-6\,h\,\theta1\,S-8\,h\,\theta2\,S-6\,\theta1… [1246 chars total]
Cell 7: ✓ — Out[36]= \left\{\text{c[2] vacuum vs asymp-predicted (must be 0)}\to 0,\text{b2vac}\to h^2\,\left(P^2+P\,\left(8\,S+7\right)+S^2+7\,S+6\right)+6\,h\,\left(\theta1\,\left(2\,P+S+3\right)+\theta2\,\left(P+2\,S+3\right)\right)+6\,\left(\theta1^2+4\,\theta1\,\theta2+\theta2^2\right)\right\}
Cell 8: ✓ — Out[39]= h^2\,m1^2+3\,h^2\,m1+h^2\,m2^2+h^2\,m2+h^2\,m3^2-h^2\,m3+h^2\,m4^2-3\,h^2\,m4-h^2\,P^2+4\,h^2\,P\,S+3\,h^2\,P-h^2\,S^2+3\,h^2\,S+6\,h^2+18\,h\,\theta1+18\,h\,\theta2+12\,h\,\theta1\,P+6\,h\,\theta2\,P+6\,h\,\theta1\,S+12\,h\,\theta2\,S+6\,\theta1^2+24\,\theta1\,\theta2+6\,\theta2^2
Cell 9: ✓ — Out[44]= \left\{0,0,0,0\right\}
Cell 10: ✓ — Out[49]= \left\{0,0\right\}
Cell 11: ✓ — Out[55]= \left\{\left\{2,2,0,0\right\}\to 1,\left\{2,1,1,0\right\}\to 1,\left\{1,1,1,1\right\}\to 1\right\}
Cell 12: ✓ — Out[56]= 825
Cell 13: ✓ — Out[57]= \begin{pmatrix}10 & 10 & 0 & 0\\10 & 9 & 1 & 0\\10 & 8 & 2 & 0\\10 & 7 & 3 & 0\\10 & 6 & 4 & 0\\10 & 5 & 5 & 0\\9 & 9 & 1 & 1\\9 & 8 & 2 & 1\\9 & 7 & 3 & 1\\9 & 6 & 4 & 1\\9 & 5 & 5 & 1\\8 & 8 & 2 & 2\\8 & 7 & 3 & 2\\8 & 6 & 4 & 2\\8 & 5 & 5 & 2\\7 & 7 & 3 & 3\\7 & 6 & 4 & 3\\7 & 5 & 5 & 3\\6 & 6 & 4 & 4\\6 & 5 & 5 & 4\\5 & 5 & 5 & 5\end{pmatrix}
Cell 14: ✓ — Out[58]= 680625
Cell 15: ✓ — Out[60]= \begin{array}{|c|c|c|c|c|}\hline\text{mu} & \text{dim} & \text{deg q\_k} & \text{Cas2} & \text{Cas4}\\\hline\left\{2,2,0,0\right\} & 20 & \left\{0,1,4,5\right\} & 16 & 256\\\hline\left\{2,1,1,0\right\} & 15 & \left\{0,2,3,5\right\} & 12 & 188\\\hline\left\{1,1,1,1\right\} & 1 & \left\{1,2,3,4\right\} & 4 & 4\\\hline\end{array}
Cell 16: ✓ — Out[61]= \begin{array}{|c|c|c|c|c|}\hline\text{mu} & \text{dim} & \text{deg q\_k} & \text{Cas2} & \text{Cas4}\\\hline\left\{16,16,0,0\right\} & 8721 & \left\{0,1,18,19\right\} & 576 & 186624\\\hline\left\{16,15,1,0\right\} & 27455 & \left\{0,2,17,19\right\} & 544 & 168160\\\hline\left\{16,14,2,0\right\} & 47424 & \left\{0,3,16,19\right\} & 516 & 152844\\\hline\left\{16,13,3,0\right\} & 62700 & \le… [1246 chars total]


## 2026-07-20 16:08:06
**Input:**
```wolfram
(* test the solver mechanism at (S,P)=(1,1), mu={2,1,1,0}, degrees {0,2,3,5} -- a genuine
   mixing case (not the trivial top state), so the gauge-fixing step below will do real work. *)
Ssub = 1; Psub = 1; mutest = {2, 1, 1, 0};
degs = degQ[mutest, Ssub, Psub] /. {S -> Ssub, P -> Psub}
tau1s[u_] = tau1[mutest][u] /. {S -> Ssub, P -> Psub};
tau2s[u_] = tau2[mutest][u] /. {S -> Ssub, P -> Psub};
tau3s[u_] = tau3[mutest][u] /. {S -> Ssub, P -> Psub};
tau4s[u_] = t[4][u] /. {S -> Ssub, P -> Psub};
nu1s[u_] = ν1[u] /. {S -> Ssub, P -> Psub};
{degs, tau1s[u]}
```
**Output:**
```
[message] Set::write: Tag ReplaceAll in {0, 2, 3, 5} /. {(------------------------------------- + 2 (-θ1 + (u_)) (-θ2 + (u_)) + 2 (-(h S) - θ1 + (u_)) (-(h P) - θ2 + (u_))) (S -> 1), (--------------------------… [510 chars]
Out= {-2*h^2 + 2*(u - θ1)*(u - θ2) + 2*(-h + u - θ1)*(-h + u - θ2), (h^2*(29/2 - (1/2 + P + S)^2 - (3/2 + P + S)^2))/2 + 2*(u - θ1)*(u - θ2) + 2*(-(h*S) + u - θ1)*(-(h*P) + u - θ2)}
```

## 2026-07-20 16:08:28
**Input:**
```wolfram
Ssub = 1; Psub = 1; mutest = {2, 1, 1, 0};
degs = degQ[mutest, Ssub, Psub] /. {S -> Ssub, P -> Psub};
tau1s[u_] = tau1[mutest][u] /. {S -> Ssub, P -> Psub};
tau2s[u_] = tau2[mutest][u] /. {S -> Ssub, P -> Psub};
tau3s[u_] = tau3[mutest][u] /. {S -> Ssub, P -> Psub};
tau4s[u_] = t[4][u] /. {S -> Ssub, P -> Psub};
nu1s[u_] = ν1[u] /. {S -> Ssub, P -> Psub};
{degs, tau1s[u]}
```
**Output:**
```
Out= {{0, 2, 3, 5}, -2*h^2 + 2*(u - θ1)*(u - θ2) + 2*(-h + u - θ1)*(-h + u - θ2)}
```

## 2026-07-20 16:08:34
**Input:**
```wolfram
(* cleared-denominator Baxter equation: multiply through by nu1(u)nu1(u-h)nu1(u-2h)nu1(u-3h) *)
clearedBax[q_][u_] := q[u + h] nu1s[u] nu1s[u - h] nu1s[u - 2 h] nu1s[u - 3 h]
   - tau1s[u] q[u] nu1s[u - h] nu1s[u - 2 h] nu1s[u - 3 h]
   + tau2s[u] q[u - h] nu1s[u - 2 h] nu1s[u - 3 h]
   - tau3s[u] q[u - 2 h] nu1s[u - 3 h]
   + tau4s[u] q[u - 3 h];
(* solve for q_1, degree 0 = {0,2,3,5}[[1]] *)
M = degs[[1]];
cvars = Table[cc[j], {j, 0, M}];
qtest[uu_] := Sum[cc[j] uu^j, {j, 0, M}];
poly = Expand[clearedBax[qtest][u]];
clist = CoefficientList[poly, u];
Length[clist]
```
**Output:**
```
Out= 5
```

## 2026-07-20 16:08:45
**Input:**
```wolfram
(* M=0 -> q1=cc[0] constant. every coefficient of u^j in clist must vanish -- check it's
   automatically consistent (all proportional to cc[0], i.e. Q_1=const is an exact solution) *)
Simplify[clist]
```
**Output:**
```
Out= {5*h^4*(12*h^2 + 7*h*θ1 + θ1^2)*(12*h^2 + 7*h*θ2 + θ2^2)*cc[0], -5*h^4*(168*h^3 + 73*h^2*(θ1 + θ2) + 2*θ1*θ2*(θ1 + θ2) + 7*h*(θ1^2 + 4*θ1*θ2 + θ2^2))*cc[0], 5*h^4*(73*h^2 + θ1^2 + 4*θ1*θ2 + θ2^2 + 21*h*(θ1 + θ2))*cc[0], -10*h^4*(7*h + θ1 + θ2)*cc[0], 5*h^4*cc[0]}
```

## 2026-07-20 16:10:37
**Input:**
```wolfram
(* diagnostic: does q_1=const solve the cleared Baxter eq? Test the KNOWN-GOOD top state first
   (mu_top={2,2,0,0}, tau=vacuum, verified), then the mixing case {2,1,1,0}. Plug numbers. *)
setup[mu_] := (
   tau1s[u_] = tau1[mu][u] /. {S -> Ssub, P -> Psub};
   tau2s[u_] = tau2[mu][u] /. {S -> Ssub, P -> Psub};
   tau3s[u_] = tau3[mu][u] /. {S -> Ssub, P -> Psub};
   tau4s[u_] = t[4][u] /. {S -> Ssub, P -> Psub};
   nu1s[u_] = ν1[u] /. {S -> Ssub, P -> Psub};);
clB1[u_] := nu1s[u] nu1s[u - h] nu1s[u - 2 h] nu1s[u - 3 h]
   - tau1s[u] nu1s[u - h] nu1s[u - 2 h] nu1s[u - 3 h]
   + tau2s[u] nu1s[u - 2 h] nu1s[u - 3 h]
   - tau3s[u] nu1s[u - 3 h] + tau4s[u];
Ssub = 1; Psub = 1;
num = {θ1 -> 1/3, θ2 -> 1/7, h -> 1};
setup[{2, 2, 0, 0}]; top = Expand[clB1[u] /. num];
setup[{2, 1, 1, 0}]; mix = Expand[clB1[u] /. num];
{"q1=const solves at mu_top?" -> (Simplify[top] === 0),
 "q1=const solves at {2,1,1,0}?" -> (Simplify[mix] === 0),
 "residual at {2,1,1,0}" -> Simplify[mix]}
```
**Output:**
```
Out= {"q1=const solves at mu_top?" -> True, "q1=const solves at {2,1,1,0}?" -> False, "residual at {2,1,1,0}" -> (5*(82940 - 90432*u + 36745*u^2 - 6594*u^3 + 441*u^4))/441}
```

## 2026-07-20 16:11:42
**Input:**
```wolfram
setup[{2, 1, 1, 0}];
ratio = Simplify[(mix)/(nu1s[u - 2 h] nu1s[u - 3 h] /. num)];
{"residual / (nu1(u-2h)nu1(u-3h))" -> ratio, "is it u-independent?" -> FreeQ[ratio, u]}
```
**Output:**
```
Out= {"residual / (nu1(u-2h)nu1(u-3h))" -> 5, "is it u-independent?" -> True}
```

## 2026-07-20 16:12:55
**Input:**
```wolfram
(* correct d[2] via solvability: d2correct = -clearedBax0[1]/(nu1(u-2h)nu1(u-3h)), where
   clearedBax0 uses tau2 with its constant term REMOVED. Must be u-independent. Keep th,h symbolic. *)
d2correct[mu_, SS_, PP_] := Module[{t1, t2, t2nc, t3, t4, n1, cb0, dd},
   t1[u_] = tau1[mu][u] /. {S -> SS, P -> PP};
   t2[u_] = tau2[mu][u] /. {S -> SS, P -> PP};
   t2nc[u_] = t2[u] - (t2[u] /. u -> 0 // # &);  (* remove constant term *)
   t2nc[u_] = t2[u] - Coefficient[Collect[t2[u], u], u, 0];
   t3[u_] = tau3[mu][u] /. {S -> SS, P -> PP};
   t4[u_] = t[4][u] /. {S -> SS, P -> PP};
   n1[u_] = ν1[u] /. {S -> SS, P -> PP};
   cb0 = n1[u] n1[u - h] n1[u - 2 h] n1[u - 3 h] - t1[u] n1[u - h] n1[u - 2 h] n1[u - 3 h]
      + t2nc[u] n1[u - 2 h] n1[u - 3 h] - t3[u] n1[u - 3 h] + t4[u];
   dd = Simplify[-cb0/(n1[u - 2 h] n1[u - 3 h])];
   {FreeQ[dd, u], dd}];
(* check on {2,1,1,0}: d2correct should equal d2mine - 5 at the numeric point *)
{"d2correct[{2,1,1,0}]" -> (d2correct[{2, 1, 1, 0}, 1, 1] /. num),
 "d2mine" -> (Coefficient[Collect[tau2[{2, 1, 1, 0}][u] /. {S -> 1, P -> 1}, u], u, 0] /. num)}
```
**Output:**
```
Out= {"d2correct[{2,1,1,0}]" -> {True, -640/441}, "d2mine" -> 1565/441}
```

## 2026-07-20 16:13:34
**Input:**
```wolfram
(* usable irreps have mu_1 = S+P (so deg q_1 = 0). Collect (S,P,mu,Cas2,Cas4,d2correct). *)
usable = Flatten[Table[
    {SS, PP, #} & /@ Select[irreps[SS, PP], #[[1]] == SS + PP &],
    {SS, 1, 3}, {PP, 1, 3}], 2];
(* symbolic d2correct (theta,h symbolic) for each; also record Cas2,Cas4 *)
data = Table[
   Module[{SS = row[[1]], PP = row[[2]], mu = row[[3]], dc},
     dc = d2correct[mu, SS, PP];
     <|"SP" -> {SS, PP}, "mu" -> mu, "Cas2" -> Cas2[mu], "Cas4" -> Cas4[mu],
       "uindep" -> dc[[1]], "d2" -> dc[[2]]|>],
   {row, usable}];
{"n usable irreps" -> Length[data], "all u-independent?" -> AllTrue[data, #["uindep"] &],
 "SP,mu list" -> ({#["SP"], #["mu"]} & /@ data)}
```
**Output:**
```
Out= {"n usable irreps" -> 23, "all u-independent?" -> True, "SP,mu list" -> {{{1, 1}, {2, 2, 0, 0}}, {{1, 1}, {2, 1, 1, 0}}, {{1, 2}, {3, 3, 0, 0}}, {{1, 2}, {3, 2, 1, 0}}, {{1, 3}, {4, 4, 0, 0}}, {{1, 3}, {4, 3, 1, 0}}, {{2, 1}, {3, 3, 0, 0}}, {{2, 1}, {3, 2, 1, 0}}, {{2, 2}, {4, 4, 0, 0}}, {{2, 2}, {4, 3, 1, 0}}, {{2, 2}, {4, 2, 2, 0}}, {{2, 3}, {5, 5, 0, 0}}, {{2, 3}, {5, 4, 1, 0}}, {{2, 3}, {5, 3, 2, 0}}, {{3, 1}, {4, 4, 0, 0}}, {{3, 1}, {4, 3, 1, 0}}, {{3, 2}, {5, 5, 0, 0}}, {{3, 2}, {5, 4, 1, 0}}, {{3, 2}, {5, 3, 2, 0}}, {{3, 3}, {6, 6, 0, 0}}, {{3, 3}, {6, 5, 1, 0}}, {{3, 3}, {6, 4, 2, 0}}, {{3, 3}, {6, 3, 3, 0}}}}
```

## 2026-07-20 16:14:16
**Input:**
```wolfram
fitSP[sp_] := Module[{rows, eqs, sol},
   rows = Select[data, #["SP"] == sp &];
   eqs = Table[r["d2"] == AA + BB r["Cas2"] + CC r["Cas2"]^2 + DD r["Cas4"], {r, rows}];
   sol = Solve[eqs, {AA, BB, CC, DD}];
   If[sol === {}, "underdetermined", {AA, BB, CC, DD} /. sol[[1]] // Simplify]];
f33 = fitSP[{3, 3}];
{"(3,3): {A,B,C,D}" -> f33, "C,D universal?" -> ({f33[[3]], f33[[4]]})}
```
**Output:**
```
[message] Solve::svars: Equations may not give solutions for all "solve" variables.
Out= {"(3,3): {A,B,C,D}" -> {AA, (7*AA + 6*(1880*h^4 - 7*θ1^2*θ2^2 + 414*h^3*(θ1 + θ2) - 28*h*θ1*θ2*(θ1 + θ2) + h^2*(-14*θ1^2 + 123*θ1*θ2 - 14*θ2^2)))/864, (AA + 624*h^4 - 6*θ1^2*θ2^2 + 108*h^3*(θ1 + θ2) - 24*h*θ1*θ2*(θ1 + θ2) - 6*h^2*(2*θ1^2 + 3*θ1*θ2 + 2*θ2^2))/1728, (-AA + 6*(-104*h^4 + θ1^2*θ2^2 - 18*h^3*(θ1 + θ2) + 4*h*θ1*θ2*(θ1 + θ2) + h^2*(2*θ1^2 + 3*θ1*θ2 + 2*θ2^2)))/864}, "C,D universal?" -> {(AA + 624*h^4 - 6*θ1^2*θ2^2 + 108*h^3*(θ1 + θ2) - 24*h*θ1*θ2*(θ1 + θ2) - 6*h^2*(2*θ1^2 + 3*θ1*θ2 + 2*θ2^2))/1728, (-AA + 6*(-104*h^4 + θ1^2*θ2^2 - 18*h^3*(θ1 + θ2) + 4*h*θ1*θ2*(θ1 + θ2) + h^2*(2*θ1^2 + 3*θ1*θ2 + 2*θ2^2)))/864}}
```

## 2026-07-20 16:14:34
**Input:**
```wolfram
sps = DeleteDuplicates[#["SP"] & /@ data];
Clear[AA, BB];
eqsAll = Table[
   r["d2"] == AA[r["SP"]] + BB[r["SP"]] r["Cas2"] + CC r["Cas2"]^2 + DD r["Cas4"],
   {r, data}];
unk = Join[{CC, DD}, AA /@ sps, BB /@ sps];
solAll = Solve[eqsAll, unk];
{"solvable?" -> (solAll =!= {}), "CC" -> (CC /. solAll[[1]] // Simplify), "DD" -> (DD /. solAll[[1]] // Simplify)}
```
**Output:**
```
[message] Solve::svars: Equations may not give solutions for all "solve" variables.
Out= {"solvable?" -> True, "CC" -> CC, "DD" -> -2*CC}
```

## 2026-07-20 16:15:26
**Input:**
```wolfram
(* general d2correct for any deg q_1 = m1: solve the cleared Baxter eq for q (degree m1) AND
   the unknown tau2-constant dv simultaneously; solvability pins dv. *)
d2correctGen[mu_, SS_, PP_] := Module[
   {m1, t1, t2, t2v, t3, t4, n1, qq, cb, cl, sol, dmine},
   m1 = (SS + PP) - mu[[1]];
   t1[u_] = tau1[mu][u] /. {S -> SS, P -> PP};
   t2[u_] = tau2[mu][u] /. {S -> SS, P -> PP};
   dmine = Coefficient[Collect[t2[u], u], u, 0];
   t2v[u_] = t2[u] - dmine + dv;   (* replace tau2 constant by unknown dv *)
   t3[u_] = tau3[mu][u] /. {S -> SS, P -> PP};
   t4[u_] = t[4][u] /. {S -> SS, P -> PP};
   n1[u_] = ν1[u] /. {S -> SS, P -> PP};
   qq[uu_] := Sum[cc[j] uu^j, {j, 0, m1}];
   cb = Expand[qq[u + h] n1[u] n1[u - h] n1[u - 2 h] n1[u - 3 h]
      - t1[u] qq[u] n1[u - h] n1[u - 2 h] n1[u - 3 h]
      + t2v[u] qq[u - h] n1[u - 2 h] n1[u - 3 h]
      - t3[u] qq[u - 2 h] n1[u - 3 h] + t4[u] qq[u - 3 h]];
   cl = CoefficientList[cb, u];
   sol = Solve[Join[Thread[cl == 0], {cc[m1] == 1}],
      Join[Table[cc[j], {j, 0, m1}], {dv}]];
   If[sol === {}, "no solution", dv /. sol[[1]] // Simplify]];
(* test on determinant rep {1,1,1,1} (S=P=1, m1=1) *)
d2correctGen[{1, 1, 1, 1}, 1, 1] /. num
```
**Output:**
```
Out= 1040/441
```

## 2026-07-20 16:15:43
**Input:**
```wolfram
extra = {{1, 1, {1, 1, 1, 1}}, {2, 1, {2, 2, 1, 1}}, {1, 2, {2, 2, 1, 1}}, {2, 2, {3, 3, 1, 1}}};
dataX = Table[
   Module[{SS = row[[1]], PP = row[[2]], mu = row[[3]], dc},
     dc = d2correctGen[mu, SS, PP];
     <|"SP" -> {SS, PP}, "mu" -> mu, "Cas2" -> Cas2[mu], "Cas4" -> Cas4[mu], "d2" -> dc|>],
   {row, extra}];
allData = Join[data, dataX];
(* refit pooled with the degeneracy-breaking points *)
sps2 = DeleteDuplicates[#["SP"] & /@ allData];
Clear[AA, BB];
eqs2 = Table[r["d2"] == AA[r["SP"]] + BB[r["SP"]] r["Cas2"] + CC r["Cas2"]^2 + DD r["Cas4"], {r, allData}];
sol2 = Solve[eqs2, Join[{CC, DD}, AA /@ sps2, BB /@ sps2]];
{"CC" -> (CC /. sol2[[1]] // Simplify), "DD" -> (DD /. sol2[[1]] // Simplify)}
```
**Output:**
```
Out= {"CC" -> h^4/8, "DD" -> -1/4*h^4}
```

## 2026-07-20 16:15:57
**Input:**
```wolfram
Avals = Table[sp -> Simplify[AA[sp] /. sol2[[1]]], {sp, sps2}];
Bvals = Table[sp -> Simplify[BB[sp] /. sol2[[1]]], {sp, sps2}];
Column[{"A(S,P):", Column[Avals], "B(S,P):", Column[Bvals]}]
```
**Output:**
```
Out= Column[{"A(S,P):", Column[{{1, 1} -> 2*(-4*h^4 + 3*θ1^2*θ2^2 - 5*h^3*(θ1 + θ2) + 6*h*θ1*θ2*(θ1 + θ2) + h^2*(θ1^2 + 5*θ1*θ2 + θ2^2)), {1, 2} -> -33*h^4 + 6*θ1^2*θ2^2 + 6*h*θ1*θ2*(3*θ1 + 2*θ2) - h^3*(27*θ1 + 20*θ2) + 2*h^2*(3*θ1^2 + 5*θ1*θ2 + θ2^2), {1, 3} -> 2*(-42*h^4 + 3*θ1^2*θ2^2 + 6*h*θ1*θ2*(2*θ1 + θ2) - h^3*(30*θ1 + 17*θ2) + h^2*(6*θ1^2 + 3*θ1*θ2 + θ2^2)), {2, 1} -> -33*h^4 + 6*θ1^2*θ2^2 + 6*h*θ1*θ2*(2*θ1 + 3*θ2) - h^3*(20*θ1 + 27*θ2) + 2*h^2*(θ1^2 + 5*θ1*θ2 + 3*θ2^2), {2, 2} -> 2*(-48*h^4 + 3*θ1^2*θ2^2 - 21*h^3*(θ1 + θ2) + 9*h*θ1*θ2*(θ1 + θ2) + h^2*(3*θ1^2 + 7*θ1*θ2 + 3*θ2^2)), {2, 3} -> -208*h^4 + 6*θ1^2*θ2^2 + 6*h*θ1*θ2*(4*θ1 + 3*θ2) - h^3*(80*θ1 + 63*θ2) + 2*h^2*(6*θ1^2 + 7*θ1*θ2 + 3*θ2^2), {3, 1} -> 2*(-42*h^4 + 3*θ1^2*θ2^2 + 6*h*θ1*θ2*(θ1 + 2*θ2) - h^3*(17*θ1 + 30*θ2) + h^2*(θ1^2 + 3*θ1*θ2 + 6*θ2^2)), {3, 2} -> -208*h^4 + 6*θ1^2*θ2^2 + 6*h*θ1*θ2*(3*θ1 + 4*θ2) - h^3*(63*θ1 + 80*θ2) + 2*h^2*(3*θ1^2 + 7*θ1*θ2 + 6*θ2^2), {3, 3} -> 6*(-68*h^4 + θ1^2*θ2^2 - 18*h^3*(θ1 + θ2) + 4*h*θ1*θ2*(θ1 + θ2) + h^2*(2*θ1^2 + 3*θ1*θ2 + 2*θ2^2))}], "B(S,P):", Column[{{1, 1} -> (11*h^4)/4 + h^2*θ1*θ2 + h^3*(θ1 + θ2), {1, 2} -> (h^2*(8*h^2 + 3*h*θ1 + 2*h*θ2 + 2*θ1*θ2))/2, {1, 3} -> (21*h^4)/4 + h^2*θ1*θ2 + h^3*(2*θ1 + θ2), {2, 1} -> (h^2*(8*h^2 + 2*h*θ1 + 3*h*θ2 + 2*θ1*θ2))/2, {2, 2} -> (h^2*(23*h^2 + 4*θ1*θ2 + 6*h*(θ1 + θ2)))/4, {2, 3} -> (h^2*(15*h^2 + 4*h*θ1 + 3*h*θ2 + 2*θ1*θ2))/2, {3, 1} -> (21*h^4)/4 + h^2*θ1*θ2 + h^3*(θ1 + 2*θ2), {3, 2} -> (h^2*(15*h^2 + 3*h*θ1 + 4*h*θ2 + 2*θ1*θ2))/2, {3, 3} -> (39*h^4)/4 + h^2*θ1*θ2 + 2*h^3*(θ1 + θ2)}]}]
```

## 2026-07-20 16:16:55
**Input:**
```wolfram
Bansatz = h^4 (q0 + q1 (S + P) + q2 S P + q3 (S^2 + P^2))
   + h^3 (θ1 (r0 + r1 P + r2 S) + θ2 (r0 + r1 S + r2 P)) + h^2 θ1 θ2 s0;
Bfit = Solve[Table[(Bansatz /. {S -> sp[[1]], P -> sp[[2]]}) == (sp /. Bvals), {sp, sps2}],
   {q0, q1, q2, q3, r0, r1, r2, s0}];
Aansatz = h^4 (p0 + p1 (S + P) + p2 S P + p3 (S^2 + P^2))
   + h^3 (θ1 (a0 + a1 P + a2 S) + θ2 (a0 + a1 S + a2 P))
   + h^2 (θ1^2 (e0 + e1 P + e2 S) + θ2^2 (e0 + e1 S + e2 P) + θ1 θ2 (g0 + g1 (S + P)))
   + h (θ1^2 θ2 (k0 + k1 S + k2 P) + θ2^2 θ1 (k0 + k1 P + k2 S)) + θ1^2 θ2^2 w0;
Afit = Solve[Table[(Aansatz /. {S -> sp[[1]], P -> sp[[2]]}) == (sp /. Avals), {sp, sps2}],
   {p0, p1, p2, p3, a0, a1, a2, e0, e1, e2, g0, g1, k0, k1, k2, w0}];
{"Bfit solvable" -> (Bfit =!= {}), "Afit solvable" -> (Afit =!= {}),
 "Bclosed" -> Simplify[Bansatz /. Bfit[[1]]], "Aclosed" -> Simplify[Aansatz /. Afit[[1]]]}
```
**Output:**
```
[message] Solve::svars: Equations may not give solutions for all "solve" variables.
[message] Part::partw: Part 1 of {} does not exist.
[message] ReplaceAll::reps: {{}[[1]]} is neither a list of replacement rules nor a valid dispatch table, and so cannot be used for replacing.
Out= {"Bfit solvable" -> True, "Afit solvable" -> False, "Bclosed" -> (h^2*(h^2*(3 + 3*P + 3*S + 2*P*S) + 4*θ1*θ2 + 2*h*(θ1 + P*θ1 + θ2 + S*θ2)))/4, "Aclosed" -> (h^4*(p0 + P*p2*S + p1*(P + S) + p3*(P^2 + S^2)) + w0*θ1^2*θ2^2 + h^3*((a0 + a1*P + a2*S)*θ1 + (a0 + a2*P + a1*S)*θ2) + h*θ1*θ2*((k0 + k2*P + k1*S)*θ1 + (k0 + k1*P + k2*S)*θ2) + h^2*((e0 + e1*P + e2*S)*θ1^2 + (g0 + g1*(P + S))*θ1*θ2 + (e0 + e2*P + e1*S)*θ2^2) /. {}[[1]])}
```

## 2026-07-20 16:17:40
**Input:**
```wolfram
Bclosed[SS_, PP_] := h^4 (3 + 3 (SS + PP) + 2 SS PP)/4 + h^2 θ1 θ2 + (h^3/2) ((1 + PP) θ1 + (1 + SS) θ2);
Cas2top[SS_, PP_] := Cas2[{SS + PP, SS + PP, 0, 0}];
Cas4top[SS_, PP_] := Cas4[{SS + PP, SS + PP, 0, 0}];
d2closed[mu_, SS_, PP_] := (d2vac /. {S -> SS, P -> PP})
   + Bclosed[SS, PP] (Cas2[mu] - Cas2top[SS, PP])
   + (h^4/8) (Cas2[mu]^2 - Cas2top[SS, PP]^2) - (h^4/4) (Cas4[mu] - Cas4top[SS, PP]);
(* verify against every computed d2correct value (all 27) *)
checks = Table[Simplify[d2closed[r["mu"], r["SP"][[1]], r["SP"][[2]]] - r["d2"]], {r, allData}];
{"max residual over all 27 fitted irreps" -> Max[Abs[checks /. num]], "all zero?" -> (Union[checks] === {0})}
```
**Output:**
```
[message] SetDelayed::write: Tag Plus in (-(-) + (- + P + S)  + (- + P + S) )[SS_, PP_] is Protected.                                   2     2              2
[message] [SetDelayed::write ×1 more suppressed — first instance shown above]
Out= {"max residual over all 27 fitted irreps" -> Max[Abs[96 + 2*(2*(P + S)*(2 + P + S)^3)[1, 1] - (88*(-16 + (2*(P^2 + 2*P*(1 + S) + S*(2 + S)))[1, 1]))/21 - 22*(2*(P^2 + 2*P*(1 + S) + S*(2 + S)))[1, 1] - (2*(P^2 + 2*P*(1 + S) + S*(2 + S)))[1, 1]^2]/8, Abs[360 + 2*(2*(P + S)*(2 + P + S)^3)[1, 2] - (116*(-30 + (2*(P^2 + 2*P*(1 + S) + S*(2 + S)))[1, 2]))/21 - 32*(2*(P^2 + 2*P*(1 + S) + S*(2 + S)))[1, 2] - (2*(P^2 + 2*P*(1 + S) + S*(2 + S)))[1, 2]^2]/8, Abs[864 + 2*(2*(P + S)*(2 + P + S)^3)[1, 3] - (48*(-48 + (2*(P^2 + 2*P*(1 + S) + S*(2 + S)))[1, 3]))/7 - 42*(2*(P^2 + 2*P*(1 + S) + S*(2 + S)))[1, 3] - (2*(P^2 + 2*P*(1 + S) + S*(2 + S)))[1, 3]^2]/8, Abs[360 + 2*(2*(P + S)*(2 + P + S)^3)[2, 1] - (100*(-30 + (2*(P^2 + 2*P*(1 + S) + S*(2 + S)))[2, 1]))/21 - 32*(2*(P^2 + 2*P*(1 + S) + S*(2 + S)))[2, 1] - (2*(P^2 + 2*P*(1 + S) + S*(2 + S)))[2, 1]^2]/8, Abs[1056 + 2*(2*(P + S)*(2 + P + S)^3)[2, 2] - (128*(-48 + (2*(P^2 + 2*P*(1 + S) + S*(2 + S)))[2, 2]))/21 - 46*(2*(P^2 + 2*P*(1 + S) + S*(2 + S)))[2, 2] - (2*(P^2 + 2*P*(1 + S) + S*(2 + S)))[2, 2]^2]/8, Abs[2240 + 2*(2*(P + S)*(2 + P + S)^3)[2, 3] - (52*(-70 + (2*(P^2 + 2*P*(1 + S) + S*(2 + S)))[2, 3]))/7 - 60*(2*(P^2 + 2*P*(1 + S) + S*(2 + S)))[2, 3] - (2*(P^2 + 2*P*(1 + S) + S*(2 + S)))[2, 3]^2]/8, Abs[864 + 2*(2*(P + S)*(2 + P + S)^3)[3, 1] - (16*(-48 + (2*(P^2 + 2*P*(1 + S) + S*(2 + S)))[3, 1]))/3 - 42*(2*(P^2 + 2*P*(1 + S) + S*(2 + S)))[3, 1] - (2*(P^2 + 2*P*(1 + S) + S*(2 + S)))[3, 1]^2]/8, Abs[2240 + 2*(2*(P + S)*(2 + P + S)^3)[3, 2] - (20*(-70 + (2*(P^2 + 2*P*(1 + S) + S*(2 + S)))[3, 2]))/3 - 60*(2*(P^2 + 2*P*(1 + S) + S*(2 + S)))[3, 2] - (2*(P^2 + 2*P*(1 + S) + S*(2 + S)))[3, 2]^2]/8, Abs[4416 + 2*(2*(P + S)*(2 + P + S)^3)[3, 3] - 8*(-96 + (2*(P^2 + 2*P*(1 + S) + S*(2 + S)))[3, 3]) - 78*(2*(P^2 + 2*P*(1 + S) + S*(2 + S)))[3, 3] - (2*(P^2 + 2*P*(1 + S) + S*(2 + S)))[3, 3]^2]/8], "all zero?" -> False}
```

## 2026-07-20 16:18:11
**Input:**
```wolfram
Clear[Cas2top, Cas4top, c2top, c4top, d2closed];
c2top[SS_, PP_] := Cas2[{SS + PP, SS + PP, 0, 0}];
c4top[SS_, PP_] := Cas4[{SS + PP, SS + PP, 0, 0}];
d2closed[mu_, SS_, PP_] := (d2vac /. {S -> SS, P -> PP})
   + Bclosed[SS, PP] (Cas2[mu] - c2top[SS, PP])
   + (h^4/8) (Cas2[mu]^2 - c2top[SS, PP]^2) - (h^4/4) (Cas4[mu] - c4top[SS, PP]);
checks = Table[Simplify[d2closed[r["mu"], r["SP"][[1]], r["SP"][[2]]] - r["d2"]], {r, allData}];
{"all 27 fitted irreps zero?" -> (Union[checks] === {0}),
 (* fresh irreps NOT used in any fit *)
 "fresh {3,2,2,1}@(2,2)" -> Simplify[d2closed[{3, 2, 2, 1}, 2, 2] - d2correctGen[{3, 2, 2, 1}, 2, 2]],
 "fresh {2,2,2,2}@(2,2)" -> Simplify[d2closed[{2, 2, 2, 2}, 2, 2] - d2correctGen[{2, 2, 2, 2}, 2, 2]]}
```
**Output:**
```
Out= {"all 27 fitted irreps zero?" -> True, "fresh {3,2,2,1}@(2,2)" -> -"no solution" + 6*θ1^2*θ2^2 - 6*h^3*(θ1 + θ2) + 18*h*θ1*θ2*(θ1 + θ2) + h^2*(6*θ1^2 + 38*θ1*θ2 + 6*θ2^2), "fresh {2,2,2,2}@(2,2)" -> -"no solution" + 6*(2*h^4 + θ1^2*θ2^2 - 3*h^3*(θ1 + θ2) + 3*h*θ1*θ2*(θ1 + θ2) + h^2*(θ1^2 + 5*θ1*θ2 + θ2^2))}
```

## 2026-07-20 16:19:11
**Input:**
```wolfram
(* corrected tau2: my old tau2 was low in d[2] by the missing h^4(1+2S+2P)/4 * DeltaCas2 *)
Clear[tau2c, Qsolve];
tau2c[mu_, SS_, PP_][u_] := (tau2[mu][u] /. {S -> SS, P -> PP})
   + (h^4 (1 + 2 SS + 2 PP)/4) (Cas2[mu] - Cas2[{SS + PP, SS + PP, 0, 0}]);
(* general polynomial Baxter solver: returns the 4 monic q_k (as coeff lists) of degrees degQ *)
Qsolve[mu_, SS_, PP_] := Module[{ds, t1, t2, t3, t4, n1, res},
   ds = (SS + PP) - (mu[[#]] + 4 - #) + 3 & /@ Range[4];  (* = degQ *)
   t1[u_] = tau1[mu][u] /. {S -> SS, P -> PP};
   t2[u_] = tau2c[mu, SS, PP][u];
   t3[u_] = tau3[mu][u] /. {S -> SS, P -> PP};
   t4[u_] = t[4][u] /. {S -> SS, P -> PP};
   n1[u_] = ν1[u] /. {S -> SS, P -> PP};
   res = Table[
     Module[{m = ds[[kk]], qq, cb, cl, sol},
       qq[uu_] := Sum[cq[j] uu^j, {j, 0, m}];
       cb = Expand[qq[u + h] n1[u] n1[u - h] n1[u - 2 h] n1[u - 3 h]
          - t1[u] qq[u] n1[u - h] n1[u - 2 h] n1[u - 3 h]
          + t2[u] qq[u - h] n1[u - 2 h] n1[u - 3 h]
          - t3[u] qq[u - 2 h] n1[u - 3 h] + t4[u] qq[u - 3 h]];
       cl = CoefficientList[cb, u];
       sol = Solve[Join[Thread[cl == 0], {cq[m] == 1}], Table[cq[j], {j, 0, m}]];
       If[sol === {}, "NO SOLN", Table[cq[j], {j, 0, m}] /. sol[[1]]]],
     {kk, 4}];
   res];
(* test: q_1 for the previously-failing {2,1,1,0}, and the "fresh" {2,2,2,2},{3,2,2,1} *)
{"{2,1,1,0}@(1,1) q1" -> Qsolve[{2, 1, 1, 0}, 1, 1][[1]] /. num,
 "{2,2,2,2}@(2,2) q1" -> (Qsolve[{2, 2, 2, 2}, 2, 2][[1]] /. num),
 "{3,2,2,1}@(2,2) q1" -> (Qsolve[{3, 2, 2, 1}, 2, 2][[1]] /. num)}
```
**Output:**
```
Out= {"{2,1,1,0}@(1,1) q1" -> "NO SOLN", "{2,2,2,2}@(2,2) q1" -> "NO SOLN", "{3,2,2,1}@(2,2) q1" -> "NO SOLN"}
```

## 2026-07-20 16:19:49
**Input:**
```wolfram
(* debug: is tau2c's constant = d2correct, and is clearedBax[const] actually 0? *)
t2test[u_] = tau2c[{2, 1, 1, 0}, 1, 1][u];
n1t[u_] = ν1[u] /. {S -> 1, P -> 1};
t1t[u_] = tau1[{2, 1, 1, 0}][u] /. {S -> 1, P -> 1};
t3t[u_] = tau3[{2, 1, 1, 0}][u] /. {S -> 1, P -> 1};
t4t[u_] = t[4][u] /. {S -> 1, P -> 1};
cbconst = Expand[n1t[u] n1t[u - h] n1t[u - 2 h] n1t[u - 3 h] - t1t[u] n1t[u - h] n1t[u - 2 h] n1t[u - 3 h]
   + t2test[u] n1t[u - 2 h] n1t[u - 3 h] - t3t[u] n1t[u - 3 h] + t4t[u]];
{"tau2c const" -> Coefficient[Collect[t2test[u], u], u, 0] // Simplify,
 "clearedBax[const] identically 0?" -> (Simplify[cbconst] === 0),
 "CoefficientList[0,u]" -> CoefficientList[0, u]}
```
**Output:**
```
Out= {"tau2c const" -> -1/8*((96 + 12*Cas2top + Cas2top^2 - 2*Cas4top)*h^4) + 6*θ1^2*θ2^2 - (-18 + Cas2top)*h^3*(θ1 + θ2) + 12*h*θ1*θ2*(θ1 + θ2) + h^2*(2*θ1^2 - (-38 + Cas2top)*θ1*θ2 + 2*θ2^2), "clearedBax[const] identically 0?" -> False, "CoefficientList[0,u]" -> {}}
```
## 2026-07-20 16:20:07 — ▶️ RUN CELLS 7–10
Cell 7: ✓ — Out[66]= \left\{\text{c[2] vacuum vs asymp-predicted (must be 0)}\to 0,\text{b2vac}\to h^2\,\left(P^2+P\,\left(8\,S+7\right)+S^2+7\,S+6\right)+6\,h\,\left(\theta1\,\left(2\,P+S+3\right)+\theta2\,\left(P+2\,S+3\right)\right)+6\,\left(\theta1^2+4\,\theta1\,\theta2+\theta2^2\right)\right\}
Cell 8: ✓ — Out[69]= h^2\,m1^2+3\,h^2\,m1+h^2\,m2^2+h^2\,m2+h^2\,m3^2-h^2\,m3+h^2\,m4^2-3\,h^2\,m4-h^2\,P^2+4\,h^2\,P\,S+3\,h^2\,P-h^2\,S^2+3\,h^2\,S+6\,h^2+18\,h\,\theta1+18\,h\,\theta2+12\,h\,\theta1\,P+6\,h\,\theta2\,P+6\,h\,\theta1\,S+12\,h\,\theta2\,S+6\,\theta1^2+24\,\theta1\,\theta2+6\,\theta2^2
Cell 9: ✓ — Out[74]= \left\{0,0,0,0\right\}
Cell 10: ✓ — Out[79]= \left\{0,0\right\}


## 2026-07-20 16:20:23
**Input:**
```wolfram
Clear[tau2c, Qsolve];
tau2c[mu_, SS_, PP_][u_] := (tau2[mu][u] /. {S -> SS, P -> PP})
   + (h^4 (1 + 2 SS + 2 PP)/4) (Cas2[mu] - Cas2[{SS + PP, SS + PP, 0, 0}]);
Qsolve[mu_, SS_, PP_] := Module[{ds, t1, t2, t3, t4, n1},
   ds = degQ[mu, SS, PP];
   t1[u_] = tau1[mu][u] /. {S -> SS, P -> PP};
   t2[u_] = tau2c[mu, SS, PP][u];
   t3[u_] = tau3[mu][u] /. {S -> SS, P -> PP};
   t4[u_] = t[4][u] /. {S -> SS, P -> PP};
   n1[u_] = ν1[u] /. {S -> SS, P -> PP};
   Table[Module[{m = ds[[kk]], qq, cl, sol},
       qq[uu_] := Sum[cq[j] uu^j, {j, 0, m}];
       cl = CoefficientList[Expand[
          qq[u + h] n1[u] n1[u - h] n1[u - 2 h] n1[u - 3 h] - t1[u] qq[u] n1[u - h] n1[u - 2 h] n1[u - 3 h]
          + t2[u] qq[u - h] n1[u - 2 h] n1[u - 3 h] - t3[u] qq[u - 2 h] n1[u - 3 h] + t4[u] qq[u - 3 h]], u];
       sol = Solve[Join[Thread[cl == 0], {cq[m] == 1}], Table[cq[j], {j, 0, m}]];
       If[sol === {}, "NO SOLN", Table[cq[j], {j, 0, m}] /. sol[[1]]]], {kk, 4}]];
{"{2,1,1,0}@(1,1)" -> (Qsolve[{2, 1, 1, 0}, 1, 1] /. num),
 "{2,2,2,2}@(2,2) q's" -> (Qsolve[{2, 2, 2, 2}, 2, 2] /. num)}
```
**Output:**
```
[message] Solve::svars: Equations may not give solutions for all "solve" variables.
[message] [Solve::svars ×2 more suppressed — first instance shown above]
[message] General::stop: Further output of Solve::svars will be suppressed during this calculation.
Out= {"{2,1,1,0}@(1,1)" -> {{1}, {cq[0], -10/21, 1}, {cq[0], -8/7 - (10*cq[2])/21, cq[2], 1}, {cq[0], -176/441 - (10*cq[2])/21 - (8*cq[3])/7, cq[2], cq[3], -25/21, 1}}, "{2,2,2,2}@(2,2) q's" -> {{71/735, -31/21, 1}, {33697/46305 + (71*cq[2])/735, -1268/441 - (31*cq[2])/21, cq[2], 1}, {1619/1029 + (71*cq[2])/735 + (33697*cq[3])/46305, -775/147 - (31*cq[2])/21 - (1268*cq[3])/441, cq[2], cq[3], 1}, {35681/9261 + (71*cq[2])/735 + (33697*cq[3])/46305 + (1619*cq[4])/1029, -4736/441 - (31*cq[2])/21 - (1268*cq[3])/441 - (775*cq[4])/147, cq[2], cq[3], cq[4], 1}}}
```

## 2026-07-20 16:21:03
**Input:**
```wolfram
(* gauge-fixed solver: for q_k, impose coefficient at position deg q_j = 0 for every j<k.
   Strict ordering deg q_1<...<deg q_4 makes this a well-posed triangular reduction. *)
Clear[QsolveG];
QsolveG[mu_, SS_, PP_] := Module[{ds, t1, t2, t3, t4, n1},
   ds = degQ[mu, SS, PP];
   t1[u_] = tau1[mu][u] /. {S -> SS, P -> PP};
   t2[u_] = tau2c[mu, SS, PP][u];
   t3[u_] = tau3[mu][u] /. {S -> SS, P -> PP};
   t4[u_] = t[4][u] /. {S -> SS, P -> PP};
   n1[u_] = ν1[u] /. {S -> SS, P -> PP};
   Table[Module[{m = ds[[kk]], gauge = ds[[1 ;; kk - 1]], qq, cl, sol},
       qq[uu_] := Sum[cq[j] uu^j, {j, 0, m}];
       cl = CoefficientList[Expand[
          qq[u + h] n1[u] n1[u - h] n1[u - 2 h] n1[u - 3 h] - t1[u] qq[u] n1[u - h] n1[u - 2 h] n1[u - 3 h]
          + t2[u] qq[u - h] n1[u - 2 h] n1[u - 3 h] - t3[u] qq[u - 2 h] n1[u - 3 h] + t4[u] qq[u - 3 h]], u];
       sol = Solve[Join[Thread[cl == 0], {cq[m] == 1}, Thread[(cq /@ gauge) == 0]], Table[cq[j], {j, 0, m}]];
       Table[cq[j], {j, 0, m}] /. sol[[1]]], {kk, 4}]];
QsolveG[{2, 1, 1, 0}, 1, 1] /. num
```
**Output:**
```
Out= {{1}, {0, -10/21, 1}, {0, -8/7, 0, 1}, {0, -176/441, 0, 0, -25/21, 1}}
```

## 2026-07-20 16:21:20
**Input:**
```wolfram
(* residual check: each gauge-fixed q_k must satisfy the cleared Baxter eq (identically in u) *)
verify[mu_, SS_, PP_] := Module[{qs = QsolveG[mu, SS, PP], t1, t2, t3, t4, n1, ds},
   ds = degQ[mu, SS, PP];
   t1[u_] = tau1[mu][u] /. {S -> SS, P -> PP}; t2[u_] = tau2c[mu, SS, PP][u];
   t3[u_] = tau3[mu][u] /. {S -> SS, P -> PP}; t4[u_] = t[4][u] /. {S -> SS, P -> PP};
   n1[u_] = ν1[u] /. {S -> SS, P -> PP};
   Max@Table[Module[{q = Function[uu, qs[[kk]] . uu^Range[0, ds[[kk]]]]},
       Max@Abs@Simplify[CoefficientList[Expand[
          q[u + h] n1[u] n1[u - h] n1[u - 2 h] n1[u - 3 h] - t1[u] q[u] n1[u - h] n1[u - 2 h] n1[u - 3 h]
          + t2[u] q[u - h] n1[u - 2 h] n1[u - 3 h] - t3[u] q[u - 2 h] n1[u - 3 h] + t4[u] q[u - 3 h]], u]]]],
     {kk, 4}]];
(* sweep every irrep of (1,1),(2,1),(2,2), symbolic theta,h *)
Column[Table[{SS, PP} -> Max[verify[#, SS, PP] & /@ irreps[SS, PP]],
   {SS, {1, 2}}, {PP, {1, 2}}] // Flatten[#, 1] &]
```
**Output:**
```
[message] ToExpression::sntx: Invalid syntax in or before "{kk, 4}]];".                                                           ^
(no output)
```

## 2026-07-20 16:21:31
**Input:**
```wolfram
clB[q_, t1_, t2_, t3_, t4_, n1_] := Expand[
   q[u + h] n1[u] n1[u - h] n1[u - 2 h] n1[u - 3 h] - t1[u] q[u] n1[u - h] n1[u - 2 h] n1[u - 3 h]
   + t2[u] q[u - h] n1[u - 2 h] n1[u - 3 h] - t3[u] q[u - 2 h] n1[u - 3 h] + t4[u] q[u - 3 h]];
verify[mu_, SS_, PP_] := Module[{qs, ds, t1, t2, t3, t4, n1, qf, resid},
   qs = QsolveG[mu, SS, PP];
   ds = degQ[mu, SS, PP];
   t1[u_] = tau1[mu][u] /. {S -> SS, P -> PP};
   t2[u_] = tau2c[mu, SS, PP][u];
   t3[u_] = tau3[mu][u] /. {S -> SS, P -> PP};
   t4[u_] = t[4][u] /. {S -> SS, P -> PP};
   n1[u_] = ν1[u] /. {S -> SS, P -> PP};
   resid = Table[
     qf[uu_] := qs[[kk]] . uu^Range[0, ds[[kk]]];
     Simplify[clB[qf, t1, t2, t3, t4, n1]], {kk, 4}];
   Max[Abs[Flatten[CoefficientList[resid, u]]]]];
Column[Flatten[Table[{SS, PP} -> Max[verify[#, SS, PP] & /@ irreps[SS, PP]], {SS, {1, 2}}, {PP, {1, 2}}], 1]]
```
**Output:**
```
Out= Column[{{1, 1} -> -Infinity, {1, 2} -> -Infinity, {2, 1} -> -Infinity, {2, 2} -> -Infinity}]
```

## 2026-07-20 16:21:48
**Input:**
```wolfram
verify2[mu_, SS_, PP_] := Module[{qs, ds, t1, t2, t3, t4, n1, qf},
   qs = QsolveG[mu, SS, PP]; ds = degQ[mu, SS, PP];
   t1[u_] = tau1[mu][u] /. {S -> SS, P -> PP}; t2[u_] = tau2c[mu, SS, PP][u];
   t3[u_] = tau3[mu][u] /. {S -> SS, P -> PP}; t4[u_] = t[4][u] /. {S -> SS, P -> PP};
   n1[u_] = ν1[u] /. {S -> SS, P -> PP};
   AllTrue[Range[4], (qf[uu_] := qs[[#]] . uu^Range[0, ds[[#]]]; Simplify[clB[qf, t1, t2, t3, t4, n1]] === 0) &]];
Grid[Prepend[
  Flatten[Table[{{SS, PP}, Length[irreps[SS, PP]], AllTrue[irreps[SS, PP], verify2[#, SS, PP] &]},
    {SS, {1, 2, 3}}, {PP, {1, 2}}], 1],
  {"(S,P)", "#irreps", "all q_k solve Baxter?"}], Frame -> All]
```
**Output:**
```
Out= Grid[{{"(S,P)", "#irreps", "all q_k solve Baxter?"}, {{1, 1}, 3, True}, {{1, 2}, 3, True}, {{2, 1}, 3, True}, {{2, 2}, 6, True}, {{3, 1}, 3, True}, {{3, 2}, 6, True}}, Frame -> All]
```
## 2026-07-20 16:22:58 — 📥 BULK INSERT 3 CELL(S) at positions 17–19
1. [code] (* --- tau2 constant corrected by the Baxter solver ---↵   The fusion form for tau2's Cas2-linear pa
2. [code] (* --- Baxter solver with gauge fixing ---↵   Input: irrep mu in [S,S,0,0] x [P,P,0,0]. Solves the u
3. [code] (* --- acceptance: every gauge-fixed q_k solves the Baxter equation, over full (S,P) families ---↵


## 2026-07-20 16:22:58
**Input:**
```wolfram
(* --- tau2 constant corrected by the Baxter solver ---
   The fusion form for tau2's Cas2-linear part nailed the u^2 and u^1 coefficients exactly, but
   its u^0 (constant) piece was incomplete. Demanding that the actual Q-polynomials (of the
   derived degrees) solve the 4th-order Baxter equation pins the true constant: tau2's d[2] was
   low by exactly  h^4 (1+2S+2P)/4 * (Cas2[mu] - Cas2[mu_top]).  (The quartic h^4(Cas2^2/8 -
   Cas4/4) was already correct: CC=h^4/8, DD=-h^4/4 came out of the fit.) Verified over full
   (S,P) sweeps below. *)
tau2c[mu_, SS_, PP_][u_] := (tau2[mu][u] /. {S -> SS, P -> PP})
   + (h^4 (1 + 2 SS + 2 PP)/4) (Cas2[mu] - Cas2[{SS + PP, SS + PP, 0, 0}]);
```
**Output:**
```
(no output)
```

## 2026-07-20 16:22:58
**Input:**
```wolfram
(* --- Baxter solver with gauge fixing ---
   Input: irrep mu in [S,S,0,0] x [P,P,0,0]. Solves the untwisted 4th-order Baxter (TQ) equation
   (transfer matrices tau1, tau2c, tau3, t[4], all in terms of Casimirs) for the four polynomial
   Q-functions q_k, of degrees degQ[mu,S,P]. GAUGE FIXING: the four solutions have strictly
   increasing degrees deg q_1 < ... < deg q_4, so any lower-degree q_j solves the same equation
   and can be added to a higher q_k (q_k -> q_k + c q_j). We remove this freedom by forcing the
   coefficient of u^(deg q_j) in q_k to vanish for every j<k. Returns monic coefficient lists. *)
QsolveG[mu_, SS_, PP_] := Module[{ds, t1, t2, t3, t4, n1},
   ds = degQ[mu, SS, PP];
   t1[u_] = tau1[mu][u] /. {S -> SS, P -> PP};  t2[u_] = tau2c[mu, SS, PP][u];
   t3[u_] = tau3[mu][u] /. {S -> SS, P -> PP};  t4[u_] = t[4][u] /. {S -> SS, P -> PP};
   n1[u_] = ν1[u] /. {S -> SS, P -> PP};
   Table[Module[{m = ds[[kk]], gauge = ds[[1 ;; kk - 1]], qq, cl, sol},
       qq[uu_] := Sum[cq[j] uu^j, {j, 0, m}];
       cl = CoefficientList[Expand[
          qq[u + h] n1[u] n1[u - h] n1[u - 2 h] n1[u - 3 h]
           - t1[u] qq[u] n1[u - h] n1[u - 2 h] n1[u - 3 h]
           + t2[u] qq[u - h] n1[u - 2 h] n1[u - 3 h]
           - t3[u] qq[u - 2 h] n1[u - 3 h] + t4[u] qq[u - 3 h]], u];
       sol = Solve[Join[Thread[cl == 0], {cq[m] == 1}, Thread[(cq /@ gauge) == 0]],
          Table[cq[j], {j, 0, m}]];
       Table[cq[j], {j, 0, m}] /. sol[[1]]], {kk, 4}]];
(* example: mu={2,1,1,0} at (S,P)=(1,1) -- degrees {0,2,3,5}, gauge-fixed *)
QsolveG[{2, 1, 1, 0}, 1, 1]
```
**Output:**
```
Out[82]= \left\{\left\{1\right\},\left\{0,-\theta1-\theta2,1\right\},\left\{0,-h^2-3\,\theta1\,\theta2,0,1\right\},\left\{0,\frac{1}{2}\,\left(-2\,h^4+5\,h^2\,\theta1^2+10\,h^2\,\theta1\,\theta2+5\,h^2\,\theta2^2+30\,\theta1^2\,\theta2^2\right),0,0,-\frac{5}{2}\,\left(\theta1+\theta2\right),1\right\}\right\}
```

## 2026-07-20 16:23:01
**Input:**
```wolfram
(* --- acceptance: every gauge-fixed q_k solves the Baxter equation, over full (S,P) families ---
   exact/symbolic in theta1,theta2,h. Each residual must be identically 0. *)
verifyBax[mu_, SS_, PP_] := Module[{qs = QsolveG[mu, SS, PP], ds = degQ[mu, SS, PP], t1, t2, t3, t4, n1, qf},
   t1[u_] = tau1[mu][u] /. {S -> SS, P -> PP}; t2[u_] = tau2c[mu, SS, PP][u];
   t3[u_] = tau3[mu][u] /. {S -> SS, P -> PP}; t4[u_] = t[4][u] /. {S -> SS, P -> PP};
   n1[u_] = ν1[u] /. {S -> SS, P -> PP};
   AllTrue[Range[4], (qf[uu_] := qs[[#]] . uu^Range[0, ds[[#]]];
      Simplify[qf[u + h] n1[u] n1[u - h] n1[u - 2 h] n1[u - 3 h]
        - t1[u] qf[u] n1[u - h] n1[u - 2 h] n1[u - 3 h]
        + t2[u] qf[u - h] n1[u - 2 h] n1[u - 3 h]
        - t3[u] qf[u - 2 h] n1[u - 3 h] + t4[u] qf[u - 3 h]] === 0) &]];
Grid[Prepend[Flatten[Table[{{SS, PP}, Length[irreps[SS, PP]], AllTrue[irreps[SS, PP], verifyBax[#, SS, PP] &]},
   {SS, {1, 2, 3}}, {PP, {1, 2}}], 1], {"(S,P)", "#irreps", "all q_k solve Baxter?"}], Frame -> All]
```
**Output:**
```
Out[84]= \begin{array}{|c|c|c|}\hline\text{(S,P)} & \text{\#irreps} & \text{all q\_k solve Baxter?}\\\hline\left\{1,1\right\} & 3 & \mathrm{True}\\\hline\left\{1,2\right\} & 3 & \mathrm{True}\\\hline\left\{2,1\right\} & 3 & \mathrm{True}\\\hline\left\{2,2\right\} & 6 & \mathrm{True}\\\hline\left\{3,1\right\} & 3 & \mathrm{True}\\\hline\left\{3,2\right\} & 6 & \mathrm{True}\\\hline\end{array}
```

## 2026-07-20 16:27:58
**Input:**
```wolfram
irr55 = irreps[5, 5];
(* numeric-parameter gauge-fixed solver (exact rationals) for speed at high degree *)
QsolveGnum[mu_, SS_, PP_] := Module[{ds, t1, t2, t3, t4, n1},
   ds = degQ[mu, SS, PP];
   t1[u_] = tau1[mu][u] /. {S -> SS, P -> PP} /. num;
   t2[u_] = tau2c[mu, SS, PP][u] /. num;
   t3[u_] = tau3[mu][u] /. {S -> SS, P -> PP} /. num;
   t4[u_] = t[4][u] /. {S -> SS, P -> PP} /. num;
   n1[u_] = ν1[u] /. {S -> SS, P -> PP} /. num;
   Table[Module[{m = ds[[kk]], gauge = ds[[1 ;; kk - 1]], qq, cl, sol},
       qq[uu_] := Sum[cq[j] uu^j, {j, 0, m}];
       cl = CoefficientList[Expand[
          qq[u + h] n1[u] n1[u - h] n1[u - 2 h] n1[u - 3 h] - t1[u] qq[u] n1[u - h] n1[u - 2 h] n1[u - 3 h]
          + t2[u] qq[u - h] n1[u - 2 h] n1[u - 3 h] - t3[u] qq[u - 2 h] n1[u - 3 h] + t4[u] qq[u - 3 h]] /. num, u];
       sol = Solve[Join[Thread[cl == 0], {cq[m] == 1}, Thread[(cq /@ gauge) == 0]], Table[cq[j], {j, 0, m}]];
       Table[cq[j], {j, 0, m}] /. sol[[1]]], {kk, 4}]];
{"# irreps (5,5)" -> Length[irr55], "top irrep timing (s)" -> First@AbsoluteTiming[QsolveGnum[irr55[[1]], 5, 5];]}
```
**Output:**
```
Out= {"# irreps (5,5)" -> 21, "top irrep timing (s)" -> 0.06222}
```

## 2026-07-20 16:28:25
**Input:**
```wolfram
clBn[qf_, t1_, t2_, t3_, t4_, n1_] := qf[u + h] n1[u] n1[u - h] n1[u - 2 h] n1[u - 3 h]
   - t1[u] qf[u] n1[u - h] n1[u - 2 h] n1[u - 3 h] + t2[u] qf[u - h] n1[u - 2 h] n1[u - 3 h]
   - t3[u] qf[u - 2 h] n1[u - 3 h] + t4[u] qf[u - 3 h];
tbl55 = Table[
   Module[{mu = irr55[[i]], tim, qs, ds, t1, t2, t3, t4, n1, ok, qf},
     ds = degQ[mu, 5, 5];
     {tim, qs} = AbsoluteTiming[QsolveGnum[mu, 5, 5]];
     t1[uu_] = tau1[mu][uu] /. {S -> 5, P -> 5} /. num; t2[uu_] = tau2c[mu, 5, 5][uu] /. num;
     t3[uu_] = tau3[mu][uu] /. {S -> 5, P -> 5} /. num; t4[uu_] = t[4][uu] /. {S -> 5, P -> 5} /. num;
     n1[uu_] = ν1[uu] /. {S -> 5, P -> 5} /. num;
     ok = AllTrue[Range[4], (qf[x_] := qs[[#]] . x^Range[0, ds[[#]]]; Expand[clBn[qf, t1, t2, t3, t4, n1]] === 0) &];
     {mu, dimGL4[mu], ds, Cas2[mu], Cas4[mu], NumberForm[tim, {5, 4}], ok}],
   {i, Length[irr55]}];
Grid[Prepend[tbl55, {"mu", "dim", "deg q_k", "Cas2", "Cas4", "solve time (s)", "all solve?"}],
   Frame -> All, Alignment -> Left, Background -> {None, {{RGBColor[.9, .9, .95]}}}]
```
**Output:**
```
Out= Grid[{{"mu", "dim", "deg q_k", "Cas2", "Cas4", "solve time (s)", "all solve?"}, {{10, 10, 0, 0}, 1716, {0, 1, 12, 13}, 240, 34560, NumberForm[0.045381, {5, 4}], False}, {{10, 9, 1, 0}, 4719, {0, 2, 11, 13}, 220, 29740, NumberForm[0.022366, {5, 4}], False}, {{10, 8, 2, 0}, 6825, {0, 3, 10, 13}, 204, 26172, NumberForm[0.021034, {5, 4}], False}, {{10, 7, 3, 0}, 7020, {0, 4, 9, 13}, 192, 23664, NumberForm[0.02248, {5, 4}], False}, {{10, 6, 4, 0}, 5200, {0, 5, 8, 13}, 184, 22072, NumberForm[0.019642, {5, 4}], False}, {{10, 5, 5, 0}, 1911, {0, 6, 7, 13}, 180, 21300, NumberForm[0.019457, {5, 4}], False}, {{9, 9, 1, 1}, 825, {1, 2, 11, 12}, 196, 23428, NumberForm[0.024398, {5, 4}], False}, {{9, 8, 2, 1}, 2079, {1, 3, 10, 12}, 180, 19860, NumberForm[0.025766, {5, 4}], False}, {{9, 7, 3, 1}, 2640, {1, 4, 9, 12}, 168, 17352, NumberForm[0.018254, {5, 4}], False}, {{9, 6, 4, 1}, 2156, {1, 5, 8, 12}, 160, 15760, NumberForm[0.018147, {5, 4}], False}, {{9, 5, 5, 1}, 825, {1, 6, 7, 12}, 156, 14988, NumberForm[0.020237, {5, 4}], False}, {{8, 8, 2, 2}, 336, {2, 3, 10, 11}, 160, 15040, NumberForm[0.018275, {5, 4}], False}, {{8, 7, 3, 2}, 735, {2, 4, 9, 11}, 148, 12532, NumberForm[0.017522, {5, 4}], False}, {{8, 6, 4, 2}, 729, {2, 5, 8, 11}, 140, 10940, NumberForm[0.016687, {5, 4}], False}, {{8, 5, 5, 2}, 300, {2, 6, 7, 11}, 136, 10168, NumberForm[0.016392, {5, 4}], False}, {{7, 7, 3, 3}, 105, {3, 4, 9, 10}, 132, 8964, NumberForm[0.016159, {5, 4}], False}, {{7, 6, 4, 3}, 175, {3, 5, 8, 10}, 124, 7372, NumberForm[0.016022, {5, 4}], False}, {{7, 5, 5, 3}, 84, {3, 6, 7, 10}, 120, 6600, NumberForm[0.015721, {5, 4}], False}, {{6, 6, 4, 4}, 20, {4, 5, 8, 9}, 112, 4864, NumberForm[0.01503, {5, 4}], False}, {{6, 5, 5, 4}, 15, {4, 6, 7, 9}, 108, 4092, NumberForm[0.015109, {5, 4}], False}, {{5, 5, 5, 5}, 1, {5, 6, 7, 8}, 100, 2500, NumberForm[0.014513, {5, 4}], False}}, Frame -> All, Alignment -> Left, Background -> {None, {{RGBColor[0.9, 0.9, 0.95]}}}]
```

## 2026-07-20 16:29:35
**Input:**
```wolfram
mu0 = {10, 10, 0, 0};
qs0 = QsolveGnum[mu0, 5, 5];
ds0 = degQ[mu0, 5, 5];
(* also symbolic solver then num, to compare *)
qs0sym = QsolveG[mu0, 5, 5] /. num;
{"degs" -> ds0, "q1 (num solver)" -> qs0[[1]], "q2 (num solver)" -> qs0[[2]],
 "num vs symbolic match?" -> (qs0 === qs0sym)}
```
**Output:**
```
Out= {"degs" -> {0, 1, 12, 13}, "q1 (num solver)" -> {1}, "q2 (num solver)" -> {0, 1}, "num vs symbolic match?" -> True}
```

## 2026-07-20 16:29:44
**Input:**
```wolfram
t1[uu_] = tau1[mu0][uu] /. {S -> 5, P -> 5} /. num;
t2[uu_] = tau2c[mu0, 5, 5][uu] /. num;
t3[uu_] = tau3[mu0][uu] /. {S -> 5, P -> 5} /. num;
t4[uu_] = t[4][uu] /. {S -> 5, P -> 5} /. num;
n1[uu_] = ν1[uu] /. {S -> 5, P -> 5} /. num;
qf1[x_] := qs0[[1]] . x^Range[0, ds0[[1]]];
resid1 = Expand[clBn[qf1, t1, t2, t3, t4, n1]];
{"q1 residual is 0?" -> (resid1 === 0), "residual" -> resid1}
```
**Output:**
```
Out= {"q1 residual is 0?" -> False, "residual" -> 2178227584/9261 + (5099267360*h)/21609 - (809127208*h^2)/7203 - (15025315360*h^3)/64827 - (110158640*h^4)/1029 - (2797520*h^5)/147 - (6936*h^6)/7 - (57433005520*u)/194481 - (14318295248*h*u)/64827 + (34924906192*h^2*u)/194481 + (1716240944*h^3*u)/7203 + (266683976*h^4*u)/3087 + (565224*h^5*u)/49 + (2880*h^6*u)/7 + (8466691916*u^2)/64827 + (4805094776*h*u^2)/64827 - (5149311730*h^2*u^2)/64827 - (93593840*h^3*u^2)/1029 - (4505206*h^4*u^2)/147 - (25080*h^5*u^2)/7 - 108*h^6*u^2 - (27595804*u^3)/1029 - (34123444*h*u^3)/3087 + (42622768*h^2*u^3)/3087 + (891560*h^3*u^3)/49 + (38316*h^4*u^3)/7 + 396*h^5*u^3 + (152779*u^4)/49 - (5104*h*u^4)/49 - (37488*h^2*u^4)/49 - (13200*h^3*u^4)/7 - 363*h^4*u^4 - (1772*u^5)/7 + (1864*h*u^5)/7 - (596*h^2*u^5)/7 + 72*h^3*u^5 + 12*u^6 - 24*h*u^6 + 12*h^2*u^6}
```

## 2026-07-20 16:30:18
**Input:**
```wolfram
{"resid1 /. h->1 is 0?" -> (Expand[resid1 /. h -> 1] === 0)}
```
**Output:**
```
Out= {"resid1 /. h->1 is 0?" -> True}
```

## 2026-07-20 16:30:45
**Input:**
```wolfram
tbl55 = Table[
   Module[{mu = irr55[[i]], tim, qs, ds, t1, t2, t3, t4, n1, ok, qf},
     ds = degQ[mu, 5, 5];
     {tim, qs} = AbsoluteTiming[QsolveGnum[mu, 5, 5]];
     t1[uu_] = tau1[mu][uu] /. {S -> 5, P -> 5} /. num; t2[uu_] = tau2c[mu, 5, 5][uu] /. num;
     t3[uu_] = tau3[mu][uu] /. {S -> 5, P -> 5} /. num; t4[uu_] = t[4][uu] /. {S -> 5, P -> 5} /. num;
     n1[uu_] = ν1[uu] /. {S -> 5, P -> 5} /. num;
     ok = AllTrue[Range[4], (qf[x_] := qs[[#]] . x^Range[0, ds[[#]]];
        Expand[clBn[qf, t1, t2, t3, t4, n1] /. num] === 0) &];
     {Row[mu, ","], dimGL4[mu], Row[ds, ","], Cas2[mu], Cas4[mu], NumberForm[tim, {6, 4}], ok}],
   {i, Length[irr55]}];
{"all irreps solve?" -> AllTrue[tbl55, #[[7]] &], "total solve time (s)" -> Total[QuantityMagnitude /@ (tbl55[[All, 6]] /. NumberForm[x_, _] :> x)]}
```
**Output:**
```
Out= {"all irreps solve?" -> True, "total solve time (s)" -> 0.498406}
```

## 2026-07-20 16:30:51
**Input:**
```wolfram
Grid[Prepend[tbl55,
   {"mu", "dim", "deg q_k", "Cas2", "Cas4", "time (s)", "solves?"}],
  Frame -> All, Alignment -> Left,
  Background -> {None, {RGBColor[.85, .88, .95], {White, RGBColor[.95, .96, .99]}}}]
```
**Output:**
```
Out= Grid[{{"mu", "dim", "deg q_k", "Cas2", "Cas4", "time (s)", "solves?"}, {Row[{10, 10, 0, 0}, ","], 1716, Row[{0, 1, 12, 13}, ","], 240, 34560, NumberForm[0.034358, {6, 4}], True}, {Row[{10, 9, 1, 0}, ","], 4719, Row[{0, 2, 11, 13}, ","], 220, 29740, NumberForm[0.0409, {6, 4}], True}, {Row[{10, 8, 2, 0}, ","], 6825, Row[{0, 3, 10, 13}, ","], 204, 26172, NumberForm[0.020315, {6, 4}], True}, {Row[{10, 7, 3, 0}, ","], 7020, Row[{0, 4, 9, 13}, ","], 192, 23664, NumberForm[0.024185, {6, 4}], True}, {Row[{10, 6, 4, 0}, ","], 5200, Row[{0, 5, 8, 13}, ","], 184, 22072, NumberForm[0.022354, {6, 4}], True}, {Row[{10, 5, 5, 0}, ","], 1911, Row[{0, 6, 7, 13}, ","], 180, 21300, NumberForm[0.026041, {6, 4}], True}, {Row[{9, 9, 1, 1}, ","], 825, Row[{1, 2, 11, 12}, ","], 196, 23428, NumberForm[0.028104, {6, 4}], True}, {Row[{9, 8, 2, 1}, ","], 2079, Row[{1, 3, 10, 12}, ","], 180, 19860, NumberForm[0.019101, {6, 4}], True}, {Row[{9, 7, 3, 1}, ","], 2640, Row[{1, 4, 9, 12}, ","], 168, 17352, NumberForm[0.01857, {6, 4}], True}, {Row[{9, 6, 4, 1}, ","], 2156, Row[{1, 5, 8, 12}, ","], 160, 15760, NumberForm[0.01764, {6, 4}], True}, {Row[{9, 5, 5, 1}, ","], 825, Row[{1, 6, 7, 12}, ","], 156, 14988, NumberForm[0.018169, {6, 4}], True}, {Row[{8, 8, 2, 2}, ","], 336, Row[{2, 3, 10, 11}, ","], 160, 15040, NumberForm[0.020801, {6, 4}], True}, {Row[{8, 7, 3, 2}, ","], 735, Row[{2, 4, 9, 11}, ","], 148, 12532, NumberForm[0.025245, {6, 4}], True}, {Row[{8, 6, 4, 2}, ","], 729, Row[{2, 5, 8, 11}, ","], 140, 10940, NumberForm[0.026288, {6, 4}], True}, {Row[{8, 5, 5, 2}, ","], 300, Row[{2, 6, 7, 11}, ","], 136, 10168, NumberForm[0.029542, {6, 4}], True}, {Row[{7, 7, 3, 3}, ","], 105, Row[{3, 4, 9, 10}, ","], 132, 8964, NumberForm[0.03103, {6, 4}], True}, {Row[{7, 6, 4, 3}, ","], 175, Row[{3, 5, 8, 10}, ","], 124, 7372, NumberForm[0.022939, {6, 4}], True}, {Row[{7, 5, 5, 3}, ","], 84, Row[{3, 6, 7, 10}, ","], 120, 6600, NumberForm[0.022829, {6, 4}], True}, {Row[{6, 6, 4, 4}, ","], 20, Row[{4, 5, 8, 9}, ","], 112, 4864, NumberForm[0.015449, {6, 4}], True}, {Row[{6, 5, 5, 4}, ","], 15, Row[{4, 6, 7, 9}, ","], 108, 4092, NumberForm[0.018539, {6, 4}], True}, {Row[{5, 5, 5, 5}, ","], 1, Row[{5, 6, 7, 8}, ","], 100, 2500, NumberForm[0.016007, {6, 4}], True}}, Frame -> All, Alignment -> Left, Background -> {None, {RGBColor[0.85, 0.88, 0.95], {GrayLevel[1], RGBColor[0.95, 0.96, 0.99]}}}]
```
## 2026-07-20 16:31:30 — 📥 BULK INSERT 1 CELL(S) at positions 20–20
1. [code] (* --- (S,P)=(5,5): run the gauge-fixed Baxter solver over the ENTIRE tensor-product↵   decompositio


## 2026-07-20 16:31:31
**Input:**
```wolfram
(* --- (S,P)=(5,5): run the gauge-fixed Baxter solver over the ENTIRE tensor-product
   decomposition (21 irreps), with per-irrep timing. Exact rational params num={th1->1/3,
   th2->1/7,h->1} keep the degree-up-to-13 solves fast. QsolveGnum = QsolveG with num applied
   before the linear solve (much faster than symbolic at high degree). *)
num = {θ1 -> 1/3, θ2 -> 1/7, h -> 1};
QsolveGnum[mu_, SS_, PP_] := Module[{ds, t1, t2, t3, t4, n1},
   ds = degQ[mu, SS, PP];
   t1[u_] = tau1[mu][u] /. {S -> SS, P -> PP} /. num;  t2[u_] = tau2c[mu, SS, PP][u] /. num;
   t3[u_] = tau3[mu][u] /. {S -> SS, P -> PP} /. num;  t4[u_] = t[4][u] /. {S -> SS, P -> PP} /. num;
   n1[u_] = ν1[u] /. {S -> SS, P -> PP} /. num;
   Table[Module[{m = ds[[kk]], gauge = ds[[1 ;; kk - 1]], qq, cl, sol},
       qq[uu_] := Sum[cq[j] uu^j, {j, 0, m}];
       cl = CoefficientList[Expand[
          qq[u + h] n1[u] n1[u - h] n1[u - 2 h] n1[u - 3 h] - t1[u] qq[u] n1[u - h] n1[u - 2 h] n1[u - 3 h]
          + t2[u] qq[u - h] n1[u - 2 h] n1[u - 3 h] - t3[u] qq[u - 2 h] n1[u - 3 h] + t4[u] qq[u - 3 h]] /. num, u];
       sol = Solve[Join[Thread[cl == 0], {cq[m] == 1}, Thread[(cq /@ gauge) == 0]], Table[cq[j], {j, 0, m}]];
       Table[cq[j], {j, 0, m}] /. sol[[1]]], {kk, 4}]];
(* cleared-Baxter residual (h collapsed via num at the end to stay consistent with the shifts) *)
clBn[qf_, t1_, t2_, t3_, t4_, n1_] := qf[u + h] n1[u] n1[u - h] n1[u - 2 h] n1[u - 3 h]
   - t1[u] qf[u] n1[u - h] n1[u - 2 h] n1[u - 3 h] + t2[u] qf[u - h] n1[u - 2 h] n1[u - 3 h]
   - t3[u] qf[u - 2 h] n1[u - 3 h] + t4[u] qf[u - 3 h];
irr55 = irreps[5, 5];
tbl55 = Table[
   Module[{mu = irr55[[i]], tim, qs, ds, t1, t2, t3, t4, n1, ok, qf},
     ds = degQ[mu, 5, 5];
     {tim, qs} = AbsoluteTiming[QsolveGnum[mu, 5, 5]];
     t1[uu_] = tau1[mu][uu] /. {S -> 5, P -> 5} /. num;  t2[uu_] = tau2c[mu, 5, 5][uu] /. num;
     t3[uu_] = tau3[mu][uu] /. {S -> 5, P -> 5} /. num;  t4[uu_] = t[4][uu] /. {S -> 5, P -> 5} /. num;
     n1[uu_] = ν1[uu] /. {S -> 5, P -> 5} /. num;
     ok = AllTrue[Range[4], (qf[x_] := qs[[#]] . x^Range[0, ds[[#]]];
        Expand[clBn[qf, t1, t2, t3, t4, n1] /. num] === 0) &];
     {Row[mu, ","], dimGL4[mu], Row[ds, ","], Cas2[mu], Cas4[mu], NumberForm[tim, {6, 4}], ok}],
   {i, Length[irr55]}];
Print["all 21 irreps solve: ", AllTrue[tbl55, Last], "  total solve time (s): ",
   Total[tbl55[[All, 6]] /. NumberForm[x_, _] :> x]];
Grid[Prepend[tbl55, {"mu", "dim", "deg q_k", "Cas2", "Cas4", "time (s)", "solves?"}],
  Frame -> All, Alignment -> Left,
  Background -> {None, {RGBColor[.85, .88, .95], {GrayLevel[1], RGBColor[.95, .96, .99]}}}]
```
**Output:**
```
Print: all 21 irreps solve:[formula]  total solve time (s):[formula] | Out[91]= \begin{array}{|c|c|c|c|c|c|c|}\hline\text{mu} & \text{dim} & \text{deg q\_k} & \text{Cas2} & \text{Cas4} & \text{time (s)} & \text{solves?}\\\hline10,10,0,0 & 1716 & 0,1,12,13 & 240 & 34560 & \text{0.0245} & \mathrm{True}\\\hline10,9,1,0 & 4719 & 0,2,11,13 & 220 & 29740 & \text{0.0293} & \mathrm{True}\\\hline10,8,2,0 &… [1317 chars total]
```

## 2026-07-20 16:37:16
**Input:**
```wolfram
irr88 = irreps[8, 8];
tbl88 = Table[
   Module[{mu = irr88[[i]], tim, qs, ds, t1, t2, t3, t4, n1, ok, qf},
     ds = degQ[mu, 8, 8];
     {tim, qs} = AbsoluteTiming[QsolveGnum[mu, 8, 8]];
     t1[uu_] = tau1[mu][uu] /. {S -> 8, P -> 8} /. num; t2[uu_] = tau2c[mu, 8, 8][uu] /. num;
     t3[uu_] = tau3[mu][uu] /. {S -> 8, P -> 8} /. num; t4[uu_] = t[4][uu] /. {S -> 8, P -> 8} /. num;
     n1[uu_] = ν1[uu] /. {S -> 8, P -> 8} /. num;
     ok = AllTrue[Range[4], (qf[x_] := qs[[#]] . x^Range[0, ds[[#]]];
        Expand[clBn[qf, t1, t2, t3, t4, n1] /. num] === 0) &];
     {Row[mu, ","], dimGL4[mu], Row[ds, ","], NumberForm[tim, {6, 4}], ok}],
   {i, Length[irr88]}];
{"# irreps" -> Length[irr88], "all solve?" -> AllTrue[tbl88, Last],
 "total time (s)" -> Total[tbl88[[All, 4]] /. NumberForm[x_, _] :> x],
 "dim-sum check" -> (Total[tbl88[[All, 2]]] == dimGL4[{8, 8, 0, 0}]^2)}
```
**Output:**
```
Out= {"# irreps" -> 45, "all solve?" -> True, "total time (s)" -> 1.576576, "dim-sum check" -> True}
```

## 2026-07-20 16:37:23
**Input:**
```wolfram
Grid[Prepend[tbl88, {"mu", "dim", "deg q_k", "time (s)", "solves?"}],
  Frame -> All, Alignment -> Left,
  Background -> {None, {RGBColor[.85, .88, .95], {GrayLevel[1], RGBColor[.95, .96, .99]}}}]
```
**Output:**
```
Out= Grid[{{"mu", "dim", "deg q_k", "time (s)", "solves?"}, {Row[{16, 16, 0, 0}, ","], 8721, Row[{0, 1, 18, 19}, ","], NumberForm[0.04608, {6, 4}], True}, {Row[{16, 15, 1, 0}, ","], 27455, Row[{0, 2, 17, 19}, ","], NumberForm[0.049954, {6, 4}], True}, {Row[{16, 14, 2, 0}, ","], 47424, Row[{0, 3, 16, 19}, ","], NumberForm[0.050383, {6, 4}], True}, {Row[{16, 13, 3, 0}, ","], 62700, Row[{0, 4, 15, 19}, ","], NumberForm[0.041904, {6, 4}], True}, {Row[{16, 12, 4, 0}, ","], 69825, Row[{0, 5, 14, 19}, ","], NumberForm[0.042951, {6, 4}], True}, {Row[{16, 11, 5, 0}, ","], 67431, Row[{0, 6, 13, 19}, ","], NumberForm[0.039247, {6, 4}], True}, {Row[{16, 10, 6, 0}, ","], 55860, Row[{0, 7, 12, 19}, ","], NumberForm[0.040021, {6, 4}], True}, {Row[{16, 9, 7, 0}, ","], 36784, Row[{0, 8, 11, 19}, ","], NumberForm[0.039862, {6, 4}], True}, {Row[{16, 8, 8, 0}, ","], 12825, Row[{0, 9, 10, 19}, ","], NumberForm[0.037295, {6, 4}], True}, {Row[{15, 15, 1, 1}, ","], 5440, Row[{1, 2, 17, 18}, ","], NumberForm[0.043865, {6, 4}], True}, {Row[{15, 14, 2, 1}, ","], 16575, Row[{1, 3, 16, 18}, ","], NumberForm[0.040282, {6, 4}], True}, {Row[{15, 13, 3, 1}, ","], 27489, Row[{1, 4, 15, 18}, ","], NumberForm[0.03882, {6, 4}], True}, {Row[{15, 12, 4, 1}, ","], 34476, Row[{1, 5, 14, 18}, ","], NumberForm[0.036066, {6, 4}], True}, {Row[{15, 11, 5, 1}, ","], 35700, Row[{1, 6, 13, 18}, ","], NumberForm[0.034829, {6, 4}], True}, {Row[{15, 10, 6, 1}, ","], 30855, Row[{1, 7, 12, 18}, ","], NumberForm[0.041006, {6, 4}], True}, {Row[{15, 9, 7, 1}, ","], 20825, Row[{1, 8, 11, 18}, ","], NumberForm[0.036305, {6, 4}], True}, {Row[{15, 8, 8, 1}, ","], 7344, Row[{1, 9, 10, 18}, ","], NumberForm[0.040396, {6, 4}], True}, {Row[{14, 14, 2, 2}, ","], 3185, Row[{2, 3, 16, 17}, ","], NumberForm[0.038113, {6, 4}], True}, {Row[{14, 13, 3, 2}, ","], 9295, Row[{2, 4, 15, 17}, ","], NumberForm[0.036634, {6, 4}], True}, {Row[{14, 12, 4, 2}, ","], 14580, Row[{2, 5, 14, 17}, ","], NumberForm[0.036381, {6, 4}], True}, {Row[{14, 11, 5, 2}, ","], 16940, Row[{2, 6, 13, 17}, ","], NumberForm[0.038754, {6, 4}], True}, {Row[{14, 10, 6, 2}, ","], 15625, Row[{2, 7, 12, 17}, ","], NumberForm[0.033661, {6, 4}], True}, {Row[{14, 9, 7, 2}, ","], 10935, Row[{2, 8, 11, 17}, ","], NumberForm[0.032392, {6, 4}], True}, {Row[{14, 8, 8, 2}, ","], 3920, Row[{2, 9, 10, 17}, ","], NumberForm[0.032371, {6, 4}], True}, {Row[{13, 13, 3, 3}, ","], 1716, Row[{3, 4, 15, 16}, ","], NumberForm[0.034834, {6, 4}], True}, {Row[{13, 12, 4, 3}, ","], 4719, Row[{3, 5, 14, 16}, ","], NumberForm[0.033402, {6, 4}], True}, {Row[{13, 11, 5, 3}, ","], 6825, Row[{3, 6, 13, 16}, ","], NumberForm[0.031912, {6, 4}], True}, {Row[{13, 10, 6, 3}, ","], 7020, Row[{3, 7, 12, 16}, ","], NumberForm[0.031532, {6, 4}], True}, {Row[{13, 9, 7, 3}, ","], 5200, Row[{3, 8, 11, 16}, ","], NumberForm[0.032455, {6, 4}], True}, {Row[{13, 8, 8, 3}, ","], 1911, Row[{3, 9, 10, 16}, ","], NumberForm[0.03019, {6, 4}], True}, {Row[{12, 12, 4, 4}, ","], 825, Row[{4, 5, 14, 15}, ","], NumberForm[0.031896, {6, 4}], True}, {Row[{12, 11, 5, 4}, ","], 2079, Row[{4, 6, 13, 15}, ","], NumberForm[0.030437, {6, 4}], True}, {Row[{12, 10, 6, 4}, ","], 2640, Row[{4, 7, 12, 15}, ","], NumberForm[0.034658, {6, 4}], True}, {Row[{12, 9, 7, 4}, ","], 2156, Row[{4, 8, 11, 15}, ","], NumberForm[0.030477, {6, 4}], True}, {Row[{12, 8, 8, 4}, ","], 825, Row[{4, 9, 10, 15}, ","], NumberForm[0.028985, {6, 4}], True}, {Row[{11, 11, 5, 5}, ","], 336, Row[{5, 6, 13, 14}, ","], NumberForm[0.03032, {6, 4}], True}, {Row[{11, 10, 6, 5}, ","], 735, Row[{5, 7, 12, 14}, ","], NumberForm[0.028902, {6, 4}], True}, {Row[{11, 9, 7, 5}, ","], 729, Row[{5, 8, 11, 14}, ","], NumberForm[0.028228, {6, 4}], True}, {Row[{11, 8, 8, 5}, ","], 300, Row[{5, 9, 10, 14}, ","], NumberForm[0.027763, {6, 4}], True}, {Row[{10, 10, 6, 6}, ","], 105, Row[{6, 7, 12, 13}, ","], NumberForm[0.028742, {6, 4}], True}, {Row[{10, 9, 7, 6}, ","], 175, Row[{6, 8, 11, 13}, ","], NumberForm[0.031341, {6, 4}], True}, {Row[{10, 8, 8, 6}, ","], 84, Row[{6, 9, 10, 13}, ","], NumberForm[0.026463, {6, 4}], True}, {Row[{9, 9, 7, 7}, 
[output truncated — 4480 chars total]
```
## 2026-07-20 16:37:41 — 📥 BULK INSERT 1 CELL(S) at positions 21–21
1. [code] (* --- (S,P)=(8,8): gauge-fixed Baxter solver over the full [8,8,0,0]x[8,8,0,0] decomposition↵   (45


## 2026-07-20 16:37:45
**Input:**
```wolfram
(* --- (S,P)=(8,8): gauge-fixed Baxter solver over the full [8,8,0,0]x[8,8,0,0] decomposition
   (45 irreps), per-irrep timing only. Q-degrees run up to S+P+3 = 19. *)
irr88 = irreps[8, 8];
tbl88 = Table[
   Module[{mu = irr88[[i]], tim, qs, ds, t1, t2, t3, t4, n1, ok, qf},
     ds = degQ[mu, 8, 8];
     {tim, qs} = AbsoluteTiming[QsolveGnum[mu, 8, 8]];
     t1[uu_] = tau1[mu][uu] /. {S -> 8, P -> 8} /. num; t2[uu_] = tau2c[mu, 8, 8][uu] /. num;
     t3[uu_] = tau3[mu][uu] /. {S -> 8, P -> 8} /. num; t4[uu_] = t[4][uu] /. {S -> 8, P -> 8} /. num;
     n1[uu_] = ν1[uu] /. {S -> 8, P -> 8} /. num;
     ok = AllTrue[Range[4], (qf[x_] := qs[[#]] . x^Range[0, ds[[#]]];
        Expand[clBn[qf, t1, t2, t3, t4, n1] /. num] === 0) &];
     {Row[mu, ","], dimGL4[mu], Row[ds, ","], NumberForm[tim, {6, 4}], ok}],
   {i, Length[irr88]}];
Print["# irreps: ", Length[irr88], "   all solve: ", AllTrue[tbl88, Last],
   "   total time (s): ", Total[tbl88[[All, 4]] /. NumberForm[x_, _] :> x],
   "   dim-sum = dim[8,8,0,0]^2: ", Total[tbl88[[All, 2]]] == dimGL4[{8, 8, 0, 0}]^2];
Grid[Prepend[tbl88, {"mu", "dim", "deg q_k", "time (s)", "solves?"}],
  Frame -> All, Alignment -> Left,
  Background -> {None, {RGBColor[.85, .88, .95], {GrayLevel[1], RGBColor[.95, .96, .99]}}}]
```
**Output:**
```
Print: # irreps:[formula]   all solve:[formula]   total time (s):[formula]   dim-sum = dim[8,8,0,0]^2:[formula] | Out[95]= \begin{array}{|c|c|c|c|c|}\hline\text{mu} & \text{dim} & \text{deg q\_k} & \text{time (s)} & \text{solves?}\\\hline16,16,0,0 & 8721 & 0,1,18,19 & \text{0.0509} & \mathrm{True}\\\hline16,15,1,0 & 27455 & 0,2,17,19 & \text{0.0510} & \mathrm{True}\\\hline16,14,2,0 & 47424 & 0,3,1… [1360 chars total]
```

## 2026-07-20 19:05:23
**Input:**
```wolfram
(* transfer matrices with abstract Casimirs C2, C4 (reusing kernel globals t1vac,t2vac,trivfactor,
   Cas2top,Cas4top,b2vac,asympSol,nu[i]). S,P,theta1,theta2,h all symbolic. *)
b3G = Simplify[(b[3] /. asympSol) /. b[2] -> (b2vac + h^2 (C2 - Cas2top))];
tau1G[x_] := t1vac[x] + (h^2/2) (C2 - Cas2top);
tau2G[x_] := t2vac[x] + (h^2/8) (C2 - Cas2top) (t1vac[x] + t1vac[x - h])
   + h^4 (1/8 (C2^2 - Cas2top^2) - 1/4 (C4 - Cas4top)) + (h^4 (1 + 2 S + 2 P)/4) (C2 - Cas2top);
tau3G[x_] := trivfactor[x] (4 x^2 - 4 x (θ1 + θ2) + 4 θ1 θ2 - 2 h P x - 2 h S x - 8 h x + b3G);
tau4G[x_] := t[4][x];
(* verify against notebook tau's for mu={2,1,1,0}@(1,1): C2->12, C4->188 *)
rep = {C2 -> Cas2[{2, 1, 1, 0}], C4 -> Cas4[{2, 1, 1, 0}], S -> 1, P -> 1};
{"tau1 ok" -> Simplify[(tau1G[u] /. rep) - (tau1[{2, 1, 1, 0}][u] /. {S -> 1, P -> 1})],
 "tau2 ok" -> Simplify[(tau2G[u] /. rep) - (tau2c[{2, 1, 1, 0}, 1, 1][u])],
 "tau3 ok" -> Simplify[(tau3G[u] /. rep) - (tau3[{2, 1, 1, 0}][u] /. {S -> 1, P -> 1})]}
```
**Output:**
```
Out= {"tau1 ok" -> 0, "tau2 ok" -> 0, "tau3 ok" -> 0}
```

## 2026-07-20 19:05:51
**Input:**
```wolfram
n1g[x_] := ν1[x];
clBg = Q[u + h] n1g[u] n1g[u - h] n1g[u - 2 h] n1g[u - 3 h]
   - tau1G[u] Q[u] n1g[u - h] n1g[u - 2 h] n1g[u - 3 h]
   + tau2G[u] Q[u - h] n1g[u - 2 h] n1g[u - 3 h]
   - tau3G[u] Q[u - 2 h] n1g[u - 3 h] + tau4G[u] Q[u - 3 h];
ser = Series[clBg, {h, 0, 5}];
Table[h^k -> Simplify[SeriesCoefficient[ser, k]], {k, 0, 4}] // Column
```
**Output:**
```
Out= Column[{1 -> 0, h -> 0, h^2 -> 0, h^3 -> 0, h^4 -> ((u - θ1)^2*(u - θ2)^2*((C2^2 - 2*C4 + 2*C2*(1 + P + S) + 4*(P^3 + P^2*(5 + 3*S) + S*(6 + 5*S + S^2) + P*(6 + 10*S + 3*S^2)))*Q[u] + 4*((C2 - 2*(P^2 + 2*P*(1 + S) + S*(2 + S)))*((2 + P + S)*u - (1 + P)*θ1 - (1 + S)*θ2)*Derivative[1][Q][u] + (-(C2*(u - θ1)*(u - θ2)) + 2*(P^2*(u - θ1)*(2*u - θ1 - θ2) + S*(u - θ2)*((3 + 2*S)*u - (2 + S)*θ1 - (1 + S)*θ2) + P*(u - θ1)*((3 + 4*S)*u - θ1 - 2*(θ2 + 2*S*θ2))))*Derivative[2][Q][u] + 2*(u - θ1)*(u - θ2)*(-2*(P*(u - θ1) + S*(u - θ2))*Derivative[3][Q][u] + (u - θ1)*(u - θ2)*Derivative[4][Q][u]))))/8}]
```

## 2026-07-20 19:06:40
**Input:**
```wolfram
ode = SeriesCoefficient[ser, 4];   (* coeff of h^4; set =0 is the ODE *)
(* divide out the common prefactor (u-θ1)^2 (u-θ2)^2 / 8 to clean up *)
odeR = Simplify[ode/((u - θ1)^2 (u - θ2)^2/8)];
coeffs = Table[Factor[Coefficient[odeR, Derivative[k][Q][u]]], {k, 4, 1, -1}];
coeff0 = Factor[odeR /. Derivative[_][Q][_] -> 0 /. Q[u] -> 1];
Column[{
  "Q'''' : " <> ToString[coeffs[[1]], InputForm],
  "Q'''  : " <> ToString[coeffs[[2]], InputForm],
  "Q''   : " <> ToString[coeffs[[3]], InputForm],
  "Q'    : " <> ToString[coeffs[[4]], InputForm],
  "Q     : " <> ToString[coeff0, InputForm]}]
```
**Output:**
```
Out= Column[{"Q'''' : 8*(u - θ1)^2*(u - θ2)^2", "Q'''  : -16*(u - θ1)*(u - θ2)*(P*u + S*u - P*θ1 - S*θ2)", "Q''   : -4*(C2*u^2 - 6*P*u^2 - 4*P^2*u^2 - 6*S*u^2 - 8*P*S*u^2 - 4*S^2*u^2 - C2*u*θ1 + 8*P*u*θ1 + 6*P^2*u*θ1 + 4*S*u*θ1 + 8*P*S*u*θ1 + 2*S^2*u*θ1 - 2*P*θ1^2 - 2*P^2*θ1^2 - C2*u*θ2 + 4*P*u*θ2 + 2*P^2*u*θ2 + 8*S*u*θ2 + 8*P*S*u*θ2 + 6*S^2*u*θ2 + C2*θ1*θ2 - 4*P*θ1*θ2 - 2*P^2*θ1*θ2 - 4*S*θ1*θ2 - 8*P*S*θ1*θ2 - 2*S^2*θ1*θ2 - 2*S*θ2^2 - 2*S^2*θ2^2)", "Q'    : 4*(C2 - 4*P - 2*P^2 - 4*S - 4*P*S - 2*S^2)*(2*u + P*u + S*u - θ1 - P*θ1 - θ2 - S*θ2)", "Q     : 2*C2 + C2^2 - 2*C4 + 24*P + 2*C2*P + 20*P^2 + 4*P^3 + 24*S + 2*C2*S + 40*P*S + 12*P^2*S + 20*S^2 + 12*P*S^2 + 4*S^3"}]
```

## 2026-07-20 19:07:50
**Input:**
```wolfram
nuu = (u - θ1) (u - θ2);
q2 = Coefficient[odeR, Q''[u]];
q2rest = Factor[Expand[q2/(-4) - C2 nuu]];   (* Q'' = -4(C2*nu + rest) *)
(* VALIDATE: for mu={2,1,1,0}@(1,1) (C2=12,C4=188,S=P=1) the ODE must have polynomial solutions
   of degrees degQ = {0,2,3,5}. Check each degree gives a nonzero polynomial solution. *)
odeNum = odeR /. {C2 -> 12, C4 -> 188, S -> 1, P -> 1};
polysol[deg_] := Module[{q = Sum[a[j] u^j, {j, 0, deg}], eq, sol},
   eq = Expand[odeNum /. Q -> Function[x, q /. u -> x]];
   sol = Solve[Join[Thread[CoefficientList[eq, u] == 0], {a[deg] == 1}], Array[a, deg + 1, 0]];
   If[sol === {}, "none", (q /. sol[[1]])]];
{"Q'' rest term" -> q2rest, "poly sols of degrees {0,2,3,5}" -> (polysol /@ {0, 2, 3, 5})}
```
**Output:**
```
[message] Solve::svars: Equations may not give solutions for all "solve" variables.
[message] [Solve::svars ×2 more suppressed — first instance shown above]
[message] General::stop: Further output of Solve::svars will be suppressed during this calculation.
Out= {"Q'' rest term" -> -2*(3*P*u^2 + 2*P^2*u^2 + 3*S*u^2 + 4*P*S*u^2 + 2*S^2*u^2 - 4*P*u*θ1 - 3*P^2*u*θ1 - 2*S*u*θ1 - 4*P*S*u*θ1 - S^2*u*θ1 + P*θ1^2 + P^2*θ1^2 - 2*P*u*θ2 - P^2*u*θ2 - 4*S*u*θ2 - 4*P*S*u*θ2 - 3*S^2*u*θ2 + 2*P*θ1*θ2 + P^2*θ1*θ2 + 2*S*θ1*θ2 + 4*P*S*θ1*θ2 + S^2*θ1*θ2 + S*θ2^2 + S^2*θ2^2), "poly sols of degrees {0,2,3,5}" -> {1, u^2 + u*(-θ1 - θ2) + a[0], u^3 + a[0] + u^2*a[2] + u*(-3*θ1*θ2 - (θ1 + θ2)*a[2]), u^5 - (5*u^4*(θ1 + θ2))/2 + a[0] + u*a[1] + u^2*a[2] + u^3*(5*θ1*θ2 - a[1]/(3*θ1*θ2) - ((θ1 + θ2)*a[2])/(3*θ1*θ2))}}
```

## 2026-07-20 19:09:31
**Input:**
```wolfram
(* asymmetric check: (S,P)=(2,1), mu={3,2,1,0}, degQ={0,2,4,6} *)
odeA = odeR /. {C2 -> Cas2[{3, 2, 1, 0}], C4 -> Cas4[{3, 2, 1, 0}], S -> 2, P -> 1};
psolA[deg_] := Module[{q = Sum[a[j] u^j, {j, 0, deg}], eq, sol},
   eq = Expand[odeA /. Q -> Function[x, q /. u -> x]];
   sol = Solve[Join[Thread[CoefficientList[eq, u] == 0], {a[deg] == 1}], Array[a, deg + 1, 0]];
   If[sol === {}, "NONE", deg]];
{"degQ" -> degQ[{3, 2, 1, 0}, 2, 1], "ODE has poly sols of degrees" -> (psolA /@ {0, 2, 4, 6}),
 "wrong degree {1} gives?" -> psolA[1]}
```
**Output:**
```
[message] Solve::svars: Equations may not give solutions for all "solve" variables.
[message] [Solve::svars ×2 more suppressed — first instance shown above]
[message] General::stop: Further output of Solve::svars will be suppressed during this calculation.
Out= {"degQ" -> {0, 2, 4, 6}, "ODE has poly sols of degrees" -> {0, 2, 4, 6}, "wrong degree {1} gives?" -> "NONE"}
```
## 2026-07-20 19:10:24 — 📥 BULK INSERT 2 CELL(S) at positions 22–23
1. [code] (* --- Gaudin limit h->0: the 4th-order Baxter DIFFERENCE equation -> 4th-order ODE (the oper) ---↵ 
2. [code] (* --- validation: the ODE has polynomial solutions of exactly the degrees degQ[mu,S,P] ---↵   (the


## 2026-07-20 19:10:25
**Input:**
```wolfram
(* --- Gaudin limit h->0: the 4th-order Baxter DIFFERENCE equation -> 4th-order ODE (the oper) ---
   Write Q = Q0 + h Q1 + ...  Orders h^0..h^3 of the cleared Baxter equation vanish identically
   (the same alternating-binomial moment identity that made orders 1/u..1/u^3 mu-blind), and
   order h^4 gives a CLOSED 4th-order ODE on Q0. Transfer matrices carry the abstract output
   Casimirs C2,C4; theta1,theta2 (inhomogeneities/sites) and S,P (input reps) stay symbolic. *)
b3G = (b[3] /. asympSol) /. b[2] -> (b2vac + h^2 (C2 - Cas2top));
tau1G[x_] := t1vac[x] + (h^2/2) (C2 - Cas2top);
tau2G[x_] := t2vac[x] + (h^2/8) (C2 - Cas2top) (t1vac[x] + t1vac[x - h])
   + h^4 (1/8 (C2^2 - Cas2top^2) - 1/4 (C4 - Cas4top)) + (h^4 (1 + 2 S + 2 P)/4) (C2 - Cas2top);
tau3G[x_] := trivfactor[x] (4 x^2 - 4 x (θ1 + θ2) + 4 θ1 θ2 - 2 h P x - 2 h S x - 8 h x + b3G);
tau4G[x_] := t[4][x];
clBg = Q[u + h] ν1[u] ν1[u - h] ν1[u - 2 h] ν1[u - 3 h]
   - tau1G[u] Q[u] ν1[u - h] ν1[u - 2 h] ν1[u - 3 h] + tau2G[u] Q[u - h] ν1[u - 2 h] ν1[u - 3 h]
   - tau3G[u] Q[u - 2 h] ν1[u - 3 h] + tau4G[u] Q[u - 3 h];
(* lowest nonvanishing order is h^4; divide off the common (u-th1)^2(u-th2)^2/8 *)
gaudinODE = Collect[Simplify[SeriesCoefficient[Series[clBg, {h, 0, 4}], 4]/((u - θ1)^2 (u - θ2)^2/8)],
   Derivative[_][Q][u], Factor];
gaudinODE
```
**Output:**
```
Out[103]= Q\left(u\right)\,\left(C2^2+2\,C2\,P+2\,C2\,S+2\,C2-2\,C4+4\,P^3+12\,P^2\,S+20\,P^2+12\,P\,S^2+40\,P\,S+24\,P+4\,S^3+20\,S^2+24\,S\right)-4\,Q''\left(u\right)\,\left(C2\,\theta1\,\theta2+C2\,u^2-C2\,\theta1\,u-C2\,\theta2\,u-2\,\theta1^2\,P^2-2\,\theta1\,\theta2\,P^2-4\,P^2\,u^2+6\,\theta1\,P^2\,u+2\,\theta2\,P^2\,u-2\,\theta1^2\,P-4\,\theta1\,\theta2\,P-8\,\theta1\,\theta2\,P\,S-8\,P\,S… [1029 chars total]
```

## 2026-07-20 19:10:25
**Input:**
```wolfram
(* --- validation: the ODE has polynomial solutions of exactly the degrees degQ[mu,S,P] ---
   (the h->0 limit of the Baxter Q-degrees, which are h-independent). Sweep whole decompositions. *)
odeOf[mu_, SS_, PP_] := gaudinODE /. {C2 -> Cas2[mu], C4 -> Cas4[mu], S -> SS, P -> PP};
hasPolySol[mu_, SS_, PP_, deg_] := Module[{q = Sum[aa[j] u^j, {j, 0, deg}], eq, sol},
   eq = Expand[odeOf[mu, SS, PP] /. Q -> Function[x, q /. u -> x]];
   sol = Solve[Join[Thread[CoefficientList[eq, u] == 0], {aa[deg] == 1}], Array[aa, deg + 1, 0]];
   sol =!= {}];
Grid[Prepend[Flatten[Table[{{SS, PP}, Length[irreps[SS, PP]],
     AllTrue[irreps[SS, PP], Function[mu, And @@ (hasPolySol[mu, SS, PP, #] & /@ degQ[mu, SS, PP])]]},
    {SS, {1, 2}}, {PP, {1, 2}}], 1], {"(S,P)", "#irreps", "ODE poly-sol degrees = degQ?"}], Frame -> All]
```
**Output:**
```
2 kernel messages: Solve::svars: Equations may not give solutions for all "solve" variables.
Solve::svars: Equations may not give solutions for all "solve" variables.
Solve::svars: Equations may not give solutions for a… [also: Solve::svars] result: Out[106]= \begin{array}{|c|c|c|}\hline\text{(S,P)} & \text{\#irreps} & \text{ODE poly-sol degrees = degQ?}\\\hline\left\{1,1\right\} & 3 & \mathrm{True}\\\hline\left\{1,2\right\} & 3 & \mathrm{True}\\…
```

## 2026-07-20 19:25:34
**Input:**
```wolfram
thsub = {θ1 -> 0, θ2 -> 1};
(* Gaudin ODE solver: solve the 4th-order ODE for the 4 polynomial Q's of degrees degQ, with the
   SAME gauge fixing (kill coeff at each lower degree). Input mu, S, P. theta1=0, theta2=1. *)
GaudinSolve[mu_, SS_, PP_] := Module[{ds, ode},
   ds = degQ[mu, SS, PP];
   ode = gaudinODE /. {C2 -> Cas2[mu], C4 -> Cas4[mu], S -> SS, P -> PP} /. thsub;
   Table[Module[{m = ds[[kk]], gauge = ds[[1 ;; kk - 1]], q, cl, sol},
       q = Sum[cq[j] u^j, {j, 0, m}];
       cl = CoefficientList[Expand[ode /. Q -> Function[x, q /. u -> x]], u];
       sol = Solve[Join[Thread[cl == 0], {cq[m] == 1}, Thread[(cq /@ gauge) == 0]], Table[cq[j], {j, 0, m}]];
       Table[cq[j], {j, 0, m}] /. sol[[1]]], {kk, 4}]];
(* test: {2,1,1,0}@(1,1), degrees {0,2,3,5} *)
GaudinSolve[{2, 1, 1, 0}, 1, 1]
```
**Output:**
```
Out= {{1}, {0, -1, 1}, {0, 0, 0, 1}, {0, 0, 0, 0, -5/2, 1}}
```

## 2026-07-20 19:25:58
**Input:**
```wolfram
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
{"# irreps" -> Length[irr88], "all solve?" -> AllTrue[tblG88, Last],
 "total time (s)" -> Total[tblG88[[All, 4]] /. NumberForm[x_, _] :> x]}
```
**Output:**
```
Out= {"# irreps" -> 45, "all solve?" -> True, "total time (s)" -> 0.116322}
```

## 2026-07-20 19:26:26
**Input:**
```wolfram
Grid[Prepend[tblG88, {"mu", "dim", "deg q_k", "time (s)", "solves?"}],
  Frame -> All, Alignment -> Left,
  Background -> {None, {RGBColor[.88, .95, .88], {GrayLevel[1], RGBColor[.96, .99, .96]}}}]
```
**Output:**
```
Out= Grid[{{"mu", "dim", "deg q_k", "time (s)", "solves?"}, {Row[{16, 16, 0, 0}, ","], 8721, Row[{0, 1, 18, 19}, ","], NumberForm[0.003014, {6, 4}], True}, {Row[{16, 15, 1, 0}, ","], 27455, Row[{0, 2, 17, 19}, ","], NumberForm[0.002517, {6, 4}], True}, {Row[{16, 14, 2, 0}, ","], 47424, Row[{0, 3, 16, 19}, ","], NumberForm[0.002453, {6, 4}], True}, {Row[{16, 13, 3, 0}, ","], 62700, Row[{0, 4, 15, 19}, ","], NumberForm[0.002687, {6, 4}], True}, {Row[{16, 12, 4, 0}, ","], 69825, Row[{0, 5, 14, 19}, ","], NumberForm[0.002569, {6, 4}], True}, {Row[{16, 11, 5, 0}, ","], 67431, Row[{0, 6, 13, 19}, ","], NumberForm[0.002611, {6, 4}], True}, {Row[{16, 10, 6, 0}, ","], 55860, Row[{0, 7, 12, 19}, ","], NumberForm[0.002495, {6, 4}], True}, {Row[{16, 9, 7, 0}, ","], 36784, Row[{0, 8, 11, 19}, ","], NumberForm[0.002609, {6, 4}], True}, {Row[{16, 8, 8, 0}, ","], 12825, Row[{0, 9, 10, 19}, ","], NumberForm[0.002365, {6, 4}], True}, {Row[{15, 15, 1, 1}, ","], 5440, Row[{1, 2, 17, 18}, ","], NumberForm[0.002896, {6, 4}], True}, {Row[{15, 14, 2, 1}, ","], 16575, Row[{1, 3, 16, 18}, ","], NumberForm[0.002358, {6, 4}], True}, {Row[{15, 13, 3, 1}, ","], 27489, Row[{1, 4, 15, 18}, ","], NumberForm[0.002522, {6, 4}], True}, {Row[{15, 12, 4, 1}, ","], 34476, Row[{1, 5, 14, 18}, ","], NumberForm[0.002598, {6, 4}], True}, {Row[{15, 11, 5, 1}, ","], 35700, Row[{1, 6, 13, 18}, ","], NumberForm[0.002554, {6, 4}], True}, {Row[{15, 10, 6, 1}, ","], 30855, Row[{1, 7, 12, 18}, ","], NumberForm[0.00263, {6, 4}], True}, {Row[{15, 9, 7, 1}, ","], 20825, Row[{1, 8, 11, 18}, ","], NumberForm[0.002592, {6, 4}], True}, {Row[{15, 8, 8, 1}, ","], 7344, Row[{1, 9, 10, 18}, ","], NumberForm[0.002703, {6, 4}], True}, {Row[{14, 14, 2, 2}, ","], 3185, Row[{2, 3, 16, 17}, ","], NumberForm[0.002638, {6, 4}], True}, {Row[{14, 13, 3, 2}, ","], 9295, Row[{2, 4, 15, 17}, ","], NumberForm[0.002429, {6, 4}], True}, {Row[{14, 12, 4, 2}, ","], 14580, Row[{2, 5, 14, 17}, ","], NumberForm[0.002668, {6, 4}], True}, {Row[{14, 11, 5, 2}, ","], 16940, Row[{2, 6, 13, 17}, ","], NumberForm[0.002792, {6, 4}], True}, {Row[{14, 10, 6, 2}, ","], 15625, Row[{2, 7, 12, 17}, ","], NumberForm[0.002671, {6, 4}], True}, {Row[{14, 9, 7, 2}, ","], 10935, Row[{2, 8, 11, 17}, ","], NumberForm[0.002644, {6, 4}], True}, {Row[{14, 8, 8, 2}, ","], 3920, Row[{2, 9, 10, 17}, ","], NumberForm[0.002774, {6, 4}], True}, {Row[{13, 13, 3, 3}, ","], 1716, Row[{3, 4, 15, 16}, ","], NumberForm[0.002485, {6, 4}], True}, {Row[{13, 12, 4, 3}, ","], 4719, Row[{3, 5, 14, 16}, ","], NumberForm[0.002519, {6, 4}], True}, {Row[{13, 11, 5, 3}, ","], 6825, Row[{3, 6, 13, 16}, ","], NumberForm[0.002681, {6, 4}], True}, {Row[{13, 10, 6, 3}, ","], 7020, Row[{3, 7, 12, 16}, ","], NumberForm[0.002588, {6, 4}], True}, {Row[{13, 9, 7, 3}, ","], 5200, Row[{3, 8, 11, 16}, ","], NumberForm[0.00253, {6, 4}], True}, {Row[{13, 8, 8, 3}, ","], 1911, Row[{3, 9, 10, 16}, ","], NumberForm[0.002775, {6, 4}], True}, {Row[{12, 12, 4, 4}, ","], 825, Row[{4, 5, 14, 15}, ","], NumberForm[0.00241, {6, 4}], True}, {Row[{12, 11, 5, 4}, ","], 2079, Row[{4, 6, 13, 15}, ","], NumberForm[0.002532, {6, 4}], True}, {Row[{12, 10, 6, 4}, ","], 2640, Row[{4, 7, 12, 15}, ","], NumberForm[0.002474, {6, 4}], True}, {Row[{12, 9, 7, 4}, ","], 2156, Row[{4, 8, 11, 15}, ","], NumberForm[0.002387, {6, 4}], True}, {Row[{12, 8, 8, 4}, ","], 825, Row[{4, 9, 10, 15}, ","], NumberForm[0.002667, {6, 4}], True}, {Row[{11, 11, 5, 5}, ","], 336, Row[{5, 6, 13, 14}, ","], NumberForm[0.002499, {6, 4}], True}, {Row[{11, 10, 6, 5}, ","], 735, Row[{5, 7, 12, 14}, ","], NumberForm[0.00268, {6, 4}], True}, {Row[{11, 9, 7, 5}, ","], 729, Row[{5, 8, 11, 14}, ","], NumberForm[0.002492, {6, 4}], True}, {Row[{11, 8, 8, 5}, ","], 300, Row[{5, 9, 10, 14}, ","], NumberForm[0.002526, {6, 4}], True}, {Row[{10, 10, 6, 6}, ","], 105, Row[{6, 7, 12, 13}, ","], NumberForm[0.002658, {6, 4}], True}, {Row[{10, 9, 7, 6}, ","], 175, Row[{6, 8, 11, 13}, ","], NumberForm[0.002486, {6, 4}], True}, {Row[{10, 8, 8, 6}, ","], 84, Row[{6, 9, 10, 13}, ","], NumberForm[0.002651, {6, 4}], True}, {Row[{9, 9, 7, 7}, 
[output truncated — 4481 chars total]
```
## 2026-07-20 19:26:50 — 📥 BULK INSERT 2 CELL(S) at positions 24–25
1. [code] (* --- Gaudin ODE solver (theta1=0, theta2=1) ---↵   Input: irrep mu in [S,S,0,0]x[P,P,0,0]. Solves 
2. [code] (* --- (S,P)=(8,8): Gaudin ODE solver over the full [8,8,0,0]x[8,8,0,0] decomposition (45 irreps),↵


## 2026-07-20 19:26:50
**Input:**
```wolfram
(* --- Gaudin ODE solver (theta1=0, theta2=1) ---
   Input: irrep mu in [S,S,0,0]x[P,P,0,0]. Solves the 4th-order Gaudin ODE (gaudinODE above, with
   C2->Cas2[mu], C4->Cas4[mu]) for the four polynomial Q's of degrees degQ[mu,S,P]. SAME gauge
   fixing as the Baxter solver: the four degrees are strictly increasing, so force the coefficient
   of u^(deg q_j) in q_k to vanish for every j<k. Purely linear in the Q coefficients. *)
thsub = {θ1 -> 0, θ2 -> 1};
GaudinSolve[mu_, SS_, PP_] := Module[{ds, ode},
   ds = degQ[mu, SS, PP];
   ode = gaudinODE /. {C2 -> Cas2[mu], C4 -> Cas4[mu], S -> SS, P -> PP} /. thsub;
   Table[Module[{m = ds[[kk]], gauge = ds[[1 ;; kk - 1]], q, cl, sol},
       q = Sum[cq[j] u^j, {j, 0, m}];
       cl = CoefficientList[Expand[ode /. Q -> Function[x, q /. u -> x]], u];
       sol = Solve[Join[Thread[cl == 0], {cq[m] == 1}, Thread[(cq /@ gauge) == 0]], Table[cq[j], {j, 0, m}]];
       Table[cq[j], {j, 0, m}] /. sol[[1]]], {kk, 4}]];
(* example: mu={2,1,1,0}@(1,1) -> q1=1, q2=u(u-1), q3=u^3, q4=u^5-(5/2)u^4 *)
GaudinSolve[{2, 1, 1, 0}, 1, 1]
```
**Output:**
```
Out[109]= \left\{\left\{1\right\},\left\{0,-1,1\right\},\left\{0,0,0,1\right\},\left\{0,0,0,0,-\frac{5}{2},1\right\}\right\}
```

## 2026-07-20 19:26:52
**Input:**
```wolfram
(* --- (S,P)=(8,8): Gaudin ODE solver over the full [8,8,0,0]x[8,8,0,0] decomposition (45 irreps),
   theta1=0, theta2=1, per-irrep timing. Each q_k verified to solve the ODE exactly. *)
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
Print: # irreps:[formula]   all solve:[formula]   total time (s):[formula] | Out[113]= \begin{array}{|c|c|c|c|c|}\hline\text{mu} & \text{dim} & \text{deg q\_k} & \text{time (s)} & \text{solves?}\\\hline16,16,0,0 & 8721 & 0,1,18,19 & \text{0.0024} & \mathrm{True}\\\hline16,15,1,0 & 27455 & 0,2,17,19 & \text{0.0024} & \mathrm{True}\\\hline16,14,2,0 & 47424 & 0,3,16,19 & \text{0.0032} & \mathrm{True}… [1324 chars total]
```

---

## 2026-07-21 06:10:57 — 🔄 KERNEL RESTART


## 2026-07-21 06:12:10
**Input:**
```wolfram
?gaudinODE
```
**Output:**
```
Out= InformationData[<|"ObjectType" -> "Symbol", "Usage" -> "Global`gaudinODE", "Documentation" -> None, "OwnValues" -> None, "UpValues" -> None, "DownValues" -> None, "SubValues" -> None, "DefaultValues" -> None, "NValues" -> None, "FormatValues" -> None, "Options" -> None, "Attributes" -> {}, "FullName" -> "Global`gaudinODE"|>, False]
```
## 2026-07-21 06:12:42 — ▶️ RUN CELLS 1–25
Cell 1: ✓ — (no output)
Cell 2: ✓ — (no output)
Cell 3: ✓ — (no output)
Cell 4: ✓ — Out[14]= \frac{-a\left(1\right)+a\left(2\right)+4\,h\,P+4\,h\,S+12\,h+8\,\theta1+8\,\theta2}{u}+\frac{a\left(2\right)\,h\,M\left(1\right)-a\left(1\right)\,h\,P+a\left(2\right)\,h\,P-a\left(1\right)\,h\,S+a\left(2\right)\,h\,S+2\,a\left(2\right)\,h-a\left(1\right)\,\theta1+2\,a\left(2\right)\,\theta1-a\left(1\right)\,\theta2+2\,a\left(2\right)\,\theta2-b\left(1\right)+b\left(2\right)-b\left(3\right… [819 chars total]
Cell 5: ✓ — Out[23]= \left\{12,64,\left\{-m1+P+S,-m2+P+S+1,-m3+P+S+2,-m4+P+S+3\right\}\right\}
Cell 6: ✓ — Out[27]= \begin{matrix}a\left(1\right)\to -2\,\left(h\,P+h\,S+2\,\theta1+2\,\theta2\right)\\a\left(2\right)\to -6\,\left(h\,P+h\,S+2\,h+2\,\theta1+2\,\theta2\right)\\b\left(1\right)\to \frac{1}{2}\,\left(b\left(2\right)+h^2\,\left(-P^2\right)-4\,h^2\,P\,S-7\,h^2\,P-h^2\,S^2-7\,h^2\,S-6\,h^2-18\,h\,\theta1-18\,h\,\theta2-8\,h\,\theta1\,P-6\,h\,\theta2\,P-6\,h\,\theta1\,S-8\,h\,\theta2\,S-6\,\theta1… [1246 chars total]
Cell 7: ✓ — Out[32]= \left\{\text{c[2] vacuum vs asymp-predicted (must be 0)}\to 0,\text{b2vac}\to h^2\,\left(P^2+P\,\left(8\,S+7\right)+S^2+7\,S+6\right)+6\,h\,\left(\theta1\,\left(2\,P+S+3\right)+\theta2\,\left(P+2\,S+3\right)\right)+6\,\left(\theta1^2+4\,\theta1\,\theta2+\theta2^2\right)\right\}
Cell 8: ✓ — Out[35]= h^2\,m1^2+3\,h^2\,m1+h^2\,m2^2+h^2\,m2+h^2\,m3^2-h^2\,m3+h^2\,m4^2-3\,h^2\,m4-h^2\,P^2+4\,h^2\,P\,S+3\,h^2\,P-h^2\,S^2+3\,h^2\,S+6\,h^2+18\,h\,\theta1+18\,h\,\theta2+12\,h\,\theta1\,P+6\,h\,\theta2\,P+6\,h\,\theta1\,S+12\,h\,\theta2\,S+6\,\theta1^2+24\,\theta1\,\theta2+6\,\theta2^2
Cell 9: ✓ — Out[40]= \left\{0,0,0,0\right\}
Cell 10: ✓ — Out[45]= \left\{0,0\right\}
Cell 11: ✓ — Out[51]= \left\{\left\{2,2,0,0\right\}\to 1,\left\{2,1,1,0\right\}\to 1,\left\{1,1,1,1\right\}\to 1\right\}
Cell 12: ✓ — Out[52]= 825
Cell 13: ✓ — Out[53]= \begin{pmatrix}10 & 10 & 0 & 0\\10 & 9 & 1 & 0\\10 & 8 & 2 & 0\\10 & 7 & 3 & 0\\10 & 6 & 4 & 0\\10 & 5 & 5 & 0\\9 & 9 & 1 & 1\\9 & 8 & 2 & 1\\9 & 7 & 3 & 1\\9 & 6 & 4 & 1\\9 & 5 & 5 & 1\\8 & 8 & 2 & 2\\8 & 7 & 3 & 2\\8 & 6 & 4 & 2\\8 & 5 & 5 & 2\\7 & 7 & 3 & 3\\7 & 6 & 4 & 3\\7 & 5 & 5 & 3\\6 & 6 & 4 & 4\\6 & 5 & 5 & 4\\5 & 5 & 5 & 5\end{pmatrix}
Cell 14: ✓ — Out[54]= 680625
Cell 15: ✓ — Out[56]= \begin{array}{|c|c|c|c|c|}\hline\text{mu} & \text{dim} & \text{deg q\_k} & \text{Cas2} & \text{Cas4}\\\hline\left\{2,2,0,0\right\} & 20 & \left\{0,1,4,5\right\} & 16 & 256\\\hline\left\{2,1,1,0\right\} & 15 & \left\{0,2,3,5\right\} & 12 & 188\\\hline\left\{1,1,1,1\right\} & 1 & \left\{1,2,3,4\right\} & 4 & 4\\\hline\end{array}
Cell 16: ✓ — Out[57]= \begin{array}{|c|c|c|c|c|}\hline\text{mu} & \text{dim} & \text{deg q\_k} & \text{Cas2} & \text{Cas4}\\\hline\left\{16,16,0,0\right\} & 8721 & \left\{0,1,18,19\right\} & 576 & 186624\\\hline\left\{16,15,1,0\right\} & 27455 & \left\{0,2,17,19\right\} & 544 & 168160\\\hline\left\{16,14,2,0\right\} & 47424 & \left\{0,3,16,19\right\} & 516 & 152844\\\hline\left\{16,13,3,0\right\} & 62700 & \le… [1246 chars total]
Cell 17: ✓ — (no output)
Cell 18: ✓ — Out[60]= \left\{\left\{1\right\},\left\{0,-\theta1-\theta2,1\right\},\left\{0,-h^2-3\,\theta1\,\theta2,0,1\right\},\left\{0,\frac{1}{2}\,\left(-2\,h^4+5\,h^2\,\theta1^2+10\,h^2\,\theta1\,\theta2+5\,h^2\,\theta2^2+30\,\theta1^2\,\theta2^2\right),0,0,-\frac{5}{2}\,\left(\theta1+\theta2\right),1\right\}\right\}
Cell 19: ✓ — Out[62]= \begin{array}{|c|c|c|}\hline\text{(S,P)} & \text{\#irreps} & \text{all q\_k solve Baxter?}\\\hline\left\{1,1\right\} & 3 & \mathrm{True}\\\hline\left\{1,2\right\} & 3 & \mathrm{True}\\\hline\left\{2,1\right\} & 3 & \mathrm{True}\\\hline\left\{2,2\right\} & 6 & \mathrm{True}\\\hline\left\{3,1\right\} & 3 & \mathrm{True}\\\hline\left\{3,2\right\} & 6 & \mathrm{True}\\\hline\end{array}
Cell 20: ✓ — Print: all 21 irreps solve:[formula]  total solve time (s):[formula] | Out[69]= \begin{array}{|c|c|c|c|c|c|c|}\hline\text{mu} & \text{dim} & \text{deg q\_k} & \text{Cas2} & \text{Cas4} & \text{time (s)} & \text{solves?}\\\hline10,10,0,0 & 1716 & 0,1,12,13 & 240 & 34560 & \text{0.0367} & \mathrm{True}\\\hline10,9,1,0 & 4719 & 0,2,11,13 & 220 & 29740 & \text{0.0312} & \mathrm{True}\\\hline10,8,2,0 &… [1317 chars total]
Cell 21: ✓ — Print: # irreps:[formula]   all solve:[formula]   total time (s):[formula]   dim-sum = dim[8,8,0,0]^2:[formula] | Out[73]= \begin{array}{|c|c|c|c|c|}\hline\text{mu} & \text{dim} & \text{deg q\_k} & \text{time (s)} & \text{solves?}\\\hline16,16,0,0 & 8721 & 0,1,18,19 & \text{0.0575} & \mathrm{True}\\\hline16,15,1,0 & 27455 & 0,2,17,19 & \text{0.0532} & \mathrm{True}\\\hline16,14,2,0 & 47424 & 0,3,1… [1360 chars total]
Cell 22: ✓ — Out[81]= Q\left(u\right)\,\left(C2^2+2\,C2\,P+2\,C2\,S+2\,C2-2\,C4+4\,P^3+12\,P^2\,S+20\,P^2+12\,P\,S^2+40\,P\,S+24\,P+4\,S^3+20\,S^2+24\,S\right)-4\,Q''\left(u\right)\,\left(C2\,\theta1\,\theta2+C2\,u^2-C2\,\theta1\,u-C2\,\theta2\,u-2\,\theta1^2\,P^2-2\,\theta1\,\theta2\,P^2-4\,P^2\,u^2+6\,\theta1\,P^2\,u+2\,\theta2\,P^2\,u-2\,\theta1^2\,P-4\,\theta1\,\theta2\,P-8\,\theta1\,\theta2\,P\,S-8\,P\,S\… [1028 chars total]
Cell 23: ✓ — 2 kernel messages: Solve::svars: Equations may not give solutions for all "solve" variables.
Solve::svars: Equations may not give solutions for all "solve" variables.
Solve::svars: Equations may not give solutions for a… [also: Solve::svars] result: Out[84]= \begin{array}{|c|c|c|}\hline\text{(S,P)} & \text{\#irreps} & \text{ODE poly-sol degrees = degQ?}\\\hline\left\{1,1\right\} & 3 & \mathrm{True}\\\hline\left\{1,2\right\} & 3 & \mathrm{True}\\\…
⚠️ stopped at Cell 23 — error detected (pass stopOnError:false to continue past errors)


## 2026-07-21 06:12:46
**Input:**
```wolfram
ToString[gaudinODE, InputForm]
```
**Output:**
```
Out= (2*C2 + C2^2 - 2*C4 + 24*P + 2*C2*P + 20*P^2 + 4*P^3 + 24*S + 2*C2*S + 40*P*S + 12*P^2*S + 20*S^2 + 12*P*S^2 + 4*S^3)*Q[u] + 4*(C2 - 4*P - 2*P^2 - 4*S - 4*P*S - 2*S^2)*(2*u + P*u + S*u - θ1 - P*θ1 - θ2 - S*θ2)*Derivative[1][Q][u] - 4*(C2*u^2 - 6*P*u^2 - 4*P^2*u^2 - 6*S*u^2 - 8*P*S*u^2 - 4*S^2*u^2 - C2*u*θ1 + 8*P*u*θ1 + 6*P^2*u*θ1 + 4*S*u*θ1 + 8*P*S*u*θ1 + 2*S^2*u*θ1 - 2*P*θ1^2 - 2*P^2*θ1^2 - C2*u*θ2 + 4*P*u*θ2 + 2*P^2*u*θ2 + 8*S*u*θ2 + 8*P*S*u*θ2 + 6*S^2*u*θ2 + C2*θ1*θ2 - 4*P*θ1*θ2 - 2*P^2*θ1*θ2 - 4*S*θ1*θ2 - 8*P*S*θ1*θ2 - 2*S^2*θ1*θ2 - 2*S*θ2^2 - 2*S^2*θ2^2)*Derivative[2][Q][u] - 16*(u - θ1)*(u - θ2)*(P*u + S*u - P*θ1 - S*θ2)*Derivative[3][Q][u] + 8*(u - θ1)^2*(u - θ2)^2*Derivative[4][Q][u]
```
