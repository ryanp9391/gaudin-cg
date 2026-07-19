# AI Action Log — Baxter_L2_XXX_SU4

## 2026-07-19 15:47:29 — 🗑️ DELETE CODE CELL 1
1+1

## 2026-07-19 15:47:41 — 📥 BULK INSERT 4 CELL(S) at positions 1–4
1. [code] (*su(4) generators*)↵ps[rep_] := ps[rep] = Block[{res = (Prepend[rep, #] & /@ (Tuples[Table[Range[re
2. [code] ClearAll[Ee,id,EE,idd];↵Ee[λ_][i_,j_]:=JJ[i,j,λ];↵id[λ_]:=IdentityMatrix[dim[λ]];↵idd[λ1_,λ2_]:=Kron
3. [code] (*central charge*)↵ClearAll[C1];↵C1[λ1_,λ2_]:=Sum[EE[λ1,λ2][i,i],{i,4}];
4. [code] (*Lax operator for a general gl(4) rep λ={λ1,λ2,λ3,λ4} at the site*)↵ClearAll[L];↵L[λ_][i_,j_][u_]:=


## 2026-07-19 15:47:41
**Input:**
```wolfram
(*su(4) generators*)
ps[rep_] := ps[rep] = Block[{res = (Prepend[rep, #] & /@ (Tuples[Table[Range[rep[[1, i + 1]], rep[[1, i]]], {i, 1, Length[rep[[1]]] - 1}]]))}, If[Length[res[[1, 1]]] == 1, Reverse@Sort@res, Reverse@Sort@Flatten[ps /@ res, 1]]];
dim[rep_] := dim[rep] = Length[ps[{rep}]];
J0[k_][pat_] := J0[k][pat] = Sum[pat[[k, i]], {i, 1, k}] - Sum[pat[[k - 1, i]], {i, 1, k - 1}];
Jp[k_][pat_] := Jp[k][pat] = -Sum[Product[pat[[k, i]] - i - pat[[k + 1, j]] + j, {j, Range[k + 1]}]/Product[pat[[k, i]] - i - pat[[k, j]] + j, {j, Complement[Range[k], {i}]}] v[ReplacePart[pat, {k, i} -> pat[[k, i]] + 1]], {i, 1, k}];
Jm[k_][pat_] := Jm[k][pat] = Sum[Product[pat[[k, i]] - i - pat[[k - 1, j]] + j, {j, Range[k - 1]}]/Product[pat[[k, i]] - i - pat[[k, j]] + j, {j, Complement[Range[k], {i}]}] v[ReplacePart[pat, {k, i} -> pat[[k, i]] - 1]], {i, 1, k}];
JJ[k_, l_, rep_] := JJ[k, l, rep] = Which[k == l, DiagonalMatrix[J0[k] /@ ps[{rep}]], l == k + 1, Transpose[Table[Coefficient[Jp[k][ps[{rep}][[i]]], v[#]] & /@ ps[{rep}], {i, 1, Length[ps[{rep}]]}]], l == k - 1, Transpose[Table[Coefficient[Jm[l][ps[{rep}][[i]]], v[#]] & /@ ps[{rep}], {i, 1, Length[ps[{rep}]]}]], l > k + 1, JJ[k, k + 1, rep] . JJ[k + 1, l, rep] - JJ[k + 1, l, rep] . JJ[k, k + 1, rep], l < k - 1, JJ[k, k - 1, rep] . JJ[k - 1, l, rep] - JJ[k - 1, l, rep] . JJ[k, k - 1, rep]];
```
**Output:**
```
(no output)
```

## 2026-07-19 15:47:42
**Input:**
```wolfram
ClearAll[Ee,id,EE,idd];
Ee[λ_][i_,j_]:=JJ[i,j,λ];
id[λ_]:=IdentityMatrix[dim[λ]];
idd[λ1_,λ2_]:=KroneckerProduct[id[λ1],id[λ2]];
EE[λ1_,λ2_][i_,j_]:=KroneckerProduct[Ee[λ1][i,j],id[λ2]]+KroneckerProduct[id[λ1],Ee[λ2][i,j]];
```
**Output:**
```
(no output)
```

## 2026-07-19 15:47:42
**Input:**
```wolfram
(*central charge*)
ClearAll[C1];
C1[λ1_,λ2_]:=Sum[EE[λ1,λ2][i,i],{i,4}];
```
**Output:**
```
(no output)
```

## 2026-07-19 15:47:42
**Input:**
```wolfram
(*Lax operator for a general gl(4) rep λ={λ1,λ2,λ3,λ4} at the site*)
ClearAll[L];
L[λ_][i_,j_][u_]:=IdentityMatrix[dim[λ]] u KroneckerDelta[i,j] - h JJ[j,i,λ];
```
**Output:**
```
(no output)
```
## 2026-07-19 15:47:46 — 📥 BULK INSERT 1 CELL(S) at positions 5–5
1. [code] (*checking central charge constraint = total boxes, on general gl(4) reps.*)↵Module[{reps={{{1,1,0,0


## 2026-07-19 15:47:46
**Input:**
```wolfram
(*checking central charge constraint = total boxes, on general gl(4) reps.*)
Module[{reps={{{1,1,0,0},{1,0,0,0}},{{2,1,1,0},{1,0,0,0}},{{1,1,0,0},{1,1,0,0}}}},
 Union@Flatten@Table[
   C1[r[[1]],r[[2]]]-(Total[r[[1]]]+Total[r[[2]]]) idd[r[[1]],r[[2]]],{r,reps}]]
```
**Output:**
```
Out[16]= \left\{0\right\}
```

## 2026-07-19 15:47:49
**Input:**
```wolfram
{dim[{1,0,0,0}],dim[{1,1,0,0}],dim[{2,1,1,0}]}
```
**Output:**
```
Out= {4, 6, 15}
```
## 2026-07-19 15:50:48 — 📥 BULK INSERT 3 CELL(S) at positions 6–8
1. [code] (*χ1..χ4: elementary symmetric functions of the four diagonal-twist eigenvalues z[1..4]*)↵χ1=z[1]+z[
2. [code] (*Numeric parameters (shared across Section B). θ,h as in the gl(3) notebook;↵  four DISTINCT unit-m
3. [code] (*monodromy with no twist and diagonal twist*)↵ClearAll[T,Tg];↵T[λ1_,λ2_][i_,j_][u_]:=Sum[KroneckerP


## 2026-07-19 15:50:49
**Input:**
```wolfram
(*χ1..χ4: elementary symmetric functions of the four diagonal-twist eigenvalues z[1..4]*)
χ1=z[1]+z[2]+z[3]+z[4];
χ2=z[1]z[2]+z[1]z[3]+z[1]z[4]+z[2]z[3]+z[2]z[4]+z[3]z[4];
χ3=z[1]z[2]z[3]+z[1]z[2]z[4]+z[1]z[3]z[4]+z[2]z[3]z[4];
χ4=z[1]z[2]z[3]z[4];
```
**Output:**
```
(no output)
```

## 2026-07-19 15:50:49
**Input:**
```wolfram
(*Numeric parameters (shared across Section B). θ,h as in the gl(3) notebook;
  four DISTINCT unit-modulus twist eigenvalues z[1..4] (generic → non-degenerate).
  χ1..χ4 were defined above as symmetric functions of z[i], so they become
  numeric automatically once z[i] are assigned.*)
θ1 = 1/3;
θ2 = 1/7;
h = 1;
z[1] = N[Exp[I Zeta[3]], 30];
z[2] = N[Exp[I Zeta[5]], 30];
z[3] = N[Exp[I Zeta[7]], 30];
z[4] = N[Exp[I Zeta[9]], 30];
{θ1, θ2, h, z[1], z[2], z[3], z[4], χ1, χ2, χ3, χ4}
```
**Output:**
```
Out[28]= \left\{\frac{1}{3},\frac{1}{7},1,0.360439875147226537266677650596+0.932782448593374716484835982071\,\mathrm{i},0.508867380990499828163321372267+0.860844927129079093188644393616\,\mathrm{i},0.533257880589599372595832345093+0.845952736734794526483381618453\,\mathrm{i},0.538611213022438940451622759023+0.842554426257732564010591177696\,\mathrm{i},1.94117634974976467847745412698+3.482134538714… [665 chars total]
```

## 2026-07-19 15:50:49
**Input:**
```wolfram
(*monodromy with no twist and diagonal twist*)
ClearAll[T,Tg];
T[λ1_,λ2_][i_,j_][u_]:=Sum[KroneckerProduct[L[λ1][i,ii][u-θ1],L[λ2][ii,j][u-θ2]],{ii,4}];
Tg[λ1_,λ2_][i_,j_][u_]:= z[i] T[λ1,λ2][i,j][u];
```
**Output:**
```
(no output)
```
## 2026-07-19 15:50:53 — 📥 BULK INSERT 1 CELL(S) at positions 9–9
1. [code] {Abs[χ1 - Total[{z[1],z[2],z[3],z[4]}]], Min[Abs[Subtract @@@ Subsets[{z[1],z[2],z[3],z[4]}, {2}]]]}


## 2026-07-19 15:50:53
**Input:**
```wolfram
{Abs[χ1 - Total[{z[1],z[2],z[3],z[4]}]], Min[Abs[Subtract @@@ Subsets[{z[1],z[2],z[3],z[4]}, {2}]]]}
```
**Output:**
```
Out[32]= \left\{0.\times 10^{-29},0.006340873933063257902352694531\right\}
```
## 2026-07-19 15:52:13 — 📥 BULK INSERT 1 CELL(S) at positions 10–10
1. [code] (*Yangian weights - eigenvalues of T[i,i] on the highest-weight state.↵  General rep: T[i,i]|hw> = (


## 2026-07-19 15:52:13
**Input:**
```wolfram
(*Yangian weights - eigenvalues of T[i,i] on the highest-weight state.
  General rep: T[i,i]|hw> = (u-θ1-h λ1[[i]])(u-θ2-h λ2[[i]]) |hw>, for i=1..4,
  from L[λ][i,i]|hw>=(u-h λ_i)|hw>.*)
ClearAll[ν,HWS];
ν[λ1_,λ2_][i_][u_]:=(u-θ1-h λ1[[i]])(u-θ2-h λ2[[i]]);
HWS[λ1_,λ2_]:=Table[{KroneckerDelta[1,ii]},{ii,1,dim[λ1] dim[λ2]}];
```
**Output:**
```
(no output)
```
## 2026-07-19 15:52:17 — 📥 BULK INSERT 1 CELL(S) at positions 11–11
1. [code] (*checking highest-weights against the general ν, on non-symmetric gl(4) reps*)↵Module[{reps={{{1,1,


## 2026-07-19 15:52:17
**Input:**
```wolfram
(*checking highest-weights against the general ν, on non-symmetric gl(4) reps*)
Module[{reps={{{1,1,0,0},{1,0,0,0}},{{2,1,1,0},{1,0,0,0}},{{1,1,0,0},{1,1,0,0}}}},
 Union@Flatten@Table[
   T[r[[1]],r[[2]]][i,i][u].HWS[r[[1]],r[[2]]]-ν[r[[1]],r[[2]]][i][u] HWS[r[[1]],r[[2]]],
   {r,reps},{i,4}]]
```
**Output:**
```
Out[36]= \left\{0\right\}
```
## 2026-07-19 15:53:05 — 📥 BULK INSERT 3 CELL(S) at positions 12–14
1. [code] (*quantum minors (q-determinants) of the monodromy matrix*)↵ClearAll[qmin,qming];↵qmin[λ1_, λ2_][II_
2. [code] (*transfer matrices with diagonal twist*)↵ClearAll[tg];↵tg[λ1_,λ2_][a_,1][u_]:= Total[qming[λ1, λ2][
3. [code] (*value of quantum determinant (χ4 carried here; general ν)*)↵ClearAll[qdetT];↵qdetT[λ1_,λ2_][u_]:=χ


## 2026-07-19 15:53:06
**Input:**
```wolfram
(*quantum minors (q-determinants) of the monodromy matrix*)
ClearAll[qmin,qming];
qmin[λ1_, λ2_][II_List, JJ_List][u_] := qmin[λ1, λ2][II, JJ][u] = Sum[Signature[σ] Dot @@ Table[T[λ1, λ2][II[[σ[[k]]]], JJ[[k]]][u - (k - 1) h], {k, Length[II]}], {σ, Permutations[Range[Length[II]]]}];
qming[λ1_, λ2_][II_List, JJ_List][u_] := qming[λ1, λ2][II, JJ][u] = Sum[Signature[σ] Dot @@ Table[Tg[λ1, λ2][II[[σ[[k]]]], JJ[[k]]][u - (k - 1) h], {k, Length[II]}], {σ, Permutations[Range[Length[II]]]}];
```
**Output:**
```
(no output)
```

## 2026-07-19 15:53:06
**Input:**
```wolfram
(*transfer matrices with diagonal twist*)
ClearAll[tg];
tg[λ1_,λ2_][a_,1][u_]:= Total[qming[λ1, λ2][#, #][u] & /@ Subsets[{1, 2, 3, 4}, {a}]]
```
**Output:**
```
(no output)
```

## 2026-07-19 15:53:06
**Input:**
```wolfram
(*value of quantum determinant (χ4 carried here; general ν)*)
ClearAll[qdetT];
qdetT[λ1_,λ2_][u_]:=χ4 ν[λ1,λ2][1][u-3h] ν[λ1,λ2][2][u-2h] ν[λ1,λ2][3][u-h] ν[λ1,λ2][4][u];
```
**Output:**
```
(no output)
```
## 2026-07-19 15:53:12 — 📥 BULK INSERT 1 CELL(S) at positions 15–15
1. [code] (*qdetT central value == tg[·][4,1]; on non-symmetric reps this is the ν shift-ordering guard.↵  Cho


## 2026-07-19 15:53:14
**Input:**
```wolfram
(*qdetT central value == tg[·][4,1]; on non-symmetric reps this is the ν shift-ordering guard.
  Chop wraps Simplify because z[k] (hence χ4) are high-precision numeric, leaving ~10^-27
  floating "0" coefficients on symbolic-u polynomials that Union never collapses to a literal 0
  even though every entry is numerically exactly 0.*)
Module[{reps={{{1,1,0,0},{1,0,0,0}},{{2,1,1,0},{1,0,0,0}},{{1,1,0,0},{1,1,0,0}}}},
 Union@Flatten@Table[
   Chop[Simplify[qdetT[r[[1]],r[[2]]][u] idd[r[[1]],r[[2]]]-tg[r[[1]],r[[2]]][4,1][u]],10^-15],{r,reps}]]
```
**Output:**
```
Out[44]= \left\{0\right\}
```
## 2026-07-19 15:54:21 — 📥 BULK INSERT 1 CELL(S) at positions 16–16
1. [code] (*coefficient matrices of tau1=tg[.][1,1] (deg 2), tau2=tg[.][2,1] (deg 4), tau3=tg[.][3,1] (deg 6),


## 2026-07-19 15:54:21
**Input:**
```wolfram
(*coefficient matrices of tau1=tg[.][1,1] (deg 2), tau2=tg[.][2,1] (deg 4), tau3=tg[.][3,1] (deg 6),
  extracted by CoefficientList in a dummy uu. c1={C0,C1,C2}, c2={D0..D4}, c3={G0..G6}.*)
ClearAll[tgCoeffs];
tgCoeffs[λ1_, λ2_] := tgCoeffs[λ1, λ2] = Module[{uu, m1, m2, m3, cl1, cl2, cl3},
   m1 = tg[λ1, λ2][1, 1][uu];
   m2 = tg[λ1, λ2][2, 1][uu];
   m3 = tg[λ1, λ2][3, 1][uu];
   cl1 = Map[PadRight[CoefficientList[#, uu], 3] &, m1, {2}];
   cl2 = Map[PadRight[CoefficientList[#, uu], 5] &, m2, {2}];
   cl3 = Map[PadRight[CoefficientList[#, uu], 7] &, m3, {2}];
   {Table[cl1[[All, All, k]], {k, 1, 3}], Table[cl2[[All, All, k]], {k, 1, 5}], Table[cl3[[All, All, k]], {k, 1, 7}]}
];
```
**Output:**
```
(no output)
```
## 2026-07-19 15:54:25 — 📥 BULK INSERT 1 CELL(S) at positions 17–17
1. [code] (*Cartan operators (gl(4) weight generators H1=E11-E22, H2=E22-E33, H3=E33-E44) on the two-site spac


## 2026-07-19 15:54:25
**Input:**
```wolfram
(*Cartan operators (gl(4) weight generators H1=E11-E22, H2=E22-E33, H3=E33-E44) on the two-site space*)
ClearAll[H1op, H2op, H3op];
H1op[λ1_, λ2_] := EE[λ1, λ2][1, 1] - EE[λ1, λ2][2, 2];
H2op[λ1_, λ2_] := EE[λ1, λ2][2, 2] - EE[λ1, λ2][3, 3];
H3op[λ1_, λ2_] := EE[λ1, λ2][3, 3] - EE[λ1, λ2][4, 4];
```
**Output:**
```
(no output)
```
## 2026-07-19 15:54:34 — 📥 BULK INSERT 1 CELL(S) at positions 18–18
1. [code] (*hard-assert: τ1,τ2,τ3 leading coeffs are χ1 Id, χ2 Id, χ3 Id; commuting family pairwise-commutes.↵


## 2026-07-19 15:54:39
**Input:**
```wolfram
(*hard-assert: τ1,τ2,τ3 leading coeffs are χ1 Id, χ2 Id, χ3 Id; commuting family pairwise-commutes.
  Tested on adjoint⊗fundamental (d=60), which has interior weight-multiplicity.*)
Module[{λ1={2,1,1,0},λ2={1,0,0,0},d,c1,c2,c3,C2,D4,G6,mats,comm},
  d=dim[λ1] dim[λ2];
  {c1,c2,c3}=tgCoeffs[λ1,λ2];
  C2=c1[[3]]; D4=c2[[5]]; G6=c3[[7]];
  If[Max[Abs[Flatten[C2-χ1 IdentityMatrix[d]]]]>10^-13,Print["C2 FAIL"];Abort[]];
  If[Max[Abs[Flatten[D4-χ2 IdentityMatrix[d]]]]>10^-13,Print["D4 FAIL"];Abort[]];
  If[Max[Abs[Flatten[G6-χ3 IdentityMatrix[d]]]]>10^-13,Print["G6 FAIL"];Abort[]];
  mats=Join[c1,c2,c3,{H1op[λ1,λ2],H2op[λ1,λ2],H3op[λ1,λ2]}];
  comm=Max[Table[Max[Abs[Flatten[mats[[i]].mats[[j]]-mats[[j]].mats[[i]]]]],
      {i,Length[mats]},{j,Length[mats]}]];
  If[comm>10^-12,Print["COMMUTATION FAIL: ",comm];Abort[]];
  {"C2=χ1 Id OK","D4=χ2 Id OK","G6=χ3 Id OK","maxCommutator"->comm}
]
```
**Output:**
```
Out[51]= \left\{\text{C2=$\chi$1 Id OK},\text{D4=$\chi$2 Id OK},\text{G6=$\chi$3 Id OK},\text{maxCommutator}\to 0.\times 10^{-23}\right\}
```
## 2026-07-19 15:55:10 — 📥 BULK INSERT 1 CELL(S) at positions 19–19
1. [code] (*gl(4) weight multiplicities in the tensor product λ1⊗λ2.↵  siteWeightMult[λ]: tally the gl(4) weig


## 2026-07-19 15:55:10
**Input:**
```wolfram
(*gl(4) weight multiplicities in the tensor product λ1⊗λ2.
  siteWeightMult[λ]: tally the gl(4) weight (J0[1..4]) of every GT pattern of λ.
  wtMult: recover (n1,n2,n3,n4) from (H1,H2,H3) with N=total boxes, then convolve the site tallies.
  H1=n1-n2, H2=n2-n3, H3=n3-n4, n1+n2+n3+n4=N ⇒ n4=(N-H1-2H2-3H3)/4.
  NB: wrap the composite weight key in Key[...] — Lookup treats a bare List as several keys.*)
ClearAll[siteWeightMult,wtMult];
siteWeightMult[λ_]:=siteWeightMult[λ]=Counts[({J0[1][#],J0[2][#],J0[3][#],J0[4][#]}&)/@ps[{λ}]];
wtMult[λ1_,λ2_][H1_,H2_,H3_]:=Module[{Ntot=Total[λ1]+Total[λ2],n1,n2,n3,n4,m1,m2},
  n4=(Ntot-H1-2 H2-3 H3)/4; n3=n4+H3; n2=n4+H2+H3; n1=n4+H1+H2+H3;
  If[!(IntegerQ[n1]&&IntegerQ[n2]&&IntegerQ[n3]&&IntegerQ[n4]),Return[0]];
  m1=siteWeightMult[λ1]; m2=siteWeightMult[λ2];
  Total[KeyValueMap[Function[{w,c}, c Lookup[m2,Key[{n1,n2,n3,n4}-w],0]],m1]]
];
```
**Output:**
```
(no output)
```
## 2026-07-19 15:55:18 — 📥 BULK INSERT 1 CELL(S) at positions 20–20
1. [code] (*hard-assert sum-rule: Σ over all (H1,H2,H3) of wtMult == dim[λ1] dim[λ2], for every sweep rep*)↵Mo


## 2026-07-19 15:55:18
**Input:**
```wolfram
(*hard-assert sum-rule: Σ over all (H1,H2,H3) of wtMult == dim[λ1] dim[λ2], for every sweep rep*)
Module[{reps={{{1,0,0,0},{1,0,0,0}},{{1,0,0,0},{1,1,1,0}},{{1,1,0,0},{1,0,0,0}},{{1,1,0,0},{1,1,0,0}},{{2,1,1,0},{1,0,0,0}}},bad},
  bad=Select[reps,Function[r,
     Module[{Ntot=Total[r[[1]]]+Total[r[[2]]],d,tot},
       d=dim[r[[1]]] dim[r[[2]]];
       tot=Total@Flatten@Table[wtMult[r[[1]],r[[2]]][H1,H2,H3],
           {H1,-Ntot,Ntot},{H2,-Ntot,Ntot},{H3,-Ntot,Ntot}];
       tot=!=d]]];
  If[bad=!={},Print["SUM-RULE FAIL: ",bad];Abort[]];
  {"all reps sum-rule OK",
   "adjoint⊗fund @(H1,H2,H3)=(1,0,0) mult"->wtMult[{2,1,1,0},{1,0,0,0}][1,0,0]}
]
```
**Output:**
```
Out[55]= \left\{\text{all reps sum-rule OK},\text{adjoint$\otimes$fund @(H1,H2,H3)=(1,0,0) mult}\to 6\right\}
```
## 2026-07-19 15:56:03 — 📥 BULK INSERT 1 CELL(S) at positions 21–21
1. [code] ClearAll[TauEigensystem4];↵TauEigensystem4::rankdef = "Generic combination failed to resolve the `1`


## 2026-07-19 15:56:03
**Input:**
```wolfram
ClearAll[TauEigensystem4];
TauEigensystem4::rankdef = "Generic combination failed to resolve the `1`-dim eigenspace for (λ1,λ2)=(`2`,`3`).";
TauEigensystem4::notgenuine = "Eigenvector #`1` for (λ1,λ2)=(`2`,`3`) failed the common-eigenvector residual check (residual=`4`).";
TauEigensystem4::noninteger = "Weight `1` not near-integer for eigenvector #`2` at (λ1,λ2)=(`3`,`4`).";
TauEigensystem4::degnull = "Sector (H1,H2,H3)=(`1`,`2`,`3`) for (λ1,λ2)=(`4`,`5`) has `6` state(s) but wtMult predicts `7`.";
TauEigensystem4[λ1_, λ2_] := TauEigensystem4[λ1, λ2] = Module[
   {d, c1, c2, c3, mats, nontriv, rSets, rset, Kmix, vals, vecs, rayleigh, recs, byHW, table},
   d = dim[λ1] dim[λ2];
   {c1, c2, c3} = tgCoeffs[λ1, λ2];
   (* full family for extraction; drop the scalar leading coeffs c1[[3]],c2[[5]],c3[[7]] from the
      splitting combination since they carry no state info *)
   mats = Join[c1, c2, c3, {H1op[λ1, λ2], H2op[λ1, λ2], H3op[λ1, λ2]}];
   nontriv = Join[c1[[1 ;; 2]], c2[[1 ;; 4]], c3[[1 ;; 6]], {H1op[λ1, λ2], H2op[λ1, λ2], H3op[λ1, λ2]}]; (* 15 matrices *)
   rSets = Table[N[Table[Prime[100 i + 7 j]^(1/3), {j, Length[nontriv]}], 30], {i, 1, 8}];
   rset = SelectFirst[rSets,
     Module[{ee, vv}, {ee, vv} = Eigensystem[#.nontriv];
        MatrixRank[vv] == d && Length[Union[Chop[ee, 10^-12]]] == d] &,
     Missing["NoGenericR"]];
   If[MissingQ[rset], Message[TauEigensystem4::rankdef, d, λ1, λ2]; Abort[]];
   Kmix = rset . nontriv;
   {vals, vecs} = Eigensystem[N[Kmix, 30]];
   rayleigh[m_, v_] := (v . (m . v))/(v . v);
   recs = Table[
     Module[{v = vecs[[i]], res, H1v, H2v, H3v, t1c, t2c, t3c},
       res = Max[Abs[#.v - rayleigh[#, v] v]] & /@ mats // Max;
       If[res > 10^-13, Message[TauEigensystem4::notgenuine, i, λ1, λ2, res]; Abort[]];
       t1c = rayleigh[#, v] & /@ c1;           (* {C0,C1,C2} quotients *)
       t2c = rayleigh[#, v] & /@ c2;           (* {D0..D4} quotients *)
       t3c = rayleigh[#, v] & /@ c3;           (* {G0..G6} quotients *)
       H1v = rayleigh[H1op[λ1, λ2], v]; H2v = rayleigh[H2op[λ1, λ2], v]; H3v = rayleigh[H3op[λ1, λ2], v];
       If[Abs[H1v - Round[Re[H1v]]] > 10^-10, Message[TauEigensystem4::noninteger, H1v, i, λ1, λ2]; Abort[]];
       If[Abs[H2v - Round[Re[H2v]]] > 10^-10, Message[TauEigensystem4::noninteger, H2v, i, λ1, λ2]; Abort[]];
       If[Abs[H3v - Round[Re[H3v]]] > 10^-10, Message[TauEigensystem4::noninteger, H3v, i, λ1, λ2]; Abort[]];
       <|"H1" -> Round[Re[H1v]], "H2" -> Round[Re[H2v]], "H3" -> Round[Re[H3v]], "t1c" -> t1c, "t2c" -> t2c, "t3c" -> t3c, "vec" -> v|>],
     {i, 1, d}];
   byHW = GroupBy[recs, {#["H1"], #["H2"], #["H3"]} &];
   table = Association @ Flatten @ KeyValueMap[
      Function[{hw, group},
        Module[{sorted, predicted},
          predicted = wtMult[λ1, λ2][hw[[1]], hw[[2]], hw[[3]]];
          (* deterministic order within a sector: full real+imag parts of all tau1/tau2/tau3 charges *)
          sorted = SortBy[group, Join[Re[#["t1c"]], Im[#["t1c"]], Re[#["t2c"]], Im[#["t2c"]], Re[#["t3c"]], Im[#["t3c"]]] &];
          If[Length[sorted] != predicted,
            Message[TauEigensystem4::degnull, hw[[1]], hw[[2]], hw[[3]], λ1, λ2, Length[sorted], predicted]; Abort[]];
          MapIndexed[
            Function[{rec, idx}, Module[{n = First[idx] - 1, vv = rec["vec"], lc},
              lc = SelectFirst[vv, Abs[#] > 10^-13 Max[Abs[vv]] &];
              {hw[[1]], hw[[2]], hw[[3]], n} -> <|"tau1c" -> rec["t1c"], "tau2c" -> rec["t2c"], "tau3c" -> rec["t3c"], "Psi" -> vv/lc|>]],
            sorted]]],
      byHW];
   If[Length[table] != d, Message[TauEigensystem4::degnull, "total", "-", "-", λ1, λ2, Length[table], d]; Abort[]];
   table
];
```
**Output:**
```
(no output)
```
## 2026-07-19 15:56:10 — 📥 BULK INSERT 1 CELL(S) at positions 22–22
1. [code] ClearAll[Psi, tau1, tau2, tau3];↵Psi::bnds = "No state at (H1,H2,H3,n)=(`1`,`2`,`3`,`4`) for (λ1,λ2)


## 2026-07-19 15:56:10
**Input:**
```wolfram
ClearAll[Psi, tau1, tau2, tau3];
Psi::bnds = "No state at (H1,H2,H3,n)=(`1`,`2`,`3`,`4`) for (λ1,λ2)=(`5`,`6`).";
Psi[λ1_, λ2_][H1_, H2_, H3_, n_] := Lookup[TauEigensystem4[λ1, λ2], Key[{H1, H2, H3, n}],
    (Message[Psi::bnds, H1, H2, H3, n, λ1, λ2]; Missing["NoState", {H1, H2, H3, n}])]["Psi"];
tau1[λ1_, λ2_][H1_, H2_, H3_, n_][u_] := TauEigensystem4[λ1, λ2][{H1, H2, H3, n}, "tau1c"] . u^Range[0, 2];
tau2[λ1_, λ2_][H1_, H2_, H3_, n_][u_] := TauEigensystem4[λ1, λ2][{H1, H2, H3, n}, "tau2c"] . u^Range[0, 4];
tau3[λ1_, λ2_][H1_, H2_, H3_, n_][u_] := TauEigensystem4[λ1, λ2][{H1, H2, H3, n}, "tau3c"] . u^Range[0, 6];
```
**Output:**
```
(no output)
```
## 2026-07-19 15:56:17 — 📥 BULK INSERT 1 CELL(S) at positions 23–23
1. [code] Module[{λ1 = {2, 1, 1, 0}, λ2 = {1, 0, 0, 0}, es, worst},↵  es = TauEigensystem4[λ1, λ2];↵  worst =


## 2026-07-19 15:56:20
**Input:**
```wolfram
Module[{λ1 = {2, 1, 1, 0}, λ2 = {1, 0, 0, 0}, es, worst},
  es = TauEigensystem4[λ1, λ2];
  worst = Max@KeyValueMap[
     Function[{k, rec},
       Module[{H1 = k[[1]], H2 = k[[2]], H3 = k[[3]], n = k[[4]], v = rec["Psi"], uu = 17/5},
         Max[
          Max[Abs[(tg[λ1, λ2][1, 1][uu] . v) - tau1[λ1, λ2][H1, H2, H3, n][uu] v]],
          Max[Abs[(tg[λ1, λ2][2, 1][uu] . v) - tau2[λ1, λ2][H1, H2, H3, n][uu] v]],
          Max[Abs[(tg[λ1, λ2][3, 1][uu] . v) - tau3[λ1, λ2][H1, H2, H3, n][uu] v]]]]],
     es];
  If[worst > 10^-11, Print["EIGEN RESIDUAL FAIL: ", worst]; Abort[]];
  {"nstates" -> Length[es], "worstResidual@u=17/5" -> worst}
]
```
**Output:**
```
Out[68]= \left\{\text{nstates}\to 60,\text{worstResidual@u=17/5}\to 4.6281415107404180764670\times 10^{-23}\right\}
```
## 2026-07-19 15:57:02 — 📥 BULK INSERT 3 CELL(S) at positions 24–26
1. [code] (*single source of truth for the Yangian weight of T11*)↵ClearAll[ν1];↵ν1[λ1_,λ2_][u_]:=ν[λ1,λ2][1][
2. [code] (*4th-order TQ residual functional (operator form). A1=τ1/ν1, A2=τ2/(ν1 ν1[-h]),↵  A3=τ3/(ν1 ν1[-h] 
3. [code] (*large-u asymptotics: alphaCoeffs = {α1,α2,α3,α4}, the 1/u coefficients of A1..A4 in BaxOp.↵  NB: f


## 2026-07-19 15:57:02
**Input:**
```wolfram
(*single source of truth for the Yangian weight of T11*)
ClearAll[ν1];
ν1[λ1_,λ2_][u_]:=ν[λ1,λ2][1][u];
```
**Output:**
```
(no output)
```

## 2026-07-19 15:57:02
**Input:**
```wolfram
(*4th-order TQ residual functional (operator form). A1=τ1/ν1, A2=τ2/(ν1 ν1[-h]),
  A3=τ3/(ν1 ν1[-h] ν1[-2h]), A4=qdetT/(ν1 ν1[-h] ν1[-2h] ν1[-3h]); qdetT carries χ4.*)
ClearAll[BaxOp];
BaxOp[λ1_, λ2_][H1_, H2_, H3_, n_][Qf_][u_] :=
  Qf[u + h]
   - tau1[λ1, λ2][H1, H2, H3, n][u]/ν1[λ1, λ2][u] Qf[u]
   + tau2[λ1, λ2][H1, H2, H3, n][u]/(ν1[λ1, λ2][u] ν1[λ1, λ2][u - h]) Qf[u - h]
   - tau3[λ1, λ2][H1, H2, H3, n][u]/(ν1[λ1, λ2][u] ν1[λ1, λ2][u - h] ν1[λ1, λ2][u - 2 h]) Qf[u - 2 h]
   + qdetT[λ1, λ2][u]/(ν1[λ1, λ2][u] ν1[λ1, λ2][u - h] ν1[λ1, λ2][u - 2 h] ν1[λ1, λ2][u - 3 h]) Qf[u - 3 h];
```
**Output:**
```
(no output)
```

## 2026-07-19 15:57:03
**Input:**
```wolfram
(*large-u asymptotics: alphaCoeffs = {α1,α2,α3,α4}, the 1/u coefficients of A1..A4 in BaxOp.
  NB: for x0=Infinity, SeriesCoefficient[f,{u,Infinity,1}] is the coeff of u^(-1) — use n=+1.*)
ClearAll[alphaCoeffs];
alphaCoeffs[λ1_, λ2_][H1_, H2_, H3_, n_] := alphaCoeffs[λ1, λ2][H1, H2, H3, n] = Module[{u, A1, A2, A3, A4},
   A1 = tau1[λ1, λ2][H1, H2, H3, n][u]/ν1[λ1, λ2][u];
   A2 = tau2[λ1, λ2][H1, H2, H3, n][u]/(ν1[λ1, λ2][u] ν1[λ1, λ2][u - h]);
   A3 = tau3[λ1, λ2][H1, H2, H3, n][u]/(ν1[λ1, λ2][u] ν1[λ1, λ2][u - h] ν1[λ1, λ2][u - 2 h]);
   A4 = qdetT[λ1, λ2][u]/(ν1[λ1, λ2][u] ν1[λ1, λ2][u - h] ν1[λ1, λ2][u - 2 h] ν1[λ1, λ2][u - 3 h]);
   SeriesCoefficient[#, {u, Infinity, 1}] & /@ {A1, A2, A3, A4}
];
```
**Output:**
```
(no output)
```
## 2026-07-19 15:57:13 — 📥 BULK INSERT 2 CELL(S) at positions 27–28
1. [code] (*asymptotic master formula (gl(4)): for each twist eigenvalue z[k] (k=1..4) the large-u balance↵  o
2. [code] (*closed-form degree: M_k = Λ_1 - n_k, Λ=λ1+λ2, Λ_1=(λ1+λ2)[[1]];↵  n4=(N-H1-2H2-3H3)/4, n3=n4+H3, n


## 2026-07-19 15:57:13
**Input:**
```wolfram
(*asymptotic master formula (gl(4)): for each twist eigenvalue z[k] (k=1..4) the large-u balance
  of the 4th-order TQ equation fixes the polynomial degree M[k]. Derived from
    M = -[Σ_j (-1)^j α_j z^{1-j}] / (h Σ_j (-1)^j χ_j (1-j) z^{1-j}), j=0..4,
  giving numerator (α1 - α2/z + α3/z^2 - α4/z^3) and
  denominator h (z - χ2/z + 2χ3/z^2 - 3χ4/z^3), the z-characteristic factor of the companion
  char poly z^4-χ1 z^3+χ2 z^2-χ3 z+χ4 differentiated and divided by z^3.*)
ClearAll[Mdeg];
Mdeg[λ1_, λ2_][H1_, H2_, H3_, n_] := Mdeg[λ1, λ2][H1, H2, H3, n] = Module[{a, M},
   a = alphaCoeffs[λ1, λ2][H1, H2, H3, n];
   M = Table[
      (a[[1]] - a[[2]]/z[k] + a[[3]]/z[k]^2 - a[[4]]/z[k]^3)/(h (z[k] - χ2/z[k] + 2 χ3/z[k]^2 - 3 χ4/z[k]^3)),
      {k, 1, 4}];
   If[Max[Abs[M - Round[Re[M]]]] > 10^-9, Print["Mdeg NON-INTEGER: ", M]; Abort[]];
   Round[Re[M]]
];
```
**Output:**
```
(no output)
```

## 2026-07-19 15:57:13
**Input:**
```wolfram
(*closed-form degree: M_k = Λ_1 - n_k, Λ=λ1+λ2, Λ_1=(λ1+λ2)[[1]];
  n4=(N-H1-2H2-3H3)/4, n3=n4+H3, n2=n4+H2+H3, n1=n4+H1+H2+H3, N=Total[λ1]+Total[λ2].*)
ClearAll[nvec, Mpred];
nvec[H1_, H2_, H3_, Ntot_] := Module[{n4 = (Ntot - H1 - 2 H2 - 3 H3)/4}, {n4 + H1 + H2 + H3, n4 + H2 + H3, n4 + H3, n4}];
Mpred[λ1_, λ2_][H1_, H2_, H3_] := (λ1 + λ2)[[1]] - nvec[H1, H2, H3, Total[λ1] + Total[λ2]];
```
**Output:**
```
(no output)
```
## 2026-07-19 15:57:23 — 📥 BULK INSERT 1 CELL(S) at positions 29–29
1. [code] (*hard-assert: closed form M_k=Λ_1-n_k matches the asymptotic Mdeg exactly over the sweep.*)↵Module[


## 2026-07-19 15:57:25
**Input:**
```wolfram
(*hard-assert: closed form M_k=Λ_1-n_k matches the asymptotic Mdeg exactly over the sweep.*)
Module[{sweep={{{1,0,0,0},{1,0,0,0}},{{1,0,0,0},{1,1,1,0}},{{1,1,0,0},{1,0,0,0}},{{1,1,0,0},{1,1,0,0}},{{2,1,1,0},{1,0,0,0}}},worst},
  worst=Max@Table[
     Module[{λ1=rep[[1]],λ2=rep[[2]],es},
       es=TauEigensystem4[λ1,λ2];
       Max@KeyValueMap[Function[{key,rec},
          Max[Abs[Mdeg[λ1,λ2][key[[1]],key[[2]],key[[3]],key[[4]]]-Mpred[λ1,λ2][key[[1]],key[[2]],key[[3]]]]]],es]],
     {rep,sweep}];
  If[worst>0,Print["DEGREE-FORMULA MISMATCH: ",worst];Abort[]];
  {"closed-form M_k=Λ_1-n_k matches Mdeg over full sweep, worst"->worst}
]
```
**Output:**
```
Out[80]= \left\{\text{closed-form M\_k=$\Lambda$\_1-n\_k matches Mdeg over full sweep, worst}\to 0\right\}
```
## 2026-07-19 15:58:29 — 📥 BULK INSERT 3 CELL(S) at positions 30–32
1. [code] (*polynomial from ascending coeff list {c0..cM}: qpoly[clist,uu]=Sum ck uu^k. x-substitution avoids 
2. [code] (*QSolve4 solves the 4th-order Baxter (TQ) equation for the four monic Q-polynomials q_1..q_4↵  (one
3. [code] (*Q accessor: Q_k(u)=z[k]^(u/h) q_k(u), q_k the monic coeff list from QSolve4.*)↵ClearAll[Qfun];↵Qfu


## 2026-07-19 15:58:29
**Input:**
```wolfram
(*polynomial from ascending coeff list {c0..cM}: qpoly[clist,uu]=Sum ck uu^k. x-substitution avoids 0^0.*)
ClearAll[qpoly];
qpoly[clist_, uu_] := Module[{x}, (clist . x^Range[0, Length[clist] - 1]) /. x -> uu];
```
**Output:**
```
(no output)
```

## 2026-07-19 15:58:29
**Input:**
```wolfram
(*QSolve4 solves the 4th-order Baxter (TQ) equation for the four monic Q-polynomials q_1..q_4
  (one per twist z[k]), degrees M_k from Mdeg. For fixed k, Q_k=z[k]^(u/h) q_k(u); substitute into
  BaxOp and multiply by ν1[u]ν1[u-h]ν1[u-2h]ν1[u-3h] to clear denominators, giving a polynomial
  identity in u linear in the coeffs of q_k. Each q=u^j is one column; NullSpace (Chop first) gives
  the coeff vector (1-dim), monic-normalized by the top coefficient.*)
ClearAll[QSolve4];
QSolve4::degnull = "Q_`5` nullspace dim `6` (expected 1) at (H1,H2,H3,n)=(`1`,`2`,`3`,`4`).";
QSolve4::lastzero = "Q_`5` top coeff ~0 at (H1,H2,H3,n)=(`1`,`2`,`3`,`4`).";
QSolve4[λ1_, λ2_][H1_, H2_, H3_, n_] := QSolve4[λ1, λ2][H1, H2, H3, n] = Module[
   {Ms, u, T1, T2, T3, QD, N1, N1m, N1mm, N1mmm, sols},
   Ms = Mdeg[λ1, λ2][H1, H2, H3, n];
   T1 = tau1[λ1, λ2][H1, H2, H3, n][u]; T2 = tau2[λ1, λ2][H1, H2, H3, n][u]; T3 = tau3[λ1, λ2][H1, H2, H3, n][u];
   QD = qdetT[λ1, λ2][u];
   N1 = ν1[λ1, λ2][u]; N1m = ν1[λ1, λ2][u - h]; N1mm = ν1[λ1, λ2][u - 2 h]; N1mmm = ν1[λ1, λ2][u - 3 h];
   sols = Table[
     Module[{M = Ms[[k]], zk = z[k], cols, mat, ns, vec, last, maxdeg},
       maxdeg = M + 8;
       cols = Table[
          PadRight[CoefficientList[Expand[
             zk (u + h)^j N1 N1m N1mm N1mmm
              - T1 u^j N1m N1mm N1mmm
              + T2/zk (u - h)^j N1mm N1mmm
              - T3/zk^2 (u - 2 h)^j N1mmm
              + QD/zk^3 (u - 3 h)^j], u], maxdeg + 1],
          {j, 0, M}];
       mat = Transpose[cols];
       ns = NullSpace[Chop[mat, 10^-12], Tolerance -> 10^-10];
       If[Length[ns] != 1, Message[QSolve4::degnull, H1, H2, H3, n, k, Length[ns]]; Abort[]];
       vec = ns[[1]]; last = vec[[-1]];
       If[Abs[last] < 10^-13 Max[Abs[vec]], Message[QSolve4::lastzero, H1, H2, H3, n, k]; Abort[]];
       vec/last],
     {k, 1, 4}];
   sols
];
```
**Output:**
```
(no output)
```

## 2026-07-19 15:58:30
**Input:**
```wolfram
(*Q accessor: Q_k(u)=z[k]^(u/h) q_k(u), q_k the monic coeff list from QSolve4.*)
ClearAll[Qfun];
Qfun[λ1_, λ2_][H1_, H2_, H3_, n_][k_][u_] := z[k]^(u/h) qpoly[QSolve4[λ1, λ2][H1, H2, H3, n][[k]], u];
```
**Output:**
```
(no output)
```
## 2026-07-19 15:58:36 — 📥 BULK INSERT 1 CELL(S) at positions 33–33
1. [code] (*single-state TQ-residual check on a non-symmetric rep: each Q_k must be annihilated by BaxOp↵  at


## 2026-07-19 15:58:36
**Input:**
```wolfram
(*single-state TQ-residual check on a non-symmetric rep: each Q_k must be annihilated by BaxOp
  at several generic numeric u. State = the first key of TauEigensystem4[{2,1,1,0},{1,0,0,0}].*)
Module[{λ1 = {2, 1, 1, 0}, λ2 = {1, 0, 0, 0}, key, H1, H2, H3, n, res},
  key = First[Keys[TauEigensystem4[λ1, λ2]]];
  {H1, H2, H3, n} = key;
  res = Table[
     Max@Table[
        Abs@Chop[N[BaxOp[λ1, λ2][H1, H2, H3, n][Function[uu, Qfun[λ1, λ2][H1, H2, H3, n][k][uu]]][u0], 30]],
        {u0, {2/3, 9/5, 4, 23/7}}],
     {k, 1, 4}];
  If[Max[res] > 10^-10, Print["TQ RESIDUAL FAIL: ", res]; Abort[]];
  {"state" -> key, "perQresidual" -> res}
]
```
**Output:**
```
Out[89]= \left\{\text{state}\to \left\{0,1,1,0\right\},\text{perQresidual}\to \left\{0,0,0,0.\times 10^{-9}\right\}\right\}
```
## 2026-07-19 15:59:32 — 📥 BULK INSERT 1 CELL(S) at positions 34–34
1. [code] (*full sweep verification: over the 5 gl(4) rep pairs and every {H1,H2,H3,n} state, check↵  (a) per-


## 2026-07-19 15:59:33
**Input:**
```wolfram
(*full sweep verification: over the 5 gl(4) rep pairs and every {H1,H2,H3,n} state, check
  (a) per-Q_k TQ residual and (b) the Casoratian reconstruction — the four Q_k span the 4th-order
  TQ solution space, so A1..A4 rebuilt from 4x4 minors over the five shift rows {u+h,u,u-h,u-2h,u-3h}
  must match tau1/tau2/tau3/qdetT over the appropriate ν1 products. Hard-asserts <10^-9.*)
Module[{sweep = {{{1,0,0,0},{1,0,0,0}},{{1,0,0,0},{1,1,1,0}},{{1,1,0,0},{1,0,0,0}},{{1,1,0,0},{1,1,0,0}},{{2,1,1,0},{1,0,0,0}}}, report},
  report = Table[
    Module[{λ1 = rep[[1]], λ2 = rep[[2]], es, worstTQ = 0, worstCas = 0},
      es = TauEigensystem4[λ1, λ2];
      KeyValueMap[
       Function[{key, rec},
        Module[{H1 = key[[1]], H2 = key[[2]], H3 = key[[3]], n = key[[4]], Qs, u0, rows, minors,
          A1rec, A2rec, A3rec, A4rec, A1, A2, A3, A4, N1, N1m, N1mm, N1mmm},
         Do[worstTQ = Max[worstTQ,
             Abs@Chop[N[BaxOp[λ1, λ2][H1, H2, H3, n][
                 Function[uu, Qfun[λ1, λ2][H1, H2, H3, n][k][uu]]][u0], 30]]],
           {k, 1, 4}, {u0, {5/4, 11/3}}];
         u0 = 13/5;
         Qs[uu_] := Table[Qfun[λ1, λ2][H1, H2, H3, n][k][uu], {k, 1, 4}];
         rows = {Qs[u0 + h], Qs[u0], Qs[u0 - h], Qs[u0 - 2 h], Qs[u0 - 3 h]}; (* 5 x 4 *)
         minors = Table[Det[Delete[rows, i]], {i, 1, 5}];
         A1rec = minors[[2]]/minors[[1]]; A2rec = minors[[3]]/minors[[1]];
         A3rec = minors[[4]]/minors[[1]]; A4rec = minors[[5]]/minors[[1]];
         N1 = ν1[λ1, λ2][u0]; N1m = ν1[λ1, λ2][u0 - h]; N1mm = ν1[λ1, λ2][u0 - 2 h]; N1mmm = ν1[λ1, λ2][u0 - 3 h];
         A1 = tau1[λ1, λ2][H1, H2, H3, n][u0]/N1;
         A2 = tau2[λ1, λ2][H1, H2, H3, n][u0]/(N1 N1m);
         A3 = tau3[λ1, λ2][H1, H2, H3, n][u0]/(N1 N1m N1mm);
         A4 = qdetT[λ1, λ2][u0]/(N1 N1m N1mm N1mmm);
         worstCas = Max[worstCas,
            Abs@Chop[N[A1rec - A1, 30]], Abs@Chop[N[A2rec - A2, 30]],
            Abs@Chop[N[A3rec - A3, 30]], Abs@Chop[N[A4rec - A4, 30]]]]],
       es];
      {rep, "nstates" -> Length[es], "worstTQ" -> worstTQ, "worstCasoratian" -> worstCas}],
    {rep, sweep}];
  If[Max[#[[3, 2]] & /@ report, #[[4, 2]] & /@ report] > 10^-9,
    Print["SWEEP FAIL: ", report]; Abort[]];
  Column[report]
]
```
**Output:**
```
Print: SWEEP FAIL:[formula] | Out[90]= \text{\$Aborted}
```

## 2026-07-19 15:59:50
**Input:**
```wolfram
Module[{sweep = {{{1,0,0,0},{1,0,0,0}},{{1,0,0,0},{1,1,1,0}},{{1,1,0,0},{1,0,0,0}},{{1,1,0,0},{1,1,0,0}},{{2,1,1,0},{1,0,0,0}}}, report},
  report = Table[
    Module[{λ1 = rep[[1]], λ2 = rep[[2]], es, worstTQ = 0, worstCas = 0},
      es = TauEigensystem4[λ1, λ2];
      KeyValueMap[
       Function[{key, rec},
        Module[{H1 = key[[1]], H2 = key[[2]], H3 = key[[3]], n = key[[4]], Qs, u0, rows, minors,
          A1rec, A2rec, A3rec, A4rec, A1, A2, A3, A4, N1, N1m, N1mm, N1mmm},
         Do[worstTQ = Max[worstTQ,
             Abs@Chop[N[BaxOp[λ1, λ2][H1, H2, H3, n][
                 Function[uu, Qfun[λ1, λ2][H1, H2, H3, n][k][uu]]][u0], 30]]],
           {k, 1, 4}, {u0, {5/4, 11/3}}];
         u0 = 13/5;
         Qs[uu_] := Table[Qfun[λ1, λ2][H1, H2, H3, n][k][uu], {k, 1, 4}];
         rows = {Qs[u0 + h], Qs[u0], Qs[u0 - h], Qs[u0 - 2 h], Qs[u0 - 3 h]};
         minors = Table[Det[Delete[rows, i]], {i, 1, 5}];
         A1rec = minors[[2]]/minors[[1]]; A2rec = minors[[3]]/minors[[1]];
         A3rec = minors[[4]]/minors[[1]]; A4rec = minors[[5]]/minors[[1]];
         N1 = ν1[λ1, λ2][u0]; N1m = ν1[λ1, λ2][u0 - h]; N1mm = ν1[λ1, λ2][u0 - 2 h]; N1mmm = ν1[λ1, λ2][u0 - 3 h];
         A1 = tau1[λ1, λ2][H1, H2, H3, n][u0]/N1;
         A2 = tau2[λ1, λ2][H1, H2, H3, n][u0]/(N1 N1m);
         A3 = tau3[λ1, λ2][H1, H2, H3, n][u0]/(N1 N1m N1mm);
         A4 = qdetT[λ1, λ2][u0]/(N1 N1m N1mm N1mmm);
         worstCas = Max[worstCas,
            Abs@Chop[N[A1rec - A1, 30]], Abs@Chop[N[A2rec - A2, 30]],
            Abs@Chop[N[A3rec - A3, 30]], Abs@Chop[N[A4rec - A4, 30]]]]],
       es];
      {rep, "nstates" -> Length[es], "worstTQ" -> worstTQ, "worstCasoratian" -> worstCas}],
    {rep, sweep}];
  report
]
```
**Output:**
```
Out= {{{{1, 0, 0, 0}, {1, 0, 0, 0}}, "nstates" -> 16, "worstTQ" -> 0, "worstCasoratian" -> 5.8240859179529194397295760861`5.41669690473683*^-10}, {{{1, 0, 0, 0}, {1, 1, 1, 0}}, "nstates" -> 16, "worstTQ" -> 0, "worstCasoratian" -> 0}, {{{1, 1, 0, 0}, {1, 0, 0, 0}}, "nstates" -> 24, "worstTQ" -> 0, "worstCasoratian" -> 7.4622301597584542807442797623`4.677218790255995*^-10}, {{{1, 1, 0, 0}, {1, 1, 0, 0}}, "nstates" -> 36, "worstTQ" -> 0, "worstCasoratian" -> 0}, {{{2, 1, 1, 0}, {1, 0, 0, 0}}, "nstates" -> 60, "worstTQ" -> 0``7.900270818900865, "worstCasoratian" -> 4.54729176095814416126440874593`5.538956360885462*^-8}}
```

## 2026-07-19 16:00:15
**Input:**
```wolfram
Module[{λ1={2,1,1,0},λ2={1,0,0,0},es,worstKey,worstVal=0},
  es=TauEigensystem4[λ1,λ2];
  KeyValueMap[Function[{key,rec},
    Module[{H1=key[[1]],H2=key[[2]],H3=key[[3]],n=key[[4]],Qs,u0,rows,minors,A1rec,A2rec,A3rec,A4rec,A1,A2,A3,A4,N1,N1m,N1mm,N1mmm,cas},
      u0=13/5;
      Qs[uu_]:=Table[Qfun[λ1,λ2][H1,H2,H3,n][k][uu],{k,1,4}];
      rows={Qs[u0+h],Qs[u0],Qs[u0-h],Qs[u0-2h],Qs[u0-3h]};
      minors=Table[Det[Delete[rows,i]],{i,1,5}];
      A1rec=minors[[2]]/minors[[1]];A2rec=minors[[3]]/minors[[1]];A3rec=minors[[4]]/minors[[1]];A4rec=minors[[5]]/minors[[1]];
      N1=ν1[λ1,λ2][u0];N1m=ν1[λ1,λ2][u0-h];N1mm=ν1[λ1,λ2][u0-2h];N1mmm=ν1[λ1,λ2][u0-3h];
      A1=tau1[λ1,λ2][H1,H2,H3,n][u0]/N1;A2=tau2[λ1,λ2][H1,H2,H3,n][u0]/(N1 N1m);
      A3=tau3[λ1,λ2][H1,H2,H3,n][u0]/(N1 N1m N1mm);A4=qdetT[λ1,λ2][u0]/(N1 N1m N1mm N1mmm);
      cas=Max[Abs@Chop[N[A1rec-A1,30]],Abs@Chop[N[A2rec-A2,30]],Abs@Chop[N[A3rec-A3,30]],Abs@Chop[N[A4rec-A4,30]]];
      If[cas>worstVal,worstVal=cas;worstKey=key]],
    es];
  {worstKey,worstVal,"Mdeg"->Mdeg@@Join[{λ1,λ2}][Sequence@@worstKey]}
]
```
**Output:**
```
[message] ToExpression::sntx: Invalid syntax in or before "]".
(no output)
```

## 2026-07-19 16:00:23
**Input:**
```wolfram
Module[{λ1={2,1,1,0},λ2={1,0,0,0},es,worstKey,worstVal=0},
  es=TauEigensystem4[λ1,λ2];
  KeyValueMap[Function[{key,rec},
    Module[{H1=key[[1]],H2=key[[2]],H3=key[[3]],n=key[[4]],Qs,u0,rows,minors,A1rec,A2rec,A3rec,A4rec,A1,A2,A3,A4,N1,N1m,N1mm,N1mmm,cas},
      u0=13/5;
      Qs[uu_]:=Table[Qfun[λ1,λ2][H1,H2,H3,n][k][uu],{k,1,4}];
      rows={Qs[u0+h],Qs[u0],Qs[u0-h],Qs[u0-2h],Qs[u0-3h]};
      minors=Table[Det[Delete[rows,i]],{i,1,5}];
      A1rec=minors[[2]]/minors[[1]];A2rec=minors[[3]]/minors[[1]];A3rec=minors[[4]]/minors[[1]];A4rec=minors[[5]]/minors[[1]];
      N1=ν1[λ1,λ2][u0];N1m=ν1[λ1,λ2][u0-h];N1mm=ν1[λ1,λ2][u0-2h];N1mmm=ν1[λ1,λ2][u0-3h];
      A1=tau1[λ1,λ2][H1,H2,H3,n][u0]/N1;A2=tau2[λ1,λ2][H1,H2,H3,n][u0]/(N1 N1m);
      A3=tau3[λ1,λ2][H1,H2,H3,n][u0]/(N1 N1m N1mm);A4=qdetT[λ1,λ2][u0]/(N1 N1m N1mm N1mmm);
      cas=Max[Abs@Chop[N[A1rec-A1,30]],Abs@Chop[N[A2rec-A2,30]],Abs@Chop[N[A3rec-A3,30]],Abs@Chop[N[A4rec-A4,30]]];
      If[cas>worstVal,worstVal=cas;worstKey=key]],
    es];
  {worstKey,worstVal,Mdeg[λ1,λ2][worstKey[[1]],worstKey[[2]],worstKey[[3]],worstKey[[4]]]}
]
```
**Output:**
```
[message] ToExpression::sntx: Invalid syntax in or before "]".
(no output)
```

## 2026-07-19 16:00:55
**Input:**
```wolfram
Module[{λ1={2,1,1,0},λ2={1,0,0,0},es,worstKey,worstVal=0},
  es=TauEigensystem4[λ1,λ2];
  KeyValueMap[Function[{key,rec},
    Module[{H1=key[[1]],H2=key[[2]],H3=key[[3]],n=key[[4]],Qs,u0,rows,minors,A1rec,A2rec,A3rec,A4rec,A1,A2,A3,A4,N1,N1m,N1mm,N1mmm,cas},
      u0=13/5;
      Qs[uu_]:=Table[Qfun[λ1,λ2][H1,H2,H3,n][k][uu],{k,1,4}];
      rows={Qs[u0+h],Qs[u0],Qs[u0-h],Qs[u0-2h],Qs[u0-3h]};
      minors=Table[Det[Delete[rows,i]],{i,1,5}];
      A1rec=minors[[2]]/minors[[1]];A2rec=minors[[3]]/minors[[1]];A3rec=minors[[4]]/minors[[1]];A4rec=minors[[5]]/minors[[1]];
      N1=ν1[λ1,λ2][u0];N1m=ν1[λ1,λ2][u0-h];N1mm=ν1[λ1,λ2][u0-2h];N1mmm=ν1[λ1,λ2][u0-3h];
      A1=tau1[λ1,λ2][H1,H2,H3,n][u0]/N1;A2=tau2[λ1,λ2][H1,H2,H3,n][u0]/(N1 N1m);
      A3=tau3[λ1,λ2][H1,H2,H3,n][u0]/(N1 N1m N1mm);A4=qdetT[λ1,λ2][u0]/(N1 N1m N1mm N1mmm);
      cas=Max[Abs@Chop[N[A1rec-A1,30]],Abs@Chop[N[A2rec-A2,30]],Abs@Chop[N[A3rec-A3,30]],Abs@Chop[N[A4rec-A4,30]]];
      If[cas>worstVal,worstVal=cas;worstKey=key]]],
    es];
  {worstKey,worstVal,Mdeg[λ1,λ2][worstKey[[1]],worstKey[[2]],worstKey[[3]],worstKey[[4]]]}
]
```
**Output:**
```
Out= {{1, 1, -2, 0}, 4.54729176095814416126440874593`5.538956360885462*^-8, {1, 2, 3, 1}}
```

## 2026-07-19 16:01:33
**Input:**
```wolfram
Module[{λ1={2,1,1,0},λ2={1,0,0,0},H1=1,H2=1,H3=-2,n=0,Ms,u,T1,T2,T3,QD,N1,N1m,N1mm,N1mmm,results},
  Ms=Mdeg[λ1,λ2][H1,H2,H3,n];
  T1=tau1[λ1,λ2][H1,H2,H3,n][u];T2=tau2[λ1,λ2][H1,H2,H3,n][u];T3=tau3[λ1,λ2][H1,H2,H3,n][u];
  QD=qdetT[λ1,λ2][u];
  N1=ν1[λ1,λ2][u];N1m=ν1[λ1,λ2][u-h];N1mm=ν1[λ1,λ2][u-2h];N1mmm=ν1[λ1,λ2][u-3h];
  results=Table[
    Module[{M=Ms[[k]],zk=z[k],cols,mat,sv,maxdeg},
      maxdeg=M+8;
      cols=Table[
        PadRight[CoefficientList[Expand[
          zk (u+h)^j N1 N1m N1mm N1mmm
           -T1 u^j N1m N1mm N1mmm
           +T2/zk (u-h)^j N1mm N1mmm
           -T3/zk^2 (u-2h)^j N1mmm
           +QD/zk^3 (u-2h)^j],u],maxdeg+1],
        {j,0,M}];
      mat=Transpose[cols];
      sv=SingularValueList[N[mat,30]];
      {k,"M"->M,"matDims"->Dimensions[mat],"singularValues"->sv}],
    {k,1,4}];
  results
]
```
**Output:**
```
Out= {{1, "M" -> 1, "matDims" -> {10, 2}, "singularValues" -> {13689.70026028494945209356471873085957924205`20.086420382419618, 1341.99411232998269793867912872873070861954`20.086420382419618}}, {2, "M" -> 2, "matDims" -> {11, 3}, "singularValues" -> {18369.50548672642258727263877907744084948799`18.770847470902318, 12209.55653656029574430251358571800740279002`18.770847470902318, 182.19692703725273814575635121075962375767`18.770847470902318}}, {3, "M" -> 3, "matDims" -> {12, 4}, "singularValues" -> {74062.2199592174100619311762803898062648534`18.271874202814583, 17192.00522827988340531834695746358610437871`18.271874202814583, 3505.86120029091766058318435671514679428957`18.271874202814583, 56.62864356477385411164053762504604863671`18.271874202814583}}, {4, "M" -> 1, "matDims" -> {10, 2}, "singularValues" -> {13779.25987203180990332216268456821459724414`17.895048385405733, 1329.59547932181320127494309506074089750542`17.895048385405733}}}
```
## 2026-07-19 16:03:31 — ✏️ EDIT CELL 7
(*Numeric parameters (shared across Section B). θ,h as in the gl(3) notebook;
  four DISTINCT unit-modulus twist eigenvalues z[1..4] (generic → non-degenerate).
  χ1..χ4 were defined above as symmetri…


---

## 2026-07-19 16:04:01 — 🔄 KERNEL RESTART

## 2026-07-19 16:04:31 — ▶️ RUN CELLS 1–34
Cell 1: ✓ — (no output)
Cell 2: ✓ — (no output)
Cell 3: ✓ — (no output)
Cell 4: ✓ — (no output)
Cell 5: ✓ — Out[16]= \left\{0\right\}
Cell 6: ✓ — (no output)
Cell 7: ✓ — Out[28]= \left\{\frac{1}{3},\frac{1}{7},1,0.36043987514722653726667765059555009129408114679289+0.93278244859337471648483598207086792945772729283224\,\mathrm{i},0.50886738099049982816332137226706861990909726908446+0.86084492712907909318864439361616424581514277367848\,\mathrm{i},0.53325788058959937259583234509312429536160862774512+0.84595273673479452648338161845314258986220080613203\,\mathrm{i},0.53… [985 chars total]
Cell 8: ✓ — (no output)
Cell 9: ✓ — Out[32]= \left\{0.\times 10^{-49},0.00634087393306325790235269453099925019937299492747\right\}
Cell 10: ✓ — (no output)
Cell 11: ✓ — Out[36]= \left\{0\right\}
Cell 12: ✓ — (no output)
Cell 13: ✓ — (no output)
Cell 14: ✓ — (no output)
Cell 15: ✓ — Out[44]= \left\{0\right\}
Cell 16: ✓ — (no output)
Cell 17: ✓ — (no output)
Cell 18: ✓ — Out[51]= \left\{\text{C2=$\chi$1 Id OK},\text{D4=$\chi$2 Id OK},\text{G6=$\chi$3 Id OK},\text{maxCommutator}\to 0.\times 10^{-43}\right\}
Cell 19: ✓ — (no output)
Cell 20: ✓ — Out[55]= \left\{\text{all reps sum-rule OK},\text{adjoint$\otimes$fund @(H1,H2,H3)=(1,0,0) mult}\to 6\right\}
Cell 21: ✓ — (no output)
Cell 22: ✓ — (no output)
Cell 23: ✓ — Out[68]= \left\{\text{nstates}\to 60,\text{worstResidual@u=17/5}\to 8.595834409382585934859838\times 10^{-26}\right\}
Cell 24: ✓ — (no output)
Cell 25: ✓ — (no output)
Cell 26: ✓ — (no output)
Cell 27: ✓ — (no output)
Cell 28: ✓ — (no output)
Cell 29: ✓ — Out[80]= \left\{\text{closed-form M\_k=$\Lambda$\_1-n\_k matches Mdeg over full sweep, worst}\to 0\right\}
Cell 30: ✓ — (no output)
Cell 31: ✓ — (no output)
Cell 32: ✓ — (no output)
Cell 33: ✓ — Out[89]= \left\{\text{state}\to \left\{0,1,1,0\right\},\text{perQresidual}\to \left\{0,0,0,0\right\}\right\}
Cell 34: ✓ — Print: SWEEP FAIL:[formula] | Out[90]= \text{\$Aborted}


## 2026-07-19 16:04:45
**Input:**
```wolfram
Module[{sweep = {{{1,0,0,0},{1,0,0,0}},{{1,0,0,0},{1,1,1,0}},{{1,1,0,0},{1,0,0,0}},{{1,1,0,0},{1,1,0,0}},{{2,1,1,0},{1,0,0,0}}}, report},
  report = Table[
    Module[{λ1 = rep[[1]], λ2 = rep[[2]], es, worstTQ = 0, worstCas = 0},
      es = TauEigensystem4[λ1, λ2];
      KeyValueMap[
       Function[{key, rec},
        Module[{H1 = key[[1]], H2 = key[[2]], H3 = key[[3]], n = key[[4]], Qs, u0, rows, minors,
          A1rec, A2rec, A3rec, A4rec, A1, A2, A3, A4, N1, N1m, N1mm, N1mmm},
         Do[worstTQ = Max[worstTQ,
             Abs@Chop[N[BaxOp[λ1, λ2][H1, H2, H3, n][
                 Function[uu, Qfun[λ1, λ2][H1, H2, H3, n][k][uu]]][u0], 30]]],
           {k, 1, 4}, {u0, {5/4, 11/3}}];
         u0 = 13/5;
         Qs[uu_] := Table[Qfun[λ1, λ2][H1, H2, H3, n][k][uu], {k, 1, 4}];
         rows = {Qs[u0 + h], Qs[u0], Qs[u0 - h], Qs[u0 - 2 h], Qs[u0 - 3 h]};
         minors = Table[Det[Delete[rows, i]], {i, 1, 5}];
         A1rec = minors[[2]]/minors[[1]]; A2rec = minors[[3]]/minors[[1]];
         A3rec = minors[[4]]/minors[[1]]; A4rec = minors[[5]]/minors[[1]];
         N1 = ν1[λ1, λ2][u0]; N1m = ν1[λ1, λ2][u0 - h]; N1mm = ν1[λ1, λ2][u0 - 2 h]; N1mmm = ν1[λ1, λ2][u0 - 3 h];
         A1 = tau1[λ1, λ2][H1, H2, H3, n][u0]/N1;
         A2 = tau2[λ1, λ2][H1, H2, H3, n][u0]/(N1 N1m);
         A3 = tau3[λ1, λ2][H1, H2, H3, n][u0]/(N1 N1m N1mm);
         A4 = qdetT[λ1, λ2][u0]/(N1 N1m N1mm N1mmm);
         worstCas = Max[worstCas,
            Abs@Chop[N[A1rec - A1, 30]], Abs@Chop[N[A2rec - A2, 30]],
            Abs@Chop[N[A3rec - A3, 30]], Abs@Chop[N[A4rec - A4, 30]]]]],
       es];
      {rep, "nstates" -> Length[es], "worstTQ" -> worstTQ, "worstCasoratian" -> worstCas}],
    {rep, sweep}];
  report
]
```
**Output:**
```
Out= {{{{1, 0, 0, 0}, {1, 0, 0, 0}}, "nstates" -> 16, "worstTQ" -> 0, "worstCasoratian" -> 0}, {{{1, 0, 0, 0}, {1, 1, 1, 0}}, "nstates" -> 16, "worstTQ" -> 0, "worstCasoratian" -> 0}, {{{1, 1, 0, 0}, {1, 0, 0, 0}}, "nstates" -> 24, "worstTQ" -> 0, "worstCasoratian" -> 0}, {{{1, 1, 0, 0}, {1, 1, 0, 0}}, "nstates" -> 36, "worstTQ" -> 0, "worstCasoratian" -> 0}, {{{2, 1, 1, 0}, {1, 0, 0, 0}}, "nstates" -> 60, "worstTQ" -> 0``9.997211473009042, "worstCasoratian" -> 1.64569535463206613073540232745`7.901462141009449*^-9}}
```
## 2026-07-19 16:05:05 — ✏️ EDIT CELL 7
(*Numeric parameters (shared across Section B). θ,h as in the gl(3) notebook;
  four DISTINCT unit-modulus twist eigenvalues z[1..4] (generic → non-degenerate).
  χ1..χ4 were defined above as symmetri…


---

## 2026-07-19 16:05:09 — 🔄 KERNEL RESTART

## 2026-07-19 16:05:40 — ▶️ RUN CELLS 1–34
Cell 1: ✓ — (no output)
Cell 2: ✓ — (no output)
Cell 3: ✓ — (no output)
Cell 4: ✓ — (no output)
Cell 5: ✓ — Out[16]= \left\{0\right\}
Cell 6: ✓ — (no output)
Cell 7: ✓ — Out[28]= \left\{\frac{1}{3},\frac{1}{7},1,0.3604398751472265372666776505955500912940811467928913386732181476418264+0.9327824485933747164848359820708679294577272928322431080055758094332927\,\mathrm{i},0.5088673809904998281633213722670686199090972690844597081081963623607767+0.8608449271290790931886443936161642458151427736784766273900288346856365\,\mathrm{i},0.5332578805895993725958323450931242953616… [1246 chars total]
Cell 8: ✓ — (no output)
Cell 9: ✓ — Out[32]= \left\{0.\times 10^{-69},0.0063408739330632579023526945309992501993729949274710916715470375206530\right\}
Cell 10: ✓ — (no output)
Cell 11: ✓ — Out[36]= \left\{0\right\}
Cell 12: ✓ — (no output)
Cell 13: ✓ — (no output)
Cell 14: ✓ — (no output)
Cell 15: ✓ — Out[44]= \left\{0\right\}
Cell 16: ✓ — (no output)
Cell 17: ✓ — (no output)
Cell 18: ✓ — Out[51]= \left\{\text{C2=$\chi$1 Id OK},\text{D4=$\chi$2 Id OK},\text{G6=$\chi$3 Id OK},\text{maxCommutator}\to 0.\times 10^{-63}\right\}
Cell 19: ✓ — (no output)
Cell 20: ✓ — Out[55]= \left\{\text{all reps sum-rule OK},\text{adjoint$\otimes$fund @(H1,H2,H3)=(1,0,0) mult}\to 6\right\}
Cell 21: ✓ — (no output)
Cell 22: ✓ — (no output)
Cell 23: ✓ — Out[68]= \left\{\text{nstates}\to 60,\text{worstResidual@u=17/5}\to 8.595834409382585934859838\times 10^{-26}\right\}
Cell 24: ✓ — (no output)
Cell 25: ✓ — (no output)
Cell 26: ✓ — (no output)
Cell 27: ✓ — (no output)
Cell 28: ✓ — (no output)
Cell 29: ✓ — Out[80]= \left\{\text{closed-form M\_k=$\Lambda$\_1-n\_k matches Mdeg over full sweep, worst}\to 0\right\}
Cell 30: ✓ — (no output)
Cell 31: ✓ — (no output)
Cell 32: ✓ — (no output)
Cell 33: ✓ — Out[89]= \left\{\text{state}\to \left\{0,1,1,0\right\},\text{perQresidual}\to \left\{0,0,0,0\right\}\right\}
Cell 34: ✓ — Print: SWEEP FAIL:[formula] | Out[90]= \text{\$Aborted}


## 2026-07-19 16:05:50
**Input:**
```wolfram
Module[{sweep = {{{2,1,1,0},{1,0,0,0}}}, report},
  report = Table[
    Module[{λ1 = rep[[1]], λ2 = rep[[2]], es, worstTQ = 0, worstCas = 0, worstKey},
      es = TauEigensystem4[λ1, λ2];
      KeyValueMap[
       Function[{key, rec},
        Module[{H1 = key[[1]], H2 = key[[2]], H3 = key[[3]], n = key[[4]], Qs, u0, rows, minors,
          A1rec, A2rec, A3rec, A4rec, A1, A2, A3, A4, N1, N1m, N1mm, N1mmm, thisCas},
         u0 = 13/5;
         Qs[uu_] := Table[Qfun[λ1, λ2][H1, H2, H3, n][k][uu], {k, 1, 4}];
         rows = {Qs[u0 + h], Qs[u0], Qs[u0 - h], Qs[u0 - 2 h], Qs[u0 - 3 h]};
         minors = Table[Det[Delete[rows, i]], {i, 1, 5}];
         A1rec = minors[[2]]/minors[[1]]; A2rec = minors[[3]]/minors[[1]];
         A3rec = minors[[4]]/minors[[1]]; A4rec = minors[[5]]/minors[[1]];
         N1 = ν1[λ1, λ2][u0]; N1m = ν1[λ1, λ2][u0 - h]; N1mm = ν1[λ1, λ2][u0 - 2 h]; N1mmm = ν1[λ1, λ2][u0 - 3 h];
         A1 = tau1[λ1, λ2][H1, H2, H3, n][u0]/N1;
         A2 = tau2[λ1, λ2][H1, H2, H3, n][u0]/(N1 N1m);
         A3 = tau3[λ1, λ2][H1, H2, H3, n][u0]/(N1 N1m N1mm);
         A4 = qdetT[λ1, λ2][u0]/(N1 N1m N1mm N1mmm);
         thisCas = Max[
            Abs@Chop[N[A1rec - A1, 30]], Abs@Chop[N[A2rec - A2, 30]],
            Abs@Chop[N[A3rec - A3, 30]], Abs@Chop[N[A4rec - A4, 30]]];
         If[thisCas > worstCas, worstCas = thisCas; worstKey = key]]],
       es];
      {rep, "nstates" -> Length[es], "worstCasoratian" -> worstCas, "worstKey" -> worstKey}],
    {rep, sweep}];
  report
]
```
**Output:**
```
Out= {{{{2, 1, 1, 0}, {1, 0, 0, 0}}, "nstates" -> 60, "worstCasoratian" -> 1.64569535463206613073540232745`7.901462141009449*^-9, "worstKey" -> {1, -1, 2, 1}}}
```

## 2026-07-19 16:06:52
**Input:**
```wolfram
Module[{λ1={2,1,1,0},λ2={1,0,0,0},H1=1,H2=-1,H3=2,n=1,Ms,u,T1,T2,T3,QD,N1,N1m,N1mm,N1mmm,results},
  Ms=Mdeg[λ1,λ2][H1,H2,H3,n];
  T1=tau1[λ1,λ2][H1,H2,H3,n][u];T2=tau2[λ1,λ2][H1,H2,H3,n][u];T3=tau3[λ1,λ2][H1,H2,H3,n][u];
  QD=qdetT[λ1,λ2][u];
  N1=ν1[λ1,λ2][u];N1m=ν1[λ1,λ2][u-h];N1mm=ν1[λ1,λ2][u-2h];N1mmm=ν1[λ1,λ2][u-3h];
  results=Table[
    Module[{M=Ms[[k]],zk=z[k],cols,mat,sv,ns1,ns2,maxdeg},
      maxdeg=M+8;
      cols=Table[
        PadRight[CoefficientList[Expand[
          zk (u+h)^j N1 N1m N1mm N1mmm
           -T1 u^j N1m N1mm N1mmm
           +T2/zk (u-h)^j N1mm N1mmm
           -T3/zk^2 (u-2h)^j N1mmm
           +QD/zk^3 (u-3h)^j],u],maxdeg+1],
        {j,0,M}];
      mat=Transpose[cols];
      sv=SingularValueList[N[mat,70]];
      ns1=NullSpace[Chop[mat,10^-12],Tolerance->10^-10];
      ns2=NullSpace[Chop[mat,10^-40],Tolerance->10^-40];
      {k,"M"->M,"singularValues"->sv,"nsLen@1e-10"->Length[ns1],"nsLen@1e-40"->Length[ns2]}],
    {k,1,4}];
  results
]
```
**Output:**
```
Out= {{1, "M" -> 1, "singularValues" -> {14900.06742161557294614415843204442243118366`21.85452141859339}, "nsLen@1e-10" -> 1, "nsLen@1e-40" -> 1}, {2, "M" -> 2, "singularValues" -> {19057.96940773073194274431042546712964366532`20.23254802824842, 3849.58057311991834183815707804235573071178`20.23254802824842}, "nsLen@1e-10" -> 1, "nsLen@1e-40" -> 1}, {3, "M" -> 1, "singularValues" -> {14789.2232693170952858894428599215011876023`19.579944845523915}, "nsLen@1e-10" -> 1, "nsLen@1e-40" -> 1}, {4, "M" -> 3, "singularValues" -> {21300.83580184066833449577145154214019850028`19.642219122973987, 17859.06895338887394426278795567684179480323`19.642219122973987, 37.77143657247830505503037866700263391156`19.642219122973987}, "nsLen@1e-10" -> 1, "nsLen@1e-40" -> 1}}
```

## 2026-07-19 16:07:25
**Input:**
```wolfram
Module[{λ1={2,1,1,0},λ2={1,0,0,0},H1=1,H2=-1,H3=2,n=1,Qs,u0,rows,minors},
  u0=13/5;
  Qs[uu_]:=Table[Qfun[λ1,λ2][H1,H2,H3,n][k][uu],{k,1,4}];
  rows={Qs[u0+h],Qs[u0],Qs[u0-h],Qs[u0-2h],Qs[u0-3h]};
  minors=Table[Det[Delete[rows,i]],{i,1,5}];
  {"minors"->minors,"rowsMagnitudes"->Map[Abs,rows,{2}]}
]
```
**Output:**
```
Out= {"minors" -> {-2.3847446095165398230398906975`18.855077804766218*^-10 - 6.25410273047821656523411864572`18.855077804766218*^-9*I, 1.6496513302730163603734706726777`18.831913491450877*^-7 - 1.0537262304573770758716455682013`18.831913491450877*^-7*I, -7.177945962860499823137704018314`18.831913491450877*^-8 - 2.540442648544253480132862922671`18.831913491450877*^-8*I, 4.96593916442891827648670694008`18.831913491450877*^-9 - 2.2825759332424531856508019910446`18.831913491450877*^-7*I, -3.248953789858219147067969424034`18.831913491450877*^-8 + 1.779804292735941265726773346413`18.831913491450877*^-8*I}, "rowsMagnitudes" -> {{10.55017199799869926838421983568050187908`21.113250577955142, 432.80613534275608029961573199750796224081`19.508656505600154, 10.55017199799869926837291841063552391883`18.8324333479095, 570430.84458588824304660517623979885312433542`18.886179241432146}, {10.37218972365750463232750324506893285908`21.221322271262252, 425.47837999148010248678213430863963772658`19.553708401014365, 10.37218972365750463231600789247476884086`18.87082062883281, 560807.46196116125805163621757156159583353054`18.92428488919217}, {10.28875649141998201759348856670919221059`21.181918411265517, 422.05429496007381980502738441254351999598`19.615433957351126, 10.28875649141998201759348856670919221123`18.9669196214492, 556296.3550506744842664103433588893884998442`19.028742239307853}, {10.30216970428692242724905972254409571188`21.188780672626034, 422.60429865094646562468272493402565644012`19.594210372310375, 10.30216970428692242724905972254409571251`18.926623113855932, 557021.58001694456016091517395576992535148857`18.98607736833343}, {10.41205508495406351888135535480180672987`21.22318253168934, 427.11691860417880132197406437061069471823`19.612888297954292, 10.4120550849540635188813553548018067305`18.949088613691522, 562962.95376940859946350523608507662346314292`19.009370562557873}}}
```

## 2026-07-19 16:08:22
**Input:**
```wolfram
Module[{λ1={2,1,1,0},λ2={1,0,0,0},H1=1,H2=-1,H3=2,n=1,Ms,u,T1,T2,T3,QD,N1,N1m,N1mm,N1mmm,k=4,M,zk,cols,mat,ns1,ns2},
  Ms=Mdeg[λ1,λ2][H1,H2,H3,n];
  T1=tau1[λ1,λ2][H1,H2,H3,n][u];T2=tau2[λ1,λ2][H1,H2,H3,n][u];T3=tau3[λ1,λ2][H1,H2,H3,n][u];
  QD=qdetT[λ1,λ2][u];
  N1=ν1[λ1,λ2][u];N1m=ν1[λ1,λ2][u-h];N1mm=ν1[λ1,λ2][u-2h];N1mmm=ν1[λ1,λ2][u-3h];
  M=Ms[[k]];zk=z[k];
  cols=Table[
    PadRight[CoefficientList[Expand[
      zk (u+h)^j N1 N1m N1mm N1mmm
       -T1 u^j N1m N1mm N1mmm
       +T2/zk (u-h)^j N1mm N1mmm
       -T3/zk^2 (u-2h)^j N1mmm
       +QD/zk^3 (u-3h)^j],u],M+9],
    {j,0,M}];
  mat=Transpose[cols];
  ns1=NullSpace[Chop[mat,10^-12],Tolerance->10^-10][[1]];
  ns2=NullSpace[Chop[mat,10^-40],Tolerance->10^-40][[1]];
  ns1=ns1/ns1[[-1]];ns2=ns2/ns2[[-1]];
  {"ns1(tol1e-10,monic)"->ns1,"ns2(tol1e-40,monic)"->ns2,"maxDiff"->Max[Abs[ns1-ns2]]}
]
```
**Output:**
```
Out= {"ns1(tol1e-10,monic)" -> {69099.74113653566381572030265900526172299331`19.057799521221774 + 555688.02358043954954790484972764503721694001`19.057799521221774*I, -57591.88531764223756429621094865325861394836`19.057799521221778 + 648.06051282574130904052672131472786104798`17.739832787124694*I, -2.30952380952380952380952380951775453476`17.494657418213944 - 354.01105865284762861714030461095799893646`19.057799521221774*I, 1.`19.057799521221774 + 0``19.668788200251427*I}, "ns2(tol1e-40,monic)" -> {69099.74113653566381572030265900526172299331`19.040159131646025 + 555688.02358043954954790484972764503721694001`19.040159131646025*I, -57591.88531764223756429621094865325861394836`19.040159131646025 + 648.06051282574130904052672131472786104798`17.72219239754894*I, -2.30952380952380952380952380951775453476`17.47701702863819 - 354.01105865284762861714030461095799893646`19.040159131646025*I, 1.`19.040159131646025 + 0``19.651147810675678*I}, "maxDiff" -> 0``12.999696705212603}
```

## 2026-07-19 16:09:10
**Input:**
```wolfram
Module[{λ1={2,1,1,0},λ2={1,0,0,0},H1=1,H2=-1,H3=2,n=1,Ms,u,T1,T2,T3,QD,N1,N1m,N1mm,N1mmm,k=4,M,zk,cols,mat,nsNoChop,nsOld},
  Ms=Mdeg[λ1,λ2][H1,H2,H3,n];
  T1=tau1[λ1,λ2][H1,H2,H3,n][u];T2=tau2[λ1,λ2][H1,H2,H3,n][u];T3=tau3[λ1,λ2][H1,H2,H3,n][u];
  QD=qdetT[λ1,λ2][u];
  N1=ν1[λ1,λ2][u];N1m=ν1[λ1,λ2][u-h];N1mm=ν1[λ1,λ2][u-2h];N1mmm=ν1[λ1,λ2][u-3h];
  M=Ms[[k]];zk=z[k];
  cols=Table[
    PadRight[CoefficientList[Expand[
      zk (u+h)^j N1 N1m N1mm N1mmm
       -T1 u^j N1m N1mm N1mmm
       +T2/zk (u-h)^j N1mm N1mmm
       -T3/zk^2 (u-2h)^j N1mmm
       +QD/zk^3 (u-3h)^j],u],M+9],
    {j,0,M}];
  mat=Transpose[cols];
  nsNoChop=NullSpace[mat][[1]];
  nsOld=NullSpace[Chop[mat,10^-12],Tolerance->10^-10][[1]];
  nsNoChop=nsNoChop/nsNoChop[[-1]];nsOld=nsOld/nsOld[[-1]];
  {"Precision[mat]entries"->Precision[mat[[1,1]]],"nsNoChop"->nsNoChop,"maxDiffVsOld"->Max[Abs[nsNoChop-nsOld]],"PrecisionNoChop"->Precision[nsNoChop[[1]]]}
]
```
**Output:**
```
Out= {"Precision[mat]entries" -> 25.433106600509664, "nsNoChop" -> {69099.74113653566381572030265900526172299331`19.040159131646025 + 555688.02358043954954790484972764503721694001`19.040159131646025*I, -57591.88531764223756429621094865325861394836`19.040159131646025 + 648.06051282574130904052672131472786104798`17.72219239754894*I, -2.30952380952380952380952380951775453476`17.47701702863819 - 354.01105865284762861714030461095799893646`19.040159131646025*I, 1.`19.040159131646025 + 0``19.651147810675678*I}, "maxDiffVsOld" -> 0``12.999696705212603, "PrecisionNoChop" -> 19.04015913164602}
```
## 2026-07-19 16:10:05 — ✏️ EDIT CELL 21
ClearAll[TauEigensystem4];
TauEigensystem4::rankdef = "Generic combination failed to resolve the `1`-dim eigenspace for (λ1,λ2)=(`2`,`3`).";
TauEigensystem4::notgenuine = "Eigenvector #`1` for (λ1,λ2)…


---

## 2026-07-19 16:10:28 — 🔄 KERNEL RESTART

## 2026-07-19 16:10:59 — ▶️ RUN CELLS 1–34
Cell 1: ✓ — (no output)
Cell 2: ✓ — (no output)
Cell 3: ✓ — (no output)
Cell 4: ✓ — (no output)
Cell 5: ✓ — Out[16]= \left\{0\right\}
Cell 6: ✓ — (no output)
Cell 7: ✓ — Out[28]= \left\{\frac{1}{3},\frac{1}{7},1,0.3604398751472265372666776505955500912940811467928913386732181476418264+0.9327824485933747164848359820708679294577272928322431080055758094332927\,\mathrm{i},0.5088673809904998281633213722670686199090972690844597081081963623607767+0.8608449271290790931886443936161642458151427736784766273900288346856365\,\mathrm{i},0.5332578805895993725958323450931242953616… [1246 chars total]
Cell 8: ✓ — (no output)
Cell 9: ✓ — Out[32]= \left\{0.\times 10^{-69},0.0063408739330632579023526945309992501993729949274710916715470375206530\right\}
Cell 10: ✓ — (no output)
Cell 11: ✓ — Out[36]= \left\{0\right\}
Cell 12: ✓ — (no output)
Cell 13: ✓ — (no output)
Cell 14: ✓ — (no output)
Cell 15: ✓ — Out[44]= \left\{0\right\}
Cell 16: ✓ — (no output)
Cell 17: ✓ — (no output)
Cell 18: ✓ — Out[51]= \left\{\text{C2=$\chi$1 Id OK},\text{D4=$\chi$2 Id OK},\text{G6=$\chi$3 Id OK},\text{maxCommutator}\to 0.\times 10^{-63}\right\}
Cell 19: ✓ — (no output)
Cell 20: ✓ — Out[55]= \left\{\text{all reps sum-rule OK},\text{adjoint$\otimes$fund @(H1,H2,H3)=(1,0,0) mult}\to 6\right\}
Cell 21: ✓ — (no output)
Cell 22: ✓ — (no output)
Cell 23: ✓ — Out[68]= \left\{\text{nstates}\to 60,\text{worstResidual@u=17/5}\to 3.4603883675228365407978639751541030187243330722550090075915650\times 10^{-64}\right\}
Cell 24: ✓ — (no output)
Cell 25: ✓ — (no output)
Cell 26: ✓ — (no output)
Cell 27: ✓ — (no output)
Cell 28: ✓ — (no output)
Cell 29: ✓ — Out[80]= \left\{\text{closed-form M\_k=$\Lambda$\_1-n\_k matches Mdeg over full sweep, worst}\to 0\right\}
Cell 30: ✓ — (no output)
Cell 31: ✓ — (no output)
Cell 32: ✓ — (no output)
Cell 33: ✓ — Out[89]= \left\{\text{state}\to \left\{0,1,1,0\right\},\text{perQresidual}\to \left\{0,0,0,0\right\}\right\}
Cell 34: ✓ — Out[90]= \begin{matrix}\left\{\begin{pmatrix}1 & 0 & 0 & 0\\1 & 0 & 0 & 0\end{pmatrix},\text{nstates}\to 16,\text{worstTQ}\to 0,\text{worstCasoratian}\to 0\right\}\\\left\{\begin{pmatrix}1 & 0 & 0 & 0\\1 & 1 & 1 & 0\end{pmatrix},\text{nstates}\to 16,\text{worstTQ}\to 0,\text{worstCasoratian}\to 0\right\}\\\left\{\begin{pmatrix}1 & 1 & 0 & 0\\1 & 0 & 0 & 0\end{pmatrix},\text{nstates}\to 24,\text{wo… [743 chars total]

