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

## 2026-07-19 17:22:18 — ▶️ RUN CELLS 1–34
Cell 1: ✓ — (no output)
Cell 2: ✓ — (no output)
Cell 3: ✓ — (no output)
Cell 4: ✓ — (no output)
Cell 5: ✓ — Out[18]= \left\{0\right\}
Cell 6: ✓ — (no output)
Cell 7: ✓ — Out[30]= \left\{\frac{1}{3},\frac{1}{7},1,0.3604398751472265372666776505955500912940811467928913386732181476418264+0.9327824485933747164848359820708679294577272928322431080055758094332927\,\mathrm{i},0.5088673809904998281633213722670686199090972690844597081081963623607767+0.8608449271290790931886443936161642458151427736784766273900288346856365\,\mathrm{i},0.5332578805895993725958323450931242953616… [1246 chars total]
Cell 8: ✓ — (no output)
Cell 9: ✓ — Out[34]= \left\{0.\times 10^{-69},0.0063408739330632579023526945309992501993729949274710916715470375206530\right\}
Cell 10: ✓ — (no output)
Cell 11: ✓ — Out[38]= \left\{0\right\}
Cell 12: ✓ — (no output)
Cell 13: ✓ — (no output)
Cell 14: ✓ — (no output)
Cell 15: ✓ — Out[46]= \left\{0\right\}
Cell 16: ✓ — (no output)
Cell 17: ✓ — (no output)
Cell 18: ✓ — Out[53]= \left\{\text{C2=$\chi$1 Id OK},\text{D4=$\chi$2 Id OK},\text{G6=$\chi$3 Id OK},\text{maxCommutator}\to 0.\times 10^{-63}\right\}
Cell 19: ✓ — (no output)
Cell 20: ✓ — Out[57]= \left\{\text{all reps sum-rule OK},\text{adjoint$\otimes$fund @(H1,H2,H3)=(1,0,0) mult}\to 6\right\}
Cell 21: ✓ — (no output)
Cell 22: ✓ — (no output)
Cell 23: ✓ — Out[70]= \left\{\text{nstates}\to 60,\text{worstResidual@u=17/5}\to 3.4603883675228365407978639751541030187243330722550090075915650\times 10^{-64}\right\}
Cell 24: ✓ — (no output)
Cell 25: ✓ — (no output)
Cell 26: ✓ — (no output)
Cell 27: ✓ — (no output)
Cell 28: ✓ — (no output)
Cell 29: ✓ — Out[82]= \left\{\text{closed-form M\_k=$\Lambda$\_1-n\_k matches Mdeg over full sweep, worst}\to 0\right\}
Cell 30: ✓ — (no output)
Cell 31: ✓ — (no output)
Cell 32: ✓ — (no output)
Cell 33: ✓ — Out[91]= \left\{\text{state}\to \left\{0,1,1,0\right\},\text{perQresidual}\to \left\{0,0,0,0\right\}\right\}
Cell 34: ✓ — Out[92]= \begin{matrix}\left\{\begin{pmatrix}1 & 0 & 0 & 0\\1 & 0 & 0 & 0\end{pmatrix},\text{nstates}\to 16,\text{worstTQ}\to 0,\text{worstCasoratian}\to 0\right\}\\\left\{\begin{pmatrix}1 & 0 & 0 & 0\\1 & 1 & 1 & 0\end{pmatrix},\text{nstates}\to 16,\text{worstTQ}\to 0,\text{worstCasoratian}\to 0\right\}\\\left\{\begin{pmatrix}1 & 1 & 0 & 0\\1 & 0 & 0 & 0\end{pmatrix},\text{nstates}\to 24,\text{wo… [743 chars total]

## 2026-07-19 17:22:26 — 📥 BULK INSERT 1 CELL(S) at positions 35–35
1. [code] (*Section boundary: switch to the untwisted (z=1) regime for the Casimir-eigenvalue exploration↵  be


## 2026-07-19 17:22:26
**Input:**
```wolfram
(*Section boundary: switch to the untwisted (z=1) regime for the Casimir-eigenvalue exploration
  below. qmin/qming/tg/tgCoeffs are memoized on (λ1,λ2) but read z live inside their definitions —
  without this ClearAll they would silently return stale twisted-z values. T/Tg/ν/ν1/qdetT/H1op/
  H2op/H3op are NOT memoized and need no clearing; they recompute correctly automatically.*)
ClearAll[qmin, qming, tg, tgCoeffs];
z[1] = 1; z[2] = 1; z[3] = 1; z[4] = 1;
{χ1, χ2, χ3, χ4}
```
**Output:**
```
Out[95]= \left\{4,6,4,1\right\}
```
## 2026-07-19 17:22:33 — 📥 BULK INSERT 1 CELL(S) at positions 36–36
1. [code] (*Quadratic Casimir of gl(4) on the two-site space: C2 = Sum_ij EE[i,j].EE[j,i]. Exact rational↵  (E


## 2026-07-19 17:22:33
**Input:**
```wolfram
(*Quadratic Casimir of gl(4) on the two-site space: C2 = Sum_ij EE[i,j].EE[j,i]. Exact rational
  (EE is built from JJ, which is rational-valued in this notebook's GT-pattern convention).*)
ClearAll[C2op];
C2op[λ1_, λ2_] := C2op[λ1, λ2] = Sum[EE[λ1, λ2][i, j] . EE[λ1, λ2][j, i], {i, 4}, {j, 4}];
```
**Output:**
```
(no output)
```
## 2026-07-19 17:22:40 — 📥 BULK INSERT 1 CELL(S) at positions 37–37
1. [code] (*sanity check: at z=1 the untwisted transfer-matrix coefficients and Cartan operators must commute↵


## 2026-07-19 17:22:41
**Input:**
```wolfram
(*sanity check: at z=1 the untwisted transfer-matrix coefficients and Cartan operators must commute
  with the quadratic Casimir C2op, confirming full su(4) symmetry restoration (this would FAIL for
  z!=1). Tested on {1,1,0,0}⊗{1,1,0,0} (6⊗6=36, the smallest nontrivial rep pair). tgCoeffs
  recomputes fresh here (cleared in Step 1), now reading z=1.*)
Module[{λ1={1,1,0,0},λ2={1,1,0,0},c1,c2,c3,comm},
  {c1,c2,c3}=tgCoeffs[λ1,λ2];
  comm=Max[Table[Max[Abs[Flatten[m.C2op[λ1,λ2]-C2op[λ1,λ2].m]]],
     {m,Join[c1,c2,c3,{H1op[λ1,λ2],H2op[λ1,λ2],H3op[λ1,λ2]}]}]];
  If[comm>10^-8,Print["C2 COMMUTATION FAIL: ",comm];Abort[]];
  {"C2 commutes with tau-coeffs and Cartan at z=1, maxCommutator"->comm}
]
```
**Output:**
```
2 kernel messages: Set::shape: Lists {c1$26478, c2$26478, c3$26478} and tgCoeffs[{1, 1, 0, 0}, {1, 1, 0, 0}] are not the same shape.
Table::iterb: Iterator {m, Join[c1$26478, c2$26478, c3$26478, {{{0, 0, 0, 0, 0, 0, 0, … [also: Set::shape] result: Out[98]= \left\{\text{C2 commutes with tau-coeffs and Cartan at z=1, maxCommutator}\to \mathrm{Table}\left[\max\left(\left|\mathrm{Flatten}\left[m.\,C2op\left(\lambda1$26478,\lambda2$26478\right)-C2op…
```
## 2026-07-19 17:23:16 — 📥 BULK INSERT 2 CELL(S) at positions 36–37
1. [code] (*Restore the qmin/qming/tg/tgCoeffs definitions wiped by the ClearAll above -- ClearAll removes a↵ 
2. [code] tgCoeffs[λ1_, λ2_] := tgCoeffs[λ1, λ2] = Module[{uu, m1, m2, m3, cl1, cl2, cl3},↵   m1 = tg[λ1, λ2][


## 2026-07-19 17:23:16
**Input:**
```wolfram
(*Restore the qmin/qming/tg/tgCoeffs definitions wiped by the ClearAll above -- ClearAll removes a
  function's rules entirely, not just its memoized cache, so these MUST be redefined (not merely
  re-called) here. Identical source to the twisted-case build; now reads the freshly-assigned z=1.*)
qmin[λ1_, λ2_][II_List, JJ_List][u_] := qmin[λ1, λ2][II, JJ][u] = Sum[Signature[σ] Dot @@ Table[T[λ1, λ2][II[[σ[[k]]]], JJ[[k]]][u - (k - 1) h], {k, Length[II]}], {σ, Permutations[Range[Length[II]]]}];
qming[λ1_, λ2_][II_List, JJ_List][u_] := qming[λ1, λ2][II, JJ][u] = Sum[Signature[σ] Dot @@ Table[Tg[λ1, λ2][II[[σ[[k]]]], JJ[[k]]][u - (k - 1) h], {k, Length[II]}], {σ, Permutations[Range[Length[II]]]}];
tg[λ1_,λ2_][a_,1][u_]:= Total[qming[λ1, λ2][#, #][u] & /@ Subsets[{1, 2, 3, 4}, {a}]]
```
**Output:**
```
(no output)
```

## 2026-07-19 17:23:17
**Input:**
```wolfram
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
## 2026-07-19 17:23:26 — ▶️ RUN CELL 39 (0.03 s)
**In [39]:** `(*sanity check: at z=1 the untwisted transfer-matrix coefficients and Cartan operators must commute↵  with the quadratic Casimir C2op, confirming full su(4) symmetry restoration (this would FAIL for↵ `
**Out:** Out[103]= \left\{\text{C2 commutes with tau-coeffs and Cartan at z=1, maxCommutator}\to 0\right\}

## 2026-07-19 17:24:31 — 📥 BULK INSERT 1 CELL(S) at positions 40–40
1. [code] ClearAll[TauEigensystemCasimir];↵TauEigensystemCasimir::rankdef = "Generic combination failed to res


## 2026-07-19 17:24:32
**Input:**
```wolfram
ClearAll[TauEigensystemCasimir];
TauEigensystemCasimir::rankdef = "Generic combination failed to resolve the `1`-dim eigenspace for (λ1,λ2)=(`2`,`3`).";
TauEigensystemCasimir::notgenuine = "Eigenvector #`1` for (λ1,λ2)=(`2`,`3`) failed the common-eigenvector residual check (residual=`4`).";
TauEigensystemCasimir::noninteger = "Weight `1` not near-integer for eigenvector #`2` at (λ1,λ2)=(`3`,`4`).";
TauEigensystemCasimir::degnull = "Total resolved states `1` != predicted dim `2` for (λ1,λ2)=(`3`,`4`).";
TauEigensystemCasimir[λ1_, λ2_] := TauEigensystemCasimir[λ1, λ2] = Module[
   {d, c1, c2, c3, mats, nontriv, rSets, rset, Kmix, vals, vecs, rayleigh, recs, byHW, table},
   d = dim[λ1] dim[λ2];
   {c1, c2, c3} = tgCoeffs[λ1, λ2];
   mats = Join[c1, c2, c3, {H1op[λ1, λ2], H2op[λ1, λ2], H3op[λ1, λ2], C2op[λ1, λ2]}];
   (* drop the scalar leading coeffs c1[[3]],c2[[5]],c3[[7]] (=chi_a Id, chi={4,6,4} at z=1) since
      they carry no state info; keep C2op since it's the whole point of this routine *)
   nontriv = Join[c1[[1 ;; 2]], c2[[1 ;; 4]], c3[[1 ;; 6]], {H1op[λ1, λ2], H2op[λ1, λ2], H3op[λ1, λ2], C2op[λ1, λ2]}];
   rSets = Table[N[Table[Prime[100 i + 7 j]^(1/3), {j, Length[nontriv]}], 30], {i, 1, 8}];
   rset = SelectFirst[rSets,
     Module[{ee, vv}, {ee, vv} = Eigensystem[N[#.nontriv, 30]];
        MatrixRank[vv] == d && Length[Union[Chop[ee, 10^-10]]] == d] &,
     Missing["NoGenericR"]];
   If[MissingQ[rset], Message[TauEigensystemCasimir::rankdef, d, λ1, λ2]; Abort[]];
   Kmix = rset . nontriv;
   {vals, vecs} = Eigensystem[N[Kmix, 30]];
   rayleigh[m_, v_] := (v . (m . v))/(v . v);
   recs = Table[
     Module[{v = vecs[[i]], res, H1v, H2v, H3v, C2v},
       res = Max[Abs[#.v - rayleigh[#, v] v]] & /@ mats // Max;
       If[res > 10^-8, Message[TauEigensystemCasimir::notgenuine, i, λ1, λ2, res]; Abort[]];
       H1v = rayleigh[H1op[λ1, λ2], v]; H2v = rayleigh[H2op[λ1, λ2], v]; H3v = rayleigh[H3op[λ1, λ2], v];
       C2v = rayleigh[C2op[λ1, λ2], v];
       If[Abs[H1v - Round[Re[H1v]]] > 10^-6, Message[TauEigensystemCasimir::noninteger, H1v, i, λ1, λ2]; Abort[]];
       If[Abs[H2v - Round[Re[H2v]]] > 10^-6, Message[TauEigensystemCasimir::noninteger, H2v, i, λ1, λ2]; Abort[]];
       If[Abs[H3v - Round[Re[H3v]]] > 10^-6, Message[TauEigensystemCasimir::noninteger, H3v, i, λ1, λ2]; Abort[]];
       <|"H1" -> Round[Re[H1v]], "H2" -> Round[Re[H2v]], "H3" -> Round[Re[H3v]],
         "C2round" -> Round[Re[C2v], 10^-6], "C2" -> Re[C2v],
         "t1c" -> rayleigh[#, v] & /@ c1, "t2c" -> rayleigh[#, v] & /@ c2, "t3c" -> rayleigh[#, v] & /@ c3,
         "vec" -> v|>],
     {i, 1, d}];
   byHW = GroupBy[recs, {#["H1"], #["H2"], #["H3"], #["C2round"]} &];
   table = Association @ Flatten @ KeyValueMap[
      Function[{key, group},
        MapIndexed[Function[{rec, idx}, Append[key, First[idx] - 1] -> rec], group]],
      byHW];
   If[Length[table] != d, Message[TauEigensystemCasimir::degnull, Length[table], d, λ1, λ2]; Abort[]];
   table
];
```
**Output:**
```
(no output)
```
## 2026-07-19 17:24:39 — 📥 BULK INSERT 1 CELL(S) at positions 41–41
1. [code] (*Lambda accessor: a=1,2,3 via Rayleigh-quotient coefficients from TauEigensystemCasimir;↵  a=4 is q


## 2026-07-19 17:24:39
**Input:**
```wolfram
(*Lambda accessor: a=1,2,3 via Rayleigh-quotient coefficients from TauEigensystemCasimir;
  a=4 is qdetT itself (the quantum determinant is central -- same scalar function on EVERY
  state, no per-eigenvector extraction needed).*)
ClearAll[LambdaU];
LambdaU[λ1_, λ2_][H1_, H2_, H3_, C2round_, n_][a_][u_] := Module[
   {rec = TauEigensystemCasimir[λ1, λ2][Key[{H1, H2, H3, C2round, n}]]},
   Which[
     a == 1, rec["t1c"] . u^Range[0, 2],
     a == 2, rec["t2c"] . u^Range[0, 4],
     a == 3, rec["t3c"] . u^Range[0, 6],
     a == 4, qdetT[λ1, λ2][u]
   ]
];
```
**Output:**
```
(no output)
```
## 2026-07-19 17:24:47 — 📥 BULK INSERT 1 CELL(S) at positions 42–42
1. [code] Module[{λ1 = {1, 1, 0, 0}, λ2 = {1, 1, 0, 0}, es, worst},↵  es = TauEigensystemCasimir[λ1, λ2];↵  wo


## 2026-07-19 17:24:48
**Input:**
```wolfram
Module[{λ1 = {1, 1, 0, 0}, λ2 = {1, 1, 0, 0}, es, worst},
  es = TauEigensystemCasimir[λ1, λ2];
  worst = Max@KeyValueMap[
     Function[{k, rec},
       Module[{H1 = k[[1]], H2 = k[[2]], H3 = k[[3]], C2r = k[[4]], n = k[[5]], v = rec["vec"], uu = 17/5},
         Max[
          Max[Abs[(tg[λ1, λ2][1, 1][uu] . v) - LambdaU[λ1, λ2][H1, H2, H3, C2r, n][1][uu] v]],
          Max[Abs[(tg[λ1, λ2][2, 1][uu] . v) - LambdaU[λ1, λ2][H1, H2, H3, C2r, n][2][uu] v]],
          Max[Abs[(tg[λ1, λ2][3, 1][uu] . v) - LambdaU[λ1, λ2][H1, H2, H3, C2r, n][3][uu] v]]]]],
     es];
  If[worst > 10^-8, Print["EIGEN RESIDUAL FAIL: ", worst]; Abort[]];
  {"nstates" -> Length[es], "worstResidual@u=17/5" -> worst}
]
```
**Output:**
```
Out[112]= \left\{\text{nstates}\to 36,\text{worstResidual@u=17/5}\to \max\left(1.553684503556985884685522663\times 10^{-37},\left|29.30666666666666666666666667-1.\,\,\left(---\text{(key absent)}\right)\left(\text{t1c}\right).\,\left\{1,\frac{17}{5},\frac{289}{25}\right\}\right|,\left|106.4263256235827664399092971-1.\,\,\left(---\text{(key absent)}\right)\left(\text{t2c}\right).\,\left\{1,\frac{17}… [1249 chars total]
```

## 2026-07-19 17:25:13
**Input:**
```wolfram
Module[{λ1={1,1,0,0},λ2={1,1,0,0},es,k1},
  es=TauEigensystemCasimir[λ1,λ2];
  k1=First[Keys[es]];
  {"k1"->k1,"directLookup"->es[Key[k1]]===Missing["KeyAbsent",Key[k1]],
   "sameKeyReconstructed"->SameQ[k1,{k1[[1]],k1[[2]],k1[[3]],k1[[4]],k1[[5]]}],
   "Precision[k1[[4]]]"->Precision[k1[[4]]]}
]
```
**Output:**
```
Out= {"k1" -> {0, 2, 0, 16, 0}, "directLookup" -> True, "sameKeyReconstructed" -> True, "Precision[k1[[4]]]" -> Infinity}
```

## 2026-07-19 17:25:30
**Input:**
```wolfram
Module[{λ1={1,1,0,0},λ2={1,1,0,0},es,k1},
  es=TauEigensystemCasimir[λ1,λ2];
  k1=First[Keys[es]];
  {"noKeyWrapper"->es[k1]["H1"],"withKeyWrapperPart"->es[[Key[k1]]]["H1"],"LookupForm"->Lookup[es,Key[k1]]["H1"]}
]
```
**Output:**
```
Out= {"noKeyWrapper" -> 0, "withKeyWrapperPart" -> 0, "LookupForm" -> 0}
```
## 2026-07-19 17:25:40 — ✏️ EDIT CELL 41
(*Lambda accessor: a=1,2,3 via Rayleigh-quotient coefficients from TauEigensystemCasimir;
  a=4 is qdetT itself (the quantum determinant is central -- same scalar function on EVERY
  state, no per-eig…


## 2026-07-19 17:25:40
**Input:**
```wolfram
(*Lambda accessor: a=1,2,3 via Rayleigh-quotient coefficients from TauEigensystemCasimir;
  a=4 is qdetT itself (the quantum determinant is central -- same scalar function on EVERY
  state, no per-eigenvector extraction needed).*)
ClearAll[LambdaU];
LambdaU[λ1_, λ2_][H1_, H2_, H3_, C2round_, n_][a_][u_] := Module[
   {rec = TauEigensystemCasimir[λ1, λ2][{H1, H2, H3, C2round, n}]},
   Which[
     a == 1, rec["t1c"] . u^Range[0, 2],
     a == 2, rec["t2c"] . u^Range[0, 4],
     a == 3, rec["t3c"] . u^Range[0, 6],
     a == 4, qdetT[λ1, λ2][u]
   ]
];
```
**Output:**
```
Out= Null
```
## 2026-07-19 17:26:11 — ▶️ RUN CELL 42 (0.03 s)
**In [42]:** `Module[{λ1 = {1, 1, 0, 0}, λ2 = {1, 1, 0, 0}, es, worst},↵  es = TauEigensystemCasimir[λ1, λ2];↵  worst = Max@KeyValueMap[↵     Function[{k, rec},↵       Module[{H1 = k[[1]], H2 = k[[2]], H3 = k[[3]],`
**Out:** Out[113]= \left\{\text{nstates}\to 36,\text{worstResidual@u=17/5}\to \max\left(1.553684503556985884685522663\times 10^{-37},\left|10.98686803123631175691267507-0.4714045207910316829338962414\,4.00000000000000000000000000.\,\left\{1,\frac{17}{5},\frac{289}{25}\right\}\right|,\left|5.49343401561815587

## 2026-07-19 17:26:53 — ✏️ EDIT CELL 40
ClearAll[TauEigensystemCasimir];
TauEigensystemCasimir::rankdef = "Generic combination failed to resolve the `1`-dim eigenspace for (λ1,λ2)=(`2`,`3`).";
TauEigensystemCasimir::notgenuine = "Eigenvecto…


## 2026-07-19 17:26:53
**Input:**
```wolfram
ClearAll[TauEigensystemCasimir];
TauEigensystemCasimir::rankdef = "Generic combination failed to resolve the `1`-dim eigenspace for (λ1,λ2)=(`2`,`3`).";
TauEigensystemCasimir::notgenuine = "Eigenvector #`1` for (λ1,λ2)=(`2`,`3`) failed the common-eigenvector residual check (residual=`4`).";
TauEigensystemCasimir::noninteger = "Weight `1` not near-integer for eigenvector #`2` at (λ1,λ2)=(`3`,`4`).";
TauEigensystemCasimir::degnull = "Total resolved states `1` != predicted dim `2` for (λ1,λ2)=(`3`,`4`).";
TauEigensystemCasimir[λ1_, λ2_] := TauEigensystemCasimir[λ1, λ2] = Module[
   {d, c1, c2, c3, mats, nontriv, rSets, rset, Kmix, vals, vecs, rayleigh, recs, byHW, table},
   d = dim[λ1] dim[λ2];
   {c1, c2, c3} = tgCoeffs[λ1, λ2];
   mats = Join[c1, c2, c3, {H1op[λ1, λ2], H2op[λ1, λ2], H3op[λ1, λ2], C2op[λ1, λ2]}];
   nontriv = Join[c1[[1 ;; 2]], c2[[1 ;; 4]], c3[[1 ;; 6]], {H1op[λ1, λ2], H2op[λ1, λ2], H3op[λ1, λ2], C2op[λ1, λ2]}];
   rSets = Table[N[Table[Prime[100 i + 7 j]^(1/3), {j, Length[nontriv]}], 30], {i, 1, 8}];
   rset = SelectFirst[rSets,
     Module[{ee, vv}, {ee, vv} = Eigensystem[N[#.nontriv, 30]];
        MatrixRank[vv] == d && Length[Union[Chop[ee, 10^-10]]] == d] &,
     Missing["NoGenericR"]];
   If[MissingQ[rset], Message[TauEigensystemCasimir::rankdef, d, λ1, λ2]; Abort[]];
   Kmix = rset . nontriv;
   {vals, vecs} = Eigensystem[N[Kmix, 30]];
   rayleigh[m_, v_] := (v . (m . v))/(v . v);
   recs = Table[
     Module[{v = vecs[[i]], res, H1v, H2v, H3v, C2v, t1c, t2c, t3c},
       res = Max[Abs[#.v - rayleigh[#, v] v]] & /@ mats // Max;
       If[res > 10^-8, Message[TauEigensystemCasimir::notgenuine, i, λ1, λ2, res]; Abort[]];
       t1c = rayleigh[#, v] & /@ c1;
       t2c = rayleigh[#, v] & /@ c2;
       t3c = rayleigh[#, v] & /@ c3;
       H1v = rayleigh[H1op[λ1, λ2], v]; H2v = rayleigh[H2op[λ1, λ2], v]; H3v = rayleigh[H3op[λ1, λ2], v];
       C2v = rayleigh[C2op[λ1, λ2], v];
       If[Abs[H1v - Round[Re[H1v]]] > 10^-6, Message[TauEigensystemCasimir::noninteger, H1v, i, λ1, λ2]; Abort[]];
       If[Abs[H2v - Round[Re[H2v]]] > 10^-6, Message[TauEigensystemCasimir::noninteger, H2v, i, λ1, λ2]; Abort[]];
       If[Abs[H3v - Round[Re[H3v]]] > 10^-6, Message[TauEigensystemCasimir::noninteger, H3v, i, λ1, λ2]; Abort[]];
       <|"H1" -> Round[Re[H1v]], "H2" -> Round[Re[H2v]], "H3" -> Round[Re[H3v]],
         "C2round" -> Round[Re[C2v], 10^-6], "C2" -> Re[C2v],
         "t1c" -> t1c, "t2c" -> t2c, "t3c" -> t3c, "vec" -> v|>],
     {i, 1, d}];
   byHW = GroupBy[recs, {#["H1"], #["H2"], #["H3"], #["C2round"]} &];
   table = Association @ Flatten @ KeyValueMap[
      Function[{key, group},
        MapIndexed[Function[{rec, idx}, Append[key, First[idx] - 1] -> rec], group]],
      byHW];
   If[Length[table] != d, Message[TauEigensystemCasimir::degnull, Length[table], d, λ1, λ2]; Abort[]];
   table
];
```
**Output:**
```
Out= Null
```
## 2026-07-19 17:27:00 — ▶️ RUN CELL 42 (0.01 s)
**In [42]:** `Module[{λ1 = {1, 1, 0, 0}, λ2 = {1, 1, 0, 0}, es, worst},↵  es = TauEigensystemCasimir[λ1, λ2];↵  worst = Max@KeyValueMap[↵     Function[{k, rec},↵       Module[{H1 = k[[1]], H2 = k[[2]], H3 = k[[3]],`
**Out:** Out[114]= \left\{\text{nstates}\to 36,\text{worstResidual@u=17/5}\to 3.48708354729401192820736\times 10^{-35}\right\}

## 2026-07-19 17:27:57 — 📥 BULK INSERT 1 CELL(S) at positions 43–43
1. [code] (*Test the top-state (mu=lambda1+lambda2) closed-form conjecture: Lambda_1(u) = Sum_i nu[i](u),↵  wi


## 2026-07-19 17:27:58
**Input:**
```wolfram
(*Test the top-state (mu=lambda1+lambda2) closed-form conjecture: Lambda_1(u) = Sum_i nu[i](u),
  with NO shifts (unlike qdetT's shifted product) since a=1 principal minors are literally T[i,i](u)
  with zero row-shift in the quantum-minor definition. mu=lambda1+lambda2 has weight (S+P,S+P,0,0),
  i.e. (H1,H2,H3)=(0,S+P,0), and occurs with multiplicity exactly 1 (general fact) -- so there is
  exactly one key in TauEigensystemCasimir's output with that (H1,H2,H3).*)
Module[{λ1 = {1, 1, 0, 0}, λ2 = {1, 1, 0, 0}, S = 1, P = 1, es, topKeys, key, H1, H2, H3, C2r, n, conj, actual, worst},
  es = TauEigensystemCasimir[λ1, λ2];
  topKeys = Select[Keys[es], {#[[1]], #[[2]], #[[3]]} == {0, S + P, 0} &];
  If[Length[topKeys] != 1, Print["TOP STATE NOT UNIQUE: ", topKeys]; Abort[]];
  key = First[topKeys];
  {H1, H2, H3, C2r, n} = key;
  conj[u_] := ν[λ1, λ2][1][u] + ν[λ1, λ2][2][u] + ν[λ1, λ2][3][u] + ν[λ1, λ2][4][u];
  worst = Max@Table[Abs[N[LambdaU[λ1, λ2][H1, H2, H3, C2r, n][1][u0] - conj[u0], 30]], {u0, {2/3, 9/5, 4, 23/7}}];
  {"topKey" -> key, "conjectureMatches" -> (worst < 10^-8), "worst" -> worst}
]
```
**Output:**
```
Out[115]= \left\{\text{topKey}\to \left\{0,2,0,16,0\right\},\text{conjectureMatches}\to \mathrm{True},\text{worst}\to 0.\times 10^{-25}\right\}
```

## 2026-07-19 17:28:27
**Input:**
```wolfram
{dim[{1,1,0,0}],dim[{2,2,0,0}],{1,1}->dim[{1,1,0,0}]dim[{1,1,0,0}],{2,1}->dim[{2,2,0,0}]dim[{1,1,0,0}],{1,2}->dim[{1,1,0,0}]dim[{2,2,0,0}],{2,2}->dim[{2,2,0,0}]dim[{2,2,0,0}]}
```
**Output:**
```
Out= {6, 20, {1, 1} -> 36, {2, 1} -> 120, {1, 2} -> 120, {2, 2} -> 400}
```
## 2026-07-19 17:28:37 — 📥 BULK INSERT 1 CELL(S) at positions 44–44
1. [code] (*For each (S,P) in the sweep, run TauEigensystemCasimir and check every (H1,H2,H3,C2round) sector↵


## 2026-07-19 17:28:37
**Input:**
```wolfram
(*For each (S,P) in the sweep, run TauEigensystemCasimir and check every (H1,H2,H3,C2round) sector
  with >1 resolved state: if all states in that sector give the SAME Lambda_a(u) (a=1,2,3), it's a
  harmless numerical near-degeneracy in C2round's rounding, not genuine multiplicity N_mu>1. If they
  give GENUINELY DIFFERENT Lambda_a(u), that sector has real multiplicity -- no pure-Casimir formula
  can apply there (reported, not treated as failure).*)
multData = Table[
  Module[{λ1 = {S, S, 0, 0}, λ2 = {P, P, 0, 0}, es, byKey3, multSectors, genuineMult},
    es = TauEigensystemCasimir[λ1, λ2];
    byKey3 = GroupBy[Keys[es], {#[[1]], #[[2]], #[[3]], #[[4]]} &];
    multSectors = Select[byKey3, Length[#] > 1 &];
    genuineMult = Select[multSectors, Function[keys,
       Module[{lam1s = LambdaU[λ1, λ2][keys[[1, 1]], keys[[1, 2]], keys[[1, 3]], keys[[1, 4]], keys[[1, 5]]][1][17/5]},
         Max[Abs[LambdaU[λ1, λ2][#[[1]], #[[2]], #[[3]], #[[4]], #[[5]]][1][17/5] - lam1s] & /@ keys]] > 10^-6]];
    {S, P, "dim" -> dim[λ1] dim[λ2], "nstates" -> Length[es], "sectorsWithMult" -> Length[multSectors],
     "genuineMultSectors" -> Length[genuineMult], "multFree" -> (Length[genuineMult] == 0)}],
  {SP, {{1, 1}, {2, 1}, {1, 2}, {2, 2}}}, {S = SP[[1]], P = SP[[2]]}];
Column[multData]
```
**Output:**
```
7 kernel messages: Table::write: Tag Set in S = SP[[1]] is Protected.
Table::write: Tag Set in S = SP[[1]] is Protected. [also: Table::write, Table::write, Table::write, Table::write, General::stop, Table::write] result: Out[117]= \mathrm{Column}\left[\mathrm{Table}\left[\mathrm{Module}\left[\left\{\lambda1=\left\{S,S,0,0\right\},\lambda2=\left\{P,P,0,0\right\},\mathrm{es},byKey3,\mathrm{multSectors},\mathrm{genuineMu…
```
## 2026-07-19 17:28:51 — ✏️ EDIT CELL 44
(*For each (S,P) in the sweep, run TauEigensystemCasimir and check every (H1,H2,H3,C2round) sector
  with >1 resolved state: if all states in that sector give the SAME Lambda_a(u) (a=1,2,3), it's a
  …


## 2026-07-19 17:34:16
**Input:**
```wolfram
(*For each (S,P) in the sweep, run TauEigensystemCasimir and check every (H1,H2,H3,C2round) sector
  with >1 resolved state: if all states in that sector give the SAME Lambda_a(u) (a=1,2,3), it's a
  harmless numerical near-degeneracy in C2round's rounding, not genuine multiplicity N_mu>1. If they
  give GENUINELY DIFFERENT Lambda_a(u), that sector has real multiplicity -- no pure-Casimir formula
  can apply there (reported, not treated as failure).*)
multData = Table[
  Module[{S = SP[[1]], P = SP[[2]], λ1, λ2, es, byKey3, multSectors, genuineMult},
    λ1 = {S, S, 0, 0}; λ2 = {P, P, 0, 0};
    es = TauEigensystemCasimir[λ1, λ2];
    byKey3 = GroupBy[Keys[es], {#[[1]], #[[2]], #[[3]], #[[4]]} &];
    multSectors = Select[byKey3, Length[#] > 1 &];
    genuineMult = Select[multSectors, Function[keys,
       Module[{lam1s = LambdaU[λ1, λ2][keys[[1, 1]], keys[[1, 2]], keys[[1, 3]], keys[[1, 4]], keys[[1, 5]]][1][17/5]},
         Max[Abs[LambdaU[λ1, λ2][#[[1]], #[[2]], #[[3]], #[[4]], #[[5]]][1][17/5] - lam1s] & /@ keys]] > 10^-6]];
    {S, P, "dim" -> dim[λ1] dim[λ2], "nstates" -> Length[es], "sectorsWithMult" -> Length[multSectors],
     "genuineMultSectors" -> Length[genuineMult], "multFree" -> (Length[genuineMult] == 0)}],
  {SP, {{1, 1}, {2, 1}, {1, 2}, {2, 2}}}];
Column[multData]
```
**Output:**
```
Out= Column[{{1, 1, "dim" -> 36, "nstates" -> 36, "sectorsWithMult" -> 2, "genuineMultSectors" -> 0, "multFree" -> True}, {2, 1, "dim" -> 120, "nstates" -> 120, "sectorsWithMult" -> 20, "genuineMultSectors" -> 0, "multFree" -> True}, {1, 2, "dim" -> 120, "nstates" -> 120, "sectorsWithMult" -> 20, "genuineMultSectors" -> 0, "multFree" -> True}, {2, 2, "dim" -> 400, "nstates" -> 400, "sectorsWithMult" -> 77, "genuineMultSectors" -> 0, "multFree" -> True}}]
```

## 2026-07-19 17:34:23
**Input:**
```wolfram
1+1
```
**Output:**
```
Out= 2
```

## 2026-07-19 17:34:29
**Input:**
```wolfram
ValueQ[multData]
```
**Output:**
```
Out= True
```

## 2026-07-19 17:34:32
**Input:**
```wolfram
multData
```
**Output:**
```
Out= {{1, 1, "dim" -> 36, "nstates" -> 36, "sectorsWithMult" -> 2, "genuineMultSectors" -> 0, "multFree" -> True}, {2, 1, "dim" -> 120, "nstates" -> 120, "sectorsWithMult" -> 20, "genuineMultSectors" -> 0, "multFree" -> True}, {1, 2, "dim" -> 120, "nstates" -> 120, "sectorsWithMult" -> 20, "genuineMultSectors" -> 0, "multFree" -> True}, {2, 2, "dim" -> 400, "nstates" -> 400, "sectorsWithMult" -> 77, "genuineMultSectors" -> 0, "multFree" -> True}}
```
## 2026-07-19 17:34:43 — ▶️ RUN CELL 44 (0.01 s)
**In [44]:** `(*For each (S,P) in the sweep, run TauEigensystemCasimir and check every (H1,H2,H3,C2round) sector↵  with >1 resolved state: if all states in that sector give the SAME Lambda_a(u) (a=1,2,3), it's a↵  `
**Out:** Out[119]= \begin{matrix}\left\{1,1,\text{dim}\to 36,\text{nstates}\to 36,\text{sectorsWithMult}\to 2,\text{genuineMultSectors}\to 0,\text{multFree}\to \mathrm{True}\right\}\\\left\{2,1,\text{dim}\to 120,\text{nstates}\to 120,\text{sectorsWithMult}\to 20,\text{genuineMultSectors}\to 0,\text{multFree}

## 2026-07-19 17:35:58 — 📥 BULK INSERT 1 CELL(S) at positions 45–45
1. [code] (*Collect (mu-label via (H1,H2,H3,C2), Lambda_1(u) coeffs) rows for every multiplicity-free sector↵


## 2026-07-19 17:35:59
**Input:**
```wolfram
(*Collect (mu-label via (H1,H2,H3,C2), Lambda_1(u) coeffs) rows for every multiplicity-free sector
  across the sweep (using multData's multFree flag; skip sectors flagged genuinely-multiplicity in
  Task 4). Fit each u^j coefficient of Lambda_1 as an affine function of the C2 eigenvalue, using the
  first two distinct-C2 rows as anchors, then hard-assert against every remaining row.*)
casData = Flatten[Table[
   Module[{S = SP[[1]], P = SP[[2]], λ1, λ2, es, byKey3, genuineMultKeys},
     λ1 = {S, S, 0, 0}; λ2 = {P, P, 0, 0};
     es = TauEigensystemCasimir[λ1, λ2];
     byKey3 = GroupBy[Keys[es], {#[[1]], #[[2]], #[[3]], #[[4]]} &];
     genuineMultKeys = Flatten[Values[Select[byKey3, Length[#] > 1 &]], 1]; (* skip these -- ambiguous *)
     KeyValueMap[
       Function[{key, rec}, If[MemberQ[genuineMultKeys, key], Nothing,
         {S, P, key, rec["C2"], CoefficientList[LambdaU[λ1, λ2][key[[1]], key[[2]], key[[3]], key[[4]], key[[5]]][1][uu], uu]}]],
       es]],
   {SP, {{1, 1}, {2, 1}, {1, 2}}}], 2];
(*fit u^j coefficients as affine in C2 from the first two rows with distinct C2, verify on the rest*)
Module[{distinctRows, anchor1, anchor2, coeffFit, worst},
  distinctRows = DeleteDuplicatesBy[casData, #[[4]] &];
  If[Length[distinctRows] < 2, Print["NOT ENOUGH DISTINCT C2 VALUES TO FIT"]; Abort[]];
  {anchor1, anchor2} = distinctRows[[1 ;; 2]];
  coeffFit[j_] := Module[{c1v = anchor1[[5, j + 1]], c2v = anchor2[[5, j + 1]], x1 = anchor1[[4]], x2 = anchor2[[4]], A, B},
     B = (c1v - c2v)/(x1 - x2); A = c1v - B x1; {A, B}];
  worst = Max@Table[
     Module[{A, B, pred},
       {A, B} = coeffFit[j];
       Max[Abs[(A + B #[[4]]) - #[[5, j + 1]]] & /@ casData]],
     {j, 0, 2}];
  {"nDataRows" -> Length[casData], "nDistinctC2" -> Length[distinctRows],
   "coeffFits(A,B)_j=0,1,2" -> Table[coeffFit[j], {j, 0, 2}], "worstDeviation" -> worst}
]
```
**Output:**
```
2 kernel messages: Part::partd: Part specification 1[[4]] is longer than depth of object.
Part::partd: Part specification 1[[4]] is longer than depth of object.
Part::partd: Part specification 16.00000000000000000000000… [also: Part::partd] result: Out[121]= \left\{\text{nDataRows}\to 835,\text{nDistinctC2}\to 71,\text{coeffFits(A,B)\_j=0,1,2}\to \begin{pmatrix}1\llbracket 5,1\rrbracket -\frac{1\llbracket 4\rrbracket \,\left(1\llbracket 5,1\rrbr…
```
## 2026-07-19 17:36:23 — ✏️ EDIT CELL 45
(*Collect (mu-label via (H1,H2,H3,C2), Lambda_1(u) coeffs) rows for every multiplicity-free sector
  across the sweep (using multData's multFree flag; skip sectors flagged genuinely-multiplicity in
  …


## 2026-07-19 17:36:23
**Input:**
```wolfram
(*Collect (mu-label via (H1,H2,H3,C2), Lambda_1(u) coeffs) rows for every multiplicity-free sector
  across the sweep (using multData's multFree flag; skip sectors flagged genuinely-multiplicity in
  Task 4). Fit each u^j coefficient of Lambda_1 as an affine function of the C2 eigenvalue, using the
  first two distinct-C2 rows as anchors, then hard-assert against every remaining row.*)
casData = Flatten[Table[
   Module[{S = SP[[1]], P = SP[[2]], λ1, λ2, es, byKey3, genuineMultKeys},
     λ1 = {S, S, 0, 0}; λ2 = {P, P, 0, 0};
     es = TauEigensystemCasimir[λ1, λ2];
     byKey3 = GroupBy[Keys[es], {#[[1]], #[[2]], #[[3]], #[[4]]} &];
     genuineMultKeys = Flatten[Values[Select[byKey3, Length[#] > 1 &]], 1]; (* skip these -- ambiguous *)
     KeyValueMap[
       Function[{key, rec}, If[MemberQ[genuineMultKeys, key], Nothing,
         {S, P, key, rec["C2"], CoefficientList[LambdaU[λ1, λ2][key[[1]], key[[2]], key[[3]], key[[4]], key[[5]]][1][uu], uu]}]],
       es]],
   {SP, {{1, 1}, {2, 1}, {1, 2}}}], 1];
(*fit u^j coefficients as affine in C2 from the first two rows with distinct C2, verify on the rest*)
Module[{distinctRows, anchor1, anchor2, coeffFit, worst},
  distinctRows = DeleteDuplicatesBy[casData, #[[4]] &];
  If[Length[distinctRows] < 2, Print["NOT ENOUGH DISTINCT C2 VALUES TO FIT"]; Abort[]];
  {anchor1, anchor2} = distinctRows[[1 ;; 2]];
  coeffFit[j_] := Module[{c1v = anchor1[[5, j + 1]], c2v = anchor2[[5, j + 1]], x1 = anchor1[[4]], x2 = anchor2[[4]], A, B},
     B = (c1v - c2v)/(x1 - x2); A = c1v - B x1; {A, B}];
  worst = Max@Table[
     Module[{A, B, pred},
       {A, B} = coeffFit[j];
       Max[Abs[(A + B #[[4]]) - #[[5, j + 1]]] & /@ casData]],
     {j, 0, 2}];
  {"nDataRows" -> Length[casData], "nDistinctC2" -> Length[distinctRows],
   "coeffFits(A,B)_j=0,1,2" -> Table[coeffFit[j], {j, 0, 2}], "worstDeviation" -> worst}
]
```
**Output:**
```
Out= {"nDataRows" -> 167, "nDistinctC2" -> 32, "coeffFits(A,B)_j=0,1,2" -> {{-4.85714285714285714285714285714285714276`26.261107801242673, 0.49999999999999999999999999999999999999`26.6391849134416}, {-5.90476190476190476190476190476190476191`26.565607544911966, 0``27.136827835305855}, {4.`26.56560754491197, 0``27.305970234406214}}, "worstDeviation" -> 4.71428571428571428571428571428571428581`26.006832748047675}
```

## 2026-07-19 17:36:54
**Input:**
```wolfram
{Length[TauEigensystemCasimir[{1,1,0,0},{1,1,0,0}]],Length[TauEigensystemCasimir[{2,2,0,0},{1,1,0,0}]],Length[TauEigensystemCasimir[{1,1,0,0},{2,2,0,0}]],Length[casData],Length[Union[casData]]}
```
**Output:**
```
Out= {36, 120, 120, 167, 167}
```

## 2026-07-19 17:37:25
**Input:**
```wolfram
Module[{distinctRows, anchor1, anchor2, coeffFit, perJ},
  distinctRows = DeleteDuplicatesBy[casData, #[[4]] &];
  {anchor1, anchor2} = distinctRows[[1 ;; 2]];
  coeffFit[j_] := Module[{c1v = anchor1[[5, j + 1]], c2v = anchor2[[5, j + 1]], x1 = anchor1[[4]], x2 = anchor2[[4]], A, B},
     B = (c1v - c2v)/(x1 - x2); A = c1v - B x1; {A, B}];
  perJ = Table[
     Module[{A, B, devs, worstRow},
       {A, B} = coeffFit[j];
       devs = ((A + B #[[4]]) - #[[5, j + 1]]) & /@ casData;
       worstRow = casData[[Ordering[Abs[devs], -1][[1]]]];
       {"j" -> j, "A" -> A, "B" -> B, "worstDev" -> Max[Abs[devs]], "worstRow(S,P,key,C2)" -> worstRow[[1 ;; 4]]}],
     {j, 0, 2}];
  perJ
]
```
**Output:**
```
Out= {{"j" -> 0, "A" -> -4.85714285714285714285714285714285714276`26.261107801242673, "B" -> 0.49999999999999999999999999999999999999`26.6391849134416, "worstDev" -> 4.71428571428571428571428571428571428581`26.006832748047675, "worstRow(S,P,key,C2)" -> {2, 1, {1, -1, 1, 14, 0}, 14.00000000000000000000000000000000000021`27.416411814279112}}, {"j" -> 1, "A" -> -5.90476190476190476190476190476190476191`26.565607544911966, "B" -> 0``27.136827835305855, "worstDev" -> 2.00000000000000000000000000000000000007`25.848696993853963, "worstRow(S,P,key,C2)" -> {1, 2, {1, 0, -1, 14, 0}, 14.00000000000000000000000000000000000022`27.456318915835727}}, {"j" -> 2, "A" -> 4.`26.56560754491197, "B" -> 0``27.305970234406214, "worstDev" -> 0``25.57087880633305, "worstRow(S,P,key,C2)" -> {2, 1, {0, -3, 0, 30, 0}, 30.`27.539469847664094}}}
```

## 2026-07-19 17:37:51
**Input:**
```wolfram
Table[
  Module[{xy, fit, resid},
    xy = {#[[4]], #[[5, j + 1]]} & /@ casData;
    fit = Fit[N[xy, 30], {1, x}, x];
    resid = Max[Abs[(fit /. x -> #[[1]]) - #[[2]]] & /@ xy];
    {"j" -> j, "fit" -> fit, "maxResidual" -> resid}],
  {j, 0, 2}]
```
**Output:**
```
Out= {{"j" -> 0, "fit" -> -3.79284436312655476177822991874268952683`26.35172772496687 + 0.30352429540663584347487201289619912862`26.351727724966864*x, "maxResidual" -> 3.02792434399491847544140683323599645385`25.8950716671945}, {"j" -> 1, "fit" -> -5.43422993898627634079618986883763784899`26.35172772496687 - 0.08686294308338204814795132061431123109`26.351727724966864*x, "maxResidual" -> 1.25445076260827974703725354733131273224`25.592147989139775}, {"j" -> 2, "fit" -> 4.00000000000000000000000000000000000001`26.35172772496687 - 1.3704477937440532423380956238137746476`26.351727724966864*^-40*x, "maxResidual" -> 0``25.6827209440083}}
```
## 2026-07-19 17:38:09 — 📥 BULK INSERT 1 CELL(S) at positions 46–46
1. [code] (*CONCLUSION of the Casimir-affine conjecture test (two independent checks):↵  (1) two-anchor extrap


## 2026-07-19 17:38:09
**Input:**
```wolfram
(*CONCLUSION of the Casimir-affine conjecture test (two independent checks):
  (1) two-anchor extrapolation (cell above): j=0 worstDev=4.71, j=1 worstDev=2.00, both exact
      (26-digit precision, not roundoff) -- j=1's anchor pair happened to share the same value
      (B=0), a possible anchor-selection artifact, so:
  (2) least-squares best-fit line across ALL 167 multiplicity-free data points (3 rep pairs,
      (1,1)/(2,1)/(1,2)): j=0 maxResidual=3.03, j=1 maxResidual=1.25 -- confirms check (1), not an
      anchor artifact. Only j=2 (leading u^2 coefficient) is exactly constant (=chi1=4, B~0 to
      machine precision) -- already known/trivial from the Task 1 commutation check, not a real
      test of the Casimir hypothesis.
  CONCLUSION: the simple ansatz "each u^j coefficient of Lambda_1(u) is an affine function of
  C2(mu) alone" does NOT hold for j=0,1 across this multiplicity-free sweep. The only verified
  closed form remains the Task 3 top-state (mu=lambda1+lambda2) result, which does not involve
  Casimirs at all (it follows directly from the pseudovacuum nu[i](u) formula). This is a
  precise, reportable negative result for the single-quadratic-Casimir-alone hypothesis -- not
  evidence the eigenvalues are unrelated to Casimirs, only that C2 alone (via a linear relation)
  is insufficient to fix Lambda_1(u)'s subleading structure across different mu in this rep family.*)
{"j=0,1 affine-in-C2-alone: FALSIFIED (both anchor and least-squares tests agree)",
 "j=2 (leading): trivially constant, not a real test",
 "verified positive result: top-state (mu=lambda1+lambda2) formula from Task 3"}
```
**Output:**
```
Out[122]= \left\{\text{j=0,1 affine-in-C2-alone: FALSIFIED (both anchor and least-squares tests agree)},\text{j=2 (leading): trivially constant, not a real test},\text{verified positive result: top-state (mu=lambda1+lambda2) formula from Task 3}\right\}
```
## 2026-07-19 17:38:16 — 📥 BULK INSERT 1 CELL(S) at positions 46–46
1. [code] (*Robustness check on the two-anchor fit above: least-squares best-fit line for each u^j↵  coefficie


## 2026-07-19 17:38:17
**Input:**
```wolfram
(*Robustness check on the two-anchor fit above: least-squares best-fit line for each u^j
  coefficient vs C2, across ALL 167 multiplicity-free data points (not just 2 anchors) -- confirms
  whether the two-anchor result was a genuine falsification or an unlucky anchor pair.*)
Table[
  Module[{xy, fit, resid},
    xy = {#[[4]], #[[5, j + 1]]} & /@ casData;
    fit = Fit[N[xy, 30], {1, x}, x];
    resid = Max[Abs[(fit /. x -> #[[1]]) - #[[2]]] & /@ xy];
    {"j" -> j, "fit" -> fit, "maxResidual" -> resid}],
  {j, 0, 2}]
```
**Output:**
```
Out[123]= \left\{\left\{\text{j}\to 0,\text{fit}\to 0.30352429540663584347487201\,x-3.7928443631265547617782299,\text{maxResidual}\to 3.0279243439949184754414068\right\},\left\{\text{j}\to 1,\text{fit}\to -0.086862943083382048147951321\,x-5.4342299389862763407961899,\text{maxResidual}\to 1.2544507626082797470372535\right\},\left\{\text{j}\to 2,\text{fit}\to 4.0000000000000000000000000-1.3704477937… [489 chars total]
```
## 2026-07-19 18:06:43 — 📥 BULK INSERT 1 CELL(S) at positions 48–48
1. [code] (*User's insight (opus): the eigenvalue of tau1 = t[1,1] is literally built from the quadratic↵  Cas


## 2026-07-19 18:06:43
**Input:**
```wolfram
(*User's insight (opus): the eigenvalue of tau1 = t[1,1] is literally built from the quadratic
  Casimir. Derivation at z=1: tau1(u) = tr_aux[L1(u-th1) L2(u-th2)], with L[λ][i,j][v]=δij v - h JJ[j,i,λ]:
    Lambda_1(u) = 4(u-th1)(u-th2) - 2P h(u-th1) - 2S h(u-th2) + (h^2/2)(C2(mu) - C2(lam1) - C2(lam2)),
  using tr(generators)=2S,2P for [S,S,0,0]/[P,P,0,0] and the split-Casimir identity
    Sum_ab E1[a,b]E2[b,a] = (C2tot - C2(lam1) - C2(lam2))/2.
  This IS affine in C2(mu) (slope h^2/2). The earlier "falsification" was a pooled-fit artifact:
  the u^1 coeff is constant in mu but varies with S+P, so pooling all (S,P) under one intercept hid
  the true per-(S,P) affine-in-C2 structure. Here we test the FULL closed form directly, per state.*)
ClearAll[Csite, Lam1pred];
Csite[λ_] := Csite[λ] = Module[{m = Sum[Ee[λ][a, b] . Ee[λ][b, a], {a, 4}, {b, 4}]},
   If[m =!= m[[1, 1]] IdentityMatrix[Length[m]], Print["Csite not scalar for ", λ]; Abort[]];
   m[[1, 1]]];
Lam1pred[S_, P_, C2mu_][u_] := 4 (u - θ1) (u - θ2) - 2 P h (u - θ1) - 2 S h (u - θ2) + (h^2/2) (C2mu - Csite[{S, S, 0, 0}] - Csite[{P, P, 0, 0}]);
Module[{sweep = {{1, 1}, {2, 1}, {1, 2}, {2, 2}}, worst},
  worst = Max@Table[
    Module[{S = SP[[1]], P = SP[[2]], λ1, λ2, es},
      λ1 = {S, S, 0, 0}; λ2 = {P, P, 0, 0};
      es = TauEigensystemCasimir[λ1, λ2];
      Max@KeyValueMap[Function[{key, rec},
         Max@Table[Abs[N[LambdaU[λ1, λ2][key[[1]], key[[2]], key[[3]], key[[4]], key[[5]]][1][u0]
                        - Lam1pred[S, P, rec["C2"]][u0], 30]], {u0, {2/3, 9/5, 4, 23/7}}]], es]],
    {SP, sweep}];
  If[worst > 10^-8, Print["CLOSED FORM FAILS: ", worst]; Abort[]];
  {"tau1 = quadratic-Casimir closed form VERIFIED over full sweep, worst" -> worst}]
```
**Output:**
```
Out[127]= \left\{\text{tau1 = quadratic-Casimir closed form VERIFIED over full sweep, worst}\to 0.\times 10^{-25}\right\}
```
## 2026-07-19 18:07:11 — ✏️ EDIT CELL 47
(*CONCLUSION -- SUPERSEDED by the closed-form verification two cells below (kept for the record).
  The two fits above (two-anchor + least-squares) show that a SINGLE affine relation
  "coeff_j(Lambda…


## 2026-07-19 18:07:11
**Input:**
```wolfram
(*CONCLUSION -- SUPERSEDED by the closed-form verification two cells below (kept for the record).
  The two fits above (two-anchor + least-squares) show that a SINGLE affine relation
  "coeff_j(Lambda_1) = A_j + B_j C2(mu)" POOLED across all (S,P) does NOT hold for j=0,1. But this
  is a POOLING ARTIFACT, not a real falsification of the Casimir connection: the correct closed form
  (derived from tau1 = tr_aux[L1 L2] and verified exactly below) IS affine in C2(mu) with slope
  h^2/2, but its u^1 coefficient is constant in mu while varying with S+P, and its u^0 intercept
  varies with (S,P) too. Pooling different (S,P) under one intercept therefore hides the true
  per-(S,P) affine-in-C2 structure. The physically-correct statement is:
    Lambda_1(u) = 4(u-th1)(u-th2) - 2P h(u-th1) - 2S h(u-th2) + (h^2/2)(C2(mu) - C2(lam1) - C2(lam2)),
  VERIFIED exactly over the full (S,P) sweep in the cell below (worst ~ 1e-25). So tau1's eigenvalue
  IS literally the quadratic Casimir (plus a fixed (S,P)-dependent pseudovacuum polynomial in u).*)
{"pooled-single-intercept affine fit: fails for j=0,1 (a pooling artifact, NOT a real negative)",
 "RESOLVED below: full closed form Lambda_1(u) affine in C2(mu), slope h^2/2, verified exactly"}
```
**Output:**
```
Out= {"pooled-single-intercept affine fit: fails for j=0,1 (a pooling artifact, NOT a real negative)", "RESOLVED below: full closed form Lambda_1(u) affine in C2(mu), slope h^2/2, verified exactly"}
```

## 2026-07-20 05:01:34
**Input:**
```wolfram
{Dimensions[EE[{1,1,0,0},{1,1,0,0}][1,2]], C2op[{1,1,0,0},{1,1,0,0}]//Dimensions, Csite[{1,1,0,0}]}
```
**Output:**
```
Out= {{36, 36}, {36, 36}, 6}
```
## 2026-07-20 05:02:21 — 📥 BULK INSERT 1 CELL(S) at positions 49–49
1. [code] (*Total cubic/quartic Casimirs on the two-site space, built via block products for speed:↵  M[i,j]=E


## 2026-07-20 05:02:21
**Input:**
```wolfram
(*Total cubic/quartic Casimirs on the two-site space, built via block products for speed:
  M[i,j]=EE[i,j]; A[i,k]=Sum_j M[i,j].M[j,k]; C3=Sum_ik A[i,k].M[k,i];
  B[i,l]=Sum_k A[i,k].M[k,l]; C4=Sum_il B[i,l].M[l,i]. Exact rational, memoized.*)
ClearAll[C3op, C4op, C3site, C4site];
C3op[λ1_, λ2_] := C3op[λ1, λ2] = Module[{M, A},
   M = Table[EE[λ1, λ2][i, j], {i, 4}, {j, 4}];
   A = Table[Sum[M[[i, j]] . M[[j, k]], {j, 4}], {i, 4}, {k, 4}];
   Sum[A[[i, k]] . M[[k, i]], {i, 4}, {k, 4}]];
C4op[λ1_, λ2_] := C4op[λ1, λ2] = Module[{M, A, B},
   M = Table[EE[λ1, λ2][i, j], {i, 4}, {j, 4}];
   A = Table[Sum[M[[i, j]] . M[[j, k]], {j, 4}], {i, 4}, {k, 4}];
   B = Table[Sum[A[[i, k]] . M[[k, l]], {k, 4}], {i, 4}, {l, 4}];
   Sum[B[[i, l]] . M[[l, i]], {i, 4}, {l, 4}]];
(*single-site cubic/quartic Casimir scalars (asserting each is a multiple of identity)*)
C3site[λ_] := C3site[λ] = Module[{m = Sum[Ee[λ][a, b] . Ee[λ][b, c] . Ee[λ][c, a], {a, 4}, {b, 4}, {c, 4}]},
   If[m =!= m[[1, 1]] IdentityMatrix[Length[m]], Print["C3site not scalar for ", λ]; Abort[]]; m[[1, 1]]];
C4site[λ_] := C4site[λ] = Module[{m = Sum[Ee[λ][a, b] . Ee[λ][b, c] . Ee[λ][c, d] . Ee[λ][d, a], {a, 4}, {b, 4}, {c, 4}, {d, 4}]},
   If[m =!= m[[1, 1]] IdentityMatrix[Length[m]], Print["C4site not scalar for ", λ]; Abort[]]; m[[1, 1]]];
(*smoke test on the small rep: dims + site scalars*)
{Dimensions[C3op[{1, 1, 0, 0}, {1, 1, 0, 0}]], Dimensions[C4op[{1, 1, 0, 0}, {1, 1, 0, 0}]],
 "C3site11" -> C3site[{1, 1, 0, 0}], "C4site11" -> C4site[{1, 1, 0, 0}],
 "C2site11" -> Csite[{1, 1, 0, 0}]}
```
**Output:**
```
Out[133]= \left\{\left\{36,36\right\},\left\{36,36\right\},\text{C3site11}\to 18,\text{C4site11}\to 54,\text{C2site11}\to 6\right\}
```
## 2026-07-20 05:02:34 — 📥 BULK INSERT 1 CELL(S) at positions 50–50
1. [code] (*Sanity: C3op, C4op must commute with the whole transfer-matrix coefficient family (tgCoeffs)↵  AND


## 2026-07-20 05:02:34
**Input:**
```wolfram
(*Sanity: C3op, C4op must commute with the whole transfer-matrix coefficient family (tgCoeffs)
  AND with the Cartan generators, on {1,1,0,0}x{1,1,0,0}. Exact rational -> commutators are 0.
  If nonzero, the "tau_a -> C_{a+1}" conjecture is ill-posed on shared eigenvectors.*)
Module[{λ1 = {1, 1, 0, 0}, λ2 = {1, 1, 0, 0}, c1, c2, c3, fam, comm3, comm4},
  {c1, c2, c3} = tgCoeffs[λ1, λ2];
  fam = Join[c1, c2, c3, {H1op[λ1, λ2], H2op[λ1, λ2], H3op[λ1, λ2], C2op[λ1, λ2]}];
  comm3 = Max[Abs[Flatten[#. C3op[λ1, λ2] - C3op[λ1, λ2] . #]]] & /@ fam // Max;
  comm4 = Max[Abs[Flatten[#. C4op[λ1, λ2] - C4op[λ1, λ2] . #]]] & /@ fam // Max;
  {"maxComm[C3op, family]" -> comm3, "maxComm[C4op, family]" -> comm4,
   "commute" -> (comm3 === 0 && comm4 === 0)}]
```
**Output:**
```
Out[134]= \left\{\text{maxComm[C3op, family]}\to 0,\text{maxComm[C4op, family]}\to 0,\text{commute}\to \mathrm{True}\right\}
```

## 2026-07-20 05:03:18
**Input:**
```wolfram
{AbsoluteTiming[Dimensions[C3op[{2,2,0,0},{2,2,0,0}]]], AbsoluteTiming[Dimensions[C4op[{2,2,0,0},{2,2,0,0}]]]}
```
**Output:**
```
Out= {{9.120411, {400, 400}}, {17.773432, {400, 400}}}
```
## 2026-07-20 05:03:33 — 📥 BULK INSERT 1 CELL(S) at positions 51–51
1. [code] (*Per-state Casimir data for the tau2/tau3 study. For every state (all sectors are multiplicity↵  fr


## 2026-07-20 05:03:42
**Input:**
```wolfram
(*Per-state Casimir data for the tau2/tau3 study. For every state (all sectors are multiplicity
  free per cell 44) collect {S,P, C2mu, C3mu, C4mu, t2c(len5), t3c(len7)} where C3mu/C4mu are
  Rayleigh quotients of C3op/C4op on the stored numeric eigenvector, and t2c/t3c are the exact
  u-coefficient lists of Lambda_2 / Lambda_3 already stored in TauEigensystemCasimir.*)
ClearAll[rayleigh];
rayleigh[m_, v_] := (v . (m . v))/(v . v);
casData23 = Flatten[Table[
   Module[{S = SP[[1]], P = SP[[2]], λ1, λ2, es, c3m, c4m},
     λ1 = {S, S, 0, 0}; λ2 = {P, P, 0, 0};
     es = TauEigensystemCasimir[λ1, λ2];
     c3m = C3op[λ1, λ2]; c4m = C4op[λ1, λ2];
     KeyValueMap[Function[{key, rec},
        {S, P, rec["C2"], Re[rayleigh[c3m, rec["vec"]]], Re[rayleigh[c4m, rec["vec"]]],
         rec["t2c"], rec["t3c"]}], es]],
   {SP, {{1, 1}, {2, 1}, {1, 2}, {2, 2}}}], 1];
{"nRows" -> Length[casData23],
 "sampleRow(S,P,C2,C3,C4)" -> N[casData23[[1, 1 ;; 5]], 12],
 "distinctC3perSP" -> (Length[DeleteDuplicatesBy[Select[casData23, {#[[1]], #[[2]]} == {2, 2} &], Round[#[[4]], 10^-6] &]])}
```
**Output:**
```
Out[138]= \left\{\text{nRows}\to 676,\text{sampleRow(S,P,C2,C3,C4)}\to \left\{1.\,,1.\,,16.\,,64.\,,256.\,\right\},\text{distinctC3perSP}\to 6\right\}
```
## 2026-07-20 05:04:36 — 📥 BULK INSERT 1 CELL(S) at positions 52–52
1. [code] (*Fitting harness. Per FIXED (S,P) (never pooled): distinct-mu rows (dedup by rounded C2,C3,C4);↵  f


## 2026-07-20 05:04:37
**Input:**
```wolfram
(*Fitting harness. Per FIXED (S,P) (never pooled): distinct-mu rows (dedup by rounded C2,C3,C4);
  for each u^k coefficient regress against escalating Casimir-monomial bases and report the FIRST
  (minimal) basis whose max residual is ~0. Data is ~30-digit exact-derived, so a correct basis
  gives residual ~1e-20; a missing term leaves O(1) residual.*)
ClearAll[distinctMu, fitResid, minimalBasis];
distinctMu[S_, P_] := DeleteDuplicatesBy[Select[casData23, {#[[1]], #[[2]]} == {S, P} &], Round[{#[[3]], #[[4]], #[[5]]}, 10^-6] &];
(* basis element = pure function of (c2,c3,c4); build design matrix, LeastSquares, max residual *)
fitResid[rows_, ys_, basis_] := Module[{X, sol},
   X = Table[Through[basis[r[[3]], r[[4]], r[[5]]]], {r, rows}];
   sol = LeastSquares[X, ys];
   {sol, Max[Abs[X . sol - ys]]}];
(* try a list of {name, basisFuncs} in order; return first with residual < tol (and its fit) *)
minimalBasis[rows_, ys_, candidates_, tol_: 10^-12] := Module[{hit},
   hit = SelectFirst[candidates, fitResid[rows, ys, #[[2]]][[2]] < tol &, Missing["none"]];
   If[MissingQ[hit], {"NONE-fit-under-tol", "worstOverAll" -> Min[fitResid[rows, ys, #[[2]]][[2]] & /@ candidates]},
     Module[{f = fitResid[rows, ys, hit[[2]]]}, {hit[[1]], "coeffs" -> Chop[f[[1]], 10^-9], "resid" -> f[[2]]}]]];
(* monomial pure-functions in c2,c3,c4 *)
b1 = Function[{c2, c3, c4}, 1]; bC2 = Function[{c2, c3, c4}, c2]; bC3 = Function[{c2, c3, c4}, c3];
bC4 = Function[{c2, c3, c4}, c4]; bC2sq = Function[{c2, c3, c4}, c2^2]; bC2C3 = Function[{c2, c3, c4}, c2 c3];
bC2cu = Function[{c2, c3, c4}, c2^3]; bC3sq = Function[{c2, c3, c4}, c3^2]; bC2sqC3 = Function[{c2, c3, c4}, c2^2 c3];
(* escalation ladder for tau2 (expect linear-in-C3) *)
tau2Cands = {
   {"{1,C2}", {b1, bC2}},
   {"{1,C2,C3}", {b1, bC2, bC3}},
   {"{1,C2,C2^2}", {b1, bC2, bC2sq}},
   {"{1,C2,C2^2,C3}", {b1, bC2, bC2sq, bC3}},
   {"{1,C2,C2^2,C3,C2*C3}", {b1, bC2, bC2sq, bC3, bC2C3}}};
(* run tau2 fit for every (S,P) and every u^k, k=0..4 (t2c is field 6) *)
tau2Fits = Table[
   Module[{S = SP[[1]], P = SP[[2]], rows}, rows = distinctMu[S, P];
     {S, P, "nMu" -> Length[rows],
      Table[k -> minimalBasis[rows, rows[[All, 6, k + 1]], tau2Cands], {k, 0, 4}]}],
   {SP, {{1, 1}, {2, 1}, {1, 2}, {2, 2}}}];
Column[tau2Fits, Dividers -> All]
```
**Output:**
```
Out[148]= \begin{matrix}\left\{1,1,\text{nMu}\to 3,\left\{0\to \left\{\text{\{1,C2,C2\^{}2\}},\text{coeffs}\to \left\{16.2630385487528344671201814,-4.47619047619047619047619048,0.25\right\},\text{resid}\to 0.\times 10^{-24}\right\},1\to \left\{\text{\{1,C2\}},\text{coeffs}\to \left\{9.1972789115646258503401361,-2.4761904761904761904761905\right\},\text{resid}\to 0.\times 10^{-24}\right\},2\to \lef… [1247 chars total]
```

## 2026-07-20 05:04:47
**Input:**
```wolfram
Grid[Prepend[Table[{SP[[1]], SP[[2]], Length[distinctMu[SP[[1]], SP[[2]]]]}, {SP, {{1,1},{2,1},{1,2},{2,2}}}], {"S","P","nDistinctMu"}], Frame->All]
```
**Output:**
```
Out= Grid[{{"S", "P", "nDistinctMu"}, {1, 1, 3}, {2, 1, 3}, {1, 2, 3}, {2, 2, 6}}, Frame -> All]
```

## 2026-07-20 05:05:30
**Input:**
```wolfram
Module[{rows=distinctMu[2,2],tab},
 tab=Table[{Round[r[[3]],1/100],Round[r[[4]],1/100],Round[r[[5]],1/100]},{r,rows}];
 {"C2,C3,C4 for 6 mu"->tab,
  "rank{1,C2}"->MatrixRank[N[({1,#[[1]]}&/@tab),30]],
  "rank{1,C2,C2^2}"->MatrixRank[N[({1,#[[1]],#[[1]]^2}&/@tab),30]],
  "rank{1,C2,C3}"->MatrixRank[N[({1,#[[1]],#[[2]]}&/@tab),30]],
  "rank{1,C2,C2^2,C3}"->MatrixRank[N[({1,#[[1]],#[[1]]^2,#[[2]]}&/@tab),30]],
  "rank{1,C2,C2^2,C3,C4}"->MatrixRank[N[({1,#[[1]],#[[1]]^2,#[[2]],#[[3]]}&/@tab),30]],
  "rank{1,C2,C2^2,C3,C2C3,C4}"->MatrixRank[N[({1,#[[1]],#[[1]]^2,#[[2]],#[[1]]#[[2]],#[[3]]}&/@tab),30]]}]
```
**Output:**
```
Out= {"C2,C3,C4 for 6 mu" -> {{48, 288, 1728}, {40, 224, 1336}, {36, 192, 1164}, {28, 128, 628}, {24, 96, 456}, {16, 32, 64}}, "rank{1,C2}" -> 2, "rank{1,C2,C2^2}" -> 3, "rank{1,C2,C3}" -> 2, "rank{1,C2,C2^2,C3}" -> 3, "rank{1,C2,C2^2,C3,C4}" -> 4, "rank{1,C2,C2^2,C3,C2C3,C4}" -> 4}
```

## 2026-07-20 05:06:29
**Input:**
```wolfram
Table[Module[{S=SP[[1]],P=SP[[2]],rows,tab,c3fit,c4fit},
  rows=distinctMu[S,P];
  tab=Table[{Round[r[[3]],1/100],Round[r[[4]],1/100],Round[r[[5]],1/100]},{r,rows}];
  c3fit=Fit[N[({#[[1]],#[[2]]}&/@tab),30],{1,x},x];  (* C3 vs C2 *)
  c4fit=Fit[N[({#[[1]],#[[3]]}&/@tab),30],{1,x,x^2},x]; (* C4 vs C2 quadratic *)
  {S,P,"nMu"->Length[tab],"C2C3C4"->tab,
   "C3=aff(C2)?resid"->Max[Abs[(c3fit/.x->#[[1]])-#[[2]]&/@tab]],"C3fit"->c3fit,
   "C4=quad(C2)?resid"->Max[Abs[(c4fit/.x->#[[1]])-#[[3]]&/@tab]],"C4fit"->Chop[c4fit,1/10^6]}],
 {SP,{{1,1},{2,1},{1,2},{2,2}}}]//Column
```
**Output:**
```
Out= Column[{{1, 1, "nMu" -> 3, "C2C3C4" -> {{16, 64, 256}, {4, 4, 4}, {12, 44, 188}}, "C3=aff(C2)?resid" -> 0``27.540607512240765, "C3fit" -> -16.`29.522878745280337 + 5.0000000000000000000000000000000000000001440833806000222844`29.52287874528033*x, "C4=quad(C2)?resid" -> 0``26.354970935278846, "C4fit" -> -112.`29.221848749616356 + 31.`29.22184874961635*x - 0.5000000000000000000000000000000000000000000000000000000018`29.221848749616356*x^2}, {2, 1, "nMu" -> 3, "C2C3C4" -> {{30, 150, 750}, {24, 111, 564}, {14, 46, 174}}, "C3=aff(C2)?resid" -> 0``27.142667503568727, "C3fit" -> -45.`29.522878745280337 + 6.5000000000000000000000000000000000000000359345675526995484`29.522878745280337*x, "C4=quad(C2)?resid" -> 0``25.7856861025756, "C4fit" -> -540.`29.221848749616356 + 58.`29.221848749616356*x - 0.4999999999999999999999999999999999999999986474472051553739`29.221848749616356*x^2}, {1, 2, "nMu" -> 3, "C2C3C4" -> {{30, 150, 750}, {24, 111, 564}, {14, 46, 174}}, "C3=aff(C2)?resid" -> 0``27.142667503568727, "C3fit" -> -45.`29.522878745280337 + 6.5000000000000000000000000000000000000000359345675526995484`29.522878745280337*x, "C4=quad(C2)?resid" -> 0``25.7856861025756, "C4fit" -> -540.`29.221848749616356 + 58.`29.221848749616356*x - 0.4999999999999999999999999999999999999999986474472051553739`29.221848749616356*x^2}, {2, 2, "nMu" -> 6, "C2C3C4" -> {{48, 288, 1728}, {40, 224, 1336}, {36, 192, 1164}, {28, 128, 628}, {24, 96, 456}, {16, 32, 64}}, "C3=aff(C2)?resid" -> 0``26.841637507904746, "C3fit" -> -96.`29.522878745280337 + 7.9999999999999999999999999999999999999999644044453804978406`29.522878745280337*x, "C4=quad(C2)?resid" -> 54.85714285714285714285714285714285714287`27.525328548289675, "C4fit" -> -809.14285714285714285714285714285714285714`29.221848749616356 + 53.28571428571428571428571428571428571429`29.221848749616356*x}}]
```
## 2026-07-20 05:07:39 — 📥 BULK INSERT 1 CELL(S) at positions 53–53
1. [code] (*KEY STRUCTURAL FACT (cell above): within each (S,P), C3mu is EXACTLY affine in C2mu↵  (C3=5C2-16,


## 2026-07-20 05:07:40
**Input:**
```wolfram
(*KEY STRUCTURAL FACT (cell above): within each (S,P), C3mu is EXACTLY affine in C2mu
  (C3=5C2-16, 6.5C2-45, 6.5C2-45, 8C2-96 for (1,1),(2,1),(1,2),(2,2)). So C3 is REDUNDANT with C2
  on this rep family -> "linear in C3" is indistinguishable from "polynomial in C2". But C4mu is
  genuinely independent (quadratic-in-C2 fit fails at (2,2), resid ~55). Decisive test uses (2,2)
  (6 distinct mu). Report max residual of each u^k coefficient of Lambda_2 / Lambda_3 vs each basis;
  bases kept <=5 cols so a 6-point fit is never trivially exact.*)
Module[{rows = distinctMu[2, 2], bases, lam2res, lam3res},
  bases = {
    {"{1,C2}", {b1, bC2}},
    {"{1,C2,C2^2}", {b1, bC2, bC2sq}},
    {"{1,C2,C2^2,C2^3}", {b1, bC2, bC2sq, bC2cu}},
    {"{1,C2,C2^2,C4}", {b1, bC2, bC2sq, bC4}},
    {"{1,C2,C2^2,C4,C2^3}", {b1, bC2, bC2sq, bC4, bC2cu}},
    {"{1,C2,C2^2,C4,C4^2}", {b1, bC2, bC2sq, bC4, Function[{c2, c3, c4}, c4^2]}}};
  lam2res = Table[{"Lam2 u^" <> ToString[k],
      Sequence @@ (Round[Log10[fitResid[rows, rows[[All, 6, k + 1]], #[[2]]][[2]] + 10^-40], 1/10] & /@ bases)}, {k, 0, 4}];
  lam3res = Table[{"Lam3 u^" <> ToString[k],
      Sequence @@ (Round[Log10[fitResid[rows, rows[[All, 7, k + 1]], #[[2]]][[2]] + 10^-40], 1/10] & /@ bases)}, {k, 0, 6}];
  Grid[Join[{Prepend[bases[[All, 1]], "log10 maxResid @(2,2)"]}, lam2res, {Table["-", {Length[bases] + 1}]}, lam3res], Frame -> All]]
```
**Output:**
```
Out[149]= \begin{array}{|c|c|c|c|c|c|c|}\hline\text{log10 maxResid @(2,2)} & \text{\{1,C2\}} & \text{\{1,C2,C2\^{}2\}} & \text{\{1,C2,C2\^{}2,C2\^{}3\}} & \text{\{1,C2,C2\^{}2,C4\}} & \text{\{1,C2,C2\^{}2,C4,C2\^{}3\}} & \text{\{1,C2,C2\^{}2,C4,C4\^{}2\}}\\\hline\text{Lam2 u\^{}0} & \frac{7}{5} & \frac{11}{10} & \frac{9}{10} & \mathrm{Indeterminate} & \mathrm{Indeterminate} & \mathrm{Indeterminate… [1247 chars total]
```

## 2026-07-20 05:07:53
**Input:**
```wolfram
Module[{rows=distinctMu[2,2],bases,fmt},
 bases={{"{1,C2}",{b1,bC2}},{"{1,C2,C2^2}",{b1,bC2,bC2sq}},{"{1,C2,C2^2,C2^3}",{b1,bC2,bC2sq,bC2cu}},{"{1,C2,C2^2,C4}",{b1,bC2,bC2sq,bC4}},{"{1,C2,C2^2,C4,C2^3}",{b1,bC2,bC2sq,bC4,bC2cu}}};
 fmt[r_]:=If[r<10^-15,0,ScientificForm[N[r],2]];
 {"Lam2"->Table[("u^"<>ToString[k])->(fmt[fitResid[rows,rows[[All,6,k+1]],#[[2]]][[2]]]&/@bases),{k,0,4}],
  "Lam3"->Table[("u^"<>ToString[k])->(fmt[fitResid[rows,rows[[All,7,k+1]],#[[2]]][[2]]]&/@bases),{k,0,6}],
  "basisOrder"->bases[[All,1]]}]//Column
```
**Output:**
```
Out= Column[{"Lam2" -> {"u^0" -> {ScientificForm[25.714285714285715, 2], ScientificForm[13.714285714285714, 2], ScientificForm[8.067226890756302, 2], 0, 0}, "u^1" -> {0, 0, 0, 0, 0}, "u^2" -> {0, 0, 0, 0, 0}, "u^3" -> {0, 0, 0, 0, 0}, "u^4" -> {0, 0, 0, 0, 0}}, "Lam3" -> {"u^0" -> {0, 0, 0, 0, 0}, "u^1" -> {0, 0, 0, 0, 0}, "u^2" -> {0, 0, 0, 0, 0}, "u^3" -> {0, 0, 0, 0, 0}, "u^4" -> {0, 0, 0, 0, 0}, "u^5" -> {0, 0, 0, 0, 0}, "u^6" -> {0, 0, 0, 0, 0}}, "basisOrder" -> {"{1,C2}", "{1,C2,C2^2}", "{1,C2,C2^2,C2^3}", "{1,C2,C2^2,C4}", "{1,C2,C2^2,C4,C2^3}"}}]
```

## 2026-07-20 05:08:55
**Input:**
```wolfram
Module[{r22=distinctMu[2,2]},
 {"tau2 u^0 fit {1,C2,C2^2,C4} @(2,2)"->Chop[fitResid[r22,r22[[All,6,1]],{b1,bC2,bC2sq,bC4}][[1]],10^-8],
  "tau2 u^0 alt {1,C2,C4} resid"->fitResid[r22,r22[[All,6,1]],{b1,bC2,bC4}][[2]],
  "tau2 u^0 {1,C2,C4} fit"->Chop[fitResid[r22,r22[[All,6,1]],{b1,bC2,bC4}][[1]],10^-8],
  "--- corroborate tau3 all-affine-in-C2 at 3-mu (S,P): resid {1,C2} per k ---"->
   Table[{SP[[1]],SP[[2]]}->Table[Chop[fitResid[distinctMu[SP[[1]],SP[[2]]],distinctMu[SP[[1]],SP[[2]]][[All,7,k+1]],{b1,bC2}][[2]]],{k,0,6}],{SP,{{1,1},{2,1},{1,2}}}],
  "--- tau2 u^0: does C4 appear at 3-mu too? test {1,C2,C4} vs {1,C2,C2^2} resid ---"->
   Table[{SP[[1]],SP[[2]]}->{"{1,C2}res"->Chop[fitResid[distinctMu@@SP,distinctMu[SP[[1]],SP[[2]]][[All,6,1]],{b1,bC2}][[2]]]},{SP,{{1,1},{2,1},{1,2}}}]}]//Column
```
**Output:**
```
Out= Column[{"tau2 u^0 fit {1,C2,C2^2,C4} @(2,2)" -> {-114.12244897959183673469387755102040816427`25.800781371936367, 6.51190476190476190476190476190476190466`25.800781371936367, 0.125`25.800781371936367, -0.24999999999999999999999999999999999999`25.800781371936367}, "tau2 u^0 alt {1,C2,C4} resid" -> 17.99999999999999999999999999999999999995`24.372503306307593, "tau2 u^0 {1,C2,C4} fit" -> {-228.12244897959183673469387755102040819427`25.800781371936367, 14.51190476190476190476190476190476190648`25.800781371936367, -0.25000000000000000000000000000000000003`25.800781371936367}, "--- corroborate tau3 all-affine-in-C2 at 3-mu (S,P): resid {1,C2} per k ---" -> {{1, 1} -> {0, 0, 0, 0, 0, 0, 0}, {2, 1} -> {0, 0, 0, 0, 0, 0, 0}, {1, 2} -> {0, 0, 0, 0, 0, 0, 0}}, "--- tau2 u^0: does C4 appear at 3-mu too? test {1,C2,C4} vs {1,C2,C2^2} resid ---" -> {{1, 1} -> {"{1,C2}res" -> 5.14285714285714285714285714285714285714`26.531293691252756}, {2, 1} -> {"{1,C2}res" -> 9.79591836734693877551020408162280208649`26.24876720199206}, {1, 2} -> {"{1,C2}res" -> 9.79591836734693877551020408163265306114`25.72818323248081}}}]
```

## 2026-07-20 05:10:12
**Input:**
```wolfram
ratfit[rows_,ys_,basis_]:=Module[{X,sol},X=Table[Through[basis[Rationalize[Round[r[[3]],10^-9],10^-6],Rationalize[Round[r[[4]],10^-9],10^-6],Rationalize[Round[r[[5]],10^-9],10^-6]]],{r,rows}];sol=LeastSquares[N[X,40],ys];Rationalize[Round[Re[sol],10^-12],10^-8]];
Module[{},
 {"Lam3: {A_k,B_k} with basis {1,C2} per (S,P), k=0..6"->
  Table[{SP[[1]],SP[[2]]}->Table[ratfit[distinctMu@@SP,distinctMu[SP[[1]],SP[[2]]][[All,7,k+1]],{b1,bC2}],{k,0,6}],{SP,{{1,1},{2,1},{1,2},{2,2}}}]}]//Column
```
**Output:**
```
Out= Column[{"Lam3: {A_k,B_k} with basis {1,C2} per (S,P), k=0..6" -> {{1, 1} -> {{2518086599719/1000000000000, 15589569161/62500000000}, {-16837058632977/500000000000, -1324263038549/500000000000}, {36084062196307/250000000000, 6803854875283/1000000000000}, {-122746571644531/500000000000, -347619047619/100000000000}, {78074829931973/500000000000, 1/2}, {-20857142857143/500000000000, 0}, {4, 0}}, {2, 1} -> {{378306878307/250000000000, 324263038549/1000000000000}, {-6569403952057/250000000000, -684126984127/200000000000}, {75799319727891/500000000000, 8613378684807/1000000000000}, {-4998564207429/15625000000, -397619047619/100000000000}, {98860544217687/500000000000, 1/2}, {-23857142857143/500000000000, 0}, {4, 0}}, {1, 2} -> {{892452218983/500000000000, 328798185941/1000000000000}, {-14637134218767/500000000000, -1732993197279/500000000000}, {160157110463233/1000000000000, 1741723356009/200000000000}, {-326329877982939/1000000000000, -397619047619/100000000000}, {198863945578231/1000000000000, 1/2}, {-23857142857143/500000000000, 0}, {4, 0}}, {2, 2} -> {{0, 427437641723/1000000000000}, {-15306338408379/1000000000000, -447619047619/100000000000}, {163709750566893/1000000000000, 11018140589569/1000000000000}, {-430363891588381/1000000000000, -447619047619/100000000000}, {12421768707483/50000000000, 1/2}, {-26857142857143/500000000000, 0}, {4, 0}}}}]
```

## 2026-07-20 05:11:25
**Input:**
```wolfram
Module[{rows=distinctMu[2,2],X,ys,sol},
 ys=rows[[All,7,2]]; (*Lam3 u^1 coeff*)
 X=Table[{1,Round[r[[3]],10^-6]},{r,rows}];
 sol=LeastSquares[N[X,50],ys];
 {"raw sol"->sol//N,"resid"->Max[Abs[X.sol-ys]],"Rationalize e-20"->Rationalize[sol,10^-20]}]
```
**Output:**
```
Out= {"raw sol" -> {-15.306338408379224, -4.476190476190476}, "resid" -> 0``24.978429111022567, "Rationalize e-20" -> {-141752/9261, -94/21}}
```
## 2026-07-20 05:13:01 — 📥 BULK INSERT 1 CELL(S) at positions 54–54
1. [code] (*Exact per-(S,P) Casimir decomposition of the Lambda_2 / Lambda_3 coefficients. Residual is↵  EXACT


## 2026-07-20 05:13:01
**Input:**
```wolfram
(*Exact per-(S,P) Casimir decomposition of the Lambda_2 / Lambda_3 coefficients. Residual is
  EXACTLY 0 over all distinct mu per (S,P) (verified in the two cells above), and the tau-coeff
  matrices are simultaneously diagonalized with C2op/C3op/C4op, so these eigenvalue fits ARE exact
  operator identities on each (S,P) block. Rationalize[.,10^-20] on the high-precision LeastSquares
  solution recovers the exact rationals.*)
ClearAll[ratCoef];
ratCoef[rows_, ys_, basis_] := Module[{X, sol},
   X = Table[Through[basis[Round[r[[3]], 10^-6], Round[r[[4]], 10^-6], Round[r[[5]], 10^-6]]], {r, rows}];
   sol = LeastSquares[N[X, 50], ys];
   {Rationalize[Re[sol], 10^-18], Max[Abs[X . sol - ys]]}];
lam3decomp = Table[Module[{S = SP[[1]], P = SP[[2]], rows = distinctMu[SP[[1]], SP[[2]]]},
    {S, P, "Lam3 u^k -> {A_k,B_k}(basis 1,C2mu)" ->
      Table[k -> ratCoef[rows, rows[[All, 7, k + 1]], {b1, bC2}][[1]], {k, 0, 6}]}],
   {SP, {{1, 1}, {2, 1}, {1, 2}, {2, 2}}}];
lam2decomp = Table[Module[{S = SP[[1]], P = SP[[2]], rows = distinctMu[SP[[1]], SP[[2]]]},
    {S, P,
     "Lam2 u^1..4 -> {A_k,B_k}(1,C2mu)" -> Table[k -> ratCoef[rows, rows[[All, 6, k + 1]], {b1, bC2}][[1]], {k, 1, 4}],
     "Lam2 u^0 -> {a,b,c,d}(1,C2mu,C2mu^2,C4mu)" -> ratCoef[rows, rows[[All, 6, 1]], {b1, bC2, bC2sq, bC4}][[1]]}],
   {SP, {{1, 1}, {2, 1}, {1, 2}, {2, 2}}}];
Column[{"=== Lam3 ===", Column[lam3decomp], "=== Lam2 ===", Column[lam2decomp]}]
```
**Output:**
```
Out[154]= \begin{matrix}\text{=== Lam3 ===}\\\begin{matrix}\left\{1,1,\text{Lam3 u\^{}k -> \{A\_k,B\_k\}(basis 1,C2mu)}\to \left\{0\to \left\{\frac{23320}{9261},\frac{110}{441}\right\},1\to \left\{-\frac{103952}{3087},-\frac{1168}{441}\right\},2\to \left\{\frac{148522}{1029},\frac{6001}{882}\right\},3\to \left\{-\frac{2273512}{9261},-\frac{73}{21}\right\},4\to \left\{\frac{22954}{147},\frac{1}{2}\… [1247 chars total]
```

## 2026-07-20 05:13:05
**Input:**
```wolfram
Column[{"Lam3:",Column[lam3decomp/.Rule[a_,b_]:>Row[{a," -> ",b}]],"Lam2:",Column[lam2decomp]}]//InputForm
```
**Output:**
```
Out= InputForm[Column[{"Lam3:", Column[{{1, 1, Row[{"Lam3 u^k -> {A_k,B_k}(basis 1,C2mu)", " -> ", {0 -> {23320/9261, 110/441}, 1 -> {-103952/3087, -1168/441}, 2 -> {148522/1029, 6001/882}, 3 -> {-2273512/9261, -73/21}, 4 -> {22954/147, 1/2}, 5 -> {-292/7, 0}, 6 -> {4, 0}}}]}, {2, 1, Row[{"Lam3 u^k -> {A_k,B_k}(basis 1,C2mu)", " -> ", {0 -> {286/189, 143/441}, 1 -> {-81119/3087, -431/126}, 2 -> {22285/147, 7597/882}, 3 -> {-2962669/9261, -167/42}, 4 -> {29065/147, 1/2}, 5 -> {-334/7, 0}, 6 -> {4, 0}}}]}, {1, 2, Row[{"Lam3 u^k -> {A_k,B_k}(basis 1,C2mu)", " -> ", {0 -> {5510/3087, 145/441}, 1 -> {-271109/9261, -1019/294}, 2 -> {494405/3087, 7681/882}, 3 -> {-3022141/9261, -167/42}, 4 -> {29233/147, 1/2}, 5 -> {-334/7, 0}, 6 -> {4, 0}}}]}, {2, 2, Row[{"Lam3 u^k -> {A_k,B_k}(basis 1,C2mu)", " -> ", {0 -> {0, 377/882}, 1 -> {-141752/9261, -94/21}, 2 -> {72196/441, 4859/441}, 3 -> {-3985600/9261, -94/21}, 4 -> {36520/147, 1/2}, 5 -> {-376/7, 0}, 6 -> {4, 0}}}]}}], "Lam2:", Column[{{1, 1, "Lam2 u^1..4 -> {A_k,B_k}(1,C2mu)" -> {1 -> {1352/147, -52/21}, 2 -> {4862/147, 1}, 3 -> {-208/7, 0}, 4 -> {6, 0}}, "Lam2 u^0 -> {a,b,c,d}(1,C2mu,C2mu^2,C4mu)" -> {3676/529445, 31673/1361430, 17254/97245, -28229/194490}}, {2, 1, "Lam2 u^1..4 -> {A_k,B_k}(1,C2mu)" -> {1 -> {5051/147, -125/42}, 2 -> {5996/147, 1}, 3 -> {-250/7, 0}, 4 -> {6, 0}}, "Lam2 u^0 -> {a,b,c,d}(1,C2mu,C2mu^2,C4mu)" -> {119017/24777081, 7489303/173439567, 28486039/173439567, -59495411/346879134}}, {1, 2, "Lam2 u^1..4 -> {A_k,B_k}(1,C2mu)" -> {1 -> {4747/147, -125/42}, 2 -> {6164/147, 1}, 3 -> {-250/7, 0}, 4 -> {6, 0}}, "Lam2 u^0 -> {a,b,c,d}(1,C2mu,C2mu^2,C4mu)" -> {36871/8259027, 6947933/173439567, 28633103/173439567, -58907155/346879134}}, {2, 2, "Lam2 u^1..4 -> {A_k,B_k}(1,C2mu)" -> {1 -> {3212/49, -73/21}, 2 -> {7886/147, 1}, 3 -> {-292/7, 0}, 4 -> {6, 0}}, "Lam2 u^0 -> {a,b,c,d}(1,C2mu,C2mu^2,C4mu)" -> {-5592/49, 547/84, 1/8, -1/4}}}]}]]
```
## 2026-07-20 05:15:05 — 📥 BULK INSERT 1 CELL(S) at positions 55–55
1. [code] (*GLOBAL HARD-ASSERT of the Casimir content over the full (S,P) sweep (all 676 states).↵  Conjecture


## 2026-07-20 05:15:05
**Input:**
```wolfram
(*GLOBAL HARD-ASSERT of the Casimir content over the full (S,P) sweep (all 676 states).
  Conjecture "tau_a introduces C_{a+1}" tested:
   - tau3 (Lambda_3): claim EACH u^k coefficient is AFFINE in C2mu (equivalently linear in C3mu,
     since C3mu=aff(C2mu) on this family) with NO C4mu -> fit {1,C2mu}, residual must be 0.
   - tau2 (Lambda_2): claim u^1..u^4 affine in C2mu, and u^0 = [affine in C2mu] + (1/8)C2mu^2
     - (1/4)C4mu with those two coefficients UNIVERSAL (S,P-independent). We SUBTRACT the fixed
     1/8, -1/4 pieces and require the remainder to be affine in C2mu. Fixing the constants makes the
     3-distinct-mu blocks (1,1),(2,1),(1,2) a genuine 1-DOF test -> confirms universality.
  affFit residual over ALL distinct mu per (S,P) must be ~0 (data is ~30-digit exact-derived).*)
Module[{sweep = {{1, 1}, {2, 1}, {1, 2}, {2, 2}}, affResid, w2, w3, w2u0, w3c4},
  affResid[rows_, ys_] := fitResid[rows, ys, {b1, bC2}][[2]];
  (* tau3: every coeff affine in C2mu *)
  w3 = Max@Table[Module[{rows = distinctMu[SP[[1]], SP[[2]]]},
      Max@Table[affResid[rows, rows[[All, 7, k + 1]]], {k, 0, 6}]], {SP, sweep}];
  (* tau3: C4 genuinely absent -- at (2,2) the {1,C2mu} residual is already 0 (checked); assert here *)
  w3c4 = Module[{rows = distinctMu[2, 2]}, Max@Table[affResid[rows, rows[[All, 7, k + 1]]], {k, 0, 6}]];
  (* tau2: u^1..u^4 affine in C2mu *)
  w2 = Max@Table[Module[{rows = distinctMu[SP[[1]], SP[[2]]]},
      Max@Table[affResid[rows, rows[[All, 6, k + 1]]], {k, 1, 4}]], {SP, sweep}];
  (* tau2 u^0: subtract UNIVERSAL (1/8)C2mu^2 - (1/4)C4mu, remainder must be affine in C2mu *)
  w2u0 = Max@Table[Module[{rows = distinctMu[SP[[1]], SP[[2]]], ys},
      ys = (#[[6, 1]] - (1/8) #[[3]]^2 + (1/4) #[[5]]) & /@ rows;
      affResid[rows, ys]], {SP, sweep}];
  If[Max[w2, w3, w2u0] > 10^-15, Print["CASIMIR-CONTENT ASSERT FAILED: ", {w2, w3, w2u0}]; Abort[]];
  {"tau3 all coeffs affine-in-C2mu (no C4), worstResid" -> w3,
   "tau2 u^1..u^4 affine-in-C2mu, worstResid" -> w2,
   "tau2 u^0 = aff(C2mu) + (1/8)C2mu^2 - (1/4)C4mu [1/8,-1/4 FIXED universal], worstResid" -> w2u0,
   "VERDICT" -> "tau3->no C4 (linear in Casimir like tau1); tau2->reaches C4 in u^0 only"}]
```
**Output:**
```
Out[155]= \left\{\text{tau3 all coeffs affine-in-C2mu (no C4), worstResid}\to 0.\times 10^{-24},\text{tau2 u\^{}1..u\^{}4 affine-in-C2mu, worstResid}\to 0.\times 10^{-24},\text{tau2 u\^{}0 = aff(C2mu) + (1/8)C2mu\^{}2 - (1/4)C4mu [1/8,-1/4 FIXED universal], worstResid}\to 0.\times 10^{-24},\text{VERDICT}\to \text{tau3->no C4 (linear in Casimir like tau1); tau2->reaches C4 in u\^{}0 only}\right\}
```
## 2026-07-20 05:16:21 — 📥 BULK INSERT 1 CELL(S) at positions 56–56
1. [code] (*Explicit closed forms with the mu-dependence carried ENTIRELY by Casimir eigenvalues:↵   Lambda_3(


## 2026-07-20 05:16:29
**Input:**
```wolfram
(*Explicit closed forms with the mu-dependence carried ENTIRELY by Casimir eigenvalues:
   Lambda_3(u) = Sum_k (A3_k + B3_k C2mu) u^k                         (affine in C2mu, NO C4)
   Lambda_2(u) = (a20 + b20 C2mu + (1/8) C2mu^2 - (1/4) C4mu)         (u^0: the ONLY C4/C2^2 term)
               + Sum_{k=1..4} (A2_k + B2_k C2mu) u^k
  The 1/8 and -1/4 are UNIVERSAL (cell 55). The per-(S,P) pseudovacuum/slope constants {A,B,a,b}
  are exact rationals (theta-dependent, S<->P asymmetric because theta1!=theta2) tabulated below.
  Hard-assert the predictors reproduce EVERY state's Lambda_2, Lambda_3 over the full sweep.*)
ClearAll[coef3, coef2, Lam3pred, Lam2pred];
Do[Module[{S = SP[[1]], P = SP[[2]], rows = distinctMu[SP[[1]], SP[[2]]], ys0},
   coef3[S, P] = Table[ratCoef[rows, rows[[All, 7, k + 1]], {b1, bC2}][[1]], {k, 0, 6}];
   coef2[S, P] = Table[ratCoef[rows, rows[[All, 6, k + 1]], {b1, bC2}][[1]], {k, 1, 4}];
   ys0 = (#[[6, 1]] - (1/8) #[[3]]^2 + (1/4) #[[5]]) & /@ rows;   (* u^0 remainder after universal terms *)
   coef2[S, P, 0] = ratCoef[rows, ys0, {b1, bC2}][[1]];           (* {a20, b20} *)
   ], {SP, {{1, 1}, {2, 1}, {1, 2}, {2, 2}}}];
Lam3pred[S_, P_, C2m_][u_] := Sum[(coef3[S, P][[k + 1, 1]] + coef3[S, P][[k + 1, 2]] C2m) u^k, {k, 0, 6}];
Lam2pred[S_, P_, C2m_, C4m_][u_] := (coef2[S, P, 0][[1]] + coef2[S, P, 0][[2]] C2m + (1/8) C2m^2 - (1/4) C4m)
   + Sum[(coef2[S, P][[k, 1]] + coef2[S, P][[k, 2]] C2m) u^k, {k, 1, 4}];
Module[{sweep = {{1, 1}, {2, 1}, {1, 2}, {2, 2}}, w2, w3},
  {w2, w3} = Transpose@Table[Module[{S = SP[[1]], P = SP[[2]], λ1, λ2, es},
      λ1 = {S, S, 0, 0}; λ2 = {P, P, 0, 0}; es = TauEigensystemCasimir[λ1, λ2];
      {Max@KeyValueMap[Function[{key, rec}, Max@Table[Abs[N[LambdaU[λ1, λ2][key[[1]], key[[2]], key[[3]], key[[4]], key[[5]]][2][u0]
              - Lam2pred[S, P, rec["C2"], Re[rayleigh[C4op[λ1, λ2], rec["vec"]]]][u0], 30]], {u0, {2/3, 9/5, 4, 23/7}}]], es],
       Max@KeyValueMap[Function[{key, rec}, Max@Table[Abs[N[LambdaU[λ1, λ2][key[[1]], key[[2]], key[[3]], key[[4]], key[[5]]][3][u0]
              - Lam3pred[S, P, rec["C2"]][u0], 30]], {u0, {2/3, 9/5, 4, 23/7}}]], es]}], {SP, sweep}];
  w2 = Max[w2]; w3 = Max[w3];
  If[Max[w2, w3] > 10^-8, Print["CLOSED FORM FAILS: ", {w2, w3}]; Abort[]];
  {"Lam3pred (affine in C2mu, no C4) VERIFIED full sweep, worst" -> w3,
   "Lam2pred (u^0 carries +1/8 C2mu^2 -1/4 C4mu) VERIFIED full sweep, worst" -> w2}]
```
**Output:**
```
Out[160]= \left\{\text{Lam3pred (affine in C2mu, no C4) VERIFIED full sweep, worst}\to 0.\times 10^{-22},\text{Lam2pred (u\^{}0 carries +1/8 C2mu\^{}2 -1/4 C4mu) VERIFIED full sweep, worst}\to 0.\times 10^{-23}\right\}
```
## 2026-07-20 05:16:55 — 📥 BULK INSERT 1 CELL(S) at positions 57–57
1. [code] (*CONCLUSION -- test of "each tau_a introduces the next Casimir C_{a+1}".↵  Structural fact on this


## 2026-07-20 05:16:56
**Input:**
```wolfram
(*CONCLUSION -- test of "each tau_a introduces the next Casimir C_{a+1}".
  Structural fact on this rep family [S,S,0,0]x[P,P,0,0]: within each (S,P), C3mu is EXACTLY affine
  in C2mu (C3=5C2-16, 6.5C2-45, 6.5C2-45, 8C2-96 for (1,1),(2,1),(1,2),(2,2)); C4mu IS independent
  of poly(C2) at (2,2). So "linear in C3" is here indistinguishable from "affine in C2"; C4 is the
  only genuinely-new Casimir direction the transfer matrices can reach.

  tau2 (Lambda_2, deg 4): u^1..u^4 coeffs are affine in C2mu (== linear in C3mu). The CONSTANT term
   u^0 is NOT any polynomial in C2 (fails even {1,C2,C2^2,C2^3}); it requires the QUARTIC Casimir,
   entering linearly: u^0 = aff(C2mu) + (1/8) C2mu^2 - (1/4) C4mu, coefficients 1/8 and -1/4
   UNIVERSAL (verified across all (S,P), cell 55). => "tau2 -> C3" is FALSE: tau2 reaches C4, and
   only in its u^0 coefficient.

  tau3 (Lambda_3, deg 6): EVERY u^k coeff is affine in C2mu (== linear in C3mu); NO C4mu anywhere
   (verified). This mirrors tau1 (Lambda^3 C^4 ~ (Lambda^1)* x det duality). => "tau3 -> C4" is
   FALSE: tau3 stays at the C2/C3 level, like tau1.

  Net: the naive ladder "tau_a -> C_{a+1}" does NOT hold. The pattern is the self-duality one:
   tau1 & tau3 are linear-in-Casimir (affine in C2mu); the self-dual middle tau2 is the unique one
   whose eigenvalue reaches the quartic Casimir C4 (constant term only).

  Explicit verified closed forms: Lam2pred[S,P,C2mu,C4mu][u], Lam3pred[S,P,C2mu][u] (cell 56),
   reproducing every state over the full sweep to ~1e-22. The universal Casimir pieces are
   {6=chi2 leading, C2mu-slope 1 in u^2, +1/8 C2mu^2, -1/4 C4mu} for tau2 and {4=chi3 leading,
   C2mu-slope 1/2 in u^4} for tau3; the remaining per-(S,P) pseudovacuum/slope constants are exact
   rationals in (S,P,theta1,theta2,h) (S<->P asymmetric since theta1!=theta2) -- collapsing THOSE
   into one global (S,P,theta,h) formula (as cell 48 did for tau1's h^2/2 slope) needs the fusion
   trace derivation or more reps; not attempted here.*)
{"tau2 -> reaches C4 (u^0 only): +1/8 C2mu^2 - 1/4 C4mu, universal; conjecture 'linear in C3' FALSE",
 "tau3 -> affine in C2mu everywhere, NO C4; conjecture 'linear in C4' FALSE (like tau1 by duality)",
 "closed forms Lam2pred/Lam3pred verified over full sweep, worst ~1e-22"}
```
**Output:**
```
Out[161]= \left\{\text{tau2 -> reaches C4 (u\^{}0 only): +1/8 C2mu\^{}2 - 1/4 C4mu, universal; conjecture 'linear in C3' FALSE},\text{tau3 -> affine in C2mu everywhere, NO C4; conjecture 'linear in C4' FALSE (like tau1 by duality)},\text{closed forms Lam2pred/Lam3pred verified over full sweep, worst \textasciitilde{}1e-22}\right\}
```
## 2026-07-20 05:18:28 — ▶️ RUN CELL 56 (0.01 s)
**In [56]:** `(*Explicit closed forms with the mu-dependence carried ENTIRELY by Casimir eigenvalues:↵   Lambda_3(u) = Sum_k (A3_k + B3_k C2mu) u^k                         (affine in C2mu, NO C4)↵   Lambda_2(u) = (`
**Out:** Out[166]= \left\{\text{Lam3pred (affine in C2mu, no C4) VERIFIED full sweep, worst}\to 0.\times 10^{-22},\text{Lam2pred (u\^{}0 carries +1/8 C2mu\^{}2 -1/4 C4mu) VERIFIED full sweep, worst}\to 0.\times 10^{-23}\right\}


## 2026-07-20 05:19:19
**Input:**
```wolfram
Module[{λ1={2,2,0,0},λ2={2,2,0,0},es,rayleigh,rows,fitC3,fitC4lin,fitC4quad,resid},
  rayleigh[m_,v_]:=(v.(m.v))/(v.v);
  es=TauEigensystemCasimir[λ1,λ2];
  rows=DeleteDuplicatesBy[
    Function[kv,Module[{key=kv[[1]],rec=kv[[2]],v},v=rec["vec"];
      {Re[rayleigh[C2op[λ1,λ2],v]],Re[rayleigh[C3op[λ1,λ2],v]],Re[rayleigh[C4op[λ1,λ2],v]]}]]/@Normal[es],
    Round[{#[[1]],#[[3]]},10^-5]&];
  resid[data_,basis_]:=Module[{f=Fit[N[data,30],basis,x]},Max[Abs[(f/.x->#[[1]])-#[[2]]]&/@data]];
  fitC3=resid[{#[[1]],#[[2]]}&/@rows,{1,x}];
  fitC4lin=resid[{#[[1]],#[[3]]}&/@rows,{1,x}];
  fitC4quad=resid[{#[[1]],#[[3]]}&/@rows,{1,x,x^2}];
  {"nDistinctIrreps"->Length[rows],"C2vals"->Sort[rows[[All,1]]],
   "C3-affine-in-C2 residual"->fitC3,"C4-affine-in-C2 residual"->fitC4lin,"C4-quadratic-in-C2 residual"->fitC4quad}]
```
**Output:**
```
Out= {"nDistinctIrreps" -> 6, "C2vals" -> {16.00000000000000000000000000000000000076`27.532260732737054, 24.00000000000000000000000000000000000027`27.626651947482788, 28.00000000000000000000000000000000000019`27.654272634178092, 36.`27.707559286138167, 39.99999999999999999999999999999999999996`27.717404220808195, 48.`27.73859351987813}, "C3-affine-in-C2 residual" -> 0``23.98923780967441, "C4-affine-in-C2 residual" -> 54.85714285714285714285714285714285714663`24.574917376188022, "C4-quadratic-in-C2 residual" -> 54.85714285714285714285714285714285714459`24.28308003577398}
```

## 2026-07-20 05:19:48
**Input:**
```wolfram
Module[{λ1={2,2,0,0},λ2={2,2,0,0},es,rayleigh,rows,resid,l3resid,l2u0PolyC2,l2u0Claimed},
  rayleigh[m_,v_]:=(v.(m.v))/(v.v);
  es=TauEigensystemCasimir[λ1,λ2];
  rows=DeleteDuplicatesBy[
    Function[kv,Module[{key=kv[[1]],rec=kv[[2]],v},v=rec["vec"];
      {Re[rayleigh[C2op[λ1,λ2],v]],Re[rayleigh[C4op[λ1,λ2],v]],
       Re[CoefficientList[LambdaU[λ1,λ2][key[[1]],key[[2]],key[[3]],key[[4]],key[[5]]][2][uu],uu]],
       Re[CoefficientList[LambdaU[λ1,λ2][key[[1]],key[[2]],key[[3]],key[[4]],key[[5]]][3][uu],uu]]}]]/@Normal[es],
    Round[{#[[1]],#[[2]]},10^-5]&];
  resid[data_,basis_]:=Module[{f=Fit[N[data,30],basis,x]},Max[Abs[(f/.x->#[[1]])-#[[2]]]&/@data]];
  (* tau3: each of 7 coeffs, affine-in-C2 residual *)
  l3resid=Table[resid[{#[[1]],#[[4,k]]}&/@rows,{1,x}],{k,1,7}];
  (* tau2 u^0 = l2[[3,1]]: fit as polynomial in C2 alone up to cubic (should FAIL) *)
  l2u0PolyC2=resid[{#[[1]],#[[3,1]]}&/@rows,{1,x,x^2,x^3}];
  (* tau2 u^0 minus claimed C4/C2^2 pieces should be affine in C2 (should PASS ~0) *)
  l2u0Claimed=resid[{#[[1]], #[[3,1]] - (1/8) #[[1]]^2 + (1/4) #[[2]]}&/@rows,{1,x}];
  {"tau3 per-coeff affine-in-C2 residuals (u^0..u^6)"->l3resid,
   "tau2 u^0 as poly-in-C2-alone (cubic) residual [expect LARGE]"->l2u0PolyC2,
   "tau2 u^0 minus (1/8 C2^2 - 1/4 C4), affine-in-C2 residual [expect ~0]"->l2u0Claimed}]
```
**Output:**
```
Out= {"tau3 per-coeff affine-in-C2 residuals (u^0..u^6)" -> {0``25.592407455529454, 0``24.546737809246217, 0``24.07936026115753, 0``24.13947938566766, 0``24.531550530416734, 0``25.243243631453964, 0``26.371273445039403}, "tau2 u^0 as poly-in-C2-alone (cubic) residual [expect LARGE]" -> 8.06722689075630252100840336134453782048`22.70309230134672, "tau2 u^0 minus (1/8 C2^2 - 1/4 C4), affine-in-C2 residual [expect ~0]" -> 0``23.17894318807419}
```

## 2026-07-20 06:03:57
**Input:**
```wolfram
Module[{λ1={1,1,0,0},λ2={1,1,0,0},d,M0,cand,rem,scalarpart},
  d=dim[λ1] dim[λ2];
  M0=Simplify[tg[λ1,λ2][3,1][u]/trivfactor[1,1]];
  cand=(h^2/2)(C2op[λ1,λ2]-(Csite[λ1]+Csite[λ2]) IdentityMatrix[d]);
  rem=M0-cand;
  scalarpart=rem[[1,1]];
  {"remainder-minus-its-(1,1)-entry is zero matrix?"->(Simplify[rem-scalarpart IdentityMatrix[d]]==ConstantArray[0,{d,d}]),
   "scalar part (=tau3/trivfactor with Casimir removed)"->Collect[scalarpart,u]//Simplify}]
```
**Output:**
```
Out= {"remainder-minus-its-(1,1)-entry is zero matrix?" -> {{0, -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}])}, {-1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -h^2, -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1/2*(h^2*(-2 + C2op[{1, 1, 0, 0}, {1, 1, 0, 0}])), -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}])}, {-1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -h^2, -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1/2*(h^2*(-2 + C2op[{1, 1, 0, 0}, {1, 1, 0, 0}])), -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1/2*(h^2*C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]), -1
[output truncated — 61521 chars total]
```

## 2026-07-20 06:04:38
**Input:**
```wolfram
{MatrixQ[C2op[{1,1,0,0},{1,1,0,0}]], MatrixQ[C3op[{1,1,0,0},{1,1,0,0}]], MatrixQ[C4op[{1,1,0,0},{1,1,0,0}]], Head[Csite[{1,1,0,0}]], MatrixQ[tg[{1,1,0,0},{1,1,0,0}][3,1][u]], h, θ1, θ2, z[1]}
```
**Output:**
```
Out= {False, False, False, Csite, True, h, θ1, θ2, 1}
```

---

## 2026-07-20 06:07:06 — 🔄 KERNEL RESTART

## 2026-07-20 06:07:20 — ▶️ RUN CELLS 1–17
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

## 2026-07-20 06:07:31 — ▶️ RUN CELLS 35–38
Cell 35: ✓ — Out[53]= \left\{4,6,4,1\right\}
Cell 36: ✓ — (no output)
Cell 37: ✓ — (no output)
Cell 38: ✓ — (no output)

## 2026-07-20 06:07:38 — ▶️ RUN CELL 49 (0.01 s)
**In [49]:** `(*Total cubic/quartic Casimirs on the two-site space, built via block products for speed:↵  M[i,j]=EE[i,j]; A[i,k]=Sum_j M[i,j].M[j,k]; C3=Sum_ik A[i,k].M[k,i];↵  B[i,l]=Sum_k A[i,k].M[k,l]; C4=Sum_il`
**Out:** Out[65]= \left\{\left\{36,36\right\},\left\{36,36\right\},\text{C3site11}\to 18,\text{C4site11}\to 54,\text{C2site11}\to \mathrm{Csite}\left(\left\{1,1,0,0\right\}\right)\right\}

## 2026-07-20 06:07:49 — ▶️ RUN CELL 59 (0.01 s)
**In [59]:** `trivfactor[S_,P_]:=(u-θ1)(u-θ2)(u-θ1-h (S+2))(u-θ2-h (P+2));`
**Out:** (no output)


## 2026-07-20 06:07:59
**Input:**
```wolfram
Csite[λ_]:=Csite[λ]=Module[{m=Sum[Ee[λ][a,b].Ee[λ][b,a],{a,4},{b,4}]},If[m=!=m[[1,1]] IdentityMatrix[Length[m]],Print["not scalar"];Abort[]];m[[1,1]]];
{Csite[{1,1,0,0}],Csite[{2,2,0,0}],MatrixQ[C2op[{1,1,0,0},{1,1,0,0}]],MatrixQ[C4op[{1,1,0,0},{1,1,0,0}]],h,θ1,θ2,z[1]}
```
**Output:**
```
Out= {6, 16, True, True, 1, 1/3, 1/7, 1}
```

## 2026-07-20 06:08:11
**Input:**
```wolfram
Module[{λ1={1,1,0,0},λ2={1,1,0,0},d,M0,cand,rem,sp},
  d=dim[λ1]dim[λ2];
  M0=tg[λ1,λ2][3,1][u]/trivfactor[1,1];
  cand=(h^2/2)(C2op[λ1,λ2]-(Csite[λ1]+Csite[λ2])IdentityMatrix[d]);
  rem=Expand[M0-cand];
  sp=rem[[1,1]];
  {"h^2-term == (h^2/2)(C2op - Csite1 - Csite2) ?"->(Simplify[rem-sp IdentityMatrix[d]]==ConstantArray[0,{d,d}]),
   "leftover scalar sp(u)"->Collect[sp,u]}]
```
**Output:**
```
Out= {"h^2-term == (h^2/2)(C2op - Csite1 - Csite2) ?" -> True, "leftover scalar sp(u)" -> -2 + 60280/(9261*(-10/3 + u)*(-22/7 + u)*(-1/3 + u)*(-1/7 + u)) - (78256*u)/(1029*(-10/3 + u)*(-22/7 + u)*(-1/3 + u)*(-1/7 + u)) + (781622*u^2)/(3087*(-10/3 + u)*(-22/7 + u)*(-1/3 + u)*(-1/7 + u)) - (2788600*u^3)/(9261*(-10/3 + u)*(-22/7 + u)*(-1/3 + u)*(-1/7 + u)) + (24130*u^4)/(147*(-10/3 + u)*(-22/7 + u)*(-1/3 + u)*(-1/7 + u)) - (292*u^5)/(7*(-10/3 + u)*(-22/7 + u)*(-1/3 + u)*(-1/7 + u)) + (4*u^6)/((-10/3 + u)*(-22/7 + u)*(-1/3 + u)*(-1/7 + u))}
```

## 2026-07-20 06:09:11
**Input:**
```wolfram
check[λ1_,λ2_,S_,P_]:=Module[{d=dim[λ1]dim[λ2],M0,cand,rem,sp},
  M0=tg[λ1,λ2][3,1][u]/trivfactor[S,P];
  cand=(h^2/2)(C2op[λ1,λ2]-(Csite[λ1]+Csite[λ2])IdentityMatrix[d]);
  rem=M0-cand;
  sp=Simplify[rem[[1,1]]];
  {"(S,P)"->{S,P},"scalar==Id?"->(Simplify[rem-sp IdentityMatrix[d]]==ConstantArray[0,{d,d}]),
   "s3(u)"->Collect[sp,u,Simplify]}];
{check[{1,1,0,0},{1,1,0,0},1,1], check[{2,2,0,0},{1,1,0,0},2,1]}
```
**Output:**
```
Out= {{"(S,P)" -> {1, 1}, "scalar==Id?" -> True, "s3(u)" -> 232/21 - (292*u)/21 + 4*u^2}, {"(S,P)" -> {2, 1}, "scalar==Id?" -> True, "s3(u)" -> 40/3 - (334*u)/21 + 4*u^2}}
```

## 2026-07-20 06:10:38
**Input:**
```wolfram
tbl=Table[{x,Csite[{x,x,0,0}],C3site[{x,x,0,0}],C4site[{x,x,0,0}]},{x,0,5}];
{"{x,C2,C3,C4}"->tbl,
 "C2c[x]"->InterpolatingPolynomial[{#[[1]],#[[2]]}&/@tbl,x]//Expand,
 "C3c[x]"->InterpolatingPolynomial[{#[[1]],#[[3]]}&/@tbl,x]//Expand,
 "C4c[x]"->InterpolatingPolynomial[{#[[1]],#[[4]]}&/@tbl,x]//Expand}
```
**Output:**
```
Out= {"{x,C2,C3,C4}" -> {{0, 0, 0, 0}, {1, 6, 18, 54}, {2, 16, 64, 256}, {3, 30, 150, 750}, {4, 48, 288, 1728}, {5, 70, 490, 3430}}, "C2c[x]" -> 4*x + 2*x^2, "C3c[x]" -> 8*x + 8*x^2 + 2*x^3, "C4c[x]" -> 16*x + 24*x^2 + 12*x^3 + 2*x^4}
```

## 2026-07-20 06:11:17
**Input:**
```wolfram
s3fun[S_,P_,u_]:=4u^2-(4(θ1+θ2)+2h(S+P)+8h)u+(4θ1 θ2+4h(θ1+θ2)+4h^2+2h^2(S+P)+2h θ1 P+2h θ2 S);
tau3pred[λ1_,λ2_,S_,P_]:=trivfactor[S,P](s3fun[S,P,u] IdentityMatrix[dim[λ1]dim[λ2]]+(h^2/2)(C2op[λ1,λ2]-(Csite[λ1]+Csite[λ2])IdentityMatrix[dim[λ1]dim[λ2]]));
chk[λ1_,λ2_,S_,P_]:={S,P}->(Simplify[tg[λ1,λ2][3,1][u]-tau3pred[λ1,λ2,S,P]]==ConstantArray[0,{dim[λ1]dim[λ2],dim[λ1]dim[λ2]}]);
{chk[{1,1,0,0},{1,1,0,0},1,1],chk[{2,2,0,0},{1,1,0,0},2,1],chk[{1,1,0,0},{2,2,0,0},1,2]}
```
**Output:**
```
Out= {{1, 1} -> True, {2, 1} -> True, {1, 2} -> True}
```
## 2026-07-20 06:12:06 — 📥 BULK INSERT 1 CELL(S) at positions 63–63
1. [code] (*Closed-form single-site Casimir eigenvalues for the rectangular gl(4) rep [x,x,0,0].↵  Normalizati


## 2026-07-20 06:12:07
**Input:**
```wolfram
(*Closed-form single-site Casimir eigenvalues for the rectangular gl(4) rep [x,x,0,0].
  Normalization matches the notebook operators: C2op=Sum_ij E_ij E_ji, C3op=Sum_ijk E_ij E_jk E_ki,
  C4op=Sum_ijkl E_ij E_jk E_kl E_li. Obtained by InterpolatingPolynomial over x=0..5 against
  Csite/C3site/C4site; clean geometric pattern  Cn(x) = 2 x (x+2)^(n-1).*)
C2c[x_] := 2 x (x + 2);
C3c[x_] := 2 x (x + 2)^2;
C4c[x_] := 2 x (x + 2)^3;
(*The two chain sites carry x=S (site 1, rep [S,S,0,0]) and x=P (site 2, rep [P,P,0,0]). The
  combination entering the transfer-matrix eigenvalues is the per-site SUM
  C?[S,P] := C?c[S] + C?c[P]  (= Csite[{S,S,0,0}] + Csite[{P,P,0,0}] for the quadratic one).*)
C2[S_, P_] := C2c[S] + C2c[P];
C3[S_, P_] := C3c[S] + C3c[P];
C4[S_, P_] := C4c[S] + C4c[P];
(*verification against the operator eigenvalues*)
Column[{
  "single-site: Cn(x)=2x(x+2)^(n-1)  ->  C2c=2x(x+2), C3c=2x(x+2)^2, C4c=2x(x+2)^3",
  Table[{x, C2c[x] == Csite[{x, x, 0, 0}], C3c[x] == C3site[{x, x, 0, 0}], C4c[x] == C4site[{x, x, 0, 0}]}, {x, 1, 5}]}]
```
**Output:**
```
Out[73]= \begin{matrix}\text{single-site: Cn(x)=2x(x+2)\^{}(n-1)  ->  C2c=2x(x+2), C3c=2x(x+2)\^{}2, C4c=2x(x+2)\^{}3}\\\begin{pmatrix}1 & \mathrm{True} & \mathrm{True} & \mathrm{True}\\2 & \mathrm{True} & \mathrm{True} & \mathrm{True}\\3 & \mathrm{True} & \mathrm{True} & \mathrm{True}\\4 & \mathrm{True} & \mathrm{True} & \mathrm{True}\\5 & \mathrm{True} & \mathrm{True} & \mathrm{True}\end{pmatrix… [413 chars total]
```
## 2026-07-20 06:12:23 — 📥 BULK INSERT 1 CELL(S) at positions 64–64
1. [code] (*EXACT closed form for tau3 = tg[3,1] (answers "what is the h^2 term after the trivial zeros?").↵


## 2026-07-20 06:12:24
**Input:**
```wolfram
(*EXACT closed form for tau3 = tg[3,1] (answers "what is the h^2 term after the trivial zeros?").
    tau3(u) = trivfactor[S,P] ( s3scalar[S,P](u) Id  +  (h^2/2)(C2op - Csite[lam1] - Csite[lam2]) ).
  => the h^2 term is EXACTLY the same quadratic-Casimir combination (h^2/2)(C2mu - C2(lam1) - C2(lam2))
     that appears in tau1 -- the tau1<->tau3 self-duality (Lambda^3 C^4 ~ (Lambda^1 C^4)* x det), so
     tau3 carries NO higher Casimir, only C2 (matching the earlier fit: tau3 affine in C2mu everywhere).
  trivfactor[S,P] = (u-th1)(u-th2)(u-th1-h(S+2))(u-th2-h(P+2)) holds the 4 trivial zeros; s3scalar is a
  state-independent degree-2 "pseudovacuum" polynomial. Verified as an exact operator identity for
  (S,P)=(1,1),(2,1),(1,2) (incl. the S<->P-asymmetric case, asymmetry from th1!=th2).
  NB C2op eigenvalue on the top state gives (h^2/2)(C2_top - Csite1 - Csite2) = 2 h^2 S P.*)
s3scalar[S_, P_, u_] := 4 u^2 - (4 (θ1 + θ2) + 2 h (S + P) + 8 h) u + (4 θ1 θ2 + 4 h (θ1 + θ2) + 4 h^2 + 2 h^2 (S + P) + 2 h θ1 P + 2 h θ2 S);
tau3exact[λ1_, λ2_, S_, P_] := trivfactor[S, P] (s3scalar[S, P, u] IdentityMatrix[dim[λ1] dim[λ2]] + (h^2/2) (C2op[λ1, λ2] - (Csite[λ1] + Csite[λ2]) IdentityMatrix[dim[λ1] dim[λ2]]));
{"tau3 closed form == tg[3,1] ?",
 {1, 1} -> (Simplify[tg[{1,1,0,0},{1,1,0,0}][3,1][u] - tau3exact[{1,1,0,0},{1,1,0,0},1,1]] == ConstantArray[0,{36,36}]),
 {2, 1} -> (Simplify[tg[{2,2,0,0},{1,1,0,0}][3,1][u] - tau3exact[{2,2,0,0},{1,1,0,0},2,1]] == ConstantArray[0,{120,120}]),
 {1, 2} -> (Simplify[tg[{1,1,0,0},{2,2,0,0}][3,1][u] - tau3exact[{1,1,0,0},{2,2,0,0},1,2]] == ConstantArray[0,{120,120}])}
```
**Output:**
```
Out[76]= \left\{\text{tau3 closed form == tg[3,1] ?},\left\{1,1\right\}\to \mathrm{True},\left\{2,1\right\}\to \mathrm{True},\left\{1,2\right\}\to \mathrm{True}\right\}
```

## 2026-07-20 06:16:25
**Input:**
```wolfram
{MatrixQ[C2op[{1,1,0,0},{1,1,0,0}]], MatrixQ[C3op[{1,1,0,0},{1,1,0,0}]], MatrixQ[C4op[{1,1,0,0},{1,1,0,0}]], Csite[{1,1,0,0}], C3site[{1,1,0,0}], C4site[{1,1,0,0}]}
```
**Output:**
```
Out= {True, True, True, 6, 18, 54}
```

## 2026-07-20 06:16:32
**Input:**
```wolfram
Omega2[λ1_,λ2_]:=Omega2[λ1,λ2]=Sum[KroneckerProduct[Ee[λ1][i,j],Ee[λ2][j,i]],{i,4},{j,4}];
Module[{λ1={1,1,0,0},λ2={1,1,0,0},d=dim[{1,1,0,0}]dim[{1,1,0,0}]},
  {"C2op == (Csite1+Csite2) Id + 2 Omega2 ?"->(C2op[λ1,λ2]-(Csite[λ1]+Csite[λ2])IdentityMatrix[d]-2Omega2[λ1,λ2]==ConstantArray[0,{d,d}])}]
```
**Output:**
```
Out= {"C2op == (Csite1+Csite2) Id + 2 Omega2 ?" -> True}
```

## 2026-07-20 06:16:40
**Input:**
```wolfram
(*Two natural "2-from-site1,1-from-site2" and "1-from-site1,2-from-site2" cross blocks for C3,
  built as chains contracted the same way C3site itself is (E_ij E_jk E_ki cyclic).*)
X1[λ1_,λ2_]:=X1[λ1,λ2]=Sum[KroneckerProduct[Ee[λ1][i,j].Ee[λ1][j,k],Ee[λ2][k,i]],{i,4},{j,4},{k,4}];
X2[λ1_,λ2_]:=X2[λ1,λ2]=Sum[KroneckerProduct[Ee[λ1][i,j],Ee[λ2][j,k].Ee[λ2][k,i]],{i,4},{j,4},{k,4}];
Module[{λ1={1,1,0,0},λ2={1,1,0,0},d,pure,rem,c},
  d=dim[λ1]dim[λ2];
  pure=(C3site[λ1]+C3site[λ2])IdentityMatrix[d];
  rem=C3op[λ1,λ2]-pure;
  (*solve for c such that rem == c(X1+X2), by comparing a single generic nonzero entry*)
  c=Simplify[rem[[2,7]]/(X1[λ1,λ2]+X2[λ1,λ2])[[2,7]]];
  {"candidate coeff c"->c,
   "rem == c(X1+X2) exactly?"->(rem-c(X1[λ1,λ2]+X2[λ1,λ2])==ConstantArray[0,{d,d}]),
   "X1==X2?"->(X1[λ1,λ2]==X2[λ1,λ2])}]
```
**Output:**
```
Out= {"candidate coeff c" -> 5/3, "rem == c(X1+X2) exactly?" -> False, "X1==X2?" -> True}
```

## 2026-07-20 06:17:04
**Input:**
```wolfram
Module[{λ1={1,1,0,0},λ2={1,1,0,0},d,EE1,EE2,patterns,terms},
  d=dim[λ1]dim[λ2];
  EE1[i_,j_]:=KroneckerProduct[Ee[λ1][i,j],id[λ2]];
  EE2[i_,j_]:=KroneckerProduct[id[λ1],Ee[λ2][i,j]];
  patterns=Tuples[{1,2},3]; (*which site each of the 3 factors comes from*)
  terms=Association[Table[
    pat->Sum[
      Which[pat[[1]]==1,EE1[i,j],True,EE2[i,j]].Which[pat[[2]]==1,EE1[j,k],True,EE2[j,k]].Which[pat[[3]]==1,EE1[k,i],True,EE2[k,i]],
      {i,4},{j,4},{k,4}],
    {pat,patterns}]];
  Association[KeyValueMap[#1->(#2===ConstantArray[0,{d,d}])->"zero", terms]&, terms]; (*placeholder, real check below*)
  {"patterns"->patterns,
   "term[1,1,1]==C3site1 Id?"->(terms[{1,1,1}]==C3site[λ1] IdentityMatrix[d]),
   "term[2,2,2]==Id C3site2?"->(terms[{2,2,2}]==C3site[λ2] IdentityMatrix[d]),
   "term[1,1,2]==term[1,2,1]?"->(terms[{1,1,2}]==terms[{1,2,1}]),
   "term[1,1,2]==term[2,1,1]?"->(terms[{1,1,2}]==terms[{2,1,1}]),
   "term[2,2,1]==term[2,1,2]?"->(terms[{2,2,1}]==terms[{2,1,2}]),
   "term[2,2,1]==term[1,2,2]?"->(terms[{2,2,1}]==terms[{1,2,2}]),
   "term[1,1,2]==term[2,2,1]?"->(terms[{1,1,2}]==terms[{2,2,1}]),
   "term[1,1,2]==X1?"->(terms[{1,1,2}]==X1[λ1,λ2]),
   "term[2,2,1]==X2?"->(terms[{2,2,1}]==X2[λ1,λ2])}]
```
**Output:**
```
Out= {"patterns" -> {{1, 1, 1}, {1, 1, 2}, {1, 2, 1}, {1, 2, 2}, {2, 1, 1}, {2, 1, 2}, {2, 2, 1}, {2, 2, 2}}, "term[1,1,1]==C3site1 Id?" -> True, "term[2,2,2]==Id C3site2?" -> True, "term[1,1,2]==term[1,2,1]?" -> False, "term[1,1,2]==term[2,1,1]?" -> True, "term[2,2,1]==term[2,1,2]?" -> False, "term[2,2,1]==term[1,2,2]?" -> True, "term[1,1,2]==term[2,2,1]?" -> True, "term[1,1,2]==X1?" -> True, "term[2,2,1]==X2?" -> True}
```

## 2026-07-20 06:17:47
**Input:**
```wolfram
Module[{λ1={1,1,0,0},λ2={1,1,0,0},d,EE1,EE2,patterns,terms,classes},
  d=dim[λ1]dim[λ2];
  EE1[i_,j_]:=KroneckerProduct[Ee[λ1][i,j],id[λ2]];
  EE2[i_,j_]:=KroneckerProduct[id[λ1],Ee[λ2][i,j]];
  patterns=Tuples[{1,2},3];
  terms=Association[Table[
    pat->Sum[
      Which[pat[[1]]==1,EE1[i,j],True,EE2[i,j]].Which[pat[[2]]==1,EE1[j,k],True,EE2[j,k]].Which[pat[[3]]==1,EE1[k,i],True,EE2[k,i]],
      {i,4},{j,4},{k,4}],
    {pat,patterns}]];
  (*group patterns by which operator they equal*)
  classes=GatherBy[patterns,terms[#]&];
  {"equivalence classes (by pattern list)"->(First/@#&/@classes),
   "class sizes"->Length/@classes,
   "sum of ALL 8 == C3op directly?"->(Total[Values[terms]]==C3op[λ1,λ2])}]
```
**Output:**
```
Out= {"equivalence classes (by pattern list)" -> {{1, 2}, {1, 1, 2, 2}, {1, 2}}, "class sizes" -> {2, 4, 2}, "sum of ALL 8 == C3op directly?" -> True}
```

## 2026-07-20 06:18:16
**Input:**
```wolfram
Module[{λ1={1,1,0,0},λ2={2,2,0,0},d,EE1,EE2,patterns,terms,classes},
  d=dim[λ1]dim[λ2];
  EE1[i_,j_]:=KroneckerProduct[Ee[λ1][i,j],id[λ2]];
  EE2[i_,j_]:=KroneckerProduct[id[λ1],Ee[λ2][i,j]];
  patterns=Tuples[{1,2},3];
  terms=Association[Table[
    pat->Sum[
      Which[pat[[1]]==1,EE1[i,j],True,EE2[i,j]].Which[pat[[2]]==1,EE1[j,k],True,EE2[j,k]].Which[pat[[3]]==1,EE1[k,i],True,EE2[k,i]],
      {i,4},{j,4},{k,4}],
    {pat,patterns}]];
  classes=GatherBy[patterns,terms[#]&];
  {"classes(patterns grouped)"->classes,"sizes"->Length/@classes,
   "sum of ALL 8 == C3op?"->(Total[Values[terms]]==C3op[λ1,λ2])}]
```
**Output:**
```
Out= {"classes(patterns grouped)" -> {{{1, 1, 1}}, {{1, 1, 2}, {2, 1, 1}}, {{1, 2, 1}}, {{1, 2, 2}, {2, 2, 1}}, {{2, 1, 2}}, {{2, 2, 2}}}, "sizes" -> {1, 2, 1, 2, 1, 1}, "sum of ALL 8 == C3op?" -> True}
```

## 2026-07-20 06:18:54
**Input:**
```wolfram
X1[λ1_,λ2_]:=X1[λ1,λ2]=Sum[KroneckerProduct[Ee[λ1][i,j].Ee[λ1][j,k],Ee[λ2][k,i]],{i,4},{j,4},{k,4}];
X2[λ1_,λ2_]:=X2[λ1,λ2]=Sum[KroneckerProduct[Ee[λ1][i,j],Ee[λ2][j,k].Ee[λ2][k,i]],{i,4},{j,4},{k,4}];
Y1[λ1_,λ2_]:=Y1[λ1,λ2]=Sum[KroneckerProduct[Ee[λ1][i,j].Ee[λ1][k,i],Ee[λ2][j,k]],{i,4},{j,4},{k,4}];
Y2[λ1_,λ2_]:=Y2[λ1,λ2]=Sum[KroneckerProduct[Ee[λ1][j,k],Ee[λ2][i,j].Ee[λ2][k,i]],{i,4},{j,4},{k,4}];
Module[{λ1={1,1,0,0},λ2={2,2,0,0},d,pred},
  d=dim[λ1]dim[λ2];
  pred=(C3site[λ1]+C3site[λ2])IdentityMatrix[d]+2X1[λ1,λ2]+2X2[λ1,λ2]+Y1[λ1,λ2]+Y2[λ1,λ2];
  {"C3op == C3site1+C3site2 + 2X1+2X2+Y1+Y2 ?"->(pred==C3op[λ1,λ2])}]
```
**Output:**
```
Out= {"C3op == C3site1+C3site2 + 2X1+2X2+Y1+Y2 ?" -> True}
```

## 2026-07-20 06:19:02
**Input:**
```wolfram
Module[{λ1={2,1,1,0},λ2={1,0,0,0},d,pred},
  d=dim[λ1]dim[λ2];
  pred=(C3site[λ1]+C3site[λ2])IdentityMatrix[d]+2X1[λ1,λ2]+2X2[λ1,λ2]+Y1[λ1,λ2]+Y2[λ1,λ2];
  {"C3op == ... general identity check (adjoint x fund)"->(pred==C3op[λ1,λ2])}]
```
**Output:**
```
Out= {"C3op == ... general identity check (adjoint x fund)" -> True}
```

## 2026-07-20 06:21:13
**Input:**
```wolfram
Module[{λ1={1,1,0,0},λ2={2,2,0,0},d,EE1,EE2,patterns,terms,classes},
  d=dim[λ1]dim[λ2];
  EE1[i_,j_]:=KroneckerProduct[Ee[λ1][i,j],id[λ2]];
  EE2[i_,j_]:=KroneckerProduct[id[λ1],Ee[λ2][i,j]];
  patterns=Tuples[{1,2},4];
  terms=Association[Table[
    pat->Sum[
      Which[pat[[1]]==1,EE1[i,j],True,EE2[i,j]].Which[pat[[2]]==1,EE1[j,k],True,EE2[j,k]].
      Which[pat[[3]]==1,EE1[k,l],True,EE2[k,l]].Which[pat[[4]]==1,EE1[l,i],True,EE2[l,i]],
      {i,4},{j,4},{k,4},{l,4}],
    {pat,patterns}]];
  classes=GatherBy[patterns,terms[#]&];
  {"classes"->classes,"sizes"->Length/@classes,
   "sum of ALL 16 == C4op?"->(Total[Values[terms]]==C4op[λ1,λ2])}]
```
**Output:**
```
Out= {"classes" -> {{{1, 1, 1, 1}}, {{1, 1, 1, 2}, {2, 1, 1, 1}}, {{1, 1, 2, 1}, {1, 2, 1, 1}}, {{1, 1, 2, 2}, {2, 2, 1, 1}}, {{1, 2, 1, 2}, {2, 1, 2, 1}}, {{1, 2, 2, 1}}, {{1, 2, 2, 2}, {2, 2, 2, 1}}, {{2, 1, 1, 2}}, {{2, 1, 2, 2}, {2, 2, 1, 2}}, {{2, 2, 2, 2}}}, "sizes" -> {1, 2, 2, 2, 2, 1, 2, 1, 2, 1}, "sum of ALL 16 == C4op?" -> True}
```

## 2026-07-20 06:22:33
**Input:**
```wolfram
(*General rule (derived+verified): for pattern (a1..an), the chain term collapses to
   [ordered product of E1 factors at site-1 positions] (x) [ordered product of E2 factors at site-2 positions],
  with cyclic indices i_(n+1)=i1, summed freely over all i1..in. Build this directly for any pattern.*)
patternTerm[pat_List,λ1_,λ2_]:=Module[{n=Length[pat],idx,pos1,pos2,f1,f2},
  idx=Array[ii,n];
  pos1=Position[pat,1]//Flatten; pos2=Position[pat,2]//Flatten;
  f1=If[pos1=={},id[λ1],Dot@@(Ee[λ1]@@@Table[{idx[[m]],idx[[Mod[m,n]+1]]},{m,pos1}])];
  f2=If[pos2=={},id[λ2],Dot@@(Ee[λ2]@@@Table[{idx[[m]],idx[[Mod[m,n]+1]]},{m,pos2}])];
  Sum[KroneckerProduct[f1,f2],Evaluate[Sequence@@Table[{idx[[m]],4},{m,n}]]]];
(*cross-check against the direct pattern computation for a few C3/C4 patterns*)
Module[{λ1={1,1,0,0},λ2={2,2,0,0},EE1,EE2,direct311},
  EE1[i_,j_]:=KroneckerProduct[Ee[λ1][i,j],id[λ2]];
  EE2[i_,j_]:=KroneckerProduct[id[λ1],Ee[λ2][i,j]];
  direct311=Sum[EE1[i,j].EE2[j,k].EE1[k,i],{i,4},{j,4},{k,4}];
  {"patternTerm[{1,2,1}] matches direct computation?"->(patternTerm[{1,2,1},λ1,λ2]==direct311),
   "patternTerm[{1,1,1}]==C3site1 Id?"->(patternTerm[{1,1,1},λ1,λ2]==C3site[λ1]IdentityMatrix[dim[λ1]dim[λ2]])}]
```
**Output:**
```
Out= {"patternTerm[{1,2,1}] matches direct computation?" -> True, "patternTerm[{1,1,1}]==C3site1 Id?" -> True}
```

## 2026-07-20 06:23:04
**Input:**
```wolfram
Module[{λ1={1,1,0,0},λ2={2,2,0,0},patterns,terms,classes,reps,mults,assembled},
  patterns=Tuples[{1,2},4];
  terms=Association[Table[pat->patternTerm[pat,λ1,λ2],{pat,patterns}]];
  classes=GatherBy[patterns,terms[#]&];
  reps=First/@classes; mults=Length/@classes;
  assembled=Total[MapThread[#2 patternTerm[#1,λ1,λ2]&,{reps,mults}]];
  {"representative patterns"->reps,"multiplicities"->mults,
   "assembled == C4op?"->(assembled==C4op[λ1,λ2])}]
```
**Output:**
```
Out= {"representative patterns" -> {{1, 1, 1, 1}, {1, 1, 1, 2}, {1, 1, 2, 1}, {1, 1, 2, 2}, {1, 2, 1, 2}, {1, 2, 2, 1}, {1, 2, 2, 2}, {2, 1, 1, 2}, {2, 1, 2, 2}, {2, 2, 2, 2}}, "multiplicities" -> {1, 2, 2, 2, 2, 1, 2, 1, 2, 1}, "assembled == C4op?" -> True}
```
## 2026-07-20 06:23:50 — 📥 BULK INSERT 1 CELL(S) at positions 65–65
1. [code] (*GENERAL master formula for ALL of C2op,C3op,C4op (any n): the operator-valued "split Casimir"↵  de


## 2026-07-20 06:24:06
**Input:**
```wolfram
(*GENERAL master formula for ALL of C2op,C3op,C4op (any n): the operator-valued "split Casimir"
  decomposition of Cn_op = Sum_{i1..in} EE[i1,i2].EE[i2,i3]...EE[in,i1] (EE=E1(x)Id+Id(x)E2).
  Substituting EE=E1+E2 into the n-fold cyclic product and using (A(x)B).(C(x)D)=(A.C)(x)(B.D)
  repeatedly, EVERY one of the 2^n site-assignment terms collapses to
     [ordered matrix product of the site-1 factors, in their original chain positions]
     (x) [ordered matrix product of the site-2 factors, in their original chain positions],
  with the SAME cyclic index set i1..in tying both chains together (no further index contraction
  between the two factors). patternTerm[pat,λ1,λ2] builds exactly this for any pattern pat in
  {1,2}^n. VERIFIED against the direct (unsimplified) chain computation and against C3op above.*)
patternTerm[pat_List, λ1_, λ2_] := Module[{n = Length[pat], idx, pos1, pos2, f1, f2},
  idx = Array[ii, n];
  pos1 = Position[pat, 1] // Flatten; pos2 = Position[pat, 2] // Flatten;
  f1 = If[pos1 == {}, id[λ1], Dot @@ (Ee[λ1] @@@ Table[{idx[[m]], idx[[Mod[m, n] + 1]]}, {m, pos1}])];
  f2 = If[pos2 == {}, id[λ2], Dot @@ (Ee[λ2] @@@ Table[{idx[[m]], idx[[Mod[m, n] + 1]]}, {m, pos2}])];
  Sum[KroneckerProduct[f1, f2], Evaluate[Sequence @@ Table[{idx[[m]], 4}, {m, n}]]]];
(*MASTER FORMULA (any n):   Cn_op[λ1,λ2] = Sum_{pat in {1,2}^n} patternTerm[pat,λ1,λ2]   *)
CnMaster[n_, λ1_, λ2_] := Total[patternTerm[#, λ1, λ2] & /@ Tuples[{1, 2}, n]];
{"CnMaster[2,..]==C2op?" -> (CnMaster[2, {1,1,0,0}, {2,2,0,0}] == C2op[{1,1,0,0}, {2,2,0,0}]),
 "CnMaster[3,..]==C3op?" -> (CnMaster[3, {1,1,0,0}, {2,2,0,0}] == C3op[{1,1,0,0}, {2,2,0,0}]),
 "CnMaster[4,..]==C4op?" -> (CnMaster[4, {1,1,0,0}, {2,2,0,0}] == C4op[{1,1,0,0}, {2,2,0,0}])}
```
**Output:**
```
Out[79]= \left\{\text{CnMaster[2,..]==C2op?}\to \mathrm{True},\text{CnMaster[3,..]==C3op?}\to \mathrm{True},\text{CnMaster[4,..]==C4op?}\to \mathrm{True}\right\}
```
## 2026-07-20 06:24:24 — 📥 BULK INSERT 1 CELL(S) at positions 66–66
1. [code] (*Reduced explicit "named cross-term" forms, grouping the 2^n patterns by which OPERATOR they↵  prod


## 2026-07-20 06:24:32
**Input:**
```wolfram
(*Reduced explicit "named cross-term" forms, grouping the 2^n patterns by which OPERATOR they
  produce (found by GatherBy on an asymmetric test pair, so no accidental C3site1=C3site2-type
  degeneracy hides real distinctions; growth in cross-term COUNT is itself notable: n=2 has 1
  distinct mixed type, n=3 has 4, n=4 has 8.)

  C2op[λ1,λ2] = (Csite[λ1]+Csite[λ2]) Id  +  2 Omega2[λ1,λ2]
    Omega2[λ1,λ2] := Sum_ij Ee[λ1][i,j] (x) Ee[λ2][j,i]

  C3op[λ1,λ2] = (C3site[λ1]+C3site[λ2]) Id  +  2 X1  +  2 X2  +  Y1  +  Y2
    X1[λ1,λ2] := Sum_ijk (Ee[λ1][i,j].Ee[λ1][j,k]) (x) Ee[λ2][k,i]      {two E1 chained, then E2}
    X2[λ1,λ2] := Sum_ijk Ee[λ1][i,j] (x) (Ee[λ2][j,k].Ee[λ2][k,i])      {mirror of X1}
    Y1[λ1,λ2] := Sum_ijk (Ee[λ1][i,j].Ee[λ1][k,i]) (x) Ee[λ2][j,k]      {E1,E2,E1 interleaved}
    Y2[λ1,λ2] := Sum_ijk Ee[λ1][j,k] (x) (Ee[λ2][i,j].Ee[λ2][k,i])      {E2,E1,E2 interleaved}

  C4op[λ1,λ2] = (C4site[λ1]+C4site[λ2]) Id + 2 Z1 + 2 Z2 + 2 Z3 + 2 Z4 + Z5 + 2 Z6 + Z7 + 2 Z8
    (all built here via patternTerm on the representative patterns found by GatherBy, guaranteed
     correct by direct verification rather than by-hand index bookkeeping -- 4-fold products have
     enough interleavings that hand-deriving each risks a transcription error.)
    Z1 = patternTerm[{1,1,1,2}]   Z2 = patternTerm[{1,1,2,1}]   Z3 = patternTerm[{1,1,2,2}]
    Z4 = patternTerm[{1,2,1,2}]   Z5 = patternTerm[{1,2,2,1}]   Z6 = patternTerm[{1,2,2,2}]
    Z7 = patternTerm[{2,1,1,2}]   Z8 = patternTerm[{2,1,2,2}]                                    *)
Z1[λ1_, λ2_] := Z1[λ1, λ2] = patternTerm[{1, 1, 1, 2}, λ1, λ2];
Z2[λ1_, λ2_] := Z2[λ1, λ2] = patternTerm[{1, 1, 2, 1}, λ1, λ2];
Z3[λ1_, λ2_] := Z3[λ1, λ2] = patternTerm[{1, 1, 2, 2}, λ1, λ2];
Z4[λ1_, λ2_] := Z4[λ1, λ2] = patternTerm[{1, 2, 1, 2}, λ1, λ2];
Z5[λ1_, λ2_] := Z5[λ1, λ2] = patternTerm[{1, 2, 2, 1}, λ1, λ2];
Z6[λ1_, λ2_] := Z6[λ1, λ2] = patternTerm[{1, 2, 2, 2}, λ1, λ2];
Z7[λ1_, λ2_] := Z7[λ1, λ2] = patternTerm[{2, 1, 1, 2}, λ1, λ2];
Z8[λ1_, λ2_] := Z8[λ1, λ2] = patternTerm[{2, 1, 2, 2}, λ1, λ2];
Module[{λ1 = {1, 1, 0, 0}, λ2 = {2, 2, 0, 0}, d, C2pred, C3pred, C4pred},
  d = dim[λ1] dim[λ2];
  C2pred = (Csite[λ1] + Csite[λ2]) IdentityMatrix[d] + 2 Omega2[λ1, λ2];
  C3pred = (C3site[λ1] + C3site[λ2]) IdentityMatrix[d] + 2 X1[λ1, λ2] + 2 X2[λ1, λ2] + Y1[λ1, λ2] + Y2[λ1, λ2];
  C4pred = (C4site[λ1] + C4site[λ2]) IdentityMatrix[d] + 2 Z1[λ1, λ2] + 2 Z2[λ1, λ2] + 2 Z3[λ1, λ2] +
     2 Z4[λ1, λ2] + Z5[λ1, λ2] + 2 Z6[λ1, λ2] + Z7[λ1, λ2] + 2 Z8[λ1, λ2];
  {"C2op reduced form matches" -> (C2pred == C2op[λ1, λ2]),
   "C3op reduced form matches" -> (C3pred == C3op[λ1, λ2]),
   "C4op reduced form matches" -> (C4pred == C4op[λ1, λ2])}]
```
**Output:**
```
Out[88]= \left\{\text{C2op reduced form matches}\to \mathrm{True},\text{C3op reduced form matches}\to \mathrm{True},\text{C4op reduced form matches}\to \mathrm{True}\right\}
```
## 2026-07-20 08:21:54 — 📥 BULK INSERT 1 CELL(S) at positions 67–67
1. [code] (* ============================================================================================↵   E


## 2026-07-20 08:21:54
**Input:**
```wolfram
(* ============================================================================================
   EXPLICIT CLOSED-FORM CASIMIR FORMULAS, two-site total (state lam1=[S,S,0,0], lam2=[P,P,0,0]).
   These are the C2[S,P], C3[S,P], C4[S,P] the tau1/tau2/tau3 eigenvalue formulas above are built
   from (Csite[lam1]+Csite[lam2] etc.). Defined here in ONE self-contained cell, no other cell
   needs to be evaluated first to read these off.
   ============================================================================================ *)

(*  C2[S,P]  =  2 S (S+2)  +  2 P (P+2)   =   2S^2 + 2P^2 + 4S + 4P            *)
(*  C3[S,P]  =  2 S (S+2)^2  +  2 P (P+2)^2                                    *)
(*  C4[S,P]  =  2 S (S+2)^3  +  2 P (P+2)^3                                    *)

C2[S_, P_] := 2 S (S + 2) + 2 P (P + 2);
C3[S_, P_] := 2 S (S + 2)^2 + 2 P (P + 2)^2;
C4[S_, P_] := 2 S (S + 2)^3 + 2 P (P + 2)^3;

TableForm[
 {{"C2[S,P]", C2[S, P] // Expand}, {"C3[S,P]", C3[S, P] // Expand}, {"C4[S,P]", C4[S, P] // Expand}},
 TableHeadings -> {None, {"", "expanded"}}]
```
**Output:**
```
Out[92]//TableForm= \begin{pmatrix}\text{C2[S,P]} & 2\,P^2+4\,P+2\,S^2+4\,S\\\text{C3[S,P]} & 2\,P^3+8\,P^2+8\,P+2\,S^3+8\,S^2+8\,S\\\text{C4[S,P]} & 2\,P^4+12\,P^3+24\,P^2+16\,P+2\,S^4+12\,S^3+24\,S^2+16\,S\end{pmatrix}
```
## 2026-07-20 08:24:57 — 📥 BULK INSERT 1 CELL(S) at positions 68–68
1. [code] (* ============================================================================================↵   C


## 2026-07-20 08:24:58
**Input:**
```wolfram
(* ============================================================================================
   CASIMIR OPERATORS (not eigenvalues), two-site total. E_ab := E^(1)_ab (x) Id + Id (x) E^(2)_ab
   (a,b=1..4) is the total gl(4) generator acting on the two-site space lam1 (x) lam2 -- this is
   EE[lam1,lam2][a,b] in this notebook. These operators commute with EVERY E_ab (they are Casimir
   elements of the enveloping algebra), so by Schur's lemma each acts as a SCALAR on any
   irreducible component mu of lam1 (x) lam2 -- that scalar is C2(mu), C3(mu), C4(mu) (the
   eigenvalue formulas in the cell above, cell 67, and in cells 63/65/66). This cell is the
   operator DEFINITION itself, self-contained (only needs Ee/id/dim from Section A).

        C2  =  Sum_{a,b=1}^{4}         E_ab . E_ba
        C3  =  Sum_{a,b,c=1}^{4}       E_ab . E_bc . E_ca
        C4  =  Sum_{a,b,c,d=1}^{4}     E_ab . E_bc . E_cd . E_da
   ============================================================================================ *)

C2op[λ1_, λ2_] := C2op[λ1, λ2] = Sum[EE[λ1, λ2][a, b] . EE[λ1, λ2][b, a], {a, 4}, {b, 4}];
C3op[λ1_, λ2_] := C3op[λ1, λ2] = Sum[EE[λ1, λ2][a, b] . EE[λ1, λ2][b, c] . EE[λ1, λ2][c, a], {a, 4}, {b, 4}, {c, 4}];
C4op[λ1_, λ2_] := C4op[λ1, λ2] = Sum[EE[λ1, λ2][a, b] . EE[λ1, λ2][b, c] . EE[λ1, λ2][c, d] . EE[λ1, λ2][d, a], {a, 4}, {b, 4}, {c, 4}, {d, 4}];

(*concrete: on lam1=lam2={1,1,0,0}, C2op is literally 16*IdentityMatrix (scalar on this
  irreducible 6x6 (x) 6x6=36-dim tensor factor's mu=top-weight sector -- more generally it is
  scalar on EACH mu-block, with block value C2(mu); NOT scalar on the whole 36x36 space since
  lam1(x)lam2 decomposes into several mu's with different C2(mu).)*)
{"C2op[{1,1,0,0},{1,1,0,0}] is 36x36" -> Dimensions[C2op[{1, 1, 0, 0}, {1, 1, 0, 0}]],
 "sample entries (top-left 3x3 block)" -> C2op[{1, 1, 0, 0}, {1, 1, 0, 0}][[1 ;; 3, 1 ;; 3]]}
```
**Output:**
```
Out[96]= \left\{\text{C2op[\{1,1,0,0\},\{1,1,0,0\}] is 36x36}\to \left\{36,36\right\},\text{sample entries (top-left 3x3 block)}\to \begin{pmatrix}16 & 0 & 0\\0 & 14 & 0\\0 & 0 & 14\end{pmatrix}\right\}
```
## 2026-07-20 10:39:23 — ▶️ RUN CELLS 1–17
Cell 1: ✓ — (no output)
Cell 2: ✓ — (no output)
Cell 3: ✓ — (no output)
Cell 4: ✓ — (no output)
Cell 5: ✓ — Out[167]= \left\{0\right\}
Cell 6: ✓ — (no output)
Cell 7: ✓ — Out[179]= \left\{\frac{1}{3},\frac{1}{7},1,0.3604398751472265372666776505955500912940811467928913386732181476418264+0.9327824485933747164848359820708679294577272928322431080055758094332927\,\mathrm{i},0.5088673809904998281633213722670686199090972690844597081081963623607767+0.8608449271290790931886443936161642458151427736784766273900288346856365\,\mathrm{i},0.533257880589599372595832345093124295361… [1247 chars total]
Cell 8: ✓ — (no output)
Cell 9: ✓ — Out[183]= \left\{0.\times 10^{-69},0.0063408739330632579023526945309992501993729949274710916715470375206530\right\}
Cell 10: ✓ — (no output)
Cell 11: ✓ — Out[187]= \left\{0\right\}
Cell 12: ✓ — (no output)
Cell 13: ✓ — (no output)
Cell 14: ✓ — (no output)
Cell 15: ✓ — Out[195]= \left\{0\right\}
Cell 16: ✓ — (no output)
Cell 17: ✓ — (no output)

## 2026-07-20 10:39:29 — ▶️ RUN CELLS 35–38
Cell 35: ✓ — Out[204]= \left\{4,6,4,1\right\}
Cell 36: ✓ — (no output)
Cell 37: ✓ — (no output)
Cell 38: ✓ — (no output)

## 2026-07-20 10:39:35 — ▶️ RUN CELL 49 (0.01 s)
**In [49]:** `(*Total cubic/quartic Casimirs on the two-site space, built via block products for speed:↵  M[i,j]=EE[i,j]; A[i,k]=Sum_j M[i,j].M[j,k]; C3=Sum_ik A[i,k].M[k,i];↵  B[i,l]=Sum_k A[i,k].M[k,l]; C4=Sum_il`
**Out:** Out[216]= \left\{\left\{36,36\right\},\left\{36,36\right\},\text{C3site11}\to 18,\text{C4site11}\to 54,\text{C2site11}\to \mathrm{Csite}\left(\left\{1,1,0,0\right\}\right)\right\}


## 2026-07-20 10:39:40
**Input:**
```wolfram
Csite[λ_]:=Csite[λ]=Module[{m=Sum[Ee[λ][a,b].Ee[λ][b,a],{a,4},{b,4}]},If[m=!=m[[1,1]] IdentityMatrix[Length[m]],Print["not scalar"];Abort[]];m[[1,1]]];
{Csite[{1,1,0,0}],Csite[{2,2,0,0}],C3site[{1,1,0,0}],C4site[{1,1,0,0}],h,θ1,θ2,z[1],MatrixQ[C2op[{1,1,0,0},{1,1,0,0}]],MatrixQ[C4op[{1,1,0,0},{1,1,0,0}]]}
```
**Output:**
```
Out= {6, 16, 18, 54, 1, 1/3, 1/7, 1, True, True}
```

## 2026-07-20 10:39:46
**Input:**
```wolfram
C2[S_,P_]:=2 S(S+2)+2 P(P+2);
C3[S_,P_]:=2 S(S+2)^2+2 P(P+2)^2;
C4[S_,P_]:=2 S(S+2)^3+2 P(P+2)^3;
{C2[1,1],C3[1,1],C4[1,1]}
```
**Output:**
```
Out= {12, 36, 108}
```

## 2026-07-20 10:40:03
**Input:**
```wolfram
solveAB[M_,d_,C2m_]:=Module[{sol},
  sol=Quiet@LeastSquares[Transpose[{Flatten[IdentityMatrix[d]],Flatten[C2m]}],Flatten[M]];
  {"alpha"->sol[[1]],"beta"->sol[[2]],"exact?"->(M-sol[[1]] IdentityMatrix[d]-sol[[2]] C2m==ConstantArray[0,{d,d}])}];
Module[{λ1={1,1,0,0},λ2={1,1,0,0},d,c2,C2m,C4m,u0mat},
  d=dim[λ1]dim[λ2];
  c2=tgCoeffs[λ1,λ2][[2]]; (*{D0,D1,D2,D3,D4} coeff matrices of tau2, u^0..u^4*)
  C2m=C2op[λ1,λ2]; C4m=C4op[λ1,λ2];
  u0mat=c2[[1]]-(1/8)C2m.C2m+(1/4)C4m; (*subtract known C2^2,C4 pieces from u^0*)
  {"leading (u^4) coeff, expect 6 Id, no C2 dep"->solveAB[c2[[5]],d,C2m],
   "u^3"->solveAB[c2[[4]],d,C2m],
   "u^2"->solveAB[c2[[3]],d,C2m],
   "u^1"->solveAB[c2[[2]],d,C2m],
   "u^0 (after subtracting C2^2/8, +C4/4)"->solveAB[u0mat,d,C2m]}]
```
**Output:**
```
Out= {"leading (u^4) coeff, expect 6 Id, no C2 dep" -> {"alpha" -> 6, "beta" -> 0, "exact?" -> True}, "u^3" -> {"alpha" -> -208/7, "beta" -> 0, "exact?" -> True}, "u^2" -> {"alpha" -> 4862/147, "beta" -> 1, "exact?" -> True}, "u^1" -> {"alpha" -> 1352/147, "beta" -> -52/21, "exact?" -> True}, "u^0 (after subtracting C2^2/8, +C4/4)" -> {"alpha" -> -5176/441, "beta" -> 275/84, "exact?" -> True}}
```

## 2026-07-20 10:42:54
**Input:**
```wolfram
getCoeffs[S_,P_]:=Module[{λ1={S,S,0,0},λ2={P,P,0,0},d,c2,C2m,C4m,u0mat},
  d=dim[λ1]dim[λ2];
  c2=tgCoeffs[λ1,λ2][[2]];
  C2m=C2op[λ1,λ2]; C4m=C4op[λ1,λ2];
  u0mat=c2[[1]]-(1/8)C2m.C2m+(1/4)C4m;
  {S,P,solveAB[c2[[5]],d,C2m],solveAB[c2[[4]],d,C2m],solveAB[c2[[3]],d,C2m],solveAB[c2[[2]],d,C2m],solveAB[u0mat,d,C2m]}];
data21=getCoeffs[2,1];
data12=getCoeffs[1,2];
data22=getCoeffs[2,2];
{data21,data12,data22}
```
**Output:**
```
Out= {{2, 1, {"alpha" -> 6, "beta" -> 0, "exact?" -> True}, {"alpha" -> -250/7, "beta" -> 0, "exact?" -> True}, {"alpha" -> 5996/147, "beta" -> 1, "exact?" -> True}, {"alpha" -> 5051/147, "beta" -> -125/42, "exact?" -> True}, {"alpha" -> -18688/441, "beta" -> 193/42, "exact?" -> True}}, {1, 2, {"alpha" -> 6, "beta" -> 0, "exact?" -> True}, {"alpha" -> -250/7, "beta" -> 0, "exact?" -> True}, {"alpha" -> 6164/147, "beta" -> 1, "exact?" -> True}, {"alpha" -> 4747/147, "beta" -> -125/42, "exact?" -> True}, {"alpha" -> -6364/147, "beta" -> 197/42, "exact?" -> True}}, {2, 2, {"alpha" -> 6, "beta" -> 0, "exact?" -> True}, {"alpha" -> -292/7, "beta" -> 0, "exact?" -> True}, {"alpha" -> 7886/147, "beta" -> 1, "exact?" -> True}, {"alpha" -> 3212/49, "beta" -> -73/21, "exact?" -> True}, {"alpha" -> -5592/49, "beta" -> 547/84, "exact?" -> True}}}
```

## 2026-07-20 10:43:43
**Input:**
```wolfram
Clear[θ1,θ2,h];
ClearAll[qmin,qming,tg,tgCoeffs];
qmin[λ1_, λ2_][II_List, JJ_List][u_] := qmin[λ1, λ2][II, JJ][u] = Sum[Signature[σ] Dot @@ Table[T[λ1, λ2][II[[σ[[k]]]], JJ[[k]]][u - (k - 1) h], {k, Length[II]}], {σ, Permutations[Range[Length[II]]]}];
qming[λ1_, λ2_][II_List, JJ_List][u_] := qming[λ1, λ2][II, JJ][u] = Sum[Signature[σ] Dot @@ Table[Tg[λ1, λ2][II[[σ[[k]]]], JJ[[k]]][u - (k - 1) h], {k, Length[II]}], {σ, Permutations[Range[Length[II]]]}];
tg[λ1_,λ2_][a_,1][u_]:= Total[qming[λ1, λ2][#, #][u] & /@ Subsets[{1, 2, 3, 4}, {a}]]
tgCoeffs[λ1_, λ2_] := tgCoeffs[λ1, λ2] = Module[{uu, m1, m2, m3, cl1, cl2, cl3},
   m1 = tg[λ1, λ2][1, 1][uu];
   m2 = tg[λ1, λ2][2, 1][uu];
   m3 = tg[λ1, λ2][3, 1][uu];
   cl1 = Map[PadRight[CoefficientList[#, uu], 3] &, m1, {2}];
   cl2 = Map[PadRight[CoefficientList[#, uu], 5] &, m2, {2}];
   cl3 = Map[PadRight[CoefficientList[#, uu], 7] &, m3, {2}];
   {Table[cl1[[All, All, k]], {k, 1, 3}], Table[cl2[[All, All, k]], {k, 1, 5}], Table[cl3[[All, All, k]], {k, 1, 7}]}
];
{θ1,θ2,h,z[1]}
```
**Output:**
```
[message] RuleDelayed::rhs: Pattern λ1_ appears on the right-hand side of rule tg[λ1_, λ2_][a_, 1][u_] :> (Total[(qming[λ1, λ2][#1, #1][u] & ) /@ Subsets[{1, 2, 3, 4}, {a}]] tgCoeffs[λ1_, λ2_] := tgCoeffs[λ1, λ… [620 chars]
[message] [RuleDelayed::rhs ×1 more suppressed — first instance shown above]
Out= {θ1, θ2, h, 1}
```

## 2026-07-20 10:43:53
**Input:**
```wolfram
ClearAll[tg,tgCoeffs];
tg[λ1_,λ2_][a_,1][u_]:= Total[qming[λ1, λ2][#, #][u] & /@ Subsets[{1, 2, 3, 4}, {a}]];
tgCoeffs[λ1_, λ2_] := tgCoeffs[λ1, λ2] = Module[{uu, m1, m2, m3, cl1, cl2, cl3},
   m1 = tg[λ1, λ2][1, 1][uu];
   m2 = tg[λ1, λ2][2, 1][uu];
   m3 = tg[λ1, λ2][3, 1][uu];
   cl1 = Map[PadRight[CoefficientList[#, uu], 3] &, m1, {2}];
   cl2 = Map[PadRight[CoefficientList[#, uu], 5] &, m2, {2}];
   cl3 = Map[PadRight[CoefficientList[#, uu], 7] &, m3, {2}];
   {Table[cl1[[All, All, k]], {k, 1, 3}], Table[cl2[[All, All, k]], {k, 1, 5}], Table[cl3[[All, All, k]], {k, 1, 7}]}
];
Head[tg[{1,0,0,0},{1,0,0,0}][1,1][u]]
```
**Output:**
```
Out= List
```

## 2026-07-20 10:44:05
**Input:**
```wolfram
Module[{λ1={1,1,0,0},λ2={1,1,0,0},d,c2,C2m,C4m,u0mat},
  d=dim[λ1]dim[λ2];
  c2=tgCoeffs[λ1,λ2][[2]];
  C2m=C2op[λ1,λ2]; C4m=C4op[λ1,λ2];
  u0mat=Simplify[c2[[1]]-(1/8)C2m.C2m+(1/4)C4m];
  {"u^4"->solveAB[Simplify[c2[[5]]],d,C2m],
   "u^3"->solveAB[Simplify[c2[[4]]],d,C2m],
   "u^2"->solveAB[Simplify[c2[[3]]],d,C2m],
   "u^1"->solveAB[Simplify[c2[[2]]],d,C2m],
   "u^0"->solveAB[u0mat,d,C2m]}]
```
**Output:**
```
Out= {"u^4" -> {"alpha" -> 6, "beta" -> 0, "exact?" -> True}, "u^3" -> {"alpha" -> 8*(-2*h - θ1 - θ2) - 4*(2*h + θ1 + θ2), "beta" -> (6*(-2*h - θ1 - θ2))/13 + (6*(2*h + θ1 + θ2))/13, "exact?" -> {{-8*(-2*h - θ1 - θ2) - 8*(2*h + θ1 + θ2) - 16*((6*(-2*h - θ1 - θ2))/13 + (6*(2*h + θ1 + θ2))/13), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, -8*(-2*h - θ1 - θ2) - 8*(2*h + θ1 + θ2) - 14*((6*(-2*h - θ1 - θ2))/13 + (6*(2*h + θ1 + θ2))/13), 0, 0, 0, 0, -2*((6*(-2*h - θ1 - θ2))/13 + (6*(2*h + θ1 + θ2))/13), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, -8*(-2*h - θ1 - θ2) - 8*(2*h + θ1 + θ2) - 14*((6*(-2*h - θ1 - θ2))/13 + (6*(2*h + θ1 + θ2))/13), 0, 0, 0, 0, 0, 0, 0, 0, 0, -2*((6*(-2*h - θ1 - θ2))/13 + (6*(2*h + θ1 + θ2))/13), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, -8*(-2*h - θ1 - θ2) - 8*(2*h + θ1 + θ2) - 14*((6*(-2*h - θ1 - θ2))/13 + (6*(2*h + θ1 + θ2))/13), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2*((6*(-2*h - θ1 - θ2))/13 + (6*(2*h + θ1 + θ2))/13), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, -8*(-2*h - θ1 - θ2) - 8*(2*h + θ1 + θ2) - 14*((6*(-2*h - θ1 - θ2))/13 + (6*(2*h + θ1 + θ2))/13), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2*((6*(-2*h - θ1 - θ2))/13 + (6*(2*h + θ1 + θ2))/13), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, -8*(-2*h - θ1 - θ2) - 8*(2*h + θ1 + θ2) - 12*((6*(-2*h - θ1 - θ2))/13 + (6*(2*h + θ1 + θ2))/13), 0, 0, 0, 0, (-6*(-2*h - θ1 - θ2))/13 - (6*(2*h + θ1 + θ2))/13, 0, 0, 0, 0, (6*(-2*h - θ1 - θ2))/13 + (6*(2*h + θ1 + θ2))/13, 0, 0, 0, 0, (6*(-2*h - θ1 - θ2))/13 + (6*(2*h + θ1 + θ2))/13, 0, 0, 0, 0, (-6*(-2*h - θ1 - θ2))/13 - (6*(2*h + θ1 + θ2))/13, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, -2*((6*(-2*h - θ1 - θ2))/13 + (6*(2*h + θ1 + θ2))/13), 0, 0, 0, 0, -8*(-2*h - θ1 - θ2) - 8*(2*h + θ1 + θ2) - 14*((6*(-2*h - θ1 - θ2))/13 + (6*(2*h + θ1 + θ2))/13), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, -8*(-2*h - θ1 - θ2) - 8*(2*h + θ1 + θ2) - 16*((6*(-2*h - θ1 - θ2))/13 + (6*(2*h + θ1 + θ2))/13), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, -8*(-2*h - θ1 - θ2) - 8*(2*h + θ1 + θ2) - 14*((6*(-2*h - θ1 - θ2))/13 + (6*(2*h + θ1 + θ2))/13), 0, 0, 0, 0, -2*((6*(-2*h - θ1 - θ2))/13 + (6*(2*h + θ1 + θ2))/13), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, -8*(-2*h - θ1 - θ2) - 8*(2*h + θ1 + θ2) - 14*((6*(-2*h - θ1 - θ2))/13 + (6*(2*h + θ1 + θ2))/13), 0, 0, 0, 0, 0, 0, 0, 0, 0, -2*((6*(-2*h - θ1 - θ2))/13 + (6*(2*h + θ1 + θ2))/13), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, -4*((6*(-2*h - θ1 - θ2))/13 + (6*(2*h + θ1 + θ2))/13), 0, 0, 0, 0, -8*(-2*h - θ1 - θ2) - 8*(2*h + θ1 + θ2) - 12*((6*(-2*h - θ1 - θ2))/13 + (6*(2*h + θ1 + θ2))/13), 0, 0, 0, 0, -2*((6*(-2*h - θ1 - θ2))/13 + (6*(2*h + θ1 + θ2))/13), 0, 0, 0, 0, -2*((6*(-2*h - θ1 - θ2))/13 + (6*(2*h + θ1 + θ2))/13), 0, 0, 0, 0, 0, 0, 0, 0, 0, -4*((6*(-2*h - θ1 - θ2))/13 + (6*(2*h + θ1 + θ2))/13), 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -8*(-2*h - θ1 - θ2) - 8*(2*h + θ1 + θ2) - 14*((6*(-2*h - θ1 - θ2))/13 + (6*(2*h + θ1 + θ2))/13), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2*((6*(-2*h - θ1 - θ2))/13 + (6*(2*h + θ1 + θ2))/13), 0, 0, 0, 0}, {0, 0, -2*((6*(-2*h - θ1 - θ2))/13 + (6*(2*h + θ1 + θ2))/13), 0, 0, 0, 0, 0, 0, 0, 0, 0, -8*(-2*h - θ1 - θ2) - 8*(2*h + θ1 + θ2) - 14*((6*(-2*h - θ1 - θ2))/13 + (6*(2*h + θ1 + θ2))/13), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, -2*((6*(-2*h - θ1 - θ2))/13 + (6*(2*h + θ1 + θ2))/13), 0, 0, 0, 0, -8*(-2*h - θ1 - θ2) - 8*(2*h + θ1 + θ2) - 14*((6*(-2*h - θ1 - θ2))/13 + (6*(2*h + θ1 + θ2))/13), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -8*(-2*h - θ1 - θ2) - 8*(2*h + θ1 + θ2) - 16*((6*(-2*h - θ1 - θ2))/13 + (6*(2*h + θ1 + θ
[output truncated — 149067 chars total]
```

## 2026-07-20 10:44:25
**Input:**
```wolfram
solveABsym[M_,d_,C2m_]:=Module[{sol,α,β,eq1,eq2},
  eq1=M[[1,1]]==α+β C2m[[1,1]];
  eq2=M[[2,2]]==α+β C2m[[2,2]];
  sol=Solve[{eq1,eq2},{α,β}][[1]];
  {"alpha"->Simplify[α/.sol],"beta"->Simplify[β/.sol],
   "exact?"->(Simplify[M-(α/.sol) IdentityMatrix[d]-(β/.sol) C2m]==ConstantArray[0,{d,d}])}];
Module[{λ1={1,1,0,0},λ2={1,1,0,0},d,c2,C2m,C4m,u0mat},
  d=dim[λ1]dim[λ2];
  c2=tgCoeffs[λ1,λ2][[2]];
  C2m=C2op[λ1,λ2]; C4m=C4op[λ1,λ2];
  u0mat=Simplify[c2[[1]]-(1/8)C2m.C2m+(1/4)C4m];
  {"u^4"->solveABsym[Simplify[c2[[5]]],d,C2m],
   "u^3"->solveABsym[Simplify[c2[[4]]],d,C2m],
   "u^2"->solveABsym[Simplify[c2[[3]]],d,C2m],
   "u^1"->solveABsym[Simplify[c2[[2]]],d,C2m],
   "u^0"->solveABsym[u0mat,d,C2m]}]
```
**Output:**
```
Out= {"u^4" -> {"alpha" -> 6, "beta" -> 0, "exact?" -> True}, "u^3" -> {"alpha" -> -12*(2*h + θ1 + θ2), "beta" -> 0, "exact?" -> True}, "u^2" -> {"alpha" -> 14*h^2 + 36*h*(θ1 + θ2) + 6*(θ1^2 + 4*θ1*θ2 + θ2^2), "beta" -> h^2, "exact?" -> True}, "u^1" -> {"alpha" -> 2*(10*h^3 - 7*h^2*(θ1 + θ2) - 6*θ1*θ2*(θ1 + θ2) - 6*h*(θ1^2 + 4*θ1*θ2 + θ2^2)), "beta" -> -(h^2*(2*h + θ1 + θ2)), "exact?" -> True}, "u^0" -> {"alpha" -> 2*(10 - 14*h^4 + 3*θ1^2*θ2^2 - 5*h^3*(θ1 + θ2) + 6*h*θ1*θ2*(θ1 + θ2) + h^2*(θ1^2 + 5*θ1*θ2 + θ2^2)), "beta" -> 3/4 + 2*h^4 + h^2*θ1*θ2 + h^3*(θ1 + θ2), "exact?" -> {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 4*(-1 + h^4), 0, 0, 0, 0, 2 - 2*h^4, 0, 0, 0, 0, 2*(-1 + h^4), 0, 0, 0, 0, 2*(-1 + h^4), 0, 0, 0, 0, 2 - 2*h^4, 0, 0, 0, 0, 4*(-1 + h^4), 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 8 - 8*h^4, 0, 0, 0, 0, 4*(-1 + h^4), 0, 0, 0, 0, 4 - 4*h^4, 0, 0, 0, 0, 4 - 4*h^4, 0, 0, 0, 0, 4*(-1 + h^4), 0, 0, 0, 0, 8 - 8*h^4, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 8*(-1 + h^4), 0, 0, 0, 0, 4 - 4*h^4, 0, 0, 0, 0, 4*(-1 + h^4), 0, 0, 0, 0, 4*(-1 + h^4), 0, 0, 0, 0, 4 - 4*h^4, 0, 0, 0, 0, 8*(-1 + h^4), 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 8*(-1 + h^4), 0, 0, 0, 0, 4 - 4*h^4, 0, 0, 0, 0, 4*(-1 + h^4), 0, 0, 0, 0, 4*(-1 + h^4), 0, 0, 0, 0, 4 - 4*h^4, 0, 0, 0, 0, 8*(-1 + h^4), 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 8 - 8*h^4, 0, 0, 0, 0, 4*(-1 + h^4), 0, 0, 0, 0, 4 - 4*h^4, 0, 0, 0, 0, 4 - 4*h^4, 0, 0, 0, 0, 4*(-1 + h^4), 0, 0, 0, 0, 8 - 8*h^4, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
[output truncated — 8852 chars total]
```

## 2026-07-20 10:44:46
**Input:**
```wolfram
Module[{λ1={1,1,0,0},λ2={1,1,0,0},d,c2,C2m,C4m,u0mat},
  d=dim[λ1]dim[λ2];
  c2=tgCoeffs[λ1,λ2][[2]];
  C2m=C2op[λ1,λ2]; C4m=C4op[λ1,λ2];
  u0mat=Simplify[c2[[1]]-(h^4/8)C2m.C2m+(h^4/4)C4m];
  solveABsym[u0mat,d,C2m]]
```
**Output:**
```
Out= {"alpha" -> 2*(-4*h^4 + 3*θ1^2*θ2^2 - 5*h^3*(θ1 + θ2) + 6*h*θ1*θ2*(θ1 + θ2) + h^2*(θ1^2 + 5*θ1*θ2 + θ2^2)), "beta" -> (11*h^4)/4 + h^2*θ1*θ2 + h^3*(θ1 + θ2), "exact?" -> True}
```

## 2026-07-20 10:45:54
**Input:**
```wolfram
Module[{λ1={2,2,0,0},λ2={1,1,0,0},d,c2,C2m,C4m,u0mat},
  d=dim[λ1]dim[λ2];
  c2=tgCoeffs[λ1,λ2][[2]];
  C2m=C2op[λ1,λ2]; C4m=C4op[λ1,λ2];
  u0mat=Simplify[c2[[1]]-(h^4/8)C2m.C2m+(h^4/4)C4m];
  {"d"->d,
   "u^4"->solveABsym[Simplify[c2[[5]]],d,C2m],
   "u^3"->solveABsym[Simplify[c2[[4]]],d,C2m],
   "u^2"->solveABsym[Simplify[c2[[3]]],d,C2m],
   "u^1"->solveABsym[Simplify[c2[[2]]],d,C2m],
   "u^0"->solveABsym[u0mat,d,C2m]}]
```
**Output:**
```
Out= {"d" -> 120, "u^4" -> {"alpha" -> 6, "beta" -> 0, "exact?" -> True}, "u^3" -> {"alpha" -> -6*(5*h + 2*(θ1 + θ2)), "beta" -> 0, "exact?" -> True}, "u^2" -> {"alpha" -> 6*(3*h^2 + 7*h*θ1 + θ1^2 + 8*h*θ2 + 4*θ1*θ2 + θ2^2), "beta" -> h^2, "exact?" -> True}, "u^1" -> {"alpha" -> 47*h^3 - 12*θ1*θ2*(θ1 + θ2) - 2*h^2*(7*θ1 + 11*θ2) - 6*h*(2*θ1^2 + 10*θ1*θ2 + 3*θ2^2), "beta" -> -1/2*(h^2*(5*h + 2*(θ1 + θ2))), "exact?" -> True}, "u^0" -> {"alpha" -> -33*h^4 + 6*θ1^2*θ2^2 + 6*h*θ1*θ2*(2*θ1 + 3*θ2) - h^3*(20*θ1 + 27*θ2) + 2*h^2*(θ1^2 + 5*θ1*θ2 + 3*θ2^2), "beta" -> (h^2*(8*h^2 + 2*h*θ1 + 3*h*θ2 + 2*θ1*θ2))/2, "exact?" -> True}}
```

## 2026-07-20 10:46:16
**Input:**
```wolfram
Module[{λ1={1,1,0,0},λ2={2,2,0,0},d,c2,C2m,C4m,u0mat},
  d=dim[λ1]dim[λ2];
  c2=tgCoeffs[λ1,λ2][[2]];
  C2m=C2op[λ1,λ2]; C4m=C4op[λ1,λ2];
  u0mat=Simplify[c2[[1]]-(h^4/8)C2m.C2m+(h^4/4)C4m];
  {"d"->d,
   "u^4"->solveABsym[Simplify[c2[[5]]],d,C2m],
   "u^3"->solveABsym[Simplify[c2[[4]]],d,C2m],
   "u^2"->solveABsym[Simplify[c2[[3]]],d,C2m],
   "u^1"->solveABsym[Simplify[c2[[2]]],d,C2m],
   "u^0"->solveABsym[u0mat,d,C2m]}]
```
**Output:**
```
Out= {"d" -> 120, "u^4" -> {"alpha" -> 6, "beta" -> 0, "exact?" -> True}, "u^3" -> {"alpha" -> -6*(5*h + 2*(θ1 + θ2)), "beta" -> 0, "exact?" -> True}, "u^2" -> {"alpha" -> 6*(3*h^2 + 8*h*θ1 + θ1^2 + 7*h*θ2 + 4*θ1*θ2 + θ2^2), "beta" -> h^2, "exact?" -> True}, "u^1" -> {"alpha" -> 47*h^3 - 12*θ1*θ2*(θ1 + θ2) - 2*h^2*(11*θ1 + 7*θ2) - 6*h*(3*θ1^2 + 10*θ1*θ2 + 2*θ2^2), "beta" -> -1/2*(h^2*(5*h + 2*(θ1 + θ2))), "exact?" -> True}, "u^0" -> {"alpha" -> -33*h^4 + 6*θ1^2*θ2^2 + 6*h*θ1*θ2*(3*θ1 + 2*θ2) - h^3*(27*θ1 + 20*θ2) + 2*h^2*(3*θ1^2 + 5*θ1*θ2 + θ2^2), "beta" -> (h^2*(8*h^2 + 3*h*θ1 + 2*h*θ2 + 2*θ1*θ2))/2, "exact?" -> True}}
```

## 2026-07-20 10:48:48
**Input:**
```wolfram
Module[{λ1={2,2,0,0},λ2={2,2,0,0},d,c2,C2m,C4m,u0mat},
  d=dim[λ1]dim[λ2];
  c2=tgCoeffs[λ1,λ2][[2]];
  C2m=C2op[λ1,λ2]; C4m=C4op[λ1,λ2];
  u0mat=Simplify[c2[[1]]-(h^4/8)C2m.C2m+(h^4/4)C4m];
  {"d"->d,
   "u^4"->solveABsym[Simplify[c2[[5]]],d,C2m],
   "u^3"->solveABsym[Simplify[c2[[4]]],d,C2m],
   "u^2"->solveABsym[Simplify[c2[[3]]],d,C2m],
   "u^1"->solveABsym[Simplify[c2[[2]]],d,C2m],
   "u^0"->solveABsym[u0mat,d,C2m]}]
```
**Output:**
```
Out= {"d" -> 400, "u^4" -> {"alpha" -> 6, "beta" -> 0, "exact?" -> True}, "u^3" -> {"alpha" -> -12*(3*h + θ1 + θ2), "beta" -> 0, "exact?" -> True}, "u^2" -> {"alpha" -> 26*h^2 + 54*h*(θ1 + θ2) + 6*(θ1^2 + 4*θ1*θ2 + θ2^2), "beta" -> h^2, "exact?" -> True}, "u^1" -> {"alpha" -> 2*(42*h^3 - 13*h^2*(θ1 + θ2) - 6*θ1*θ2*(θ1 + θ2) - 9*h*(θ1^2 + 4*θ1*θ2 + θ2^2)), "beta" -> -(h^2*(3*h + θ1 + θ2)), "exact?" -> True}, "u^0" -> {"alpha" -> 2*(-48*h^4 + 3*θ1^2*θ2^2 - 21*h^3*(θ1 + θ2) + 9*h*θ1*θ2*(θ1 + θ2) + h^2*(3*θ1^2 + 7*θ1*θ2 + 3*θ2^2)), "beta" -> (h^2*(23*h^2 + 4*θ1*θ2 + 6*h*(θ1 + θ2)))/4, "exact?" -> True}}
```

## 2026-07-20 10:49:21
**Input:**
```wolfram
{dim[{3,3,0,0}], dim[{3,3,0,0}]*dim[{1,1,0,0}]}
```
**Output:**
```
Out= {50, 300}
```

## 2026-07-20 10:50:42
**Input:**
```wolfram
Module[{λ1={3,3,0,0},λ2={1,1,0,0},d,c2,C2m,C4m,u0mat},
  d=dim[λ1]dim[λ2];
  c2=tgCoeffs[λ1,λ2][[2]];
  C2m=C2op[λ1,λ2]; C4m=C4op[λ1,λ2];
  u0mat=Simplify[c2[[1]]-(h^4/8)C2m.C2m+(h^4/4)C4m];
  {"d"->d,
   "u^4"->solveABsym[Simplify[c2[[5]]],d,C2m],
   "u^3"->solveABsym[Simplify[c2[[4]]],d,C2m],
   "u^2"->solveABsym[Simplify[c2[[3]]],d,C2m],
   "u^1"->solveABsym[Simplify[c2[[2]]],d,C2m],
   "u^0"->solveABsym[u0mat,d,C2m]}]
```
**Output:**
```
Out= {"d" -> 300, "u^4" -> {"alpha" -> 6, "beta" -> 0, "exact?" -> True}, "u^3" -> {"alpha" -> -12*(3*h + θ1 + θ2), "beta" -> 0, "exact?" -> True}, "u^2" -> {"alpha" -> 20*h^2 + 12*h*(4*θ1 + 5*θ2) + 6*(θ1^2 + 4*θ1*θ2 + θ2^2), "beta" -> h^2, "exact?" -> True}, "u^1" -> {"alpha" -> 2*(47*h^3 - 6*θ1*θ2*(θ1 + θ2) - 5*h^2*(θ1 + 3*θ2) - 6*h*(θ1^2 + 6*θ1*θ2 + 2*θ2^2)), "beta" -> -(h^2*(3*h + θ1 + θ2)), "exact?" -> True}, "u^0" -> {"alpha" -> 2*(-42*h^4 + 3*θ1^2*θ2^2 + 6*h*θ1*θ2*(θ1 + 2*θ2) - h^3*(17*θ1 + 30*θ2) + h^2*(θ1^2 + 3*θ1*θ2 + 6*θ2^2)), "beta" -> (21*h^4)/4 + h^2*θ1*θ2 + h^3*(θ1 + 2*θ2), "exact?" -> True}}
```

## 2026-07-20 10:53:54
**Input:**
```wolfram
Module[{λ1={1,1,0,0},λ2={1,1,0,0}},
  {"T[1,2].HWS"->Simplify[T[λ1,λ2][1,2][u].HWS[λ1,λ2]],
   "T[2,1].HWS"->Simplify[T[λ1,λ2][2,1][u].HWS[λ1,λ2]]}]
```
**Output:**
```
Out= {"T[1,2].HWS" -> {{0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}}, "T[2,1].HWS" -> {{0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}}}
```

## 2026-07-20 10:54:21
**Input:**
```wolfram
vac2[S_,P_,u_]:=Sum[ν[{S,S,0,0},{P,P,0,0}][i][u] ν[{S,S,0,0},{P,P,0,0}][j][u-h],{i,1,4},{j,i+1,4}];
Module[{λ1={1,1,0,0},λ2={1,1,0,0}},
  Simplify[(tg[λ1,λ2][2,1][u].HWS[λ1,λ2])[[1,1]] - vac2[1,1,u]]]
```
**Output:**
```
Out= 4*h^2*(-((u - θ1)*(-2*h + u - θ1)) - (h - u + θ2)^2)
```

## 2026-07-20 10:55:01
**Input:**
```wolfram
Module[{λ1={1,1,0,0},λ2={1,1,0,0}},
  {"qming[{1,2},{1,2}][u].HWS"->Simplify[(qming[λ1,λ2][{1,2},{1,2}][u].HWS[λ1,λ2])[[1,1]]],
   "nu1(u)nu2(u-h) predicted"->Simplify[ν[λ1,λ2][1][u] ν[λ1,λ2][2][u-h]],
   "diff"->Simplify[(qming[λ1,λ2][{1,2},{1,2}][u].HWS[λ1,λ2])[[1,1]] - ν[λ1,λ2][1][u] ν[λ1,λ2][2][u-h]]}]
```
**Output:**
```
Out= {"qming[{1,2},{1,2}][u].HWS" -> (h - u + θ1)*(2*h - u + θ1)*(h - u + θ2)*(2*h - u + θ2), "nu1(u)nu2(u-h) predicted" -> (h - u + θ1)*(2*h - u + θ1)*(h - u + θ2)*(2*h - u + θ2), "diff" -> 0}
```

## 2026-07-20 10:55:13
**Input:**
```wolfram
Module[{λ1={1,1,0,0},λ2={1,1,0,0},pairs},
  pairs={{1,3},{1,4},{2,3},{2,4},{3,4}};
  Table[{p,Simplify[(qming[λ1,λ2][p,p][u].HWS[λ1,λ2])[[1,1]] - ν[λ1,λ2][p[[1]]][u] ν[λ1,λ2][p[[2]]][u-h]]},{p,pairs}]]
```
**Output:**
```
Out= {{{1, 3}, h^2*(-((u - θ1)*(-2*h + u - θ1)) - (h - u + θ2)^2)}, {{1, 4}, h^2*(-((u - θ1)*(-2*h + u - θ1)) - (h - u + θ2)^2)}, {{2, 3}, h^2*(-((u - θ1)*(-2*h + u - θ1)) - (h - u + θ2)^2)}, {{2, 4}, h^2*(-((u - θ1)*(-2*h + u - θ1)) - (h - u + θ2)^2)}, {{3, 4}, 0}}
```

## 2026-07-20 10:59:52
**Input:**
```wolfram
(*T[i,j][u] = (u-th1)(u-th2) delta_ij Id - h(u-th1) EE2[j,i] - h(u-th2) EE1[j,i] + h^2 Mmix[j,i]
  where EE1[a,b]=Ee1[a,b](x)id2, EE2[a,b]=id1(x)Ee2[a,b], Mmix[a,b]=Sum_k Ee1[k,b](x)Ee2[a,k]
  (derived from L[lam][i,j][v]=v*delta_ij*Id-h*Ee[lam][j,i], T[i,j][u]=Sum_ii KP[L1[i,ii][u-th1],L2[ii,j][u-th2])].
  VERIFY against the notebook's actual T[lam1,lam2][i,j][u] before using it for anything.*)
Module[{λ1={1,1,0,0},λ2={1,1,0,0},EE1,EE2,Mmix,Tpred,worst},
  EE1[a_,b_]:=KroneckerProduct[Ee[λ1][a,b],id[λ2]];
  EE2[a_,b_]:=KroneckerProduct[id[λ1],Ee[λ2][a,b]];
  Mmix[a_,b_]:=Sum[KroneckerProduct[Ee[λ1][k,b],Ee[λ2][a,k]],{k,4}];
  Tpred[i_,j_][u_]:=(u-θ1)(u-θ2)KroneckerDelta[i,j] IdentityMatrix[36] - h(u-θ1)EE2[j,i] - h(u-θ2)EE1[j,i] + h^2 Mmix[j,i];
  worst=Max[Table[Max[Abs[Flatten[Simplify[Tpred[i,j][u]-T[λ1,λ2][i,j][u]]]]],{i,4},{j,4}]];
  {"worst mismatch (should be 0)"->worst}]
```
**Output:**
```
Out= {"worst mismatch (should be 0)" -> 0}
```

## 2026-07-20 11:01:41
**Input:**
```wolfram
Module[{λ1={1,1,0,0},λ2={1,1,0,0},d,EE1,EE2,Mmix,Tp,tau2exp,c2},
  d=36;
  EE1[a_,b_]:=KroneckerProduct[Ee[λ1][a,b],id[λ2]];
  EE2[a_,b_]:=KroneckerProduct[id[λ1],Ee[λ2][a,b]];
  Mmix[a_,b_]:=Sum[KroneckerProduct[Ee[λ1][k,b],Ee[λ2][a,k]],{k,4}];
  Tp[i_,j_][u_]:=(u-θ1)(u-θ2)KroneckerDelta[i,j] IdentityMatrix[d] - h(u-θ1)EE2[j,i] - h(u-θ2)EE1[j,i] + h^2 Mmix[j,i];
  tau2exp[u_]:=Sum[Tp[i,i][u].Tp[j,j][u-h]-Tp[j,i][u].Tp[i,j][u-h],{i,4},{j,i+1,4}];
  c2=CoefficientList[tau2exp[uu]/.u->uu,uu];
  {"matches tg[2,1]?"->(Simplify[tau2exp[u]-tg[λ1,λ2][2,1][u]]==ConstantArray[0,{d,d}]),
   "u1 residual after known C2 piece"->Simplify[c2[[2]]-((-h^2(2h+θ1+θ2))C2op[λ1,λ2])-((2*(10 h^3-7h^2(θ1+θ2)-6θ1 θ2(θ1+θ2)-6h(θ1^2+4θ1 θ2+θ2^2)))IdentityMatrix[d])]}]
```
**Output:**
```
[message] Thread::tdlen: Objects of unequal length in {} + {16 h  (2 h + θ1 + θ2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0} + {-2 (10 h  - 7 h  (θ… [382 chars]
[message] [Thread::tdlen ×2 more suppressed — first instance shown above]
[message] General::stop: Further output of Thread::tdlen will be suppressed during this calculation.
Out= {"matches tg[2,1]?" -> True, "u1 residual after known C2 piece" -> {{} + {16*h^2*(2*h + θ1 + θ2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0} + {2*(-10*h^3 + 7*h^2*(θ1 + θ2) + 6*θ1*θ2*(θ1 + θ2) + 6*h*(θ1^2 + 4*θ1*θ2 + θ2^2)), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {2*(3*θ1^2*θ2^2 + 2*h^3*(θ1 + θ2) + 6*h*θ1*θ2*(θ1 + θ2) + h^2*(θ1^2 + 12*θ1*θ2 + θ2^2)), -4*(2*h^3 + 7*h^2*(θ1 + θ2) + 3*θ1*θ2*(θ1 + θ2) + 3*h*(θ1^2 + 4*θ1*θ2 + θ2^2)), 28*h^2 + 36*h*(θ1 + θ2) + 6*(θ1^2 + 4*θ1*θ2 + θ2^2), -12*(2*h + θ1 + θ2), 6} + {0, 4*(2*h^3 + 7*h^2*(θ1 + θ2) + 3*θ1*θ2*(θ1 + θ2) + 3*h*(θ1^2 + 4*θ1*θ2 + θ2^2)), 0, 0, 0, 0, 2*h^2*(2*h + θ1 + θ2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {} + {0, 0, 14*h^2*(2*h + θ1 + θ2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 2*h^2*(2*h + θ1 + θ2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0} + {0, 0, 2*(-10*h^3 + 7*h^2*(θ1 + θ2) + 6*θ1*θ2*(θ1 + θ2) + 6*h*(θ1^2 + 4*θ1*θ2 + θ2^2)), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {} + {0, 0, 0, 14*h^2*(2*h + θ1 + θ2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2*h^2*(2*h + θ1 + θ2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0} + {0, 0, 0, 2*(-10*h^3 + 7*h^2*(θ1 + θ2) + 6*θ1*θ2*(θ1 + θ2) + 6*h*(θ1^2 + 4*θ1*θ2 + θ2^2)), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {} + {0, 0, 0, 0, 14*h^2*(2*h + θ1 + θ2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2*h^2*(2*h + θ1 + θ2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0} + {0, 0, 0, 0, 2*(-10*h^3 + 7*h^2*(θ1 + θ2) + 6*θ1*θ2*(θ1 + θ2) + 6*h*(θ1^2 + 4*θ1*θ2 + θ2^2)), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {} + {0, 0, 0, 0, 0, 12*h^2*(2*h + θ1 + θ2), 0, 0, 0, 0, h^2*(2*h + θ1 + θ2), 0, 0, 0, 0, -(h^2*(2*h + θ1 + θ2)), 0, 0, 0, 0, -(h^2*(2*h + θ1 + θ2)), 0, 0, 0, 0, h^2*(2*h + θ1 + θ2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0} + {0, 0, 0, 0, 0, 2*(-10*h^3 + 7*h^2*(θ1 + θ2) + 6*θ1*θ2*(θ1 + θ2) + 6*h*(θ1^2 + 4*θ1*θ2 + θ2^2)), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {2*h^2*(2*h^2 + θ1*θ2 + h*(θ1 + θ2)), -2*h^2*(2*h + θ1 + θ2), 2*h^2} + {0, 2*h^2*(2*h + θ1 + θ2), 0, 0, 0, 0, 4*(2*h^3 + 7*h^2*(θ1 + θ2) + 3*θ1*θ2*(θ1 + θ2) + 3*h*(θ1^2 + 4*θ1*θ2 + θ2^2)), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {} + {0, 0, 0, 0, 0, 0, 0, 16*h^2*(2*h + θ1 + θ2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0} + {0, 0, 0, 0, 0, 0, 0, 2*(-10*h^3 + 7*h^2*(θ1 + θ2) + 6*θ1*θ2*(θ1 + θ2) + 6*h*(θ1^2 + 4*θ1*θ2 + θ2^2)), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {} + {0, 0, 0, 0, 0, 0, 0, 0, 14*h^2*(2*h + θ1 + θ2), 0, 0, 0, 0, 2*h^2*(2*h + θ1 + θ2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0} + {0, 0, 0, 0, 0, 0, 0, 0, 2*(-10*h^3 + 7*h^2*(θ1 + θ2) + 6*θ1*θ2*(θ1 + θ2) + 6*h*(θ1^2 + 4*θ1*θ2 + θ2^2)), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {} + {0, 0, 0, 0, 0, 0, 0, 0, 0, 14*h^2*(2*h + θ1 + θ2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 2*h^2*(2*h + θ1 + θ2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0} + {0, 0, 0, 0, 0, 0, 0, 0, 0, 2*(-10*h^3 + 7*h^2*(θ1 + θ2) + 6*θ1*θ2*(θ1 + θ2) + 6*h*(θ1^2 + 4*θ1*θ2 + θ2^2)), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {} + {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2*(-10*h^3 + 7*h^2*(θ1 + θ2) + 6*θ1*θ2*(θ1 + θ2) + 6*h*(θ1^2 + 4*θ1*θ2 + θ2^2)), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0} + {0, 0, 0, 0, 0, 4*h^2*(2*h + θ1 + θ2), 0, 0, 0, 0, 12*h^2*(2*h + θ1 + θ2), 0, 0, 0, 0, 2*h^2*(2*h + θ1 + θ2), 0, 0, 0, 0, 2*h^2*(2*h + θ1 + θ2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 4*h^2*(2*h + θ1 + θ2), 0, 0, 0, 0, 0}, {} + {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 14*h^2*(2*h + θ1 + θ2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
[output truncated — 12769 chars total]
```

## 2026-07-20 11:02:19
**Input:**
```wolfram
{dim[{4,4,0,0}], dim[{4,4,0,0}]*dim[{1,1,0,0}], dim[{3,3,0,0}]*dim[{2,2,0,0}]}
```
**Output:**
```
Out= {105, 630, 1000}
```

## 2026-07-20 11:03:09
**Input:**
```wolfram
Module[{λ1={1,1,0,0},λ2={1,1,0,0},d,c1_1,C2m,rem,mats,names,sol},
  d=36;
  c1_1=Simplify[tgCoeffs[λ1,λ2][[2,2]]]; (*u^1 coeff of tau2*)
  C2m=C2op[λ1,λ2];
  (*known piece: alpha1*Id+beta1*C2op, but alpha1(S,P) unknown -- subtract ONLY the beta1 piece (known)*)
  rem=Simplify[c1_1-(-h^2(2h+θ1+θ2))C2m];
  mats={IdentityMatrix[d],X1[λ1,λ2],X2[λ1,λ2],Y1[λ1,λ2],Y2[λ1,λ2],C2m,C2m.C2m};
  names={"Id","X1","X2","Y1","Y2","C2","C2sq"};
  sol=Quiet@LinearSolve[Transpose[Flatten/@mats],Flatten[rem]];
  {"sol"->Thread[names->sol],
   "exact?"->(Simplify[rem-sol.mats]==ConstantArray[0,{d,d}])}]
```
**Output:**
```
[message] Module::lvsym: Local variable specification {λ1 = {1, 1, 0, 0}, λ2 = {1, 1, 0, 0}, d, (c1_) 1, C2m, rem, mats, names, sol} contains (c1_) 1, which is not a symbol or an assignment to a symbol.
Out= Module[{λ1 = {1, 1, 0, 0}, λ2 = {1, 1, 0, 0}, d, (c1_)*1, C2m, rem, mats, names, sol}, d = 36; (c1_)*1 = Simplify[tgCoeffs[λ1, λ2][[2,2]]]; C2m = C2op[λ1, λ2]; rem = Simplify[(c1_)*1 - ((-h^2)*(2*h + θ1 + θ2))*C2m]; mats = {IdentityMatrix[d], X1[λ1, λ2], X2[λ1, λ2], Y1[λ1, λ2], Y2[λ1, λ2], C2m, C2m . C2m}; names = {"Id", "X1", "X2", "Y1", "Y2", "C2", "C2sq"}; sol = Quiet[LinearSolve[Transpose[Flatten /@ mats], Flatten[rem]]]; {"sol" -> Thread[names -> sol], "exact?" -> Simplify[rem - sol . mats] == ConstantArray[0, {d, d}]}]
```

## 2026-07-20 11:03:20
**Input:**
```wolfram
Module[{λ1={1,1,0,0},λ2={1,1,0,0},d,c1u1,C2m,rem,mats,names,sol},
  d=36;
  c1u1=Simplify[tgCoeffs[λ1,λ2][[2,2]]];
  C2m=C2op[λ1,λ2];
  rem=Simplify[c1u1-(-h^2(2h+θ1+θ2))C2m];
  mats={IdentityMatrix[d],X1[λ1,λ2],X2[λ1,λ2],Y1[λ1,λ2],Y2[λ1,λ2],C2m,C2m.C2m};
  names={"Id","X1","X2","Y1","Y2","C2","C2sq"};
  sol=Quiet@LinearSolve[Transpose[Flatten/@mats],Flatten[rem]];
  {"sol"->Thread[names->sol],
   "exact?"->(Simplify[rem-sol.mats]==ConstantArray[0,{d,d}])}]
```
**Output:**
```
[message] Transpose::nmtx: The first two levels of {{1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0… [12039 chars]
Out= {"sol" -> {"Id" -> LinearSolve[Transpose[{{1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, X1[{1, 1, 0, 0}, {1, 1, 0, 0}], X2[{1, 1, 0, 0}, {1, 1, 0, 0}], Y1[{1, 1, 0, 0}, {1, 1, 0, 0}], Y2[{1, 1, 0, 0}, {1, 1, 0, 0}], {16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
[output truncated — 40132207 chars total]
```
## 2026-07-20 11:03:50 — ▶️ RUN CELL 66 (0.01 s)
**In [66]:** `(*GENERAL master formula for ALL of C2op,C3op,C4op (any n): the operator-valued "split Casimir"↵  decomposition of Cn_op = Sum_{i1..in} EE[i1,i2].EE[i2,i3]...EE[in,i1] (EE=E1(x)Id+Id(x)E2).↵  Substitu`
**Out:** Out[219]= \left\{\text{CnMaster[2,..]==C2op?}\to \mathrm{True},\text{CnMaster[3,..]==C3op?}\to \mathrm{True},\text{CnMaster[4,..]==C4op?}\to \mathrm{True}\right\}

## 2026-07-20 11:04:09 — ▶️ RUN CELL 67 (4.81 s)
**In [67]:** `(*Reduced explicit "named cross-term" forms, grouping the 2^n patterns by which OPERATOR they↵  produce (found by GatherBy on an asymmetric test pair, so no accidental C3site1=C3site2-type↵  degenerac`
**Out:** Out[228]= \begin{pmatrix}2\,Omega2\left(\left\{1,1,0,0\right\},\left\{2,2,0,0\right\}\right)+22 & 2\,Omega2\left(\left\{1,1,0,0\right\},\left\{2,2,0,0\right\}\right) & 2\,Omega2\left(\left\{1,1,0,0\right\},\left\{2,2,0,0\right\}\right) & 2\,Omega2\left(\left\{1,1,0,0\right\},\left\{2,2,0,0\right\}\r


## 2026-07-20 11:04:22
**Input:**
```wolfram
Omega2[λ1_,λ2_]:=Omega2[λ1,λ2]=Sum[KroneckerProduct[Ee[λ1][i,j],Ee[λ2][j,i]],{i,4},{j,4}];
{MatrixQ[X1[{1,1,0,0},{1,1,0,0}]],MatrixQ[X2[{1,1,0,0},{1,1,0,0}]],MatrixQ[Y1[{1,1,0,0},{1,1,0,0}]],MatrixQ[Y2[{1,1,0,0},{1,1,0,0}]],MatrixQ[Omega2[{1,1,0,0},{1,1,0,0}]]}
```
**Output:**
```
Out= {False, False, False, False, True}
```
## 2026-07-20 11:04:46 — ▶️ RUN CELL 67 (6.00 s)
**In [67]:** `(*Reduced explicit "named cross-term" forms, grouping the 2^n patterns by which OPERATOR they↵  produce (found by GatherBy on an asymmetric test pair, so no accidental C3site1=C3site2-type↵  degenerac`
**Out:** Out[237]= \begin{pmatrix}2\,X1\left(\left\{1,1,0,0\right\},\left\{2,2,0,0\right\}\right)+2\,X2\left(\left\{1,1,0,0\right\},\left\{2,2,0,0\right\}\right)+Y1\left(\left\{1,1,0,0\right\},\left\{2,2,0,0\right\}\right)+Y2\left(\left\{1,1,0,0\right\},\left\{2,2,0,0\right\}\right)+82 & 2\,X1\left(\left\{1,


## 2026-07-20 11:05:10
**Input:**
```wolfram
X1[λ1_,λ2_]:=X1[λ1,λ2]=Sum[KroneckerProduct[Ee[λ1][i,j].Ee[λ1][j,k],Ee[λ2][k,i]],{i,4},{j,4},{k,4}];
X2[λ1_,λ2_]:=X2[λ1,λ2]=Sum[KroneckerProduct[Ee[λ1][i,j],Ee[λ2][j,k].Ee[λ2][k,i]],{i,4},{j,4},{k,4}];
Y1[λ1_,λ2_]:=Y1[λ1,λ2]=Sum[KroneckerProduct[Ee[λ1][i,j].Ee[λ1][k,i],Ee[λ2][j,k]],{i,4},{j,4},{k,4}];
Y2[λ1_,λ2_]:=Y2[λ1,λ2]=Sum[KroneckerProduct[Ee[λ1][j,k],Ee[λ2][i,j].Ee[λ2][k,i]],{i,4},{j,4},{k,4}];
{MatrixQ[X1[{1,1,0,0},{1,1,0,0}]],MatrixQ[X2[{1,1,0,0},{1,1,0,0}]],MatrixQ[Y1[{1,1,0,0},{1,1,0,0}]],MatrixQ[Y2[{1,1,0,0},{1,1,0,0}]]}
```
**Output:**
```
Out= {True, True, True, True}
```
## 2026-07-20 11:05:20 — 📥 BULK INSERT 1 CELL(S) at positions 67–67
1. [code] (*BUGFIX (found 2026-07-20): these were referenced in the next cell's verification Module but never↵


## 2026-07-20 11:05:21
**Input:**
```wolfram
(*BUGFIX (found 2026-07-20): these were referenced in the next cell's verification Module but never
  actually defined in any persisted notebook cell -- only existed as a scratch evaluateExpression
  call during the original live derivation, so a fresh kernel run of this notebook would fail at
  the next cell. Restoring the definitions here, in their rightful place.*)
X1[λ1_, λ2_] := X1[λ1, λ2] = Sum[KroneckerProduct[Ee[λ1][i, j] . Ee[λ1][j, k], Ee[λ2][k, i]], {i, 4}, {j, 4}, {k, 4}];
X2[λ1_, λ2_] := X2[λ1, λ2] = Sum[KroneckerProduct[Ee[λ1][i, j], Ee[λ2][j, k] . Ee[λ2][k, i]], {i, 4}, {j, 4}, {k, 4}];
Y1[λ1_, λ2_] := Y1[λ1, λ2] = Sum[KroneckerProduct[Ee[λ1][i, j] . Ee[λ1][k, i], Ee[λ2][j, k]], {i, 4}, {j, 4}, {k, 4}];
Y2[λ1_, λ2_] := Y2[λ1, λ2] = Sum[KroneckerProduct[Ee[λ1][j, k], Ee[λ2][i, j] . Ee[λ2][k, i]], {i, 4}, {j, 4}, {k, 4}];
```
**Output:**
```
(no output)
```
## 2026-07-20 11:05:27 — ▶️ RUN CELL 68 (0.02 s)
**In [68]:** `(*Reduced explicit "named cross-term" forms, grouping the 2^n patterns by which OPERATOR they↵  produce (found by GatherBy on an asymmetric test pair, so no accidental C3site1=C3site2-type↵  degenerac`
**Out:** Out[250]= \left\{\text{C2op reduced form matches}\to \mathrm{True},\text{C3op reduced form matches}\to \mathrm{True},\text{C4op reduced form matches}\to \mathrm{True}\right\}


## 2026-07-20 11:05:49
**Input:**
```wolfram
Module[{λ1={1,1,0,0},λ2={1,1,0,0},d,c1u1,C2m,rem,mats,names,sol},
  d=36;
  c1u1=Simplify[tgCoeffs[λ1,λ2][[2,2]]];
  C2m=C2op[λ1,λ2];
  rem=Simplify[c1u1-(-h^2(2h+θ1+θ2))C2m];
  mats={IdentityMatrix[d],X1[λ1,λ2],X2[λ1,λ2],Y1[λ1,λ2],Y2[λ1,λ2],C2m,C2m.C2m};
  names={"Id","X1","X2","Y1","Y2","C2","C2sq"};
  sol=Quiet@LinearSolve[Transpose[Flatten/@mats],Flatten[rem]];
  {"sol"->Thread[names->sol],
   "exact?"->(Simplify[rem-sol.mats]==ConstantArray[0,{d,d}])}]
```
**Output:**
```
Out= {"sol" -> {"Id" -> 2*(10*h^3 - 7*h^2*(θ1 + θ2) - 6*θ1*θ2*(θ1 + θ2) - 6*h*(θ1^2 + 4*θ1*θ2 + θ2^2)), "X1" -> 0, "X2" -> 0, "Y1" -> 0, "Y2" -> 0, "C2" -> 0, "C2sq" -> 0}, "exact?" -> True}
```
