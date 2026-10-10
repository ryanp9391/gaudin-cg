# AI Action Log — su3_V2_general

## 2026-07-17 15:22:35 — ✏️ EDIT CELL 2 [batch]
ClearAll[Ee,id,EE,idd];
Ee[λ_][i_,j_]:=JJ[i,j,λ];
id[λ_]:=IdentityMatrix[dim[λ]];
idd[λ1_,λ2_]:=KroneckerProduct[id[λ1],id[λ2]];
EE[λ1_,λ2_][i_,j_]:=KroneckerProduct[Ee[λ1][i,j],id[λ2]]+KroneckerProdu…


## 2026-07-17 15:22:35
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
Out= (no output)
```
## 2026-07-17 15:22:35 — ✏️ EDIT CELL 3 [batch]
(*central charge*)
ClearAll[C1];
C1[λ1_,λ2_]:=Sum[EE[λ1,λ2][i,i],{i,3}];


## 2026-07-17 15:22:36
**Input:**
```wolfram
(*central charge*)
ClearAll[C1];
C1[λ1_,λ2_]:=Sum[EE[λ1,λ2][i,i],{i,3}];
```
**Output:**
```
Out= (no output)
```
## 2026-07-17 15:22:36 — ✏️ EDIT CELL 5 [batch]
(*Lax operator for a general gl(3) rep λ={λ1,λ2,λ3} at the site*)
ClearAll[L];
L[λ_][i_,j_][u_]:=IdentityMatrix[dim[λ]] u KroneckerDelta[i,j] - h JJ[j,i,λ];


## 2026-07-17 15:22:36
**Input:**
```wolfram
(*Lax operator for a general gl(3) rep λ={λ1,λ2,λ3} at the site*)
ClearAll[L];
L[λ_][i_,j_][u_]:=IdentityMatrix[dim[λ]] u KroneckerDelta[i,j] - h JJ[j,i,λ];
```
**Output:**
```
Out= (no output)
```
## 2026-07-17 15:22:36 — ✏️ EDIT CELL 8 [batch]
(*monodromy with no twist, diagonal twist, and companion twist*)
ClearAll[T,Tg,TG];
T[λ1_,λ2_][i_,j_][u_]:=Sum[KroneckerProduct[L[λ1][i,ii][u-θ1],L[λ2][ii,j][u-θ2]],{ii,3}];
Tg[λ1_,λ2_][i_,j_][u_]:= z…


## 2026-07-17 15:22:36
**Input:**
```wolfram
(*monodromy with no twist, diagonal twist, and companion twist*)
ClearAll[T,Tg,TG];
T[λ1_,λ2_][i_,j_][u_]:=Sum[KroneckerProduct[L[λ1][i,ii][u-θ1],L[λ2][ii,j][u-θ2]],{ii,3}];
Tg[λ1_,λ2_][i_,j_][u_]:= z[i] T[λ1,λ2][i,j][u];
TG[λ1_,λ2_][i_,j_][u_]:=Sum[T[λ1,λ2][i,ii][u] G[[ii,j]],{ii,3}];
```
**Output:**
```
Out= (no output)
```
## 2026-07-17 15:22:36 — ✏️ EDIT CELL 9 [batch]
(*Yangian weights - eigenvalues of T[i,i] on the highest-weight state.
  General rep: T[i,i]|hw> = (u-θ1-h λ1[[i]])(u-θ2-h λ2[[i]]) |hw>, from L[λ][i,i]|hw>=(u-h λ_i)|hw>.
  Symmetric case λ=(S,0,0) r…


## 2026-07-17 15:22:36
**Input:**
```wolfram
(*Yangian weights - eigenvalues of T[i,i] on the highest-weight state.
  General rep: T[i,i]|hw> = (u-θ1-h λ1[[i]])(u-θ2-h λ2[[i]]) |hw>, from L[λ][i,i]|hw>=(u-h λ_i)|hw>.
  Symmetric case λ=(S,0,0) reduces to ν[1]=(u-θ1-hS1)(u-θ2-hS2), ν[2]=ν[3]=(u-θ1)(u-θ2).*)
ClearAll[ν,HWS];
ν[λ1_,λ2_][i_][u_]:=(u-θ1-h λ1[[i]])(u-θ2-h λ2[[i]]);
HWS[λ1_,λ2_]:=Table[{KroneckerDelta[1,ii]},{ii,1,dim[λ1] dim[λ2]}];
```
**Output:**
```
Out= (no output)
```
## 2026-07-17 15:22:49 — ✏️ EDIT CELL 4 [batch]
(*checking central charge constraint = total boxes, on general reps*)
Module[{reps={{{2,1,0},{1,0,0}},{{1,1,0},{1,1,0}},{{2,2,0},{1,0,0}},{{2,1,0},{1,1,0}}}},
 Union@Flatten@Table[
   C1[r[[1]],r[[2]]…


## 2026-07-17 15:22:49
**Input:**
```wolfram
(*checking central charge constraint = total boxes, on general reps*)
Module[{reps={{{2,1,0},{1,0,0}},{{1,1,0},{1,1,0}},{{2,2,0},{1,0,0}},{{2,1,0},{1,1,0}}}},
 Union@Flatten@Table[
   C1[r[[1]],r[[2]]]-(Total[r[[1]]]+Total[r[[2]]]) idd[r[[1]],r[[2]]],{r,reps}]]
```
**Output:**
```
Out= Out[14]= \left\{0\right\}
```
## 2026-07-17 15:22:49 — ✏️ EDIT CELL 10 [batch]
(*checking highest-weights against the general ν, on non-symmetric reps*)
Module[{reps={{{2,1,0},{1,0,0}},{{1,1,0},{1,1,0}},{{2,2,0},{1,0,0}},{{2,1,0},{1,1,0}}}},
 Union@Flatten@Table[
   T[r[[1]],r[[…


## 2026-07-17 15:22:49
**Input:**
```wolfram
(*checking highest-weights against the general ν, on non-symmetric reps*)
Module[{reps={{{2,1,0},{1,0,0}},{{1,1,0},{1,1,0}},{{2,2,0},{1,0,0}},{{2,1,0},{1,1,0}}}},
 Union@Flatten@Table[
   T[r[[1]],r[[2]]][i,i][u].HWS[r[[1]],r[[2]]]-ν[r[[1]],r[[2]]][i][u] HWS[r[[1]],r[[2]]],
   {r,reps},{i,3}]]
```
**Output:**
```
Out= Out[37]= \left\{0\right\}
```
## 2026-07-17 15:22:59 — ▶️ RUN CELLS 2–10
Cell 2: ✓ — (no output)
Cell 3: ✓ — (no output)
Cell 4: ✓ — Out[14]= \left\{0\right\}
Cell 5: ✓ — (no output)
Cell 6: ✓ — (no output)
Cell 7: ✓ — Out[27]= \left\{\frac{1}{3},\frac{1}{7},1,0.360439875147226537266677650596+0.932782448593374716484835982071\,\mathrm{i},0.508867380990499828163321372267+0.860844927129079093188644393616\,\mathrm{i},0.533257880589599372595832345093+0.845952736734794526483381618453\,\mathrm{i},1.40256513672732573802583136796+2.63958011245724833615686199414\,\mathrm{i},-1.6733239859016290991182472189+2.47680418504351… [511 chars total]
Cell 8: ✓ — (no output)
Cell 9: ✓ — (no output)
Cell 10: ✓ — Out[37]= \left\{0\right\}

## 2026-07-17 15:29:34 — ✏️ EDIT CELL 11 [batch]
(*quantum minors (q-determinants) of the monodromy matrix*)
ClearAll[qmin,qming,qminG];
qmin[λ1_, λ2_][II_List, JJ_List][u_] := 
    qmin[λ1, λ2][II, JJ][u] = 
     Sum[Signature[σ] 
        Dot @@ Ta…


## 2026-07-17 15:29:35
**Input:**
```wolfram
(*quantum minors (q-determinants) of the monodromy matrix*)
ClearAll[qmin,qming,qminG];
qmin[λ1_, λ2_][II_List, JJ_List][u_] := 
    qmin[λ1, λ2][II, JJ][u] = 
     Sum[Signature[σ] 
        Dot @@ Table[T[λ1, λ2][II[[σ[[k]]]], JJ[[k]]][u - (k - 1) h], {k, 
            Length[II]}], {σ, Permutations[Range[Length[II]]]}];

qming[λ1_, λ2_][II_List, JJ_List][u_] := 
    qming[λ1, λ2][II, JJ][u] = 
     Sum[Signature[σ] 
        Dot @@ Table[Tg[λ1, λ2][II[[σ[[k]]]], JJ[[k]]][u - (k - 1) h], {k, 
            Length[II]}], {σ, Permutations[Range[Length[II]]]}];

qminG[λ1_, λ2_][II_List, JJ_List][u_] := 
    qminG[λ1, λ2][II, JJ][u] = 
     Sum[Signature[σ] 
        Dot @@ Table[TG[λ1, λ2][II[[σ[[k]]]], JJ[[k]]][u - (k - 1) h], {k, 
            Length[II]}], {σ, Permutations[Range[Length[II]]]}];
```
**Output:**
```
Out= (no output)
```
## 2026-07-17 15:29:35 — ✏️ EDIT CELL 12 [batch]
(*transfer matrices with diagonal twist*)
ClearAll[tg];
tg[λ1_,λ2_][a_,1][u_]:= Total[qming[λ1, λ2][#, #][u] & /@ Subsets[{1, 2, 3}, {a}]]

(*transfer matrices with companion twist*)
ClearAll[tG];
tG[…


## 2026-07-17 15:29:35
**Input:**
```wolfram
(*transfer matrices with diagonal twist*)
ClearAll[tg];
tg[λ1_,λ2_][a_,1][u_]:= Total[qming[λ1, λ2][#, #][u] & /@ Subsets[{1, 2, 3}, {a}]]

(*transfer matrices with companion twist*)
ClearAll[tG];
tG[λ1_,λ2_][a_,1][u_]:= Total[qminG[λ1, λ2][#, #][u] & /@ Subsets[{1, 2, 3}, {a}]]
```
**Output:**
```
Out= (no output)
```
## 2026-07-17 15:29:35 — ✏️ EDIT CELL 13 [batch]
(*value of quantum determinant (χ3 carried here; general ν)*)
ClearAll[qdetT];
qdetT[λ1_,λ2_][u_]:=χ3 ν[λ1,λ2][1][u-2h] ν[λ1,λ2][2][u-h] ν[λ1,λ2][3][u];


## 2026-07-17 15:29:35
**Input:**
```wolfram
(*value of quantum determinant (χ3 carried here; general ν)*)
ClearAll[qdetT];
qdetT[λ1_,λ2_][u_]:=χ3 ν[λ1,λ2][1][u-2h] ν[λ1,λ2][2][u-h] ν[λ1,λ2][3][u];
```
**Output:**
```
Out= (no output)
```
## 2026-07-17 15:29:35 — ✏️ EDIT CELL 16 [batch]
(*single source of truth for the Yangian weight of T11*)
ClearAll[ν1];
ν1[λ1_,λ2_][u_]:=ν[λ1,λ2][1][u];
ν1[{1,0,0},{1,0,0}][u]


## 2026-07-17 15:29:35
**Input:**
```wolfram
(*single source of truth for the Yangian weight of T11*)
ClearAll[ν1];
ν1[λ1_,λ2_][u_]:=ν[λ1,λ2][1][u];
ν1[{1,0,0},{1,0,0}][u]
```
**Output:**
```
Out= Out[52]= \left(u-\frac{4}{3}\right)\,\left(u-\frac{8}{7}\right)
```
## 2026-07-17 15:29:35 — ✏️ EDIT CELL 18 [batch]
(*3rd-order Baxter/TQ equation (reference form). Diagonal-twist transfer matrices τ1=tg[·][1,1],
  τ2=tg[·][2,1]; qdetT already carries χ3, so NO explicit χ3 here (would double-count).
  A1=τ1/ν1[u], …


## 2026-07-17 15:29:35
**Input:**
```wolfram
(*3rd-order Baxter/TQ equation (reference form). Diagonal-twist transfer matrices τ1=tg[·][1,1],
  τ2=tg[·][2,1]; qdetT already carries χ3, so NO explicit χ3 here (would double-count).
  A1=τ1/ν1[u], A2=τ2/(ν1[u]ν1[u-h]), A3=qdetT/(ν1[u]ν1[u-h]ν1[u-2h]).*)
ClearAll[Bax];
Bax[λ1_, λ2_] :=
  f[u + h]
   - tg[λ1, λ2][1, 1][u]/ν1[λ1, λ2][u] f[u]
   + tg[λ1, λ2][2, 1][u]/(ν1[λ1, λ2][u] ν1[λ1, λ2][u - h]) f[u - h]
   - qdetT[λ1, λ2][u]/(ν1[λ1, λ2][u] ν1[λ1, λ2][u - h] ν1[λ1, λ2][u - 2 h]) f[u - 2 h];
Head[Bax[{1,0,0}, {1,0,0}]]
```
**Output:**
```
Out= Out[57]= \mathrm{List}
```
## 2026-07-17 15:29:48 — ✏️ EDIT CELL 14
(*qdetT central value == tg[·][3,1]; on non-symmetric reps this is the ν shift-ordering guard*)
Module[{reps={{{2,1,0},{1,0,0}},{{1,1,0},{1,1,0}},{{2,2,0},{1,0,0}},{{2,1,0},{1,1,0}}}},
 Union@Flatten@…


## 2026-07-17 15:29:48
**Input:**
```wolfram
(*qdetT central value == tg[·][3,1]; on non-symmetric reps this is the ν shift-ordering guard*)
Module[{reps={{{2,1,0},{1,0,0}},{{1,1,0},{1,1,0}},{{2,2,0},{1,0,0}},{{2,1,0},{1,1,0}}}},
 Union@Flatten@Table[
   Simplify[qdetT[r[[1]],r[[2]]][u] idd[r[[1]],r[[2]]]-tg[r[[1]],r[[2]]][3,1][u]],{r,reps}]]
```
**Output:**
```
Out= {0, (0``28.792425882065558 + 0``28.79242588206556*I)*(-10/3 + u)*(-22/7 + u)*(-1/3 + u)*(-1/7 + u), (0``29.143855764391372 + 0``29.143855764391372*I)*(-13/3 + u)*(-22/7 + u)*(-7/3 + u)*(-15/7 + u)*(-1/3 + u)*(-1/7 + u), (0``29.143855764391372 + 0``29.143855764391372*I)*(-10/3 + u)*(-22/7 + u)*(-7/3 + u)*(-15/7 + u)*(-1/3 + u)*(-1/7 + u), (0``29.143855764391372 + 0``29.143855764391372*I)*(-13/3 + u)*(-10/3 + u)*(-22/7 + u)*(-8/7 + u)*(-1/3 + u)*(-1/7 + u), (0``29.143855764391372 + 0``29.143855764391372*I)*(-13/3 + u)*(-22/7 + u)*(-7/3 + u)*(-8/7 + u)*(-1/3 + u)*(-1/7 + u), (0``29.368197005148932 + 0``29.36819700514893*I) + (0``27.34995413925883 + 0``27.349954139258827*I)*u + (0``27.18085934747624 + 0``27.18085934747624*I)*u^2 + (0``27.267897180569527 + 0``27.26789718056952*I)*u^3, (0``29.278122895797562 + 0``29.278122895797555*I) + (0``27.34995413925883 + 0``27.349954139258827*I)*u + (0``27.18085934747624 + 0``27.18085934747624*I)*u^2 + (0``27.267897180569527 + 0``27.26789718056952*I)*u^3, (0``28.506940455091872 + 0``28.506940455091872*I) + (0``27.34995413925883 + 0``27.349954139258827*I)*u + (0``27.18085934747624 + 0``27.18085934747624*I)*u^2 + (0``27.267897180569527 + 0``27.26789718056952*I)*u^3, (0``28.497862220919053 + 0``28.497862220919053*I) + (0``27.34995413925883 + 0``27.349954139258827*I)*u + (0``27.18085934747624 + 0``27.18085934747624*I)*u^2 + (0``27.267897180569527 + 0``27.26789718056952*I)*u^3, (0``28.49723127761503 + 0``28.49723127761502*I) + (0``27.34995413925883 + 0``27.349954139258827*I)*u + (0``27.18085934747624 + 0``27.18085934747624*I)*u^2 + (0``27.267897180569527 + 0``27.26789718056952*I)*u^3, (0``28.3654352186077 + 0``28.3654352186077*I) + (0``27.34995413925883 + 0``27.349954139258827*I)*u + (0``27.18085934747624 + 0``27.18085934747624*I)*u^2 + (0``27.267897180569527 + 0``27.26789718056952*I)*u^3, (0``27.629847720265417 + 0``27.62984772026541*I) + (0``27.34995413925883 + 0``27.349954139258827*I)*u + (0``27.18085934747624 + 0``27.18085934747624*I)*u^2 + (0``27.267897180569527 + 0``27.26789718056952*I)*u^3, (0``27.627691509041714 + 0``27.627691509041714*I) + (0``27.34995413925883 + 0``27.349954139258827*I)*u + (0``27.18085934747624 + 0``27.18085934747624*I)*u^2 + (0``27.267897180569527 + 0``27.26789718056952*I)*u^3, (0``30.010895739261123 + 0``30.010895739261112*I) + (0``27.34995413925883 + 0``27.349954139258827*I)*u + (0``27.18085934747624 + 0``27.18085934747624*I)*u^2 + (0``27.267897180569527 + 0``27.26789718056952*I)*u^3 + (0``27.58380203700932 + 0``27.583802037009317*I)*u^4, (0``29.375044513638386 + 0``29.375044513638386*I) + (0``27.34995413925883 + 0``27.349954139258827*I)*u + (0``27.18085934747624 + 0``27.18085934747624*I)*u^2 + (0``27.267897180569527 + 0``27.26789718056952*I)*u^3 + (0``27.58380203700932 + 0``27.583802037009317*I)*u^4, (0``29.327863279905728 + 0``29.327863279905728*I) + (0``27.34995413925883 + 0``27.349954139258827*I)*u + (0``27.18085934747624 + 0``27.18085934747624*I)*u^2 + (0``27.267897180569527 + 0``27.26789718056952*I)*u^3 + (0``27.58380203700932 + 0``27.583802037009317*I)*u^4, (0``29.323941146057887 + 0``29.32394114605788*I) + (0``27.34995413925883 + 0``27.349954139258827*I)*u + (0``27.18085934747624 + 0``27.18085934747624*I)*u^2 + (0``27.267897180569527 + 0``27.26789718056952*I)*u^3 + (0``27.58380203700932 + 0``27.583802037009317*I)*u^4, (0``29.29214450880943 + 0``29.292144508809425*I) + (0``27.34995413925883 + 0``27.349954139258827*I)*u + (0``27.18085934747624 + 0``27.18085934747624*I)*u^2 + (0``27.267897180569527 + 0``27.26789718056952*I)*u^3 + (0``27.58380203700932 + 0``27.583802037009317*I)*u^4, (0``29.191531884612953 + 0``29.19153188461295*I) + (0``27.34995413925883 + 0``27.349954139258827*I)*u + (0``27.18085934747624 + 0``27.18085934747624*I)*u^2 + (0``27.267897180569527 + 0``27.26789718056952*I)*u^3 + (0``27.58380203700932 + 0``27.583802037009317*I)*u^4, (0``29.18902507100624 + 0``29.18902507100624*I) + (0``27.34995413925883 + 0``27.349954139258827*I)*u + (0``27.18085934747624 + 0``27.18085934747624*I)*u^2 + (0``27.267897180569527 + 0``27.26789718056952*I)*u^3 + (0``
[output truncated — 55207 chars total]
```

## 2026-07-17 15:30:29
**Input:**
```wolfram
reps={{{2,1,0},{1,0,0}},{{1,1,0},{1,1,0}},{{2,2,0},{1,0,0}},{{2,1,0},{1,1,0}}};
Table[Chop[N[qdetT[r[[1]],r[[2]]][7/2] idd[r[[1]],r[[2]]]-tg[r[[1]],r[[2]]][3,1][7/2]],10^-15]//Max@*Abs@*Flatten,{r,reps}]
```
**Output:**
```
Out[1]= (suppressed)
Out[2]= {0, 0, 0, 0}
```

## 2026-07-17 15:30:42
**Input:**
```wolfram
reps={{{2,1,0},{1,0,0}},{{1,1,0},{1,1,0}},{{2,2,0},{1,0,0}},{{2,1,0},{1,1,0}}};
raw = Table[Simplify[qdetT[r[[1]],r[[2]]][u] idd[r[[1]],r[[2]]]-tg[r[[1]],r[[2]]][3,1][u]],{r,reps}];
Union@Flatten@Map[Chop[#,10^-15]&,raw]
```
**Output:**
```
Out[1]= (suppressed)
Out[2]= (suppressed)
Out[3]= {0}
```

## 2026-07-17 15:31:00
**Input:**
```wolfram
reps={{{2,1,0},{1,0,0}},{{1,1,0},{1,1,0}},{{2,2,0},{1,0,0}},{{2,1,0},{1,1,0}}};
(* test reversed ordering nu[3][u-2h] nu[2][u-h] nu[1][u] as a control, should NOT vanish if default ordering is the unique correct one *)
testRev = Table[Chop[N[χ3 ν[r[[1]],r[[2]]][3][7/2-2h] ν[r[[1]],r[[2]]][2][7/2-h] ν[r[[1]],r[[2]]][1][7/2] idd[r[[1]],r[[2]]] - tg[r[[1]],r[[2]]][3,1][7/2]],10^-15]//Max@*Abs@*Flatten,{r,reps}]
```
**Output:**
```
Out[1]= (suppressed)
Out[2]= {20.674886621315192, 11.801398337112621, 2.953555231616456, 11.903722600151172}
```
## 2026-07-17 15:31:09 — ✏️ EDIT CELL 14
(*qdetT central value == tg[·][3,1]; on non-symmetric reps this is the ν shift-ordering guard.
  NB: Chop added around Simplify — z[k] (hence χ3) are high-precision numeric (not exact rationals),
  so…


## 2026-07-17 15:31:09
**Input:**
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
**Output:**
```
Out= {0}
```
## 2026-07-17 15:31:17 — ✏️ EDIT CELL 17
(* NB: exploratory, not used by alphaCoeffs — ee[λ1,λ2][i,i] was never defined even in V1 (stayed
   symbolic), so this cell has always been a scratch expression, not a validated formula. *)
tasympt[1…


## 2026-07-17 15:31:17
**Input:**
```wolfram
(* NB: exploratory, not used by alphaCoeffs — ee[λ1,λ2][i,i] was never defined even in V1 (stayed
   symbolic), so this cell has always been a scratch expression, not a validated formula. *)
tasympt[1,1]=u^2 χ1 - 
         u (θ1 + θ2) χ1 - 
         h u Sum[z[ii] ee[λ1, λ2][ii, ii],{ii,3}]
```
**Output:**
```
Out= (-0.6678881603463455895361101752170204793164128942279036824046`29.820316782640784 - 1.256942910693927779122315235304845126254785777906856664183`30.094928590684432*I)*u + (1.4025651367273257380258313679557430065644670778785977330496`29.820316782640784 + 2.639580112457248336156861994140174765135195253625596338774`30.094928590684432*I)*u^2 - u*((0.3604398751472265372666776505955500912937645232622850024938`29.707347828414424 + 0.9327824485933747164848359820708679294578496404879183613065`30.120295363558725*I)*ee[λ1, λ2][1, 1] + (0.5088673809904998281633213722670686199090953963505721183772`29.857119610800275 + 0.8608449271290790931886443936161642458151438806992171086196`30.085439922389305*I)*ee[λ1, λ2][2, 2] + (0.5332578805895993725958323450931242953616071582657406121786`29.877452280089955 + 0.8459527367347945264833816184531425898622017324384608688479`30.077861097574722*I)*ee[λ1, λ2][3, 3])
```
## 2026-07-17 15:31:38 — ▶️ RUN CELL 14
**In [14]:** `(*qdetT central value == tg[·][3,1]; on non-symmetric reps this is the ν shift-ordering guard.↵  NB: Chop added around Simplify — z[k] (hence χ3) are high-precision numeric (not exact rationals),↵  so`
**Out:** Out[48]= \left\{0,\left(0.\times 10^{-29}+0.\times 10^{-29}\,\mathrm{i}\right)\,\left(u-\frac{36}{7}\right)\,\left(u-\frac{10}{3}\right)\,\left(u-\frac{1}{3}\right)\,\left(u-\frac{1}{7}\right),\left(0.\times 10^{-29}+0.\times 10^{-29}\,\mathrm{i}\right)\,\left(u-\frac{29}{7}\right)\,\left(u-\frac{10

## 2026-07-17 15:31:38 — ▶️ RUN CELL 17
**In [17]:** `(* NB: exploratory, not used by alphaCoeffs — ee[λ1,λ2][i,i] was never defined even in V1 (stayed↵   symbolic), so this cell has always been a scratch expression, not a validated formula. *)↵tasympt[1`
**Out:** Out[124]= -h\,u\,\left(z\left(1\right)\,\mathrm{ee}\left(S1,S2\right)\left(1,1\right)+z\left(2\right)\,\mathrm{ee}\left(S1,S2\right)\left(2,2\right)+z\left(3\right)\,\mathrm{ee}\left(S1,S2\right)\left(3,3\right)\right)+u^2\,\left(z\left(1\right)+z\left(2\right)+z\left(3\right)\right)-u\,\left(z\left


## 2026-07-17 15:32:24
**Input:**
```wolfram
$Line
```
**Output:**
```
Out= 175
```
## 2026-07-17 15:32:59 — 🗑️ DELETE CODE CELL 14
(*qdetT central value == tg[·][3,1]; on non-symmetric reps this is the ν shift-ordering guard.
  NB: Chop added around Simplify — z[k] (hence χ3) are high-precision numeric (not exact rationals),
  so…

## 2026-07-17 15:33:11 — 📥 BULK INSERT 1 CELL(S) at positions 14–14
1. [code] (*qdetT central value == tg[·][3,1]; on non-symmetric reps this is the ν shift-ordering guard.↵  NB:


## 2026-07-17 15:33:11
**Input:**
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
**Output:**
```
(no output)
```
## 2026-07-17 15:33:36 — ▶️ RUN CELL 14
**In [14]:** `(*qdetT central value == tg[·][3,1]; on non-symmetric reps this is the ν shift-ordering guard.↵  NB: Chop added around Simplify — z[k] (hence χ3) are high-precision numeric (not exact rationals),↵  so`
**Out:** (no output)

## 2026-07-17 15:34:17 — ▶️ RUN CELL 14
**In [14]:** `(*qdetT central value == tg[·][3,1]; on non-symmetric reps this is the ν shift-ordering guard.↵  NB: Chop added around Simplify — z[k] (hence χ3) are high-precision numeric (not exact rationals),↵  so`
**Out:** (no output)

## 2026-07-17 15:34:33 — ✏️ EDIT CELL 14
(*qdetT central value == tg[·][3,1]; on non-symmetric reps this is the ν shift-ordering guard.
  NB: Chop added around Simplify — z[k] (hence χ3) are high-precision numeric (not exact rationals),
  so…


## 2026-07-17 15:34:33
**Input:**
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
**Output:**
```
Out= {0}
```
## 2026-07-17 17:40:23 — ✏️ EDIT CELL 14
(*qdetT central value == tg[·][3,1]; on non-symmetric reps this is the ν shift-ordering guard*)
Module[{reps={{{2,1,0},{1,0,0}},{{1,1,0},{1,1,0}},{{2,2,0},{1,0,0}},{{2,1,0},{1,1,0}}}},
 Union@Flatten@…

## 2026-07-17 17:40:27 — ▶️ RUN CELL 11
**In [11]:** `(*quantum minors (q-determinants) of the monodromy matrix*)↵ClearAll[qmin,qming,qminG];↵qmin[λ1_, λ2_][II_List, JJ_List][u_] := ↵    qmin[λ1, λ2][II, JJ][u] = ↵     Sum[Signature[σ] ↵        Dot @@ Ta`
**Out:** (no output)

## 2026-07-17 17:40:27 — ▶️ RUN CELL 12
**In [12]:** `(*transfer matrices with diagonal twist*)↵ClearAll[tg];↵tg[λ1_,λ2_][a_,1][u_]:= Total[qming[λ1, λ2][#, #][u] & /@ Subsets[{1, 2, 3}, {a}]]↵↵(*transfer matrices with companion twist*)↵ClearAll[tG];↵tG[`
**Out:** (no output)

## 2026-07-17 17:40:28 — ▶️ RUN CELL 13
**In [13]:** `(*value of quantum determinant (χ3 carried here; general ν)*)↵ClearAll[qdetT];↵qdetT[λ1_,λ2_][u_]:=χ3 ν[λ1,λ2][1][u-2h] ν[λ1,λ2][2][u-h] ν[λ1,λ2][3][u];`
**Out:** (no output)

## 2026-07-17 17:40:30 — ▶️ RUN CELL 16
**In [16]:** `(*single source of truth for the Yangian weight of T11*)↵ClearAll[ν1];↵ν1[λ1_,λ2_][u_]:=ν[λ1,λ2][1][u];↵ν1[{1,0,0},{1,0,0}][u]`
**Out:** Out[52]= \left(u-\frac{4}{3}\right)\,\left(u-\frac{8}{7}\right)

## 2026-07-17 17:40:31 — ▶️ RUN CELL 18
**In [18]:** `(*3rd-order Baxter/TQ equation (reference form). Diagonal-twist transfer matrices τ1=tg[·][1,1],↵  τ2=tg[·][2,1]; qdetT already carries χ3, so NO explicit χ3 here (would double-count).↵  A1=τ1/ν1[u], `
**Out:** Out[57]= \mathrm{List}

## 2026-07-17 17:40:34 — ▶️ RUN CELL 14
**In [14]:** `(*qdetT central value == tg[·][3,1]; on non-symmetric reps this is the ν shift-ordering guard*)↵Module[{reps={{{2,1,0},{1,0,0}},{{1,1,0},{1,1,0}},{{2,2,0},{1,0,0}},{{2,1,0},{1,1,0}}}},↵ Union@Flatten@`
**Out:** (no output)


## 2026-07-17 17:40:41
**Input:**
```wolfram
Module[{reps={{{2,1,0},{1,0,0}},{{1,1,0},{1,1,0}},{{2,2,0},{1,0,0}},{{2,1,0},{1,1,0}}}},
 Union@Flatten@Table[
   Simplify[qdetT[r[[1]],r[[2]]][u] idd[r[[1]],r[[2]]]-tg[r[[1]],r[[2]]][3,1][u]],{r,reps}]]
```
**Output:**
```
Out= {0, (0``28.792425882065558 + 0``28.79242588206556*I)*(-10/3 + u)*(-22/7 + u)*(-1/3 + u)*(-1/7 + u), (0``29.143855764391372 + 0``29.143855764391372*I)*(-13/3 + u)*(-22/7 + u)*(-7/3 + u)*(-15/7 + u)*(-1/3 + u)*(-1/7 + u), (0``29.143855764391372 + 0``29.143855764391372*I)*(-10/3 + u)*(-22/7 + u)*(-7/3 + u)*(-15/7 + u)*(-1/3 + u)*(-1/7 + u), (0``29.143855764391372 + 0``29.143855764391372*I)*(-13/3 + u)*(-10/3 + u)*(-22/7 + u)*(-8/7 + u)*(-1/3 + u)*(-1/7 + u), (0``29.143855764391372 + 0``29.143855764391372*I)*(-13/3 + u)*(-22/7 + u)*(-7/3 + u)*(-8/7 + u)*(-1/3 + u)*(-1/7 + u), (0``29.368197005148932 + 0``29.36819700514893*I) + (0``27.34995413925883 + 0``27.349954139258827*I)*u + (0``27.18085934747624 + 0``27.18085934747624*I)*u^2 + (0``27.267897180569527 + 0``27.26789718056952*I)*u^3, (0``29.278122895797562 + 0``29.278122895797555*I) + (0``27.34995413925883 + 0``27.349954139258827*I)*u + (0``27.18085934747624 + 0``27.18085934747624*I)*u^2 + (0``27.267897180569527 + 0``27.26789718056952*I)*u^3, (0``28.506940455091872 + 0``28.506940455091872*I) + (0``27.34995413925883 + 0``27.349954139258827*I)*u + (0``27.18085934747624 + 0``27.18085934747624*I)*u^2 + (0``27.267897180569527 + 0``27.26789718056952*I)*u^3, (0``28.497862220919053 + 0``28.497862220919053*I) + (0``27.34995413925883 + 0``27.349954139258827*I)*u + (0``27.18085934747624 + 0``27.18085934747624*I)*u^2 + (0``27.267897180569527 + 0``27.26789718056952*I)*u^3, (0``28.49723127761503 + 0``28.49723127761502*I) + (0``27.34995413925883 + 0``27.349954139258827*I)*u + (0``27.18085934747624 + 0``27.18085934747624*I)*u^2 + (0``27.267897180569527 + 0``27.26789718056952*I)*u^3, (0``28.3654352186077 + 0``28.3654352186077*I) + (0``27.34995413925883 + 0``27.349954139258827*I)*u + (0``27.18085934747624 + 0``27.18085934747624*I)*u^2 + (0``27.267897180569527 + 0``27.26789718056952*I)*u^3, (0``27.629847720265417 + 0``27.62984772026541*I) + (0``27.34995413925883 + 0``27.349954139258827*I)*u + (0``27.18085934747624 + 0``27.18085934747624*I)*u^2 + (0``27.267897180569527 + 0``27.26789718056952*I)*u^3, (0``27.627691509041714 + 0``27.627691509041714*I) + (0``27.34995413925883 + 0``27.349954139258827*I)*u + (0``27.18085934747624 + 0``27.18085934747624*I)*u^2 + (0``27.267897180569527 + 0``27.26789718056952*I)*u^3, (0``30.010895739261123 + 0``30.010895739261112*I) + (0``27.34995413925883 + 0``27.349954139258827*I)*u + (0``27.18085934747624 + 0``27.18085934747624*I)*u^2 + (0``27.267897180569527 + 0``27.26789718056952*I)*u^3 + (0``27.58380203700932 + 0``27.583802037009317*I)*u^4, (0``29.375044513638386 + 0``29.375044513638386*I) + (0``27.34995413925883 + 0``27.349954139258827*I)*u + (0``27.18085934747624 + 0``27.18085934747624*I)*u^2 + (0``27.267897180569527 + 0``27.26789718056952*I)*u^3 + (0``27.58380203700932 + 0``27.583802037009317*I)*u^4, (0``29.327863279905728 + 0``29.327863279905728*I) + (0``27.34995413925883 + 0``27.349954139258827*I)*u + (0``27.18085934747624 + 0``27.18085934747624*I)*u^2 + (0``27.267897180569527 + 0``27.26789718056952*I)*u^3 + (0``27.58380203700932 + 0``27.583802037009317*I)*u^4, (0``29.323941146057887 + 0``29.32394114605788*I) + (0``27.34995413925883 + 0``27.349954139258827*I)*u + (0``27.18085934747624 + 0``27.18085934747624*I)*u^2 + (0``27.267897180569527 + 0``27.26789718056952*I)*u^3 + (0``27.58380203700932 + 0``27.583802037009317*I)*u^4, (0``29.29214450880943 + 0``29.292144508809425*I) + (0``27.34995413925883 + 0``27.349954139258827*I)*u + (0``27.18085934747624 + 0``27.18085934747624*I)*u^2 + (0``27.267897180569527 + 0``27.26789718056952*I)*u^3 + (0``27.58380203700932 + 0``27.583802037009317*I)*u^4, (0``29.191531884612953 + 0``29.19153188461295*I) + (0``27.34995413925883 + 0``27.349954139258827*I)*u + (0``27.18085934747624 + 0``27.18085934747624*I)*u^2 + (0``27.267897180569527 + 0``27.26789718056952*I)*u^3 + (0``27.58380203700932 + 0``27.583802037009317*I)*u^4, (0``29.18902507100624 + 0``29.18902507100624*I) + (0``27.34995413925883 + 0``27.349954139258827*I)*u + (0``27.18085934747624 + 0``27.18085934747624*I)*u^2 + (0``27.267897180569527 + 0``27.26789718056952*I)*u^3 + (0``
[output truncated — 55207 chars total]
```

## 2026-07-17 17:41:07
**Input:**
```wolfram
(* Confirm the raw (non-Chopped) verbatim-check residual is numerical zero at high precision, not a real discrepancy: *)
verbatimCheck = Module[{reps={{{2,1,0},{1,0,0}},{{1,1,0},{1,1,0}},{{2,2,0},{1,0,0}},{{2,1,0},{1,1,0}}}},
 Union@Flatten@Table[
   Simplify[qdetT[r[[1]],r[[2]]][u] idd[r[[1]],r[[2]]]-tg[r[[1]],r[[2]]][3,1][u]],{r,reps}]];
Length[verbatimCheck]
Chop[verbatimCheck /. u->RandomReal[{-5,5}]]
Chop[verbatimCheck /. u -> 0]
(* magnitude sanity: how small are these "zero" coefficients really, at u=1? *)
Max[Abs[N[verbatimCheck /. u->1]]]
```
**Output:**
```
Out= {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}
```

## 2026-07-17 17:41:11
**Input:**
```wolfram
Max[Abs[N[verbatimCheck /. u->1]]]
```
**Output:**
```
Out= 0.
```

## 2026-07-17 17:41:14
**Input:**
```wolfram
(* sanity: test that a DIFFERENT (wrong) shift ordering does NOT vanish, confirming default u-2h,u-h,u ordering is the unique correct one *)
qdetTrev[λ1_,λ2_][u_]:=χ3 ν[λ1,λ2][3][u-2h] ν[λ1,λ2][2][u-h] ν[λ1,λ2][1][u];
revCheck = Module[{reps={{{2,1,0},{1,0,0}},{{1,1,0},{1,1,0}},{{2,2,0},{1,0,0}},{{2,1,0},{1,1,0}}}},
 Union@Flatten@Table[
   N[qdetTrev[r[[1]],r[[2]]][1] idd[r[[1]],r[[2]]]-tg[r[[1]],r[[2]]][3,1][1]],{r,reps}]];
Max[Abs[revCheck]]
```
**Output:**
```
Out= 5.777345858978512
```
## 2026-07-17 17:41:39 — ✏️ EDIT CELL 14
(*qdetT central value == tg[·][3,1]; on non-symmetric reps this is the ν shift-ordering guard.
  NB: Chop wraps Simplify because z[k] (hence χ3) are high-precision *numeric* (not exact
  symbolic) her…


## 2026-07-17 17:41:39
**Input:**
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
**Output:**
```
Out= {0}
```
## 2026-07-17 17:41:52 — ▶️ RUN CELL 14
**In [14]:** `(*qdetT central value == tg[·][3,1]; on non-symmetric reps this is the ν shift-ordering guard.↵  NB: Chop wraps Simplify because z[k] (hence χ3) are high-precision *numeric* (not exact↵  symbolic) her`
**Out:** (no output)

## 2026-07-17 17:42:09 — 🗑️ DELETE CODE CELL 14
(*qdetT central value == tg[·][3,1]; on non-symmetric reps this is the ν shift-ordering guard.
  NB: Chop wraps Simplify because z[k] (hence χ3) are high-precision *numeric* (not exact
  symbolic) her…

## 2026-07-17 17:42:16 — 📥 BULK INSERT 1 CELL(S) at positions 14–14
1. [code] (*qdetT central value == tg[·][3,1]; on non-symmetric reps this is the ν shift-ordering guard.↵  NB:


## 2026-07-17 17:42:16
**Input:**
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
**Output:**
```
(no output)
```
## 2026-07-17 17:42:33 — ▶️ RUN CELL 14
**In [14]:** `(*qdetT central value == tg[·][3,1]; on non-symmetric reps this is the ν shift-ordering guard.↵  NB: Chop wraps Simplify because z[k] (hence χ3) are high-precision *numeric* (not exact↵  symbolic) her`
**Out:** (no output)

## 2026-07-17 17:42:39 — ▶️ RUN CELL 4
**In [4]:** `(*checking central charge constraint = total boxes, on general reps*)↵Module[{reps={{{2,1,0},{1,0,0}},{{1,1,0},{1,1,0}},{{2,2,0},{1,0,0}},{{2,1,0},{1,1,0}}}},↵ Union@Flatten@Table[↵   C1[r[[1]],r[[2]]`
**Out:** Out[14]= \left\{0\right\}


## 2026-07-17 17:42:45
**Input:**
```wolfram
$MessageList
```
**Output:**
```
Out= {}
```
## 2026-07-17 17:42:46 — 📥 BULK INSERT 1 CELL(S) at positions 15–15
1. [code] (*scratch: does a trivial {0} result render via runCell?*)↵{0}


## 2026-07-17 17:42:47
**Input:**
```wolfram
(*scratch: does a trivial {0} result render via runCell?*)
{0}
```
**Output:**
```
(no output)
```
## 2026-07-17 17:43:02 — ✏️ EDIT CELL 15
{1,2,3}


## 2026-07-17 17:43:02
**Input:**
```wolfram
{1,2,3}
```
**Output:**
```
Out= {1, 2, 3}
```
## 2026-07-17 17:43:05 — ✏️ EDIT CELL 15
{0}


## 2026-07-17 17:43:05
**Input:**
```wolfram
{0}
```
**Output:**
```
Out= {0}
```
## 2026-07-17 17:43:17 — ✏️ EDIT CELL 14
(*qdetT central value == tg[·][3,1]; on non-symmetric reps this is the ν shift-ordering guard.
  NB: Chop wraps Simplify because z[k] (hence χ3) are high-precision *numeric* (not exact
  symbolic) her…


## 2026-07-17 17:43:17
**Input:**
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
**Output:**
```
Out= {0}
```
## 2026-07-17 17:43:31 — 🗑️ DELETE CODE CELL 15
{0}

## 2026-07-17 17:43:40 — ▶️ RUN CELL 14
**In [14]:** `(*qdetT central value == tg[·][3,1]; on non-symmetric reps this is the ν shift-ordering guard.↵  NB: Chop wraps Simplify because z[k] (hence χ3) are high-precision *numeric* (not exact↵  symbolic) her`
**Out:** (no output)

## 2026-07-17 17:43:52 — ✏️ EDIT CELL 14
(*qdetT check, plain ascii comment test*)
Module[{reps={{{2,1,0},{1,0,0}},{{1,1,0},{1,1,0}},{{2,2,0},{1,0,0}},{{2,1,0},{1,1,0}}}},
 Union@Flatten@Table[
   Chop[Simplify[qdetT[r[[1]],r[[2]]][u] idd[r[…


## 2026-07-17 17:43:52
**Input:**
```wolfram
(*qdetT check, plain ascii comment test*)
Module[{reps={{{2,1,0},{1,0,0}},{{1,1,0},{1,1,0}},{{2,2,0},{1,0,0}},{{2,1,0},{1,1,0}}}},
 Union@Flatten@Table[
   Chop[Simplify[qdetT[r[[1]],r[[2]]][u] idd[r[[1]],r[[2]]]-tg[r[[1]],r[[2]]][3,1][u]],10^-15],{r,reps}]]
```
**Output:**
```
Out= {0}
```
## 2026-07-17 17:44:00 — ✏️ EDIT CELL 14
(*qdetT central value == tg[·][3,1]; on non-symmetric reps this is the ν shift-ordering guard.
  NB: Chop wraps Simplify because z[k] (hence χ3) are high-precision *numeric* (not exact
  symbolic) her…


## 2026-07-17 17:44:00
**Input:**
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
**Output:**
```
Out= {0}
```

## 2026-07-17 17:44:18
**Input:**
```wolfram
$Line
```
**Output:**
```
Out= 201
```
## 2026-07-17 17:44:18 — ▶️ RUN CELL 14
**In [14]:** `(*qdetT central value == tg[·][3,1]; on non-symmetric reps this is the ν shift-ordering guard.↵  NB: Chop wraps Simplify because z[k] (hence χ3) are high-precision *numeric* (not exact↵  symbolic) her`
**Out:** (no output)


## 2026-07-17 17:44:22
**Input:**
```wolfram
{$Line, Out[$Line-1]}
```
**Output:**
```
Out= {202, {0}}
```
## 2026-07-17 17:44:56 — ▶️ RUN CELL 17
**In [17]:** `(* NB: exploratory, not used by alphaCoeffs — ee[λ1,λ2][i,i] was never defined even in V1 (stayed↵   symbolic), so this cell has always been a scratch expression, not a validated formula. *)↵tasympt[1`
**Out:** Out[124]= -h\,u\,\left(z\left(1\right)\,\mathrm{ee}\left(S1,S2\right)\left(1,1\right)+z\left(2\right)\,\mathrm{ee}\left(S1,S2\right)\left(2,2\right)+z\left(3\right)\,\mathrm{ee}\left(S1,S2\right)\left(3,3\right)\right)+u^2\,\left(z\left(1\right)+z\left(2\right)+z\left(3\right)\right)-u\,\left(z\left


## 2026-07-17 17:45:01
**Input:**
```wolfram
{λ1, λ2, ee[λ1,λ2][1,1], FullForm[Hold[ee[λ1, λ2][1,1]]]}
```
**Output:**
```
Out= {λ1, λ2, ee[λ1, λ2][1, 1], FullForm[Hold[ee[λ1, λ2][1, 1]]]}
```

## 2026-07-17 17:45:13
**Input:**
```wolfram
tasympt[1,1]
```
**Output:**
```
Out= (-0.6678881603463455895361101752170204793164128942279036824046`29.820316782640784 - 1.256942910693927779122315235304845126254785777906856664183`30.094928590684432*I)*u + (1.4025651367273257380258313679557430065644670778785977330496`29.820316782640784 + 2.639580112457248336156861994140174765135195253625596338774`30.094928590684432*I)*u^2 - u*((0.3604398751472265372666776505955500912937645232622850024938`29.707347828414424 + 0.9327824485933747164848359820708679294578496404879183613065`30.120295363558725*I)*ee[λ1, λ2][1, 1] + (0.5088673809904998281633213722670686199090953963505721183772`29.857119610800275 + 0.8608449271290790931886443936161642458151438806992171086196`30.085439922389305*I)*ee[λ1, λ2][2, 2] + (0.5332578805895993725958323450931242953616071582657406121786`29.877452280089955 + 0.8459527367347945264833816184531425898622017324384608688479`30.077861097574722*I)*ee[λ1, λ2][3, 3])
```
## 2026-07-17 17:46:32 — 📥 BULK INSERT 1 CELL(S) at positions 41–41
1. [code] (*scratch: test if a cell at the END of the notebook persists its output correctly*)↵Module[{reps={{


## 2026-07-17 17:46:32
**Input:**
```wolfram
(*scratch: test if a cell at the END of the notebook persists its output correctly*)
Module[{reps={{{2,1,0},{1,0,0}},{{1,1,0},{1,1,0}},{{2,2,0},{1,0,0}},{{2,1,0},{1,1,0}}}},
 Union@Flatten@Table[
   Chop[Simplify[qdetT[r[[1]],r[[2]]][u] idd[r[[1]],r[[2]]]-tg[r[[1]],r[[2]]][3,1][u]],10^-15],{r,reps}]]
```
**Output:**
```
(no output)
```
## 2026-07-17 17:46:43 — ✏️ EDIT CELL 41
Chop[{0.,0.,0.},10^-15]


## 2026-07-17 17:46:43
**Input:**
```wolfram
Chop[{0.,0.,0.},10^-15]
```
**Output:**
```
Out= {0, 0, 0}
```
## 2026-07-17 17:46:49 — ✏️ EDIT CELL 41
{0.``28, 1.5}


## 2026-07-17 17:46:49
**Input:**
```wolfram
{0.``28, 1.5}
```
**Output:**
```
Out= {0``28., 1.5}
```
## 2026-07-17 17:47:01 — ✏️ EDIT CELL 41
Chop[5]


## 2026-07-17 17:47:01
**Input:**
```wolfram
Chop[5]
```
**Output:**
```
Out= 5
```
## 2026-07-17 17:47:07 — ✏️ EDIT CELL 41
myChop = Chop; myChop[5]


## 2026-07-17 17:47:07
**Input:**
```wolfram
myChop = Chop; myChop[5]
```
**Output:**
```
Out= 5
```
## 2026-07-17 17:47:17 — ✏️ EDIT CELL 41
(* the word Chop appears only in this comment *)
5


## 2026-07-17 17:47:17
**Input:**
```wolfram
(* the word Chop appears only in this comment *)
5
```
**Output:**
```
Out= 5
```
## 2026-07-17 17:47:24 — 📥 BULK INSERT 1 CELL(S) at positions 42–42
1. [code] 42


## 2026-07-17 17:47:24
**Input:**
```wolfram
42
```
**Output:**
```
(no output)
```
## 2026-07-17 17:47:30 — ▶️ RUN CELL 4
**In [4]:** `(*checking central charge constraint = total boxes, on general reps*)↵Module[{reps={{{2,1,0},{1,0,0}},{{1,1,0},{1,1,0}},{{2,2,0},{1,0,0}},{{2,1,0},{1,1,0}}}},↵ Union@Flatten@Table[↵   C1[r[[1]],r[[2]]`
**Out:** Out[14]= \left\{0\right\}

## 2026-07-17 17:50:48 — ▶️ RUN CELL 42
**In [42]:** `42`
**Out:** (no output)

## 2026-07-17 17:51:18 — 🗑️ DELETE CODE CELL 41
(* the word Chop appears only in this comment *)
5

## 2026-07-17 17:51:18 — 🗑️ DELETE CODE CELL 42
42


## 2026-07-17 17:51:31
**Input:**
```wolfram
Print["hello test"]
```
**Output:**
```
Print:
*STR*hello test
Out= Null
```
## 2026-07-17 17:51:35 — 📥 BULK INSERT 1 CELL(S) at positions 14–14
1. [code] Print["print channel test"]; result = {0,0}


## 2026-07-17 17:51:36
**Input:**
```wolfram
Print["print channel test"]; result = {0,0}
```
**Output:**
```
(no output)
```
## 2026-07-17 17:51:42 — 🗑️ DELETE CODE CELL 14
Print["print channel test"]; result = {0,0}

## 2026-07-17 17:51:55 — ▶️ RUN CELL 11
**In [11]:** `(*quantum minors (q-determinants) of the monodromy matrix*)↵ClearAll[qmin,qming,qminG];↵qmin[λ1_, λ2_][II_List, JJ_List][u_] := ↵    qmin[λ1, λ2][II, JJ][u] = ↵     Sum[Signature[σ] ↵        Dot @@ Ta`
**Out:** (no output)

## 2026-07-17 17:51:55 — ▶️ RUN CELL 12
**In [12]:** `(*transfer matrices with diagonal twist*)↵ClearAll[tg];↵tg[λ1_,λ2_][a_,1][u_]:= Total[qming[λ1, λ2][#, #][u] & /@ Subsets[{1, 2, 3}, {a}]]↵↵(*transfer matrices with companion twist*)↵ClearAll[tG];↵tG[`
**Out:** (no output)

## 2026-07-17 17:51:56 — ▶️ RUN CELL 13
**In [13]:** `(*value of quantum determinant (χ3 carried here; general ν)*)↵ClearAll[qdetT];↵qdetT[λ1_,λ2_][u_]:=χ3 ν[λ1,λ2][1][u-2h] ν[λ1,λ2][2][u-h] ν[λ1,λ2][3][u];`
**Out:** (no output)


## 2026-07-17 17:52:00
**Input:**
```wolfram
finalCheck = Module[{reps={{{2,1,0},{1,0,0}},{{1,1,0},{1,1,0}},{{2,2,0},{1,0,0}},{{2,1,0},{1,1,0}}}},
 Union@Flatten@Table[
   Chop[Simplify[qdetT[r[[1]],r[[2]]][u] idd[r[[1]],r[[2]]]-tg[r[[1]],r[[2]]][3,1][u]],10^-15],{r,reps}]]
```
**Output:**
```
Out= {0}
```

## 2026-07-17 17:59:28
**Input:**
```wolfram
{J0[1], dim[{2,1,0}], dim[{1,1,0}], dim[{2,0,0}], dim[{2,2,0}]}
```
**Output:**
```
Out= {J0[1], 8, 3, 6, 6}
```
## 2026-07-17 17:59:40 — ✏️ EDIT CELL 22 [batch]
(*gl(3) weight multiplicities in the tensor product of two general reps λ1⊗λ2.
  siteWeightMult[λ]: tally the gl(3) weight (J0[1],J0[2],J0[3]) of every GT pattern of λ.
  wtMult: recover (n1,n2,n3) fr…


## 2026-07-17 17:59:40
**Input:**
```wolfram
(*gl(3) weight multiplicities in the tensor product of two general reps λ1⊗λ2.
  siteWeightMult[λ]: tally the gl(3) weight (J0[1],J0[2],J0[3]) of every GT pattern of λ.
  wtMult: recover (n1,n2,n3) from (H1,H2) with N=total boxes, then convolve the two site tallies.*)
ClearAll[siteWeightMult,wtMult];
siteWeightMult[λ_]:=siteWeightMult[λ]=Counts[({J0[1][#],J0[2][#],J0[3][#]}&)/@ps[{λ}]];
wtMult[λ1_,λ2_][H1_,H2_]:=Module[{Ntot=Total[λ1]+Total[λ2],n1,n2,n3,m1,m2},
  n3=(Ntot-H1-2 H2)/3; n2=n3+H2; n1=n3+H1+H2;
  If[!(IntegerQ[n1]&&IntegerQ[n2]&&IntegerQ[n3]),Return[0]];
  m1=siteWeightMult[λ1]; m2=siteWeightMult[λ2];
  Total[KeyValueMap[Function[{w,c}, c Lookup[m2,{n1,n2,n3}-w,0]],m1]]
];
```
**Output:**
```
Out= (no output)
```
## 2026-07-17 17:59:40 — ✏️ EDIT CELL 23 [batch]
(*hard-assert sum-rule: Σ over all (H1,H2) of wtMult == dim[λ1] dim[λ2], for every test rep*)
Module[{reps={{{2,1,0},{1,0,0}},{{2,0,0},{1,1,0}},{{1,1,0},{1,1,0}},{{2,1,0},{1,1,0}},{{2,2,0},{1,0,0}}},b…


## 2026-07-17 17:59:40
**Input:**
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
**Output:**
```
Out= Out[66]= \left\{\text{total weights}\to 18,\text{dim}\to 18,\text{spot (2,1)@(H1,H2)=(1,1)}\to 2\right\}
```
## 2026-07-17 17:59:40 — ✏️ EDIT CELL 24 [batch]
(*guard on {2,1,0}⊗{2,1,0} (d=64) — a rep with genuine interior multiplicity on both sites*)
Module[{λ1={2,1,0},λ2={2,1,0},Ntot,d,tot},
  Ntot=Total[λ1]+Total[λ2]; d=dim[λ1] dim[λ2];
  tot=Total@Flatt…


## 2026-07-17 17:59:41
**Input:**
```wolfram
(*guard on {2,1,0}⊗{2,1,0} (d=64) — a rep with genuine interior multiplicity on both sites*)
Module[{λ1={2,1,0},λ2={2,1,0},Ntot,d,tot},
  Ntot=Total[λ1]+Total[λ2]; d=dim[λ1] dim[λ2];
  tot=Total@Flatten@Table[wtMult[λ1,λ2][H1,H2],{H1,-3 Ntot,3 Ntot},{H2,-3 Ntot,3 Ntot}];
  If[tot!=d,Print["SUM-RULE FAIL: ",tot," != ",d];Abort[]];
  {"total weights"->tot,"dim"->d}
]
```
**Output:**
```
Out= Out[67]= \left\{\text{total weights}\to 36,\text{dim}\to 36\right\}
```
## 2026-07-17 17:59:49 — ▶️ RUN CELL 22
**In [22]:** `(*gl(3) weight multiplicities in the tensor product of two general reps λ1⊗λ2.↵  siteWeightMult[λ]: tally the gl(3) weight (J0[1],J0[2],J0[3]) of every GT pattern of λ.↵  wtMult: recover (n1,n2,n3) fr`
**Out:** (no output)

## 2026-07-17 17:59:51 — ▶️ RUN CELL 23
**In [23]:** `(*hard-assert sum-rule: Σ over all (H1,H2) of wtMult == dim[λ1] dim[λ2], for every test rep*)↵Module[{reps={{{2,1,0},{1,0,0}},{{2,0,0},{1,1,0}},{{1,1,0},{1,1,0}},{{2,1,0},{1,1,0}},{{2,2,0},{1,0,0}}},b`
**Out:** Out[66]= \left\{\text{total weights}\to 18,\text{dim}\to 18,\text{spot (2,1)@(H1,H2)=(1,1)}\to 2\right\}

## 2026-07-17 17:59:52 — ▶️ RUN CELL 24
**In [24]:** `(*guard on {2,1,0}⊗{2,1,0} (d=64) — a rep with genuine interior multiplicity on both sites*)↵Module[{λ1={2,1,0},λ2={2,1,0},Ntot,d,tot},↵  Ntot=Total[λ1]+Total[λ2]; d=dim[λ1] dim[λ2];↵  tot=Total@Flatt`
**Out:** Out[67]= \left\{\text{total weights}\to 36,\text{dim}\to 36\right\}


## 2026-07-17 18:00:06
**Input:**
```wolfram
Module[{reps={{{2,1,0},{1,0,0}},{{2,0,0},{1,1,0}},{{1,1,0},{1,1,0}},{{2,1,0},{1,1,0}},{{2,2,0},{1,0,0}}},bad},
  bad=Select[reps,Function[r,
     Module[{Ntot=Total[r[[1]]]+Total[r[[2]]],d,tot},
       d=dim[r[[1]]] dim[r[[2]]];
       tot=Total@Flatten@Table[wtMult[r[[1]],r[[2]]][H1,H2],
           {H1,-3 Ntot,3 Ntot},{H2,-3 Ntot,3 Ntot}];
       tot=!=d]]];
  If[bad=!={},Print["SUM-RULE FAIL: ",bad];Abort[]];
  {"all reps sum-rule OK",
   "adjoint\[CircleTimes]fund @(H1,H2)=(0,0) mult"->wtMult[{2,1,0},{1,0,0}][0,0]}
]
```
**Output:**
```
Print:
*STR*SUM-RULE FAIL: \034BoxData[RowBox[{"(", "\342\201\240", GridBox[{{RowBox[{"{", RowBox[{"2", ",", "1", ",", "0"}], "}"}], RowBox[{"{", RowBox[{"1", ",", "0", ",", "0"}], "}"}]}, {RowBox[{"{", RowBox[{"2", ",", "0", ",", "0"}], "}"}], RowBox[{"{", RowBox[{"1", ",", "1", ",", "0"}], "}"}]}, {RowBox[{"{", RowBox[{"1", ",", "1", ",", "0"}], "}"}], RowBox[{"{", RowBox[{"1", ",", "1", ",", "0"}], "}"}]}, {RowBox[{"{", RowBox[{"2", ",", "1", ",", "0"}], "}"}], RowBox[{"{", RowBox[{"1", ",", "1", ",", "0"}], "}"}]}, {RowBox[{"{", RowBox[{"2", ",", "2", ",", "0"}], "}"}], RowBox[{"{", RowBox[{"1", ",", "0", ",", "0"}], "}"}]}}, RowSpacings -> 1, ColumnSpacings -> 1, RowAlignments -> Baseline, ColumnAlignments -> Center], "\342\201\240", ")"}]]
(no output)
```

## 2026-07-17 18:00:11
**Input:**
```wolfram
siteWeightMult[{1,0,0}]
siteWeightMult[{2,1,0}]
J0[1][{{1},{1,0},{1,0,0}}]
J0[2][{{1},{1,0},{1,0,0}}]
J0[3][{{1},{1,0},{1,0,0}}]
```
**Output:**
```
Out[1]= <|{1, 0, 0} -> 1, {0, 1, 0} -> 1, {0, 0, 1} -> 1|>
Out[2]= <|{2, 1, 0} -> 1, {2, 0, 1} -> 1, {1, 2, 0} -> 1, {1, 1, 1} -> 2, {1, 0, 2} -> 1, {0, 2, 1} -> 1, {0, 1, 2} -> 1|>
Out[3]= 1
Out[4]= 0
Out[5]= 0
```

## 2026-07-17 18:00:25
**Input:**
```wolfram
Table[{H1,H2,wtMult[{1,0,0},{1,0,0}][H1,H2]},{H1,-2,2},{H2,-2,2}]//Flatten[#,1]&//Select[#,#[[3]]!=0&]&
Total[Flatten[Table[wtMult[{1,0,0},{1,0,0}][H1,H2],{H1,-6,6},{H2,-6,6}]]]
dim[{1,0,0}]^2
```
**Output:**
```
Out[1]= {}
Out[2]= 0
Out[3]= 9
```

## 2026-07-17 18:00:31
**Input:**
```wolfram
wtMult[{1,0,0},{1,0,0}][2,0]
wtMult[{1,0,0},{1,0,0}][0,2]
wtMult[{1,0,0},{1,0,0}][1,-1]
wtMult[{1,0,0},{1,0,0}][-1,1]
wtMult[{1,0,0},{1,0,0}][1,1]
```
**Output:**
```
Out[1]= {0, 0, 0}
Out[2]= 0
Out[3]= {0, 0, 0}
Out[4]= 0
Out[5]= 0
```

## 2026-07-17 18:00:45
**Input:**
```wolfram
m1=siteWeightMult[{1,0,0}];
m2=siteWeightMult[{1,0,0}];
n1=2;n2=0;n3=0;
KeyValueMap[Function[{w,c}, c Lookup[m2,{n1,n2,n3}-w,0]],m1]
Total[%]
```
**Output:**
```
Out[1]= (suppressed)
Out[2]= (suppressed)
Out[3]= (suppressed)
Out[4]= {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}}
Out[5]= Total[$Aborted]
```

## 2026-07-17 18:00:54
**Input:**
```wolfram
Lookup[<|{1,0,0}->5|>, {1,0,0}, 0]
Lookup[<|{1,0,0}->5|>, Key[{1,0,0}], 0]
```
**Output:**
```
Out[1]= {0, 0, 0}
Out[2]= 5
```

## 2026-07-17 18:01:06
**Input:**
```wolfram
{Lookup[<|{1,0,0}->5|>, {1,0,0}, 0], Lookup[<|{1,0,0}->5|>, Key[{1,0,0}], 0]}
```
**Output:**
```
Out= {{0, 0, 0}, 5}
```
## 2026-07-17 18:01:13 — ✏️ EDIT CELL 22
(*gl(3) weight multiplicities in the tensor product of two general reps λ1⊗λ2.
  siteWeightMult[λ]: tally the gl(3) weight (J0[1],J0[2],J0[3]) of every GT pattern of λ.
  wtMult: recover (n1,n2,n3) fr…


## 2026-07-17 18:01:13
**Input:**
```wolfram
(*gl(3) weight multiplicities in the tensor product of two general reps λ1⊗λ2.
  siteWeightMult[λ]: tally the gl(3) weight (J0[1],J0[2],J0[3]) of every GT pattern of λ.
  wtMult: recover (n1,n2,n3) from (H1,H2) with N=total boxes, then convolve the two site tallies.
  NB: Lookup's second argument, when itself a List, is treated as a list of *separate* keys
  (not one list-valued key) — must wrap the composite weight key in Key[...] or Lookup silently
  returns a per-component {0,0,0} instead of the single association value.*)
ClearAll[siteWeightMult,wtMult];
siteWeightMult[λ_]:=siteWeightMult[λ]=Counts[({J0[1][#],J0[2][#],J0[3][#]}&)/@ps[{λ}]];
wtMult[λ1_,λ2_][H1_,H2_]:=Module[{Ntot=Total[λ1]+Total[λ2],n1,n2,n3,m1,m2},
  n3=(Ntot-H1-2 H2)/3; n2=n3+H2; n1=n3+H1+H2;
  If[!(IntegerQ[n1]&&IntegerQ[n2]&&IntegerQ[n3]),Return[0]];
  m1=siteWeightMult[λ1]; m2=siteWeightMult[λ2];
  Total[KeyValueMap[Function[{w,c}, c Lookup[m2,Key[{n1,n2,n3}-w],0]],m1]]
];
```
**Output:**
```
Out= Null
```

## 2026-07-17 18:01:18
**Input:**
```wolfram
Total[Flatten[Table[wtMult[{1,0,0},{1,0,0}][H1,H2],{H1,-6,6},{H2,-6,6}]]]
dim[{1,0,0}]^2
wtMult[{1,0,0},{1,0,0}][2,0]
wtMult[{1,0,0},{1,0,0}][1,1]
wtMult[{2,1,0},{1,0,0}][0,0]
```
**Output:**
```
Out[1]= 9
Out[2]= 9
Out[3]= 1
Out[4]= 0
Out[5]= 0
```
## 2026-07-17 18:01:49 — ▶️ RUN CELL 23
**In [23]:** `(*hard-assert sum-rule: Σ over all (H1,H2) of wtMult == dim[λ1] dim[λ2], for every test rep*)↵Module[{reps={{{2,1,0},{1,0,0}},{{2,0,0},{1,1,0}},{{1,1,0},{1,1,0}},{{2,1,0},{1,1,0}},{{2,2,0},{1,0,0}}},b`
**Out:** Out[66]= \left\{\text{total weights}\to 18,\text{dim}\to 18,\text{spot (2,1)@(H1,H2)=(1,1)}\to 2\right\}

## 2026-07-17 18:01:49 — ▶️ RUN CELL 24
**In [24]:** `(*guard on {2,1,0}⊗{2,1,0} (d=64) — a rep with genuine interior multiplicity on both sites*)↵Module[{λ1={2,1,0},λ2={2,1,0},Ntot,d,tot},↵  Ntot=Total[λ1]+Total[λ2]; d=dim[λ1] dim[λ2];↵  tot=Total@Flatt`
**Out:** Out[67]= \left\{\text{total weights}\to 36,\text{dim}\to 36\right\}


## 2026-07-17 18:02:01
**Input:**
```wolfram
Module[{reps={{{2,1,0},{1,0,0}},{{2,0,0},{1,1,0}},{{1,1,0},{1,1,0}},{{2,1,0},{1,1,0}},{{2,2,0},{1,0,0}}},bad},
  bad=Select[reps,Function[r,
     Module[{Ntot=Total[r[[1]]]+Total[r[[2]]],d,tot},
       d=dim[r[[1]]] dim[r[[2]]];
       tot=Total@Flatten@Table[wtMult[r[[1]],r[[2]]][H1,H2],
           {H1,-3 Ntot,3 Ntot},{H2,-3 Ntot,3 Ntot}];
       tot=!=d]]];
  {bad,
   "all reps sum-rule OK",
   "adjoint\[CircleTimes]fund @(H1,H2)=(0,0) mult"->wtMult[{2,1,0},{1,0,0}][0,0]}
]
```
**Output:**
```
Out= {{}, "all reps sum-rule OK", "adjoint⊗fund @(H1,H2)=(0,0) mult" -> 0}
```

## 2026-07-17 18:02:08
**Input:**
```wolfram
Module[{λ1={2,1,0},λ2={2,1,0},Ntot,d,tot},
  Ntot=Total[λ1]+Total[λ2]; d=dim[λ1] dim[λ2];
  tot=Total@Flatten@Table[wtMult[λ1,λ2][H1,H2],{H1,-3 Ntot,3 Ntot},{H2,-3 Ntot,3 Ntot}];
  {"total weights"->tot,"dim"->d}
]
```
**Output:**
```
Out= {"total weights" -> 64, "dim" -> 64}
```
## 2026-07-17 18:02:13 — ▶️ RUN CELL 23
**In [23]:** `(*hard-assert sum-rule: Σ over all (H1,H2) of wtMult == dim[λ1] dim[λ2], for every test rep*)↵Module[{reps={{{2,1,0},{1,0,0}},{{2,0,0},{1,1,0}},{{1,1,0},{1,1,0}},{{2,1,0},{1,1,0}},{{2,2,0},{1,0,0}}},b`
**Out:** Out[66]= \left\{\text{total weights}\to 18,\text{dim}\to 18,\text{spot (2,1)@(H1,H2)=(1,1)}\to 2\right\}

## 2026-07-17 18:02:13 — ▶️ RUN CELL 24
**In [24]:** `(*guard on {2,1,0}⊗{2,1,0} (d=64) — a rep with genuine interior multiplicity on both sites*)↵Module[{λ1={2,1,0},λ2={2,1,0},Ntot,d,tot},↵  Ntot=Total[λ1]+Total[λ2]; d=dim[λ1] dim[λ2];↵  tot=Total@Flatt`
**Out:** Out[67]= \left\{\text{total weights}\to 36,\text{dim}\to 36\right\}

## 2026-07-17 18:02:26 — ▶️ RUN CELLS 23–24
Cell 23: ✓ — Out[66]= \left\{\text{total weights}\to 18,\text{dim}\to 18,\text{spot (2,1)@(H1,H2)=(1,1)}\to 2\right\}
Cell 24: ✓ — Out[67]= \left\{\text{total weights}\to 36,\text{dim}\to 36\right\}


## 2026-07-17 18:02:31
**Input:**
```wolfram
{Out[66], Out[67], $Line}
```
**Output:**
```
Out= {{"total weights" -> 18, "dim" -> 18, "spot (2,1)@(H1,H2)=(1,1)" -> 2}, {"total weights" -> 36, "dim" -> 36}, 234}
```
## 2026-07-17 18:02:39 — 🗑️ DELETE CODE CELL 23
(*hard-assert sum-rule: Σ over all (H1,H2) of wtMult == dim[λ1] dim[λ2], for every test rep*)
Module[{reps={{{2,1,0},{1,0,0}},{{2,0,0},{1,1,0}},{{1,1,0},{1,1,0}},{{2,1,0},{1,1,0}},{{2,2,0},{1,0,0}}},b…

## 2026-07-17 18:02:39 — 🗑️ DELETE CODE CELL 24
(*guard on {2,1,0}⊗{2,1,0} (d=64) — a rep with genuine interior multiplicity on both sites*)
Module[{λ1={2,1,0},λ2={2,1,0},Ntot,d,tot},
  Ntot=Total[λ1]+Total[λ2]; d=dim[λ1] dim[λ2];
  tot=Total@Flatt…

## 2026-07-17 18:02:48 — 📥 BULK INSERT 2 CELL(S) at positions 23–24
1. [code] (*hard-assert sum-rule: Σ over all (H1,H2) of wtMult == dim[λ1] dim[λ2], for every test rep*)↵Module
2. [code] (*guard on {2,1,0}⊗{2,1,0} (d=64) — a rep with genuine interior multiplicity on both sites*)↵Module[


## 2026-07-17 18:02:48
**Input:**
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
**Output:**
```
(no output)
```

## 2026-07-17 18:02:49
**Input:**
```wolfram
(*guard on {2,1,0}⊗{2,1,0} (d=64) — a rep with genuine interior multiplicity on both sites*)
Module[{λ1={2,1,0},λ2={2,1,0},Ntot,d,tot},
  Ntot=Total[λ1]+Total[λ2]; d=dim[λ1] dim[λ2];
  tot=Total@Flatten@Table[wtMult[λ1,λ2][H1,H2],{H1,-3 Ntot,3 Ntot},{H2,-3 Ntot,3 Ntot}];
  If[tot!=d,Print["SUM-RULE FAIL: ",tot," != ",d];Abort[]];
  {"total weights"->tot,"dim"->d}
]
```
**Output:**
```
(no output)
```
## 2026-07-17 18:02:55 — ▶️ RUN CELL 23
**In [23]:** `(*hard-assert sum-rule: Σ over all (H1,H2) of wtMult == dim[λ1] dim[λ2], for every test rep*)↵Module[{reps={{{2,1,0},{1,0,0}},{{2,0,0},{1,1,0}},{{1,1,0},{1,1,0}},{{2,1,0},{1,1,0}},{{2,2,0},{1,0,0}}},b`
**Out:** (no output)

## 2026-07-17 18:02:56 — ▶️ RUN CELL 24
**In [24]:** `(*guard on {2,1,0}⊗{2,1,0} (d=64) — a rep with genuine interior multiplicity on both sites*)↵Module[{λ1={2,1,0},λ2={2,1,0},Ntot,d,tot},↵  Ntot=Total[λ1]+Total[λ2]; d=dim[λ1] dim[λ2];↵  tot=Total@Flatt`
**Out:** (no output)

## 2026-07-17 18:03:07 — ▶️ RUN CELL 23
**In [23]:** `(*hard-assert sum-rule: Σ over all (H1,H2) of wtMult == dim[λ1] dim[λ2], for every test rep*)↵Module[{reps={{{2,1,0},{1,0,0}},{{2,0,0},{1,1,0}},{{1,1,0},{1,1,0}},{{2,1,0},{1,1,0}},{{2,2,0},{1,0,0}}},b`
**Out:** (no output)


## 2026-07-17 18:07:07
**Input:**
```wolfram
$Line
```
**Output:**
```
Out= 239
```

## 2026-07-17 18:07:12
**Input:**
```wolfram
Table[{n, Quiet[Out[n]]}, {n, 230, 240}]
```
**Output:**
```
Out= {{230, {"all reps sum-rule OK", "adjoint⊗fund @(H1,H2)=(0,0) mult" -> 0}}, {231, {"total weights" -> 64, "dim" -> 64}}, {232, {"all reps sum-rule OK", "adjoint⊗fund @(H1,H2)=(0,0) mult" -> 0}}, {233, {"total weights" -> 64, "dim" -> 64}}, {234, {"all reps sum-rule OK", "adjoint⊗fund @(H1,H2)=(0,0) mult" -> 0}}, {235, {"total weights" -> 64, "dim" -> 64}}, {236, {"all reps sum-rule OK", "adjoint⊗fund @(H1,H2)=(0,0) mult" -> 0}}, {237, {"total weights" -> 64, "dim" -> 64}}, {238, {"all reps sum-rule OK", "adjoint⊗fund @(H1,H2)=(0,0) mult" -> 0}}, {239, Out[239]}, {240, Out[240]}}
```
## 2026-07-17 18:08:10 — 📥 BULK INSERT 1 CELL(S) at positions 25–25
1. [code] {"probe test"->1,"probe test2"->2}


## 2026-07-17 18:08:10
**Input:**
```wolfram
{"probe test"->1,"probe test2"->2}
```
**Output:**
```
(no output)
```

## 2026-07-17 18:08:15
**Input:**
```wolfram
Out[$Line-1]
```
**Output:**
```
Out= {"probe test" -> 1, "probe test2" -> 2}
```
## 2026-07-17 18:08:21 — 🗑️ DELETE CODE CELL 25
{"probe test"->1,"probe test2"->2}

## 2026-07-17 18:08:22 — ▶️ RUN CELL 4
**In [4]:** `(*checking central charge constraint = total boxes, on general reps*)↵Module[{reps={{{2,1,0},{1,0,0}},{{1,1,0},{1,1,0}},{{2,2,0},{1,0,0}},{{2,1,0},{1,1,0}}}},↵ Union@Flatten@Table[↵   C1[r[[1]],r[[2]]`
**Out:** Out[14]= \left\{0\right\}

## 2026-07-17 18:08:46 — ▶️ RUN CELL 23
**In [23]:** `(*hard-assert sum-rule: Σ over all (H1,H2) of wtMult == dim[λ1] dim[λ2], for every test rep*)↵Module[{reps={{{2,1,0},{1,0,0}},{{2,0,0},{1,1,0}},{{1,1,0},{1,1,0}},{{2,1,0},{1,1,0}},{{2,2,0},{1,0,0}}},b`
**Out:** (no output)

## 2026-07-17 18:08:58 — ✏️ EDIT CELL 23
(*hard-assert sum-rule: Σ over all (H1,H2) of wtMult == dim[λ1] dim[λ2], for every test rep.*)
Module[{reps={{{2,1,0},{1,0,0}},{{2,0,0},{1,1,0}},{{1,1,0},{1,1,0}},{{2,1,0},{1,1,0}},{{2,2,0},{1,0,0}}},…


## 2026-07-17 18:08:58
**Input:**
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
**Output:**
```
Out= {"all reps sum-rule OK", "adjoint⊗fund @(H1,H2)=(0,0) mult" -> 0}
```
## 2026-07-17 18:09:06 — ✏️ EDIT CELL 24
(*guard on {2,1,0}⊗{2,1,0} (d=64) — a rep with genuine interior multiplicity on both sites.*)
Module[{λ1={2,1,0},λ2={2,1,0},Ntot,d,tot,result},
  Ntot=Total[λ1]+Total[λ2]; d=dim[λ1] dim[λ2];
  tot=Tot…


## 2026-07-17 18:09:06
**Input:**
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
**Output:**
```
Out= {"total weights" -> 64, "dim" -> 64}
```
## 2026-07-17 18:09:42 — ▶️ RUN CELLS 23–24
Cell 23: ✓ — (no output)
Cell 24: ✓ — (no output)

## 2026-07-17 18:11:17 — ✏️ EDIT CELL 4
(*checking central charge constraint = total boxes, on general reps.*)
Module[{reps={{{2,1,0},{1,0,0}},{{1,1,0},{1,1,0}},{{2,2,0},{1,0,0}},{{2,1,0},{1,1,0}}}},
 Union@Flatten@Table[
   C1[r[[1]],r[[2]…


## 2026-07-17 18:11:17
**Input:**
```wolfram
(*checking central charge constraint = total boxes, on general reps.*)
Module[{reps={{{2,1,0},{1,0,0}},{{1,1,0},{1,1,0}},{{2,2,0},{1,0,0}},{{2,1,0},{1,1,0}}}},
 Union@Flatten@Table[
   C1[r[[1]],r[[2]]]-(Total[r[[1]]]+Total[r[[2]]]) idd[r[[1]],r[[2]]],{r,reps}]]
```
**Output:**
```
Out= {0}
```
## 2026-07-17 18:11:40 — 🗑️ DELETE CODE CELL 23
(*hard-assert sum-rule: Σ over all (H1,H2) of wtMult == dim[λ1] dim[λ2], for every test rep.*)
Module[{reps={{{2,1,0},{1,0,0}},{{2,0,0},{1,1,0}},{{1,1,0},{1,1,0}},{{2,1,0},{1,1,0}},{{2,2,0},{1,0,0}}},…

## 2026-07-17 18:11:40 — 🗑️ DELETE CODE CELL 24
(*guard on {2,1,0}⊗{2,1,0} (d=64) — a rep with genuine interior multiplicity on both sites.*)
Module[{λ1={2,1,0},λ2={2,1,0},Ntot,d,tot,result},
  Ntot=Total[λ1]+Total[λ2]; d=dim[λ1] dim[λ2];
  tot=Tot…

## 2026-07-17 18:11:49 — 📥 BULK INSERT 2 CELL(S) at positions 23–24
1. [code] (*hard-assert sum-rule: Σ over all (H1,H2) of wtMult == dim[λ1] dim[λ2], for every test rep*)↵Module
2. [code] (*guard on {2,1,0}⊗{2,1,0} (d=64) — a rep with genuine interior multiplicity on both sites*)↵Module[

## 2026-07-17 18:11:53 — ▶️ RUN CELL 23
**In [23]:** `(*hard-assert sum-rule: Σ over all (H1,H2) of wtMult == dim[λ1] dim[λ2], for every test rep*)↵Module[{reps={{{2,1,0},{1,0,0}},{{2,0,0},{1,1,0}},{{1,1,0},{1,1,0}},{{2,1,0},{1,1,0}},{{2,2,0},{1,0,0}}},b`
**Out:** (no output)

## 2026-07-17 18:11:59 — ✏️ EDIT CELL 23
(*hard-assert sum-rule: Σ over all (H1,H2) of wtMult == dim[λ1] dim[λ2], for every test rep.*)
Module[{reps={{{2,1,0},{1,0,0}},{{2,0,0},{1,1,0}},{{1,1,0},{1,1,0}},{{2,1,0},{1,1,0}},{{2,2,0},{1,0,0}}},…


## 2026-07-17 18:12:00
**Input:**
```wolfram
(*hard-assert sum-rule: Σ over all (H1,H2) of wtMult == dim[λ1] dim[λ2], for every test rep.*)
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
**Output:**
```
Out= {"all reps sum-rule OK", "adjoint⊗fund @(H1,H2)=(0,0) mult" -> 0}
```
## 2026-07-17 18:14:11 — ✏️ EDIT CELL 23 [batch]
(*hard-assert sum-rule: Σ over all (H1,H2) of wtMult == dim[λ1] dim[λ2], for every test rep*)
Module[{reps={{{2,1,0},{1,0,0}},{{2,0,0},{1,1,0}},{{1,1,0},{1,1,0}},{{2,1,0},{1,1,0}},{{2,2,0},{1,0,0}}},b…


## 2026-07-17 18:14:11
**Input:**
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
**Output:**
```
Out= (no output)
```
## 2026-07-17 18:14:11 — ✏️ EDIT CELL 24 [batch]
(*guard on {2,1,0}⊗{2,1,0} (d=64) — a rep with genuine interior multiplicity on both sites*)
Module[{λ1={2,1,0},λ2={2,1,0},Ntot,d,tot},
  Ntot=Total[λ1]+Total[λ2]; d=dim[λ1] dim[λ2];
  tot=Total@Flatt…


## 2026-07-17 18:14:11
**Input:**
```wolfram
(*guard on {2,1,0}⊗{2,1,0} (d=64) — a rep with genuine interior multiplicity on both sites*)
Module[{λ1={2,1,0},λ2={2,1,0},Ntot,d,tot},
  Ntot=Total[λ1]+Total[λ2]; d=dim[λ1] dim[λ2];
  tot=Total@Flatten@Table[wtMult[λ1,λ2][H1,H2],{H1,-3 Ntot,3 Ntot},{H2,-3 Ntot,3 Ntot}];
  If[tot!=d,Print["SUM-RULE FAIL: ",tot," != ",d];Abort[]];
  {"total weights"->tot,"dim"->d}
]
```
**Output:**
```
Out= (no output)
```

## 2026-07-17 18:14:21
**Input:**
```wolfram
check22 = Module[{reps={{{2,1,0},{1,0,0}},{{2,0,0},{1,1,0}},{{1,1,0},{1,1,0}},{{2,1,0},{1,1,0}},{{2,2,0},{1,0,0}}},bad},
  bad=Select[reps,Function[r,
     Module[{Ntot=Total[r[[1]]]+Total[r[[2]]],d,tot},
       d=dim[r[[1]]] dim[r[[2]]];
       tot=Total@Flatten@Table[wtMult[r[[1]],r[[2]]][H1,H2],
           {H1,-3 Ntot,3 Ntot},{H2,-3 Ntot,3 Ntot}];
       tot=!=d]]];
  If[bad=!={},{"SUM-RULE FAIL",bad},
  {"all reps sum-rule OK",
   "adjoint\[CircleTimes]fund @(H1,H2)=(0,0) mult"->wtMult[{2,1,0},{1,0,0}][0,0]}]
]
check23 = Module[{\[Lambda]1={2,1,0},\[Lambda]2={2,1,0},Ntot,d,tot},
  Ntot=Total[\[Lambda]1]+Total[\[Lambda]2]; d=dim[\[Lambda]1] dim[\[Lambda]2];
  tot=Total@Flatten@Table[wtMult[\[Lambda]1,\[Lambda]2][H1,H2],{H1,-3 Ntot,3 Ntot},{H2,-3 Ntot,3 Ntot}];
  If[tot!=d,{"SUM-RULE FAIL",tot,d},{"total weights"->tot,"dim"->d}]
]
{check22, check23}
```
**Output:**
```
Out[1]= {"all reps sum-rule OK", "adjoint⊗fund @(H1,H2)=(0,0) mult" -> 0}
Out[2]= {"total weights" -> 64, "dim" -> 64}
Out[3]= {{"all reps sum-rule OK", "adjoint⊗fund @(H1,H2)=(0,0) mult" -> 0}, {"total weights" -> 64, "dim" -> 64}}
```

## 2026-07-17 18:14:24
**Input:**
```wolfram
{"dim{2,1,0}"->dim[{2,1,0}], "dim{1,1,0}"->dim[{1,1,0}], "dim{2,0,0}"->dim[{2,0,0}], "dim{2,2,0}"->dim[{2,2,0}]}
siteWeightMult[{2,1,0}]
```
**Output:**
```
Out[1]= {"dim{2,1,0}" -> 8, "dim{1,1,0}" -> 3, "dim{2,0,0}" -> 6, "dim{2,2,0}" -> 6}
Out[2]= <|{2, 1, 0} -> 1, {2, 0, 1} -> 1, {1, 2, 0} -> 1, {1, 1, 1} -> 2, {1, 0, 2} -> 1, {0, 2, 1} -> 1, {0, 1, 2} -> 1|>
```
## 2026-07-17 18:14:37 — ✏️ EDIT CELL 23
(*hard-assert sum-rule: Σ over all (H1,H2) of wtMult == dim[λ1] dim[λ2], for every test rep*)
Module[{reps={{{2,1,0},{1,0,0}},{{2,0,0},{1,1,0}},{{1,1,0},{1,1,0}},{{2,1,0},{1,1,0}},{{2,2,0},{1,0,0}}},b…


## 2026-07-17 18:14:37
**Input:**
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
**Output:**
```
Out= {"all reps sum-rule OK", "adjoint⊗fund @(H1,H2)=(0,0) mult" -> 0}
```
## 2026-07-17 18:14:44 — ✏️ EDIT CELL 24
(*guard on {2,1,0}⊗{2,1,0} (d=64) — a rep with genuine interior multiplicity on both sites*)
Module[{λ1={2,1,0},λ2={2,1,0},Ntot,d,tot},
  Ntot=Total[λ1]+Total[λ2]; d=dim[λ1] dim[λ2];
  tot=Total@Flatt…


## 2026-07-17 18:14:44
**Input:**
```wolfram
(*guard on {2,1,0}⊗{2,1,0} (d=64) — a rep with genuine interior multiplicity on both sites*)
Module[{λ1={2,1,0},λ2={2,1,0},Ntot,d,tot},
  Ntot=Total[λ1]+Total[λ2]; d=dim[λ1] dim[λ2];
  tot=Total@Flatten@Table[wtMult[λ1,λ2][H1,H2],{H1,-3 Ntot,3 Ntot},{H2,-3 Ntot,3 Ntot}];
  If[tot!=d,Print["SUM-RULE FAIL: ",tot," != ",d];Abort[]];
  {"total weights"->tot,"dim"->d}
]
```
**Output:**
```
Out= {"total weights" -> 64, "dim" -> 64}
```
## 2026-07-17 18:57:09 — ✏️ EDIT CELL 19 [batch]
(*coefficient matrices of tau1=tg[.][1,1] (deg 2 in u) and tau2=tg[.][2,1] (deg 4 in u),
  extracted by CoefficientList in a dummy variable uu so tg[..][..][uu] is treated as a polynomial in uu.
  c1=…


## 2026-07-17 18:57:09
**Input:**
```wolfram
(*coefficient matrices of tau1=tg[.][1,1] (deg 2 in u) and tau2=tg[.][2,1] (deg 4 in u),
  extracted by CoefficientList in a dummy variable uu so tg[..][..][uu] is treated as a polynomial in uu.
  c1={C0,C1,C2} with tau1(u)=Sum_k C_k u^k (index k+1 in the list); c2={D0,...,D4} similarly for tau2.*)
ClearAll[tgCoeffs];
tgCoeffs[λ1_, λ2_] := tgCoeffs[λ1, λ2] = Module[{uu, m1, m2, cl1, cl2},
   m1 = tg[λ1, λ2][1, 1][uu];
   m2 = tg[λ1, λ2][2, 1][uu];
   cl1 = Map[PadRight[CoefficientList[#, uu], 3] &, m1, {2}];
   cl2 = Map[PadRight[CoefficientList[#, uu], 5] &, m2, {2}];
   {Table[cl1[[All, All, k]], {k, 1, 3}], Table[cl2[[All, All, k]], {k, 1, 5}]}
];
```
**Output:**
```
Out= (no output)
```
## 2026-07-17 18:57:09 — ✏️ EDIT CELL 20 [batch]
(*Cartan operators (gl(3) weight generators H1=E11-E22, H2=E22-E33) acting on the two-site space*)
ClearAll[H1op, H2op];
H1op[λ1_, λ2_] := EE[λ1, λ2][1, 1] - EE[λ1, λ2][2, 2];
H2op[λ1_, λ2_] := EE[λ1,…


## 2026-07-17 18:57:09
**Input:**
```wolfram
(*Cartan operators (gl(3) weight generators H1=E11-E22, H2=E22-E33) acting on the two-site space*)
ClearAll[H1op, H2op];
H1op[λ1_, λ2_] := EE[λ1, λ2][1, 1] - EE[λ1, λ2][2, 2];
H2op[λ1_, λ2_] := EE[λ1, λ2][2, 2] - EE[λ1, λ2][3, 3];
```
**Output:**
```
Out= (no output)
```
## 2026-07-17 18:57:09 — ✏️ EDIT CELL 25 [batch]
ClearAll[TauEigensystem3];
TauEigensystem3::rankdef = "Generic combination failed to resolve the `1`-dim eigenspace for (S1,S2)=(`2`,`3`) at every tried coefficient set.";
TauEigensystem3::notgenuine …


## 2026-07-17 18:57:09
**Input:**
```wolfram
ClearAll[TauEigensystem3];
TauEigensystem3::rankdef = "Generic combination failed to resolve the `1`-dim eigenspace for (S1,S2)=(`2`,`3`) at every tried coefficient set.";
TauEigensystem3::notgenuine = "Eigenvector #`1` for (S1,S2)=(`2`,`3`) failed the common-eigenvector residual check (residual=`4`).";
TauEigensystem3::noninteger = "Weight (H1 or H2)=`1` not near-integer for eigenvector #`2` at (S1,S2)=(`3`,`4`).";
TauEigensystem3::degnull = "Sector (H1,H2)=(`1`,`2`) for (S1,S2)=(`3`,`4`) has `5` state(s) but wtMult predicts `6`.";
TauEigensystem3[λ1_, λ2_] := TauEigensystem3[λ1, λ2] = Module[
   {d, c1, c2, mats, nontriv, rSets, rset, Kmix, vals, vecs, rayleigh, recs, byHW, table},
   d = dim[λ1] dim[λ2];
   {c1, c2} = tgCoeffs[λ1, λ2];
   (* full family for extraction; drop the two scalar leading coeffs (c1[[3]],c2[[5]]) from the
      splitting combination since they carry no state info *)
   mats = Join[c1, c2, {H1op[λ1, λ2], H2op[λ1, λ2]}];
   nontriv = Join[c1[[1 ;; 2]], c2[[1 ;; 4]], {H1op[λ1, λ2], H2op[λ1, λ2]}]; (* 8 matrices *)
   rSets = Table[N[Table[Prime[100 i + 7 j]^(1/3), {j, Length[nontriv]}], 30], {i, 1, 8}];
   rset = SelectFirst[rSets,
     Module[{ee, vv}, {ee, vv} = Eigensystem[#.nontriv];
        MatrixRank[vv] == d && Length[Union[Chop[ee, 10^-12]]] == d] &,
     Missing["NoGenericR"]];
   If[MissingQ[rset], Message[TauEigensystem3::rankdef, d, λ1, λ2]; Abort[]];
   Kmix = rset . nontriv;
   {vals, vecs} = Eigensystem[N[Kmix, 30]];
   rayleigh[m_, v_] := (v . (m . v))/(v . v);
   recs = Table[
     Module[{v = vecs[[i]], quotients, res, H1v, H2v, t1c, t2c},
       quotients = rayleigh[#, v] & /@ mats;
       res = Max[Abs[#.v - rayleigh[#, v] v]] & /@ mats // Max;
       If[res > 10^-13, Message[TauEigensystem3::notgenuine, i, λ1, λ2, res]; Abort[]];
       t1c = rayleigh[#, v] & /@ c1;           (* {C0,C1,C2} quotients *)
       t2c = rayleigh[#, v] & /@ c2;           (* {D0..D4} quotients *)
       H1v = rayleigh[H1op[λ1, λ2], v]; H2v = rayleigh[H2op[λ1, λ2], v];
       If[Abs[H1v - Round[Re[H1v]]] > 10^-10, Message[TauEigensystem3::noninteger, H1v, i, λ1, λ2]; Abort[]];
       If[Abs[H2v - Round[Re[H2v]]] > 10^-10, Message[TauEigensystem3::noninteger, H2v, i, λ1, λ2]; Abort[]];
       <|"H1" -> Round[Re[H1v]], "H2" -> Round[Re[H2v]], "t1c" -> t1c, "t2c" -> t2c, "vec" -> v|>],
     {i, 1, d}];
   byHW = GroupBy[recs, {#["H1"], #["H2"]} &];
   table = Association @ Flatten @ KeyValueMap[
      Function[{hw, group},
        Module[{sorted, predicted, last},
          predicted = wtMult[λ1, λ2][hw[[1]], hw[[2]]];
          (* order n by the full real+imaginary parts of ALL tau1/tau2 Rayleigh quotients:
             within a fixed (H1,H2) sector, states are separated by their tau1/tau2 charges
             since the generic combination Kmix has simple spectrum there, so this tuple gives
             a deterministic total order (a single C0 coefficient is not guaranteed distinct) *)
          sorted = SortBy[group, Join[Re[#["t1c"]], Im[#["t1c"]], Re[#["t2c"]], Im[#["t2c"]]] &];
          If[Length[sorted] != predicted,
            Message[TauEigensystem3::degnull, hw[[1]], hw[[2]], λ1, λ2, Length[sorted], predicted]; Abort[]];
          MapIndexed[
            Function[{rec, idx}, Module[{n = First[idx] - 1, vv = rec["vec"], lc},
              lc = SelectFirst[vv, Abs[#] > 10^-13 Max[Abs[vv]] &];
              {hw[[1]], hw[[2]], n} -> <|"tau1c" -> rec["t1c"], "tau2c" -> rec["t2c"], "Psi" -> vv/lc|>]],
            sorted]]],
      byHW];
   If[Length[table] != d, Message[TauEigensystem3::degnull, "total", "-", λ1, λ2, Length[table], d]; Abort[]];
   table
];
```
**Output:**
```
Out= (no output)
```
## 2026-07-17 18:57:09 — ✏️ EDIT CELL 26 [batch]
ClearAll[Psi, tau1, tau2];
Psi::bnds = "No state at (H1,H2,n)=(`1`,`2`,`3`) for (S1,S2)=(`4`,`5`).";
Psi[λ1_, λ2_][H1_, H2_, n_] := Lookup[TauEigensystem3[λ1, λ2], Key[{H1, H2, n}],
    (Message[Psi::…


## 2026-07-17 18:57:09
**Input:**
```wolfram
ClearAll[Psi, tau1, tau2];
Psi::bnds = "No state at (H1,H2,n)=(`1`,`2`,`3`) for (S1,S2)=(`4`,`5`).";
Psi[λ1_, λ2_][H1_, H2_, n_] := Lookup[TauEigensystem3[λ1, λ2], Key[{H1, H2, n}],
    (Message[Psi::bnds, H1, H2, n, λ1, λ2]; Missing["NoState", {H1, H2, n}])]["Psi"];
tau1[λ1_, λ2_][H1_, H2_, n_][u_] := TauEigensystem3[λ1, λ2][{H1, H2, n}, "tau1c"] . u^Range[0, 2];
tau2[λ1_, λ2_][H1_, H2_, n_][u_] := TauEigensystem3[λ1, λ2][{H1, H2, n}, "tau2c"] . u^Range[0, 4];
```
**Output:**
```
Out= (no output)
```
## 2026-07-17 18:57:09 — ✏️ EDIT CELL 21 [batch]
(*hard-assert: τ1,τ2 leading coeffs are χ1 Id, χ2 Id; full commuting family pairwise-commutes.
  Tested on the adjoint⊗fundamental (d=24), which has an interior multiplicity-2 weight.*)
Module[{λ1={2,…


## 2026-07-17 18:57:10
**Input:**
```wolfram
(*hard-assert: τ1,τ2 leading coeffs are χ1 Id, χ2 Id; full commuting family pairwise-commutes.
  Tested on the adjoint⊗fundamental (d=24), which has an interior multiplicity-2 weight.*)
Module[{λ1={2,1,0},λ2={1,0,0},d,c1,c2,C2,D4,mats,comm},
  d=dim[λ1] dim[λ2];
  {c1,c2}=tgCoeffs[λ1,λ2];
  C2=c1[[3]]; D4=c2[[5]];
  If[Max[Abs[Flatten[C2-χ1 IdentityMatrix[d]]]]>10^-13,Print["C2 FAIL"];Abort[]];
  If[Max[Abs[Flatten[D4-χ2 IdentityMatrix[d]]]]>10^-13,Print["D4 FAIL"];Abort[]];
  mats=Join[c1,c2,{H1op[λ1,λ2],H2op[λ1,λ2]}];
  comm=Max[Table[Max[Abs[Flatten[mats[[i]].mats[[j]]-mats[[j]].mats[[i]]]]],
      {i,Length[mats]},{j,Length[mats]}]];
  If[comm>10^-12,Print["COMMUTATION FAIL: ",comm];Abort[]];
  {"C2=χ1 Id OK","D4=χ2 Id OK","maxCommutator"->comm}
]
```
**Output:**
```
Out= Out[63]= \left\{\text{C2=$\chi$1 Id OK},\text{D4=$\chi$2 Id OK},\text{maxCommutator}\to 0.\times 10^{-26}\right\}
```
## 2026-07-17 18:57:10 — ✏️ EDIT CELL 27 [batch]
Module[{λ1 = {2, 1, 0}, λ2 = {1, 0, 0}, es, worst},
  es = TauEigensystem3[λ1, λ2];
  worst = Max@KeyValueMap[
     Function[{k, rec},
       Module[{H1 = k[[1]], H2 = k[[2]], n = k[[3]], v = rec["Psi…


## 2026-07-17 18:57:10
**Input:**
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
  {"nstates" -> Length[es], "worstResidual@u=17/5" -> worst}
]
```
**Output:**
```
Out= Out[85]= \left\{\text{nstates}\to 18,\text{worstResidual@u=17/5}\to 1.359292769213167893921303\times 10^{-27}\right\}
```
## 2026-07-17 18:57:18 — ▶️ RUN CELL 21
**In [21]:** `(*hard-assert: τ1,τ2 leading coeffs are χ1 Id, χ2 Id; full commuting family pairwise-commutes.↵  Tested on the adjoint⊗fundamental (d=24), which has an interior multiplicity-2 weight.*)↵Module[{λ1={2,`
**Out:** Out[63]= \left\{\text{C2=$\chi$1 Id OK},\text{D4=$\chi$2 Id OK},\text{maxCommutator}\to 0.\times 10^{-26}\right\}


## 2026-07-17 18:57:23
**Input:**
```wolfram
{dim[{2,1,0}], dim[{1,0,0}], dim[{2,1,0}] dim[{1,0,0}]}
```
**Output:**
```
Out= {8, 3, 24}
```

## 2026-07-17 18:57:23
**Input:**
```wolfram
{c1,c2}=tgCoeffs[{2,1,0},{1,0,0}]; Dimensions[c1[[3]]]
```
**Output:**
```
Out= {24, 24}
```
## 2026-07-17 18:57:28 — ▶️ RUN CELL 27
**In [27]:** `Module[{λ1 = {2, 1, 0}, λ2 = {1, 0, 0}, es, worst},↵  es = TauEigensystem3[λ1, λ2];↵  worst = Max@KeyValueMap[↵     Function[{k, rec},↵       Module[{H1 = k[[1]], H2 = k[[2]], n = k[[3]], v = rec["Psi`
**Out:** Out[85]= \left\{\text{nstates}\to 18,\text{worstResidual@u=17/5}\to 1.359292769213167893921303\times 10^{-27}\right\}


## 2026-07-17 18:57:36
**Input:**
```wolfram
ClearAll[TauEigensystem3]; DownValues[TauEigensystem3]
```
**Output:**
```
Out= {}
```
## 2026-07-17 18:57:41 — ▶️ RUN CELL 25
**In [25]:** `ClearAll[TauEigensystem3];↵TauEigensystem3::rankdef = "Generic combination failed to resolve the `1`-dim eigenspace for (S1,S2)=(`2`,`3`) at every tried coefficient set.";↵TauEigensystem3::notgenuine `
**Out:** (no output)

## 2026-07-17 18:57:41 — ▶️ RUN CELL 26
**In [26]:** `ClearAll[Psi, tau1, tau2];↵Psi::bnds = "No state at (H1,H2,n)=(`1`,`2`,`3`) for (S1,S2)=(`4`,`5`).";↵Psi[λ1_, λ2_][H1_, H2_, n_] := Lookup[TauEigensystem3[λ1, λ2], Key[{H1, H2, n}],↵    (Message[Psi::`
**Out:** (no output)


## 2026-07-17 18:57:43
**Input:**
```wolfram
Length[TauEigensystem3[{2,1,0},{1,0,0}]]
```
**Output:**
```
Out= 24
```
## 2026-07-17 18:57:50 — ▶️ RUN CELL 27
**In [27]:** `Module[{λ1 = {2, 1, 0}, λ2 = {1, 0, 0}, es, worst},↵  es = TauEigensystem3[λ1, λ2];↵  worst = Max@KeyValueMap[↵     Function[{k, rec},↵       Module[{H1 = k[[1]], H2 = k[[2]], n = k[[3]], v = rec["Psi`
**Out:** Out[85]= \left\{\text{nstates}\to 18,\text{worstResidual@u=17/5}\to 1.359292769213167893921303\times 10^{-27}\right\}


## 2026-07-17 18:57:56
**Input:**
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
  {"nstates" -> Length[es], "worstResidual@u=17/5" -> worst}
]
```
**Output:**
```
Out= {"nstates" -> 24, "worstResidual@u=17/5" -> 1.48498530494543347294990299017607961167850701950923`25.011916239636957*^-27}
```
## 2026-07-17 18:58:09 — ✏️ EDIT CELL 27
Module[{λ1 = {2, 1, 0}, λ2 = {1, 0, 0}, es, worst},
  es = TauEigensystem3[λ1, λ2];
  worst = Max@KeyValueMap[
     Function[{k, rec},
       Module[{H1 = k[[1]], H2 = k[[2]], n = k[[3]], v = rec["Psi…


## 2026-07-17 18:58:09
**Input:**
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
**Output:**
```
Out= {"nstates" -> 24, "worstResidual@u=17/5" -> 1.48498530494543347294990299017607961167850701950923`25.011916239636957*^-27}
```

## 2026-07-17 18:58:12
**Input:**
```wolfram
$MessageList
```
**Output:**
```
Out= {}
```
## 2026-07-17 18:58:18 — ▶️ RUN CELL 27
**In [27]:** `Module[{λ1 = {2, 1, 0}, λ2 = {1, 0, 0}, es, worst},↵  es = TauEigensystem3[λ1, λ2];↵  worst = Max@KeyValueMap[↵     Function[{k, rec},↵       Module[{H1 = k[[1]], H2 = k[[2]], n = k[[3]], v = rec["Psi`
**Out:** Out[85]= \left\{\text{nstates}\to 18,\text{worstResidual@u=17/5}\to 1.359292769213167893921303\times 10^{-27}\right\}


## 2026-07-17 18:58:26
**Input:**
```wolfram
Out[85]
```
**Output:**
```
Out= {"nstates" -> 18, "worstResidual@u=17/5" -> 1.35929276921316789392130261643327320417130403047792`24.938614497066567*^-27}
```

## 2026-07-17 18:58:26
**Input:**
```wolfram
$Line
```
**Output:**
```
Out= 280
```

## 2026-07-17 18:58:37
**Input:**
```wolfram
Length[TauEigensystem3[{2,1,0},{1,0,0}]]
```
**Output:**
```
Out= 24
```

## 2026-07-17 18:58:37
**Input:**
```wolfram
Length[DownValues[TauEigensystem3]]
```
**Output:**
```
Out= 2
```
## 2026-07-17 18:59:34 — 🗑️ DELETE CODE CELL 27
Module[{λ1 = {2, 1, 0}, λ2 = {1, 0, 0}, es, worst},
  es = TauEigensystem3[λ1, λ2];
  worst = Max@KeyValueMap[
     Function[{k, rec},
       Module[{H1 = k[[1]], H2 = k[[2]], n = k[[3]], v = rec["Psi…

## 2026-07-17 18:59:40 — 📥 BULK INSERT 1 CELL(S) at positions 27–27
1. [code] Module[{λ1 = {2, 1, 0}, λ2 = {1, 0, 0}, es, worst},↵  es = TauEigensystem3[λ1, λ2];↵  worst = Max@Ke


## 2026-07-17 18:59:40
**Input:**
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
**Output:**
```
(no output)
```
## 2026-07-17 18:59:47 — ▶️ RUN CELL 27
**In [27]:** `Module[{λ1 = {2, 1, 0}, λ2 = {1, 0, 0}, es, worst},↵  es = TauEigensystem3[λ1, λ2];↵  worst = Max@KeyValueMap[↵     Function[{k, rec},↵       Module[{H1 = k[[1]], H2 = k[[2]], n = k[[3]], v = rec["Psi`
**Out:** (no output)


## 2026-07-17 18:59:50
**Input:**
```wolfram
Length[TauEigensystem3[{2,1,0},{1,0,0}]]
```
**Output:**
```
Out= 24
```
## 2026-07-17 18:59:55 — ▶️ RUN CELL 27
**In [27]:** `Module[{λ1 = {2, 1, 0}, λ2 = {1, 0, 0}, es, worst},↵  es = TauEigensystem3[λ1, λ2];↵  worst = Max@KeyValueMap[↵     Function[{k, rec},↵       Module[{H1 = k[[1]], H2 = k[[2]], n = k[[3]], v = rec["Psi`
**Out:** (no output)

## 2026-07-17 19:00:15 — ✏️ EDIT CELL 27
Module[{λ1 = {2, 1, 0}, λ2 = {1, 0, 0}, es, worst},
  es = TauEigensystem3[λ1, λ2];
  worst = Max@KeyValueMap[
     Function[{k, rec},
       Module[{H1 = k[[1]], H2 = k[[2]], n = k[[3]], v = rec["Psi…


## 2026-07-17 19:00:15
**Input:**
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
**Output:**
```
Out= {"nstates" -> 24, "worstResidual@u=17/5" -> 1.48498530494543347294990299017607961167850701950923`25.011916239636957*^-27}
```

## 2026-07-17 19:00:25
**Input:**
```wolfram
ClearAll[TauEigensystem3]
```
**Output:**
```
Out= Null
```
## 2026-07-17 19:00:25 — ▶️ RUN CELL 25
**In [25]:** `ClearAll[TauEigensystem3];↵TauEigensystem3::rankdef = "Generic combination failed to resolve the `1`-dim eigenspace for (S1,S2)=(`2`,`3`) at every tried coefficient set.";↵TauEigensystem3::notgenuine `
**Out:** (no output)


## 2026-07-17 19:00:28
**Input:**
```wolfram
Quiet[Off[]]; ClearAll[msgs]; msgs = {}; Internal`WithLocalSettings[Null, res = Check[TauEigensystem3[{2,1,0},{1,0,0}], $Failed], Null]; {Length[res], res===$Failed}
```
**Output:**
```
Out= {24, False}
```

## 2026-07-17 19:00:34
**Input:**
```wolfram
ClearAll[TauEigensystem3]; res = Check[TauEigensystem3[{2,1,0},{1,0,0}], "MSGFIRED", {TauEigensystem3::rankdef, TauEigensystem3::notgenuine, TauEigensystem3::noninteger, TauEigensystem3::degnull}]; {res===$Failed, res==="MSGFIRED", Length[res]}
```
**Output:**
```
Out= {False, False, 2}
```

## 2026-07-17 19:00:42
**Input:**
```wolfram
Head[res]
Length[res]
Keys[res]//Length
```
**Output:**
```
Out[1]= TauEigensystem3
Out[2]= 2
[message] Keys::invrl: The argument TauEigensystem3[{2, 1, 0}, {1, 0, 0}] is not a valid Association or a list of rules.
Out[3]= 1
```
## 2026-07-17 19:00:49 — ▶️ RUN CELL 25
**In [25]:** `ClearAll[TauEigensystem3];↵TauEigensystem3::rankdef = "Generic combination failed to resolve the `1`-dim eigenspace for (S1,S2)=(`2`,`3`) at every tried coefficient set.";↵TauEigensystem3::notgenuine `
**Out:** (no output)

## 2026-07-17 19:00:49 — ▶️ RUN CELL 26
**In [26]:** `ClearAll[Psi, tau1, tau2];↵Psi::bnds = "No state at (H1,H2,n)=(`1`,`2`,`3`) for (S1,S2)=(`4`,`5`).";↵Psi[λ1_, λ2_][H1_, H2_, n_] := Lookup[TauEigensystem3[λ1, λ2], Key[{H1, H2, n}],↵    (Message[Psi::`
**Out:** (no output)


## 2026-07-17 19:00:53
**Input:**
```wolfram
res2 = Check[TauEigensystem3[{2,1,0},{1,0,0}], "MSGFIRED", {TauEigensystem3::rankdef, TauEigensystem3::notgenuine, TauEigensystem3::noninteger, TauEigensystem3::degnull}]
{res2 === "MSGFIRED", Head[res2], Length[res2]}
```
**Output:**
```
Out[1]= <|{-2, 0, 0} -> <|"tau1c" -> {0.46381111429156251648688025744937764605`26.88749351941445 + 0.77593928190604048669556527022060552902`26.889208877815733*I, -2.75213868350654399105441760993740630986`26.995676300283108 - 4.67053823842167501846636725944345879762`26.995465816592894*I, 1.40256513672732573802583136795574300657`26.995566725475317 + 2.63958011245724833615686199414017476514`26.995571924717474*I}, "tau2c" -> {-0.84778952360710218023484598604636440273`26.975088827073936 + 1.50538603715106001670329455404266524107`26.985331713596658*I, 7.81514217268905581437290108262205659057`26.98650936178515 - 13.59699411294372192820981042805691576295`26.991692034391836*I, -14.44454351108724654331452688325230007773`26.987614808999254 + 23.37094424011514985850803801856920445867`26.99177390198123*I, 9.20069099678327144792717223660316574221`26.99272879523839 - 14.04513805866919027293107472602002663704`26.99562003693214*I, -1.67332398590162909911824721890393600537`26.992842675051804 + 2.47680418504351011497841273933203086295`26.99560680057734*I}, "Psi" -> {0, 1.495795369633247384233108152393900977832`27.297104956853076*^-38 - 6.78640934011821406956277471806890883441`27.29816989002037*^-39*I, -4.82407083178249438669760754428865842118`27.29670219683249*^-39 + 3.01286353206107880707972524208102320805`27.29816989002038*^-39*I, 0, 3.499139236798184683854659404808697649131`27.297743177675432*^-38 - 6.41367872811600429943006027417493639641`27.298169890020368*^-39*I, -3.6525966630370871360265991949576715180702`27.297318164290054*^-37 + 1.3270671760938563283711706830181042153687`27.298169890020368*^-37*I, -6.73680668764023242159822301078752664182`27.29816989002037*^-39 - 1.129030206034075665495516346272652518434`27.296767844570617*^-38*I, 0, -1.34409788947338412200370429972829641438`27.29758188147997*^-39 + 3.3816785000345043147577156961703027384`27.298169890020365*^-40*I, 1.256579933315460121886525315592948529058`27.296278448512656*^-38 - 1.010882252560031082489893724135982398632`27.298169890020368*^-38*I, -6.3965264612479498118955782215675804083`27.29816989002037*^-40 - 6.7826521683468818613086570114119610609542767134116524342`27.27643258016508*^-41*I, 0, 6.50864049001534999449652129881411957495`27.29583368346157*^-39 - 6.46637912372236177976268031514859835158`27.29816989002038*^-39*I, -9.3929994855880844540453150773107534445`27.29732667136141*^-40 + 3.3788159149627337555735782782518291394`27.298169890020365*^-40*I, 0, 2.9176901849028391132013036718276717116325`27.297007278260622*^-37 - 1.4446047006124378981148925998279680154435`27.29816989002037*^-37*I, 0, 0, 1.1893319816302674979701575943281639295`27.298169890020382*^-40 + 1.8474480879434709163308226710449803697`27.29665690562218*^-40*I, 0, 1.00000000000000000000000000000000000001`27.29816989002037 + 0``29.56225370354431*I, 0, -1.33332839380882229262156015420535226882`27.29816989002037 + 0.003629330649565128314083549247410045244559536918279656335`27.29816989002037*I, 0}|>, {-2, 0, 1} -> <|"tau1c" -> {1.50589197865692900275789832722865206244`26.962031489628593 + 2.48266423180242716093199122667676125435`26.990551081457607*I, -2.75213868350654399105441760993740630985`26.961135724440545 - 4.67053823842167501846636725944345879761`26.990809871934182*I, 1.40256513672732573802583136795574300656`26.95774557791067 + 2.63958011245724833615686199414017476513`26.99102577575576*I}, "tau2c" -> {-0.90571312271160645327094608212017138338`26.972535458346 + 1.58096729753898082790241087765140503646`26.974368075976283*I, 8.39437816373409854473390204336012639709`26.98390749792578 - 14.35280671682293004020097366414431371683`26.980538665971117*I, -15.66093909228183627707262890080224667143`26.98572224226365 + 24.95815070826148689368948081435274016185`26.98046971995756*I, 9.20069099678327144792717223660316574217`26.986279662278818 - 14.04513805866919027293107472602002663698`26.980073004771423*I, -1.67332398590162909911824721890393600536`26.986480976366956 + 2.47680418504351011497841273933203086294`26.97967761907439*I}, "Psi" -> {0, 2.5935407412980038746320839783501772912`26.824250255978107*^-40 + 8.2081320085585782634
[output truncated — 58738 chars total]
Out[2]= {False, Association, 24}
```

## 2026-07-17 19:01:01
**Input:**
```wolfram
worst2 = Max@KeyValueMap[
     Function[{k, rec},
       Module[{H1 = k[[1]], H2 = k[[2]], n = k[[3]], v = rec["Psi"], uu = 17/5},
         Max[
          Max[Abs[(tg[{2,1,0}, {1,0,0}][1, 1][uu] . v) - (rec["tau1c"].uu^Range[0,2]) v]],
          Max[Abs[(tg[{2,1,0}, {1,0,0}][2, 1][uu] . v) - (rec["tau2c"].uu^Range[0,4]) v]]]]],
     res2];
{"nstates"->Length[res2],"worstResidual@u=17/5"->worst2}
```
**Output:**
```
Out= {"nstates" -> 24, "worstResidual@u=17/5" -> 1.48498530494543347294990299017607961167850701950923`25.011916239636957*^-27}
```

## 2026-07-17 19:01:07
**Input:**
```wolfram
(* confirm the interior multiplicity>1 sector(s) actually got exercised, not just trivially all-mult-1 *)
sectorCounts = Counts[Keys[res2][[All,1;;2]]];
Select[sectorCounts, #>1&]
```
**Output:**
```
Out= <|{-2, 0} -> 2, {0, 2} -> 2, {1, 0} -> 4, {0, -1} -> 4, {2, -2} -> 2, {-1, 1} -> 4|>
```

## 2026-07-17 19:01:12
**Input:**
```wolfram
Table[{s, wtMult[{2,1,0},{1,0,0}][s[[1]],s[[2]]]}, {s, Keys[sectorCounts]}]
```
**Output:**
```
Out= {{{-2, 0}, 2}, {{0, 2}, 2}, {{1, 0}, 4}, {{0, -1}, 4}, {{2, -2}, 2}, {{-1, 1}, 4}, {{2, 1}, 1}, {{-1, -2}, 1}, {{3, -1}, 1}, {{1, -3}, 1}, {{-3, 2}, 1}, {{-2, 3}, 1}}
```

## 2026-07-17 19:01:24
**Input:**
```wolfram
check21 = Check[
Module[{λ1={2,1,0},λ2={1,0,0},d,c1,c2,C2,D4,mats,comm},
  d=dim[λ1] dim[λ2];
  {c1,c2}=tgCoeffs[λ1,λ2];
  C2=c1[[3]]; D4=c2[[5]];
  If[Max[Abs[Flatten[C2-χ1 IdentityMatrix[d]]]]>10^-13,Print["C2 FAIL"];Abort[]];
  If[Max[Abs[Flatten[D4-χ2 IdentityMatrix[d]]]]>10^-13,Print["D4 FAIL"];Abort[]];
  mats=Join[c1,c2,{H1op[λ1,λ2],H2op[λ1,λ2]}];
  comm=Max[Table[Max[Abs[Flatten[mats[[i]].mats[[j]]-mats[[j]].mats[[i]]]]],
      {i,Length[mats]},{j,Length[mats]}]];
  If[comm>10^-12,Print["COMMUTATION FAIL: ",comm];Abort[]];
  {"C2=χ1 Id OK","D4=χ2 Id OK","maxCommutator"->comm}
], "ABORTED"]
```
**Output:**
```
Out= {"C2=χ1 Id OK", "D4=χ2 Id OK", "maxCommutator" -> 0``25.701675213462764}
```

## 2026-07-17 19:07:18
**Input:**
```wolfram
{Head[TauEigensystem3], Keys[TauEigensystem3[{2,1,0},{1,0,0}]]}
```
**Output:**
```
Out= {Symbol, {{-2, 0, 0}, {-2, 0, 1}, {0, 2, 0}, {0, 2, 1}, {1, 0, 0}, {1, 0, 1}, {1, 0, 2}, {1, 0, 3}, {0, -1, 0}, {0, -1, 1}, {0, -1, 2}, {0, -1, 3}, {2, -2, 0}, {2, -2, 1}, {-1, 1, 0}, {-1, 1, 1}, {-1, 1, 2}, {-1, 1, 3}, {2, 1, 0}, {-1, -2, 0}, {3, -1, 0}, {1, -3, 0}, {-3, 2, 0}, {-2, 3, 0}}}
```
## 2026-07-17 19:07:59 — ✏️ EDIT CELL 28 [batch]
(*3rd-order TQ residual functional (operator form used for verification below): Qf is a pure
  function of u; residual should vanish on the true Baxter Q-function once M[k] is known.
  Same A1/A2/A3 s…


## 2026-07-17 19:07:59
**Input:**
```wolfram
(*3rd-order TQ residual functional (operator form used for verification below): Qf is a pure
  function of u; residual should vanish on the true Baxter Q-function once M[k] is known.
  Same A1/A2/A3 structure as Bax[λ1,λ2] above but exposed as a callable per-state functional.*)
ClearAll[BaxOp];
BaxOp[λ1_, λ2_][H1_, H2_, n_][Qf_][u_] :=
  Qf[u + h]
   - tau1[λ1, λ2][H1, H2, n][u]/ν1[λ1, λ2][u] Qf[u]
   + tau2[λ1, λ2][H1, H2, n][u]/(ν1[λ1, λ2][u] ν1[λ1, λ2][u - h]) Qf[u - h]
   - qdetT[λ1, λ2][u]/(ν1[λ1, λ2][u] ν1[λ1, λ2][u - h] ν1[λ1, λ2][u - 2 h]) Qf[u - 2 h];
```
**Output:**
```
Out= (no output)
```
## 2026-07-17 19:07:59 — ✏️ EDIT CELL 29 [batch]
(*large-u asymptotics: alphaCoeffs[λ1,λ2][H1,H2,n] = {α1,α2,α3}, the 1/u coefficients of the
  three rational combinations A1=τ1/ν1, A2=τ2/(ν1[u]ν1[u-h]), A3=qdetT/(ν1[u]ν1[u-h]ν1[u-2h])
  that appear…


## 2026-07-17 19:07:59
**Input:**
```wolfram
(*large-u asymptotics: alphaCoeffs[λ1,λ2][H1,H2,n] = {α1,α2,α3}, the 1/u coefficients of the
  three rational combinations A1=τ1/ν1, A2=τ2/(ν1[u]ν1[u-h]), A3=qdetT/(ν1[u]ν1[u-h]ν1[u-2h])
  that appear in BaxOp. These feed the asymptotic master formula for the polynomial degree Mdeg.
  NB (deviation from task-5-brief.md verbatim code, flagged not silently patched): for x0=Infinity,
  SeriesCoefficient[f,{x,Infinity,n}] returns the coefficient of x^(-n), so the 1/u term needs
  n=+1, NOT n=-1. Brief's literal `{u,Infinity,-1}` asks for the (nonexistent, since A1,A2,A3->
  finite constants) u^{+1} term and silently returns 0 for every state -- verified directly:
  SeriesCoefficient[5+7/u-11/u^2,{u,Infinity,-1}]=0 but {u,Infinity,1}=7 (the true 1/u coeff).
  Using n=-1 made Mdeg trivially {0,0,0} everywhere (a false pass, not a NON-INTEGER abort).*)
ClearAll[alphaCoeffs];
alphaCoeffs[λ1_, λ2_][H1_, H2_, n_] := alphaCoeffs[λ1, λ2][H1, H2, n] = Module[{u, A1, A2, A3},
   A1 = tau1[λ1, λ2][H1, H2, n][u]/ν1[λ1, λ2][u];
   A2 = tau2[λ1, λ2][H1, H2, n][u]/(ν1[λ1, λ2][u] ν1[λ1, λ2][u - h]);
   A3 = qdetT[λ1, λ2][u]/(ν1[λ1, λ2][u] ν1[λ1, λ2][u - h] ν1[λ1, λ2][u - 2 h]);
   SeriesCoefficient[#, {u, Infinity, 1}] & /@ {A1, A2, A3}
];
```
**Output:**
```
Out= (no output)
```
## 2026-07-17 19:07:59 — ✏️ EDIT CELL 30 [batch]
(*asymptotic master formula: for each companion-twist eigenvalue z[k] (k=1,2,3), the large-u
  balance of the TQ equation fixes the polynomial degree M[k] of the Q-function branch built on
  z[k]. Den…


## 2026-07-17 19:07:59
**Input:**
```wolfram
(*asymptotic master formula: for each companion-twist eigenvalue z[k] (k=1,2,3), the large-u
  balance of the TQ equation fixes the polynomial degree M[k] of the Q-function branch built on
  z[k]. Denominator h(z[k]-χ2/z[k]+2χ3/z[k]^2) is the z[k]-characteristic factor of z^3-χ1 z^2+χ2 z-χ3
  divided by z^2 (companion twist has eigenvalues z[1],z[2],z[3] with χ1=Σz,χ2=Σzz,χ3=z1z2z3).
  Mdeg is rep-agnostic: only the signature is λ-vector, the body is unchanged from V1.*)
ClearAll[Mdeg];
Mdeg[λ1_, λ2_][H1_, H2_, n_] := Mdeg[λ1, λ2][H1, H2, n] = Module[{a, M},
   a = alphaCoeffs[λ1, λ2][H1, H2, n];
   M = Table[
      (a[[1]] - a[[2]]/z[k] + a[[3]]/z[k]^2)/(h (z[k] - χ2/z[k] + 2 χ3/z[k]^2)),
      {k, 1, 3}];
   If[Max[Abs[M - Round[Re[M]]]] > 10^-9, Print["Mdeg NON-INTEGER: ", M]; Abort[]];
   Round[Re[M]]
];
```
**Output:**
```
Out= (no output)
```
## 2026-07-17 19:07:59 — ✏️ EDIT CELL 36 [batch]
(*polynomial from an ascending coeff list {c0,c1,...,cM}: qpoly[clist,uu]=Sum ck uu^k.
  Dummy-symbol trick (build in x, substitute uu at the end) avoids 0^0 issues when uu=0.*)
ClearAll[qpoly];
qpoly…


## 2026-07-17 19:08:00
**Input:**
```wolfram
(*polynomial from an ascending coeff list {c0,c1,...,cM}: qpoly[clist,uu]=Sum ck uu^k.
  Dummy-symbol trick (build in x, substitute uu at the end) avoids 0^0 issues when uu=0.*)
ClearAll[qpoly];
qpoly[clist_, uu_] := Module[{x}, (clist . x^Range[0, Length[clist] - 1]) /. x -> uu];
```
**Output:**
```
Out= (no output)
```
## 2026-07-17 19:08:00 — ✏️ EDIT CELL 37 [batch]
(*QSolve3[λ1,λ2][H1,H2,n] solves the 3rd-order Baxter (TQ) equation for the three monic
  Q-polynomials q_1,q_2,q_3 (one per companion-twist eigenvalue z[k]), degrees M_k from Mdeg.
  For fixed k, wri…


## 2026-07-17 19:08:00
**Input:**
```wolfram
(*QSolve3[λ1,λ2][H1,H2,n] solves the 3rd-order Baxter (TQ) equation for the three monic
  Q-polynomials q_1,q_2,q_3 (one per companion-twist eigenvalue z[k]), degrees M_k from Mdeg.
  For fixed k, write Q_k(u)=z[k]^(u/h) q_k(u) and substitute into BaxOp[λ1,λ2][H1,H2,n][Q_k][u]=0,
  then multiply through by ν1[u]ν1[u-h]ν1[u-2h] to clear denominators -- this turns the TQ
  equation into a polynomial identity in u, linear in the unknown coefficients of q_k. Each
  monomial q=u^j contributes one column (coefficient list in u); NullSpace of the resulting
  matrix gives the coefficient vector (expected 1-dim nullspace), monic-normalized by dividing
  by the top coefficient. Chop the (high-precision) matrix before NullSpace -- same pitfall as
  su(2)'s QSolve: residual roundoff in the numeric tau1/tau2/qdetT/ν1 data can otherwise mask
  the true nullspace.*)
ClearAll[QSolve3];
QSolve3::degnull = "Q_`4` nullspace dim `5` (expected 1) at (H1,H2,n)=(`1`,`2`,`3`).";
QSolve3::lastzero = "Q_`4` top coeff ~0 at (H1,H2,n)=(`1`,`2`,`3`).";
QSolve3[λ1_, λ2_][H1_, H2_, n_] := QSolve3[λ1, λ2][H1, H2, n] = Module[
   {Ms, u, T1, T2, QD, N1, N1m, N1mm, sols},
   Ms = Mdeg[λ1, λ2][H1, H2, n];
   T1 = tau1[λ1, λ2][H1, H2, n][u]; T2 = tau2[λ1, λ2][H1, H2, n][u];
   QD = qdetT[λ1, λ2][u];
   N1 = ν1[λ1, λ2][u]; N1m = ν1[λ1, λ2][u - h]; N1mm = ν1[λ1, λ2][u - 2 h];
   sols = Table[
     Module[{M = Ms[[k]], zk = z[k], cols, mat, ns, vec, last, maxdeg},
       maxdeg = M + 6;
       (* column j (j=0..M): coefficient list of the polynomial identity for q=u^j *)
       cols = Table[
          PadRight[CoefficientList[Expand[
             zk (u + h)^j N1 N1m N1mm
              - T1 u^j N1m N1mm
              + T2/zk (u - h)^j N1mm
              - QD/zk^2 (u - 2 h)^j], u], maxdeg + 1],
          {j, 0, M}];
       mat = Transpose[cols];                 (* (maxdeg+1) x (M+1) *)
       ns = NullSpace[Chop[mat, 10^-12], Tolerance -> 10^-10];
       If[Length[ns] != 1, Message[QSolve3::degnull, H1, H2, n, k, Length[ns]]; Abort[]];
       vec = ns[[1]]; last = vec[[-1]];
       If[Abs[last] < 10^-13 Max[Abs[vec]], Message[QSolve3::lastzero, H1, H2, n, k]; Abort[]];
       vec/last],
     {k, 1, 3}];
   sols
];
```
**Output:**
```
Out= (no output)
```
## 2026-07-17 19:08:00 — ✏️ EDIT CELL 38 [batch]
(*Q accessor: Q_k(u)=z[k]^(u/h) q_k(u), q_k the monic coeff list from QSolve3.*)
ClearAll[Qfun];
Qfun[λ1_, λ2_][H1_, H2_, n_][k_][u_] := z[k]^(u/h) qpoly[QSolve3[λ1, λ2][H1, H2, n][[k]], u];


## 2026-07-17 19:08:00
**Input:**
```wolfram
(*Q accessor: Q_k(u)=z[k]^(u/h) q_k(u), q_k the monic coeff list from QSolve3.*)
ClearAll[Qfun];
Qfun[λ1_, λ2_][H1_, H2_, n_][k_][u_] := z[k]^(u/h) qpoly[QSolve3[λ1, λ2][H1, H2, n][[k]], u];
```
**Output:**
```
Out= (no output)
```
## 2026-07-17 19:08:13 — ✏️ EDIT CELL 39
(*single-state TQ-residual check on a non-symmetric rep: each Q_k (built from QSolve3) must be
  annihilated by BaxOp, tested at several generic numeric u. State: λ1={2,1,0}, λ2={1,0,0}
  (adjoint⊗fun…


## 2026-07-17 19:08:13
**Input:**
```wolfram
(*single-state TQ-residual check on a non-symmetric rep: each Q_k (built from QSolve3) must be
  annihilated by BaxOp, tested at several generic numeric u. State: λ1={2,1,0}, λ2={1,0,0}
  (adjoint⊗fund), (H1,H2,n)=(1,0,0) (confirmed a valid n=0 key of TauEigensystem3[{2,1,0},{1,0,0}]).*)
Module[{λ1 = {2, 1, 0}, λ2 = {1, 0, 0}, H1 = 1, H2 = 0, n = 0, res},
  res = Table[
     Max@Table[
        Abs@Chop[N[BaxOp[λ1, λ2][H1, H2, n][Function[uu, Qfun[λ1, λ2][H1, H2, n][k][uu]]][u0], 30]],
        {u0, {2/3, 9/5, 4, 23/7}}],
     {k, 1, 3}];
  If[Max[res] > 10^-10, Print["TQ RESIDUAL FAIL: ", res]; Abort[]];
  {"perQresidual" -> res}
]
```
**Output:**
```
Out= {"perQresidual" -> {0, 0, 0}}
```
## 2026-07-17 19:08:23 — ✏️ EDIT CELL 31
(*collect Mdeg over the whole sweep as {λ1,λ2,H1,H2,n,{M1,M2,M3}} rows for closed-form fitting.
  Sanity gate: every {M1,M2,M3} must be a triple of non-negative integers (Mdeg already asserts
  near-i…


## 2026-07-17 19:08:24
**Input:**
```wolfram
(*collect Mdeg over the whole sweep as {λ1,λ2,H1,H2,n,{M1,M2,M3}} rows for closed-form fitting.
  Sanity gate: every {M1,M2,M3} must be a triple of non-negative integers (Mdeg already asserts
  near-integer). M_k depends only on (λ1,λ2,H1,H2), NOT on the branch index n (verified below).*)
degData=Flatten[Table[
   Module[{λ1=rep[[1]],λ2=rep[[2]],es},
     es=TauEigensystem3[λ1,λ2];
     KeyValueMap[Function[{key,rec},
        {λ1,λ2,key[[1]],key[[2]],key[[3]],Mdeg[λ1,λ2][key[[1]],key[[2]],key[[3]]]}],es]],
   {rep,{{{1,0,0},{1,0,0}},{{2,0,0},{1,1,0}},{{1,1,0},{1,1,0}},
         {{2,1,0},{1,0,0}},{{2,1,0},{1,1,0}},{{2,2,0},{1,0,0}}}}],1];
Column[degData]
```
**Output:**
```
Out= Column[{{{1, 0, 0}, {1, 0, 0}, 0, 1, 0, {1, 1, 2}}, {{1, 0, 0}, {1, 0, 0}, 0, 1, 1, {1, 1, 2}}, {{1, 0, 0}, {1, 0, 0}, 1, -1, 0, {1, 2, 1}}, {{1, 0, 0}, {1, 0, 0}, 1, -1, 1, {1, 2, 1}}, {{1, 0, 0}, {1, 0, 0}, -1, 0, 0, {2, 1, 1}}, {{1, 0, 0}, {1, 0, 0}, -1, 0, 1, {2, 1, 1}}, {{1, 0, 0}, {1, 0, 0}, 2, 0, 0, {0, 2, 2}}, {{1, 0, 0}, {1, 0, 0}, 0, -2, 0, {2, 2, 0}}, {{1, 0, 0}, {1, 0, 0}, -2, 2, 0, {2, 0, 2}}, {{2, 0, 0}, {1, 1, 0}, 0, -1, 0, {2, 2, 1}}, {{2, 0, 0}, {1, 1, 0}, 0, -1, 1, {2, 2, 1}}, {{2, 0, 0}, {1, 1, 0}, 0, -1, 2, {2, 2, 1}}, {{2, 0, 0}, {1, 1, 0}, -1, 1, 0, {2, 1, 2}}, {{2, 0, 0}, {1, 1, 0}, -1, 1, 1, {2, 1, 2}}, {{2, 0, 0}, {1, 1, 0}, -1, 1, 2, {2, 1, 2}}, {{2, 0, 0}, {1, 1, 0}, 1, 0, 0, {1, 2, 2}}, {{2, 0, 0}, {1, 1, 0}, 1, 0, 1, {1, 2, 2}}, {{2, 0, 0}, {1, 1, 0}, 1, 0, 2, {1, 2, 2}}, {{2, 0, 0}, {1, 1, 0}, -1, -2, 0, {3, 2, 0}}, {{2, 0, 0}, {1, 1, 0}, 2, 1, 0, {0, 2, 3}}, {{2, 0, 0}, {1, 1, 0}, -2, 0, 0, {3, 1, 1}}, {{2, 0, 0}, {1, 1, 0}, 0, 2, 0, {1, 1, 3}}, {{2, 0, 0}, {1, 1, 0}, 3, -1, 0, {0, 3, 2}}, {{2, 0, 0}, {1, 1, 0}, 1, -3, 0, {2, 3, 0}}, {{2, 0, 0}, {1, 1, 0}, -3, 2, 0, {3, 0, 2}}, {{2, 0, 0}, {1, 1, 0}, -2, 3, 0, {2, 0, 3}}, {{2, 0, 0}, {1, 1, 0}, 2, -2, 0, {1, 3, 1}}, {{1, 1, 0}, {1, 1, 0}, 0, -1, 0, {1, 1, 0}}, {{1, 1, 0}, {1, 1, 0}, 0, -1, 1, {1, 1, 0}}, {{1, 1, 0}, {1, 1, 0}, -1, 1, 0, {1, 0, 1}}, {{1, 1, 0}, {1, 1, 0}, -1, 1, 1, {1, 0, 1}}, {{1, 1, 0}, {1, 1, 0}, 1, 0, 0, {0, 1, 1}}, {{1, 1, 0}, {1, 1, 0}, 1, 0, 1, {0, 1, 1}}, {{1, 1, 0}, {1, 1, 0}, -2, 0, 0, {2, 0, 0}}, {{1, 1, 0}, {1, 1, 0}, 0, 2, 0, {0, 0, 2}}, {{1, 1, 0}, {1, 1, 0}, 2, -2, 0, {0, 2, 0}}, {{2, 1, 0}, {1, 0, 0}, -2, 0, 0, {3, 1, 1}}, {{2, 1, 0}, {1, 0, 0}, -2, 0, 1, {3, 1, 1}}, {{2, 1, 0}, {1, 0, 0}, 0, 2, 0, {1, 1, 3}}, {{2, 1, 0}, {1, 0, 0}, 0, 2, 1, {1, 1, 3}}, {{2, 1, 0}, {1, 0, 0}, 1, 0, 0, {1, 2, 2}}, {{2, 1, 0}, {1, 0, 0}, 1, 0, 1, {1, 2, 2}}, {{2, 1, 0}, {1, 0, 0}, 1, 0, 2, {1, 2, 2}}, {{2, 1, 0}, {1, 0, 0}, 1, 0, 3, {1, 2, 2}}, {{2, 1, 0}, {1, 0, 0}, 0, -1, 0, {2, 2, 1}}, {{2, 1, 0}, {1, 0, 0}, 0, -1, 1, {2, 2, 1}}, {{2, 1, 0}, {1, 0, 0}, 0, -1, 2, {2, 2, 1}}, {{2, 1, 0}, {1, 0, 0}, 0, -1, 3, {2, 2, 1}}, {{2, 1, 0}, {1, 0, 0}, 2, -2, 0, {1, 3, 1}}, {{2, 1, 0}, {1, 0, 0}, 2, -2, 1, {1, 3, 1}}, {{2, 1, 0}, {1, 0, 0}, -1, 1, 0, {2, 1, 2}}, {{2, 1, 0}, {1, 0, 0}, -1, 1, 1, {2, 1, 2}}, {{2, 1, 0}, {1, 0, 0}, -1, 1, 2, {2, 1, 2}}, {{2, 1, 0}, {1, 0, 0}, -1, 1, 3, {2, 1, 2}}, {{2, 1, 0}, {1, 0, 0}, 2, 1, 0, {0, 2, 3}}, {{2, 1, 0}, {1, 0, 0}, -1, -2, 0, {3, 2, 0}}, {{2, 1, 0}, {1, 0, 0}, 3, -1, 0, {0, 3, 2}}, {{2, 1, 0}, {1, 0, 0}, 1, -3, 0, {2, 3, 0}}, {{2, 1, 0}, {1, 0, 0}, -3, 2, 0, {3, 0, 2}}, {{2, 1, 0}, {1, 0, 0}, -2, 3, 0, {2, 0, 3}}, {{2, 1, 0}, {1, 1, 0}, 2, 0, 0, {0, 2, 2}}, {{2, 1, 0}, {1, 1, 0}, 2, 0, 1, {0, 2, 2}}, {{2, 1, 0}, {1, 1, 0}, 0, 1, 0, {1, 1, 2}}, {{2, 1, 0}, {1, 1, 0}, 0, 1, 1, {1, 1, 2}}, {{2, 1, 0}, {1, 1, 0}, 0, 1, 2, {1, 1, 2}}, {{2, 1, 0}, {1, 1, 0}, 0, 1, 3, {1, 1, 2}}, {{2, 1, 0}, {1, 1, 0}, -2, 2, 0, {2, 0, 2}}, {{2, 1, 0}, {1, 1, 0}, -2, 2, 1, {2, 0, 2}}, {{2, 1, 0}, {1, 1, 0}, 0, -2, 0, {2, 2, 0}}, {{2, 1, 0}, {1, 1, 0}, 0, -2, 1, {2, 2, 0}}, {{2, 1, 0}, {1, 1, 0}, -1, 0, 0, {2, 1, 1}}, {{2, 1, 0}, {1, 1, 0}, -1, 0, 1, {2, 1, 1}}, {{2, 1, 0}, {1, 1, 0}, -1, 0, 2, {2, 1, 1}}, {{2, 1, 0}, {1, 1, 0}, -1, 0, 3, {2, 1, 1}}, {{2, 1, 0}, {1, 1, 0}, 1, -1, 0, {1, 2, 1}}, {{2, 1, 0}, {1, 1, 0}, 1, -1, 1, {1, 2, 1}}, {{2, 1, 0}, {1, 1, 0}, 1, -1, 2, {1, 2, 1}}, {{2, 1, 0}, {1, 1, 0}, 1, -1, 3, {1, 2, 1}}, {{2, 1, 0}, {1, 1, 0}, -2, -1, 0, {3, 1, 0}}, {{2, 1, 0}, {1, 1, 0}, -3, 1, 0, {3, 0, 1}}, {{2, 1, 0}, {1, 1, 0}, 1, 2, 0, {0, 1, 3}}, {{2, 1, 0}, {1, 1, 0}, -1, 3, 0, {1, 0, 3}}, {{2, 1, 0}, {1, 1, 0}, 2, -3, 0, {1, 3, 0}}, {{2, 1, 0}, {1, 1, 0}, 3, -2, 0, {0, 3, 1}}, {{2, 2, 0}, {1, 0, 0}, 1, 2, 0, {0, 1, 3}}, {{2, 2, 0}, {1, 0, 0}, -1, 3, 0, {1, 0, 3}}, {{2, 2, 0}, {1, 0, 0}, -2, -1, 0, {3, 1, 0}}, {{2, 2, 0}, {1, 0, 0}, 2, 0, 0, {0, 2, 2}}, {{2, 2, 0}, {1, 0, 0}, -1, 0, 0, {2, 1, 1}}, {{2, 2, 0}, {1, 0, 0}, -1, 0, 1, {2, 1, 1}}, {{2, 2, 0}, {1, 0, 0}, -1, 0, 2, {2, 1, 1}}, {{2, 2, 0}, {1, 
[output truncated — 4572 chars total]
```
## 2026-07-17 19:08:57 — ✏️ EDIT CELL 32
(*test candidate closed forms for M_k against degData. n_k from (H1,H2) with N=total boxes.
  Candidates (each gives M_k in index k, matching z[k]/n_k pairing):
    cand["N-n"]     : M_k = N - n_k,   …


## 2026-07-17 19:08:57
**Input:**
```wolfram
(*test candidate closed forms for M_k against degData. n_k from (H1,H2) with N=total boxes.
  Candidates (each gives M_k in index k, matching z[k]/n_k pairing):
    cand["N-n"]     : M_k = N - n_k,          N = total boxes = Total[λ1]+Total[λ2]
    cand["Lam1-n"]  : M_k = Λ_1 - n_k,        Λ = λ1 + λ2 (pseudovacuum hw), Λ_1 = (λ1+λ2)[[1]]
    cand["Lamk-n"]  : M_k = Λ_k - n_k,        per-index Λ component (partial-sum-family probe)
  For symmetric reps Λ=(N,0,0) so "N-n" and "Lam1-n" coincide; the non-symmetric sweep
  distinguishes them. WINNER (zero max-deviation over ALL rows): cand["Lam1-n"], i.e.
  M_k = (λ1+λ2)[[1]] - n_k.  (sum rule then reads Σ_k M_k = 3 Λ_1 - N, not 2N.)*)
nvec[H1_,H2_,Ntot_]:=Module[{n3=(Ntot-H1-2 H2)/3},{n3+H1+H2,n3+H2,n3}];
cand["N-n"][λ1_,λ2_,H1_,H2_]:=(Total[λ1]+Total[λ2])-nvec[H1,H2,Total[λ1]+Total[λ2]];
cand["Lam1-n"][λ1_,λ2_,H1_,H2_]:=(λ1+λ2)[[1]]-nvec[H1,H2,Total[λ1]+Total[λ2]];
cand["Lamk-n"][λ1_,λ2_,H1_,H2_]:=(λ1+λ2)-nvec[H1,H2,Total[λ1]+Total[λ2]];
results=Association@Table[
   name->Max[Abs[Flatten[
      (cand[name][#[[1]],#[[2]],#[[3]],#[[4]]]-#[[6]])&/@degData]]],
   {name,{"N-n","Lam1-n","Lamk-n"}}];
results
```
**Output:**
```
Out= <|"N-n" -> 2, "Lam1-n" -> 0, "Lamk-n" -> 3|>
```
## 2026-07-17 19:09:12 — ✏️ EDIT CELL 33 [batch]
(* superseded by degData tabulation + closed-form fit above *)


## 2026-07-17 19:09:12
**Input:**
```wolfram
(* superseded by degData tabulation + closed-form fit above *)
```
**Output:**
```
Out= Out[97]= \begin{matrix}\left\{\text{H1H2n}\to \left\{-4,4,0\right\},\text{M}\to \left\{4,0,4\right\},\text{sumM}\to 8\right\}\\\left\{\text{H1H2n}\to \left\{-3,2,0\right\},\text{M}\to \left\{4,1,3\right\},\text{sumM}\to 8\right\}\\\left\{\text{H1H2n}\to \left\{-3,2,1\right\},\text{M}\to \left\{4,1,3\right\},\text{sumM}\to 8\right\}\\\left\{\text{H1H2n}\to \left\{-2,0,0\right\},\text{M}\to \left\{4… [1246 chars total]
```
## 2026-07-17 19:09:12 — ✏️ EDIT CELL 34 [batch]
(* superseded by degData tabulation + closed-form fit above *)


## 2026-07-17 19:09:13
**Input:**
```wolfram
(* superseded by degData tabulation + closed-form fit above *)
```
**Output:**
```
Out= Out[98]= \begin{matrix}\left\{\text{H1H2n}\to \left\{-4,4,0\right\},\text{M}\to \left\{4,0,4\right\},\text{sumM}\to 8\right\}\\\left\{\text{H1H2n}\to \left\{-3,2,0\right\},\text{M}\to \left\{4,1,3\right\},\text{sumM}\to 8\right\}\\\left\{\text{H1H2n}\to \left\{-3,2,1\right\},\text{M}\to \left\{4,1,3\right\},\text{sumM}\to 8\right\}\\\left\{\text{H1H2n}\to \left\{-2,0,0\right\},\text{M}\to \left\{4… [1246 chars total]
```
## 2026-07-17 19:09:13 — ✏️ EDIT CELL 35 [batch]
(*hard-assert: extracted closed form for M_k matches Mdeg exactly (zero residual) over the sweep.
  WINNER (fit above): M_k = Λ_1 - n_k, with Λ = λ1+λ2 (pseudovacuum highest weight), Λ_1=(λ1+λ2)[[1]],…


## 2026-07-17 19:09:13
**Input:**
```wolfram
(*hard-assert: extracted closed form for M_k matches Mdeg exactly (zero residual) over the sweep.
  WINNER (fit above): M_k = Λ_1 - n_k, with Λ = λ1+λ2 (pseudovacuum highest weight), Λ_1=(λ1+λ2)[[1]],
  and (n1,n2,n3) the weight components: n3=(N-H1-2H2)/3, n2=n3+H2, n1=n3+H1+H2, N=Total[λ1]+Total[λ2].
  For symmetric reps λ=(S,0,0) this reduces to V1's M_k=N-n_k since then Λ_1=N.*)
nvec[H1_,H2_,Ntot_]:=Module[{n3=(Ntot-H1-2 H2)/3},{n3+H1+H2,n3+H2,n3}];
Mpred[λ1_,λ2_][H1_,H2_]:=(λ1+λ2)[[1]]-nvec[H1,H2,Total[λ1]+Total[λ2]];
Module[{sweep={{{1,0,0},{1,0,0}},{{2,0,0},{1,1,0}},{{1,1,0},{1,1,0}},
               {{2,1,0},{1,0,0}},{{2,1,0},{1,1,0}},{{2,2,0},{1,0,0}}},worst},
  worst=Max@Table[
     Module[{λ1=rep[[1]],λ2=rep[[2]],es},
       es=TauEigensystem3[λ1,λ2];
       Max@KeyValueMap[Function[{key,rec},
          Max[Abs[Mdeg[λ1,λ2][key[[1]],key[[2]],key[[3]]]-Mpred[λ1,λ2][key[[1]],key[[2]]]]]],es]],
     {rep,sweep}];
  If[worst>0,Print["DEGREE-FORMULA MISMATCH: ",worst];Abort[]];
  {"closed-form M_k=Λ_1-n_k matches Mdeg over full sweep, worst"->worst}
]
```
**Output:**
```
Out= Out[99]= \left\{\text{worst |Mdeg - (N-n\_k formula)| over (2,1),(1,1),(2,2),(1,3)}\to 0\right\}
```
## 2026-07-17 19:09:21 — ▶️ RUN CELL 35
**In [35]:** `(*hard-assert: extracted closed form for M_k matches Mdeg exactly (zero residual) over the sweep.↵  WINNER (fit above): M_k = Λ_1 - n_k, with Λ = λ1+λ2 (pseudovacuum highest weight), Λ_1=(λ1+λ2)[[1]],`
**Out:** Out[99]= \left\{\text{worst |Mdeg - (N-n\_k formula)| over (2,1),(1,1),(2,2),(1,3)}\to 0\right\}


## 2026-07-17 19:09:33
**Input:**
```wolfram
nvec[H1_,H2_,Ntot_]:=Module[{n3=(Ntot-H1-2 H2)/3},{n3+H1+H2,n3+H2,n3}];
Mpred[λ1_,λ2_][H1_,H2_]:=(λ1+λ2)[[1]]-nvec[H1,H2,Total[λ1]+Total[λ2]];
Module[{sweep={{{1,0,0},{1,0,0}},{{2,0,0},{1,1,0}},{{1,1,0},{1,1,0}},{{2,1,0},{1,0,0}},{{2,1,0},{1,1,0}},{{2,2,0},{1,0,0}}},worst},worst=Max@Table[Module[{λ1=rep[[1]],λ2=rep[[2]],es},es=TauEigensystem3[λ1,λ2];Max@KeyValueMap[Function[{key,rec},Max[Abs[Mdeg[λ1,λ2][key[[1]],key[[2]],key[[3]]]-Mpred[λ1,λ2][key[[1]],key[[2]]]]]],es]],{rep,sweep}];{"worst"->worst}]
```
**Output:**
```
Out= {"worst" -> 0}
```
## 2026-07-17 19:18:01 — ✏️ EDIT CELL 40
(*full sweep verification: over the six general gl(3) rep pairs (λ1,λ2) and every
  {H1,H2,n} state, check (a) per-Q_k TQ residual and (b) the 4x4 Casoratian reconstruction --
  the three Q_k span the…


---

## 2026-07-17 19:18:07 — 🔄 KERNEL RESTART

## 2026-07-17 19:18:21 — ▶️ RUN CELLS 1–14
Cell 1: ✓ — (no output)
Cell 2: ✓ — (no output)
Cell 3: ✓ — (no output)
Cell 4: ✓ — Out[14]= \left\{0\right\}
Cell 5: ✓ — (no output)
Cell 6: ✓ — (no output)
Cell 7: ✓ — Out[27]= \left\{\frac{1}{3},\frac{1}{7},1,0.360439875147226537266677650596+0.932782448593374716484835982071\,\mathrm{i},0.508867380990499828163321372267+0.860844927129079093188644393616\,\mathrm{i},0.533257880589599372595832345093+0.845952736734794526483381618453\,\mathrm{i},1.40256513672732573802583136796+2.63958011245724833615686199414\,\mathrm{i},-1.6733239859016290991182472189+2.47680418504351… [511 chars total]
Cell 8: ✓ — (no output)
Cell 9: ✓ — (no output)
Cell 10: ✓ — Out[37]= \left\{0\right\}
Cell 11: ✓ — (no output)
Cell 12: ✓ — (no output)
Cell 13: ✓ — (no output)
Cell 14: ✓ — (no output)

## 2026-07-17 19:18:34 — ▶️ RUN CELLS 15–27
Cell 15: ✓ — Out[49]= \left\{29.998415624999733,-0.9944145978415286-0.10554433949611333\,\mathrm{i},0\right\}
Cell 16: ✓ — Out[52]= \left(u-\frac{4}{3}\right)\,\left(u-\frac{8}{7}\right)
Cell 17: ✓ — Out[124]= -h\,u\,\left(z\left(1\right)\,\mathrm{ee}\left(S1,S2\right)\left(1,1\right)+z\left(2\right)\,\mathrm{ee}\left(S1,S2\right)\left(2,2\right)+z\left(3\right)\,\mathrm{ee}\left(S1,S2\right)\left(3,3\right)\right)+u^2\,\left(z\left(1\right)+z\left(2\right)+z\left(3\right)\right)-u\,\left(z\left(1\right)+z\left(2\right)+z\left(3\right)\right)\,\left(\theta1+\theta2\right)
Cell 18: ✓ — Out[57]= \mathrm{List}
Cell 19: ✓ — (no output)
Cell 20: ✓ — (no output)
Cell 21: ✓ — Out[63]= \left\{\text{C2=$\chi$1 Id OK},\text{D4=$\chi$2 Id OK},\text{maxCommutator}\to 0.\times 10^{-26}\right\}
Cell 22: ✓ — (no output)
Cell 23: ✓ — (no output)
Cell 24: ✓ — (no output)
Cell 25: ✓ — (no output)
Cell 26: ✓ — (no output)
Cell 27: ✓ — (no output)

## 2026-07-17 19:18:49 — ▶️ RUN CELLS 28–40
Cell 28: ✓ — (no output)
Cell 29: ✓ — (no output)
Cell 30: ✓ — (no output)
Cell 31: ✓ — Out[95]= \begin{matrix}\left\{\text{H1H2n}\to \left\{-3,3,0\right\},\text{M}\to \left\{3,0,3\right\},\text{sumM}\to 6\right\}\\\left\{\text{H1H2n}\to \left\{-2,1,0\right\},\text{M}\to \left\{3,1,2\right\},\text{sumM}\to 6\right\}\\\left\{\text{H1H2n}\to \left\{-2,1,1\right\},\text{M}\to \left\{3,1,2\right\},\text{sumM}\to 6\right\}\\\left\{\text{H1H2n}\to \left\{-1,-1,0\right\},\text{M}\to \left\{… [1246 chars total]
Cell 32: ✓ — Out[96]= \begin{matrix}\left\{\text{H1H2n}\to \left\{-2,2,0\right\},\text{M}\to \left\{2,0,2\right\},\text{sumM}\to 4\right\}\\\left\{\text{H1H2n}\to \left\{-1,0,0\right\},\text{M}\to \left\{2,1,1\right\},\text{sumM}\to 4\right\}\\\left\{\text{H1H2n}\to \left\{-1,0,1\right\},\text{M}\to \left\{2,1,1\right\},\text{sumM}\to 4\right\}\\\left\{\text{H1H2n}\to \left\{0,-2,0\right\},\text{M}\to \left\{2… [966 chars total]
Cell 33: ✓ — Out[97]= \begin{matrix}\left\{\text{H1H2n}\to \left\{-4,4,0\right\},\text{M}\to \left\{4,0,4\right\},\text{sumM}\to 8\right\}\\\left\{\text{H1H2n}\to \left\{-3,2,0\right\},\text{M}\to \left\{4,1,3\right\},\text{sumM}\to 8\right\}\\\left\{\text{H1H2n}\to \left\{-3,2,1\right\},\text{M}\to \left\{4,1,3\right\},\text{sumM}\to 8\right\}\\\left\{\text{H1H2n}\to \left\{-2,0,0\right\},\text{M}\to \left\{4… [1246 chars total]
Cell 34: ✓ — Out[98]= \begin{matrix}\left\{\text{H1H2n}\to \left\{-4,4,0\right\},\text{M}\to \left\{4,0,4\right\},\text{sumM}\to 8\right\}\\\left\{\text{H1H2n}\to \left\{-3,2,0\right\},\text{M}\to \left\{4,1,3\right\},\text{sumM}\to 8\right\}\\\left\{\text{H1H2n}\to \left\{-3,2,1\right\},\text{M}\to \left\{4,1,3\right\},\text{sumM}\to 8\right\}\\\left\{\text{H1H2n}\to \left\{-2,0,0\right\},\text{M}\to \left\{4… [1246 chars total]
Cell 35: ✓ — Out[99]= \left\{\text{worst |Mdeg - (N-n\_k formula)| over (2,1),(1,1),(2,2),(1,3)}\to 0\right\}
Cell 36: ✓ — (no output)
Cell 37: ✓ — (no output)
Cell 38: ✓ — (no output)
Cell 39: ✓ — Out[108]= \left\{\text{perQresidual}\to \left\{0,0,0\right\}\right\}
Cell 40: ✓ — Out[109]= \begin{matrix}\left\{\left\{1,1\right\},\text{nstates}\to 9,\text{worstTQ}\to 0,\text{worstCasoratian}\to 0\right\}\\\left\{\left\{2,1\right\},\text{nstates}\to 18,\text{worstTQ}\to 0,\text{worstCasoratian}\to 0\right\}\\\left\{\left\{1,2\right\},\text{nstates}\to 18,\text{worstTQ}\to 0,\text{worstCasoratian}\to 0\right\}\\\left\{\left\{2,2\right\},\text{nstates}\to 36,\text{worstTQ}\to … [657 chars total]

## 2026-07-17 19:19:19 — ▶️ RUN CELL 40
**In [40]:** `(*full sweep verification: over the six general gl(3) rep pairs (λ1,λ2) and every↵  {H1,H2,n} state, check (a) per-Q_k TQ residual and (b) the 4x4 Casoratian reconstruction --↵  the three Q_k span the`
**Out:** Out[109]= \begin{matrix}\left\{\left\{1,1\right\},\text{nstates}\to 9,\text{worstTQ}\to 0,\text{worstCasoratian}\to 0\right\}\\\left\{\left\{2,1\right\},\text{nstates}\to 18,\text{worstTQ}\to 0,\text{worstCasoratian}\to 0\right\}\\\left\{\left\{1,2\right\},\text{nstates}\to 18,\text{worstTQ}\to 0,\t


## 2026-07-17 19:19:35
**Input:**
```wolfram
Module[{sweep = {{{1,0,0},{1,0,0}},{{2,0,0},{1,1,0}},{{1,1,0},{1,1,0}},{{2,1,0},{1,0,0}},{{2,1,0},{1,1,0}},{{2,2,0},{1,0,0}}}, report},
  report = Table[
    Module[{λ1 = rep[[1]], λ2 = rep[[2]], es, worstTQ = 0, worstCas = 0},
      es = TauEigensystem3[λ1, λ2];
      KeyValueMap[Function[{key, rec},
        Module[{H1 = key[[1]], H2 = key[[2]], n = key[[3]], Qs, u0, rows, minors, A1rec, A2rec, A3rec, A1, A2, A3, N1, N1m, N1mm},
         Do[worstTQ = Max[worstTQ, Abs@Chop[N[BaxOp[λ1,λ2][H1,H2,n][Function[uu, Qfun[λ1,λ2][H1,H2,n][k][uu]]][u0],30]]], {k,1,3},{u0,{5/4,11/3}}];
         u0 = 13/5;
         Qs[uu_] := Table[Qfun[λ1,λ2][H1,H2,n][k][uu],{k,1,3}];
         rows = {Qs[u0+h],Qs[u0],Qs[u0-h],Qs[u0-2 h]};
         minors = Table[Det[Delete[rows,i]],{i,1,4}];
         A1rec = minors[[2]]/minors[[1]]; A2rec = minors[[3]]/minors[[1]]; A3rec = minors[[4]]/minors[[1]];
         N1 = ν1[λ1,λ2][u0]; N1m = ν1[λ1,λ2][u0-h]; N1mm = ν1[λ1,λ2][u0-2 h];
         A1 = tau1[λ1,λ2][H1,H2,n][u0]/N1;
         A2 = tau2[λ1,λ2][H1,H2,n][u0]/(N1 N1m);
         A3 = qdetT[λ1,λ2][u0]/(N1 N1m N1mm);
         worstCas = Max[worstCas, Abs@Chop[N[A1rec-A1,30]], Abs@Chop[N[A2rec-A2,30]], Abs@Chop[N[A3rec-A3,30]]]]], es];
      {rep, Length[es], worstTQ, worstCas}],
    {rep, sweep}];
  {report, "totalStates" -> Total[report[[All,2]]], "worstTQall" -> Max[report[[All,3]]], "worstCasAll" -> Max[report[[All,4]]]}]
```
**Output:**
```
Out= {{{{{1, 0, 0}, {1, 0, 0}}, 9, 0, 0}, {{{2, 0, 0}, {1, 1, 0}}, 18, 0, 0}, {{{1, 1, 0}, {1, 1, 0}}, 9, 0, 0}, {{{2, 1, 0}, {1, 0, 0}}, 24, 0, 0}, {{{2, 1, 0}, {1, 1, 0}}, 24, 0, 0}, {{{2, 2, 0}, {1, 0, 0}}, 18, 0, 0}}, "totalStates" -> 102, "worstTQall" -> 0, "worstCasAll" -> 0}
```

## 2026-07-17 19:20:28
**Input:**
```wolfram
{"c14qdet" -> Module[{reps={{{2,1,0},{1,0,0}},{{1,1,0},{1,1,0}},{{2,2,0},{1,0,0}},{{2,1,0},{1,1,0}}}}, Union@Flatten@Table[Chop[Simplify[qdetT[r[[1]],r[[2]]][u] idd[r[[1]],r[[2]]]-tg[r[[1]],r[[2]]][3,1][u]],10^-15],{r,reps}]],
"c23sumrule" -> Module[{reps={{{2,1,0},{1,0,0}},{{2,0,0},{1,1,0}},{{1,1,0},{1,1,0}},{{2,1,0},{1,1,0}},{{2,2,0},{1,0,0}}},bad}, bad=Select[reps,Function[r, Module[{Ntot=Total[r[[1]]]+Total[r[[2]]],d,tot}, d=dim[r[[1]]] dim[r[[2]]]; tot=Total@Flatten@Table[wtMult[r[[1]],r[[2]]][H1,H2],{H1,-3 Ntot,3 Ntot},{H2,-3 Ntot,3 Ntot}]; tot=!=d]]]; {"badReps"->bad,"adjFundMult"->wtMult[{2,1,0},{1,0,0}][0,0]}],
"c24guard" -> Module[{λ1={2,1,0},λ2={2,1,0},Ntot,d,tot}, Ntot=Total[λ1]+Total[λ2]; d=dim[λ1] dim[λ2]; tot=Total@Flatten@Table[wtMult[λ1,λ2][H1,H2],{H1,-3 Ntot,3 Ntot},{H2,-3 Ntot,3 Ntot}]; {"tot"->tot,"dim"->d}],
"c27eigen" -> Module[{λ1={2,1,0},λ2={1,0,0},es,worst}, es=TauEigensystem3[λ1,λ2]; worst=Max@KeyValueMap[Function[{k,rec}, Module[{H1=k[[1]],H2=k[[2]],n=k[[3]],v=rec["Psi"],uu=17/5}, Max[Max[Abs[(tg[λ1,λ2][1,1][uu].v)-tau1[λ1,λ2][H1,H2,n][uu] v]], Max[Abs[(tg[λ1,λ2][2,1][uu].v)-tau2[λ1,λ2][H1,H2,n][uu] v]]]]],es]; {"nstates"->Length[es],"worst"->worst}]}
```
**Output:**
```
Out= {"c14qdet" -> {0}, "c23sumrule" -> {"badReps" -> {}, "adjFundMult" -> 0}, "c24guard" -> {"tot" -> 64, "dim" -> 64}, "c27eigen" -> {"nstates" -> 24, "worst" -> 1.48498530494543347294990299017607961167850701950923`25.011916239636957*^-27}}
```
## 2026-07-17 19:21:48 — ▶️ RUN CELL 24
**In [24]:** `(*guard on {2,1,0}⊗{2,1,0} (d=64) — a rep with genuine interior multiplicity on both sites*)↵Module[{λ1={2,1,0},λ2={2,1,0},Ntot,d,tot},↵  Ntot=Total[λ1]+Total[λ2]; d=dim[λ1] dim[λ2];↵  tot=Total@Flatt`
**Out:** (no output)

## 2026-07-17 19:21:48 — ▶️ RUN CELL 14
**In [14]:** `(*qdetT central value == tg[·][3,1]; on non-symmetric reps this is the ν shift-ordering guard.↵  NB: Chop wraps Simplify because z[k] (hence χ3) are high-precision *numeric* (not exact↵  symbolic) her`
**Out:** (no output)


---

## 2026-07-17 19:30:25 — 🔄 KERNEL RESTART

