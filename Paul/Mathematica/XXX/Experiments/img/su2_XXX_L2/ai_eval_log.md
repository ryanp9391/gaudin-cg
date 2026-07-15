# AI Action Log — su2_XXX_L2

## 2026-07-14 13:45:23 — 📥 BULK INSERT 5 CELL(S) at positions 1–5
1. [markdown] (*Two-site XXX_{su(2)} chain (L=2). Sites carry independent gl(2) representations [λ1,0], [λ2,0], in
2. [code] ps[rep_] := ↵  ps[rep] = ↵   Block[{res = (Prepend[rep, #] & /@ (Tuples[↵          Table[Range[rep[[
3. [code] (*spin-s generators*)↵ClearAll[Ee];↵Ee[λ_][i_,j_]:=JJ[i,j,{λ,0}]
4. [code] (*Lax matrix for spin-s representation (no θ baked in — shift applied when a site is instantiated)*)
5. [code] (*companion twist matrix — shared by the whole chain*)↵G={{χ1,-χ2},{1,0}};


## 2026-07-14 13:45:24
**Input:**
```wolfram
ps[rep_] := 
  ps[rep] = 
   Block[{res = (Prepend[rep, #] & /@ (Tuples[
          Table[Range[rep[[1, i + 1]], rep[[1, i]]], {i, 1, 
            Length[rep[[1]]] - 1}]]))}, 
    If[Length[res[[1, 1]]] == 1, Reverse@Sort@res, 
     Reverse@Sort@Flatten[ps /@ res, 1]]];
dim[rep_] := dim[rep] = Length[ps[{rep}]];
J0[k_][pat_] := 
  J0[k][pat] = 
   Sum[pat[[k, i]], {i, 1, k}] - Sum[pat[[k - 1, i]], {i, 1, k - 1}];
Jp[k_][pat_] := 
  Jp[k][pat] = -Sum[
     Product[pat[[k, i]] - i - pat[[k + 1, j]] + j, {j, 
         Range[k + 1]}]/
       Product[pat[[k, i]] - i - pat[[k, j]] + j, {j, 
         Complement[Range[k], {i}]}] v[
       ReplacePart[pat, {k, i} -> pat[[k, i]] + 1]], {i, 1, k}];
Jm[k_][pat_] := 
  Jm[k][pat] = 
   Sum[Product[
       pat[[k, i]] - i - pat[[k - 1, j]] + j, {j, Range[k - 1]}]/
      Product[pat[[k, i]] - i - pat[[k, j]] + j, {j, 
        Complement[Range[k], {i}]}] v[
      ReplacePart[pat, {k, i} -> pat[[k, i]] - 1]], {i, 1, k}];
JJ[k_, l_, rep_] := 
  JJ[k, l, rep] = 
   Which[k == l, DiagonalMatrix[J0[k] /@ ps[{rep}]], l == k + 1, 
    Transpose[
     Table[Coefficient[Jp[k][ps[{rep}][[i]]], v[#]] & /@ 
       ps[{rep}], {i, 1, Length[ps[{rep}]]}]], l == k - 1, 
    Transpose[
     Table[Coefficient[Jm[l][ps[{rep}][[i]]], v[#]] & /@ 
       ps[{rep}], {i, 1, Length[ps[{rep}]]}]], l > k + 1, 
    JJ[k, k + 1, rep] . JJ[k + 1, l, rep] - 
     JJ[k + 1, l, rep] . JJ[k, k + 1, rep], l < k - 1, 
    JJ[k, k - 1, rep] . JJ[k - 1, l, rep] - 
     JJ[k - 1, l, rep] . JJ[k, k - 1, rep]];
```
**Output:**
```
(no output)
```

## 2026-07-14 13:45:24
**Input:**
```wolfram
(*spin-s generators*)
ClearAll[Ee];
Ee[λ_][i_,j_]:=JJ[i,j,{λ,0}]
```
**Output:**
```
(no output)
```

## 2026-07-14 13:45:24
**Input:**
```wolfram
(*Lax matrix for spin-s representation (no θ baked in — shift applied when a site is instantiated)*)
ClearAll[L];
L[λ_][i_,j_][u_]:=KroneckerDelta[i,j]*u*IdentityMatrix[dim[{λ,0}]]-h Ee[λ][j,i];
```
**Output:**
```
(no output)
```

## 2026-07-14 13:45:24
**Input:**
```wolfram
(*companion twist matrix — shared by the whole chain*)
G={{χ1,-χ2},{1,0}};
```
**Output:**
```
(no output)
```
## 2026-07-14 13:45:35 — 📥 BULK INSERT 3 CELL(S) at positions 6–8
1. [markdown] (*Two-site Lax operators on V_λ1⊗V_λ2, built as KroneckerProduct of the single-site Lax matrix (eval
2. [code] ClearAll[L1,L2];↵L1[λ1_,λ2_][i_,j_][u_]:=KroneckerProduct[L[λ1][i,j][u-θ1],IdentityMatrix[dim[{λ2,0}
3. [code] (*dimension check: both site Lax matrices must be (dimλ1*dimλ2)x(dimλ1*dimλ2) for every aux index pa


## 2026-07-14 13:45:36
**Input:**
```wolfram
ClearAll[L1,L2];
L1[λ1_,λ2_][i_,j_][u_]:=KroneckerProduct[L[λ1][i,j][u-θ1],IdentityMatrix[dim[{λ2,0}]]];
L2[λ1_,λ2_][i_,j_][u_]:=KroneckerProduct[IdentityMatrix[dim[{λ1,0}]],L[λ2][i,j][u-θ2]];
```
**Output:**
```
(no output)
```

## 2026-07-14 13:45:36
**Input:**
```wolfram
(*dimension check: both site Lax matrices must be (dimλ1*dimλ2)x(dimλ1*dimλ2) for every aux index pair*)
Table[{λ1,λ2,dim[{λ1,0}]*dim[{λ2,0}],Union@Flatten[Table[Dimensions[L1[λ1,λ2][i,j][u]],{i,2},{j,2}],1],Union@Flatten[Table[Dimensions[L2[λ1,λ2][i,j][u]],{i,2},{j,2}],1]},{λ1,1,3},{λ2,1,3}]
```
**Output:**
```
(no output)
```
## 2026-07-14 13:45:45 — ▶️ RUN CELL 8
**In [8]:** `(*dimension check: both site Lax matrices must be (dimλ1*dimλ2)x(dimλ1*dimλ2) for every aux index pair*)↵Table[{λ1,λ2,dim[{λ1,0}]*dim[{λ2,0}],Union@Flatten[Table[Dimensions[L1[λ1,λ2][i,j][u]],{i,2},{j`
**Out:** (no output)


## 2026-07-14 13:45:49
**Input:**
```wolfram
Table[{λ1,λ2,dim[{λ1,0}]*dim[{λ2,0}],Union@Flatten[Table[Dimensions[L1[λ1,λ2][i,j][u]],{i,2},{j,2}],1],Union@Flatten[Table[Dimensions[L2[λ1,λ2][i,j][u]],{i,2},{j,2}],1]},{λ1,1,3},{λ2,1,3}]
```
**Output:**
```
Out= {{{1, 1, 4, {{4, 4}}, {{4, 4}}}, {1, 2, 6, {{6, 6}}, {{6, 6}}}, {1, 3, 8, {{8, 8}}, {{8, 8}}}}, {{2, 1, 6, {{6, 6}}, {{6, 6}}}, {2, 2, 9, {{9, 9}}, {{9, 9}}}, {2, 3, 12, {{12, 12}}, {{12, 12}}}}, {{3, 1, 8, {{8, 8}}, {{8, 8}}}, {3, 2, 12, {{12, 12}}, {{12, 12}}}, {3, 3, 16, {{16, 16}}, {{16, 16}}}}}
```

## 2026-07-14 13:46:08
**Input:**
```wolfram
(*dimension check: both site Lax matrices must be (dimλ1*dimλ2)x(dimλ1*dimλ2) for every aux index pair*)\nTable[{λ1,λ2,dim[{λ1,0}]*dim[{λ2,0}],Union@Flatten[Table[Dimensions[L1[λ1,λ2][i,j][u]],{i,2},{j,2}],1],Union@Flatten[Table[Dimensions[L2[λ1,λ2][i,j][u]],{i,2},{j,2}],1]},{λ1,1,3},{λ2,1,3}]
```
**Output:**
```
[message] ToExpression::sntx: Invalid syntax in or before "(*dimension check: both site Lax matrices must be (dim\\:03bb1*dim\\:03bb2)x(dim\\:03bb1*dim\\:03bb2) for every aux index pair*)\\nTable[{\\:03bb1,\\:0… [430 chars]
Out[1]= (no output)
```
## 2026-07-14 13:46:22 — ✏️ EDIT CELL 8
(*dimension check: both site Lax matrices must be (dim\[Lambda]1*dim\[Lambda]2)x(dim\[Lambda]1*dim\[Lambda]2) for every aux index pair*)
dimCheck=Table[{\[Lambda]1,\[Lambda]2,dim[{\[Lambda]1,0}]*dim[{…


## 2026-07-14 13:46:22
**Input:**
```wolfram
(*dimension check: both site Lax matrices must be (dim\[Lambda]1*dim\[Lambda]2)x(dim\[Lambda]1*dim\[Lambda]2) for every aux index pair*)
dimCheck=Table[{\[Lambda]1,\[Lambda]2,dim[{\[Lambda]1,0}]*dim[{\[Lambda]2,0}],Union@Flatten[Table[Dimensions[L1[\[Lambda]1,\[Lambda]2][i,j][u]],{i,2},{j,2}],1],Union@Flatten[Table[Dimensions[L2[\[Lambda]1,\[Lambda]2][i,j][u]],{i,2},{j,2}],1]},{\[Lambda]1,1,3},{\[Lambda]2,1,3}]
```
**Output:**
```
Out= {{{1, 1, 4, {{4, 4}}, {{4, 4}}}, {1, 2, 6, {{6, 6}}, {{6, 6}}}, {1, 3, 8, {{8, 8}}, {{8, 8}}}}, {{2, 1, 6, {{6, 6}}, {{6, 6}}}, {2, 2, 9, {{9, 9}}, {{9, 9}}}, {2, 3, 12, {{12, 12}}, {{12, 12}}}}, {{3, 1, 8, {{8, 8}}, {{8, 8}}}, {3, 2, 12, {{12, 12}}, {{12, 12}}}, {3, 3, 16, {{16, 16}}, {{16, 16}}}}}
```
## 2026-07-14 13:46:32 — 📥 BULK INSERT 4 CELL(S) at positions 9–12
1. [markdown] (*Monodromy matrix T[λ1,λ2][i,j][u] = product of the two site Lax operators in auxiliary space, T(u)
2. [code] ClearAll[T];↵T[λ1_,λ2_][i_,j_][u_]:=T[λ1,λ2][i,j][u]=Sum[L1[λ1,λ2][i,k][u].L2[λ1,λ2][k,j][u],{k,2}];
3. [code] ClearAll[t2];↵t2[λ1_,λ2_][1,1][u_]:=Sum[T[λ1,λ2][i,j][u]G[[j,i]],{i,2},{j,2}];
4. [code] (*sanity check: if site 2 carries the trivial rep λ2=0 (dim 1, Ee[0][i,j]=0), the second Lax factor


## 2026-07-14 13:46:33
**Input:**
```wolfram
ClearAll[T];
T[λ1_,λ2_][i_,j_][u_]:=T[λ1,λ2][i,j][u]=Sum[L1[λ1,λ2][i,k][u].L2[λ1,λ2][k,j][u],{k,2}];
```
**Output:**
```
(no output)
```

## 2026-07-14 13:46:33
**Input:**
```wolfram
ClearAll[t2];
t2[λ1_,λ2_][1,1][u_]:=Sum[T[λ1,λ2][i,j][u]G[[j,i]],{i,2},{j,2}];
```
**Output:**
```
(no output)
```

## 2026-07-14 13:46:33
**Input:**
```wolfram
(*sanity check: if site 2 carries the trivial rep λ2=0 (dim 1, Ee[0][i,j]=0), the second Lax factor is just the scalar (u-θ2) times the identity, so it should decouple: t2[λ1,0][1,1][u] == (u-θ2) t[λ1][1,1][u] with θ->θ1. First check Ee[0][i,j] vanishes, then check the reduction.*)
check0=Table[Ee[0][i,j],{i,1},{j,1}];
check0
reduction=Simplify[t2[1,0][1,1][u]-(u-θ2)(t[1][1,1][u]/.θ->θ1)]
```
**Output:**
```
(no output)
```
## 2026-07-14 13:46:43 — ▶️ RUN CELL 12
**In [12]:** `(*sanity check: if site 2 carries the trivial rep λ2=0 (dim 1, Ee[0][i,j]=0), the second Lax factor is just the scalar (u-θ2) times the identity, so it should decouple: t2[λ1,0][1,1][u] == (u-θ2) t[λ1`
**Out:** (no output)

## 2026-07-14 13:46:51 — ✏️ EDIT CELL 12
(*sanity check: if site 2 carries the trivial rep \[Lambda]2=0 (dim 1, Ee[0][i,j]=0), the second Lax factor is just the scalar (u-\[Theta]2) times the identity, so it should decouple: t2[\[Lambda]1,0]…


## 2026-07-14 13:46:51
**Input:**
```wolfram
(*sanity check: if site 2 carries the trivial rep \[Lambda]2=0 (dim 1, Ee[0][i,j]=0), the second Lax factor is just the scalar (u-\[Theta]2) times the identity, so it should decouple: t2[\[Lambda]1,0][1,1][u] == (u-\[Theta]2) t1site[\[Lambda]1][u] where t1site is the single-site L=1 transfer matrix built the same way but with shift \[Theta]1.*)
ClearAll[t1site];
t1site[\[Lambda]1_][u_]:=Sum[L[\[Lambda]1][i,j][u-\[Theta]1]G[[j,i]],{i,2},{j,2}];
check0=Table[Ee[0][i,j],{i,1},{j,1}];
reduction=Simplify[t2[1,0][1,1][u]-(u-\[Theta]2) t1site[1][u]];
{check0,reduction}
```
**Output:**
```
Out= {{{{{0}}}}, {{0, 0}, {0, 0}}}
```
## 2026-07-14 13:47:03 — 📥 BULK INSERT 3 CELL(S) at positions 13–15
1. [markdown] (*Two-site quantum determinant qdetT[λ1,λ2][u] built the same way as the single-site one in XXX_su2.
2. [code] ClearAll[qdetT,qdet1,qdet2];↵qdetT[λ1_,λ2_][u_]:=T[λ1,λ2][1,1][u].T[λ1,λ2][2,2][u-h]-T[λ1,λ2][2,1][u
3. [code] (*test factorization qdetT == qdet1 (x) qdet2(.-h) for a few small (λ1,λ2)*)↵factorTest=Table[{λ1,λ2


## 2026-07-14 13:47:03
**Input:**
```wolfram
ClearAll[qdetT,qdet1,qdet2];
qdetT[λ1_,λ2_][u_]:=T[λ1,λ2][1,1][u].T[λ1,λ2][2,2][u-h]-T[λ1,λ2][2,1][u].T[λ1,λ2][1,2][u-h];
qdet1[λ1_][u_]:=L[λ1][1,1][u-θ1].L[λ1][2,2][u-θ1-h]-L[λ1][2,1][u-θ1].L[λ1][1,2][u-θ1-h];
qdet2[λ2_][u_]:=L[λ2][1,1][u-θ2].L[λ2][2,2][u-θ2-h]-L[λ2][2,1][u-θ2].L[λ2][1,2][u-θ2-h];
```
**Output:**
```
(no output)
```

## 2026-07-14 13:47:04
**Input:**
```wolfram
(*test factorization qdetT == qdet1 (x) qdet2(.-h) for a few small (λ1,λ2)*)
factorTest=Table[{λ1,λ2,Simplify[qdetT[λ1,λ2][u]-KroneckerProduct[qdet1[λ1][u],qdet2[λ2][u-h]]]//Flatten//Union},{λ1,1,2},{λ2,1,2}];
factorTest
```
**Output:**
```
(no output)
```

## 2026-07-14 13:47:12
**Input:**
```wolfram
factorTest
```
**Output:**
```
Out= {{{1, 1, {0, h*(u - θ1)*(2*h - u + θ1)*(3*h - 2*u + 2*θ2)}}, {1, 2, {0, 2*h*(u - θ1)*(2*h - u + θ1)*(2*h - u + θ2)}}}, {{2, 1, {0, h*(u - θ1)*(3*h - u + θ1)*(3*h - 2*u + 2*θ2)}}, {2, 2, {0, 2*h*(u - θ1)*(3*h - u + θ1)*(2*h - u + θ2)}}}}
```

## 2026-07-14 13:47:33
**Input:**
```wolfram
Table[{a,b,Simplify[qdetT[1,1][u]-KroneckerProduct[qdet1[1][u+a*h],qdet2[1][u+b*h]]]//Flatten//Union},{a,-2,2},{b,-2,2}]
```
**Output:**
```
Out= {{{-2, -2, {0, -4*h*(2*h - u + θ1)*(2*h - u + θ2)*(4*h - 2*u + θ1 + θ2)}}, {-2, -1, {0, h*(-2*h + u - θ1)*(12*h^2 + 2*(3*u - θ1 - 2*θ2)*(u - θ2) + h*(-19*u + 3*θ1 + 16*θ2)), -(h*(2*h - u + θ1)*(12*h^2 + 2*(3*u - θ1 - 2*θ2)*(u - θ2) + h*(-19*u + 3*θ1 + 16*θ2)))}}, {-2, 0, {0, 4*h*(2*h - u + θ1)*(u - θ2)*(2*h - u + θ2)}}, {-2, 1, {0, h*(2*h - u + θ1)*(4*h^2 + h*(-u + θ1) - 2*(u + θ1 - 2*θ2)*(u - θ2))}}, {-2, 2, {0, -4*h*(2*h - u + θ1)*(u - θ2)*(2*h + θ1 - θ2)}}}, {{-1, -2, {0, h*(-2*h + u - θ2)*(12*h^2 + 2*(u - θ1)*(3*u - 2*θ1 - θ2) + h*(-19*u + 16*θ1 + 3*θ2)), -(h*(2*h - u + θ2)*(12*h^2 + 2*(u - θ1)*(3*u - 2*θ1 - θ2) + h*(-19*u + 16*θ1 + 3*θ2)))}}, {-1, -1, {0, (u - θ1)*(-2*h + u - θ1)*(u - θ2)*(-2*h + u - θ2) - (h - u + θ1)*(3*h - u + θ1)*(h - u + θ2)*(3*h - u + θ2), h*(-9*h^3 + 12*h^2*(2*u - θ1 - θ2) + 2*(u - θ1)*(u - θ2)*(2*u - θ1 - θ2) - 3*h*(6*u^2 + θ1^2 + 4*θ1*θ2 + θ2^2 - 6*u*(θ1 + θ2)))}}, {-1, 0, {0, h*(3*h - 2*u + 2*θ1)*(u - θ2)*(2*h - u + θ2)}}, {-1, 1, {0, h*(h + θ1 - θ2)*(3*h^2 + 2*(u - θ1)*(u - θ2) + h*(-4*u + θ1 + 3*θ2))}}, {-1, 2, {0, -(h*(u - θ2)*(6*h^2 - 3*h*(3*u - 4*θ1 + θ2) + 2*(u - θ1)*(u - 2*θ1 + θ2)))}}}, {{0, -2, {0, 4*h*(u - θ1)*(2*h - u + θ1)*(2*h - u + θ2)}}, {0, -1, {0, h*(u - θ1)*(2*h - u + θ1)*(3*h - 2*u + 2*θ2)}}, {0, 0, {0}}, {0, 1, {0, h*(u - θ1)*(-2*h + u - θ1)*(h - 2*u + 2*θ2), -(h*(u - θ1)*(2*h - u + θ1)*(h - 2*u + 2*θ2))}}, {0, 2, {0, 4*h*(u - θ1)*(2*h - u + θ1)*(u - θ2)}}}, {{1, -2, {0, h*(2*h - u + θ2)*(4*h^2 + h*(-u + θ2) - 2*(u - θ1)*(u - 2*θ1 + θ2))}}, {1, -1, {0, h*(h - θ1 + θ2)*(3*h^2 + 2*(u - θ1)*(u - θ2) + h*(-4*u + 3*θ1 + θ2))}}, {1, 0, {0, -(h*(h - 2*u + 2*θ1)*(u - θ2)*(2*h - u + θ2))}}, {1, 1, {0, -(h*(h^3 + 2*(u - θ1)*(u - θ2)*(2*u - θ1 - θ2) - h*(6*u^2 + θ1^2 + 4*θ1*θ2 + θ2^2 - 6*u*(θ1 + θ2)))), h*(-h^3 - 2*(u - θ1)*(u - θ2)*(2*u - θ1 - θ2) + h*(6*u^2 + θ1^2 + 4*θ1*θ2 + θ2^2 - 6*u*(θ1 + θ2)))}}, {1, 2, {0, h*(2*h^2 + h*(5*u - 4*θ1 - θ2) - 2*(u - θ1)*(3*u - 2*θ1 - θ2))*(u - θ2), -(h*(u - θ2)*(-2*h^2 + 2*(u - θ1)*(3*u - 2*θ1 - θ2) + h*(-5*u + 4*θ1 + θ2)))}}}, {{2, -2, {0, -4*h*(u - θ1)*(2*h - u + θ2)*(2*h - θ1 + θ2)}}, {2, -1, {0, -(h*(u - θ1)*(6*h^2 - 3*h*(3*u + θ1 - 4*θ2) + 2*(u + θ1 - 2*θ2)*(u - θ2)))}}, {2, 0, {0, 4*h*(u - θ1)*(u - θ2)*(2*h - u + θ2)}}, {2, 1, {0, h*(u - θ1)*(2*h^2 + h*(5*u - θ1 - 4*θ2) - 2*(3*u - θ1 - 2*θ2)*(u - θ2)), -(h*(u - θ1)*(-2*h^2 + 2*(3*u - θ1 - 2*θ2)*(u - θ2) + h*(-5*u + θ1 + 4*θ2)))}}, {2, 2, {0, -4*h*(u - θ1)*(u - θ2)*(2*u - θ1 - θ2)}}}}
```

## 2026-07-14 13:47:38
**Input:**
```wolfram
factorTestNoShift=Table[{λ1,λ2,Simplify[qdetT[λ1,λ2][u]-KroneckerProduct[qdet1[λ1][u],qdet2[λ2][u]]]//Flatten//Union},{λ1,1,3},{λ2,1,3}]
```
**Output:**
```
Out= {{{1, 1, {0}}, {1, 2, {0}}, {1, 3, {0}}}, {{2, 1, {0}}, {2, 2, {0}}, {2, 3, {0}}}, {{3, 1, {0}}, {3, 2, {0}}, {3, 3, {0}}}}
```
## 2026-07-14 13:47:48 — ✏️ EDIT CELL 13 [batch]
(*Two-site quantum determinant qdetT[λ1,λ2][u], built the same way as the single-site one in XXX_su2.wb (using T instead of L, with the second row shifted by h). We also build the single-site quantum …

## 2026-07-14 13:47:48 — ✏️ EDIT CELL 15 [batch]
(*confirm factorization qdetT == qdet1(u) (x) qdet2(u), no shift, for λ1,λ2∈{1,2,3}*)
factorTest=Table[{λ1,λ2,Simplify[qdetT[λ1,λ2][u]-KroneckerProduct[qdet1[λ1][u],qdet2[λ2][u]]]//Flatten//Union},{λ1…


## 2026-07-14 13:47:48
**Input:**
```wolfram
(*confirm factorization qdetT == qdet1(u) (x) qdet2(u), no shift, for λ1,λ2∈{1,2,3}*)
factorTest=Table[{λ1,λ2,Simplify[qdetT[λ1,λ2][u]-KroneckerProduct[qdet1[λ1][u],qdet2[λ2][u]]]//Flatten//Union},{λ1,1,3},{λ2,1,3}];
factorTest
```
**Output:**
```
Out= (no output)
```

## 2026-07-14 13:47:52
**Input:**
```wolfram
factorTest
```
**Output:**
```
Out= {{{1, 1, {0}}, {1, 2, {0}}, {1, 3, {0}}}, {{2, 1, {0}}, {2, 2, {0}}, {2, 3, {0}}}, {{3, 1, {0}}, {3, 2, {0}}, {3, 3, {0}}}}
```
## 2026-07-14 13:48:02 — 📥 BULK INSERT 5 CELL(S) at positions 16–20
1. [markdown] (*Fusion/Hirota hierarchy for the two-site transfer matrices t2[λ1,λ2][1,s][u], built exactly as in 
2. [code] t2[λ1_,λ2_][0,s_][u_]:=IdentityMatrix[dim[{λ1,0}]*dim[{λ2,0}]];↵t2[λ1_,λ2_][a_,0][u_]:=IdentityMatri
3. [code] (*Qθ-type factors for the two sites*)↵ClearAll[Qθ1,Qθ2];↵Qθ1[u_]:=(u-θ1);↵Qθ2[u_]:=(u-θ2);
4. [markdown] (*SoV B-operator for L=2: the bare (1,1) monodromy element (no G, no trace) — same convention as B[λ
5. [code] ClearAll[B2];↵B2[λ1_,λ2_][u_]:=T[λ1,λ2][1,1][u];


## 2026-07-14 13:48:02
**Input:**
```wolfram
t2[λ1_,λ2_][0,s_][u_]:=IdentityMatrix[dim[{λ1,0}]*dim[{λ2,0}]];
t2[λ1_,λ2_][a_,0][u_]:=IdentityMatrix[dim[{λ1,0}]*dim[{λ2,0}]];
t2[λ1_,λ2_][2,1][u_]:=χ2 qdetT[λ1,λ2][u];
t2[λ1_,λ2_][1,s_/;s>1][u_]:=t2[λ1,λ2][1,s-1][u].t2[λ1,λ2][1,1][u+(s-1)h]-t2[λ1,λ2][2,1][u+(s-1)h].t2[λ1,λ2][1,s-2][u];
```
**Output:**
```
(no output)
```

## 2026-07-14 13:48:02
**Input:**
```wolfram
(*Qθ-type factors for the two sites*)
ClearAll[Qθ1,Qθ2];
Qθ1[u_]:=(u-θ1);
Qθ2[u_]:=(u-θ2);
```
**Output:**
```
(no output)
```

## 2026-07-14 13:48:03
**Input:**
```wolfram
ClearAll[B2];
B2[λ1_,λ2_][u_]:=T[λ1,λ2][1,1][u];
```
**Output:**
```
(no output)
```
## 2026-07-14 13:48:23 — 📥 BULK INSERT 3 CELL(S) at positions 21–23
1. [markdown] (*Step 8/9 setup: single-site fused transfer matrices and SoV covectors for an ARBITRARY inhomogenei
2. [code] ClearAll[tSingle,xSingle];↵tSingle[λ_,θ0_][0,s_][u_]:=IdentityMatrix[dim[{λ,0}]];↵tSingle[λ_,θ0_][a_
3. [code] (*cross-check: xSingle at θ0→θ-independent instance must reproduce the same B-eigenvalue property as


## 2026-07-14 13:48:23
**Input:**
```wolfram
ClearAll[tSingle,xSingle];
tSingle[λ_,θ0_][0,s_][u_]:=IdentityMatrix[dim[{λ,0}]];
tSingle[λ_,θ0_][a_,0][u_]:=IdentityMatrix[dim[{λ,0}]];
tSingle[λ_,θ0_][1,1][u_]:=Sum[L[λ][i,j][u-θ0]G[[j,i]],{i,2},{j,2}];
tSingle[λ_,θ0_][2,1][u_]:=χ2(L[λ][1,1][u-θ0].L[λ][2,2][u-θ0-h]-L[λ][2,1][u-θ0].L[λ][1,2][u-θ0-h]//Simplify);
tSingle[λ_,θ0_][1,s_/;s>1][u_]:=tSingle[λ,θ0][1,s-1][u].tSingle[λ,θ0][1,1][u+(s-1)h]-tSingle[λ,θ0][2,1][u+(s-1)h].tSingle[λ,θ0][1,s-2][u];
xSingle[λ_,θ0_][{0}]:={Table[KroneckerDelta[i,dim[{λ,0}]],{i,dim[{λ,0}]}]};
xSingle[λ_,θ0_][{s_}]:=xSingle[λ,θ0][{0}].tSingle[λ,θ0][1,s][θ0]/Product[(θ0-λ h+(p-1)h),{p,1,s}]//Simplify;
```
**Output:**
```
(no output)
```

## 2026-07-14 13:48:23
**Input:**
```wolfram
(*cross-check: xSingle at θ0→θ-independent instance must reproduce the same B-eigenvalue property as the L=1 notebook, using B built the single-site way*)
checkSingle=Table[xSingle[λ,θ1][{s}].L[λ][1,1][u-θ1]-(u-θ1-s h) xSingle[λ,θ1][{s}],{λ,1,4},{s,0,λ}]//Simplify//Flatten//Union;
checkSingle
```
**Output:**
```
(no output)
```

## 2026-07-14 13:48:25
**Input:**
```wolfram
checkSingle
```
**Output:**
```
Out= {0}
```
## 2026-07-14 13:48:34 — 📥 BULK INSERT 5 CELL(S) at positions 24–28
1. [markdown] (*Step 8: SoV pseudovacuum for L=2 as the tensor product of the two single-site lowest-weight covect
2. [code] ClearAll[x2];↵x2[λ1_,λ2_][{0,0}]:=KroneckerProduct[xSingle[λ1,θ1][{0}],xSingle[λ2,θ2][{0}]];
3. [markdown] (*Step 9 — THE KEY CHECK: candidate SoV covector basis x2[λ1,λ2][{s1,s2}] built as the naive tensor 
4. [code] ClearAll[x2];↵x2[λ1_,λ2_][{s1_,s2_}]:=KroneckerProduct[xSingle[λ1,θ1][{s1}],xSingle[λ2,θ2][{s2}]];
5. [code] (*direct eigenvector check for λ1=λ2=1: all 4 candidate states {s1,s2}∈{0,1}x{0,1}*)↵case11=Table[{s


## 2026-07-14 13:48:34
**Input:**
```wolfram
ClearAll[x2];
x2[λ1_,λ2_][{0,0}]:=KroneckerProduct[xSingle[λ1,θ1][{0}],xSingle[λ2,θ2][{0}]];
```
**Output:**
```
(no output)
```

## 2026-07-14 13:48:34
**Input:**
```wolfram
ClearAll[x2];
x2[λ1_,λ2_][{s1_,s2_}]:=KroneckerProduct[xSingle[λ1,θ1][{s1}],xSingle[λ2,θ2][{s2}]];
```
**Output:**
```
(no output)
```

## 2026-07-14 13:48:35
**Input:**
```wolfram
(*direct eigenvector check for λ1=λ2=1: all 4 candidate states {s1,s2}∈{0,1}x{0,1}*)
case11=Table[{s1,s2,Simplify[x2[1,1][{s1,s2}].B2[1,1][u]-(u-θ1-s1 h)(u-θ2-s2 h) x2[1,1][{s1,s2}]]//Flatten//Union},{s1,0,1},{s2,0,1}];
case11
```
**Output:**
```
(no output)
```

## 2026-07-14 13:48:40
**Input:**
```wolfram
case11
```
**Output:**
```
Out= {{{0, 0, {0}}, {0, 1, {0, h^3/(h - θ2)}}}, {{1, 0, {0}}, {1, 1, {0}}}}
```

## 2026-07-14 13:48:47
**Input:**
```wolfram
v01 = x2[1,1][{0,1}];
result = Simplify[v01.B2[1,1][u]];
ratio = Simplify[result / v01];
{v01, result, ratio}
```
**Output:**
```
[message] Power::infy: Infinite expression - encountered.                                  0
[message] [Power::infy ×2 more suppressed — first instance shown above]
[message] General::stop: Further output of Power::infy will be suppressed during this calculation.
[message] Infinity::indet: Indeterminate expression 0 ComplexInfinity encountered.
[message] [Infinity::indet ×1 more suppressed — first instance shown above]
Out= {{{0, 0, h/(h - θ2), 0}}, {{0, h^3/(h - θ2), -((h*(u - θ1)*(h - u + θ2))/(h - θ2)), 0}}, {{Indeterminate, ComplexInfinity, -((u - θ1)*(h - u + θ2)), Indeterminate}}}
```

## 2026-07-14 13:48:56
**Input:**
```wolfram
xSingle[1,θ1][{0}]
xSingle[1,θ1][{1}]
xSingle[1,θ2][{0}]
xSingle[1,θ2][{1}]
x2[1,1][{0,1}]
x2[1,1][{1,0}]
B2[1,1][u]//MatrixForm
```
**Output:**
```
Out[1]= {{0, 1}}
Out[2]= {{h/(h - θ1), 0}}
Out[3]= {{0, 1}}
Out[4]= {{h/(h - θ2), 0}}
Out[5]= {{0, 0, h/(h - θ2), 0}}
Out[6]= {{0, h/(h - θ1), 0, 0}}
Out[7]= MatrixForm[{{(-h + u - θ1)*(-h + u - θ2), 0, 0, 0}, {0, (-h + u - θ1)*(u - θ2), 0, 0}, {0, h^2, (u - θ1)*(-h + u - θ2), 0}, {0, 0, 0, (u - θ1)*(u - θ2)}}]
```

## 2026-07-14 13:50:03
**Input:**
```wolfram
M12 = B2[1,2][u];
Dimensions[M12]
MatrixForm[Simplify[M12]]
```
**Output:**
```
Out[1]= (suppressed)
Out[2]= {6, 6}
Out[3]= MatrixForm[{{(h - u + θ1)*(2*h - u + θ2), 0, 0, 0, 0, 0}, {0, (h - u + θ1)*(h - u + θ2), 0, 0, 0, 0}, {0, 0, -((h - u + θ1)*(u - θ2)), 0, 0, 0}, {0, 2*h^2, 0, (u - θ1)*(-2*h + u - θ2), 0, 0}, {0, 0, 2*h^2, 0, (u - θ1)*(-h + u - θ2), 0}, {0, 0, 0, 0, 0, (u - θ1)*(u - θ2)}}]
```

## 2026-07-14 13:50:46
**Input:**
```wolfram
M22 = B2[2,2][u];
Dimensions[M22]
Simplify[M22]//MatrixForm
```
**Output:**
```
Out[1]= (suppressed)
Out[2]= {9, 9}
Out[3]= MatrixForm[{{(2*h - u + θ1)*(2*h - u + θ2), 0, 0, 0, 0, 0, 0, 0, 0}, {0, (2*h - u + θ1)*(h - u + θ2), 0, 0, 0, 0, 0, 0, 0}, {0, 0, (-2*h + u - θ1)*(u - θ2), 0, 0, 0, 0, 0, 0}, {0, 2*h^2, 0, (h - u + θ1)*(2*h - u + θ2), 0, 0, 0, 0, 0}, {0, 0, 2*h^2, 0, (h - u + θ1)*(h - u + θ2), 0, 0, 0, 0}, {0, 0, 0, 0, 0, -((h - u + θ1)*(u - θ2)), 0, 0, 0}, {0, 0, 0, 0, 2*h^2, 0, (u - θ1)*(-2*h + u - θ2), 0, 0}, {0, 0, 0, 0, 0, 2*h^2, 0, (u - θ1)*(-h + u - θ2), 0}, {0, 0, 0, 0, 0, 0, 0, 0, (u - θ1)*(u - θ2)}}]
```
## 2026-07-14 13:53:40 — 📥 BULK INSERT 5 CELL(S) at positions 29–33
1. [markdown] (*RESULT of the key check: the naive tensor-product ansatz x2[λ1,λ2][{s1,s2}] = xSingle[λ1,θ1][{s1}]
2. [code] (*direct diagonalization confirms this for λ1=λ2=1*)↵MT=Transpose[B2[1,1][u]];↵{evals,evecs}=Simplif
3. [markdown] (*General fix: build the exact SoV covector xCorrected[λ1,λ2][{s1,s2}] by an ansatz X = x2[{s1,s2}] 
4. [code] ClearAll[xCorrected];↵xCorrected[λ1_,λ2_][{s1_,s2_}]:=xCorrected[λ1,λ2][{s1,s2}]=Module[{kmax,a,ansa
5. [code] (*final verification pass: xCorrected is an EXACT B2(u) eigencovector, with eigenvalue (u-θ1-s1h)(u-


## 2026-07-14 13:53:41
**Input:**
```wolfram
(*direct diagonalization confirms this for λ1=λ2=1*)
MT=Transpose[B2[1,1][u]];
{evals,evecs}=Simplify[{Eigenvalues[MT],Eigenvectors[MT]}]
```
**Output:**
```
(no output)
```

## 2026-07-14 13:53:41
**Input:**
```wolfram
ClearAll[xCorrected];
xCorrected[λ1_,λ2_][{s1_,s2_}]:=xCorrected[λ1,λ2][{s1,s2}]=Module[{kmax,a,ansatz,eigenvalue,diff,eqSet,sol},
  If[s2==0||s1==λ1,Return[x2[λ1,λ2][{s1,s2}]]];
  kmax=Min[λ1-s1,s2];
  eigenvalue=(u-θ1-s1 h)(u-θ2-s2 h);
  ansatz=x2[λ1,λ2][{s1,s2}]+Sum[a[k] x2[λ1,λ2][{s1+k,s2-k}],{k,1,kmax}];
  diff=Simplify[ansatz.B2[λ1,λ2][u]-eigenvalue ansatz][[1]];
  eqSet=Flatten[Table[Coefficient[diff[[j]],u,p]==0,{j,Length[diff]},{p,0,2}]];
  sol=Solve[eqSet,Table[a[k],{k,1,kmax}]]//Flatten;
  ansatz/.sol
];
```
**Output:**
```
(no output)
```

## 2026-07-14 13:53:41
**Input:**
```wolfram
(*final verification pass: xCorrected is an EXACT B2(u) eigencovector, with eigenvalue (u-θ1-s1h)(u-θ2-s2h), for all (s1,s2) at (λ1,λ2)∈{(1,1),(1,2),(2,1),(2,2)}*)
finalCheck=Table[
  {λ1,λ2,s1,s2,Simplify[xCorrected[λ1,λ2][{s1,s2}].B2[λ1,λ2][u]-(u-θ1-s1 h)(u-θ2-s2 h) xCorrected[λ1,λ2][{s1,s2}]]//Flatten//Union},
  {λ1,1,2},{λ2,1,2},{s1,0,λ1},{s2,0,λ2}]//Flatten[#,3]&;
Union[finalCheck[[All,-1]]]
```
**Output:**
```
(no output)
```

## 2026-07-14 13:53:45
**Input:**
```wolfram
Union[finalCheck[[All,-1]]]
```
**Output:**
```
Out= {{0}}
```
## 2026-07-14 13:54:00 — ▶️ RUN CELL 33
**In [33]:** `(*final verification pass: xCorrected is an EXACT B2(u) eigencovector, with eigenvalue (u-θ1-s1h)(u-θ2-s2h), for all (s1,s2) at (λ1,λ2)∈{(1,1),(1,2),(2,1),(2,2)}*)↵finalCheck=Table[↵  {λ1,λ2,s1,s2,Sim`
**Out:** (no output)

## 2026-07-14 13:54:09 — ✏️ EDIT CELL 33
(*final verification pass: xCorrected is an EXACT B2(u) eigencovector, with eigenvalue (u-\[Theta]1-s1 h)(u-\[Theta]2-s2 h), for all (s1,s2) at (\[Lambda]1,\[Lambda]2)\[Element]{(1,1),(1,2),(2,1),(2,2…


## 2026-07-14 13:54:09
**Input:**
```wolfram
(*final verification pass: xCorrected is an EXACT B2(u) eigencovector, with eigenvalue (u-\[Theta]1-s1 h)(u-\[Theta]2-s2 h), for all (s1,s2) at (\[Lambda]1,\[Lambda]2)\[Element]{(1,1),(1,2),(2,1),(2,2)}*)
allResiduals=Flatten[Table[
  Simplify[xCorrected[\[Lambda]1,\[Lambda]2][{s1,s2}].B2[\[Lambda]1,\[Lambda]2][u]-(u-\[Theta]1-s1 h)(u-\[Theta]2-s2 h) xCorrected[\[Lambda]1,\[Lambda]2][{s1,s2}]],
  {\[Lambda]1,1,2},{\[Lambda]2,1,2},{s1,0,\[Lambda]1},{s2,0,\[Lambda]2}]];
summaryResult=Union[allResiduals]
```
**Output:**
```
Out= {0}
```

---

## 2026-07-14 13:54:27 — 🔄 KERNEL RESTART

## 2026-07-14 13:54:40 — ▶️ RUN CELLS 1–33
Cell 2: ✓ — (no output)
Cell 3: ✓ — (no output)
Cell 4: ✓ — (no output)
Cell 5: ✓ — (no output)
Cell 7: ✓ — (no output)
Cell 8: ✓ — (no output)
Cell 10: ✓ — (no output)
Cell 11: ✓ — (no output)
Cell 12: ✓ — (no output)
Cell 14: ✓ — (no output)
Cell 15: ✓ — (no output)
Cell 17: ✓ — (no output)
Cell 18: ✓ — (no output)
Cell 20: ✓ — (no output)
Cell 22: ✓ — (no output)
Cell 23: ✓ — (no output)
Cell 25: ✓ — (no output)
Cell 27: ✓ — (no output)
Cell 28: ✓ — (no output)
Cell 30: ✓ — (no output)
Cell 32: ✓ — (no output)
Cell 33: ✓ — (no output)


## 2026-07-14 13:54:50
**Input:**
```wolfram
summaryResult
checkSingle
factorTest
```
**Output:**
```
Out[1]= {0}
Out[2]= {0}
Out[3]= {{{1, 1, {0}}, {1, 2, {0}}, {1, 3, {0}}}, {{2, 1, {0}}, {2, 2, {0}}, {2, 3, {0}}}, {{3, 1, {0}}, {3, 2, {0}}, {3, 3, {0}}}}
```

## 2026-07-14 17:53:58
**Input:**
```wolfram
ClearAll[κ1,κ2];
twistSubstL2 = {χ1 -> κ1+κ2, χ2 -> κ1 κ2};
tm11 = t2[1,1][1,1][u]/.twistSubstL2;
polyDeg = Exponent[tm11[[1,1]],u]
tm11 // MatrixForm
```
**Output:**
```
Out[1]= (suppressed)
Out[2]= (suppressed)
Out[3]= (suppressed)
Out[4]= 2
Out[5]= MatrixForm[{{(-h + u - θ1)*(-h + u - θ2)*(κ1 + κ2), h*(u - θ1)*κ1*κ2, h*(-h + u - θ2)*κ1*κ2, 0}, {-(h*(-h + u - θ1)), (-h + u - θ1)*(u - θ2)*(κ1 + κ2), 0, h*(u - θ2)*κ1*κ2}, {-(h*(u - θ2)), h^2*(κ1 + κ2), (u - θ1)*(-h + u - θ2)*(κ1 + κ2), h*(-h + u - θ1)*κ1*κ2}, {0, -(h*(-h + u - θ2)), -(h*(u - θ1)), (u - θ1)*(u - θ2)*(κ1 + κ2)}}]
```

## 2026-07-14 17:54:04
**Input:**
```wolfram
MTt = Transpose[tm11];
{evalsT, evecsT} = Eigensystem[MTt];
Simplify[evalsT]
```
**Output:**
```
Out[1]= (suppressed)
Out[2]= (suppressed)
Out[3]= {h^2*κ1 + h*(-2*u + θ1 + θ2)*κ1 + (u - θ1)*(u - θ2)*(κ1 + κ2), h*θ1*κ2 + h*(h + θ2)*κ2 + u^2*(κ1 + κ2) + θ1*θ2*(κ1 + κ2) - u*(2*h*κ2 + θ1*(κ1 + κ2) + θ2*(κ1 + κ2)), (u - θ1)*(u - θ2)*(κ1 + κ2) + (h*(-2*u*κ1 + θ1*κ1 + θ2*κ1 - 2*u*κ2 + θ1*κ2 + θ2*κ2 - Sqrt[θ1^2*(κ1 - κ2)^2 - 2*θ1*θ2*(κ1 - κ2)^2 + θ2^2*(κ1 - κ2)^2 + 4*h^2*κ1*κ2]))/2, (u - θ1)*(u - θ2)*(κ1 + κ2) + (h*(-2*u*κ1 + θ1*κ1 + θ2*κ1 - 2*u*κ2 + θ1*κ2 + θ2*κ2 + Sqrt[θ1^2*(κ1 - κ2)^2 - 2*θ1*θ2*(κ1 - κ2)^2 + θ2^2*(κ1 - κ2)^2 + 4*h^2*κ1*κ2]))/2}
```

## 2026-07-14 17:55:03
**Input:**
```wolfram
simpEvecs = Simplify[evecsT];
FreeQ[simpEvecs, u]
simpEvecs
```
**Output:**
```
Out[1]= (suppressed)
Out[2]= False
Out[3]= {{κ2^(-2), -κ2^(-1), -κ2^(-1), 1}, {κ1^(-2), -κ1^(-1), -κ1^(-1), 1}, {1/(κ1*κ2), (2*h^3*κ1*κ2 + h^2*(2*u*κ1^2 - θ1*κ1^2 - θ2*κ1^2 - 2*θ1*κ1*κ2 + 2*θ2*κ1*κ2 + 2*u*κ2^2 - θ1*κ2^2 - θ2*κ2^2 + κ1*Sqrt[θ1^2*(κ1 - κ2)^2 - 2*θ1*θ2*(κ1 - κ2)^2 + θ2^2*(κ1 - κ2)^2 + 4*h^2*κ1*κ2] + κ2*Sqrt[θ1^2*(κ1 - κ2)^2 - 2*θ1*θ2*(κ1 - κ2)^2 + θ2^2*(κ1 - κ2)^2 + 4*h^2*κ1*κ2]) + h*(-2*u^2*κ1^2 + u*θ1*κ1^2 + 3*u*θ2*κ1^2 - θ1*θ2*κ1^2 - θ2^2*κ1^2 + 2*u*θ1*κ1*κ2 - 2*θ1^2*κ1*κ2 - 2*u*θ2*κ1*κ2 + 2*θ1*θ2*κ1*κ2 - 2*u^2*κ2^2 + u*θ1*κ2^2 + 3*u*θ2*κ2^2 - θ1*θ2*κ2^2 - θ2^2*κ2^2 - u*κ1*Sqrt[θ1^2*(κ1 - κ2)^2 - 2*θ1*θ2*(κ1 - κ2)^2 + θ2^2*(κ1 - κ2)^2 + 4*h^2*κ1*κ2] + θ2*κ1*Sqrt[θ1^2*(κ1 - κ2)^2 - 2*θ1*θ2*(κ1 - κ2)^2 + θ2^2*(κ1 - κ2)^2 + 4*h^2*κ1*κ2] - u*κ2*Sqrt[θ1^2*(κ1 - κ2)^2 - 2*θ1*θ2*(κ1 - κ2)^2 + θ2^2*(κ1 - κ2)^2 + 4*h^2*κ1*κ2] + θ2*κ2*Sqrt[θ1^2*(κ1 - κ2)^2 - 2*θ1*θ2*(κ1 - κ2)^2 + θ2^2*(κ1 - κ2)^2 + 4*h^2*κ1*κ2]) + (u - θ1)*(u - θ2)*(θ1*(κ1 - κ2)^2 - θ2*(κ1 - κ2)^2 + (κ1 + κ2)*Sqrt[θ1^2*(κ1 - κ2)^2 - 2*θ1*θ2*(κ1 - κ2)^2 + θ2^2*(κ1 - κ2)^2 + 4*h^2*κ1*κ2]))/(κ1*κ2*(-2*(-u + θ1)*(-u + θ2)*Sqrt[θ1^2*(κ1 - κ2)^2 - 2*θ1*θ2*(κ1 - κ2)^2 + θ2^2*(κ1 - κ2)^2 + 4*h^2*κ1*κ2] + h^2*(-2*u*κ1 + θ1*κ1 + θ2*κ1 - 2*u*κ2 + θ1*κ2 + θ2*κ2 - Sqrt[θ1^2*(κ1 - κ2)^2 - 2*θ1*θ2*(κ1 - κ2)^2 + θ2^2*(κ1 - κ2)^2 + 4*h^2*κ1*κ2]) + h*(2*u^2*κ1 - 2*u*θ1*κ1 + θ1^2*κ1 - 2*u*θ2*κ1 + θ2^2*κ1 + 2*u^2*κ2 - 2*u*θ1*κ2 + θ1^2*κ2 - 2*u*θ2*κ2 + θ2^2*κ2 + 2*u*Sqrt[θ1^2*(κ1 - κ2)^2 - 2*θ1*θ2*(κ1 - κ2)^2 + θ2^2*(κ1 - κ2)^2 + 4*h^2*κ1*κ2] - θ1*Sqrt[θ1^2*(κ1 - κ2)^2 - 2*θ1*θ2*(κ1 - κ2)^2 + θ2^2*(κ1 - κ2)^2 + 4*h^2*κ1*κ2] - θ2*Sqrt[θ1^2*(κ1 - κ2)^2 - 2*θ1*θ2*(κ1 - κ2)^2 + θ2^2*(κ1 - κ2)^2 + 4*h^2*κ1*κ2]))), (-2*h^3*κ1*κ2 + 4*h^2*(u - θ2)*κ1*κ2 + h*(u*θ1*κ1^2 - θ1^2*κ1^2 - u*θ2*κ1^2 + θ1*θ2*κ1^2 - 4*u^2*κ1*κ2 + 2*u*θ1*κ1*κ2 + 6*u*θ2*κ1*κ2 - 2*θ1*θ2*κ1*κ2 - 2*θ2^2*κ1*κ2 + u*θ1*κ2^2 - θ1^2*κ2^2 - u*θ2*κ2^2 + θ1*θ2*κ2^2 - u*κ1*Sqrt[θ1^2*(κ1 - κ2)^2 - 2*θ1*θ2*(κ1 - κ2)^2 + θ2^2*(κ1 - κ2)^2 + 4*h^2*κ1*κ2] + θ1*κ1*Sqrt[θ1^2*(κ1 - κ2)^2 - 2*θ1*θ2*(κ1 - κ2)^2 + θ2^2*(κ1 - κ2)^2 + 4*h^2*κ1*κ2] - u*κ2*Sqrt[θ1^2*(κ1 - κ2)^2 - 2*θ1*θ2*(κ1 - κ2)^2 + θ2^2*(κ1 - κ2)^2 + 4*h^2*κ1*κ2] + θ1*κ2*Sqrt[θ1^2*(κ1 - κ2)^2 - 2*θ1*θ2*(κ1 - κ2)^2 + θ2^2*(κ1 - κ2)^2 + 4*h^2*κ1*κ2]) + (-u + θ1)*(-u + θ2)*(-(θ1*(κ1 - κ2)^2) + θ2*(κ1 - κ2)^2 + (κ1 + κ2)*Sqrt[θ1^2*(κ1 - κ2)^2 - 2*θ1*θ2*(κ1 - κ2)^2 + θ2^2*(κ1 - κ2)^2 + 4*h^2*κ1*κ2]))/(κ1*κ2*(-2*(-u + θ1)*(-u + θ2)*Sqrt[θ1^2*(κ1 - κ2)^2 - 2*θ1*θ2*(κ1 - κ2)^2 + θ2^2*(κ1 - κ2)^2 + 4*h^2*κ1*κ2] + h^2*(-2*u*κ1 + θ1*κ1 + θ2*κ1 - 2*u*κ2 + θ1*κ2 + θ2*κ2 - Sqrt[θ1^2*(κ1 - κ2)^2 - 2*θ1*θ2*(κ1 - κ2)^2 + θ2^2*(κ1 - κ2)^2 + 4*h^2*κ1*κ2]) + h*(2*u^2*κ1 - 2*u*θ1*κ1 + θ1^2*κ1 - 2*u*θ2*κ1 + θ2^2*κ1 + 2*u^2*κ2 - 2*u*θ1*κ2 + θ1^2*κ2 - 2*u*θ2*κ2 + θ2^2*κ2 + 2*u*Sqrt[θ1^2*(κ1 - κ2)^2 - 2*θ1*θ2*(κ1 - κ2)^2 + θ2^2*(κ1 - κ2)^2 + 4*h^2*κ1*κ2] - θ1*Sqrt[θ1^2*(κ1 - κ2)^2 - 2*θ1*θ2*(κ1 - κ2)^2 + θ2^2*(κ1 - κ2)^2 + 4*h^2*κ1*κ2] - θ2*Sqrt[θ1^2*(κ1 - κ2)^2 - 2*θ1*θ2*(κ1 - κ2)^2 + θ2^2*(κ1 - κ2)^2 + 4*h^2*κ1*κ2]))), 1}, {1/(κ1*κ2), -((-2*h^3*κ1*κ2 + h^2*(-2*u*κ1^2 + θ1*κ1^2 + θ2*κ1^2 + 2*θ1*κ1*κ2 - 2*θ2*κ1*κ2 - 2*u*κ2^2 + θ1*κ2^2 + θ2*κ2^2 + κ1*Sqrt[θ1^2*(κ1 - κ2)^2 - 2*θ1*θ2*(κ1 - κ2)^2 + θ2^2*(κ1 - κ2)^2 + 4*h^2*κ1*κ2] + κ2*Sqrt[θ1^2*(κ1 - κ2)^2 - 2*θ1*θ2*(κ1 - κ2)^2 + θ2^2*(κ1 - κ2)^2 + 4*h^2*κ1*κ2]) + h*(2*u^2*κ1^2 - u*θ1*κ1^2 - 3*u*θ2*κ1^2 + θ1*θ2*κ1^2 + θ2^2*κ1^2 - 2*u*θ1*κ1*κ2 + 2*θ1^2*κ1*κ2 + 2*u*θ2*κ1*κ2 - 2*θ1*θ2*κ1*κ2 + 2*u^2*κ2^2 - u*θ1*κ2^2 - 3*u*θ2*κ2^2 + θ1*θ2*κ2^2 + θ2^2*κ2^2 - u*κ1*Sqrt[θ1^2*(κ1 - κ2)^2 - 2*θ1*θ2*(κ1 - κ2)^2 + θ2^2*(κ1 - κ2)^2 + 4*h^2*κ1*κ2] + θ2*κ1*Sqrt[θ1^2*(κ1 - κ2)^2 - 2*θ1*θ2*(κ1 - κ2)^2 + θ2^2*(κ1 - κ2)^2 + 4*h^2*κ1*κ2] - u*κ2*Sqrt[θ1^2*(κ1 - κ2)^2 - 2*θ1*θ2*(κ1 - κ2)^2 + θ2^2*(κ1 - κ2)^2 + 4*h^2*κ1*κ2] + θ2*κ2*Sqrt[θ1^2*(κ1 - κ2)^2 - 2*θ1*θ2*(κ1 - κ2)^2 + θ2^2*(κ1 - κ2)^2 + 4*h^2*κ1*κ2]) + (-u + θ1)*(-u + θ2)*(-(θ1*(κ1 - κ2)^2) + θ2*(κ1 - κ2)^2 + (κ1 + κ2)*Sqrt[θ1^2*(κ1 - κ2)^2 - 2*θ1*θ2*(κ1 - κ2)^2 + θ2^2*(κ1 - κ2)^2 + 4*h^2*κ1*κ2]))/(κ1*κ2*(2*(-u + θ1)*(-u + θ2)*Sqrt[θ1^2*(κ1 - κ2)^2 - 2*θ1*θ2*(κ1 - κ2)^2 + θ2^2*(κ1 - κ2)^2 + 4*h^2*κ1*κ2] + h^2*(-2*u*κ1 + θ1*κ1 + θ2*κ1 - 2*u*κ2 + θ
[output truncated — 5874 chars total]
```

## 2026-07-14 17:55:15
**Input:**
```wolfram
ClearAll[xNew];
xNew[λ1_,λ2_][{s1_,s2_}] := x2[λ1,λ2][{0,0}] . t2[λ1,λ2][1,s1][θ1] . t2[λ1,λ2][1,s2][θ2];
(* quick test λ1=λ2=1, all 4 states *)
test11 = Table[{s1,s2,Simplify[xNew[1,1][{s1,s2}].B2[1,1][u]-(u-θ1-s1 h)(u-θ2-s2 h) xNew[1,1][{s1,s2}]]//Flatten//Union},{s1,0,1},{s2,0,1}];
test11
```
**Output:**
```
Out[1]= (suppressed)
Out[2]= (suppressed)
Out[3]= (suppressed)
Out[4]= {{{0, 0, {0}}, {0, 1, {0}}}, {{1, 0, {0}}, {1, 1, {0}}}}
```

## 2026-07-14 17:55:19
**Input:**
```wolfram
testAll = Flatten[Table[
  {λ1,λ2,s1,s2,Simplify[xNew[λ1,λ2][{s1,s2}].B2[λ1,λ2][u]-(u-θ1-s1 h)(u-θ2-s2 h) xNew[λ1,λ2][{s1,s2}]]//Flatten//Union},
  {λ1,1,2},{λ2,1,2},{s1,0,λ1},{s2,0,λ2}],3];
Union[testAll[[All,-1]]]
```
**Output:**
```
Out= {{0}}
```

## 2026-07-14 17:55:23
**Input:**
```wolfram
xNewSwapped[λ1_,λ2_][{s1_,s2_}] := x2[λ1,λ2][{0,0}] . t2[λ1,λ2][1,s2][θ2] . t2[λ1,λ2][1,s1][θ1];
orderTest = Table[Simplify[xNew[λ1,λ2][{s1,s2}]-xNewSwapped[λ1,λ2][{s1,s2}]]//Flatten//Union,{λ1,1,2},{λ2,1,2},{s1,0,λ1},{s2,0,λ2}]//Flatten//Union
```
**Output:**
```
Out= {0}
```
## 2026-07-14 17:56:04 — ✏️ EDIT CELL 29 [batch]
(*RESULT of the naive-tensor-product test: diagonalizing B2[1,1][u] directly (Eigensystem of Transpose[B2[u]], symbolic in u) confirms the eigenvalue guess (u-θ1-s1 h)(u-θ2-s2 h) is exactly right for …

## 2026-07-14 17:56:04 — ✏️ EDIT CELL 31 [batch]
(*Correct SoV covector construction (per direct verification): ⟨x_{0,0}| = x[λ1]{0}⊗x[λ2]{0} (already x2[λ1,λ2][{0,0}]), and

⟨x_{s1,s2}| := ⟨x_{0,0}| . t2[λ1,λ2][1,s1][θ1] . t2[λ1,λ2][1,s2][θ2]

i.e.…

## 2026-07-14 17:56:04 — ✏️ EDIT CELL 32 [batch]
ClearAll[xCorrected];
xCorrected[λ1_,λ2_][{s1_,s2_}] := xCorrected[λ1,λ2][{s1,s2}] = x2[λ1,λ2][{0,0}] . t2[λ1,λ2][1,s1][θ1] . t2[λ1,λ2][1,s2][θ2];


## 2026-07-14 17:56:04
**Input:**
```wolfram
ClearAll[xCorrected];
xCorrected[λ1_,λ2_][{s1_,s2_}] := xCorrected[λ1,λ2][{s1,s2}] = x2[λ1,λ2][{0,0}] . t2[λ1,λ2][1,s1][θ1] . t2[λ1,λ2][1,s2][θ2];
```
**Output:**
```
Out= (no output)
```
## 2026-07-14 17:56:04 — ✏️ EDIT CELL 33 [batch]
(*final verification pass: xCorrected is an EXACT B2(u) eigencovector, with eigenvalue (u-\[Theta]1-s1 h)(u-\[Theta]2-s2 h), for all (s1,s2) at (\[Lambda]1,\[Lambda]2)\[Element]{(1,1),(1,2),(2,1),(2,2…


## 2026-07-14 17:56:05
**Input:**
```wolfram
(*final verification pass: xCorrected is an EXACT B2(u) eigencovector, with eigenvalue (u-\[Theta]1-s1 h)(u-\[Theta]2-s2 h), for all (s1,s2) at (\[Lambda]1,\[Lambda]2)\[Element]{(1,1),(1,2),(2,1),(2,2)} — no ad hoc correction term needed.*)
allResiduals=Flatten[Table[
  Simplify[xCorrected[\[Lambda]1,\[Lambda]2][{s1,s2}].B2[\[Lambda]1,\[Lambda]2][u]-(u-\[Theta]1-s1 h)(u-\[Theta]2-s2 h) xCorrected[\[Lambda]1,\[Lambda]2][{s1,s2}]],
  {\[Lambda]1,1,2},{\[Lambda]2,1,2},{s1,0,\[Lambda]1},{s2,0,\[Lambda]2}]];
summaryResult=Union[allResiduals]
```
**Output:**
```
Out= (no output)
```

## 2026-07-14 17:56:09
**Input:**
```wolfram
summaryResult
```
**Output:**
```
Out= {0}
```

## 2026-07-14 17:56:15
**Input:**
```wolfram
(* re-derive eigenvalues/eigenvectors of t2[1,1][u] with the twist diagonalized *)
ClearAll[κ1,κ2];
twistSubstL2 = {χ1 -> κ1+κ2, χ2 -> κ1 κ2};
tm11 = t2[1,1][1,1][u]/.twistSubstL2;
MTt = Transpose[tm11];
{evalsT,evecsT} = Eigensystem[MTt];
evalsSimplified = Simplify[evalsT]
```
**Output:**
```
Out[1]= (suppressed)
Out[2]= (suppressed)
Out[3]= (suppressed)
Out[4]= (suppressed)
Out[5]= (suppressed)
Out[6]= {h^2*κ1 + h*(-2*u + θ1 + θ2)*κ1 + (u - θ1)*(u - θ2)*(κ1 + κ2), h*θ1*κ2 + h*(h + θ2)*κ2 + u^2*(κ1 + κ2) + θ1*θ2*(κ1 + κ2) - u*(2*h*κ2 + θ1*(κ1 + κ2) + θ2*(κ1 + κ2)), (u - θ1)*(u - θ2)*(κ1 + κ2) + (h*(-2*u*κ1 + θ1*κ1 + θ2*κ1 - 2*u*κ2 + θ1*κ2 + θ2*κ2 - Sqrt[θ1^2*(κ1 - κ2)^2 - 2*θ1*θ2*(κ1 - κ2)^2 + θ2^2*(κ1 - κ2)^2 + 4*h^2*κ1*κ2]))/2, (u - θ1)*(u - θ2)*(κ1 + κ2) + (h*(-2*u*κ1 + θ1*κ1 + θ2*κ1 - 2*u*κ2 + θ1*κ2 + θ2*κ2 + Sqrt[θ1^2*(κ1 - κ2)^2 - 2*θ1*θ2*(κ1 - κ2)^2 + θ2^2*(κ1 - κ2)^2 + 4*h^2*κ1*κ2]))/2}
```

## 2026-07-14 17:56:38
**Input:**
```wolfram
(* check whether the sqrt-branch eigenvectors are genuinely u-independent by evaluating at two different u *)
numSubst = {κ1->3,κ2->7/2,h->11/10,θ1->13/10,θ2->17/10};
vecAtU0 = evecsT[[3]]/.numSubst/.u->2;
vecAtU1 = evecsT[[3]]/.numSubst/.u->5;
{N[vecAtU0], N[vecAtU1], N[vecAtU0]/N[vecAtU0[[1]]], N[vecAtU1]/N[vecAtU1[[1]]]}
```
**Output:**
```
Out= {{0.09523809523809523, -1.2810625667938889, 0.6620149477462707, 1.}, {0.09523809523809523, -1.2810625667938904, 0.6620149477462723, 1.}, {1., -13.451156951335832, 6.951156951335842, 10.5}, {1., -13.45115695133585, 6.95115695133586, 10.5}}
```

## 2026-07-14 17:57:11
**Input:**
```wolfram
(* candidate formula generalizing L=1's tEigenvalue to L=2, per-site independent M's *)
ClearAll[tEigenvalue2];
tEigenvalue2[λ1_,λ2_,M1_,M2_][u_] := κ1 (u-θ1-(λ1-M1)h)(u-θ2-(λ2-M2)h) + κ2 (u-θ1-M1 h)(u-θ2-M2 h);
numSubst2 = {κ1->3,κ2->7/2,h->11/10,θ1->13/10,θ2->17/10,u->17/10};
actualEvals = Sort[N[evalsT/.twistSubstL2/.numSubst2]];
candidateEvals = Sort[N[Flatten[Table[tEigenvalue2[1,1,M1,M2][u],{M1,0,1},{M2,0,1}]]/.numSubst2]];
{actualEvals, candidateEvals, Chop[actualEvals-candidateEvals]}
```
**Output:**
```
Out= {{-5.352390852528595, 2.31, 2.4923908525285947, 2.695}, {-1.54, -1.32, 2.31, 2.695}, {-3.812390852528595, 3.63, 0.18239085252859466, 0}}
```

## 2026-07-14 17:58:44
**Input:**
```wolfram
(* redo properly: RIGHT eigenvectors (matching L=1 convention: t[λ][1,1][u].Ψ[λ,M]==eigenvalue Ψ[λ,M]) *)
tm11 = t2[1,1][1,1][u]/.twistSubstL2;
esR = Eigensystem[tm11];
valsR = esR[[1]]; vecsR = esR[[2]];
(* check u-independence numerically at two points *)
numSubst2 = {κ1->3,κ2->7/2,h->11/10,θ1->13/10,θ2->17/10};
v0 = vecsR/.numSubst2/.u->2;
v1 = vecsR/.numSubst2/.u->5;
Table[N[v0[[k]]]/N[v0[[k,1]]] - N[v1[[k]]]/N[v1[[k,1]]],{k,4}]//Chop
```
**Output:**
```
Out[1]= (suppressed)
Out[2]= (suppressed)
Out[3]= (suppressed)
Out[4]= (suppressed)
Out[5]= (suppressed)
Out[6]= (suppressed)
Out[7]= {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}}
```

## 2026-07-14 17:58:59
**Input:**
```wolfram
(* test small-h continuity matching for λ1=λ2=1 *)
numSmallH = {κ1->3,κ2->7/2,h->1/1000,θ1->13/10,θ2->17/10,u->17/10};
actualEvalsSmallH = N[valsR/.twistSubstL2/.numSmallH]//Chop;
Sort[actualEvalsSmallH]
candList = Table[{M1,M2}->N[tEigenvalue2[1,1,M1,M2][u]/.numSmallH],{M1,0,1},{M2,0,1}]//Flatten
Sort[candList[[All,2]]]
```
**Output:**
```
Out[1]= (suppressed)
Out[2]= (suppressed)
Out[3]= {-0.0014000524862259804, -0.0013965, -0.0011999475137740195, -0.001197}
Out[4]= {{0, 0} -> -0.001197, {0, 1} -> -0.0014, {1, 0} -> -0.0012, {1, 1} -> -0.0013965}
Out[5]= {-0.0014, -0.0013965, -0.0012, -0.001197}
```

## 2026-07-14 17:59:25
**Input:**
```wolfram
ClearAll[buildPsi2Raw];
buildPsi2Raw[λ1_,λ2_] := Module[{tm,es,vals,vecs,pairs,numSmall,actualNum,candNum,matchIdx,assoc},
  tm = t2[λ1,λ2][1,1][u]/.twistSubstL2;
  es = Eigensystem[tm];
  vals = es[[1]]; vecs = es[[2]];
  pairs = Flatten[Table[{M1,M2},{M1,0,λ1},{M2,0,λ2}],1];
  numSmall = {κ1->3,κ2->7/2,h->1/1000,θ1->13/10,θ2->17/10,u->17/10};
  actualNum = N[vals/.numSmall];
  candNum = N[(tEigenvalue2[λ1,λ2,#[[1]],#[[2]]][u]&/@pairs)/.numSmall];
  matchIdx = Table[First@Ordering[Abs[candNum[[p]]-actualNum],1],{p,Length[pairs]}];
  If[!(Sort[matchIdx]===Range[Length[pairs]]), Print["WARNING: matching not a bijection for λ1=",λ1,",λ2=",λ2]];
  assoc = AssociationThread[pairs -> Transpose[{vals[[matchIdx]],vecs[[matchIdx]]}]];
  assoc
];
data11 = buildPsi2Raw[1,1];
Keys[data11]
data11[{0,1}][[1]] (* eigenvalue for (0,1) *)
```
**Output:**
```
Out[1]= (suppressed)
Out[2]= (suppressed)
Out[3]= (suppressed)
Out[4]= {{0, 0}, {0, 1}, {1, 0}, {1, 1}}
Out[5]= (-2*h*u*κ1 + 2*u^2*κ1 + h*θ1*κ1 - 2*u*θ1*κ1 + h*θ2*κ1 - 2*u*θ2*κ1 + 2*θ1*θ2*κ1 - 2*h*u*κ2 + 2*u^2*κ2 + h*θ1*κ2 - 2*u*θ1*κ2 + h*θ2*κ2 - 2*u*θ2*κ2 + 2*θ1*θ2*κ2 - h*Sqrt[θ1^2*κ1^2 - 2*θ1*θ2*κ1^2 + θ2^2*κ1^2 + 4*h^2*κ1*κ2 - 2*θ1^2*κ1*κ2 + 4*θ1*θ2*κ1*κ2 - 2*θ2^2*κ1*κ2 + θ1^2*κ2^2 - 2*θ1*θ2*κ2^2 + θ2^2*κ2^2])/2
```

## 2026-07-14 17:59:51
**Input:**
```wolfram
ClearAll[buildPsi2Raw];
buildPsi2Raw[λ1_,λ2_] := buildPsi2Raw[λ1,λ2] = Module[{tm0,es0,vals0,vecs0,pairs,numSmall,candNum0,matchIdx,tmU,eigenvalueOfU,assoc},
  tm0 = t2[λ1,λ2][1,1][u]/.twistSubstL2/.u->0;
  es0 = Eigensystem[tm0];
  vals0 = es0[[1]]; vecs0 = es0[[2]];
  tmU = t2[λ1,λ2][1,1][u]/.twistSubstL2;
  pairs = Flatten[Table[{M1,M2},{M1,0,λ1},{M2,0,λ2}],1];
  numSmall = {κ1->3,κ2->7/2,h->1/1000,θ1->13/10,θ2->17/10};
  candNum0 = N[(tEigenvalue2[λ1,λ2,#[[1]],#[[2]]][0]&/@pairs)/.numSmall];
  matchIdx = Table[First@Ordering[Abs[candNum0[[p]]-N[vals0/.numSmall]],1],{p,Length[pairs]}];
  If[!(Sort[matchIdx]===Range[Length[pairs]]), Print["WARNING: non-bijective match for ",λ1,λ2]];
  assoc = AssociationThread[pairs -> Table[
     Module[{v=vecs0[[matchIdx[[p]]]],nz,eigU,prod},
       prod = tmU.v; (* vector as function of u *)
       nz = First@FirstPosition[v, x_/;x=!=0];
       eigU = Simplify[prod[[nz]]/v[[nz]]];
       {eigU, v}
     ],
  {p,Length[pairs]}]];
  assoc
];
data11 = buildPsi2Raw[1,1];
Table[{pk, data11[pk][[1]]}, {pk, Keys[data11]}]
```
**Output:**
```
Out[1]= (suppressed)
Out[2]= (suppressed)
Out[3]= (suppressed)
Out[4]= {{{0, 0}, 1}, {{0, 1}, 1}, {{1, 0}, 1}, {{1, 1}, 1}}
```

## 2026-07-14 17:59:57
**Input:**
```wolfram
tm0test = t2[1,1][1,1][u]/.twistSubstL2/.u->0;
es0test = Eigensystem[tm0test];
vals0test = es0test[[1]];
vecs0test = es0test[[2]];
vals0test
vecs0test
```
**Output:**
```
Out[1]= (suppressed)
Out[2]= (suppressed)
Out[3]= (suppressed)
Out[4]= (suppressed)
Out[5]= {h^2*κ1 + h*θ1*κ1 + h*θ2*κ1 + θ1*θ2*κ1 + θ1*θ2*κ2, θ1*θ2*κ1 + h^2*κ2 + h*θ1*κ2 + h*θ2*κ2 + θ1*θ2*κ2, (h*θ1*κ1 + h*θ2*κ1 + 2*θ1*θ2*κ1 + h*θ1*κ2 + h*θ2*κ2 + 2*θ1*θ2*κ2 - h*Sqrt[θ1^2*κ1^2 - 2*θ1*θ2*κ1^2 + θ2^2*κ1^2 + 4*h^2*κ1*κ2 - 2*θ1^2*κ1*κ2 + 4*θ1*θ2*κ1*κ2 - 2*θ2^2*κ1*κ2 + θ1^2*κ2^2 - 2*θ1*θ2*κ2^2 + θ2^2*κ2^2])/2, (h*θ1*κ1 + h*θ2*κ1 + 2*θ1*θ2*κ1 + h*θ1*κ2 + h*θ2*κ2 + 2*θ1*θ2*κ2 + h*Sqrt[θ1^2*κ1^2 - 2*θ1*θ2*κ1^2 + θ2^2*κ1^2 + 4*h^2*κ1*κ2 - 2*θ1^2*κ1*κ2 + 4*θ1*θ2*κ1*κ2 - 2*θ2^2*κ1*κ2 + θ1^2*κ2^2 - 2*θ1*θ2*κ2^2 + θ2^2*κ2^2])/2}
Out[6]= {{κ1^2, κ1, κ1, 1}, {κ2^2, κ2, κ2, 1}, {κ1*κ2, -((h*θ1*θ2*κ1^2 + θ1^2*θ2*κ1^2 - h*θ2^2*κ1^2 - θ1*θ2^2*κ1^2 - 2*h^3*κ1*κ2 - 4*h^2*θ1*κ1*κ2 - 2*h*θ1^2*κ1*κ2 - 2*h*θ1*θ2*κ1*κ2 - 2*θ1^2*θ2*κ1*κ2 + 2*θ1*θ2^2*κ1*κ2 + h*θ1*θ2*κ2^2 + θ1^2*θ2*κ2^2 - h*θ2^2*κ2^2 - θ1*θ2^2*κ2^2 + h*θ2*κ1*Sqrt[θ1^2*κ1^2 - 2*θ1*θ2*κ1^2 + θ2^2*κ1^2 + 4*h^2*κ1*κ2 - 2*θ1^2*κ1*κ2 + 4*θ1*θ2*κ1*κ2 - 2*θ2^2*κ1*κ2 + θ1^2*κ2^2 - 2*θ1*θ2*κ2^2 + θ2^2*κ2^2] + θ1*θ2*κ1*Sqrt[θ1^2*κ1^2 - 2*θ1*θ2*κ1^2 + θ2^2*κ1^2 + 4*h^2*κ1*κ2 - 2*θ1^2*κ1*κ2 + 4*θ1*θ2*κ1*κ2 - 2*θ2^2*κ1*κ2 + θ1^2*κ2^2 - 2*θ1*θ2*κ2^2 + θ2^2*κ2^2] + h*θ2*κ2*Sqrt[θ1^2*κ1^2 - 2*θ1*θ2*κ1^2 + θ2^2*κ1^2 + 4*h^2*κ1*κ2 - 2*θ1^2*κ1*κ2 + 4*θ1*θ2*κ1*κ2 - 2*θ2^2*κ1*κ2 + θ1^2*κ2^2 - 2*θ1*θ2*κ2^2 + θ2^2*κ2^2] + θ1*θ2*κ2*Sqrt[θ1^2*κ1^2 - 2*θ1*θ2*κ1^2 + θ2^2*κ1^2 + 4*h^2*κ1*κ2 - 2*θ1^2*κ1*κ2 + 4*θ1*θ2*κ1*κ2 - 2*θ2^2*κ1*κ2 + θ1^2*κ2^2 - 2*θ1*θ2*κ2^2 + θ2^2*κ2^2])/(h^2*θ1*κ1 + h*θ1^2*κ1 + h^2*θ2*κ1 + h*θ2^2*κ1 + h^2*θ1*κ2 + h*θ1^2*κ2 + h^2*θ2*κ2 + h*θ2^2*κ2 - h^2*Sqrt[θ1^2*κ1^2 - 2*θ1*θ2*κ1^2 + θ2^2*κ1^2 + 4*h^2*κ1*κ2 - 2*θ1^2*κ1*κ2 + 4*θ1*θ2*κ1*κ2 - 2*θ2^2*κ1*κ2 + θ1^2*κ2^2 - 2*θ1*θ2*κ2^2 + θ2^2*κ2^2] - h*θ1*Sqrt[θ1^2*κ1^2 - 2*θ1*θ2*κ1^2 + θ2^2*κ1^2 + 4*h^2*κ1*κ2 - 2*θ1^2*κ1*κ2 + 4*θ1*θ2*κ1*κ2 - 2*θ2^2*κ1*κ2 + θ1^2*κ2^2 - 2*θ1*θ2*κ2^2 + θ2^2*κ2^2] - h*θ2*Sqrt[θ1^2*κ1^2 - 2*θ1*θ2*κ1^2 + θ2^2*κ1^2 + 4*h^2*κ1*κ2 - 2*θ1^2*κ1*κ2 + 4*θ1*θ2*κ1*κ2 - 2*θ2^2*κ1*κ2 + θ1^2*κ2^2 - 2*θ1*θ2*κ2^2 + θ2^2*κ2^2] - 2*θ1*θ2*Sqrt[θ1^2*κ1^2 - 2*θ1*θ2*κ1^2 + θ2^2*κ1^2 + 4*h^2*κ1*κ2 - 2*θ1^2*κ1*κ2 + 4*θ1*θ2*κ1*κ2 - 2*θ2^2*κ1*κ2 + θ1^2*κ2^2 - 2*θ1*θ2*κ2^2 + θ2^2*κ2^2])), -((-4*h^2*θ1*κ1^2 - 4*h*θ1^2*κ1^2 - 4*h^2*θ2*κ1^2 - 4*h*θ1*θ2*κ1^2 - 2*θ1^2*θ2*κ1^2 + 2*θ2^3*κ1^2 + 8*h^3*κ1*κ2 + 8*h^2*θ1*κ1*κ2 + 8*h*θ1*θ2*κ1*κ2 + 4*θ1^2*θ2*κ1*κ2 - 8*h*θ2^2*κ1*κ2 - 4*θ2^3*κ1*κ2 - 4*h^2*θ1*κ2^2 - 4*h*θ1^2*κ2^2 - 4*h^2*θ2*κ2^2 - 4*h*θ1*θ2*κ2^2 - 2*θ1^2*θ2*κ2^2 + 2*θ2^3*κ2^2 + 4*h^2*κ1*Sqrt[θ1^2*κ1^2 - 2*θ1*θ2*κ1^2 + θ2^2*κ1^2 + 4*h^2*κ1*κ2 - 2*θ1^2*κ1*κ2 + 4*θ1*θ2*κ1*κ2 - 2*θ2^2*κ1*κ2 + θ1^2*κ2^2 - 2*θ1*θ2*κ2^2 + θ2^2*κ2^2] + 4*h*θ1*κ1*Sqrt[θ1^2*κ1^2 - 2*θ1*θ2*κ1^2 + θ2^2*κ1^2 + 4*h^2*κ1*κ2 - 2*θ1^2*κ1*κ2 + 4*θ1*θ2*κ1*κ2 - 2*θ2^2*κ1*κ2 + θ1^2*κ2^2 - 2*θ1*θ2*κ2^2 + θ2^2*κ2^2] + 4*θ1*θ2*κ1*Sqrt[θ1^2*κ1^2 - 2*θ1*θ2*κ1^2 + θ2^2*κ1^2 + 4*h^2*κ1*κ2 - 2*θ1^2*κ1*κ2 + 4*θ1*θ2*κ1*κ2 - 2*θ2^2*κ1*κ2 + θ1^2*κ2^2 - 2*θ1*θ2*κ2^2 + θ2^2*κ2^2] + 4*h^2*κ2*Sqrt[θ1^2*κ1^2 - 2*θ1*θ2*κ1^2 + θ2^2*κ1^2 + 4*h^2*κ1*κ2 - 2*θ1^2*κ1*κ2 + 4*θ1*θ2*κ1*κ2 - 2*θ2^2*κ1*κ2 + θ1^2*κ2^2 - 2*θ1*θ2*κ2^2 + θ2^2*κ2^2] + 4*h*θ1*κ2*Sqrt[θ1^2*κ1^2 - 2*θ1*θ2*κ1^2 + θ2^2*κ1^2 + 4*h^2*κ1*κ2 - 2*θ1^2*κ1*κ2 + 4*θ1*θ2*κ1*κ2 - 2*θ2^2*κ1*κ2 + θ1^2*κ2^2 - 2*θ1*θ2*κ2^2 + θ2^2*κ2^2] + 4*θ1*θ2*κ2*Sqrt[θ1^2*κ1^2 - 2*θ1*θ2*κ1^2 + θ2^2*κ1^2 + 4*h^2*κ1*κ2 - 2*θ1^2*κ1*κ2 + 4*θ1*θ2*κ1*κ2 - 2*θ2^2*κ1*κ2 + θ1^2*κ2^2 - 2*θ1*θ2*κ2^2 + θ2^2*κ2^2] - 2*θ2*(θ1^2*κ1^2 - 2*θ1*θ2*κ1^2 + θ2^2*κ1^2 + 4*h^2*κ1*κ2 - 2*θ1^2*κ1*κ2 + 4*θ1*θ2*κ1*κ2 - 2*θ2^2*κ1*κ2 + θ1^2*κ2^2 - 2*θ1*θ2*κ2^2 + θ2^2*κ2^2))/(4*h^2*θ1*κ1 + 4*h*θ1^2*κ1 + 4*h^2*θ2*κ1 + 4*h*θ2^2*κ1 + 4*h^2*θ1*κ2 + 4*h*θ1^2*κ2 + 4*h^2*θ2*κ2 + 4*h*θ2^2*κ2 - 4*h^2*Sqrt[θ1^2*κ1^2 - 2*θ1*θ2*κ1^2 + θ2^2*κ1^2 + 4*h^2*κ1*κ2 - 2*θ1^2*κ1*κ2 + 4*θ1*θ2*κ1*κ2 - 2*θ2^2*κ1*κ2 + θ1^2*κ2^2 - 2*θ1*θ2*κ2^2 + θ2^2*κ2^2] - 4*h*θ1*Sqrt[θ1^2*κ1^2 - 2*θ1*θ2*κ1^2 + θ2^2*κ1^2 + 4*h^2*κ1*κ2 - 2*θ1^2*κ1*κ2 + 4*θ1*θ2*κ1*κ2 - 2*θ2^2*κ1*κ2 + θ1^2*κ2^2 - 2*θ1*θ2*κ2^2 + θ2^2*κ2^2] - 4*h*θ2*Sqrt[θ1^2*κ1^2 - 2*θ1*θ2*κ1^2 + θ2^2*κ1^2 + 4*h^2*κ1*κ2 - 2*θ1^2*κ1*κ2 + 4*θ1*θ2*κ1*κ2 - 2*θ2^2*κ1*κ2 + θ1^2*κ2^2 - 2*θ1*θ2*κ2^2 + θ2^2*κ2^2] - 8*θ1*θ2*Sqrt[θ1^2*κ1^2 - 2*θ1*θ2*κ1^2 + θ2^2*κ1^2 + 4*h^2*κ1*κ2 - 2*θ1^2*κ1*κ2 + 4*θ1*θ2*κ1*κ2 - 2*θ2^2*κ1*κ2 + θ1^2*κ2^2 - 2*θ1*θ2*κ2^2 + θ2^2*κ2^2])), 1}, {κ1*κ2, -((h*θ1*θ2*κ1^2 + θ1^2*θ2*κ1^2 - h*θ2^2*κ1^2 - θ1*θ2^2*κ1^2 - 2*h^3*κ1*κ2 - 4*h^2*θ1*κ1*κ2 - 2*h*θ1^2*κ1*κ2 - 2*h*θ1*θ2*κ1*κ2 - 2*θ1^2*θ2*κ1*κ2 + 2*θ1*θ2^2*κ1*κ2 + h*θ1*θ2*κ2^2 + θ1^2*θ2*κ2^2 - h*θ2^2*κ2^2 - θ1*θ2^2*κ2^2 - h*θ2*κ1*Sqrt[θ1^2*κ1^2 - 2*θ1*θ2*κ1^2 + θ2^2*κ1^2 + 4*h^2*κ1*κ2 - 2*θ1^2*κ1*κ2 + 4*θ1*θ2*κ1*κ2 - 2*θ2^2*κ1*κ2 + θ1^2*κ2^2 - 2*θ1*θ2*κ2^2 + θ2^2*κ2^2] - θ1*θ2*κ1*Sqrt[θ1^2*κ1^2 - 2*θ1*θ2*κ1^2 + θ2^2*κ1^2 + 
[output truncated — 7280 chars total]
```

## 2026-07-14 18:00:12
**Input:**
```wolfram
v00 = vecs0test[[1]]; (* the {κ1^2,κ1,κ1,1} eigenvector, corresponds to eigenvalue h^2κ1+hθ1κ1+hθ2κ1+θ1θ2κ1+θ1θ2κ2 *)
tmUtest = t2[1,1][1,1][u]/.twistSubstL2;
prodTest = tmUtest.v00;
prodTest
v00
```
**Output:**
```
Out[1]= Null
Out[2]= (suppressed)
Out[3]= (suppressed)
Out[4]= {h*(u - θ1)*κ1^2*κ2 + h*(-h + u - θ2)*κ1^2*κ2 + (-h + u - θ1)*(-h + u - θ2)*κ1^2*(κ1 + κ2), -(h*(-h + u - θ1)*κ1^2) + h*(u - θ2)*κ1*κ2 + (-h + u - θ1)*(u - θ2)*κ1*(κ1 + κ2), -(h*(u - θ2)*κ1^2) + h*(-h + u - θ1)*κ1*κ2 + h^2*κ1*(κ1 + κ2) + (u - θ1)*(-h + u - θ2)*κ1*(κ1 + κ2), -(h*(u - θ1)*κ1) - h*(-h + u - θ2)*κ1 + (u - θ1)*(u - θ2)*(κ1 + κ2)}
Out[5]= {κ1^2, κ1, κ1, 1}
```

## 2026-07-14 18:00:39
**Input:**
```wolfram
v = vecs0test[[1]];
nzTest = First@FirstPosition[v, x_/;x=!=0];
nzTest
prodT = tmUtest.v;
prodT[[nzTest]]
v[[nzTest]]
Simplify[prodT[[nzTest]]/v[[nzTest]]]
```
**Output:**
```
Out[1]= (suppressed)
Out[2]= (suppressed)
Out[3]= 0
Out[4]= (suppressed)
Out[5]= List
Out[6]= List
Out[7]= 1
```

## 2026-07-14 18:00:47
**Input:**
```wolfram
FirstPosition[{κ1^2,κ1,κ1,1}, x_/;x=!=0]
FirstPosition[{0,κ1,κ1,1}, x_/;x=!=0]
```
**Output:**
```
Out[1]= {0}
Out[2]= {0}
```

## 2026-07-14 18:01:00
**Input:**
```wolfram
ClearAll[buildPsi2Raw];
buildPsi2Raw[λ1_,λ2_] := buildPsi2Raw[λ1,λ2] = Module[{tm0,es0,vals0,vecs0,pairs,numSmall,candNum0,matchIdx,tmU,assoc},
  tm0 = t2[λ1,λ2][1,1][u]/.twistSubstL2/.u->0;
  es0 = Eigensystem[tm0];
  vals0 = es0[[1]]; vecs0 = es0[[2]];
  tmU = t2[λ1,λ2][1,1][u]/.twistSubstL2;
  pairs = Flatten[Table[{M1,M2},{M1,0,λ1},{M2,0,λ2}],1];
  numSmall = {κ1->3,κ2->7/2,h->1/1000,θ1->13/10,θ2->17/10};
  candNum0 = N[(tEigenvalue2[λ1,λ2,#[[1]],#[[2]]][0]&/@pairs)/.numSmall];
  matchIdx = Table[First@Ordering[Abs[candNum0[[p]]-N[vals0/.numSmall]],1],{p,Length[pairs]}];
  If[!(Sort[matchIdx]===Range[Length[pairs]]), Print["WARNING: non-bijective match for ",λ1,λ2]];
  assoc = AssociationThread[pairs -> Table[
     Module[{v=vecs0[[matchIdx[[p]]]],nz,eigU,prod},
       prod = tmU.v;
       nz = First@FirstPosition[v, x_/;x=!=0, {}, {1}];
       eigU = Simplify[prod[[nz]]/v[[nz]]];
       {eigU, v}
     ],
  {p,Length[pairs]}]];
  assoc
];
data11 = buildPsi2Raw[1,1];
Table[{pk, data11[pk][[1]]}, {pk, Keys[data11]}]
```
**Output:**
```
Out[1]= (suppressed)
Out[2]= (suppressed)
Out[3]= (suppressed)
Out[4]= {{{0, 0}, 1}, {{0, 1}, 1}, {{1, 0}, 1}, {{1, 1}, 1}}
```

## 2026-07-14 18:01:06
**Input:**
```wolfram
FirstPosition[{κ1^2,κ1,κ1,1}, x_/;x=!=0, {}, {1}]
ClearAll[buildPsi2Raw];
?buildPsi2Raw
```
**Output:**
```
Out[1]= {0}
Out[2]= (suppressed)
Out[3]= InformationData[<|"ObjectType" -> "Symbol", "Usage" -> "Global`buildPsi2Raw", "Documentation" -> None, "OwnValues" -> None, "UpValues" -> None, "DownValues" -> None, "SubValues" -> None, "DefaultValues" -> None, "NValues" -> None, "FormatValues" -> None, "Options" -> None, "Attributes" -> {}, "FullName" -> "Global`buildPsi2Raw"|>, False]
```

## 2026-07-14 18:01:11
**Input:**
```wolfram
FirstPosition[{0,5,3}, x_/;x=!=0, {}, {1}]
FirstPosition[{0,5,3}, x_/;x=!=0]
Position[{0,5,3},x_/;x=!=0,{1}]
```
**Output:**
```
Out[1]= {0}
Out[2]= {0}
Out[3]= {{0}, {2}, {3}}
```

## 2026-07-14 18:01:23
**Input:**
```wolfram
SelectFirst[Range[Length[v]], v[[#]]=!=0 &]
```
**Output:**
```
Out= 1
```

## 2026-07-14 18:01:33
**Input:**
```wolfram
ClearAll[buildPsi2Raw];
buildPsi2Raw[λ1_,λ2_] := buildPsi2Raw[λ1,λ2] = Module[{tm0,es0,vals0,vecs0,pairs,numSmall,candNum0,matchIdx,tmU,assoc},
  tm0 = t2[λ1,λ2][1,1][u]/.twistSubstL2/.u->0;
  es0 = Eigensystem[tm0];
  vals0 = es0[[1]]; vecs0 = es0[[2]];
  tmU = t2[λ1,λ2][1,1][u]/.twistSubstL2;
  pairs = Flatten[Table[{M1,M2},{M1,0,λ1},{M2,0,λ2}],1];
  numSmall = {κ1->3,κ2->7/2,h->1/1000,θ1->13/10,θ2->17/10};
  candNum0 = N[(tEigenvalue2[λ1,λ2,#[[1]],#[[2]]][0]&/@pairs)/.numSmall];
  matchIdx = Table[First@Ordering[Abs[candNum0[[p]]-N[vals0/.numSmall]],1],{p,Length[pairs]}];
  If[!(Sort[matchIdx]===Range[Length[pairs]]), Print["WARNING: non-bijective match for ",λ1,λ2]];
  assoc = AssociationThread[pairs -> Table[
     Module[{v=vecs0[[matchIdx[[p]]]],nz,eigU,prod},
       prod = tmU.v;
       nz = SelectFirst[Range[Length[v]], v[[#]]=!=0 &];
       eigU = Simplify[prod[[nz]]/v[[nz]]];
       {eigU, v}
     ],
  {p,Length[pairs]}]];
  assoc
];
data11 = buildPsi2Raw[1,1];
Table[{pk, data11[pk][[1]]}, {pk, Keys[data11]}]
```
**Output:**
```
Out[1]= (suppressed)
Out[2]= (suppressed)
Out[3]= (suppressed)
Out[4]= {{{0, 0}, h^2*κ1 + h*(-2*u + θ1 + θ2)*κ1 + (u - θ1)*(u - θ2)*(κ1 + κ2)}, {{0, 1}, ((h - u + θ1)*(h - u + θ2)*(κ1 + κ2)*(h^2*θ1*κ1 + h*θ1^2*κ1 + h^2*θ2*κ1 + h*θ2^2*κ1 + h^2*θ1*κ2 + h*θ1^2*κ2 + h^2*θ2*κ2 + h*θ2^2*κ2 - h^2*Sqrt[θ1^2*(κ1 - κ2)^2 - 2*θ1*θ2*(κ1 - κ2)^2 + θ2^2*(κ1 - κ2)^2 + 4*h^2*κ1*κ2] - h*θ1*Sqrt[θ1^2*(κ1 - κ2)^2 - 2*θ1*θ2*(κ1 - κ2)^2 + θ2^2*(κ1 - κ2)^2 + 4*h^2*κ1*κ2] - h*θ2*Sqrt[θ1^2*(κ1 - κ2)^2 - 2*θ1*θ2*(κ1 - κ2)^2 + θ2^2*(κ1 - κ2)^2 + 4*h^2*κ1*κ2] - 2*θ1*θ2*Sqrt[θ1^2*(κ1 - κ2)^2 - 2*θ1*θ2*(κ1 - κ2)^2 + θ2^2*(κ1 - κ2)^2 + 4*h^2*κ1*κ2]) - h*(u - θ1)*(h*θ1*θ2*κ1^2 + θ1^2*θ2*κ1^2 - h*θ2^2*κ1^2 - θ1*θ2^2*κ1^2 - 2*h^3*κ1*κ2 - 4*h^2*θ1*κ1*κ2 - 2*h*θ1^2*κ1*κ2 - 2*h*θ1*θ2*κ1*κ2 - 2*θ1^2*θ2*κ1*κ2 + 2*θ1*θ2^2*κ1*κ2 + h*θ1*θ2*κ2^2 + θ1^2*θ2*κ2^2 - h*θ2^2*κ2^2 - θ1*θ2^2*κ2^2 + h*θ2*κ1*Sqrt[θ1^2*(κ1 - κ2)^2 - 2*θ1*θ2*(κ1 - κ2)^2 + θ2^2*(κ1 - κ2)^2 + 4*h^2*κ1*κ2] + θ1*θ2*κ1*Sqrt[θ1^2*(κ1 - κ2)^2 - 2*θ1*θ2*(κ1 - κ2)^2 + θ2^2*(κ1 - κ2)^2 + 4*h^2*κ1*κ2] + h*θ2*κ2*Sqrt[θ1^2*(κ1 - κ2)^2 - 2*θ1*θ2*(κ1 - κ2)^2 + θ2^2*(κ1 - κ2)^2 + 4*h^2*κ1*κ2] + θ1*θ2*κ2*Sqrt[θ1^2*(κ1 - κ2)^2 - 2*θ1*θ2*(κ1 - κ2)^2 + θ2^2*(κ1 - κ2)^2 + 4*h^2*κ1*κ2]) + h*(h - u + θ2)*(2*h^3*κ1*κ2 + θ1*θ2*(-(θ1*(κ1 - κ2)^2) + θ2*(κ1 - κ2)^2 + (κ1 + κ2)*Sqrt[θ1^2*(κ1 - κ2)^2 - 2*θ1*θ2*(κ1 - κ2)^2 + θ2^2*(κ1 - κ2)^2 + 4*h^2*κ1*κ2]) + h*(-(θ1*θ2*(κ1 - κ2)^2) - 2*θ2^2*κ1*κ2 + θ1*(κ1 + κ2)*Sqrt[θ1^2*(κ1 - κ2)^2 - 2*θ1*θ2*(κ1 - κ2)^2 + θ2^2*(κ1 - κ2)^2 + 4*h^2*κ1*κ2] - θ1^2*(κ1^2 + κ2^2)) + h^2*(-(θ1*(κ1 - κ2)^2) + (κ1 + κ2)*(-(θ2*κ1) - θ2*κ2 + Sqrt[θ1^2*(κ1 - κ2)^2 - 2*θ1*θ2*(κ1 - κ2)^2 + θ2^2*(κ1 - κ2)^2 + 4*h^2*κ1*κ2]))))/(h^2*θ1*κ1 + h*θ1^2*κ1 + h^2*θ2*κ1 + h*θ2^2*κ1 + h^2*θ1*κ2 + h*θ1^2*κ2 + h^2*θ2*κ2 + h*θ2^2*κ2 - h^2*Sqrt[θ1^2*(κ1 - κ2)^2 - 2*θ1*θ2*(κ1 - κ2)^2 + θ2^2*(κ1 - κ2)^2 + 4*h^2*κ1*κ2] - h*θ1*Sqrt[θ1^2*(κ1 - κ2)^2 - 2*θ1*θ2*(κ1 - κ2)^2 + θ2^2*(κ1 - κ2)^2 + 4*h^2*κ1*κ2] - h*θ2*Sqrt[θ1^2*(κ1 - κ2)^2 - 2*θ1*θ2*(κ1 - κ2)^2 + θ2^2*(κ1 - κ2)^2 + 4*h^2*κ1*κ2] - 2*θ1*θ2*Sqrt[θ1^2*(κ1 - κ2)^2 - 2*θ1*θ2*(κ1 - κ2)^2 + θ2^2*(κ1 - κ2)^2 + 4*h^2*κ1*κ2])}, {{1, 0}, (h - u + θ1)*(h - u + θ2)*(κ1 + κ2) + (h*(u - θ1)*(-(h*θ1*θ2*κ1^2) - θ1^2*θ2*κ1^2 + h*θ2^2*κ1^2 + θ1*θ2^2*κ1^2 + 2*h^3*κ1*κ2 + 4*h^2*θ1*κ1*κ2 + 2*h*θ1^2*κ1*κ2 + 2*h*θ1*θ2*κ1*κ2 + 2*θ1^2*θ2*κ1*κ2 - 2*θ1*θ2^2*κ1*κ2 - h*θ1*θ2*κ2^2 - θ1^2*θ2*κ2^2 + h*θ2^2*κ2^2 + θ1*θ2^2*κ2^2 + h*θ2*κ1*Sqrt[θ1^2*(κ1 - κ2)^2 - 2*θ1*θ2*(κ1 - κ2)^2 + θ2^2*(κ1 - κ2)^2 + 4*h^2*κ1*κ2] + θ1*θ2*κ1*Sqrt[θ1^2*(κ1 - κ2)^2 - 2*θ1*θ2*(κ1 - κ2)^2 + θ2^2*(κ1 - κ2)^2 + 4*h^2*κ1*κ2] + h*θ2*κ2*Sqrt[θ1^2*(κ1 - κ2)^2 - 2*θ1*θ2*(κ1 - κ2)^2 + θ2^2*(κ1 - κ2)^2 + 4*h^2*κ1*κ2] + θ1*θ2*κ2*Sqrt[θ1^2*(κ1 - κ2)^2 - 2*θ1*θ2*(κ1 - κ2)^2 + θ2^2*(κ1 - κ2)^2 + 4*h^2*κ1*κ2]))/(h^2*θ1*κ1 + h*θ1^2*κ1 + h^2*θ2*κ1 + h*θ2^2*κ1 + h^2*θ1*κ2 + h*θ1^2*κ2 + h^2*θ2*κ2 + h*θ2^2*κ2 + h^2*Sqrt[θ1^2*(κ1 - κ2)^2 - 2*θ1*θ2*(κ1 - κ2)^2 + θ2^2*(κ1 - κ2)^2 + 4*h^2*κ1*κ2] + h*θ1*Sqrt[θ1^2*(κ1 - κ2)^2 - 2*θ1*θ2*(κ1 - κ2)^2 + θ2^2*(κ1 - κ2)^2 + 4*h^2*κ1*κ2] + h*θ2*Sqrt[θ1^2*(κ1 - κ2)^2 - 2*θ1*θ2*(κ1 - κ2)^2 + θ2^2*(κ1 - κ2)^2 + 4*h^2*κ1*κ2] + 2*θ1*θ2*Sqrt[θ1^2*(κ1 - κ2)^2 - 2*θ1*θ2*(κ1 - κ2)^2 + θ2^2*(κ1 - κ2)^2 + 4*h^2*κ1*κ2]) + (h*(h - u + θ2)*(2*h^3*κ1*κ2 - θ1*θ2*(θ1*(κ1 - κ2)^2 - θ2*(κ1 - κ2)^2 + (κ1 + κ2)*Sqrt[θ1^2*(κ1 - κ2)^2 - 2*θ1*θ2*(κ1 - κ2)^2 + θ2^2*(κ1 - κ2)^2 + 4*h^2*κ1*κ2]) - h*(θ1*θ2*(κ1 - κ2)^2 + 2*θ2^2*κ1*κ2 + θ1*(κ1 + κ2)*Sqrt[θ1^2*(κ1 - κ2)^2 - 2*θ1*θ2*(κ1 - κ2)^2 + θ2^2*(κ1 - κ2)^2 + 4*h^2*κ1*κ2] + θ1^2*(κ1^2 + κ2^2)) - h^2*(θ1*(κ1 - κ2)^2 + (κ1 + κ2)*(θ2*κ1 + θ2*κ2 + Sqrt[θ1^2*(κ1 - κ2)^2 - 2*θ1*θ2*(κ1 - κ2)^2 + θ2^2*(κ1 - κ2)^2 + 4*h^2*κ1*κ2]))))/(h^2*θ1*κ1 + h*θ1^2*κ1 + h^2*θ2*κ1 + h*θ2^2*κ1 + h^2*θ1*κ2 + h*θ1^2*κ2 + h^2*θ2*κ2 + h*θ2^2*κ2 + h^2*Sqrt[θ1^2*(κ1 - κ2)^2 - 2*θ1*θ2*(κ1 - κ2)^2 + θ2^2*(κ1 - κ2)^2 + 4*h^2*κ1*κ2] + h*θ1*Sqrt[θ1^2*(κ1 - κ2)^2 - 2*θ1*θ2*(κ1 - κ2)^2 + θ2^2*(κ1 - κ2)^2 + 4*h^2*κ1*κ2] + h*θ2*Sqrt[θ1^2*(κ1 - κ2)^2 - 2*θ1*θ2*(κ1 - κ2)^2 + θ2^2*(κ1 - κ2)^2 + 4*h^2*κ1*κ2] + 2*θ1*θ2*Sqrt[θ1^2*(κ1 - κ2)^2 - 2*θ1*θ2*(κ1 - κ2)^2 + θ2^2*(κ1 - κ2)^2 + 4*h^2*κ1*κ2])}, {{1, 1}, h*(u - θ1)*κ1 - h*(h - u + θ2)*κ1 + (h - u + θ1)*(h - u + θ2)*(κ1 + κ2)}}
```

## 2026-07-14 18:01:56
**Input:**
```wolfram
ClearAll[Ψ2];
Ψ2[λ1_,λ2_][{M1_,M2_}] := Module[{data,eigU,raw,overlap},
  data = buildPsi2Raw[λ1,λ2][{M1,M2}];
  eigU = data[[1]]; raw = data[[2]];
  overlap = Simplify[First[xCorrected[λ1,λ2][{0,0}].raw]];
  If[overlap===0, $Failed, raw/overlap]
];
(* sanity check λ1=λ2=1: t2[1,1][u].Ψ2 == eigenvalue Ψ2, identically in u *)
tmU11 = t2[1,1][1,1][u]/.twistSubstL2;
check11 = Table[
  {M1,M2, Simplify[tmU11.Ψ2[1,1][{M1,M2}] - buildPsi2Raw[1,1][{M1,M2}][[1]] Ψ2[1,1][{M1,M2}]]//Flatten//Union},
{M1,0,1},{M2,0,1}];
check11
```
**Output:**
```
Out[1]= (suppressed)
Out[2]= (suppressed)
Out[3]= (suppressed)
Out[4]= (suppressed)
Out[5]= {{{0, 0, {0}}, {0, 1, {0}}}, {{1, 0, {0}}, {1, 1, {0}}}}
```

## 2026-07-14 18:03:16
**Input:**
```wolfram
(* normalization check *)
normCheck11 = Table[Simplify[First[xCorrected[1,1][{0,0}].Ψ2[1,1][{M1,M2}]]-1],{M1,0,1},{M2,0,1}]//Flatten//Union
(* now λ1=1,λ2=2 *)
data12 = buildPsi2Raw[1,2];
Keys[data12]
```
**Output:**
```
Out[1]= {0}
Out[2]= (suppressed)
Out[3]= {{0, 0}, {0, 1}, {0, 2}, {1, 0}, {1, 1}, {1, 2}}
```

## 2026-07-14 18:04:32
**Input:**
```wolfram
tmU12 = t2[1,2][1,1][u]/.twistSubstL2;
check12 = Table[
  {M1,M2, Simplify[tmU12.Ψ2[1,2][{M1,M2}] - buildPsi2Raw[1,2][{M1,M2}][[1]] Ψ2[1,2][{M1,M2}]]//Flatten//Union},
{M1,0,1},{M2,0,2}];
check12
normCheck12 = Table[Simplify[First[xCorrected[1,2][{0,0}].Ψ2[1,2][{M1,M2}]]-1],{M1,0,1},{M2,0,2}]//Flatten//Union
```
**Output:**
```
Out[1]= (suppressed)
Out[2]= (suppressed)
Out[3]= {{{0, 0, {0}}, {0, 1, {0}}, {0, 2, {0}}}, {{1, 0, {0}}, {1, 1, {0}}, {1, 2, {0}}}}
Out[4]= {0}
```

## 2026-07-14 18:04:51
**Input:**
```wolfram
(* which states are "clean" (match naive per-site product formula exactly) vs mixed, for (1,2) *)
cleanCheck12 = Table[
  {M1,M2, Simplify[buildPsi2Raw[1,2][{M1,M2}][[1]] - tEigenvalue2[1,2,M1,M2][u]]===0},
{M1,0,1},{M2,0,2}]
```
**Output:**
```
Out= {{{0, 0, True}, {0, 1, False}, {0, 2, False}}, {{1, 0, False}, {1, 1, False}, {1, 2, True}}}
```
## 2026-07-14 18:05:21 — 📥 BULK INSERT 5 CELL(S) at positions 34–38
1. [markdown] (*Ψ2[λ1,λ2][{M1,M2}]: L=2 analog of Ψ[λ,M] from XXX_su2.wb — the genuine eigenvector of the PHYSICAL
2. [code] ClearAll[κ1,κ2];↵twistSubstL2 = {χ1 -> κ1+κ2, χ2 -> κ1 κ2};↵(* naive per-site product candidate, gen
3. [markdown] (*STEP 2: the physical eigenvalues of t2[1,1][u] are NOT simply tEigenvalue2[λ1,λ2,M1,M2][u] for eve
4. [code] (* diagonalize t2[1,1][u] once at u=0 (eigenvectors are u-independent, matching L=1's construction p
5. [code] (*confirm which (M1,M2) states are "clean" (exact naive product formula) vs genuinely mixed, for (λ1


## 2026-07-14 18:05:22
**Input:**
```wolfram
ClearAll[κ1,κ2];
twistSubstL2 = {χ1 -> κ1+κ2, χ2 -> κ1 κ2};
(* naive per-site product candidate, generalizing L=1's tEigenvalue[λ,M][u] = κ1(u-θ-(λ-M)h)+κ2(u-θ-Mh) *)
ClearAll[tEigenvalue2];
tEigenvalue2[λ1_,λ2_,M1_,M2_][u_] := κ1 (u-θ1-(λ1-M1)h)(u-θ2-(λ2-M2)h) + κ2 (u-θ1-M1 h)(u-θ2-M2 h);
```
**Output:**
```
(no output)
```

## 2026-07-14 18:05:22
**Input:**
```wolfram
(* diagonalize t2[1,1][u] once at u=0 (eigenvectors are u-independent, matching L=1's construction point), label
branches (M1,M2) via small-h nearest-candidate matching, then recover the exact eigenvalue as a function of u
by taking the ratio of a nonzero component of t2[1,1][u].v to the same component of v (u-independent v). *)
ClearAll[buildPsi2Raw];
buildPsi2Raw[λ1_,λ2_] := buildPsi2Raw[λ1,λ2] = Module[{tm0,es0,vals0,vecs0,pairs,numSmall,candNum0,matchIdx,tmU,assoc},
  tm0 = t2[λ1,λ2][1,1][u]/.twistSubstL2/.u->0;
  es0 = Eigensystem[tm0];
  vals0 = es0[[1]]; vecs0 = es0[[2]];
  tmU = t2[λ1,λ2][1,1][u]/.twistSubstL2;
  pairs = Flatten[Table[{M1,M2},{M1,0,λ1},{M2,0,λ2}],1];
  numSmall = {κ1->3,κ2->7/2,h->1/1000,θ1->13/10,θ2->17/10};
  candNum0 = N[(tEigenvalue2[λ1,λ2,#[[1]],#[[2]]][0]&/@pairs)/.numSmall];
  matchIdx = Table[First@Ordering[Abs[candNum0[[p]]-N[vals0/.numSmall]],1],{p,Length[pairs]}];
  If[!(Sort[matchIdx]===Range[Length[pairs]]), Print["WARNING: non-bijective match for λ1=",λ1,", λ2=",λ2]];
  assoc = AssociationThread[pairs -> Table[
     Module[{v=vecs0[[matchIdx[[p]]]],nz,eigU,prod},
       prod = tmU.v;
       nz = SelectFirst[Range[Length[v]], v[[#]]=!=0 &];
       eigU = Simplify[prod[[nz]]/v[[nz]]];
       {eigU, v}
     ],
  {p,Length[pairs]}]];
  assoc
];
```
**Output:**
```
(no output)
```

## 2026-07-14 18:05:51
**Input:**
```wolfram
(*confirm which (M1,M2) states are "clean" (exact naive product formula) vs genuinely mixed, for (λ1,λ2)=(1,1) and (1,2)*)
cleanCheck11 = Table[{M1,M2,Simplify[buildPsi2Raw[1,1][{M1,M2}][[1]]-tEigenvalue2[1,1,M1,M2][u]]===0},{M1,0,1},{M2,0,1}];
cleanCheck12 = Table[{M1,M2,Simplify[buildPsi2Raw[1,2][{M1,M2}][[1]]-tEigenvalue2[1,2,M1,M2][u]]===0},{M1,0,1},{M2,0,2}];
{cleanCheck11, cleanCheck12}
```
**Output:**
```
(no output)
```
## 2026-07-14 18:06:42 — ▶️ RUN CELL 38
**In [38]:** `(*confirm which (M1,M2) states are "clean" (exact naive product formula) vs genuinely mixed, for (λ1,λ2)=(1,1) and (1,2)*)↵cleanCheck11 = Table[{M1,M2,Simplify[buildPsi2Raw[1,1][{M1,M2}][[1]]-tEigenva`
**Out:** (no output)


## 2026-07-14 18:06:44
**Input:**
```wolfram
{cleanCheck11, cleanCheck12}
```
**Output:**
```
Out= {{{{0, 0, True}, {0, 1, False}}, {{1, 0, False}, {1, 1, True}}}, {{{0, 0, True}, {0, 1, False}, {0, 2, False}}, {{1, 0, False}, {1, 1, False}, {1, 2, True}}}}
```
## 2026-07-14 18:06:54 — 📥 BULK INSERT 5 CELL(S) at positions 39–43
1. [markdown] (*STEP 3: normalize Ψ2[λ1,λ2][{M1,M2}] so that xCorrected[λ1,λ2][{0,0}].Ψ2[λ1,λ2][{M1,M2}]==1, mirro
2. [code] ClearAll[Ψ2];↵Ψ2[λ1_,λ2_][{M1_,M2_}] := Module[{data,eigU,raw,overlap},↵  data = buildPsi2Raw[λ1,λ2]
3. [markdown] (*STEP 4: sanity checks. (i) Ψ2 is a genuine eigenvector of t2[1,1][u] (after twistSubstL2), identic
4. [code] tmU11 = t2[1,1][1,1][u]/.twistSubstL2;↵eigenvectorCheck11 = Table[↵  Simplify[tmU11.Ψ2[1,1][{M1,M2}]
5. [code] tmU12 = t2[1,2][1,1][u]/.twistSubstL2;↵eigenvectorCheck12 = Table[↵  Simplify[tmU12.Ψ2[1,2][{M1,M2}]


## 2026-07-14 18:06:54
**Input:**
```wolfram
ClearAll[Ψ2];
Ψ2[λ1_,λ2_][{M1_,M2_}] := Module[{data,eigU,raw,overlap},
  data = buildPsi2Raw[λ1,λ2][{M1,M2}];
  eigU = data[[1]]; raw = data[[2]];
  overlap = Simplify[First[xCorrected[λ1,λ2][{0,0}].raw]];
  If[overlap===0, $Failed, raw/overlap]
];
```
**Output:**
```
(no output)
```

## 2026-07-14 18:06:55
**Input:**
```wolfram
tmU11 = t2[1,1][1,1][u]/.twistSubstL2;
eigenvectorCheck11 = Table[
  Simplify[tmU11.Ψ2[1,1][{M1,M2}] - buildPsi2Raw[1,1][{M1,M2}][[1]] Ψ2[1,1][{M1,M2}]],
{M1,0,1},{M2,0,1}]//Flatten//Union;
normCheck11 = Table[Simplify[First[xCorrected[1,1][{0,0}].Ψ2[1,1][{M1,M2}]]-1],{M1,0,1},{M2,0,1}]//Flatten//Union;
{eigenvectorCheck11, normCheck11}
```
**Output:**
```
(no output)
```

## 2026-07-14 18:07:24
**Input:**
```wolfram
tmU12 = t2[1,2][1,1][u]/.twistSubstL2;
eigenvectorCheck12 = Table[
  Simplify[tmU12.Ψ2[1,2][{M1,M2}] - buildPsi2Raw[1,2][{M1,M2}][[1]] Ψ2[1,2][{M1,M2}]],
{M1,0,1},{M2,0,2}]//Flatten//Union;
normCheck12 = Table[Simplify[First[xCorrected[1,2][{0,0}].Ψ2[1,2][{M1,M2}]]-1],{M1,0,1},{M2,0,2}]//Flatten//Union;
{eigenvectorCheck12, normCheck12}
```
**Output:**
```
(no output)
```
## 2026-07-14 18:08:10 — ▶️ RUN CELL 43
**In [43]:** `tmU12 = t2[1,2][1,1][u]/.twistSubstL2;↵eigenvectorCheck12 = Table[↵  Simplify[tmU12.Ψ2[1,2][{M1,M2}] - buildPsi2Raw[1,2][{M1,M2}][[1]] Ψ2[1,2][{M1,M2}]],↵{M1,0,1},{M2,0,2}]//Flatten//Union;↵normCheck1`
**Out:** (no output)


## 2026-07-14 18:08:12
**Input:**
```wolfram
{eigenvectorCheck11, normCheck11}
{eigenvectorCheck12, normCheck12}
```
**Output:**
```
Out[1]= {{0}, {0}}
Out[2]= {{0}, {0}}
```

---

## 2026-07-14 18:08:15 — 🔄 KERNEL RESTART

## 2026-07-14 18:11:10 — ▶️ RUN CELLS 1–43
Cell 2: ✓ — (no output)
Cell 3: ✓ — (no output)
Cell 4: ✓ — (no output)
Cell 5: ✓ — (no output)
Cell 7: ✓ — (no output)
Cell 8: ✓ — (no output)
Cell 10: ✓ — (no output)
Cell 11: ✓ — (no output)
Cell 12: ✓ — (no output)
Cell 14: ✓ — (no output)
Cell 15: ✓ — (no output)
Cell 17: ✓ — (no output)
Cell 18: ✓ — (no output)
Cell 20: ✓ — (no output)
Cell 22: ✓ — (no output)
Cell 23: ✓ — (no output)
Cell 25: ✓ — (no output)
Cell 27: ✓ — (no output)
Cell 28: ✓ — (no output)
Cell 30: ✓ — (no output)
Cell 32: ✓ — (no output)
Cell 33: ✓ — (no output)
Cell 35: ✓ — (no output)
Cell 37: ✓ — (no output)
Cell 38: ✓ — (no output)
Cell 40: ✓ — (no output)
Cell 42: ✓ — (no output)
Cell 43: ✓ — (no output)


## 2026-07-14 18:11:13
**Input:**
```wolfram
summaryResult
{eigenvectorCheck11, normCheck11}
{eigenvectorCheck12, normCheck12}
```
**Output:**
```
Out[1]= {0}
Out[2]= {{0}, {0}}
Out[3]= {{0}, {0}}
```
## 2026-07-14 18:22:20 — 🗑️ DELETE MARKDOWN CELL 26
(*Step 9 — THE KEY CHECK: candidate SoV covector basis x2[λ1,λ2][{s1,s2}] built as the naive tensor product of the single-site fused-transfer-matrix covectors, x2[λ1,λ2][{s1,s2}] = xSingle[λ1,θ1][{s1}…

## 2026-07-14 18:22:20 — 🗑️ DELETE CODE CELL 27
ClearAll[x2];
x2[λ1_,λ2_][{s1_,s2_}]:=KroneckerProduct[xSingle[λ1,θ1][{s1}],xSingle[λ2,θ2][{s2}]];

## 2026-07-14 18:22:20 — 🗑️ DELETE CODE CELL 28
(*direct eigenvector check for λ1=λ2=1: all 4 candidate states {s1,s2}∈{0,1}x{0,1}*)
case11=Table[{s1,s2,Simplify[x2[1,1][{s1,s2}].B2[1,1][u]-(u-θ1-s1 h)(u-θ2-s2 h) x2[1,1][{s1,s2}]]//Flatten//Union},…

## 2026-07-14 18:22:20 — 🗑️ DELETE MARKDOWN CELL 29
(*RESULT of the naive-tensor-product test: diagonalizing B2[1,1][u] directly (Eigensystem of Transpose[B2[u]], symbolic in u) confirms the eigenvalue guess (u-θ1-s1 h)(u-θ2-s2 h) is exactly right for …

## 2026-07-14 18:22:20 — 🗑️ DELETE CODE CELL 30
(*direct diagonalization confirms this for λ1=λ2=1*)
MT=Transpose[B2[1,1][u]];
{evals,evecs}=Simplify[{Eigenvalues[MT],Eigenvectors[MT]}]

## 2026-07-14 18:22:20 — 🗑️ DELETE CODE CELL 38
(*confirm which (M1,M2) states are "clean" (exact naive product formula) vs genuinely mixed, for (λ1,λ2)=(1,1) and (1,2)*)
cleanCheck11 = Table[{M1,M2,Simplify[buildPsi2Raw[1,1][{M1,M2}][[1]]-tEigenva…

## 2026-07-14 18:22:20 — ✏️ EDIT CELL 26
(*The naive tensor-product ansatz \[LeftAngleBracket]s1,s2| = xSingle[\[Lambda]1,\[Theta]1][{s1}]\[CircleTimes]xSingle[\[Lambda]2,\[Theta]2][{s2}] is an exact B2(u)-eigencovector (eigenvalue (u-\[Thet…


---

## 2026-07-14 18:22:27 — 🔄 KERNEL RESTART

## 2026-07-14 18:25:09 — ▶️ RUN CELLS 1–37
Cell 2: ✓ — (no output)
Cell 3: ✓ — (no output)
Cell 4: ✓ — (no output)
Cell 5: ✓ — (no output)
Cell 7: ✓ — (no output)
Cell 8: ✓ — (no output)
Cell 10: ✓ — (no output)
Cell 11: ✓ — (no output)
Cell 12: ✓ — (no output)
Cell 14: ✓ — (no output)
Cell 15: ✓ — (no output)
Cell 17: ✓ — (no output)
Cell 18: ✓ — (no output)
Cell 20: ✓ — (no output)
Cell 22: ✓ — (no output)
Cell 23: ✓ — (no output)
Cell 25: ✓ — (no output)
Cell 27: ✓ — (no output)
Cell 28: ✓ — (no output)
Cell 30: ✓ — (no output)
Cell 32: ✓ — (no output)
Cell 34: ✓ — (no output)
Cell 36: ✓ — (no output)
Cell 37: ✓ — (no output)


## 2026-07-14 18:25:12
**Input:**
```wolfram
summaryResult
{eigenvectorCheck11, normCheck11}
{eigenvectorCheck12, normCheck12}
```
**Output:**
```
Out[1]= {0}
Out[2]= {{0}, {0}}
Out[3]= {{0}, {0}}
```
## 2026-07-14 18:30:58 — ✏️ EDIT CELL 27 [batch]
ClearAll[x2];
x2[λ1_,λ2_][{0,0}] := KroneckerProduct[xSingle[λ1,θ1][{0}],xSingle[λ2,θ2][{0}]];
x2[λ1_,λ2_][{s1_,s2_}] := x2[λ1,λ2][{s1,s2}] = x2[λ1,λ2][{0,0}] . t2[λ1,λ2][1,s1][θ1] . t2[λ1,λ2][1,s2][θ…

## 2026-07-14 18:30:58 — ✏️ EDIT CELL 28 [batch]
(*final verification pass: x2 is an EXACT B2(u) eigencovector, with eigenvalue (u-\[Theta]1-s1 h)(u-\[Theta]2-s2 h), for all (s1,s2) at (\[Lambda]1,\[Lambda]2)\[Element]{(1,1),(1,2),(2,1),(2,2)} — no …

## 2026-07-14 18:30:58 — ✏️ EDIT CELL 33 [batch]
(*STEP 3: normalize Ψ2[λ1,λ2][{M1,M2}] so that x2[λ1,λ2][{0,0}].Ψ2[λ1,λ2][{M1,M2}]==1, mirroring x[λ][{0}].Ψ[λ,M]==1 at L=1.*)

## 2026-07-14 18:30:58 — ✏️ EDIT CELL 34 [batch]
ClearAll[Ψ2];
Ψ2[λ1_,λ2_][{M1_,M2_}] := Module[{data,eigU,raw,overlap},
  data = buildPsi2Raw[λ1,λ2][{M1,M2}];
  eigU = data[[1]]; raw = data[[2]];
  overlap = Simplify[First[x2[λ1,λ2][{0,0}].raw]];
 …

## 2026-07-14 18:30:58 — ✏️ EDIT CELL 35 [batch]
(*STEP 4: sanity checks. (i) Ψ2 is a genuine eigenvector of t2[1,1][u] (after twistSubstL2), identically in u, with the eigenvalue read off in buildPsi2Raw — checked for (λ1,λ2)=(1,1) and (1,2), all s…

## 2026-07-14 18:30:58 — ✏️ EDIT CELL 36 [batch]
tmU11 = t2[1,1][1,1][u]/.twistSubstL2;
eigenvectorCheck11 = Table[
  Simplify[tmU11.Ψ2[1,1][{M1,M2}] - buildPsi2Raw[1,1][{M1,M2}][[1]] Ψ2[1,1][{M1,M2}]],
{M1,0,1},{M2,0,1}]//Flatten//Union;
normCheck1…

## 2026-07-14 18:30:58 — ✏️ EDIT CELL 37 [batch]
tmU12 = t2[1,2][1,1][u]/.twistSubstL2;
eigenvectorCheck12 = Table[
  Simplify[tmU12.Ψ2[1,2][{M1,M2}] - buildPsi2Raw[1,2][{M1,M2}][[1]] Ψ2[1,2][{M1,M2}]],
{M1,0,1},{M2,0,2}]//Flatten//Union;
normCheck1…

## 2026-07-14 18:31:15 — ✏️ EDIT CELL 27
x2[λ1_,λ2_][{s1_,s2_}] := x2[λ1,λ2][{s1,s2}] = x2[λ1,λ2][{0,0}] . t2[λ1,λ2][1,s1][θ1] . t2[λ1,λ2][1,s2][θ2];


---

## 2026-07-14 18:31:21 — 🔄 KERNEL RESTART

## 2026-07-14 18:34:03 — ▶️ RUN CELLS 1–37
Cell 2: ✓ — (no output)
Cell 3: ✓ — (no output)
Cell 4: ✓ — (no output)
Cell 5: ✓ — (no output)
Cell 7: ✓ — (no output)
Cell 8: ✓ — (no output)
Cell 10: ✓ — (no output)
Cell 11: ✓ — (no output)
Cell 12: ✓ — (no output)
Cell 14: ✓ — (no output)
Cell 15: ✓ — (no output)
Cell 17: ✓ — (no output)
Cell 18: ✓ — (no output)
Cell 20: ✓ — (no output)
Cell 22: ✓ — (no output)
Cell 23: ✓ — (no output)
Cell 25: ✓ — (no output)
Cell 27: ✓ — (no output)
Cell 28: ✓ — (no output)
Cell 30: ✓ — (no output)
Cell 32: ✓ — (no output)
Cell 34: ✓ — (no output)
Cell 36: ✓ — (no output)
Cell 37: ✓ — (no output)


## 2026-07-14 18:34:06
**Input:**
```wolfram
summaryResult
{eigenvectorCheck11, normCheck11}
{eigenvectorCheck12, normCheck12}
```
**Output:**
```
Out[1]= {0}
Out[2]= {{0}, {0}}
Out[3]= {{0}, {0}}
```

## 2026-07-14 19:56:34
**Input:**
```wolfram
(* symbolic sanity check of the M+1 and M+2 degree cancellations, BEFORE assigning numerics *)
Block[{θ1,θ2,h,κ1,κ2,λ1,λ2,M,ee1,ee2,u,k,a,d,τ0,O,cl},
  ee1=λ1+λ2-M; ee2=M;
  a=(u-θ1-h λ1)(u-θ2-h λ2);
  d=(u-θ1)(u-θ2);
  τ0=(κ1+κ2)u^2-(θ1+θ2)(κ1+κ2)u-u h(κ1 ee1+κ2 ee2);
  O=κ1 a (u+h)^k+κ2 d (u-h)^k-τ0 u^k;
  (* coefficient of u^(k+2) and u^(k+1) in O, general k *)
  {Simplify[Coefficient[Expand[O],u,k+2]],
   Simplify[Coefficient[Expand[O],u,k+1]]}
]
```
**Output:**
```
Out= {-κ1 - κ2, h*M*(-κ1 + κ2) + (θ1 + θ2)*(κ1 + κ2) + h*κ1*(λ1 + λ2)}
```

## 2026-07-14 19:56:51
**Input:**
```wolfram
Block[{θ1,θ2,h,κ1,κ2,λ1=1,λ2=1,u},
 Table[
  Module[{ee1=λ1+λ2-M,ee2=M,a,d,τ0,cl},
   a=(u-θ1-h λ1)(u-θ2-h λ2);
   d=(u-θ1)(u-θ2);
   τ0=(κ1+κ2)u^2-(θ1+θ2)(κ1+κ2)u-u h(κ1 ee1+κ2 ee2);
   (* for each k, coeffs at degree M+1 and M+2 *)
   Table[
    cl=PadRight[CoefficientList[Expand[κ1 a (u+h)^k+κ2 d (u-h)^k-τ0 u^k],u],M+3];
    {M,k,Simplify[cl[[M+2]]],Simplify[cl[[M+3]]]},
    {k,0,M}]
  ],{M,0,3}]//Column]
```
**Output:**
```
Out= Column[{{{0, 0, 0, 0}}, {{1, 0, 0, 0}, {1, 1, 0, 0}}, {{2, 0, 0, 0}, {2, 1, 0, 0}, {2, 2, 0, 0}}, {{3, 0, 0, 0}, {3, 1, 0, 0}, {3, 2, 0, 0}, {3, 3, 0, 0}}}]
```
## 2026-07-14 19:57:14 — 📥 BULK INSERT 3 CELL(S) at positions 25–27
1. [markdown] (*Numerical solution of the L=2 Baxter (TQ) equation.↵↵Baxter equation:  a(u) Q(u+h) - τ(u) Q(u) + d
2. [code] (*Numeric parameters. NOTE: this assigns global numeric values to θ1,θ2,h,κ1,κ2 (previously symbolic
3. [code] (*a(u); d(u) is Qθ[u] from Cell 19 (=(u-θ1)(u-θ2)) — reused, no redundant dFun. τ0 = τ minus the J t


## 2026-07-14 19:57:14
**Input:**
```wolfram
(*Numeric parameters. NOTE: this assigns global numeric values to θ1,θ2,h,κ1,κ2 (previously symbolic). Cells 1-24 are already verified and are not re-run.*)
θ1 = 1/3;
θ2 = 1/7;
h = 1;
κ1 = N[Exp[I Zeta[3]], 30];
κ2 = 1/κ1;
{θ1, θ2, h, κ1, κ2}
```
**Output:**
```
Out[228]= \left\{\frac{1}{3},\frac{1}{7},1,0.360439875147226537266677650596+0.932782448593374716484835982071\,\mathrm{i},0.36043987514722653726667765060-0.932782448593374716484835982071\,\mathrm{i}\right\}
```

## 2026-07-14 19:57:15
**Input:**
```wolfram
(*a(u); d(u) is Qθ[u] from Cell 19 (=(u-θ1)(u-θ2)) — reused, no redundant dFun. τ0 = τ minus the J term, using ee1=λ1+λ2-M, ee2=M.*)
ClearAll[aFun, τ0];
aFun[λ1_, λ2_][u_] := (u - θ1 - h λ1) (u - θ2 - h λ2);
τ0[λ1_, λ2_, M_][u_] := (κ1 + κ2) u^2 - (θ1 + θ2) (κ1 + κ2) u - u h (κ1 (λ1 + λ2 - M) + κ2 M);
{aFun[1, 1][u], Qθ[u], τ0[1, 1, 1][u]}
```
**Output:**
```
Out[232]= \left\{\left(u-\frac{4}{3}\right)\,\left(u-\frac{8}{7}\right),\left(u-\frac{1}{3}\right)\,\left(u-\frac{1}{7}\right),\left(0.720879750294453074533355301191+0.\times 10^{-30}\,\mathrm{i}\right)\,u^2-\left(1.06415582186324025288257211128+0.\times 10^{-29}\,\mathrm{i}\right)\,u\right\}
```
## 2026-07-14 19:57:34 — 📥 BULK INSERT 2 CELL(S) at positions 28–29
1. [code] (*BaxterMatrix: build O[u^k]=κ1 a(u)(u+h)^k + κ2 Qθ[u](u-h)^k - τ0(u) u^k for k=0..M.↵  HARD ASSERT 
2. [code] (*Explicit display of the high-degree (u^(M+1), u^(M+2)) coefficients for λ1=λ2=1, M=0,1,2 —↵  confi


## 2026-07-14 19:57:35
**Input:**
```wolfram
(*BaxterMatrix: build O[u^k]=κ1 a(u)(u+h)^k + κ2 Qθ[u](u-h)^k - τ0(u) u^k for k=0..M.
  HARD ASSERT that the degree M+1 and M+2 coefficients vanish (this is where a wrong a/d/τ or a
  wrong ee1=λ1+λ2-M convention would show up). Then columns = degree-0..M coefficients.*)
ClearAll[BaxterMatrix];
BaxterMatrix[λ1_, λ2_, M_] := Module[{clists, high, maxhigh, cols},
  clists = Table[
     PadRight[
      CoefficientList[
       Expand[κ1 aFun[λ1, λ2][u] (u + h)^k + κ2 Qθ[u] (u - h)^k - τ0[λ1, λ2, M][u] u^k],
       u],
      M + 3],
     {k, 0, M}];
  high = clists[[All, {M + 2, M + 3}]];   (* degree M+1 and M+2 coeffs for each column k *)
  maxhigh = Max[Abs[Flatten[high]]];
  If[maxhigh > 10^-18,
    Print["BaxterMatrix ASSERTION FAILED (λ1,λ2,M)=", {λ1, λ2, M},
          "  max|high-degree coeff|=", maxhigh, "  high=", high];
    Abort[]];
  cols = clists[[All, 1 ;; M + 1]];       (* degrees 0..M, one column per k *)
  Transpose[cols]
];
Dimensions /@ (BaxterMatrix[1, 1, #] & /@ {0, 1, 2})
```
**Output:**
```
Out[235]= \begin{pmatrix}1 & 1\\2 & 2\\3 & 3\end{pmatrix}
```

## 2026-07-14 19:57:35
**Input:**
```wolfram
(*Explicit display of the high-degree (u^(M+1), u^(M+2)) coefficients for λ1=λ2=1, M=0,1,2 —
  confirms the cancellation the assertion relies on (should all be ~10^-30).*)
Table[
 Module[{clists},
  clists = Table[
     PadRight[CoefficientList[Expand[κ1 aFun[1, 1][u] (u + h)^k + κ2 Qθ[u] (u - h)^k - τ0[1, 1, M][u] u^k], u], M + 3],
     {k, 0, M}];
  {M, "max|deg M+1,M+2 coeff|" -> Max[Abs[Flatten[clists[[All, {M + 2, M + 3}]]]]]}],
 {M, 0, 2}] // Column
```
**Output:**
```
Out[236]= \begin{matrix}\left\{0,\text{max|deg M+1,M+2 coeff|}\to 0.\times 10^{-29}\right\}\\\left\{1,\text{max|deg M+1,M+2 coeff|}\to 0.\times 10^{-29}\right\}\\\left\{2,\text{max|deg M+1,M+2 coeff|}\to 0.\times 10^{-29}\right\}\end{matrix}
```
## 2026-07-14 19:57:58 — 📥 BULK INSERT 3 CELL(S) at positions 30–32
1. [code] (*BaxterSolve: diagonalize BaxterMatrix. Eigenvalues are the admissible J_n; eigenvectors are↵  the 
2. [code] (*INDEPENDENT residual check: re-evaluate the FULL Baxter equation from scratch (using τ=τ0+J_n,↵  N
3. [code] (*Test run: full (M,n,J_n,q_n(u)) table for λ1=λ2=1 (so λ1+λ2=2, M=0,1,2). Values rounded to 10↵  di


## 2026-07-14 19:57:59
**Input:**
```wolfram
(*BaxterSolve: diagonalize BaxterMatrix. Eigenvalues are the admissible J_n; eigenvectors are
  the ascending coefficient vectors of q_n, renormalized so the u^M (top) component is 1 (monic).
  Returns {{J_0,q_0},...,{J_M,q_M}} with q_n = {c0,...,cM} (cM=1), n=0..M.*)
ClearAll[BaxterSolve, qpoly];
qpoly[clist_, uu_] := clist . uu^Range[0, Length[clist] - 1];   (* build polynomial from ascending coeff list *)
BaxterSolve[λ1_, λ2_, M_] := Module[{mat, vals, vecs},
  mat = BaxterMatrix[λ1, λ2, M];
  {vals, vecs} = Eigensystem[mat];
  Table[{vals[[n]], vecs[[n]]/vecs[[n, M + 1]]}, {n, 1, M + 1}]
];
(*quick look at the M=1 sector*)
BaxterSolve[1, 1, 1]
```
**Output:**
```
Out[240]= \begin{pmatrix}1.1900552523337025695566871451+6.5302155526567322175297423340\times 10^{-40}\,\mathrm{i} & \left\{-0.73809523809523809523809523810+0.72070903915869767537485500365\,\mathrm{i},1.\,+0.\times 10^{-29}\,\mathrm{i}\right\}\\-0.77812396645115795553762697294-6.5302155526567322127504666006\times 10^{-40}\,\mathrm{i} & \left\{-0.73809523809523809523809523810-0.334295384301913272877… [473 chars total]
```

## 2026-07-14 19:57:59
**Input:**
```wolfram
(*INDEPENDENT residual check: re-evaluate the FULL Baxter equation from scratch (using τ=τ0+J_n,
  NOT any BaxterMatrix internals) at several u that are not roots of a(u) or d(u).
  a(u) roots: 4/3, 8/7 ; d(u) roots: 1/3, 1/7 ; sample u in {0,1,3/2,-2,10} avoid all of these.
  Exact rationals are used for u so only κ1 carries finite (30-digit) precision.*)
sampleU = {0, 1, 3/2, -2, 10};
resTable = Flatten[
   Table[
    Module[{sols = BaxterSolve[1, 1, M]},
     Table[
      With[{Jn = sols[[n, 1]], cl = sols[[n, 2]]},
       {M, n - 1,
        Max[Abs[Table[
            κ1 aFun[1, 1][uu] qpoly[cl, uu + h]
            - (τ0[1, 1, M][uu] + Jn) qpoly[cl, uu]
            + κ2 Qθ[uu] qpoly[cl, uu - h],
            {uu, sampleU}]]]}],
      {n, 1, M + 1}]],
    {M, 0, 2}],
   1];
maxResidual = Max[resTable[[All, 3]]];
Print["max residual over all (M,n) and sample u = ", maxResidual];
Print[If[maxResidual < 10^-18, "RESIDUAL CHECK PASSED", "RESIDUAL CHECK FAILED"]];
TableForm[resTable, TableHeadings -> {None, {"M", "n", "max|residual|"}}]
```
**Output:**
```
2 kernel messages: Power::indet: Indeterminate expression 0  encountered.
Power::indet: Indeterminate expression 0  encountered.
Power::indet: Indeterminate expression 0  encountered.
General::stop: Further output of Po… [also: Power::indet] result: Print: max residual over all (M,n) and sample u =[formula] | Print: [formula] | Out[246]//TableForm= \begin{pmatrix}0 & 0 & \mathrm{Indeterminate}\\1 & 0 & \mathrm{Indeterminate}\\1 & 1 & \mathrm{Inde…
```

## 2026-07-14 19:58:00
**Input:**
```wolfram
(*Test run: full (M,n,J_n,q_n(u)) table for λ1=λ2=1 (so λ1+λ2=2, M=0,1,2). Values rounded to 10
  digits for display; full 30-digit values are available from BaxterSolve[1,1,M].*)
resultTable = Flatten[
   Table[
    Module[{sols = BaxterSolve[1, 1, M]},
     Table[{M, n - 1, N[sols[[n, 1]], 10], N[qpoly[sols[[n, 2]], u], 10]}, {n, 1, M + 1}]],
    {M, 0, 2}],
   1];
TableForm[resultTable, TableHeadings -> {None, {"M", "n", "J_n", "q_n(u)"}}]
```
**Output:**
```
Out[248]//TableForm= \begin{pmatrix}0 & 0 & 0.5664055181+1.376964567\,\mathrm{i} & 1.\,\\1 & 0 & 1.190055252+0.\times 10^{-39}\,\mathrm{i} & \left(1.\,+0.\times 10^{-10}\,\mathrm{i}\right)\,u-\left(0.7380952381-0.7207090392\,\mathrm{i}\right)\\1 & 1 & -0.7781239665+0.\times 10^{-39}\,\mathrm{i} & 1.\,\,u-\left(0.7380952381+0.3342953843\,\mathrm{i}\right)\\2 & 0 & -0.154474232-1.732310262\,\mathrm{… [1004 chars total]
```
## 2026-07-14 19:58:31 — ✏️ EDIT CELL 30
(*BaxterSolve: diagonalize BaxterMatrix. Eigenvalues are the admissible J_n; eigenvectors are
  the ascending coefficient vectors of q_n, renormalized so the u^M (top) component is 1 (monic).
  Return…


## 2026-07-14 19:58:31
**Input:**
```wolfram
(*BaxterSolve: diagonalize BaxterMatrix. Eigenvalues are the admissible J_n; eigenvectors are
  the ascending coefficient vectors of q_n, renormalized so the u^M (top) component is 1 (monic).
  Returns {{J_0,q_0},...,{J_M,q_M}} with q_n = {c0,...,cM} (cM=1), n=0..M.*)
ClearAll[BaxterSolve, qpoly];
(*build polynomial from ascending coeff list; form it in a dummy symbol first so that
  substituting uu=0 gives the constant term rather than an indeterminate 0^0.*)
qpoly[clist_, uu_] := Module[{x}, (clist . x^Range[0, Length[clist] - 1]) /. x -> uu];
BaxterSolve[λ1_, λ2_, M_] := Module[{mat, vals, vecs},
  mat = BaxterMatrix[λ1, λ2, M];
  {vals, vecs} = Eigensystem[mat];
  Table[{vals[[n]], vecs[[n]]/vecs[[n, M + 1]]}, {n, 1, M + 1}]
];
(*quick look at the M=1 sector*)
BaxterSolve[1, 1, 1]
```
**Output:**
```
Out= {{1.1900552523337025695566871450502370524528514334087876680712`29.272055931994913 + 6.53021555265673221752974233400125296228677658315650892484`29.272055931994913*^-40*I, {-0.7380952380952380952380952380952380952405638156481974163032`28.66999594066695 + 0.7207090391586976753748550036528131217877401778330615040709`28.66999594066695*I, 1.`28.66999594066695 + 0``28.670119316342326*I}}, {-0.7781239664511579555376269729410369481171205496804619509347`29.272055931994913 - 6.53021555265673221275046660060389613120013538985337350807`29.272055931994913*^-40*I, {-0.7380952380952380952380952380952380952388370373758838536617`28.971025936330932 - 0.3342953843019132728770027485554114090640161353343750603849`28.971025936330932*I, 0.9999999999999999999999999999999999999999999999999999999998`28.971025936330932}}}
```
## 2026-07-14 19:58:38 — ▶️ RUN CELL 31 (0.12 s)
**In [31]:** `(*INDEPENDENT residual check: re-evaluate the FULL Baxter equation from scratch (using τ=τ0+J_n,↵  NOT any BaxterMatrix internals) at several u that are not roots of a(u) or d(u).↵  a(u) roots: 4/3, 8`
**Out:** Print: max residual over all (M,n) and sample u =[formula] | Print: RESIDUAL CHECK PASSED | Out[254]//TableForm= \begin{pmatrix}0 & 0 & 0.\times 10^{-27}\\1 & 0 & 0.\times 10^{-25}\\1 & 1 & 0.\times 10^{-26}\\2 & 0 & 0.\times 10^{-24}\\2 & 1 & 0.\times 10^{-24}\\2 & 2 & 0.\times 10^{-24}\end{pmatrix

## 2026-07-14 19:58:46 — ▶️ RUN CELL 32 (0.06 s)
**In [32]:** `(*Test run: full (M,n,J_n,q_n(u)) table for λ1=λ2=1 (so λ1+λ2=2, M=0,1,2). Values rounded to 10↵  digits for display; full 30-digit values are available from BaxterSolve[1,1,M].*)↵resultTable = Flatte`
**Out:** Out[256]//TableForm= \begin{pmatrix}0 & 0 & 0.5664055181+1.376964567\,\mathrm{i} & 1.\,\\1 & 0 & 1.190055252+0.\times 10^{-39}\,\mathrm{i} & \left(1.\,+0.\times 10^{-10}\,\mathrm{i}\right)\,u-\left(0.7380952381-0.7207090392\,\mathrm{i}\right)\\1 & 1 & -0.7781239665+0.\times 10^{-39}\,\mathrm{i} & 1.


## 2026-07-14 19:59:03
**Input:**
```wolfram
maxResidual
```
**Output:**
```
Out= 0``23.88453537724287
```
## 2026-07-14 20:07:37 — ✏️ EDIT CELL 25
(*Numerical solution of the L=2 Baxter (TQ) equation.

Baxter equation:  a(u) Q(u+h) - τ(u) Q(u) + d(u) Q(u-h) = 0, with
  a(u) = (u-θ1-h λ1)(u-θ2-h λ2)   [B2 eigenvalue at s1=λ1,s2=λ2, top of SoV lat…

## 2026-07-14 20:07:41 — ✏️ EDIT CELL 30
(*BaxterSolve: diagonalize BaxterMatrix. Eigenvalues are the admissible J_n; eigenvectors are
  the ascending coefficient vectors of q_n, renormalized so the u^M (top) component is 1 (monic).
  FIXED …


## 2026-07-14 20:07:41
**Input:**
```wolfram
(*BaxterSolve: diagonalize BaxterMatrix. Eigenvalues are the admissible J_n; eigenvectors are
  the ascending coefficient vectors of q_n, renormalized so the u^M (top) component is 1 (monic).
  FIXED n-ordering: the M+1 solutions are sorted by ascending Re[J] (primary key) then ascending
  Im[J] (tie-break); n=0..M is assigned in that sorted order — deterministic and reproducible.
  MEMOIZED per (λ1,λ2,M) so that Q1 and τ always read the identical cached eigensystem.
  Returns {{J_0,q_0},...,{J_M,q_M}} with q_n = {c0,...,cM} (cM=1), n=0..M.*)
ClearAll[BaxterSolve, qpoly];
(*build polynomial from ascending coeff list; form it in a dummy symbol first so that
  substituting uu=0 gives the constant term rather than an indeterminate 0^0.*)
qpoly[clist_, uu_] := Module[{x}, (clist . x^Range[0, Length[clist] - 1]) /. x -> uu];
BaxterSolve[λ1_, λ2_, M_] := BaxterSolve[λ1, λ2, M] = Module[{mat, vals, vecs, pairs},
  mat = BaxterMatrix[λ1, λ2, M];
  {vals, vecs} = Eigensystem[mat];
  pairs = Table[{vals[[n]], vecs[[n]]/vecs[[n, M + 1]]}, {n, 1, M + 1}];
  (*sort into the fixed (Re[J] asc, then Im[J] asc) convention*)
  SortBy[pairs, {Re[First[#]], Im[First[#]]} &]
];
(*quick look at the M=1 sector*)
BaxterSolve[1, 1, 1]
```
**Output:**
```
Out= {{-0.7781239664511579555376269729410369481171205496804619509347`29.272055931994913 - 6.53021555265673221275046660060389613120013538985337350807`29.272055931994913*^-40*I, {-0.7380952380952380952380952380952380952388370373758838536617`28.971025936330932 - 0.3342953843019132728770027485554114090640161353343750603849`28.971025936330932*I, 0.9999999999999999999999999999999999999999999999999999999998`28.971025936330932}}, {1.1900552523337025695566871450502370524528514334087876680712`29.272055931994913 + 6.53021555265673221752974233400125296228677658315650892484`29.272055931994913*^-40*I, {-0.7380952380952380952380952380952380952405638156481974163032`28.66999594066695 + 0.7207090391586976753748550036528131217877401778330615040709`28.66999594066695*I, 1.`28.66999594066695 + 0``28.670119316342326*I}}}
```
## 2026-07-14 20:09:14 — 📥 BULK INSERT 4 CELL(S) at positions 33–36
1. [code] (*On-demand wrappers over the memoized BaxterSolve cache.↵  Q1[λ1,λ2][M,n][u] = κ1^(u/h) q_n(u)     
2. [code] (*VALIDATION 5a: end-to-end check of the ORIGINAL (pre-ansatz-substitution) Baxter equation↵    a(u)
3. [code] (*VALIDATION 5b: genuine generalization to a NEW, previously-untested representation λ1=2, λ2=1↵  (a
4. [code] (*VALIDATION 5c: confirm the fixed n-ordering actually holds — Re[J_n] must be NON-DECREASING in n↵


## 2026-07-14 20:09:14
**Input:**
```wolfram
(*On-demand wrappers over the memoized BaxterSolve cache.
  Q1[λ1,λ2][M,n][u] = κ1^(u/h) q_n(u)      — the actual Baxter Q-function value at u.
  τ[λ1,λ2][M,n][u]  = τ0[λ1,λ2,M][u] + J_n — the FULL τ eigenvalue for that (M,n) solution.
  Both pull from the SAME cached BaxterSolve[λ1,λ2,M] (n is 0-indexed → list entry n+1), so they
  are guaranteed mutually consistent — neither recomputes independently.
  Bounds: M an integer in 0..λ1+λ2, n an integer in 0..M; violation → clear message + $Failed.*)
ClearAll[Q1, τ];
Q1::bnds = "Invalid indices (M,n)=(`1`,`2`) for (λ1,λ2)=(`3`,`4`): need M an integer in 0..λ1+λ2 and n an integer in 0..M.";
τ::bnds = "Invalid indices (M,n)=(`1`,`2`) for (λ1,λ2)=(`3`,`4`): need M an integer in 0..λ1+λ2 and n an integer in 0..M.";
Q1[λ1_, λ2_][M_, n_][u_] := If[
   IntegerQ[M] && IntegerQ[n] && 0 <= M <= λ1 + λ2 && 0 <= n <= M,
   κ1^(u/h) qpoly[BaxterSolve[λ1, λ2, M][[n + 1, 2]], u],
   Message[Q1::bnds, M, n, λ1, λ2]; $Failed];
τ[λ1_, λ2_][M_, n_][u_] := If[
   IntegerQ[M] && IntegerQ[n] && 0 <= M <= λ1 + λ2 && 0 <= n <= M,
   τ0[λ1, λ2, M][u] + BaxterSolve[λ1, λ2, M][[n + 1, 1]],
   Message[τ::bnds, M, n, λ1, λ2]; $Failed];
(*smoke test: two valid values, plus a deliberately out-of-range call to show the guard fires*)
Print["bounds-guard test (expect a Q1::bnds message and $Failed): ", Q1[1, 1][3, 0][2]];
{Q1[1, 1][1, 0][2], τ[1, 1][1, 0][2]}
```
**Output:**
```
2 kernel messages: Q1::bnds: Invalid indices (M,n)=(3,0) for (λ1,λ2)=(1,1): need M an integer in 0..λ1+λ2 and n an integer in 0..M. [also: Q1::bnds] result: Print: bounds-guard test (expect a Q1::bnds message and $Failed):[formula] | Out[263]= \left\{-0.7092310109580169265372251422+1.0959691625127720670675272805\,\mathrm{i},-0.0229166089998261631693499907…
```

## 2026-07-14 20:09:15
**Input:**
```wolfram
(*VALIDATION 5a: end-to-end check of the ORIGINAL (pre-ansatz-substitution) Baxter equation
    a(u) Q(u+h) - τ(u) Q(u) + d(u) Q(u-h) = 0
  assembled ENTIRELY from the wrapper functions Q1 and τ (NOT the raw q_n/J_n as Cell 31 did),
  for λ1=λ2=1, every (M,n), at several sample u. Exact-rational u so only κ1 carries finite
  (30-digit) precision. This exercises the wrappers — memoization, sorting, and the κ1^(u/h)
  prefactor — top to bottom, independently of the internal matrix-building logic.*)
sampleU5a = {0, 1, 3/2, -2, 10};
res5a = Flatten[Table[
    Table[{M, n,
      Max[Abs[Table[
         aFun[1, 1][uu] Q1[1, 1][M, n][uu + h]
         - τ[1, 1][M, n][uu] Q1[1, 1][M, n][uu]
         + Qθ[uu] Q1[1, 1][M, n][uu - h],
         {uu, sampleU5a}]]]},
     {n, 0, M}],
    {M, 0, 2}], 1];
max5a = Max[res5a[[All, 3]]];
Print["5a max residual (via Q1/τ wrappers, λ1=λ2=1) = ", max5a];
Print[If[max5a < 10^-18, "5a PASSED", "5a FAILED"]];
TableForm[res5a, TableHeadings -> {None, {"M", "n", "max|residual|"}}]
```
**Output:**
```
Print: 5a max residual (via Q1/τ wrappers, λ1=λ2=1) =[formula] | Print: 5a PASSED | Out[269]//TableForm= \begin{pmatrix}0 & 0 & 0.\times 10^{-26}\\1 & 0 & 0.\times 10^{-25}\\1 & 1 & 0.\times 10^{-25}\\2 & 0 & 0.\times 10^{-24}\\2 & 1 & 0.\times 10^{-24}\\2 & 2 & 0.\times 10^{-24}\end{pmatrix}
```

## 2026-07-14 20:09:15
**Input:**
```wolfram
(*VALIDATION 5b: genuine generalization to a NEW, previously-untested representation λ1=2, λ2=1
  (asymmetric; λ1+λ2=3, so M=0..3). Same end-to-end wrapper residual check, plus the full
  (M,n,J_n) table for this case.*)
sampleU5b = {0, 1, 3/2, -2, 10};
res5b = Flatten[Table[
    Table[{M, n,
      Max[Abs[Table[
         aFun[2, 1][uu] Q1[2, 1][M, n][uu + h]
         - τ[2, 1][M, n][uu] Q1[2, 1][M, n][uu]
         + Qθ[uu] Q1[2, 1][M, n][uu - h],
         {uu, sampleU5b}]]]},
     {n, 0, M}],
    {M, 0, 3}], 1];
max5b = Max[res5b[[All, 3]]];
Print["5b max residual (via Q1/τ wrappers, λ1=2,λ2=1) = ", max5b];
Print[If[max5b < 10^-18, "5b PASSED", "5b FAILED"]];
JnTable5b = Flatten[Table[
    Table[{M, n, N[BaxterSolve[2, 1, M][[n + 1, 1]], 10]}, {n, 0, M}],
    {M, 0, 3}], 1];
Column[{
  TableForm[res5b, TableHeadings -> {None, {"M", "n", "max|residual|"}}],
  "(M,n,J_n) table for λ1=2, λ2=1:",
  TableForm[JnTable5b, TableHeadings -> {None, {"M", "n", "J_n"}}]}]
```
**Output:**
```
Print: 5b max residual (via Q1/τ wrappers, λ1=2,λ2=1) =[formula] | Print: 5b PASSED | Out[276]= \begin{matrix}\begin{aligned}\text{M} & \text{n} & \text{max|residual|}\\0 & 0 & 0.\times 10^{-26}\\1 & 0 & 0.\times 10^{-25}\\1 & 1 & 0.\times 10^{-25}\\2 & 0 & 0.\times 10^{-24}\\2 & 1 & 0.\times 10^{-24}\\2 & 2 & 0.\times 10^{-24}\\3 & 0 & 0.\times 10^{-22}\\3 & 1 & 0.\times 10^{-23}\\3 & 2 & 0.\time… [1046 chars total]
```

## 2026-07-14 20:09:15
**Input:**
```wolfram
(*VALIDATION 5c: confirm the fixed n-ordering actually holds — Re[J_n] must be NON-DECREASING in n
  within each sector. Spot-checked on sectors with M>=2 (where >1 solution exists to order):
  (λ1,λ2,M) = (1,1,2), (2,1,2), (2,1,3).*)
order5c = Table[
   With[{sols = BaxterSolve[lp[[1]], lp[[2]], lp[[3]]]},
    {lp, N[sols[[All, 1]], 6], OrderedQ[Re[sols[[All, 1]]]]}],
   {lp, {{1, 1, 2}, {2, 1, 2}, {2, 1, 3}}}];
Print["5c: Re[J_n] non-decreasing in n for every checked sector? ", And @@ (Last /@ order5c)];
Print[If[And @@ (Last /@ order5c), "5c PASSED", "5c FAILED"]];
TableForm[order5c, TableHeadings -> {None, {"(λ1,λ2,M)", "J_n  (n=0..M)", "Re non-decreasing?"}}]
```
**Output:**
```
Print: 5c: Re[J_n] non-decreasing in n for every checked sector?[formula] | Print: 5c PASSED | Out[280]//TableForm= \begin{pmatrix}\left\{1,1,2\right\} & \left\{-0.15447-1.73231\,\mathrm{i},-0.15447-1.02162\,\mathrm{i},0.566406-1.37696\,\mathrm{i}\right\} & \mathrm{True}\\\left\{2,1,2\right\} & \left\{-0.837449-0.508617\,\mathrm{i},-0.46342-2.53184\,\mathrm{i},1.71280-0.690675\,\mathrm{i}\right\} … [599 chars total]
```

## 2026-07-14 20:10:58
**Input:**
```wolfram
(* (a) memoization is live: a cached DownValue now exists for the sectors we solved *)
memoCount = Length[DownValues[BaxterSolve]];
(* (b) determinism: two independent calls return byte-identical results *)
detOK = BaxterSolve[2,1,3] === BaxterSolve[2,1,3];
(* (c) Q1 and τ read the SAME cached J/q: reconstruct the pre-ansatz polynomial identity residual
       using J from τ and q from Q1, for a generic (2,1,2,1) at u=5/2 *)
Jfromτ = τ[2,1][2,1][5/2] - τ0[2,1,2][5/2];  (* should equal cached J_{2,1} *)
JcacheDirect = BaxterSolve[2,1,2][[2,1]];
consistency = Abs[Jfromτ - JcacheDirect];
{memoCount, detOK, consistency}
```
**Output:**
```
Out= {8, True, 0``27.86666234283269}
```

---

## 2026-07-14 20:23:03 — 🔄 KERNEL RESTART

## 2026-07-14 20:24:00 — ▶️ RUN CELLS 1–12
Cell 2: ✓ — (no output)
Cell 3: ✓ — (no output)
Cell 4: ✓ — (no output)
Cell 5: ✓ — (no output)
Cell 6: ✓ — (no output)
Cell 7: ✓ — (no output)
Cell 8: ✓ — (no output)
Cell 9: ✓ — (no output)
Cell 10: ✓ — (no output)
Cell 11: ✓ — (no output)
Cell 12: ✓ — (no output)


## 2026-07-14 20:24:37
**Input:**
```wolfram
{"defs present" -> {ValueQ[θ1], Head[aFun[1,1][u]]=!=aFun, ValueQ[Qθ[u]], MatchQ[BaxterSolve[1,1,1],{_,_}]},
 "5a max res" -> Max[Flatten[Table[Table[Abs[aFun[1,1][uu] Q1[1,1][M,n][uu+h]-τ[1,1][M,n][uu] Q1[1,1][M,n][uu]+Qθ[uu] Q1[1,1][M,n][uu-h]],{uu,{0,1,3/2,-2,10}}],{M,0,2},{n,0,M}]]],
 "5b max res" -> Max[Flatten[Table[Table[Abs[aFun[2,1][uu] Q1[2,1][M,n][uu+h]-τ[2,1][M,n][uu] Q1[2,1][M,n][uu]+Qθ[uu] Q1[2,1][M,n][uu-h]],{uu,{0,1,3/2,-2,10}}],{M,0,3},{n,0,M}]]],
 "5c ordered" -> (And @@ Table[OrderedQ[Re[BaxterSolve[lp[[1]],lp[[2]],lp[[3]]][[All,1]]]],{lp,{{1,1,2},{2,1,2},{2,1,3}}}]),
 "Qθ from local" -> (Qθ[5] === (5-θ1)(5-θ2))}
```
**Output:**
```
Out= {"defs present" -> {True, True, True, True}, "5a max res" -> 0``23.61899743977343, "5b max res" -> 0``22.47920114117539, "5c ordered" -> True, "Qθ from local" -> True}
```
## 2026-07-14 20:25:00 — ▶️ RUN CELL 2
**In [2]:** `(*Numeric parameters (standalone solver — this file defines everything it needs).↵  Assigns global numeric values to θ1,θ2,h,κ1,κ2.*)↵θ1 = 1/3;↵θ2 = 1/7;↵h = 1;↵κ1 = N[Exp[I Zeta[3]], 30];↵κ2 = 1/κ1;↵`
**Out:** (no output)

## 2026-07-14 20:25:27 — ✏️ EDIT CELL 2
(*Numeric parameters (standalone solver — this file defines everything it needs).
  Assigns global numeric values to θ1,θ2,h,κ1,κ2.*)
θ1 = 1/3;
θ2 = 1/7;
h = 1;
κ1 = N[Exp[I Zeta[3]], 30];
κ2 = 1/κ1;
…


## 2026-07-14 20:25:27
**Input:**
```wolfram
(*Numeric parameters (standalone solver — this file defines everything it needs).
  Assigns global numeric values to θ1,θ2,h,κ1,κ2.*)
θ1 = 1/3;
θ2 = 1/7;
h = 1;
κ1 = N[Exp[I Zeta[3]], 30];
κ2 = 1/κ1;
{θ1, θ2, h, κ1, κ2}
```
**Output:**
```
Out= {1/3, 1/7, 1, 0.3604398751472265372666776505955500912937645232622850024938`29.707347828414424 + 0.9327824485933747164848359820708679294578496404879183613065`30.120295363558725*I, 0.360439875147226537266677650595550091293764523262285002494`29.421282762500656 - 0.9327824485933747164848359820708679294578496404879183613068`29.58721122792411*I}
```
## 2026-07-14 20:27:30 — ✏️ EDIT CELL 3 [batch]
(*a(u), d(u)=Qθ(u), and τ0 = τ minus the J term (using ee1=λ1+λ2-M, ee2=M).
  Qθ (=d) is defined LOCALLY here so this file is self-contained; in su2_XXX_L2.wb it came from
  the spin-chain machinery (…


## 2026-07-14 20:27:30
**Input:**
```wolfram
(*a(u), d(u)=Qθ(u), and τ0 = τ minus the J term (using ee1=λ1+λ2-M, ee2=M).
  Qθ (=d) is defined LOCALLY here so this file is self-contained; in su2_XXX_L2.wb it came from
  the spin-chain machinery (Cell 19).*)
ClearAll[aFun, τ0, Qθ];
aFun[λ1_, λ2_][u_] := (u - θ1 - h λ1) (u - θ2 - h λ2);
Qθ[u_] := (u - θ1) (u - θ2);
τ0[λ1_, λ2_, M_][u_] := (κ1 + κ2) u^2 - (θ1 + θ2) (κ1 + κ2) u - u h (κ1 (λ1 + λ2 - M) + κ2 M);
{aFun[1, 1][u], Qθ[u], τ0[1, 1, 1][u]}
```
**Output:**
```
Out= (no output)
```
## 2026-07-14 20:27:30 — ✏️ EDIT CELL 4 [batch]
(*BaxterMatrix: build O[u^k]=κ1 a(u)(u+h)^k + κ2 Qθ[u](u-h)^k - τ0(u) u^k for k=0..M.
  HARD ASSERT that the degree M+1 and M+2 coefficients vanish (this is where a wrong a/d/τ or a
  wrong ee1=λ1+λ2-…


## 2026-07-14 20:27:30
**Input:**
```wolfram
(*BaxterMatrix: build O[u^k]=κ1 a(u)(u+h)^k + κ2 Qθ[u](u-h)^k - τ0(u) u^k for k=0..M.
  HARD ASSERT that the degree M+1 and M+2 coefficients vanish (this is where a wrong a/d/τ or a
  wrong ee1=λ1+λ2-M convention would show up). Then columns = degree-0..M coefficients.*)
ClearAll[BaxterMatrix];
BaxterMatrix[λ1_, λ2_, M_] := Module[{clists, high, maxhigh, cols},
  clists = Table[
     PadRight[
      CoefficientList[
       Expand[κ1 aFun[λ1, λ2][u] (u + h)^k + κ2 Qθ[u] (u - h)^k - τ0[λ1, λ2, M][u] u^k],
       u],
      M + 3],
     {k, 0, M}];
  high = clists[[All, {M + 2, M + 3}]];   (* degree M+1 and M+2 coeffs for each column k *)
  maxhigh = Max[Abs[Flatten[high]]];
  If[maxhigh > 10^-18,
    Print["BaxterMatrix ASSERTION FAILED (λ1,λ2,M)=", {λ1, λ2, M},
          "  max|high-degree coeff|=", maxhigh, "  high=", high];
    Abort[]];
  cols = clists[[All, 1 ;; M + 1]];       (* degrees 0..M, one column per k *)
  Transpose[cols]
];
Dimensions /@ (BaxterMatrix[1, 1, #] & /@ {0, 1, 2})
```
**Output:**
```
Out= (no output)
```
## 2026-07-14 20:27:30 — ✏️ EDIT CELL 5 [batch]
(*Explicit display of the high-degree (u^(M+1), u^(M+2)) coefficients for λ1=λ2=1, M=0,1,2 —
  confirms the cancellation the assertion relies on (should all be ~10^-30).*)
Table[
 Module[{clists},
  c…


## 2026-07-14 20:27:30
**Input:**
```wolfram
(*Explicit display of the high-degree (u^(M+1), u^(M+2)) coefficients for λ1=λ2=1, M=0,1,2 —
  confirms the cancellation the assertion relies on (should all be ~10^-30).*)
Table[
 Module[{clists},
  clists = Table[
     PadRight[CoefficientList[Expand[κ1 aFun[1, 1][u] (u + h)^k + κ2 Qθ[u] (u - h)^k - τ0[1, 1, M][u] u^k], u], M + 3],
     {k, 0, M}];
  {M, "max|deg M+1,M+2 coeff|" -> Max[Abs[Flatten[clists[[All, {M + 2, M + 3}]]]]]}],
 {M, 0, 2}] // Column
```
**Output:**
```
Out= (no output)
```
## 2026-07-14 20:27:30 — ✏️ EDIT CELL 6 [batch]
(*BaxterSolve: diagonalize BaxterMatrix. Eigenvalues are the admissible J_n; eigenvectors are
  the ascending coefficient vectors of q_n, renormalized so the u^M (top) component is 1 (monic).
  FIXED …


## 2026-07-14 20:27:31
**Input:**
```wolfram
(*BaxterSolve: diagonalize BaxterMatrix. Eigenvalues are the admissible J_n; eigenvectors are
  the ascending coefficient vectors of q_n, renormalized so the u^M (top) component is 1 (monic).
  FIXED n-ordering: the M+1 solutions are sorted by ascending Re[J] (primary key) then ascending
  Im[J] (tie-break); n=0..M is assigned in that sorted order — deterministic and reproducible.
  MEMOIZED per (λ1,λ2,M) so that Q1 and τ always read the identical cached eigensystem.
  Returns {{J_0,q_0},...,{J_M,q_M}} with q_n = {c0,...,cM} (cM=1), n=0..M.*)
ClearAll[BaxterSolve, qpoly];
(*build polynomial from ascending coeff list; form it in a dummy symbol first so that
  substituting uu=0 gives the constant term rather than an indeterminate 0^0.*)
qpoly[clist_, uu_] := Module[{x}, (clist . x^Range[0, Length[clist] - 1]) /. x -> uu];
BaxterSolve[λ1_, λ2_, M_] := BaxterSolve[λ1, λ2, M] = Module[{mat, vals, vecs, pairs},
  mat = BaxterMatrix[λ1, λ2, M];
  {vals, vecs} = Eigensystem[mat];
  pairs = Table[{vals[[n]], vecs[[n]]/vecs[[n, M + 1]]}, {n, 1, M + 1}];
  (*sort into the fixed (Re[J] asc, then Im[J] asc) convention*)
  SortBy[pairs, {Re[First[#]], Im[First[#]]} &]
];
(*quick look at the M=1 sector*)
BaxterSolve[1, 1, 1]
```
**Output:**
```
Out= (no output)
```
## 2026-07-14 20:27:31 — ✏️ EDIT CELL 7 [batch]
(*INDEPENDENT residual check: re-evaluate the FULL Baxter equation from scratch (using τ=τ0+J_n,
  NOT any BaxterMatrix internals) at several u that are not roots of a(u) or d(u).
  a(u) roots: 4/3, 8…


## 2026-07-14 20:27:31
**Input:**
```wolfram
(*INDEPENDENT residual check: re-evaluate the FULL Baxter equation from scratch (using τ=τ0+J_n,
  NOT any BaxterMatrix internals) at several u that are not roots of a(u) or d(u).
  a(u) roots: 4/3, 8/7 ; d(u) roots: 1/3, 1/7 ; sample u in {0,1,3/2,-2,10} avoid all of these.
  Exact rationals are used for u so only κ1 carries finite (30-digit) precision.*)
sampleU = {0, 1, 3/2, -2, 10};
resTable = Flatten[
   Table[
    Module[{sols = BaxterSolve[1, 1, M]},
     Table[
      With[{Jn = sols[[n, 1]], cl = sols[[n, 2]]},
       {M, n - 1,
        Max[Abs[Table[
            κ1 aFun[1, 1][uu] qpoly[cl, uu + h]
            - (τ0[1, 1, M][uu] + Jn) qpoly[cl, uu]
            + κ2 Qθ[uu] qpoly[cl, uu - h],
            {uu, sampleU}]]]}],
      {n, 1, M + 1}]],
    {M, 0, 2}],
   1];
maxResidual = Max[resTable[[All, 3]]];
Print["max residual over all (M,n) and sample u = ", maxResidual];
Print[If[maxResidual < 10^-18, "RESIDUAL CHECK PASSED", "RESIDUAL CHECK FAILED"]];
TableForm[resTable, TableHeadings -> {None, {"M", "n", "max|residual|"}}]
```
**Output:**
```
Out= (no output)
```
## 2026-07-14 20:27:31 — ✏️ EDIT CELL 8 [batch]
(*Test run: full (M,n,J_n,q_n(u)) table for λ1=λ2=1 (so λ1+λ2=2, M=0,1,2). Values rounded to 10
  digits for display; full 30-digit values are available from BaxterSolve[1,1,M].*)
resultTable = Flatte…


## 2026-07-14 20:27:31
**Input:**
```wolfram
(*Test run: full (M,n,J_n,q_n(u)) table for λ1=λ2=1 (so λ1+λ2=2, M=0,1,2). Values rounded to 10
  digits for display; full 30-digit values are available from BaxterSolve[1,1,M].*)
resultTable = Flatten[
   Table[
    Module[{sols = BaxterSolve[1, 1, M]},
     Table[{M, n - 1, N[sols[[n, 1]], 10], N[qpoly[sols[[n, 2]], u], 10]}, {n, 1, M + 1}]],
    {M, 0, 2}],
   1];
TableForm[resultTable, TableHeadings -> {None, {"M", "n", "J_n", "q_n(u)"}}]
```
**Output:**
```
Out= (no output)
```
## 2026-07-14 20:27:31 — ✏️ EDIT CELL 9 [batch]
(*On-demand wrappers over the memoized BaxterSolve cache.
  Q1[λ1,λ2][M,n][u] = κ1^(u/h) q_n(u)      — the actual Baxter Q-function value at u.
  τ[λ1,λ2][M,n][u]  = τ0[λ1,λ2,M][u] + J_n — the FULL τ …


## 2026-07-14 20:27:31
**Input:**
```wolfram
(*On-demand wrappers over the memoized BaxterSolve cache.
  Q1[λ1,λ2][M,n][u] = κ1^(u/h) q_n(u)      — the actual Baxter Q-function value at u.
  τ[λ1,λ2][M,n][u]  = τ0[λ1,λ2,M][u] + J_n — the FULL τ eigenvalue for that (M,n) solution.
  Both pull from the SAME cached BaxterSolve[λ1,λ2,M] (n is 0-indexed → list entry n+1), so they
  are guaranteed mutually consistent — neither recomputes independently.
  Bounds: M an integer in 0..λ1+λ2, n an integer in 0..M; violation → clear message + $Failed.*)
ClearAll[Q1, τ];
Q1::bnds = "Invalid indices (M,n)=(`1`,`2`) for (λ1,λ2)=(`3`,`4`): need M an integer in 0..λ1+λ2 and n an integer in 0..M.";
τ::bnds = "Invalid indices (M,n)=(`1`,`2`) for (λ1,λ2)=(`3`,`4`): need M an integer in 0..λ1+λ2 and n an integer in 0..M.";
Q1[λ1_, λ2_][M_, n_][u_] := If[
   IntegerQ[M] && IntegerQ[n] && 0 <= M <= λ1 + λ2 && 0 <= n <= M,
   κ1^(u/h) qpoly[BaxterSolve[λ1, λ2, M][[n + 1, 2]], u],
   Message[Q1::bnds, M, n, λ1, λ2]; $Failed];
τ[λ1_, λ2_][M_, n_][u_] := If[
   IntegerQ[M] && IntegerQ[n] && 0 <= M <= λ1 + λ2 && 0 <= n <= M,
   τ0[λ1, λ2, M][u] + BaxterSolve[λ1, λ2, M][[n + 1, 1]],
   Message[τ::bnds, M, n, λ1, λ2]; $Failed];
(*smoke test: two valid values, plus a deliberately out-of-range call to show the guard fires*)
Print["bounds-guard test (expect a Q1::bnds message and $Failed): ", Q1[1, 1][3, 0][2]];
{Q1[1, 1][1, 0][2], τ[1, 1][1, 0][2]}
```
**Output:**
```
Out= (no output)
```
## 2026-07-14 20:27:31 — ✏️ EDIT CELL 10 [batch]
(*VALIDATION 5a: end-to-end check of the ORIGINAL (pre-ansatz-substitution) Baxter equation
    a(u) Q(u+h) - τ(u) Q(u) + d(u) Q(u-h) = 0
  assembled ENTIRELY from the wrapper functions Q1 and τ (NOT …


## 2026-07-14 20:27:31
**Input:**
```wolfram
(*VALIDATION 5a: end-to-end check of the ORIGINAL (pre-ansatz-substitution) Baxter equation
    a(u) Q(u+h) - τ(u) Q(u) + d(u) Q(u-h) = 0
  assembled ENTIRELY from the wrapper functions Q1 and τ (NOT the raw q_n/J_n), for λ1=λ2=1,
  every (M,n), at several sample u. Exact-rational u so only κ1 carries finite (30-digit)
  precision. Exercises the wrappers — memoization, sorting, and the κ1^(u/h) prefactor — top to
  bottom, independently of the internal matrix-building logic.*)
sampleU5a = {0, 1, 3/2, -2, 10};
res5a = Flatten[Table[
    Table[{M, n,
      Max[Abs[Table[
         aFun[1, 1][uu] Q1[1, 1][M, n][uu + h]
         - τ[1, 1][M, n][uu] Q1[1, 1][M, n][uu]
         + Qθ[uu] Q1[1, 1][M, n][uu - h],
         {uu, sampleU5a}]]]},
     {n, 0, M}],
    {M, 0, 2}], 1];
max5a = Max[res5a[[All, 3]]];
Print["5a max residual (via Q1/τ wrappers, λ1=λ2=1) = ", max5a];
Print[If[max5a < 10^-18, "5a PASSED", "5a FAILED"]];
TableForm[res5a, TableHeadings -> {None, {"M", "n", "max|residual|"}}]
```
**Output:**
```
Out= (no output)
```
## 2026-07-14 20:27:31 — ✏️ EDIT CELL 11 [batch]
(*VALIDATION 5b: genuine generalization to a NEW representation λ1=2, λ2=1 (asymmetric; λ1+λ2=3,
  so M=0..3). Same end-to-end wrapper residual check, plus the full (M,n,J_n) table for this case.*)
sa…


## 2026-07-14 20:27:32
**Input:**
```wolfram
(*VALIDATION 5b: genuine generalization to a NEW representation λ1=2, λ2=1 (asymmetric; λ1+λ2=3,
  so M=0..3). Same end-to-end wrapper residual check, plus the full (M,n,J_n) table for this case.*)
sampleU5b = {0, 1, 3/2, -2, 10};
res5b = Flatten[Table[
    Table[{M, n,
      Max[Abs[Table[
         aFun[2, 1][uu] Q1[2, 1][M, n][uu + h]
         - τ[2, 1][M, n][uu] Q1[2, 1][M, n][uu]
         + Qθ[uu] Q1[2, 1][M, n][uu - h],
         {uu, sampleU5b}]]]},
     {n, 0, M}],
    {M, 0, 3}], 1];
max5b = Max[res5b[[All, 3]]];
Print["5b max residual (via Q1/τ wrappers, λ1=2,λ2=1) = ", max5b];
Print[If[max5b < 10^-18, "5b PASSED", "5b FAILED"]];
JnTable5b = Flatten[Table[
    Table[{M, n, N[BaxterSolve[2, 1, M][[n + 1, 1]], 10]}, {n, 0, M}],
    {M, 0, 3}], 1];
Column[{
  TableForm[res5b, TableHeadings -> {None, {"M", "n", "max|residual|"}}],
  "(M,n,J_n) table for λ1=2, λ2=1:",
  TableForm[JnTable5b, TableHeadings -> {None, {"M", "n", "J_n"}}]}]
```
**Output:**
```
Out= (no output)
```
## 2026-07-14 20:27:32 — ✏️ EDIT CELL 12 [batch]
(*VALIDATION 5c: confirm the fixed n-ordering actually holds — Re[J_n] must be NON-DECREASING in n
  within each sector. Spot-checked on sectors with M>=2 (where >1 solution exists to order):
  (λ1,λ2…


## 2026-07-14 20:27:32
**Input:**
```wolfram
(*VALIDATION 5c: confirm the fixed n-ordering actually holds — Re[J_n] must be NON-DECREASING in n
  within each sector. Spot-checked on sectors with M>=2 (where >1 solution exists to order):
  (λ1,λ2,M) = (1,1,2), (2,1,2), (2,1,3).*)
order5c = Table[
   With[{sols = BaxterSolve[lp[[1]], lp[[2]], lp[[3]]]},
    {lp, N[sols[[All, 1]], 6], OrderedQ[Re[sols[[All, 1]]]]}],
   {lp, {{1, 1, 2}, {2, 1, 2}, {2, 1, 3}}}];
Print["5c: Re[J_n] non-decreasing in n for every checked sector? ", And @@ (Last /@ order5c)];
Print[If[And @@ (Last /@ order5c), "5c PASSED", "5c FAILED"]];
TableForm[order5c, TableHeadings -> {None, {"(λ1,λ2,M)", "J_n  (n=0..M)", "Re non-decreasing?"}}]
```
**Output:**
```
Out= (no output)
```
## 2026-07-14 20:28:22 — 🗑️ DELETE CODE CELL 29
(*Explicit display of the high-degree (u^(M+1), u^(M+2)) coefficients for λ1=λ2=1, M=0,1,2 —
  confirms the cancellation the assertion relies on (should all be ~10^-30).*)
Table[
 Module[{clists},
  c…

## 2026-07-14 20:28:22 — 🗑️ DELETE CODE CELL 31
(*INDEPENDENT residual check: re-evaluate the FULL Baxter equation from scratch (using τ=τ0+J_n,
  NOT any BaxterMatrix internals) at several u that are not roots of a(u) or d(u).
  a(u) roots: 4/3, 8…

## 2026-07-14 20:28:22 — 🗑️ DELETE CODE CELL 32
(*Test run: full (M,n,J_n,q_n(u)) table for λ1=λ2=1 (so λ1+λ2=2, M=0,1,2). Values rounded to 10
  digits for display; full 30-digit values are available from BaxterSolve[1,1,M].*)
resultTable = Flatte…

## 2026-07-14 20:28:22 — 🗑️ DELETE CODE CELL 34
(*VALIDATION 5a: end-to-end check of the ORIGINAL (pre-ansatz-substitution) Baxter equation
    a(u) Q(u+h) - τ(u) Q(u) + d(u) Q(u-h) = 0
  assembled ENTIRELY from the wrapper functions Q1 and τ (NOT …

## 2026-07-14 20:28:22 — 🗑️ DELETE CODE CELL 35
(*VALIDATION 5b: genuine generalization to a NEW, previously-untested representation λ1=2, λ2=1
  (asymmetric; λ1+λ2=3, so M=0..3). Same end-to-end wrapper residual check, plus the full
  (M,n,J_n) ta…

## 2026-07-14 20:28:22 — 🗑️ DELETE CODE CELL 36
(*VALIDATION 5c: confirm the fixed n-ordering actually holds — Re[J_n] must be NON-DECREASING in n
  within each sector. Spot-checked on sectors with M>=2 (where >1 solution exists to order):
  (λ1,λ2…


---

## 2026-07-14 20:28:32 — 🔄 KERNEL RESTART

## 2026-07-14 20:28:52 — ▶️ RUN CELLS 1–30
Cell 2: ✓ — (no output)
Cell 3: ✓ — (no output)
Cell 4: ✓ — (no output)
Cell 5: ✓ — (no output)
Cell 7: ✓ — (no output)
Cell 9: ✓ — (no output)
Cell 10: ✓ — (no output)
Cell 11: ✓ — (no output)
Cell 12: ✓ — (no output)
Cell 13: ✓ — Out[22]= \left\{0\right\}
Cell 14: ✓ — (no output)
Cell 15: ✓ — Out[25]= \begin{pmatrix}0 & 0 & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\end{pmatrix}
Cell 16: ✓ — (no output)
Cell 18: ✓ — (no output)
Cell 19: ✓ — (no output)
Cell 21: ✓ — (no output)
Cell 23: ✓ — (no output)
Cell 24: ✓ — Out[41]= \left\{0\right\}
Cell 26: ✓ — Out[47]= \left\{\frac{1}{3},\frac{1}{7},1,0.360439875147226537266677650596+0.932782448593374716484835982071\,\mathrm{i},0.36043987514722653726667765060-0.932782448593374716484835982071\,\mathrm{i}\right\}
Cell 27: ✓ — Out[51]= \left\{\left(u-\frac{4}{3}\right)\,\left(u-\frac{8}{7}\right),\left(u-\frac{1}{3}\right)\,\left(u-\frac{1}{7}\right),\left(0.720879750294453074533355301191+0.\times 10^{-30}\,\mathrm{i}\right)\,u^2-\left(1.06415582186324025288257211128+0.\times 10^{-29}\,\mathrm{i}\right)\,u\right\}
Cell 28: ✓ — Out[54]= \begin{pmatrix}1 & 1\\2 & 2\\3 & 3\end{pmatrix}
Cell 29: ✓ — Out[58]= \begin{pmatrix}-0.77812396645115795553762697294-6.5302155526567322127504666006\times 10^{-40}\,\mathrm{i} & \left\{-0.73809523809523809523809523810-0.33429538430191327287700274856\,\mathrm{i},1.00000000000000000000000000000\right\}\\1.1900552523337025695566871451+6.5302155526567322175297423340\times 10^{-40}\,\mathrm{i} & \left\{-0.73809523809523809523809523810+0.7207090391586976753748550… [472 chars total]
Cell 30: ✓ — 2 kernel messages: Q1::bnds: Invalid indices (M,n)=(3,0) for (λ1,λ2)=(1,1): need M an integer in 0..λ1+λ2 and n an integer in 0..M. [also: Q1::bnds] result: Print: bounds-guard test (expect a Q1::bnds message and $Failed):[formula] | Out[65]= \left\{-0.7092310109580169265372251422+1.0959691625127720670675272805\,\mathrm{i},-0.0229166089998261631693499907+…
⚠️ stopped at Cell 30 — error detected (pass stopOnError:false to continue past errors)


## 2026-07-14 20:29:19
**Input:**
```wolfram
{"su2 Qθ from Cell19" -> (Qθ[5] === (5-θ1)(5-θ2)),
 "su2 5a-style res" -> Max[Flatten[Table[Table[Abs[aFun[1,1][uu] Q1[1,1][M,n][uu+h]-τ[1,1][M,n][uu] Q1[1,1][M,n][uu]+Qθ[uu] Q1[1,1][M,n][uu-h]],{uu,{0,1,3/2,-2,10}}],{M,0,2},{n,0,M}]]],
 "su2 2,1 res" -> Max[Flatten[Table[Table[Abs[aFun[2,1][uu] Q1[2,1][M,n][uu+h]-τ[2,1][M,n][uu] Q1[2,1][M,n][uu]+Qθ[uu] Q1[2,1][M,n][uu-h]],{uu,{0,1,3/2,-2,10}}],{M,0,3},{n,0,M}]]],
 "cellcount defs" -> {ValueQ[θ1], ValueQ[Qθ[u]], MatchQ[BaxterSolve[2,1,3],{Repeated[{_,_}]}]}}
```
**Output:**
```
Out= {"su2 Qθ from Cell19" -> True, "su2 5a-style res" -> 0``23.61899743977343, "su2 2,1 res" -> 0``22.47920114117539, "cellcount defs" -> {True, True, True}}
```
## 2026-07-14 20:30:59 — ▶️ RUN CELLS 1–12
Cell 2: ✓ — (no output)
Cell 3: ✓ — (no output)
Cell 4: ✓ — (no output)
Cell 5: ✓ — (no output)
Cell 6: ✓ — (no output)
Cell 7: ✓ — (no output)
Cell 8: ✓ — (no output)
Cell 9: ✓ — (no output)
Cell 10: ✓ — (no output)
Cell 11: ✓ — (no output)
Cell 12: ✓ — (no output)


---

## 2026-07-15 07:34:54 — 🔄 KERNEL RESTART

