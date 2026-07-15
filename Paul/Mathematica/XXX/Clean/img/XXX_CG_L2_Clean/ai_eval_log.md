# AI Action Log — XXX_CG_L2_Clean

## 2026-07-15 10:21:28 — ✏️ EDIT CELL 1
(*Numeric parameters (shared across the whole notebook). NOTE: this assigns global numeric values to θ1,θ2,h,κ1,κ2.
  χ1,χ2 are the companion twist G's eigenvalues (Section A's G={{χ1,-χ2},{1,0}} has …


## 2026-07-15 10:21:28
**Input:**
```wolfram
(*Numeric parameters (shared across the whole notebook). NOTE: this assigns global numeric values to θ1,θ2,h,κ1,κ2.
  χ1,χ2 are the companion twist G's eigenvalues (Section A's G={{χ1,-χ2},{1,0}} has char. poly
  λ²-χ1 λ+χ2=(λ-κ1)(λ-κ2), so χ1=κ1+κ2, χ2=κ1κ2=1 here since κ2=1/κ1) — fixed here, numerically,
  once for the whole merged notebook (both Section A and Section B rely on G/t2 being numeric).*)
θ1 = 1/3;
θ2 = 1/7;
h = 1;
κ1 = N[Exp[I Zeta[3]], 30];
κ2 = 1/κ1;
χ1 = κ1 + κ2;
χ2 = 1;
{θ1, θ2, h, κ1, κ2, χ1, χ2}
```
**Output:**
```
Out= {1/3, 1/7, 1, 0.3604398751472265372666776505955500912937645232622850024938`29.707347828414424 + 0.9327824485933747164848359820708679294578496404879183613065`30.120295363558725*I, 0.360439875147226537266677650595550091293764523262285002494`29.421282762500656 - 0.9327824485933747164848359820708679294578496404879183613068`29.58721122792411*I, 0.7208797502944530745333553011911001825875290465245700049878`29.54117567596917 + 0``29.50582140452307*I, 1}
```
## 2026-07-15 10:22:04 — 📥 BULK INSERT 25 CELL(S) at positions 2–26
1. [markdown] (*Two-site XXX_{su(2)} chain (L=2). Sites carry independent gl(2) representations [λ1,0], [λ2,0], in
2. [code] ps[rep_] := ↵  ps[rep] = ↵   Block[{res = (Prepend[rep, #] & /@ (Tuples[↵          Table[Range[rep[[
3. [code] (*spin-s generators*)↵ClearAll[Ee];↵Ee[λ_][i_,j_]:=JJ[i,j,{λ,0}]
4. [code] (*Lax matrix for spin-s representation (no θ baked in — shift applied when a site is instantiated)*)
5. [code] (*companion twist matrix — shared by the whole chain*)↵G={{χ1,-χ2},{1,0}};
6. [markdown] (*Two-site Lax operators on V_λ1⊗V_λ2, built as KroneckerProduct of the single-site Lax matrix (eval
7. [code] ClearAll[L1,L2];↵L1[λ1_,λ2_][i_,j_][u_]:=KroneckerProduct[L[λ1][i,j][u-θ1],IdentityMatrix[dim[{λ2,0}
8. [markdown] (*Monodromy matrix T[λ1,λ2][i,j][u] = product of the two site Lax operators in auxiliary space, T(u)
9. [code] ClearAll[T];↵T[λ1_,λ2_][i_,j_][u_]:=T[λ1,λ2][i,j][u]=Sum[L1[λ1,λ2][i,k][u].L2[λ1,λ2][k,j][u],{k,2}];
10. [code] (*transfer matrix in fundamental representation*)↵ClearAll[t2];↵t2[λ1_,λ2_][1,1][u_]:=Sum[T[λ1,λ2][i
11. [code] (*transfer matrices commute at different values of spectral parameter*)↵(*so coefficients generate a
12. [code] (*transfer matrix in fundamental representation with diagonal twist*)↵t2diag[λ1_,λ2_][1,1][u_]:=Sum[
13. [code] (*global Lie algebra generators*)↵ClearAll[EE]↵EE[λ1_,λ2_][i_,j_]:=KroneckerProduct[Ee[λ1][i,j],Iden
14. [code] (*central charge*)↵Table[EE[λ1,λ2][1,1]+EE[λ1,λ2][2,2]-(λ1+λ2)IdentityMatrix[dim[{λ1,0}]dim[{λ2,0}]]
15. [code] (*asymptotics of transfer matrix in diagonal frame*)↵tdiagasympt[λ1_,λ2_][u_]:=(κ[1]+κ[2])u^2 Identi
16. [code] (*checking asymptotics matches*)↵t2diag[1,1][1,1][u]-tdiagasympt[1,1][u];↵Coefficient[%,u,1]//Simpli
17. [code] (*quantum determinant*)↵ClearAll[qdetT];↵qdetT[λ1_,λ2_][u_]:=T[λ1,λ2][1,1][u].T[λ1,λ2][2,2][u-h]-T[λ
18. [markdown] (*Fusion/Hirota hierarchy for the two-site transfer matrices t2[λ1,λ2][1,s][u], built exactly as in 
19. [code] (*fused transfer matrices*)↵t2[λ1_,λ2_][0,s_][u_]:=IdentityMatrix[dim[{λ1,0}]*dim[{λ2,0}]];↵t2[λ1_,λ
20. [code] (*Qθ-type factors for the two sites*)↵ClearAll[Qθ1,Qθ];↵Qθ1[u_]:=(u-θ1);↵Qθ[u_]:=(u-θ1)(u-θ2);
21. [markdown] (*SoV B-operator for L=2: the bare (1,1) monodromy element (no G, no trace) — same convention as B[λ
22. [code] ClearAll[B2];↵B2[λ1_,λ2_][u_]:=T[λ1,λ2][1,1][u];
23. [markdown] (*Step 8: SoV pseudovacuum for L=2 as the tensor product of the two single-site lowest-weight covect
24. [code] ClearAll[xSingle,x2];↵xSingle[λ_][{0}]:={Table[KroneckerDelta[ii,dim[{λ,0}]],{ii,1,dim[{λ,0}]}]};↵x2
25. [code] (*checking x2 diagonalizes B*)↵Table[Table[x2[λ1,λ2][{s1,s2}].B2[λ1,λ2][u]-(u-θ1-h s1)(u-θ2-h s2)x2[


## 2026-07-15 10:22:04
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

## 2026-07-15 10:22:05
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

## 2026-07-15 10:22:05
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

## 2026-07-15 10:22:05
**Input:**
```wolfram
(*companion twist matrix — shared by the whole chain*)
G={{χ1,-χ2},{1,0}};
```
**Output:**
```
(no output)
```

## 2026-07-15 10:22:05
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

## 2026-07-15 10:22:06
**Input:**
```wolfram
ClearAll[T];
T[λ1_,λ2_][i_,j_][u_]:=T[λ1,λ2][i,j][u]=Sum[L1[λ1,λ2][i,k][u].L2[λ1,λ2][k,j][u],{k,2}];
```
**Output:**
```
(no output)
```

## 2026-07-15 10:22:06
**Input:**
```wolfram
(*transfer matrix in fundamental representation*)
ClearAll[t2];
t2[λ1_,λ2_][1,1][u_]:=Sum[T[λ1,λ2][i,j][u]G[[j,i]],{i,2},{j,2}];
```
**Output:**
```
(no output)
```

## 2026-07-15 10:22:07
**Input:**
```wolfram
(*transfer matrices commute at different values of spectral parameter*)
(*so coefficients generate a complete family of integrals of motion*)
Table[t2[λ1,λ2][1,1][u].t2[λ1,λ2][1,1][v]-t2[λ1,λ2][1,1][v].t2[λ1,λ2][1,1][u],{λ1,1,4},{λ2,1,4}]//Simplify//Flatten//Union
```
**Output:**
```
(no output)
```

## 2026-07-15 10:22:07
**Input:**
```wolfram
(*transfer matrix in fundamental representation with diagonal twist*)
t2diag[λ1_,λ2_][1,1][u_]:=Sum[T[λ1,λ2][i,i][u]κ[i],{i,2}];
```
**Output:**
```
(no output)
```

## 2026-07-15 10:22:07
**Input:**
```wolfram
(*global Lie algebra generators*)
ClearAll[EE]
EE[λ1_,λ2_][i_,j_]:=KroneckerProduct[Ee[λ1][i,j],IdentityMatrix[dim[{λ2,0}]]]+KroneckerProduct[IdentityMatrix[dim[{λ1,0}]],Ee[λ2][i,j]];
```
**Output:**
```
(no output)
```

## 2026-07-15 10:22:07
**Input:**
```wolfram
(*central charge*)
Table[EE[λ1,λ2][1,1]+EE[λ1,λ2][2,2]-(λ1+λ2)IdentityMatrix[dim[{λ1,0}]dim[{λ2,0}]],{λ1,1,4},{λ2,1,4}]//Flatten//Union
```
**Output:**
```
(no output)
```

## 2026-07-15 10:22:08
**Input:**
```wolfram
(*asymptotics of transfer matrix in diagonal frame*)
tdiagasympt[λ1_,λ2_][u_]:=(κ[1]+κ[2])u^2 IdentityMatrix[dim[{λ1,0}]*dim[{λ2,0}]]-u(κ[1]+κ[2])(θ1+θ2) IdentityMatrix[dim[{λ1,0}]*dim[{λ2,0}]]- u h(κ[1] EE[λ1,λ2][1,1]+κ[2] EE[λ1,λ2][2,2])
```
**Output:**
```
(no output)
```

## 2026-07-15 10:22:08
**Input:**
```wolfram
(*checking asymptotics matches*)
t2diag[1,1][1,1][u]-tdiagasympt[1,1][u];
Coefficient[%,u,1]//Simplify
```
**Output:**
```
(no output)
```

## 2026-07-15 10:22:08
**Input:**
```wolfram
(*quantum determinant*)
ClearAll[qdetT];
qdetT[λ1_,λ2_][u_]:=T[λ1,λ2][1,1][u].T[λ1,λ2][2,2][u-h]-T[λ1,λ2][2,1][u].T[λ1,λ2][1,2][u-h];
```
**Output:**
```
(no output)
```

## 2026-07-15 10:22:08
**Input:**
```wolfram
(*fused transfer matrices*)
t2[λ1_,λ2_][0,s_][u_]:=IdentityMatrix[dim[{λ1,0}]*dim[{λ2,0}]];
t2[λ1_,λ2_][a_,0][u_]:=IdentityMatrix[dim[{λ1,0}]*dim[{λ2,0}]];
t2[λ1_,λ2_][2,1][u_]:=χ2 qdetT[λ1,λ2][u];
t2[λ1_,λ2_][1,s_/;s>1][u_]:=t2[λ1,λ2][1,s-1][u].t2[λ1,λ2][1,1][u+(s-1)h]-t2[λ1,λ2][2,1][u+(s-1)h].t2[λ1,λ2][1,s-2][u];
```
**Output:**
```
(no output)
```

## 2026-07-15 10:22:09
**Input:**
```wolfram
(*Qθ-type factors for the two sites*)
ClearAll[Qθ1,Qθ];
Qθ1[u_]:=(u-θ1);
Qθ[u_]:=(u-θ1)(u-θ2);
```
**Output:**
```
(no output)
```

## 2026-07-15 10:22:09
**Input:**
```wolfram
ClearAll[B2];
B2[λ1_,λ2_][u_]:=T[λ1,λ2][1,1][u];
```
**Output:**
```
(no output)
```

## 2026-07-15 10:22:09
**Input:**
```wolfram
ClearAll[xSingle,x2];
xSingle[λ_][{0}]:={Table[KroneckerDelta[ii,dim[{λ,0}]],{ii,1,dim[{λ,0}]}]};
x2[λ1_,λ2_][{0,0}]:=KroneckerProduct[xSingle[λ1][{0}],xSingle[λ2][{0}]];

x2[λ1_,λ2_][{s1_,s2_}] := x2[λ1,λ2][{s1,s2}] = x2[λ1,λ2][{0,0}].t2[λ1,λ2][1,s1][θ1].t2[λ1,λ2][1,s2][θ2]/(Product[Qθ[θ1 + (-λ1 + k - 1) h], {k, 1, s1}] Product[Qθ[θ2 + (-λ2 + k - 1) h], {k, 1, s2}]);
```
**Output:**
```
(no output)
```

## 2026-07-15 10:22:10
**Input:**
```wolfram
(*checking x2 diagonalizes B*)
Table[Table[x2[λ1,λ2][{s1,s2}].B2[λ1,λ2][u]-(u-θ1-h s1)(u-θ2-h s2)x2[λ1,λ2][{s1,s2}],{s1,0,λ1},{s2,0,λ2}],{λ1,1,3},{λ2,1,3}]//Simplify//Flatten//Union
```
**Output:**
```
(no output)
```
## 2026-07-15 10:23:03 — ▶️ RUN CELL 12
**In [12]:** `(*transfer matrices commute at different values of spectral parameter*)↵(*so coefficients generate a complete family of integrals of motion*)↵Table[t2[λ1,λ2][1,1][u].t2[λ1,λ2][1,1][v]-t2[λ1,λ2][1,1][v`
**Out:** (no output)

## 2026-07-15 10:23:24 — ▶️ RUN CELL 12
**In [12]:** `(*transfer matrices commute at different values of spectral parameter*)↵(*so coefficients generate a complete family of integrals of motion*)↵Table[t2[λ1,λ2][1,1][u].t2[λ1,λ2][1,1][v]-t2[λ1,λ2][1,1][v`
**Out:** (no output)


## 2026-07-15 10:23:36
**Input:**
```wolfram
Table[t2[λ1,λ2][1,1][u].t2[λ1,λ2][1,1][v]-t2[λ1,λ2][1,1][v].t2[λ1,λ2][1,1][u],{λ1,1,4},{λ2,1,4}]//Simplify//Flatten//Union
```
**Output:**
```
Out= {0, 0``28.22254937886658 + 0``28.44383859562073*I, 0``28.14420878610075 + 0``27.96671734090107*I, 0``27.84317879043677 + 0``27.66568734523709*I, 0``27.268306869427246 + 0``27.489596086181407*I, (0``29.85940410877844 + 0``29.681912663578753*I) + (0``28.65546697109418 + 0``28.87675618784834*I)*u + (0``29.382282854058772 + 0``29.204791408859087*I)*v, (0``28.82460937019454 + 0``29.04589858694869*I)*(-7/3 + u)*(-15/7 + u)*(-7/3 + v)*(-15/7 + v), (0``28.82460937019454 + 0``29.04589858694869*I)*(-15/7 + u)*(-4/3 + u)*(-15/7 + v)*(-4/3 + v), (0``28.82460937019454 + 0``29.04589858694869*I)*(-7/3 + u)*(-8/7 + u)*(-7/3 + v)*(-8/7 + v), (0``28.82460937019454 + 0``29.04589858694869*I)*(-4/3 + u)*(-8/7 + u)*(-4/3 + v)*(-8/7 + v), (0``28.641678686608554 + 0``28.862967903362705*I) + (0``28.65546697109418 + 0``28.87675618784834*I)*u + (0``28.82460937019454 + 0``29.04589858694869*I)*u^2 + (0``28.43082532129366 + 0``28.65211453804781*I)*v + (0``28.82460937019454 + 0``29.04589858694869*I)*v^2, (-7/3 + u)*(7/3 - v) + (0``28.82460937019454 + 0``29.04589858694869*I)*(-7/3 + u)*(-15/7 + u)*(-7/3 + v)*(-15/7 + v) + ((-7 + 3*u)*(-7 + 3*v))/9, (-7/3 + u)*(7/3 - v) + (0``28.82460937019454 + 0``29.04589858694869*I)*(-7/3 + u)*(-8/7 + u)*(-7/3 + v)*(-8/7 + v) + ((-7 + 3*u)*(-7 + 3*v))/9, (7/3 - u)*(-7/3 + v) + (0``28.82460937019454 + 0``29.04589858694869*I)*(-7/3 + u)*(-8/7 + u)*(-7/3 + v)*(-8/7 + v) + ((-7 + 3*u)*(-7 + 3*v))/9, (7/3 - u)*(-7/3 + v) + (0``28.82460937019454 + 0``29.04589858694869*I)*(-7/3 + u)*(-1/7 + u)*(-7/3 + v)*(-1/7 + v) + ((-7 + 3*u)*(-7 + 3*v))/9, 2*(-7/3 + u)*(7/3 - v) + (0``28.82460937019454 + 0``29.04589858694869*I)*(-22/7 + u)*(-7/3 + u)*(-22/7 + v)*(-7/3 + v) + (2*(-7 + 3*u)*(-7 + 3*v))/9, 2*(-7/3 + u)*(7/3 - v) + (0``28.82460937019454 + 0``29.04589858694869*I)*(-7/3 + u)*(-15/7 + u)*(-7/3 + v)*(-15/7 + v) + (2*(-7 + 3*u)*(-7 + 3*v))/9, 2*(7/3 - u)*(-7/3 + v) + (0``28.82460937019454 + 0``29.04589858694869*I)*(-7/3 + u)*(-8/7 + u)*(-7/3 + v)*(-8/7 + v) + (2*(-7 + 3*u)*(-7 + 3*v))/9, 2*(7/3 - u)*(-7/3 + v) + (0``28.82460937019454 + 0``29.04589858694869*I)*(-7/3 + u)*(-1/7 + u)*(-7/3 + v)*(-1/7 + v) + (2*(-7 + 3*u)*(-7 + 3*v))/9, 3*(-7/3 + u)*(7/3 - v) + (0``28.82460937019454 + 0``29.04589858694869*I)*(-22/7 + u)*(-7/3 + u)*(-22/7 + v)*(-7/3 + v) + ((-7 + 3*u)*(-7 + 3*v))/3, 3*(7/3 - u)*(-7/3 + v) + (0``28.82460937019454 + 0``29.04589858694869*I)*(-7/3 + u)*(-1/7 + u)*(-7/3 + v)*(-1/7 + v) + ((-7 + 3*u)*(-7 + 3*v))/3, 4*(-7/3 + u)*(7/3 - v) + (0``28.82460937019454 + 0``29.04589858694869*I)*(-29/7 + u)*(-7/3 + u)*(-29/7 + v)*(-7/3 + v) + (4*(-7 + 3*u)*(-7 + 3*v))/9, 4*(7/3 - u)*(-7/3 + v) + (0``28.82460937019454 + 0``29.04589858694869*I)*(-7/3 + u)*(-1/7 + u)*(-7/3 + v)*(-1/7 + v) + (4*(-7 + 3*u)*(-7 + 3*v))/9, (-4/3 + u)*(4/3 - v) + (0``28.82460937019454 + 0``29.04589858694869*I)*(-15/7 + u)*(-4/3 + u)*(-15/7 + v)*(-4/3 + v) + ((-4 + 3*u)*(-4 + 3*v))/9, (-4/3 + u)*(4/3 - v) + (0``28.82460937019454 + 0``29.04589858694869*I)*(-4/3 + u)*(-8/7 + u)*(-4/3 + v)*(-8/7 + v) + ((-4 + 3*u)*(-4 + 3*v))/9, (4/3 - u)*(-4/3 + v) + (0``28.82460937019454 + 0``29.04589858694869*I)*(-4/3 + u)*(-8/7 + u)*(-4/3 + v)*(-8/7 + v) + ((-4 + 3*u)*(-4 + 3*v))/9, (4/3 - u)*(-4/3 + v) + (0``28.82460937019454 + 0``29.04589858694869*I)*(-4/3 + u)*(-1/7 + u)*(-4/3 + v)*(-1/7 + v) + ((-4 + 3*u)*(-4 + 3*v))/9, 2*(-4/3 + u)*(4/3 - v) + (0``28.82460937019454 + 0``29.04589858694869*I)*(-22/7 + u)*(-4/3 + u)*(-22/7 + v)*(-4/3 + v) + (2*(-4 + 3*u)*(-4 + 3*v))/9, 2*(-4/3 + u)*(4/3 - v) + (0``28.82460937019454 + 0``29.04589858694869*I)*(-15/7 + u)*(-4/3 + u)*(-15/7 + v)*(-4/3 + v) + (2*(-4 + 3*u)*(-4 + 3*v))/9, 2*(4/3 - u)*(-4/3 + v) + (0``28.82460937019454 + 0``29.04589858694869*I)*(-4/3 + u)*(-8/7 + u)*(-4/3 + v)*(-8/7 + v) + (2*(-4 + 3*u)*(-4 + 3*v))/9, 2*(4/3 - u)*(-4/3 + v) + (0``28.82460937019454 + 0``29.04589858694869*I)*(-4/3 + u)*(-1/7 + u)*(-4/3 + v)*(-1/7 + v) + (2*(-4 + 3*u)*(-4 + 3*v))/9, 3*(-4/3 + u)*(4/3 - v) + (0``28.82460937019454 + 0``29.04589858694869*I)*(-22/7 + u)*(-4/3 + u)*(-22/7 + v)*(-4/3 + v) + ((-4 + 3*u)*(-4 + 3*v))/3, 3*(4/3 - u)*(-4/3 +
[output truncated — 8592 chars total]
```

## 2026-07-15 10:23:59
**Input:**
```wolfram
Chop[Table[t2[λ1,λ2][1,1][u].t2[λ1,λ2][1,1][v]-t2[λ1,λ2][1,1][v].t2[λ1,λ2][1,1][u],{λ1,1,4},{λ2,1,4}]/.{u->1.3,v->2.7},10^-15]//Flatten//Union
```
**Output:**
```
Out= {-2.1316282072803006*^-14, -1.7763568394002505*^-14, -1.5987211554602254*^-14, -1.0658141036401503*^-14, -7.105427357601002*^-15, -6.217248937900877*^-15, -5.329070518200751*^-15, -4.884981308350689*^-15, -4.440892098500626*^-15, -3.9968028886505635*^-15, -3.552713678800501*^-15, -3.4416913763379853*^-15, -3.3306690738754696*^-15, -2.6645352591003757*^-15, -2.220446049250313*^-15, -1.7763568394002505*^-15, -1.7486012637846216*^-15, -1.5543122344752192*^-15, -1.5126788710517758*^-15, -1.4988010832439613*^-15, -1.4710455076283324*^-15, -1.3322676295501878*^-15, -1.1102230246251565*^-15, 0, 1.1102230246251565*^-15, 1.3322676295501878*^-15, 1.5543122344752192*^-15, 1.7208456881689926*^-15, 1.7763568394002505*^-15, 1.9984014443252818*^-15, 2.220446049250313*^-15, 2.4424906541753444*^-15, 2.6645352591003757*^-15, 3.4416913763379853*^-15, 3.552713678800501*^-15, 3.9968028886505635*^-15, 4.440892098500626*^-15, 5.329070518200751*^-15, 7.105427357601002*^-15, 7.993605777301127*^-15, 8.881784197001252*^-15, 1.0658141036401503*^-14}
```
## 2026-07-15 10:24:50 — ▶️ RUN CELL 15
**In [15]:** `(*central charge*)↵Table[EE[λ1,λ2][1,1]+EE[λ1,λ2][2,2]-(λ1+λ2)IdentityMatrix[dim[{λ1,0}]dim[{λ2,0}]],{λ1,1,4},{λ2,1,4}]//Flatten//Union`
**Out:** (no output)

## 2026-07-15 10:24:54 — ▶️ RUN CELL 17
**In [17]:** `(*checking asymptotics matches*)↵t2diag[1,1][1,1][u]-tdiagasympt[1,1][u];↵Coefficient[%,u,1]//Simplify`
**Out:** (no output)

## 2026-07-15 10:24:58 — ▶️ RUN CELL 26
**In [26]:** `(*checking x2 diagonalizes B*)↵Table[Table[x2[λ1,λ2][{s1,s2}].B2[λ1,λ2][u]-(u-θ1-h s1)(u-θ2-h s2)x2[λ1,λ2][{s1,s2}],{s1,0,λ1},{s2,0,λ2}],{λ1,1,3},{λ2,1,3}]//Simplify//Flatten//Union`
**Out:** (no output)


## 2026-07-15 10:25:18
**Input:**
```wolfram
Table[EE[λ1,λ2][1,1]+EE[λ1,λ2][2,2]-(λ1+λ2)IdentityMatrix[dim[{λ1,0}]dim[{λ2,0}]],{λ1,1,4},{λ2,1,4}]//Flatten//Union
```
**Output:**
```
Out= {0}
```

## 2026-07-15 10:25:21
**Input:**
```wolfram
t2diag[1,1][1,1][u]-tdiagasympt[1,1][u];Coefficient[%,u,1]//Simplify
```
**Output:**
```
Out[1]= {0, 0``28.65546697109418 + 0``28.87675618784834*I, 0``28.65546697109418 + 0``28.87675618784834*I, 0``28.65546697109418 + 0``28.87675618784834*I, 0``28.65546697109418 + 0``28.87675618784834*I, 0``28.65546697109418 + 0``28.87675618784834*I, 0``28.65546697109418 + 0``28.87675618784834*I, 0``28.65546697109418 + 0``28.87675618784834*I, 0``28.65546697109418 + 0``28.87675618784834*I, 0``28.65546697109418 + 0``28.87675618784834*I, 0``28.65546697109418 + 0``28.87675618784834*I, 0``28.65546697109418 + 0``28.87675618784834*I, 0``28.65546697109418 + 0``28.87675618784834*I, 0``28.65546697109418 + 0``28.87675618784834*I, 0``28.65546697109418 + 0``28.87675618784834*I}
```

## 2026-07-15 10:25:28
**Input:**
```wolfram
Table[Table[x2[λ1,λ2][{s1,s2}].B2[λ1,λ2][u]-(u-θ1-h s1)(u-θ2-h s2)x2[λ1,λ2][{s1,s2}],{s1,0,λ1},{s2,0,λ2}],{λ1,1,3},{λ2,1,3}]//Simplify//Flatten//Union
```
**Output:**
```
Out= {0, (0``31.16371041018973 + 0``30.986218964990048*I) + (0``28.65546697109418 + 0``28.87675618784834*I)*u + (0``28.82460937019454 + 0``29.04589858694869*I)*u^2, (0``30.98304837503374 + 0``30.805556929834047*I) + (0``28.65546697109418 + 0``28.87675618784834*I)*u + (0``28.82460937019454 + 0``29.04589858694869*I)*u^2, (0``30.63720603495339 + 0``30.459714589753702*I) + (0``28.65546697109418 + 0``28.87675618784834*I)*u + (0``28.82460937019454 + 0``29.04589858694869*I)*u^2, (0``30.562167102131593 + 0``30.384675656931908*I) + (0``28.65546697109418 + 0``28.87675618784834*I)*u + (0``28.82460937019454 + 0``29.04589858694869*I)*u^2, (0``30.543430646216596 + 0``30.365939201016907*I) + (0``28.65546697109418 + 0``28.87675618784834*I)*u + (0``28.82460937019454 + 0``29.04589858694869*I)*u^2, (0``30.22420227984548 + 0``30.04671083464579*I) + (0``28.65546697109418 + 0``28.87675618784834*I)*u + (0``28.82460937019454 + 0``29.04589858694869*I)*u^2, (0``30.15370779602599 + 0``29.97621635082631*I) + (0``28.65546697109418 + 0``28.87675618784834*I)*u + (0``28.82460937019454 + 0``29.04589858694869*I)*u^2, (0``30.038079363024714 + 0``29.860587917825033*I) + (0``28.65546697109418 + 0``28.87675618784834*I)*u + (0``28.82460937019454 + 0``29.04589858694869*I)*u^2, (0``29.913802232849132 + 0``29.736310787649444*I) + (0``28.65546697109418 + 0``28.87675618784834*I)*u + (0``28.82460937019454 + 0``29.04589858694869*I)*u^2, (0``29.612717684345665 + 0``29.435226239145983*I) + (0``28.65546697109418 + 0``28.87675618784834*I)*u + (0``28.82460937019454 + 0``29.04589858694869*I)*u^2, (0``29.55034626010249 + 0``29.372854814902798*I) + (0``28.65546697109418 + 0``28.87675618784834*I)*u + (0``28.82460937019454 + 0``29.04589858694869*I)*u^2, (0``29.5199958787035 + 0``29.34250443350382*I) + (0``28.65546697109418 + 0``28.87675618784834*I)*u + (0``28.82460937019454 + 0``29.04589858694869*I)*u^2, (0``29.394039951336364 + 0``29.216548506136682*I) + (0``28.65546697109418 + 0``28.87675618784834*I)*u + (0``28.82460937019454 + 0``29.04589858694869*I)*u^2, (0``29.37828249032985 + 0``29.200791045130167*I) + (0``28.65546697109418 + 0``28.87675618784834*I)*u + (0``28.82460937019454 + 0``29.04589858694869*I)*u^2}
```
## 2026-07-15 10:25:44 — ▶️ RUN CELL 26
**In [26]:** `(*checking x2 diagonalizes B*)↵Table[Table[x2[λ1,λ2][{s1,s2}].B2[λ1,λ2][u]-(u-θ1-h s1)(u-θ2-h s2)x2[λ1,λ2][{s1,s2}],{s1,0,λ1},{s2,0,λ2}],{λ1,1,3},{λ2,1,3}]//Simplify//Flatten//Union`
**Out:** (no output)

## 2026-07-15 10:25:52 — ▶️ RUN CELL 17
**In [17]:** `(*checking asymptotics matches*)↵t2diag[1,1][1,1][u]-tdiagasympt[1,1][u];↵Coefficient[%,u,1]//Simplify`
**Out:** (no output)

## 2026-07-15 10:26:39 — ▶️ RUN CELL 12
**In [12]:** `(*transfer matrices commute at different values of spectral parameter*)↵(*so coefficients generate a complete family of integrals of motion*)↵Table[t2[λ1,λ2][1,1][u].t2[λ1,λ2][1,1][v]-t2[λ1,λ2][1,1][v`
**Out:** (no output)

## 2026-07-15 10:27:19 — ▶️ RUN CELL 12
**In [12]:** `(*transfer matrices commute at different values of spectral parameter*)↵(*so coefficients generate a complete family of integrals of motion*)↵Table[t2[λ1,λ2][1,1][u].t2[λ1,λ2][1,1][v]-t2[λ1,λ2][1,1][v`
**Out:** (no output)

## 2026-07-15 10:28:48 — 📥 BULK INSERT 8 CELL(S) at positions 27–34
1. [markdown] (*Section B — Baxter (TQ) solver: minimal reusable routine set (BaxterMatrix, TauEigensystem/Psi/τ, 
2. [code] (*a(u), d(u)=Qθ(u), and τ0 = τ minus the J term (using ee1=λ1+λ2-M, ee2=M).↵  Qθ (=d) is defined LOC
3. [code] (*BaxterMatrix: build O[u^k]=κ1 a(u)(u+h)^k + κ2 Qθ[u](u-h)^k - τ0(u) u^k for k=0..M.↵  HARD ASSERT 
4. [code] (* Coefficient matrices C0,C1,C2 with t2[λ1,λ2][1,1][u] = C0 + C1 u + C2 u^2.↵   Extracted entrywise
5. [code] (* dPred[M,λ1,λ2]: predicted dimension of the weight-M subspace of Vλ1⊗Vλ2 (number of (m1,m2)↵   pai
6. [code] (* Psi[λ1,λ2][M,n]: on-demand accessor over the memoized TauEigensystem cache. Returns the fixed↵   
7. [code] (* QSolve[λ1,λ2][M,n]: the TQ equation is now a LINEAR solve, not an eigenvalue problem, because↵   
8. [code] (* Q1[λ1,λ2][M,n][u] = κ1^(u/h) q_n(u): the actual Baxter Q-function value, built from the SAME↵   Q


## 2026-07-15 10:28:49
**Input:**
```wolfram
(*a(u), d(u)=Qθ(u), and τ0 = τ minus the J term (using ee1=λ1+λ2-M, ee2=M).
  Qθ (=d) is defined LOCALLY here so this file is self-contained.
  qpoly builds a polynomial from an ASCENDING coefficient list {c0,...,cM} (cM=1 for monic q_n);
  formed in a dummy symbol first so substituting uu=0 gives the constant term rather than an
  indeterminate 0^0.*)
ClearAll[aFun, τ0, Qθ, qpoly];
aFun[λ1_, λ2_][u_] := (u - θ1 - h λ1) (u - θ2 - h λ2);
Qθ[u_] := (u - θ1) (u - θ2);
τ0[λ1_, λ2_, M_][u_] := (κ1 + κ2) u^2 - (θ1 + θ2) (κ1 + κ2) u - u h (κ1 (λ1 + λ2 - M) + κ2 M);
qpoly[clist_, uu_] := Module[{x}, (clist . x^Range[0, Length[clist] - 1]) /. x -> uu];
{aFun[1, 1][u], Qθ[u], τ0[1, 1, 1][u]}
```
**Output:**
```
(no output)
```

## 2026-07-15 10:28:49
**Input:**
```wolfram
(*BaxterMatrix: build O[u^k]=κ1 a(u)(u+h)^k + κ2 Qθ[u](u-h)^k - τ0(u) u^k for k=0..M.
  HARD ASSERT that the degree M+1 and M+2 coefficients vanish (this is where a wrong a/d/τ or a
  wrong ee1=λ1+λ2-M convention would show up). Then columns = degree-0..M coefficients.
  Used TWICE downstream: (old role) never diagonalized directly anymore; (new role) as the FIXED
  linear operator whose nullspace at a KNOWN genuine J gives the Q-function coefficients
  (see QSolve below) — this is now a linear solve, not an eigenvalue problem.*)
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
(no output)
```

## 2026-07-15 10:28:49
**Input:**
```wolfram
(* Coefficient matrices C0,C1,C2 with t2[λ1,λ2][1,1][u] = C0 + C1 u + C2 u^2.
   Extracted entrywise via CoefficientList in a dummy variable (NOT by evaluating at numeric
   points). Memoized per (λ1,λ2) — reusable across all (M,n) of that representation. *)
ClearAll[t2Coeffs];
t2Coeffs[λ1_, λ2_] := t2Coeffs[λ1, λ2] = Module[{uu, mat, cl},
   mat = t2[λ1, λ2][1, 1][uu];
   cl = Map[PadRight[CoefficientList[#, uu], 3] &, mat, {2}];
   Table[cl[[All, All, k]], {k, 1, 3}]
];
(* check: C2 must equal (κ1+κ2) Id — leading-order asymptotics match τ by construction *)
With[{Cs = t2Coeffs[1, 1]}, "max|C2-(κ1+κ2)Id|" -> Max[Abs[Flatten[Cs[[3]] - (κ1 + κ2) IdentityMatrix[4]]]]]
```
**Output:**
```
(no output)
```

## 2026-07-15 10:28:50
**Input:**
```wolfram
(* dPred[M,λ1,λ2]: predicted dimension of the weight-M subspace of Vλ1⊗Vλ2 (number of (m1,m2)
   pairs with m1+m2 fixed) = Min[M,λ1,λ2,λ1+λ2-M]+1. This is what the GENUINE t2-diagonalization
   below must reproduce exactly for every M — no more, no less — since it is now an honest
   d=(λ1+1)(λ2+1)-dimensional diagonalization of the actual transfer matrix, not an oversized
   (M+1)-dimensional Baxter-matrix proxy that can manufacture spurious extra roots. *)
ClearAll[dPred];
dPred[M_, λ1_, λ2_] := Min[M, λ1, λ2, λ1 + λ2 - M] + 1;

(* TauEigensystem[λ1,λ2]: diagonalize the GENUINE transfer matrix t2[λ1,λ2][1,1][u] directly.
   C0,C1,C2 pairwise commute (all built from the same commuting-family monodromy), so a GENERIC
   linear combination C0+r C1 has the full set of common eigenvectors of C0 and C1 (C2 is already
   known to be the scalar (κ1+κ2) Id, so it carries no extra information). We try a short list of
   generic r candidates and keep the first that fully resolves the d-dim eigenspace (rank d,
   d distinct combo-eigenvalues); this guards against an accidental coincidental degeneracy in
   C0+r C1 for a specific r.
   For each eigenvector v: extract t0=(v.C0.v)/(v.v), t1=(v.C1.v)/(v.v) (valid for ANY genuine
   right eigenvector, no conjugation needed since C0,C1 are not Hermitian). The residual checks
   ||C0.v-t0 v|| and ||C1.v-t1 v|| being ~0 IS the genuineness test (replaces the old
   Baxter-matrix + post-hoc Psi::spurious filtering entirely: this diagonalization is exact, so a
   real failure here means a genuine construction bug, hence a hard Abort).
   M is read off ALGEBRAICALLY from t1 via τ0's u^1 coefficient (τ=τ0+J only shifts u^0):
     M = [-(θ1+θ2)(κ1+κ2) - hκ1(λ1+λ2) - t1] / (h(κ2-κ1))
   and must land within 10^-10 of an integer (hard assert). n within each M-sector: J=t0 directly
   (τ0 has zero constant term), sorted ascending by (Re[J],Im[J]) — same tie-break convention as
   the old BaxterSolve. Each M-sector's count is asserted to equal dPred[M,λ1,λ2] exactly.
   Returns an Association keyed by {M,n} -> <|"J"->J_n,"Psi"->eigenvector normalized to last
   component 1|>. Memoized per (λ1,λ2). *)
ClearAll[TauEigensystem];
TauEigensystem::rankdef = "Combination C0+r C1 failed to resolve the `1`-dim eigenspace for (λ1,λ2)=(`2`,`3`) at every tried r in `4`.";
TauEigensystem::noninteger = "M = `1` is not within tolerance of an integer for eigenvector #`2` at (λ1,λ2)=(`3`,`4`) (combo r=`5`) — indicates a bug in the τ0 u^1-coefficient formula or the t2 construction.";
Psi::notgenuine = "Eigenvector #`1` for (λ1,λ2)=(`2`,`3`) failed the common-eigenvector residual check against C0/C1 (residual=`4`, combo r=`5`) — the generic-combination diagonalization should make this impossible; investigate before trusting results.";
Psi::degnull = "Weight sector M=`1` for (λ1,λ2)=(`2`,`3`) has `4` genuine state(s) but the predicted weight multiplicity dPred is `5`.";
Psi::lastzero = "Genuine eigenvector for (M,n)=(`1`,`2`) at (λ1,λ2)=(`3`,`4`) has last component ~0; cannot normalize.";
TauEigensystem[λ1_, λ2_] := TauEigensystem[λ1, λ2] = Module[
   {d, Cs, C0, C1, C2, rCandidates, r, vals, vecs, recs, byM, table},
   d = (λ1 + 1) (λ2 + 1);
   Cs = t2Coeffs[λ1, λ2]; {C0, C1, C2} = Cs;
   rCandidates = {7/3, N[Pi, 30], N[E, 30], 11/5, 13/4, N[Sqrt[2], 30]};
   r = SelectFirst[rCandidates,
     Module[{ee, vv}, {ee, vv} = Eigensystem[C0 + # C1];
       MatrixRank[vv] == d && Length[Union[Chop[ee, 10^-12]]] == d] &,
     Missing["NoGenericR"]];
   If[MissingQ[r], Message[TauEigensystem::rankdef, d, λ1, λ2, rCandidates]; Abort[]];
   {vals, vecs} = Eigensystem[N[C0 + r C1, 30]];
   recs = Table[
     Module[{v = vecs[[i]], t0, t1, res0, res1, Mi},
       t0 = (v . (C0 . v))/(v . v);
       t1 = (v . (C1 . v))/(v . v);
       res0 = Max[Abs[C0 . v - t0 v]];
       res1 = Max[Abs[C1 . v - t1 v]];
       If[Max[res0, res1] > 10^-15,
         Message[Psi::notgenuine, i, λ1, λ2, Max[res0, res1], r]; Abort[]];
       Mi = (-(θ1 + θ2) (κ1 + κ2) - h κ1 (λ1 + λ2) - t1)/(h (κ2 - κ1));
       If[Abs[Mi - Round[Re[Mi]]] > 10^-10,
         Message[TauEigensystem::noninteger, Mi, i, λ1, λ2, r]; Abort[]];
       <|"M" -> Round[Re[Mi]], "J" -> t0, "vec" -> v|>],
     {i, 1, d}];
   byM = GroupBy[recs, #["M"] &];
   table = Association @ Flatten @ KeyValueMap[
     Function[{Mval, group},
       Module[{sorted = SortBy[group, {Re[#["J"]], Im[#["J"]]} &], predicted = dPred[Mval, λ1, λ2]},
         If[Length[sorted] != predicted,
           Message[Psi::degnull, Mval, λ1, λ2, Length[sorted], predicted]; Abort[]];
         MapIndexed[
           Function[{rec, idx}, Module[{n = First[idx] - 1, vv = rec["vec"], last},
             last = vv[[-1]];
             If[Abs[last] < 10^-15 Max[Abs[vv]],
               Message[Psi::lastzero, Mval, n, λ1, λ2]; Abort[]];
             {Mval, n} -> <|"J" -> rec["J"], "Psi" -> vv/last|>]],
           sorted]]],
     byM];
   If[Length[table] != d, Message[Psi::degnull, "total", λ1, λ2, Length[table], d]; Abort[]];
   table
];
(* smoke test: (λ1,λ2)=(1,1) table keys and J values *)
KeySort[TauEigensystem[1, 1]][[All, "J"]] // N
```
**Output:**
```
(no output)
```

## 2026-07-15 10:28:50
**Input:**
```wolfram
(* Psi[λ1,λ2][M,n]: on-demand accessor over the memoized TauEigensystem cache. Returns the fixed
   (last-component-normalized) eigenvector for a genuine state.
   τ[λ1,λ2][M,n][u]: the FULL τ eigenvalue polynomial for that state, τ0[λ1,λ2,M][u] + J_n.
   Bounds/failure modes (Message + $Failed, matching the existing bounds-checked style):
     Psi::bnds     — gross bounds violation: M not an integer in 0..λ1+λ2, or n not an integer in 0..M.
     Psi::spurious — REPURPOSED from the old Baxter-matrix-filtering meaning: now fires when n is
                     within the gross 0..M range but exceeds the ACTUAL weight multiplicity
                     dPred[M,λ1,λ2]-1 for that M-sector (e.g. (λ1,λ2,M,n)=(2,1,2,2) — dPred[2,2,1]=2
                     so only n=0,1 exist even though the gross bound allows n up to M=2). There is no
                     other "spurious" failure mode left: every state actually produced by
                     TauEigensystem is, by construction, a genuine common eigenvector of C0/C1/C2.
   Both pull from the SAME cached TauEigensystem[λ1,λ2] table, so Psi and τ are mutually consistent
   by construction. *)
ClearAll[Psi, τ];
Psi::bnds = "Invalid indices (M,n)=(`1`,`2`) for (λ1,λ2)=(`3`,`4`): need M an integer in 0..λ1+λ2 and n an integer in 0..M.";
τ::bnds = "Invalid indices (M,n)=(`1`,`2`) for (λ1,λ2)=(`3`,`4`): need M an integer in 0..λ1+λ2 and n an integer in 0..M.";
Psi::spurious = "No genuine state at (M,n)=(`1`,`2`) for (λ1,λ2)=(`3`,`4`): weight sector M=`1` has only `5` genuine state(s) (n=0..`6`), so n=`2` does not exist.";
τ::spurious = "No genuine state at (M,n)=(`1`,`2`) for (λ1,λ2)=(`3`,`4`): weight sector M=`1` has only `5` genuine state(s) (n=0..`6`), so n=`2` does not exist.";
Psi[λ1_, λ2_][M_, n_] := Module[{tab, mult},
  If[! (IntegerQ[M] && IntegerQ[n] && 0 <= M <= λ1 + λ2 && 0 <= n <= M),
    Message[Psi::bnds, M, n, λ1, λ2]; Return[$Failed]];
  tab = TauEigensystem[λ1, λ2];
  mult = dPred[M, λ1, λ2];
  If[n > mult - 1,
    Message[Psi::spurious, M, n, λ1, λ2, mult, mult - 1]; Return[Missing["NoState", {λ1, λ2, M, n}]]];
  tab[{M, n}, "Psi"]];
τ[λ1_, λ2_][M_, n_][u_] := Module[{tab, mult},
  If[! (IntegerQ[M] && IntegerQ[n] && 0 <= M <= λ1 + λ2 && 0 <= n <= M),
    Message[τ::bnds, M, n, λ1, λ2]; Return[$Failed]];
  mult = dPred[M, λ1, λ2];
  If[n > mult - 1,
    Message[τ::spurious, M, n, λ1, λ2, mult, mult - 1]; Return[Missing["NoState", {λ1, λ2, M, n}]]];
  τ0[λ1, λ2, M][u] + TauEigensystem[λ1, λ2][{M, n}, "J"]];
(* smoke test: a valid state, plus a deliberately out-of-range call (n exceeds actual multiplicity
   for M=2 at (λ1,λ2)=(2,1), which is dPred[2,2,1]=2, so n=2 should fire Psi::spurious) *)
Print["spurious-guard test (expect a Psi::spurious message and Missing[...]): ", Psi[2, 1][2, 2]];
{Psi[1, 1][1, 0], N[τ[1, 1][1, 0][2], 10]}
```
**Output:**
```
(no output)
```

## 2026-07-15 10:28:50
**Input:**
```wolfram
(* QSolve[λ1,λ2][M,n]: the TQ equation is now a LINEAR solve, not an eigenvalue problem, because
   J_n (hence τ) is already KNOWN from the genuine TauEigensystem diagonalization above. The
   monic-degree-M polynomial q_n(u) satisfies O[q_n]=0 with O the SAME BaxterMatrix operator as
   before, now evaluated at the fixed known J_n: (BaxterMatrix[λ1,λ2,M] - J_n·Id).q_n = 0.
   Reuses BaxterMatrix (built above) and the bounds/multiplicity-check style from Psi/τ.
   CRITICAL numerical-precision pitfall (same theme as flagged for the old Psi construction, but a
   DIFFERENT concrete manifestation): BaxterMatrix[λ1,λ2,M]-J_n·Id is singular by construction, but
   its (M+1)-th (would-be-zero) entries/singular-values come out as precision-TRACKED zeros, e.g.
   literally `0``26` (value exactly 0, but carrying a finite Accuracy/Precision tag rather than
   Infinite precision). NullSpace treats such an entry's Precision (which evaluates to ~0, since
   Precision[0]-style values are accuracy-driven) as UNKNOWN rather than confirmed-zero, and
   conservatively returns an EMPTY nullspace even though the matrix is genuinely singular — this
   bites the M=0 sector (1x1 matrix, entry identically the precision-tracked zero) with EITHER
   default OR explicit Tolerance->10^-10 alone. The fix is to Chop the matrix (threshold 10^-12,
   comfortably above the ~10^-26..10^-29 noise floor of this construction but far below any
   genuine O(1) entry) BEFORE calling NullSpace, converting those precision-tracked zeros to exact
   (infinite-precision) 0 so NullSpace's linear-algebra can recognize the singularity. The explicit
   Tolerance->10^-10 is ALSO required per the design brief (matches the historical Psi-construction
   bug at (λ1,λ2,M,n)=(2,1,2,0) in the old stacked-C0/C1/C2 nullspace, and provides headroom for any
   genuinely-small-but-nonzero near-null direction that Chop alone would not catch in general).
   Both are used together defensively. Messages/behavior on failure mirror Psi/τ (bnds/spurious)
   plus new degnull/lastzero variants for this linear solve. Memoized per (λ1,λ2,M,n). *)
ClearAll[QSolve];
QSolve::bnds = "Invalid indices (M,n)=(`1`,`2`) for (λ1,λ2)=(`3`,`4`): need M an integer in 0..λ1+λ2 and n an integer in 0..M.";
QSolve::spurious = "No genuine state at (M,n)=(`1`,`2`) for (λ1,λ2)=(`3`,`4`): weight sector M=`1` has only `5` genuine state(s) (n=0..`6`), so n=`2` does not exist.";
QSolve::degnull = "Unexpected Q-function nullspace dimension `5` (expected 1) for (M,n)=(`1`,`2`) at (λ1,λ2)=(`3`,`4`).";
QSolve::lastzero = "Q-function nullspace vector found but its top (u^M) coefficient is ~0 for (M,n)=(`1`,`2`) at (λ1,λ2)=(`3`,`4`); cannot normalize to monic.";
QSolve[λ1_, λ2_][M_, n_] := QSolve[λ1, λ2][M, n] = Module[{mult, Jn, mat, ns, vec, last},
  If[! (IntegerQ[M] && IntegerQ[n] && 0 <= M <= λ1 + λ2 && 0 <= n <= M),
    Message[QSolve::bnds, M, n, λ1, λ2]; Return[$Failed]];
  mult = dPred[M, λ1, λ2];
  If[n > mult - 1,
    Message[QSolve::spurious, M, n, λ1, λ2, mult, mult - 1]; Return[$Failed]];
  Jn = TauEigensystem[λ1, λ2][{M, n}, "J"];
  mat = BaxterMatrix[λ1, λ2, M] - Jn IdentityMatrix[M + 1];
  ns = NullSpace[Chop[mat, 10^-12], Tolerance -> 10^-10];
  Which[
    Length[ns] != 1, Message[QSolve::degnull, M, n, λ1, λ2, Length[ns]]; Return[$Failed],
    True,
      vec = ns[[1]]; last = vec[[-1]];
      If[Abs[last] < 10^-15 Max[Abs[vec]],
        Message[QSolve::lastzero, M, n, λ1, λ2]; Return[$Failed]];
      vec/last]
];
(* smoke test: the regression-critical case (λ1,λ2,M,n)=(2,1,2,0) *)
QSolve[2, 1][2, 0]
```
**Output:**
```
(no output)
```

## 2026-07-15 10:28:51
**Input:**
```wolfram
(* Q1[λ1,λ2][M,n][u] = κ1^(u/h) q_n(u): the actual Baxter Q-function value, built from the SAME
   QSolve cache (mutually consistent with τ by construction, since both trace back to the same
   TauEigensystem[λ1,λ2][{M,n}] entry). *)
ClearAll[Q1];
Q1[λ1_, λ2_][M_, n_][u_] := Module[{q = QSolve[λ1, λ2][M, n]},
  If[q === $Failed, $Failed, κ1^(u/h) qpoly[q, u]]];
(* smoke test *)
{Q1[1, 1][1, 0][2], N[τ[1, 1][1, 0][2], 10]}
```
**Output:**
```
(no output)
```

---

## 2026-07-15 10:29:13 — 🔄 KERNEL RESTART

## 2026-07-15 10:29:32 — ▶️ RUN CELLS 1–34
Cell 1: ✓ — (no output)
Cell 3: ✓ — (no output)
Cell 4: ✓ — (no output)
Cell 5: ✓ — (no output)
Cell 6: ✓ — (no output)
Cell 8: ✓ — (no output)
Cell 10: ✓ — (no output)
Cell 11: ✓ — (no output)
Cell 12: ✓ — (no output)
Cell 13: ✓ — (no output)
Cell 14: ✓ — (no output)
Cell 15: ✓ — (no output)
Cell 16: ✓ — (no output)
Cell 17: ✓ — (no output)
Cell 18: ✓ — (no output)
Cell 20: ✓ — (no output)
Cell 21: ✓ — (no output)
Cell 23: ✓ — (no output)
Cell 25: ✓ — (no output)
Cell 26: ✓ — (no output)
Cell 28: ✓ — (no output)
Cell 29: ✓ — (no output)
Cell 30: ✓ — (no output)
Cell 31: ✓ — (no output)
Cell 32: ✓ — (no output)
Cell 33: ✓ — (no output)
Cell 34: ✓ — (no output)


## 2026-07-15 10:30:12
**Input:**
```wolfram
{Q1[1, 1][1, 0][2], N[τ[1, 1][1, 0][2], 10]}
```
**Output:**
```
Out= {-0.70923101095801692653722514224657913666`26.84522127451995 + 1.09596916251277206706752728049497901496`27.042725760166295*I, -0.02291660899982616316934999074083675684`10.15051499783199 + 0``11.790364642896865*I}
```

## 2026-07-15 10:30:20
**Input:**
```wolfram
QSolve[2, 1][2, 0]
```
**Output:**
```
Out= {1.03464899078135229594035120859018010109`26.919487563830202 + 0.38186956838094628389957594115425242959`26.919487563830202*I, -2.37022278949001677248172618524167157572`26.919487563830202 - 0.20048970475590097847093779292961086368`26.07936750871613*I, 1.00000000000000000000000000000000000001`26.919487563830202 + 0``27.106798363243854*I}
```
## 2026-07-15 10:31:03 — ▶️ RUN CELL 1
**In [1]:** `(*Numeric parameters (shared across the whole notebook). NOTE: this assigns global numeric values to θ1,θ2,h,κ1,κ2.↵  χ1,χ2 are the companion twist G's eigenvalues (Section A's G={{χ1,-χ2},{1,0}} has `
**Out:** (no output)


## 2026-07-15 10:31:32
**Input:**
```wolfram
Print["spurious-guard test (expect a Psi::spurious message and Missing[...]): ", Psi[2, 1][2, 2]]
```
**Output:**
```
Print:
*STR*spurious-guard test (expect a Psi::spurious message and Missing[...]): \034BoxData[InterpretationBox[TooltipBox["\342\200\224", "(missing)"], Missing["NoState", {2, 1, 2, 2}]]]
[message] Psi::spurious: No genuine state at (M,n)=(2,2) for (λ1,λ2)=(2,1): weight sector M=2 has only 2 genuine state(s) (n=0..1), so n=2 does not exist.
Out= Null
```

## 2026-07-15 10:31:33
**Input:**
```wolfram
{dPred[2,2,1], Length[TauEigensystem[2,1]], (λ1=2;λ2=1;(λ1+1)(λ2+1))}
```
**Output:**
```
Out= {2, 6, 6}
```
## 2026-07-15 10:32:28 — ✏️ EDIT CELL 1
(*Numeric parameters (shared across the whole notebook). NOTE: this assigns global numeric values to θ1,θ2,h,κ1,κ2.
  χ1,χ2 are the companion twist G's eigenvalues (Section A's G={{χ1,-χ2},{1,0}} has …


## 2026-07-15 10:32:28
**Input:**
```wolfram
(*Numeric parameters (shared across the whole notebook). NOTE: this assigns global numeric values to θ1,θ2,h,κ1,κ2.
  χ1,χ2 are the companion twist G's eigenvalues (Section A's G={{χ1,-χ2},{1,0}} has char. poly
  λ²-χ1 λ+χ2=(λ-κ1)(λ-κ2), so χ1=κ1+κ2, χ2=κ1κ2=1 here since κ2=1/κ1) — fixed here, numerically,
  once for the whole merged notebook (both Section A and Section B rely on G/t2 being numeric).*)
θ1 = 1/3;
θ2 = 1/7;
h = 1;
κ1 = N[Exp[I Zeta[3]], 30];
κ2 = 1/κ1;
χ1 = κ1 + κ2;
χ2 = 1;
{θ1, θ2, h, κ1, κ2, χ1, χ2}
```
**Output:**
```
Out= {1/3, 1/7, 1, 0.3604398751472265372666776505955500912937645232622850024938`29.707347828414424 + 0.9327824485933747164848359820708679294578496404879183613065`30.120295363558725*I, 0.360439875147226537266677650595550091293764523262285002494`29.421282762500656 - 0.9327824485933747164848359820708679294578496404879183613068`29.58721122792411*I, 0.7208797502944530745333553011911001825875290465245700049878`29.54117567596917 + 0``29.50582140452307*I, 1}
```

## 2026-07-15 10:33:13
**Input:**
```wolfram
{"Q1[1,1][1,0][2] match" -> Q1[1,1][1,0][2]-(-0.709231010958016926537225142+1.09596916251277206706752728 I) // Chop[#,10^-20]&, "QSolve[2,1][2,0] match" -> QSolve[2,1][2,0]-{1.03464899078135229594035121+0.381869568380946283899575941 I,-2.37022278949001677248172619-0.20048970475590097847093779 I,1.} // Chop[#,10^-20]&}
```
**Output:**
```
Out= {"Q1[1,1][1,0][2] match" -> 0, "QSolve[2,1][2,0] match" -> {0, 0, 0}}
```
