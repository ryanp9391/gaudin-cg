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
## 2026-07-16 14:43:59 — 📥 BULK INSERT 11 CELL(S) at positions 35–45
1. [markdown] (*Section C — FSoV scalar product / orthogonality check. mu/Nn/nn/br/det ported verbatim from↵  Expe
2. [code] (*FSoV measures*)↵ClearAll[mu];↵mu[λ1_,λ2_][1][u_]:=Nn[λ1,λ2][1]/(Product[(u-θ1-h k),{k,0,λ1}]Produc
3. [code] (*Normalization factors in measure*)↵Nn[λ1_,λ2_][1]:=1/(2 π I)/(1-Exp[2 π I/h(θ1-θ2)]) nn[λ1,λ2][1];
4. [code] (*brackets of function f, given by contour integral around all poles in measure*)↵ClearAll[br];↵br[λ
5. [code] ClearAll[det];↵det[λ1_,λ2_][f_]:=1/(f/.u:>θ1)/(f/.u:>θ2)/(θ2-θ1)({{br[λ1,λ2][1][f],br[λ1,λ2][1][u f]
6. [code] (*Port-fidelity check: fsovL2CachedFormula is FSoV_L2.wb's own cached Out[215] for↵  det[1,1][f[u] z
7. [code] (*statesOf: the sorted {M,n} key list for a representation, read off TauEigensystem's memoized↵  cac
8. [code] states21 = statesOf[2, 1]
9. [code] G21 = GramMatrix[2, 1];↵Dimensions[G21]
10. [code] G21//Chop[#,10^(-30)]&
11. [code] (*Verification: G21 should be numerically diagonal. maxOffDiag is the largest |off-diagonal|↵  entry


## 2026-07-16 14:43:59
**Input:**
```wolfram
(*FSoV measures*)
ClearAll[mu];
mu[λ1_,λ2_][1][u_]:=Nn[λ1,λ2][1]/(Product[(u-θ1-h k),{k,0,λ1}]Product[(u-θ2-h k),{k,0,λ2}])(1-Exp[2 π I/h (u-θ2)]);
mu[λ1_,λ2_][2][u_]:=Nn[λ1,λ2][2]/(Product[(u-θ1-h k),{k,0,λ1}]Product[(u-θ2-h k),{k,0,λ2}])(1-Exp[2 π I/h (u-θ1)]);
```
**Output:**
```
(no output)
```

## 2026-07-16 14:43:59
**Input:**
```wolfram
(*Normalization factors in measure*)
Nn[λ1_,λ2_][1]:=1/(2 π I)/(1-Exp[2 π I/h(θ1-θ2)]) nn[λ1,λ2][1];
Nn[λ1_,λ2_][2]:=1/(2 π I)/(1-Exp[2 π I/h(θ2-θ1)]) nn[λ1,λ2][2];

nn[λ1_,λ2_][1]:= (-1)^(λ1+λ2-1)h^(λ1)Gamma[λ1+1] Product[θ2+k h-θ1,{k,0,λ2}]//Simplify;
nn[λ1_,λ2_][2]:= (-1)^(λ1+λ2-1)h^(λ2)Gamma[λ2+1] Product[θ1+k h-θ2,{k,0,λ1}]//Simplify;
```
**Output:**
```
(no output)
```

## 2026-07-16 14:44:00
**Input:**
```wolfram
(*brackets of function f, given by contour integral around all poles in measure*)
ClearAll[br];
br[λ1_,λ2_][1][f_]:=2 π I Plus@@(Residue[f mu[λ1,λ2][1][u],{u,#}]&/@Table[θ1+h k,{k,0,λ1}])//Expand//Simplify;
br[λ1_,λ2_][2][f_]:=2 π I Plus@@(Residue[f mu[λ1,λ2][2][u],{u,#}]&/@Table[θ2+h k,{k,0,λ2}])//Expand//Simplify;
```
**Output:**
```
(no output)
```

## 2026-07-16 14:44:00
**Input:**
```wolfram
ClearAll[det];
det[λ1_,λ2_][f_]:=1/(f/.u:>θ1)/(f/.u:>θ2)/(θ2-θ1)({{br[λ1,λ2][1][f],br[λ1,λ2][1][u f]},{br[λ1,λ2][2][f],br[λ1,λ2][2][u f]}}//Det//Expand//Simplify)
```
**Output:**
```
(no output)
```

## 2026-07-16 14:44:00
**Input:**
```wolfram
(*Port-fidelity check: fsovL2CachedFormula is FSoV_L2.wb's own cached Out[215] for
  det[1,1][f[u] z^(u/h)], transcribed verbatim with θ1,θ2,h,z,f still symbolic (as they were in
  that standalone notebook). Here θ1,θ2,h are already bound (cell 0), so comparing the freshly
  computed det[...] against this formula with the SAME substitution applied is a direct
  transcription check, independent of re-deriving the algebra.*)
ClearAll[fsovL2CachedFormula];
fsovL2CachedFormula[θ1_,θ2_,h_,z_,fn_] := -(z^((θ1+θ2)/h-θ1/h-θ2/h) (z fn[h+θ1] (fn[θ2] (h-θ1+θ2)+z (θ1-θ2) fn[h+θ2])+fn[θ1] (fn[θ2] (θ1-θ2)-z fn[h+θ2] (h+θ1-θ2))))/(fn[θ1] fn[θ2] (θ2-θ1));
portCheckDiff = Simplify[det[1,1][f[u] z^(u/h)] - fsovL2CachedFormula[θ1,θ2,h,z,f]]
```
**Output:**
```
(no output)
```

## 2026-07-16 14:44:01
**Input:**
```wolfram
(*statesOf: the sorted {M,n} key list for a representation, read off TauEigensystem's memoized
  cache (built once in Section B) rather than re-deriving dPred here.
  GramMatrix: the full Gram matrix G_{(M,n),(M',n')} = det[λ1,λ2][Q1(M,n)(u) Q1(M',n')(u)] over
  ALL states of the representation at once — this covers same-M degeneracy AND cross-M pairs in
  a single object, per the design.*)
ClearAll[statesOf, GramMatrix];
statesOf[λ1_,λ2_] := SortBy[Keys[TauEigensystem[λ1,λ2]], {First,Last}];
GramMatrix[λ1_,λ2_] := Module[{s = statesOf[λ1,λ2]},
  Table[
    det[λ1,λ2][Q1[λ1,λ2][s1[[1]],s1[[2]]][u] Q1[λ1,λ2][s2[[1]],s2[[2]]][u]],
    {s1, s}, {s2, s}]
];
```
**Output:**
```
(no output)
```

## 2026-07-16 14:44:01
**Input:**
```wolfram
states21 = statesOf[2, 1]
```
**Output:**
```
(no output)
```

## 2026-07-16 14:44:02
**Input:**
```wolfram
G21 = GramMatrix[2, 1];
Dimensions[G21]
```
**Output:**
```
(no output)
```

## 2026-07-16 14:44:02
**Input:**
```wolfram
G21//Chop[#,10^(-30)]&
```
**Output:**
```
(no output)
```

## 2026-07-16 14:44:02
**Input:**
```wolfram
(*Verification: G21 should be numerically diagonal. maxOffDiag is the largest |off-diagonal|
  entry over all 30 off-diagonal pairs (same-M AND cross-M, since states21 spans all of M=0..3);
  minDiagAbs is the smallest |diagonal| entry, confirming non-degeneracy (a genuine scalar
  product, not one that degenerates to 0 on some state). Threshold 10^-15 matches this notebook's
  existing hard-assert convention (e.g. Section B's TauEigensystem residual check, cell 30) — well
  above the ~10^-20..10^-27 noise floor expected at this notebook's 30-digit working precision,
  and far below any genuine O(1)-scale diagonal entry.*)
FSoVGram::notorthogonal = "Off-diagonal Gram entries do not vanish for (λ1,λ2)=(2,1): max|off-diag|=`1`, min|diag|=`2`.";
FSoVGram::degenerate = "A diagonal Gram entry is suspiciously small for (λ1,λ2)=(2,1) (min|diag|=`1`) — scalar product may be degenerate.";
offDiagVals = Flatten[Table[If[i != j, G21[[i, j]], Nothing], {i, 1, 6}, {j, 1, 6}]];
diagVals = Table[G21[[i, i]], {i, 1, 6}];
maxOffDiag = Max[Abs[N[offDiagVals]]];
minDiagAbs = Min[Abs[N[diagVals]]];
If[maxOffDiag > 10^-15,
  Message[FSoVGram::notorthogonal, maxOffDiag, minDiagAbs];
  Abort[]];
If[minDiagAbs < 10^-15,
  Message[FSoVGram::degenerate, minDiagAbs];
  Abort[]];
<|
  "states" -> states21,
  "max|off-diagonal|" -> maxOffDiag,
  "digits vanished (-Log10)" -> N[-Log10[maxOffDiag]],
  "min|diagonal|" -> minDiagAbs,
  "diagonal values" -> diagVals
|>
```
**Output:**
```
(no output)
```
## 2026-07-16 14:45:22 — 📥 BULK INSERT 7 CELL(S) at positions 46–52
1. [markdown] (*Section D — left eigenvectors of the companion-twisted transfer matrix, paired with the↵  existing
2. [code] ClearAll[TauEigensystemLeft];↵TauEigensystemLeft::rankdef = "Left combination C0+r C1 failed to reso
3. [code] (*smoke test: (λ1,λ2)=(2,1) must produce the same 6 keys as statesOf[2,1] (already verified in↵  Sec
4. [code] (*PsiLeft[λ1,λ2][M,n]: on-demand accessor over the memoized TauEigensystemLeft cache, same↵  bounds/
5. [code] (*Biorthogonality: a genuine left/right eigenvector pair of a diagonalizable operator with↵  simple 
6. [code] (*braket = <Psi|Psi>: plain bilinear pairing (no conjugation), matching this notebook's existing↵  n
7. [code] (*Final report: braket, G21 diagonal entry, and normRatio for all 6 states of (λ1,λ2)=(2,1).↵  Hard-


## 2026-07-16 14:45:22
**Input:**
```wolfram
ClearAll[TauEigensystemLeft];
TauEigensystemLeft::rankdef = "Left combination C0+r C1 failed to resolve the `1`-dim eigenspace for (λ1,λ2)=(`2`,`3`) at every tried r in `4`.";
TauEigensystemLeft::nomatch = "No left eigenvalue matched target for state `1` at (λ1,λ2)=(`2`,`3`) within tolerance (closest distance `4`).";
TauEigensystemLeft::dupmatch = "Left eigenvector index `4` matched more than one state (latest: `1`) at (λ1,λ2)=(`2`,`3`) — combo eigenvalues are not resolving distinctly.";
TauEigensystemLeft::lastzero = "Matched left eigenvector for state `1` at (λ1,λ2)=(`2`,`3`) has last component ~0; cannot normalize.";
TauEigensystemLeft::countmismatch = "TauEigensystemLeft produced `3` states for (λ1,λ2)=(`1`,`2`), expected `4`.";
TauEigensystemLeft[λ1_, λ2_] := TauEigensystemLeft[λ1, λ2] = Module[
   {d, Cs, C0, C1, C2, rCandidatesL, r, valsL, vecsL, states, usedIdx, table},
   d = (λ1 + 1) (λ2 + 1);
   Cs = t2Coeffs[λ1, λ2]; {C0, C1, C2} = Cs;
   rCandidatesL = {17/6, N[Pi, 30] + 1, N[E, 30] + 1, 9/4, 19/6, N[Sqrt[3], 30]};
   r = SelectFirst[rCandidatesL,
     Module[{ee, vv}, {ee, vv} = Eigensystem[Transpose[C0 + # C1]];
       MatrixRank[vv] == d && Length[Union[Chop[ee, 10^-12]]] == d] &,
     Missing["NoGenericR"]];
   If[MissingQ[r], Message[TauEigensystemLeft::rankdef, d, λ1, λ2, rCandidatesL]; Abort[]];
   {valsL, vecsL} = Eigensystem[Transpose[N[C0 + r C1, 30]]];
   states = statesOf[λ1, λ2];
   usedIdx = {};
   table = Association @ Table[
      Module[{Mn = states[[k]], target, dists, idx, w, last},
        target = τ[λ1, λ2][Mn[[1]], Mn[[2]]][r] - (κ1 + κ2) r^2;
        dists = Abs[valsL - target];
        idx = First[Ordering[dists]];
        If[dists[[idx]] > 10^-10,
          Message[TauEigensystemLeft::nomatch, Mn, λ1, λ2, dists[[idx]]]; Abort[]];
        If[MemberQ[usedIdx, idx],
          Message[TauEigensystemLeft::dupmatch, Mn, λ1, λ2, idx]; Abort[]];
        AppendTo[usedIdx, idx];
        w = vecsL[[idx]];
        last = w[[-1]];
        If[Abs[last] < 10^-15 Max[Abs[w]],
          Message[TauEigensystemLeft::lastzero, Mn, λ1, λ2]; Abort[]];
        Mn -> w/last],
      {k, 1, Length[states]}];
   If[Length[table] != d,
     Message[TauEigensystemLeft::countmismatch, λ1, λ2, Length[table], d]; Abort[]];
   table
];
```
**Output:**
```
(no output)
```

## 2026-07-16 14:45:23
**Input:**
```wolfram
(*smoke test: (λ1,λ2)=(2,1) must produce the same 6 keys as statesOf[2,1] (already verified in
  Section C), just via the independent left-eigenvector construction.*)
KeySort[TauEigensystemLeft[2, 1]] // Keys
```
**Output:**
```
(no output)
```

## 2026-07-16 14:45:23
**Input:**
```wolfram
(*PsiLeft[λ1,λ2][M,n]: on-demand accessor over the memoized TauEigensystemLeft cache, same
  bounds/spurious-guard pattern as Psi (Section B).*)
ClearAll[PsiLeft];
PsiLeft::bnds = "Invalid indices (M,n)=(`1`,`2`) for (λ1,λ2)=(`3`,`4`): need M an integer in 0..λ1+λ2 and n an integer in 0..M.";
PsiLeft::spurious = "No genuine state at (M,n)=(`1`,`2`) for (λ1,λ2)=(`3`,`4`): weight sector M=`1` has only `5` genuine state(s) (n=0..`6`), so n=`2` does not exist.";
PsiLeft[λ1_, λ2_][M_, n_] := Module[{tab, mult},
  If[! (IntegerQ[M] && IntegerQ[n] && 0 <= M <= λ1 + λ2 && 0 <= n <= M),
    Message[PsiLeft::bnds, M, n, λ1, λ2]; Return[$Failed]];
  tab = TauEigensystemLeft[λ1, λ2];
  mult = dPred[M, λ1, λ2];
  If[n > mult - 1,
    Message[PsiLeft::spurious, M, n, λ1, λ2, mult, mult - 1]; Return[Missing["NoState", {λ1, λ2, M, n}]]];
  tab[{M, n}]];
(*smoke test: a valid state's left eigenvector, plus confirming its last component is exactly 1*)
{PsiLeft[2, 1][1, 0][[-1]], Length[PsiLeft[2, 1][1, 0]]}
```
**Output:**
```
(no output)
```

## 2026-07-16 14:45:23
**Input:**
```wolfram
(*Biorthogonality: a genuine left/right eigenvector pair of a diagonalizable operator with
  simple spectrum must satisfy PsiLeft[M,n].Psi[M',n']=0 for {M,n}!={M',n'}. This is the
  correctness check on TauEigensystemLeft's construction itself (independent of Section C's
  det[Q,Q] machinery) — the natural companion to Section C's off-diagonal Gram check. Threshold
  10^-15 matches this notebook's existing hard-assert convention.*)
BiorthoCheck::fail = "Off-diagonal <PsiLeft|Psi> pairing does not vanish for (λ1,λ2)=(2,1): max|off-diag|=`1`.";
offBiortho = Flatten[Table[
   If[i != j, PsiLeft[2, 1] @@ states21[[i]] . (Psi[2, 1] @@ states21[[j]]), Nothing],
   {i, 1, 6}, {j, 1, 6}]];
maxOffBiortho = Max[Abs[N[offBiortho]]];
If[maxOffBiortho > 10^-15,
  Message[BiorthoCheck::fail, maxOffBiortho];
  Abort[]];
<|"max|<PsiLeft|Psi>| off-diagonal" -> maxOffBiortho,
  "digits vanished (-Log10)" -> N[-Log10[maxOffBiortho]]|>
```
**Output:**
```
(no output)
```

## 2026-07-16 14:45:24
**Input:**
```wolfram
(*braket = <Psi|Psi>: plain bilinear pairing (no conjugation), matching this notebook's existing
  non-Hermitian convention. normRatio divides by G21's diagonal entry for the SAME state (G21 is
  Section C's already-verified Gram matrix over states21) — only valid for (λ1,λ2)=(2,1) since
  G21/states21 are fixed globals for that representation, matching this plan's fixed test case.*)
ClearAll[braket, normRatio];
braket[λ1_, λ2_][M_, n_] := PsiLeft[λ1, λ2][M, n] . Psi[λ1, λ2][M, n];
normRatio[λ1_, λ2_][M_, n_] := Module[{i = First[FirstPosition[states21, {M, n}]]},
  braket[λ1, λ2][M, n]/G21[[i, i]]];
(*smoke test*)
{braket[2, 1][1, 0], normRatio[2, 1][1, 0]}
```
**Output:**
```
(no output)
```

## 2026-07-16 14:45:24
**Input:**
```wolfram
(*Final report: braket, G21 diagonal entry, and normRatio for all 6 states of (λ1,λ2)=(2,1).
  Hard-assert every braket is nonzero (non-degenerate pairing) before reporting — a zero braket
  here would mean PsiLeft and Psi are somehow orthogonal even on the diagonal, which should be
  impossible given Task 2's biorthogonality check already passed.*)
BraketReport::degenerate = "braket vanished for state `1` at (λ1,λ2)=(2,1) (value `2`) — pairing may be degenerate.";
reportRows = Table[
   Module[{Mn = states21[[i]], b, g, ratio},
     b = braket[2, 1] @@ Mn;
     g = G21[[i, i]];
     ratio = normRatio[2, 1] @@ Mn;
     If[Abs[N[b]] < 10^-15,
       Message[BraketReport::degenerate, Mn, b];
       Abort[]];
     <|"M" -> Mn[[1]], "n" -> Mn[[2]], "braket" -> N[b, 10], "G21 diagonal" -> N[g, 10],
       "normRatio" -> N[ratio, 10]|>],
   {i, 1, 6}];
reportRows // Column
```
**Output:**
```
(no output)
```
## 2026-07-16 14:46:10 — 📥 BULK INSERT 8 CELL(S) at positions 53–60
1. [markdown] (*Section E — operatorial factorization of the L=2 SoV basis. Builds a genuine STANDALONE↵  length-1
2. [code] (*Single-site (L=1) fused transfer matrices — the L=1 analogue of Section A's two-site t2 block,↵  b
3. [code] (*smoke: shape check for the fused hierarchy*)↵{Dimensions[t1[2][1,1][u]], Dimensions[t1[2][1,2][u]]
4. [code] (*Genuine length-1 SoV covector (paper normalization, 1-site Qθ1(u)=u-θ1). xSingle[λ][{0}] is the↵  
5. [code] (*Internal check (independent of the factorization): xL1[λ][{s}] must be a LEFT eigencovector of↵  B
6. [code] (*Operatorial factorization: x2[λ1,λ2][{s1,0}] = c (xL1[λ1][{s1}] ⊗ xSingle[λ2][{0}]). Flatten both↵
7. [code] (*Summary: worst ratio-constancy spread across the whole sweep (the hard-asserted factorization↵  co
8. [code] (*Physics-level closure: the site-2 pseudovacuum xSingle[λ2][{0}] is a LEFT eigencovector of the↵  s


## 2026-07-16 14:46:10
**Input:**
```wolfram
(*Single-site (L=1) fused transfer matrices — the L=1 analogue of Section A's two-site t2 block,
  built from the SAME single-site Lax L[λ] (no θ baked in) and shared companion twist G. qdet1 is
  the single-site quantum determinant: both factors are this site's L evaluated at u-θ1 and u-θ1-h
  (the internal -h shift), unlike the two-site qdetT where θ is already in T. t1[λ][1,1](u)=tr_a[G
  L(u-θ1)] matches paper eq. t_{1,1}(u)=tr_a[G L(u-θ)].*)
ClearAll[t1, qdet1];
t1[λ_][0, s_][u_] := IdentityMatrix[dim[{λ, 0}]];
t1[λ_][a_, 0][u_] := IdentityMatrix[dim[{λ, 0}]];
t1[λ_][1, 1][u_] := Sum[L[λ][i, j][u - θ1] G[[j, i]], {i, 2}, {j, 2}];
qdet1[λ_][u_] := L[λ][1, 1][u - θ1] . L[λ][2, 2][u - θ1 - h] - L[λ][2, 1][u - θ1] . L[λ][1, 2][u - θ1 - h];
t1[λ_][2, 1][u_] := χ2 qdet1[λ][u];
t1[λ_][1, s_ /; s > 1][u_] := t1[λ][1, s - 1][u] . t1[λ][1, 1][u + (s - 1) h] - t1[λ][2, 1][u + (s - 1) h] . t1[λ][1, s - 2][u];
```
**Output:**
```
(no output)
```

## 2026-07-16 14:46:10
**Input:**
```wolfram
(*smoke: shape check for the fused hierarchy*)
{Dimensions[t1[2][1,1][u]], Dimensions[t1[2][1,2][u]]}
```
**Output:**
```
(no output)
```

## 2026-07-16 14:46:10
**Input:**
```wolfram
(*Genuine length-1 SoV covector (paper normalization, 1-site Qθ1(u)=u-θ1). xSingle[λ][{0}] is the
  site lowest-weight covector (Section A) — a 1×(λ+1) row. xL1 stays that shape: (1×n).(n×n)=1×n.
  B1 is the L=1 SoV B-operator.*)
ClearAll[B1, xL1];
B1[λ_][u_] := L[λ][1, 1][u - θ1];
xL1[λ_][{0}] := xSingle[λ][{0}];
xL1[λ_][{s_}] := xL1[λ][{s}] = xSingle[λ][{0}] . t1[λ][1, s][θ1] / Product[Qθ1[θ1 - λ h + (p - 1) h], {p, 1, s}];
```
**Output:**
```
(no output)
```

## 2026-07-16 14:46:11
**Input:**
```wolfram
(*Internal check (independent of the factorization): xL1[λ][{s}] must be a LEFT eigencovector of
  B1 with eigenvalue u-θ1-s h — the paper's L=1 SoV property (verified there symbolically for
  λ=1..5). This validates the whole t1/xL1 construction before it is used downstream. Symbolic in
  u; residual Chopped at the 10^-20 noise floor.*)
xL1BCheck = Table[
   xL1[λ][{s}] . B1[λ][u] - (u - θ1 - s h) xL1[λ][{s}],
   {λ, 1, 3}, {s, 0, λ}] // ExpandAll // Chop[#, 10^-20] & // Flatten // Union;
xL1BCheck
```
**Output:**
```
(no output)
```

## 2026-07-16 14:46:11
**Input:**
```wolfram
(*Operatorial factorization: x2[λ1,λ2][{s1,0}] = c (xL1[λ1][{s1}] ⊗ xSingle[λ2][{0}]). Flatten both
  to length-D vectors, take the component ratio over the nonzero components of the tensor product.
  HARD requirement #1: the ratio lhs/tp must be a CONSTANT c (spread below the 10^-15 floor) — the
  actual factorization content ("t2[1,s1][θ1] acts only through a rank-1 site-2 eigenvalue"). HARD
  requirement #2 (hardened after the μEff≠1 finding was pinned down): c must equal the CORRECTED
  scalar cPredCorr[λ1,λ2,s1] = Product[(θ1-θ2+(k-1)h-λ2 h)/(θ1-θ2+(k-1)h-λ1 h),{k,1,s1}], which
  accounts for the genuine site-2 passthrough eigenvalue μSite2[λ2,s1] (defined/verified at the
  operator level in the cell below) — NOT the naive 1/Product[(θ1-λ1 h+(k-1)h)-θ2,{k,1,s1}] that
  assumed a trivial (identity) site-2 passthrough. cPredCorr reduces to 1 at λ1=λ2 and to
  f(λ2)/f(λ1), f(λ)=21λ-4, at s1=1 — matching the earlier empirical finding exactly. mueff=cval/cpred
  is now a closure check (must be exactly 1 to the noise floor) rather than a soft warning.*)
ClearAll[cPredCorr];
cPredCorr[λ1_, λ2_, s1_] := Product[(θ1 - θ2 + (k - 1) h - λ2 h)/(θ1 - θ2 + (k - 1) h - λ1 h), {k, 1, s1}];
SoVFact::notconst = "Ratio lhs/tp not constant for (λ1,λ2,s1)=(`1`,`2`,`3`): spread=`4`.";
SoVFact::mu = "c ≠ cPredCorr for (λ1,λ2,s1)=(`1`,`2`,`3`): μEff=`4` (corrected scalar formula failed — a construction bug, not a finding).";
factResults = Flatten[Table[
    Module[{lhs, tp, nz, rl, cval, cpred, mueff, spread},
     lhs = Flatten[x2[λ1, λ2][{s1, 0}]];
     tp = Flatten[KroneckerProduct[xL1[λ1][{s1}], xSingle[λ2][{0}]]];
     nz = Position[tp, x_ /; Abs[Chop[x, 10^-20]] > 10^-12];
     rl = Extract[lhs, nz]/Extract[tp, nz];
     cval = Mean[rl];
     spread = Max[Abs[rl - cval]];
     cpred = cPredCorr[λ1, λ2, s1];
     mueff = cval/cpred;
     If[spread > 10^-15, Message[SoVFact::notconst, λ1, λ2, s1, N[spread]]; Abort[]];
     If[Abs[N[mueff - 1]] > 10^-15, Message[SoVFact::mu, λ1, λ2, s1, N[mueff]]; Abort[]];  (* now HARD *)
     {λ1, λ2, s1, cval, cpred, N[mueff], N[spread]}],
    {λ1, 1, 3}, {λ2, 1, 3}, {s1, 0, λ1}], 2];
Column[factResults]
```
**Output:**
```
(no output)
```

## 2026-07-16 14:46:11
**Input:**
```wolfram
(*Summary: worst ratio-constancy spread across the whole sweep (the hard-asserted factorization
  content), and the range of μEff=c/cPred (should be exactly 1 everywhere if the site-2 passthrough
  is trivial). Columns of factResults: {λ1,λ2,s1,c,cPred,μEff,spread}.*)
<|
  "max spread (ratio non-constancy)" -> Max[factResults[[All, 7]]],
  "max |μEff - 1| over sweep" -> Max[Abs[factResults[[All, 6]] - 1]],
  "num (λ1,λ2,s1) cases checked" -> Length[factResults]
|>
```
**Output:**
```
(no output)
```

## 2026-07-16 14:46:12
**Input:**
```wolfram
(*Physics-level closure: the site-2 pseudovacuum xSingle[λ2][{0}] is a LEFT eigencovector of the
  site-2 Lax factor appearing in t2[λ1,λ2][1,s1][θ1] at each fused step, with eigenvalue
  μSite2[λ2,s1] = Product[(θ1-θ2+(k-1)h-λ2 h),{k,1,s1}]. This is WHY the naive "passthrough is the
  identity" guess failed: L[λ2][1,1](θ1-θ2) has a nontrivial GT lowest-weight eigenvalue at each of
  the s1 fused steps (shifted by (k-1)h from the internal Hirota recursion), not 1. Verified directly
  at the NUMERATOR level (i.e. before any xL1-normalization is applied — the cleanest statement of the
  effect): x2[λ1,λ2][{0,0}].t2[λ1,λ2][1,s1][θ1] must equal μSite2[λ2,s1] times
  (xSingle[λ1][{0}].t1[λ1][1,s1][θ1]) ⊗ xSingle[λ2][{0}], componentwise constant over the nonzero
  tensor-product entries. HARD assert on both constancy (spread) and the closed form (deviation
  from μSite2).*)
ClearAll[μSite2];
μSite2[λ2_, s1_] := Product[(θ1 - θ2 + (k - 1) h - λ2 h), {k, 1, s1}];
Site2Mu::notconst = "μ-ratio not constant at operator level for (λ1,λ2,s1)=(`1`,`2`,`3`): spread=`4`.";
Site2Mu::mismatch = "μ ≠ μSite2 at operator level for (λ1,λ2,s1)=(`1`,`2`,`3`): μ=`4`, μSite2=`5`.";
opMuResults = Flatten[Table[
    Module[{numLhs, numTp, nz, rl, muval, mupred, spread},
     numLhs = Flatten[x2[λ1, λ2][{0, 0}] . t2[λ1, λ2][1, s1][θ1]];
     numTp = Flatten[KroneckerProduct[xSingle[λ1][{0}] . t1[λ1][1, s1][θ1], xSingle[λ2][{0}]]];
     nz = Position[numTp, x_ /; Abs[Chop[x, 10^-20]] > 10^-12];
     rl = Extract[numLhs, nz]/Extract[numTp, nz];
     muval = Mean[rl];
     spread = Max[Abs[rl - muval]];
     mupred = μSite2[λ2, s1];
     If[spread > 10^-15, Message[Site2Mu::notconst, λ1, λ2, s1, N[spread]]; Abort[]];
     If[Abs[N[muval - mupred]] > 10^-15, Message[Site2Mu::mismatch, λ1, λ2, s1, N[muval], N[mupred]]; Abort[]];
     {λ1, λ2, s1, N[muval], N[mupred], N[spread]}],
    {λ1, 1, 3}, {λ2, 1, 3}, {s1, 0, λ1}], 2];
<|"max operator-level spread" -> Max[opMuResults[[All, 6]]],
  "max |μ - μSite2|" -> Max[Abs[opMuResults[[All, 4]] - opMuResults[[All, 5]]]],
  "num (λ1,λ2,s1) cases" -> Length[opMuResults]|>
```
**Output:**
```
(no output)
```
## 2026-07-16 14:47:08 — 📥 BULK INSERT 11 CELL(S) at positions 61–71
1. [markdown] (*Section F — functional CG overlap. Sets up the length-1 FSoV measure/bracket (Section C's↵  mu[...
2. [code] (*Length-1 FSoV measure and bracket. mu1 is a single-site rational measure with simple poles at↵  th
3. [code] (*Internal check: the leading (θ1-node) weight must be exactly 1 for λ=1,2,3 (this IS the↵  br1[f]=f
4. [code] (*Sanity: br1 acting on a generic analytic f reproduces the weighted tower sum Σ_j w_j f(θ1+h j),↵  
5. [code] (*RHS operatorial CG overlap: rotate the site-1 pseudovacuum covector by g=exp[φ Ee[λ1][a,b]]↵  (Ee[
6. [code] (*Which single-site generator is nontrivial? The lowest-weight covector ⟨x0| annihilates one of↵  E1
7. [code] (*OPg fixes the nontrivial generator so downstream code reads OPg[λ1,λ2][M,n][φ]. Degree in φ is ≤λ1
8. [code] (*Functional object: the 1-site FSoV DETERMINANT det1[f] = br1[f]/f(θ1) — the 1-site reduction of↵  
9. [code] (*λ1=1 universality test (the sharp discovery step), rep (1,2), all 6 states: det1 is 1+D_state k;↵ 
10. [code] (*Validate the universal φ(k)=αVal k across the FULL rep sweep and ALL their states. This is the↵  r
11. [code] (*Final dictionary report: the functional CG overlap and the operatorial CG overlap match EXACTLY↵


## 2026-07-16 14:47:08
**Input:**
```wolfram
(*Length-1 FSoV measure and bracket. mu1 is a single-site rational measure with simple poles at
  the tower θ1+h k, k=0..λ. Nn1 is DERIVED (not guessed) by forcing the θ1-node weight
  w_0 = 2πI Res[mu1, θ1] to equal 1, so that br1[λ][f] = f(θ1) + (higher-tower). br1 is the residue
  sum over the whole tower. (For λ≥1 the residues of mu1 sum to 0, so br1[λ][1]=0 — the correct
  normalization statement is w_0=1, checked next.)*)
ClearAll[Nn1, mu1, br1];
Nn1[λ_] := Nn1[λ] = Module[{uu}, 1/(2 π I Residue[1/Product[(uu - θ1 - h k), {k, 0, λ}], {uu, θ1}])];
mu1[λ_][u_] := Nn1[λ]/Product[(u - θ1 - h k), {k, 0, λ}];
br1[λ_][f_] := 2 π I Plus @@ (Residue[f mu1[λ][u], {u, #}] & /@ Table[θ1 + h k, {k, 0, λ}]);
```
**Output:**
```
(no output)
```

## 2026-07-16 14:47:09
**Input:**
```wolfram
(*Internal check: the leading (θ1-node) weight must be exactly 1 for λ=1,2,3 (this IS the
  br1[f]=f(θ1)+... normalization), and br1[λ][1] must be 0 (rational-measure residues sum to
  zero). Also list the full weight vectors w_j for the record.*)
w1[λ_] := Table[2 π I Residue[mu1[λ][u], {u, θ1 + h j}], {j, 0, λ}];
Nn1Check::badlead = "Leading weight w_0 ≠ 1 for λ=`1`: w_0=`2`.";
Nn1Check::badsum = "br1[λ][1] ≠ 0 for λ=`1`: got `2`.";
Do[
  Module[{w = w1[λ], s},
   If[Abs[N[w[[1]] - 1]] > 10^-15, Message[Nn1Check::badlead, λ, w[[1]]]; Abort[]];
   s = br1[λ][1];
   If[Abs[N[s]] > 10^-15, Message[Nn1Check::badsum, λ, s]; Abort[]]],
  {λ, 1, 3}];
<|"w_0 (λ=1,2,3)" -> Table[w1[λ][[1]], {λ, 1, 3}],
  "full weights w_j" -> Table[λ -> w1[λ], {λ, 1, 3}],
  "br1[λ][1] (λ=1,2,3)" -> Table[br1[λ][1], {λ, 1, 3}]|>
```
**Output:**
```
(no output)
```

## 2026-07-16 14:47:09
**Input:**
```wolfram
(*Sanity: br1 acting on a generic analytic f reproduces the weighted tower sum Σ_j w_j f(θ1+h j),
  and picks out f(θ1) as the leading term. Test with f=g[u] (undefined head, stays symbolic).*)
Clear[g];
{br1[1][g[u]] - Sum[w1[1][[j + 1]] g[θ1 + h j], {j, 0, 1}] // Simplify,
 br1[2][g[u]] - Sum[w1[2][[j + 1]] g[θ1 + h j], {j, 0, 2}] // Simplify}
```
**Output:**
```
(no output)
```

## 2026-07-16 14:47:09
**Input:**
```wolfram
(*RHS operatorial CG overlap: rotate the site-1 pseudovacuum covector by g=exp[φ Ee[λ1][a,b]]
  (Ee[a,b] nilpotent of order ≤λ1+1 on the finite-dim rep). Build g as the explicit truncated
  series Σ_{m=0}^{λ} φ^m/m! Ee[a,b]^m — this is EXACT (higher powers vanish) and guaranteed to be a
  polynomial in φ, avoiding MatrixExp returning an unevaluated symbolic exponential. Tensor with the
  site-2 pseudovacuum covector, contract with the two-site right eigenvector Psi. Shapes:
  (1×n1).(n1×n1)=1×n1 → KroneckerProduct with 1×n2 → Flatten to length D → .Psi = scalar.
  NOTE (kernel quirk, not a typo): MatrixPower[singularMatrix,0] does NOT auto-reduce to
  IdentityMatrix in this kernel (fires MatrixPower::sing and stays unevaluated) even though m=0 is
  a concrete nonnegative integer, while m>=1 evaluates cleanly. Ee[λ][a,b] is nilpotent/singular for
  the generators used here, so the m=0 identity term is split off explicitly and only
  MatrixPower is summed over m=1..λ, which evaluates without any message.*)
ClearAll[gRot, OP];
gRot[λ_][a_, b_][φ_] := IdentityMatrix[dim[{λ, 0}]] + Sum[φ^m/m! MatrixPower[Ee[λ][a, b], m], {m, 1, λ}];
OP[λ1_, λ2_][a_, b_][M_, n_][φ_] := Flatten[KroneckerProduct[xSingle[λ1][{0}] . gRot[λ1][a, b][φ], xSingle[λ2][{0}]]] . Psi[λ1, λ2][M, n];
```
**Output:**
```
(no output)
```

## 2026-07-16 14:47:10
**Input:**
```wolfram
(*Which single-site generator is nontrivial? The lowest-weight covector ⟨x0| annihilates one of
  E12/E21 (giving φ-independent OP); keep the other. Test on (λ1,λ2)=(2,1), state (1,0).*)
op12 = OP[2, 1][1, 2][1, 0][φ] // Simplify;
op21 = OP[2, 1][2, 1][1, 0][φ] // Simplify;
genNontrivial = Which[
   !FreeQ[op12, φ] && FreeQ[op21, φ], {1, 2},
   !FreeQ[op21, φ] && FreeQ[op12, φ], {2, 1},
   True, {1, 2}  (* fallback: if both depend on φ, default to {1,2} and flag in the report *)];
<|"OP(E12)" -> op12, "OP(E21)" -> op21, "nontrivial generator {a,b}" -> genNontrivial,
  "both φ-dependent?" -> (!FreeQ[op12, φ] && !FreeQ[op21, φ])|>
```
**Output:**
```
(no output)
```

## 2026-07-16 14:47:10
**Input:**
```wolfram
(*OPg fixes the nontrivial generator so downstream code reads OPg[λ1,λ2][M,n][φ]. Degree in φ is ≤λ1.*)
ClearAll[OPg];
OPg[λ1_, λ2_][M_, n_][φ_] := OP[λ1, λ2][genNontrivial[[1]], genNontrivial[[2]]][M, n][φ];
(*smoke: degree in φ for a λ1=1 and a λ1=2 case*)
{Exponent[OPg[1, 1][1, 0][φ] // Simplify, φ], Exponent[OPg[2, 1][1, 0][φ] // Simplify, φ]}
```
**Output:**
```
(no output)
```

## 2026-07-16 14:47:10
**Input:**
```wolfram
(*Functional object: the 1-site FSoV DETERMINANT det1[f] = br1[f]/f(θ1) — the 1-site reduction of
  Section C's det[λ1,λ2][f] (which divides by f at BOTH nodes θ1,θ2; only the θ1-node survives at
  1 site). This 1/f(θ1) normalization cancels the (k κ1)^(θ1/h) fractional prefactor that the bare
  bracket br1[f] would otherwise carry, leaving det1[k^(u/h) Q1] = Σ_j w_j k^j Q1(θ1+jh)/Q1(θ1), a
  CLEAN degree-λ1 polynomial in k with leading term exactly 1 (since w_0=1 by Nn1's construction).*)
ClearAll[det1, FLdet];
det1[λ1_][f_] := br1[λ1][f]/(f /. u -> θ1);
FLdet[λ1_, λ2_][M_, n_][k_] := det1[λ1][k^(u/h) Q1[λ1, λ2][M, n][u]];
(*smoke: FLdet[1,2][0,0][k] must be a degree-1 polynomial in k, constant term exactly 1, NO
  fractional power of k.*)
FLdet[1, 2][0, 0][k] // Simplify
```
**Output:**
```
(no output)
```

## 2026-07-16 14:47:11
**Input:**
```wolfram
(*λ1=1 universality test (the sharp discovery step), rep (1,2), all 6 states: det1 is 1+D_state k;
  OPg is 1+c1(state) φ. HARD assert #1: both constant terms are exactly 1 (a structural fact of the
  normalization, not a coincidence). HARD assert #2: α_state := D_state/c1(state) is the SAME
  constant across all 6 states — this is the actual universality content.*)
sts12 = statesOf[1, 2];
DictDet::badconst = "Constant term ≠ 1 for state `1`: FLdet const=`2`, OPg const=`3`.";
DictDet::notuniv = "α_state not universal: spread=`1` (see alphaTab).";
alphaTab = Table[
   Module[{fl, op, D0, Dc, c0, c1, alpha},
    fl = FLdet[1, 2][s[[1]], s[[2]]][k] // Simplify;
    op = OPg[1, 2][s[[1]], s[[2]]][φ] // Simplify;
    D0 = Coefficient[fl, k, 0]; Dc = Coefficient[fl, k, 1];
    c0 = Coefficient[op, φ, 0]; c1 = Coefficient[op, φ, 1];
    If[Abs[N[D0 - 1]] > 10^-15 || Abs[N[c0 - 1]] > 10^-15,
     Message[DictDet::badconst, s, N[D0], N[c0]]; Abort[]];
    alpha = Dc/c1;
    {s, N[D0], N[c0], N[Dc], N[c1], N[alpha]}],
   {s, sts12}];
alphaSpread = Max[Abs[alphaTab[[All, 6]] - alphaTab[[1, 6]]]];
If[alphaSpread > 10^-12, Message[DictDet::notuniv, alphaSpread]; Abort[]];
αVal = alphaTab[[1, 6]];
<|"alphaTab {state,D0,c0,Dc,c1,α}" -> alphaTab, "α spread" -> alphaSpread, "α (universal value)" -> αVal|>
```
**Output:**
```
(no output)
```

## 2026-07-16 14:47:11
**Input:**
```wolfram
(*Validate the universal φ(k)=αVal k across the FULL rep sweep and ALL their states. This is the
  real content: it tests that the binomial tower weights w_m (from Nn1/mu1) match OPg's φ^m Taylor
  coefficients (built from Ee[a,b]^m via the fixed gRot) under the single rescaling α, for every
  (λ1,λ2,M,n). HARD assert: residual FLdet-OPg[αVal k] vanishes (Simplify reduces it to literal 0,
  or Chop below 10^-15) for every state of every rep in the sweep.*)
ClearAll[φuniversal];
φuniversal[k_] := αVal k;
DictDet::residual = "Nonzero residual for (λ1,λ2,M,n)=(`1`,`2`,`3`,`4`): resid=`5`.";
repSweepDet = {{1, 1}, {2, 1}, {1, 2}, {2, 2}, {3, 1}, {1, 3}};
validationDet = Table[
   Module[{λ1 = rep[[1]], λ2 = rep[[2]], sts, resids},
    sts = statesOf[λ1, λ2];
    resids = Table[
      Module[{lhs, rhs, r},
       lhs = FLdet[λ1, λ2][st[[1]], st[[2]]][k];
       rhs = OPg[λ1, λ2][st[[1]], st[[2]]][φuniversal[k]];
       r = Chop[Simplify[lhs - rhs], 10^-15];
       If[r =!= 0, Message[DictDet::residual, λ1, λ2, st[[1]], st[[2]], r]; Abort[]];
       {st, r}],
      {st, sts}];
    {λ1, λ2} -> resids],
   {rep, repSweepDet}];
allResidsDet = Flatten[Table[validationDet[[i, 2]][[All, 2]], {i, Length[validationDet]}]];
maxResidDet = Max[Abs[allResidsDet]];
<|"validation" -> validationDet, "max residual over full sweep" -> maxResidDet,
  "num states checked" -> Length[allResidsDet]|>
```
**Output:**
```
(no output)
```

## 2026-07-16 14:47:11
**Input:**
```wolfram
(*Final dictionary report: the functional CG overlap and the operatorial CG overlap match EXACTLY
  (symbolically zero after Simplify, not merely below a noise floor) once the LHS uses the 1-site
  FSoV DETERMINANT det1 (not the bare bracket br1) and the parameter map is the single LINEAR
  rescaling φ(k)=α k. This is the Section F dictionary: generator, φ(k), and normalization.*)
dictReportDet = <|
  "generator {a,b}" -> genNontrivial,
  "φ(k)" -> "α k",
  "α (universal value)" -> αVal,
  "normalization" -> "1-site FSoV determinant det1[f] = br1[f]/(f/.u->θ1); Section C's det[λ1,λ2][f] restricted to 1 site",
  "max residual over rep sweep" -> maxResidDet,
  "reps validated" -> repSweepDet,
  "num states validated" -> Length[allResidsDet]
|>;
dictReportDet
```
**Output:**
```
(no output)
```

---

## 2026-07-16 14:47:32 — 🔄 KERNEL RESTART

## 2026-07-16 14:47:49 — ▶️ RUN CELLS 1–34
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


## 2026-07-16 14:48:05
**Input:**
```wolfram
{θ1, θ2, h, κ1, TauEigensystem[1,1][{0,0}], Psi[1,1][1,0]}
```
**Output:**
```
Out= {1/3, 1/7, 1, 0.3604398751472265372666776505955500912937645232622850024938`29.707347828414424 + 0.9327824485933747164848359820708679294578496404879183613065`30.120295363558725*I, <|"J" -> 0.56640551808849884427620773665015014352`26.563908186381784 + 1.37696456697117220052523406877128122917`26.915783392193127*I, "Psi" -> {-0.74016619280770349120031215726511215204`28.25736179271083 + 0.67242397862104047713498006666438147933`28.24691480308365*I, 0.36043987514722653726667765059555009129`27.843010048988198 + 0.93278244859337471648483598207086792946`28.246958303210636*I, 0.36043987514722653726667765059555009129`27.843010048988198 + 0.93278244859337471648483598207086792946`28.246958303210636*I, 1.`28.25273175422156 + 0``28.253729389534094*I}|>, {1.00000000000000000000000000000000000001`28.402747934397283 + 0``28.703777930061257*I, -1.30044948810764362720923581419356337479`28.402747934397286 + 0``28.589684442012377*I, 2.02132923840209670174259111538466355738`28.359954584956583 + 0``28.69917086764009*I, 1.`28.402747934397283 + 0``28.703777930061257*I}}
```
## 2026-07-16 14:48:20 — ▶️ RUN CELLS 35–52
Cell 36: ✓ — (no output)
Cell 37: ✓ — (no output)
Cell 38: ✓ — (no output)
Cell 39: ✓ — (no output)
Cell 40: ✓ — (no output)
Cell 41: ✓ — (no output)
Cell 42: ✓ — (no output)
Cell 43: ✓ — (no output)
Cell 44: ✓ — (no output)
Cell 45: ✓ — (no output)
Cell 47: ✓ — (no output)
Cell 48: ✓ — (no output)
Cell 49: ✓ — (no output)
Cell 50: ✓ — (no output)
Cell 51: ✓ — (no output)
Cell 52: ✓ — (no output)

## 2026-07-16 14:48:33 — ▶️ RUN CELLS 53–71
Cell 54: ✓ — (no output)
Cell 55: ✓ — (no output)
Cell 56: ✓ — (no output)
Cell 57: ✓ — (no output)
Cell 58: ✓ — (no output)
Cell 59: ✓ — (no output)
Cell 60: ✓ — (no output)
Cell 62: ✓ — (no output)
Cell 63: ✓ — (no output)
Cell 64: ✓ — (no output)
Cell 65: ✓ — (no output)
Cell 66: ✓ — (no output)
Cell 67: ✓ — (no output)
Cell 68: ✓ — (no output)
Cell 69: ✓ — (no output)
Cell 70: ✓ — (no output)
Cell 71: ✓ — (no output)


## 2026-07-16 14:48:53
**Input:**
```wolfram
{
 "xL1BCheck" -> xL1BCheck,
 "factResults max spread/mueff" -> {Max[factResults[[All,7]]], Max[Abs[factResults[[All,6]]-1]]},
 "opMuResults check" -> {Max[opMuResults[[All,6]]], Max[Abs[opMuResults[[All,4]]-opMuResults[[All,5]]]]},
 "w1 check" -> {Table[w1[λ][[1]],{λ,1,3}], Table[br1[λ][1],{λ,1,3}]},
 "genNontrivial" -> genNontrivial,
 "alphaSpread,αVal" -> {alphaSpread, αVal},
 "maxResidDet" -> maxResidDet,
 "dictReportDet" -> dictReportDet,
 "G21 offdiag/diag" -> {maxOffDiag, minDiagAbs},
 "biortho" -> maxOffBiortho,
 "portCheckDiff" -> portCheckDiff
}
```
**Output:**
```
Out= {"xL1BCheck" -> {0}, "factResults max spread/mueff" -> {0., 0.}, "opMuResults check" -> {0., 0.}, "w1 check" -> {{1, 1, 1}, {0, 0, 0}}, "genNontrivial" -> {2, 1}, "alphaSpread,αVal" -> {0., -1. + 0.*I}, "maxResidDet" -> 0, "dictReportDet" -> <|"generator {a,b}" -> {2, 1}, "φ(k)" -> "α k", "α (universal value)" -> -1. + 0.*I, "normalization" -> "1-site FSoV determinant det1[f] = br1[f]/(f/.u->θ1); Section C's det[λ1,λ2][f] restricted to 1 site", "max residual over rep sweep" -> 0, "reps validated" -> {{1, 1}, {2, 1}, {1, 2}, {2, 2}, {3, 1}, {1, 3}}, "num states validated" -> 41|>, "G21 offdiag/diag" -> {5.065342385207116*^-34, 3.0051232400647923}, "biortho" -> 0., "portCheckDiff" -> 0}
```
## 2026-07-16 14:49:08 — ▶️ RUN CELL 32
**In [32]:** `(* Psi[λ1,λ2][M,n]: on-demand accessor over the memoized TauEigensystem cache. Returns the fixed↵   (last-component-normalized) eigenvector for a genuine state.↵   τ[λ1,λ2][M,n][u]: the FULL τ eigenva`
**Out:** (no output)


## 2026-07-16 14:49:15
**Input:**
```wolfram
Print["spurious-guard test (expect a Psi::spurious message and Missing[...]): ", Psi[2, 1][2, 2]];
{Psi[1, 1][1, 0], N[τ[1, 1][1, 0][2], 10]}
```
**Output:**
```
Print:
*STR*spurious-guard test (expect a Psi::spurious message and Missing[...]): \034BoxData[InterpretationBox[TooltipBox["\342\200\224", "(missing)"], Missing["NoState", {2, 1, 2, 2}]]]
[message] Psi::spurious: No genuine state at (M,n)=(2,2) for (λ1,λ2)=(2,1): weight sector M=2 has only 2 genuine state(s) (n=0..1), so n=2 does not exist.
Out[1]= (suppressed)
Out[2]= {{1.00000000000000000000000000000000000001`28.402747934397283 + 0``28.703777930061257*I, -1.30044948810764362720923581419356337479`28.402747934397286 + 0``28.589684442012377*I, 2.02132923840209670174259111538466355738`28.359954584956583 + 0``28.69917086764009*I, 1.`28.402747934397283 + 0``28.703777930061257*I}, -0.02291660899982616316934999074083675684`10.15051499783199 + 0``11.790364642896865*I}
```

---

## 2026-07-16 17:04:33 — 🔄 KERNEL RESTART

## 2026-07-16 17:24:20 — ✏️ EDIT CELL 30
(*a(u)=aFun is now an ALIAS for Section A's two-site ν1 (Yangian weight eigenvalue of T[1,1] on
  the HWS) — same formula (u-θ1-hλ1)(u-θ2-hλ2), deduped so there is a single source of truth
  instead o…


## 2026-07-16 17:24:20
**Input:**
```wolfram
(*a(u)=aFun is now an ALIAS for Section A's two-site ν1 (Yangian weight eigenvalue of T[1,1] on
  the HWS) — same formula (u-θ1-hλ1)(u-θ2-hλ2), deduped so there is a single source of truth
  instead of two copies that could drift. d(u)=Qθ(u), and τ0 = τ minus the J term (using
  ee1=λ1+λ2-M, ee2=M). Qθ (=d) is defined LOCALLY here so this file is self-contained.
  qpoly builds a polynomial from an ASCENDING coefficient list {c0,...,cM} (cM=1 for monic q_n);
  formed in a dummy symbol first so substituting uu=0 gives the constant term rather than an
  indeterminate 0^0.*)
ClearAll[aFun, τ0, Qθ, qpoly];
aFun[λ1_, λ2_][u_] := ν1[λ1, λ2][u];
Qθ[u_] := (u - θ1) (u - θ2);
τ0[λ1_, λ2_, M_][u_] := (κ1 + κ2) u^2 - (θ1 + θ2) (κ1 + κ2) u - u h (κ1 (λ1 + λ2 - M) + κ2 M);
qpoly[clist_, uu_] := Module[{x}, (clist . x^Range[0, Length[clist] - 1]) /. x -> uu];
{aFun[1, 1][u], Qθ[u], τ0[1, 1, 1][u]}
```
**Output:**
```
Out= {(-4/3 + u)*(-8/7 + u), (-1/3 + u)*(-1/7 + u), (-1.064155821863240252882572111282100269533971449631508102601`29.54117567596917 + 0``29.336679005422717*I)*u + (0.7208797502944530745333553011911001825875290465245700049878`29.54117567596917 + 0``29.50582140452307*I)*u^2}
```
## 2026-07-16 17:24:24 — ✏️ EDIT CELL 58
(*Genuine length-1 SoV covector (paper normalization, 1-site Qθ1(u)=u-θ1). xSingle[λ][{0}] is the
  site lowest-weight covector (Section A) — a 1×(λ+1) row. xL1 stays that shape: (1×n).(n×n)=1×n.
  B1…


## 2026-07-16 17:24:25
**Input:**
```wolfram
(*Genuine length-1 SoV covector (paper normalization, 1-site Qθ1(u)=u-θ1). xSingle[λ][{0}] is the
  site lowest-weight covector (Section A) — a 1×(λ+1) row. xL1 stays that shape: (1×n).(n×n)=1×n.
  B1 is the L=1 SoV B-operator. The single-site ν1[λ][u] defined here is the arity-1 companion to
  Section A's two-site ν1[λ1,λ2][u] (same symbol, different arity — both DownValues coexist; NOT
  ClearAll'd here so Section A's two-site pattern survives). Qθ1[θ1-λh+(p-1)h] and
  ν1[λ][θ1+(p-1)h] give IDENTICAL numeric values (Qθ1 is linear, so shifting the argument vs. the
  function is equivalent) — xL1 is rewritten here purely for notational parity with the Section A
  x2 fix pattern; xL1 itself never had a bug.*)
ClearAll[B1, xL1];
B1[λ_][u_] := L[λ][1, 1][u - θ1];
ν1[λ_][u_] := (u - θ1 - h λ);
xL1[λ_][{0}] := xSingle[λ][{0}];
xL1[λ_][{s_}] := xL1[λ][{s}] = xSingle[λ][{0}] . t1[λ][1, s][θ1] / Product[ν1[λ][θ1 + (p - 1) h], {p, 1, s}];
```
**Output:**
```
Out= Null
```
## 2026-07-16 17:24:34 — ▶️ RUN CELL 59 (0.05 s)
**In [59]:** `(*Internal check (independent of the factorization): xL1[λ][{s}] must be a LEFT eigencovector of↵  B1 with eigenvalue u-θ1-s h — the paper's L=1 SoV property (verified there symbolically for↵  λ=1..5)`
**Out:** Out[126]= \left\{0\right\}


## 2026-07-16 17:24:36
**Input:**
```wolfram
testRatio1 = Flatten[Table[
   Module[{lhs, tp, nz, rl},
    lhs = Flatten[x2[λ1, λ2][{s1, 0}]];
    tp = Flatten[KroneckerProduct[xL1[λ1][{s1}], xSingle[λ2][{0}]]];
    nz = Position[tp, x_ /; Abs[Chop[x, 10^-20]] > 10^-12];
    rl = Extract[lhs, nz]/Extract[tp, nz];
    {λ1, λ2, s1, N[Max[Abs[rl - 1]]]}],
   {λ1, 1, 3}, {λ2, 1, 3}, {s1, 0, λ1}], 2];
{Max[testRatio1[[All,4]]], testRatio1}
```
**Output:**
```
Out= {0., {{1, 1, 0, 0.}, {1, 1, 1, 0.}, {1, 2, 0, 0.}, {1, 2, 1, 0.}, {1, 3, 0, 0.}, {1, 3, 1, 0.}, {2, 1, 0, 0.}, {2, 1, 1, 0.}, {2, 1, 2, 0.}, {2, 2, 0, 0.}, {2, 2, 1, 0.}, {2, 2, 2, 0.}, {2, 3, 0, 0.}, {2, 3, 1, 0.}, {2, 3, 2, 0.}, {3, 1, 0, 0.}, {3, 1, 1, 0.}, {3, 1, 2, 0.}, {3, 1, 3, 0.}, {3, 2, 0, 0.}, {3, 2, 1, 0.}, {3, 2, 2, 0.}, {3, 2, 3, 0.}, {3, 3, 0, 0.}, {3, 3, 1, 0.}, {3, 3, 2, 0.}, {3, 3, 3, 0.}}}
```

## 2026-07-16 17:25:23
**Input:**
```wolfram
muSite2Formula[λ2_,s1_] := Product[(θ1-θ2+(k-1)h-λ2 h),{k,1,s1}];
testOpMu = Flatten[Table[
   Module[{numLhs,numTp,nz,rl,muval},
    numLhs = Flatten[x2[λ1,λ2][{0,0}] . t2[λ1,λ2][1,s1][θ1]];
    numTp = Flatten[KroneckerProduct[xSingle[λ1][{0}] . t1[λ1][1,s1][θ1], xSingle[λ2][{0}]]];
    nz = Position[numTp, x_ /; Abs[Chop[x,10^-20]] > 10^-12];
    rl = Extract[numLhs,nz]/Extract[numTp,nz];
    muval = Mean[rl];
    {λ1,λ2,s1, N[muval - muSite2Formula[λ2,s1]]}],
   {λ1,1,3},{λ2,1,3},{s1,0,λ1}], 2];
Max[Abs[testOpMu[[All,4]]]]
```
**Output:**
```
Out= 0.
```
## 2026-07-16 17:25:52 — ✏️ EDIT CELL 60 [batch]
(*Operatorial factorization: x2[λ1,λ2][{s1,0}] = xL1[λ1][{s1}] ⊗ xSingle[λ2][{0}] — a PURE tensor
  factorization with NO leftover scalar. (An earlier version of this notebook found an apparent
  nont…


## 2026-07-16 17:25:52
**Input:**
```wolfram
(*Operatorial factorization: x2[λ1,λ2][{s1,0}] = xL1[λ1][{s1}] ⊗ xSingle[λ2][{0}] — a PURE tensor
  factorization with NO leftover scalar. (An earlier version of this notebook found an apparent
  nontrivial "site-2 passthrough eigenvalue" μEff here — that was an ARTIFACT of a Qθ normalization
  bug in x2's original definition, since fixed in Section A: x2 divided by Qθ(u)=(u-θ1)(u-θ2)
  evaluated at a point shifted by -λ1h in BOTH factors, incorrectly putting site-1's λ1 shift into
  the θ2 factor that should have carried site 2's OWN λ2 shift. The corrected x2 divides by
  ν1[λ1,λ2][u]=(u-θ1-hλ1)(u-θ2-hλ2) instead — each factor shifted by its own site's λ — and with
  that fix the ratio below comes out EXACTLY 1 for every (λ1,λ2,s1), not just a constant. There is
  no genuine extra site-2 eigenvalue beyond what ν1 already encodes: the old cPredCorr/μSite2
  "finding" was bookkeeping fallout from the bug, not new physics — see the cell below for the
  honest restatement of that operator identity.)
  HARD requirement: the ratio lhs/tp must be IDENTICALLY 1 (both non-constancy spread AND
  deviation-from-1 below the 10^-15 floor) for all 27 (λ1,λ2,s1) cases.*)
SoVFact::notconst = "Ratio lhs/tp not constant for (λ1,λ2,s1)=(`1`,`2`,`3`): spread=`4`.";
SoVFact::notone = "Ratio lhs/tp ≠ 1 for (λ1,λ2,s1)=(`1`,`2`,`3`): ratio=`4`.";
factResults = Flatten[Table[
    Module[{lhs, tp, nz, rl, cval, spread},
     lhs = Flatten[x2[λ1, λ2][{s1, 0}]];
     tp = Flatten[KroneckerProduct[xL1[λ1][{s1}], xSingle[λ2][{0}]]];
     nz = Position[tp, x_ /; Abs[Chop[x, 10^-20]] > 10^-12];
     rl = Extract[lhs, nz]/Extract[tp, nz];
     cval = Mean[rl];
     spread = Max[Abs[rl - cval]];
     If[spread > 10^-15, Message[SoVFact::notconst, λ1, λ2, s1, N[spread]]; Abort[]];
     If[Abs[N[cval - 1]] > 10^-15, Message[SoVFact::notone, λ1, λ2, s1, N[cval]]; Abort[]];
     {λ1, λ2, s1, N[cval], N[spread]}],
    {λ1, 1, 3}, {λ2, 1, 3}, {s1, 0, λ1}], 2];
Column[factResults]
```
**Output:**
```
Out= Out[130]= \begin{matrix}\left\{1,1,0,1.\,,0.\,\right\}\\\left\{1,1,1,1.\,,0.\,\right\}\\\left\{1,2,0,1.\,,0.\,\right\}\\\left\{1,2,1,1.\,,0.\,\right\}\\\left\{1,3,0,1.\,,0.\,\right\}\\\left\{1,3,1,1.\,,0.\,\right\}\\\left\{2,1,0,1.\,,0.\,\right\}\\\left\{2,1,1,1.\,,0.\,\right\}\\\left\{2,1,2,1.\,,0.\,\right\}\\\left\{2,2,0,1.\,,0.\,\right\}\\\left\{2,2,1,1.\,,0.\,\right\}\\\left\{2,2,2,1.\,,0.\,\r… [898 chars total]
```
## 2026-07-16 17:25:52 — ✏️ EDIT CELL 61 [batch]
(*Summary: worst deviation from the now-EXACT tensor factorization (ratio ≡ 1, no leftover
  scalar) across the whole sweep. Columns of factResults: {λ1,λ2,s1,ratio,spread}.*)
<|
  "max spread (ratio …


## 2026-07-16 17:25:52
**Input:**
```wolfram
(*Summary: worst deviation from the now-EXACT tensor factorization (ratio ≡ 1, no leftover
  scalar) across the whole sweep. Columns of factResults: {λ1,λ2,s1,ratio,spread}.*)
<|
  "max spread (ratio non-constancy)" -> Max[factResults[[All, 5]]],
  "max |ratio - 1| over sweep" -> Max[Abs[factResults[[All, 4]] - 1]],
  "num (λ1,λ2,s1) cases checked" -> Length[factResults]
|>
```
**Output:**
```
Out= Out[131]= \langle\!\lvert \text{max spread (ratio non-constancy)}\to 0.\,,\text{max |ratio - 1| over sweep}\to 0.\,,\text{num ($\lambda$1,$\lambda$2,s1) cases checked}\to 27\rvert\!\rangle
```
## 2026-07-16 17:25:52 — ✏️ EDIT CELL 62 [batch]
(*Operator-level identity (STILL TRUE, and unaffected by the Qθ normalization bug — it never used
  x2[{s1,0}]'s buggy denominator, only the raw x2[{0,0}] pseudovacuum and t2 itself): the site-2
  pse…


## 2026-07-16 17:25:52
**Input:**
```wolfram
(*Operator-level identity (STILL TRUE, and unaffected by the Qθ normalization bug — it never used
  x2[{s1,0}]'s buggy denominator, only the raw x2[{0,0}] pseudovacuum and t2 itself): the site-2
  pseudovacuum xSingle[λ2][{0}] passes through t2[λ1,λ2][1,s1][θ1] picking up the factor
  μSite2[λ2,s1] = Product[(θ1-θ2+(k-1)h-λ2 h),{k,1,s1}] relative to the bare tensor product
  (xSingle[λ1][{0}].t1[λ1][1,s1][θ1]) ⊗ xSingle[λ2][{0}]. CORRECTED READING (see the cell above):
  this is NOT an extra, unexplained site-2 eigenvalue — it is exactly the λ2-dependent factor of
  Section A's two-site ν1[λ1,λ2][u]=(u-θ1-hλ1)(u-θ2-hλ2) evaluated on the site-1 tower, i.e.
  ν1[λ1,λ2][θ1+(k-1)h] factors as ((k-1-λ1)h)·(θ1-θ2+(k-1-λ2)h), and the second factor here IS
  μSite2[λ2,s1]. Once x2's denominator uses the FULL two-site ν1 (Section A's fix) instead of a
  λ1-only Qθ factor, this piece is already built into the normalization, which is why the ratio in
  the cell above is exactly 1 — μSite2 was always just part of the two-site Yangian weight
  function's structure, not a surprising new eigenvalue. Kept here as a genuine, still-verified
  operator identity (hard-asserted below), relabeled honestly.*)
ClearAll[μSite2];
μSite2[λ2_, s1_] := Product[(θ1 - θ2 + (k - 1) h - λ2 h), {k, 1, s1}];
Site2Mu::notconst = "μ-ratio not constant at operator level for (λ1,λ2,s1)=(`1`,`2`,`3`): spread=`4`.";
Site2Mu::mismatch = "μ ≠ μSite2 at operator level for (λ1,λ2,s1)=(`1`,`2`,`3`): μ=`4`, μSite2=`5`.";
opMuResults = Flatten[Table[
    Module[{numLhs, numTp, nz, rl, muval, mupred, spread},
     numLhs = Flatten[x2[λ1, λ2][{0, 0}] . t2[λ1, λ2][1, s1][θ1]];
     numTp = Flatten[KroneckerProduct[xSingle[λ1][{0}] . t1[λ1][1, s1][θ1], xSingle[λ2][{0}]]];
     nz = Position[numTp, x_ /; Abs[Chop[x, 10^-20]] > 10^-12];
     rl = Extract[numLhs, nz]/Extract[numTp, nz];
     muval = Mean[rl];
     spread = Max[Abs[rl - muval]];
     mupred = μSite2[λ2, s1];
     If[spread > 10^-15, Message[Site2Mu::notconst, λ1, λ2, s1, N[spread]]; Abort[]];
     If[Abs[N[muval - mupred]] > 10^-15, Message[Site2Mu::mismatch, λ1, λ2, s1, N[muval], N[mupred]]; Abort[]];
     {λ1, λ2, s1, N[muval], N[mupred], N[spread]}],
    {λ1, 1, 3}, {λ2, 1, 3}, {s1, 0, λ1}], 2];
<|"max operator-level spread" -> Max[opMuResults[[All, 6]]],
  "max |μ - μSite2|" -> Max[Abs[opMuResults[[All, 4]] - opMuResults[[All, 5]]]],
  "num (λ1,λ2,s1) cases" -> Length[opMuResults]|>
```
**Output:**
```
Out= Out[137]= \langle\!\lvert \text{max operator-level spread}\to 0.\,,\text{max |$\mu$ - $\mu$Site2|}\to 0.\,,\text{num ($\lambda$1,$\lambda$2,s1) cases}\to 27\rvert\!\rangle
```

## 2026-07-16 17:25:59
**Input:**
```wolfram
dictReportDet
```
**Output:**
```
Out= dictReportDet
```

---

## 2026-07-16 17:26:09 — 🔄 KERNEL RESTART

## 2026-07-16 17:26:27 — ▶️ RUN CELLS 1–36
Cell 1: ✓ — Out[8]= \left\{\frac{1}{3},\frac{1}{7},1,0.360439875147226537266677650596+0.932782448593374716484835982071\,\mathrm{i},0.36043987514722653726667765060-0.932782448593374716484835982071\,\mathrm{i},0.720879750294453074533355301191+0.\times 10^{-30}\,\mathrm{i},1\right\}
Cell 3: ✓ — (no output)
Cell 4: ✓ — (no output)
Cell 5: ✓ — (no output)
Cell 6: ✓ — (no output)
Cell 8: ✓ — (no output)
Cell 10: ✓ — (no output)
Cell 11: ✓ — (no output)
Cell 12: ✓ — Out[27]= \left\{0\right\}
Cell 13: ✓ — (no output)
Cell 14: ✓ — (no output)
Cell 15: ✓ — Out[31]= \left\{0\right\}
Cell 16: ✓ — (no output)
Cell 17: ✓ — Out[34]= \begin{pmatrix}0 & 0 & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\end{pmatrix}
Cell 18: ✓ — (no output)
Cell 20: ✓ — (no output)
Cell 21: ✓ — (no output)
Cell 23: ✓ — (no output)
Cell 25: ✓ — (no output)
Cell 26: ✓ — (no output)
Cell 27: ✓ — Out[52]= \begin{pmatrix}0 & \frac{1}{2} & 0 & 0 & 0 & 0\end{pmatrix}
Cell 28: ✓ — Out[53]= \left\{0\right\}
Cell 30: ✓ — Out[59]= \left\{\left(u-\frac{4}{3}\right)\,\left(u-\frac{8}{7}\right),\left(u-\frac{1}{3}\right)\,\left(u-\frac{1}{7}\right),\left(0.720879750294453074533355301191+0.\times 10^{-30}\,\mathrm{i}\right)\,u^2-\left(1.06415582186324025288257211128+0.\times 10^{-29}\,\mathrm{i}\right)\,u\right\}
Cell 31: ✓ — Out[62]= \begin{pmatrix}1 & 1\\2 & 2\\3 & 3\end{pmatrix}
Cell 32: ✓ — Out[65]= \text{max|C2-($\kappa$1+$\kappa$2)Id|}\to 0.\times 10^{-29}
Cell 33: ✓ — Out[75]= \langle\!\lvert \left\{0,0\right\}\to 0.5664055180884988 +1.3769645669711723\,\mathrm{i},\left\{1,0\right\}\to -0.778123966451158+0.\,\,\mathrm{i},\left\{1,1\right\}\to 1.1900552523337025 +0.\,\,\mathrm{i},\left\{2,0\right\}\to 0.5664055180884988 -1.3769645669711723\,\mathrm{i}\rvert\!\rangle
Cell 34: ✓ — 2 kernel messages: Psi::spurious: No genuine state at (M,n)=(2,2) for (λ1,λ2)=(2,1): weight sector M=2 has only 2 genuine state(s) (n=0..1), so n=2 does not exist. [also: Psi::spurious] result: Print: spurious-guard test (expect a Psi::spurious message and Missing[...]):[formula] | Out[84]= \left\{\left\{1.000000000000000000000000000+0.\times 10^{-29}\,\mathrm{i},-1.3004494881076436272092358…
⚠️ stopped at Cell 34 — error detected (pass stopOnError:false to continue past errors)

## 2026-07-16 17:26:50 — ▶️ RUN CELLS 35–73
Cell 35: ✓ — Out[91]= \left\{1.03464899078135229594035121+0.381869568380946283899575941\,\mathrm{i},-2.37022278949001677248172619-0.20048970475590097847093779\,\mathrm{i},1.00000000000000000000000000+0.\times 10^{-27}\,\mathrm{i}\right\}
Cell 36: ✓ — Out[94]= \left\{-0.709231010958016926537225142+1.09596916251277206706752728\,\mathrm{i},-0.02291660900+0.\times 10^{-12}\,\mathrm{i}\right\}
Cell 38: ✓ — (no output)
Cell 39: ✓ — (no output)
Cell 40: ✓ — (no output)
Cell 41: ✓ — (no output)
Cell 42: ✓ — Out[109]= 0
Cell 43: ✓ — (no output)
Cell 44: ✓ — Out[113]= \begin{pmatrix}0 & 0\\1 & 0\\1 & 1\\2 & 0\\2 & 1\\3 & 0\end{pmatrix}
Cell 45: ✓ — Out[115]= \left\{6,6\right\}
Cell 46: ✓ — Out[116]= \begin{pmatrix}2.9090642893150832630837933001-5.8046200635659338430425538562\,\mathrm{i} & 0 & 0 & 0 & 0 & 0\\0 & 7.36608403016431791181-0.968046299794461400266\,\mathrm{i} & 0 & 0 & 0 & 0\\0 & 0 & 2.68702160063169134012-1.3456153260557858257\,\mathrm{i} & 0 & 0 & 0\\0 & 0 & 0 & 7.3660840301643179118107+0.9680462997944614002659\,\mathrm{i} & 0 & 0\\0 & 0 & 0 & 0 & 2.687021600631691340118… [537 chars total]
Cell 47: ✓ — Out[125]= \langle\!\lvert \text{states}\to \begin{pmatrix}0 & 0\\1 & 0\\1 & 1\\2 & 0\\2 & 1\\3 & 0\end{pmatrix},\text{max|off-diagonal|}\to 5.065342385207116\times 10^{-34},\text{digits vanished (-Log10)}\to 33.29539119374288,\text{min|diagonal|}\to 3.0051232400647923,\text{diagonal values}\to \left\{2.9090642893150832630837933001-5.8046200635659338430425538562\,\mathrm{i},7.36608403016431791181-0… [699 chars total]
Cell 49: ✓ — (no output)
Cell 50: ✓ — Out[133]= \begin{pmatrix}0 & 0\\1 & 0\\1 & 1\\2 & 0\\2 & 1\\3 & 0\end{pmatrix}
Cell 51: ✓ — Out[138]= \left\{1.00000000000000000000000000000,6\right\}
Cell 52: ✓ — Out[143]= \langle\!\lvert \text{max|<PsiLeft|Psi>| off-diagonal}\to 0.\,,\text{digits vanished (-Log10)}\to \mathrm{Indeterminate}\rvert\!\rangle
Cell 53: ✓ — Out[147]= \left\{7.366084030164317911810727273-0.9680462997944614002658874268\,\mathrm{i},1.00000000000000000000+0.\times 10^{-21}\,\mathrm{i}\right\}
Cell 54: ✓ — Out[150]= \begin{matrix}\langle\!\lvert \text{M}\to 0,\text{n}\to 0,\text{braket}\to 2.909064289-5.804620064\,\mathrm{i},\text{G21 diagonal}\to 2.909064289-5.804620064\,\mathrm{i},\text{normRatio}\to 1.000000000+0.\times 10^{-10}\,\mathrm{i}\rvert\!\rangle \\\langle\!\lvert \text{M}\to 1,\text{n}\to 0,\text{braket}\to 7.366084030-0.968046300\,\mathrm{i},\text{G21 diagonal}\to 7.366084030-0.9680463… [1247 chars total]
Cell 56: ✓ — (no output)
Cell 57: ✓ — Out[158]= \begin{pmatrix}3 & 3\\3 & 3\end{pmatrix}
Cell 58: ✓ — (no output)
Cell 59: ✓ — Out[165]= \left\{0\right\}
Cell 60: ✓ — Out[169]= \begin{matrix}\left\{1,1,0,1.\,,0.\,\right\}\\\left\{1,1,1,1.\,,0.\,\right\}\\\left\{1,2,0,1.\,,0.\,\right\}\\\left\{1,2,1,1.\,,0.\,\right\}\\\left\{1,3,0,1.\,,0.\,\right\}\\\left\{1,3,1,1.\,,0.\,\right\}\\\left\{2,1,0,1.\,,0.\,\right\}\\\left\{2,1,1,1.\,,0.\,\right\}\\\left\{2,1,2,1.\,,0.\,\right\}\\\left\{2,2,0,1.\,,0.\,\right\}\\\left\{2,2,1,1.\,,0.\,\right\}\\\left\{2,2,2,1.\,,0.\,\r… [898 chars total]
Cell 61: ✓ — Out[170]= \langle\!\lvert \text{max spread (ratio non-constancy)}\to 0.\,,\text{max |ratio - 1| over sweep}\to 0.\,,\text{num ($\lambda$1,$\lambda$2,s1) cases checked}\to 27\rvert\!\rangle
Cell 62: ✓ — Out[176]= \langle\!\lvert \text{max operator-level spread}\to 0.\,,\text{max |$\mu$ - $\mu$Site2|}\to 0.\,,\text{num ($\lambda$1,$\lambda$2,s1) cases}\to 27\rvert\!\rangle
Cell 64: ✓ — (no output)
Cell 65: ✓ — Out[185]= \langle\!\lvert \text{w\_0 ($\lambda$=1,2,3)}\to \left\{1,1,1\right\},\text{full weights w\_j}\to \left\{1\to \left\{1,-1\right\},2\to \left\{1,-2,1\right\},3\to \left\{1,-3,3,-1\right\}\right\},\text{br1[$\lambda$][1] ($\lambda$=1,2,3)}\to \left\{0,0,0\right\}\rvert\!\rangle
Cell 66: ✓ — Out[187]= \left\{0,0\right\}
Cell 67: ✓ — (no output)
Cell 68: ✓ — Out[194]= \langle\!\lvert \text{OP(E12)}\to 1.\,+0.\times 10^{-30}\,\mathrm{i},\text{OP(E21)}\to \left(-0.03626734060947795145020244939-2.004237714666237867354326095\,\mathrm{i}\right)\,\varphi ^2-\left(1.522150083623538480912456055-0.2442047957313920727733393355\,\mathrm{i}\right)\,\varphi +\left(1.\,+0.\times 10^{-30}\,\mathrm{i}\right),\text{nontrivial generator \{a,b\}}\to \left\{2,1\right\},\… [464 chars total]
Cell 69: ✓ — Out[197]= \left\{1,2\right\}
Cell 70: ✓ — Out[201]= \left(1.000000000000000000000000000+0.\times 10^{-27}\,\mathrm{i}\right)-\left(0.360439875147226537266677651+0.932782448593374716484835982\,\mathrm{i}\right)\,k
Cell 71: ✓ — Out[209]= \langle\!\lvert \text{alphaTab \{state,D0,c0,Dc,c1,$\alpha$\}}\to \begin{pmatrix}\left\{0,0\right\} & 1.\, +0.\,\,\mathrm{i} & 1.\, +0.\,\,\mathrm{i} & -0.36043987514722653-0.9327824485933747\,\mathrm{i} & 0.36043987514722653 +0.9327824485933747\,\mathrm{i} & -1.\,+0.\,\,\mathrm{i}\\\left\{1,0\right\} & 1.\, +0.\,\,\mathrm{i} & 1.\, +0.\,\,\mathrm{i} & 0.7101545203859353 -0.2371533418544… [1247 chars total]
Cell 72: ✓ — Out[217]= \langle\!\lvert \text{validation}\to \left\{\left\{1,1\right\}\to \begin{pmatrix}\left\{0,0\right\} & 0\\\left\{1,0\right\} & 0\\\left\{1,1\right\} & 0\\\left\{2,0\right\} & 0\end{pmatrix},\left\{2,1\right\}\to \begin{pmatrix}\left\{0,0\right\} & 0\\\left\{1,0\right\} & 0\\\left\{1,1\right\} & 0\\\left\{2,0\right\} & 0\\\left\{2,1\right\} & 0\\\left\{3,0\right\} & 0\end{pmatrix},\left\{1… [1247 chars total]
Cell 73: ✓ — Out[219]= \langle\!\lvert \text{generator \{a,b\}}\to \left\{2,1\right\},\text{$\varphi$(k)}\to \text{$\alpha$ k},\text{$\alpha$ (universal value)}\to -1.\,+0.\,\,\mathrm{i},\text{normalization}\to \text{1-site FSoV determinant det1[f] = br1[f]/(f/.u->$\theta$1); Section C's det[$\lambda$1,$\lambda$2][f] restricted to 1 site},\text{max residual over rep sweep}\to 0,\text{reps validated}\to \begin{… [510 chars total]

