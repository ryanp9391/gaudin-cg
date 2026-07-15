# AI Action Log — FastQ

## 2026-07-12 04:40:57 — 📥 BULK INSERT 5 CELL(S) at positions 1–5
1. [markdown] # FastQ↵↵Clean, validated residue-evaluation primitives for su(4) rectangular-rep Clebsch-Gordan coe
2. [code] ClearAll[f, res1, res0];↵f[u_, n_, S1_, S2_, w_] := Exp[w u] u^n / ((u - 1)^S1 u^S2);↵res1[n_, S1_, 
3. [code] doRes0[S1_,S2_,f_]:=Module[{nmax,coeffs,res},↵  nmax=Exponent[f,u];↵  coeffs=CoefficientList[f,u];↵ 
4. [markdown] ## Validation↵↵`res1`/`res0` vs. direct `Residue[f[u,n,S1,S2,w], {u,1}]` / `{u,0}`, for `n ∈ {0,...,
5. [code] mismatches1 = Reap[↵   Do[↵    diff = Simplify[res1[n, S1, S2] - Residue[f[u, n, S1, S2, w], {u, 1}]


## 2026-07-12 04:40:58
**Input:**
```wolfram
ClearAll[f, res1, res0];
f[u_, n_, S1_, S2_, w_] := Exp[w u] u^n / ((u - 1)^S1 u^S2);
res1[n_, S1_, S2_] := Exp[w] LaguerreL[S1 - 1, n - S2 - S1 + 1, -w];
res0[n_, S1_, S2_] := If[S2 > n, (-1)^(S1 + S2 - n + 1) LaguerreL[S2 - n - 1, 1 - S1 - S2 + n, w], 0];
```
**Output:**
```
(no output)
```

## 2026-07-12 04:40:58
**Input:**
```wolfram
doRes0[S1_,S2_,f_]:=Module[{nmax,coeffs,res},
  nmax=Exponent[f,u];
  coeffs=CoefficientList[f,u];
  res=Plus@@(coeffs Table[res0[n,S1,S2],{n,0,nmax}]);
  res
  ]

doRes1[S1_,S2_,f_]:=Module[{nmax,coeffs,res},
  nmax=Exponent[f,u];
  coeffs=CoefficientList[f,u];
  res=Plus@@(coeffs Table[res1[n,S1,S2],{n,0,nmax}]);
  res
  ]
```
**Output:**
```
(no output)
```

## 2026-07-12 04:40:58
**Input:**
```wolfram
mismatches1 = Reap[
   Do[
    diff = Simplify[res1[n, S1, S2] - Residue[f[u, n, S1, S2, w], {u, 1}]];
    If[diff =!= 0, Sow[{n, S1, S2, diff}]],
    {n, 0, 4}, {S1, 1, 4}, {S2, 0, 4}
    ]
   ][[2]] // Flatten[#, 1] &;
mismatches0 = Reap[
   Do[
    diff = Simplify[res0[n, S1, S2] - Residue[f[u, n, S1, S2, w], {u, 0}]];
    If[diff =!= 0, Sow[{n, S1, S2, diff}]],
    {n, 0, 4}, {S1, 0, 4}, {S2, 0, 4}
    ]
   ][[2]] // Flatten[#, 1] &;
{Length[mismatches1], Length[mismatches0]}
```
**Output:**
```
(no output)
```
## 2026-07-12 04:41:06 — ▶️ RUN CELL 5
**In [5]:** `mismatches1 = Reap[↵   Do[↵    diff = Simplify[res1[n, S1, S2] - Residue[f[u, n, S1, S2, w], {u, 1}]];↵    If[diff =!= 0, Sow[{n, S1, S2, diff}]],↵    {n, 0, 4}, {S1, 1, 4}, {S2, 0, 4}↵    ]↵   ][[2]]`
**Out:** (no output)


## 2026-07-12 04:41:12
**Input:**
```wolfram
{Length[mismatches1], Length[mismatches0]}
```
**Output:**
```
Out= {0, 0}
```
## 2026-07-12 04:41:24 — ▶️ RUN CELLS 2–5
Cell 2: ✓ — (no output)
Cell 3: ✓ — (no output)
Cell 5: ✓ — (no output)

## 2026-07-12 04:41:41 — ✏️ EDIT CELL 5
mismatches1 = Reap[
   Do[
    diff = Simplify[res1[n, S1, S2] - Residue[f[u, n, S1, S2, w], {u, 1}]];
    If[diff =!= 0, Sow[{n, S1, S2, diff}]],
    {n, 0, 4}, {S1, 1, 4}, {S2, 0, 4}
    ]
   ][[2]]…


## 2026-07-12 04:41:41
**Input:**
```wolfram
mismatches1 = Reap[
   Do[
    diff = Simplify[res1[n, S1, S2] - Residue[f[u, n, S1, S2, w], {u, 1}]];
    If[diff =!= 0, Sow[{n, S1, S2, diff}]],
    {n, 0, 4}, {S1, 1, 4}, {S2, 0, 4}
    ]
   ][[2]] // Flatten[#, 1] &;
mismatches0 = Reap[
   Do[
    diff = Simplify[res0[n, S1, S2] - Residue[f[u, n, S1, S2, w], {u, 0}]];
    If[diff =!= 0, Sow[{n, S1, S2, diff}]],
    {n, 0, 4}, {S1, 0, 4}, {S2, 0, 4}
    ]
   ][[2]] // Flatten[#, 1] &;
Print["mismatches (res1, res0) = ", {Length[mismatches1], Length[mismatches0]}];
```
**Output:**
```
Out= Null
```

## 2026-07-12 04:42:00
**Input:**
```wolfram
{Length[mismatches1], Length[mismatches0]}
```
**Output:**
```
Out= {0, 0}
```
