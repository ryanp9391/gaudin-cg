# AI Action Log — Hamiltonians_SU4

## 2026-07-12 10:15:02 — ✏️ EDIT CELL 1 [batch]
## Representations of su(4)

Constructs finite-dimensional irreps of $\mathfrak{gl}(4)$ (equivalently $\mathfrak{su}(4)$ irreps, realized as $\mathfrak{gl}(4)$-modules) via Gelfand-Tsetlin patterns. `…

## 2026-07-12 10:15:02 — ✏️ EDIT CELL 4 [batch]
## Gaudin Lax matrix

Builds the single-site Gaudin Lax matrix $L(z)$ from the $\mathfrak{su}(4)$ generators $e_{ij} = J_{ij}$ acting on the representation `rep`: $L(z)_{ij} = p_j\,\delta_{ij} + \dfra…

## 2026-07-12 10:24:23 — ✏️ EDIT CELL 6 [batch]
M[{i_, j_}][f_] := D[f, z] KroneckerDelta[i, j] + L[i, j][z] . f;

## 2026-07-12 10:24:23 — ✏️ EDIT CELL 7 [batch]
applyMs = {
   MChain[{a__, b__, c__, d__}] :> M[a][M[b][M[c][M[d][f[z]]]]] //.
     Dot[a___, b_ + c_, d___] :> Dot[a, b, d] + Dot[a, c, d]
  };

## 2026-07-12 10:24:23 — ✏️ EDIT CELL 8 [batch]
EE[II_, JJ_] :=
  1/4! Sum[
     Signature[\[Sigma]] Times @@ (e[#] & /@ Transpose[{II[[\[Sigma]]], JJ}]),
     {\[Sigma], Permutations[{1, 2, 3, 4}]}
    ] MChain[Transpose[{II, JJ}]] /. e[{a_, b_}] …

## 2026-07-12 10:24:23 — ✏️ EDIT CELL 9 [batch]
detOp =
  (Sum[EE[II, JJ], {II, Tuples[Range[4], 4]}, {JJ, Tuples[Range[4], 4]}] //. applyMs // ExpandAll) //.
    Dot[a___, b_ + c_, d___] :> Dot[a, b, d] + Dot[a, c, d] //.
    Dot[a___, b_, d___] :…

## 2026-07-12 10:24:23 — ✏️ EDIT CELL 10 [batch]
H[1][z_] = Coefficient[detOp, f'''[z]];
H[2][z_] = Coefficient[detOp, f''[z]] // ExpandAll;
H[3][z_] = Coefficient[detOp, f'[z]] // ExpandAll;
H[4][z_] = Coefficient[detOp, f[z]] // ExpandAll;

## 2026-07-12 10:28:04 — 📥 BULK INSERT 1 CELL(S) at positions 6–6
1. [markdown] ## Quantum determinant and Talalaev's construction↵↵Constructs the quantum (Capelli) determinant $\d

## 2026-07-12 10:30:37 — 📥 BULK INSERT 1 CELL(S) at positions 3–3
1. [markdown] Sanity check: evaluates the diagonal (Cartan) generator $J_{1,1}$ on the $\mathfrak{su}(4)$ represen

## 2026-07-12 10:30:39 — 📥 BULK INSERT 1 CELL(S) at positions 8–8
1. [markdown] $M_{ij}[f] = \partial_z f\,\delta_{ij} + L(z)_{ij}\cdot f$, the $(i,j)$ entry of the operator $(\par

## 2026-07-12 10:30:40 — 📥 BULK INSERT 1 CELL(S) at positions 10–10
1. [markdown] Expands a length-4 index chain `MChain[{a,b,c,d}]` into the nested operator composition $M_a\circ M_

## 2026-07-12 10:30:42 — 📥 BULK INSERT 1 CELL(S) at positions 12–12
1. [markdown] The $\epsilon$-$\epsilon$ contraction term: sums $\mathrm{sign}(\sigma)$ over permutations $\sigma$

## 2026-07-12 10:30:43 — 📥 BULK INSERT 1 CELL(S) at positions 14–14
1. [markdown] Sums `EE[I,J]` over all index tuples $I,J\in\{1,\dots,4\}^4$, i.e. the identity $\det A = \frac{1}{N

## 2026-07-12 10:30:44 — 📥 BULK INSERT 1 CELL(S) at positions 16–16
1. [markdown] Reads off the coefficients of $f'''(z), f''(z), f'(z), f(z)$ to give $H_1(z),\dots,H_4(z)$, the Tala

## 2026-07-12 10:32:50 — ✏️ EDIT CELL 15
detOp =
  (Sum[EE[II, JJ], {II, Tuples[Range[4], 4]}, {JJ, Tuples[Range[4], 4]}] //. applyMs // ExpandAll) //.
    Dot[a___, b_ + c_, d___] :> Dot[a, b, d] + Dot[a, c, d] //.
    Dot[a___, b_, d___] :…

## 2026-07-12 10:54:01 — 📥 BULK INSERT 8 CELL(S) at positions 7–14
1. [markdown] ## SU(3) determinant limit: contour-integral building block↵↵Before specializing the quantum determi
2. [markdown] Defines the measure `mu[u,theta,h]` with simple poles at $u=\theta$ and $u=\theta+h$, and `gaudinMea
3. [code] mu[u_, theta_, h_] := 1/((u - theta) (u - theta - h));↵↵gaudinMeasureResidue[f_, theta_, h_] :=↵  Si
4. [markdown] Checks the reduction symbolically for a generic black-box function $f_0$: since both poles are simpl
5. [code] finiteDiff[f_, theta_, h_] := (f[theta + h] - f[theta])/h;↵↵genericResult = gaudinMeasureResidue[Fun
6. [markdown] Validates on concrete pole-free test functions (polynomial, exponential, trig — no test function may
7. [code] testFunctions = {Function[u, u^3], Function[u, Exp[u]], Function[u, Sin[u]]};↵↵(* symbolic check: ge
8. [markdown] Both checks vanish identically: the contour integral of a pole-free $f$ against the two-point Gaudin


## 2026-07-12 10:54:01
**Input:**
```wolfram
mu[u_, theta_, h_] := 1/((u - theta) (u - theta - h));

gaudinMeasureResidue[f_, theta_, h_] :=
  Simplify[
   Residue[f[u] mu[u, theta, h], {u, theta}] +
    Residue[f[u] mu[u, theta, h], {u, theta + h}]];
```
**Output:**
```
(no output)
```

## 2026-07-12 10:54:01
**Input:**
```wolfram
finiteDiff[f_, theta_, h_] := (f[theta + h] - f[theta])/h;

genericResult = gaudinMeasureResidue[Function[u, f0[u]], theta, h];
genericTarget = finiteDiff[Function[u, f0[u]], theta, h];

{genericResult, genericTarget, Simplify[genericResult - genericTarget] === 0}
```
**Output:**
```
Out[6]= \left\{\frac{f0\left(h+\mathrm{theta}\right)-f0\left(\mathrm{theta}\right)}{h},\frac{f0\left(h+\mathrm{theta}\right)-f0\left(\mathrm{theta}\right)}{h},\mathrm{True}\right\}
```

## 2026-07-12 10:54:02
**Input:**
```wolfram
testFunctions = {Function[u, u^3], Function[u, Exp[u]], Function[u, Sin[u]]};

(* symbolic check: generic theta, h *)
symbolicChecks = Simplify[gaudinMeasureResidue[#, theta, h] - finiteDiff[#, theta, h]] & /@ testFunctions

(* numeric check: specific theta, h *)
numericChecks = Chop[
   gaudinMeasureResidue[#, 7/3, 11/5] - finiteDiff[#, 7/3, 11/5]] & /@ testFunctions
```
**Output:**
```
Out[8]= \left\{0,0,0\right\} | Out[9]= \left\{0,\frac{5}{11}\,\mathrm{e}^{7/3}\,\left(\mathrm{e}^{11/5}-1\right)-\frac{5}{11}\,\left(\mathrm{e}^{68/15}-\mathrm{e}^{7/3}\right),0\right\}
```
## 2026-07-12 10:54:20 — ✏️ EDIT CELL 13
testFunctions = {Function[u, u^3], Function[u, Exp[u]], Function[u, Sin[u]]};

(* symbolic check: generic theta, h *)
symbolicChecks = Simplify[gaudinMeasureResidue[#, theta, h] - finiteDiff[#, theta,…


## 2026-07-12 10:54:20
**Input:**
```wolfram
testFunctions = {Function[u, u^3], Function[u, Exp[u]], Function[u, Sin[u]]};

(* symbolic check: generic theta, h *)
symbolicChecks = Simplify[gaudinMeasureResidue[#, theta, h] - finiteDiff[#, theta, h]] & /@ testFunctions

(* exact-rational check at specific theta, h: Simplify to confirm the exact difference vanishes *)
exactChecks = Simplify[
   gaudinMeasureResidue[#, 7/3, 11/5] - finiteDiff[#, 7/3, 11/5]] & /@ testFunctions

(* numeric (floating-point) check at specific theta, h *)
numericChecks = Chop[
   N[gaudinMeasureResidue[#, 7/3, 11/5] - finiteDiff[#, 7/3, 11/5]]] & /@ testFunctions
```
**Output:**
```
[message] Set::write: Tag Times in {0, 0, 0} {0, ------------------- - ------------------, 0} is Protected.                                                11                    11
[message] [Set::write ×1 more suppressed — first instance shown above]
Out= {0, 0, 0}
```
## 2026-07-12 10:54:45 — ✏️ EDIT CELL 13
testFunctions = {Function[u, u^3], Function[u, Exp[u]], Function[u, Sin[u]]};

(* symbolic check: generic theta, h *)
symbolicChecks = Simplify[gaudinMeasureResidue[#, theta, h] - finiteDiff[#, theta,…


## 2026-07-12 10:54:45
**Input:**
```wolfram
testFunctions = {Function[u, u^3], Function[u, Exp[u]], Function[u, Sin[u]]};

(* symbolic check: generic theta, h *)
symbolicChecks = Simplify[gaudinMeasureResidue[#, theta, h] - finiteDiff[#, theta, h]] & /@ testFunctions;

(* exact-rational check at specific theta, h: Simplify to confirm the exact difference vanishes *)
exactChecks = Simplify[gaudinMeasureResidue[#, 7/3, 11/5] - finiteDiff[#, 7/3, 11/5]] & /@ testFunctions;

(* numeric (floating-point) check at specific theta, h *)
numericChecks = Chop[N[gaudinMeasureResidue[#, 7/3, 11/5] - finiteDiff[#, 7/3, 11/5]]] & /@ testFunctions;

{symbolicChecks, exactChecks, numericChecks}
```
**Output:**
```
Out= {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}}
```
## 2026-07-12 10:55:05 — ↕️ MOVE MARKDOWN CELL 7 → position 25
## SU(3) determinant limit: contour-integral building block

Before specializing the quantum determi

## 2026-07-12 10:55:14 — ↕️ MOVE MARKDOWN CELL 7 → position 25
Defines the measure `mu[u,theta,h]` with simple poles at $u=\theta$ and $u=\theta+h$, and `gaudinMea

## 2026-07-12 10:55:21 — ↕️ MOVE CODE CELL 7 → position 25
mu[u_, theta_, h_] := 1/((u - theta) (u - theta - h));

gaudinMeasureResidue[f_, theta_, h_] :=
  Si

## 2026-07-12 10:55:26 — ↕️ MOVE MARKDOWN CELL 7 → position 25
Checks the reduction symbolically for a generic black-box function $f_0$: since both poles are simpl

## 2026-07-12 10:55:35 — ↕️ MOVE CODE CELL 7 → position 25
finiteDiff[f_, theta_, h_] := (f[theta + h] - f[theta])/h;

genericResult = gaudinMeasureResidue[Fun

## 2026-07-12 10:55:47 — ↕️ MOVE MARKDOWN CELL 7 → position 25
Validates on concrete pole-free test functions (polynomial, exponential, trig — no test function may

## 2026-07-12 10:55:53 — ↕️ MOVE CODE CELL 7 → position 25
testFunctions = {Function[u, u^3], Function[u, Exp[u]], Function[u, Sin[u]]};

(* symbolic check: ge

## 2026-07-12 10:55:59 — ↕️ MOVE MARKDOWN CELL 7 → position 25
Both checks vanish identically: the contour integral of a pole-free $f$ against the two-point Gaudin

## 2026-07-12 10:56:12 — ▶️ RUN CELLS 20–24
Cell 20: ✓ — (no output)
Cell 22: ✓ — Out[15]= \left\{\frac{f0\left(h+\mathrm{theta}\right)-f0\left(\mathrm{theta}\right)}{h},\frac{f0\left(h+\mathrm{theta}\right)-f0\left(\mathrm{theta}\right)}{h},\mathrm{True}\right\}
Cell 24: ✓ — Out[20]= \begin{pmatrix}0 & 0 & 0\\0 & 0 & 0\\0 & 0 & 0\end{pmatrix}

## 2026-07-12 11:03:25 — ✏️ EDIT CELL 19 [batch]
Defines the measure `mu[u,theta,h]` with simple poles at $u=\theta$ and $u=\theta+h$, and `gaudinMeasureResidue[expr,theta,h]`, the sum of the residues of $\mathrm{expr}(u)\,\mu(u)$ at those two poles…

## 2026-07-12 11:03:25 — ✏️ EDIT CELL 20 [batch]
ClearAll[mu, gaudinMeasureResidue];

mu[u_, theta_, h_] := 1/((u - theta) (u - theta - h));

gaudinMeasureResidue[expr_, theta_, h_] :=
  Simplify[
   Residue[expr mu[u, theta, h], {u, theta}] +
    R…


## 2026-07-12 11:03:25
**Input:**
```wolfram
ClearAll[mu, gaudinMeasureResidue];

mu[u_, theta_, h_] := 1/((u - theta) (u - theta - h));

gaudinMeasureResidue[expr_, theta_, h_] :=
  Simplify[
   Residue[expr mu[u, theta, h], {u, theta}] +
    Residue[expr mu[u, theta, h], {u, theta + h}]];
```
**Output:**
```
Out= (no output)
```
## 2026-07-12 11:03:25 — ✏️ EDIT CELL 21 [batch]
Checks the reduction symbolically for a generic black-box expression $f_0(u)$ (i.e. `f0[u]`, plain in the variable `u`): since both poles are simple, $\mathrm{Res}_{u=\theta} = -f_0(\theta)/h$ and $\m…

## 2026-07-12 11:03:25 — ✏️ EDIT CELL 22 [batch]
ClearAll[finiteDiff];

(* substitution-based finite difference: expr is a bare expression in u, not a function head *)
finiteDiff[expr_, theta_, h_] := ((expr /. u -> theta + h) - (expr /. u -> theta)…


## 2026-07-12 11:03:26
**Input:**
```wolfram
ClearAll[finiteDiff];

(* substitution-based finite difference: expr is a bare expression in u, not a function head *)
finiteDiff[expr_, theta_, h_] := ((expr /. u -> theta + h) - (expr /. u -> theta))/h;

genericResult = gaudinMeasureResidue[f0[u], theta, h];
genericTarget = finiteDiff[f0[u], theta, h];

{genericResult, genericTarget, Simplify[genericResult - genericTarget] === 0}
```
**Output:**
```
Out= Out[28]= \left\{\frac{f0\left(h+\mathrm{theta}\right)-f0\left(\mathrm{theta}\right)}{h},\frac{f0\left(h+\mathrm{theta}\right)-f0\left(\mathrm{theta}\right)}{h},\mathrm{True}\right\}
```
## 2026-07-12 11:03:26 — ✏️ EDIT CELL 23 [batch]
Validates on concrete pole-free test expressions in `u` (polynomial, exponential, trig, and a symbolic-rate exponential $e^{k u}$ — no test expression may have a pole near $\theta,\theta+h$, unlike e.…

## 2026-07-12 11:03:26 — ✏️ EDIT CELL 24 [batch]
testExprs = {u^3, Exp[u], Sin[u], Exp[k u]};

(* symbolic check: generic theta, h *)
symbolicChecks = Simplify[gaudinMeasureResidue[#, theta, h] - finiteDiff[#, theta, h]] & /@ testExprs;

(* exact-ra…


## 2026-07-12 11:03:26
**Input:**
```wolfram
testExprs = {u^3, Exp[u], Sin[u], Exp[k u]};

(* symbolic check: generic theta, h *)
symbolicChecks = Simplify[gaudinMeasureResidue[#, theta, h] - finiteDiff[#, theta, h]] & /@ testExprs;

(* exact-rational check at specific theta, h: Simplify to confirm the exact difference vanishes *)
exactChecks = Simplify[gaudinMeasureResidue[#, 7/3, 11/5] - finiteDiff[#, 7/3, 11/5]] & /@ testExprs;

(* numeric (floating-point) check at specific theta, h; also substitute a numeric rate k for the Exp[k u] case *)
numericChecks = Chop[N[(gaudinMeasureResidue[#, 7/3, 11/5] - finiteDiff[#, 7/3, 11/5]) /. k -> 3/2]] & /@ testExprs;

{symbolicChecks, exactChecks, numericChecks}
```
**Output:**
```
Out= Out[33]= \begin{pmatrix}0 & 0 & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\end{pmatrix}
```
## 2026-07-12 11:03:26 — ✏️ EDIT CELL 25 [batch]
Both checks vanish identically: the contour integral of a pole-free expression $\mathrm{expr}(u)$ against the two-point Gaudin measure $\mu(u)$ is exactly the finite-difference operator $\mathrm{expr}…


## 2026-07-12 11:03:46
**Input:**
```wolfram
gaudinMeasureResidue[Exp[k u], theta, h]
```
**Output:**
```
Out= (E^(k*theta)*(-1 + E^(h*k)))/h
```

## 2026-07-24 05:48:52
**Input:**
```wolfram
ClearAll[Mg, Lg, fg, MChaing, applyMsg, EEg, eg, cleanupg];
Mg[{i_, j_}][ff_] := D[ff, z] KroneckerDelta[i, j] + Lg[i, j][z] . ff;
cleanupg[x_] := x //. Dot[a___, b_ + c_, d___] :> Dot[a, b, d] + Dot[a, c, d] //. Dot[a___, b_, d___] :> b Dot[a, d] /; FreeQ[b, fg | Lg] //. Dot[a___, b_*c_, d___] :> b Dot[a, c, d] /; FreeQ[b, Lg] //. Dot[a___, fg'''[z]] :> Dot[a] fg'''[z] //. Dot[a___, fg''[z]] :> Dot[a] fg''[z] //. Dot[a___, fg'[z]] :> Dot[a] fg'[z] //. Dot[a___, fg[z]] :> Dot[a] fg[z];
applyMsg = {MChaing[{a__, b__, c__, d__}] :> Mg[a][Mg[b][Mg[c][Mg[d][fg[z]]]]] //. Dot[a___, b_ + c_, d___] :> Dot[a, b, d] + Dot[a, c, d]};
EEg[II_, JJ_] := 1/4! Sum[Signature[s] Times @@ (eg[#] & /@ Transpose[{II[[s]], JJ}]), {s, Permutations[{1, 2, 3, 4}]}] MChaing[Transpose[{II, JJ}]] /. eg[{a_, b_}] :> KroneckerDelta[a, b];
detOpOld = cleanupg[ExpandAll[Sum[EEg[II, JJ], {II, Tuples[Range[4], 4]}, {JJ, Tuples[Range[4], 4]}] //. applyMsg]];
detOpNew = cleanupg[ExpandAll[Sum[Signature[s] Mg[{s[[1]], 1}][Mg[{s[[2]], 2}][Mg[{s[[3]], 3}][Mg[{s[[4]], 4}][fg[z]]]]], {s, Permutations[{1, 2, 3, 4}]}]]];
Table[ExpandAll[Coefficient[detOpOld, d] - Coefficient[detOpNew, d]], {d, {fg'''[z], fg''[z], fg'[z], fg[z]}}]
```
**Output:**
```
Out[1]= (suppressed)
Out[2]= (suppressed)
Out[3]= (suppressed)
Out[4]= (suppressed)
Out[5]= (suppressed)
Out[6]= (suppressed)
Out[7]= (suppressed)
Out[8]= {0, -1/2*Lg[1, 1][z] . Lg[2, 2][z] - Lg[1, 1][z] . Lg[3, 3][z]/2 - Lg[1, 1][z] . Lg[4, 4][z]/2 - Lg[1, 2][z] . Lg[2, 1][z]/2 - Lg[1, 3][z] . Lg[3, 1][z]/2 - Lg[1, 4][z] . Lg[4, 1][z]/2 + Lg[2, 1][z] . Lg[1, 2][z]/2 + Lg[2, 2][z] . Lg[1, 1][z]/2 - Lg[2, 2][z] . Lg[3, 3][z]/2 - Lg[2, 2][z] . Lg[4, 4][z]/2 - Lg[2, 3][z] . Lg[3, 2][z]/2 - Lg[2, 4][z] . Lg[4, 2][z]/2 + Lg[3, 1][z] . Lg[1, 3][z]/2 + Lg[3, 2][z] . Lg[2, 3][z]/2 + Lg[3, 3][z] . Lg[1, 1][z]/2 + Lg[3, 3][z] . Lg[2, 2][z]/2 - Lg[3, 3][z] . Lg[4, 4][z]/2 - Lg[3, 4][z] . Lg[4, 3][z]/2 + Lg[4, 1][z] . Lg[1, 4][z]/2 + Lg[4, 2][z] . Lg[2, 4][z]/2 + Lg[4, 3][z] . Lg[3, 4][z]/2 + Lg[4, 4][z] . Lg[1, 1][z]/2 + Lg[4, 4][z] . Lg[2, 2][z]/2 + Lg[4, 4][z] . Lg[3, 3][z]/2 + (3*Derivative[1][Lg[1, 1]][z])/2 + Derivative[1][Lg[2, 2]][z]/2 - Derivative[1][Lg[3, 3]][z]/2 - (3*Derivative[1][Lg[4, 4]][z])/2, (2*Lg[1, 1][z] . Derivative[1][Lg[2, 2]][z])/3 - Lg[1, 1][z] . Derivative[1][Lg[3, 3]][z]/3 - (4*Lg[1, 1][z] . Derivative[1][Lg[4, 4]][z])/3 - (2*Lg[1, 2][z] . Derivative[1][Lg[2, 1]][z])/3 - (2*Lg[1, 3][z] . Derivative[1][Lg[3, 1]][z])/3 - (2*Lg[1, 4][z] . Derivative[1][Lg[4, 1]][z])/3 - (2*Lg[2, 1][z] . Derivative[1][Lg[1, 2]][z])/3 + (2*Lg[2, 2][z] . Derivative[1][Lg[1, 1]][z])/3 - Lg[2, 2][z] . Derivative[1][Lg[3, 3]][z]/3 - (4*Lg[2, 2][z] . Derivative[1][Lg[4, 4]][z])/3 - (2*Lg[2, 3][z] . Derivative[1][Lg[3, 2]][z])/3 - (2*Lg[2, 4][z] . Derivative[1][Lg[4, 2]][z])/3 + Lg[3, 1][z] . Derivative[1][Lg[1, 3]][z]/3 + Lg[3, 2][z] . Derivative[1][Lg[2, 3]][z]/3 + (2*Lg[3, 3][z] . Derivative[1][Lg[1, 1]][z])/3 + (2*Lg[3, 3][z] . Derivative[1][Lg[2, 2]][z])/3 - (4*Lg[3, 3][z] . Derivative[1][Lg[4, 4]][z])/3 - (2*Lg[3, 4][z] . Derivative[1][Lg[4, 3]][z])/3 + (4*Lg[4, 1][z] . Derivative[1][Lg[1, 4]][z])/3 + (4*Lg[4, 2][z] . Derivative[1][Lg[2, 4]][z])/3 + (4*Lg[4, 3][z] . Derivative[1][Lg[3, 4]][z])/3 + (2*Lg[4, 4][z] . Derivative[1][Lg[1, 1]][z])/3 + (2*Lg[4, 4][z] . Derivative[1][Lg[2, 2]][z])/3 + (2*Lg[4, 4][z] . Derivative[1][Lg[3, 3]][z])/3 + Derivative[1][Lg[1, 1]][z] . Lg[2, 2][z]/3 + Derivative[1][Lg[1, 1]][z] . Lg[3, 3][z]/3 + Derivative[1][Lg[1, 1]][z] . Lg[4, 4][z]/3 - Derivative[1][Lg[1, 2]][z] . Lg[2, 1][z]/3 - Derivative[1][Lg[1, 3]][z] . Lg[3, 1][z]/3 - Derivative[1][Lg[1, 4]][z] . Lg[4, 1][z]/3 - Derivative[1][Lg[2, 1]][z] . Lg[1, 2][z]/3 + Derivative[1][Lg[2, 2]][z] . Lg[1, 1][z]/3 - (2*Derivative[1][Lg[2, 2]][z] . Lg[3, 3][z])/3 - (2*Derivative[1][Lg[2, 2]][z] . Lg[4, 4][z])/3 - Derivative[1][Lg[2, 3]][z] . Lg[3, 2][z]/3 - Derivative[1][Lg[2, 4]][z] . Lg[4, 2][z]/3 - Derivative[1][Lg[3, 1]][z] . Lg[1, 3][z]/3 + (2*Derivative[1][Lg[3, 2]][z] . Lg[2, 3][z])/3 + Derivative[1][Lg[3, 3]][z] . Lg[1, 1][z]/3 + Derivative[1][Lg[3, 3]][z] . Lg[2, 2][z]/3 - (5*Derivative[1][Lg[3, 3]][z] . Lg[4, 4][z])/3 - Derivative[1][Lg[3, 4]][z] . Lg[4, 3][z]/3 - Derivative[1][Lg[4, 1]][z] . Lg[1, 4][z]/3 + (2*Derivative[1][Lg[4, 2]][z] . Lg[2, 4][z])/3 + (5*Derivative[1][Lg[4, 3]][z] . Lg[3, 4][z])/3 + Derivative[1][Lg[4, 4]][z] . Lg[1, 1][z]/3 + Derivative[1][Lg[4, 4]][z] . Lg[2, 2][z]/3 + Derivative[1][Lg[4, 4]][z] . Lg[3, 3][z]/3 - (5*Lg[1, 1][z] . Lg[2, 2][z] . Lg[3, 3][z])/6 - (5*Lg[1, 1][z] . Lg[2, 2][z] . Lg[4, 4][z])/6 - Lg[1, 1][z] . Lg[2, 3][z] . Lg[3, 2][z]/6 - Lg[1, 1][z] . Lg[2, 4][z] . Lg[4, 2][z]/6 + (5*Lg[1, 1][z] . Lg[3, 2][z] . Lg[2, 3][z])/6 + Lg[1, 1][z] . Lg[3, 3][z] . Lg[2, 2][z]/6 - (5*Lg[1, 1][z] . Lg[3, 3][z] . Lg[4, 4][z])/6 - Lg[1, 1][z] . Lg[3, 4][z] . Lg[4, 3][z]/6 + (5*Lg[1, 1][z] . Lg[4, 2][z] . Lg[2, 4][z])/6 + (5*Lg[1, 1][z] . Lg[4, 3][z] . Lg[3, 4][z])/6 + Lg[1, 1][z] . Lg[4, 4][z] . Lg[2, 2][z]/6 + Lg[1, 1][z] . Lg[4, 4][z] . Lg[3, 3][z]/6 - Lg[1, 2][z] . Lg[2, 1][z] . Lg[3, 3][z]/6 - Lg[1, 2][z] . Lg[2, 1][z] . Lg[4, 4][z]/6 + Lg[1, 2][z] . Lg[2, 3][z] . Lg[3, 1][z]/6 + Lg[1, 2][z] . Lg[2, 4][z] . Lg[4, 1][z]/6 + Lg[1, 2][z] . Lg[3, 1][z] . Lg[2, 3][z]/6 - Lg[1, 2][z] . Lg[3, 3][z] . Lg[2, 1][z]/6 + Lg[1, 2][z] . Lg[4, 1][z] . Lg[2, 4][z]/6 - Lg[1, 2][z] . Lg[4, 4][z] . Lg[2, 1][z]/6 + Lg[1, 3][z] . Lg[2, 1][z] . Lg[3, 2][z]/6 - Lg[1, 3][z] . Lg[2, 2][z] 
[output truncated — 73484 chars total]
```

## 2026-07-24 05:49:33
**Input:**
```wolfram
coeffs[x_] := Table[ExpandAll[Coefficient[x, d]], {d, {fg'''[z], fg''[z], fg'[z], fg[z]}}];
ref = coeffs[detOpOld];
candB = cleanupg[ExpandAll[Sum[Signature[s] Mg[{1, s[[1]]}][Mg[{2, s[[2]]}][Mg[{3, s[[3]]}][Mg[{4, s[[4]]}][fg[z]]]]], {s, Permutations[{1, 2, 3, 4}]}]]];
candC = cleanupg[ExpandAll[Sum[Signature[s] Mg[{s[[4]], 4}][Mg[{s[[3]], 3}][Mg[{s[[2]], 2}][Mg[{s[[1]], 1}][fg[z]]]]], {s, Permutations[{1, 2, 3, 4}]}]]];
candD = cleanupg[ExpandAll[Sum[Signature[s] Mg[{4, s[[4]]}][Mg[{3, s[[3]]}][Mg[{2, s[[2]]}][Mg[{1, s[[1]]}][fg[z]]]]], {s, Permutations[{1, 2, 3, 4}]}]]];
{"rowFwd(B)" -> (coeffs[candB] === ref), "colRev(C)" -> (coeffs[candC] === ref), "rowRev(D)" -> (coeffs[candD] === ref)}
```
**Output:**
```
Out[1]= (suppressed)
Out[2]= (suppressed)
Out[3]= (suppressed)
Out[4]= (suppressed)
Out[5]= (suppressed)
Out[6]= {"rowFwd(B)" -> False, "colRev(C)" -> False, "rowRev(D)" -> False}
```

## 2026-07-24 05:50:27
**Input:**
```wolfram
detOpFast = 1/4! cleanupg[ExpandAll[Sum[
     Signature[iperm] Signature[jperm] Mg[{iperm[[1]], jperm[[1]]}][Mg[{iperm[[2]], jperm[[2]]}][Mg[{iperm[[3]], jperm[[3]]}][Mg[{iperm[[4]], jperm[[4]]}][fg[z]]]]],
     {iperm, Permutations[{1, 2, 3, 4}]}, {jperm, Permutations[{1, 2, 3, 4}]}]]];
coeffs[detOpFast] === ref
```
**Output:**
```
Out= True
```

## 2026-07-24 05:50:59
**Input:**
```wolfram
tOld = First@AbsoluteTiming[cleanupg[ExpandAll[Sum[EEg[II, JJ], {II, Tuples[Range[4], 4]}, {JJ, Tuples[Range[4], 4]}] //. applyMsg]]];
tNew = First@AbsoluteTiming[1/4! cleanupg[ExpandAll[Sum[Signature[iperm] Signature[jperm] Mg[{iperm[[1]], jperm[[1]]}][Mg[{iperm[[2]], jperm[[2]]}][Mg[{iperm[[3]], jperm[[3]]}][Mg[{iperm[[4]], jperm[[4]]}][fg[z]]]]], {iperm, Permutations[{1, 2, 3, 4}]}, {jperm, Permutations[{1, 2, 3, 4}]}]]]];
{"old (s)" -> tOld, "new (s)" -> tNew, "speedup" -> tOld/tNew}
```
**Output:**
```
Out= {"old (s)" -> 9.90374, "new (s)" -> 0.321433, "speedup" -> 30.811211045536705}
```
## 2026-07-24 05:52:18 — ✏️ EDIT CELL 15
(* Symmetrized quantum determinant det_q(\[PartialD]_z + L). The \[Epsilon]-\[Epsilon] sum over all
   index tuples collapses: only permutations I,J\[Element]S_4 survive, with weight
   Signature[I] S…

## 2026-07-24 05:52:34 — ✏️ EDIT CELL 10
The determinant is assembled by directly nesting `M[{i,j}]` into the operator composition $M_{i_1 j_1}\circ M_{i_2 j_2}\circ M_{i_3 j_3}\circ M_{i_4 j_4}$ acting on $f(z)$, summed over row/column inde…

## 2026-07-24 05:52:37 — 🗑️ DELETE CODE CELL 11
applyMs = {
   MChain[{a__, b__, c__, d__}] :> M[a][M[b][M[c][M[d][f[z]]]]] //.
     Dot[a___, b_ + c_, d___] :> Dot[a, b, d] + Dot[a, c, d]
  };

## 2026-07-24 05:52:37 — 🗑️ DELETE MARKDOWN CELL 12
The $\epsilon$-$\epsilon$ contraction term: sums $\mathrm{sign}(\sigma)$ over permutations $\sigma$ matching $I_\sigma$ to $J$ (zero unless $I$ is a permutation of $J$), weighting the ordered operator…

## 2026-07-24 05:52:37 — 🗑️ DELETE CODE CELL 13
EE[II_, JJ_] :=
  1/4! Sum[
     Signature[\[Sigma]] Times @@ (e[#] & /@ Transpose[{II[[\[Sigma]]], JJ}]),
     {\[Sigma], Permutations[{1, 2, 3, 4}]}
    ] MChain[Transpose[{II, JJ}]] /. e[{a_, b_}] …

## 2026-07-24 05:52:37 — 🗑️ DELETE MARKDOWN CELL 14
Sums `EE[I,J]` over all index tuples $I,J\in\{1,\dots,4\}^4$, i.e. the identity $\det A = \frac{1}{N!}\sum_{I,J} \epsilon_I \epsilon_J\, A_{i_1 j_1}\cdots A_{i_N j_N}$ applied noncommutatively (scalar…

## 2026-07-24 05:53:17 — ▶️ RUN CELL 9
**In [9]:** `M[{i_, j_}][f_] := D[f, z] KroneckerDelta[i, j] + L[i, j][z] . f;`
**Out:** (no output)

## 2026-07-24 05:53:20 — ▶️ RUN CELL 11
**In [11]:** `(* Symmetrized quantum determinant det_q(\[PartialD]_z + L). The \[Epsilon]-\[Epsilon] sum over all↵   index tuples collapses: only permutations I,J\[Element]S_4 survive, with weight↵   Signature[I] S`
**Out:** (no output)


## 2026-07-24 05:53:24
**Input:**
```wolfram
refFromScratch = coeffs[detOpFast] /. Lg -> L /. fg -> f;
notebookH = Table[ExpandAll[Coefficient[detOp, d]], {d, {f'''[z], f''[z], f'[z], f[z]}}];
ExpandAll[notebookH - refFromScratch]
```
**Output:**
```
Out= {0, 0, 0, 0}
```

## 2026-07-24 05:53:35
**Input:**
```wolfram
ClearAll[Mg, Lg, fg, MChaing, applyMsg, EEg, eg, cleanupg, detOpOld, detOpNew, detOpFast, candB, candC, candD, coeffs, ref, refFromScratch, notebookH, tOld, tNew]; "scratch cleared"
```
**Output:**
```
Out= scratch cleared
```

## 2026-07-24 06:52:16
**Input:**
```wolfram
ClearAll[Mg, Lg, fg];
Mg[{i_, j_}][ff_] := D[ff, z] KroneckerDelta[i, j] + Lg[i, j][z] . ff;
raw = 1/4! Sum[Signature[ip] Signature[jp] Mg[{ip[[1]], jp[[1]]}][Mg[{ip[[2]], jp[[2]]}][Mg[{ip[[3]], jp[[3]]}][Mg[{ip[[4]], jp[[4]]}][fg[z]]]]], {ip, Permutations[Range[4]]}, {jp, Permutations[Range[4]]}];
old7[e_] := e //. Dot[a___, b_ + c_, d___] :> Dot[a, b, d] + Dot[a, c, d] //. Dot[a___, b_, d___] :> b Dot[a, d] /; FreeQ[b, fg | Lg] //. Dot[a___, b_*c_, d___] :> b Dot[a, c, d] /; FreeQ[b, Lg] //. Dot[a___, fg'''[z]] :> Dot[a] fg'''[z] //. Dot[a___, fg''[z]] :> Dot[a] fg''[z] //. Dot[a___, fg'[z]] :> Dot[a] fg'[z] //. Dot[a___, fg[z]] :> Dot[a] fg[z];
new2[e_] := Expand[e] //. Dot[x___, p_Plus, y___] :> Total[Dot[x, #, y] & /@ (List @@ p)] //. Dot[x___, s_, y___] /; FreeQ[s, Lg] :> s Dot[x, y];
newExpandOnly[e_] := Expand[e] //. Dot[x___, s_, y___] /; FreeQ[s, Lg] :> s Dot[x, y];
r7 = old7[ExpandAll[raw]];
Table[ExpandAll[Coefficient[r7, d] - Coefficient[new2[raw], d]], {d, {fg'''[z], fg''[z], fg'[z], fg[z]}}] -> Table[ExpandAll[Coefficient[r7, d] - Coefficient[newExpandOnly[raw], d]], {d, {fg'''[z], fg''[z], fg'[z], fg[z]}}]
```
**Output:**
```
Out[1]= (suppressed)
Out[2]= (suppressed)
Out[3]= (suppressed)
Out[4]= (suppressed)
Out[5]= (suppressed)
Out[6]= (suppressed)
Out[7]= (suppressed)
Out[8]= {0, 0, Lg[1, 1][z] . Derivative[1][Lg[2, 2]][z]/2 + Lg[1, 1][z] . Derivative[1][Lg[3, 3]][z]/2 + Lg[1, 1][z] . Derivative[1][Lg[4, 4]][z]/2 - Lg[1, 2][z] . Derivative[1][Lg[2, 1]][z]/3 - Lg[1, 3][z] . Derivative[1][Lg[3, 1]][z]/3 - Lg[1, 4][z] . Derivative[1][Lg[4, 1]][z]/3 - Lg[2, 1][z] . Derivative[1][Lg[1, 2]][z]/3 + Lg[2, 2][z] . Derivative[1][Lg[1, 1]][z]/2 + Lg[2, 2][z] . Derivative[1][Lg[3, 3]][z]/2 + Lg[2, 2][z] . Derivative[1][Lg[4, 4]][z]/2 - Lg[2, 3][z] . Derivative[1][Lg[3, 2]][z]/3 - Lg[2, 4][z] . Derivative[1][Lg[4, 2]][z]/3 - Lg[3, 1][z] . Derivative[1][Lg[1, 3]][z]/3 - Lg[3, 2][z] . Derivative[1][Lg[2, 3]][z]/3 + Lg[3, 3][z] . Derivative[1][Lg[1, 1]][z]/2 + Lg[3, 3][z] . Derivative[1][Lg[2, 2]][z]/2 + Lg[3, 3][z] . Derivative[1][Lg[4, 4]][z]/2 - Lg[3, 4][z] . Derivative[1][Lg[4, 3]][z]/3 - Lg[4, 1][z] . Derivative[1][Lg[1, 4]][z]/3 - Lg[4, 2][z] . Derivative[1][Lg[2, 4]][z]/3 - Lg[4, 3][z] . Derivative[1][Lg[3, 4]][z]/3 + Lg[4, 4][z] . Derivative[1][Lg[1, 1]][z]/2 + Lg[4, 4][z] . Derivative[1][Lg[2, 2]][z]/2 + Lg[4, 4][z] . Derivative[1][Lg[3, 3]][z]/2, 0} -> {(3*Lg[1, 1][z])/4 + (3*Lg[2, 2][z])/4 + (3*Lg[3, 3][z])/4 + (3*Lg[4, 4][z])/4, Lg[1, 1][z] . Lg[2, 2][z]/2 + Lg[1, 1][z] . Lg[3, 3][z]/2 + Lg[1, 1][z] . Lg[4, 4][z]/2 - (5*Lg[1, 2][z] . Lg[2, 1][z])/12 - (5*Lg[1, 3][z] . Lg[3, 1][z])/12 - (5*Lg[1, 4][z] . Lg[4, 1][z])/12 - (5*Lg[2, 1][z] . Lg[1, 2][z])/12 + Lg[2, 2][z] . Lg[1, 1][z]/2 + Lg[2, 2][z] . Lg[3, 3][z]/2 + Lg[2, 2][z] . Lg[4, 4][z]/2 - (5*Lg[2, 3][z] . Lg[3, 2][z])/12 - (5*Lg[2, 4][z] . Lg[4, 2][z])/12 - (5*Lg[3, 1][z] . Lg[1, 3][z])/12 - (5*Lg[3, 2][z] . Lg[2, 3][z])/12 + Lg[3, 3][z] . Lg[1, 1][z]/2 + Lg[3, 3][z] . Lg[2, 2][z]/2 + Lg[3, 3][z] . Lg[4, 4][z]/2 - (5*Lg[3, 4][z] . Lg[4, 3][z])/12 - (5*Lg[4, 1][z] . Lg[1, 4][z])/12 - (5*Lg[4, 2][z] . Lg[2, 4][z])/12 - (5*Lg[4, 3][z] . Lg[3, 4][z])/12 + Lg[4, 4][z] . Lg[1, 1][z]/2 + Lg[4, 4][z] . Lg[2, 2][z]/2 + Lg[4, 4][z] . Lg[3, 3][z]/2 + (3*Derivative[1][Lg[1, 1]][z])/4 + (3*Derivative[1][Lg[2, 2]][z])/4 + (3*Derivative[1][Lg[3, 3]][z])/4 + (3*Derivative[1][Lg[4, 4]][z])/4, (2*Lg[1, 1][z] . Derivative[1][Lg[2, 2]][z])/3 + (2*Lg[1, 1][z] . Derivative[1][Lg[3, 3]][z])/3 + (2*Lg[1, 1][z] . Derivative[1][Lg[4, 4]][z])/3 - Lg[1, 2][z] . Derivative[1][Lg[2, 1]][z]/2 - Lg[1, 3][z] . Derivative[1][Lg[3, 1]][z]/2 - Lg[1, 4][z] . Derivative[1][Lg[4, 1]][z]/2 - Lg[2, 1][z] . Derivative[1][Lg[1, 2]][z]/2 + (2*Lg[2, 2][z] . Derivative[1][Lg[1, 1]][z])/3 + (2*Lg[2, 2][z] . Derivative[1][Lg[3, 3]][z])/3 + (2*Lg[2, 2][z] . Derivative[1][Lg[4, 4]][z])/3 - Lg[2, 3][z] . Derivative[1][Lg[3, 2]][z]/2 - Lg[2, 4][z] . Derivative[1][Lg[4, 2]][z]/2 - Lg[3, 1][z] . Derivative[1][Lg[1, 3]][z]/2 - Lg[3, 2][z] . Derivative[1][Lg[2, 3]][z]/2 + (2*Lg[3, 3][z] . Derivative[1][Lg[1, 1]][z])/3 + (2*Lg[3, 3][z] . Derivative[1][Lg[2, 2]][z])/3 + (2*Lg[3, 3][z] . Derivative[1][Lg[4, 4]][z])/3 - Lg[3, 4][z] . Derivative[1][Lg[4, 3]][z]/2 - Lg[4, 1][z] . Derivative[1][Lg[1, 4]][z]/2 - Lg[4, 2][z] . Derivative[1][Lg[2, 4]][z]/2 - Lg[4, 3][z] . Derivative[1][Lg[3, 4]][z]/2 + (2*Lg[4, 4][z] . Derivative[1][Lg[1, 1]][z])/3 + (2*Lg[4, 4][z] . Derivative[1][Lg[2, 2]][z])/3 + (2*Lg[4, 4][z] . Derivative[1][Lg[3, 3]][z])/3 + Derivative[1][Lg[1, 1]][z] . Lg[2, 2][z]/3 + Derivative[1][Lg[1, 1]][z] . Lg[3, 3][z]/3 + Derivative[1][Lg[1, 1]][z] . Lg[4, 4][z]/3 - Derivative[1][Lg[1, 2]][z] . Lg[2, 1][z]/6 - Derivative[1][Lg[1, 3]][z] . Lg[3, 1][z]/6 - Derivative[1][Lg[1, 4]][z] . Lg[4, 1][z]/6 - Derivative[1][Lg[2, 1]][z] . Lg[1, 2][z]/6 + Derivative[1][Lg[2, 2]][z] . Lg[1, 1][z]/3 + Derivative[1][Lg[2, 2]][z] . Lg[3, 3][z]/3 + Derivative[1][Lg[2, 2]][z] . Lg[4, 4][z]/3 - Derivative[1][Lg[2, 3]][z] . Lg[3, 2][z]/6 - Derivative[1][Lg[2, 4]][z] . Lg[4, 2][z]/6 - Derivative[1][Lg[3, 1]][z] . Lg[1, 3][z]/6 - Derivative[1][Lg[3, 2]][z] . Lg[2, 3][z]/6 + Derivative[1][Lg[3, 3]][z] . Lg[1, 1][z]/3 + Derivative[1][Lg[3, 3]][z] . Lg[2, 2][z]/3 + Derivative[1][Lg[3, 3]][z] . Lg[4, 4][z]/3 - Derivative[1][Lg[3, 4]][z] . Lg[4, 3][z]/6 - Derivative[1][Lg[4, 1]][z] . Lg[1, 4][z]/6 - Derivative[1][Lg[4, 2]
[output truncated — 53776 chars total]
```

## 2026-07-24 06:53:38
**Input:**
```wolfram
strip = Dot[a___, last_] :> Dot[a] last /; FreeQ[last, Lg];
r1 = Dot[a___, b_ + c_, d___] :> Dot[a, b, d] + Dot[a, c, d];
r2 = Dot[a___, b_, d___] :> b Dot[a, d] /; FreeQ[b, fg | Lg];
r3 = Dot[a___, b_*c_, d___] :> b Dot[a, c, d] /; FreeQ[b, Lg];
Bver[e_] := e //. r1 //. r2 //. r3 //. strip;
Cver[e_] := e //. r1 //. r2 //. strip;
re = ExpandAll[raw];
{"B (merge 4-7, keep r3)" -> (Table[ExpandAll[Coefficient[r7, d] - Coefficient[Bver[re], d]], {d, {fg'''[z], fg''[z], fg'[z], fg[z]}}] === {0, 0, 0, 0}),
 "C (also drop r3)" -> (Table[ExpandAll[Coefficient[r7, d] - Coefficient[Cver[re], d]], {d, {fg'''[z], fg''[z], fg'[z], fg[z]}}] === {0, 0, 0, 0})}
```
**Output:**
```
Out[1]= (suppressed)
Out[2]= (suppressed)
Out[3]= (suppressed)
Out[4]= (suppressed)
Out[5]= (suppressed)
Out[6]= (suppressed)
Out[7]= (suppressed)
Out[8]= {"B (merge 4-7, keep r3)" -> True, "C (also drop r3)" -> False}
```

## 2026-07-24 06:54:03
**Input:**
```wolfram
ExpandAll[Bver[re] - r7] === 0
```
**Output:**
```
Out= True
```
## 2026-07-24 06:54:33 — ✏️ EDIT CELL 11
(* Symmetrized quantum determinant det_q(\[PartialD]_z + L). The \[Epsilon]-\[Epsilon] sum over all
   index tuples collapses: only permutations I,J\[Element]S_4 survive, with weight
   Signature[I] S…


## 2026-07-24 06:54:34
**Input:**
```wolfram
(* Symmetrized quantum determinant det_q(\[PartialD]_z + L). The \[Epsilon]-\[Epsilon] sum over all
   index tuples collapses: only permutations I,J\[Element]S_4 survive, with weight
   Signature[I] Signature[J], so det_q = 1/4! \[Sum]_{I,J\[Element]S_4} sgn(I) sgn(J) M_{I_1 J_1}\[SmallCircle]\[CenterDot]\[SmallCircle]M_{I_4 J_4}.
   Cleanup: distribute Dot over sums, pull scalar summands and scalar factors out
   (L stays noncommutative, never crossing \[PartialD]/f), then strip the trailing
   f-derivative to the right so its coefficient is a plain ordered product of L's. *)
detOp =
  (1/4! Sum[
       Signature[iperm] Signature[jperm]
        M[{iperm[[1]], jperm[[1]]}][M[{iperm[[2]], jperm[[2]]}][
          M[{iperm[[3]], jperm[[3]]}][M[{iperm[[4]], jperm[[4]]}][f[z]]]]],
       {iperm, Permutations[Range[4]]}, {jperm, Permutations[Range[4]]}] // ExpandAll) //.
    Dot[a___, b_ + c_, d___] :> Dot[a, b, d] + Dot[a, c, d] //.
    Dot[a___, b_, d___] :> b Dot[a, d] /; FreeQ[b, f | L] //.
    Dot[a___, b_*c_, d___] :> b Dot[a, c, d] /; FreeQ[b, L] //.
    Dot[a___, last_] :> Dot[a] last /; FreeQ[last, L]
```
**Output:**
```
Out= M[{1, 1}][M[{2, 2}][M[{3, 3}][M[{4, 4}][f[z]]]]]/24 - M[{1, 1}][M[{2, 2}][M[{3, 4}][M[{4, 3}][f[z]]]]]/24 - M[{1, 1}][M[{2, 2}][M[{4, 3}][M[{3, 4}][f[z]]]]]/24 + M[{1, 1}][M[{2, 2}][M[{4, 4}][M[{3, 3}][f[z]]]]]/24 - M[{1, 1}][M[{2, 3}][M[{3, 2}][M[{4, 4}][f[z]]]]]/24 + M[{1, 1}][M[{2, 3}][M[{3, 4}][M[{4, 2}][f[z]]]]]/24 + M[{1, 1}][M[{2, 3}][M[{4, 2}][M[{3, 4}][f[z]]]]]/24 - M[{1, 1}][M[{2, 3}][M[{4, 4}][M[{3, 2}][f[z]]]]]/24 + M[{1, 1}][M[{2, 4}][M[{3, 2}][M[{4, 3}][f[z]]]]]/24 - M[{1, 1}][M[{2, 4}][M[{3, 3}][M[{4, 2}][f[z]]]]]/24 - M[{1, 1}][M[{2, 4}][M[{4, 2}][M[{3, 3}][f[z]]]]]/24 + M[{1, 1}][M[{2, 4}][M[{4, 3}][M[{3, 2}][f[z]]]]]/24 - M[{1, 1}][M[{3, 2}][M[{2, 3}][M[{4, 4}][f[z]]]]]/24 + M[{1, 1}][M[{3, 2}][M[{2, 4}][M[{4, 3}][f[z]]]]]/24 + M[{1, 1}][M[{3, 2}][M[{4, 3}][M[{2, 4}][f[z]]]]]/24 - M[{1, 1}][M[{3, 2}][M[{4, 4}][M[{2, 3}][f[z]]]]]/24 + M[{1, 1}][M[{3, 3}][M[{2, 2}][M[{4, 4}][f[z]]]]]/24 - M[{1, 1}][M[{3, 3}][M[{2, 4}][M[{4, 2}][f[z]]]]]/24 - M[{1, 1}][M[{3, 3}][M[{4, 2}][M[{2, 4}][f[z]]]]]/24 + M[{1, 1}][M[{3, 3}][M[{4, 4}][M[{2, 2}][f[z]]]]]/24 - M[{1, 1}][M[{3, 4}][M[{2, 2}][M[{4, 3}][f[z]]]]]/24 + M[{1, 1}][M[{3, 4}][M[{2, 3}][M[{4, 2}][f[z]]]]]/24 + M[{1, 1}][M[{3, 4}][M[{4, 2}][M[{2, 3}][f[z]]]]]/24 - M[{1, 1}][M[{3, 4}][M[{4, 3}][M[{2, 2}][f[z]]]]]/24 + M[{1, 1}][M[{4, 2}][M[{2, 3}][M[{3, 4}][f[z]]]]]/24 - M[{1, 1}][M[{4, 2}][M[{2, 4}][M[{3, 3}][f[z]]]]]/24 - M[{1, 1}][M[{4, 2}][M[{3, 3}][M[{2, 4}][f[z]]]]]/24 + M[{1, 1}][M[{4, 2}][M[{3, 4}][M[{2, 3}][f[z]]]]]/24 - M[{1, 1}][M[{4, 3}][M[{2, 2}][M[{3, 4}][f[z]]]]]/24 + M[{1, 1}][M[{4, 3}][M[{2, 4}][M[{3, 2}][f[z]]]]]/24 + M[{1, 1}][M[{4, 3}][M[{3, 2}][M[{2, 4}][f[z]]]]]/24 - M[{1, 1}][M[{4, 3}][M[{3, 4}][M[{2, 2}][f[z]]]]]/24 + M[{1, 1}][M[{4, 4}][M[{2, 2}][M[{3, 3}][f[z]]]]]/24 - M[{1, 1}][M[{4, 4}][M[{2, 3}][M[{3, 2}][f[z]]]]]/24 - M[{1, 1}][M[{4, 4}][M[{3, 2}][M[{2, 3}][f[z]]]]]/24 + M[{1, 1}][M[{4, 4}][M[{3, 3}][M[{2, 2}][f[z]]]]]/24 - M[{1, 2}][M[{2, 1}][M[{3, 3}][M[{4, 4}][f[z]]]]]/24 + M[{1, 2}][M[{2, 1}][M[{3, 4}][M[{4, 3}][f[z]]]]]/24 + M[{1, 2}][M[{2, 1}][M[{4, 3}][M[{3, 4}][f[z]]]]]/24 - M[{1, 2}][M[{2, 1}][M[{4, 4}][M[{3, 3}][f[z]]]]]/24 + M[{1, 2}][M[{2, 3}][M[{3, 1}][M[{4, 4}][f[z]]]]]/24 - M[{1, 2}][M[{2, 3}][M[{3, 4}][M[{4, 1}][f[z]]]]]/24 - M[{1, 2}][M[{2, 3}][M[{4, 1}][M[{3, 4}][f[z]]]]]/24 + M[{1, 2}][M[{2, 3}][M[{4, 4}][M[{3, 1}][f[z]]]]]/24 - M[{1, 2}][M[{2, 4}][M[{3, 1}][M[{4, 3}][f[z]]]]]/24 + M[{1, 2}][M[{2, 4}][M[{3, 3}][M[{4, 1}][f[z]]]]]/24 + M[{1, 2}][M[{2, 4}][M[{4, 1}][M[{3, 3}][f[z]]]]]/24 - M[{1, 2}][M[{2, 4}][M[{4, 3}][M[{3, 1}][f[z]]]]]/24 + M[{1, 2}][M[{3, 1}][M[{2, 3}][M[{4, 4}][f[z]]]]]/24 - M[{1, 2}][M[{3, 1}][M[{2, 4}][M[{4, 3}][f[z]]]]]/24 - M[{1, 2}][M[{3, 1}][M[{4, 3}][M[{2, 4}][f[z]]]]]/24 + M[{1, 2}][M[{3, 1}][M[{4, 4}][M[{2, 3}][f[z]]]]]/24 - M[{1, 2}][M[{3, 3}][M[{2, 1}][M[{4, 4}][f[z]]]]]/24 + M[{1, 2}][M[{3, 3}][M[{2, 4}][M[{4, 1}][f[z]]]]]/24 + M[{1, 2}][M[{3, 3}][M[{4, 1}][M[{2, 4}][f[z]]]]]/24 - M[{1, 2}][M[{3, 3}][M[{4, 4}][M[{2, 1}][f[z]]]]]/24 + M[{1, 2}][M[{3, 4}][M[{2, 1}][M[{4, 3}][f[z]]]]]/24 - M[{1, 2}][M[{3, 4}][M[{2, 3}][M[{4, 1}][f[z]]]]]/24 - M[{1, 2}][M[{3, 4}][M[{4, 1}][M[{2, 3}][f[z]]]]]/24 + M[{1, 2}][M[{3, 4}][M[{4, 3}][M[{2, 1}][f[z]]]]]/24 - M[{1, 2}][M[{4, 1}][M[{2, 3}][M[{3, 4}][f[z]]]]]/24 + M[{1, 2}][M[{4, 1}][M[{2, 4}][M[{3, 3}][f[z]]]]]/24 + M[{1, 2}][M[{4, 1}][M[{3, 3}][M[{2, 4}][f[z]]]]]/24 - M[{1, 2}][M[{4, 1}][M[{3, 4}][M[{2, 3}][f[z]]]]]/24 + M[{1, 2}][M[{4, 3}][M[{2, 1}][M[{3, 4}][f[z]]]]]/24 - M[{1, 2}][M[{4, 3}][M[{2, 4}][M[{3, 1}][f[z]]]]]/24 - M[{1, 2}][M[{4, 3}][M[{3, 1}][M[{2, 4}][f[z]]]]]/24 + M[{1, 2}][M[{4, 3}][M[{3, 4}][M[{2, 1}][f[z]]]]]/24 - M[{1, 2}][M[{4, 4}][M[{2, 1}][M[{3, 3}][f[z]]]]]/24 + M[{1, 2}][M[{4, 4}][M[{2, 3}][M[{3, 1}][f[z]]]]]/24 + M[{1, 2}][M[{4, 4}][M[{3, 1}][M[{2, 3}][f[z]]]]]/24 - M[{1, 2}][M[{4, 4}][M[{3, 3}][M[{2, 1}][f[z]]]]]/24 + M[{1, 3}][M[{2, 1}][M[{3, 2}][M[{4, 4}][f[z]]]]]/24 - M[{1, 3}][M[{2, 1}][M[{3, 4}][M[{4, 2}][f[z]]]]]/24 - M[{1, 3}][M[{2, 1}][M[{4, 2}][M[{3, 4}][f[z]]]]]/24 + M[{1, 3}][M[{2, 1}][M[{4, 4}][M[{3, 2}][f[z]]]
[output truncated — 31101 chars total]
```
## 2026-07-24 06:54:47 — ▶️ RUN CELL 9
**In [9]:** `M[{i_, j_}][f_] := D[f, z] KroneckerDelta[i, j] + L[i, j][z] . f;`
**Out:** (no output)

## 2026-07-24 06:54:52 — ▶️ RUN CELL 11
**In [11]:** `(* Symmetrized quantum determinant det_q(\[PartialD]_z + L). The \[Epsilon]-\[Epsilon] sum over all↵   index tuples collapses: only permutations I,J\[Element]S_4 survive, with weight↵   Signature[I] S`
**Out:** (no output)


## 2026-07-24 06:54:55
**Input:**
```wolfram
Table[ExpandAll[Coefficient[detOp, d] - (Coefficient[r7, d] /. Lg -> L /. fg -> f)], {d, {f'''[z], f''[z], f'[z], f[z]}}]
```
**Output:**
```
Out= {L[1, 1][z] + L[2, 2][z] + L[3, 3][z] + L[4, 4][z], L[1, 1][z] . L[2, 2][z]/2 + L[1, 1][z] . L[3, 3][z]/2 + L[1, 1][z] . L[4, 4][z]/2 - L[1, 2][z] . L[2, 1][z]/2 - L[1, 3][z] . L[3, 1][z]/2 - L[1, 4][z] . L[4, 1][z]/2 - L[2, 1][z] . L[1, 2][z]/2 + L[2, 2][z] . L[1, 1][z]/2 + L[2, 2][z] . L[3, 3][z]/2 + L[2, 2][z] . L[4, 4][z]/2 - L[2, 3][z] . L[3, 2][z]/2 - L[2, 4][z] . L[4, 2][z]/2 - L[3, 1][z] . L[1, 3][z]/2 - L[3, 2][z] . L[2, 3][z]/2 + L[3, 3][z] . L[1, 1][z]/2 + L[3, 3][z] . L[2, 2][z]/2 + L[3, 3][z] . L[4, 4][z]/2 - L[3, 4][z] . L[4, 3][z]/2 - L[4, 1][z] . L[1, 4][z]/2 - L[4, 2][z] . L[2, 4][z]/2 - L[4, 3][z] . L[3, 4][z]/2 + L[4, 4][z] . L[1, 1][z]/2 + L[4, 4][z] . L[2, 2][z]/2 + L[4, 4][z] . L[3, 3][z]/2 + (3*Derivative[1][L[1, 1]][z])/2 + (3*Derivative[1][L[2, 2]][z])/2 + (3*Derivative[1][L[3, 3]][z])/2 + (3*Derivative[1][L[4, 4]][z])/2, (2*L[1, 1][z] . Derivative[1][L[2, 2]][z])/3 + (2*L[1, 1][z] . Derivative[1][L[3, 3]][z])/3 + (2*L[1, 1][z] . Derivative[1][L[4, 4]][z])/3 - (2*L[1, 2][z] . Derivative[1][L[2, 1]][z])/3 - (2*L[1, 3][z] . Derivative[1][L[3, 1]][z])/3 - (2*L[1, 4][z] . Derivative[1][L[4, 1]][z])/3 - (2*L[2, 1][z] . Derivative[1][L[1, 2]][z])/3 + (2*L[2, 2][z] . Derivative[1][L[1, 1]][z])/3 + (2*L[2, 2][z] . Derivative[1][L[3, 3]][z])/3 + (2*L[2, 2][z] . Derivative[1][L[4, 4]][z])/3 - (2*L[2, 3][z] . Derivative[1][L[3, 2]][z])/3 - (2*L[2, 4][z] . Derivative[1][L[4, 2]][z])/3 - (2*L[3, 1][z] . Derivative[1][L[1, 3]][z])/3 - (2*L[3, 2][z] . Derivative[1][L[2, 3]][z])/3 + (2*L[3, 3][z] . Derivative[1][L[1, 1]][z])/3 + (2*L[3, 3][z] . Derivative[1][L[2, 2]][z])/3 + (2*L[3, 3][z] . Derivative[1][L[4, 4]][z])/3 - (2*L[3, 4][z] . Derivative[1][L[4, 3]][z])/3 - (2*L[4, 1][z] . Derivative[1][L[1, 4]][z])/3 - (2*L[4, 2][z] . Derivative[1][L[2, 4]][z])/3 - (2*L[4, 3][z] . Derivative[1][L[3, 4]][z])/3 + (2*L[4, 4][z] . Derivative[1][L[1, 1]][z])/3 + (2*L[4, 4][z] . Derivative[1][L[2, 2]][z])/3 + (2*L[4, 4][z] . Derivative[1][L[3, 3]][z])/3 + Derivative[1][L[1, 1]][z] . L[2, 2][z]/3 + Derivative[1][L[1, 1]][z] . L[3, 3][z]/3 + Derivative[1][L[1, 1]][z] . L[4, 4][z]/3 - Derivative[1][L[1, 2]][z] . L[2, 1][z]/3 - Derivative[1][L[1, 3]][z] . L[3, 1][z]/3 - Derivative[1][L[1, 4]][z] . L[4, 1][z]/3 - Derivative[1][L[2, 1]][z] . L[1, 2][z]/3 + Derivative[1][L[2, 2]][z] . L[1, 1][z]/3 + Derivative[1][L[2, 2]][z] . L[3, 3][z]/3 + Derivative[1][L[2, 2]][z] . L[4, 4][z]/3 - Derivative[1][L[2, 3]][z] . L[3, 2][z]/3 - Derivative[1][L[2, 4]][z] . L[4, 2][z]/3 - Derivative[1][L[3, 1]][z] . L[1, 3][z]/3 - Derivative[1][L[3, 2]][z] . L[2, 3][z]/3 + Derivative[1][L[3, 3]][z] . L[1, 1][z]/3 + Derivative[1][L[3, 3]][z] . L[2, 2][z]/3 + Derivative[1][L[3, 3]][z] . L[4, 4][z]/3 - Derivative[1][L[3, 4]][z] . L[4, 3][z]/3 - Derivative[1][L[4, 1]][z] . L[1, 4][z]/3 - Derivative[1][L[4, 2]][z] . L[2, 4][z]/3 - Derivative[1][L[4, 3]][z] . L[3, 4][z]/3 + Derivative[1][L[4, 4]][z] . L[1, 1][z]/3 + Derivative[1][L[4, 4]][z] . L[2, 2][z]/3 + Derivative[1][L[4, 4]][z] . L[3, 3][z]/3 + L[1, 1][z] . L[2, 2][z] . L[3, 3][z]/6 + L[1, 1][z] . L[2, 2][z] . L[4, 4][z]/6 - L[1, 1][z] . L[2, 3][z] . L[3, 2][z]/6 - L[1, 1][z] . L[2, 4][z] . L[4, 2][z]/6 - L[1, 1][z] . L[3, 2][z] . L[2, 3][z]/6 + L[1, 1][z] . L[3, 3][z] . L[2, 2][z]/6 + L[1, 1][z] . L[3, 3][z] . L[4, 4][z]/6 - L[1, 1][z] . L[3, 4][z] . L[4, 3][z]/6 - L[1, 1][z] . L[4, 2][z] . L[2, 4][z]/6 - L[1, 1][z] . L[4, 3][z] . L[3, 4][z]/6 + L[1, 1][z] . L[4, 4][z] . L[2, 2][z]/6 + L[1, 1][z] . L[4, 4][z] . L[3, 3][z]/6 - L[1, 2][z] . L[2, 1][z] . L[3, 3][z]/6 - L[1, 2][z] . L[2, 1][z] . L[4, 4][z]/6 + L[1, 2][z] . L[2, 3][z] . L[3, 1][z]/6 + L[1, 2][z] . L[2, 4][z] . L[4, 1][z]/6 + L[1, 2][z] . L[3, 1][z] . L[2, 3][z]/6 - L[1, 2][z] . L[3, 3][z] . L[2, 1][z]/6 + L[1, 2][z] . L[4, 1][z] . L[2, 4][z]/6 - L[1, 2][z] . L[4, 4][z] . L[2, 1][z]/6 + L[1, 3][z] . L[2, 1][z] . L[3, 2][z]/6 - L[1, 3][z] . L[2, 2][z] . L[3, 1][z]/6 - L[1, 3][z] . L[3, 1][z] . L[2, 2][z]/6 - L[1, 3][z] . L[3, 1][z] . L[4, 4][z]/6 + L[1, 3][z] . L[3, 2][z] . L[2, 1][z]/6 + L[1, 3][z] . L[3, 4][z] . L[4, 1][z]/6 + L[1, 3][z] 
[output truncated — 68792 chars total]
```

## 2026-07-24 06:56:08
**Input:**
```wolfram
detOp7 = (1/4! Sum[Signature[ip] Signature[jp] M[{ip[[1]], jp[[1]]}][M[{ip[[2]], jp[[2]]}][M[{ip[[3]], jp[[3]]}][M[{ip[[4]], jp[[4]]}][f[z]]]]], {ip, Permutations[Range[4]]}, {jp, Permutations[Range[4]]}] // ExpandAll) //.
   Dot[a___, b_ + c_, d___] :> Dot[a, b, d] + Dot[a, c, d] //.
   Dot[a___, b_, d___] :> b Dot[a, d] /; FreeQ[b, f | L] //.
   Dot[a___, b_*c_, d___] :> b Dot[a, c, d] /; FreeQ[b, L] //.
   Dot[a___, f'''[z]] :> Dot[a] f'''[z] //. Dot[a___, f''[z]] :> Dot[a] f''[z] //.
   Dot[a___, f'[z]] :> Dot[a] f'[z] //. Dot[a___, f[z]] :> Dot[a] f[z];
ExpandAll[detOp - detOp7]
```
**Output:**
```
Out= 0
```

## 2026-07-24 06:56:18
**Input:**
```wolfram
ClearAll[detOp7, Mg, Lg, fg, raw, old7, new2, newExpandOnly, r7, strip, r1, r2, r3, Bver, Cver, re]; "cleared"
```
**Output:**
```
Out= cleared
```
