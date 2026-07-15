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
