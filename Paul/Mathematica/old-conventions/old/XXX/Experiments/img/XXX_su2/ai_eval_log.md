# AI Action Log — XXX_su2

## 2026-07-14 10:43:54 — 📥 BULK INSERT 1 CELL(S) at positions 13–13
1. [code] (*confirming c[λ][s]=Binomial[λ,s]ζ^s solves the defining equation exactly, one step further to λ=6*


## 2026-07-14 10:43:55
**Input:**
```wolfram
(*confirming c[λ][s]=Binomial[λ,s]ζ^s solves the defining equation exactly, one step further to λ=6*)
Table[x[λ][{0}]+Sum[x[λ][{s}]Binomial[λ,s]ζ^s,{s,1,λ}]-x[λ][{0}].MatrixExp[Ee[λ][2,1]ζ],{λ,1,6}]//Simplify//Flatten//Union
```
**Output:**
```
Out[258]= \left\{0\right\}
```
## 2026-07-14 10:44:06 — 📥 BULK INSERT 1 CELL(S) at positions 14–14
1. [markdown] (*Closed-form result: the global-rotation coefficients solving↵x[λ][{0}]+Sum[x[λ][{s}] c[λ][s],{s,1,

## 2026-07-14 12:38:29 — ▶️ RUN CELLS 1–14
Cell 2: ✓ — (no output)
Cell 3: ✓ — (no output)
Cell 4: ✓ — (no output)
Cell 5: ✓ — (no output)
Cell 7: ✓ — (no output)
Cell 8: ✓ — (no output)
Cell 9: ✓ — (no output)
Cell 10: ✓ — (no output)
Cell 11: ✓ — Out[25]= \left\{0\right\}
Cell 12: ✓ — Out[27]= \left\{c\left(1\right)\left(1\right)\to \zeta ,c\left(2\right)\left(1\right)\to 2\,\zeta ,c\left(2\right)\left(2\right)\to \zeta ^2,c\left(3\right)\left(1\right)\to 3\,\zeta ,c\left(3\right)\left(2\right)\to 3\,\zeta ^2,c\left(3\right)\left(3\right)\to \zeta ^3,c\left(4\right)\left(1\right)\to 4\,\zeta ,c\left(4\right)\left(2\right)\to 6\,\zeta ^2,c\left(4\right)\left(3\right)\to 4\,\zeta… [675 chars total]
Cell 13: ✓ — Out[28]= \left\{0\right\}


## 2026-07-14 12:38:36
**Input:**
```wolfram
Table[
  Module[{lhs, ratios},
    lhs = x[λ][{M}].t[λ][1,1][u]//Simplify;
    ratios = Union[Simplify[lhs/x[λ][{M}]]];
    {λ, M, ratios}
  ],
  {λ,1,3},{M,0,λ}
]
```
**Output:**
```
[message] Power::infy: Infinite expression - encountered.                                  0
[message] [Power::infy ×2 more suppressed — first instance shown above]
[message] General::stop: Further output of Power::infy will be suppressed during this calculation.
[message] Infinity::indet: Indeterminate expression 0 ComplexInfinity encountered.
[message] [Infinity::indet ×2 more suppressed — first instance shown above]
[message] [General::stop ×1 more suppressed — first instance shown above]
Out= {{{1, 0, {{ComplexInfinity, (u - θ)*χ1}}}, {1, 1, {{-((h - u + θ)*χ1), ComplexInfinity}}}}, {{2, 0, {{Indeterminate, ComplexInfinity, (u - θ)*χ1}}}, {2, 1, {{ComplexInfinity, -((h - u + θ)*χ1), ComplexInfinity}}}, {2, 2, {{-((2*h - u + θ)*χ1), ComplexInfinity, Indeterminate}}}}, {{3, 0, {{Indeterminate, Indeterminate, ComplexInfinity, (u - θ)*χ1}}}, {3, 1, {{Indeterminate, ComplexInfinity, -((h - u + θ)*χ1), ComplexInfinity}}}, {3, 2, {{ComplexInfinity, -((2*h - u + θ)*χ1), ComplexInfinity, Indeterminate}}}, {3, 3, {{-((3*h - u + θ)*χ1), ComplexInfinity, Indeterminate, Indeterminate}}}}}
```

## 2026-07-14 12:38:41
**Input:**
```wolfram
Table[
  Module[{xv, lhs, nz, ratios},
    xv = x[λ][{M}];
    lhs = xv.t[λ][1,1][u]//Simplify;
    nz = Select[Range[Length[xv[[1]]]], xv[[1,#]]=!=0 &];
    ratios = Union[Simplify[lhs[[1,#]]/xv[[1,#]] & /@ nz]];
    {λ, M, ratios}
  ],
  {λ,1,3},{M,0,λ}
]
```
**Output:**
```
Out= {{{1, 0, {(u - θ)*χ1}}, {1, 1, {-((h - u + θ)*χ1)}}}, {{2, 0, {(u - θ)*χ1}}, {2, 1, {-((h - u + θ)*χ1)}}, {2, 2, {-((2*h - u + θ)*χ1)}}}, {{3, 0, {(u - θ)*χ1}}, {3, 1, {-((h - u + θ)*χ1)}}, {3, 2, {-((2*h - u + θ)*χ1)}}, {3, 3, {-((3*h - u + θ)*χ1)}}}}
```

## 2026-07-14 12:39:20
**Input:**
```wolfram
ClearAll[tEigenvalue];
tEigenvalue[λ_,M_][u_] := χ1 (u - θ - M h);
qAnsatz[λ_,M_][u_] := Sum[a[k] u^k, {k,0,M}];
tqEq[λ_,M_][u_] := κ1 Qθ[u - λ h] qAnsatz[λ,M][u+h] - (κ1+κ2)(u-θ-M h) qAnsatz[λ,M][u] + κ2 Qθ[u] qAnsatz[λ,M][u-h] // Expand;
sol1 = CoefficientList[tqEq[1,0][u], u]
sol2 = CoefficientList[tqEq[1,1][u], u]
```
**Output:**
```
[message] Set::write: Tag Times in sol2 {-(h κ1 a[0])} is Protected.
Out= {h*κ2*a[0] - h^2*κ1*a[1] - h*θ*κ1*a[1] + h*θ*κ2*a[1], h*κ1*a[1]}
```

## 2026-07-14 12:39:27
**Input:**
```wolfram
Clear[cf1,cf2];
cf1 = CoefficientList[tqEq[1,0][u], u];
cf2 = CoefficientList[tqEq[1,1][u], u];
{cf1, cf2}
```
**Output:**
```
Out= {{-(h*κ1*a[0])}, {h*κ2*a[0] - h^2*κ1*a[1] - h*θ*κ1*a[1] + h*θ*κ2*a[1], h*κ1*a[1]}}
```

## 2026-07-14 12:39:56
**Input:**
```wolfram
Eigenvalues[t[1][1,1][u]]//Simplify
```
**Output:**
```
Out= {(u - θ)*χ1 - (h*(χ1 + Sqrt[χ1^2 - 4*χ2]))/2, (u - θ)*χ1 + (h*(-χ1 + Sqrt[χ1^2 - 4*χ2]))/2}
```

## 2026-07-14 12:40:09
**Input:**
```wolfram
x[1][{0}]
t[1][1,1][u]
x[1][{0}].t[1][1,1][u]//Simplify
x[1][{1}]
x[1][{1}].t[1][1,1][u]//Simplify
```
**Output:**
```
Out[1]= {{0, 1}}
Out[2]= {{(-h + u - θ)*χ1, h*χ2}, {-h, (u - θ)*χ1}}
Out[3]= {{-h, (u - θ)*χ1}}
Out[4]= {{1, 0}}
Out[5]= {{-((h - u + θ)*χ1), h*χ2}}
```

## 2026-07-14 12:40:22
**Input:**
```wolfram
Table[
  Module[{xv, lhs, eigCandidate},
    xv = x[λ][{M}];
    lhs = xv.t[λ][1,1][u]//Simplify;
    eigCandidate = χ1(u-θ-M h);
    {λ, M, Simplify[lhs - eigCandidate xv]}
  ],
  {λ,1,3},{M,0,λ}
]
```
**Output:**
```
Out= {{{1, 0, {{-h, 0}}}, {1, 1, {{0, h*χ2}}}}, {{2, 0, {{0, -h, 0}}}, {2, 1, {{-1/2*h, 0, h*χ2}}}, {2, 2, {{0, h*χ2, 0}}}}, {{3, 0, {{0, 0, -h, 0}}}, {3, 1, {{0, -1/3*h, 0, h*χ2}}}, {3, 2, {{-1/6*h, 0, (2*h*χ2)/3, 0}}}, {3, 3, {{0, (h*χ2)/2, 0, 0}}}}}
```

## 2026-07-14 12:42:20
**Input:**
```wolfram
ClearAll[κ1,κ2];
subst = {χ1 -> κ1+κ2, χ2 -> κ1 κ2};
ev1 = Eigenvalues[t[1][1,1][u]] /. subst // Simplify
ev2 = Eigenvalues[t[2][1,1][u]] /. subst // Simplify
ev3 = Eigenvalues[t[3][1,1][u]] /. subst // Simplify
```
**Output:**
```
Out[1]= (suppressed)
Out[2]= (suppressed)
Out[3]= {(u - θ)*(κ1 + κ2) - (h*(κ1 + Sqrt[(κ1 - κ2)^2] + κ2))/2, (h*(-κ1 + Sqrt[(κ1 - κ2)^2] - κ2))/2 + (u - θ)*(κ1 + κ2)}
Out[4]= {-((h - u + θ)*(κ1 + κ2)), -Sqrt[h^2*(κ1 - κ2)^2] - h*(κ1 + κ2) + u*(κ1 + κ2) - θ*(κ1 + κ2), Sqrt[h^2*(κ1 - κ2)^2] - h*(κ1 + κ2) + u*(κ1 + κ2) - θ*(κ1 + κ2)}
Out[5]= {(u - θ)*(κ1 + κ2) - (h*(3*κ1 + Sqrt[(κ1 - κ2)^2] + 3*κ2))/2, (h*(-3*κ1 + Sqrt[(κ1 - κ2)^2] - 3*κ2))/2 + (u - θ)*(κ1 + κ2), -(θ*κ1) - (3*Sqrt[h^2*(κ1 - κ2)^2])/2 - θ*κ2 - (3*h*(κ1 + κ2))/2 + u*(κ1 + κ2), -(θ*κ1) + (3*Sqrt[h^2*(κ1 - κ2)^2])/2 - θ*κ2 - (3*h*(κ1 + κ2))/2 + u*(κ1 + κ2)}
```

## 2026-07-14 12:42:26
**Input:**
```wolfram
simp[e_]:=Simplify[e, κ1>κ2>0];
Table[simp[e],{e,ev1}]
Table[simp[e],{e,ev2}]
Table[simp[e],{e,ev3}]
```
**Output:**
```
Out[1]= (suppressed)
Out[2]= {-(h*κ1) + (u - θ)*(κ1 + κ2), -(h*κ2) + (u - θ)*(κ1 + κ2)}
Out[3]= {-((h - u + θ)*(κ1 + κ2)), Sqrt[h^2]*(-κ1 + κ2) - h*(κ1 + κ2) + u*(κ1 + κ2) - θ*(κ1 + κ2), Sqrt[h^2]*(κ1 - κ2) - h*(κ1 + κ2) + u*(κ1 + κ2) - θ*(κ1 + κ2)}
Out[4]= {(u - θ)*(κ1 + κ2) - h*(2*κ1 + κ2), (u - θ)*(κ1 + κ2) - h*(κ1 + 2*κ2), -(θ*κ1) - (3*Sqrt[h^2]*(κ1 - κ2))/2 - θ*κ2 - (3*h*(κ1 + κ2))/2 + u*(κ1 + κ2), -(θ*κ1) + (3*Sqrt[h^2]*(κ1 - κ2))/2 - θ*κ2 - (3*h*(κ1 + κ2))/2 + u*(κ1 + κ2)}
```

## 2026-07-14 12:42:59
**Input:**
```wolfram
tEigenvalueGuess[λ_,M_][u_] := κ1 (u-θ-M h) + κ2 (u-θ-(λ-M) h);
(* verify against Eigenvalues[t[λ][1,1][u]] for λ=1..3, matching by sorting *)
Table[
  Sort[Table[tEigenvalueGuess[λ,M][u],{M,0,λ}]] === Sort[Simplify[Eigenvalues[t[λ][1,1][u]]/.subst,Assumptions->{h>0}]],
  {λ,1,3}
]
```
**Output:**
```
Out= {False, False, False}
```

## 2026-07-14 12:43:05
**Input:**
```wolfram
Table[
  Module[{guess, actual, diffs},
    guess = Sort[Table[tEigenvalueGuess[λ,M][u],{M,0,λ}]];
    actual = Sort[Eigenvalues[t[λ][1,1][u]]/.subst];
    diffs = Simplify[guess-actual, Assumptions->{h>0,κ1>κ2>0}];
    {λ, diffs}
  ],
  {λ,1,3}
]
```
**Output:**
```
Out= {{1, {0, 0}}, {2, {h*(-κ1 + κ2), 2*h*(κ1 - κ2), h*(-κ1 + κ2)}}, {3, {h*(-κ1 + κ2), h*(κ1 - κ2), 2*h*(κ1 - κ2), 2*h*(-κ1 + κ2)}}}
```

## 2026-07-14 12:43:10
**Input:**
```wolfram
numSubst = {κ1->3, κ2->7/2, h->11/10, θ->13/10, u->17/10};
Table[
  Module[{guess, actual},
    guess = Sort[N[Table[tEigenvalueGuess[λ,M][u],{M,0,λ}] /. numSubst]];
    actual = Sort[N[(Eigenvalues[t[λ][1,1][u]]/.subst) /. numSubst]];
    {λ, Chop[guess-actual]}
  ],
  {λ,1,3}
]
```
**Output:**
```
Out= {{1, {0, 0}}, {2, {0, 0, 0}}, {3, {0, 0, 0, 0}}}
```

## 2026-07-14 12:43:28
**Input:**
```wolfram
ClearAll[tEigenvalue,qAnsatz,tqEq,a];
tEigenvalue[λ_,M_][u_] := κ1 (u - θ - M h) + κ2 (u - θ - (λ-M) h);
qAnsatz[λ_,M_][u_] := Sum[a[k] u^k, {k,0,M}];
tqEq[λ_,M_][u_] := κ1 Qθ[u - λ h] qAnsatz[λ,M][u+h] - tEigenvalue[λ,M][u] qAnsatz[λ,M][u] + κ2 Qθ[u] qAnsatz[λ,M][u-h] // Expand;
(* test λ=1,M=1 *)
cflist = CoefficientList[tqEq[1,1][u], u]
```
**Output:**
```
Out= {-(h^2*κ1*a[1]) - h*θ*κ1*a[1] + h*θ*κ2*a[1], h*κ1*a[1] - h*κ2*a[1]}
```

## 2026-07-14 12:43:35
**Input:**
```wolfram
tqEq[1,1][u]
```
**Output:**
```
Out= -(h^2*κ1*a[1]) + h*u*κ1*a[1] - h*θ*κ1*a[1] - h*u*κ2*a[1] + h*θ*κ2*a[1]
```

## 2026-07-14 12:44:13
**Input:**
```wolfram
ClearAll[tEigNaive,qAns,eqNaive,a];
tEigNaive[λ_,M_][u_] := χ1 (u - θ - M h);
qAns[λ_,M_][u_] := Sum[a[k] u^k, {k,0,M}];
eqNaive[λ_,M_][u_] := κ1 Qθ[u - λ h] qAns[λ,M][u+h] - tEigNaive[λ,M][u] qAns[λ,M][u] + κ2 Qθ[u] qAns[λ,M][u-h] // Expand;
eqNaive[1,1][u] /. χ1 -> κ1+κ2 // Expand
```
**Output:**
```
Out= h*κ2*a[0] - h^2*κ1*a[1] + h*u*κ1*a[1] - h*θ*κ1*a[1] + h*θ*κ2*a[1]
```

## 2026-07-14 12:45:20
**Input:**
```wolfram
ClearAll[tEigenvalue,qAnsatz,tqEq,a];
tEigenvalue[λ_,M_][u_] := κ1 (u-θ-(λ-M) h) + κ2 (u-θ-M h);
qAnsatz[λ_,M_][u_] := Sum[a[k] u^k, {k,0,M}];
tqEq[λ_,M_][u_] := κ1 Qθ[u - λ h] qAnsatz[λ,M][u+h] - tEigenvalue[λ,M][u] qAnsatz[λ,M][u] + κ2 Qθ[u] qAnsatz[λ,M][u-h] // Expand;
Table[CoefficientList[tqEq[1,M][u],u],{M,0,1}]
```
**Output:**
```
Out= {{}, {-(h*κ1*a[0]) + h*κ2*a[0] - h^2*κ1*a[1] - h*θ*κ1*a[1] + h*θ*κ2*a[1]}}
```

## 2026-07-14 12:45:36
**Input:**
```wolfram
solveQ[λ_,M_] := Module[{eq, cf, vars, sol},
   ClearAll[a];
   eq = tqEq[λ,M][u];
   cf = CoefficientList[eq, u];
   vars = Table[a[k],{k,0,M-1}]; (* solve for a0..a[M-1] in terms of a[M], normalize a[M]->1 *)
   sol = Solve[(cf/.a[M]->1)==0, vars];
   {cf, sol}
];
Table[{λ,M,solveQ[λ,M]},{λ,1,3},{M,0,λ}]
```
**Output:**
```
Out= {{{1, 0, {{}, {{}}}}, {1, 1, {{-(h*κ1*a[0]) + h*κ2*a[0] - h^2*κ1*a[1] - h*θ*κ1*a[1] + h*θ*κ2*a[1]}, {{a[0] -> (-(h*κ1) - θ*κ1 + θ*κ2)/(κ1 - κ2)}}}}}, {{2, 0, {{}, {{}}}}, {2, 1, {{-(h*κ1*a[0]) + h*κ2*a[0] - 2*h^2*κ1*a[1] - h*θ*κ1*a[1] + h*θ*κ2*a[1]}, {{a[0] -> (-2*h*κ1 - θ*κ1 + θ*κ2)/(κ1 - κ2)}}}}, {2, 2, {{-2*h*κ1*a[0] + 2*h*κ2*a[0] - 2*h^2*κ1*a[1] - h*θ*κ1*a[1] + h*θ*κ2*a[1] - 2*h^3*κ1*a[2] - h^2*θ*κ1*a[2] - h^2*θ*κ2*a[2], -(h*κ1*a[1]) + h*κ2*a[1] - 3*h^2*κ1*a[2] - 2*h*θ*κ1*a[2] + h^2*κ2*a[2] + 2*h*θ*κ2*a[2]}, {{a[0] -> -((-2*h^2*κ1^2 - 3*h*θ*κ1^2 - θ^2*κ1^2 + 4*h*θ*κ1*κ2 + 2*θ^2*κ1*κ2 - h*θ*κ2^2 - θ^2*κ2^2)/(κ1 - κ2)^2), a[1] -> -((3*h*κ1 + 2*θ*κ1 - h*κ2 - 2*θ*κ2)/(κ1 - κ2))}}}}}, {{3, 0, {{}, {{}}}}, {3, 1, {{-(h*κ1*a[0]) + h*κ2*a[0] - 3*h^2*κ1*a[1] - h*θ*κ1*a[1] + h*θ*κ2*a[1]}, {{a[0] -> (-3*h*κ1 - θ*κ1 + θ*κ2)/(κ1 - κ2)}}}}, {3, 2, {{-2*h*κ1*a[0] + 2*h*κ2*a[0] - 3*h^2*κ1*a[1] - h*θ*κ1*a[1] + h*θ*κ2*a[1] - 3*h^3*κ1*a[2] - h^2*θ*κ1*a[2] - h^2*θ*κ2*a[2], -(h*κ1*a[1]) + h*κ2*a[1] - 5*h^2*κ1*a[2] - 2*h*θ*κ1*a[2] + h^2*κ2*a[2] + 2*h*θ*κ2*a[2]}, {{a[0] -> -((-6*h^2*κ1^2 - 5*h*θ*κ1^2 - θ^2*κ1^2 + 6*h*θ*κ1*κ2 + 2*θ^2*κ1*κ2 - h*θ*κ2^2 - θ^2*κ2^2)/(κ1 - κ2)^2), a[1] -> -((5*h*κ1 + 2*θ*κ1 - h*κ2 - 2*θ*κ2)/(κ1 - κ2))}}}}, {3, 3, {{-3*h*κ1*a[0] + 3*h*κ2*a[0] - 3*h^2*κ1*a[1] - h*θ*κ1*a[1] + h*θ*κ2*a[1] - 3*h^3*κ1*a[2] - h^2*θ*κ1*a[2] - h^2*θ*κ2*a[2] - 3*h^4*κ1*a[3] - h^3*θ*κ1*a[3] + h^3*θ*κ2*a[3], -2*h*κ1*a[1] + 2*h*κ2*a[1] - 5*h^2*κ1*a[2] - 2*h*θ*κ1*a[2] + h^2*κ2*a[2] + 2*h*θ*κ2*a[2] - 8*h^3*κ1*a[3] - 3*h^2*θ*κ1*a[3] - h^3*κ2*a[3] - 3*h^2*θ*κ2*a[3], -(h*κ1*a[2]) + h*κ2*a[2] - 6*h^2*κ1*a[3] - 3*h*θ*κ1*a[3] + 3*h^2*κ2*a[3] + 3*h*θ*κ2*a[3]}, {{a[0] -> -((6*h^3*κ1^3 + 11*h^2*θ*κ1^3 + 6*h*θ^2*κ1^3 + θ^3*κ1^3 - 18*h^2*θ*κ1^2*κ2 - 15*h*θ^2*κ1^2*κ2 - 3*θ^3*κ1^2*κ2 + 9*h^2*θ*κ1*κ2^2 + 12*h*θ^2*κ1*κ2^2 + 3*θ^3*κ1*κ2^2 - 2*h^2*θ*κ2^3 - 3*h*θ^2*κ2^3 - θ^3*κ2^3)/(κ1 - κ2)^3), a[1] -> -((-11*h^2*κ1^2 - 12*h*θ*κ1^2 - 3*θ^2*κ1^2 + 7*h^2*κ1*κ2 + 18*h*θ*κ1*κ2 + 6*θ^2*κ1*κ2 - 2*h^2*κ2^2 - 6*h*θ*κ2^2 - 3*θ^2*κ2^2)/(κ1 - κ2)^2), a[2] -> (-3*(2*h*κ1 + θ*κ1 - h*κ2 - θ*κ2))/(κ1 - κ2)}}}}}}
```

## 2026-07-14 12:46:50
**Input:**
```wolfram
(* re-derive q_M(u) explicitly in terms of n=(u-θ)/h and ζ=κ2/κ1, normalized a[M]=1 *)
getQ[λ_,M_] := Module[{eq,cf,vars,sol,qpoly},
   ClearAll[a];
   eq = tqEq[λ,M][u];
   cf = CoefficientList[eq, u];
   vars = Table[a[k],{k,0,M-1}];
   sol = If[M==0,{},First@Solve[(cf/.a[M]->1)==0, vars]];
   qpoly = (qAnsatz[λ,M][u]/.a[M]->1)/.sol//Expand;
   qpoly /. u -> θ + n h
];
qResults = Table[{λ,M, Simplify[getQ[λ,M]/.κ2->κ1 ζ]}, {λ,1,3},{M,0,λ}];
qResults//TableForm
```
**Output:**
```
Out= TableForm[{{{1, 0, 1}, {1, 1, h*(n + (-1 + ζ)^(-1))}}, {{2, 0, 1}, {2, 1, h*(n + 2/(-1 + ζ))}, {2, 2, (h^2*(2 + n^2*(-1 + ζ)^2 - n*(3 - 4*ζ + ζ^2)))/(-1 + ζ)^2}}, {{3, 0, 1}, {3, 1, h*(n + 3/(-1 + ζ))}, {3, 2, (h^2*(6 + n^2*(-1 + ζ)^2 - n*(5 - 6*ζ + ζ^2)))/(-1 + ζ)^2}, {3, 3, (h^3*(6 - 3*n^2*(-2 + ζ)*(-1 + ζ)^2 + n^3*(-1 + ζ)^3 + n*(-11 + 18*ζ - 9*ζ^2 + 2*ζ^3)))/(-1 + ζ)^3}}}]
```

## 2026-07-14 12:46:55
**Input:**
```wolfram
Table[{λ,M,Factor[getQ[λ,M]/.κ2->κ1 ζ]},{λ,1,3},{M,0,λ}]//TableForm
```
**Output:**
```
Out= TableForm[{{{1, 0, 1}, {1, 1, (h*(1 - n + n*ζ))/(-1 + ζ)}}, {{2, 0, 1}, {2, 1, (h*(2 - n + n*ζ))/(-1 + ζ)}, {2, 2, (h^2*(2 - 3*n + n^2 + 4*n*ζ - 2*n^2*ζ - n*ζ^2 + n^2*ζ^2))/(-1 + ζ)^2}}, {{3, 0, 1}, {3, 1, (h*(3 - n + n*ζ))/(-1 + ζ)}, {3, 2, (h^2*(6 - 5*n + n^2 + 6*n*ζ - 2*n^2*ζ - n*ζ^2 + n^2*ζ^2))/(-1 + ζ)^2}, {3, 3, (h^3*(6 - 11*n + 6*n^2 - n^3 + 18*n*ζ - 15*n^2*ζ + 3*n^3*ζ - 9*n*ζ^2 + 12*n^2*ζ^2 - 3*n^3*ζ^2 + 2*n*ζ^3 - 3*n^2*ζ^3 + n^3*ζ^3))/(-1 + ζ)^3}}}]
```

## 2026-07-14 12:47:27
**Input:**
```wolfram
qGuess[λ_,M_][n_] := Sum[Binomial[M,k] Pochhammer[-n,k] Pochhammer[n-λ,M-k] (-ζ)^k,{k,0,M}]//Expand;
(* normalize to leading coeff 1 in n *)
normQ[expr_,M_] := Expand[expr/Coefficient[expr,n,M]];
Table[{λ,M, Simplify[normQ[qGuess[λ,M][n],M] - (getQ[λ,M]/.κ2->κ1 ζ)/h^M]},{λ,1,3},{M,0,λ}]
```
**Output:**
```
Out= {{{1, 0, 0}, {1, 1, (-2*ζ)/(-1 + ζ^2)}}, {{2, 0, 0}, {2, 1, (-4*ζ)/(-1 + ζ^2)}, {2, 2, (4*ζ*(-2 + n - n*ζ^2))/(-1 + ζ^2)^2}}, {{3, 0, 0}, {3, 1, (-6*ζ)/(-1 + ζ^2)}, {3, 2, (-8*ζ*(3 + n*(-1 + ζ^2)))/(-1 + ζ^2)^2}, {3, 3, (-6*ζ*(n^2*(-1 + ζ^2)^2 + 2*(3 + ζ^2) - n*(5 - 6*ζ^2 + ζ^4)))/(-1 + ζ^2)^3}}}
```

## 2026-07-14 12:47:50
**Input:**
```wolfram
qGuess2[λ_,M_][n_] := Sum[Pochhammer[-M,k] Pochhammer[-n,k]/(Pochhammer[-λ,k] k!) (1-ζ)^k, {k,0,M}]//Expand;
Table[{λ,M, Together[Simplify[normQ[qGuess2[λ,M][n],M] - (getQ[λ,M]/.κ2->κ1 ζ)/h^M]]},{λ,1,3},{M,0,λ}]
```
**Output:**
```
Out= {{{1, 0, 0}, {1, 1, 0}}, {{2, 0, 0}, {2, 1, 0}, {2, 2, 0}}, {{3, 0, 0}, {3, 1, 0}, {3, 2, 0}, {3, 3, 0}}}
```

## 2026-07-14 12:48:00
**Input:**
```wolfram
(* Direct verification: plug the closed-form Krawtchouk-type q_M into the ORIGINAL TQ relation (before normalization) *)
qClosed[λ_,M_][u_] := (Sum[Pochhammer[-M,k] Pochhammer[-(u-θ)/h,k]/(Pochhammer[-λ,k] k!) (1-κ2/κ1)^k, {k,0,M}]);
tqEqClosed[λ_,M_][u_] := κ1 Qθ[u - λ h] qClosed[λ,M][u+h] - tEigenvalue[λ,M][u] qClosed[λ,M][u] + κ2 Qθ[u] qClosed[λ,M][u-h] // Simplify;
Table[tqEqClosed[λ,M][u],{λ,1,4},{M,0,λ}]//Flatten//Union
```
**Output:**
```
Out= {0}
```

## 2026-07-14 12:48:49
**Input:**
```wolfram
Qfull[λ_,M_][u_] := κ1^(u/h) * (qClosed[λ,M][u] /. u -> u);
fullCheck[λ_,M_][u_] := Module[{expr},
  expr = Qθ[u-λ h] Qfull[λ,M][u+h] - tEigenvalue[λ,M][u] Qfull[λ,M][u] + (κ1 κ2) Qθ[u] Qfull[λ,M][u-h];
  expr = expr /. {κ1^(u/h+1) -> κ1 κ1^(u/h), κ1^(u/h-1) -> κ1^(u/h)/κ1};
  Simplify[expr / κ1^(u/h)]
];
Table[fullCheck[λ,M][u], {λ,1,3},{M,0,λ}]//Flatten//Union
```
**Output:**
```
Out= {0}
```
## 2026-07-14 12:49:19 — 📥 BULK INSERT 10 CELL(S) at positions 15–24
1. [markdown] (* Baxter TQ-equation solver, part 1: eigenvector sanity check.↵↵Introduce κ1,κ2 (eigenvalues of the
2. [code] ClearAll[κ1,κ2];↵twistSubst = {χ1 -> κ1+κ2, χ2 -> κ1 κ2};↵↵(* weak check: ratio at nonzero-matching 
3. [code] (* strict check: full residual against the candidate eigenvalue χ1(u-θ-Mh) suggested by the weak che
4. [markdown] (* Sanity-check result: the weak check passes trivially (a single well-defined ratio↵χ1(u-θ-Mh) at t
5. [code] (* the genuine eigenvalues of t[λ][1,1][u], compared (numerically, several sample points) against↵th
6. [markdown] (* Baxter TQ-equation solver, part 2: solve for q(u).↵↵Relation: Qθ[u-hλ] Q(u+h) - t(u) Q(u) + χ2 Qθ
7. [code] ClearAll[qAnsatz,tqEq,a,solveQ,qSolution];↵qAnsatz[λ_,M_][u_] := Sum[a[k] u^k, {k,0,M}];↵tqEq[λ_,M_]
8. [code] (* print q(u) explicitly for each (λ,M) *)↵Table[{λ,M, qSolution[λ,M][u]}, {λ,1,3},{M,0,λ}] // Table
9. [markdown] (* Closed-form pattern.↵↵Writing n=(u-θ)/h and ζ=κ2/κ1, the explicit q_M(u) computed above (for each
10. [code] ClearAll[qClosed];↵qClosed[λ_,M_][u_] := Sum[Pochhammer[-M,k] Pochhammer[-(u-θ)/h,k]/(Pochhammer[-λ,


## 2026-07-14 12:49:19
**Input:**
```wolfram
ClearAll[κ1,κ2];
twistSubst = {χ1 -> κ1+κ2, χ2 -> κ1 κ2};

(* weak check: ratio at nonzero-matching components only *)
weakCheck = Table[
   Module[{xv,lhs,nz,ratios},
     xv = x[λ][{M}];
     lhs = xv.t[λ][1,1][u]//Simplify;
     nz = Select[Range[Length[xv[[1]]]], xv[[1,#]]=!=0 &];
     ratios = Union[(Simplify[lhs[[1,#]]/xv[[1,#]]])& /@ nz];
     {λ,M,ratios}
   ],{λ,1,3},{M,0,λ}];
weakCheck
```
**Output:**
```
Out[32]= \left\{\begin{pmatrix}1 & 0 & \left\{\chi1\,\left(u-\theta \right)\right\}\\1 & 1 & \left\{-\left(\chi1\,\left(h+\theta -u\right)\right)\right\}\end{pmatrix},\begin{pmatrix}2 & 0 & \left\{\chi1\,\left(u-\theta \right)\right\}\\2 & 1 & \left\{-\left(\chi1\,\left(h+\theta -u\right)\right)\right\}\\2 & 2 & \left\{-\left(\chi1\,\left(2\,h+\theta -u\right)\right)\right\}\end{pmatrix},\begin{pm… [695 chars total]
```

## 2026-07-14 12:49:19
**Input:**
```wolfram
(* strict check: full residual against the candidate eigenvalue χ1(u-θ-Mh) suggested by the weak check *)
strictCheck = Table[
   Module[{xv,lhs,ratio,resid},
     xv = x[λ][{M}];
     lhs = xv.t[λ][1,1][u]//Simplify;
     ratio = χ1(u-θ-M h);
     resid = Simplify[lhs - ratio xv];
     {λ,M,resid}
   ],{λ,1,3},{M,0,λ}];
strictCheck
```
**Output:**
```
Out[34]= \left\{\begin{pmatrix}1 & 0 & \begin{pmatrix}-h & 0\end{pmatrix}\\1 & 1 & \begin{pmatrix}0 & h\,\chi2\end{pmatrix}\end{pmatrix},\begin{pmatrix}2 & 0 & \begin{pmatrix}0 & -h & 0\end{pmatrix}\\2 & 1 & \begin{pmatrix}-\frac{h}{2} & 0 & h\,\chi2\end{pmatrix}\\2 & 2 & \begin{pmatrix}0 & h\,\chi2 & 0\end{pmatrix}\end{pmatrix},\begin{pmatrix}3 & 0 & \begin{pmatrix}0 & 0 & -h & 0\end{pmatrix}\\3 … [636 chars total]
```

## 2026-07-14 12:49:20
**Input:**
```wolfram
(* the genuine eigenvalues of t[λ][1,1][u], compared (numerically, several sample points) against
the candidate closed form κ1(u-θ-(λ-M)h)+κ2(u-θ-Mh); matching verified as a SET for λ=1..4 *)
ClearAll[tEigenvalue];
tEigenvalue[λ_,M_][u_] := κ1 (u-θ-(λ-M) h) + κ2 (u-θ-M h);

numSubst = {κ1->3, κ2->7/2, h->11/10, θ->13/10, u->17/10};
Table[
  Module[{guess,actual},
    guess = Sort[N[Table[tEigenvalue[λ,M][u],{M,0,λ}] /. numSubst]];
    actual = Sort[N[(Eigenvalues[t[λ][1,1][u]]/.twistSubst) /. numSubst]];
    {λ, Chop[guess-actual]}
  ],
{λ,1,4}]
```
**Output:**
```
Out[38]= \begin{pmatrix}1 & \left\{0,0\right\}\\2 & \left\{0,0,0\right\}\\3 & \left\{0,0,0,0\right\}\\4 & \left\{0,0,0,0,0\right\}\end{pmatrix}
```

## 2026-07-14 12:49:20
**Input:**
```wolfram
ClearAll[qAnsatz,tqEq,a,solveQ,qSolution];
qAnsatz[λ_,M_][u_] := Sum[a[k] u^k, {k,0,M}];
tqEq[λ_,M_][u_] := κ1 Qθ[u - λ h] qAnsatz[λ,M][u+h] - tEigenvalue[λ,M][u] qAnsatz[λ,M][u] + κ2 Qθ[u] qAnsatz[λ,M][u-h] // Expand;

solveQ[λ_,M_] := Module[{eq,cf,vars,sol},
   eq = tqEq[λ,M][u] /. a[M]->1;
   cf = CoefficientList[eq /. a[M]->1, u];
   vars = Table[a[k],{k,0,M-1}];
   sol = If[M==0, {}, First@Solve[cf==0, vars]];
   qAnsatz[λ,M][u] /. a[M]->1 /. sol // Expand
];
qSolution[λ_,M_][u_] := qSolution[λ,M][u] = solveQ[λ,M];

(* run for λ=1,2,3, all M=0..λ, and verify the ORIGINAL (unstripped) TQ relation holds
identically once Q(u)=κ1^(u/h) qSolution(u) is substituted back *)
fullVerify[λ_,M_][u_] := Module[{Qf, expr},
   Qf[v_] := κ1^(v/h) qSolution[λ,M][v];
   expr = Qθ[u-λ h] Qf[u+h] - tEigenvalue[λ,M][u] Qf[u] + (κ1 κ2) Qθ[u] Qf[u-h];
   expr = expr /. {κ1^(u/h+1) -> κ1 κ1^(u/h), κ1^(u/h-1) -> κ1^(u/h)/κ1};
   Simplify[expr / κ1^(u/h)]
];
Table[fullVerify[λ,M][u], {λ,1,3},{M,0,λ}]//Flatten//Union
```
**Output:**
```
Out[45]= \left\{0,h\,\left(h\,\kappa1-\left(\kappa1-\kappa2\right)\,\left(u-\theta \right)\right),h\,\left(2\,h\,\kappa1-\left(\kappa1-\kappa2\right)\,\left(u-\theta \right)\right),h\,\left(3\,h\,\kappa1-\left(\kappa1-\kappa2\right)\,\left(u-\theta \right)\right),-\frac{2\,h\,\left(6\,h^2\,\kappa1^2-h\,\left(5\,\kappa1^2-6\,\kappa1\,\kappa2+\kappa2^2\right)\,\left(u-\theta \right)+\left(\kappa1-\k… [1018 chars total]
```

## 2026-07-14 12:49:21
**Input:**
```wolfram
(* print q(u) explicitly for each (λ,M) *)
Table[{λ,M, qSolution[λ,M][u]}, {λ,1,3},{M,0,λ}] // TableForm[#, TableHeadings->{None,{"λ","M","q(u)"}}]&
```
**Output:**
```
Out[46]//TableForm= \left\{\begin{pmatrix}1 & 0 & 1\\1 & 1 & u+\frac{\theta \,\kappa2}{\kappa1-\kappa2}-\frac{h\,\kappa1}{\kappa1-\kappa2}-\frac{\theta \,\kappa1}{\kappa1-\kappa2}\end{pmatrix},\begin{pmatrix}2 & 0 & 1\\2 & 1 & u+\frac{\theta \,\kappa2}{\kappa1-\kappa2}-\frac{2\,h\,\kappa1}{\kappa1-\kappa2}-\frac{\theta \,\kappa1}{\kappa1-\kappa2}\\2 & 2 & u^2+\frac{h\,\kappa2\,u}{\kappa1-\kappa2}+… [1257 chars total]
```

## 2026-07-14 12:49:21
**Input:**
```wolfram
ClearAll[qClosed];
qClosed[λ_,M_][u_] := Sum[Pochhammer[-M,k] Pochhammer[-(u-θ)/h,k]/(Pochhammer[-λ,k] k!) (1-κ2/κ1)^k, {k,0,M}];

tqEqClosed[λ_,M_][u_] := κ1 Qθ[u - λ h] qClosed[λ,M][u+h] - tEigenvalue[λ,M][u] qClosed[λ,M][u] + κ2 Qθ[u] qClosed[λ,M][u-h] // Simplify;
Table[tqEqClosed[λ,M][u], {λ,1,4},{M,0,λ}] // Flatten // Union
```
**Output:**
```
Out[50]= \left\{0\right\}
```
## 2026-07-14 12:49:55 — ✏️ EDIT CELL 21
ClearAll[qAnsatz,tqEq,a,solveQ,qPoly,qSolution];
qAnsatz[λ_,M_][u_] := Sum[a[k] u^k, {k,0,M}];
tqEq[λ_,M_][u_] := κ1 Qθ[u - λ h] qAnsatz[λ,M][u+h] - tEigenvalue[λ,M][u] qAnsatz[λ,M][u] + κ2 Qθ[u] qAns…


## 2026-07-14 12:49:55
**Input:**
```wolfram
ClearAll[qAnsatz,tqEq,a,solveQ,qPoly,qSolution];
qAnsatz[λ_,M_][u_] := Sum[a[k] u^k, {k,0,M}];
tqEq[λ_,M_][u_] := κ1 Qθ[u - λ h] qAnsatz[λ,M][u+h] - tEigenvalue[λ,M][u] qAnsatz[λ,M][u] + κ2 Qθ[u] qAnsatz[λ,M][u-h] // Expand;

solveQ[λ_,M_] := Module[{cf,vars,sol},
   cf = CoefficientList[tqEq[λ,M][u] /. a[M]->1, u];
   vars = Table[a[k],{k,0,M-1}];
   sol = If[M==0, {}, First@Solve[cf==0, vars]];
   qAnsatz[λ,M][u] /. a[M]->1 /. sol // Expand
];
(* memoize the u-polynomial itself, then substitute the actual argument on every call
   (memoizing qSolution[λ,M][u]:=... directly would freeze the shift variable to the
   literal global symbol u, breaking qSolution[λ,M][u+h] etc. - substitute explicitly instead) *)
qPoly[λ_,M_] := qPoly[λ,M] = solveQ[λ,M];
qSolution[λ_,M_][uArg_] := qPoly[λ,M] /. u -> uArg;

(* run for λ=1,2,3, all M=0..λ, and verify the ORIGINAL (unstripped) TQ relation holds
identically once Q(u)=κ1^(u/h) qSolution(u) is substituted back *)
fullVerify[λ_,M_][u_] := Module[{Qf, expr},
   Qf[v_] := κ1^(v/h) qSolution[λ,M][v];
   expr = Qθ[u-λ h] Qf[u+h] - tEigenvalue[λ,M][u] Qf[u] + (κ1 κ2) Qθ[u] Qf[u-h];
   expr = expr /. {κ1^(u/h+1) -> κ1 κ1^(u/h), κ1^(u/h-1) -> κ1^(u/h)/κ1};
   Simplify[expr / κ1^(u/h)]
];
Table[fullVerify[λ,M][u], {λ,1,3},{M,0,λ}]//Flatten//Union
```
**Output:**
```
Out= {0}
```
## 2026-07-14 12:50:00 — ▶️ RUN CELLS 22–24
Cell 22: ✓ — Out[51]//TableForm= \left\{\begin{pmatrix}1 & 0 & 1\\1 & 1 & u+\frac{\theta \,\kappa2}{\kappa1-\kappa2}-\frac{h\,\kappa1}{\kappa1-\kappa2}-\frac{\theta \,\kappa1}{\kappa1-\kappa2}\end{pmatrix},\begin{pmatrix}2 & 0 & 1\\2 & 1 & u+\frac{\theta \,\kappa2}{\kappa1-\kappa2}-\frac{2\,h\,\kappa1}{\kappa1-\kappa2}-\frac{\theta \,\kappa1}{\kappa1-\kappa2}\\2 & 2 & u^2+\frac{h\,\kappa2\,u}{\kappa1-\kappa2}+… [1257 chars total]
Cell 24: ✓ — Out[55]= \left\{0\right\}


## 2026-07-14 12:56:33
**Input:**
```wolfram
t[1][1,1][u]
```
**Output:**
```
Out= {{(-h + u - θ)*χ1, h*χ2}, {-h, (u - θ)*χ1}}
```

## 2026-07-14 12:56:33
**Input:**
```wolfram
dim[{1,0}]
```
**Output:**
```
Out= 2
```

## 2026-07-14 12:56:37
**Input:**
```wolfram
Eigensystem[t[1][1,1][u]] // Simplify
```
**Output:**
```
Out= {{(u - θ)*χ1 - (h*(χ1 + Sqrt[χ1^2 - 4*χ2]))/2, (u - θ)*χ1 + (h*(-χ1 + Sqrt[χ1^2 - 4*χ2]))/2}, {{(χ1 + Sqrt[χ1^2 - 4*χ2])/2, 1}, {(χ1 - Sqrt[χ1^2 - 4*χ2])/2, 1}}}
```

## 2026-07-14 12:56:37
**Input:**
```wolfram
Eigensystem[t[2][1,1][u]] // Simplify
```
**Output:**
```
Out= {{-((h - u + θ)*χ1), -(h*χ1) + u*χ1 - θ*χ1 - Sqrt[h^2*(χ1^2 - 4*χ2)], -(h*χ1) + u*χ1 - θ*χ1 + Sqrt[h^2*(χ1^2 - 4*χ2)]}, {{2*χ2, χ1, 1}, {χ1^2 + (χ1*Sqrt[h^2*(χ1^2 - 4*χ2)])/h - 2*χ2, (h*χ1 + Sqrt[h^2*(χ1^2 - 4*χ2)])/h, 1}, {χ1^2 - (χ1*Sqrt[h^2*(χ1^2 - 4*χ2)])/h - 2*χ2, (h*χ1 - Sqrt[h^2*(χ1^2 - 4*χ2)])/h, 1}}}
```

## 2026-07-14 12:56:43
**Input:**
```wolfram
Eigensystem[t[2][1,1][u]] /. twistSubst // PowerExpand // Simplify
```
**Output:**
```
Out= {{-((h - u + θ)*(κ1 + κ2)), (u - θ)*(κ1 + κ2) - h*(κ1 + Sqrt[(κ1 - κ2)^2] + κ2), h*(-κ1 + Sqrt[(κ1 - κ2)^2] - κ2) + (u - θ)*(κ1 + κ2)}, {{2*κ1*κ2, κ1 + κ2, 1}, {κ1^2 + κ1*Sqrt[(κ1 - κ2)^2] + κ2*(Sqrt[(κ1 - κ2)^2] + κ2), κ1 + Sqrt[(κ1 - κ2)^2] + κ2, 1}, {κ1^2 - κ1*Sqrt[(κ1 - κ2)^2] + κ2*(-Sqrt[(κ1 - κ2)^2] + κ2), κ1 - Sqrt[(κ1 - κ2)^2] + κ2, 1}}}
```

## 2026-07-14 12:56:50
**Input:**
```wolfram
Eigensystem[t[2][1,1][u]] /. twistSubst /. Sqrt[(κ1-κ2)^2]->κ1-κ2 // Simplify
```
**Output:**
```
Out= {{-((h - u + θ)*(κ1 + κ2)), -Sqrt[h^2*(κ1 - κ2)^2] - h*(κ1 + κ2) + u*(κ1 + κ2) - θ*(κ1 + κ2), Sqrt[h^2*(κ1 - κ2)^2] - h*(κ1 + κ2) + u*(κ1 + κ2) - θ*(κ1 + κ2)}, {{2*κ1*κ2, κ1 + κ2, 1}, {κ1^2 + κ2^2 + (Sqrt[h^2*(κ1 - κ2)^2]*(κ1 + κ2))/h, κ1 + Sqrt[h^2*(κ1 - κ2)^2]/h + κ2, 1}, {κ1^2 + κ2^2 - (Sqrt[h^2*(κ1 - κ2)^2]*(κ1 + κ2))/h, κ1 - Sqrt[h^2*(κ1 - κ2)^2]/h + κ2, 1}}}
```

## 2026-07-14 12:56:50
**Input:**
```wolfram
Eigensystem[t[3][1,1][u]] /. twistSubst // Simplify
```
**Output:**
```
Out= {{(u - θ)*(κ1 + κ2) - (h*(3*κ1 + Sqrt[(κ1 - κ2)^2] + 3*κ2))/2, (h*(-3*κ1 + Sqrt[(κ1 - κ2)^2] - 3*κ2))/2 + (u - θ)*(κ1 + κ2), (-3*Sqrt[h^2*(κ1 - κ2)^2])/2 - (3*h*(κ1 + κ2))/2 + u*(κ1 + κ2) - θ*(κ1 + κ2), (3*Sqrt[h^2*(κ1 - κ2)^2])/2 - (3*h*(κ1 + κ2))/2 + u*(κ1 + κ2) - θ*(κ1 + κ2)}, {{3*κ1*κ2*(κ1 + Sqrt[(κ1 - κ2)^2] + κ2), 2*κ1*κ2 + Sqrt[(κ1 - κ2)^2]*(κ1 + κ2) + (κ1 + κ2)^2, (3*κ1 + Sqrt[(κ1 - κ2)^2] + 3*κ2)/2, 1}, {3*κ1*κ2*(κ1 - Sqrt[(κ1 - κ2)^2] + κ2), 2*κ1*κ2 - Sqrt[(κ1 - κ2)^2]*(κ1 + κ2) + (κ1 + κ2)^2, (3*κ1 - Sqrt[(κ1 - κ2)^2] + 3*κ2)/2, 1}, {3*(κ1^3 + κ2^3 + (Sqrt[h^2*(κ1 - κ2)^2]*(κ1^2 + κ1*κ2 + κ2^2))/h), 3*(κ1^2 + κ2^2 + (Sqrt[h^2*(κ1 - κ2)^2]*(κ1 + κ2))/h), (3*(Sqrt[h^2*(κ1 - κ2)^2] + h*(κ1 + κ2)))/(2*h), 1}, {3*(κ1^3 + κ2^3 - (Sqrt[h^2*(κ1 - κ2)^2]*(κ1^2 + κ1*κ2 + κ2^2))/h), 3*(κ1^2 + κ2^2 - (Sqrt[h^2*(κ1 - κ2)^2]*(κ1 + κ2))/h), (3*(-Sqrt[h^2*(κ1 - κ2)^2] + h*(κ1 + κ2)))/(2*h), 1}}}
```

## 2026-07-14 12:57:00
**Input:**
```wolfram
Eigensystem[t[3][1,1][u]] /. twistSubst // Simplify[#, κ1>κ2>0 && h>0]&
```
**Output:**
```
Out= {{(u - θ)*(κ1 + κ2) - h*(2*κ1 + κ2), (u - θ)*(κ1 + κ2) - h*(κ1 + 2*κ2), -3*h*κ1 + (u - θ)*(κ1 + κ2), -3*h*κ2 + u*(κ1 + κ2) - θ*(κ1 + κ2)}, {{6*κ1^2*κ2, 2*κ1*(κ1 + 2*κ2), 2*κ1 + κ2, 1}, {6*κ1*κ2^2, 2*κ2*(2*κ1 + κ2), κ1 + 2*κ2, 1}, {6*κ1^3, 6*κ1^2, 3*κ1, 1}, {6*κ2^3, 6*κ2^2, 3*κ2, 1}}}
```

## 2026-07-14 12:58:32
**Input:**
```wolfram
Simplify[(t[2][1,1][u] /. twistSubst) - (t[2][1,1][0] /. twistSubst) - u (κ1+κ2) IdentityMatrix[dim[{2,0}]]]
```
**Output:**
```
Out= {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}}
```

## 2026-07-14 12:58:50
**Input:**
```wolfram
ClearAll[testΨ];
testΨ[λ_,M_] := Module[{tm, es, vals, vecs, pos},
   tm = t[λ][1,1][u] /. twistSubst /. u -> 0;
   es = Eigensystem[tm];
   vals = es[[1]];
   vecs = es[[2]];
   pos = First @ FirstPosition[vals, v_ /; Simplify[v - (tEigenvalue[λ,M][u] /. u -> 0)] === 0];
   vecs[[pos]]
];
Table[Simplify[(t[λ][1,1][u]/.twistSubst).testΨ[λ,M] - tEigenvalue[λ,M][u] testΨ[λ,M]], {λ,2,3},{M,0,λ}]//Flatten//Union
```
**Output:**
```
Out[1]= (suppressed)
Out[2]= (suppressed)
Out[3]= {0}
```
## 2026-07-14 12:59:01 — 📥 BULK INSERT 4 CELL(S) at positions 25–28
1. [markdown] (* Baxter TQ-equation solver, part 3: eigenvector Ψ[λ,M] of t[λ][1,1][u].↵↵Key observation: t[λ][1,1
2. [code] ClearAll[Ψ];↵Ψ[λ_,M_] := Module[{tm,es,vals,vecs,pos},↵   tm = t[λ][1,1][u] /. twistSubst /. u -> 0;
3. [markdown] (* Sanity check: Ψ[λ,M] must satisfy t[λ][1,1][u].Ψ[λ,M] == tEigenvalue[λ,M][u] Ψ[λ,M]↵identically i
4. [code] Table[↵  Simplify[(t[λ][1,1][u]/.twistSubst).Ψ[λ,M] - tEigenvalue[λ,M][u] Ψ[λ,M]],↵{λ,2,3},{M,0,λ}]


## 2026-07-14 12:59:02
**Input:**
```wolfram
ClearAll[Ψ];
Ψ[λ_,M_] := Module[{tm,es,vals,vecs,pos},
   tm = t[λ][1,1][u] /. twistSubst /. u -> 0;
   es = Eigensystem[tm];
   vals = es[[1]];
   vecs = es[[2]];
   pos = First @ FirstPosition[vals, v_ /; Simplify[v - (tEigenvalue[λ,M][u] /. u -> 0)] === 0];
   vecs[[pos]]
];
```
**Output:**
```
(no output)
```

## 2026-07-14 12:59:02
**Input:**
```wolfram
Table[
  Simplify[(t[λ][1,1][u]/.twistSubst).Ψ[λ,M] - tEigenvalue[λ,M][u] Ψ[λ,M]],
{λ,2,3},{M,0,λ}] // Flatten // Union
```
**Output:**
```
Out[58]= \left\{0\right\}
```

## 2026-07-14 13:05:00
**Input:**
```wolfram
Table[Simplify[First[x[λ][{0}].Ψ[λ,M]]], {λ,1,4},{M,0,λ}]
```
**Output:**
```
Out= {{1, 1}, {1, 1, 1}, {1, 1, 1, 1}, {1, 1, 1, 1, 1}}
```

## 2026-07-14 13:05:10
**Input:**
```wolfram
Module[{tm,es,vals,vecs},tm = t[2][1,1][u] /. twistSubst /. u -> 0; es = Eigensystem[tm]; es[[2]]]//Simplify
```
**Output:**
```
Out= {{2*κ1*κ2, κ1 + κ2, 1}, {2*κ1^2, 2*κ1, 1}, {2*κ2^2, 2*κ2, 1}}
```
## 2026-07-14 13:05:20 — ✏️ EDIT CELL 26
ClearAll[Ψ];
Ψ[λ_,M_] := Module[{tm,es,vals,vecs,pos,raw,overlap},
   tm = t[λ][1,1][u] /. twistSubst /. u -> 0;
   es = Eigensystem[tm];
   vals = es[[1]];
   vecs = es[[2]];
   pos = First @ FirstPo…


## 2026-07-14 13:05:20
**Input:**
```wolfram
ClearAll[Ψ];
Ψ[λ_,M_] := Module[{tm,es,vals,vecs,pos,raw,overlap},
   tm = t[λ][1,1][u] /. twistSubst /. u -> 0;
   es = Eigensystem[tm];
   vals = es[[1]];
   vecs = es[[2]];
   pos = First @ FirstPosition[vals, v_ /; Simplify[v - (tEigenvalue[λ,M][u] /. u -> 0)] === 0];
   raw = vecs[[pos]];
   (* normalize so that x[λ][{0}].Ψ[λ,M] == 1; overlap is proportional to the bottom
      component of raw (x[λ][{0}] is a KroneckerDelta pick-out of the last entry) *)
   overlap = Simplify[First[x[λ][{0}].raw]];
   If[overlap === 0, Message[Ψ::zeroOverlap,λ,M]; $Failed, raw/overlap]
];
Ψ::zeroOverlap = "Raw eigenvector overlap with x[`1`][{0}] vanishes for λ=`1`, M=`2` - cannot normalize.";
```
**Output:**
```
Out= Null
```
## 2026-07-14 13:05:27 — ▶️ RUN CELL 28 (0.03 s)
**In [28]:** `Table[↵  Simplify[(t[λ][1,1][u]/.twistSubst).Ψ[λ,M] - tEigenvalue[λ,M][u] Ψ[λ,M]],↵{λ,2,3},{M,0,λ}] // Flatten // Union`
**Out:** Out[59]= \left\{0\right\}

## 2026-07-14 13:05:30 — 📥 BULK INSERT 2 CELL(S) at positions 29–30
1. [markdown] (* Normalization check: x[λ][{0}].Ψ[λ,M] == 1 identically, for λ=1,...,4, all valid M. *)
2. [code] Table[↵  Simplify[First[x[λ][{0}].Ψ[λ,M]] - 1],↵{λ,1,4},{M,0,λ}] // Flatten // Union


## 2026-07-14 13:05:31
**Input:**
```wolfram
Table[
  Simplify[First[x[λ][{0}].Ψ[λ,M]] - 1],
{λ,1,4},{M,0,λ}] // Flatten // Union
```
**Output:**
```
Out[60]= \left\{0\right\}
```

## 2026-07-14 13:08:19
**Input:**
```wolfram
ClearAll[QFull]; QFull[λ_,M_][v_] := κ1^(v/h) qSolution[λ,M][v];
Table[Simplify[First[x[λ][{s}].Ψ[λ,M]] - QFull[λ,M][θ+h s]/QFull[λ,M][θ]], {λ,1,4},{M,0,λ},{s,0,λ}] // Flatten // Union
```
**Output:**
```
Out[1]= (suppressed)
Out[2]= {0}
```

## 2026-07-14 13:08:25
**Input:**
```wolfram
Table[Simplify[First[x[λ][{s}].Ψ[λ,M]] - qSolution[λ,M][θ+h s]/qSolution[λ,M][θ]], {λ,1,4},{M,1,λ},{s,1,λ}] // Flatten // Union
```
**Output:**
```
Out= {((-1 + κ1)*κ2)/κ1, ((-1 + κ1^2)*κ2)/κ1, ((-1 + κ1^3)*κ2)/κ1, ((-1 + κ1^4)*κ2)/κ1, ((-1 + κ1^2)*κ2^2)/κ1^2, ((-1 + κ1^3)*κ2^2)/κ1^2, ((-1 + κ1^4)*κ2^2)/κ1^2, ((-1 + κ1^3)*κ2^3)/κ1^3, ((-1 + κ1^4)*κ2^3)/κ1^3, ((-1 + κ1^4)*κ2^4)/κ1^4, ((-1 + κ1)*(κ1 + κ2))/(2*κ1), ((-1 + κ1^2)*(κ1 + κ2))/(2*κ1), ((-1 + κ1^2)*κ2*(κ1 + κ2))/(2*κ1^2), ((-1 + κ1^3)*κ2*(κ1 + κ2))/(2*κ1^2), ((-1 + κ1)*(2*κ1 + κ2))/(3*κ1), ((-1 + κ1^2)*κ2*(2*κ1 + κ2))/(3*κ1^2), ((-1 + κ1)*(3*κ1 + κ2))/(4*κ1), ((-1 + κ1^3)*κ2^2*(3*κ1 + κ2))/(4*κ1^3), ((-1 + κ1)*(κ1 + 2*κ2))/(3*κ1), ((-1 + κ1^2)*(κ1 + 2*κ2))/(3*κ1), ((-1 + κ1)*(κ1 + 3*κ2))/(4*κ1), ((-1 + κ1^3)*(κ1 + 3*κ2))/(4*κ1), ((-1 + κ1^2)*(κ1^2 + 4*κ1*κ2 + κ2^2))/(6*κ1^2)}
```
## 2026-07-14 13:08:39 — 📥 BULK INSERT 4 CELL(S) at positions 31–34
1. [markdown] (* SoV wavefunction reconstruction check.↵↵Claim: x[λ][{s}].Ψ[λ,M] == Q_M(θ+hs)/Q_M(θ), s=0,...,λ, w
2. [code] ClearAll[QFull];↵QFull[λ_,M_][v_] := κ1^(v/h) qSolution[λ,M][v];↵↵reconstructionCheck = Table[↵   Si
3. [code] (* control: same check but with the κ1^(u/h) prefactor OMITTED (bare polynomial ratio↵qSolution[λ,M]
4. [markdown] (* Result: reconstructionCheck == {0} identically for all λ=1,...,4, M=0,...,λ, s=0,...,λ -↵the SoV


## 2026-07-14 13:08:40
**Input:**
```wolfram
ClearAll[QFull];
QFull[λ_,M_][v_] := κ1^(v/h) qSolution[λ,M][v];

reconstructionCheck = Table[
   Simplify[First[x[λ][{s}].Ψ[λ,M]] - QFull[λ,M][θ+h s]/QFull[λ,M][θ]],
{λ,1,4},{M,0,λ},{s,0,λ}] // Flatten // Union
```
**Output:**
```
Out[63]= \left\{0\right\}
```

## 2026-07-14 13:08:41
**Input:**
```wolfram
(* control: same check but with the κ1^(u/h) prefactor OMITTED (bare polynomial ratio
qSolution[λ,M][θ+hs]/qSolution[λ,M][θ] only) - confirms the prefactor is NOT redundant,
i.e. it does not cancel in the ratio and is genuinely required for the identity to hold *)
Table[
   Simplify[First[x[λ][{s}].Ψ[λ,M]] - qSolution[λ,M][θ+h s]/qSolution[λ,M][θ]],
{λ,1,4},{M,1,λ},{s,1,λ}] // Flatten // Union // Length
```
**Output:**
```
Out[64]= 23
```

## 2026-07-14 13:16:57
**Input:**
```wolfram
ClearAll[μ]; μ[λ_][u_] := Product[1/(u - θ - h(k-1)), {k, 1, λ}]
ClearAll[residueSum]; residueSum[λ_][f_] := 2 Pi I Sum[Residue[μ[λ][u] f[u], {u, θ + h(k-1)}], {k, 1, λ}]
residueSum[1][f]
residueSum[2][f]
residueSum[3][f]
residueSum[4][f]
```
**Output:**
```
Out[1]= Null
Out[2]= Null
Out[3]= (2*I)*Pi*f[θ]
Out[4]= (2*I)*Pi*(-(f[θ]/h) + f[h + θ]/h)
Out[5]= (2*I)*Pi*(f[θ]/(2*h^2) - f[h + θ]/h^2 + f[2*h + θ]/(2*h^2))
Out[6]= (2*I)*Pi*(-1/6*f[θ]/h^3 + f[h + θ]/(2*h^3) - f[2*h + θ]/(2*h^3) + f[3*h + θ]/(6*h^3))
```

## 2026-07-14 13:17:03
**Input:**
```wolfram
ClearAll[normConst]; normConst[λ_] := (-1)^(λ-1) h^(λ-1) (λ-1)! / (2 Pi I)
ClearAll[normalizedSum]; normalizedSum[λ_][f_] := Simplify[normConst[λ] residueSum[λ][f]]
ClearAll[claimedClosed]; claimedClosed[λ_][f_] := Sum[(-1)^i Binomial[λ-1,i] f[θ+h i], {i,0,λ-1}]
Table[Simplify[normalizedSum[λ][f] - claimedClosed[λ][f]], {λ,1,6}]
```
**Output:**
```
Out[1]= Null
Out[2]= Null
Out[3]= Null
Out[4]= {0, 0, 0, 0, 0, 0}
```

## 2026-07-14 13:17:07
**Input:**
```wolfram
Table[
  Simplify[
    (normConst[λ] * 2 Pi I Sum[Residue[μ[λ][u] u^n, {u, θ + h(k-1)}], {k,1,λ}])
      - Sum[(-1)^i Binomial[λ-1,i] (θ+h i)^n, {i,0,λ-1}]
  ],
  {λ,1,4},{n,0,λ-1}
] // Flatten // Union
```
**Output:**
```
Out= {0}
```

## 2026-07-14 13:17:11
**Input:**
```wolfram
ClearAll[SoVMeasure];
SoVMeasure[λ_][f_] := Sum[(-1)^i Binomial[λ-1,i] f[θ+h i], {i,0,λ-1}];
(* cross-check SoVMeasure agrees with the raw normalized contour-integral definition, symbolic f and monomials *)
Table[Simplify[SoVMeasure[λ][f] - normalizedSum[λ][f]], {λ,1,6}]
```
**Output:**
```
Out[1]= (suppressed)
Out[2]= (suppressed)
Out[3]= {0, 0, 0, 0, 0, 0}
```
## 2026-07-14 13:17:48 — 📥 BULK INSERT 7 CELL(S) at positions 35–41
1. [md] (* SoV measure: normalized contour-integral functional.↵↵Define μ[λ](u) = Product[1/(u - θ - h(k-1))
2. [code] (* measure with λ simple poles at u = θ+h(k-1), k=1,...,λ *)↵ClearAll[μ];↵μ[λ_][u_] := Product[1/(u 
3. [code] (* raw contour functional = 2πi * sum of residues, computed directly via Residue *)↵ClearAll[residue
4. [code] (* normalization constant fixing the coefficient of f[θ] to 1, and the claimed alternating-binomial↵
5. [code] (* full check on a monomial basis f(u)=u^n, n=0,...,λ-1 (legitimate since the functional is↵linear):
6. [code] (* final definition: normalized SoV measure functional *)↵ClearAll[SoVMeasure];↵SoVMeasure[λ_][f_] :
7. [md] (* Result: normConst[λ] = (-1)^{λ-1} h^{λ-1} (λ-1)! / (2πi) exactly cancels the f(θ)-term↵coefficien


## 2026-07-14 13:17:48
**Input:**
```wolfram
(* SoV measure: normalized contour-integral functional.

Define μ[λ](u) = Product[1/(u - θ - h(k-1)), {k,1,λ}] (λ simple poles at u = θ, θ+h, ..., θ+h(λ-1)),
and the contour functional ⟨f⟩ = 2πi Σ_{k=1}^λ Residue[μ[λ](u) f(u), u=θ+h(k-1)] (contour just
encircling the λ poles, f assumed entire).

The residue at pole j (u=θ+h(j-1)) is f(θ+h(j-1)) (-1)^{λ-j} / [h^{λ-1} (j-1)! (λ-j)!], so the
f(θ) term (j=1) has coefficient 2πi (-1)^{λ-1} / [h^{λ-1} (λ-1)!]. Normalizing by

  normConst[λ] = (-1)^{λ-1} h^{λ-1} (λ-1)! / (2πi)

so that the f(θ) coefficient is exactly 1, the residue sum collapses (verified below, both
symbolically for generic f and on monomials u^n) to the alternating-binomial finite-difference form

  ⟨f⟩_normalized = Sum[(-1)^i Binomial[λ-1,i] f(θ+h i), {i,0,λ-1}]

i.e. exactly (-1)^{λ-1} times the (λ-1)-th forward finite difference of f at θ with step h. *)
```
**Output:**
```
(no output)
```

## 2026-07-14 13:17:48
**Input:**
```wolfram
(* measure with λ simple poles at u = θ+h(k-1), k=1,...,λ *)
ClearAll[μ];
μ[λ_][u_] := Product[1/(u - θ - h(k-1)), {k, 1, λ}];
```
**Output:**
```
(no output)
```

## 2026-07-14 13:17:49
**Input:**
```wolfram
(* raw contour functional = 2πi * sum of residues, computed directly via Residue *)
ClearAll[residueSum];
residueSum[λ_][f_] := 2 Pi I Sum[Residue[μ[λ][u] f[u], {u, θ + h(k-1)}], {k, 1, λ}];
Table[residueSum[λ][f], {λ,1,4}]
```
**Output:**
```
Out[69]= \left\{2\,\mathrm{i}\,\pi \,f\left(\theta \right),2\,\mathrm{i}\,\pi \,\left(\frac{f\left(h+\theta \right)}{h}-\frac{f\left(\theta \right)}{h}\right),2\,\mathrm{i}\,\pi \,\left(\frac{f\left(\theta \right)}{2\,h^2}-\frac{f\left(h+\theta \right)}{h^2}+\frac{f\left(2\,h+\theta \right)}{2\,h^2}\right),2\,\mathrm{i}\,\pi \,\left(-\frac{f\left(\theta \right)}{6\,h^3}+\frac{f\left(h+\theta \righ… [510 chars total]
```

## 2026-07-14 13:17:49
**Input:**
```wolfram
(* normalization constant fixing the coefficient of f[θ] to 1, and the claimed alternating-binomial
closed form; check they agree with the raw residue sum for symbolic f, λ=1,...,6 *)
ClearAll[normConst,normalizedSum,claimedClosed];
normConst[λ_] := (-1)^(λ-1) h^(λ-1) (λ-1)! / (2 Pi I);
normalizedSum[λ_][f_] := Simplify[normConst[λ] residueSum[λ][f]];
claimedClosed[λ_][f_] := Sum[(-1)^i Binomial[λ-1,i] f[θ+h i], {i,0,λ-1}];
Table[Simplify[normalizedSum[λ][f] - claimedClosed[λ][f]], {λ,1,6}]
```
**Output:**
```
Out[74]= \left\{0,0,0,0,0,0\right\}
```

## 2026-07-14 13:17:49
**Input:**
```wolfram
(* full check on a monomial basis f(u)=u^n, n=0,...,λ-1 (legitimate since the functional is
linear): normalized raw residue sum vs claimed closed form, λ=1,...,4 *)
Table[
  Simplify[
    (normConst[λ] * 2 Pi I Sum[Residue[μ[λ][u] u^n, {u, θ + h(k-1)}], {k,1,λ}])
      - Sum[(-1)^i Binomial[λ-1,i] (θ+h i)^n, {i,0,λ-1}]
  ],
  {λ,1,4},{n,0,λ-1}
] // Flatten // Union
```
**Output:**
```
Out[75]= \left\{0\right\}
```

## 2026-07-14 13:17:50
**Input:**
```wolfram
(* final definition: normalized SoV measure functional *)
ClearAll[SoVMeasure];
SoVMeasure[λ_][f_] := Sum[(-1)^i Binomial[λ-1,i] f[θ+h i], {i,0,λ-1}];
(* cross-check against the normalized raw contour-integral form, symbolic f, λ=1,...,6 *)
Table[Simplify[SoVMeasure[λ][f] - normalizedSum[λ][f]], {λ,1,6}]
```
**Output:**
```
Out[78]= \left\{0,0,0,0,0,0\right\}
```

## 2026-07-14 13:17:50
**Input:**
```wolfram
(* Result: normConst[λ] = (-1)^{λ-1} h^{λ-1} (λ-1)! / (2πi) exactly cancels the f(θ)-term
coefficient of the raw residue sum, and the normalized functional SoVMeasure[λ][f] =
Sum[(-1)^i Binomial[λ-1,i] f(θ+h i), {i,0,λ-1}] reproduces the normalized contour integral
identically for λ=1,...,6 (symbolic f) and on the monomial basis f(u)=u^n, n=0,...,λ-1, for
λ=1,...,4 (all differences {0}). So ⟨f⟩_normalized = f(θ) - (λ-1) f(θ+h) + Binomial[λ-1,2] f(θ+2h)
- ... + (-1)^{λ-1} f(θ+h(λ-1)), matching the anticipated alternating finite-difference pattern. *)
```
**Output:**
```
(no output)
```

## 2026-07-14 13:21:49
**Input:**
```wolfram
ClearAll[R,c];
R[λ_] := SoVMeasure[λ][Function[u, w^(u/h) Q[u]]] / (Q[θ] w^(θ/h));
c[s_,λ_] := (-1)^s Binomial[λ-1,s] w^s;
claimedR[λ_] := Sum[c[s,λ] (Q[θ+h s]/Q[θ]), {s,0,λ-1}];
Table[Simplify[R[λ] - claimedR[λ]], {λ,1,4}]
```
**Output:**
```
Out= {0, 0, 0, 0}
```
## 2026-07-14 13:22:04 — 📥 BULK INSERT 5 CELL(S) at positions 42–46
1. [markdown] (* Generic-Q algebraic check: apply SoVMeasure to u |-> w^(u/h) Q(u) for an ABSTRACT↵function Q (not
2. [code] ClearAll[R,c,claimedR];↵R[λ_] := SoVMeasure[λ][Function[u, w^(u/h) Q[u]]] / (Q[θ] w^(θ/h));↵(* claim
3. [markdown] (* Verification: R[λ] - claimedR[λ] must vanish identically for λ=1,...,4, with Q an↵opaque symbolic
4. [code] genericCheck = Table[Simplify[R[λ] - claimedR[λ]], {λ,1,4}]
5. [markdown] (* Result: genericCheck == {0,0,0,0} - the identity↵↵  R[λ] = <w^(u/h) Q(u)> / (Q(θ) w^(θ/h)) = Sum[


## 2026-07-14 13:22:05
**Input:**
```wolfram
ClearAll[R,c,claimedR];
R[λ_] := SoVMeasure[λ][Function[u, w^(u/h) Q[u]]] / (Q[θ] w^(θ/h));
(* claimed closed-form coefficients c_s(w,λ) and the resulting sum over SoV wavefunction values *)
c[s_,λ_] := (-1)^s Binomial[λ-1,s] w^s;
claimedR[λ_] := Sum[c[s,λ] (Q[θ+h s]/Q[θ]), {s,0,λ-1}];
R[1]
```
**Output:**
```
Out[83]= 1
```

## 2026-07-14 13:22:05
**Input:**
```wolfram
genericCheck = Table[Simplify[R[λ] - claimedR[λ]], {λ,1,4}]
```
**Output:**
```
Out[84]= \left\{0,0,0,0\right\}
```

## 2026-07-14 13:25:14
**Input:**
```wolfram
phiTest[λ_,w_] := Sum[(-1)^s Binomial[λ-1,s] w^s x[λ][{s}], {s,0,λ-1}]
phiCheckTest = Table[Simplify[First[phiTest[λ,w].Ψ[λ,M]] - SoVMeasure[λ][Function[u, w^(u/h) QFull[λ,M][u]]]/(QFull[λ,M][θ] w^(θ/h))], {λ,1,4},{M,0,λ}] // Flatten // Union
```
**Output:**
```
Out[1]= Null
Out[2]= {0}
```
## 2026-07-14 13:25:27 — 📥 BULK INSERT 4 CELL(S) at positions 47–50
1. [markdown] (* Identifying the covector state realizing R[λ] via the genuine SoV pairing.↵↵Combining Q[θ+hs]/Q[θ
2. [code] ClearAll[φ];↵φ[λ_,w_] := Sum[(-1)^s Binomial[λ-1,s] w^s x[λ][{s}], {s,0,λ-1}];↵↵phiCheck = Table[↵  
3. [markdown] (* Result: phiCheck == {0} identically for λ=1,...,4, all valid M=0,...,λ - the covector↵φ[λ,w] genu
4. [markdown] (* State identification: the covector↵↵  <φ_w| := Sum[(-1)^s Binomial[λ-1,s] w^s <x_s|, {s,0,λ-1}]↵↵


## 2026-07-14 13:25:28
**Input:**
```wolfram
ClearAll[φ];
φ[λ_,w_] := Sum[(-1)^s Binomial[λ-1,s] w^s x[λ][{s}], {s,0,λ-1}];

phiCheck = Table[
   Simplify[First[φ[λ,w].Ψ[λ,M]] - SoVMeasure[λ][Function[u, w^(u/h) QFull[λ,M][u]]]/(QFull[λ,M][θ] w^(θ/h))],
{λ,1,4},{M,0,λ}] // Flatten // Union
```
**Output:**
```
Out[87]= \left\{0\right\}
```
## 2026-07-14 13:32:14 — ✏️ EDIT CELL 35 [batch]
(* SoV measure: normalized contour-integral functional.

Define μ[λ](u) = Product[1/(u - θ - h(k-1)), {k,1,λ+1}] (λ+1 simple poles at u = θ, θ+h, ..., θ+hλ,
matching the FULL SoV covector basis x[λ][{…


## 2026-07-14 13:32:14
**Input:**
```wolfram
(* SoV measure: normalized contour-integral functional.

Define μ[λ](u) = Product[1/(u - θ - h(k-1)), {k,1,λ+1}] (λ+1 simple poles at u = θ, θ+h, ..., θ+hλ,
matching the FULL SoV covector basis x[λ][{s}], s=0,...,λ), and the contour functional
⟨f⟩ = 2πi Σ_{k=1}^{λ+1} Residue[μ[λ](u) f(u), u=θ+h(k-1)] (contour just encircling the λ+1
poles, f assumed entire).

The residue at pole j (u=θ+h(j-1)) is f(θ+h(j-1)) (-1)^{λ+1-j} / [h^λ (j-1)! (λ+1-j)!], so the
f(θ) term (j=1) has coefficient 2πi (-1)^λ / [h^λ λ!]. Normalizing by

  normConst[λ] = (-1)^λ h^λ λ! / (2πi)

so that the f(θ) coefficient is exactly 1, the residue sum collapses (verified below, both
symbolically for generic f and on monomials u^n) to the alternating-binomial finite-difference form

  ⟨f⟩_normalized = Sum[(-1)^i Binomial[λ,i] f(θ+h i), {i,0,λ}]

i.e. exactly (-1)^λ times the λ-th forward finite difference of f at θ with step h. *)
```
**Output:**
```
Out= (no output)
```
## 2026-07-14 13:32:14 — ✏️ EDIT CELL 36 [batch]
(* measure with λ+1 simple poles at u = θ+h(k-1), k=1,...,λ+1 - matches the FULL SoV
covector basis x[λ][{s}], s=0,...,λ (λ+1 basis states total) *)
ClearAll[μ];
μ[λ_][u_] := Product[1/(u - θ - h(k-1)…


## 2026-07-14 13:32:14
**Input:**
```wolfram
(* measure with λ+1 simple poles at u = θ+h(k-1), k=1,...,λ+1 - matches the FULL SoV
covector basis x[λ][{s}], s=0,...,λ (λ+1 basis states total) *)
ClearAll[μ];
μ[λ_][u_] := Product[1/(u - θ - h(k-1)), {k, 1, λ+1}];
```
**Output:**
```
Out= (no output)
```
## 2026-07-14 13:32:14 — ✏️ EDIT CELL 37 [batch]
(* raw contour functional = 2πi * sum of residues, computed directly via Residue *)
ClearAll[residueSum];
residueSum[λ_][f_] := 2 Pi I Sum[Residue[μ[λ][u] f[u], {u, θ + h(k-1)}], {k, 1, λ+1}];
Table[r…


## 2026-07-14 13:32:14
**Input:**
```wolfram
(* raw contour functional = 2πi * sum of residues, computed directly via Residue *)
ClearAll[residueSum];
residueSum[λ_][f_] := 2 Pi I Sum[Residue[μ[λ][u] f[u], {u, θ + h(k-1)}], {k, 1, λ+1}];
Table[residueSum[λ][f], {λ,1,4}]
```
**Output:**
```
Out= Out[92]= \left\{2\,\mathrm{i}\,\pi \,\left(\frac{f\left(h+\theta \right)}{h}-\frac{f\left(\theta \right)}{h}\right),2\,\mathrm{i}\,\pi \,\left(\frac{f\left(\theta \right)}{2\,h^2}-\frac{f\left(h+\theta \right)}{h^2}+\frac{f\left(2\,h+\theta \right)}{2\,h^2}\right),2\,\mathrm{i}\,\pi \,\left(-\frac{f\left(\theta \right)}{6\,h^3}+\frac{f\left(h+\theta \right)}{2\,h^3}-\frac{f\left(2\,h+\theta \right… [705 chars total]
```
## 2026-07-14 13:32:14 — ✏️ EDIT CELL 38 [batch]
(* normalization constant fixing the coefficient of f[θ] to 1, and the claimed alternating-binomial
closed form; check they agree with the raw residue sum for symbolic f, λ=1,...,6 *)
ClearAll[normCon…


## 2026-07-14 13:32:14
**Input:**
```wolfram
(* normalization constant fixing the coefficient of f[θ] to 1, and the claimed alternating-binomial
closed form; check they agree with the raw residue sum for symbolic f, λ=1,...,6 *)
ClearAll[normConst,normalizedSum,claimedClosed];
normConst[λ_] := (-1)^λ h^λ λ! / (2 Pi I);
normalizedSum[λ_][f_] := Simplify[normConst[λ] residueSum[λ][f]];
claimedClosed[λ_][f_] := Sum[(-1)^i Binomial[λ,i] f[θ+h i], {i,0,λ}];
Table[Simplify[normalizedSum[λ][f] - claimedClosed[λ][f]], {λ,1,6}]
```
**Output:**
```
Out= Out[97]= \left\{0,0,0,0,0,0\right\}
```
## 2026-07-14 13:32:14 — ✏️ EDIT CELL 39 [batch]
(* full check on a monomial basis f(u)=u^n, n=0,...,λ (legitimate since the functional is
linear): normalized raw residue sum vs claimed closed form, λ=1,...,4 *)
Table[
  Simplify[
    (normConst[λ] …


## 2026-07-14 13:32:15
**Input:**
```wolfram
(* full check on a monomial basis f(u)=u^n, n=0,...,λ (legitimate since the functional is
linear): normalized raw residue sum vs claimed closed form, λ=1,...,4 *)
Table[
  Simplify[
    (normConst[λ] * 2 Pi I Sum[Residue[μ[λ][u] u^n, {u, θ + h(k-1)}], {k,1,λ+1}])
      - Sum[(-1)^i Binomial[λ,i] (θ+h i)^n, {i,0,λ}]
  ],
  {λ,1,4},{n,0,λ}
] // Flatten // Union
```
**Output:**
```
Out= Out[98]= \left\{0\right\}
```
## 2026-07-14 13:32:15 — ✏️ EDIT CELL 40 [batch]
(* final definition: normalized SoV measure functional *)
ClearAll[SoVMeasure];
SoVMeasure[λ_][f_] := Sum[(-1)^i Binomial[λ,i] f[θ+h i], {i,0,λ}];
(* cross-check against the normalized raw contour-int…


## 2026-07-14 13:32:15
**Input:**
```wolfram
(* final definition: normalized SoV measure functional *)
ClearAll[SoVMeasure];
SoVMeasure[λ_][f_] := Sum[(-1)^i Binomial[λ,i] f[θ+h i], {i,0,λ}];
(* cross-check against the normalized raw contour-integral form, symbolic f, λ=1,...,6 *)
Table[Simplify[SoVMeasure[λ][f] - normalizedSum[λ][f]], {λ,1,6}]
```
**Output:**
```
Out= Out[101]= \left\{0,0,0,0,0,0\right\}
```
## 2026-07-14 13:32:15 — ✏️ EDIT CELL 41 [batch]
(* Result: normConst[λ] = (-1)^λ h^λ λ! / (2πi) exactly cancels the f(θ)-term
coefficient of the raw residue sum, and the normalized functional SoVMeasure[λ][f] =
Sum[(-1)^i Binomial[λ,i] f(θ+h i), {i…


## 2026-07-14 13:32:15
**Input:**
```wolfram
(* Result: normConst[λ] = (-1)^λ h^λ λ! / (2πi) exactly cancels the f(θ)-term
coefficient of the raw residue sum, and the normalized functional SoVMeasure[λ][f] =
Sum[(-1)^i Binomial[λ,i] f(θ+h i), {i,0,λ}] reproduces the normalized contour integral
identically for λ=1,...,6 (symbolic f) and on the monomial basis f(u)=u^n, n=0,...,λ, for
λ=1,...,4 (all differences {0}). So ⟨f⟩_normalized = f(θ) - λ f(θ+h) + Binomial[λ,2] f(θ+2h)
- ... + (-1)^λ f(θ+hλ), matching the anticipated alternating finite-difference pattern over
the FULL (λ+1)-point SoV basis. *)
```
**Output:**
```
Out= (no output)
```
## 2026-07-14 13:32:15 — ✏️ EDIT CELL 42 [batch]
(* Generic-Q algebraic check: apply SoVMeasure to u |-> w^(u/h) Q(u) for an ABSTRACT
function Q (not assumed to solve the TQ equation) and a free formal parameter w (new
symbol, unused elsewhere in th…

## 2026-07-14 13:32:15 — ✏️ EDIT CELL 43 [batch]
ClearAll[R,c,claimedR];
R[λ_] := SoVMeasure[λ][Function[u, w^(u/h) Q[u]]] / (Q[θ] w^(θ/h));
(* claimed closed-form coefficients c_s(w,λ) and the resulting sum over SoV wavefunction values *)
c[s_,λ_] …


## 2026-07-14 13:32:15
**Input:**
```wolfram
ClearAll[R,c,claimedR];
R[λ_] := SoVMeasure[λ][Function[u, w^(u/h) Q[u]]] / (Q[θ] w^(θ/h));
(* claimed closed-form coefficients c_s(w,λ) and the resulting sum over SoV wavefunction values *)
c[s_,λ_] := (-1)^s Binomial[λ,s] w^s;
claimedR[λ_] := Sum[c[s,λ] (Q[θ+h s]/Q[θ]), {s,0,λ}];
R[1]
```
**Output:**
```
Out= Out[106]= \frac{w^{-\frac{\theta }{h}}\,\left(Q\left(\theta \right)\,w^{\frac{\theta }{h}}-Q\left(h+\theta \right)\,w^{\frac{h+\theta }{h}}\right)}{Q\left(\theta \right)}
```
## 2026-07-14 13:32:15 — ✏️ EDIT CELL 44 [batch]
(* Verification: R[λ] - claimedR[λ] must vanish identically for λ=1,...,4, with Q an
opaque symbolic function (Q[u_] left undefined) and w a free formal parameter. *)

## 2026-07-14 13:32:15 — ✏️ EDIT CELL 45 [batch]
genericCheck = Table[Simplify[R[λ] - claimedR[λ]], {λ,1,4}]


## 2026-07-14 13:32:16
**Input:**
```wolfram
genericCheck = Table[Simplify[R[λ] - claimedR[λ]], {λ,1,4}]
```
**Output:**
```
Out= Out[107]= \left\{0,0,0,0\right\}
```
## 2026-07-14 13:32:16 — ✏️ EDIT CELL 46 [batch]
(* Result: genericCheck == {0,0,0,0} - the identity

  R[λ] = <w^(u/h) Q(u)> / (Q(θ) w^(θ/h)) = Sum[(-1)^s Binomial[λ,s] w^s Q[θ+h s]/Q[θ], {s,0,λ}]

holds EXACTLY for λ=1,...,4 with Q completely abst…

## 2026-07-14 13:32:16 — ✏️ EDIT CELL 47 [batch]
(* Identifying the covector state realizing R[λ] via the genuine SoV pairing.

Combining Q[θ+hs]/Q[θ] == x[λ][{s}].Ψ[λ,M] (cells 25-34, valid for Q=Q_M a genuine
TQ-solution) with the cell 42-46 resul…

## 2026-07-14 13:32:16 — ✏️ EDIT CELL 48 [batch]
ClearAll[φ];
φ[λ_,w_] := Sum[(-1)^s Binomial[λ,s] w^s x[λ][{s}], {s,0,λ}];

phiCheck = Table[
   Simplify[First[φ[λ,w].Ψ[λ,M]] - SoVMeasure[λ][Function[u, w^(u/h) QFull[λ,M][u]]]/(QFull[λ,M][θ] w^(θ/h…


## 2026-07-14 13:32:16
**Input:**
```wolfram
ClearAll[φ];
φ[λ_,w_] := Sum[(-1)^s Binomial[λ,s] w^s x[λ][{s}], {s,0,λ}];

phiCheck = Table[
   Simplify[First[φ[λ,w].Ψ[λ,M]] - SoVMeasure[λ][Function[u, w^(u/h) QFull[λ,M][u]]]/(QFull[λ,M][θ] w^(θ/h))],
{λ,1,4},{M,0,λ}] // Flatten // Union
```
**Output:**
```
Out= Out[110]= \left\{0\right\}
```
## 2026-07-14 13:32:16 — ✏️ EDIT CELL 49 [batch]
(* Result: phiCheck == {0} identically for λ=1,...,4, all valid M=0,...,λ - the covector
φ[λ,w] genuinely reproduces the SoV-measure ratio when paired with the true eigenvector
Ψ[λ,M], not just algebr…

## 2026-07-14 13:32:16 — ✏️ EDIT CELL 50 [batch]
(* State identification: the covector

  <φ_w| := Sum[(-1)^s Binomial[λ,s] w^s <x_s|, {s,0,λ}]

i.e. φ[λ,w] = Sum[(-1)^s Binomial[λ,s] w^s x[λ][{s}], {s,0,λ}], is the state satisfying

  <φ_w|Ψ[λ,M]> …


## 2026-07-14 13:34:58
**Input:**
```wolfram
Table[φ[λ,w] - (x[λ][{0}].MatrixExp[Ee[λ][2,1] ζ] /. ζ->-w), {λ,1,4}] // Simplify
```
**Output:**
```
Out= {{{0, 0}}, {{0, 0, 0}}, {{0, 0, 0, 0}}, {{0, 0, 0, 0, 0}}}
```

## 2026-07-14 13:35:04
**Input:**
```wolfram
Table[φ[λ,w] - (x[λ][{0}].MatrixExp[Ee[λ][1,2] ζ] /. ζ->-w), {λ,1,4}] // Simplify
Table[φ[λ,w] - (x[λ][{0}].MatrixExp[Ee[λ][1,2] ζ] /. ζ->w), {λ,1,4}] // Simplify
Table[φ[λ,w] - (x[λ][{0}].MatrixExp[Ee[λ][1,2] ζ] /. ζ->1/w), {λ,1,4}] // Simplify
Table[φ[λ,w] - (x[λ][{0}].MatrixExp[Ee[λ][1,2] ζ] /. ζ->-1/w), {λ,1,4}] // Simplify
Table[φ[λ,w] - (x[λ][{0}].MatrixExp[Ee[λ][2,1] ζ] /. ζ->w), {λ,1,4}] // Simplify
Table[φ[λ,w] - (x[λ][{0}].MatrixExp[Ee[λ][2,1] ζ] /. ζ->1/w), {λ,1,4}] // Simplify
Table[φ[λ,w] - (x[λ][{0}].MatrixExp[Ee[λ][2,1] ζ] /. ζ->-1/w), {λ,1,4}] // Simplify
```
**Output:**
```
Out[1]= {{{-w, 0}}, {{w^2/2, -w, 0}}, {{-1/6*w^3, w^2/2, -w, 0}}, {{w^4/24, -1/6*w^3, w^2/2, -w, 0}}}
Out[2]= {{{-w, 0}}, {{w^2/2, -w, 0}}, {{-1/6*w^3, w^2/2, -w, 0}}, {{w^4/24, -1/6*w^3, w^2/2, -w, 0}}}
Out[3]= {{{-w, 0}}, {{w^2/2, -w, 0}}, {{-1/6*w^3, w^2/2, -w, 0}}, {{w^4/24, -1/6*w^3, w^2/2, -w, 0}}}
Out[4]= {{{-w, 0}}, {{w^2/2, -w, 0}}, {{-1/6*w^3, w^2/2, -w, 0}}, {{w^4/24, -1/6*w^3, w^2/2, -w, 0}}}
Out[5]= {{{-2*w, 0}}, {{0, -2*w, 0}}, {{-1/3*w^3, 0, -2*w, 0}}, {{0, -1/3*w^3, 0, -2*w, 0}}}
Out[6]= {{{-((1 + w^2)/w), 0}}, {{(-1 + w^4)/(2*w^2), -((1 + w^2)/w), 0}}, {{-1/6*(1 + w^6)/w^3, (-1 + w^4)/(2*w^2), -((1 + w^2)/w), 0}}, {{(-1 + w^8)/(24*w^4), -1/6*(1 + w^6)/w^3, (-1 + w^4)/(2*w^2), -((1 + w^2)/w), 0}}}
Out[7]= {{{w^(-1) - w, 0}}, {{(-1 + w^4)/(2*w^2), w^(-1) - w, 0}}, {{-1/6*(-1 + w^6)/w^3, (-1 + w^4)/(2*w^2), w^(-1) - w, 0}}, {{(-1 + w^8)/(24*w^4), -1/6*(-1 + w^6)/w^3, (-1 + w^4)/(2*w^2), w^(-1) - w, 0}}}
```
## 2026-07-14 13:35:45 — 📥 BULK INSERT 2 CELL(S) at positions 51–52
1. [code] (* Pinning down the state identification: compare φ[λ,w] (cells 47-50) directly against↵x[λ][{0}].Ma
2. [markdown] (* Result: of the 8 (generator, ζ-relation) combinations tested (all against φ[λ,w],↵λ=1,...,4), onl


## 2026-07-14 13:35:46
**Input:**
```wolfram
(* Pinning down the state identification: compare φ[λ,w] (cells 47-50) directly against
x[λ][{0}].MatrixExp[Ee[λ][i,j] ζ] (the global-rotation covector of cells 12-14) for the
candidate generator (i,j)=(2,1) [natural guess from matching Binomial[λ,s]ζ^s vs
(-1)^s Binomial[λ,s]w^s] and the alternative (i,j)=(1,2), each against several candidate
relations ζ(w). λ=1,...,4 throughout. *)
stateCheck = <|
   "(2,1), ζ=-w"   -> (Table[φ[λ,w] - (x[λ][{0}].MatrixExp[Ee[λ][2,1] ζ] /. ζ -> -w), {λ,1,4}] // Simplify // Flatten // Union),
   "(2,1), ζ=w"    -> (Table[φ[λ,w] - (x[λ][{0}].MatrixExp[Ee[λ][2,1] ζ] /. ζ -> w), {λ,1,4}] // Simplify // Flatten // Union),
   "(2,1), ζ=1/w"  -> (Table[φ[λ,w] - (x[λ][{0}].MatrixExp[Ee[λ][2,1] ζ] /. ζ -> 1/w), {λ,1,4}] // Simplify // Flatten // Union),
   "(2,1), ζ=-1/w" -> (Table[φ[λ,w] - (x[λ][{0}].MatrixExp[Ee[λ][2,1] ζ] /. ζ -> -1/w), {λ,1,4}] // Simplify // Flatten // Union),
   "(1,2), ζ=-w"   -> (Table[φ[λ,w] - (x[λ][{0}].MatrixExp[Ee[λ][1,2] ζ] /. ζ -> -w), {λ,1,4}] // Simplify // Flatten // Union),
   "(1,2), ζ=w"    -> (Table[φ[λ,w] - (x[λ][{0}].MatrixExp[Ee[λ][1,2] ζ] /. ζ -> w), {λ,1,4}] // Simplify // Flatten // Union),
   "(1,2), ζ=1/w"  -> (Table[φ[λ,w] - (x[λ][{0}].MatrixExp[Ee[λ][1,2] ζ] /. ζ -> 1/w), {λ,1,4}] // Simplify // Flatten // Union),
   "(1,2), ζ=-1/w" -> (Table[φ[λ,w] - (x[λ][{0}].MatrixExp[Ee[λ][1,2] ζ] /. ζ -> -1/w), {λ,1,4}] // Simplify // Flatten // Union)
|>
```
**Output:**
```
Out[111]= \langle\!\lvert \text{(2,1), $\zeta$=-w}\to \left\{0\right\},\text{(2,1), $\zeta$=w}\to \left\{0,-2\,w,-\frac{w^3}{3}\right\},\text{(2,1), $\zeta$=1/w}\to \left\{0,-\frac{w^2+1}{w},\frac{w^4-1}{2\,w^2},-\frac{w^6+1}{6\,w^3},\frac{w^8-1}{24\,w^4}\right\},\text{(2,1), $\zeta$=-1/w}\to \left\{0,\frac{1}{w}-w,\frac{w^4-1}{2\,w^2},-\frac{w^6-1}{6\,w^3},\frac{w^8-1}{24\,w^4}\right\},\text{(1,2… [774 chars total]
```

## 2026-07-14 13:35:49
**Input:**
```wolfram
Table[x[λ][{0}].Ee[λ][1,2], {λ,1,4}]
```
**Output:**
```
Out= {{{0, 0}}, {{0, 0, 0}}, {{0, 0, 0, 0}}, {{0, 0, 0, 0, 0}}}
```

## 2026-07-14 13:52:18
**Input:**
```wolfram
MT = Transpose[B2[1,1][u]];
evals = Eigenvalues[MT];
evecs = Eigenvectors[MT];
{evals, evecs}//Simplify
```
**Output:**
```
Out[1]= (suppressed)
Out[2]= (suppressed)
Out[3]= (suppressed)
Out[4]= {{(u - θ1)*(u - θ2), -((h - u + θ1)*(u - θ2)), (u - θ1)*(-h + u - θ2), (h - u + θ1)*(h - u + θ2)}, {{0, 0, 0, 1}, {0, 1, 0, 0}, {0, h/(θ1 - θ2), 1, 0}, {1, 0, 0, 0}}}
```

## 2026-07-14 13:53:00
**Input:**
```wolfram
(* general solver: find correction coefficients a[k] for interior state (s1,s2) *)
buildXCorrected[λ1_,λ2_,s1_,s2_]:=Module[{kmax,a,ansatz,eigenvalue,eqSet,sol,diff},
  If[s2==0||s1==λ1,Return[x2[λ1,λ2][{s1,s2}]]];
  kmax=Min[λ1-s1,s2];
  eigenvalue=(u-θ1-s1 h)(u-θ2-s2 h);
  ansatz=x2[λ1,λ2][{s1,s2}]+Sum[a[k] x2[λ1,λ2][{s1+k,s2-k}],{k,1,kmax}];
  diff=Simplify[ansatz.B2[λ1,λ2][u]-eigenvalue ansatz][[1]];
  eqSet=Flatten[Table[Coefficient[diff[[j]],u,p]==0,{j,Length[diff]},{p,0,2}]];
  sol=Solve[eqSet,Table[a[k],{k,1,kmax}]]//Flatten;
  ansatz/.sol
];
(* test on λ1=1,λ2=1 interior state (0,1) *)
buildXCorrected[1,1,0,1]
```
**Output:**
```
Out= {{0, (h*(-(h^4/((h - θ1)*(h - θ2))) + (h^3*θ1)/((h - θ1)*(h - θ2))))/((h - θ1)*(-((h^3*θ1)/((h - θ1)*(h - θ2))) + (h^3*θ2)/((h - θ1)*(h - θ2)) + (h^2*(θ1 - θ2)*θ2)/((h - θ1)*(h - θ2)))), h/(h - θ2), 0}}
```

## 2026-07-14 13:53:05
**Input:**
```wolfram
FullSimplify[buildXCorrected[1,1,0,1]]
```
**Output:**
```
Out= {{0, h^2/((h - θ2)*(θ1 - θ2)), h/(h - θ2), 0}}
```

## 2026-07-14 13:53:13
**Input:**
```wolfram
(* verify residual exactly zero for λ1=1,λ2=2, all interior states *)
verify12=Table[
  Module[{X,eigenvalue,resid},
    X = buildXCorrected[1,2,s1,s2];
    eigenvalue=(u-θ1-s1 h)(u-θ2-s2 h);
    resid = Simplify[X.B2[1,2][u]-eigenvalue X];
    {s1,s2,Flatten[resid]//Union}
  ],
  {s1,0,1},{s2,0,2}];
verify12
```
**Output:**
```
Out= {{{0, 0, {0}}, {0, 1, {0}}, {0, 2, {0}}}, {{1, 0, {0}}, {1, 1, {0}}, {1, 2, {0}}}}
```

## 2026-07-14 13:53:23
**Input:**
```wolfram
verify22=Table[
  Module[{X,eigenvalue,resid},
    X = buildXCorrected[2,2,s1,s2];
    eigenvalue=(u-θ1-s1 h)(u-θ2-s2 h);
    resid = Simplify[X.B2[2,2][u]-eigenvalue X];
    {s1,s2,Flatten[resid]//Union}
  ],
  {s1,0,2},{s2,0,2}];
verify22
```
**Output:**
```
Out= {{{0, 0, {0}}, {0, 1, {0}}, {0, 2, {0}}}, {{1, 0, {0}}, {1, 1, {0}}, {1, 2, {0}}}, {{2, 0, {0}}, {2, 1, {0}}, {2, 2, {0}}}}
```
