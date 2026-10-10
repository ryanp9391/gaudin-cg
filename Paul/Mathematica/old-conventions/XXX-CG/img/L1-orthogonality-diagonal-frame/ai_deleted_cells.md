## Deleted Cell 14 — 2026-09-23 06:00:41
**Kind:** code
```wolfram
(* ---- 5. Diagnostic: which reference points work with THIS vacL/vacR? Old notebook's scan: all (x, y) in a box with
   <Psi_A|Psi_A> Psi_A(x) D_A(y) / det N^{AA} the same for every state A. Labels (a;b,c) are shifts of theta in units of h. ---- *)
ClearAll[boxO, psiLabO, dualLabO, matchPairsO];
boxO = Tuples[Range[lamO[[3]] - 1, lamO[[1]] + 1], 3];
psiLabO[B_, {a_, b_, c_}] := QDO[B, 1, thO + a hO] (QDO[B, 1, thO + b hO] QDO[B, 2, thO + (c - 1) hO] - QDO[B, 2, thO + b hO] QDO[B, 1, thO + (c - 1) hO]);
dualLabO[A_, {a_, b_, c_}] := QUO[A, 3, thO + a hO] (QUO[A, 3, thO + b hO] QUO[A, 2, thO + (c + 1) hO] - QUO[A, 2, thO + b hO] QUO[A, 3, thO + (c + 1) hO]);

matchPairsO = Module[{nS = Length[statesO], g, P, D, pa, db, res = {}},
    g = N[Table[ovMatO[[A, A]] / detMatO[[A, A]], {A, nS}], 30];
    P = N[Table[psiLabO[A, x], {x, boxO}, {A, nS}], 30];
    D = N[Table[dualLabO[A, y], {y, boxO}, {A, nS}], 30];
    (* a pair matches iff the two state-profiles, each normalized to its Omega entry, multiply to 1 componentwise *)
    pa = Table[If[Min[Abs[P[[k]]]] < 10^-25, Nothing, {boxO[[k]], g P[[k]] / (g[[iOmO]] P[[k, iOmO]])}], {k, Length[boxO]}];
    db = Table[If[Min[Abs[D[[k]]]] < 10^-25, Nothing, {boxO[[k]], D[[k]] / D[[k, iOmO]]}], {k, Length[boxO]}];
    Do[If[Max[Abs[x[[2]] y[[2]] - 1]] < 10^-20, AppendTo[res, {x[[1]], y[[1]]}]], {x, pa}, {y, db}];
    res
];

Grid[
    Prepend[
        If[matchPairsO === {}, {{"none", "none"}}, {Row[#[[1]], ","], Row[#[[2]], ","]} & /@ matchPairsO],
        {"x = (a;b,c) for Psi_B", "y = (a;b,c) for c_A"}
    ],
    Frame -> All, Alignment -> Left
]
```

## Deleted Cell 8 — 2026-09-23 06:00:41
**Kind:** code
```wolfram
(3 h^4 (z[1] - z[2]) (z[1] - z[3])^2 (z[2] - z[3]))
```

