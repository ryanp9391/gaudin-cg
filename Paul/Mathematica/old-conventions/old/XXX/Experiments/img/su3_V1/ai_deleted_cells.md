## Deleted Cell 18 — 2026-07-17 13:28:29
**Kind:** code
```wolfram
(*Baxter equation for S1, S2 reps*)
ClearAll[Bax]
Bax[S1_,S2_]:=f[u+h] - τg[S1,S2][1,1][u]/ν1[S1,S2][u] f[u] + τg[S1,S2][2,1][u]/ν1[S1,S2][u]/ν1[S1,S2][u-h]f[u-h] + χ3 qdet[S1,S2][3,1][u]/ν1[S1,S2][u]/ν1[S1,S2][u-h]/ν1[S1,S2][u-2h]f[u-2h]
```

