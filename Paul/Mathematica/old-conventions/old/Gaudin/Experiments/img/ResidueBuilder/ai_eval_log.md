# AI Action Log — ResidueBuilder

## 2026-07-12 04:46:12 — ✏️ EDIT CELL 1
# ResidueBuilder — residue-evaluation experiments

Experiments folder companion to `Clean/FastQ.wb`. This notebook keeps the exploratory work behind the `res1`/`res0` closed-form residue formulas for …

## 2026-07-12 04:46:15 — 🗑️ DELETE MARKDOWN CELL 19
*(empty)*

## 2026-07-12 04:46:15 — 🗑️ DELETE CODE CELL 22
*(empty)*

## 2026-07-12 04:46:26 — 📥 BULK INSERT 1 CELL(S) at positions 18–18
1. [markdown] ## Generalized wrapper: sum residues over an arbitrary polynomial's coefficients↵↵`doRes0`/`doRes1`

## 2026-07-12 04:46:33 — ✏️ EDIT CELL 21
### Next steps

Need a Baxter equation and solver for the `[S1,S1,0,0] x [S2,S2,0,0]` representation.

## 2026-07-12 04:46:57 — ▶️ RUN CELLS 2–21
Cell 2: ✓ — (no output)
Cell 4: ✓ — Out[9]= 0
Cell 5: ✓ — Out[11]= 0
Cell 6: ✓ — Out[30]= \frac{7\,u^{20}}{2}-\frac{13\,u^{18}}{10}-\frac{u^{17}}{5}-\frac{u^{16}}{10}+\frac{3\,u^{15}}{2}-2\,u^{14}-\frac{7\,u^{13}}{3}+\frac{2\,u^{12}}{5}-5\,u^{11}-\frac{11\,u^{10}}{2}+6\,u^9+\frac{2\,u^8}{5}+\frac{15\,u^7}{7}-\frac{20\,u^6}{9}+\frac{13\,u^5}{9}-\frac{4\,u^4}{7}-\frac{20\,u^3}{7}-\frac{17\,u^2}{4}-4\,u+\frac{7}{5}
Cell 8: ✓ — Out[48]= 21
Cell 9: ✓ — Out[49]= \left\{0.030851,\mathrm{Null}\right\}
Cell 10: ✓ — Out[50]= \left\{0.021656,\mathrm{Null}\right\}
Cell 11: ✓ — Out[52]= \left\{0.077055,\mathrm{Null}\right\}
Cell 12: ✓ — Out[54]= \left\{0,0\right\}
Cell 14: ✓ — Out[56]= \begin{array}{|c|c|c|c|}\hline\text{S1=S2} & \text{direct Residue (s)} & \text{closed-form Laguerre (s)} & \text{fastResidue (s)}\\\hline8 & 0.011856 & 0.008005 & 0.004026\\\hline16 & 0.020836 & 0.003989 & 0.003732\\\hline24 & 0.024365 & 0.008816 & 0.005524\\\hline32 & 0.032679 & 0.019103 & 0.00955\\\hline40 & 0.020031 & 0.013879 & 0.008478\\\hline\end{array}
Cell 16: ✓ — Out[58]= \begin{array}{|c|c|c|c|}\hline\text{S1=S2} & \text{direct Residue (s)} & \text{closed-form Laguerre (s)} & \text{fastResidue (s)}\\\hline50 & 0.052777 & 0.025125 & 0.017995\\\hline100 & 0.109625 & 0.040383 & 0.044044\\\hline200 & 0.737406 & 0.092564 & 0.249711\\\hline400 & 5.342746 & 0.246564 & 1.75029\\\hline\end{array}
Cell 19: ✓ — 2 kernel messages: General::munfl: ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------… [also: General::munfl] result: Out[60]= \left\{5.393482,\mathrm{Null}\right\} | Out[61]= \left\{0.318719,\mathrm{Null}\right\} | Out[62]= \left\{0.004109,\mathrm{Null}\right\}
Cell 20: ✓ — (no output)


## 2026-07-12 04:47:11
**Input:**
```wolfram
{Length[mismatches1], Length[mismatches0], ValueQ[doRes0], ValueQ[doRes1], Length[DownValues[res1]]}
```
**Output:**
```
Out= {0, 0, False, False, 290}
```

## 2026-07-12 04:47:19
**Input:**
```wolfram
{Length[DownValues[doRes0]], Length[DownValues[doRes1]], doRes0[8, 8, poly20] // AbsoluteTiming}
```
**Output:**
```
Out= {1, 1, {0.000309, 15/7 - (20*(8 + w))/9 + (13*(72 + 16*w + w^2))/18 - (2*(720 + 216*w + 24*w^2 + w^3))/21 - (5*(7920 + 2880*w + 432*w^2 + 32*w^3 + w^4))/42 - (17*(95040 + 39600*w + 7200*w^2 + 720*w^3 + 40*w^4 + w^5))/480 + (-1235520 - 570240*w - 118800*w^2 - 14400*w^3 - 1080*w^4 - 48*w^5 - w^6)/180 + (17297280 + 8648640*w + 1995840*w^2 + 277200*w^3 + 25200*w^4 + 1512*w^5 + 56*w^6 + w^7)/3600}}
```
## 2026-07-12 04:48:14 — ✏️ EDIT CELL 6
(* poly20: an arbitrary degree-20 test polynomial in u with random rational coefficients,
   used only as a stand-in workload to validate/benchmark the residue-evaluation
   approaches below (n plays …


## 2026-07-12 04:48:14
**Input:**
```wolfram
(* poly20: an arbitrary degree-20 test polynomial in u with random rational coefficients,
   used only as a stand-in workload to validate/benchmark the residue-evaluation
   approaches below (n plays the role of a generic monomial degree, not a physical
   quantum number) *)
SeedRandom[42];
poly20 = Sum[RandomInteger[{-20, 20}]/RandomInteger[{1, 10}] u^k, {k, 0, 20}]
```
**Output:**
```
Out= 7/5 - 4*u - (17*u^2)/4 - (20*u^3)/7 - (4*u^4)/7 + (13*u^5)/9 - (20*u^6)/9 + (15*u^7)/7 + (2*u^8)/5 + 6*u^9 - (11*u^10)/2 - 5*u^11 + (2*u^12)/5 - (7*u^13)/3 - 2*u^14 + (3*u^15)/2 - u^16/10 - u^17/5 - (13*u^18)/10 + (7*u^20)/2
```
## 2026-07-12 04:50:38 — ↕️ MOVE MARKDOWN CELL 18 → position 19
## Generalized wrapper: sum residues over an arbitrary polynomial's coefficients

`doRes0`/`doRes1`

## 2026-07-12 04:50:47 — 📥 BULK INSERT 1 CELL(S) at positions 19–19
1. [markdown] ### What these timings show↵↵- At small pole order (S1=S2=8) all three evaluation strategies are wit

