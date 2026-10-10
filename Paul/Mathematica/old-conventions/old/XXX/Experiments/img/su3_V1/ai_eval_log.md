# AI Action Log — su3_V1

## 2026-07-17 08:58:20 — 🌟 KERNEL START


---

## 2026-07-17 13:14:33 — 🔄 KERNEL RESTART

## 2026-07-17 13:14:42 — ▶️ RUN CELLS 1–6
Cell 1: ✓ — (no output)
Cell 2: ✓ — (no output)
Cell 3: ✓ — (no output)
Cell 4: ✓ — Out[14]= \left\{0\right\}
Cell 5: ✓ — (no output)
Cell 6: ✓ — (no output)

## 2026-07-17 13:14:50 — 📥 BULK INSERT 1 CELL(S) at positions 7–7
1. [code] (*Numeric parameters (shared across Section B). θ,h copied from the su(2) Clean notebook;↵  three DI


## 2026-07-17 13:14:50
**Input:**
```wolfram
(*Numeric parameters (shared across Section B). θ,h copied from the su(2) Clean notebook;
  three DISTINCT unit-modulus twist eigenvalues z[1],z[2],z[3] (generic → non-degenerate;
  no ∏z=1 constraint since gl(3) is not unimodular-normalized). χ1,χ2,χ3 and G were defined
  in Section A as symmetric functions of z[i], so they become numeric automatically.*)
θ1 = 1/3;
θ2 = 1/7;
h = 1;
z[1] = N[Exp[I Zeta[3]], 30];
z[2] = N[Exp[I Zeta[5]], 30];
z[3] = N[Exp[I Zeta[7]], 30];
{θ1, θ2, h, z[1], z[2], z[3], χ1, χ2, χ3}
```
**Output:**
```
Out[27]= \left\{\frac{1}{3},\frac{1}{7},1,0.360439875147226537266677650596+0.932782448593374716484835982071\,\mathrm{i},0.508867380990499828163321372267+0.860844927129079093188644393616\,\mathrm{i},0.533257880589599372595832345093+0.845952736734794526483381618453\,\mathrm{i},1.40256513672732573802583136796+2.63958011245724833615686199414\,\mathrm{i},-1.6733239859016290991182472189+2.47680418504351… [511 chars total]
```
## 2026-07-17 13:15:34 — ▶️ RUN CELLS 8–14
Cell 8: ✓ — (no output)
Cell 9: ✓ — (no output)
Cell 10: ✓ — Out[37]= \left\{0\right\}
Cell 11: ✓ — (no output)
Cell 12: ✓ — (no output)
Cell 13: ✓ — (no output)
Cell 14: ✓ — Out[48]= \left\{0,\left(0.\times 10^{-29}+0.\times 10^{-29}\,\mathrm{i}\right)\,\left(u-\frac{36}{7}\right)\,\left(u-\frac{10}{3}\right)\,\left(u-\frac{1}{3}\right)\,\left(u-\frac{1}{7}\right),\left(0.\times 10^{-29}+0.\times 10^{-29}\,\mathrm{i}\right)\,\left(u-\frac{29}{7}\right)\,\left(u-\frac{10}{3}\right)\,\left(u-\frac{1}{3}\right)\,\left(u-\frac{1}{7}\right),\left(0.\times 10^{-29}+0.\times… [1248 chars total]

## 2026-07-17 13:15:48 — 📥 BULK INSERT 1 CELL(S) at positions 15–15
1. [code] (*verify numeric parameters took hold before building Section B on top of them*)↵{Precision[χ1], N[χ


## 2026-07-17 13:15:48
**Input:**
```wolfram
(*verify numeric parameters took hold before building Section B on top of them*)
{Precision[χ1], N[χ3], Chip = χ1 - (z[1] + z[2] + z[3])} /. Chip -> Chop[Chip]
```
**Output:**
```
Out[49]= \left\{29.998415624999733,-0.9944145978415286-0.10554433949611333\,\mathrm{i},0\right\}
```
## 2026-07-17 13:15:54 — 📥 BULK INSERT 1 CELL(S) at positions 16–16
1. [code] (*single source of truth for the Yangian weight of T11, mirroring su(2)'s aFun→ν1 dedup*)↵ClearAll[ν


## 2026-07-17 13:15:55
**Input:**
```wolfram
(*single source of truth for the Yangian weight of T11, mirroring su(2)'s aFun→ν1 dedup*)
ClearAll[ν1];
ν1[S1_, S2_][u_] := ν[S1, S2][1][u];   (* = (u-θ1-h S1)(u-θ2-h S2) *)
ν1[1, 1][u]
```
**Output:**
```
Out[52]= \left(u-\frac{4}{3}\right)\,\left(u-\frac{8}{7}\right)
```
## 2026-07-17 13:16:05 — ✏️ EDIT CELL 18
(*3rd-order Baxter/TQ equation (reference form). Diagonal-twist transfer matrices τ1=tg[·][1,1],
  τ2=tg[·][2,1]; qdetT already carries χ3, so NO explicit χ3 here (would double-count).
  A1=τ1/ν1[u], …


## 2026-07-17 13:16:05
**Input:**
```wolfram
(*3rd-order Baxter/TQ equation (reference form). Diagonal-twist transfer matrices τ1=tg[·][1,1],
  τ2=tg[·][2,1]; qdetT already carries χ3, so NO explicit χ3 here (would double-count).
  A1=τ1/ν1[u], A2=τ2/(ν1[u]ν1[u-h]), A3=qdetT/(ν1[u]ν1[u-h]ν1[u-2h]).*)
ClearAll[Bax];
Bax[S1_,S2_] :=
  f[u + h]
   - tg[S1,S2][1,1][u]/ν1[S1,S2][u] f[u]
   + tg[S1,S2][2,1][u]/(ν1[S1,S2][u] ν1[S1,S2][u - h]) f[u - h]
   - qdetT[S1,S2][u]/(ν1[S1,S2][u] ν1[S1,S2][u - h] ν1[S1,S2][u - 2 h]) f[u - 2 h];
Head[Bax[1, 1]]
```
**Output:**
```
Out= List
```
## 2026-07-17 13:16:19 — ▶️ RUN CELL 18
**In [18]:** `(*Baxter equation for S1, S2 reps*)↵ClearAll[Bax]↵Bax[S1_,S2_]:=f[u+h] - τg[S1,S2][1,1][u]/ν1[S1,S2][u] f[u] + τg[S1,S2][2,1][u]/ν1[S1,S2][u]/ν1[S1,S2][u-h]f[u-h] + χ3 qdet[S1,S2][3,1][u]/ν1[S1,S2][u]`
**Out:** (no output)


## 2026-07-17 13:16:29
**Input:**
```wolfram
Head[Bax[1,1]]
```
**Output:**
```
Out= Plus
```
## 2026-07-17 13:17:17 — ✏️ EDIT CELL 18
(*3rd-order Baxter/TQ equation (reference form). Diagonal-twist transfer matrices τ1=tg[·][1,1],
  τ2=tg[·][2,1]; qdetT already carries χ3, so NO explicit χ3 here (would double-count).
  A1=τ1/ν1[u], …


## 2026-07-17 13:17:17
**Input:**
```wolfram
(*3rd-order Baxter/TQ equation (reference form). Diagonal-twist transfer matrices τ1=tg[·][1,1],
  τ2=tg[·][2,1]; qdetT already carries χ3, so NO explicit χ3 here (would double-count).
  A1=τ1/ν1[u], A2=τ2/(ν1[u]ν1[u-h]), A3=qdetT/(ν1[u]ν1[u-h]ν1[u-2h]).*)
ClearAll[Bax];
Bax[S1_,S2_] :=
  f[u + h]
   - tg[S1,S2][1,1][u]/ν1[S1,S2][u] f[u]
   + tg[S1,S2][2,1][u]/(ν1[S1,S2][u] ν1[S1,S2][u - h]) f[u - h]
   - qdetT[S1,S2][u]/(ν1[S1,S2][u] ν1[S1,S2][u - h] ν1[S1,S2][u - 2 h]) f[u - 2 h];
Head[Bax[1, 1]]
```
**Output:**
```
Out= List
```

## 2026-07-17 13:17:30
**Input:**
```wolfram
Head[Bax[1,1]]
```
**Output:**
```
Out= List
```

## 2026-07-17 13:17:59
**Input:**
```wolfram
tg[1,1][1,1][u]
tg[1,1][2,1][u]
ν1[1,1][u]
qdetT[1,1][u]
Bax[1,1]
```
**Output:**
```
Out[1]= {{(0.3604398751472265372666776505955500912937645232622850024938`29.707347828414424 + 0.9327824485933747164848359820708679294578496404879183613065`30.120295363558725*I)*(-4/3 + u)*(-8/7 + u) + (1.0421252615800992007591537173601929152707025546163127305558`29.867404925683804 + 1.7067976638638736196720260120693068356773456131376779774675`30.081667041926103*I)*(-1/3 + u)*(-1/7 + u), 0, 0, 0, 0, 0, 0, 0, 0}, {0, (0.5088673809904998281633213722670686199090953963505721183772`29.857119610800275 + 0.8608449271290790931886443936161642458151438806992171086196`30.085439922389305*I)*(-8/7 + u)*(-1/3 + u) + (0.3604398751472265372666776505955500912937645232622850024938`29.707347828414424 + 0.9327824485933747164848359820708679294578496404879183613065`30.120295363558725*I)*(-4/3 + u)*(-1/7 + u) + (0.5332578805895993725958323450931242953616071582657406121786`29.877452280089955 + 0.8459527367347945264833816184531425898622017324384608688479`30.077861097574722*I)*(-1/3 + u)*(-1/7 + u), 0, 0.5088673809904998281633213722670686199090953963505721183772`29.857119610800275 + 0.8608449271290790931886443936161642458151438806992171086196`30.085439922389305*I, 0, 0, 0, 0, 0}, {0, 0, (0.5332578805895993725958323450931242953616071582657406121786`29.877452280089955 + 0.8459527367347945264833816184531425898622017324384608688479`30.077861097574722*I)*(-8/7 + u)*(-1/3 + u) + (0.3604398751472265372666776505955500912937645232622850024938`29.707347828414424 + 0.9327824485933747164848359820708679294578496404879183613065`30.120295363558725*I)*(-4/3 + u)*(-1/7 + u) + (0.5088673809904998281633213722670686199090953963505721183772`29.857119610800275 + 0.8608449271290790931886443936161642458151438806992171086196`30.085439922389305*I)*(-1/3 + u)*(-1/7 + u), 0, 0, 0, 0.5332578805895993725958323450931242953616071582657406121786`29.877452280089955 + 0.8459527367347945264833816184531425898622017324384608688479`30.077861097574722*I, 0, 0}, {0, 0.3604398751472265372666776505955500912937645232622850024938`29.707347828414424 + 0.9327824485933747164848359820708679294578496404879183613065`30.120295363558725*I, 0, (0.3604398751472265372666776505955500912937645232622850024938`29.707347828414424 + 0.9327824485933747164848359820708679294578496404879183613065`30.120295363558725*I)*(-8/7 + u)*(-1/3 + u) + (0.5088673809904998281633213722670686199090953963505721183772`29.857119610800275 + 0.8608449271290790931886443936161642458151438806992171086196`30.085439922389305*I)*(-4/3 + u)*(-1/7 + u) + (0.5332578805895993725958323450931242953616071582657406121786`29.877452280089955 + 0.8459527367347945264833816184531425898622017324384608688479`30.077861097574722*I)*(-1/3 + u)*(-1/7 + u), 0, 0, 0, 0, 0}, {0, 0, 0, 0, (0.5088673809904998281633213722670686199090953963505721183772`29.857119610800275 + 0.8608449271290790931886443936161642458151438806992171086196`30.085439922389305*I)*(-4/3 + u)*(-8/7 + u) + (0.8936977557368259098625099956886743866553716815280256146724`29.800675669499228 + 1.7787351853281692429682176005240105193200513729263792301544`30.099596298137637*I)*(-1/3 + u)*(-1/7 + u), 0, 0, 0, 0}, {0, 0, 0, 0, 0, (0.5332578805895993725958323450931242953616071582657406121786`29.877452280089955 + 0.8459527367347945264833816184531425898622017324384608688479`30.077861097574722*I)*(-8/7 + u)*(-1/3 + u) + (0.5088673809904998281633213722670686199090953963505721183772`29.857119610800275 + 0.8608449271290790931886443936161642458151438806992171086196`30.085439922389305*I)*(-4/3 + u)*(-1/7 + u) + (0.3604398751472265372666776505955500912937645232622850024938`29.707347828414424 + 0.9327824485933747164848359820708679294578496404879183613065`30.120295363558725*I)*(-1/3 + u)*(-1/7 + u), 0, 0.5332578805895993725958323450931242953616071582657406121786`29.877452280089955 + 0.8459527367347945264833816184531425898622017324384608688479`30.077861097574722*I, 0}, {0, 0, 0.3604398751472265372666776505955500912937645232622850024938`29.707347828414424 + 0.9327824485933747164848359820708679294578496404879183613065`30.120295363558725*I, 0, 0, 0, (0.3604398751472265372666776505955500912937645232622850024938`29.70734782841
[output truncated — 5755 chars total]
Out[2]= {{(0.6195649437159357064208469289648537218014899499180425647488`29.517241328555162 - 0.784945399705144727340209210312420467181148591100079349728`29.61999394061785*I)*(8/7 - u)^2 + (0.6195649437159357064208469289648537218014899499180425647488`29.517241328555162 - 0.784945399705144727340209210312420467181148591100079349728`29.61999394061785*I)*(1/3 - u)*(7/3 - u) - (0.5968824612647582951218549041615314257475333055057271216695`29.499496999309883 - 0.8023286903972238692144500919746258807419595316683941809551`29.627960503873116*I)*(8/7 - u)*(-8/7 + u) - (1.2164474049806940015427018331263851475490232554237696864186`29.508443989177355 - 1.5872740901023685965546593022870463479231081227684735306831`29.62400258000614*I)*(-4/3 + u)^2*(-8/7 + u)^2 - (0.5968824612647582951218549041615314257475333055057271216695`29.499496999309883 - 0.8023286903972238692144500919746258807419595316683941809551`29.627960503873116*I)*(7/3 - u)*(-1/3 + u) - (0.4568765809209350975755453857775508578129657999197682896188`29.371151331691877 - 0.8895300949411415184237534370449845150193994358697773187929`29.66051308050147*I)*(-4/3 + u)*(-8/7 + u)*(-1/3 + u)*(-1/7 + u), 0, 0, 0, 0, 0, 0, 0, 0}, {0, (0.6195649437159357064208469289648537218014899499180425647488`29.517241328555162 - 0.784945399705144727340209210312420467181148591100079349728`29.61999394061785*I)*(1/7 - u)*(15/7 - u) - (0.4568765809209350975755453857775508578129657999197682896188`29.371151331691877 - 0.8895300949411415184237534370449845150193994358697773187931`29.66051308050147*I)*(4/3 - u)*(-1/3 + u) - (0.4568765809209350975755453857775508578129657999197682896188`29.371151331691877 - 0.8895300949411415184237534370449845150193994358697773187929`29.66051308050147*I)*(-4/3 + u)*(-8/7 + u)^2*(-1/3 + u) - (0.5968824612647582951218549041615314257475333055057271216695`29.499496999309883 - 0.8023286903972238692144500919746258807419595316683941809551`29.627960503873116*I)*(8/7 - u)*(-1/7 + u) - (0.6195649437159357064208469289648537218014899499180425647488`29.517241328555162 - 0.784945399705144727340209210312420467181148591100079349728`29.61999394061785*I)*(-15/7 + u)*(-4/3 + u)^2*(-1/7 + u) - (0.5968824612647582951218549041615314257475333055057271216698`29.499496999309883 - 0.8023286903972238692144500919746258807419595316683941809551`29.627960503873116*I)*(-4/3 + u)^2*(-8/7 + u)*(-1/7 + u), 0, (0.4568765809209350975755453857775508578129657999197682896188`29.371151331691877 - 0.8895300949411415184237534370449845150193994358697773187931`29.66051308050147*I)*(8/7 - u) + (0.6195649437159357064208469289648537218014899499180425647488`29.517241328555162 - 0.784945399705144727340209210312420467181148591100079349728`29.61999394061785*I)*(1/7 - u)*(4/3 - u) - (0.4568765809209350975755453857775508578129657999197682896188`29.371151331691877 - 0.8895300949411415184237534370449845150193994358697773187929`29.66051308050147*I)*(-4/3 + u)*(-8/7 + u) - (0.4568765809209350975755453857775508578129657999197682896188`29.371151331691877 - 0.8895300949411415184237534370449845150193994358697773187931`29.66051308050147*I)*(-1/3 + u) - (0.6195649437159357064208469289648537218014899499180425647488`29.517241328555162 - 0.784945399705144727340209210312420467181148591100079349728`29.61999394061785*I)*(-4/3 + u)*(-1/7 + u), 0, 0, 0, 0, 0}, {0, 0, (0.6195649437159357064208469289648537218014899499180425647488`29.517241328555162 - 0.784945399705144727340209210312420467181148591100079349728`29.61999394061785*I)*(1/7 - u)*(8/7 - u) - (0.5968824612647582951218549041615314257475333055057271216695`29.499496999309883 - 0.8023286903972238692144500919746258807419595316683941809551`29.627960503873116*I)*(15/7 - u)*(-1/7 + u) - (0.5968824612647582951218549041615314257475333055057271216698`29.499496999309883 - 0.8023286903972238692144500919746258807419595316683941809551`29.627960503873116*I)*(-15/7 + u)*(-4/3 + u)^2*(-1/7 + u) - (0.6195649437159357064208469289648537218014899499180425647488`29.517241328555162 - 0.784945399705144727340209210312420467181148591100079349728`29.61999394061785*I)*(-4/3 + u)^2*(-8/7 + u)*(-1/7 + u) - (0.45687658092093509757554
[output truncated — 13678 chars total]
Out[3]= (-4/3 + u)*(-8/7 + u)
Out[4]= (-0.9944145978415286577806456150645332908393735038958268340201`29.44245325110677 - 0.1055443394961133295075060113382859913751619075649164082601`28.46832070644821*I)*(-10/3 + u)*(-22/7 + u)*(-4/3 + u)*(-8/7 + u)*(-1/3 + u)*(-1/7 + u)
Out[5]= {{((0.9944145978415286577806456150645332908393735038958268340201`29.44245325110677 + 0.1055443394961133295075060113382859913751619075649164082601`28.46832070644821*I)*(-1/3 + u)*(-1/7 + u)*f[-2 + u])/((-7/3 + u)*(-15/7 + u)) + (((0.6195649437159357064208469289648537218014899499180425647488`29.517241328555162 - 0.784945399705144727340209210312420467181148591100079349728`29.61999394061785*I)*(8/7 - u)^2 + (0.6195649437159357064208469289648537218014899499180425647488`29.517241328555162 - 0.784945399705144727340209210312420467181148591100079349728`29.61999394061785*I)*(1/3 - u)*(7/3 - u) - (0.5968824612647582951218549041615314257475333055057271216695`29.499496999309883 - 0.8023286903972238692144500919746258807419595316683941809551`29.627960503873116*I)*(8/7 - u)*(-8/7 + u) - (1.2164474049806940015427018331263851475490232554237696864186`29.508443989177355 - 1.5872740901023685965546593022870463479231081227684735306831`29.62400258000614*I)*(-4/3 + u)^2*(-8/7 + u)^2 - (0.5968824612647582951218549041615314257475333055057271216695`29.499496999309883 - 0.8023286903972238692144500919746258807419595316683941809551`29.627960503873116*I)*(7/3 - u)*(-1/3 + u) - (0.4568765809209350975755453857775508578129657999197682896188`29.371151331691877 - 0.8895300949411415184237534370449845150193994358697773187929`29.66051308050147*I)*(-4/3 + u)*(-8/7 + u)*(-1/3 + u)*(-1/7 + u))*f[-1 + u])/((-7/3 + u)*(-15/7 + u)*(-4/3 + u)*(-8/7 + u)) - (((0.3604398751472265372666776505955500912937645232622850024938`29.707347828414424 + 0.9327824485933747164848359820708679294578496404879183613065`30.120295363558725*I)*(-4/3 + u)*(-8/7 + u) + (1.0421252615800992007591537173601929152707025546163127305558`29.867404925683804 + 1.7067976638638736196720260120693068356773456131376779774675`30.081667041926103*I)*(-1/3 + u)*(-1/7 + u))*f[u])/((-4/3 + u)*(-8/7 + u)) + f[1 + u], ((0.9944145978415286577806456150645332908393735038958268340201`29.44245325110677 + 0.1055443394961133295075060113382859913751619075649164082601`28.46832070644821*I)*(-1/3 + u)*(-1/7 + u)*f[-2 + u])/((-7/3 + u)*(-15/7 + u)) + f[1 + u], ((0.9944145978415286577806456150645332908393735038958268340201`29.44245325110677 + 0.1055443394961133295075060113382859913751619075649164082601`28.46832070644821*I)*(-1/3 + u)*(-1/7 + u)*f[-2 + u])/((-7/3 + u)*(-15/7 + u)) + f[1 + u], ((0.9944145978415286577806456150645332908393735038958268340201`29.44245325110677 + 0.1055443394961133295075060113382859913751619075649164082601`28.46832070644821*I)*(-1/3 + u)*(-1/7 + u)*f[-2 + u])/((-7/3 + u)*(-15/7 + u)) + f[1 + u], ((0.9944145978415286577806456150645332908393735038958268340201`29.44245325110677 + 0.1055443394961133295075060113382859913751619075649164082601`28.46832070644821*I)*(-1/3 + u)*(-1/7 + u)*f[-2 + u])/((-7/3 + u)*(-15/7 + u)) + f[1 + u], ((0.9944145978415286577806456150645332908393735038958268340201`29.44245325110677 + 0.1055443394961133295075060113382859913751619075649164082601`28.46832070644821*I)*(-1/3 + u)*(-1/7 + u)*f[-2 + u])/((-7/3 + u)*(-15/7 + u)) + f[1 + u], ((0.9944145978415286577806456150645332908393735038958268340201`29.44245325110677 + 0.1055443394961133295075060113382859913751619075649164082601`28.46832070644821*I)*(-1/3 + u)*(-1/7 + u)*f[-2 + u])/((-7/3 + u)*(-15/7 + u)) + f[1 + u], ((0.9944145978415286577806456150645332908393735038958268340201`29.44245325110677 + 0.1055443394961133295075060113382859913751619075649164082601`28.46832070644821*I)*(-1/3 + u)*(-1/7 + u)*f[-2 + u])/((-7/3 + u)*(-15/7 + u)) + f[1 + u], ((0.9944145978415286577806456150645332908393735038958268340201`29.44245325110677 + 0.1055443394961133295075060113382859913751619075649164082601`28.46832070644821*I)*(-1/3 + u)*(-1/7 + u)*f[-2 + u])/((-7/3 + u)*(-15/7 + u)) + f[1 + u]}, {((0.9944145978415286577806456150645332908393735038958268340201`29.44245325110677 + 0.1055443394961133295075060113382859913751619075649164082601`28.46832070644821*I)*(-1/3 + u)*(-1/7 + u)*f[-2 + u])/((-7/3 + u)*(-15/7 + u)) + f[1 + u], ((0.9944145978415286577806456150645332908393735038958268340201`29.44245325110677 + 0.1055443394961133295075060113382859913
[output truncated — 39494 chars total]
```

## 2026-07-17 13:19:10
**Input:**
```wolfram
{Head[tg[1,1][1,1][u]], Dimensions[tg[1,1][1,1][u]], Head[qdetT[1,1][u]], Dimensions[Bax[1,1]], Head[Bax[1,1]]}
```
**Output:**
```
Out= {List, {9, 9}, Times, {9, 9}, List}
```
## 2026-07-17 13:28:05 — ✏️ EDIT CELL 18
(*3rd-order Baxter/TQ equation (reference form). Diagonal-twist transfer matrices τ1=tg[·][1,1],
  τ2=tg[·][2,1]; qdetT already carries χ3, so NO explicit χ3 here (would double-count).
  A1=τ1/ν1[u], …

## 2026-07-17 13:28:29 — 🗑️ DELETE CODE CELL 18
(*Baxter equation for S1, S2 reps*)
ClearAll[Bax]
Bax[S1_,S2_]:=f[u+h] - τg[S1,S2][1,1][u]/ν1[S1,S2][u] f[u] + τg[S1,S2][2,1][u]/ν1[S1,S2][u]/ν1[S1,S2][u-h]f[u-h] + χ3 qdet[S1,S2][3,1][u]/ν1[S1,S2][u]…

## 2026-07-17 13:28:38 — 📥 BULK INSERT 1 CELL(S) at positions 18–18
1. [code] (*3rd-order Baxter/TQ equation (reference form). Diagonal-twist transfer matrices τ1=tg[·][1,1],↵  τ


## 2026-07-17 13:28:38
**Input:**
```wolfram
(*3rd-order Baxter/TQ equation (reference form). Diagonal-twist transfer matrices τ1=tg[·][1,1],
  τ2=tg[·][2,1]; qdetT already carries χ3, so NO explicit χ3 here (would double-count).
  A1=τ1/ν1[u], A2=τ2/(ν1[u]ν1[u-h]), A3=qdetT/(ν1[u]ν1[u-h]ν1[u-2h]).*)
ClearAll[Bax];
Bax[S1_, S2_] :=
  f[u + h]
   - tg[S1, S2][1, 1][u]/ν1[S1, S2][u] f[u]
   + tg[S1, S2][2, 1][u]/(ν1[S1, S2][u] ν1[S1, S2][u - h]) f[u - h]
   - qdetT[S1, S2][u]/(ν1[S1, S2][u] ν1[S1, S2][u - h] ν1[S1, S2][u - 2 h]) f[u - 2 h];
Head[Bax[1, 1]]
```
**Output:**
```
Out[57]= \mathrm{List}
```
## 2026-07-17 13:32:18 — 📥 BULK INSERT 3 CELL(S) at positions 19–21
1. [code] (*coefficient matrices of tau1=tg[.][1,1] (deg 2 in u) and tau2=tg[.][2,1] (deg 4 in u),↵  extracted
2. [code] (*Cartan operators (gl(3) weight generators H1=E11-E22, H2=E22-E33) acting on the two-site space*)↵C
3. [code] (*hard-assert check: leading coefficients of tau1,tau2 are chi1*Id, chi2*Id (quantum-determinant↵  n


## 2026-07-17 13:32:19
**Input:**
```wolfram
(*coefficient matrices of tau1=tg[.][1,1] (deg 2 in u) and tau2=tg[.][2,1] (deg 4 in u),
  extracted by CoefficientList in a dummy variable uu so tg[..][..][uu] is treated as a polynomial in uu.
  c1={C0,C1,C2} with tau1(u)=Sum_k C_k u^k (index k+1 in the list); c2={D0,...,D4} similarly for tau2.*)
ClearAll[tgCoeffs];
tgCoeffs[S1_, S2_] := tgCoeffs[S1, S2] = Module[{uu, m1, m2, cl1, cl2},
   m1 = tg[S1, S2][1, 1][uu];
   m2 = tg[S1, S2][2, 1][uu];
   cl1 = Map[PadRight[CoefficientList[#, uu], 3] &, m1, {2}];
   cl2 = Map[PadRight[CoefficientList[#, uu], 5] &, m2, {2}];
   {Table[cl1[[All, All, k]], {k, 1, 3}], Table[cl2[[All, All, k]], {k, 1, 5}]}
];
```
**Output:**
```
(no output)
```

## 2026-07-17 13:32:19
**Input:**
```wolfram
(*Cartan operators (gl(3) weight generators H1=E11-E22, H2=E22-E33) acting on the two-site space*)
ClearAll[H1op, H2op];
H1op[S1_, S2_] := EE[S1, S2][1, 1] - EE[S1, S2][2, 2];
H2op[S1_, S2_] := EE[S1, S2][2, 2] - EE[S1, S2][3, 3];
```
**Output:**
```
(no output)
```

## 2026-07-17 13:32:19
**Input:**
```wolfram
(*hard-assert check: leading coefficients of tau1,tau2 are chi1*Id, chi2*Id (quantum-determinant
  normalization), and the full commuting family {C0,C1,C2,D0,...,D4,H1op,H2op} pairwise-commutes*)
Module[{S1 = 2, S2 = 1, d, c1, c2, C2, D4, mats, comm},
  d = dim[{S1, 0, 0}] dim[{S2, 0, 0}];
  {c1, c2} = tgCoeffs[S1, S2];
  C2 = c1[[3]]; D4 = c2[[5]];
  (* leading coeffs: τ1 ~ χ1 u^2, τ2 ~ χ2 u^4 *)
  If[Max[Abs[Flatten[C2 - χ1 IdentityMatrix[d]]]] > 10^-13, Print["C2 FAIL"]; Abort[]];
  If[Max[Abs[Flatten[D4 - χ2 IdentityMatrix[d]]]] > 10^-13, Print["D4 FAIL"]; Abort[]];
  (* full commuting family pairwise-commutes *)
  mats = Join[c1, c2, {H1op[S1, S2], H2op[S1, S2]}];
  comm = Max[Table[Max[Abs[Flatten[mats[[i]].mats[[j]] - mats[[j]].mats[[i]]]]],
      {i, Length[mats]}, {j, Length[mats]}]];
  If[comm > 10^-12, Print["COMMUTATION FAIL: ", comm]; Abort[]];
  {"C2=χ1 Id OK", "D4=χ2 Id OK", "maxCommutator" -> comm}
]
```
**Output:**
```
Out[63]= \left\{\text{C2=$\chi$1 Id OK},\text{D4=$\chi$2 Id OK},\text{maxCommutator}\to 0.\times 10^{-26}\right\}
```

## 2026-07-17 13:32:30
**Input:**
```wolfram
Module[{S1 = 1, S2 = 1, d, c1, c2, C2, D4, mats, comm},
  d = dim[{S1, 0, 0}] dim[{S2, 0, 0}];
  {c1, c2} = tgCoeffs[S1, S2];
  C2 = c1[[3]]; D4 = c2[[5]];
  If[Max[Abs[Flatten[C2 - χ1 IdentityMatrix[d]]]] > 10^-13, Print["C2 FAIL"]; Abort[]];
  If[Max[Abs[Flatten[D4 - χ2 IdentityMatrix[d]]]] > 10^-13, Print["D4 FAIL"]; Abort[]];
  mats = Join[c1, c2, {H1op[S1, S2], H2op[S1, S2]}];
  comm = Max[Table[Max[Abs[Flatten[mats[[i]].mats[[j]] - mats[[j]].mats[[i]]]]],
      {i, Length[mats]}, {j, Length[mats]}]];
  If[comm > 10^-12, Print["COMMUTATION FAIL: ", comm]; Abort[]];
  {"C2=χ1 Id OK", "D4=χ2 Id OK", "maxCommutator" -> comm}
]
```
**Output:**
```
Out= {"C2=χ1 Id OK", "D4=χ2 Id OK", "maxCommutator" -> 0``26.15841266270934}
```

## 2026-07-17 13:32:48
**Input:**
```wolfram
Module[{S1 = 2, S2 = 2, d, c1, c2, C2, D4, mats, comm},
  d = dim[{S1, 0, 0}] dim[{S2, 0, 0}];
  {c1, c2} = tgCoeffs[S1, S2];
  C2 = c1[[3]]; D4 = c2[[5]];
  If[Max[Abs[Flatten[C2 - χ1 IdentityMatrix[d]]]] > 10^-13, Print["C2 FAIL"]; Abort[]];
  If[Max[Abs[Flatten[D4 - χ2 IdentityMatrix[d]]]] > 10^-13, Print["D4 FAIL"]; Abort[]];
  mats = Join[c1, c2, {H1op[S1, S2], H2op[S1, S2]}];
  comm = Max[Table[Max[Abs[Flatten[mats[[i]].mats[[j]] - mats[[j]].mats[[i]]]]],
      {i, Length[mats]}, {j, Length[mats]}]];
  If[comm > 10^-12, Print["COMMUTATION FAIL: ", comm]; Abort[]];
  {"C2=χ1 Id OK", "D4=χ2 Id OK", "maxCommutator" -> comm}
]
```
**Output:**
```
Out= {"C2=χ1 Id OK", "D4=χ2 Id OK", "maxCommutator" -> 0``25.6545786594317}
```
## 2026-07-17 13:42:22 — 📥 BULK INSERT 2 CELL(S) at positions 22–23
1. [code] (*combinatorial predictor for gl(3) weight multiplicities in [S1,0,0]\[CircleTimes][S2,0,0]:↵  recov
2. [code] (*hard-assert sum-rule check: total weight count over all (H1,H2) must equal dim[S1]*dim[S2]*)↵Modul


## 2026-07-17 13:42:22
**Input:**
```wolfram
(*combinatorial predictor for gl(3) weight multiplicities in [S1,0,0]\[CircleTimes][S2,0,0]:
  recover (n1,n2,n3) from (H1,H2) via n1+n2+n3=S1+S2, then count splittings
  (a1,a2,a3) of [S1,0,0] with a componentwise complement (n_i-a_i>=0) landing in [S2,0,0]*)
ClearAll[wtMult];
wtMult[S1_, S2_][H1_, H2_] := Module[{Ntot = S1 + S2, n1, n2, n3},
   n3 = (Ntot - H1 - 2 H2)/3; n2 = n3 + H2; n1 = n3 + H1 + H2;
   If[! (IntegerQ[n1] && IntegerQ[n2] && IntegerQ[n3] && n1 >= 0 && n2 >= 0 && n3 >= 0),
     Return[0]];
   Length@Select[
     Flatten[Table[{a1, a2, S1 - a1 - a2}, {a1, 0, S1}, {a2, 0, S1 - a1}], 1],
     (n1 - #[[1]] >= 0 && n2 - #[[2]] >= 0 && n3 - #[[3]] >= 0) &]
];
```
**Output:**
```
(no output)
```

## 2026-07-17 13:42:23
**Input:**
```wolfram
(*hard-assert sum-rule check: total weight count over all (H1,H2) must equal dim[S1]*dim[S2]*)
Module[{S1 = 2, S2 = 1, d, tot},
  d = dim[{S1, 0, 0}] dim[{S2, 0, 0}];
  tot = Total@Flatten@Table[
      wtMult[S1, S2][H1, H2], {H1, -3 (S1 + S2), 3 (S1 + S2)}, {H2, -3 (S1 + S2), 3 (S1 + S2)}];
  If[tot != d, Print["SUM-RULE FAIL: ", tot, " != ", d]; Abort[]];
  {"total weights" -> tot, "dim" -> d,
   "spot (2,1)@(H1,H2)=(1,1)" -> wtMult[2, 1][1, 1]}
]
```
**Output:**
```
Out[66]= \left\{\text{total weights}\to 18,\text{dim}\to 18,\text{spot (2,1)@(H1,H2)=(1,1)}\to 2\right\}
```
## 2026-07-17 13:42:30 — 📥 BULK INSERT 1 CELL(S) at positions 24–24
1. [code] (*guard the combinatorics on a second rep: (S1,S2)=(2,2), d=36*)\nModule[{S1 = 2, S2 = 2, d, tot},\n


## 2026-07-17 13:42:30
**Input:**
```wolfram
(*guard the combinatorics on a second rep: (S1,S2)=(2,2), d=36*)\nModule[{S1 = 2, S2 = 2, d, tot},\n  d = dim[{S1, 0, 0}] dim[{S2, 0, 0}];\n  tot = Total@Flatten@Table[\n      wtMult[S1, S2][H1, H2], {H1, -3 (S1 + S2), 3 (S1 + S2)}, {H2, -3 (S1 + S2), 3 (S1 + S2)}];\n  If[tot != d, Print[\"SUM-RULE FAIL: \", tot, \" != \", d]; Abort[]];\n  {\"total weights\" -> tot, \"dim\" -> d}\n]
```
**Output:**
```
Out[67]= \left\{\text{total weights}\to 36,\text{dim}\to 36\right\}
```
## 2026-07-17 13:48:14 — 📥 BULK INSERT 3 CELL(S) at positions 25–27
1. [code] ClearAll[TauEigensystem3];↵TauEigensystem3::rankdef = "Generic combination failed to resolve the `1`
2. [code] ClearAll[Psi, tau1, tau2];↵Psi::bnds = "No state at (H1,H2,n)=(`1`,`2`,`3`) for (S1,S2)=(`4`,`5`).";
3. [code] Module[{S1 = 2, S2 = 1, es, worst},↵  es = TauEigensystem3[S1, S2];↵  worst = Max@KeyValueMap[↵


## 2026-07-17 13:48:14
**Input:**
```wolfram
ClearAll[TauEigensystem3];
TauEigensystem3::rankdef = "Generic combination failed to resolve the `1`-dim eigenspace for (S1,S2)=(`2`,`3`) at every tried coefficient set.";
TauEigensystem3::notgenuine = "Eigenvector #`1` for (S1,S2)=(`2`,`3`) failed the common-eigenvector residual check (residual=`4`).";
TauEigensystem3::noninteger = "Weight (H1 or H2)=`1` not near-integer for eigenvector #`2` at (S1,S2)=(`3`,`4`).";
TauEigensystem3::degnull = "Sector (H1,H2)=(`1`,`2`) for (S1,S2)=(`3`,`4`) has `5` state(s) but wtMult predicts `6`.";
TauEigensystem3[S1_, S2_] := TauEigensystem3[S1, S2] = Module[
   {d, c1, c2, mats, nontriv, rSets, rset, Kmix, vals, vecs, rayleigh, recs, byHW, table},
   d = dim[{S1, 0, 0}] dim[{S2, 0, 0}];
   {c1, c2} = tgCoeffs[S1, S2];
   (* full family for extraction; drop the two scalar leading coeffs (c1[[3]],c2[[5]]) from the
      splitting combination since they carry no state info *)
   mats = Join[c1, c2, {H1op[S1, S2], H2op[S1, S2]}];
   nontriv = Join[c1[[1 ;; 2]], c2[[1 ;; 4]], {H1op[S1, S2], H2op[S1, S2]}]; (* 8 matrices *)
   rSets = Table[N[Table[Prime[100 i + 7 j]^(1/3), {j, Length[nontriv]}], 30], {i, 1, 8}];
   rset = SelectFirst[rSets,
     Module[{ee, vv}, {ee, vv} = Eigensystem[#.nontriv];
        MatrixRank[vv] == d && Length[Union[Chop[ee, 10^-12]]] == d] &,
     Missing["NoGenericR"]];
   If[MissingQ[rset], Message[TauEigensystem3::rankdef, d, S1, S2]; Abort[]];
   Kmix = rset . nontriv;
   {vals, vecs} = Eigensystem[N[Kmix, 30]];
   rayleigh[m_, v_] := (v . (m . v))/(v . v);
   recs = Table[
     Module[{v = vecs[[i]], quotients, res, H1v, H2v, t1c, t2c},
       quotients = rayleigh[#, v] & /@ mats;
       res = Max[Abs[#.v - rayleigh[#, v] v]] & /@ mats // Max;
       If[res > 10^-13, Message[TauEigensystem3::notgenuine, i, S1, S2, res]; Abort[]];
       t1c = rayleigh[#, v] & /@ c1;           (* {C0,C1,C2} quotients *)
       t2c = rayleigh[#, v] & /@ c2;           (* {D0..D4} quotients *)
       H1v = rayleigh[H1op[S1, S2], v]; H2v = rayleigh[H2op[S1, S2], v];
       If[Abs[H1v - Round[Re[H1v]]] > 10^-10, Message[TauEigensystem3::noninteger, H1v, i, S1, S2]; Abort[]];
       If[Abs[H2v - Round[Re[H2v]]] > 10^-10, Message[TauEigensystem3::noninteger, H2v, i, S1, S2]; Abort[]];
       <|"H1" -> Round[Re[H1v]], "H2" -> Round[Re[H2v]], "t1c" -> t1c, "t2c" -> t2c, "vec" -> v|>],
     {i, 1, d}];
   byHW = GroupBy[recs, {#["H1"], #["H2"]} &];
   table = Association @ Flatten @ KeyValueMap[
      Function[{hw, group},
        Module[{sorted, predicted, last},
          predicted = wtMult[S1, S2][hw[[1]], hw[[2]]];
          (* order n by (Re,Im) of the C0 Rayleigh quotient (energy-like scalar) *)
          sorted = SortBy[group, {Re[#["t1c"][[1]]], Im[#["t1c"][[1]]]} &];
          If[Length[sorted] != predicted,
            Message[TauEigensystem3::degnull, hw[[1]], hw[[2]], S1, S2, Length[sorted], predicted]; Abort[]];
          MapIndexed[
            Function[{rec, idx}, Module[{n = First[idx] - 1, vv = rec["vec"], lc},
              lc = SelectFirst[vv, Abs[#] > 10^-13 Max[Abs[vv]] &];
              {hw[[1]], hw[[2]], n} -> <|"tau1c" -> rec["t1c"], "tau2c" -> rec["t2c"], "Psi" -> vv/lc|>]],
            sorted]]],
      byHW];
   If[Length[table] != d, Message[TauEigensystem3::degnull, "total", "-", S1, S2, Length[table], d]; Abort[]];
   table
];
```
**Output:**
```
(no output)
```

## 2026-07-17 13:48:14
**Input:**
```wolfram
ClearAll[Psi, tau1, tau2];
Psi::bnds = "No state at (H1,H2,n)=(`1`,`2`,`3`) for (S1,S2)=(`4`,`5`).";
Psi[S1_, S2_][H1_, H2_, n_] := Lookup[TauEigensystem3[S1, S2], Key[{H1, H2, n}],
    (Message[Psi::bnds, H1, H2, n, S1, S2]; Missing["NoState", {H1, H2, n}])]["Psi"];
tau1[S1_, S2_][H1_, H2_, n_][u_] := TauEigensystem3[S1, S2][{H1, H2, n}, "tau1c"] . u^Range[0, 2];
tau2[S1_, S2_][H1_, H2_, n_][u_] := TauEigensystem3[S1, S2][{H1, H2, n}, "tau2c"] . u^Range[0, 4];
```
**Output:**
```
(no output)
```

## 2026-07-17 13:48:15
**Input:**
```wolfram
Module[{S1 = 2, S2 = 1, es, worst},
  es = TauEigensystem3[S1, S2];
  worst = Max@KeyValueMap[
     Function[{k, rec},
       Module[{H1 = k[[1]], H2 = k[[2]], n = k[[3]], v = rec["Psi"], uu = 17/5},
         Max[
          Max[Abs[(tg[S1, S2][1, 1][uu] . v) - tau1[S1, S2][H1, H2, n][uu] v]],
          Max[Abs[(tg[S1, S2][2, 1][uu] . v) - tau2[S1, S2][H1, H2, n][uu] v]]]]],
     es];
  If[worst > 10^-11, Print["EIGEN RESIDUAL FAIL: ", worst]; Abort[]];
  {"nstates" -> Length[es], "worstResidual@u=17/5" -> worst}
]
```
**Output:**
```
Out[79]= \left\{\text{nstates}\to 18,\text{worstResidual@u=17/5}\to 1.359292769213167893921303\times 10^{-27}\right\}
```

## 2026-07-17 14:02:18
**Input:**
```wolfram
Module[{S1 = 1, S2 = 1, es, worst},
  es = TauEigensystem3[S1, S2];
  worst = Max@KeyValueMap[
     Function[{k, rec},
       Module[{H1 = k[[1]], H2 = k[[2]], n = k[[3]], v = rec["Psi"], uu = 17/5},
         Max[
          Max[Abs[(tg[S1, S2][1, 1][uu] . v) - tau1[S1, S2][H1, H2, n][uu] v]],
          Max[Abs[(tg[S1, S2][2, 1][uu] . v) - tau2[S1, S2][H1, H2, n][uu] v]]]]],
     es];
  If[worst > 10^-11, Print["EIGEN RESIDUAL FAIL: ", worst]; Abort[]];
  {"nstates" -> Length[es], "worstResidual@u=17/5" -> worst, "expected" -> dim[{S1,0,0}] dim[{S2,0,0}]}
]
```
**Output:**
```
Out= {"nstates" -> 9, "worstResidual@u=17/5" -> 4.351882429728977936066513775178630089543`24.81639060791746*^-38, "expected" -> 9}
```

## 2026-07-17 14:02:21
**Input:**
```wolfram
Module[{S1 = 2, S2 = 2, es, worst},
  es = TauEigensystem3[S1, S2];
  worst = Max@KeyValueMap[
     Function[{k, rec},
       Module[{H1 = k[[1]], H2 = k[[2]], n = k[[3]], v = rec["Psi"], uu = 17/5},
         Max[
          Max[Abs[(tg[S1, S2][1, 1][uu] . v) - tau1[S1, S2][H1, H2, n][uu] v]],
          Max[Abs[(tg[S1, S2][2, 1][uu] . v) - tau2[S1, S2][H1, H2, n][uu] v]]]]],
     es];
  If[worst > 10^-11, Print["EIGEN RESIDUAL FAIL: ", worst]; Abort[]];
  {"nstates" -> Length[es], "worstResidual@u=17/5" -> worst, "expected" -> dim[{S1,0,0}] dim[{S2,0,0}]}
]
```
**Output:**
```
Out= {"nstates" -> 36, "worstResidual@u=17/5" -> 1.798311988555260034005821300982471455924115020616055`25.012439748753433*^-26, "expected" -> 36}
```

## 2026-07-17 14:02:27
**Input:**
```wolfram
Keys[TauEigensystem3[2,1]]
```
**Output:**
```
Out= {{-1, -1, 0}, {-1, -1, 1}, {1, -2, 0}, {1, -2, 1}, {1, 1, 0}, {1, 1, 1}, {-2, 1, 0}, {-2, 1, 1}, {0, 0, 0}, {0, 0, 1}, {0, 0, 2}, {2, -1, 0}, {2, -1, 1}, {-1, 2, 0}, {-1, 2, 1}, {0, -3, 0}, {3, 0, 0}, {-3, 3, 0}}
```

## 2026-07-17 14:02:28
**Input:**
```wolfram
SortBy[Keys[TauEigensystem3[2,1]], {#[[1]],#[[2]],#[[3]]}&] // TableForm
```
**Output:**
```
Out= TableForm[{{-3, 3, 0}, {-2, 1, 0}, {-2, 1, 1}, {-1, -1, 0}, {-1, -1, 1}, {-1, 2, 0}, {-1, 2, 1}, {0, -3, 0}, {0, 0, 0}, {0, 0, 1}, {0, 0, 2}, {1, -2, 0}, {1, -2, 1}, {1, 1, 0}, {1, 1, 1}, {2, -1, 0}, {2, -1, 1}, {3, 0, 0}}]
```
## 2026-07-17 14:11:36 — ✏️ EDIT CELL 25
ClearAll[TauEigensystem3];
TauEigensystem3::rankdef = "Generic combination failed to resolve the `1`-dim eigenspace for (S1,S2)=(`2`,`3`) at every tried coefficient set.";
TauEigensystem3::notgenuine …


## 2026-07-17 14:11:36
**Input:**
```wolfram
ClearAll[TauEigensystem3];
TauEigensystem3::rankdef = "Generic combination failed to resolve the `1`-dim eigenspace for (S1,S2)=(`2`,`3`) at every tried coefficient set.";
TauEigensystem3::notgenuine = "Eigenvector #`1` for (S1,S2)=(`2`,`3`) failed the common-eigenvector residual check (residual=`4`).";
TauEigensystem3::noninteger = "Weight (H1 or H2)=`1` not near-integer for eigenvector #`2` at (S1,S2)=(`3`,`4`).";
TauEigensystem3::degnull = "Sector (H1,H2)=(`1`,`2`) for (S1,S2)=(`3`,`4`) has `5` state(s) but wtMult predicts `6`.";
TauEigensystem3[S1_, S2_] := TauEigensystem3[S1, S2] = Module[
   {d, c1, c2, mats, nontriv, rSets, rset, Kmix, vals, vecs, rayleigh, recs, byHW, table},
   d = dim[{S1, 0, 0}] dim[{S2, 0, 0}];
   {c1, c2} = tgCoeffs[S1, S2];
   (* full family for extraction; drop the two scalar leading coeffs (c1[[3]],c2[[5]]) from the
      splitting combination since they carry no state info *)
   mats = Join[c1, c2, {H1op[S1, S2], H2op[S1, S2]}];
   nontriv = Join[c1[[1 ;; 2]], c2[[1 ;; 4]], {H1op[S1, S2], H2op[S1, S2]}]; (* 8 matrices *)
   rSets = Table[N[Table[Prime[100 i + 7 j]^(1/3), {j, Length[nontriv]}], 30], {i, 1, 8}];
   rset = SelectFirst[rSets,
     Module[{ee, vv}, {ee, vv} = Eigensystem[#.nontriv];
        MatrixRank[vv] == d && Length[Union[Chop[ee, 10^-12]]] == d] &,
     Missing["NoGenericR"]];
   If[MissingQ[rset], Message[TauEigensystem3::rankdef, d, S1, S2]; Abort[]];
   Kmix = rset . nontriv;
   {vals, vecs} = Eigensystem[N[Kmix, 30]];
   rayleigh[m_, v_] := (v . (m . v))/(v . v);
   recs = Table[
     Module[{v = vecs[[i]], quotients, res, H1v, H2v, t1c, t2c},
       quotients = rayleigh[#, v] & /@ mats;
       res = Max[Abs[#.v - rayleigh[#, v] v]] & /@ mats // Max;
       If[res > 10^-13, Message[TauEigensystem3::notgenuine, i, S1, S2, res]; Abort[]];
       t1c = rayleigh[#, v] & /@ c1;           (* {C0,C1,C2} quotients *)
       t2c = rayleigh[#, v] & /@ c2;           (* {D0..D4} quotients *)
       H1v = rayleigh[H1op[S1, S2], v]; H2v = rayleigh[H2op[S1, S2], v];
       If[Abs[H1v - Round[Re[H1v]]] > 10^-10, Message[TauEigensystem3::noninteger, H1v, i, S1, S2]; Abort[]];
       If[Abs[H2v - Round[Re[H2v]]] > 10^-10, Message[TauEigensystem3::noninteger, H2v, i, S1, S2]; Abort[]];
       <|"H1" -> Round[Re[H1v]], "H2" -> Round[Re[H2v]], "t1c" -> t1c, "t2c" -> t2c, "vec" -> v|>],
     {i, 1, d}];
   byHW = GroupBy[recs, {#["H1"], #["H2"]} &];
   table = Association @ Flatten @ KeyValueMap[
      Function[{hw, group},
        Module[{sorted, predicted, last},
          predicted = wtMult[S1, S2][hw[[1]], hw[[2]]];
          (* order n by the full real+imaginary parts of ALL tau1/tau2 Rayleigh quotients:
             within a fixed (H1,H2) sector, states are separated by their tau1/tau2 charges
             since the generic combination Kmix has simple spectrum there, so this tuple gives
             a deterministic total order (a single C0 coefficient is not guaranteed distinct) *)
          sorted = SortBy[group, Join[Re[#["t1c"]], Im[#["t1c"]], Re[#["t2c"]], Im[#["t2c"]]] &];
          If[Length[sorted] != predicted,
            Message[TauEigensystem3::degnull, hw[[1]], hw[[2]], S1, S2, Length[sorted], predicted]; Abort[]];
          MapIndexed[
            Function[{rec, idx}, Module[{n = First[idx] - 1, vv = rec["vec"], lc},
              lc = SelectFirst[vv, Abs[#] > 10^-13 Max[Abs[vv]] &];
              {hw[[1]], hw[[2]], n} -> <|"tau1c" -> rec["t1c"], "tau2c" -> rec["t2c"], "Psi" -> vv/lc|>]],
            sorted]]],
      byHW];
   If[Length[table] != d, Message[TauEigensystem3::degnull, "total", "-", S1, S2, Length[table], d]; Abort[]];
   table
];
```
**Output:**
```
Out= Null
```
## 2026-07-17 14:11:48 — ▶️ RUN CELL 26
**In [26]:** `ClearAll[Psi, tau1, tau2];↵Psi::bnds = "No state at (H1,H2,n)=(`1`,`2`,`3`) for (S1,S2)=(`4`,`5`).";↵Psi[S1_, S2_][H1_, H2_, n_] := Lookup[TauEigensystem3[S1, S2], Key[{H1, H2, n}],↵    (Message[Psi::`
**Out:** (no output)

## 2026-07-17 14:11:48 — ▶️ RUN CELL 27 (0.02 s)
**In [27]:** `Module[{S1 = 2, S2 = 1, es, worst},↵  es = TauEigensystem3[S1, S2];↵  worst = Max@KeyValueMap[↵     Function[{k, rec},↵       Module[{H1 = k[[1]], H2 = k[[2]], n = k[[3]], v = rec["Psi"], uu = 17/5},↵`
**Out:** Out[85]= \left\{\text{nstates}\to 18,\text{worstResidual@u=17/5}\to 1.359292769213167893921303\times 10^{-27}\right\}


## 2026-07-17 14:12:11
**Input:**
```wolfram
DownValues[TauEigensystem3] // Length
```
**Output:**
```
Out= 2
```

## 2026-07-17 14:12:11
**Input:**
```wolfram
StringContainsQ[ToString[DownValues[TauEigensystem3][[1]], InputForm], "Join[Re"]
```
**Output:**
```
Out= False
```

## 2026-07-17 14:12:27
**Input:**
```wolfram
DownValues[TauEigensystem3][[2]] // InputForm // Short
```
**Output:**
```
Out= Short[InputForm[HoldPattern[TauEigensystem3[S1_, S2_]] :> (TauEigensystem3[S1, S2] = Module[{d, c1, c2, mats, nontriv, rSets, rset, Kmix, vals, vecs, rayleigh, recs, byHW, table}, d = dim[{S1, 0, 0}]*dim[{S2, 0, 0}]; {c1, c2} = tgCoeffs[S1, S2]; mats = Join[c1, c2, {H1op[S1, S2], H2op[S1, S2]}]; nontriv = Join[c1[[1 ;; 2]], c2[[1 ;; 4]], {H1op[S1, S2], H2op[S1, S2]}]; rSets = Table[N[Table[Prime[100*i + 7*j]^(1/3), {j, Length[nontriv]}], 30], {i, 1, 8}]; rset = SelectFirst[rSets, Module[{ee, vv}, {ee, vv} = Eigensystem[#1 . nontriv]; MatrixRank[vv] == d && Length[Union[Chop[ee, 10^(-12)]]] == d] & , Missing["NoGenericR"]]; If[MissingQ[rset], Message[TauEigensystem3::rankdef, d, S1, S2]; Abort[]]; Kmix = rset . nontriv; {vals, vecs} = Eigensystem[N[Kmix, 30]]; rayleigh[m_, v_] := v . (m . v)/v . v; recs = Table[Module[{v = vecs[[i]], quotients, res, H1v, H2v, t1c, t2c}, quotients = (rayleigh[#1, v] & ) /@ mats; res = Max[(Max[Abs[#1 . v - rayleigh[#1, v]*v]] & ) /@ mats]; If[res > 10^(-13), Message[TauEigensystem3::notgenuine, i, S1, S2, res]; Abort[]]; t1c = (rayleigh[#1, v] & ) /@ c1; t2c = (rayleigh[#1, v] & ) /@ c2; H1v = rayleigh[H1op[S1, S2], v]; H2v = rayleigh[H2op[S1, S2], v]; If[Abs[H1v - Round[Re[H1v]]] > 10^(-10), Message[TauEigensystem3::noninteger, H1v, i, S1, S2]; Abort[]]; If[Abs[H2v - Round[Re[H2v]]] > 10^(-10), Message[TauEigensystem3::noninteger, H2v, i, S1, S2]; Abort[]]; Association["H1" -> Round[Re[H1v]], "H2" -> Round[Re[H2v]], "t1c" -> t1c, "t2c" -> t2c, "vec" -> v]], {i, 1, d}]; byHW = GroupBy[recs, {#1["H1"], #1["H2"]} & ]; table = Association[Flatten[KeyValueMap[Function[{hw, group}, Module[{sorted, predicted, last}, predicted = wtMult[S1, S2][hw[[1]], hw[[2]]]; sorted = SortBy[group, Join[Re[#1["t1c"]], Im[#1["t1c"]], Re[#1["t2c"]], Im[#1["t2c"]]] & ]; If[Length[sorted] != predicted, Message[TauEigensystem3::degnull, hw[[1]], hw[[2]], S1, S2, Length[sorted], predicted]; Abort[]]; MapIndexed[Function[{rec, idx}, Module[{n = First[idx] - 1, vv = rec["vec"], lc}, lc = SelectFirst[vv, Abs[#1] > Max[Abs[vv]]/10^13 & ]; {hw[[1]], hw[[2]], n} -> Association["tau1c" -> rec["t1c"], "tau2c" -> rec["t2c"], "Psi" -> vv/lc]]], sorted]]], byHW]]]; If[Length[table] != d, Message[TauEigensystem3::degnull, "total", "-", S1, S2, Length[table], d]; Abort[]]; table])]]
```

## 2026-07-17 14:12:27
**Input:**
```wolfram
ToString[DownValues[TauEigensystem3][[2]], InputForm] // StringLength
```
**Output:**
```
Out= 2297
```

## 2026-07-17 14:12:37
**Input:**
```wolfram
Module[{S1 = 1, S2 = 1, es, worst},
  es = TauEigensystem3[S1, S2];
  worst = Max@KeyValueMap[
     Function[{k, rec},
       Module[{H1 = k[[1]], H2 = k[[2]], n = k[[3]], v = rec["Psi"], uu = 17/5},
         Max[
          Max[Abs[(tg[S1, S2][1, 1][uu] . v) - tau1[S1, S2][H1, H2, n][uu] v]],
          Max[Abs[(tg[S1, S2][2, 1][uu] . v) - tau2[S1, S2][H1, H2, n][uu] v]]]]],
     es];
  If[worst > 10^-11, Print["EIGEN RESIDUAL FAIL: ", worst]; Abort[]];
  {"nstates" -> Length[es], "worstResidual@u=17/5" -> worst, "expected" -> dim[{S1,0,0}] dim[{S2,0,0}]}
]
```
**Output:**
```
Out= {"nstates" -> 9, "worstResidual@u=17/5" -> 4.351882429728977936066513775178630089543`24.81639060791746*^-38, "expected" -> 9}
```

## 2026-07-17 14:12:40
**Input:**
```wolfram
Module[{S1 = 2, S2 = 2, es, worst},
  es = TauEigensystem3[S1, S2];
  worst = Max@KeyValueMap[
     Function[{k, rec},
       Module[{H1 = k[[1]], H2 = k[[2]], n = k[[3]], v = rec["Psi"], uu = 17/5},
         Max[
          Max[Abs[(tg[S1, S2][1, 1][uu] . v) - tau1[S1, S2][H1, H2, n][uu] v]],
          Max[Abs[(tg[S1, S2][2, 1][uu] . v) - tau2[S1, S2][H1, H2, n][uu] v]]]]],
     es];
  If[worst > 10^-11, Print["EIGEN RESIDUAL FAIL: ", worst]; Abort[]];
  {"nstates" -> Length[es], "worstResidual@u=17/5" -> worst, "expected" -> dim[{S1,0,0}] dim[{S2,0,0}]}
]
```
**Output:**
```
Out= {"nstates" -> 36, "worstResidual@u=17/5" -> 1.798311988555260034005821300982471455924115020616055`25.012439748753433*^-26, "expected" -> 36}
```
## 2026-07-17 14:23:58 — 📥 BULK INSERT 3 CELL(S) at positions 28–30
1. [code] (*3rd-order TQ residual functional (operator form used for verification below): Qf is a pure↵  funct
2. [code] (*large-u asymptotics: alphaCoeffs[S1,S2][H1,H2,n] = {α1,α2,α3}, the 1/u coefficients of the↵  three
3. [code] (*asymptotic master formula: for each companion-twist eigenvalue z[k] (k=1,2,3), the large-u↵  balan


## 2026-07-17 14:23:58
**Input:**
```wolfram
(*3rd-order TQ residual functional (operator form used for verification below): Qf is a pure
  function of u; residual should vanish on the true Baxter Q-function once M[k] is known.
  Same A1/A2/A3 structure as Bax[S1,S2] above but exposed as a callable per-state functional.*)
ClearAll[BaxOp];
BaxOp[S1_, S2_][H1_, H2_, n_][Qf_][u_] :=
  Qf[u + h]
   - tau1[S1, S2][H1, H2, n][u]/ν1[S1, S2][u] Qf[u]
   + tau2[S1, S2][H1, H2, n][u]/(ν1[S1, S2][u] ν1[S1, S2][u - h]) Qf[u - h]
   - qdetT[S1, S2][u]/(ν1[S1, S2][u] ν1[S1, S2][u - h] ν1[S1, S2][u - 2 h]) Qf[u - 2 h];
```
**Output:**
```
(no output)
```

## 2026-07-17 14:23:58
**Input:**
```wolfram
(*large-u asymptotics: alphaCoeffs[S1,S2][H1,H2,n] = {α1,α2,α3}, the 1/u coefficients of the
  three rational combinations A1=τ1/ν1, A2=τ2/(ν1[u]ν1[u-h]), A3=qdetT/(ν1[u]ν1[u-h]ν1[u-2h])
  that appear in BaxOp. These feed the asymptotic master formula for the polynomial degree Mdeg.*)
ClearAll[alphaCoeffs];
alphaCoeffs[S1_, S2_][H1_, H2_, n_] := alphaCoeffs[S1, S2][H1, H2, n] = Module[{u, A1, A2, A3},
   A1 = tau1[S1, S2][H1, H2, n][u]/ν1[S1, S2][u];
   A2 = tau2[S1, S2][H1, H2, n][u]/(ν1[S1, S2][u] ν1[S1, S2][u - h]);
   A3 = qdetT[S1, S2][u]/(ν1[S1, S2][u] ν1[S1, S2][u - h] ν1[S1, S2][u - 2 h]);
   SeriesCoefficient[#, {u, Infinity, -1}] & /@ {A1, A2, A3}
];
```
**Output:**
```
(no output)
```

## 2026-07-17 14:23:58
**Input:**
```wolfram
(*asymptotic master formula: for each companion-twist eigenvalue z[k] (k=1,2,3), the large-u
  balance of the TQ equation fixes the polynomial degree M[k] of the Q-function branch built on
  z[k]. Denominator h(z[k]-χ2/z[k]+2χ3/z[k]^2) is the z[k]-characteristic factor of z^3-χ1 z^2+χ2 z-χ3
  divided by z^2 (companion twist has eigenvalues z[1],z[2],z[3] with χ1=Σz,χ2=Σzz,χ3=z1z2z3).*)
ClearAll[Mdeg];
Mdeg[S1_, S2_][H1_, H2_, n_] := Mdeg[S1, S2][H1, H2, n] = Module[{a, M},
   a = alphaCoeffs[S1, S2][H1, H2, n];
   M = Table[
      (a[[1]] - a[[2]]/z[k] + a[[3]]/z[k]^2)/(h (z[k] - χ2/z[k] + 2 χ3/z[k]^2)),
      {k, 1, 3}];
   If[Max[Abs[M - Round[Re[M]]]] > 10^-9, Print["Mdeg NON-INTEGER: ", M]; Abort[]];
   Round[Re[M]]
];
```
**Output:**
```
(no output)
```
## 2026-07-17 14:24:09 — 📥 BULK INSERT 1 CELL(S) at positions 31–31
1. [code] (*degree tabulation for (S1,S2)=(2,1): tabulate Mdeg over every (H1,H2,n) sector and check↵  M1+M2+M


## 2026-07-17 14:24:09
**Input:**
```wolfram
(*degree tabulation for (S1,S2)=(2,1): tabulate Mdeg over every (H1,H2,n) sector and check
  M1+M2+M3 against the sum rule (placeholder comment -- filled in once the pattern is inspected).*)
Module[{S1 = 2, S2 = 1, es, tab},
  es = TauEigensystem3[S1, S2];
  tab = KeyValueMap[
     Function[{k, rec},
       {"H1H2n" -> k, "M" -> Mdeg[S1, S2][k[[1]], k[[2]], k[[3]]],
        "sumM" -> Total[Mdeg[S1, S2][k[[1]], k[[2]], k[[3]]]]}],
     KeySort[es]];
  Column[tab]
]
```
**Output:**
```
Out[92]= \begin{matrix}\left\{\text{H1H2n}\to \left\{-3,3,0\right\},\text{M}\to \left\{0,0,0\right\},\text{sumM}\to 0\right\}\\\left\{\text{H1H2n}\to \left\{-2,1,0\right\},\text{M}\to \left\{0,0,0\right\},\text{sumM}\to 0\right\}\\\left\{\text{H1H2n}\to \left\{-2,1,1\right\},\text{M}\to \left\{0,0,0\right\},\text{sumM}\to 0\right\}\\\left\{\text{H1H2n}\to \left\{-1,-1,0\right\},\text{M}\to \left\{… [1246 chars total]
```

## 2026-07-17 14:24:15
**Input:**
```wolfram
Module[{S1 = 2, S2 = 1, es}, es = TauEigensystem3[S1, S2]; KeyValueMap[{#1, Mdeg[S1, S2][#1[[1]], #1[[2]], #1[[3]]], Total[Mdeg[S1, S2][#1[[1]], #1[[2]], #1[[3]]]]} &, KeySort[es]]]
```
**Output:**
```
Out= {{{-3, 3, 0}, {0, 0, 0}, 0}, {{-2, 1, 0}, {0, 0, 0}, 0}, {{-2, 1, 1}, {0, 0, 0}, 0}, {{-1, -1, 0}, {0, 0, 0}, 0}, {{-1, -1, 1}, {0, 0, 0}, 0}, {{-1, 2, 0}, {0, 0, 0}, 0}, {{-1, 2, 1}, {0, 0, 0}, 0}, {{0, -3, 0}, {0, 0, 0}, 0}, {{0, 0, 0}, {0, 0, 0}, 0}, {{0, 0, 1}, {0, 0, 0}, 0}, {{0, 0, 2}, {0, 0, 0}, 0}, {{1, -2, 0}, {0, 0, 0}, 0}, {{1, -2, 1}, {0, 0, 0}, 0}, {{1, 1, 0}, {0, 0, 0}, 0}, {{1, 1, 1}, {0, 0, 0}, 0}, {{2, -1, 0}, {0, 0, 0}, 0}, {{2, -1, 1}, {0, 0, 0}, 0}, {{3, 0, 0}, {0, 0, 0}, 0}}
```

## 2026-07-17 14:24:27
**Input:**
```wolfram
alphaCoeffs[2, 1][0, 0, 0]
tau1[2,1][0,0,0][u]
tau2[2,1][0,0,0][u]
qdetT[2,1][u]
ν1[2,1][u]
```
**Output:**
```
Out[1]= {0, 0, 0}
Out[2]= (-0.15512013376208799846104511365444896446`26.500385257974646 - 0.18794996225177864310346104731711533909`26.493151198980833*I) - (2.07045329707367132756194154317276348587`26.747410459548707 + 3.8965230231511761152791772294450198914`26.79227529613108*I)*u + (1.40256513672732573802583136795574300656`26.747410459548707 + 2.63958011245724833615686199414017476514`26.79227529613108*I)*u^2
Out[3]= (-0.33491884715956547424630537673157289559`26.737309510338484 + 0.5069654674892599678957945274667658799`26.773465994626267*I) + (3.70586070441822421030231734685534343324`26.74380186484431 - 5.5975902608202412907674028653642821575`26.77807385239469*I)*u - (10.67970021314331388038020333623079158101`26.749019890651734 - 16.04357373298151878271475682849841377097`26.781712999420648*I)*u^2 + (8.28693783494140125277608146504806402659`26.757498397181568 - 12.26607786878690723608356785193005760694`26.78761486212262*I)*u^3 - (1.67332398590162909911824721890393600537`26.757498397181568 - 2.47680418504351011497841273933203086294`26.787614862122627*I)*u^4
Out[4]= (-0.9944145978415286577806456150645332908393735038958268340201`29.44245325110677 - 0.1055443394961133295075060113382859913751619075649164082601`28.46832070644821*I)*(-13/3 + u)*(-22/7 + u)*(-4/3 + u)*(-8/7 + u)*(-1/3 + u)*(-1/7 + u)
Out[5]= (-7/3 + u)*(-8/7 + u)
```

## 2026-07-17 14:24:47
**Input:**
```wolfram
Module[{S1=2,S2=1,H1=1,H2=1,n=0,u,A1,A2,A3,a},
  A1 = tau1[S1,S2][H1,H2,n][u]/ν1[S1,S2][u];
  A2 = tau2[S1,S2][H1,H2,n][u]/(ν1[S1,S2][u] ν1[S1,S2][u-h]);
  A3 = qdetT[S1,S2][u]/(ν1[S1,S2][u] ν1[S1,S2][u-h] ν1[S1,S2][u-2h]);
  a = SeriesCoefficient[#,{u,Infinity,-1}]&/@{A1,A2,A3};
  {a, Series[A1,{u,Infinity,2}], Series[A2,{u,Infinity,2}]}
]
```
**Output:**
```
Out= {{0, 0, 0}, SeriesData[u$7539, Infinity, {1.40256513672732573802583136795574300656`26.775788737879044 + 2.63958011245724833615686199414017476514`26.806178429897088*I, 2.97794827889702431138081743040906021718`26.41736829742432 + 5.19233051305591648231226962466262419068`26.402972864910733*I, 6.48071775978741516877975030881149336448`26.282567120425124 + 10.80634225323392026440255505740526579581`26.247996502131244*I}, 0, 3, 1], SeriesData[u$7539, Infinity, {-1.67332398590162909911824721890393600537`26.76969982785933 + 2.47680418504351011497841273933203086293`26.807113153392038*I, -6.53060758081151578762768733242844115749`26.215423888654712 + 10.01180143541003725099719518406068749956`26.271705252089504*I, -21.42813752278974898322441925002127622476`25.86365383607618 + 33.86724124761183916771377041769875355945`25.93535501778722*I}, 0, 3, 1]}
```

## 2026-07-17 14:25:02
**Input:**
```wolfram
Module[{S1=2,S2=1,H1=1,H2=1,n=0,u,A1,s},
  A1 = tau1[S1,S2][H1,H2,n][u]/ν1[S1,S2][u];
  s = Series[A1,{u,Infinity,2}];
  {SeriesCoefficient[A1,{u,Infinity,-1}], SeriesCoefficient[s,{u,Infinity,-1}], SeriesCoefficient[A1,{u,Infinity,0}], SeriesCoefficient[A1,{u,Infinity,-2}]}
]
```
**Output:**
```
Out= {0, 0, 1.40256513672732573802583136795574300656`26.775788737879044 + 2.63958011245724833615686199414017476514`26.806178429897088*I, 0}
```

## 2026-07-17 14:25:18
**Input:**
```wolfram
f = 1 + 2/u + 3/u^2;
{SeriesCoefficient[f,{u,Infinity,0}], SeriesCoefficient[f,{u,Infinity,-1}], SeriesCoefficient[f,{u,Infinity,-2}], SeriesCoefficient[f,{u,Infinity,1}]}
```
**Output:**
```
Out= {1, 0, 0, 2}
```

## 2026-07-17 14:25:41
**Input:**
```wolfram
(* Confirm SeriesCoefficient[f,{x,Infinity,n}] convention: for x0=Infinity, WMA returns coefficient of x^(-n), NOT x^n.
   So to get the coefficient of the 1/u term (i.e. x^(-1)) we need n=+1, not n=-1 as literally written in the brief. *)
g = 5 + 7/u - 11/u^2 + 13/u^3;
{"coeff of u^0 (n=0)" -> SeriesCoefficient[g,{u,Infinity,0}],
 "coeff of u^-1, requested n=-1 (BRIEF LITERAL)" -> SeriesCoefficient[g,{u,Infinity,-1}],
 "coeff of u^-1, requested n=+1 (CORRECT)" -> SeriesCoefficient[g,{u,Infinity,1}],
 "coeff of u^-2, requested n=+2" -> SeriesCoefficient[g,{u,Infinity,2}]}
```
**Output:**
```
Out= {"coeff of u^0 (n=0)" -> 5, "coeff of u^-1, requested n=-1 (BRIEF LITERAL)" -> 0, "coeff of u^-1, requested n=+1 (CORRECT)" -> 7, "coeff of u^-2, requested n=+2" -> -11}
```
## 2026-07-17 14:26:05 — ✏️ EDIT CELL 29
(*large-u asymptotics: alphaCoeffs[S1,S2][H1,H2,n] = {α1,α2,α3}, the 1/u coefficients of the
  three rational combinations A1=τ1/ν1, A2=τ2/(ν1[u]ν1[u-h]), A3=qdetT/(ν1[u]ν1[u-h]ν1[u-2h])
  that appear…


## 2026-07-17 14:26:05
**Input:**
```wolfram
(*large-u asymptotics: alphaCoeffs[S1,S2][H1,H2,n] = {α1,α2,α3}, the 1/u coefficients of the
  three rational combinations A1=τ1/ν1, A2=τ2/(ν1[u]ν1[u-h]), A3=qdetT/(ν1[u]ν1[u-h]ν1[u-2h])
  that appear in BaxOp. These feed the asymptotic master formula for the polynomial degree Mdeg.
  NB (deviation from task-5-brief.md verbatim code, flagged not silently patched): for x0=Infinity,
  SeriesCoefficient[f,{x,Infinity,n}] returns the coefficient of x^(-n), so the 1/u term needs
  n=+1, NOT n=-1. Brief's literal `{u,Infinity,-1}` asks for the (nonexistent, since A1,A2,A3->
  finite constants) u^{+1} term and silently returns 0 for every state -- verified directly:
  SeriesCoefficient[5+7/u-11/u^2,{u,Infinity,-1}]=0 but {u,Infinity,1}=7 (the true 1/u coeff).
  Using n=-1 made Mdeg trivially {0,0,0} everywhere (a false pass, not a NON-INTEGER abort).*)
ClearAll[alphaCoeffs];
alphaCoeffs[S1_, S2_][H1_, H2_, n_] := alphaCoeffs[S1, S2][H1, H2, n] = Module[{u, A1, A2, A3},
   A1 = tau1[S1, S2][H1, H2, n][u]/ν1[S1, S2][u];
   A2 = tau2[S1, S2][H1, H2, n][u]/(ν1[S1, S2][u] ν1[S1, S2][u - h]);
   A3 = qdetT[S1, S2][u]/(ν1[S1, S2][u] ν1[S1, S2][u - h] ν1[S1, S2][u - 2 h]);
   SeriesCoefficient[#, {u, Infinity, 1}] & /@ {A1, A2, A3}
];
```
**Output:**
```
Out= Null
```
## 2026-07-17 14:26:16 — ▶️ RUN CELL 30
**In [30]:** `(*asymptotic master formula: for each companion-twist eigenvalue z[k] (k=1,2,3), the large-u↵  balance of the TQ equation fixes the polynomial degree M[k] of the Q-function branch built on↵  z[k]. Den`
**Out:** (no output)

## 2026-07-17 14:26:17 — ▶️ RUN CELL 31 (0.04 s)
**In [31]:** `(*degree tabulation for (S1,S2)=(2,1): tabulate Mdeg over every (H1,H2,n) sector and check↵  M1+M2+M3 against the sum rule (placeholder comment -- filled in once the pattern is inspected).*)↵Module[{S`
**Out:** Out[95]= \begin{matrix}\left\{\text{H1H2n}\to \left\{-3,3,0\right\},\text{M}\to \left\{3,0,3\right\},\text{sumM}\to 6\right\}\\\left\{\text{H1H2n}\to \left\{-2,1,0\right\},\text{M}\to \left\{3,1,2\right\},\text{sumM}\to 6\right\}\\\left\{\text{H1H2n}\to \left\{-2,1,1\right\},\text{M}\to \left\{3,1,2


## 2026-07-17 14:26:23
**Input:**
```wolfram
Module[{S1=2,S2=1,es},es=TauEigensystem3[S1,S2];KeyValueMap[{#1,Mdeg[S1,S2][#1[[1]],#1[[2]],#1[[3]]],Total[Mdeg[S1,S2][#1[[1]],#1[[2]],#1[[3]]]]}&,KeySort[es]]]
```
**Output:**
```
Out= {{{-3, 3, 0}, {3, 0, 3}, 6}, {{-2, 1, 0}, {3, 1, 2}, 6}, {{-2, 1, 1}, {3, 1, 2}, 6}, {{-1, -1, 0}, {3, 2, 1}, 6}, {{-1, -1, 1}, {3, 2, 1}, 6}, {{-1, 2, 0}, {2, 1, 3}, 6}, {{-1, 2, 1}, {2, 1, 3}, 6}, {{0, -3, 0}, {3, 3, 0}, 6}, {{0, 0, 0}, {2, 2, 2}, 6}, {{0, 0, 1}, {2, 2, 2}, 6}, {{0, 0, 2}, {2, 2, 2}, 6}, {{1, -2, 0}, {2, 3, 1}, 6}, {{1, -2, 1}, {2, 3, 1}, 6}, {{1, 1, 0}, {1, 2, 3}, 6}, {{1, 1, 1}, {1, 2, 3}, 6}, {{2, -1, 0}, {1, 3, 2}, 6}, {{2, -1, 1}, {1, 3, 2}, 6}, {{3, 0, 0}, {0, 3, 3}, 6}}
```

## 2026-07-17 14:26:24
**Input:**
```wolfram
Module[{S1=1,S2=1,es},es=TauEigensystem3[S1,S2];KeyValueMap[{#1,Mdeg[S1,S2][#1[[1]],#1[[2]],#1[[3]]],Total[Mdeg[S1,S2][#1[[1]],#1[[2]],#1[[3]]]]}&,KeySort[es]]]
```
**Output:**
```
Out= {{{-2, 2, 0}, {2, 0, 2}, 4}, {{-1, 0, 0}, {2, 1, 1}, 4}, {{-1, 0, 1}, {2, 1, 1}, 4}, {{0, -2, 0}, {2, 2, 0}, 4}, {{0, 1, 0}, {1, 1, 2}, 4}, {{0, 1, 1}, {1, 1, 2}, 4}, {{1, -1, 0}, {1, 2, 1}, 4}, {{1, -1, 1}, {1, 2, 1}, 4}, {{2, 0, 0}, {0, 2, 2}, 4}}
```

## 2026-07-17 14:26:25
**Input:**
```wolfram
Module[{S1=2,S2=2,es},es=TauEigensystem3[S1,S2];KeyValueMap[{#1,Mdeg[S1,S2][#1[[1]],#1[[2]],#1[[3]]],Total[Mdeg[S1,S2][#1[[1]],#1[[2]],#1[[3]]]]}&,KeySort[es]]]
```
**Output:**
```
Out= {{{-4, 4, 0}, {4, 0, 4}, 8}, {{-3, 2, 0}, {4, 1, 3}, 8}, {{-3, 2, 1}, {4, 1, 3}, 8}, {{-2, 0, 0}, {4, 2, 2}, 8}, {{-2, 0, 1}, {4, 2, 2}, 8}, {{-2, 0, 2}, {4, 2, 2}, 8}, {{-2, 3, 0}, {3, 1, 4}, 8}, {{-2, 3, 1}, {3, 1, 4}, 8}, {{-1, -2, 0}, {4, 3, 1}, 8}, {{-1, -2, 1}, {4, 3, 1}, 8}, {{-1, 1, 0}, {3, 2, 3}, 8}, {{-1, 1, 1}, {3, 2, 3}, 8}, {{-1, 1, 2}, {3, 2, 3}, 8}, {{-1, 1, 3}, {3, 2, 3}, 8}, {{0, -4, 0}, {4, 4, 0}, 8}, {{0, -1, 0}, {3, 3, 2}, 8}, {{0, -1, 1}, {3, 3, 2}, 8}, {{0, -1, 2}, {3, 3, 2}, 8}, {{0, -1, 3}, {3, 3, 2}, 8}, {{0, 2, 0}, {2, 2, 4}, 8}, {{0, 2, 1}, {2, 2, 4}, 8}, {{0, 2, 2}, {2, 2, 4}, 8}, {{1, -3, 0}, {3, 4, 1}, 8}, {{1, -3, 1}, {3, 4, 1}, 8}, {{1, 0, 0}, {2, 3, 3}, 8}, {{1, 0, 1}, {2, 3, 3}, 8}, {{1, 0, 2}, {2, 3, 3}, 8}, {{1, 0, 3}, {2, 3, 3}, 8}, {{2, -2, 0}, {2, 4, 2}, 8}, {{2, -2, 1}, {2, 4, 2}, 8}, {{2, -2, 2}, {2, 4, 2}, 8}, {{2, 1, 0}, {1, 3, 4}, 8}, {{2, 1, 1}, {1, 3, 4}, 8}, {{3, -1, 0}, {1, 4, 3}, 8}, {{3, -1, 1}, {1, 4, 3}, 8}, {{4, 0, 0}, {0, 4, 4}, 8}}
```

## 2026-07-17 14:26:27
**Input:**
```wolfram
Module[{S1=1,S2=3,es},es=TauEigensystem3[S1,S2];KeyValueMap[{#1,Mdeg[S1,S2][#1[[1]],#1[[2]],#1[[3]]],Total[Mdeg[S1,S2][#1[[1]],#1[[2]],#1[[3]]]]}&,KeySort[es]]]
```
**Output:**
```
Out= {{{-4, 4, 0}, {4, 0, 4}, 8}, {{-3, 2, 0}, {4, 1, 3}, 8}, {{-3, 2, 1}, {4, 1, 3}, 8}, {{-2, 0, 0}, {4, 2, 2}, 8}, {{-2, 0, 1}, {4, 2, 2}, 8}, {{-2, 3, 0}, {3, 1, 4}, 8}, {{-2, 3, 1}, {3, 1, 4}, 8}, {{-1, -2, 0}, {4, 3, 1}, 8}, {{-1, -2, 1}, {4, 3, 1}, 8}, {{-1, 1, 0}, {3, 2, 3}, 8}, {{-1, 1, 1}, {3, 2, 3}, 8}, {{-1, 1, 2}, {3, 2, 3}, 8}, {{0, -4, 0}, {4, 4, 0}, 8}, {{0, -1, 0}, {3, 3, 2}, 8}, {{0, -1, 1}, {3, 3, 2}, 8}, {{0, -1, 2}, {3, 3, 2}, 8}, {{0, 2, 0}, {2, 2, 4}, 8}, {{0, 2, 1}, {2, 2, 4}, 8}, {{1, -3, 0}, {3, 4, 1}, 8}, {{1, -3, 1}, {3, 4, 1}, 8}, {{1, 0, 0}, {2, 3, 3}, 8}, {{1, 0, 1}, {2, 3, 3}, 8}, {{1, 0, 2}, {2, 3, 3}, 8}, {{2, -2, 0}, {2, 4, 2}, 8}, {{2, -2, 1}, {2, 4, 2}, 8}, {{2, 1, 0}, {1, 3, 4}, 8}, {{2, 1, 1}, {1, 3, 4}, 8}, {{3, -1, 0}, {1, 4, 3}, 8}, {{3, -1, 1}, {1, 4, 3}, 8}, {{4, 0, 0}, {0, 4, 4}, 8}}
```

## 2026-07-17 14:27:05
**Input:**
```wolfram
(* verify M_k = N - n_k, i.e. M1=(2N-2H1-H2)/3, M2=(2N+H1-H2)/3, M3=(2N+H1+2H2)/3, across all four reps *)
Module[{reps = {{2,1},{1,1},{2,2},{1,3}}, worst},
 worst = Max@Table[
   Module[{S1=rep[[1]],S2=rep[[2]],N,es,maxdiff},
     N = S1+S2; es = TauEigensystem3[S1,S2];
     maxdiff = Max@KeyValueMap[Function[{k,rec},
       Module[{H1=k[[1]],H2=k[[2]],n=k[[3]],M,Mpred},
         M = Mdeg[S1,S2][H1,H2,n];
         Mpred = {(2N-2H1-H2)/3,(2N+H1-H2)/3,(2N+H1+2H2)/3};
         Max[Abs[M-Mpred]]]], es];
     maxdiff],
   {rep, reps}];
 {"worst |M - Mpredicted| over all four reps" -> worst}
]
```
**Output:**
```
Out= {"worst |M - Mpredicted| over all four reps" -> 0}
```
## 2026-07-17 14:27:20 — ✏️ EDIT CELL 31
(*degree tabulation for (S1,S2)=(2,1): tabulate Mdeg over every (H1,H2,n) sector and check
  M1+M2+M3=2N (N=S1+S2) against the sum rule.
  DEGREE RELATION FOUND (verified exactly, zero residual, acros…


## 2026-07-17 14:27:20
**Input:**
```wolfram
(*degree tabulation for (S1,S2)=(2,1): tabulate Mdeg over every (H1,H2,n) sector and check
  M1+M2+M3=2N (N=S1+S2) against the sum rule.
  DEGREE RELATION FOUND (verified exactly, zero residual, across (2,1),(1,1),(2,2),(1,3)):
    M_k = N - n_k,  N=S1+S2,  (n1,n2,n3) the same weight components used in wtMult:
      n3=(N-H1-2H2)/3,  n2=n3+H2,  n1=n3+H1+H2
    equivalently, directly in (H1,H2):
      M1=(2N-2H1-H2)/3,  M2=(2N+H1-H2)/3,  M3=(2N+H1+2H2)/3
    Ordering: M[k] pairs directly with z[k] (companion-twist eigenvalue k) and with n_k in the
    SAME index k used by wtMult's (n1,n2,n3) -- no permutation needed. Mdeg depends only on
    (H1,H2), NOT on the branch index n (all n at fixed (H1,H2) give identical M), consistent
    with M_k being a weight-only (not multiplicity-resolving) quantity.*)
Module[{S1 = 2, S2 = 1, es, tab},
  es = TauEigensystem3[S1, S2];
  tab = KeyValueMap[
     Function[{k, rec},
       {"H1H2n" -> k, "M" -> Mdeg[S1, S2][k[[1]], k[[2]], k[[3]]],
        "sumM" -> Total[Mdeg[S1, S2][k[[1]], k[[2]], k[[3]]]]}],
     KeySort[es]];
  Column[tab]
]
```
**Output:**
```
Out= Column[{{"H1H2n" -> {-3, 3, 0}, "M" -> {3, 0, 3}, "sumM" -> 6}, {"H1H2n" -> {-2, 1, 0}, "M" -> {3, 1, 2}, "sumM" -> 6}, {"H1H2n" -> {-2, 1, 1}, "M" -> {3, 1, 2}, "sumM" -> 6}, {"H1H2n" -> {-1, -1, 0}, "M" -> {3, 2, 1}, "sumM" -> 6}, {"H1H2n" -> {-1, -1, 1}, "M" -> {3, 2, 1}, "sumM" -> 6}, {"H1H2n" -> {-1, 2, 0}, "M" -> {2, 1, 3}, "sumM" -> 6}, {"H1H2n" -> {-1, 2, 1}, "M" -> {2, 1, 3}, "sumM" -> 6}, {"H1H2n" -> {0, -3, 0}, "M" -> {3, 3, 0}, "sumM" -> 6}, {"H1H2n" -> {0, 0, 0}, "M" -> {2, 2, 2}, "sumM" -> 6}, {"H1H2n" -> {0, 0, 1}, "M" -> {2, 2, 2}, "sumM" -> 6}, {"H1H2n" -> {0, 0, 2}, "M" -> {2, 2, 2}, "sumM" -> 6}, {"H1H2n" -> {1, -2, 0}, "M" -> {2, 3, 1}, "sumM" -> 6}, {"H1H2n" -> {1, -2, 1}, "M" -> {2, 3, 1}, "sumM" -> 6}, {"H1H2n" -> {1, 1, 0}, "M" -> {1, 2, 3}, "sumM" -> 6}, {"H1H2n" -> {1, 1, 1}, "M" -> {1, 2, 3}, "sumM" -> 6}, {"H1H2n" -> {2, -1, 0}, "M" -> {1, 3, 2}, "sumM" -> 6}, {"H1H2n" -> {2, -1, 1}, "M" -> {1, 3, 2}, "sumM" -> 6}, {"H1H2n" -> {3, 0, 0}, "M" -> {0, 3, 3}, "sumM" -> 6}}]
```
## 2026-07-17 14:27:42 — 📥 BULK INSERT 4 CELL(S) at positions 32–35
1. [code] (*degree tabulation for (S1,S2)=(1,1): same M_k=N-n_k relation holds, N=S1+S2=2.*)↵Module[{S1 = 1, S
2. [code] (*degree tabulation for (S1,S2)=(2,2): same M_k=N-n_k relation holds, N=S1+S2=4.*)↵Module[{S1 = 2, S
3. [code] (*degree tabulation for (S1,S2)=(1,3): same M_k=N-n_k relation holds, N=S1+S2=4.*)↵Module[{S1 = 1, S
4. [code] (*hard-assert: closed-form M_k=N-n_k (equivalently M1=(2N-2H1-H2)/3, M2=(2N+H1-H2)/3,↵  M3=(2N+H1+2H


## 2026-07-17 14:27:42
**Input:**
```wolfram
(*degree tabulation for (S1,S2)=(1,1): same M_k=N-n_k relation holds, N=S1+S2=2.*)
Module[{S1 = 1, S2 = 1, es, tab},
  es = TauEigensystem3[S1, S2];
  tab = KeyValueMap[
     Function[{k, rec},
       {"H1H2n" -> k, "M" -> Mdeg[S1, S2][k[[1]], k[[2]], k[[3]]],
        "sumM" -> Total[Mdeg[S1, S2][k[[1]], k[[2]], k[[3]]]]}],
     KeySort[es]];
  Column[tab]
]
```
**Output:**
```
Out[96]= \begin{matrix}\left\{\text{H1H2n}\to \left\{-2,2,0\right\},\text{M}\to \left\{2,0,2\right\},\text{sumM}\to 4\right\}\\\left\{\text{H1H2n}\to \left\{-1,0,0\right\},\text{M}\to \left\{2,1,1\right\},\text{sumM}\to 4\right\}\\\left\{\text{H1H2n}\to \left\{-1,0,1\right\},\text{M}\to \left\{2,1,1\right\},\text{sumM}\to 4\right\}\\\left\{\text{H1H2n}\to \left\{0,-2,0\right\},\text{M}\to \left\{2… [966 chars total]
```

## 2026-07-17 14:27:42
**Input:**
```wolfram
(*degree tabulation for (S1,S2)=(2,2): same M_k=N-n_k relation holds, N=S1+S2=4.*)
Module[{S1 = 2, S2 = 2, es, tab},
  es = TauEigensystem3[S1, S2];
  tab = KeyValueMap[
     Function[{k, rec},
       {"H1H2n" -> k, "M" -> Mdeg[S1, S2][k[[1]], k[[2]], k[[3]]],
        "sumM" -> Total[Mdeg[S1, S2][k[[1]], k[[2]], k[[3]]]]}],
     KeySort[es]];
  Column[tab]
]
```
**Output:**
```
Out[97]= \begin{matrix}\left\{\text{H1H2n}\to \left\{-4,4,0\right\},\text{M}\to \left\{4,0,4\right\},\text{sumM}\to 8\right\}\\\left\{\text{H1H2n}\to \left\{-3,2,0\right\},\text{M}\to \left\{4,1,3\right\},\text{sumM}\to 8\right\}\\\left\{\text{H1H2n}\to \left\{-3,2,1\right\},\text{M}\to \left\{4,1,3\right\},\text{sumM}\to 8\right\}\\\left\{\text{H1H2n}\to \left\{-2,0,0\right\},\text{M}\to \left\{4… [1246 chars total]
```

## 2026-07-17 14:27:43
**Input:**
```wolfram
(*degree tabulation for (S1,S2)=(1,3): same M_k=N-n_k relation holds, N=S1+S2=4.*)
Module[{S1 = 1, S2 = 3, es, tab},
  es = TauEigensystem3[S1, S2];
  tab = KeyValueMap[
     Function[{k, rec},
       {"H1H2n" -> k, "M" -> Mdeg[S1, S2][k[[1]], k[[2]], k[[3]]],
        "sumM" -> Total[Mdeg[S1, S2][k[[1]], k[[2]], k[[3]]]]}],
     KeySort[es]];
  Column[tab]
]
```
**Output:**
```
Out[98]= \begin{matrix}\left\{\text{H1H2n}\to \left\{-4,4,0\right\},\text{M}\to \left\{4,0,4\right\},\text{sumM}\to 8\right\}\\\left\{\text{H1H2n}\to \left\{-3,2,0\right\},\text{M}\to \left\{4,1,3\right\},\text{sumM}\to 8\right\}\\\left\{\text{H1H2n}\to \left\{-3,2,1\right\},\text{M}\to \left\{4,1,3\right\},\text{sumM}\to 8\right\}\\\left\{\text{H1H2n}\to \left\{-2,0,0\right\},\text{M}\to \left\{4… [1246 chars total]
```

## 2026-07-17 14:27:43
**Input:**
```wolfram
(*hard-assert: closed-form M_k=N-n_k (equivalently M1=(2N-2H1-H2)/3, M2=(2N+H1-H2)/3,
  M3=(2N+H1+2H2)/3, N=S1+S2) matches Mdeg exactly (zero residual) across all four test reps.*)
Module[{reps = {{2, 1}, {1, 1}, {2, 2}, {1, 3}}, worst},
  worst = Max@Table[
     Module[{S1 = rep[[1]], S2 = rep[[2]], Nt, es, maxdiff},
       Nt = S1 + S2; es = TauEigensystem3[S1, S2];
       maxdiff = Max@KeyValueMap[
          Function[{k, rec},
            Module[{H1 = k[[1]], H2 = k[[2]], n = k[[3]], M, Mpred},
              M = Mdeg[S1, S2][H1, H2, n];
              Mpred = {(2 Nt - 2 H1 - H2)/3, (2 Nt + H1 - H2)/3, (2 Nt + H1 + 2 H2)/3};
              Max[Abs[M - Mpred]]]],
          es];
       maxdiff],
     {rep, reps}];
  If[worst > 0, Print["DEGREE-FORMULA MISMATCH: ", worst]; Abort[]];
  {"worst |Mdeg - (N-n_k formula)| over (2,1),(1,1),(2,2),(1,3)" -> worst}
]
```
**Output:**
```
Out[99]= \left\{\text{worst |Mdeg - (N-n\_k formula)| over (2,1),(1,1),(2,2),(1,3)}\to 0\right\}
```

## 2026-07-17 14:27:49
**Input:**
```wolfram
{alphaCoeffs[2,1][0,0,0], alphaCoeffs[2,1][1,1,0], alphaCoeffs[1,1][-1,0,0]}
```
**Output:**
```
Out= {{2.80513027345465147605166273591148601313`26.35362641064783 + 5.27916022491449667231372398828034953027`26.398491247230204*I, -6.69329594360651639647298887561574402151`26.21639483179503 + 9.90721674017404045991365095732812345174`26.24651129673609*I, -5.96648758704917194668387369038719974504`29.44245325110677 - 0.63326603697667997704503606802971594825`28.46832070644821*I}, {2.97794827889702431138081743040906021718`26.41736829742432 + 5.19233051305591648231226962466262419068`26.402972864910733*I, -6.53060758081151578762768733242844115749`26.215423888654712 + 10.01180143541003725099719518406068749956`26.271705252089504*I, -5.96648758704917194668387369038719974504`29.44245325110677 - 0.63326603697667997704503606802971594825`28.46832070644821*I}, {1.76300501187455227529250901855129309784`26.559813714906316 + 3.5723625610506230526416979762110426946`26.667709131123708*I, -4.56309537678395219977919627093425715828`26.421878061669982 + 6.54088246018938882651148478095110807376`26.45797232495877*I, -3.9776583913661146311225824602581331633574940155833073360804`29.44245325110677 - 0.4221773579844533180300240453531439655006476302596656330406`28.46832070644821*I}}
```

## 2026-07-17 14:35:40
**Input:**
```wolfram
KeyExistsQ[TauEigensystem3[2,1], {1,1,0}]
Keys[TauEigensystem3[2,1]][[1]]
```
**Output:**
```
Out[1]= True
Out[2]= {-1, -1, 0}
```
## 2026-07-17 14:36:00 — 📥 BULK INSERT 4 CELL(S) at positions 36–39
1. [code] (*polynomial from an ascending coeff list {c0,c1,...,cM}: qpoly[clist,uu]=Sum ck uu^k.↵  Dummy-symbo
2. [code] (*QSolve3[S1,S2][H1,H2,n] solves the 3rd-order Baxter (TQ) equation for the three monic↵  Q-polynomi
3. [code] (*Q accessor: Q_k(u)=z[k]^(u/h) q_k(u), q_k the monic coeff list from QSolve3.*)↵ClearAll[Qfun];↵Qfu
4. [code] (*single-state TQ-residual check: each Q_k (built from QSolve3) must be annihilated by BaxOp,↵  test


## 2026-07-17 14:36:01
**Input:**
```wolfram
(*polynomial from an ascending coeff list {c0,c1,...,cM}: qpoly[clist,uu]=Sum ck uu^k.
  Dummy-symbol trick (build in x, substitute uu at the end) avoids 0^0 issues when uu=0.*)
ClearAll[qpoly];
qpoly[clist_, uu_] := Module[{x}, (clist . x^Range[0, Length[clist] - 1]) /. x -> uu];
```
**Output:**
```
(no output)
```

## 2026-07-17 14:36:01
**Input:**
```wolfram
(*QSolve3[S1,S2][H1,H2,n] solves the 3rd-order Baxter (TQ) equation for the three monic
  Q-polynomials q_1,q_2,q_3 (one per companion-twist eigenvalue z[k]), degrees M_k from Mdeg.
  For fixed k, write Q_k(u)=z[k]^(u/h) q_k(u) and substitute into BaxOp[S1,S2][H1,H2,n][Q_k][u]=0,
  then multiply through by ν1[u]ν1[u-h]ν1[u-2h] to clear denominators -- this turns the TQ
  equation into a polynomial identity in u, linear in the unknown coefficients of q_k. Each
  monomial q=u^j contributes one column (coefficient list in u); NullSpace of the resulting
  matrix gives the coefficient vector (expected 1-dim nullspace), monic-normalized by dividing
  by the top coefficient. Chop the (high-precision) matrix before NullSpace -- same pitfall as
  su(2)'s QSolve: residual roundoff in the numeric tau1/tau2/qdetT/ν1 data can otherwise mask
  the true nullspace.*)
ClearAll[QSolve3];
QSolve3::degnull = "Q_`4` nullspace dim `5` (expected 1) at (H1,H2,n)=(`1`,`2`,`3`).";
QSolve3::lastzero = "Q_`4` top coeff ~0 at (H1,H2,n)=(`1`,`2`,`3`).";
QSolve3[S1_, S2_][H1_, H2_, n_] := QSolve3[S1, S2][H1, H2, n] = Module[
   {Ms, u, T1, T2, QD, N1, N1m, N1mm, sols},
   Ms = Mdeg[S1, S2][H1, H2, n];
   T1 = tau1[S1, S2][H1, H2, n][u]; T2 = tau2[S1, S2][H1, H2, n][u];
   QD = qdetT[S1, S2][u];
   N1 = ν1[S1, S2][u]; N1m = ν1[S1, S2][u - h]; N1mm = ν1[S1, S2][u - 2 h];
   sols = Table[
     Module[{M = Ms[[k]], zk = z[k], cols, mat, ns, vec, last, maxdeg},
       maxdeg = M + 6;
       (* column j (j=0..M): coefficient list of the polynomial identity for q=u^j *)
       cols = Table[
          PadRight[CoefficientList[Expand[
             zk (u + h)^j N1 N1m N1mm
              - T1 u^j N1m N1mm
              + T2/zk (u - h)^j N1mm
              - QD/zk^2 (u - 2 h)^j], u], maxdeg + 1],
          {j, 0, M}];
       mat = Transpose[cols];                 (* (maxdeg+1) x (M+1) *)
       ns = NullSpace[Chop[mat, 10^-12], Tolerance -> 10^-10];
       If[Length[ns] != 1, Message[QSolve3::degnull, H1, H2, n, k, Length[ns]]; Abort[]];
       vec = ns[[1]]; last = vec[[-1]];
       If[Abs[last] < 10^-13 Max[Abs[vec]], Message[QSolve3::lastzero, H1, H2, n, k]; Abort[]];
       vec/last],
     {k, 1, 3}];
   sols
];
```
**Output:**
```
(no output)
```

## 2026-07-17 14:36:01
**Input:**
```wolfram
(*Q accessor: Q_k(u)=z[k]^(u/h) q_k(u), q_k the monic coeff list from QSolve3.*)
ClearAll[Qfun];
Qfun[S1_, S2_][H1_, H2_, n_][k_][u_] := z[k]^(u/h) qpoly[QSolve3[S1, S2][H1, H2, n][[k]], u];
```
**Output:**
```
(no output)
```

## 2026-07-17 14:36:01
**Input:**
```wolfram
(*single-state TQ-residual check: each Q_k (built from QSolve3) must be annihilated by BaxOp,
  tested at several generic numeric u. State: (S1,S2)=(2,1), (H1,H2,n)=(1,1,0) (confirmed to
  exist as a key of TauEigensystem3[2,1]).*)
Module[{S1 = 2, S2 = 1, H1 = 1, H2 = 1, n = 0, res},
  res = Table[
     Max@Table[
        Abs@Chop[N[BaxOp[S1, S2][H1, H2, n][Function[uu, Qfun[S1, S2][H1, H2, n][k][uu]]][u0], 30]],
        {u0, {2/3, 9/5, 4, 23/7}}],
     {k, 1, 3}];
  If[Max[res] > 10^-10, Print["TQ RESIDUAL FAIL: ", res]; Abort[]];
  {"perQresidual" -> res}
]
```
**Output:**
```
Out[108]= \left\{\text{perQresidual}\to \left\{0,0,0\right\}\right\}
```

## 2026-07-17 14:36:10
**Input:**
```wolfram
Mdeg[2,1][1,1,0]
QSolve3[2,1][1,1,0]
Length /@ QSolve3[2,1][1,1,0]
(* raw (unchopped) residual magnitudes at the check's sample points, for the report *)
Table[Max@Table[Abs@N[BaxOp[2,1][1,1,0][Function[uu,Qfun[2,1][1,1,0][k][uu]]][u0],30],{u0,{2/3,9/5,4,23/7}}],{k,1,3}]
```
**Output:**
```
Out[1]= {1, 2, 3}
Out[2]= {{-0.87329087531727873212628924736267922064`21.93822674206486 - 0.02170986053280489736756284158544904972`21.93822674206486*I, 1.`21.93822674206486 + 0``23.241966951186356*I}, {1.2328790553408441735579573982047936162`21.882740047994726 + 5.2425786893665903094800670764527553777`21.882740047994726*I, -2.43623293420653079168323456216113030321`21.882740047994726 - 6.02038916124948659109338937230062644996`21.851112521058038*I, 0.99999999999999999999999999999999999996`21.882740047994726 + 0``22.23371183098432*I}, {201.33650312334816926186024785952498046408`21.932167092538627 - 53.60276779392818897919178394830966184851`21.932167092538627*I, -229.14084188932420322841758855062085225257`21.89105029669843 + 106.7671993628138167764748607844613960598`21.932167092538624*I, -3.30952380952380952380952380952381040863`21.15870889762202 - 45.28151923343296681469469675223827980061`21.932167092538624*I, 1.0000000000000000000000000000000002916`21.932167092538627 + 0``22.23561358659119*I}}
Out[3]= {2, 3, 4}
Out[4]= {0``18.87513251073512, 0``17.825131129807, 0``16.257496427896893}
```
## 2026-07-17 14:41:40 — 📥 BULK INSERT 1 CELL(S) at positions 40–40
1. [code] (*full sweep verification: over reps {(1,1),(2,1),(1,2),(2,2),(3,1),(1,3)} and every↵  {H1,H2,n} sta


## 2026-07-17 14:41:42
**Input:**
```wolfram
(*full sweep verification: over reps {(1,1),(2,1),(1,2),(2,2),(3,1),(1,3)} and every
  {H1,H2,n} state, check (a) per-Q_k TQ residual and (b) the 4x4 Casoratian reconstruction --
  the three Q_k span the TQ solution space, so A1,A2,A3 rebuilt from 3x3 minors over the four
  shift rows {u+h,u,u-h,u-2h} must match the coefficients built from tau1,tau2,qdetT (divided
  by the appropriate ν1 products). Hard-asserts worstTQ, worstCasoratian < 10^-9.*)
Module[{sweep = {{1, 1}, {2, 1}, {1, 2}, {2, 2}, {3, 1}, {1, 3}}, report},
  report = Table[
    Module[{S1 = rep[[1]], S2 = rep[[2]], es, worstTQ = 0, worstCas = 0},
      es = TauEigensystem3[S1, S2];
      KeyValueMap[
       Function[{key, rec},
        Module[{H1 = key[[1]], H2 = key[[2]], n = key[[3]], Qs, u0, rows, minors,
          A1rec, A2rec, A3rec, A1, A2, A3, N1, N1m, N1mm},
         (* (a) TQ residual per Q_k at sample points *)
         Do[worstTQ = Max[worstTQ,
             Abs@Chop[N[BaxOp[S1, S2][H1, H2, n][
                 Function[uu, Qfun[S1, S2][H1, H2, n][k][uu]]][u0], 30]]],
           {k, 1, 3}, {u0, {5/4, 11/3}}];
         (* (b) Casoratian reconstruction of A1,A2,A3 at a sample u0 *)
         u0 = 13/5;
         Qs[uu_] := Table[Qfun[S1, S2][H1, H2, n][k][uu], {k, 1, 3}];
         rows = {Qs[u0 + h], Qs[u0], Qs[u0 - h], Qs[u0 - 2 h]}; (* 4 x 3 *)
         (* minor_i = det of the 3 rows OTHER than row i (i=1..4) *)
         minors = Table[Det[Delete[rows, i]], {i, 1, 4}];
         (* TQ from Casoratian: Q(u+h) - A1 Q(u) + A2 Q(u-h) - A3 Q(u-2h) = 0
            ⇒ A1 = minor2/minor1, A2 = minor3/minor1, A3 = minor4/minor1 (signs from cofactor) *)
         A1rec = minors[[2]]/minors[[1]];
         A2rec = minors[[3]]/minors[[1]];
         A3rec = minors[[4]]/minors[[1]];
         N1 = ν1[S1, S2][u0]; N1m = ν1[S1, S2][u0 - h]; N1mm = ν1[S1, S2][u0 - 2 h];
         A1 = tau1[S1, S2][H1, H2, n][u0]/N1;
         A2 = tau2[S1, S2][H1, H2, n][u0]/(N1 N1m);
         A3 = qdetT[S1, S2][u0]/(N1 N1m N1mm);
         worstCas = Max[worstCas,
            Abs@Chop[N[A1rec - A1, 30]], Abs@Chop[N[A2rec - A2, 30]], Abs@Chop[N[A3rec - A3, 30]]]]],
       es];
      {rep, "nstates" -> Length[es], "worstTQ" -> worstTQ, "worstCasoratian" -> worstCas}],
    {rep, sweep}];
  If[Max[#[[3, 2]] & /@ report, #[[4, 2]] & /@ report] > 10^-9,
    Print["SWEEP FAIL: ", report]; Abort[]];
  Column[report]
]
```
**Output:**
```
Out[109]= \begin{matrix}\left\{\left\{1,1\right\},\text{nstates}\to 9,\text{worstTQ}\to 0,\text{worstCasoratian}\to 0\right\}\\\left\{\left\{2,1\right\},\text{nstates}\to 18,\text{worstTQ}\to 0,\text{worstCasoratian}\to 0\right\}\\\left\{\left\{1,2\right\},\text{nstates}\to 18,\text{worstTQ}\to 0,\text{worstCasoratian}\to 0\right\}\\\left\{\left\{2,2\right\},\text{nstates}\to 36,\text{worstTQ}\to … [657 chars total]
```

## 2026-07-17 14:41:48
**Input:**
```wolfram
Out[109]
```
**Output:**
```
Out= {{1, 1}, nstates -> 9, worstTQ -> 0, worstCasoratian -> 0}
{{2, 1}, nstates -> 18, worstTQ -> 0, worstCasoratian -> 0}
{{1, 2}, nstates -> 18, worstTQ -> 0, worstCasoratian -> 0}
{{2, 2}, nstates -> 36, worstTQ -> 0, worstCasoratian -> 0}
{{3, 1}, nstates -> 30, worstTQ -> 0, worstCasoratian -> 0}
{{1, 3}, nstates -> 30, worstTQ -> 0, worstCasoratian -> 0}
```
