## Deleted Cell 41 — 2026-07-16 08:27:46
**Kind:** code
```wolfram
1+1
```

## Deleted Cell 45 — 2026-07-16 08:42:52
**Kind:** code
```wolfram
1+1
```

## Deleted Cell 45 — 2026-07-16 08:46:52
**Kind:** code
```wolfram
1+1
```

## Deleted Cell 45 — 2026-07-16 09:00:25
**Kind:** code
```wolfram
1+1
```

## Deleted Cell 45 — 2026-07-16 09:02:16
**Kind:** code
```wolfram
1+1
```

## Deleted Cell 47 — 2026-07-16 09:58:22
**Kind:** code
```wolfram

```

## Deleted Cell 45 — 2026-07-16 09:58:22
**Kind:** code
```wolfram

```

## Deleted Cell 91 — 2026-07-18 11:28:31
**Kind:** code
```wolfram
(*Discovery step (sharp, λ2=1 rep, mirrors Cell 71): rep (2,1), all 6 states — solve φ per state
  assuming exact match against OPg2 (both linear in φ since λ2=1), then check the ratio
  α2_state = D_state/c1_state is the SAME constant across all 6 states.*)
statesR21 = statesOf[2,1];
alphaTab2 = Table[
  Fl = FLdet2[2,1][st[[1]],st[[2]]][k]//Simplify;
  c0 = Fl/.k->0;
  Dcoef = Coefficient[Fl,k,1];
  Og = OPg2[2,1][st[[1]],st[[2]]][φ]//Simplify;
  cc0 = Og/.φ->0;
  cc1 = Coefficient[Og,φ,1];
  If[Chop[c0-1,10^-12]=!=0 || Chop[cc0-1,10^-12]=!=0, Abort[]];
  {st, Dcoef, cc1, Dcoef/cc1}
,{st,statesR21}];
alpha2Vals = alphaTab2[[All,4]];
alpha2Spread = Max[Abs[N[alpha2Vals-alpha2Vals[[1]]]]];
If[alpha2Spread > 10^-12, Abort[]];
alpha2Val = N[alpha2Vals[[1]]];
{"alphaTab2 {state,D,c1,alpha}",alphaTab2,"alpha2Val",alpha2Val,"spread",alpha2Spread}
```

