## Deleted Cell 14 — 2026-07-24 05:52:37
**Kind:** markdown
```markdown
Sums `EE[I,J]` over all index tuples $I,J\in\{1,\dots,4\}^4$, i.e. the identity $\det A = \frac{1}{N!}\sum_{I,J} \epsilon_I \epsilon_J\, A_{i_1 j_1}\cdots A_{i_N j_N}$ applied noncommutatively (scalars are pulled through `Dot`, but $L$ never crosses past $\partial_z$/$f$), collapsing the result into a differential operator in $z$ acting on $f(z)$.
```

## Deleted Cell 13 — 2026-07-24 05:52:37
**Kind:** code
```wolfram
EE[II_, JJ_] :=
  1/4! Sum[
     Signature[\[Sigma]] Times @@ (e[#] & /@ Transpose[{II[[\[Sigma]]], JJ}]),
     {\[Sigma], Permutations[{1, 2, 3, 4}]}
    ] MChain[Transpose[{II, JJ}]] /. e[{a_, b_}] :> KroneckerDelta[a, b];
```

## Deleted Cell 12 — 2026-07-24 05:52:37
**Kind:** markdown
```markdown
The $\epsilon$-$\epsilon$ contraction term: sums $\mathrm{sign}(\sigma)$ over permutations $\sigma$ matching $I_\sigma$ to $J$ (zero unless $I$ is a permutation of $J$), weighting the ordered operator chain `MChain[Transpose[{I,J}]]`.
```

## Deleted Cell 11 — 2026-07-24 05:52:37
**Kind:** code
```wolfram
applyMs = {
   MChain[{a__, b__, c__, d__}] :> M[a][M[b][M[c][M[d][f[z]]]]] //.
     Dot[a___, b_ + c_, d___] :> Dot[a, b, d] + Dot[a, c, d]
  };
```

