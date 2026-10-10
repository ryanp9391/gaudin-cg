# Do the t's commute with their Hermitian conjugates? An abstract computation

These notes cover `Hermitian-structure-abstract.wb`, from 2026-10-01. Everything is abstract, with no representation: the gl(3) generators E⁽ˢ⁾_ij satisfy [E_ij, E_kl] = δ_jk E_il − δ_li E_kj on each site, different sites commute, and E_ij† = E_ji (the EH convention). The notebook:
- normal-orders words by (site, i, j);
- defines † as reversing products, sending E_ij → E_ji, and conjugating scalars;
- keeps the conjugated scalars as separate symbols (ub, θb[s], hb, zb[i]), so any reality conditions come out explicitly.

It uses the new conventions: ℒ_ij(x) = x δ_ij + h E_ij, T = ℒ⁽¹⁾(u−θ₁)⋯ℒ⁽ᴸ⁾(u−θ_L), the twist diagonal on the right, t₁ = Σ z_i T_ii and t₂ = Σ_{i<j} z_i z_j [T_ii(u)T_jj(u−h) − T_ji(u)T_ij(u−h)].

**Sanity checks.** [t₁(u),t₁(v)] = 0 and [t₁(u),t₂(v)] = 0 hold at L=1 and L=2, which tests the normal ordering and the quantum minor together.

## Results

Write t̄_a for t_a with every scalar conjugated.

1. **L=1.** t₁† = t̄₁ and t₂† = t̄₂ exactly. Off-diagonal generators enter t₂ only as E_jiE_ij, which is self-conjugate. So for real parameters the t's are Hermitian, and they commute with their conjugates.
2. **L=2, the conjugate of t₁.**

   t₁† − t̄₁ = h̄² Σ_{i≠k} (z̄_k − z̄_i) E⁽¹⁾_ik E⁽²⁾_ki

   Conjugation moves the twist weight of the two-site hopping term from one index to the other. (t₂† − t̄₂ has 240 terms.)
3. **L=2, the commutator.** For real parameters, write:
   - H_s = Σ z_i E⁽ˢ⁾_ii;
   - C_z = Σ_{i,k} z_i E⁽¹⁾_ik E⁽²⁾_ki and C_z† = Σ_{i,k} z_k E⁽¹⁾_ik E⁽²⁾_ki.

   Then t₁(u) = c(u) + hu(H₁+H₂) − h(θ₂H₁+θ₁H₂) + h²C_z, and

   [t₁(u), t₁(v)†] = −h³ θ₁₂ Σ_{i,k} (z_i−z_k)² E⁽¹⁾_ik E⁽²⁾_ki + h⁴ [C_z, C_z†]

   This is independent of u and v and checked exactly against the brute-force commutator. Two facts make it work: H₁+H₂ commutes with everything here, and [H₂, E⁽¹⁾_ik E⁽²⁾_ki] = (z_k−z_i) E⁽¹⁾_ik E⁽²⁾_ki.
4. **When does it vanish?**
   - [C_z, C_z†] contains genuinely cubic words, and these carry no θ dependence. So no choice of θ's can cancel the commutator.
   - It is zero only for an **untwisted** chain (all z_i equal), where C_z = C_z†.
   - It is not zero for z₁ = z₂ ≠ z₃, for a unitary twist z̄_i = 1/z_i, or for imaginary h, u and θ.

## Interpretation

With a generic diagonal twist at L ≥ 2, the EH structure does not make t₁ normal. What holds instead is t_a(u)ᵀ = t_a of the site-reversed chain; it was checked numerically for a = 1, 2, 3 in `L2-diagonal-frame-overlap-EH.wb`, Section 1b. So the natural pairing is † combined with site reversal, which brings in the R-matrix that swaps the two sites. That has not been explored.

## Not done

- [t₂, t₂†] at L=2. The conclusion already follows from t₁.
- The "† ∘ site reversal" structure.
