# gl(3) single-site generators — Python port design

**Date:** 2026-07-18
**Subproject:** XXX (gl(3))
**Source:** `Paul/Mathematica/XXX/Experiments/su3_V2_general.wb`, cell 0 (`ps`, `dim`, `J0`, `Jp`,
`Jm`, `JJ`) and cell 1 (`Ee[λ_][i_,j_] := JJ[i,j,λ]`)
**Status:** design approved 2026-07-18

## Goal

Port the single-site Gelfand-Tsetlin (GT) construction of the gl(3) generators `E_ij` from the
Mathematica notebook into an idiomatic Python module, as a first worked example of doing
cross-check math in Python in this repo (see `HANDOFF.md`/session history for the broader
"build a Python side for cross-checks" context). The two-site tensor-product layer (`EE`, Lax
matrix, monodromy, Baxter equation) is explicitly **out of scope for this pass** — it will reuse
this module once it exists, per the user, but is not part of this deliverable.

The Mathematica source for this machinery lives in `Experiments/` (not yet promoted to `Clean/`),
but per `HANDOFF.md` it has already passed a full validation sweep (Tasks 1–6, hard-asserts,
102-state full-sweep check) — it is trusted as ground truth for this port.

## Fixed decisions (confirmed with the user)

- **Scope:** single-site generators only (`ps`/`dim`/`J0`/`Jp`/`Jm`/`JJ`/`Ee` from cell 0–1). No
  two-site embedding (`EE`), no Lax/monodromy/Baxter machinery.
- **Interface:** a small class, not a literal function-for-function mirror of the Mathematica
  names — this is a "port," not a transliteration:
  ```python
  rep = GL3Rep((2, 1, 0))
  rep.dim            # 8
  rep.patterns        # list of GT patterns, basis order
  rep.generator(1, 2)  # sympy.Matrix, E_12, dim x dim, exact rationals
  ```
- **Arithmetic:** `sympy.Rational` / `sympy.Matrix` throughout — GT matrix elements are
  generically rational, not integer, and exactness matters for comparison against the source.
- **Verification:** internal algebraic checks only (no live Mathematica cross-check for this
  pass), run via `pytest`.
- **Generality:** implemented specifically for gl(3), not a general gl(n) engine. The one
  non-adjacent generator pair (`E_13`, `E_31`) is built via a single commutator step
  (`E_13 = [E_12, E_23]`) rather than a general recursive-commutator routine — gl(3) never needs
  more than one step, so a fully general recursion would be unused generality.

## Data representation

A GT pattern for `λ = (λ1, λ2, λ3)` (integers, `λ1 ≥ λ2 ≥ λ3 ≥ 0`) is a tuple of 3 tuples:
`(row1, row2, row3)` where `row_k` has `k` integer entries, satisfying the interlacing condition
`row_{k}[i] ≤ row_{k+1}[i] ≤ row_{k}[i-1]` (1-indexed) between consecutive rows, and `row3 == λ`.

**Row ordering matches the Mathematica source exactly**: index 0 is the bottom singleton row,
index 2 (last) is the top row `λ`. This is the opposite of what a fresh reader might expect
(reading "top-down"), but re-ordering rows relative to the source is exactly the kind of silent
change that introduces an indexing bug in a port — so the convention is kept identical to
`pat[[k]]` in the notebook (`pat[[1]]` = bottom, `pat[[3]]` = top = λ), just re-based to Python's
0-indexing. A one-line comment on the pattern-generating function documents this explicitly.

Patterns are generated once per `GL3Rep` instance (in `__init__`), via the same recursive
interlacing construction as `ps` in the notebook: given the current top row, enumerate all valid
next-rows-down via `range` intersections, recurse until row length 1.

## Algorithm

- **`dim`**: `len(patterns)`.
- **Weight / diagonal generator** (`J0` in the source): for row `k`, weight component
  `w_k(pat) = sum(row_k) - sum(row_{k-1})` (with `row_0` treated as empty, sum 0). `E_kk` is the
  diagonal matrix of `w_k` over the pattern basis.
- **Raising `E_{k,k+1}`** (`Jp` in the source): for each pattern and each row-`k` entry `i`, the
  GT raising formula gives a coefficient (product over adjacent-row entries, divided by a product
  over same-row entries) and a target pattern (row `k`, entry `i`, incremented by 1). Represented
  directly as a `{target_pattern: coefficient}` dict per source pattern — this replaces the
  Mathematica source's `v[...]` formal-symbol/`Coefficient[...]` bookkeeping trick, which exists
  in Mathematica to extract matrix entries from a symbolic expression; Python builds the dict
  directly, which is the equivalent computation without a symbolic-expression round-trip.
- **Lowering `E_{k,k-1}`** (`Jm`): same shape, decrementing row `k` entry `i`, coefficients from
  the corresponding source formula.
- **Non-adjacent (`E_13`, `E_31`)**: `E_13 = E_12 @ E_23 - E_23 @ E_12` (matrix commutator), and
  symmetrically `E_31 = E_32 @ E_21 - E_21 @ E_32`.
- **Memoization**: `generator(i, j)` caches results in an instance dict, mirroring the source's
  `JJ[k,l,rep] := JJ[k,l,rep] = ...` self-memoizing pattern — avoids rebuilding `E_12`/`E_23`
  every time `E_13` is requested.

## Testing (pytest, `Paul/Python/test_gl3_reps.py`)

All checks are internal/algebraic — no comparison against a live Mathematica evaluation in this
pass:

1. **Dimension formula**: `rep.dim` matches the closed-form Weyl dimension formula for gl(3),
   `((λ1−λ2+1)(λ2−λ3+1)(λ1−λ3+2))/2`, for a handful of reps (fundamental `(1,0,0)`, symmetric
   `(2,0,0)`, antisymmetric `(1,1,0)`, adjoint `(2,1,0)`). This is an independent formula, not
   derived from the same GT code, so it's a real cross-check rather than circular.
2. **Commutation relations**: `[E_ij, E_kl] = δ_jk E_il − δ_li E_kj` for all `i,j,k,l ∈ {1,2,3}`,
   across the same set of reps.
3. **Highest-weight eigenvalues**: `E_ii` acting on the highest-weight basis pattern (the first
   pattern in basis order, which is `λ` in every row) gives eigenvalue `λ_i`.

## Files

Flat layout in `Paul/Python/` (no subproject split yet, per the user):

- `gl3_reps.py` — the `GL3Rep` class.
- `test_gl3_reps.py` — pytest checks above.
- `pytest` added as a dev dependency via `uv add --dev pytest`.

## Out of scope

- Two-site tensor-product layer (`EE`, Lax matrix `L`, monodromy `T`/`Tg`/`TG`, quantum minors,
  Baxter/TQ machinery) — planned for a follow-up once this module exists.
- General gl(n) engine — this is gl(3)-specific.
- Live cross-check against a running Mathematica kernel for this pass (available locally via
  `WolframKernel`/`math` CLI if needed later, but not used here per the user's choice of
  internal-checks-only verification).

## Open questions

None — scope, interface, arithmetic, verification strategy, and file layout all confirmed with
the user during brainstorming.
