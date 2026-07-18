# gl(3) Single-Site Generators Python Port Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Port the Gelfand-Tsetlin (GT) construction of gl(3) generators `E_ij` from
`Paul/Mathematica/XXX/Experiments/su3_V2_general.wb` (cell 0–1) into a Python class `GL3Rep`, as
the first worked example of a cross-check module in `Paul/Python/`.

**Architecture:** A GT pattern is a tuple of 3 row-tuples (row `k` has `k` entries), generated
recursively via the same interlacing-range construction as the notebook's `ps`. `GL3Rep(lam)`
builds the pattern list once, then computes generator matrices on demand: diagonal entries from a
weight formula (`J0`), the two adjacent raising/lowering pairs `(1,2)/(2,1)` and `(2,3)/(3,2)`
from the GT raising/lowering formulas (`Jp`/`Jm`), and the two non-adjacent generators `(1,3)/(3,1)`
from a single commutator step, matching the source's `JJ` dispatch exactly.

**Tech Stack:** Python 3.10, `sympy` (exact `Rational`/`Matrix` arithmetic), `pytest`, managed via
`uv` in `Paul/Python/` (existing `pyproject.toml`/`uv.lock`, see prior session).

## Global Constraints

- Scope: single-site generators only (`ps`/`dim`/`J0`/`Jp`/`Jm`/`JJ`/`Ee` from the notebook). No
  two-site tensor-product layer (`EE`, Lax, monodromy, Baxter) in this plan.
- gl(3)-specific implementation, not a general gl(n) engine — the one non-adjacent pair per
  direction (`E_13`, `E_31`) is a single hardcoded commutator step, not a recursive routine.
- All matrix arithmetic via `sympy.Rational` / `sympy.Matrix` — no floats, no plain Python `Fraction`.
- GT pattern row ordering matches the Mathematica source exactly: row index 0 (Python) is the
  bottom singleton row, the last row index is `lam` itself (top row) — this mirrors the
  notebook's `pat[[k]]` convention (row `k` has `k` entries) and must not be reversed.
- Verification is internal/algebraic only in this pass (dimension formula, commutation relations,
  highest-weight eigenvalues) — no live Mathematica cross-check.
- Files live flat in `Paul/Python/` (no subproject split): `gl3_reps.py`, `test_gl3_reps.py`.
- Design spec: `docs/superpowers/specs/2026-07-18-gl3-python-port-design.md`.

---

### Task 1: GT pattern enumeration + dimension check

**Files:**
- Create: `Paul/Python/gl3_reps.py`
- Create: `Paul/Python/test_gl3_reps.py`

**Interfaces:**
- Produces: `gt_patterns(lam: tuple[int, int, int]) -> list[tuple[tuple[int, ...], ...]]` — all GT
  patterns for top row `lam`, each pattern a 3-tuple of rows (lengths 1, 2, 3), row 0 = bottom.

- [ ] **Step 1: Add `pytest` as a dev dependency**

Run (from `Paul/Python/`):
```bash
cd Paul/Python && uv add --dev pytest
```
Expected: `pyproject.toml` gains a `[dependency-groups]` (or `[tool.uv]`) `dev = ["pytest>=..."]`
entry, `uv.lock` updates, command exits 0.

- [ ] **Step 2: Write the failing test**

Create `Paul/Python/test_gl3_reps.py`:

```python
import pytest

from gl3_reps import gt_patterns

REPS = [(1, 0, 0), (2, 0, 0), (1, 1, 0), (2, 1, 0)]


def weyl_dim(lam):
    l1, l2, l3 = lam
    return (l1 - l2 + 1) * (l2 - l3 + 1) * (l1 - l3 + 2) // 2


@pytest.mark.parametrize("lam", REPS)
def test_dimension_matches_weyl_formula(lam):
    patterns = gt_patterns(lam)
    assert len(patterns) == weyl_dim(lam)
    assert len(set(patterns)) == len(patterns)  # no duplicate patterns
```

- [ ] **Step 3: Run test to verify it fails**

Run:
```bash
cd Paul/Python && uv run pytest test_gl3_reps.py -v
```
Expected: FAIL / ERROR — `ModuleNotFoundError: No module named 'gl3_reps'` (the module doesn't
exist yet).

- [ ] **Step 4: Write the implementation**

Create `Paul/Python/gl3_reps.py`:

```python
"""GT-basis construction of gl(3) generators E_ij in an irrep lam.

Ported from Paul/Mathematica/XXX/Experiments/su3_V2_general.wb (cell 0-1):
ps/dim/J0/Jp/Jm/JJ (single-site GT machinery) and Ee[lam][i,j] := JJ[i,j,lam].
"""
from itertools import product


def gt_patterns(lam):
    """All Gelfand-Tsetlin patterns with top row `lam`.

    A pattern is a tuple of rows (row of length 1, ..., row of length len(lam)).
    Row ordering matches the Mathematica source: index 0 is the bottom
    singleton row, the last index is `lam` itself -- kept identical to the
    notebook's pat[[k]] convention (row k has k entries) rather than reversed,
    since re-ordering rows relative to the source is exactly the kind of
    silent change that introduces an indexing bug in a port.
    """
    return _extend((tuple(lam),))


def _extend(rows):
    top = rows[0]
    if len(top) == 1:
        return [rows]
    ranges = [range(top[i + 1], top[i] + 1) for i in range(len(top) - 1)]
    patterns = []
    for combo in product(*ranges):
        patterns.extend(_extend((combo,) + rows))
    return patterns
```

- [ ] **Step 5: Run test to verify it passes**

Run:
```bash
cd Paul/Python && uv run pytest test_gl3_reps.py -v
```
Expected: `test_dimension_matches_weyl_formula[lam0]` through `[lam3]` all PASS (4 passed).

- [ ] **Step 6: Commit**

```bash
cd /home/paul/git_projects/gaudin-CG
git add Paul/Python/gl3_reps.py Paul/Python/test_gl3_reps.py Paul/Python/pyproject.toml Paul/Python/uv.lock
git commit -m "$(cat <<'EOF'
Add GT pattern enumeration for gl(3) reps (Python port, Task 1)

Ports ps[rep] from su3_V2_general.wb: recursive interlacing
construction of Gelfand-Tsetlin patterns, verified against the
independent Weyl dimension formula.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>
EOF
)"
```

---

### Task 2: `GL3Rep` class skeleton + diagonal generator

**Files:**
- Modify: `Paul/Python/gl3_reps.py`
- Modify: `Paul/Python/test_gl3_reps.py`

**Interfaces:**
- Consumes: `gt_patterns(lam)` from Task 1.
- Produces: `_weight(k: int, pattern) -> int`; `class GL3Rep` with `__init__(self, lam)`,
  `.lam` (tuple), `.patterns` (list), `.dim` (property, int), `.generator(i: int, j: int) ->
  sympy.Matrix` (diagonal case `i == j` only — other `(i, j)` raise `ValueError` until Tasks 3–4).

- [ ] **Step 1: Write the failing test**

Append to `Paul/Python/test_gl3_reps.py` (add import at top, then the test):

```python
import sympy

from gl3_reps import GL3Rep, gt_patterns  # update the existing import line


def highest_weight_pattern(lam):
    return tuple(tuple(lam[:k]) for k in range(1, len(lam) + 1))


@pytest.mark.parametrize("lam", REPS)
def test_diagonal_generator_highest_weight_eigenvalue(lam):
    rep = GL3Rep(lam)
    hw = highest_weight_pattern(lam)
    hw_index = rep.patterns.index(hw)
    for k in range(1, 4):
        E_kk = rep.generator(k, k)
        assert E_kk[hw_index, hw_index] == lam[k - 1]
```

(Replace the file's `from gl3_reps import gt_patterns` line with
`from gl3_reps import GL3Rep, gt_patterns`, and add `import sympy` near the top with the other
imports — `sympy` isn't used directly in this test yet but the next tasks will need it in this
file.)

- [ ] **Step 2: Run test to verify it fails**

Run:
```bash
cd Paul/Python && uv run pytest test_gl3_reps.py -v
```
Expected: the new `test_diagonal_generator_highest_weight_eigenvalue` cases FAIL/ERROR with
`ImportError: cannot import name 'GL3Rep' from 'gl3_reps'`; the Task 1 dimension test still PASSES.

- [ ] **Step 3: Write the implementation**

In `Paul/Python/gl3_reps.py`, add `import sympy` directly below the existing
`from itertools import product` line at the top of the file. Then append the following after
`_extend`:

```python
def _weight(k, pattern):
    """w_k(pattern) = sum(row_k) - sum(row_{k-1}), 1-indexed rows; row_0 is empty."""
    row_k = pattern[k - 1]
    row_prev = pattern[k - 2] if k > 1 else ()
    return sum(row_k) - sum(row_prev)


class GL3Rep:
    """Irreducible gl(3) representation in the Gelfand-Tsetlin basis.

    Ported from Paul/Mathematica/XXX/Experiments/su3_V2_general.wb
    (ps/dim/J0/Jp/Jm/JJ, cell 0-1).
    """

    def __init__(self, lam):
        lam = tuple(lam)
        if len(lam) != 3 or list(lam) != sorted(lam, reverse=True) or lam[-1] < 0:
            raise ValueError(
                f"lam must be a weakly decreasing triple of non-negative integers, got {lam}"
            )
        self.lam = lam
        self.patterns = gt_patterns(lam)
        self._index = {pattern: i for i, pattern in enumerate(self.patterns)}
        self._cache = {}

    @property
    def dim(self):
        return len(self.patterns)

    def generator(self, i, j):
        key = (i, j)
        if key not in self._cache:
            self._cache[key] = self._build_generator(i, j)
        return self._cache[key]

    def _build_generator(self, i, j):
        if i == j:
            return sympy.diag(*[_weight(i, pat) for pat in self.patterns])
        raise ValueError(f"generator({i}, {j}) not yet implemented")
```

- [ ] **Step 4: Run test to verify it passes**

Run:
```bash
cd Paul/Python && uv run pytest test_gl3_reps.py -v
```
Expected: all tests PASS (Task 1's dimension test + Task 2's diagonal/highest-weight test, 8 total
across both parametrized tests).

- [ ] **Step 5: Commit**

```bash
cd /home/paul/git_projects/gaudin-CG
git add Paul/Python/gl3_reps.py Paul/Python/test_gl3_reps.py
git commit -m "$(cat <<'EOF'
Add GL3Rep class + diagonal generator (Python port, Task 2)

Ports J0[k] (GT weight formula) as _weight, and wires up the
GL3Rep.generator(i,j) dispatch for the diagonal case. Verified via
highest-weight eigenvalues E_ii|hw> = lam_i|hw>.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>
EOF
)"
```

---

### Task 3: Adjacent raising/lowering generators

**Files:**
- Modify: `Paul/Python/gl3_reps.py`
- Modify: `Paul/Python/test_gl3_reps.py`

**Interfaces:**
- Consumes: `GL3Rep`, `_weight` from Task 2.
- Produces: `_raise(k: int, pattern) -> dict` (E_{k,k+1} matrix elements from `pattern`),
  `_lower(k: int, pattern) -> dict` (E_{k+1,k} matrix elements from `pattern`);
  `GL3Rep.generator(i, j)` now also handles `j == i + 1` and `j == i - 1`.

- [ ] **Step 1: Write the failing test**

Append to `Paul/Python/test_gl3_reps.py`:

```python
@pytest.mark.parametrize("lam", REPS)
def test_adjacent_commutators(lam):
    rep = GL3Rep(lam)
    comm12 = rep.generator(1, 2) * rep.generator(2, 1) - rep.generator(2, 1) * rep.generator(1, 2)
    assert comm12 == rep.generator(1, 1) - rep.generator(2, 2)
    comm23 = rep.generator(2, 3) * rep.generator(3, 2) - rep.generator(3, 2) * rep.generator(2, 3)
    assert comm23 == rep.generator(2, 2) - rep.generator(3, 3)
```

- [ ] **Step 2: Run test to verify it fails**

Run:
```bash
cd Paul/Python && uv run pytest test_gl3_reps.py -v
```
Expected: `test_adjacent_commutators` FAILS with `ValueError: generator(1, 2) not yet
implemented`; earlier tests still PASS.

- [ ] **Step 3: Write the implementation**

Add to `Paul/Python/gl3_reps.py` (after `_weight`, before `class GL3Rep`):

```python
def _raise(k, pattern):
    """Coefficients of E_{k,k+1} acting on `pattern`: {new_pattern: coeff}. Ports Jp[k][pat]."""
    row_k = pattern[k - 1]
    row_next = pattern[k]
    terms = {}
    for i in range(1, k + 1):
        num = 1
        for j in range(1, k + 2):
            num *= row_k[i - 1] - i - row_next[j - 1] + j
        den = 1
        for j in range(1, k + 1):
            if j == i:
                continue
            den *= row_k[i - 1] - i - row_k[j - 1] + j
        coeff = -sympy.Rational(num, den)
        if coeff == 0:
            continue
        new_row_k = list(row_k)
        new_row_k[i - 1] += 1
        new_pattern = pattern[: k - 1] + (tuple(new_row_k),) + pattern[k:]
        terms[new_pattern] = terms.get(new_pattern, 0) + coeff
    return terms


def _lower(k, pattern):
    """Coefficients of E_{k+1,k} acting on `pattern`: {new_pattern: coeff}. Ports Jm[k][pat]."""
    row_k = pattern[k - 1]
    row_prev = pattern[k - 2] if k > 1 else ()
    terms = {}
    for i in range(1, k + 1):
        num = 1
        for j in range(1, k):
            num *= row_k[i - 1] - i - row_prev[j - 1] + j
        den = 1
        for j in range(1, k + 1):
            if j == i:
                continue
            den *= row_k[i - 1] - i - row_k[j - 1] + j
        coeff = sympy.Rational(num, den)
        if coeff == 0:
            continue
        new_row_k = list(row_k)
        new_row_k[i - 1] -= 1
        new_pattern = pattern[: k - 1] + (tuple(new_row_k),) + pattern[k:]
        terms[new_pattern] = terms.get(new_pattern, 0) + coeff
    return terms
```

Replace the `GL3Rep._build_generator` method with:

```python
    def _build_generator(self, i, j):
        if i == j:
            return sympy.diag(*[_weight(i, pat) for pat in self.patterns])
        if j == i + 1:
            return self._adjacent_up(i)
        if j == i - 1:
            return self._adjacent_down(j)
        raise ValueError(f"generator({i}, {j}) not yet implemented")

    def _adjacent_up(self, k):
        M = sympy.zeros(self.dim, self.dim)
        for col, pat in enumerate(self.patterns):
            for target, coeff in _raise(k, pat).items():
                M[self._index[target], col] = coeff
        return M

    def _adjacent_down(self, k):
        M = sympy.zeros(self.dim, self.dim)
        for col, pat in enumerate(self.patterns):
            for target, coeff in _lower(k, pat).items():
                M[self._index[target], col] = coeff
        return M
```

- [ ] **Step 4: Run test to verify it passes**

Run:
```bash
cd Paul/Python && uv run pytest test_gl3_reps.py -v
```
Expected: all tests PASS, including `test_adjacent_commutators[lam0]` through `[lam3]`.

- [ ] **Step 5: Commit**

```bash
cd /home/paul/git_projects/gaudin-CG
git add Paul/Python/gl3_reps.py Paul/Python/test_gl3_reps.py
git commit -m "$(cat <<'EOF'
Add adjacent raising/lowering generators (Python port, Task 3)

Ports Jp[k]/Jm[k] (GT raising/lowering matrix-element formulas) as
_raise/_lower, replacing the notebook's symbolic v[...]/Coefficient
bookkeeping with a direct {pattern: coefficient} dict. Wires E_12,
E_21, E_23, E_32 into GL3Rep.generator. Verified via
[E_{k,k+1}, E_{k+1,k}] = E_kk - E_{k+1,k+1}.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>
EOF
)"
```

---

### Task 4: Non-adjacent generators + full commutation relation suite

**Files:**
- Modify: `Paul/Python/gl3_reps.py`
- Modify: `Paul/Python/test_gl3_reps.py`

**Interfaces:**
- Consumes: `GL3Rep.generator` (diagonal + adjacent cases) from Tasks 2–3.
- Produces: `GL3Rep.generator(i, j)` fully defined for all `i, j` in `1..3` (adds `(1, 3)` and
  `(3, 1)`).

- [ ] **Step 1: Write the failing test**

Append to `Paul/Python/test_gl3_reps.py`:

```python
@pytest.mark.parametrize("lam", REPS)
def test_full_commutation_relations(lam):
    rep = GL3Rep(lam)
    n = rep.dim
    for i in range(1, 4):
        for j in range(1, 4):
            for k in range(1, 4):
                for l in range(1, 4):
                    lhs = rep.generator(i, j) * rep.generator(k, l) - rep.generator(
                        k, l
                    ) * rep.generator(i, j)
                    rhs = sympy.zeros(n, n)
                    if j == k:
                        rhs += rep.generator(i, l)
                    if l == i:
                        rhs -= rep.generator(k, j)
                    assert lhs == rhs, f"[E_{i}{j}, E_{k}{l}] failed for lam={lam}"
```

- [ ] **Step 2: Run test to verify it fails**

Run:
```bash
cd Paul/Python && uv run pytest test_gl3_reps.py -v
```
Expected: `test_full_commutation_relations` FAILS with `ValueError: generator(1, 3) not yet
implemented` (or `(3, 1)`); earlier tests still PASS.

- [ ] **Step 3: Write the implementation**

In `Paul/Python/gl3_reps.py`, replace the `_build_generator` method's final `raise` line with the
two non-adjacent cases:

```python
    def _build_generator(self, i, j):
        if i == j:
            return sympy.diag(*[_weight(i, pat) for pat in self.patterns])
        if j == i + 1:
            return self._adjacent_up(i)
        if j == i - 1:
            return self._adjacent_down(j)
        if (i, j) == (1, 3):
            return self.generator(1, 2) * self.generator(2, 3) - self.generator(
                2, 3
            ) * self.generator(1, 2)
        if (i, j) == (3, 1):
            return self.generator(3, 2) * self.generator(2, 1) - self.generator(
                2, 1
            ) * self.generator(3, 2)
        raise ValueError(f"i, j must be in 1..3, got ({i}, {j})")
```

- [ ] **Step 4: Run test to verify it passes**

Run:
```bash
cd Paul/Python && uv run pytest test_gl3_reps.py -v
```
Expected: all tests PASS (dimension, highest-weight, adjacent-commutator, and full
commutation-relation tests, across all 4 reps).

- [ ] **Step 5: Commit**

```bash
cd /home/paul/git_projects/gaudin-CG
git add Paul/Python/gl3_reps.py Paul/Python/test_gl3_reps.py
git commit -m "$(cat <<'EOF'
Add non-adjacent generators, complete GL3Rep (Python port, Task 4)

Ports the JJ commutator-recursion cases (E_13 = [E_12,E_23], E_31 =
[E_32,E_21]), completing generator(i,j) for all i,j in 1..3. Verified
against the full gl(3) commutation relations
[E_ij,E_kl] = delta_jk E_il - delta_li E_kj across all i,j,k,l and all
four test reps -- the acceptance check from the design spec.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>
EOF
)"
```
