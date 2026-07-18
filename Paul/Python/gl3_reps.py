"""GT-basis construction of gl(3) generators E_ij in an irrep lam.

Ported from Paul/Mathematica/XXX/Experiments/su3_V2_general.wb (cell 0-1):
ps/dim/J0/Jp/Jm/JJ (single-site GT machinery) and Ee[lam][i,j] := JJ[i,j,lam].
"""
from itertools import product

import sympy


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


def _weight(k, pattern):
    """w_k(pattern) = sum(row_k) - sum(row_{k-1}), 1-indexed rows; row_0 is empty."""
    row_k = pattern[k - 1]
    row_prev = pattern[k - 2] if k > 1 else ()
    return sum(row_k) - sum(row_prev)


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
        if j == i + 1:
            return self._adjacent_up(i)
        if j == i - 1:
            return self._adjacent_down(j)
        raise ValueError(f"generator({i}, {j}) not yet implemented")

    def _adjacent_up(self, k):
        M = sympy.zeros(self.dim, self.dim)
        for col, pat in enumerate(self.patterns):
            for target, coeff in _raise(k, pat).items():
                # _raise's raw formula can propose a target row outside the
                # interlacing range (e.g. off the top of the GT pattern) with
                # a nonzero coefficient. The Mathematica source relies on
                # Coefficient[...] silently discarding such terms because
                # they never match a v[pat] for pat in ps[{rep}]; the
                # equivalent here is to skip targets absent from self._index.
                row_index = self._index.get(target)
                if row_index is not None:
                    M[row_index, col] = coeff
        return M

    def _adjacent_down(self, k):
        M = sympy.zeros(self.dim, self.dim)
        for col, pat in enumerate(self.patterns):
            for target, coeff in _lower(k, pat).items():
                # Same boundary-term filtering as _adjacent_up, mirroring the
                # Mathematica Coefficient[] extraction (see comment above).
                row_index = self._index.get(target)
                if row_index is not None:
                    M[row_index, col] = coeff
        return M
