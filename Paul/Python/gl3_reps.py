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
