# -*- coding: utf-8 -*-
"""Maximum loneliness of an integer speed tuple, computed exactly.

ML(v_1,...,v_n) = max_t min_i ||t v_i||   over t in [0,1),  ||x|| = dist(x, Z).

Each ||t v_i|| is piecewise linear in t with kinks at the multiples of
1/(2 v_i).  Between consecutive kinks every ||t v_i|| is affine, so their
minimum is concave there and its maximum on that interval is attained at an
endpoint or at a crossing of two of the affine pieces.  Enumerating those
candidates gives the exact value as a Fraction.

If any speed is zero the corresponding runner never leaves the observer and
ML = 0 by definition.
"""
from fractions import Fraction as F
from itertools import combinations

__all__ = ["ml_tuple", "is_tight_quintuple", "TIGHT_QUINTUPLES"]

TIGHT_QUINTUPLES = [(1, 2, 3, 4, 5), (1, 3, 4, 5, 9)]


def _nrm(x):
    y = x - int(x)
    if y < 0:
        y += 1
    return min(y, 1 - y)


def ml_tuple(speeds):
    """Exact ML of a tuple of integers (repeats and signs are irrelevant)."""
    if any(s == 0 for s in speeds):
        return F(0)
    vs = sorted({abs(s) for s in speeds})
    kinks = {F(0), F(1)}
    for v in vs:
        for m in range(2 * v + 1):
            kinks.add(F(m, 2 * v))
    kinks = sorted(kinks)
    best = F(0)
    for a, b in zip(kinks, kinks[1:]):
        mid = (a + b) / 2
        lines = []                      # (slope, intercept) of ||t v|| on [a,b]
        for v in vs:
            y = mid * v - int(mid * v)
            if y < 0:
                y += 1
            s = 1 if y < F(1, 2) else -1
            lines.append((s * v, _nrm(mid * v) - s * v * mid))
        cands = [a, b]
        for i, j in combinations(range(len(lines)), 2):
            m1, c1 = lines[i]
            m2, c2 = lines[j]
            if m1 != m2:
                t = F(c2 - c1) / (m1 - m2)
                if a <= t <= b:
                    cands.append(t)
        for t in cands:
            val = min(m * t + c for m, c in lines)
            if val > best:
                best = val
    return best


def is_tight_quintuple(speeds):
    """True iff the five distinct positive speeds are tight, i.e. ML = 1/6."""
    vs = sorted({abs(s) for s in speeds})
    return len(vs) == 5 and ml_tuple(vs) == F(1, 6)
