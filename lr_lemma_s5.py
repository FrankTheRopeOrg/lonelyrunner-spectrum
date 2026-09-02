# -*- coding: utf-8 -*-
"""Lemma 6.5: if U has at least one free repeat direction, it has at least two.

Let w be the tight repeat vector on that direction; up to a signed permutation
w = (a,a,b,c,d,e) with 0 < b < c < d < e -- ten possibilities.  Take any second
generator z, so that the coordinate forms are phi_r = (w_r, z_r).  Since the
|w_r| are pairwise distinct apart from |w_1| = |w_2|, the only triple with
w_p - eps*w_q = 0 is (1,2,+), the direction of w itself.  "At most one free
direction" therefore forces, for each of the other 29 triples,

      (w_p - eps*w_q)*z_r - w_r*z_p + eps*w_r*z_q = 0    for some r,

a union of six hyperplanes in z.  We intersect the 29 unions exactly, branching
on r, keeping solution spaces in RREF so duplicates collapse, and discarding
any branch whose space lies inside <w> (there U would not be 2-dimensional).

Result: empty for all ten w.
"""
import itertools
from fractions import Fraction as F
from lr_classify import canonical_first_vectors, rref

__all__ = ["solutions_with_at_most_one_free_direction", "run"]


def _intersect_hyperplane(B, a):
    """Basis of {z in span(B) : a.z = 0}."""
    k = [sum(F(a[i]) * F(b[i]) for i in range(6)) for b in B]
    nz = [i for i, x in enumerate(k) if x != 0]
    if not nz:
        return B
    j = nz[0]
    return [[F(B[i][t]) - k[i] / k[j] * F(B[j][t]) for t in range(6)]
            for i in range(len(B)) if i != j]


def _inside_span_w(B, w):
    r0 = next(i for i in range(6) if w[i] != 0)
    for b in B:
        lam = F(b[r0], w[r0])
        if any(F(b[i]) != lam * w[i] for i in range(6)):
            return False
    return True


def solutions_with_at_most_one_free_direction(w):
    conds = []
    for p, q in itertools.combinations(range(6), 2):
        for eps in (1, -1):
            if (p, q, eps) == (0, 1, 1):
                continue                      # the direction of w itself
            hyps = set()
            for r in range(6):
                a = [0] * 6
                a[r] += w[p] - eps * w[q]
                a[p] -= w[r]
                a[q] += eps * w[r]
                if any(a):
                    hyps.add(tuple(a))
            conds.append(sorted(hyps))
    conds.sort(key=len)                       # most restrictive first
    I6 = [[1 if i == j else 0 for i in range(6)] for j in range(6)]
    level = {rref(I6): I6}
    for hyps in conds:
        nxt = {}
        for B in level.values():
            for a in hyps:
                B2 = _intersect_hyperplane(B, a)
                if not B2 or _inside_span_w(B2, w):
                    continue
                k2 = rref(B2)
                nxt.setdefault(k2, [list(map(F, r)) for r in k2])
        level = nxt
        if not level:
            break
    return list(level)


def run(verbose=True):
    total = 0
    for w in canonical_first_vectors():
        sols = solutions_with_at_most_one_free_direction(w)
        total += len(sols)
        if verbose:
            print("   w = %-24s solutions with <=1 free direction: %d"
                  % (str(w), len(sols)))
    return total
