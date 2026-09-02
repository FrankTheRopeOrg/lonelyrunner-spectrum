# -*- coding: utf-8 -*-
"""Lemma 6.6 (closes Hypothesis (H)).  No 2-dimensional subtorus of (R/Z)^6
whose six coordinate forms are pairwise non-proportional can have all thirty
repeat directions carrying a vanishing coordinate.

Proof.  Suppose it does.  By the three-term Pluecker relation the condition is
    phi_p - eps*phi_q  proportional to some phi_m   for all p<q, eps = +-1.
Normalise by GL_2(Q) -- legitimate, every D_ij just scales by det g -- so that
    phi_1 = (1,0),  phi_2 = (0,1).
The pair {1,2} forces two of the others to be proportional to (1,1) and (1,-1);
relabel phi_3 = c(1,1), phi_4 = d(1,-1) with c, d nonzero.  Only phi_5, phi_6
are then free: the set S of the six directions is
    { 0, oo, 1, -1 }  together with at most TWO further points.
(Affine coordinate xi = second/first; 0 = [1:0], oo = [0:1], 1 = [1:1],
-1 = [1:-1].)

Pairs {1,3} and {2,3} give the four directions
    t, 1/t   with t = (c-1)/c,        s, 1/s   with s = (1+c)/c.
Pairs {1,4} and {2,4} give
    t', 1/t' with t' = (1-d)/d,       s', 1/s' with s' = -(1+d)/d.

Each quadruple consists of four DISTINCT points (t = s, t = 1/s, t = 1/t and
s = 1/s are each impossible or force an excluded value of the parameter), and
each lies in {0, oo, 1, -1} only for the four special values of the parameter.
Hence c and d must both lie in {1, 1/2, -1, -1/2}, and then

    c = +-1   -> new directions { 2, 1/2 }        (positive)
    c = +-1/2 -> new directions { 3, 1/3 }        (positive)
    d = +-1   -> new directions { -2, -1/2 }      (negative)
    d = +-1/2 -> new directions { -3, -1/3 }      (negative)

The two sets are disjoint -- one positive, one negative -- so FOUR distinct new
directions are required, and only two slots exist.  Contradiction. QED
"""
from fractions import Fraction as F
from itertools import product

INF = None
KNOWN_DIRS = {INF, F(0), F(1), F(-1)}
SPECIAL = [F(1), F(1, 2), F(-1), F(-1, 2)]


def inv(x):
    return F(0) if x is INF else (INF if x == 0 else 1 / x)


def quad_c(c):
    t, s = (c - 1) / c, (1 + c) / c
    return [t, inv(t), s, inv(s)]


def quad_d(d):
    t, s = (1 - d) / d, -(1 + d) / d
    return [t, inv(t), s, inv(s)]


def check_step1(height=80):
    """Every parameter outside {1,1/2,-1,-1/2} yields four distinct new points."""
    vals = sorted({F(p, q) for p in range(-height, height + 1)
                   for q in range(1, height + 1) if p != 0})
    for quad, name in ((quad_c, "c"), (quad_d, "d")):
        for x in vals:
            q = quad(x)
            new = {y for y in q if y not in KNOWN_DIRS}
            if x in SPECIAL:
                assert len(new) == 2, (name, x, new)
            else:
                # outside the special set: four points, distinct, all new
                assert len(set(q)) == 4, ("not distinct", name, x, q)
                assert len(new) == 4, (name, x, new)
    return len(vals)


def check_step2():
    """All sixteen surviving (c,d) need four distinct new directions."""
    rows = []
    for c, d in product(SPECIAL, SPECIAL):
        nc = {y for y in quad_c(c) if y not in KNOWN_DIRS}
        nd = {y for y in quad_d(d) if y not in KNOWN_DIRS}
        rows.append((c, d, nc, nd, len(nc | nd)))
    return rows


def run(verbose=True):
    n = check_step1()
    if verbose:
        print("   step 1: %d parameter values; outside {1,1/2,-1,-1/2} the quadruple" % n)
        print("           is always four distinct NEW directions")
    worst = 99
    for c, d, nc, nd, tot in check_step2():
        worst = min(worst, tot)
    if verbose:
        print("   step 2: over the sixteen surviving (c,d), new directions required =",
              worst, "  (slots available: 2)")
    return worst
