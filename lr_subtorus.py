# -*- coding: utf-8 -*-
"""Two-dimensional subtori of (R/Z)^n: coordinate forms, repeat directions,
and an exact decision procedure for D(U) = 1/3 when n = 6.

A 2-dimensional subtorus is U = <u, v>_R with u, v in Z^n independent.  Its
points are (s u + t v) mod 1.  The r-th coordinate function on the parameter
plane is the linear form  phi_r(s,t) = s*u_r + t*v_r,  identified with the
integer pair (u_r, v_r).

  * U is proper  iff  no phi_r vanishes identically.
  * The integer vectors of U with a vanishing r-th entry are those on the
    line ker phi_r  ("vanishing line").
  * The integer vectors with |w_p| = |w_q| are those on one of the at most
    thirty lines ker(phi_p -+ phi_q)  ("repeat lines").
  * A repeat line that is not a vanishing line is a FREE repeat direction;
    the primitive vector on it is a zero-free vector with a repeat.

D(U) = 1/2 - ML(U) with ML(U) = max_{s,t} min_r ||s u_r + t v_r||.
"""
from fractions import Fraction as F
from math import gcd, floor, ceil
from itertools import combinations

__all__ = ["forms", "primitive", "vanishing_directions", "repeat_directions",
           "free_repeat_directions", "good_pairs", "projective_classes",
           "exceeds_sixth", "vector_on_direction"]


# ---------------------------------------------------------------- geometry

def forms(u, v):
    return [(u[r], v[r]) for r in range(len(u))]


def primitive(f):
    """Canonical representative of the line ker(f), f = (alpha,beta) != 0."""
    a, b = f
    if a == 0 and b == 0:
        return None
    g = gcd(abs(a), abs(b))
    a, b = a // g, b // g
    if a < 0 or (a == 0 and b < 0):
        a, b = -a, -b
    return (a, b)


def vanishing_directions(u, v):
    return {primitive(f) for f in forms(u, v)}


def repeat_directions(u, v):
    """{direction: [(p,q,eps), ...]} over all pairs and both signs.
    A pair with phi_p == eps*phi_q identically is degenerate (U lies inside
    {x_p = eps x_q}) and is reported separately."""
    phi = forms(u, v)
    n = len(phi)
    out, degenerate = {}, []
    for p, q in combinations(range(n), 2):
        for eps in (1, -1):
            f = (phi[p][0] - eps * phi[q][0], phi[p][1] - eps * phi[q][1])
            if f == (0, 0):
                degenerate.append((p, q, eps))
            else:
                out.setdefault(primitive(f), []).append((p, q, eps))
    return out, degenerate


def free_repeat_directions(u, v):
    R, _ = repeat_directions(u, v)
    Z = vanishing_directions(u, v)
    return {d: w for d, w in R.items() if d not in Z}


def good_pairs(u, v):
    """Pairs (p,q) for which BOTH generators of the Jain-Kravitz normal form
    are free of vanishing entries.  This is their criterion; it is strictly
    stronger than having two free repeat directions."""
    phi = forms(u, v)
    Z = vanishing_directions(u, v)
    out = []
    for p, q in combinations(range(len(phi)), 2):
        fp, fq = phi[p], phi[q]
        if fp[0] * fq[1] - fp[1] * fq[0] == 0:
            continue
        dm = primitive((fp[0] - fq[0], fp[1] - fq[1]))
        dp = primitive((fp[0] + fq[0], fp[1] + fq[1]))
        if dm not in Z and dp not in Z:
            out.append((p, q))
    return out


def projective_classes(u, v):
    """Partition of {0..n-1} by proportionality of the coordinate forms."""
    phi = forms(u, v)
    cls = {}
    for r, f in enumerate(phi):
        cls.setdefault(primitive(f), []).append(r)
    return cls


def vector_on_direction(u, v, d):
    """Primitive integer vector of U spanning the line ker(d)."""
    a, b = d
    s, t = -b, a
    w = [s * u[r] + t * v[r] for r in range(len(u))]
    g = 0
    for x in w:
        g = gcd(g, abs(x))
    return tuple(x // g for x in w) if g else tuple(w)


# ------------------------------------------------- exact decision D(U)=1/3

def _t_intervals(s, ui, vi, lo, hi):
    """Open t-intervals in [0,1) with frac(s*ui + t*vi) in (lo,hi); vi != 0."""
    out = []
    a, b = (s * ui, s * ui + vi) if vi > 0 else (s * ui + vi, s * ui)
    for m in range(floor(a - 1), ceil(b + 1) + 1):
        x0 = (m + lo - s * ui) / vi
        x1 = (m + hi - s * ui) / vi
        p, q = (x0, x1) if x0 < x1 else (x1, x0)
        p, q = max(p, F(0)), min(q, F(1))
        if p < q:
            out.append((p, q))
    out.sort()
    return out


def _intersect(A, B):
    out, i, j = [], 0, 0
    while i < len(A) and j < len(B):
        lo, hi = max(A[i][0], B[j][0]), min(A[i][1], B[j][1])
        if lo < hi:
            out.append((lo, hi))
        if A[i][1] < B[j][1]:
            i += 1
        else:
            j += 1
    return out


def _feasible_at(s, u, v, lo, hi):
    cur = [(F(0), F(1))]
    for ui, vi in zip(u, v):
        if vi == 0:
            if not (lo < (s * ui) % 1 < hi):
                return False
            continue
        cur = _intersect(cur, _t_intervals(s, ui, vi, lo, hi))
        if not cur:
            return False
    return True


def _critical_s(u, v):
    """All s where the ordering of the t-interval endpoints can change.
    Endpoints satisfy 6(s u_i + t v_i) = c with c = +-1 mod 6, so two of them
    coincide at s = (c' v_i - c v_j) / (6 (u_j v_i - u_i v_j))."""
    n = len(u)
    S = {F(0)}
    for i in range(n):
        if v[i] == 0 and u[i] != 0:
            for c in range(-6 * abs(u[i]) - 6, 6 * abs(u[i]) + 7):
                if c % 6 in (1, 5):
                    s = F(c, 6 * u[i])
                    if 0 <= s < 1:
                        S.add(s)
    for i, j in combinations(range(n), 2):
        det = u[j] * v[i] - u[i] * v[j]
        if det == 0:
            continue
        Ri = 6 * (abs(u[i]) + abs(v[i])) + 12
        Rj = 6 * (abs(u[j]) + abs(v[j])) + 12
        ci = [c for c in range(-Ri, Ri + 1) if c % 6 in (1, 5)]
        cj = [c for c in range(-Rj, Rj + 1) if c % 6 in (1, 5)]
        for c in ci:
            for cp in cj:
                s = F(cp * v[i] - c * v[j], 6 * det)
                if 0 <= s < 1:
                    S.add(s)
    return sorted(S)


def exceeds_sixth(u, v):
    """True iff ML(U) > 1/6, i.e. iff D(U) < 1/3.  Exact.

    Sweep on s: for fixed s the admissible t form an intersection of finite
    unions of open intervals with rational endpoints, and non-emptiness is
    constant between consecutive critical values of s, so testing midpoints
    decides the question."""
    lo, hi = F(1, 6), F(5, 6)
    S = _critical_s(u, v) + [F(1)]
    for a, b in zip(S, S[1:]):
        if a != b and _feasible_at((a + b) / 2, u, v, lo, hi):
            return True
    return False
