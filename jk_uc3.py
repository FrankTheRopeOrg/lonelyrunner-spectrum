#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Stage 2: the far region.  For every residue class of (A,B) modulo 6 and every
cone cut by the nine lines L_G = 0, subdivide by the rays where two of the
homogeneous forms F_j = s_j L_j / (K_j gamma_j) agree; on each subcone a single
F_j is maximal, ML = 1/6 - 1/F_j, P = F_j/6.  Then read off P mod 6 on classes
modulo 36 and hence k.  Also proves the period M = 6 exhaustively.
"""
import sys, itertools
from fractions import Fraction as F
from math import gcd, floor, ceil
from jk_uc import U, SIXTH, build_quadruples, coset_data
from jk_uc2 import improved_q0, gamma_of, gamma_table

def prove_period(Y):
    """gamma_j depends on (q mod x_h, c mod K, sign L); x_h | 6 and K | 6, so on
    ((A,B) mod 6, sign).  Exhaustive check: over classes mod 36 the value only
    depends on the class mod 6."""
    tab36 = gamma_table(Y, 36)
    for Q in Y:
        d = tab36[id(Q)]
        seen = {}
        for (A0, B0, s), g in d.items():
            key = (A0 % 6, B0 % 6, s)
            if key in seen and seen[key] != g: return False
            seen[key] = g
    return True

def lin(form, d):   # form = (alpha, beta) ; d = (dA, dB)
    return form[0] * d[0] + form[1] * d[1]

def far_region_analysis(Y, gtab, verbose=True):
    # the nine primitive forms L_G and their ray slopes in the half plane A > 0
    forms = {}
    for Q in Y:
        g = gcd(abs(Q.a), abs(Q.b)); prim = (Q.a // g, Q.b // g)
        if prim[0] < 0 or (prim[0] == 0 and prim[1] < 0): prim = (-prim[0], -prim[1])
        forms.setdefault(prim, []).append(Q)
    slopes = set()
    for (a, b) in forms:
        if b == 0: continue                    # the line A = 0: boundary of the half plane
        # a A + b B = 0  -> slope B/A = -a/b  (b != 0)
        slopes.add(F(-a, b))
    slopes = sorted(slopes)
    if verbose: print("ray slopes in A>0:", [str(s) for s in slopes])
    bounds = [None] + slopes + [None]
    cones = []
    for lo, hi in zip(bounds, bounds[1:]):
        if lo is None: mid = hi - 1
        elif hi is None: mid = lo + 1
        else: mid = (lo + hi) / 2
        cones.append((lo, hi, (mid.denominator, mid.numerator)))   # test direction (A,B)
    results = []     # (class36, cone index, slope interval, winner, P form, P mod 6)
    non_near_tight = []
    for ci, (lo, hi, dtest) in enumerate(cones):
        for A0 in range(6):
            for B0 in range(6):
                if gcd(gcd(A0, B0), 6) != 1: continue
                # signs and gammas
                Fj = []
                zero = False
                for Q in Y:
                    s = 1 if lin((Q.a, Q.b), dtest) > 0 else -1
                    g = gtab[id(Q)][(A0, B0, s)]
                    if g == 0: zero = True; break
                    Fj.append(((F(s * Q.a) / (Q.K * g), F(s * Q.b) / (Q.K * g)), Q))
                if zero:
                    non_near_tight.append(((A0, B0), ci)); continue
                # tie rays inside the cone
                rays = set()
                for (f1, _), (f2, _) in itertools.combinations(Fj, 2):
                    a, b = f1[0] - f2[0], f1[1] - f2[1]
                    if a == 0 and b == 0: continue
                    if b == 0: continue          # ray is A = 0, not in open half plane
                    s = F(-a, b)
                    if (lo is None or s > lo) and (hi is None or s < hi): rays.add(s)
                rays = sorted(rays)
                sub = [lo] + rays + [hi]
                for k in range(len(sub) - 1):
                    l2, h2 = sub[k], sub[k + 1]
                    if l2 is None: mid = h2 - 1
                    elif h2 is None: mid = l2 + 1
                    else: mid = (l2 + h2) / 2
                    d = (mid.denominator, mid.numerator)
                    vals = [(lin(f, d), f, Q) for f, Q in Fj]
                    best = max(v[0] for v in vals)
                    winners = [(f, Q) for v, f, Q in vals if v == best]
                    f, Q = winners[0]
                    Pform = (f[0] / 6, f[1] / 6)
                    results.append(((A0, B0), ci, (l2, h2), Q, Pform, [w[1] for w in winners]))
                # points on the tie rays themselves: same value from both forms; include as
                # degenerate "subcones" (they are handled by evaluating either form)
                for s in rays:
                    d = (s.denominator, s.numerator)
                    vals = [(lin(f, d), f, Q) for f, Q in Fj]
                    best = max(v[0] for v in vals)
                    f, Q = [(f, Q) for v, f, Q in vals if v == best][0]
                    results.append(((A0, B0), ci, (s, s), Q, (f[0] / 6, f[1] / 6), None))
    return cones, results, non_near_tight

def residue_analysis(results):
    """P = alpha A + beta B on lattice points of the class; check 6alpha,6beta integers,
    integrality, and P mod 6 on classes modulo 36."""
    out = []      # (class36, cone, interval, Pmod6, k, Pform)
    bad = []
    for (A0, B0), ci, (l2, h2), Q, (al, be), ties in results:
        if l2 != h2 and ((6 * al).denominator != 1 or (6 * be).denominator != 1):
            bad.append(("non-sixth-integral coefficients", (A0, B0), ci, l2, h2, al, be)); continue
        if l2 == h2:
            d = (l2.denominator, l2.numerator)
            cc = al * d[0] + be * d[1]
            if (6 * cc).denominator != 1:
                bad.append(("ray coefficient not sixth-integral", (A0, B0), ci, l2, al, be)); continue
        for A1 in range(A0, 36, 6):
            for B1 in range(B0, 36, 6):
                if l2 == h2:
                    # tie ray: lattice points t*(dA,dB); need some t with t*d = (A1,B1) mod 36
                    d = (l2.denominator, l2.numerator)
                    if not any(((t * d[0] - A1) % 36 == 0 and (t * d[1] - B1) % 36 == 0) for t in range(36)):
                        continue
                    # P(t d) = t (al dA + be dB): find t mod 36 with t d = (A1,B1) mod 36, any such t
                    t = next(t for t in range(36) if (t * d[0] - A1) % 36 == 0 and (t * d[1] - B1) % 36 == 0)
                    Pv = t * (al * d[0] + be * d[1])
                    # P(t+36 s) - P(t) = 36 s (al dA + be dB) in 6Z since 6al,6be integers
                else:
                    Pv = al * A1 + be * B1
                if Pv.denominator != 1:
                    bad.append(("P not integral", (A1, B1), ci, l2, h2, al, be)); continue
                r = int(Pv) % 6
                k = 6 // gcd(r - 1, 6) if r != 1 else 1
                out.append(((A1, B1), ci, (l2, h2), r, k, (al, be)))
    return out, bad

if __name__ == "__main__":
    quads = build_quadruples()
    Y = [Q for Q in quads if Q.inY]
    for Q in Y: Q.q0 = improved_q0(Q)
    print("period 6 proven exhaustively:", prove_period(Y))
    gtab = gamma_table(Y, 6)
    cones, results, nnt = far_region_analysis(Y, gtab)
    print("cones:", len(cones), " (class, cone) pairs that are never near-tight:", len(nnt))
    print("subcones/rays with a winner:", len(results))
    out, bad = residue_analysis(results)
    print("problems:", len(bad))
    for b in bad[:10]: print("  ", b)
    from collections import Counter
    print("P mod 6 distribution over (class mod 36, piece):", Counter(r for _, _, _, r, _, _ in out))
    print("k distribution:", Counter(k for _, _, _, _, k, _ in out))
    # aggregate the k=3 pieces by class mod 6 and slope interval
    k3 = {}
    for (A1, B1), ci, (l2, h2), r, k, Pf in out:
        if k == 3:
            k3.setdefault((A1 % 6, B1 % 6), set()).add((str(l2), str(h2)))
    for c, ivs in sorted(k3.items()):
        print("  k = 3 on class", c, ":", sorted(ivs, key=lambda x: (float(F(x[0])) if x[0] != 'None' else -1e9)))
    # forms P for the k=3 pieces
    k3forms = {}
    for (A1, B1), ci, (l2, h2), r, k, Pf in out:
        k3forms.setdefault(((A1 % 6, B1 % 6), k), set()).add((str(Pf[0]), str(Pf[1])))
    for key, fs in sorted(k3forms.items()): print("  class", key[0], "k =", key[1], "P forms:", sorted(fs))
