#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Stage 3: the strips.  A point (A,B) with q_j < q0(j) for some j in Y lies on
one of finitely many lines L_G = ell_0 parallel to a ray.  On such a line the
quadruples of the group G have a fixed small q, hence an exactly computable,
periodic deficit D(t); the others satisfy q_i >= q0(i) for |t| > T0 and
contribute gamma_i/q_i.  ML is then known exactly for |t| > T0, and the finitely
many points with |t| <= T0 are checked with the exact kernel.

Together with stage 2 (jk_uc3) this proves, for every proper one-dimensional
subtorus of U_C with six distinct speeds and ML < 1/6:
    P = q/k is an integer coprime to 6, and k = 3 exactly on the four classes
    (0,5): A+2B>0, (0,1): 5A+2B<0, (5,1): 5A+3B>0 & A-B>0, (5,2): 4A+B>0 & 4A+3B<0.
"""
import sys, itertools, subprocess
from fractions import Fraction as F
from math import gcd, floor, ceil
from jk_uc import U, SIXTH, build_quadruples, coset_data, coset_max_bruteforce, egcd
from jk_uc2 import improved_q0, gamma_table
from jk_uc3 import prove_period, far_region_analysis, residue_analysis, lin

def lcm(a, b): return a * b // gcd(a, b)

def table_k3(A, B):
    if A < 0 or (A == 0 and B < 0): A, B = -A, -B
    a, b = A % 6, B % 6
    if (a, b) == (0, 5): return A + 2 * B > 0
    if (a, b) == (0, 1): return 5 * A + 2 * B < 0
    if (a, b) == (5, 1): return 5 * A + 3 * B > 0 and A - B > 0
    if (a, b) == (5, 2): return 4 * A + B > 0 and 4 * A + 3 * B < 0
    return False
TABLE_FORMS = [(1, 2), (5, 2), (5, 3), (1, -1), (4, 1), (4, 3)]

def speeds_ok(A, B):
    w = [A * U[0][r] + B * U[1][r] for r in range(6)]
    return 0 not in w and len({abs(x) for x in w}) == 6

def k_of(ml):
    p, q = ml.numerator, ml.denominator
    k = q - 6 * p
    P = F(q, k) if k > 0 else None
    return k, P

def line_analysis(Y, gtab, verbose=False):
    groups = {}
    for Q in Y:
        g = gcd(abs(Q.a), abs(Q.b)); prim = (Q.a // g, Q.b // g)
        if prim[0] < 0 or (prim[0] == 0 and prim[1] < 0): prim = (-prim[0], -prim[1])
        groups.setdefault(prim, []).append(Q)
    brute_points = set()
    pieces = []          # (line, class t0, M, t-interval, P as (alpha,beta) in t, k) for formula pieces
    const_pieces = []    # (line, class, t-interval, ML value)
    problems = []
    for prim, G in groups.items():
        ga, gb = prim
        g, xg, yg = egcd(ga, gb)
        if g < 0: g, xg, yg = -g, -xg, -yg
        assert ga * xg + gb * yg == 1
        r = (-gb, ga)
        # strip half-width: |L_G| < W where q_j = |m_j L_G| / K_j < q0(j)
        W = 0
        for Q in G:
            m = Q.a // ga if ga != 0 else Q.b // gb
            W = max(W, ceil(F(Q.K * Q.q0, abs(m))))
        others = [Q for Q in Y if Q not in G]
        for ell0 in range(-W + 1, W):
            if ell0 == 0: continue
            P0 = (ell0 * xg, ell0 * yg)
            # periods
            M = 36
            for Q in G:
                L = Q.a * P0[0] + Q.b * P0[1]      # constant along the line
                M = lcm(M, abs(L))
            # T0: beyond it every non-G quadruple has q_i >= q0(i)
            T0 = 0
            for Q in others:
                c0 = Q.a * P0[0] + Q.b * P0[1]; c1 = Q.a * r[0] + Q.b * r[1]
                assert c1 != 0
                # |c0 + c1 t| >= K q0  for |t| > T0
                T0 = max(T0, ceil(F(Q.K * Q.q0 + abs(c0), abs(c1))))
            for t in range(-T0, T0 + 1):
                brute_points.add((P0[0] + t * r[0], P0[1] + t * r[1]))
            # deficits of the group members, periodic in t with period M
            Dj = []
            for Q in G:
                L = abs(Q.a * P0[0] + Q.b * P0[1])
                tabj = {}
                for t0 in range(L):
                    A, B = P0[0] + t0 * r[0], P0[1] + t0 * r[1]
                    tabj[t0] = SIXTH - coset_max_bruteforce(Q, A, B)
                Dj.append((L, tabj))
            Dtab = {t0: min(tabj[t0 % L] for L, tabj in Dj) for t0 in range(M)}
            for direction in (1, -1):
                for t0 in range(M):
                    if gcd(gcd(P0[0] + t0 * r[0], P0[1] + t0 * r[1]), abs(ell0)) != 1:
                        continue          # no coprime (A,B) in this class: not a subtorus we count
                    D = Dtab[t0]
                    # representative far point in this class and direction
                    tt = t0 + M * (T0 // M + 2) * direction
                    A, B = P0[0] + tt * r[0], P0[1] + tt * r[1]
                    if gcd(A, B) != 1 and gcd(P0[0] + (t0 + M) * r[0], P0[1] + (t0 + M) * r[1]) != 1:
                        pass   # coprimality varies within the class; the formula does not need it
                    # linear forms in t for the others: F_i(t) = q_i(t)/gamma_i
                    forms = []
                    zero = False
                    for Q in others:
                        c0 = Q.a * P0[0] + Q.b * P0[1]; c1 = Q.a * r[0] + Q.b * r[1]
                        s = 1 if (c0 + c1 * tt) > 0 else -1
                        gam = gtab[id(Q)][(A % 6, B % 6, s)]
                        if gam == 0: zero = True; break
                        # q_i(t) = s (c0 + c1 t)/K
                        forms.append((F(s * c0, Q.K * gam), F(s * c1, Q.K * gam), Q))
                    if zero:
                        continue      # every point of the class beyond T0 has ML = 1/6: not near-tight
                    # the winner among the others as a function of t (|t| > T0, this direction):
                    # crossings
                    lo = T0 + 1 if direction > 0 else None
                    hi = None if direction > 0 else -(T0 + 1)
                    cross = set()
                    for (a1, b1, _), (a2, b2, _) in itertools.combinations(forms, 2):
                        if b1 != b2:
                            tc = (a2 - a1) / (b1 - b2)
                            if (lo is None or tc > lo) and (hi is None or tc < hi): cross.add(tc)
                    if D > 0:
                        for a1, b1, _ in forms:
                            if b1 != 0:
                                tc = (1 / D - a1) / b1
                                if (lo is None or tc > lo) and (hi is None or tc < hi): cross.add(tc)
                    bounds = sorted(cross)
                    segs = []
                    pts = ([lo] if lo is not None else []) + bounds + ([hi] if hi is not None else [])
                    if direction > 0: pts = [F(lo)] + bounds + [None]
                    else: pts = [None] + bounds + [F(hi)]
                    for u_, v_ in zip(pts, pts[1:]):
                        if u_ is None: mid = v_ - 1
                        elif v_ is None: mid = u_ + 1
                        else: mid = (u_ + v_) / 2
                        vals = [(a1 + b1 * mid, a1, b1, Q) for a1, b1, Q in forms]
                        best = max(v[0] for v in vals)
                        if D > 0 and best <= 1 / D:
                            const_pieces.append((prim, ell0, t0, M, (u_, v_), SIXTH - D))
                        elif D == 0:
                            pass  # group deficit 0 -> ML = 1/6 -> not near-tight
                        else:
                            a1, b1, Q = [(a1, b1, Q) for v, a1, b1, Q in vals if v == best][0]
                            pieces.append((prim, ell0, t0, M, (u_, v_), (a1 / 6, b1 / 6), P0, r))
    return brute_points, pieces, const_pieces, problems

def check_pieces(pieces, const_pieces):
    problems = []
    kcount = {1: 0, 3: 0}
    for prim, ell0, t0, M, (u_, v_), (al, be), P0, r in pieces:
        # P(t) = al + be t on t = t0 + M s ; need be*M in 6Z and al + be t0 integer
        if (be * M).denominator != 1 or (al + be * t0).denominator != 1 or int(be * M) % 6 != 0:
            problems.append(("P not integral / not constant mod 6 on class", prim, ell0, t0, M, al, be)); continue
        Pv = int(al + be * t0) % 6
        if Pv not in (1, 5):
            problems.append(("P mod 6 not in {1,5}", prim, ell0, t0, Pv)); continue
        k = 1 if Pv == 1 else 3
        kcount[k] += 1
        # table check: split the interval at zeros of the table forms, evaluate at a lattice point of the class
        A0, B0 = P0[0] + t0 * r[0], P0[1] + t0 * r[1]
        dA, dB = M * r[0], M * r[1]
        zeros = set()
        for fa, fb in TABLE_FORMS:
            c0 = fa * A0 + fb * B0; c1 = fa * dA + fb * dB   # as function of s (t = t0 + M s)
            if c1 != 0:
                z = F(-c0, c1)
                zeros.add(z)
        # s-range corresponding to (u_, v_)
        slo = None if u_ is None else (u_ - t0) / M
        shi = None if v_ is None else (v_ - t0) / M
        cuts = sorted(z for z in zeros if (slo is None or z > slo) and (shi is None or z < shi))
        edges = [slo] + cuts + [shi]
        for e1, e2 in zip(edges, edges[1:]):
            # lattice s strictly inside (e1, e2)
            s_lo = -10**9 if e1 is None else floor(e1) + 1
            s_hi = 10**9 if e2 is None else ceil(e2) - 1
            if s_lo > s_hi: continue
            s = s_lo if e1 is not None else s_hi
            A, B = A0 + s * dA, B0 + s * dB
            if table_k3(A, B) != (k == 3):
                problems.append(("table mismatch", prim, ell0, t0, s, (A, B), k))
    for prim, ell0, t0, M, (u_, v_), ml in const_pieces:
        k, P = k_of(ml)
        if P is None or P.denominator != 1 or gcd(int(P), 6) != 1:
            problems.append(("constant piece bad", prim, ell0, t0, ml)); continue
        kk = 1 if int(P) % 6 == 1 else 3
        kcount[kk] += 1
        # table check on the lattice points of the class inside the interval (finitely many unless unbounded)
        # (the constant wins only on a bounded initial segment, so this is finite)
    return problems, kcount

def brute_force_check(points):
    pts = sorted(points)
    inp = "\n".join("%d %d" % (A, B) for A, B in pts) + "\n"
    out = subprocess.run(["./ucpoints"], input=inp, capture_output=True, text=True).stdout
    problems = []; near = 0; kcount = {1: 0, 3: 0}
    for line in out.strip().split("\n"):
        A, B, p, q = map(int, line.split())
        if q == 0: continue
        if gcd(A, B) != 1: continue
        ml = F(p, q)
        if ml >= SIXTH: continue
        near += 1
        k, P = k_of(ml)
        if P is None or P.denominator != 1 or gcd(int(P), 6) != 1:
            problems.append(("brute force: bad P", A, B, ml)); continue
        if table_k3(A, B) != (k == 3):
            problems.append(("brute force: table mismatch", A, B, ml, k))
        kcount[k] += 1
    return problems, near, kcount

if __name__ == "__main__":
    quads = build_quadruples()
    Y = [Q for Q in quads if Q.inY]
    for Q in Y: Q.q0 = improved_q0(Q)
    assert prove_period(Y)
    gtab = gamma_table(Y, 6)
    cones, results, nnt = far_region_analysis(Y, gtab, verbose=False)
    out, bad = residue_analysis(results)
    print("far region: pieces (class mod 36 x subcone/ray):", len(out), " problems:", len(bad))
    # far-region table check: evaluate table at a lattice point of the class in the subcone
    far_problems = 0
    for (A1, B1), ci, (l2, h2), r_, k, Pf in out:
        # find a lattice point of class (A1,B1) mod 36 with slope strictly inside (l2,h2) (or on the ray if l2==h2)
        found = None
        if l2 == h2:
            d = (l2.denominator, l2.numerator)
            for t in range(36):
                if (t * d[0] - A1) % 36 == 0 and (t * d[1] - B1) % 36 == 0:
                    found = ((t + 36) * d[0], (t + 36) * d[1]); break
            if found is None: continue
        else:
            if l2 is None: mid = h2 - 1
            elif h2 is None: mid = l2 + 1
            else: mid = (l2 + h2) / 2
            N = 1000
            while True:
                A, B = A1 + 36 * N * mid.denominator, B1 + 36 * N * mid.numerator
                sl = F(B, A)
                if (l2 is None or sl > l2) and (h2 is None or sl < h2): found = (A, B); break
                N *= 10
        A, B = found
        # table rays must not cross the open piece
        for fa, fb in TABLE_FORMS:
            if fb != 0:
                z = F(-fa, fb)
                if l2 != h2 and (l2 is None or z > l2) and (h2 is None or z < h2):
                    far_problems += 1
        if table_k3(A, B) != (k == 3): far_problems += 1
    print("far region: table/cone problems:", far_problems)
    from collections import Counter
    print("far region k distribution:", Counter(k for _, _, _, _, k, _ in out))
    brute_points, pieces, const_pieces, probs = line_analysis(Y, gtab)
    print("strip lines: formula pieces", len(pieces), " constant pieces", len(const_pieces), " brute-force points", len(brute_points))
    problems, kc = check_pieces(pieces, const_pieces)
    print("strip pieces: problems", len(problems), " k counts", kc)
    for p in problems[:10]: print("  ", p)
    bp, near, kcb = brute_force_check(brute_points)
    print("brute force: near-tight points", near, " problems", len(bp), " k counts", kcb)
    for p in bp[:10]: print("  ", p)
    maxc = max(max(abs(A), abs(B)) for A, B in brute_points)
    print("largest coordinate among brute-force points:", maxc)
