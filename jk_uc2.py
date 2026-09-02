#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Exact, explicit Jain-Kravitz analysis of the relative spectrum of U_C.

Stage 1 (this file): quadruples, thresholds q0, the offset gamma as a function
of the residue class of (A,B) modulo M and the sign of L, and validation of
   ML(T) = 1/6 - min_{j in Y} gamma_j / q_j
against the exact kernel on every (A,B) with |A|,|B| <= H for which all
q_j >= q0(j) ("far region").
"""
import sys, itertools, random
from fractions import Fraction as F
from math import gcd, floor, ceil
from jk_uc import (U, SIXTH, build_quadruples, coset_data, coset_max_bruteforce,
                   f_eval, egcd)

def improved_q0(Q):
    lam_s = min(min(z[1], z[3]) for z in Q.zones)
    q1 = ceil(1 / Q.rho_min)
    q2 = floor(lam_s / (SIXTH - Q.f_out)) + 1
    return max(q1, q2)

def gamma_of(Q, A, B):
    """gamma such that deficit = gamma/q  (formula, no validity check)."""
    q, th0 = coset_data(Q, A, B)
    best = None
    for tau, lm, rm, lp, rp in Q.zones:
        am = (tau - th0) % F(1, q); ap = (th0 - tau) % F(1, q)
        for val in (lm * am, lp * ap):
            if best is None or val < best: best = val
    return best * q

def find_period(Y, Mcands=(6, 12, 18, 36, 72, 108, 216), samples=4000, seed=3):
    """smallest M in the list such that, for every j in Y, gamma_j depends only on
    ((A,B) mod M, sign L_j).  Verified by sampling; the proof is Lemma 2.4 of JK."""
    rng = random.Random(seed)
    pts = []
    while len(pts) < samples:
        A = rng.randint(0, 3000); B = rng.randint(-3000, 3000)
        if gcd(A, B) == 1: pts.append((A, B))
    for M in Mcands:
        ok = True
        for Q in Y:
            table = {}
            for A, B in pts:
                L = Q.a * A + Q.b * B
                if L == 0: continue
                key = (A % M, B % M, L > 0)
                g = gamma_of(Q, A, B)
                if key in table and table[key] != g: ok = False; break
                table[key] = g
            if not ok: break
        if ok: return M
    return None

def gamma_table(Y, M):
    """gamma_j[(A0,B0,sign)] for every class and sign, computed at a representative."""
    tab = {}
    for Q in Y:
        d = {}
        for A0 in range(M):
            for B0 in range(M):
                for sgn in (1, -1):
                    # find a representative (A,B) in the class with sign(L) = sgn and L != 0
                    rep = None
                    for t in range(1, 60):
                        for (A, B) in ((A0 + M * t, B0), (A0, B0 + M * t), (A0 + M * t, B0 + M * t),
                                       (A0 + M * t, B0 - M * t), (A0, B0 - M * t), (A0 - M * t, B0)):
                            L = Q.a * A + Q.b * B
                            if L != 0 and (L > 0) == (sgn > 0): rep = (A, B); break
                        if rep: break
                    if rep is None: continue
                    d[(A0, B0, sgn)] = gamma_of(Q, rep[0], rep[1])
        tab[id(Q)] = d
    return tab

def predicted_ml(Y, M, gtab, A, B):
    """1/6 - min_j gamma_j/q_j, or None if some q_j < q0(j) or L_j = 0."""
    best = None
    for Q in Y:
        L = Q.a * A + Q.b * B
        if L == 0: return None
        q = abs(L) // Q.K
        if q < Q.q0: return None
        g = gtab[id(Q)][(A % M, B % M, 1 if L > 0 else -1)]
        d = g / q
        if best is None or d < best: best = d
    return SIXTH - best

if __name__ == "__main__":
    H = int(sys.argv[1]) if len(sys.argv) > 1 else 300
    quads = build_quadruples()
    Y = [Q for Q in quads if Q.inY]
    for Q in Y: Q.q0 = improved_q0(Q)
    print("q0 per quadruple:", [(Q.i + 1, Q.j + 1, Q.eps, Q.ell, Q.q0) for Q in Y])
    # re-validate the formula with the improved thresholds
    rng = random.Random(5); bad = 0; tested = 0
    while tested < 20000:
        A = rng.randint(0, 400); B = rng.randint(-400, 400)
        if gcd(A, B) != 1: continue
        for Q in Y:
            r = coset_data(Q, A, B)
            if r is None or r[0] < Q.q0: continue
            tested += 1
            if SIXTH - gamma_of(Q, A, B) / r[0] != coset_max_bruteforce(Q, A, B):
                bad += 1; print("MISMATCH", (Q.i, Q.j, Q.eps, Q.ell), A, B)
    print("coset formula with improved q0: tested", tested, "mismatches", bad)
    M = find_period(Y)
    print("period M =", M)
    gtab = gamma_table(Y, M)
    # validation against the exact dump
    exact = {}
    with open("ucdump%d.txt" % H) as fh:
        for line in fh:
            A, B, p, q, Qd = map(int, line.split()); exact[(A, B)] = F(p, q)
    tested = bad = far = 0
    for (A, B), ml in exact.items():
        pred = predicted_ml(Y, M, gtab, A, B)
        if pred is None: continue
        far += 1
        # non-Y quadruples have coset max 0, so the prediction is exact whenever ML < 1/6,
        # and when the prediction is 1/6 the true ML is 1/6 (both directions: ML >= coset max)
        if pred != ml:
            bad += 1
            if bad < 10: print("ML MISMATCH", A, B, pred, ml)
    print("far-region points up to H=%d: %d, mismatches %d (of %d points)" % (H, far, bad, len(exact)))
