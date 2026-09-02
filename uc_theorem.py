#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Runs the complete verification for the subtorus U_C and prints the numbers
quoted in the paper:

  (1) the 18 quadruples of Y with L, K, the maximisers tau and the threshold q0;
  (2) the period 6 of the offsets gamma;
  (3) the far region: cones, pieces, P mod 6 always in {1,5}, k = 3 table;
  (4) the strips: lines, pieces, brute-force points, all checks;
  (5) the supremum of max_i |v_i| / q over the near-tight subtori of U_C.

Usage: python3 uc_theorem.py      (needs ./ucpoints compiled from ucpoints.c)
"""
import sys, itertools, subprocess, time
from fractions import Fraction as F
from math import gcd, floor, ceil
from collections import Counter
from jk_uc import U, SIXTH, build_quadruples
from jk_uc2 import improved_q0, gamma_table
from jk_uc3 import prove_period, far_region_analysis, residue_analysis
from jk_uc4 import line_analysis, check_pieces, brute_force_check, table_k3, TABLE_FORMS, k_of

def speeds_form(d):      # the six linear forms w_i evaluated at direction/point d
    return [d[0] * U[0][r] + d[1] * U[1][r] for r in range(6)]

def ratio_sup_far(out):
    """sup of max|w_i| / q over each far piece; q = P (k=1) or 3P (k=3)."""
    best = F(0); where = None; attained = []
    for (A1, B1), ci, (l2, h2), r_, k, (al, be) in out:
        c = 1 if k == 1 else 3
        if l2 == h2:
            d = (l2.denominator, l2.numerator)
            P = al * d[0] + be * d[1]
            val = max(abs(x) for x in speeds_form(d)) / (c * P)
            if val > best: best, where = val, ((A1, B1), l2, h2)
            if val == 3: attained.append(((A1, B1), l2))
            continue
        # subdivide at the zeros of the speed forms: slopes 0, -1, -2, -3, -3/2
        cuts = sorted({F(0), F(-1), F(-2), F(-3), F(-3, 2)})
        cuts = [s for s in cuts if (l2 is None or s > l2) and (h2 is None or s < h2)]
        edges = [l2] + cuts + [h2]
        for e in edges:
            if e is None: continue
            d = (e.denominator, e.numerator)
            P = al * d[0] + be * d[1]
            if P <= 0: continue
            val = max(abs(x) for x in speeds_form(d)) / (c * P)
            if val > best: best, where = val, ((A1, B1), l2, h2, e)
        # limits at vertical directions
        if l2 is None or h2 is None:
            d = (0, 1) if h2 is None else (0, -1)
            P = al * d[0] + be * d[1]
            if P > 0:
                val = max(abs(x) for x in speeds_form(d)) / (c * P)
                if val > best: best, where = val, ((A1, B1), l2, h2, 'inf')
    return best, where, attained

def ratio_sup_lines(pieces):
    best = F(0); where = None
    for prim, ell0, t0, M, (u_, v_), (al, be), P0, r in pieces:
        Pv0 = al + be * t0
        k = 1 if int(Pv0) % 6 == 1 else 3
        c = 1 if k == 1 else 3
        # w_i(t) = w_i(P0) + t w_i(r);  P(t) = al + be t
        wP = speeds_form(P0); wr = speeds_form(r)
        cand = set()
        if u_ is not None: cand.add(u_)
        if v_ is not None: cand.add(v_)
        for i in range(6):
            if wr[i] != 0:
                z = F(-wP[i], wr[i])
                if (u_ is None or z > u_) and (v_ is None or z < v_): cand.add(z)
        for t in cand:
            P = al + be * t
            if P <= 0: continue
            val = max(abs(wP[i] + wr[i] * t) for i in range(6)) / (c * P)
            if val > best: best, where = val, (prim, ell0, t0, t)
        if u_ is None or v_ is None:
            if be != 0:
                val = max(abs(wr[i]) for i in range(6)) / (c * abs(be))
                if val > best: best, where = val, (prim, ell0, t0, 'inf')
    return best, where

if __name__ == "__main__":
    t_start = time.time()
    quads = build_quadruples()
    Y = [Q for Q in quads if Q.inY]
    for Q in Y: Q.q0 = improved_q0(Q)
    print("(1) components of the subgroups U ∩ {x_i = eps x_j}: %d, of which %d with maximum 1/6 (the set Y);"
          % (len(quads), len(Y)))
    print("    the other %d components have maximum 0 (D = 1/2)." % (len(quads) - len(Y)))
    print("    (i,j,eps,ell)   L = w_i - eps w_j       K   maximisers tau          q0")
    for Q in Y:
        print("    (%d,%d,%s,%d)      %3d A %+3d B          %d   %-22s  %d" % (Q.i + 1, Q.j + 1, '+' if Q.eps > 0 else '-', Q.ell, Q.a, Q.b, Q.K, ','.join(str(t) for t in Q.taus), Q.q0))
    print("(2) the offsets gamma depend only on (A,B) mod 6 and the sign of L:", prove_period(Y))
    gtab = gamma_table(Y, 6)
    cones, results, nnt = far_region_analysis(Y, gtab, verbose=True)
    out, bad = residue_analysis(results)
    print("(3) far region: %d cones in the half-plane A > 0; %d (class mod 6, cone) pairs on which every subtorus has ML = 1/6;"
          % (len(cones), len(nnt)))
    print("    %d subcones and tie rays with a winner; %d pieces (class mod 36 x subcone/ray); problems: %d"
          % (len(results), len(out), len(bad)))
    print("    P mod 6 distribution:", dict(Counter(r for _, _, _, r, _, _ in out)), " k:", dict(Counter(k for _, _, _, _, k, _ in out)))
    # winners
    print("    winning forms P by class mod 6 and k:")
    forms = {}
    for (A1, B1), ci, (l2, h2), r_, k, (al, be) in out:
        forms.setdefault(((A1 % 6, B1 % 6), k), set()).add((al, be))
    for key in sorted(forms):
        print("      class %s, k = %d: %s" % (key[0], key[1], ", ".join("%sA%+sB" % (a, b) for a, b in sorted(forms[key]))))
    brute_points, pieces, const_pieces, probs = line_analysis(Y, gtab)
    nlines = len({(p, e) for p, e, *_ in pieces})
    print("(4) strips: %d lines; %d formula pieces, %d pieces where a strip quadruple wins; %d points checked directly"
          % (nlines, len(pieces), len(const_pieces), len(brute_points)))
    problems, kc = check_pieces(pieces, const_pieces)
    print("    piece checks: problems %d; k counts %s" % (len(problems), kc))
    bp, near, kcb = brute_force_check(brute_points)
    maxc = max(max(abs(A), abs(B)) for A, B in brute_points)
    print("    direct check: %d near-tight among them, problems %d, k counts %s, largest coordinate %d" % (near, len(bp), kcb, maxc))
    bf, wf, att = ratio_sup_far(out)
    bl, wl = ratio_sup_lines(pieces)
    print("(5) sup of max|v_i|/q: far pieces %s (at %s), attained on rays: %d; strip pieces %s (at %s)" % (bf, wf, len(att), bl, wl))
    # brute force points ratio
    pts = sorted(brute_points)
    inp = "\n".join("%d %d" % (A, B) for A, B in pts) + "\n"
    outp = subprocess.run(["./ucpoints"], input=inp, capture_output=True, text=True).stdout
    bb = F(0); wb = None
    for line in outp.strip().split("\n"):
        A, B, p, q = map(int, line.split())
        if q == 0 or gcd(A, B) != 1 or F(p, q) >= SIXTH: continue
        val = F(max(abs(x) for x in speeds_form((A, B))), q)
        if val > bb: bb, wb = val, (A, B)
    print("    directly checked points: max ratio %s at %s" % (bb, wb))
    print("total time %.0f s" % (time.time() - t_start))
