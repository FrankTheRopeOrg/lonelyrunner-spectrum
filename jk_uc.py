#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""The Jain-Kravitz apparatus on U_C = <u, v>, u = (1,0,1,2,3,3), v = (0,1,1,1,1,2),
implemented from scratch in exact rational arithmetic, with explicit thresholds.

For T = <A u + B v>, w = A u + B v, the maximum loneliness is attained at a time
where two runners coincide or are antipodal (Kravitz, Prop. 2.1), i.e.

    ML(w) = max over (i<j, eps) of  max_{m} min_r || m w_r / L ||,   L = w_i - eps w_j.

The points m w / L lie in the 1-dimensional subgroup U_{ij eps} = U ∩ {x_i = eps x_j},
which has K = gcd(u_i - eps u_j, v_i - eps v_j) components (circles).  On the
component ell the points of T form a coset theta_0 + <1/q>, q = |L|/K, and the
restriction f_ell of min_r ||x_r|| to the circle is piecewise linear.  Its maximum
is <= 1/6; the components with maximum exactly 1/6 form the set Y (18 quadruples).

For a quadruple in Y with maximisers tau_h, slopes lambda_h^-, lambda_h^+ and
q >= q_0 (explicit), the coset maximum is  1/6 - gamma/q  with gamma depending
only on the residue class of (A,B) modulo M and the sign of L (Lemma 2.4-2.5 of
Jain-Kravitz).  Everything below is verified against brute force.
"""
import sys, itertools
from fractions import Fraction as F
from math import gcd, floor, ceil

U = ([1, 0, 1, 2, 3, 3], [0, 1, 1, 1, 1, 2])
SIXTH = F(1, 6)

def egcd(a, b):
    if b == 0: return (a, 1, 0)
    g, x, y = egcd(b, a % b)
    return (g, y, x - (a // b) * y)

def nrm(x):
    x = x - floor(x)
    return min(x, 1 - x)

# ------------------------------------------------------------------ pieces
def pieces_of_min_tents(c, d):
    """f(theta) = min_r ||c_r + theta d_r|| on [0,1): list of (a, b, slope, icpt)
    with f = slope*theta + icpt on [a,b].  Exact."""
    kinks = {F(0), F(1)}
    for cr, dr in zip(c, d):
        if dr == 0: continue
        # c + theta d in (1/2) Z  ->  theta = (n/2 - c)/d
        lo = floor(2 * min(cr, cr + dr)) - 1; hi = ceil(2 * max(cr, cr + dr)) + 1
        for n in range(lo, hi + 1):
            th = (F(n, 2) - cr) / dr
            if 0 < th < 1: kinks.add(th)
    kinks = sorted(kinks)
    out = []
    for a, b in zip(kinks, kinks[1:]):
        mid = (a + b) / 2
        lines = []
        for cr, dr in zip(c, d):
            y = cr + mid * dr; y = y - floor(y)
            s = 1 if y < F(1, 2) else -1     # ||.|| increasing or decreasing here
            lines.append((s * dr, nrm(cr + mid * dr) - s * dr * mid))
        # subdivide [a,b] at crossings of the lines, keep the min piece
        cuts = {a, b}
        for (m1, c1), (m2, c2) in itertools.combinations(lines, 2):
            if m1 != m2:
                t = (c2 - c1) / (m1 - m2)
                if a < t < b: cuts.add(t)
        cuts = sorted(cuts)
        for p, q in zip(cuts, cuts[1:]):
            m = (p + q) / 2
            sl, ic = min(lines, key=lambda l: l[0] * m + l[1])
            out.append((p, q, sl, ic))
    # merge consecutive identical affine pieces
    merged = []
    for pc in out:
        if merged and merged[-1][2] == pc[2] and merged[-1][3] == pc[3] and merged[-1][1] == pc[0]:
            merged[-1] = (merged[-1][0], pc[1], pc[2], pc[3])
        else:
            merged.append(pc)
    return merged

def f_eval(pieces, th):
    th = th - floor(th)
    for a, b, sl, ic in pieces:
        if a <= th <= b: return sl * th + ic
    raise ValueError

# --------------------------------------------------------------- components
class Quad:
    pass

def build_quadruples(u=U[0], v=U[1]):
    """All (i,j,eps,ell) with the data of the component; flag in_Y."""
    n = len(u); quads = []
    for i, j in itertools.combinations(range(n), 2):
        for eps in (1, -1):
            a, b = u[i] - eps * u[j], v[i] - eps * v[j]
            assert (a, b) != (0, 0)
            K = gcd(abs(a), abs(b)); ap, bp = a // K, b // K
            g, x0, y0 = egcd(ap, bp)
            if g < 0: g, x0, y0 = -g, -x0, -y0
            assert ap * x0 + bp * y0 == 1
            for ell in range(K):
                Q = Quad()
                Q.i, Q.j, Q.eps, Q.ell, Q.K = i, j, eps, ell, K
                Q.ap, Q.bp, Q.x0, Q.y0 = ap, bp, x0, y0
                Q.a, Q.b = a, b
                Q.c = [F(ell, K) * (x0 * u[r] + y0 * v[r]) for r in range(n)]
                Q.d = [-bp * u[r] + ap * v[r] for r in range(n)]
                Q.pieces = pieces_of_min_tents(Q.c, Q.d)
                Q.fmax = max(max(sl * a_ + ic, sl * b_ + ic) for a_, b_, sl, ic in Q.pieces)
                Q.inY = (Q.fmax == SIXTH)
                assert Q.fmax <= SIXTH
                if Q.inY: analyse_Y(Q)
                quads.append(Q)
    return quads

def analyse_Y(Q):
    """maximisers tau_h, slopes, linear zones, f_out, threshold q_0."""
    P = Q.pieces
    taus = set()
    for a, b, sl, ic in P:
        if sl * a + ic == SIXTH: taus.add(a % 1)
        if sl * b + ic == SIXTH: taus.add(b % 1)
        assert not (sl == 0 and ic == SIXTH), "flat maximum"
    Q.taus = sorted(taus)
    Q.zones = []      # (tau, lam_minus, rho_minus, lam_plus, rho_plus)
    for tau in Q.taus:
        # right zone: pieces starting at tau going right (cyclically)
        def walk(direction):
            # returns (lambda, rho): f(tau + dir*t) = 1/6 - lambda t for 0<=t<=rho, maximal rho
            lam = None; rho = F(0); pos = tau
            for _ in range(2 * len(P) + 2):
                # find piece adjacent to pos in the given direction
                pp = pos - floor(pos)
                found = None
                for a, b, sl, ic in P:
                    if direction > 0 and a == pp: found = (a, b, sl, ic); break
                    if direction < 0 and b == pp: found = (a, b, sl, ic); break
                    if direction > 0 and pp == 1 and a == 0: found = (a, b, sl, ic); break
                    if direction < 0 and pp == 0 and b == 1: found = (a, b, sl, ic); break
                a, b, sl, ic = found
                L = -direction * sl
                if lam is None: lam = L
                if L != lam: break
                rho += (b - a); pos = pos + direction * (b - a)
            return lam, rho
        lp, rp = walk(+1); lm, rm = walk(-1)
        assert lp > 0 and lm > 0
        Q.zones.append((tau, lm, rm, lp, rp))
    Q.rho_min = min(min(z[2], z[4]) for z in Q.zones)
    Q.lam_max = max(max(z[1], z[3]) for z in Q.zones)
    # f_out: max of f on the complement of the open zones
    # evaluate at all piece endpoints not strictly inside a zone, plus zone boundary points
    def in_open_zone(th):
        for tau, lm, rm, lp, rp in Q.zones:
            dl = (tau - th) % 1; dr = (th - tau) % 1
            if 0 < dr < rp or 0 < dl < rm or th % 1 == tau: return True
        return False
    cands = set()
    for a, b, sl, ic in P: cands.add(a % 1); cands.add(b % 1)
    for tau, lm, rm, lp, rp in Q.zones: cands.add((tau + rp) % 1); cands.add((tau - rm) % 1)
    Q.f_out = max(f_eval(P, th) for th in cands if not in_open_zone(th))
    assert Q.f_out < SIXTH
    q1 = ceil(1 / Q.rho_min)
    q2 = floor(Q.lam_max / (2 * (SIXTH - Q.f_out))) + 1
    Q.q0 = max(q1, q2)

# ----------------------------------------------------- coset for given (A,B)
def coset_data(Q, A, B):
    """q, theta_0 for T = <Au+Bv> on this component; None if L = 0."""
    L = Q.a * A + Q.b * B
    if L == 0: return None
    c = Q.x0 * B - Q.y0 * A
    q = abs(L) // Q.K
    th0 = (F(Q.ell * c, L)) % 1
    return q, th0

def coset_max_bruteforce(Q, A, B):
    r = coset_data(Q, A, B)
    if r is None: return None
    q, th0 = r
    return max(f_eval(Q.pieces, th0 + F(m, q)) for m in range(q))

def deficit_formula(Q, A, B):
    """1/6 - coset max, by the Approx formula (valid when q >= q0)."""
    r = coset_data(Q, A, B)
    if r is None: return None
    q, th0 = r
    best = None
    for tau, lm, rm, lp, rp in Q.zones:
        am = (tau - th0) % F(1, q)      # Approx^-
        ap = (th0 - tau) % F(1, q)      # Approx^+
        for val in (lm * am, lp * ap):
            if best is None or val < best: best = val
    return best, q

def ml_bruteforce(A, B):
    w = [A * U[0][r] + B * U[1][r] for r in range(6)]
    vs = sorted({abs(x) for x in w})
    dens = set()
    for x in vs: dens.add(2 * x)
    for x, y in itertools.combinations(vs, 2): dens.add(x + y); dens.add(abs(x - y))
    best = F(0)
    for bden in dens:
        if bden == 0: continue
        for a in range(1, bden):
            r = min(min((a * x) % bden, bden - (a * x) % bden) for x in vs)
            if F(r, bden) > best: best = F(r, bden)
    return best

if __name__ == "__main__":
    quads = build_quadruples()
    Y = [Q for Q in quads if Q.inY]
    print("quadruples:", len(quads), " in Y:", len(Y))
    for Q in Y:
        print(" (%d,%d,%+d,%d) K=%d  L = %d A + %d B   taus=%s  q0=%d  rho_min=%s lam_max=%s f_out=%s" %
              (Q.i + 1, Q.j + 1, Q.eps, Q.ell, Q.K, Q.a, Q.b, [str(t) for t in Q.taus], Q.q0, Q.rho_min, Q.lam_max, Q.f_out))
    print("non-Y gaps eta = 1/6 - fmax:", sorted({str(SIXTH - Q.fmax) for Q in quads if not Q.inY}))
    # validation 1: formula vs brute force on the coset
    import random
    random.seed(1); bad = 0; tested = 0
    for _ in range(3000):
        A = random.randint(0, 60); B = random.randint(-60, 60)
        if gcd(A, B) != 1: continue
        for Q in Y:
            r = coset_data(Q, A, B)
            if r is None or r[0] < Q.q0: continue
            tested += 1
            d, q = deficit_formula(Q, A, B)
            if SIXTH - d != coset_max_bruteforce(Q, A, B):
                bad += 1; print("MISMATCH", (Q.i, Q.j, Q.eps, Q.ell), A, B)
    print("formula vs brute force on cosets: tested", tested, "mismatches", bad)
    # validation 2: predicted ML vs brute force ML for all (A,B) up to 40, outside thresholds
    bad = 0; tested = 0; skipped = 0
    for A in range(0, 41):
        for B in range(-40, 41):
            if gcd(A, B) != 1 or (A == 0 and B != 1): continue
            w = [A * U[0][r] + B * U[1][r] for r in range(6)]
            if 0 in w or len({abs(x) for x in w}) < 6: continue
            defs = []
            ok = True
            for Q in Y:
                r = coset_data(Q, A, B)
                if r is None: continue
                if r[0] < Q.q0: ok = False; break
                d, q = deficit_formula(Q, A, B); defs.append(d)
            if not ok: skipped += 1; continue
            pred = SIXTH - min(defs)
            true = ml_bruteforce(A, B)
            eta_min = min(SIXTH - Q.fmax for Q in quads if not Q.inY)
            if true < SIXTH:
                tested += 1
                if min(defs) < eta_min and pred != true:
                    bad += 1; print("ML MISMATCH", A, B, pred, true)
    print("predicted ML vs brute force (near-tight, all q >= q0):", tested, "mismatches", bad, " skipped", skipped)
