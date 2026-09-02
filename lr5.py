#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Classification of the two-dimensional proper subtori U of (R/Z)^5 with
D(U) = 3/10, i.e. ML(U) = 1/5 -- the n = 5 analogue of Theorem 6.7 of the note.

Input from the literature: the tight quadruples are (1,2,3,4) and (1,3,4,7)
up to dilation (Y.G. Chen, J. Number Theory 37 (1991); see Perarnau-Serra,
"The Lonely Runner Conjecture turns 60", Section on tight instances).

Argument (same skeleton as Lemmas 6.1-6.6 of the note, one dimension down):

  L1  no proper 2-dim subtorus W of (R/Z)^4 has ML(W) = 1/5: a zero-free
      primitive w in W has ML(w) <= 1/5; with <= 3 distinct absolute values
      ML(w) >= 1/4, so it has 4 distinct values, ML(w) >= 1/5 by LRC(4), hence
      tight, hence entries <= 7; infinitely many primitive directions, at most
      4 vanishing lines, bounded height: contradiction.
  C2  hence U is in no {x_p = +-x_q}: the five coordinate forms are pairwise
      distinct up to sign.
  L3  a zero-free w in U with |w_p| = |w_q| has exactly 4 distinct absolute
      values forming a tight quadruple ("tight repeat vector").
  L4  d = number of projective classes of the forms.  For a pair of classes
      (k,l) the repeat directions are >= 2 max(m_k, m_l) distinct lines, at most
      d - 2 of which are vanishing lines.  d = 2: >= 6 free; d = 3: >= 3;
      d = 4: >= 2.  So d <= 4 gives two free repeat directions outright.
  L6  d = 5 (all forms pairwise non-proportional): the Pluecker closure
      argument of Lemma 6.6 gives >= 1 free direction (with only ONE slot left
      the contradiction is immediate: c not in {+-1, +-1/2} needs four new
      directions, c in that set needs two, one slot is available).
  L5  one free direction implies two: branch and prune, verified below.

  So U is spanned by two non-parallel tight repeat vectors; enumerate all pairs
  and decide ML(U) = 1/5 exactly by a sweep.
"""
import itertools, sys
from fractions import Fraction as F
from math import gcd, floor, ceil

N0 = 5                        # threshold 1/N0
TIGHT4 = [(1, 2, 3, 4), (1, 3, 4, 7)]
DIM = 5

# ------------------------------------------------------------ exact ML kernel
def ml_tuple(speeds):
    vs = sorted({abs(s) for s in speeds})
    if 0 in vs:
        return F(0)
    best = F(0)
    dens = set()
    for v in vs:
        dens.add(2 * v)
    for a, b in itertools.combinations(vs, 2):
        dens.add(a + b); dens.add(abs(a - b))
    for b in dens:
        if b == 0: continue
        for a in range(1, b):
            r = min(min((a * v) % b, b - (a * v) % b) for v in vs)
            if r * best.denominator > best.numerator * b:
                best = F(r, b)
    return best

# ------------------------------------------------- exact decision ML(U) > 1/N0
def _t_intervals(s, ui, vi, lo, hi):
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

def _critical_s(u, v, n0):
    """s where two interval endpoints coincide: n0 (s u_i + t v_i) = c, c = +-1 mod n0."""
    n = len(u)
    S = {F(0)}
    good = lambda c: c % n0 in (1, n0 - 1)
    for i in range(n):
        if v[i] == 0 and u[i] != 0:
            for c in range(-n0 * abs(u[i]) - n0, n0 * abs(u[i]) + n0 + 1):
                if good(c):
                    s = F(c, n0 * u[i])
                    if 0 <= s < 1: S.add(s)
    for i, j in itertools.combinations(range(n), 2):
        det = u[j] * v[i] - u[i] * v[j]
        if det == 0: continue
        Ri = n0 * (abs(u[i]) + abs(v[i])) + 2 * n0
        Rj = n0 * (abs(u[j]) + abs(v[j])) + 2 * n0
        ci = [c for c in range(-Ri, Ri + 1) if good(c)]
        cj = [c for c in range(-Rj, Rj + 1) if good(c)]
        for c in ci:
            for cp in cj:
                s = F(cp * v[i] - c * v[j], n0 * det)
                if 0 <= s < 1: S.add(s)
    return sorted(S)

def exceeds(u, v, n0=N0):
    """True iff ML(U) > 1/n0 for U = <u, v>: some (s,t) with all frac(s u_r + t v_r) in (1/n0, 1-1/n0)."""
    lo, hi = F(1, n0), F(n0 - 1, n0)
    S = _critical_s(u, v, n0) + [F(1)]
    for a, b in zip(S, S[1:]):
        if a != b and _feasible_at((a + b) / 2, u, v, lo, hi):
            return True
    return False

# ------------------------------------------------------------- enumeration
def canonical_first_vectors():
    out = []
    for T in TIGHT4:
        for a in T:
            out.append((a, a) + tuple(sorted(x for x in T if x != a)))
    return out

def all_repeat_vectors():
    seen = set()
    for T in TIGHT4:
        for rep in T:
            others = [x for x in T if x != rep]
            for p, q in itertools.combinations(range(DIM), 2):
                free = [i for i in range(DIM) if i not in (p, q)]
                for perm in itertools.permutations(others):
                    base = [0] * DIM
                    base[p] = base[q] = rep
                    for i, val in zip(free, perm):
                        base[i] = val
                    for signs in itertools.product((1, -1), repeat=DIM):
                        w = tuple(s * b for s, b in zip(signs, base))
                        if w[0] < 0: w = tuple(-x for x in w)
                        seen.add(w)
    return sorted(seen)

def rref(rows):
    M = [[F(x) for x in r] for r in rows]
    piv = 0
    for col in range(len(M[0])):
        r = next((i for i in range(piv, len(M)) if M[i][col] != 0), None)
        if r is None: continue
        M[piv], M[r] = M[r], M[piv]
        pv = M[piv][col]
        M[piv] = [x / pv for x in M[piv]]
        for i in range(len(M)):
            if i != piv and M[i][col] != 0:
                f = M[i][col]
                M[i] = [x - f * y for x, y in zip(M[i], M[piv])]
        piv += 1
        if piv == len(M): break
    return tuple(tuple(r) for r in M if any(x != 0 for x in r))

_PERMS = list(itertools.permutations(range(DIM)))
_SIGNS = [s for s in itertools.product((1, -1), repeat=DIM) if s[0] == 1]

def orbit_canonical(basis):
    best = None
    for p in _PERMS:
        b0 = [[r[p[i]] for i in range(DIM)] for r in basis]
        for s in _SIGNS:
            k = rref([[s[i] * r[i] for i in range(DIM)] for r in b0])
            if best is None or k < best: best = k
    return best

def witness_filter(w1, cands, N=60):
    """cheap rational-witness rejection: a grid point (p/N, q/N) with all
    coordinates strictly inside (1/N0, 1-1/N0) certifies ML(U) > 1/N0."""
    lo, hi = N // N0, N - N // N0
    alive = []
    grid = [(p, q) for p in range(N) for q in range(N) if (p, q) != (0, 0)]
    for w2 in cands:
        if all(w1[0] * w2[j] == w2[0] * w1[j] for j in range(1, DIM)):
            continue  # parallel
        killed = False
        for p, q in grid:
            ok = True
            for r in range(DIM):
                x = (p * w1[r] + q * w2[r]) % N
                if not (lo < x < hi): ok = False; break
            if ok: killed = True; break
        if not killed: alive.append(w2)
    return alive

# ---------------------------------------------- Lemma L5: one free -> two
def _intersect_hyperplane(B, a):
    """B: basis (list of rows) of a subspace of Q^DIM; return basis of B ∩ {a.z = 0}."""
    coeffs = [sum(F(a[i]) * b[i] for i in range(DIM)) for b in B]
    piv = next((i for i, c in enumerate(coeffs) if c != 0), None)
    if piv is None: return B
    out = []
    for i, b in enumerate(B):
        if i == piv: continue
        f = coeffs[i] / coeffs[piv]
        out.append([b[j] - f * B[piv][j] for j in range(DIM)])
    return out

def _inside_span_w(B, w):
    for b in B:
        if any(b[i] * w[j] != b[j] * w[i] for i in range(DIM) for j in range(DIM)):
            return False
    return True

def solutions_with_at_most_one_free_direction(w):
    """Exact branch and prune: z such that every repeat line other than that of
    w carries a vanishing coordinate.  Returns number of surviving spaces (0 expected)."""
    conds = []
    for p, q in itertools.combinations(range(DIM), 2):
        for eps in (1, -1):
            if w[p] - eps * w[q] == 0: continue
            hyper = []
            for r in range(DIM):
                a = [0] * DIM
                a[r] += (w[p] - eps * w[q]); a[p] -= w[r]; a[q] += eps * w[r]
                hyper.append(a)
            conds.append(hyper)
    start = [[F(int(i == j)) for j in range(DIM)] for i in range(DIM)]
    frontier = {rref(start): start}
    for hyper in conds:
        nxt = {}
        for key, B in frontier.items():
            for a in hyper:
                B2 = _intersect_hyperplane(B, a)
                if not B2 or _inside_span_w(B2, w): continue
                nxt.setdefault(rref(B2), B2)
        frontier = nxt
        if not frontier: break
    return len(frontier)

def main():
    print("Tight quadruples:", [(T, ml_tuple(T)) for T in TIGHT4])
    print("L5 (one free direction implies two), branch and prune:")
    tot = 0
    for w in canonical_first_vectors():
        s = solutions_with_at_most_one_free_direction(w); tot += s
        print("   w =", w, " surviving spaces:", s)
    print("  total:", tot)
    W1 = canonical_first_vectors(); W2 = all_repeat_vectors()
    print("canonical first vectors:", len(W1), " tight repeat vectors:", len(W2))
    survivors = []
    for w1 in W1:
        alive = witness_filter(w1, W2)
        survivors += [(w1, w2) for w2 in alive]
        print("   w1 = %-18s survivors %d" % (str(w1), len(alive)))
    print("pairs surviving the witness filter:", len(survivors))
    exact = [(a, b) for a, b in survivors if not exceeds(list(a), list(b))]
    print("confirmed exactly to have ML(U) = 1/5:", len(exact))
    planes = {}
    for a, b in exact:
        planes.setdefault(rref([a, b]), []).append((a, b))
    print("distinct planes:", len(planes))
    classes = {}
    for _, reps in planes.items():
        classes.setdefault(orbit_canonical([list(reps[0][0]), list(reps[0][1])]), []).append(reps[0])
    print("orbits under signed permutations:", len(classes))
    for key, reps in classes.items():
        a, b = reps[0]
        print("   class spanned by", a, ",", b, "   rref:", [[str(x) for x in r] for r in key])

if __name__ == "__main__":
    main()
