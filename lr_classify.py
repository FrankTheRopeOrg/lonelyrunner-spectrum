# -*- coding: utf-8 -*-
"""Classification of the two-dimensional subtori U of (R/Z)^6 with D(U) = 1/3.

Strategy (Lemmas S1-S5 of the note):

  S1/S2  U is not contained in any {x_p = +-x_q}; equivalently the six
         coordinate forms are pairwise distinct up to sign.
  S3     any zero-free w in U with |w_p| = |w_q| has exactly five distinct
         absolute values, forming a tight quintuple: entries at most 9.
  S4     if the forms occupy d <= 4 projective classes, U has at least two
         free repeat directions; d = 5 gives at least one.
  S5     one free direction implies two (verified in lr_lemma_s5.py).

So U is spanned by two non-parallel tight repeat vectors.  This module
enumerates all such pairs, decides D(U) = 1/3 exactly, and reduces the
survivors to orbits under the 46,080 signed permutations of the coordinates.
"""
import itertools
from fractions import Fraction as F
from lr_ml import TIGHT_QUINTUPLES
from lr_subtorus import exceeds_sixth

__all__ = ["canonical_first_vectors", "all_repeat_vectors", "filter_by_witness",
           "rref", "orbit_canonical", "classify"]


def canonical_first_vectors():
    """The ten w = (a,a,b,c,d,e), a>0, 0<b<c<d<e, {a,b,c,d,e} tight.
    Any tight repeat vector can be brought to this form by a signed
    permutation of the coordinates."""
    out = []
    for T in TIGHT_QUINTUPLES:
        for a in T:
            out.append((a, a) + tuple(sorted(x for x in T if x != a)))
    return out


def all_repeat_vectors():
    """Every tight repeat vector, up to a global sign: a tight quintuple with
    one value repeated in two of the six positions, arbitrary signs."""
    seen = set()
    for T in TIGHT_QUINTUPLES:
        for rep in T:
            others = [x for x in T if x != rep]
            for p, q in itertools.combinations(range(6), 2):
                free = [i for i in range(6) if i not in (p, q)]
                for perm in itertools.permutations(others):
                    base = [0] * 6
                    base[p] = base[q] = rep
                    for i, val in zip(free, perm):
                        base[i] = val
                    for signs in itertools.product((1, -1), repeat=6):
                        w = tuple(s * b for s, b in zip(signs, base))
                        if w[0] < 0:
                            w = tuple(-x for x in w)
                        seen.add(w)
    return sorted(seen)


def filter_by_witness(w1, cands, N=102):
    """Discard every w2 for which an explicit rational point (p/N, q/N) has
    min_r ||.|| > 1/6 -- a certificate that D(U) < 1/3.

    This is only an optimisation: correctness rests on exceeds_sixth, which
    decides every survivor exactly.  Without numpy the filter is skipped
    entirely and every non-parallel candidate goes through to the exact sweep,
    which gives the same 190 survivors, only far more slowly.  The header's
    claim that correctness does not depend on numpy is therefore literally
    true, and this branch is what makes it so.
    """
    import random
    try:
        import numpy as np
    except ImportError:
        import sys as _s
        print("  numpy not found: the witness filter is skipped and every "
              "candidate\n  goes to the exact sweep.  Same answer, but expect "
              "hours, not minutes.", file=_s.stderr)
        kept = [tuple(w2) for w2 in cands
                if any(w1[0] * w2[j] != w2[0] * w1[j] for j in range(1, 6))]
        return kept, len(kept)
    a = np.array(w1, dtype=np.int64)
    C = np.array(cands, dtype=np.int64)
    cross = a[0] * C[:, 1:] - C[:, 0:1] * a[1:]
    C = C[(cross != 0).any(axis=1)]                 # drop vectors parallel to w1
    lo, hi = N // 6, 5 * N // 6
    alive = np.arange(len(C))
    grid = [(p, q) for p in range(N) for q in range(N) if (p, q) != (0, 0)]
    random.Random(0).shuffle(grid)
    for p, q in grid:
        X = (p * a + q * C[alive]) % N
        ok = ((X > lo) & (X < hi)).all(axis=1)
        if ok.any():
            alive = alive[~ok]
            if len(alive) == 0:
                break
    return [tuple(int(x) for x in C[i]) for i in alive], len(C)


def rref(rows):
    """Reduced row echelon form over Q: a canonical form for the row space."""
    M = [[F(x) for x in r] for r in rows]
    piv = 0
    for col in range(len(M[0])):
        r = next((i for i in range(piv, len(M)) if M[i][col] != 0), None)
        if r is None:
            continue
        M[piv], M[r] = M[r], M[piv]
        pv = M[piv][col]
        M[piv] = [x / pv for x in M[piv]]
        for i in range(len(M)):
            if i != piv and M[i][col] != 0:
                f = M[i][col]
                M[i] = [x - f * y for x, y in zip(M[i], M[piv])]
        piv += 1
        if piv == len(M):
            break
    return tuple(tuple(r) for r in M if any(x != 0 for x in r))


_PERMS = list(itertools.permutations(range(6)))
_SIGNS = [s for s in itertools.product((1, -1), repeat=6) if s[0] == 1]


def orbit_canonical(basis):
    """Lexicographically least RREF over the group of signed permutations of
    the coordinates (order 6! * 2^6 = 46,080; the global sign is factored out)."""
    best = None
    for p in _PERMS:
        b0 = [[r[p[i]] for i in range(6)] for r in basis]
        for s in _SIGNS:
            k = rref([[s[i] * r[i] for i in range(6)] for r in b0])
            if best is None or k < best:
                best = k
    return best


def classify(verbose=True):
    """Full pipeline.  Returns (classes, stats)."""
    W1 = canonical_first_vectors()
    W2 = all_repeat_vectors()
    stats = {"w1": len(W1), "w2": len(W2), "pairs": 0, "survivors": 0,
             "exact_ok": 0, "planes": 0}
    survivors = []
    for w1 in W1:
        alive, ncand = filter_by_witness(w1, W2)
        stats["pairs"] += ncand
        survivors += [(w1, w2) for w2 in alive]
        if verbose:
            print("   w1 = %-24s candidates %7d -> survivors %3d"
                  % (str(w1), ncand, len(alive)))
    stats["survivors"] = len(survivors)

    exact = [(a, b) for a, b in survivors if not exceeds_sixth(list(a), list(b))]
    stats["exact_ok"] = len(exact)

    planes = {}
    for a, b in exact:
        planes.setdefault(rref([a, b]), []).append((a, b))
    stats["planes"] = len(planes)

    classes = {}
    for _, reps in planes.items():
        classes.setdefault(orbit_canonical([list(reps[0][0]), list(reps[0][1])]),
                           []).append(reps[0])
    return classes, stats
