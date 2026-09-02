#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Theorem 3.7 (the three critical subtori of (R/Z)^6): reproduces every number
quoted in Section 3 of the paper.  About one minute with numpy.
"""
import sys, os, time
from fractions import Fraction as F
from lr_ml import ml_tuple, TIGHT_QUINTUPLES
from lr_subtorus import exceeds_sixth
from lr_classify import classify, orbit_canonical
from lr_lemma_s5 import run as lemma_s5_run
from lr_lemma_s6 import run as lemma_s6_run

KNOWN = {
    "U_A": ([1, 2, 3, 4, 5, 0], [0, 0, 0, 0, 0, 1]),
    "U_B": ([1, 3, 4, 5, 9, 0], [0, 0, 0, 0, 0, 1]),
    "U_C": ([1, 0, 1, 2, 3, 3], [0, 1, 1, 1, 1, 2]),
}
FAIL = []
def check(label, got, want):
    ok = got == want
    print("  [%s] %-58s %s" % ("ok" if ok else "FAIL", label, got))
    if not ok: FAIL.append((label, got, want))

def main():
    t0 = time.time()
    print("(1) tight quintuples")
    for T in TIGHT_QUINTUPLES: check("ML%s" % (T,), ml_tuple(T), F(1, 6))
    print("(2) D(U) = 1/3 for U_A, U_B, U_C")
    for nm, (u, v) in KNOWN.items(): check("D(%s) = 1/3" % nm, not exceeds_sixth(u, v), True)
    print("(3) Lemma 3.6 (one free repeat line implies two): branch and prune")
    check("surviving solution spaces over the ten canonical w", lemma_s5_run(), 0)
    print("(4) Lemma 3.5 (d = 6): the sixteen-case check")
    check("new directions required (two slots available)", lemma_s6_run() >= 4, True)
    print("(5) enumeration of pairs of tight repeat vectors")
    classes, st = classify(verbose=True)
    check("canonical first vectors", st["w1"], 10)
    check("tight repeat vectors up to global sign", st["w2"], 115200)
    check("non-parallel pairs examined", st["pairs"], 1151990)
    check("survivors of the rational-witness filter", st["survivors"], 190)
    check("survivors with D(U) = 1/3 (exact sweep)", st["exact_ok"], 190)
    check("distinct planes", st["planes"], 22)
    check("orbits under the 46,080 signed permutations", len(classes), 3)
    kc = {nm: orbit_canonical([list(u), list(v)]) for nm, (u, v) in KNOWN.items()}
    found = sorted(nm for key in classes for nm, k in kc.items() if k == key)
    check("the orbits are U_A, U_B, U_C", found, ["U_A", "U_B", "U_C"])
    print("%s   (%.0f s)" % ("ALL CHECKS PASSED" if not FAIL else "FAILURES: %s" % FAIL, time.time() - t0))
    sys.exit(1 if FAIL else 0)

if __name__ == "__main__":
    main()
