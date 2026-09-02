/* Independent exact maximum-loneliness kernel (integer arithmetic only).
 *
 * ML(v) = max_t min_i ||t v_i||.  The maximum of the piecewise-linear
 * function f(t) = min_i ||t v_i|| is attained at a time where either a single
 * active runner sits at the antipode (t v_i = 1/2 mod 1) or two active runners
 * cross (t (v_i +- v_j) in Z).  So it suffices to scan t = a/b with
 * b in { 2 v_i, v_i + v_j, |v_i - v_j| } and 0 <= a < b.
 * For t = a/b, ||t v_i|| = r/b with r = min(a v_i mod b, b - a v_i mod b).
 * Fractions are compared by cross multiplication in 128-bit integers.
 */
#ifndef LRK_H
#define LRK_H
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef __int128 i128;

static inline int64_t gcd64(int64_t a, int64_t b) {
    if (a < 0) a = -a; if (b < 0) b = -b;
    while (b) { int64_t t = a % b; a = b; b = t; }
    return a;
}

/* min_i ||a v_i / b|| as numerator over b */
static inline int64_t minr(const int64_t *v, int n, int64_t a, int64_t b) {
    int64_t best = b;
    for (int i = 0; i < n; i++) {
        int64_t r = (a * v[i]) % b;
        if (r < 0) r += b;
        if (b - r < r) r = b - r;
        if (r < best) { best = r; if (best == 0) return 0; }
    }
    return best;
}

/* collect candidate denominators (distinct) */
static inline int cand_denoms(const int64_t *v, int n, int64_t *den) {
    int m = 0;
    for (int i = 0; i < n; i++) den[m++] = 2 * v[i];
    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) {
        den[m++] = v[i] + v[j];
        int64_t d = v[i] - v[j]; if (d < 0) d = -d; if (d > 0) den[m++] = d;
    }
    /* dedupe */
    int k = 0;
    for (int i = 0; i < m; i++) {
        int dup = 0;
        for (int j = 0; j < k; j++) if (den[j] == den[i]) { dup = 1; break; }
        if (!dup) den[k++] = den[i];
    }
    return k;
}

/* Exact ML: returns numerator p and denominator q (lowest terms) and the
 * smallest denominator Q of an optimal time (over the candidate set, which
 * contains every optimal time up to the representation a/b in lowest terms;
 * we reduce a/b to lowest terms before taking the minimum). */
static void ml_exact(const int64_t *v, int n, int64_t *p, int64_t *q, int64_t *Q) {
    int64_t den[128]; int nd = cand_denoms(v, n, den);
    int64_t bp = 0, bq = 1, bQ = 0;
    for (int d = 0; d < nd; d++) {
        int64_t b = den[d];
        for (int64_t a = 1; a < b; a++) {
            int64_t r = minr(v, n, a, b);
            if (r == 0) continue;
            /* compare r/b with bp/bq */
            i128 lhs = (i128)r * bq, rhs = (i128)bp * b;
            if (lhs > rhs) {
                int64_t g = gcd64(r, b); bp = r / g; bq = b / g;
                int64_t g2 = gcd64(a, b); bQ = b / g2;
            } else if (lhs == rhs) {
                int64_t g2 = gcd64(a, b); int64_t Qc = b / g2;
                if (Qc < bQ) bQ = Qc;
            }
        }
    }
    *p = bp; *q = bq; *Q = bQ;
}

/* 1 iff ML(v) < 1/n0, i.e. no time with min_i ||t v_i|| >= 1/n0.
 * Exact: we need  n0 * r < b  for every candidate. */
static inline int near_tight(const int64_t *v, int n, int n0) {
    int64_t den[128]; int nd = cand_denoms(v, n, den);
    for (int d = 0; d < nd; d++) {
        int64_t b = den[d];
        for (int64_t a = 1; a < b; a++) {
            int64_t r = minr(v, n, a, b);
            if ((i128)n0 * r >= b) return 0;
        }
    }
    return 1;
}
#endif
