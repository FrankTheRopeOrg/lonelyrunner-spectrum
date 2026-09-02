/* Exhaustive survey of near-tight sextuples (ML < 1/6) with distinct speeds
 * v_1 < ... < v_6 <= N and gcd 1.  For each: exact ML = p/q, k = q - 6p,
 * P = q/k, Q, and membership in U_A, U_B, U_C up to signed permutations.
 * usage: survey N [n]   (n = number of speeds, default 6; near-tight = ML < 1/n)
 * Output: one line per near-tight tuple, then a summary. */
#include <stdio.h>
#include <omp.h>
#include "lrk.h"

static int same_multiset4(int64_t *a, int64_t *b) {
    /* a, b arrays of 4 nonneg ints */
    int64_t x[4], y[4]; memcpy(x, a, sizeof x); memcpy(y, b, sizeof y);
    for (int i = 0; i < 4; i++) for (int j = i + 1; j < 4; j++) { if (x[j] < x[i]) { int64_t t = x[i]; x[i] = x[j]; x[j] = t; } if (y[j] < y[i]) { int64_t t = y[i]; y[i] = y[j]; y[j] = t; } }
    for (int i = 0; i < 4; i++) if (x[i] != y[i]) return 0;
    return 1;
}
static int in_UC(const int64_t *v) {  /* v sorted positive, 6 entries */
    for (int i = 0; i < 6; i++) for (int j = 0; j < 6; j++) if (i != j) for (int s = -1; s <= 1; s += 2) {
        int64_t A = v[i], B = s * v[j];
        int64_t need[4] = { llabs(A + B), llabs(2*A + B), llabs(3*A + B), llabs(3*A + 2*B) };
        int64_t rest[4]; int m = 0; for (int t = 0; t < 6; t++) if (t != i && t != j) rest[m++] = v[t];
        if (same_multiset4(need, rest)) return 1;
    }
    return 0;
}
static int in_UV(const int64_t *v, const int64_t *V) { /* five of the speeds are b*V */
    for (int drop = 0; drop < 6; drop++) {
        int64_t r[5]; int m = 0; for (int t = 0; t < 6; t++) if (t != drop) r[m++] = v[t];
        int64_t b = r[0] / V[0]; if (b * V[0] != r[0]) continue;
        int ok = 1; for (int t = 1; t < 5; t++) if (r[t] != b * V[t]) ok = 0;
        if (ok) return 1;
    }
    return 0;
}
static int in_UV_n5(const int64_t *v, const int64_t *V) { /* four of the five speeds are b*V */
    for (int drop = 0; drop < 5; drop++) {
        int64_t r[4]; int m = 0; for (int t = 0; t < 5; t++) if (t != drop) r[m++] = v[t];
        int64_t b = r[0] / V[0]; if (b * V[0] != r[0]) continue;
        int ok = 1; for (int t = 1; t < 4; t++) if (r[t] != b * V[t]) ok = 0;
        if (ok) return 1;
    }
    return 0;
}
/* n = 4: U^1 = <(0,1,2,3),(1,0,0,0)>: three speeds b*(1,2,3); U^2 = <(1,0,1,1),(1,1,0,2)>: {|A|,|B|,|A+B|,|A+2B|} */
static int in_U1_n4(const int64_t *v) {
    for (int drop = 0; drop < 4; drop++) {
        int64_t r[3]; int m = 0; for (int t = 0; t < 4; t++) if (t != drop) r[m++] = v[t];
        if (r[1] == 2*r[0] && r[2] == 3*r[0]) return 1;
    }
    return 0;
}
static int in_U2_n4(const int64_t *v) {
    for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) if (i != j) for (int s = -1; s <= 1; s += 2) {
        int64_t A = v[i], B = s * v[j];
        int64_t need[2] = { llabs(A + B), llabs(A + 2*B) };
        int64_t rest[2]; int m = 0; for (int t = 0; t < 4; t++) if (t != i && t != j) rest[m++] = v[t];
        if ((need[0] == rest[0] && need[1] == rest[1]) || (need[0] == rest[1] && need[1] == rest[0])) return 1;
    }
    return 0;
}

int main(int argc, char **argv) {
    int64_t N = argc > 1 ? atoll(argv[1]) : 40;
    int n = argc > 2 ? atoi(argv[2]) : 6;
    int64_t VA[5] = {1,2,3,4,5}, VB[5] = {1,3,4,5,9};
    long total = 0, nt = 0, cA = 0, cB = 0, cC = 0, cNone = 0, kcount[16] = {0};
    /* enumerate by first two speeds in parallel */
    #pragma omp parallel for schedule(dynamic,1) reduction(+:total,nt,cA,cB,cC,cNone)
    for (int64_t v1 = 1; v1 <= N - n + 1; v1++) {
        int64_t v[8];
        v[0] = v1;
        /* iterative enumeration of remaining n-1 speeds */
        int64_t idx[8]; int d = 1; idx[1] = v1 + 1;
        while (d >= 1) {
            if (idx[d] > N - (n - 1 - d)) { d--; if (d >= 1) idx[d]++; continue; }
            v[d] = idx[d];
            if (d == n - 1) {
                total++;
                int64_t g = 0; for (int i = 0; i < n; i++) g = gcd64(g, v[i]);
                if (g == 1 && near_tight(v, n, n)) {
                    int64_t p, q, Q; ml_exact(v, n, &p, &q, &Q);
                    int64_t k = q - n * p; int64_t P = (k > 0 && q % k == 0) ? q / k : -1;
                    int A = 0, B = 0, C = 0;
                    if (n == 6) { A = in_UV(v, VA); B = in_UV(v, VB); C = in_UC(v); }
                    if (n == 4) { A = in_U1_n4(v); B = in_U2_n4(v); C = 0; }
                    if (n == 5) { int64_t V1[4] = {1,2,3,4}, V2[4] = {1,3,4,7}; A = in_UV_n5(v, V1); B = in_UV_n5(v, V2); C = 0; }
                    nt++; if (A) cA++; if (B) cB++; if (C) cC++; if (!A && !B && !C) cNone++;
                    #pragma omp critical
                    { if (k >= 0 && k < 16) kcount[k]++;
                      printf("NT");
                      for (int i = 0; i < n; i++) printf(" %lld", (long long)v[i]);
                      printf("  ML=%lld/%lld k=%lld P=%lld Q=%lld maxv/Q=%.3f tori=%s%s%s%s\n", (long long)p, (long long)q, (long long)k, (long long)P, (long long)Q,
                          (double)v[n-1] / (double)Q, A ? "A" : "", B ? "B" : "", C ? "C" : "", (!A && !B && !C) ? "-" : "");
                      fflush(stdout); }
                }
                idx[d]++;
            } else { d++; idx[d] = idx[d-1] + 1; }
        }
    }
    printf("N=%lld n=%d: tuples %ld, near-tight primitive %ld; in U_A %ld, U_B %ld, U_C %ld, in none %ld\n", (long long)N, n, total, nt, cA, cB, cC, cNone);
    printf("k distribution:"); for (int k = 0; k < 16; k++) if (kcount[k]) printf(" k=%d:%ld", k, kcount[k]); printf("\n");
    return 0;
}
