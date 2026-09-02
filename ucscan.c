/* Scan the one-dimensional subtori <A u + B v> of U_C = <(1,0,1,2,3,3),(0,1,1,1,1,2)>
 * (the basis of Jain-Kravitz's U^7) for 0 <= A <= H, |B| <= H, gcd(A,B)=1.
 * For every near-tight proper subtorus with six distinct speeds: exact ML = p/q,
 * k = q - 6p, P = q/k; check P integer, odd, coprime to 6; q == Q.
 * Print per (A mod 6, B mod 6) the slope range where k = 3 occurs, and the
 * paper's U^7-basis table as a biconditional.
 * usage: ucscan H */
#include <stdio.h>
#include <omp.h>
#include "lrk.h"

static int k3_table_u7(int64_t A, int64_t B) {  /* paper's table, U^7 basis, A > 0 */
    int64_t a = ((A % 6) + 6) % 6, b = ((B % 6) + 6) % 6;
    if (a == 0 && b == 1) return 2*B + A < 0;
    if (a == 0 && b == 5) return 2*B - 3*A > 0;
    if (a == 5 && b == 0) return (B + 2*A > 0) && (3*B - 2*A < 0);
    if (a == 5 && b == 5) return (3*B - A > 0) && (B - 3*A < 0);
    return 0;
}

/* the same table rewritten directly in the Jain-Kravitz basis (A,B), A >= 0 */
static int k3_table_jk(int64_t A, int64_t B) {
    int64_t a = ((A % 6) + 6) % 6, b = ((B % 6) + 6) % 6;
    if (a == 0 && b == 5) return A + 2*B > 0;                       /* B/A > -1/2 */
    if (a == 0 && b == 1) return 5*A + 2*B < 0;                     /* B/A < -5/2 */
    if (a == 5 && b == 1) return (5*A + 3*B > 0) && (A - B > 0);    /* -5/3 < B/A < 1 */
    if (a == 5 && b == 2) return (4*A + B > 0) && (4*A + 3*B < 0);  /* -4 < B/A < -4/3 */
    return 0;
}
typedef struct { int64_t A, B, p, q, Q, k, P; } rec;

int main(int argc, char **argv) {
    int64_t H = argc > 1 ? atoll(argv[1]) : 100;
    int64_t u[6] = {1,0,1,2,3,3}, vv[6] = {0,1,1,1,1,2};
    rec *R = malloc(sizeof(rec) * (size_t)(2*H+2) * (H+2)); long nR = 0;
    long nprop = 0;
    #pragma omp parallel for schedule(dynamic,4) reduction(+:nprop)
    for (int64_t A = 0; A <= H; A++) {
        for (int64_t B = -H; B <= H; B++) {
            if (A == 0 && B != 1) continue;
            if (gcd64(A, B) != 1) continue;
            int64_t w[6], aw[6];
            for (int i = 0; i < 6; i++) { w[i] = A*u[i] + B*vv[i]; aw[i] = w[i] < 0 ? -w[i] : w[i]; }
            int ok = 1;
            for (int i = 0; i < 6 && ok; i++) { if (aw[i] == 0) ok = 0; for (int j = 0; j < i; j++) if (aw[i] == aw[j]) ok = 0; }
            if (!ok) continue;
            nprop++;
            if (!near_tight(aw, 6, 6)) continue;
            int64_t p, q, Q; ml_exact(aw, 6, &p, &q, &Q);
            int64_t k = q - 6*p;
            rec r = {A, B, p, q, Q, k, (k > 0 && q % k == 0) ? q / k : -1};
            #pragma omp critical
            R[nR++] = r;
        }
    }
    printf("H = %lld: proper subtori with six distinct speeds: %ld, near-tight: %ld\n", (long long)H, nprop, nR);
    long badP = 0, badOdd = 0, bad3 = 0, badQ = 0, fp = 0, fn = 0, fpj = 0, fnj = 0, cnt[6][6][2] = {{{0}}};
    /* slope ranges of k=3 per class in U_C basis: track min/max of B/A over k=3 and k=1 */
    double k3min[6][6], k3max[6][6], k1min[6][6], k1max[6][6];
    for (int a = 0; a < 6; a++) for (int b = 0; b < 6; b++) { k3min[a][b] = 1e9; k3max[a][b] = -1e9; k1min[a][b] = 1e9; k1max[a][b] = -1e9; }
    for (long i = 0; i < nR; i++) {
        rec *r = &R[i];
        if (r->P < 0) { badP++; printf("  P not integer: A=%lld B=%lld ML=%lld/%lld k=%lld\n",(long long)r->A,(long long)r->B,(long long)r->p,(long long)r->q,(long long)r->k); continue; }
        if (r->P % 2 == 0) { badOdd++; printf("  P even: A=%lld B=%lld P=%lld\n",(long long)r->A,(long long)r->B,(long long)r->P); }
        if (r->P % 3 == 0) { bad3++; printf("  3|P: A=%lld B=%lld P=%lld\n",(long long)r->A,(long long)r->B,(long long)r->P); }
        if (r->q != r->Q) { badQ++; printf("  q != Q: A=%lld B=%lld q=%lld Q=%lld\n",(long long)r->A,(long long)r->B,(long long)r->q,(long long)r->Q); }
        int a = (int)(((r->A % 6) + 6) % 6), b = (int)(((r->B % 6) + 6) % 6);
        cnt[a][b][r->k == 3]++;
        double s = r->A ? (double)r->B / (double)r->A : 1e9;
        if (r->k == 3) { if (s < k3min[a][b]) k3min[a][b] = s; if (s > k3max[a][b]) k3max[a][b] = s; }
        else           { if (s < k1min[a][b]) k1min[a][b] = s; if (s > k1max[a][b]) k1max[a][b] = s; }
        /* paper's table in U^7 coordinates */
        int64_t A7 = r->A + r->B, B7 = -r->A; if (A7 < 0) { A7 = -A7; B7 = -B7; }
        int pred = k3_table_u7(A7, B7);
        if (pred && r->k != 3) fp++;
        if (!pred && r->k == 3) fn++;
        int pj = k3_table_jk(r->A, r->B);
        if (pj && r->k != 3) fpj++;
        if (!pj && r->k == 3) fnj++;
    }
    printf("P non-integer: %ld, P even: %ld, 3|P: %ld, q!=Q: %ld\n", badP, badOdd, bad3, badQ);
    printf("paper's k=3 table (U^7 basis): false positives %ld, false negatives %ld\n", fp, fn);
    printf("k=3 table in the Jain-Kravitz basis: false positives %ld, false negatives %ld\n", fpj, fnj);
    printf("classes (A mod 6, B mod 6) in the U_C = U^7(JK) basis, A >= 0:\n");
    for (int a = 0; a < 6; a++) for (int b = 0; b < 6; b++) if (cnt[a][b][0] + cnt[a][b][1]) {
        printf("  (%d,%d): k=1: %5ld  k=3: %5ld", a, b, cnt[a][b][0], cnt[a][b][1]);
        if (cnt[a][b][1]) printf("   k=3 slopes B/A in [%.4f, %.4f]", k3min[a][b], k3max[a][b]);
        if (cnt[a][b][0]) printf("   k=1 slopes in [%.4f, %.4f]", k1min[a][b], k1max[a][b]);
        printf("\n");
    }
    /* max ratio max|w|/Q */
    double best = 0; long bi = -1;
    for (long i = 0; i < nR; i++) {
        rec *r = &R[i]; int64_t w[6], m = 0;
        for (int j = 0; j < 6; j++) { w[j] = r->A*u[j] + r->B*vv[j]; if (w[j] < 0) w[j] = -w[j]; if (w[j] > m) m = w[j]; }
        double ratio = (double)m / (double)r->Q; if (ratio > best) { best = ratio; bi = i; }
    }
    if (bi >= 0) printf("largest max v / Q = %.5f at (A,B)=(%lld,%lld), Q=%lld\n", best, (long long)R[bi].A, (long long)R[bi].B, (long long)R[bi].Q);
    return 0;
}
