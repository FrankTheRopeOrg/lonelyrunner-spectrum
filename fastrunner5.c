/* One-very-fast-runner tori U_A, U_B: speeds {c} u V, V tight quintuple.
 * For every c up to CMAX (c not in V) decide near-tightness exactly and, when
 * near-tight, check k = q - 6p == 1 and report P = q/k.
 * usage: fastrunner CMAX  (runs both quintuples) */
#include <stdio.h>
#include <omp.h>
#include "lrk.h"
int main(int argc, char **argv) {
    int64_t CMAX = argc > 1 ? atoll(argv[1]) : 2500;
    int64_t V[2][4] = {{1,2,3,4},{1,3,4,7}};
    for (int t = 0; t < 2; t++) {
        long nt = 0, bad = 0; int64_t maxc = 0;
        #pragma omp parallel for schedule(dynamic,16) reduction(+:nt,bad)
        for (int64_t c = 1; c <= CMAX; c++) {
            int skip = 0; for (int i = 0; i < 4; i++) if (V[t][i] == c) skip = 1;
            if (skip) continue;
            int64_t v[5]; for (int i = 0; i < 4; i++) v[i] = V[t][i]; v[4] = c;
            if (!near_tight(v, 5, 5)) continue;
            int64_t p, q, Q; ml_exact(v, 5, &p, &q, &Q);
            int64_t k = q - 5 * p;
            nt++;
            if (k != 1 || p * 5 != c || q != c + 1) { bad++;
                #pragma omp critical
                printf("  k != 1: c=%lld ML=%lld/%lld k=%lld\n", (long long)c,(long long)p,(long long)q,(long long)k); }
            #pragma omp critical
            { if (c > maxc) maxc = c; }
        }
        printf("V=(%lld,%lld,%lld,%lld): c <= %lld: near-tight quintuples %ld, with k != 1: %ld, largest near-tight c %lld\n",
            (long long)V[t][0],(long long)V[t][1],(long long)V[t][2],(long long)V[t][3],
            (long long)CMAX, nt, bad, (long long)maxc);
    }
    return 0;
}
