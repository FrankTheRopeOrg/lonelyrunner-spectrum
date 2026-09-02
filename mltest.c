/* usage: mltest v1 v2 ... : prints ML p/q and Q, and near-tight flag (1/n) */
#include <stdio.h>
#include "lrk.h"
int main(int argc, char **argv) {
    int64_t v[16]; int n = 0;
    for (int i = 1; i < argc; i++) v[n++] = atoll(argv[i]);
    int64_t p, q, Q; ml_exact(v, n, &p, &q, &Q);
    printf("%lld/%lld Q=%lld near_tight=%d\n", (long long)p, (long long)q, (long long)Q, near_tight(v, n, n));
    return 0;
}
