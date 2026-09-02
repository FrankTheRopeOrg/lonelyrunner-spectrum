#include <stdio.h>
#include "lrk.h"
int main(){ int64_t V[2][5]={{1,2,3,4,5},{1,3,4,5,9}}; long bad=0, tested=0;
 for(int t=0;t<2;t++) for(int64_t b=2;b<=12;b++) for(int64_t c=1;c<=400;c++){ if(gcd64(b,c)!=1) continue; int skip=0; for(int i=0;i<5;i++) if(b*V[t][i]==c) skip=1; if(skip) continue;
   int64_t v[6]; for(int i=0;i<5;i++) v[i]=b*V[t][i]; v[5]=c; tested++; if(near_tight(v,6,6)){bad++; printf("near-tight: b=%lld c=%lld t=%d\n",(long long)b,(long long)c,t);} }
 printf("pre-jump check: tested %ld, near-tight %ld\n",tested,bad); return 0; }
