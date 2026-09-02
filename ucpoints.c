/* read "A B" pairs from stdin, print "A B p q" with exact ML of A u + B v on U_C (absolute speeds).
   prints p q = 0 0 if not proper / speeds not distinct */
#include <stdio.h>
#include "lrk.h"
int main(){ int64_t u[6]={1,0,1,2,3,3}, vv[6]={0,1,1,1,1,2}; long long A,B;
 while(scanf("%lld %lld",&A,&B)==2){ int64_t aw[6]; int ok=1;
   for(int i=0;i<6;i++){ int64_t w=A*u[i]+B*vv[i]; aw[i]=w<0?-w:w; }
   for(int i=0;i<6&&ok;i++){ if(aw[i]==0) ok=0; for(int j=0;j<i;j++) if(aw[i]==aw[j]) ok=0; }
   if(!ok){ printf("%lld %lld 0 0\n",A,B); continue; }
   int64_t p,q,Q; ml_exact(aw,6,&p,&q,&Q); printf("%lld %lld %lld %lld\n",A,B,(long long)p,(long long)q); }
 return 0; }
