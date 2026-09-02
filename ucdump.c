/* dump exact ML for all proper six-distinct-speed subtori of U_C with 0<=A<=H, |B|<=H: "A B p q Q" */
#include <stdio.h>
#include <omp.h>
#include "lrk.h"
int main(int argc,char**argv){ int64_t H=atoll(argv[1]); int64_t u[6]={1,0,1,2,3,3}, vv[6]={0,1,1,1,1,2};
 #pragma omp parallel for schedule(dynamic,4)
 for(int64_t A=0;A<=H;A++){ for(int64_t B=-H;B<=H;B++){ if(A==0&&B!=1) continue; if(gcd64(A,B)!=1) continue;
   int64_t w[6],aw[6]; for(int i=0;i<6;i++){w[i]=A*u[i]+B*vv[i]; aw[i]=w[i]<0?-w[i]:w[i];}
   int ok=1; for(int i=0;i<6&&ok;i++){ if(aw[i]==0) ok=0; for(int j=0;j<i;j++) if(aw[i]==aw[j]) ok=0;} if(!ok) continue;
   int64_t p,q,Q; ml_exact(aw,6,&p,&q,&Q);
   #pragma omp critical
   printf("%lld %lld %lld %lld %lld\n",(long long)A,(long long)B,(long long)p,(long long)q,(long long)Q); } }
 return 0; }
