/* GPL>=2. Exact x5=3*triangle count, using degree-oriented sparse adjacency.
 * Replaces only existing skew-correction graph count; no new statistic.
 */
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <stdint.h>
#include <limits.h>
#include <string.h>
SEXP gseg_graph_x5(SEXP ns,SEXP es){
 int n=asInteger(ns);if(n<2||TYPEOF(es)!=INTSXP||XLENGTH(es)%2)error("Invalid graph-count input");
 R_xlen_t ml=XLENGTH(es)/2;if(ml>INT_MAX)return R_NilValue;int m=(int)ml;
 const int*e=INTEGER(es);
 int*deg=(int*)R_alloc(n,sizeof(int));memset(deg,0,(size_t)n*sizeof(int));
 for(int i=0;i<m;i++){int a=e[i]-1,b=e[i+m]-1;if(a<0||a>=n||b<0||b>=n||a==b)error("Invalid graph endpoint");deg[a]++;deg[b]++;}
 int*offset=(int*)R_alloc((size_t)n+1,sizeof(int));memset(offset,0,((size_t)n+1)*sizeof(int));
 for(int i=0;i<m;i++){int a=e[i]-1,b=e[i+m]-1;int u=(deg[a]<deg[b]||(deg[a]==deg[b]&&a<b))?a:b;offset[u+1]++;}
 for(int u=0;u<n;u++)offset[u+1]+=offset[u];
 int*cursor=(int*)R_alloc(n,sizeof(int));memcpy(cursor,offset,(size_t)n*sizeof(int));
 int*to=(int*)R_alloc(m,sizeof(int));
 for(int i=0;i<m;i++){int a=e[i]-1,b=e[i+m]-1;if(deg[a]>deg[b]||(deg[a]==deg[b]&&a>b)){int t=a;a=b;b=t;}to[cursor[a]++]=b;}
 int*mark=(int*)R_alloc(n,sizeof(int));memset(mark,0,(size_t)n*sizeof(int));uint64_t tri=0;
 for(int u=0;u<n;u++){
  if((u&255)==0)R_CheckUserInterrupt();
  for(int j=offset[u];j<offset[u+1];j++)mark[to[j]]=u+1;
  for(int j=offset[u];j<offset[u+1];j++){
   int v=to[j];for(int k=offset[v];k<offset[v+1];k++)if(mark[to[k]]==u+1){if(tri>=3002399751580330ULL)return R_NilValue;tri++;}
  }
 }
 return ScalarReal(3.0*(double)tri);
}
