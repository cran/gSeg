/* Exact-order implementation based on of ade4 1.7-22 MSTgraph.
 * Reference: CRAN ade4_1.7-22.tar.gz, src/tests.c, Algorithm 422.
 * GPL (>= 2), same as ade4/gSeg. Preserve origin N, strict update >,
 * last-minimum <= and swap-with-last active order. Excluded edges use
 * infinity so they cannot compete with finite distances at any scale.
 * Optimizations: symmetric row access; fused update/min scan; compact
 * union flags and direct column-major upper-triangle integer edge output.
 */
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <string.h>
#include <stdint.h>
#include <limits.h>
SEXP gseg_mst_exact(SEXP distances, SEXP size, SEXP layers) {
 int n=asInteger(size), kmax=asInteger(layers);
 if((TYPEOF(distances)!=REALSXP && TYPEOF(distances)!=INTSXP) || n<2 || kmax<1 || kmax>=n)
  error("Invalid prototype MST input");
 R_xlen_t nn=(R_xlen_t)n*n;
 if(XLENGTH(distances)!=(R_xlen_t)n*(n-1)/2)error("Invalid dist length");
 int nprotect=0;
 if(TYPEOF(distances)==INTSXP){distances=PROTECT(coerceVector(distances,REALSXP));nprotect++;}
 double *dm=(double*)R_alloc((size_t)nn,sizeof(double));
 unsigned char *used=(unsigned char*)R_alloc((size_t)nn,sizeof(unsigned char));
 memset(used,0,(size_t)nn);
 const double *d=REAL(distances); R_xlen_t at=0;
 for(int j=0;j<n;j++){
  dm[(R_xlen_t)j*n+j]=R_PosInf;
  for(int i=j+1;i<n;i++){
   double v=d[at++];if(!R_FINITE(v))error("NA/NaN/Inf in foreign function call (arg 1)");
   dm[(R_xlen_t)j*n+i]=v;dm[(R_xlen_t)i*n+j]=v;
  }
 }
 double *ui=(double*)R_alloc(n,sizeof(double));
 int *nit=(int*)R_alloc(n,sizeof(int)),*ji=(int*)R_alloc(n,sizeof(int));
 int *a=(int*)R_alloc(n,sizeof(int)),*b=(int*)R_alloc(n,sizeof(int));
 for(int layer=0;layer<kmax;layer++){
  R_CheckUserInterrupt();int np=n-1,kp=n-1,ne=0;
  const double *row=dm+(R_xlen_t)kp*n;
  for(int i=0;i<np;i++){nit[i]=i;ui[i]=row[i];ji[i]=kp;}
  while(np>0){
   row=dm+(R_xlen_t)kp*n;
   double lane_best[8]={R_PosInf,R_PosInf,R_PosInf,R_PosInf,R_PosInf,R_PosInf,R_PosInf,R_PosInf};
   int lane_ix[8]={-1,-1,-1,-1,-1,-1,-1,-1};int i=0;
#define STEP(L) {int at=i+(L);double v=row[nit[at]];if(ui[at]>v){ui[at]=v;ji[at]=kp;}if(ui[at]<=lane_best[L]){lane_best[L]=ui[at];lane_ix[L]=at;}}
   for(;i+7<np;i+=8){STEP(0);STEP(1);STEP(2);STEP(3);STEP(4);STEP(5);STEP(6);STEP(7);}
#undef STEP
   for(;i<np;i++){double v=row[nit[i]];if(ui[i]>v){ui[i]=v;ji[i]=kp;}if(ui[i]<=lane_best[0]){lane_best[0]=ui[i];lane_ix[0]=i;}}
   double best=lane_best[0];int ix=lane_ix[0];
   for(int l=1;l<8;l++)if(lane_best[l]<best||(lane_best[l]==best&&lane_ix[l]>ix)){best=lane_best[l];ix=lane_ix[l];}
   if(ix<0 || !R_FINITE(best))
    error("The graph constructor could not produce `k` complete edge-disjoint MST layers; set a smaller `k` or supply your own graph explicitly.");
   a[ne]=nit[ix];b[ne]=ji[ix];ne++;kp=nit[ix];
   np--;ui[ix]=ui[np];nit[ix]=nit[np];ji[ix]=ji[np];
  }
  for(int i=0;i<ne;i++){
   R_xlen_t ab=(R_xlen_t)a[i]*n+b[i],ba=(R_xlen_t)b[i]*n+a[i];
   used[ab]=used[ba]=1;dm[ab]=dm[ba]=R_PosInf;
  }
 }
 R_xlen_t count=0;
 for(int j=0;j<n;j++)for(int i=0;i<j;i++)count+=used[(R_xlen_t)j*n+i];
 if(count>INT_MAX)error("Too many edges for an R matrix");
 SEXP ans=PROTECT(allocMatrix(INTSXP,(int)count,2));int *edge=INTEGER(ans);R_xlen_t q=0;
 for(int j=0;j<n;j++)for(int i=0;i<j;i++)if(used[(R_xlen_t)j*n+i]){
  edge[q]=i+1;edge[q+count]=j+1;q++;
 }
 UNPROTECT(nprotect+1);return ans;
}
