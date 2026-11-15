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
SEXP gseg_mst_matrix(SEXP distances, SEXP size, SEXP layers) {
 int n=asInteger(size), kmax=asInteger(layers);
 if(TYPEOF(distances)!=REALSXP || n<2 || kmax<1 || kmax>=n)
  error("Invalid prototype MST input");
 R_xlen_t nn=(R_xlen_t)n*n;
 if(XLENGTH(distances)!=nn)error("Invalid matrix distance length");
 double *dm=(double*)R_alloc((size_t)nn,sizeof(double));
 unsigned char *used=(unsigned char*)R_alloc((size_t)nn,sizeof(unsigned char));
 memset(used,0,(size_t)nn);
 const double *d=REAL(distances);
 /* ade4 copies R column-major values to its transposed DM, then reads
  * DM[NI][KP]. This row-major view therefore stores R[KP,NI]. */
 for(int a=0;a<n;a++)for(int b=0;b<n;b++){
  double v=d[(R_xlen_t)b*n+a];
  if(!R_FINITE(v))error("NA/NaN/Inf in foreign function call (arg 1)");
  dm[(R_xlen_t)a*n+b]=(a==b)?R_PosInf:v;
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
 SEXP ans=PROTECT(allocMatrix(REALSXP,n,n));double *out=REAL(ans);
 for(R_xlen_t i=0;i<nn;i++)out[i]=used[i]?1.0:0.0;
 UNPROTECT(1);return ans;
}
