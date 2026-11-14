/* GPL>=2. Finite-matrix Euclidean distance, based on R4.4.0 distance.c.
 * Four independent pairs share coordinate loads; per-pair dimension order
 * and subtraction/multiplication/addition/sqrt are unchanged. Compile with
 * -ffp-contract=off, without fast-math. Not a general missing-data dist API.
 */
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <math.h>
/* Strict arithmetic without non-portable Makevars flags. */
#if defined(__clang__)
#pragma clang fp contract(off)
#define GSEG_SQUARE(v) ((v)*(v))
#elif defined(__GNUC__) && __GNUC__ >= 5
#pragma GCC optimize ("fp-contract=off")
#define GSEG_SQUARE(v) ((v)*(v))
#else
static double gseg_separate_square(double v) { volatile double p=v*v; return p; }
#define GSEG_SQUARE(v) gseg_separate_square(v)
#endif
SEXP gseg_dist_exact(SEXP xs,SEXP nr,SEXP nc){
 int n=asInteger(nr),p=asInteger(nc);
 if(TYPEOF(xs)!=REALSXP||n<2||p<1||XLENGTH(xs)!=(R_xlen_t)n*p)error("Invalid finite distance input");
 const double*x=REAL(xs);for(R_xlen_t i=0;i<XLENGTH(xs);i++)if(!R_FINITE(x[i]))error("Distance fast path requires finite data");
 R_xlen_t m=(R_xlen_t)n*(n-1)/2;SEXP ans=PROTECT(allocVector(REALSXP,m));double*out=REAL(ans);R_xlen_t at=0;
 for(int j=0;j<n-1;j++){
  R_CheckUserInterrupt();int i=j+1;
  for(;i+3<n;i+=4){
   double s0=0,s1=0,s2=0,s3=0;
   for(int k=0;k<p;k++){
    R_xlen_t c=(R_xlen_t)n*k;double v=x[c+j];
    double d0=x[c+i]-v,d1=x[c+i+1]-v,d2=x[c+i+2]-v,d3=x[c+i+3]-v;
    s0+=GSEG_SQUARE(d0);s1+=GSEG_SQUARE(d1);s2+=GSEG_SQUARE(d2);s3+=GSEG_SQUARE(d3);
   }
   out[at++]=sqrt(s0);out[at++]=sqrt(s1);out[at++]=sqrt(s2);out[at++]=sqrt(s3);
  }
  for(;i<n;i++){
   double s=0;for(int k=0;k<p;k++){R_xlen_t c=(R_xlen_t)n*k;double v=x[c+i]-x[c+j];s+=GSEG_SQUARE(v);}out[at++]=sqrt(s);
  }
 }
 UNPROTECT(1);return ans;
}
