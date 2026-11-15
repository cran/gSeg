/* GPL (>=2). Finite Euclidean kernels derived from R stats distance.c.
 * Preserve coordinate accumulation order; explicit fused or separate arithmetic.
 * Optional row-major copy changes access only, never the reduction order.
 */
#include <R.h>
#include <Rinternals.h>
#include <math.h>
#if defined(__clang__)
#pragma clang fp contract(off)
#pragma clang fp reassociate(off)
#elif defined(__GNUC__) && __GNUC__ >= 5
#pragma GCC optimize ("no-fast-math", "fp-contract=off")
#endif
#if defined(__clang__) || (defined(__GNUC__) && __GNUC__ >= 5)
/* The TU pragmas forbid contraction/reassociation. Independent pair lanes
 * may vectorize, but every pair retains its exact scalar reduction order. */
static double square_separate(double d) { return d*d; }
#else
static double square_separate(double d) { volatile double v=d*d; return v; }
#endif
#define ACC_SEPARATE(S,D) ((S)+square_separate(D))
#define ACC_FUSED(S,D) fma((D),(D),(S))
#define DEFINE_KERNEL(NAME,ACC,ROW) \
static void NAME(const double *x,double *out,int n,int p) { \
 R_xlen_t at=0; \
 for(int j=0;j<n-1;j++){ \
  R_CheckUserInterrupt();int i=j+1; \
  for(;i+3<n;i+=4){double s0=0,s1=0,s2=0,s3=0; \
   for(int k=0;k<p;k++){ \
    double v=ROW?x[(R_xlen_t)j*p+k]:x[(R_xlen_t)n*k+j]; \
    double d0=(ROW?x[(R_xlen_t)i*p+k]:x[(R_xlen_t)n*k+i])-v; \
    double d1=(ROW?x[(R_xlen_t)(i+1)*p+k]:x[(R_xlen_t)n*k+i+1])-v; \
    double d2=(ROW?x[(R_xlen_t)(i+2)*p+k]:x[(R_xlen_t)n*k+i+2])-v; \
    double d3=(ROW?x[(R_xlen_t)(i+3)*p+k]:x[(R_xlen_t)n*k+i+3])-v; \
    s0=ACC(s0,d0);s1=ACC(s1,d1);s2=ACC(s2,d2);s3=ACC(s3,d3); \
   } \
   out[at++]=sqrt(s0);out[at++]=sqrt(s1);out[at++]=sqrt(s2);out[at++]=sqrt(s3); \
  } \
  for(;i<n;i++){double s=0;for(int k=0;k<p;k++){ \
   double d=(ROW?x[(R_xlen_t)i*p+k]:x[(R_xlen_t)n*k+i])-(ROW?x[(R_xlen_t)j*p+k]:x[(R_xlen_t)n*k+j]);s=ACC(s,d); \
  }out[at++]=sqrt(s);} \
 } \
}
DEFINE_KERNEL(separate_col,ACC_SEPARATE,0)
DEFINE_KERNEL(separate_row,ACC_SEPARATE,1)
DEFINE_KERNEL(fused_col,ACC_FUSED,0)
DEFINE_KERNEL(fused_row,ACC_FUSED,1)
SEXP gseg_dist_reference(SEXP xs,SEXP nr,SEXP nc,SEXP fused,SEXP rowmajor){
 int n=asInteger(nr),p=asInteger(nc),fm=asLogical(fused),row=asLogical(rowmajor);
 if((TYPEOF(xs)!=REALSXP && TYPEOF(xs)!=INTSXP)||n<2||p<1||XLENGTH(xs)!=(R_xlen_t)n*p||fm==NA_LOGICAL||row==NA_LOGICAL)error("Invalid finite distance input");
 int nprotect=0;
 if(TYPEOF(xs)==INTSXP){xs=PROTECT(coerceVector(xs,REALSXP));nprotect++;}
 const double *x=REAL(xs);
 for(R_xlen_t i=0;i<XLENGTH(xs);i++)if(!R_FINITE(x[i])){UNPROTECT(nprotect);return R_NilValue;}
 if(row){double *t=(double*)R_alloc((size_t)n*p,sizeof(double));
  for(int i=0;i<n;i++){if((i&1023)==0)R_CheckUserInterrupt();for(int k=0;k<p;k++)t[(R_xlen_t)i*p+k]=x[(R_xlen_t)n*k+i];}x=t;}
 SEXP ans=PROTECT(allocVector(REALSXP,(R_xlen_t)n*(n-1)/2));
 if(fm){if(row)fused_row(x,REAL(ans),n,p);else fused_col(x,REAL(ans),n,p);}
 else{if(row)separate_row(x,REAL(ans),n,p);else separate_col(x,REAL(ans),n,p);}
 UNPROTECT(nprotect+1);return ans;
}
