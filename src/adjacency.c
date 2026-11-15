/* GPL>=2. Build exactly the original scalar-loop adjacency, in input-edge order. */
#include <R.h>
#include <Rinternals.h>
#include <limits.h>
#include <string.h>
SEXP gseg_adjacency(SEXP edges, SEXP size) {
 int n=asInteger(size), type=TYPEOF(edges);
 if(n<1 || (type!=INTSXP && type!=REALSXP) || XLENGTH(edges)%2)error("Invalid adjacency input");
 R_xlen_t m=XLENGTH(edges)/2;
 int *deg=(int*)R_alloc(n,sizeof(int)),*at=(int*)R_alloc(n,sizeof(int));
 memset(deg,0,(size_t)n*sizeof(int));memset(at,0,(size_t)n*sizeof(int));
 for(R_xlen_t i=0;i<m;i++){
  double da=type==INTSXP?INTEGER(edges)[i]:REAL(edges)[i],db=type==INTSXP?INTEGER(edges)[i+m]:REAL(edges)[i+m];
  if(!R_FINITE(da)||!R_FINITE(db)||da<1||db<1||da>n||db>n||da!=(int)da||db!=(int)db)error("Invalid adjacency endpoint");
  int a=(int)da-1,b=(int)db-1;
  if(deg[a]==INT_MAX)error("Adjacency degree too large");
  deg[a]++;
  if(deg[b]==INT_MAX)error("Adjacency degree too large");
  deg[b]++;
 }
 SEXP ans=PROTECT(allocVector(VECSXP,n));
 for(int i=0;i<n;i++)SET_VECTOR_ELT(ans,i,allocVector(REALSXP,deg[i]));
 for(R_xlen_t i=0;i<m;i++){
  if((i&65535)==0)R_CheckUserInterrupt();
  int a=(type==INTSXP?INTEGER(edges)[i]:(int)REAL(edges)[i])-1;
  int b=(type==INTSXP?INTEGER(edges)[i+m]:(int)REAL(edges)[i+m])-1;
  REAL(VECTOR_ELT(ans,a))[at[a]++]=b+1;
  REAL(VECTOR_ELT(ans,b))[at[b]++]=a+1;
 }
 UNPROTECT(1);return ans;
}
