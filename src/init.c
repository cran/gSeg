#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
extern SEXP gseg_mst_exact(SEXP,SEXP,SEXP);
extern SEXP gseg_graph_x5(SEXP,SEXP);
extern SEXP gseg_dist_exact(SEXP,SEXP,SEXP);
static const R_CallMethodDef entries[]={
 {"gseg_mst_exact",(DL_FUNC)&gseg_mst_exact,3},
 {"gseg_graph_x5",(DL_FUNC)&gseg_graph_x5,2},
 {"gseg_dist_exact",(DL_FUNC)&gseg_dist_exact,3},
 {NULL,NULL,0}};
void R_init_gSeg(DllInfo*dll){R_registerRoutines(dll,NULL,entries,NULL,NULL);R_useDynamicSymbols(dll,FALSE);R_forceSymbols(dll,TRUE);}
