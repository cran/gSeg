library(gSeg)
ns<-asNamespace("gSeg");fd<-get(".gseg_fast_dist",ns);policy<-get(".gseg_dist_policy",ns)
set.seed(314);rng<-.Random.seed;mode<-policy();stopifnot(identical(rng,.Random.seed))
cases<-0L
for(n in c(2,3,4,5,8,9,17,50))for(d in c(1,2,3,4,7,16,31,40,200,257,1023))for(kind in c("normal","integer","tiny","huge","offset","duplicate","binary_scale","near","zero","extreme")){
 set.seed(53+n+d);x<-matrix(rnorm(n*d),n)
 if(kind=="integer"){x<-round(x*3);storage.mode(x)<-"integer"}
 if(kind=="tiny")x<-x*1e-160
 if(kind=="huge")x<-x*1e160
 if(kind=="offset")x<-1e150+x*1e140
 if(kind=="duplicate")x[2,]<-x[1,]
 if(kind=="binary_scale")x<-matrix(sample(c(-1,1),n*d,TRUE)*2^sample(-1000:1000,n*d,TRUE),n)
 if(kind=="near")x[2,]<-x[1,]+.Machine$double.eps*abs(x[1,])
 if(kind=="zero")x[]<-rep(c(0,-0),length.out=length(x))
 if(kind=="extreme")x[]<-rep(c(.Machine$double.xmax,-.Machine$double.xmax,0),length.out=length(x))
 if(n%%2==0)rownames(x)<-paste0("row_",seq_len(n))
 ref<-stats::dist(x);stopifnot(identical(fd(x),ref,num.eq=FALSE,single.NA=FALSE))
 if(!is.na(mode))for(row in c(FALSE,TRUE)){
  raw<-.Call(get("C_gseg_dist_reference",ns),as.double(x),as.integer(n),as.integer(d),as.logical(mode),row)
  stopifnot(identical(as.double(ref),raw,num.eq=FALSE,single.NA=FALSE))
 }
 cases<-cases+1L
}
# Missing/nonfinite and tiny matrix behavior remains entirely in stats.
for(x in list(matrix(NA_real_,20,2),matrix(c(Inf,0,1),20,3),matrix(1:3,1),matrix(numeric(),0,2)))
 stopifnot(identical(fd(x),stats::dist(x),num.eq=FALSE,single.NA=FALSE))
match_build<-get(".gseg_dist_match_build",ns);known<-get(".gseg_dist_builds",ns)[[1]]
stopifnot(is.na(match_build(NULL)),!is.na(match_build(known)))
bad<-known;bad$stats_md5<-paste0("unverified-",bad$stats_md5);stopifnot(is.na(match_build(bad)))
bad<-known;bad$body_text<-paste0(bad$body_text,"changed");stopifnot(is.na(match_build(bad)))
local({
 state<-get(".gseg_dist_state",ns);previous<-state$mode;on.exit(state$mode<-previous)
 state$mode<-NA_integer_;set.seed(5);x<-matrix(rnorm(50*31),50)
 stopifnot(identical(fd(x),stats::dist(x),num.eq=FALSE,single.NA=FALSE))
})
cat("PASS reference-exact finite distances, full attributes, safe unknown-build fallback and unchanged RNG; cases",cases,"mode",mode,"\n")
