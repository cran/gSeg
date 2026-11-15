library(gSeg)
ns<-asNamespace("gSeg");fd<-get(".gseg_fast_dist",ns);policy<-get(".gseg_dist_policy",ns)
mode<-policy();if(!is.na(mode))stopifnot(get(".gseg_dist_conversion_ok",ns)())
counter<-new.env();counter$method<-0L;counter$getter<-0L
reference<-function(x)stats::dist(x)
clone_namespace<-function(){
 e<-new.env(parent=ns)
 for(nm in ls(ns,all.names=TRUE)){f<-get(nm,ns);if(is.function(f)&&identical(environment(f),ns)){environment(f)<-e;assign(nm,f,e)}}
 stopifnot(identical(environment(get("gseg1",ns)),ns));e
}
fresh_fd<-function(){e<-clone_namespace();e$.gseg_dist_state<-new.env(parent=emptyenv());e$.gseg_fast_dist}
reference_engine<-clone_namespace();reference_engine$.gseg_fast_dist<-compiler::cmpfun(reference)
observe<-function(f,x){
 set.seed(901);counter$method<-0L;counter$getter<-0L;w<-character();ok<-TRUE
 console<-capture.output(value<-tryCatch(withCallingHandlers(f(x),warning=function(e){w<<-c(w,conditionMessage(e));invokeRestart("muffleWarning")}),error=function(e){ok<<-FALSE;list(error=conditionMessage(e))}))
 list(value=value,ok=ok,console=console,warnings=w,rng=.Random.seed,methods=counter$method,getters=counter$getter)
}
checks<-0L
check<-function(label,x,public=FALSE,allow_reference_error=FALSE){
 for(f in list(fd,fresh_fd())){
  a<-observe(f,x);b<-observe(reference,x);stopifnot(identical(a,b,num.eq=FALSE,single.NA=FALSE));checks<<-checks+1L
 }
 if(public)for(j in 1:2){
  new<-get(paste0("gseg",j),ns);old<-reference_engine[[paste0("gseg",j)]]
  args<-c(list(statistics="m",pval.appr=TRUE,pval.perm=TRUE,B=19),if(j==1)list(n0=3,n1=37)else list(l0=3,l1=37))
  a<-observe(function(x)do.call(new,c(list(x=x),args)),x)
  b<-observe(function(x)do.call(old,c(list(x=x),args)),x)
  stopifnot(identical(a,b,num.eq=FALSE,single.NA=FALSE));if(!allow_reference_error)stopifnot(a$ok);checks<<-checks+1L
 }
 cat("PASS dispatch",label,"\n")
}
with_global<-function(name,value,expr){
 had<-exists(name,.GlobalEnv,inherits=FALSE);if(had){stopifnot(!bindingIsActive(name,.GlobalEnv));previous<-get(name,.GlobalEnv,inherits=FALSE)}
 assign(name,value,.GlobalEnv)
 on.exit(if(had)assign(name,previous,.GlobalEnv)else rm(list=name,envir=.GlobalEnv))
 force(expr)
}
with_registered<-function(class,value,expr){
 table<-get(".__S3MethodsTable__.",asNamespace("base"));name<-paste0("as.matrix.",class)
 had<-exists(name,table,inherits=FALSE);if(had)previous<-get(name,table,inherits=FALSE)
 registerS3method("as.matrix",class,value,envir=asNamespace("base"))
 on.exit(if(had)assign(name,previous,table)else rm(list=name,envir=table))
 force(expr)
}
set.seed(730);plain<-matrix(rnorm(40*20),40)
for(type in c("double","integer")){
 x<-plain;if(type=="integer")storage.mode(x)<-"integer"
 check(paste("plain",type),x)
}
with_global("as.double.gseg_review_scale",function(x,...){counter$method<-counter$method+1L;z<-unclass(x);z[,1]<-20*z[,1];as.vector(z)},{
 x<-plain;class(x)<-c("gseg_review_scale","matrix","array");check("S3 as.double scale",x,TRUE)
})
with_global("as.double.gseg_review_error",function(x,...)stop("unexpected as.double dispatch"),{
 x<-plain;class(x)<-c("gseg_review_error","matrix","array");check("S3 as.double error",x)
})
with_global("as.double.dist",function(x,...){counter$method<-counter$method+1L;runif(1);stop("unexpected dist coercion")},{
 check("as.double.dist ignored, including fresh compiler self-check",plain,TRUE)
})
transform<-function(x,...){counter$method<-counter$method+1L;u<-runif(1);x[,1]<-x[,1]*(1+u);x}
for(cl in c("matrix","array","double","numeric","integer")){
 x<-if(cl=="integer")matrix(seq_len(40*20),40)else plain
 with_global(paste0("as.matrix.",cl),transform,check(paste("visible implicit",cl),x,cl=="matrix",TRUE))
}
with_global("as.matrix.default",function(x,...){counter$method<-counter$method+1L;runif(1);base::as.matrix.default(x)},check("visible default",plain))
with_registered("matrix",transform,check("registered matrix",plain,TRUE,TRUE))
with_registered("default",transform,check("registered default",plain))
with_global("as.matrix.gseg_review_retain",transform,{
 x<-plain;class(x)<-c("gseg_review_retain","matrix","array");check("classed conversion retaining class",x,TRUE)
})
with_global("as.matrix.gseg_review_drop",function(x,...){counter$method<-counter$method+1L;runif(1);unclass(x)},{
 x<-plain;class(x)<-c("gseg_review_drop","matrix","array");check("classed conversion dropping class",x,TRUE)
})
for(name in c("as.matrix.matrix","as.matrix.default"))local({
 stopifnot(!exists(name,.GlobalEnv,inherits=FALSE))
 makeActiveBinding(name,function(value){counter$getter<-counter$getter+1L;transform},.GlobalEnv)
 on.exit(rm(list=name,envir=.GlobalEnv))
 check(paste("active binding",name),plain)
})
methods::setClass("gsegDispatchReviewMatrix",contains="matrix")
x<-methods::new("gsegDispatchReviewMatrix",plain);check("S4 matrix",x)
for(x in list(matrix(seq_len(40*20),40),matrix(c(TRUE,FALSE),40,4),
              matrix(NA_integer_,40,3),matrix(c(Inf,0,1),40,3),
              matrix(complex(real=1:120,imaginary=1),40),
              matrix(as.character(1:120),40),as.data.frame(plain)))
 check(paste("fallback or ALTREP",typeof(x)),x)

# Compare the optimized constructor to its original-conversion bridge,
# including neig class behavior before the canonical integer edge output.
canonical<-function(E){E<-as.matrix(E[,1:2,drop=FALSE]);storage.mode(E)<-"integer";E}
mst_new<-function(x)canonical(get(".gseg_native_mst",ns)(x,6L))
mst_reference<-function(x)canonical(get(".gseg_mst_legacy_conversion",ns)(x,6L))
check_mst<-function(label,d){
 a<-observe(mst_new,d);b<-observe(mst_reference,d)
 stopifnot(identical(a,b,num.eq=FALSE,single.NA=FALSE));checks<<-checks+1L
 cat("PASS constructor",label,"\n")
}
d<-stats::dist(plain)
cat("STANDARD_CONSTRUCTOR_FAST_PATH",get(".gseg_mst_conversion_ok",ns)(d),"\n")
with_global("as.double.dist",function(x,...)stop("unexpected dist coercion"),check_mst("no added as.double.dist",d))
with_global("is.finite.dist",function(x){counter$method<-counter$method+1L;runif(1);rep(TRUE,length(unclass(x)))},{
 check_mst("no added is.finite.dist",d);check("is.finite.dist public path",plain,TRUE)
})
matrix_dist<-get("as.matrix.dist",asNamespace("stats"))
with_global("as.matrix.dist",function(x,...){counter$method<-counter$method+1L;u<-runif(1);matrix_dist(x)*(1+u)},check_mst("custom as.matrix.dist",d))
with_global("t.matrix",function(x){counter$method<-counter$method+1L;runif(1);base::t.default(x)},check_mst("custom transpose",d))
with_global("[.neig",function(x,...){
 counter$method<-counter$method+1L;runif(1)
 stopifnot(!is.null(attr(x,"degrees")),identical(attr(x,"call"),quote(neig(mat01=mst))))
 NextMethod("[")
},check_mst("neig subset dispatch and metadata",d))
bad<-d;attr(bad,"Labels")<-letters[1:3];check_mst("malformed labels preserve error",bad)
sub<-d;class(sub)<-c("gseg_review_dist","dist")
with_global("as.matrix.gseg_review_dist",function(x,...){counter$method<-counter$method+1L;runif(1);matrix_dist(structure(unclass(x),class="dist"))},check_mst("dist subclass conversion",sub))
# Classed graph matrices use the original indexing/adjacency/count order.
E<-rbind(cbind(1:9,2:10),c(1,4),c(2,6),c(3,8));rawE<-E
Eb<-get(".gseg_reference_adjacency",ns)(rawE,10L);deg<-vapply(Eb,length,integer(1))
class(E)<-c("gsegReviewGraph","matrix","array")
with_global("[.gsegReviewGraph",function(x,...){counter$method<-counter$method+1L;runif(1);NextMethod("[")},
 with_global("as.integer.gsegReviewGraph",function(x,...)stop("unexpected graph coercion"),{
  adj_ref<-function(E){out<-vector("list",10);for(i in 1:10)out[[i]]<-rep(0,0);for(i in 1:nrow(E)){out[[E[i,1]]]<-c(out[[E[i,1]]],E[i,2]);out[[E[i,2]]]<-c(out[[E[i,2]]],E[i,1])};out}
  x3_ref<-function(E){v<-0;for(i in 1:nrow(E))v<-v+(deg[E[i,1]]-1)*(deg[E[i,2]]-1);v}
  x5_ref<-function(E){v<-0;for(i in 1:nrow(E)){j<-E[i,1];k<-E[i,2];v<-v+length(which(!is.na(match(Eb[[j]],Eb[[k]]))))};v}
  for(pair in list(list(function(E)get(".gseg_adjacency",ns)(E,10L),adj_ref),list(function(E)get(".gseg_graph_x3",ns)(deg,E),x3_ref),list(function(E)get(".gseg_native_x5",ns)(10L,E,Eb),x5_ref))){
   a<-observe(pair[[1]],E);b<-observe(pair[[2]],E);stopifnot(identical(a,b,num.eq=FALSE,single.NA=FALSE),a$ok);checks<-checks+1L
  }
 }))
cat("PASS classed graph adjacency and moment-count dispatch\n")

cat("CONVERSION_STATE_AFTER_TESTS",get(".gseg_dist_conversion_ok",ns)(),"\n")
# Only the native routine is stressed here, avoiding unrelated R allocation
# work under torture. Integer coercion and transient layout buffers must survive GC.
if(!is.na(mode))local({
 ptr<-get("C_gseg_dist_reference",ns);x<-matrix(seq_len(16*7),16)
 ref<-stats::dist(x);attributes(ref)<-NULL
 previous<-gctorture(TRUE);on.exit(gctorture(previous))
 for(row in c(FALSE,TRUE)){
  value<-.Call(ptr,x,16L,7L,mode,row)
  stopifnot(identical(value,ref,num.eq=FALSE,single.NA=FALSE))
 }
})
cat("PASS dispatch/coercion regression assertions",checks,"mode",mode,"\n")
