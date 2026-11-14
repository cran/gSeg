library(gSeg)
fd<-getFromNamespace(".gseg_fast_dist","gSeg");x5<-getFromNamespace(".gseg_native_x5","gSeg");x3<-getFromNamespace(".gseg_graph_x3","gSeg")
bits<-function(z)writeBin(as.double(z),raw(),size=8)
ncases<-0L
for(n in c(2:9,17,50))for(d in c(1,2,3,4,5,7,16,31,200))for(kind in c("random","integer","tiny","huge","offset","duplicate")){
 set.seed(53+n+d);x<-matrix(rnorm(n*d),n)
 if(kind=="integer")x<-round(x*3)
 if(kind=="tiny")x<-x*1e-160
 if(kind=="huge")x<-x*1e160
 if(kind=="offset")x<-1e150+x*1e140
 if(kind=="duplicate")x[2,]<-x[1,]
 stopifnot(identical(bits(stats::dist(x)),bits(fd(x))));ncases<-ncases+1L
}
cat("PASS",ncases,"bitwise distance cases\n")
count<-0L
for(n in c(3,8,20,50))for(seed in 1:10)for(prob in c(0,.2,.7,1)){
 set.seed(53+n+seed);E<-t(combn(n,2));E<-E[runif(nrow(E))<prob,,drop=FALSE]
 A<-matrix(FALSE,n,n);if(nrow(E)){A[E]<-TRUE;A[cbind(E[,2],E[,1])]<-TRUE}
 t<-combn(n,3);tri<-sum(A[t(t[c(1,2),,drop=FALSE])]&A[t(t[c(1,3),,drop=FALSE])]&A[t(t[c(2,3),,drop=FALSE])])
 stopifnot(identical(x5(n,E),3*as.numeric(tri)));count<-count+1L
}
cat("PASS",count,"independent graph-count cases\n")
# Verify scalar-order fallback above the exact-in-double bound.
deg<-c(1e8,1e8,2);E<-rbind(c(1,2),c(1,3),c(2,3));ref<-0
for(i in seq_len(nrow(E)))ref<-ref+(deg[E[i,1]]-1)*(deg[E[i,2]]-1)
stopifnot(identical(x3(deg,E),ref))
cat("PASS large-count scalar-order safeguard\n")
set.seed(42);x<-matrix(rnorm(200),40);rownames(x)<-paste0("id",1:40)
stopifnot(identical(stats::dist(x),fd(x)))
cat("PASS distance metadata\nNATIVE CHECKS PASSED\n")
