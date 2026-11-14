library(gSeg)
quiet<-function(expr){invisible(capture.output(z<-force(expr)));z}
# This fails with recursive DFS even on a small chain, without heavy CRAN computation.
with_limited_stack<-function(){
 saved<-options(expressions=100);on.exit(options(saved))
 n<-220L;E<-nnl(as.matrix(dist(seq_len(n))),K=1)
 stopifnot(identical(unname(E),unname(cbind(as.numeric(1:(n-1)),as.numeric(2:n)))))
}
with_limited_stack()
# Normal and large finite scales must have identical ordered edge unions.
set.seed(1);D<-as.matrix(dist(matrix(rnorm(80),20)))
for(k in c(1,2,4)){
 ref<-gseg_kmst(dissimilarity=D,k=k)
 for(scale in c(2^70,2^900))stopifnot(identical(ref,gseg_kmst(dissimilarity=D*scale,k=k)))
}
# Independent exhaustive enumeration of all 5-node labelled trees via Pruefer codes.
n<-5L;D<-outer(1:n,1:n,function(i,j)abs(i-j)+.01*(i+j));diag(D)<-0
codes<-as.matrix(expand.grid(rep(list(1:n),n-2)))
decode<-function(code){
 deg<-1L+tabulate(code,nbins=n);E<-matrix(0L,n-1,2)
 for(j in seq_along(code)){v<-code[j];u<-which(deg==1L)[1];E[j,]<-c(u,v);deg[c(u,v)]<-deg[c(u,v)]-1L}
 E[n-1,]<-which(deg==1L);E
}
keys<-function(E)paste(pmin(E[,1],E[,2]),pmax(E[,1],E[,2]),sep=":")
trees<-lapply(seq_len(nrow(codes)),function(i)decode(codes[i,]))
used<-character()
for(layer in 1:2){
 costs<-vapply(trees,function(E)if(any(keys(E)%in%used))Inf else sum(D[E]),numeric(1))
 stopifnot(is.finite(min(costs)));chosen<-trees[[which.min(costs)]];used<-c(used,keys(chosen))
}
for(scale in c(1,2^70))stopifnot(setequal(keys(gseg_kmst(dissimilarity=D*scale,k=2)),used))
# A unique star consumes all incident edges of its centre: another full layer is impossible.
S<-matrix(2,5,5);diag(S)<-0;S[5,1:4]<-S[1:4,5]<-1
stopifnot(inherits(try(gseg_kmst(dissimilarity=S,k=2),silent=TRUE),"try-error"))
# Nonpositive maxima return 1 independently for average and union components.
id<-rep(1:10,each=2);E<-as.matrix(expand.grid(1:5,6:10))
for(stat in c("o","w"))for(skew in c(FALSE,TRUE)){
 z<-quiet(gseg1_discrete(20,E,id,statistics=stat,n0=9,n1=11,skew.corr=skew))
 stopifnot(all(unlist(z$pval.appr)==1))
}
# Interval case has positive average maximum and negative union maximum.
E<-cbind(1:5,6:10)
for(stat in c("o","w"))for(skew in c(FALSE,TRUE)){
 z<-quiet(gseg2_discrete(20,E,id,statistics=stat,l0=9,l1=11,skew.corr=skew))
 key<-if(stat=="o")"ori" else "weighted"
 stopifnot(z$pval.appr[[paste0(key,"_u")]]==1,is.finite(z$pval.appr[[paste0(key,"_a")]]),z$pval.appr[[paste0(key,"_a")]]>=0,z$pval.appr[[paste0(key,"_a")]]<=1)
}
# Exercise exact zero as well as negative boundaries at both analytic entry points.
for(mode in 1:2)for(stat in c("o","w"))for(skew in c(FALSE,TRUE)){
 raw<-quiet(if(mode==1)gseg1_discrete(20,E,id,statistics=stat,n0=9,n1=11,pval.appr=FALSE) else gseg2_discrete(20,E,id,statistics=stat,l0=9,l1=11,pval.appr=FALSE))
 key<-if(stat=="o")"ori" else "weighted";field<-if(stat=="o")"Zo" else "Zw"
 raw$scanZ[[key]][[paste0(field,"_a_max")]]<-0;raw$scanZ[[key]][[paste0(field,"_u_max")]]<--1
 f<-getFromNamespace(paste0("pval",mode,"_discrete"),"gSeg")
 z<-quiet(f(20,E,id,raw$scanZ,statistics=stat,skew.corr=skew,lower=9,upper=11))
 stopifnot(all(unlist(z)==1))
}
cat("PASS iterative NNL, scale-invariant MST, exhaustive tree oracle, exhausted-layer rejection, repeated OET/WET boundary and mixed-component cases\n")

# Equivalent raw, matrix and dist inputs enforce the same minimum N.
D2 <- matrix(c(0,1,1,0),2)
for (args in list(list(x=matrix(c(0,1),2)),list(dissimilarity=D2),list(dissimilarity=as.dist(D2)))) {
  stopifnot(inherits(try(do.call(gseg_kmst,args),silent=TRUE),"try-error"))
}
x3 <- matrix(c(0,1,3),3);D3 <- as.matrix(dist(x3))
stopifnot(identical(gseg_kmst(x=x3),gseg_kmst(dissimilarity=D3)),
          identical(gseg_kmst(x=x3),gseg_kmst(dissimilarity=as.dist(D3))))
cat("PASS consistent three-observation minimum for raw, matrix and dist inputs\n")
