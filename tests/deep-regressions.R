library(gSeg)
fail <- character(); checks <- 0L
check <- function(name, code) {
  checks <<- checks+1L
  ans <- tryCatch({force(code); TRUE},error=function(e)conditionMessage(e))
  cat(if(isTRUE(ans)) "PASS" else "FAIL",name,if(!isTRUE(ans)) ans else "","\n")
  if(!isTRUE(ans)) fail <<- c(fail,name)
}
quiet <- function(expr) {invisible(capture.output(ans<-force(expr)));ans}
must_error <- function(expr) inherits(try(force(expr),silent=TRUE),"try-error")
# Independent exact finite-population moments: no package moment functions.
n <- 8L; E <- rbind(cbind(1:7,2:8),c(1,4),c(2,6),c(3,8))
counts <- function(A) {
 inside <- seq_len(n) %in% A
 c(cross=sum(xor(inside[E[,1]],inside[E[,2]])),
   left=sum(inside[E[,1]] & inside[E[,2]]),
   right=sum(!inside[E[,1]] & !inside[E[,2]]))
}
mom <- lapply(2:(n-2),function(t){
 cs <- t(combn(n,t,FUN=counts)); colnames(cs)<-c("cross","left","right"); rw<-((n-t-1)*cs[,"left"]+(t-1)*cs[,"right"])/(n-2)
 rd<-cs[,"left"]-cs[,"right"]
 list(cross=mean(cs[,"cross"]), vcross=mean((cs[,"cross"]-mean(cs[,"cross"]))^2),
      rw=mean(rw),vrw=mean((rw-mean(rw))^2),rd=mean(rd),vrd=mean((rd-mean(rd))^2))
})
reference <- function(A) {
 t<-length(A);c<-counts(A);m<-mom[[t-1]]
 rw<-((n-t-1)*c['left']+(t-1)*c['right'])/(n-2)
 zw<-as.numeric((rw-m$rw)/sqrt(m$vrw));zd<-as.numeric((c['left']-c['right']-m$rd)/sqrt(m$vrd))
 c(ori=as.numeric((m$cross-c['cross'])/sqrt(m$vcross)),weighted=zw,
   generalized=zw^2+zd^2,max.type=max(zw,abs(zd)))
}
z<-quiet(gseg1_graph(n,E,statistics='all',pval.appr=FALSE))
check('single-scan exact subset moments',{
 for(t in 2:(n-2)) {
   ref<-reference(seq_len(t))
   vals<-c(ori=z$scanZ$ori$Z[t],weighted=z$scanZ$weighted$Zw[t],
           generalized=z$scanZ$generalized$S[t],max.type=z$scanZ$max.type$M[t])
   stopifnot(max(abs(vals-ref))<1e-10)
 }
})
z<-quiet(gseg2_graph(n,E,statistics='all',l0=2,l1=n-2,pval.appr=FALSE))
check('interval-scan exact subset moments',{
 for(i in 1:(n-2))for(j in seq.int(i+2,min(n,i+n-2))) {
   ref<-reference(seq.int(i+1,j));vals<-c(ori=z$scanZ$ori$Z[i,j],
     weighted=z$scanZ$weighted$Zw[i,j],generalized=z$scanZ$generalized$S[i,j],max.type=z$scanZ$max.type$M[i,j])
   stopifnot(max(abs(vals-ref))<1e-10)
 }
})
# The domain restriction applies to returned maximizers, including ties.
set.seed(6); x<-matrix(rnorm(40),10);edges<-gseg_kmst(x=x,k=1)
z<-quiet(gseg2_graph(10,edges,statistics='all',l0=2,l1=2,pval.appr=FALSE))
check('continuous interval maximizers stay in requested window',{
 for(s in z$scanZ){t<-s$tauhat;half<-length(t)/2;stopifnot(half>=1,half==floor(half),
   all(t[seq_len(half)]>=1),all(t[half+seq_len(half)]<=10),
   all(t[half+seq_len(half)]-t[seq_len(half)]==2))}
})
id<-rep(1:6,each=2); ed<-rbind(cbind(1:5,2:6),c(1,4),c(2,5))
z<-quiet(gseg2_discrete(12,ed,id,statistics='all',l0=3,l1=3,pval.appr=FALSE))
check('repeated-observation interval maximizers stay in requested window',{
 for(s in z$scanZ)for(suffix in c('a','u')){
   t<-s[[paste0('tauhat_',suffix)]];half<-length(t)/2
   stopifnot(half>=1,half==floor(half),all(t[half+seq_len(half)]-t[seq_len(half)]==3))
 }
})
# Recover every seeded permutation independently by relabeling E.
set.seed(44);z<-quiet(gseg1_graph(10,edges,statistics='g',pval.appr=FALSE,pval.perm=TRUE,B=9))
check('GET returns the actual permutation scan paths',{
 set.seed(44)
 for(b in 1:9){perm<-sample(10);ep<-matrix(perm[edges],ncol=2)
  expected<-quiet(gseg1_graph(10,ep,statistics='g',pval.appr=FALSE))$scanZ$generalized$S
  stopifnot(isTRUE(all.equal(z$pval.perm$generalized$Z[b,],expected,tolerance=0)))
 }
})
set.seed(1);z<-quiet(gseg1_graph(10,edges,statistics='all',pval.appr=FALSE,pval.perm=TRUE,B=9))
check('permutation tail curves count ties consistently',{
 for(s in z$pval.perm){v<-s$curve[,1];p<-(1+vapply(v,function(t)sum(v>=t),numeric(1)))/10
  stopifnot(identical(as.numeric(s$curve[,2]),as.numeric(p)))}
})
check('undefined scan locations cannot become small p-values',{
 stopifnot(must_error(quiet(gseg1_graph(10,cbind(1,2:10),statistics='o',pval.appr=FALSE,pval.perm=TRUE,B=19))))
})
check('zero-width analytic integration is rejected',{
 stopifnot(must_error(quiet(gseg1_graph(10,edges,n0=5,n1=5))),
           must_error(quiet(gseg2_graph(10,edges,l0=5,l1=5))),
           must_error(quiet(gseg1_discrete(12,ed,id,n0=6,n1=6))),
           must_error(quiet(gseg2_discrete(12,ed,id,l0=4,l1=4))))
})
check('NNL zero-distance ties exclude self-links',{
 ee<-nnl(as.dist(matrix(0,3,3)),1)
 stopifnot(nrow(ee)==3,all(ee[,1]!=ee[,2]))
})
check('NNL rejects invalid layer counts explicitly',{
 msg<-tryCatch({nnl(dist(1:4),0);''},error=function(e)conditionMessage(e))
 stopifnot(grepl('K',msg,fixed=TRUE))
})
# Independent Kruskal reference for the default k-MST union, unique weights.
check('default graph equals independent successive-MST union',{
 set.seed(77);x<-matrix(rnorm(60),20); D<-as.matrix(dist(x));pairs<-t(combn(20,2));keep<-rep(TRUE,nrow(pairs));ref<-matrix(integer(),0,2)
 for(layer in 1:floor(sqrt(20))){
  parent<-1:20
  find<-function(i){while(parent[i]!=i)i<-parent[i];i}
  ord<-which(keep)[order(D[pairs[keep,,drop=FALSE]])];taken<-integer()
  for(ix in ord){a<-find(pairs[ix,1]);b<-find(pairs[ix,2]);if(a!=b){parent[b]<-a;taken<-c(taken,ix);if(length(taken)==19)break}}
  stopifnot(length(taken)==19);keep[taken]<-FALSE;ref<-rbind(ref,pairs[taken,,drop=FALSE])
 }
 actual<-gseg_kmst(x=x);key<-function(E)sort(paste(pmin(E[,1],E[,2]),pmax(E[,1],E[,2]),sep=':'))
 stopifnot(identical(key(actual),key(ref)))
})
check('repeated-observation permutation tail curves count ties', {
  for (which in 1:2) {
    set.seed(82)
    fun <- get(paste0("gseg",which,"_discrete"))
    z <- quiet(fun(12,ed,id,statistics="all",pval.appr=FALSE,pval.perm=TRUE,B=7))
    for (s in z$pval.perm) {
      v<-s$curve[,1];prob<-(1+vapply(v,function(t)sum(v>=t),numeric(1)))/8
      stopifnot(identical(as.numeric(s$curve[,2]),as.numeric(prob)))
    }
  }
})
check('NNL exhaustion is explicit instead of recycling used edges', {
  msg<-tryCatch({nnl(as.dist(matrix(0,4,4)),2);''},error=function(e)conditionMessage(e))
  stopifnot(grepl('unused',msg,fixed=TRUE))
})
check('explicit data names cannot silently become graph calls', {
  edges<-cbind(1:9,2:10)
  for (fun in list(gseg1,gseg2)) {
    stopifnot(must_error(quiet(fun(x=edges,n=10,pval.appr=FALSE))),
              must_error(quiet(fun(x=10,E=edges,pval.appr=FALSE))))
  }
})
check('sparse category IDs fail before allocating their maximum', {
  stopifnot(must_error(quiet(gseg1_discrete(6,cbind(1:2,2:3),c(1,2,3,4,5,1e9),pval.appr=FALSE))))
})
check('incomplete MST layers are not reported as complete', {
  # A Euclidean star: centre at the origin, leaves on five coordinate axes.
  x<-rbind(rep(0,5),diag(5))
  stopifnot(nrow(gseg_kmst(x=x,k=1))==5,
            must_error(gseg_kmst(x=x,k=2)),
            must_error(gseg_kmst(dissimilarity=dist(x),k=2)))
})
cat('RESULT:',checks-length(fail),'/',checks,'passed\n')
if(length(fail)) stop(paste("Deep regression failures:", paste(fail, collapse=", ")))
