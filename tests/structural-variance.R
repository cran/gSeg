library(gSeg)
ns<-asNamespace('gSeg');check<-get('.gseg_check_discrete_structure',ns);assertions<-0L
# Exact integer-weight subset enumeration, independent of structural predicates.
verify<-function(id,E){n<-length(id);mk<-tabulate(id);degree<-tabulate(as.integer(E),nbins=length(mk));Wu<-Wa<-matrix(0,n,n);L<-3600
 for(i in 1:(n-1))for(j in (i+1):n){u<-id[i];v<-id[j];linked<-any((E[,1]==u&E[,2]==v)|(E[,1]==v&E[,2]==u));w<-if(u==v)2*L/mk[u]else if(linked)L/(mk[u]*mk[v])else 0;stopifnot(w==floor(w));Wa[i,j]<-Wa[j,i]<-w;Wu[i,j]<-Wu[j,i]<-as.numeric(u==v||linked)}
 union.degree<-vapply(seq_along(mk),function(i)rowSums(Wu)[match(i,id)],numeric(1))
 for(t in c(2L,3L,as.integer(n/2))){zeros<-lapply(list(Wa,Wu),function(W){C<-t(combn(n,t,FUN=function(A){inside<-seq_len(n)%in%A;R1<-sum(W[upper.tri(W)&outer(inside,inside,'&')]);R2<-sum(W[upper.tri(W)&outer(!inside,!inside,'&')]);c(w=(n-t-1)*R1+(t-1)*R2,d=R1-R2,o=sum(W[upper.tri(W)])-R1-R2)}));colnames(C)<-c('w','d','o');apply(C,2,function(z)all(z==z[1]))})
  for(st in c('o','w','g','m','all')){need<-switch(st,o='o',w='w',g=c('w','d'),m=c('w','d'),all=c('w','d','o'));singular<-any(vapply(zeros,function(z)any(z[need]),logical(1)));err<-inherits(try(check(mk,degree,union.degree,st,t),silent=TRUE),'try-error');stopifnot(identical(err,singular));assertions<<-assertions+1L}
 }
}
pairs<-t(combn(4,2))
for(id in list(c(1,2,3,4,1,2,3,4),c(1,2,3,4,3,4,4,4)))for(mask in 1:63)verify(id,pairs[as.logical(intToBits(mask)[1:6]),,drop=FALSE])
pairs<-t(combn(3,2));for(mask in 1:7)verify(c(rep(1,6),2,3),pairs[as.logical(intToBits(mask)[1:3]),,drop=FALSE])
allpairs<-t(combn(8,2));for(E in list(allpairs,matrix(integer(),0,2),cbind(1L,2:8),allpairs[allpairs[,1]!=1,,drop=FALSE],rbind(cbind(1:7,2:8),c(1,4))))verify(1:8,E)
# The reported counterexample must fail before consuming RNG or emitting sqrt warnings.
E<-t(combn(4,2));id<-rep(1:4,3)
for(mode in 1:2)for(st in c('o','w','g','m','all')){ar<-list(n=12,E=E,id=id,statistics=st,pval.appr=FALSE,pval.perm=TRUE,B=19);ar[[if(mode==1)'n0'else'l0']]<-3;ar[[if(mode==1)'n1'else'l1']]<-3;set.seed(44);seed<-.Random.seed;warnings<-0L;msg<-tryCatch(withCallingHandlers(do.call(get(paste0('gseg',mode,'_discrete')),ar),warning=function(w){warnings<<-warnings+1L;invokeRestart('muffleWarning')}),error=function(e)conditionMessage(e));stopifnot(is.character(msg),grepl('exactly zero variance',msg,fixed=TRUE),identical(seed,.Random.seed),warnings==0L)}
# A regular graph invalidates Zd, but must not disable the valid weighted scan.
E<-rbind(cbind(1:5,2:6),c(1,6));id<-rep(1:6,2)
for(mode in 1:2){ar<-list(n=12,E=E,id=id,statistics='w',pval.appr=FALSE,pval.perm=TRUE,B=1);ar[[if(mode==1)'n0'else'l0']]<-3;ar[[if(mode==1)'n1'else'l1']]<-5;invisible(capture.output(z<-suppressWarnings(do.call(get(paste0('gseg',mode,'_discrete')),ar))));stopifnot(is.finite(z$pval.perm$weighted_a$pval),is.finite(z$pval.perm$weighted_u$pval))}
cat('PASS',assertions,'exact enumerated projection-degeneracy checks; 10 early rejections; valid weighted controls\n')
