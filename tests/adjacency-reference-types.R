library(gSeg)
ns<-asNamespace('gSeg');adj<-get('.gseg_adjacency',ns)
original<-function(E,n){out<-vector('list',n);for(i in 1:n)out[[i]]<-rep(0,0);for(i in 1:nrow(E)){out[[E[i,1]]]<-c(out[[E[i,1]]],E[i,2]);out[[E[i,2]]]<-c(out[[E[i,2]]],E[i,1])};out}
checks<-0L
for(n in c(3L,8L,19L,40L))for(style in c('path','complete','disconnected','irregular')){
 pairs<-t(combn(n,2));E<-switch(style,path=cbind(1:(n-1),2:n),complete=pairs,disconnected=pairs[pairs[,2]<=ceiling(n/2),,drop=FALSE],irregular=pairs[((seq_len(nrow(pairs))*7)%%11)<4,,drop=FALSE]);if(!nrow(E))next
 E<-E[rev(seq_len(nrow(E))),2:1,drop=FALSE]
 for(type in c('integer','double')){storage.mode(E)<-type;a<-adj(E,n);b<-original(E,n);stopifnot(identical(a,b,num.eq=FALSE,single.NA=FALSE),all(vapply(a,is.double,logical(1))));checks<-checks+1L}
}
E<-matrix(c(1L,2L),1,2);a<-adj(E,6L);stopifnot(identical(a,original(E,6L)),identical(a[[6]],numeric(0)))
local({old<-gctorture(TRUE);on.exit(gctorture(old));stopifnot(identical(adj(E,6L),original(E,6L)))})
cat('PASS',checks,'strict original-loop adjacency cases; isolated vertices and GC protection\n')
