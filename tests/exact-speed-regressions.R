library(gSeg)
ns<-asNamespace("gSeg");keys<-get(".gseg_edge_keys",ns);adj<-get(".gseg_adjacency",ns)
set.seed(911)
for(n in c(6,5000,94906265,94906266,2147483647))for(type in 1:4){
 E<-cbind(sample.int(n,1000,TRUE),sample.int(n,1000,TRUE));if(type%%2==0)storage.mode(E)<-"double";E<-rbind(E,E[1:50,2:1])
 stopifnot(identical(anyDuplicated(paste(pmin(E[,1],E[,2]),pmax(E[,1],E[,2]),sep=":")),anyDuplicated(keys(E,n))))
}
for(n in c(1,2,7,100))for(type in 1:4){
 E<-cbind(sample.int(n,500,TRUE),sample.int(n,500,TRUE));if(type%%2==0)storage.mode(E)<-"double";if(type>=3)E<-E[1:2,,drop=FALSE]
 reference<-unname(split(as.double(as.vector(t(E[,2:1,drop=FALSE]))),factor(as.vector(t(E)),levels=seq_len(n))))
 stopifnot(identical(reference,adj(E,n)))
}
cat("PASS exact duplicate detection at integer bound and native adjacency order/types/empty nodes\n")
