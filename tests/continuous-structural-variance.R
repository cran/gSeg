library(gSeg);rows<-list()
observe<-function(fun,ar){warnings<-character();ok<-TRUE;console<-capture.output(z<-tryCatch(withCallingHandlers(do.call(fun,ar),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e){ok<<-FALSE;conditionMessage(e)}));list(ok=ok,z=z,warnings=warnings)}
for(n in c(6L,8L,9L,10L,12L,15L,20L))for(shape in c('complete','cycle','star','complement_star')){
 pairs<-t(combn(n,2));E<-switch(shape,complete=pairs,cycle=rbind(cbind(1:(n-1),2:n),c(1,n)),star=cbind(1,2:n),complement_star=pairs[pairs[,1]!=1,,drop=FALSE])
 for(mode in 1:2)for(t in unique(c(2L,3L,as.integer(n%/%2))))for(st in c('o','w','g','m')){
  zero.w<-shape!='cycle';zero.d<-shape%in%c('complete','cycle');singular<-switch(st,o=zero.w&&(zero.d||2*t==n),w=zero.w,g=zero.w||zero.d,m=zero.w||zero.d)
  ar<-list(n=n,E=E,statistics=st,pval.appr=FALSE,pval.perm=FALSE);ar[[if(mode==1)'n0'else'l0']]<-t;ar[[if(mode==1)'n1'else'l1']]<-t;v<-observe(get(paste0('gseg',mode,'_graph')),ar)
  value<-if(v$ok)v$z$scanZ[[switch(st,o='ori',w='weighted',g='generalized',m='max.type')]]$Zmax else NA_real_
  rows[[length(rows)+1L]]<-data.frame(n=n,shape=shape,mode=mode,t=t,statistic=st,exactly_singular=singular,accepted=v$ok,max=value,warnings=length(v$warnings),error=if(v$ok)''else v$z)
 }
}
z<-do.call(rbind,rows)
stopifnot(nrow(z)==640,sum(z$exactly_singular)==500,sum(!z$exactly_singular)==140,
          all(!z$accepted[z$exactly_singular]),all(z$accepted[!z$exactly_singular]))
# Unique-ID discrete and graph entry points agree on validity; rejection must
# precede sampling and warning-producing normalizations, with graph wording.
E<-rbind(cbind(1:7,2:8),c(1,8));n<-8L
for(mode in 1:2)for(kind in c('graph','alias','discrete')){
 ar<-list(n=n,E=E,statistics='m',pval.appr=FALSE,pval.perm=TRUE,B=19)
 ar[[if(mode==1)'n0'else'l0']]<-2;ar[[if(mode==1)'n1'else'l1']]<-2
 if(kind=='discrete')ar$id<-1:n
 fn<-get(paste0('gseg',mode,switch(kind,graph='_graph',alias='',discrete='_discrete')))
 set.seed(44);seed<-.Random.seed;v<-observe(fn,ar)
 stopifnot(!v$ok,grepl('exactly zero variance',v$z,fixed=TRUE),!length(v$warnings),identical(seed,.Random.seed))
 if(kind!='discrete')stopifnot(grepl('graph channel',v$z,fixed=TRUE))
}
# Weighted-only cycle scans remain valid with both permutation and (n>=9)
# analytic calibration; a singular difference component is not requested.
n<-12;E<-rbind(cbind(1:11,2:12),c(1,12))
for(mode in 1:2){ar<-list(n=n,E=E,statistics='w',pval.appr=FALSE,pval.perm=TRUE,B=3);ar[[if(mode==1)'n0'else'l0']]<-3;ar[[if(mode==1)'n1'else'l1']]<-5;stopifnot(observe(get(paste0('gseg',mode,'_graph')),ar)$ok)}
cat('PASS 640 graph boundary cases: 500 singular rejections, 140 valid controls; unique-ID/alias/RNG rejection controls\n')
