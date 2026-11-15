library(gSeg)
ns<-asNamespace('gSeg');calls<-0L
with_method<-function(class,fun,expr){tab<-get('.__S3MethodsTable__.',asNamespace('base'));key<-paste0('as.matrix.',class);had<-exists(key,tab,inherits=FALSE);if(had)old<-get(key,tab,inherits=FALSE);registerS3method('as.matrix',class,fun,envir=asNamespace('base'));on.exit(if(had)registerS3method('as.matrix',class,old,envir=asNamespace('base'))else rm(list=key,envir=tab));force(expr)}
observe<-function(fun){set.seed(44);calls<<-0L;warnings<-character();console<-capture.output(z<-withCallingHandlers(fun(),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}));list(value=z,rng=.Random.seed,calls=calls,warnings=warnings,console=console)}
set.seed(221);x<-matrix(rnorm(24*3),24);d<-dist(x);matrix.dist<-get('as.matrix.dist',asNamespace('stats'))
with_method('dist',function(x,...){calls<<-calls+1L;runif(1);matrix.dist(x,...)},{
 # Dispatch alone must not evaluate any matrix method.
 z<-get('.gseg_graph_call',ns)(NULL,d,NULL,NULL,1,c('dissimilarity','k'));stopifnot(identical(z,list(graph=FALSE)),calls==0L)
 for(mode in 1:2){fn<-get(paste0('gseg',mode));fg<-get(paste0('gseg',mode,'_graph'));ar<-list(statistics='m',pval.appr=FALSE,pval.perm=TRUE,B=19)
  a<-observe(function()do.call(fn,c(list(dissimilarity=d,k=1),ar)))
  b<-observe(function(){E<-gseg_kmst(dissimilarity=d,k=1);z<-do.call(fg,c(list(n=24,E=E),ar));z$graph<-E;z$k<-attr(E,'k');z})
  stopifnot(a$calls==1L,identical(a,b,num.eq=FALSE,single.NA=FALSE))
 }
})
# Stateful row conversion: metadata and graph use the same conversion.
xx<-x;class(xx)<-c('gsegNormalizeRows','matrix','array')
with_method('gsegNormalizeRows',function(x,...){calls<<-calls+1L;z<-unclass(x);if(calls%%2==1L)z<-z[-nrow(z),,drop=FALSE];z},{
 for(mode in 1:2){fn<-get(paste0('gseg',mode));calls<-0L;E<-gseg_kmst(xx,k=1);stopifnot(calls==1L,max(E)==23L)
  a<-observe(function()fn(xx,k=1,pval.appr=FALSE));stopifnot(a$calls==1L,identical(a$value$graph,E),length(a$value$scanZ$max.type[[if(mode==1)'M'else'Mv']])==if(mode==1)23 else 23^2)
 }
})
# Retained classes still get the legitimate stats::dist conversion. Use its
# actual Size, rather than the preliminary data matrix row count.
with_method('gsegNormalizeRows',function(x,...){calls<<-calls+1L;z<-unclass(x)[-nrow(x),,drop=FALSE];class(z)<-class(x);z},{
 for(mode in 1:2){a<-observe(function()get(paste0('gseg',mode))(xx,k=1,pval.appr=FALSE));stopifnot(a$calls==2L,max(a$value$graph)==22L)}
})
# Matrix dissimilarity conversion also happens once.
DD<-as.matrix(d);class(DD)<-c('gsegNormalizeDistance','matrix','array')
with_method('gsegNormalizeDistance',function(x,...){calls<<-calls+1L;z<-unclass(x);if(calls%%2==1L)z<-z[-nrow(z),-ncol(z),drop=FALSE];z},{
 for(mode in 1:2){a<-observe(function()get(paste0('gseg',mode))(dissimilarity=DD,k=1,pval.appr=FALSE));stopifnot(a$calls==1L,max(a$value$graph)==23L)}
})
# Old positional and mixed named graph calls normalize E exactly once, and
# agree with the explicit graph API, including user-method RNG side effects.
E<-gseg_kmst(x,k=1);class(E)<-c('gsegNormalizeEdges','matrix','array')
with_method('gsegNormalizeEdges',function(x,...){calls<<-calls+1L;runif(1);unclass(x)},{
 for(mode in 1:2){fg<-get(paste0('gseg',mode,'_graph'));fn<-get(paste0('gseg',mode));ar<-list(statistics='all',pval.appr=FALSE,pval.perm=TRUE,B=3);b<-observe(function()do.call(fg,c(list(n=24,E=E),ar)))
  for(input in list(list(24,E),list(24,E=E),list(E,n=24),list(n=24,E=E))){a<-observe(function()do.call(fn,c(input,ar)));stopifnot(a$calls==1L,identical(a,b,num.eq=FALSE,single.NA=FALSE))}
 }
})
# Standard dataframes retain constructor metadata and scan results.
for(mode in 1:2){a<-observe(function()get(paste0('gseg',mode))(x,k=1,pval.appr=FALSE));b<-observe(function()get(paste0('gseg',mode))(as.data.frame(x),k=1,pval.appr=FALSE));stopifnot(identical(a$value,b$value,num.eq=FALSE,single.NA=FALSE))}
# Equivalent raw, matrix and dist forms reject scans below N=6.
for(mode in 1:2)for(n in c(2L,3L,5L)){z<-matrix(seq_len(n*2),n);for(input in list(list(x=z),list(dissimilarity=as.matrix(dist(z))),list(dissimilarity=dist(z)))){msg<-tryCatch({do.call(get(paste0('gseg',mode)),c(input,list(k=1,pval.appr=FALSE)));NULL},error=function(e)conditionMessage(e));stopifnot(is.character(msg))}}
cat('PASS dispatch without coercion; registered-method/RNG parity; stateful and retained-class rows; distance matrices; legacy calls; dataframe and minimum-N controls\n')
