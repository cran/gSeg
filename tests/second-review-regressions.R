library(gSeg)
.gseg_analytic <- getFromNamespace(".gseg_analytic", "gSeg")
.gseg_integrate <- getFromNamespace(".gseg_integrate", "gSeg")
quiet<-function(expr){invisible(capture.output(z<-force(expr)));z}
fails<-character();checks<-0L
check<-function(label,expr){checks<<-checks+1L;z<-tryCatch({force(expr);TRUE},error=function(e)conditionMessage(e));cat(if(isTRUE(z))'PASS' else 'FAIL',label,if(!isTRUE(z))z else '', '\n');if(!isTRUE(z))fails<<-c(fails,label)}
must_error<-function(expr)inherits(try(force(expr),silent=TRUE),'try-error')
set.seed(1);x<-matrix(rnorm(80),20)
check('max alias matches m including analytic p-values',{
 for(fun in list(gseg1,gseg2)){
  a<-quiet(fun(x,statistics='m'));b<-quiet(fun(x,statistics='max'))
  stopifnot(isTRUE(all.equal(a,b,tolerance=0)))
 }
})
check('dist objects cannot silently become observation matrices',{
 stopifnot(must_error(gseg_kmst(x=dist(x))),must_error(quiet(gseg1(dist(x),pval.appr=FALSE))))
})
# Realistic replication weights, independently expanded to observation pairs.
id<-c(1,1,2,2,2,3,4,4);n<-length(id);E<-rbind(c(1,2),c(2,3),c(3,4),c(1,3));mk<-tabulate(id)
pairs<-t(combn(n,2));wa<-wu<-numeric(nrow(pairs))
for(q in seq_len(nrow(pairs))){a<-id[pairs[q,1]];b<-id[pairs[q,2]]
 if(a==b){wa[q]<-2/mk[a];wu[q]<-1}else if(any((E[,1]==a&E[,2]==b)|(E[,1]==b&E[,2]==a))){wa[q]<-1/(mk[a]*mk[b]);wu[q]<-1}
}
counts<-function(A,w){inside<-seq_len(n)%in%A;c(cross=sum(w[xor(inside[pairs[,1]],inside[pairs[,2]])]),left=sum(w[inside[pairs[,1]]&inside[pairs[,2]]]),right=sum(w[!inside[pairs[,1]]&!inside[pairs[,2]]]))}
expected<-function(A,w){t<-length(A);cs<-t(combn(n,t,FUN=function(a)counts(a,w)));colnames(cs)<-c('cross','left','right');obs<-counts(A,w)
 rw<-((n-t-1)*cs[,'left']+(t-1)*cs[,'right'])/(n-2);rd<-cs[,'left']-cs[,'right'];ow<-((n-t-1)*obs['left']+(t-1)*obs['right'])/(n-2);od<-obs['left']-obs['right']
 z<-function(v,a)as.numeric((a-mean(v))/sqrt(mean((v-mean(v))^2)))
 zw<-z(rw,ow);zd<-z(rd,od);c(ori=-z(cs[,'cross'],obs['cross']),weighted=zw,generalized=zw^2+zd^2,max.type=max(zw,abs(zd)))}
a<-quiet(gseg1_discrete(n,E,id,statistics='all',pval.appr=FALSE));b<-quiet(gseg2_discrete(n,E,id,statistics='all',l0=2,l1=n-2,pval.appr=FALSE))
check('repeated single-scan scores match independent weighted-pair enumeration',{
 for(suffix in c('a','u'))for(t in 2:(n-2)){
  w<-if(suffix=='a')wa else wu;ref<-expected(seq_len(t),w)
  val<-c(a$scanZ$ori[[paste0('Zo_',suffix)]][t],a$scanZ$weighted[[paste0('Zw_',suffix)]][t],a$scanZ$generalized[[paste0('S_',suffix)]][t],a$scanZ$max.type[[paste0('M_',suffix)]][t]);stopifnot(max(abs(val-ref))<1e-9)
 }
})
check('repeated interval scores match independent weighted-pair enumeration',{
 for(suffix in c('a','u'))for(i in 1:(n-2))for(j in seq.int(i+2,min(n,i+n-2))){
  w<-if(suffix=='a')wa else wu;ref<-expected(seq.int(i+1,j),w)
  val<-c(b$scanZ$ori[[paste0('Zo_',suffix)]][i,j],b$scanZ$weighted[[paste0('Zw_',suffix)]][i,j],b$scanZ$generalized[[paste0('S_',suffix)]][i,j],b$scanZ$max.type[[paste0('M_',suffix)]][i,j]);stopifnot(max(abs(val-ref))<1e-9)
 }
})
check('repeated interval scans reject undefined component variances',{
 stopifnot(must_error(quiet(gseg2_discrete(10,cbind(1,2:10),1:10,pval.appr=FALSE))))
})
check('zero-column observations are explicitly rejected',{
 stopifnot(must_error(gseg_kmst(x=matrix(numeric(),20,0))))
})
check('p-value output rejects nonfinite or out-of-range results',{
 stopifnot(isTRUE(all.equal(.gseg_analytic(function()list(x=.2)),list(x=.2))))
 for(v in list(NA_real_,NaN,Inf,-.1,1.1,numeric()))stopifnot(must_error(.gseg_analytic(function()list(x=v))))
})
check('non-converged quadrature cannot silently supply a p-value',{
 stopifnot(abs(.gseg_integrate(function(x)x,0,1)$value-.5)<1e-12)
 stopifnot(must_error(.gseg_integrate(function(x)1/sqrt(x),0,1,subdivisions=1,stop.on.error=FALSE)))
})
cat('RESULT:',checks-length(fails),'/',checks,'passed\n');if(length(fails))stop(paste(fails,collapse=', '))
