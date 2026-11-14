library(gSeg)
ns <- asNamespace("gSeg")
root <- get(".gseg_saddlepoint_root",ns)
combine <- get(".gseg_combine_tails",ns)
Nu <- get("Nu",ns)
check <- function(x,y)stopifnot(isTRUE(all.equal(unname(x),unname(y),tolerance=1e-11)))
check(root(4,1),4)
# Verify the saddlepoint equation, independently of the explicit root formula.
r<-c(-.1,-1e-18,0,1e-18,.1);theta<-root(4,1+2*r*4)
check(theta+r*theta^2/2,rep(4,length(r)))
stopifnot(all(is.finite(theta)),root(.Machine$double.xmax,1)==.Machine$double.xmax)
check(combine(1e-20,2e-20),3e-20)
stopifnot(combine(1e-20,2e-20)>0,combine(0,0)==0,combine(1,.2)==1,combine(.2,1)==1,combine(2,.2)==1)
check(combine(.2,.3),.44)
# Exact constant-integrand zero-skew limits in all eight correction helpers.
n<-20;b<-4;lower<-2;upper<-18;xx<-.1
for(kind in 1:2)for(discrete in c(FALSE,TRUE))for(component in 1:2){
 name<-paste0("pval",kind,if(discrete)"_discrete" else "","_sub_",component)
 f<-get(name,ns);len<-if(discrete)n else n-1;factor<-if(component==1)2 else 1
 expected<-if(kind==1)factor*dnorm(b)*b*xx*Nu(sqrt(2*b*b*xx))*(upper-lower) else factor*dnorm(b)/b*(b*b*xx*Nu(sqrt(2*b*b*xx)))^2*(n*(upper-lower)-(upper^2-lower^2)/2)
 for(r0 in c(-1e-18,0,1e-18))check(f(n,b,rep(r0,len),rep(xx,len),lower,upper),expected)
}
# Public-interface regression for positive, nontruncated strong-signal tails.
set.seed(901);x<-rbind(matrix(rnorm(60),20,3),matrix(rnorm(60)+3,20,3))
for(f in list(gseg1,gseg2)){
 invisible(capture.output(ans<-f(x)))
 stopifnot(is.finite(ans$pval.appr$max.type),ans$pval.appr$max.type>0,ans$pval.appr$max.type<1e-10)
}
cat("PASS saddlepoint equation; all eight zero/near-zero helper limits; tiny-tail combination; public MET regression\n")
