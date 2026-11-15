library(gSeg)
ns<-asNamespace("gSeg");root<-get(".gseg_saddlepoint_root",ns);combine<-get(".gseg_combine_tails",ns)
# Reliable original arithmetic is preserved exactly, not merely within tolerance.
for(b in c(1,4,12))for(r in c(.1,.25,1)){
 disc<-1+2*r*b;legacy<-(sqrt(disc)-1)/r
 stopifnot(identical(root(b,disc,r),legacy))
}
for(p in list(c(.2,.3),c(.5,.5),c(.8,.1),c(0,0),c(1,.2))){
 old<-1-(1-min(p[1],1))*(1-min(p[2],1))
 stopifnot(identical(combine(p[1],p[2]),old))
}
# Counterexamples require an actual repair, checked against defining equations.
for(r in c(-1e-18,0,1e-18)){
 theta<-root(4,1+2*r*4,r)
 stopifnot(identical(theta,4),abs(theta+r*theta^2/2-4)<=.Machine$double.eps*4)
}
stopifnot(combine(1e-20,2e-20)>0,
 abs(combine(1e-12,2e-12)-(3e-12-2e-24))<=1e-27)
# Keep the original double return type while correcting boundary decoding.
f<-get(".gseg_interval_argmax",ns)
stopifnot(identical(f(c(0,1,rep(0,7)),2L,3L),c(1,2)),
 identical(f(c(rep(0,5),1,rep(0,3)),6L,3L),c(2,3)))
# Exact original public distances, including rounding-sensitive geometry.
set.seed(53+3+31);x<-matrix(rnorm(3*31),3)
fd<-get(".gseg_fast_dist",ns)
stopifnot(identical(fd(x),stats::dist(x)))
cat("PASS strict reliable-reference arithmetic, singular/cancellation repairs, original location type and public distances\n")
