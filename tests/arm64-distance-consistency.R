library(gSeg)
# A Linux ARM64/GCC rounding witness. Separate and fused multiply-adds
# differ by one ulp on two nearly tied edges; the resulting graphs can
# give different scan statistics and analytic/permutation p-values.
# Require exact reference distances, not an enlarged numerical tolerance.
v <- c(0.270870490107157, 0.810618837050124, 1.91883954684686, 0.550899107103773, 
1.76951289776982, -1.22908355082537, -1.06844372800939, -0.0365802083408006, 
1.52883370717159, 1.56322493215698, -0.847546905575772, 0.749838372845859, 
-1.58182656527391, -1.84244822892253, 1.94683358604506, -2.59436374146425, 
0.692942363583685, -0.857298144829627, 0.585843868347192, -1.13508613594088, 
0.85867848450079, 0.340396444072215, 0.934919568437006, 0.157201590378431, 
0.757265020745389, -1.49818227953382, 0.393378006247957, -2.9397554017498, 
0.169921160765807, -0.354361436096665, -2.15163445793965)
w <- c(0.270870490107157, 0.810618837050124, 1.91883954684686, 0.169921160765807, 
1.76951289776982, -1.22908355082537, -1.06844372800939, -0.0365802083408006, 
1.52883370717159, 1.56322493215698, -0.847546905575772, 0.749838372845859, 
-1.58182656527391, -1.84244822892253, 1.94683358604506, -2.59436374146425, 
0.692942363583685, -0.857298144829627, 0.585843868347192, -1.13508613594088, 
0.85867848450079, 0.340396444072215, 0.934919568437006, 0.157201590378431, 
0.757265020745389, -1.49818227953382, 0.393378006247957, -2.9397554017498, 
0.550899107103773, -0.354361436096665, -2.15163445793965)
order <- c(27L, 14L, 30L, 10L, 22L, 12L, 7L, 29L, 9L, 26L, 2L, 4L, 1L, 
8L, 11L, 28L, 18L, 24L, 19L, 15L, 3L, 20L, 23L, 17L, 25L, 13L, 
16L, 6L, 5L, 21L)

x <- cbind(rbind(rep(0,length(v)),v,w),matrix(0,3,9))
r <- as.double(stats::dist(rbind(0,v)))[1]
for(i in 1:27){z<-c(w,rep(0,9));z[length(v)+1L]<-r*.03*i;x<-rbind(x,z)}
x <- x[order,,drop=FALSE]
rownames(x) <- paste0("witness_",seq_len(nrow(x)))
fd <- getFromNamespace(".gseg_fast_dist","gSeg")
bits <- function(z)writeBin(as.double(z),raw(),size=8)
stopifnot(identical(bits(fd(x)),bits(stats::dist(x))),
          identical(fd(x),stats::dist(x)))
quiet <- function(expr){invisible(capture.output(value<-force(expr)));value}
for(k in c(1L,2L)){
 E <- gseg_kmst(x=x,k=k)
 stopifnot(identical(E,gseg_kmst(dissimilarity=stats::dist(x),k=k)),
           identical(E,gseg_kmst(dissimilarity=as.matrix(stats::dist(x)),k=k)))
 for(mode in 1:2){
  fun <- if(mode==1)gseg1 else gseg2
  bounds <- if(mode==1)list(n0=3,n1=27)else list(l0=3,l1=27)
  args <- c(list(statistics="all",k=k,pval.appr=TRUE,pval.perm=TRUE,B=39),bounds)
  set.seed(901);a <- quiet(do.call(fun,c(list(x=x),args)))
  set.seed(901);b <- quiet(do.call(fun,c(list(dissimilarity=stats::dist(x)),args)))
  set.seed(901);c <- quiet(do.call(fun,c(list(dissimilarity=as.matrix(stats::dist(x))),args)))
  stopifnot(identical(a,b),identical(a,c))
 }
}
cat("PASS rounding-sensitive raw/dist/matrix graphs, complete scans, analytic p-values and seeded permutations\n")
