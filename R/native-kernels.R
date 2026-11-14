.gseg_native_mst <- function(d,k) {
  if(any(!is.finite(d))) stop("NA/NaN/Inf in foreign function call (arg 1)")
  .Call(C_gseg_mst_exact,as.double(d),as.integer(attr(d,"Size")),as.integer(k))
}
.gseg_native_x5 <- function(n,E) .Call(C_gseg_graph_x5,as.integer(n),as.integer(E))
.gseg_graph_x3 <- function(deg,E) {
  v <- (deg[E[,1]]-1)*(deg[E[,2]]-1)
  if(!length(v) || max(v)<=2^53/length(v)) return(sum(v))
  value <- 0
  for(i in seq_len(nrow(E))) value <- value+(deg[E[i,1]]-1)*(deg[E[i,2]]-1)
  value
}
.gseg_fast_dist <- function(x) {
  n <- nrow(x)
  z <- .Call(C_gseg_dist_exact,as.double(x),as.integer(n),as.integer(ncol(x)))
  attributes(z) <- list(Size=n,Labels=dimnames(x)[[1L]],Diag=FALSE,Upper=FALSE,
                       method="euclidean",call=quote(stats::dist(x=x)),class="dist")
  z
}
