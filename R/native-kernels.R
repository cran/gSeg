.gseg_native_mst <- function(d,k) {
  if (!.gseg_mst_conversion_ok(d)) return(.gseg_mst_legacy_conversion(d,k))
  .Call(C_gseg_mst_exact,d,attr(d,"Size"),k)
}
.gseg_native_x5 <- function(n,E,Ebynode=NULL) {
  if(is.object(E)||is.object(n)){
    if(is.null(Ebynode))Ebynode<-.gseg_reference_adjacency(E,n)
    value<-0
    for(i in 1:nrow(E)){
      j<-E[i,1];k<-E[i,2]
      value<-value+length(which(!is.na(match(Ebynode[[j]],Ebynode[[k]]))))
    }
    return(value)
  }
  .Call(C_gseg_graph_x5,n,as.integer(E))
}
.gseg_graph_x3 <- function(deg,E) {
  if(is.object(E)||is.object(deg)){
    value<-0
    for(i in 1:nrow(E))value<-value+(deg[E[i,1]]-1)*(deg[E[i,2]]-1)
    return(value)
  }
  v <- (deg[E[,1]]-1)*(deg[E[,2]]-1)
  if(!length(v) || max(v)<=2^53/length(v)) return(sum(v))
  value <- 0
  for(i in seq_len(nrow(E))) value <- value+(deg[E[i,1]]-1)*(deg[E[i,2]]-1)
  value
}
# Follow the original package's public distance implementation on every
# platform; compiler-specific contraction can otherwise change MST ties.
.gseg_fast_dist <- function(x) .gseg_reference_dist(x)

.gseg_edge_keys <- function(E, n) {
  if(is.object(E)||is.object(n))
    return(paste(pmin(E[,1],E[,2]),pmax(E[,1],E[,2]),sep=":"))
  # n <= floor(sqrt(2^53)): all integer products/sums are exactly
  # representable. The mapping is injective for valid 1 <= endpoints <= n.
  # Keep the original representation beyond this exact-integer bound.
  a <- pmin(E[,1], E[,2])
  b <- pmax(E[,1], E[,2])
  if (n <= 94906265) a + as.double(n) * (b - 1) else
    paste(a, b, sep = ":")
}

.gseg_adjacency <- function(E,n) {
  if(is.object(E)||is.object(n))return(.gseg_reference_adjacency(E,n))
  .Call(C_gseg_adjacency,E,n)
}
