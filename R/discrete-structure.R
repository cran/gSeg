# Exact zero-variance predicates for the expanded repeated-observation graph.
# For 2 <= t <= n-2, Rd has zero variance iff weighted degrees are constant;
# Rw has zero variance iff off-diagonal weights are additive (a_i + a_j).
# Ro = const - 2 Rw - (2t-n)/(n-2) Rd; its two projections are orthogonal.
# Use category counts/degrees already computed by getMV_discrete1: no extra
# graph conversion/indexing and no expanded N-by-N matrix or variance tolerance.
.gseg_check_discrete_structure <- function(mk, degree, union.degree, statistics, sizes,
                                           channel.names = c("average", "union")) {
  n <- sum(mk)
  simple.additive <- function(d) {
    all(d == 0) || all(d == n-1) ||
      (sum(mk[d == n-1]) == 1 && all(d == n-1 | d == 1)) ||
      (sum(mk[d == 0]) == 1 && all(d == 0 | d == n-2))
  }
  # Reduce (degree - 2)/multiplicity without potentially inexact products.
  rational.constant <- function(a,b) {
    u <- abs(a); v <- b
    while(any(v != 0)) {
      ix <- which(v != 0); z <- u[ix] %% v[ix]
      u[ix] <- v[ix]; v[ix] <- z
    }
    numerator <- a/u; denominator <- b/u
    all(numerator == numerator[1L] & denominator == denominator[1L])
  }
  duplicated <- which(mk > 1)
  average.additive <- if(!length(duplicated)) simple.additive(degree) else
    if(length(duplicated) > 1L) FALSE else {
      center <- duplicated[1L]; others <- setdiff(seq_along(mk),center)
      length(others) <= 1L ||
        (degree[center] == length(others) && all(degree[others] == 1))
    }
  zero.w <- c(average=average.additive, union=simple.additive(union.degree))
  zero.d <- c(average=rational.constant(degree-2,mk),
              union=all(union.degree == union.degree[1L]))
  need.w <- any(statistics %in% c('w','weighted','m','max','g','generalized','all'))
  need.d <- any(statistics %in% c('m','max','g','generalized','all'))
  need.o <- any(statistics %in% c('o','ori','original','all'))
  invalid <- (need.w & zero.w) | (need.d & zero.d) |
    (need.o & zero.w & (zero.d | any(2*sizes == n)))
  if(any(invalid)) stop('The requested scan statistic has exactly zero variance in the ',
    paste(unique(channel.names[invalid]),collapse=' and '),
    ' channel; use a nondegenerate similarity graph or a statistic with positive variance.',
    call.=FALSE)
  invisible(NULL)
}
