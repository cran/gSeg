# Preserve reference arithmetic when it is numerically reliable. Repair
# cancellation/zero-skew branches without changing the statistical formula.
.gseg_saddlepoint_root <- function(b, discriminant, r = NULL) {
  stable <- b / (0.5 + 0.5 * sqrt(discriminant))
  if (is.null(r)) return(stable)
  legacy <- (sqrt(discriminant) - 1) / r
  term <- (r * legacy) * (legacy / 2)
  residual <- legacy + term - b
  # Eight rounding units cover the short root/residual arithmetic chains.
  # Agreement with the stable root also guards near-multiple-root conditioning.
  close <- abs(legacy - stable) <= 8 * .Machine$double.eps *
    pmax(abs(legacy), abs(stable), .Machine$double.xmin)
  consistent <- abs(residual) <= 8 * .Machine$double.eps *
    pmax(abs(legacy) + abs(term) + abs(b), .Machine$double.xmin)
  keep <- r != 0 & is.finite(legacy) & close & consistent
  keep[is.na(keep)] <- FALSE
  stable[keep] <- legacy[keep]
  stable
}

.gseg_combine_tails <- function(p1, p2) {
  p1 <- min(p1, 1)
  p2 <- min(p2, 1)
  legacy <- 1 - (1 - p1) * (1 - p2)
  stable <- p1 + (1 - p1) * p2
  # Retain the original value unless cancellation exceeds rounding-level
  # agreement; a positive tail rounded to zero is always repaired.
  keep <- is.finite(legacy) && is.finite(stable) &&
    abs(legacy - stable) <= 4 * .Machine$double.eps *
      max(abs(stable), .Machine$double.xmin)
  if (keep) legacy else stable
}
