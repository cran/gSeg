# Stable algebra for analytic p-values. Callers retain the existing
# discriminant domain and extrapolation/fallback decisions.
.gseg_saddlepoint_root <- function(b, discriminant) {
  # Equal to (sqrt(1+2*r*b)-1)/r, with continuous limit b at r=0.
  # Avoids both cancellation and forming 2*b.
  b / (0.5 + 0.5 * sqrt(discriminant))
}

.gseg_combine_tails <- function(p1, p2) {
  # Preserve the existing upper clipping of each approximate tail.
  p1 <- min(p1, 1)
  p2 <- min(p2, 1)
  p1 + (1 - p1) * p2
}
