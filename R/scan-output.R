# Decode all tied maximizers from the requested interval domain. Scan surfaces
# are flattened row-wise and (i, j) denotes observations (i + 1), ..., j.
.gseg_interval_argmax <- function(values, eligible, n) {
  maxima <- eligible[which(values[eligible] == max(values[eligible]))]
  c((maxima - 1L) %/% n + 1L, (maxima - 1L) %% n + 1L)
}

# The same inclusive-tail, plus-one convention is used for both the reported
# scalar p-value and the plotted permutation tail curve, including ties.
.gseg_perm_curve <- function(maxZs, value.name = "maxZs") {
  if (!length(maxZs) || any(!is.finite(maxZs))) {
    stop("Permutation scan maxima must be finite; the graph may be degenerate.")
  }
  p <- (length(maxZs) + 2 - match(maxZs, sort(maxZs))) / (length(maxZs) + 1)
  out <- cbind(maxZs, p)
  colnames(out) <- c(value.name, "p")
  out
}

.gseg_finite_maxima <- function(z) {
  if (any(!is.finite(z))) {
    stop("A permutation scan is undefined at a requested location; use a graph with positive scan variance throughout the requested range.")
  }
  z
}

.gseg_check_single_paths <- function(scanZ, n0, n1, repeated = FALSE) {
  fields <- if (repeated) {
    list(ori = c("Zo_a", "Zo_u"), weighted = c("Zw_a", "Zw_u"),
         generalized = c("S_a", "S_u"), max.type = c("M_a", "M_u"))
  } else {
    list(ori = "Z", weighted = "Zw", generalized = "S", max.type = "M")
  }
  for (statistic in names(scanZ)) {
    for (field in fields[[statistic]]) {
      if (any(!is.finite(scanZ[[statistic]][[field]][n0:n1]))) {
        stop("The requested scan statistic is undefined at one or more candidate locations; use a graph with positive scan variance throughout the requested range.")
      }
    }
  }
  invisible(NULL)
}

# A numerical value is not a usable approximation if quadrature failed.
.gseg_integrate <- function(...) {
  result <- stats::integrate(...)
  if (!identical(result$message, "OK") || !is.finite(result$value)) {
    stop("Analytic p-value integration did not converge: ", result$message)
  }
  result
}

.gseg_analytic <- function(compute) {
  result <- tryCatch(compute(), error = function(e) {
    stop("Analytic p-value approximation could not be computed: ",
         conditionMessage(e),
         ". Set `pval.appr = FALSE` and use justified permutation calibration.",
         call. = FALSE)
  })
  valid <- is.list(result) && length(result) > 0L &&
    all(vapply(result, function(v) is.numeric(v) && length(v) == 1L &&
                 is.finite(v) && v >= 0 && v <= 1, logical(1)))
  if (!valid) {
    stop("Analytic p-value approximation did not produce finite probabilities in [0, 1]; set `pval.appr = FALSE` and use justified permutation calibration.",
         call. = FALSE)
  }
  result
}
