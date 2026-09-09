#' Construct a k-MST similarity graph
#'
#' @param x Optional data matrix with observations in rows.
#' @param dissimilarity Optional square dissimilarity matrix or a \code{dist}
#'   object. Supply exactly one of \code{x} and \code{dissimilarity}.
#' @param k Number of complete MST layers. The default is
#'   \code{floor(sqrt(N))}, where \code{N} is the number of observations;
#'   the largest allowed value is \code{floor(N/2)}.
#' @return A two-column edge matrix containing the union of the first
#'   \code{k} MST layers.
#' @examples
#' set.seed(1)
#' x <- matrix(rnorm(60), nrow = 20)
#' E <- gseg_kmst(x = x)
#' attr(E, "k")
#' @export
gseg_kmst <- function(x = NULL, dissimilarity = NULL, k = NULL) {
  if (is.null(x) == is.null(dissimilarity)) {
    stop("Supply exactly one of `x` and `dissimilarity`.")
  }

  if (!is.null(x)) {
    x <- as.matrix(x)
    if (nrow(x) < 3L) stop("At least three observations are required.")
    if (!is.numeric(x) || any(!is.finite(x))) {
      stop("`x` must be a finite numeric matrix.")
    }
    d <- stats::dist(x)
    N <- nrow(x)
  } else {
    if (inherits(dissimilarity, "dist")) {
      d <- dissimilarity
      N <- attr(d, "Size")
      if (is.null(N) || length(N) != 1L || !is.numeric(N) || !is.finite(N) ||
          N != as.integer(N) || N < 3L || !is.numeric(d) ||
          length(d) != N * (N - 1) / 2 || any(!is.finite(d)) || any(d < 0)) {
        stop("`dissimilarity` must be a finite nonnegative `dist` object with at least three observations.")
      }
    } else {
      D <- as.matrix(dissimilarity)
      if (nrow(D) != ncol(D)) stop("`dissimilarity` must be square.")
      if (!is.numeric(D)) stop("`dissimilarity` must be numeric.")
      if (any(!is.finite(D))) stop("`dissimilarity` must contain finite values.")
      if (any(D < 0)) stop("`dissimilarity` must be nonnegative.")
      if (any(abs(diag(D)) > sqrt(.Machine$double.eps))) {
        stop("`dissimilarity` must have a zero diagonal.")
      }
      if (!isTRUE(all.equal(D, t(D), tolerance = sqrt(.Machine$double.eps)))) {
        stop("`dissimilarity` must be symmetric.")
      }
      d <- stats::as.dist(D)
      N <- nrow(D)
    }
  }

  if (is.null(k)) k <- floor(sqrt(N))
  if (!is.numeric(k) || length(k) != 1L || !is.finite(k) || k < 1 ||
      k > floor(N / 2) || k != as.integer(k)) {
    stop("`k` must be an integer between 1 and floor(N/2).")
  }
  k <- as.integer(k)
  E <- ade4::mstree(d, k)
  E <- as.matrix(E[, 1:2, drop = FALSE])
  storage.mode(E) <- "integer"
  attr(E, "k") <- k
  E
}

#' Single-change graph scan from data or dissimilarities
#'
#' Constructs a k-MST with \code{k=floor(sqrt(N))} by default and reports
#' the max-type edge-count test (MET).
#' The original and weighted scans remain available through \code{statistics}.
#'
#' @inheritParams gseg_kmst
#' @param statistics Statistics passed to \code{gseg1}. The default
#'   \code{"m"} reports MET only; use \code{"g"} to request GET or
#'   \code{"all"} for all historical statistics.
#' @param ... Further arguments passed to \code{gseg1}.
#' @return The result from \code{gseg1}, with the constructed graph and k added.
#' @examples
#' set.seed(1)
#' x <- matrix(rnorm(60), nrow = 20)
#' ans <- gseg1_data(x = x, pval.appr = FALSE)
#' names(ans$scanZ)
#' @export
gseg1_data <- function(x = NULL, dissimilarity = NULL, k = NULL,
                       statistics = "m", ...) {
  if (is.null(x) == is.null(dissimilarity)) {
    stop("Supply exactly one of `x` and `dissimilarity`.")
  }
  N <- if (!is.null(x)) nrow(as.matrix(x)) else if (inherits(dissimilarity, "dist")) {
    attr(dissimilarity, "Size")
  } else {
    nrow(as.matrix(dissimilarity))
  }
  if (is.null(N) || N < 6L) stop("At least six observations are required for the scan.")
  E <- gseg_kmst(x = x, dissimilarity = dissimilarity, k = k)
  out <- gseg1(N, E, statistics = statistics, ...)
  out$graph <- E
  out$k <- attr(E, "k")
  out
}

#' Changed-interval graph scan from data or dissimilarities
#'
#' Constructs a k-MST with \code{k=floor(sqrt(N))} by default and reports
#' the max-type edge-count test (MET).
#'
#' @inheritParams gseg1_data
#' @param ... Further arguments passed to \code{gseg2}.
#' @return The result from \code{gseg2}, with the constructed graph and k added.
#' @examples
#' set.seed(1)
#' x <- matrix(rnorm(60), nrow = 20)
#' ans <- gseg2_data(x = x, pval.appr = FALSE)
#' names(ans$scanZ)
#' @export
gseg2_data <- function(x = NULL, dissimilarity = NULL, k = NULL,
                       statistics = "m", ...) {
  if (is.null(x) == is.null(dissimilarity)) {
    stop("Supply exactly one of `x` and `dissimilarity`.")
  }
  N <- if (!is.null(x)) nrow(as.matrix(x)) else if (inherits(dissimilarity, "dist")) {
    attr(dissimilarity, "Size")
  } else {
    nrow(as.matrix(dissimilarity))
  }
  if (is.null(N) || N < 6L) stop("At least six observations are required for the scan.")
  E <- gseg_kmst(x = x, dissimilarity = dissimilarity, k = k)
  out <- gseg2(N, E, statistics = statistics, ...)
  out$graph <- E
  out$k <- attr(E, "k")
  out
}
