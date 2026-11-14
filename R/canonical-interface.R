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
    if (inherits(x, "dist")) {
      stop("A `dist` object contains dissimilarities, not observations; supply it as `dissimilarity = ...`.")
    }
    x <- as.matrix(x)
    if (ncol(x) < 1L) stop("`x` must have at least one observation coordinate.")
    if (nrow(x) < 3L) stop("At least three observations are required.")
    if (!is.numeric(x) || any(!is.finite(x))) {
      stop("`x` must be a finite numeric matrix.")
    }
    d <- .gseg_fast_dist(x)
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
      if (N < 3L) stop("At least three observations are required.")
    }
  }

  if (is.null(k)) k <- floor(sqrt(N))
  if (!is.numeric(k) || length(k) != 1L || !is.finite(k) || k < 1 ||
      k > floor(N / 2) || k != as.integer(k)) {
    stop("`k` must be an integer between 1 and floor(N/2).")
  }
  k <- as.integer(k)
  E <- .gseg_native_mst(d, k)
  E <- as.matrix(E[, 1:2, drop = FALSE])
  storage.mode(E) <- "integer"
  if (nrow(E) != k * (N - 1L) ||
      anyDuplicated(paste(pmin(E[, 1L], E[, 2L]),
                          pmax(E[, 1L], E[, 2L]), sep = ":"))) {
    stop("The graph constructor could not produce `k` complete edge-disjoint MST layers; set a smaller `k` or supply your own graph explicitly.")
  }
  attr(E, "k") <- k
  E
}

.gseg_data_input <- function(x, dissimilarity, k) {
  if (inherits(x, "dist")) {
    stop("A `dist` object contains dissimilarities, not observations; supply it as `dissimilarity = ...`.")
  }
  if (is.null(x) == is.null(dissimilarity)) {
    stop("Supply exactly one of `x` and `dissimilarity`.")
  }
  N <- if (!is.null(x)) nrow(as.matrix(x)) else if (inherits(dissimilarity, "dist")) {
    attr(dissimilarity, "Size")
  } else {
    nrow(as.matrix(dissimilarity))
  }
  if (is.null(N) || length(N) != 1L || N < 6L) {
    stop("At least six observations are required for the scan.")
  }
  graph <- gseg_kmst(x = x, dissimilarity = dissimilarity, k = k)
  list(N = N, graph = graph, k = attr(graph, "k"))
}

.gseg_graph_call <- function(x, dissimilarity, n, E, k, argument.names = character()) {
  # R matches named arguments before positional ones. Retain the earlier
  # gseg1(n, E) / gseg2(n, E) API even when just one argument is named.
  named.graph <- !is.null(n) || !is.null(E)
  explicit.data <- any(c("x", "dissimilarity") %in% argument.names)
  if (named.graph) {
    if (!explicit.data && is.null(n) && !is.null(E) && is.numeric(x) &&
        length(x) == 1L && is.null(dissimilarity)) {
      n <- x
      x <- NULL
    } else if (!explicit.data && !is.null(n) && is.null(E) && !is.null(x) &&
               is.null(dissimilarity)) {
      E <- x
      x <- NULL
    }
    if (is.null(n) || is.null(E)) {
      stop("Supply both `n` and `E` for a graph-level call.")
    }
    if (!is.null(x) || !is.null(dissimilarity)) {
      stop("Do not combine `x` or `dissimilarity` with graph arguments `n` and `E`.")
    }
    if (!is.null(k)) stop("`k` is only used when constructing a graph from data or dissimilarities.")
    return(list(graph = TRUE, N = n, E = E))
  }

  positional.edges <- if (is.null(dissimilarity)) NULL else
    try(as.matrix(dissimilarity), silent = TRUE)
  positional.graph <- !explicit.data && is.numeric(x) && length(x) == 1L &&
    !is.null(dissimilarity) && !inherits(positional.edges, "try-error") &&
    ncol(positional.edges) == 2L
  if (positional.graph) {
    if (!is.null(k)) stop("`k` is only used when constructing a graph from data or dissimilarities.")
    return(list(graph = TRUE, N = x, E = dissimilarity))
  }
  list(graph = FALSE)
}

#' Single-change graph scan
#'
#' With data or a dissimilarity matrix, constructs a k-MST with
#' \code{k=floor(sqrt(N))} by default and reports the max-type edge-count test.
#' Calls of the earlier form \code{gseg1(n, E, ...)} remain supported.
#'
#' @param x Optional data matrix with observations in rows. In the earlier
#'   graph-level interface, the first positional argument may instead be the
#'   number of observations.
#' @param dissimilarity Optional square dissimilarity matrix or a \code{dist}
#'   object. In the earlier graph-level interface, the second positional
#'   argument may instead be a two-column edge matrix.
#' @param statistics Scan statistic. The default \code{"m"} reports the
#'   max-type edge-count test; use \code{"g"} for the generalized edge-count
#'   test or \code{"all"} to request all available statistics.
#' @param n0,n1 Beginning and end of the candidate change-point range.
#' @param pval.appr Whether to calculate the analytic p-value approximation.
#' @param skew.corr Whether to use the skewness correction.
#' @param pval.perm Whether to calculate a permutation p-value.
#' @param B Number of permutations.
#' @param k Number of MST layers. The data-level default is
#'   \code{floor(sqrt(N))}.
#' @param n,E Optional named graph-level arguments retained for compatibility.
#' @return A list containing scan statistics and requested p-values. Data-level
#'   calls also return the constructed graph and the selected value of \code{k}.
#' @examples
#' set.seed(1)
#' x <- matrix(rnorm(60), nrow = 20)
#' ans <- gseg1(x, pval.appr = FALSE)
#' names(ans$scanZ)
#' @export
gseg1 <- function(x = NULL, dissimilarity = NULL, statistics = "m",
                  n0 = NULL, n1 = NULL, pval.appr = TRUE,
                  skew.corr = TRUE, pval.perm = FALSE, B = 100,
                  k = NULL, n = NULL, E = NULL) {
  graph.call <- .gseg_graph_call(x, dissimilarity, n, E, k,
                                  names(sys.call())[-1L])
  if (graph.call$graph) {
    N <- graph.call$N
    if (is.null(n0)) n0 <- 0.05 * N
    if (is.null(n1)) n1 <- 0.95 * N
    return(.gseg1_graph_impl(N, graph.call$E, statistics, n0, n1,
                             pval.appr, skew.corr, pval.perm, B))
  }

  input <- .gseg_data_input(x, dissimilarity, k)
  if (is.null(n0)) n0 <- 0.05 * input$N
  if (is.null(n1)) n1 <- 0.95 * input$N
  out <- .gseg1_graph_impl(input$N, input$graph, statistics, n0, n1,
                           pval.appr, skew.corr, pval.perm, B)
  out$graph <- input$graph
  out$k <- input$k
  out
}

#' Changed-interval graph scan
#'
#' With data or a dissimilarity matrix, constructs a k-MST with
#' \code{k=floor(sqrt(N))} by default and reports the max-type edge-count test.
#' Calls of the earlier form \code{gseg2(n, E, ...)} remain supported.
#'
#' @inheritParams gseg1
#' @param l0,l1 Minimum and maximum candidate interval lengths.
#' @return A list containing scan statistics and requested p-values. Data-level
#'   calls also return the constructed graph and the selected value of \code{k}.
#' @examples
#' set.seed(1)
#' x <- matrix(rnorm(60), nrow = 20)
#' ans <- gseg2(x, pval.appr = FALSE)
#' names(ans$scanZ)
#' @export
gseg2 <- function(x = NULL, dissimilarity = NULL, statistics = "m",
                  l0 = NULL, l1 = NULL, pval.appr = TRUE,
                  skew.corr = TRUE, pval.perm = FALSE, B = 100,
                  k = NULL, n = NULL, E = NULL) {
  graph.call <- .gseg_graph_call(x, dissimilarity, n, E, k,
                                  names(sys.call())[-1L])
  if (graph.call$graph) {
    N <- graph.call$N
    if (is.null(l0)) l0 <- 0.05 * N
    if (is.null(l1)) l1 <- 0.95 * N
    return(.gseg2_graph_impl(N, graph.call$E, statistics, l0, l1,
                             pval.appr, skew.corr, pval.perm, B))
  }

  input <- .gseg_data_input(x, dissimilarity, k)
  if (is.null(l0)) l0 <- 0.05 * input$N
  if (is.null(l1)) l1 <- 0.95 * input$N
  out <- .gseg2_graph_impl(input$N, input$graph, statistics, l0, l1,
                           pval.appr, skew.corr, pval.perm, B)
  out$graph <- input$graph
  out$k <- input$k
  out
}

#' Single-change graph scan with a supplied graph
#'
#' @param n Number of observations.
#' @param E Two-column edge matrix.
#' @inheritParams gseg1
#' @return A list containing scan statistics and requested p-values.
#' @export
gseg1_graph <- function(n, E, statistics = "m", n0 = 0.05 * n,
                        n1 = 0.95 * n, pval.appr = TRUE,
                        skew.corr = TRUE, pval.perm = FALSE, B = 100) {
  .gseg1_graph_impl(n, E, statistics, n0, n1, pval.appr, skew.corr,
                     pval.perm, B)
}

#' Changed-interval graph scan with a supplied graph
#'
#' @param n Number of observations.
#' @param E Two-column edge matrix.
#' @inheritParams gseg2
#' @return A list containing scan statistics and requested p-values.
#' @export
gseg2_graph <- function(n, E, statistics = "m", l0 = 0.05 * n,
                        l1 = 0.95 * n, pval.appr = TRUE,
                        skew.corr = TRUE, pval.perm = FALSE, B = 100) {
  .gseg2_graph_impl(n, E, statistics, l0, l1, pval.appr, skew.corr,
                     pval.perm, B)
}

#' Data-level alias for \code{gseg1}
#'
#' @inheritParams gseg1
#' @return The result from \code{gseg1}.
#' @export
gseg1_data <- function(x = NULL, dissimilarity = NULL, k = NULL,
                       statistics = "m", ...) {
  gseg1(x = x, dissimilarity = dissimilarity, statistics = statistics,
        k = k, ...)
}

#' Data-level alias for \code{gseg2}
#'
#' @inheritParams gseg2
#' @return The result from \code{gseg2}.
#' @export
gseg2_data <- function(x = NULL, dissimilarity = NULL, k = NULL,
                       statistics = "m", ...) {
  gseg2(x = x, dissimilarity = dissimilarity, statistics = statistics,
        k = k, ...)
}
