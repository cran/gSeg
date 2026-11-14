library(gSeg)

must_error <- function(expr) inherits(try(force(expr), silent = TRUE), "try-error")

set.seed(1)
x <- matrix(rnorm(80), nrow = 20)
E <- gseg_kmst(x = x)
stopifnot(attr(E, "k") == floor(sqrt(nrow(x))))
stopifnot(must_error(gseg_kmst(x = x, k = 11L)))
stopifnot(must_error(gseg_kmst(x = x, k = "2")))
stopifnot(must_error(gseg_kmst(x = replace(x, 1, NA_real_))))
stopifnot(must_error(gseg1_data(x = matrix(rnorm(10), nrow = 5))))
small_x <- matrix(rnorm(24), nrow = 8)
stopifnot(must_error(gseg1_data(x = small_x)))
stopifnot(!must_error(gseg1_data(x = small_x, pval.appr = FALSE)))
stopifnot(must_error(gseg1_data(x = x, statistics = "typo")))
stopifnot(must_error(gseg1_data(dissimilarity = matrix("x", 20, 20))))

out <- gseg1_data(x = x, pval.appr = FALSE)
stopifnot(identical(names(out$scanZ), "max.type"))
invisible(capture.output(out_main <- gseg1(x, pval.appr = FALSE)))
stopifnot(identical(names(out_main$scanZ), "max.type"),
          out_main$k == floor(sqrt(nrow(x))),
          isTRUE(all.equal(out_main$scanZ, out$scanZ, tolerance = 0)))
invisible(capture.output(out_graph <-
  gseg1(nrow(x), out$graph, statistics = "m", pval.appr = FALSE)))
stopifnot(isTRUE(all.equal(out$scanZ, out_graph$scanZ, tolerance = 0)))
invisible(capture.output(out_graph_named <-
  gseg1(n = nrow(x), E = out$graph, statistics = "m", pval.appr = FALSE)))
invisible(capture.output(out_graph_explicit <-
  gseg1_graph(nrow(x), out$graph, statistics = "m", pval.appr = FALSE)))
stopifnot(isTRUE(all.equal(out_graph$scanZ, out_graph_named$scanZ, tolerance = 0)),
          isTRUE(all.equal(out_graph$scanZ, out_graph_explicit$scanZ, tolerance = 0)))

# The complete earlier positional interface remains valid.
invisible(capture.output(out_graph_positional <-
  gseg1(nrow(x), out$graph, "m", 4, 16, FALSE, TRUE, FALSE, 100)))
invisible(capture.output(out_graph_positional_ref <-
  gseg1_graph(nrow(x), out$graph, "m", 4, 16, FALSE, TRUE, FALSE, 100)))
stopifnot(isTRUE(all.equal(out_graph_positional$scanZ,
                          out_graph_positional_ref$scanZ, tolerance = 0)))

out_dist <- gseg1_data(dissimilarity = as.matrix(dist(x)), pval.appr = FALSE)
stopifnot(isTRUE(all.equal(out$scanZ, out_dist$scanZ, tolerance = 0)))
invisible(capture.output(out_dist_main <-
  gseg1(dissimilarity = as.matrix(dist(x)), pval.appr = FALSE)))
stopifnot(isTRUE(all.equal(out_dist$scanZ, out_dist_main$scanZ, tolerance = 0)))

out_get <- gseg1_data(x = x, statistics = "g", pval.appr = FALSE)
stopifnot(identical(names(out_get$scanZ), "generalized"))

out_all <- gseg1_data(x = x, statistics = "all", pval.appr = FALSE)
stopifnot(all(c("ori", "weighted", "generalized", "max.type") %in% names(out_all$scanZ)))

set.seed(2)
out_perm <- gseg1_data(x = x, pval.appr = FALSE, pval.perm = TRUE, B = 9)
perm_p <- vapply(out_perm$pval.perm, function(z) z$pval, numeric(1))
stopifnot(all(perm_p >= 0.1), all(perm_p <= 1))

set.seed(21)
invisible(capture.output(out_perm_one <-
  gseg1_data(x = x, pval.appr = FALSE, pval.perm = TRUE, B = 1)))
stopifnot(length(out_perm_one$pval.perm$max.type$maxZs) == 1L,
          out_perm_one$pval.perm$max.type$pval %in% c(0.5, 1))
set.seed(211)
invisible(capture.output(out_single_location <-
  gseg1_data(x = x, n0 = 10, n1 = 10, pval.appr = FALSE,
             pval.perm = TRUE, B = 3)))
stopifnot(length(out_single_location$pval.perm$max.type$maxZs) == 3L)

set.seed(22)
invisible(capture.output(out_interval <-
  gseg2_data(x = x, pval.appr = FALSE, pval.perm = TRUE, B = 7)))
stopifnot(length(out_interval$pval.perm$max.type$Zmax) == 7L)
stopifnot(nrow(out_interval$pval.perm$max.type$curve) == 7L)
stopifnot(out_interval$pval.perm$max.type$pval ==
            (1 + sum(out_interval$pval.perm$max.type$Zmax >=
                       out_interval$scanZ$max.type$Zmax)) / 8)
invisible(capture.output(out_interval_main <-
  gseg2(x, pval.appr = FALSE)))
invisible(capture.output(out_interval_graph <-
  gseg2_graph(nrow(x), out_interval_main$graph, pval.appr = FALSE)))
invisible(capture.output(out_interval_old <-
  gseg2(nrow(x), out_interval_main$graph, pval.appr = FALSE)))
invisible(capture.output(out_interval_named <-
  gseg2(n = nrow(x), E = out_interval_main$graph, pval.appr = FALSE)))
stopifnot(identical(names(out_interval_main$scanZ), "max.type"),
          isTRUE(all.equal(out_interval_graph$scanZ,
                           out_interval_old$scanZ, tolerance = 0)),
          isTRUE(all.equal(out_interval_graph$scanZ,
                           out_interval_named$scanZ, tolerance = 0)))
stopifnot(must_error(gseg1_data(x = x, pval.perm = TRUE, B = 0)))

set.seed(3)
n_discrete <- 30L
d_discrete <- 5L
tau <- n_discrete / 2L
y1_pool <- matrix(rnorm(d_discrete * tau), tau)
y2_pool <- matrix(rnorm(d_discrete * tau, 1 / sqrt(d_discrete)), tau)
y <- rbind(y1_pool[sample.int(tau, tau, replace = TRUE), , drop = FALSE],
           y2_pool[sample.int(tau, tau, replace = TRUE), , drop = FALSE])
y_unique <- unique(y)
E_discrete <- nnl(dist(y_unique), 1)
keys <- do.call(paste, as.data.frame(y))
id <- match(keys, unique(keys))
invisible(capture.output(out_discrete <-
  gseg1_discrete(n_discrete, E_discrete, id, pval.appr = FALSE)))
stopifnot(identical(names(out_discrete$scanZ), "max.type"))
stopifnot(must_error(gseg1_discrete(n_discrete, E_discrete, id,
                                    statistics = "typo", pval.appr = FALSE)))
stopifnot(must_error(gseg1(nrow(x), matrix(c(1, 21), nrow = 1),
                           pval.appr = FALSE)))
stopifnot(must_error(gseg1(10, t(combn(10, 2)), pval.appr = FALSE)))
stopifnot(must_error(gseg1(10, rbind(c(1, 2), c(2, 1)), pval.appr = FALSE)))
stopifnot(must_error(gseg2(10, rbind(c(1, 2), c(2, 1)), pval.appr = FALSE)))
stopifnot(must_error(gseg1_discrete(n_discrete, rbind(c(1, 2), c(2, 1)), id,
                                    pval.appr = FALSE)))
stopifnot(must_error(gseg2_discrete(n_discrete, rbind(c(1, 2), c(2, 1)), id,
                                    pval.appr = FALSE)))

set.seed(23)
invisible(capture.output(out_interval_discrete <-
  gseg2_discrete(n_discrete, E_discrete, id, pval.appr = FALSE,
                 pval.perm = TRUE, B = 7)))
stopifnot(all(vapply(out_interval_discrete$pval.perm,
                     function(z) length(z$Z) == 7L, logical(1))))

set.seed(24)
invisible(capture.output(out_single_discrete <-
  gseg1_discrete(n_discrete, E_discrete, id, pval.appr = FALSE,
                 pval.perm = TRUE, B = 1)))
stopifnot(all(vapply(out_single_discrete$pval.perm,
                     function(z) length(z$maxZs_a) == 1L ||
                       length(z$maxZs_u) == 1L, logical(1))))

# R's named-argument matching must not break formerly valid graph calls.
for (fun in c("gseg1", "gseg2")) {
  f <- get(fun)
  invisible(capture.output(reference <- f(n = nrow(x), E = E,
                                          statistics = "all", pval.appr = FALSE)))
  calls <- list(list(nrow(x), E = E), list(n = nrow(x), E),
                list(nrow(x), E), list(n = nrow(x), E = E))
  for (args in calls) {
    invisible(capture.output(actual <- do.call(f, c(args,
      list(statistics = "all", pval.appr = FALSE)))))
    stopifnot(isTRUE(all.equal(actual, reference, tolerance = 0)))
  }
  stopifnot(must_error(f(n = nrow(x), E = E, k = 1, pval.appr = FALSE)),
            must_error(f(nrow(x), E, k = 1, pval.appr = FALSE)),
            must_error(f(x = x, n = nrow(x), E = E, pval.appr = FALSE)))
}
# A vector is a one-coordinate sequence; dist input retains the same geometry.
invisible(capture.output(vector_result <- gseg1(x[, 1], pval.appr = FALSE)))
invisible(capture.output(matrix_result <- gseg1(x[, 1, drop = FALSE], pval.appr = FALSE)))
stopifnot(isTRUE(all.equal(vector_result, matrix_result, tolerance = 0)))
invisible(capture.output(dist_result <- gseg1(dissimilarity = dist(x), pval.appr = FALSE)))
invisible(capture.output(data_result <- gseg1(x, pval.appr = FALSE)))
stopifnot(isTRUE(all.equal(dist_result, data_result, tolerance = 0)))
