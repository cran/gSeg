# gSeg

Graph-based scans for a single change-point or one changed interval in a
sequence of multivariate or object-valued observations.

## Default use

```r
library(gSeg)
set.seed(1)
x <- rbind(matrix(rnorm(200), 40), matrix(rnorm(200, 0.5), 40))
fit <- gseg1(x)
fit$k                       # floor(sqrt(nrow(x))) = 8
fit$scanZ$max.type$tauhat
fit$pval.appr$max.type
```

`gseg1()` accepts observations in rows; a numeric vector is a one-dimensional
sequence. `gseg2()` uses the same inputs to scan for one changed interval.
Both construct a k-MST with `k = floor(sqrt(N))` and report the max-type scan
(MET) by default. `k` may be set explicitly. The graph and realized `k` are
returned, so the representation used in an analysis can be recorded.

For non-Euclidean observations, supply pairwise dissimilarities:

```r
fit <- gseg1(dissimilarity = dist(x))
```

`statistics = "g"`, `"w"`, `"o"`, or `"all"` selects GET, WET, OET, or all
four scans. Analytic p-values are approximations, not exact finite-sample
calibration. For exchangeable observations, Monte Carlo permutation p-values
can be requested with `pval.appr = FALSE, pval.perm = TRUE, B = 999`.
The permutation test compares scan maxima and uses the plus-one correction
with ties included. Arbitrarily serially dependent observations are not made
exchangeable by this interface; ordinary permutations are not justified merely
by calling these functions. Input validation does not certify an analytic
approximation for every supplied graph.

## Existing graph-based code

```r
E <- gseg_kmst(x = x, k = 3)
fit <- gseg1_graph(nrow(x), E)
fit <- gseg1(nrow(x), E)       # earlier API still works
fit <- gseg1(nrow(x), E = E)   # mixed named/positional calls also work
```

An explicitly supplied graph is used as supplied. `k` cannot be specified in
that call because no graph is being constructed. `gseg1_data()` and
`gseg2_data()` remain supported aliases; new code can call `gseg1()` and
`gseg2()` directly. Repeated-observation scans remain available through
`gseg1_discrete()` and `gseg2_discrete()`.

The requested MST union is verified to have `k * (N - 1)` distinct edges.
The counting bound `k <= floor(N/2)` alone does not guarantee greedy tree
packing: a previous tree can exhaust a vertex. Such inputs now raise a clear
construction error; choose a smaller explicit `k` or supply a graph. No silent
fallback changes `k` or the similarity representation.

Pass a `dist` object through `dissimilarity = ...`; using it as `x` is rejected
so that its geometry cannot silently change. Analytic integrations must report
successful convergence and finite probabilities in [0, 1]. A failed analytic
approximation raises an informative error; set `pval.appr = FALSE` and use
justified permutation calibration rather than treating a failed approximation
as a probability.

The package uses registered native kernels for exact MST construction, finite-data Euclidean distances and existing skew-correction graph counts. These implementations retain the same graph selection, tie handling and statistical formulas; floating-point contraction is disabled during compilation.
