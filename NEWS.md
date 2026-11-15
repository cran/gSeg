# gSeg 1.2-1

- Preserve platform-reference stats::dist arithmetic, including ARM64 rounding,
  conversion semantics and original graph tie ordering. Use optimized finite
  distance kernels only for verified R builds; otherwise use stats::dist.
- Preserve original adjacency storage, neighbor ordering and graph-constructor
  dispatch while accelerating distance, graph, scan and moment calculations.
- Correct numerical boundary cases and reject exactly singular scan channels
  using structural zero-variance checks, including repeated observations.
- Avoid unused dense distance conversions during dispatch; normalize graph
  inputs once and use the observation count of the actual distance geometry.
- Add reference, near-tie, permutation, conversion and structural-variance
  regression tests. Internal candidate revisions are retained in the release
  evidence and were not separate published CRAN versions.

# gSeg 1.2

- Fixed the documented `statistics = "max"` alias in analytic scan calls.
- Reject `dist` objects in `x` instead of silently recomputing distances
  between distance-matrix rows; use `dissimilarity = ...`.
- Reject non-converged analytic quadrature and invalid analytic probability
  outputs explicitly. Existing skew-correction fallbacks remain available.
- Repeated-observation interval MET no longer silently substitutes zero for
  an undefined weighted component.

- Validate that requested k-MST layers produce exactly k*(N-1) distinct
  edges; fail explicitly if the underlying constructor reuses exhausted edges.

- Fixed changed-interval maximizer indexing and restricted all tied maximizers
  to the requested interval-length range, including repeated-observation scans.
- Corrected returned GET permutation paths and tied permutation tail curves.
- Reject undefined candidate scan locations instead of silently dropping
  non-finite permutation maxima, which could yield spuriously small p-values.
- Reject zero-width ranges for analytic scan approximations; permutation
  calibration remains available for those restricted ranges.
- NNL construction excludes self-links under zero-distance ties, validates
  inputs/layer counts, and detects exhausted layers.

- Fixed mixed named/positional legacy graph calls, including `gseg1(n, E=E)`
  and `gseg1(n=n, E)` (and the corresponding `gseg2` calls).
- Supplying `k` with an explicit graph now raises an error rather than silently
  ignoring it. Set `k` only when constructing a graph from data/dissimilarities.
- Reworked package documentation and examples around the main data-level calls,
  with explicit graph entry points and a migration guide.

- `gseg1()` and `gseg2()` now accept observations or dissimilarities directly.
- Added explicit `gseg1_graph()` and `gseg2_graph()` functions for supplied
  similarity graphs.
- Calls using the earlier `gseg1(n, E, ...)` and `gseg2(n, E, ...)` forms
  remain supported.
- `gseg1_data()` and `gseg2_data()` remain available as aliases.
- The main data-level calls construct a k-MST with `k = floor(sqrt(N))` by
  default and report MET.

# gSeg 1.1

- Added data-level `gseg1_data()` and `gseg2_data()` interfaces.
- The data-level interfaces use a k-MST with `k = floor(sqrt(N))` by default.
- Default calls report MET. GET, OET, and WET remain available through
  the explicit `statistics = "g"`, `"o"`, or `"w"` options, or together
  through `statistics = "all"`.
- Existing graph-level interfaces remain available.
- Monte Carlo permutation p-values use the standard plus-one correction.
- Corrected changed-interval permutation output so that it retains one scan
  maximum per permutation and uses all `B` maxima to compute the p-value.
- Added explicit checks for analytic small-sample limits, degenerate graphs,
  duplicate edges, and malformed inputs.
