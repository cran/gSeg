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
