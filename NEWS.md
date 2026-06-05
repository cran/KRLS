# KRLS 1.7-0

Adds the two convention items deferred from Chad Hazlett's review of
the 1.6 line: a `cat_columns` argument for kbal/GPSS-style one-hot
rescaling of categorical inputs, and a `b_maxvarK` bandwidth selector
that becomes the new default for `sigma`.

## New arguments

* `krls(..., cat_columns = NULL)` — character vector of column names
  (or integer indices) flagging categorical inputs. When supplied,
  those columns are one-hot encoded with all levels (no reference
  cell) and multiplied by `sqrt(0.5)`, matching `kbal::one_hot` and
  `gpss:::one_hot`. The `sqrt(0.5)` factor compensates for the
  doubled squared-Euclidean distance from carrying both "in level k"
  and "not in level k" indicators, restoring a Hamming-like distance
  between two observations differing in one category. Continuous
  columns are still standardized to `sd = 1`.

  When `cat_columns` is left at `NULL`, `krls()` scans the input and
  emits a single warning if it finds columns that look categorical
  (factor / character / logical, or numeric with `<= 10` unique
  values) but were not declared. There is no autodetection; the
  warning is a friendly nudge that mirrors how `kbal` and `gpss`
  handle the same ambiguity. Pass `cat_columns = integer(0)` (or
  `character(0)`) to acknowledge the situation and silence the
  warning.

## Changed defaults

* `sigma = NULL` (the default) now triggers `b_maxvarK()` — bandwidth
  chosen by maximizing the variance of the off-diagonal kernel
  entries. The pre-1.7 default `sigma = ncol(X)` was a dimensional
  heuristic; the maxvarK choice picks the bandwidth that makes the
  columns of K most informative and matches the convention in `kbal`
  and `gpss`. A one-line `message()` notes the change at first call
  and tells users how to pin pre-1.7 behavior
  (`sigma = ncol(X_processed)`).

  Bit-identical reproduction of v1.5-2 and v1.6-x fits remains
  available by pinning `sigma` and `approx = "none"` explicitly.

## New exports

* `b_maxvarK(X_processed)` — exact-path bandwidth selector. Returns
  the sigma that maximizes off-diagonal `var(K)` over the search
  interval.
* `b_maxvarK_nystrom(X_processed, Z_processed)` — Nystrom-path
  bandwidth selector that operates on the n × m cross-kernel
  `C = K(X, Z)`.

## Fit object additions

* `$X_proc` — the kernel-space (preprocessed) matrix used internally.
  Identical column count to `$X` when `cat_columns` is not used.
* `$prep` — preprocessing metadata (centers, scales, factor levels)
  for `predict()` to re-apply the same transformation to newdata.

## Compatibility

* All 188 tests from 1.6-1 still pass.
* Fits produced by 1.6-x without `$prep` continue to work in
  `predict()` via a legacy code path.

# KRLS 1.6-1

Bug-fix patch on top of 1.6-0.

## Fixes

* `krls(derivative = TRUE, vcov = FALSE)` under the default
  `approx = "auto"` no longer crashes with the cryptic
  `object 'vcovmatc' not found`. The derivative+vcov guard was
  previously only enforced under `approx = "none"` (pre-dispatch);
  it now re-runs after the auto-dispatch resolves and produces a
  clear error when the exact path is chosen. The same combination
  remains supported under `approx = "nystrom"`, which returns
  point-estimate derivatives without SEs.
* Auto-dispatch now validates `nystrom_m` before using it as the
  dispatch threshold or printing the switchover message. Previously
  `nystrom_m = Inf` produced an opaque "missing value where
  TRUE/FALSE needed", and `nystrom_m = 1.5` announced
  `nystrom_m = 1` before the real validator caught it. Both now
  error cleanly with the same "must be a finite integer" message.

## Documentation

* `README.md` and `vignette("krls-nystrom-scaling")` now describe
  the default `nystrom_m` as `min(500, N)` (was still listed as the
  pre-1.6-0 `min(500, ceiling(sqrt(N)))` heuristic). The scaling
  vignette's benchmark code now uses the actual default rather than
  `ceiling(sqrt(n))`.

# KRLS 1.6-0

Refinements from Chad Hazlett's review of the Nystrom path. Three
user-visible additions; one default-behavior change at large samples.

## Auto-dispatch (`approx = "auto"` is now the default)

* `krls()` now defaults to `approx = "auto"`, which uses the exact
  path when `N <= nystrom_m` and switches to `approx = "nystrom"`
  with a one-line `message()` when `N > nystrom_m`. The new behavior
  benefits casual users who do not realize the Nystrom path exists,
  while preserving the call-site honesty of an explicit
  approximation flag.
* `approx = "none"` forces the exact path regardless of `N` and
  reproduces every KRLS <= 1.5-2 fit bit-for-bit. Users who need the
  historical default for reproducibility should set this explicitly.
* Conditions that aren't supported under Nystrom (non-Gaussian
  kernel, non-null `eigtrunc`, `derivative = TRUE` with
  `vcov = FALSE`) fall through to the exact path silently under
  `approx = "auto"`.

## Larger default `nystrom_m`

* Default `nystrom_m` is now `min(500, N)` (was `min(500, ceiling(sqrt(N)))`).
  The old `sqrt(N)` floor gave overly coarse approximations at moderate
  `N` (e.g. `m = 100` at `N = 10000`) given how cheap the cached-SVD
  Nystrom path is. The new default caps at 500 regardless and uses
  every row as a landmark for small samples.

## New `landmark_seed` argument with RNG-hygiene

* `krls(..., landmark_seed = NULL)` lets users pin Nystrom landmark
  selection without touching their global RNG state. When non-NULL,
  `.Random.seed` is saved before landmark selection and restored on
  exit, so the seed is consumed locally and the caller's downstream
  draws are unaffected. Two krls() calls with the same `landmark_seed`
  produce bit-identical landmarks regardless of the surrounding
  `set.seed()` context.

## Other

* `man/krls.Rd` now states the Nystrom approximation explicitly as
  `K ≈ C W_reg^{-1} C'` after the relative-ridge floor, matching the
  internal implementation.

# KRLS 1.5-2

## Fixes

* Nystrom fits with an underflowed cross-kernel (very small `sigma`
  for the data scale, or landmarks placed far from the observations)
  now signal a clear error pointing at sigma/landmark scale, instead
  of failing later inside the lambda-bound search with a cryptic
  `missing value where TRUE/FALSE needed`. The guard runs whether
  `lambda` is supplied directly or auto-selected.

# KRLS 1.5-1

Bug-fix patch for the Nystrom path that landed in 1.4-0 / 1.4-1.

## Fixes

* Nystrom lambda search at small `nystrom_m` (notably `m = 1`) no
  longer collapses. The previous bound heuristic could only satisfy
  `EDF >= 1` at `U = 0`, leaving the golden-section search running
  over a near-zero interval and selecting an essentially unregularized
  lambda. The bound logic now anchors at the dominant landmark
  eigenvalue and grows U as needed; the search runs over a usable
  window for every supported m.
* `landmark_method = "kmeans"` now handles `nystrom_m == nrow(X)`
  by short-circuiting to one cluster per row, instead of letting
  `stats::kmeans()` error on `centers >= n`.
* When `lambda_method = "gcv"` and `print.level > 1`, both backends
  now correctly label the printed scalar as the GCV minimum (was
  hard-coded as "Loo-Loss").
* `man/krls.Rd` no longer claims `landmark_method = "kmeans"` is
  reserved for a future release; it documents the actual behavior.
* `dev/03_nystrom_bootstrap_validation.R` assigns column names to its
  simulated X and prefers `pkgload::load_all()` when available so the
  script validates the source tree rather than a stale installed
  package.

# KRLS 1.5-0

## New: `lambda_method = "gcv"`

* `krls(..., lambda_method = c("loo", "gcv"))` adds a generalized
  cross-validation alternative to the historical leave-one-out
  criterion. GCV minimizes
  `RSS(lambda) / (1 - tr(S(lambda))/n)^2` and is computed in closed
  form from the kernel eigendecomposition (exact path) or from the
  cached SVD of Phi (Nystrom path). Both methods evaluate at the same
  per-iteration cost, so the choice is essentially free. Default is
  `"loo"` -- backward-compatible with all earlier KRLS releases.
* The chosen objective is recorded on the fit object as
  `lambda_method`.

## New: Nystrom scaling vignette

* `vignette("krls-nystrom-scaling", package = "KRLS")` walks through
  the exact-vs-Nystrom runtime comparison, the landmark-reuse pattern
  via `get_landmarks()`, and the LOO-vs-GCV trade-off. Covers what to
  expect when moving past the exact path's ~5,000-row ceiling.

# KRLS 1.4-1

Polish patch for the Nystrom path that landed in 1.4-0. All changes
are additive and non-breaking; existing code paths are unaffected.

## New: `landmark_method = "kmeans"`

* `krls(..., approx = "nystrom", landmark_method = "kmeans")` selects
  landmarks as the cluster centers of a `stats::kmeans()` run on the
  standardized predictors. More robust than random sampling when the
  predictor distribution is imbalanced. Random sampling remains the
  default because it is simpler and bit-reproducible under
  `set.seed()`; k-means initialization is reproducible within an R
  version but not bit-stable across versions.

## New: `get_landmarks()` accessor

* `get_landmarks(fit, scale = c("original", "standardized"))` returns
  the landmark coordinates from a Nystrom fit. The default
  `"original"` scale un-standardizes the stored matrix so it can be
  passed back through `krls(..., landmarks = ...)` without the
  standardize-twice bug. Useful for sensitivity-check refits and
  for inspecting which points the approximation is anchored at.

## Improved diagnostics

* `summary()` on a Nystrom fit now prints the approximation
  configuration (m, landmark method, inference type) and the
  landmark-kernel condition number alongside the count of eigenvalues
  that landed at the relative-ridge floor. When more than half of the
  landmark eigenvalues are floored, an explanatory note is emitted
  pointing at sigma / m / landmark distinctness as the likely cause.
* `glance()` gains `approx`, `nystrom_m`, and `inference` columns.
  Under Nystrom, `eff_df` is now computed correctly from the cached
  Phi spectrum (`Sigma2`) rather than returning `NA`.

## Better input validation

* `krls()` now rejects user-supplied lambda-search windows where
  `L >= U` up front rather than producing a silently invalid
  golden-section sweep.
* User-supplied `landmarks` matrices are checked for NAs, non-finite
  entries, and duplicate rows (which would make the landmark kernel
  rank-deficient).

## Fit object additions (Nystrom only)

* `landmark_method`: `"random"`, `"kmeans"`, `"user_indices"`, or
  `"user_matrix"` recording which path produced the landmarks.
* `Sigma2`, `floored_count`, `D_min_raw`, `D_max_raw`: diagnostic
  scalars from the W eigendecomposition; consumed by `summary()` and
  `glance()` and available to downstream code.

# KRLS 1.4-0

## New: Nystrom approximation with conditional approximate inference

* `krls(..., approx = "nystrom")` adds an explicitly approximate
  low-rank fit path for larger data. The approximation uses landmark
  points, a stabilized Nystrom feature map, and an `m`-space ridge
  solve. Time becomes O(N m^2 + m^3) and memory O(N m), making KRLS
  feasible at sample sizes well beyond the exact path's ~5000-row
  ceiling.
* Five new arguments: `approx`, `nystrom_m` (default
  `min(500, ceiling(sqrt(N)))`), `landmarks`, `landmark_method`
  (`"random"` implemented; `"kmeans"` reserved for a future release),
  `nystrom_eps` (relative-ridge stabilization of W).
* `landmarks` accepts three forms: `NULL` (auto-select),
  an integer vector of row indices into `X`, or an `m` by `D`
  numeric matrix in the original `X` units (standardized internally).
* When `vcov = TRUE`, Nystrom fits report conditional approximate
  standard errors for coefficient weights, predictions, and average
  derivatives. These standard errors condition on the selected
  landmarks, fixed `lambda`, and the low-rank feature approximation
  -- they are not equivalent to exact KRLS standard errors. The full
  `N x N` fitted-value covariance matrix is not stored, preserving
  the memory benefit of the approximation. A new `inference` field
  on the fit object reports `"conditional_nystrom"` or `"none"`.
* `predict(fit, ..., se.fit = TRUE)` works under Nystrom when the fit
  was built with `vcov = TRUE`.
* `summary()`, `tidy()`, and `fdskrls()` (binary first-differences)
  handle Nystrom fits cleanly, with or without standard errors.
* Bootstrap calibration of the AME standard-error formula is included
  in `dev/03_nystrom_bootstrap_validation.R` (developer-side check,
  not shipped to CRAN).

## Faster: O(n^2) AME-variance computation

* The variance of the average marginal effect for each predictor is
  now computed via the algebraic identity `sum(L' V L) = (L 1)' V (L 1)`,
  reducing the per-predictor work from O(n^3) to O(n^2) without
  materializing the full `L' V L` product. The default `krls()` call
  is roughly 1.2x to 3x faster on moderate-to-large fits (n in the
  few hundred to a few thousand range), with the largest gains where
  the variance computation dominated the runtime.
* Numerical results are bit-identical to KRLS 1.2-0 for every fit
  quantity except `var.avgderivatives`, which differs by ~1e-15
  (machine precision; an irreducible side effect of the change in
  summation order). All other quantities — `coeffs`, `fitted`,
  `lambda`, `R2`, `Looe`, `derivatives`, `avgderivatives`, `vcov.c` —
  are bit-for-bit identical.

## Bug fixes

* `krls()` now correctly propagates the user-supplied `tol` argument
  to `lambdasearch()`. Previously, `krls(X, y, tol = ...)` silently
  ignored `tol` and used the default convergence tolerance; calling
  `lambdasearch()` directly worked as documented.
* `solveforc()` and the variance-covariance construction in `krls()`
  no longer error when `eigtrunc` is aggressive enough to keep only
  a single eigenvector. The single-column subset of `Eigenobject$vectors`
  is now extracted with `drop = FALSE`, preserving its matrix shape
  for downstream `multdiag()` / `tcrossprod()`.

# KRLS 1.2-0

## New: formula interface

* `krls(y ~ x1 + x2, data = df)` is now supported alongside the
  matrix interface. The original `krls(X, y)` call continues to work
  unchanged. Drops the intercept column automatically; rejects NAs
  in either side with the same error messages the matrix interface
  uses.

## New: tidyverse-friendly extractors

* `tidy(fit)` returns a per-predictor data frame: `estimate` (the
  average marginal effect), `std.error`, `statistic`, `p.value`,
  `conf.low`/`conf.high`, plus pointwise quartiles
  (`q25`, `median`, `q75`) so users can see the heterogeneity that
  the AME averages over.
* `glance(fit)` is a one-row fit summary: `nobs`, `n_predictors`,
  `r.squared`, leave-one-out MSE, `lambda`, `sigma`, and effective
  degrees of freedom from the kernel ridge.
* `augment(fit, data = ...)` joins `.fitted`, `.resid`, and one
  `.dy_d_<predictor>` column per predictor (the pointwise
  derivatives) back to the data.

  Methods are registered against the `generics` package generics, so
  `library(broom)` makes them discoverable.

## New: `ggplot2` autoplot

* `autoplot(fit)` returns a faceted ggplot of the pointwise
  marginal-effect distribution, one panel per predictor, with the
  AME overlaid in blue. Discoverable via `library(ggplot2)`;
  `ggplot2` is in `Suggests:`.

## New: vignette

* `vignette("krls-quickstart", package = "KRLS")` walks through a
  small simulated example showing the AME-vs-pointwise heterogeneity
  story end-to-end.

# KRLS 1.1-0

## Bug fixes

* `lambdasearch()` now returns a plain scalar instead of a 1x1 matrix
  (it previously used `ifelse()` which inherits the shape of its test
  argument). This eliminates the R 4.4+ deprecation warning
  "Recycling array of length 1 in vector-array arithmetic is
  deprecated" that fired during `Eigenobject$values + lambda` in both
  `solveforc()` and `krls()`'s vcov calculation. The numerical results
  are identical to 1.0-0 within rounding.

* `plot.krls()` no longer dispatches `UseMethod("summary")` on a wrong-
  class input — that was a copy-paste mistake. All three S3 methods
  (`predict`, `summary`, `plot`) now `stop()` cleanly when the input is
  not a `krls` object, with a clear message. The class check itself is
  now `inherits(x, "krls")` rather than `class(x) != "krls"` (the old
  pattern misbehaves when the object has multiple classes).

## Internal cleanups (no user-visible change)

* `krls()`: removed dead code — pre-allocation of `Eigenobject$values`
  and `$vectors` immediately overwritten by `eigen()`, and several
  large commented-out alternative implementations.
* `solveforc()`: trimmed dead alternatives, modernized formatting.
* `predict.krls()`: removed dead commented `if(is.vector(newdata))`
  branch; fixed typo "standart errors" -> "standard errors".
* Redundant `sd(y) == 0` check in `krls()` removed (`var(y) == 0`
  already covers it).

## DESCRIPTION

* Switched to `Authors@R` (Jens Hainmueller as cre+aut, Chad Hazlett
  as aut) per current CRAN style.
* `Imports: grDevices, graphics, stats` added explicitly (these were
  already used via the NAMESPACE; making the dependency declaration
  explicit).
* `Suggests: testthat (>= 3.0.0)` added; `Config/testthat/edition: 3`
  added.
* `URL` field gains the GitHub repo and corrects the Stanford URL.
  `BugReports` field added.
* `Encoding: UTF-8` declared explicitly.
* Description text adds a DOI to the 2014 Political Analysis paper.
