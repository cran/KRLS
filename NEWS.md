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
