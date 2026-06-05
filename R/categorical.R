## Categorical preprocessing for krls()
##
## Convention matches kbal / GPSS:
##
##   - Each declared categorical column is one-hot encoded with ALL
##     levels (contrasts = FALSE, no reference cell), then multiplied
##     by sqrt(0.5). The factor of sqrt(1/2) compensates for the
##     doubled squared-Euclidean distance from carrying both the
##     "in level k" and "not in level k" indicators: between two
##     observations differing in one category, the raw one-hot
##     contributes 2 to ||x_i - x_j||^2; sqrt(0.5) brings that back
##     to a Hamming-like 1.
##
##   - Continuous columns are standardized to sd = 1 (as in the
##     existing krls() default).
##
## There is no autodetection. Users must list categorical columns
## explicitly in `cat_columns`. .warn_unmarked_categoricals() scans
## the supplied X for columns that look categorical but were not
## listed, and emits a warning -- mirroring the friendly nudge in
## kbal/GPSS while keeping the API explicit.

.ONEHOT_SCALE <- sqrt(0.5)

## resolve_cat_columns(): turn user-supplied character or numeric
## indices into a canonical integer index vector against ncol(X).
.resolve_cat_columns <- function(cat_columns, X) {
  if (is.null(cat_columns) || length(cat_columns) == 0L)
    return(integer(0))
  d <- ncol(X)
  cn <- colnames(X)
  if (is.character(cat_columns)) {
    if (is.null(cn))
      stop("cat_columns is character but X has no column names")
    bad <- setdiff(cat_columns, cn)
    if (length(bad))
      stop("cat_columns not found in X: ", paste(bad, collapse = ", "))
    idx <- match(cat_columns, cn)
  } else if (is.numeric(cat_columns)) {
    if (any(cat_columns != as.integer(cat_columns)) ||
        any(cat_columns < 1L) || any(cat_columns > d))
      stop("numeric cat_columns must be integer indices in 1:ncol(X)")
    idx <- as.integer(cat_columns)
  } else {
    stop("cat_columns must be a character vector of column names ",
         "or an integer vector of column indices")
  }
  if (anyDuplicated(idx))
    stop("cat_columns contains duplicate entries")
  sort(idx)
}

## warn_unmarked_categoricals(): heuristic check for likely
## categorical columns the user did not declare. Triggers on:
##   - factor or character columns
##   - logical columns
##   - numeric columns with <= max_unique distinct values
## Documented as a friendly nudge, not a guarantee.
.warn_unmarked_categoricals <- function(X, cat_idx, max_unique = 10L) {
  d <- ncol(X)
  cn <- colnames(X)
  flagged <- character(0)
  reasons <- character(0)
  for (j in seq_len(d)) {
    if (j %in% cat_idx) next
    col <- X[, j]
    nm <- if (!is.null(cn)) cn[j] else paste0("col ", j)
    if (is.factor(col) || is.character(col)) {
      flagged <- c(flagged, nm)
      reasons <- c(reasons, "factor/character")
    } else if (is.logical(col)) {
      flagged <- c(flagged, nm)
      reasons <- c(reasons, "logical")
    } else if (is.numeric(col)) {
      u <- length(unique(stats::na.omit(col)))
      if (u <= max_unique) {
        flagged <- c(flagged, nm)
        reasons <- c(reasons, sprintf("%d unique values", u))
      }
    }
  }
  if (length(flagged)) {
    msg <- paste0(
      "krls: the following column(s) look categorical but were not ",
      "listed in cat_columns:\n  ",
      paste(sprintf("%s (%s)", flagged, reasons), collapse = ", "),
      "\nKRLS does NOT autodetect categoricals. Pass these in ",
      "cat_columns to get the sqrt(0.5)-scaled one-hot encoding ",
      "(matching kbal / GPSS), or pass cat_columns = character(0) ",
      "/ integer(0) to silence this warning.")
    warning(msg, call. = FALSE)
  }
  invisible(NULL)
}

## one_hot(): kbal-style one-hot expansion. Keeps every level
## (contrasts = FALSE) and strips the intercept column. Returns a
## numeric matrix; preserves column names that encode level identity.
.one_hot <- function(df) {
  df <- as.data.frame(lapply(as.data.frame(df, stringsAsFactors = FALSE),
                             as.factor),
                      stringsAsFactors = FALSE)
  if (any(vapply(df, function(v) length(levels(v)) < 2L, logical(1)))) {
    bad <- names(df)[vapply(df, function(v) length(levels(v)) < 2L,
                            logical(1))]
    stop("categorical column(s) with fewer than 2 levels: ",
         paste(bad, collapse = ", "))
  }
  mm <- stats::model.matrix(
    ~ ., data = df,
    contrasts.arg = lapply(df, stats::contrasts, contrasts = FALSE))
  mm[, colnames(mm) != "(Intercept)", drop = FALSE]
}

## preprocess_X(): the workhorse. Splits X into categorical and
## continuous columns, one-hot encodes the categoricals with the
## sqrt(0.5) scale, standardizes the continuous columns to sd=1, and
## returns the concatenated numeric matrix plus the metadata needed to
## apply the same transformation at predict time.
##
## Returns:
##   list(
##     X_proc           = numeric matrix (n x d_proc),
##     cont_centers     = named numeric vector (for predict),
##     cont_scales      = named numeric vector (for predict),
##     cat_levels       = list of factor levels per original cat column,
##     cat_columns_orig = original cat column indices (for predict),
##     col_origin       = integer vector mapping each processed column
##                        back to its original column (for cat_columns
##                        bookkeeping and AME aggregation),
##     d_orig           = ncol(X) before processing
##   )
##
## With cat_columns = integer(0) this is a near-identity passthrough:
## continuous standardize-to-sd=1 on every column (matching the
## existing krls() default). The categorical machinery only kicks in
## when the user opts into it.
.preprocess_X <- function(X, cat_columns = NULL,
                          warn_unmarked = TRUE) {
  if (is.data.frame(X)) {
    # Allow data.frame input so that factor/character cat columns can
    # come through untouched. We keep them as columns until one_hot()
    # converts them.
    X_df_in <- X
  } else {
    X_df_in <- as.data.frame(X, stringsAsFactors = FALSE)
  }
  d_orig <- ncol(X_df_in)
  cn <- colnames(X_df_in)
  cat_idx <- .resolve_cat_columns(cat_columns, X_df_in)

  if (warn_unmarked)
    .warn_unmarked_categoricals(X_df_in, cat_idx)

  cont_idx <- setdiff(seq_len(d_orig), cat_idx)

  # Continuous block: numeric matrix, standardize to sd=1 (kbal default).
  if (length(cont_idx)) {
    cont_block_raw <- X_df_in[, cont_idx, drop = FALSE]
    if (any(vapply(cont_block_raw, function(v) !is.numeric(v),
                   logical(1)))) {
      bad <- names(cont_block_raw)[
        vapply(cont_block_raw, function(v) !is.numeric(v), logical(1))]
      stop("continuous column(s) are not numeric (declare them in ",
           "cat_columns?): ", paste(bad, collapse = ", "))
    }
    cont_mat <- as.matrix(cont_block_raw)
    storage.mode(cont_mat) <- "double"
    cont_centers <- colMeans(cont_mat)
    cont_scales  <- apply(cont_mat, 2L, stats::sd)
    if (any(cont_scales == 0))
      stop("continuous column(s) with zero variance: ",
           paste(names(cont_scales)[cont_scales == 0], collapse = ", "))
    cont_std <- scale(cont_mat, center = cont_centers, scale = cont_scales)
    attr(cont_std, "scaled:center") <- NULL
    attr(cont_std, "scaled:scale")  <- NULL
    cont_origin <- cont_idx
  } else {
    cont_std <- matrix(numeric(0), nrow = nrow(X_df_in), ncol = 0)
    cont_centers <- numeric(0)
    cont_scales  <- numeric(0)
    cont_origin  <- integer(0)
  }

  # Categorical block: one-hot with all levels, scaled by sqrt(0.5).
  if (length(cat_idx)) {
    cat_block_raw <- X_df_in[, cat_idx, drop = FALSE]
    cat_levels <- lapply(as.data.frame(cat_block_raw, stringsAsFactors = FALSE),
                         function(v) levels(as.factor(v)))
    names(cat_levels) <- if (!is.null(cn)) cn[cat_idx] else
      paste0("V", cat_idx)
    cat_mat <- .one_hot(cat_block_raw) * .ONEHOT_SCALE
    # Track which original column each expanded one-hot column
    # came from. .one_hot() concatenates blocks in cat_idx order; we
    # rebuild the mapping from level counts.
    n_levels <- vapply(cat_levels, length, integer(1))
    cat_origin_expanded <- rep(cat_idx, times = n_levels)
  } else {
    cat_mat <- matrix(numeric(0), nrow = nrow(X_df_in), ncol = 0)
    cat_levels <- list()
    cat_origin_expanded <- integer(0)
  }

  X_proc <- cbind(cont_std, cat_mat)
  col_origin <- c(cont_origin, cat_origin_expanded)
  storage.mode(X_proc) <- "double"

  list(X_proc            = X_proc,
       cont_centers      = cont_centers,
       cont_scales       = cont_scales,
       cont_idx          = cont_idx,
       cat_idx           = cat_idx,
       cat_levels        = cat_levels,
       col_origin        = col_origin,
       d_orig            = d_orig,
       onehot_scale      = .ONEHOT_SCALE)
}

## apply_preprocess_X(): re-apply the training-time transformation to
## a new X (predict path). Uses the cached centers/scales for
## continuous columns and the cached level set for one-hots. Unseen
## levels in newdata raise a hard error -- silently dropping a level
## or coercing it to NA at predict time would shift kernel distances
## without the user noticing.
.apply_preprocess_X <- function(X_new, prep) {
  if (is.data.frame(X_new)) {
    X_df <- X_new
  } else {
    X_df <- as.data.frame(X_new, stringsAsFactors = FALSE)
  }
  if (ncol(X_df) != prep$d_orig)
    stop("newdata has ", ncol(X_df), " columns but training X had ",
         prep$d_orig)

  # Continuous part
  if (length(prep$cont_idx)) {
    cont_block <- as.matrix(X_df[, prep$cont_idx, drop = FALSE])
    storage.mode(cont_block) <- "double"
    cont_std <- scale(cont_block, center = prep$cont_centers,
                      scale = prep$cont_scales)
    attr(cont_std, "scaled:center") <- NULL
    attr(cont_std, "scaled:scale")  <- NULL
  } else {
    cont_std <- matrix(numeric(0), nrow = nrow(X_df), ncol = 0)
  }

  # Categorical part: factor with frozen levels, then one-hot.
  if (length(prep$cat_idx)) {
    cat_block_raw <- X_df[, prep$cat_idx, drop = FALSE]
    cat_levels_train <- prep$cat_levels
    cat_block_factored <- as.data.frame(
      Map(function(v, lev) {
        v_chr <- as.character(v)
        bad <- setdiff(unique(v_chr), lev)
        if (length(bad))
          stop("newdata contains unseen level(s) in categorical column: ",
               paste(bad, collapse = ", "))
        factor(v_chr, levels = lev)
      }, cat_block_raw, cat_levels_train),
      stringsAsFactors = FALSE)
    names(cat_block_factored) <- names(cat_levels_train)
    cat_mat <- .one_hot(cat_block_factored) * prep$onehot_scale
  } else {
    cat_mat <- matrix(numeric(0), nrow = nrow(X_df), ncol = 0)
  }

  out <- cbind(cont_std, cat_mat)
  storage.mode(out) <- "double"
  out
}
