################################################################################
# gelman_normalization.R — R-side Gelman (2 SD) normalization for HPC jobs
################################################################################
#
# WHY THIS FILE EXISTS
# --------------------
# Until 2026-08-22 the master CSV shipped to the cluster was normalized in
# Python by generate_hpc_jobs.py, but the row-level subsetting keywords
# ("spusers", "iqr", ...) are applied HERE, in R, after the CSV is read.  Every
# subset model was therefore fitted on predictors centred and scaled on a
# LARGER sample than the one being fitted: mean was not 0 and sd was not 0.5 on
# the rows that actually entered the likelihood.
#
# The fix is to ship raw values and normalize in R, last, so the normalization
# sample and the estimation sample are the same set of rows.  The pipeline order
# in both nonsp_predictors.R and hpc_mediation.R is now:
#
#     read raw CSV
#       -> row-level subsetting keywords          (spusers, iqr, ...)
#       -> drop rows incomplete on any model term (drop_incomplete_model_rows)
#       -> gelman_normalize_df()                  <- this file
#       -> verify_gelman_normalization()
#       -> as.factor() / as.ordered() conversions
#       -> brm()
#
# WHAT IT IS A PORT OF
# --------------------
# modules/data_prep.py :: normalize_analysis_df() / verify_normalization().
# The four rules, their order, the "skip columns already at mean 0 / sd 0.5"
# guard and the NaN handling are deliberately identical, so that running
# gelman_normalize_df() on the FULL unsubset dataframe reproduces the Python
# output column for column.  That equivalence is the regression test for this
# file.  If you change a rule here, change it in the Python half too, or the two
# halves of the pipeline diverge.
#
# WHERE THE VARIABLE LISTS COME FROM
# ----------------------------------
# NOT from this file.  They arrive in normalization_vars.R, which
# generate_hpc_jobs.py auto-generates from VARIABLE_REGISTRY and
# BASE_COVARIATE_SETS in modules/master_config.py.  This file contains only the
# transform; that file contains only the membership.  Do not hardcode a variable
# name here.
#
# GOTCHA
# ------
# "Gelman normalization" divides by 2 SD, not 1 SD.  A correctly normalized
# column has mean 0 and sd 0.5 — NOT sd 1.
################################################################################


# A Gelman-normalized column has mean 0 and sd 0.5.  That pair is exactly the
# fixed point of the transform, so "normalize only if not already at it"
# produces the same values as normalizing unconditionally, without re-running
# the transform and accumulating floating-point drift.  Matches GELMAN_ATOL in
# modules/data_prep.py.
GELMAN_ATOL <- 1e-9


#' Is this column already Gelman 2SD normalized (mean 0, sd 0.5)?
#'
#' Returns FALSE for all-NA and zero-variance columns: those cannot be
#' normalized, so reporting them as normalized would hide the reason they were
#' skipped.  Port of data_prep.py :: is_gelman_normalized().
is_gelman_normalized <- function(x, atol = GELMAN_ATOL) {
  s <- x[!is.na(x)]
  if (length(s) < 2L) return(FALSE)
  s_sd <- stats::sd(s)
  if (is.na(s_sd) || s_sd == 0) return(FALSE)
  abs(mean(s)) < atol && abs(s_sd - 0.5) < atol
}


#' Gelman 2SD normalize the named columns of `df`, in place.
#'
#' Port of data_prep.py :: gelman_standardize().  Columns are skipped — with a
#' message, never silently — when they are absent, when they are in `skip_vars`
#' (brms factor-style binaries, kept unscaled), or when mean/sd is NA or the
#' variance is zero.
#'
#' `stats::sd()` and pandas `Series.std()` both use the n-1 denominator, and
#' both statistics here are computed over non-missing values only, so the two
#' implementations agree.
gelman_standardize_cols <- function(df, cols, skip_vars = character(0), verbose = TRUE) {
  for (col in cols) {
    if (!(col %in% names(df))) {
      if (verbose) message(paste0("  [gelman] skipping ", col, " (missing column)"))
      next
    }
    if (col %in% skip_vars) {
      if (verbose) message(paste0("  [gelman] leaving binary column unscaled: ", col))
      next
    }
    x <- df[[col]]
    mean_val <- mean(x, na.rm = TRUE)
    sd_val   <- stats::sd(x, na.rm = TRUE)
    if (is.na(mean_val) || is.na(sd_val) || sd_val == 0) {
      if (verbose) message(paste0("  [gelman] skipping ", col, " (NA mean/sd or zero variance)"))
      next
    }
    df[[col]] <- (x - mean_val) / (2 * sd_val)
  }
  df
}


#' Gelman-normalize only the columns not already at mean 0 / sd 0.5.
#' Returns list(df = <df>, done = <character vector of columns transformed>).
.norm_if_needed <- function(df, cols, skip_vars, verbose) {
  todo <- cols[cols %in% names(df)]
  todo <- todo[!vapply(todo, function(c) is_gelman_normalized(df[[c]]), logical(1))]
  if (length(todo) > 0) df <- gelman_standardize_cols(df, todo, skip_vars, verbose)
  list(df = df, done = todo)
}


#' Apply the project's normalization conventions to an analysis dataframe.
#'
#' Port of data_prep.py :: normalize_analysis_df().  Four rules, applied in
#' order — the order matters, because rules 3 and 4 rely on rules 1 and 2 having
#' already created and normalized their columns:
#'
#'   1. `need_non_normalized` — the brms family constrains the raw scale (gamma,
#'      beta, zero_inflated_beta, hurdle/negbinomial).  The raw column is left
#'      untouched and a separate `{col}_normalized` column is created.
#'   2. `inplace_normalized` — the family is an unbounded real (student_t /
#'      gaussian).  The raw column is normalized and `{col}_normalized` is an
#'      alias of it, so the same name is valid in the mediator formula and the
#'      DV formula.
#'   3. `normalize_in_place` — remaining continuous covariates and main
#'      predictors are normalized in place.  `categorical` variables are
#'      excluded (R applies as.factor()), as are `monotonic` ones (mo() needs
#'      raw integers).
#'   4. `ensure_normalized_copy` — every job predictor is guaranteed a
#'      `{col}_normalized` column.  Predictors already normalized by rules 1-3
#'      are aliased; any others keep their raw column and receive a normalized
#'      copy.
#'
#' @param skip_vars Columns never scaled even when named by a rule.  This
#'   mirrors BINARY_FACTOR_VARS in data_prep.py, which is a SEPARATE and shorter
#'   list than CATEGORICAL_FACTOR_VARS in master_config.py.  The difference is
#'   inert today (no categorical variable appears in any rule's list) but the
#'   shorter list is passed through deliberately so this port cannot diverge
#'   from the Python it is checked against.
gelman_normalize_df <- function(df,
                                normalize_in_place     = character(0),
                                ensure_normalized_copy = character(0),
                                need_non_normalized    = character(0),
                                inplace_normalized     = character(0),
                                categorical            = character(0),
                                monotonic              = character(0),
                                skip_vars              = character(0),
                                verbose                = TRUE) {

  # Rule 1 — raw column preserved, normalized values in a separate column.
  need_present <- need_non_normalized[need_non_normalized %in% names(df)]
  for (col in need_present) {
    norm_col <- paste0(col, "_normalized")
    if (!(norm_col %in% names(df))) df[[norm_col]] <- df[[col]]
  }
  r1 <- .norm_if_needed(df, paste0(need_present, "_normalized"), skip_vars, verbose)
  df <- r1$df

  # Rule 2 — raw column normalized; `{col}_normalized` aliases it.
  inplace_present <- inplace_normalized[inplace_normalized %in% names(df)]
  r2 <- .norm_if_needed(df, inplace_present, skip_vars, verbose)
  df <- r2$df
  for (col in inplace_present) df[[paste0(col, "_normalized")]] <- df[[col]]

  # Rule 3 — everything else brms will see as a continuous term.
  in_place <- sort(setdiff(unique(normalize_in_place),
                           c(categorical, monotonic,
                             need_non_normalized, inplace_normalized)))
  r3 <- .norm_if_needed(df, in_place, skip_vars, verbose)
  df <- r3$df

  # Rule 4 — guarantee a `{col}_normalized` column for every job predictor.
  copies_made <- character(0)
  for (col in ensure_normalized_copy) {
    if (!(col %in% names(df))) {
      stop(paste0("gelman_normalize_df(): predictor '", col, "' not found in dataframe"))
    }
    if (!is.numeric(df[[col]])) {
      stop(paste0("gelman_normalize_df(): non-numeric predictor '", col, "'"))
    }
    norm_col <- paste0(col, "_normalized")
    if (norm_col %in% names(df)) next
    df[[norm_col]] <- df[[col]]
    if (!(col %in% in_place)) {
      r4 <- .norm_if_needed(df, norm_col, skip_vars, verbose)
      df <- r4$df
    }
    copies_made <- c(copies_made, col)
  }

  if (verbose) {
    message(sprintf("  [normalize] rule 1 need_non_normalized : %d cols (%d normalized, rest already were)",
                    length(need_present), length(r1$done)))
    message(sprintf("  [normalize] rule 2 inplace_normalized  : %d cols (%d normalized, rest already were)",
                    length(inplace_present), length(r2$done)))
    message(sprintf("  [normalize] rule 3 in-place            : %d cols (%d normalized, rest already were)",
                    length(in_place), length(r3$done)))
    message(sprintf("  [normalize] rule 4 predictor copies    : %d cols", length(copies_made)))
  }

  df
}


#' Assert the normalization invariants declared in VARIABLE_REGISTRY.
#'
#' Port of data_prep.py :: verify_normalization().  Stops, listing every
#' violation, rather than fitting a model on a mis-scaled dataframe.
#'
#' A raw need_non_normalized column that looks normalized is the serious case:
#' those variables carry a gamma, beta, or hurdle/negbinomial family whose
#' support excludes the negative values normalization introduces, so the column
#' is unusable as a brms response and the original values are no longer
#' recoverable from the dataframe.
verify_gelman_normalization <- function(df,
                                        need_non_normalized = character(0),
                                        inplace_normalized  = character(0),
                                        atol    = GELMAN_ATOL,
                                        verbose = TRUE) {
  problems <- character(0)
  checked  <- 0L

  for (col in need_non_normalized) {
    if (!(col %in% names(df))) next
    checked  <- checked + 1L
    norm_col <- paste0(col, "_normalized")
    if (!(norm_col %in% names(df))) {
      problems <- c(problems, paste0(norm_col, ": missing (need_non_normalized requires it)"))
    } else if (!is_gelman_normalized(df[[norm_col]], atol)) {
      s <- df[[norm_col]][!is.na(df[[norm_col]])]
      problems <- c(problems, sprintf("%s: not normalized (mean=%.6g, sd=%.6g; expected mean 0, sd 0.5)",
                                      norm_col, mean(s), stats::sd(s)))
    }
    if (is_gelman_normalized(df[[col]], atol)) {
      problems <- c(problems, paste0(
        col, ": RAW COLUMN IS NORMALIZED — its brms family requires the raw scale, ",
        "and the raw values are no longer present in the dataframe"))
    }
  }

  for (col in inplace_normalized) {
    if (!(col %in% names(df))) next
    checked <- checked + 1L
    for (pair in list(c(col, "raw"), c(paste0(col, "_normalized"), "alias"))) {
      cname <- pair[1]; role <- pair[2]
      if (!(cname %in% names(df))) {
        problems <- c(problems, paste0(cname, ": missing (inplace_normalized requires it)"))
        next
      }
      if (!is_gelman_normalized(df[[cname]], atol)) {
        s <- df[[cname]][!is.na(df[[cname]])]
        problems <- c(problems, sprintf("%s: %s column not normalized (mean=%.6g, sd=%.6g; expected mean 0, sd 0.5)",
                                        cname, role, mean(s), stats::sd(s)))
      }
    }
  }

  if (length(problems) > 0) {
    stop(paste0("Normalization verification failed:\n  - ",
                paste(problems, collapse = "\n  - ")))
  }
  if (verbose) {
    message(paste0("  [normalize] verification passed on ", checked,
                   " registry-declared variables"))
  }
  invisible(TRUE)
}


#' Split brms formula RHS strings into the bare variable names they reference.
#'
#' Handles the term syntax this pipeline actually emits: "+" separated terms,
#' interactions written with "*" or ":", and monotonic covariates wrapped in
#' mo().  Intercept literals ("1", "0") and empty strings are dropped.
#'
#' Used to work out which columns a model depends on, so rows incomplete on any
#' of them can be dropped BEFORE normalization — see drop_incomplete_model_rows().
model_formula_vars <- function(...) {
  parts <- unlist(list(...), use.names = FALSE)
  parts <- parts[!is.na(parts) & nzchar(parts)]
  if (length(parts) == 0) return(character(0))
  terms <- unlist(strsplit(parts, "[+*:]"), use.names = FALSE)
  terms <- trimws(terms)
  # Unwrap mo(x) — and any other single-argument wrapper of the same shape.
  terms <- sub("^[A-Za-z_][A-Za-z0-9_.]*\\((.*)\\)$", "\\1", terms)
  terms <- trimws(terms)
  terms <- terms[nzchar(terms) & !(terms %in% c("0", "1"))]
  unique(terms)
}


#' Drop rows that are missing any variable the model will use.
#'
#' brms drops incomplete cases itself, but it does so AFTER this script has
#' normalized — which would leave the normalization sample a strict superset of
#' the estimation sample and reintroduce the very off-centring this reordering
#' exists to fix.  Doing it here makes the two sets identical, so mean 0 / sd
#' 0.5 holds on exactly the rows that enter the likelihood.
#'
#' `_normalized` columns do not exist yet at this point in the script (they are
#' created by gelman_normalize_df() further down), so a variable named
#' `{x}_normalized` is resolved to its raw column `{x}`.  This is safe because
#' normalization is row-wise: `{x}_normalized` is missing exactly where `{x}` is.
#' The mapping is reported per variable, never applied silently, and a name that
#' resolves to neither column is a hard stop rather than a skipped filter.
drop_incomplete_model_rows <- function(df, vars, verbose = TRUE) {
  vars <- unique(vars[!is.na(vars) & nzchar(vars)])
  resolved <- character(0)
  missing_cols <- character(0)

  for (v in vars) {
    if (v %in% names(df)) {
      resolved <- c(resolved, v)
      next
    }
    base <- sub("_normalized$", "", v)
    if (base != v && base %in% names(df)) {
      if (verbose) {
        message(paste0("  [complete-cases] '", v, "' not yet created; ",
                       "using raw column '", base, "' (identical missingness)"))
      }
      resolved <- c(resolved, base)
      next
    }
    missing_cols <- c(missing_cols, v)
  }

  if (length(missing_cols) > 0) {
    stop(paste0(
      "drop_incomplete_model_rows(): model variable(s) not found in the dataframe: ",
      paste(missing_cols, collapse = ", "),
      "\nThe master CSV does not contain every term this model's formula references. ",
      "Re-run generate_hpc_jobs.py and re-transfer rather than proceeding."))
  }

  resolved <- unique(resolved)
  n_before <- nrow(df)
  keep     <- stats::complete.cases(df[, resolved, drop = FALSE])
  df       <- df[keep, , drop = FALSE]

  if (verbose) {
    message(paste0("  [complete-cases] model variables (", length(resolved), "): ",
                   paste(sort(resolved), collapse = ", ")))
    message(paste0("  [complete-cases] ", nrow(df), " rows (removed ",
                   n_before - nrow(df), " incomplete of ", n_before, ")"))
  }
  if (nrow(df) == 0) {
    stop("drop_incomplete_model_rows(): no rows remain after dropping incomplete cases.")
  }
  df
}


#' Resolve a column name to the column that actually exists in a raw dataframe.
#'
#' The job files address the focal predictor by its `{x}_normalized` name, but
#' the row filters run BEFORE normalization, when only `{x}` exists.  This maps
#' the one to the other, reporting the substitution rather than making it
#' silently, and stopping if neither name is present.
#'
#' Only safe for filters that are invariant to a positive affine transform of the
#' column — the IQR fence is (quantiles are equivariant, so the same rows fall
#' inside it either way).  A filter comparing against an absolute threshold on
#' the normalized scale is NOT invariant, which is why the "nooutlier" keyword
#' was retired rather than ported.
resolve_raw_column <- function(df, col, verbose = TRUE) {
  if (col %in% names(df)) return(col)
  base <- sub("_normalized$", "", col)
  if (base != col && base %in% names(df)) {
    if (verbose) {
      message(paste0("  [resolve] '", col, "' not yet created; using raw column '",
                     base, "'"))
    }
    return(base)
  }
  stop(paste0("resolve_raw_column(): neither '", col, "' nor '", base,
              "' is present in the dataframe."))
}
