#!/usr/bin/env Rscript
###############################################################################
# nonsp_diagnostic_worker.R
#
# HPC SLURM worker: generate a 2×2 diagnostic compilation PNG for one
# nonsp single-path brms model.  Called from a dSQ job array produced by
# generate_nonsp_diagnostic_jobs.py.  Must be run from HPC_BASE
# (the nonsp_predictor_analyses/ directory), so that relative paths
# file.path(predictor, model, ...) resolve correctly.
#
# Arguments (positional, matching nonsp_predictors.R conventions):
#   1  settings   — brms family string (e.g. hurdle_negbinom_huvary, bernoulli)
#   2  model      — covariate-set / model-type key (e.g. nice_covariates_spusers)
#   3  dv         — DV column name (e.g. caps_vision)
#   4  predictor  — normalized spvar column name (e.g. avg_life_dose_normalized)
#
# Output:
#   {predictor}/{model}/results/diagnostics/{dv}_diagnostic_compilation.png
#
# Panel layout (2 × 2):
#   (a) top-left    — posterior predictive check
#   (b) top-right   — MCMC trace (spvar)
#   (c) bottom-left — DHARMa comprehensive: QQ-uniform + residuals vs fitted
#   (d) bottom-right— DHARMa residuals vs spvar (heteroscedasticity check)
#
# The DHARMa comprehensive (panel c) is generated fresh from the fit; all other
# panels load existing PNGs saved by nonsp_predictors.R.
#
# Staleness check: skips if the compilation PNG already exists AND is newer
# than all inputs (fit .RData, pp_check, trace, heteroscedasticity PNGs).
# If any input was regenerated (e.g. after re-running nonsp_predictors.R),
# the compilation is rebuilt automatically.
#
# Requires: brms, DHARMa, DHARMa.helpers, ggplot2, cowplot, png
###############################################################################

# ── R library path (mirrors nonsp_predictors.R setup) ────────────────────────
r_minor    <- sub("\\..*$", "", R.version$minor)
user_r_lib <- path.expand(file.path("~", "R", paste0(R.version$major, ".", r_minor)))
if (!dir.exists(user_r_lib)) dir.create(user_r_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(user_r_lib, .libPaths()))

suppressPackageStartupMessages({
  library(brms)
  library(DHARMa)
  library(DHARMa.helpers)
  library(ggplot2)
  library(cowplot)
  library(png)
})

# ── Parse arguments ───────────────────────────────────────────────────────────
args      <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4L) {
  stop("Usage: Rscript nonsp_diagnostic_worker.R <settings> <model> <dv> <predictor>")
}
settings  <- args[1]
model     <- args[2]
dv        <- args[3]
predictor <- args[4]  # normalized spvar column name

cat(sprintf("\n[nonsp_diagnostic_worker] settings=%s | model=%s | dv=%s | predictor=%s\n",
            settings, model, dv, predictor))

# ── Paths (all relative to HPC_BASE, i.e. the submission directory) ──────────
base_dir         <- file.path(predictor, model)
pp_path          <- file.path(base_dir, "pp_checks",
                               paste0(dv, "_posterior_check_plot.png"))
trace_path       <- file.path(base_dir, "convergence_tests",
                               paste0(dv, "_convergence_plot.png"))
hetero_path      <- file.path(base_dir, "results", "diagnostics",
                               paste0(dv, "_heteroscedasticity_check.png"))
dharma_comp_path <- file.path(base_dir, "results", "diagnostics",
                               paste0(dv, "_dharma_comprehensive.png"))
out_path         <- file.path(base_dir, "results", "diagnostics",
                               paste0(dv, "_diagnostic_compilation.png"))
fit_path         <- file.path(base_dir, "results", "fits",
                               paste0(dv, "_fit.RData"))

# Skip if compilation already exists AND no input file is newer (stale check).
# If the fit .RData or any source PNG was regenerated after the compilation,
# the compilation is stale and must be rebuilt.
if (file.exists(out_path)) {
  out_mtime   <- file.mtime(out_path)
  input_files <- c(fit_path, pp_path, trace_path, hetero_path)
  existing    <- input_files[file.exists(input_files)]
  if (length(existing) > 0 && any(file.mtime(existing) > out_mtime)) {
    cat("  Output exists but inputs are newer — regenerating.\n")
  } else {
    cat("  Output already exists and is up to date — skipping.\n")
    quit(status = 0)
  }
}

# Verify fit file exists
if (!file.exists(fit_path)) {
  cat(sprintf("  ERROR: Fit file not found: %s\n", fit_path))
  quit(status = 1)
}

# ── Load fit ─────────────────────────────────────────────────────────────────
cat(sprintf("  Loading: %s\n", fit_path))
e    <- new.env(parent = emptyenv())
vars <- tryCatch(load(fit_path, envir = e), error = function(e) character(0))
var  <- intersect(vars, c("fit1", "fit"))
if (length(var) == 0L) {
  cat("  ERROR: RData does not contain fit1 or fit.\n")
  quit(status = 1)
}
fit1 <- get(var[1], envir = e)
rm(e); invisible(gc(verbose = FALSE))
cat("  Fit loaded.\n")

# ── DHARMa: integer flag ─────────────────────────────────────────────────────
# These settings produce integer-valued responses where integer=TRUE is required.
INTEGER_SETTINGS <- c(
  "hurdle_negbinom_huvary", "hurdle_negbinomial",
  "negbinomial", "zero_negbinomial", "zero_inflated_poisson",
  "bernoulli", "logistic"
)
integer_flag <- settings %in% INTEGER_SETTINGS
cat(sprintf("  DHARMa integer=%s (settings=%s)\n", integer_flag, settings))

# ── Generate DHARMa residuals ─────────────────────────────────────────────────
sim_res <- NULL
tryCatch({
  sim_res <- DHARMa.helpers::dh_check_brms(fit1, integer = integer_flag, plot = FALSE)
  cat("  DHARMa residuals generated.\n")
}, error = function(e) {
  cat(sprintf("  WARNING: DHARMa residual generation failed: %s\n", e$message))
})

# ── Save DHARMa comprehensive figure (QQ-uniform + residuals vs fitted) ──────
if (!is.null(sim_res)) {
  tryCatch({
    png(dharma_comp_path, width = 10, height = 5, units = "in", res = 150)
    plot(sim_res)
    dev.off()
    cat(sprintf("  Saved DHARMa comprehensive: %s\n", basename(dharma_comp_path)))
  }, error = function(e) {
    if (dev.cur() > 1) dev.off()
    cat(sprintf("  WARNING: DHARMa comprehensive plot failed: %s\n", e$message))
    if (file.exists(dharma_comp_path)) file.remove(dharma_comp_path)
  })
}

# ── Helper: load PNG as ggplot panel, with placeholder on failure ─────────────
# Uses annotation_raster() to embed the PNG into a blank ggplot coordinate system.
load_png_panel <- function(path) {
  if (!file.exists(path)) {
    return(
      ggplot() +
        annotate("text", x = 0.5, y = 0.5, hjust = 0.5, vjust = 0.5,
                 label = paste0("Not found:\n", basename(path)),
                 size = 3.5, color = "firebrick") +
        theme_void() +
        theme(plot.background = element_rect(fill = "white",
                                             color = "grey70", linewidth = 0.5))
    )
  }
  img <- tryCatch(png::readPNG(path), error = function(e) NULL)
  if (is.null(img)) {
    return(
      ggplot() +
        annotate("text", x = 0.5, y = 0.5, hjust = 0.5, vjust = 0.5,
                 label = paste0("Read error:\n", basename(path)),
                 size = 3.5, color = "firebrick") +
        theme_void()
    )
  }
  ggplot() +
    xlim(0, 1) + ylim(0, 1) +
    annotation_raster(img, xmin = 0, xmax = 1, ymin = 0, ymax = 1) +
    theme_void() +
    theme(plot.background = element_rect(fill = "white", color = NA))
}

# ── Assemble 2 × 2 grid ───────────────────────────────────────────────────────
panel_a <- load_png_panel(pp_path)
panel_b <- load_png_panel(trace_path)
panel_c <- load_png_panel(dharma_comp_path)
panel_d <- load_png_panel(hetero_path)

grid_4 <- cowplot::plot_grid(
  panel_a, panel_b,
  panel_c, panel_d,
  ncol   = 2,
  nrow   = 2,
  labels = c("a  PP check", "b  Trace", "c  DHARMa", "d  DHARMa vs predictor"),
  label_size       = 10,
  label_fontfamily = "sans",
  label_x     = 0.01,
  label_y     = 0.99,
  hjust       = 0,
  vjust       = 1
)

ggplot2::ggsave(out_path, grid_4, width = 16, height = 12, dpi = 120)
cat(sprintf("  Saved: %s\n", out_path))
cat("[nonsp_diagnostic_worker] Done.\n")
