################################################################################
# _compile_diagnostics_helper.R
#
# Generates all diagnostic figures for create_mediation_diagnostic_compilation.py.
# Called as a subprocess with command-line arguments; saves "compiled_*" PNGs
# to the specified output directory.
#
# All figures use white/light backgrounds so they render correctly in the
# compiled supplementary figure.
#
# Arguments (positional):
#   [1]  fit_path       — path to .RData file containing brms 'fit' object
#   [2]  out_dir        — directory in which to save compiled_*.png files
#   [3]  spvar          — spvar column name as it appears in fit$data
#   [4]  spvar_label    — human-readable label (from VARIABLE_REGISTRY)
#   [5]  mediator       — raw mediator column (used as response in mediator formula)
#   [6]  mediator_label — human-readable label
#   [7]  mediator_in_dv — mediator column used as *predictor* in DV formula
#                         (same as mediator if inplace_normalized; + "_normalized" if not)
#   [8]  dv             — DV column name
#   [9]  dv_label       — human-readable label
#   [10] W2             — figure width  for 2-column panels (default 8.5 inches)
#   [11] H2             — figure height for 2-column panels (default 5.2 inches)
#   [12] W3             — figure width  for 3-column panels (default 5.5 inches)
#   [13] H3             — figure height for 3-column panels at full height (default 4.1 inches)
#   [14] DPI            — output DPI (default 150)
#   [15] H3_HALF        — half height for split trace rows (default H3/2 = 2.05 inches)
#
# Trace splitting logic:
#   mu traces (compiled_traces_mu_*): match ^b_{resp_clean}_ (exact, no submodel prefix).
#     Generated at H3_HALF for dv=="caps_vision" (displayed in split top row);
#     at H3 for all other DVs (single full-height trace row).
#   hu traces (compiled_traces_hu_*): match ^b_[a-z]+_{resp_clean}_ (non-empty submodel
#     prefix, e.g. "hu_").  Always generated at H3_HALF.  For DVs without a hu submodel
#     (e.g. hppd_binary/Bernoulli), a blank white placeholder PNG is saved instead so
#     the downstream cache check in Python always passes.
#
# Outputs (all in out_dir/):
#   compiled_pp_check_dv.png         — posterior predictive check, DV
#   compiled_pp_check_med.png        — posterior predictive check, mediator
#   compiled_dharma_dv.png           — DHARMa QQ + residuals vs. fitted, DV
#   compiled_dharma_med.png          — DHARMa QQ + residuals vs. fitted, mediator
#   compiled_traces_mu_spvar_dv.png  — MCMC traces, spvar   in DV  mu submodel
#   compiled_traces_mu_med_dv.png    — MCMC traces, mediator in DV  mu submodel
#   compiled_traces_mu_spvar_med.png — MCMC traces, spvar   in med mu submodel
#   compiled_traces_hu_spvar_dv.png  — MCMC traces, spvar   in DV  hu submodel (or blank)
#   compiled_traces_hu_med_dv.png    — MCMC traces, mediator in DV  hu submodel (or blank)
#   compiled_traces_hu_spvar_med.png — MCMC traces, spvar   in med hu submodel (or blank)
#   compiled_resid_dv_sp.png         — DV DHARMa residuals vs. spvar
#   compiled_resid_dv_med.png        — DV DHARMa residuals vs. mediator predictor
#   compiled_resid_med_sp.png        — mediator DHARMa residuals vs. spvar
################################################################################

# ── User library (avoid HPC/cluster package conflicts) ────────────────────────
r_minor    <- sub("\\..*$", "", R.version$minor)
user_r_lib <- path.expand(file.path("~", "R", paste0(R.version$major, ".", r_minor)))
if (!dir.exists(user_r_lib)) dir.create(user_r_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(user_r_lib, .libPaths()))

suppressPackageStartupMessages({
  library(brms)
  library(bayesplot)
  library(ggplot2)
  library(posterior)
  library(DHARMa)
  library(DHARMa.helpers)
})

# ── Parse arguments ────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
fit_path       <- args[1]
out_dir        <- args[2]
spvar          <- args[3]
spvar_label    <- args[4]
mediator       <- args[5]
mediator_label <- args[6]
mediator_in_dv <- args[7]
dv             <- args[8]
dv_label       <- args[9]
W2      <- if (length(args) >= 10 && nchar(args[10]) > 0) as.numeric(args[10]) else 8.5
H2      <- if (length(args) >= 11 && nchar(args[11]) > 0) as.numeric(args[11]) else 5.2
W3      <- if (length(args) >= 12 && nchar(args[12]) > 0) as.numeric(args[12]) else 5.5
H3      <- if (length(args) >= 13 && nchar(args[13]) > 0) as.numeric(args[13]) else 4.1
DPI     <- if (length(args) >= 14 && nchar(args[14]) > 0) as.numeric(args[14]) else 150
H3_HALF <- if (length(args) >= 15 && nchar(args[15]) > 0) as.numeric(args[15]) else H3 / 2

message(paste(rep("=", 70), collapse = ""))
message("compile_diagnostics_helper.R")
message("  fit      : ", fit_path)
message("  out_dir  : ", out_dir)
message("  dv       : ", dv, "  (", dv_label, ")")
message("  spvar    : ", spvar, "  (", spvar_label, ")")
message("  mediator : ", mediator, "  (", mediator_label, ")")
message("  med_in_dv: ", mediator_in_dv)
message("  2-col    : ", W2, "\" x ", H2, "\"  DPI=", DPI)
message("  3-col    : ", W3, "\" x ", H3, "\"  DPI=", DPI)
message(paste(rep("=", 70), collapse = ""))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── brms strips underscores from response names ───────────────────────────────
dv_clean  <- gsub("_", "", dv)
med_clean <- gsub("_", "", mediator)

# ── Load fit ──────────────────────────────────────────────────────────────────
message("\nLoading fit from: ", fit_path)
load(fit_path)   # loads object named 'fit'
message("  Family: ", family(fit)$family, "  N=", nobs(fit))

# ── Canonical white compilation theme (ggplot2) ───────────────────────────────
COMP_THEME <- theme_bw(base_size = 11) +
  theme(
    plot.background   = element_rect(fill = "white", colour = NA),
    panel.background  = element_rect(fill = "white", colour = NA),
    legend.background = element_rect(fill = "white", colour = NA),
    strip.background  = element_rect(fill = "#f2f2f2", colour = NA),
    legend.position   = "none",
    plot.title        = element_blank(),
    plot.margin       = margin(5, 5, 5, 5)
  )

bayesplot_theme_set(theme_bw(base_size = 10) +
  theme(plot.background = element_rect(fill = "white", colour = NA),
        strip.background = element_rect(fill = "#f2f2f2", colour = NA)))

# ── Helper: save ggplot figure ────────────────────────────────────────────────
sg <- function(p, fname, w = W2, h = H2) {
  path <- file.path(out_dir, fname)
  tryCatch(
    ggsave(path, plot = p, width = w, height = h, dpi = DPI, bg = "white"),
    error = function(e) message("  [WARN] ggsave failed for ", fname, ": ", e$message)
  )
  if (file.exists(path)) message("  Saved: ", fname)
}

# ── Helper: save base-R figure ────────────────────────────────────────────────
# plot_fn is an anonymous function (closure) capturing local variables correctly.
spng <- function(fname, w = W2, h = H2, plot_fn) {
  path <- file.path(out_dir, fname)
  tryCatch({
    grDevices::png(path,
                   width  = round(w * DPI),
                   height = round(h * DPI),
                   res    = DPI,
                   bg     = "white")
    par(bg = "white")
    plot_fn()
    dev.off()
    if (file.exists(path)) message("  Saved: ", fname)
  }, error = function(e) {
    message("  [WARN] base-R png failed for ", fname, ": ", e$message)
    if (dev.cur() > 1L) dev.off()
  })
}

################################################################################
# 1. POSTERIOR PREDICTIVE CHECKS
################################################################################
message("\n[1/4] Posterior predictive checks...")

tryCatch({
  p_pp_dv <- pp_check(fit, resp = dv_clean, ndraws = 100, type = "dens_overlay") +
    COMP_THEME +
    labs(x = dv_label, y = "Density")
  sg(p_pp_dv, "compiled_pp_check_dv.png")
}, error = function(e) message("  [WARN] DV pp_check failed: ", e$message))

tryCatch({
  p_pp_med <- pp_check(fit, resp = med_clean, ndraws = 100, type = "dens_overlay") +
    COMP_THEME +
    labs(x = mediator_label, y = "Density")
  sg(p_pp_med, "compiled_pp_check_med.png")
}, error = function(e) message("  [WARN] Med pp_check failed: ", e$message))

################################################################################
# 2. MCMC TRACES — split into mu and hu submodel panels
#
# mu traces (compiled_traces_mu_*):
#   Match ^b_{resp_clean}_ exactly — the mu (count/mean) submodel only.
#   Generated at H3_HALF for dv=="caps_vision" (displayed in top half of split
#   trace row) and at H3 for all other DVs (single full-height trace row).
#
# hu traces (compiled_traces_hu_*):
#   Match ^b_[a-z]+_{resp_clean}_ — any non-empty submodel prefix (e.g. "hu_").
#   Always generated at H3_HALF.  If no hu params exist (e.g. hppd_binary /
#   Bernoulli), a blank white placeholder is saved so Python's cache check passes.
#   The spvar→mediator hu panel (col 3) is almost always blank because mediator
#   submodels (student, Gamma) carry no hu component.
################################################################################
message("\n[2/4] MCMC traces (mu and hu)...")

post_arr <- as.array(fit)
all_pars <- dimnames(post_arr)$variable

# Height for mu traces: H3_HALF for caps_vision (split display), H3 for others.
mu_trace_h   <- if (dv == "caps_vision") H3_HALF else H3
max_pars_mu  <- max(3L, floor(mu_trace_h  / 0.55))   # rows that fit at mu height
max_pars_hu  <- max(2L, floor(H3_HALF     / 0.55))   # rows that fit at H3_HALF

# mu-only: parameter names match ^b_{resp_clean}_ exactly (no submodel prefix).
gen_trace_path_mu <- function(resp_clean, var_col, out_name) {
  mu_pars <- grep(paste0("^b_", resp_clean, "_"), all_pars, value = TRUE, perl = TRUE)
  pars    <- mu_pars[grepl(var_col, mu_pars, fixed = TRUE)]
  if (length(pars) == 0L) {
    message("  [WARN] No mu params: resp=", resp_clean, "  predictor=", var_col,
            "  (skipping ", out_name, ")")
    return(invisible(NULL))
  }
  if (length(pars) > max_pars_mu) {
    message("  Limiting mu to first ", max_pars_mu, " of ", length(pars), " params.")
    pars <- pars[seq_len(max_pars_mu)]
  }
  p <- mcmc_trace(post_arr, pars = pars) + COMP_THEME + ggtitle(NULL)
  sg(p, out_name, w = W3, h = mu_trace_h)
}

# hu-only: parameter names match ^b_[a-z]+_{resp_clean}_ (non-empty submodel
# prefix separated by underscore, e.g. "b_hu_capsvision_*").
# Saves a blank white placeholder when no hu params exist so the Python
# _compiled_paths() cache check always passes for this model.
gen_trace_path_hu <- function(resp_clean, var_col, out_name) {
  hu_pars <- grep(paste0("^b_[a-z]+_", resp_clean, "_"), all_pars, value = TRUE, perl = TRUE)
  pars    <- hu_pars[grepl(var_col, hu_pars, fixed = TRUE)]
  if (length(pars) == 0L) {
    message("  No hu params for resp=", resp_clean, ", var=", var_col,
            " — saving blank placeholder for ", out_name)
    spng(out_name, w = W3, h = H3_HALF, plot_fn = function() {
      par(bg = "white", mar = rep(0, 4))
      plot.new()
    })
    return(invisible(NULL))
  }
  if (length(pars) > max_pars_hu) {
    message("  Limiting hu to first ", max_pars_hu, " of ", length(pars), " params.")
    pars <- pars[seq_len(max_pars_hu)]
  }
  p <- mcmc_trace(post_arr, pars = pars) + COMP_THEME + ggtitle(NULL)
  sg(p, out_name, w = W3, h = H3_HALF)
}

# mu traces
tryCatch(gen_trace_path_mu(dv_clean,  spvar,          "compiled_traces_mu_spvar_dv.png"),
         error = function(e) message("  [WARN] mu spvar→DV: ", e$message))
tryCatch(gen_trace_path_mu(dv_clean,  mediator_in_dv, "compiled_traces_mu_med_dv.png"),
         error = function(e) message("  [WARN] mu med→DV: ", e$message))
tryCatch(gen_trace_path_mu(med_clean, spvar,          "compiled_traces_mu_spvar_med.png"),
         error = function(e) message("  [WARN] mu spvar→med: ", e$message))

# hu traces (blank placeholders when no hu submodel)
tryCatch(gen_trace_path_hu(dv_clean,  spvar,          "compiled_traces_hu_spvar_dv.png"),
         error = function(e) message("  [WARN] hu spvar→DV: ", e$message))
tryCatch(gen_trace_path_hu(dv_clean,  mediator_in_dv, "compiled_traces_hu_med_dv.png"),
         error = function(e) message("  [WARN] hu med→DV: ", e$message))
tryCatch(gen_trace_path_hu(med_clean, spvar,          "compiled_traces_hu_spvar_med.png"),
         error = function(e) message("  [WARN] hu spvar→med: ", e$message))

################################################################################
# 3. DHARMA COMPREHENSIVE (QQ + residuals vs. fitted)
# Uses plotQQunif() and plotResiduals() side by side with proper axis labels.
################################################################################
message("\n[3/4] DHARMa comprehensive...")

# Determine whether each response family is discrete so DHARMa applies the
# correct randomized quantile residual computation. Without integer = TRUE,
# DHARMa treats count/binary responses as continuous, producing incorrect
# QQ plots and p-values. Logic mirrors process_and_save_mediation_model_function.R
# and nonsp_predictors.R.
INTEGER_FAMILIES <- c("hurdle_negbinomial", "bernoulli", "negbinomial",
                      "zero_inflated_negbinomial", "zero_inflated_poisson",
                      "poisson", "ordinal")

dv_fam_str <- tryCatch(
  fit$family[[dv_clean]]$family,
  error = function(e) tryCatch(family(fit)$family, error = function(e2) "unknown")
)
integer_dv <- isTRUE(any(sapply(INTEGER_FAMILIES,
                                function(f) grepl(f, dv_fam_str, ignore.case = TRUE))))

med_fam_str <- tryCatch(
  fit$family[[med_clean]]$family,
  error = function(e) tryCatch(family(fit)$family, error = function(e2) "unknown")
)
integer_med <- isTRUE(any(sapply(INTEGER_FAMILIES,
                                 function(f) grepl(f, med_fam_str, ignore.case = TRUE))))

message("  DV family : ", dv_fam_str, "  →  integer_dv  = ", integer_dv)
message("  Med family: ", med_fam_str, "  →  integer_med = ", integer_med)

sim_res_dv  <- NULL
sim_res_med <- NULL

tryCatch({
  message("  Computing DV residuals (resp = ", dv_clean, ")...")
  sim_res_dv <- dh_check_brms(fit, resp = dv_clean, integer = integer_dv)
  message("  Done.")
}, error = function(e) message("  [WARN] DV DHARMa residuals: ", e$message))

tryCatch({
  message("  Computing mediator residuals (resp = ", med_clean, ")...")
  sim_res_med <- dh_check_brms(fit, resp = med_clean, integer = integer_med)
  message("  Done.")
}, error = function(e) message("  [WARN] Med DHARMa residuals: ", e$message))

if (!is.null(sim_res_dv)) {
  local({
    sr <- sim_res_dv
    lbl <- dv_label
    spng("compiled_dharma_dv.png", w = W2, h = H2, plot_fn = function() {
      par(mfrow = c(1, 2),
          bg   = "white",
          mar  = c(4.5, 4.2, 3.0, 1.0),
          oma  = c(0, 0, 0, 0))
      plotQQunif(sr,
                 testUniformity = TRUE,
                 testOutliers   = TRUE,
                 testDispersion = TRUE)
      plotResiduals(sr,
                    xlab = "Fitted value (rank transformed)",
                    ylab = paste0(lbl, "  [DHARMa residual]"))
    })
  })
}

if (!is.null(sim_res_med)) {
  local({
    sr  <- sim_res_med
    lbl <- mediator_label
    spng("compiled_dharma_med.png", w = W2, h = H2, plot_fn = function() {
      par(mfrow = c(1, 2),
          bg   = "white",
          mar  = c(4.5, 4.2, 3.0, 1.0),
          oma  = c(0, 0, 0, 0))
      plotQQunif(sr,
                 testUniformity = TRUE,
                 testOutliers   = TRUE,
                 testDispersion = TRUE)
      plotResiduals(sr,
                    xlab = "Fitted value (rank transformed)",
                    ylab = paste0(lbl, "  [DHARMa residual]"))
    })
  })
}

################################################################################
# 4. DHARMA RESIDUALS vs. SPECIFIC PREDICTORS
################################################################################
message("\n[4/4] DHARMa residuals vs. predictors...")

# Resolve which mediator column to use when plotting DV residuals vs. mediator.
# mediator_in_dv may be "vch_bl_yes_0_normalized" even though mediator = "vch_bl_yes_0".
med_dv_col <- if (mediator_in_dv %in% colnames(fit$data)) mediator_in_dv else mediator
message("  mediator column in DV formula: ", med_dv_col)

if (!is.null(sim_res_dv)) {

  local({
    sr  <- sim_res_dv
    xd  <- fit$data[[spvar]]
    xl  <- spvar_label
    yl  <- paste0(dv_label, "  [DHARMa residual]")
    spng("compiled_resid_dv_sp.png", w = W3, h = H3, plot_fn = function() {
      par(bg = "white", mar = c(5.5, 5.0, 2.0, 1.0))
      plotResiduals(sr, form = xd, xlab = xl, ylab = yl)
    })
  })

  local({
    sr  <- sim_res_dv
    xd  <- fit$data[[med_dv_col]]
    xl  <- mediator_label
    yl  <- paste0(dv_label, "  [DHARMa residual]")
    spng("compiled_resid_dv_med.png", w = W3, h = H3, plot_fn = function() {
      par(bg = "white", mar = c(5.5, 5.0, 2.0, 1.0))
      plotResiduals(sr, form = xd, xlab = xl, ylab = yl)
    })
  })
}

if (!is.null(sim_res_med)) {
  local({
    sr  <- sim_res_med
    xd  <- fit$data[[spvar]]
    xl  <- spvar_label
    yl  <- paste0(mediator_label, "  [DHARMa residual]")
    spng("compiled_resid_med_sp.png", w = W3, h = H3, plot_fn = function() {
      par(bg = "white", mar = c(5.5, 5.0, 2.0, 1.0))
      plotResiduals(sr, form = xd, xlab = xl, ylab = yl)
    })
  })
}

message("\n", paste(rep("=", 70), collapse = ""))
message("Compilation figure generation complete.")
message(paste(rep("=", 70), collapse = ""))
