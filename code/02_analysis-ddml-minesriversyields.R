#!/usr/bin/env Rscript
# Adversarial & PTD Bias-Correction Report  (Gordon et al. 2026)
#
# Structure:
#   Part 1 – Apples-to-apples: fe-OLS vs PTD vs Adversarial, all using the
#             binary downstream indicator (dist_order ∈ {-1, +1}).
#   Part 2 – Disaggregation: fe-OLS vs PTD-corrected across all dist_order
#             levels, showing how the bias correction shifts estimates.
#
# Requires (working directory = project root):
#   data/ptd_labelled.RDS   -- built by 01_prepare-data-ddml-minesriversyields.R
#   data/ptd_unlabelled.RDS -- built by 01_prepare-data-ddml-minesriversyields.R
#   code/helpers_ddml.R

library("tidyverse")
library("fixest")
library("ggplot2")

source("code/helpers_ddml.R")

# ── Output folders ────────────────────────────────────────────────────────────
p_folder <- "output/plots/"
t_folder <- "output/tables/"
dir.create(p_folder, showWarnings = FALSE, recursive = TRUE)
dir.create(t_folder, showWarnings = FALSE, recursive = TRUE)

# ── Shared specification (mirrors 31_regression_main.R) ──────────────────────
covs    <- c("elevation", "slope", "soilgrid_grouped",
             "tmax_tc", "pre_chirps",
             "accessibility_to_cities_2015", "pop_2015")
cluster <- "mine_basin"
fe      <- c("year", "as.factor(mine_basin)")

# W for adversarial: continuous covariates only.
# soilgrid_grouped (~15 dummies) is excluded because it makes the adversarial
# optimisation underdetermined on ~567 labelled obs.
adv_W <- c("max_EVI_16_cci_c_broad", "elevation", "slope",
           "tmax_tc", "pre_chirps",
           "accessibility_to_cities_2015", "pop_2015")

# ── Load data ─────────────────────────────────────────────────────────────────
lab_raw   <- readRDS("data/ptd_labelled.RDS")
unlab_raw <- readRDS("data/ptd_unlabelled.RDS")

# ── Complete-case filter ──────────────────────────────────────────────────────
req_lab   <- c("log_yield_usd", "y_pred_theta", "downstream", "mine_basin",
               intersect(covs, names(lab_raw)))
req_unlab <- c("y_pred_theta", "downstream", "mine_basin", "year", covs)

lab_complete   <- lab_raw[  complete.cases(lab_raw[,   intersect(req_lab,   names(lab_raw))]),   ]
unlab_complete <- unlab_raw[complete.cases(unlab_raw[, intersect(req_unlab, names(unlab_raw))]), ]

message(sprintf("Labelled rows: %d | Unlabelled rows: %d",
                nrow(lab_complete), nrow(unlab_complete)))

fe_str   <- paste(fe,   collapse = " + ")
cov_str  <- paste(covs, collapse = " + ")
ypred_col <- "y_pred_theta"   # proxy outcome used by all three estimators


# ════════════════════════════════════════════════════════════════════════════
# PART 1 – Apples-to-apples comparison on binary downstream
#   All three estimators use downstream ∈ {0,1} (dist_order ∈ {-1,+1} only)
#   and the same outcome: y_pred_theta (EVI-based predicted log yield).
# ════════════════════════════════════════════════════════════════════════════

message("\n── Part 1: fe-OLS vs PTD vs Adversarial (binary downstream) ──")

# ── (a) fe-OLS: raw y_pred_theta ~ downstream + covs | FE ────────────────────
# This is the uncorrected baseline: regressing the proxy prediction on
# downstream exactly as PTD and adversarial do, but without bias correction.
# beta_raw here should equal ptd_res$beta_raw by construction.
unlab_rdd <- unlab_complete |> dplyr::filter(dist_order %in% c(-1, 1))

ols_fml <- as.formula(paste0(ypred_col, " ~ downstream + ", cov_str, " | ", fe_str))
ols_fit <- fixest::feols(ols_fml, data = unlab_rdd,
                         vcov = as.formula(paste0("~", cluster)))
ols_est <- coef(ols_fit)["downstream"]
ols_se  <- sqrt(diag(vcov(ols_fit)))["downstream"]

message(sprintf("fe-OLS:  beta = %.4f (SE %.4f)", ols_est, ols_se))

# ── (b) PTD-corrected estimate ────────────────────────────────────────────────
message("Running PTD debias...")
ptd_res <- ptd_debias(
  y_true_name     = "log_yield_usd",
  y_pred_name     = "y_pred_theta",
  x_name          = "downstream",
  labelled_data   = lab_complete,
  unlabelled_data = unlab_complete,
  W_names         = covs,
  use_feols       = TRUE,
  cluster         = cluster,
  fe              = fe,
  fe_lab          = "year"   # mine_basin FE collinear with treatment in LSMS data
)
message(sprintf("PTD:     beta_raw = %.4f | gamma = %.4f (SE %.4f) | beta_corrected = %.4f",
                ptd_res$beta_raw, ptd_res$gamma_bias, ptd_res$gamma_se, ptd_res$beta_corrected))

message("Bootstrapping PTD SEs (B=100)...")
boot_ptd <- bootstrap_ptd_se(
  y_true_name     = "log_yield_usd",
  y_pred_name     = "y_pred_theta",
  x_name          = "downstream",
  labelled_data   = lab_complete,
  unlabelled_data = unlab_complete,
  W_names         = covs,
  use_feols       = TRUE,
  cluster         = cluster,
  fe              = fe,
  fe_lab          = "year",
  B               = 100,
  seed            = 42
)
message(sprintf("PTD bootstrap: se = %.4f | CI [%.4f, %.4f]",
                boot_ptd$se_beta_corrected,
                boot_ptd$ci_beta_corrected[1], boot_ptd$ci_beta_corrected[2]))

# ── (c) Adversarial debiasing ─────────────────────────────────────────────────
# The adversarial model (Gordon et al. 2026, eq. 9) trains a prediction model
# with a modified loss that penalises Cov(prediction error, downstream).
# beta_D is then the coefficient from regressing debiased predictions on downstream.
message("Running adversarial debias...")
adv_res <- adversarial_debias(
  y_true_name     = "log_yield_usd",
  y_pred_name     = "y_pred_theta",
  x_name          = "downstream",
  W_names         = adv_W,
  labelled_data   = lab_complete,
  unlabelled_data = unlab_complete,
  alpha           = 1.0,
  use_feols       = TRUE,
  cluster         = cluster,
  fe              = fe
)
message(sprintf("Adversarial: beta_D = %.4f | convergence = %d | cov_nu_x_debiased = %.6f",
                adv_res$beta_D, adv_res$convergence, adv_res$cov_nu_x_debiased))

# Bootstrapping propagates uncertainty from the adversarial optimisation itself
message("Bootstrapping adversarial SEs (B=100)...")
boot_adv <- bootstrap_adversarial_se(
  y_true_name     = "log_yield_usd",
  y_pred_name     = "y_pred_theta",
  x_name          = "downstream",
  W_names         = adv_W,
  labelled_data   = lab_complete,
  unlabelled_data = unlab_complete,
  alpha           = 1.0,
  use_feols       = TRUE,
  cluster         = cluster,
  fe              = fe,
  B               = 100,
  seed            = 42
)
message(sprintf("Adversarial bootstrap: se = %.4f | CI [%.4f, %.4f]",
                boot_adv$se_beta_D, boot_adv$ci_beta_D[1], boot_adv$ci_beta_D[2]))

# ── Part 1 plot ───────────────────────────────────────────────────────────────
comparison_df <- tibble::tibble(
  method   = c("fe-OLS (raw)", "PTD-corrected", "Adversarial"),
  estimate = c(ols_est, ptd_res$beta_corrected, adv_res$beta_D),
  se       = c(ols_se,  boot_ptd$se_beta_corrected, boot_adv$se_beta_D)
) |>
  dplyr::mutate(
    ci_lo  = estimate - 1.96 * se,
    ci_hi  = estimate + 1.96 * se,
    method = factor(method, levels = c("fe-OLS (raw)", "PTD-corrected", "Adversarial"))
  )

p1 <- ggplot(comparison_df,
             aes(x = method, y = estimate, ymin = ci_lo, ymax = ci_hi)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_pointrange(aes(colour = method), linewidth = 0.8, size = 0.9,
                  show.legend = FALSE) +
  scale_colour_manual(values = c(
    "fe-OLS (raw)"   = "#d30303d1",
    "PTD-corrected"  = "#2ca02c",
    "Adversarial"    = "#1f77b4"
  )) +
  labs(
    title    = "Effect of mine downstream exposure on EVI-predicted yield: method comparison",
    subtitle = sprintf("PTD bias \u03b3 = %.3f (SE %.3f) | PTD corrected = %.3f | Adversarial = %.3f",
                       ptd_res$gamma_bias, ptd_res$gamma_se,
                       ptd_res$beta_corrected, adv_res$beta_D),
    x = NULL,
    y = "Downstream coefficient (log EVI)"
  ) +
  theme_bw(base_size = 12)

ggsave(paste0(p_folder, "01_method_comparison.pdf"), p1, width = 7, height = 5)
message("Plot saved: output/plots/01_method_comparison.pdf")


# ════════════════════════════════════════════════════════════════════════════
# PART 2 – Disaggregation: fe-OLS vs PTD-corrected across dist_order levels
# ════════════════════════════════════════════════════════════════════════════

message("\n── Part 2: fe-OLS vs PTD vs Adversarial disaggregated by dist_order ──")

# All three estimators use the same joint feols with i(dist_order, ref = -1),
# mirroring 31_regression_main.R exactly.  dist_order = 0 (mine's own basin)
# is kept in the data but dropped from the plot.
dist_orders_plot <- sort(unique(unlab_complete$dist_order))
dist_orders_plot <- dist_orders_plot[dist_orders_plot != -1]  # ref level; coef = 0 by def.

# ── (a) fe-OLS raw: y_pred_theta ~ i(dist_order, ref=-1) + covs | FE ─────────
raw_fit <- fixest::feols(
  as.formula(paste0(ypred_col, " ~ i(dist_order, ref = -1) + ", cov_str, " | ", fe_str)),
  data = unlab_complete,
  vcov = as.formula(paste0("~", cluster))
)

raw_ests <- tibble::tibble(
  dist_order = dist_orders_plot,
  beta_raw   = coef(raw_fit)[paste0("dist_order::", dist_orders_plot)],
  se_raw     = sqrt(diag(vcov(raw_fit)))[paste0("dist_order::", dist_orders_plot)]
)

# ── (b) PTD-corrected: shift every dist_order coef by the global gamma ────────
# Gamma is a single scalar estimated from labelled data independently of the
# unlabelled feols, so Var(beta_raw - gamma) = Var(beta_raw) + Var(gamma).
gamma_global <- ptd_res$gamma_bias
se_gamma     <- boot_ptd$se_gamma

ptd_ests <- raw_ests |>
  dplyr::mutate(
    beta_ptd = beta_raw - gamma_global,
    se_ptd   = sqrt(se_raw^2 + se_gamma^2)
  )

# ── (c) Adversarial: re-run final regression on debiased predictions ──────────
# adversarial_debias() trains omega on labelled data to zero cov(nu, downstream)
# and stores debiased predictions (y_pred_D) for every unlabelled row.
# We take those predictions — already computed in Part 1 (adv_res) — and run
# the same i(dist_order, ref=-1) regression on them to get per-order coefficients.
# This is valid because the debiasing happens at the prediction level; the final
# regression is just reading off the corrected signal at each stream order.
unlab_with_pred <- unlab_complete
unlab_with_pred$y_pred_D <- adv_res$predictions_unlab

adv_fit <- fixest::feols(
  as.formula(paste0("y_pred_D ~ i(dist_order, ref = -1) + ", cov_str, " | ", fe_str)),
  data = unlab_with_pred,
  vcov = as.formula(paste0("~", cluster))
)

# Bootstrap SEs for adversarial dist_order coefficients by resampling the
# full adversarial optimisation B times and re-running the joint regression.
message("Bootstrapping adversarial dist_order SEs (B=100)...")
set.seed(42)
B <- 100
boot_adv_order <- purrr::map(seq_len(B), function(b) {
  lab_b   <- lab_complete[  sample(nrow(lab_complete),   replace = TRUE), ]
  unlab_b <- unlab_complete[sample(nrow(unlab_complete), replace = TRUE), ]
  tryCatch({
    res_b <- adversarial_debias(
      y_true_name     = "log_yield_usd",
      y_pred_name     = "y_pred_theta",
      x_name          = "downstream",
      W_names         = adv_W,
      labelled_data   = lab_b,
      unlabelled_data = unlab_b,
      alpha           = 1.0,
      use_feols       = TRUE,
      cluster         = cluster,
      fe              = fe
    )
    unlab_b$y_pred_D <- res_b$predictions_unlab
    fit_b <- fixest::feols(
      as.formula(paste0("y_pred_D ~ i(dist_order, ref = -1) + ", cov_str, " | ", fe_str)),
      data = unlab_b,
      vcov = as.formula(paste0("~", cluster))
    )
    coef(fit_b)[paste0("dist_order::", dist_orders_plot)]
  }, error = function(e) rep(NA_real_, length(dist_orders_plot)))
})

boot_adv_order <- do.call(rbind, boot_adv_order)  # B x length(dist_orders_plot) matrix

adv_ests <- tibble::tibble(
  dist_order = dist_orders_plot,
  beta_adv   = coef(adv_fit)[paste0("dist_order::", dist_orders_plot)],
  se_adv     = apply(boot_adv_order, 2, sd, na.rm = TRUE)
)

# ── Part 2 plot ───────────────────────────────────────────────────────────────
plot2_df <- raw_ests |>
  dplyr::left_join(adv_ests, by = "dist_order") |>
  dplyr::filter(dist_order != 0) |>   # drop mine-basin level; not meaningful
  tidyr::pivot_longer(
    cols      = c(beta_raw, beta_adv),
    names_to  = "method",
    values_to = "estimate"
  ) |>
  dplyr::mutate(
    se    = dplyr::case_when(
      method == "beta_raw" ~ se_raw,
      method == "beta_adv" ~ se_adv
    ),
    ci_lo  = estimate - 1.96 * se,
    ci_hi  = estimate + 1.96 * se,
    method = dplyr::recode(method,
      "beta_raw" = "fe-OLS (raw)",
      "beta_adv" = "Adversarial ML correction"
    ),
    method = factor(method, levels = c("fe-OLS (raw)", "Adversarial ML correction"))
  )

dist_orders_axis <- sort(unique(plot2_df$dist_order))

p2 <- ggplot(plot2_df,
             aes(x = dist_order, y = estimate,
                 ymin = ci_lo, ymax = ci_hi, colour = method)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_pointrange(position = position_dodge(width = 0.4),
                  linewidth = 0.7, size = 0.6) +
  scale_colour_manual(values = c(
    "fe-OLS (raw)" = "#d30303d1",
    "Adversarial ML correction"  = "#1f77b4"
  )) +
  scale_x_continuous(breaks = dist_orders_axis) +
  labs(
    title    = "Effect of downstream mine exposure on EVI-predicted yield by stream order",
    subtitle = "Adversarial debiasing (Gordon et al. 2026); ref = dist_order -1 (upstream)",
    x        = "Stream order relative to mine (negative = upstream, positive = downstream)",
    y        = "Coefficient vs. upstream (dist_order = -1)",
    colour   = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")
p2
ggsave(paste0(p_folder, "02_dist_order_comparison.pdf"), p2, width = 9, height = 5)
message("Plot saved: output/plots/02_dist_order_comparison.pdf")


# ── Summary table ─────────────────────────────────────────────────────────────
summary_tbl <- ptd_ests |>
  dplyr::left_join(adv_ests, by = "dist_order") |>
  dplyr::mutate(gamma_global = gamma_global, se_gamma = se_gamma)
print(summary_tbl)

invisible(list(
  ols_est    = ols_est,
  ptd_res    = ptd_res,
  boot_ptd   = boot_ptd,
  adv_res    = adv_res,
  boot_adv   = boot_adv,
  summary    = summary_tbl
))
