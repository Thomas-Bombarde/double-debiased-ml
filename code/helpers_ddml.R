# ── Approach 1: PTD / Bias-correction (simpler, what your code approximates) ──
# Estimate bias γ from labelled data, subtract from full-sample β
#
# Arguments
#   y_true_name   column name of true outcome in labelled_data
#   y_pred_name   column name of proxy prediction in both datasets
#   x_name        column name of treatment (e.g. "downstream")
#   labelled_data data.frame with true outcome available
#   unlabelled_data data.frame for the full-sample β regression
#   W_names       optional character vector of covariate names to condition on
#   use_feols     logical; if TRUE use fixest::feols with fe/cluster
#   cluster       column name (string) for clustered SEs (feols only)
#   fe            character vector of fixed-effect terms (feols only)
get_ols_estimate <- function(data) {
  mod <- lm(sattelite_cov ~ roads, data = data)
  coefs <- summary(mod)$coefficients
  beta_ols <- coefs["roads", "Estimate"]
  std_ols <- coefs["roads", "Std. Error"]
  return(data.frame("beta_ols" = beta_ols, "std_ols" = std_ols))
}

ptd_debias <- function(y_true_name, y_pred_name, x_name,
                       labelled_data, unlabelled_data,
                       W_names   = NULL,
                       use_feols = FALSE,
                       cluster   = NULL,
                       fe        = NULL,
                       fe_lab    = NULL) {
  # fe_lab: FE terms for the gamma (labelled-data) regression.
  # Defaults to fe, but callers can pass a simpler spec (e.g. just "year")
  # when the full fe would be collinear in the smaller labelled dataset.
  if (is.null(fe_lab)) fe_lab <- fe
  lab   <- labelled_data
  unlab <- unlabelled_data

  # Prediction errors on labelled points
  lab$nu <- lab[[y_pred_name]] - lab[[y_true_name]]

  # Helper: build and fit a regression (lm or feols) returning a fitted model
  # `covariates` is the full RHS vector: c(x_name, W_names) or just c(x_name)
  # `fe` terms referencing columns absent from `data` are silently dropped.
  .fit_reg <- function(response, covariates, data, use_feols, cluster, fe) {
    if (use_feols) {
      xterms <- paste(covariates, collapse = " + ")
      # Drop any FE term whose base column does not exist in data
      # (e.g. "as.factor(mine_basin)" when mine_basin is absent, or
      #  terms that would be perfectly collinear because every unit maps
      #  1-to-1 with treatment — detected by checking column presence)
      fe_valid <- fe[sapply(fe, function(term) {
        # extract bare column name from expressions like "as.factor(col)"
        col <- gsub(".*\\((.+)\\).*", "\\1", term)
        col %in% names(data)
      })]
      fe_str <- if (length(fe_valid) > 0) paste(fe_valid, collapse = " + ") else "0"
      fml    <- as.formula(paste(response, "~", xterms, "|", fe_str))
      fit    <- fixest::feols(fml, data = data,
                              vcov = if (!is.null(cluster) && cluster %in% names(data))
                                       as.formula(paste0("~", cluster)) else "iid")
    } else {
      rhs_terms <- paste(covariates, collapse = " + ")
      fml        <- as.formula(paste(response, "~", rhs_terms))
      fit        <- lm(fml, data = data)
    }
    fit
  }

  all_x_lab   <- c(x_name, intersect(W_names, names(lab)))    # only covs present in labelled data
  all_x_unlab <- c(x_name, intersect(W_names, names(unlab)))  # only covs present in unlabelled data

  # Estimate bias γ: regress errors on treatment (+ covs) in labelled data
  bias_mod  <- .fit_reg("nu",        all_x_lab,   lab,   use_feols, cluster, fe_lab)
  gamma_hat <- coef(bias_mod)[x_name]

  # Estimate β: regress proxy on treatment (+ covs) in unlabelled data
  full_mod  <- .fit_reg(y_pred_name, all_x_unlab, unlab, use_feols, cluster, fe)
  beta_hat  <- coef(full_mod)[x_name]

  # SE of gamma (from model)
  gamma_se  <- tryCatch(
    sqrt(diag(vcov(bias_mod)))[x_name],
    error = function(e) NA_real_
  )

  # Bias-corrected estimate (equation 6 in paper)
  beta_corrected <- beta_hat - gamma_hat

  list(
    beta_raw       = beta_hat,
    gamma_bias     = gamma_hat,
    gamma_se       = gamma_se,
    beta_corrected = beta_corrected,
    bias_model     = bias_mod,
    full_model     = full_mod
  )
}


# ── Bootstrap standard errors for PTD-corrected estimate ──────────────────────
# Same arguments as ptd_debias, plus B (number of replicates) and seed.
# Returns a list with:
#   boot_beta_corrected  numeric vector of length B (bias-corrected beta per replicate)
#   boot_gamma           numeric vector of length B (gamma per replicate)
#   boot_beta_raw        numeric vector of length B (raw beta on unlabelled per replicate)
#   se_beta_corrected    bootstrap SE of beta_corrected
#   se_gamma             bootstrap SE of gamma
#   ci_beta_corrected    named numeric[2]: 2.5% and 97.5% quantile CI

bootstrap_ptd_se <- function(y_true_name, y_pred_name, x_name,
                              labelled_data, unlabelled_data,
                              W_names   = NULL,
                              use_feols = FALSE,
                              cluster   = NULL,
                              fe        = NULL,
                              fe_lab    = NULL,
                              B         = 500,
                              seed      = NULL) {

  if (!is.null(seed)) set.seed(seed)

  n_lab   <- nrow(labelled_data)
  n_unlab <- nrow(unlabelled_data)

  boot_gamma  <- numeric(B)
  boot_beta   <- numeric(B)
  boot_bc     <- numeric(B)

  for (b in seq_len(B)) {
    # Resample both datasets independently
    lab_b   <- labelled_data[  sample(n_lab,   replace = TRUE), , drop = FALSE]
    unlab_b <- unlabelled_data[sample(n_unlab, replace = TRUE), , drop = FALSE]

    res <- tryCatch(
      ptd_debias(
        y_true_name     = y_true_name,
        y_pred_name     = y_pred_name,
        x_name          = x_name,
        labelled_data   = lab_b,
        unlabelled_data = unlab_b,
        W_names         = W_names,
        use_feols       = use_feols,
        cluster         = cluster,
        fe              = fe,
        fe_lab          = fe_lab
      ),
      error = function(e) NULL
    )

    boot_gamma[b] <- if (!is.null(res)) res$gamma_bias   else NA_real_
    boot_beta[b]  <- if (!is.null(res)) res$beta_raw     else NA_real_
    boot_bc[b]    <- if (!is.null(res)) res$beta_corrected else NA_real_
  }

  # Drop failed replicates
  ok <- is.finite(boot_gamma) & is.finite(boot_beta)
  boot_gamma <- boot_gamma[ok]
  boot_beta  <- boot_beta[ok]
  boot_bc    <- boot_bc[ok]

  list(
    boot_beta_corrected = boot_bc,
    boot_gamma          = boot_gamma,
    boot_beta_raw       = boot_beta,
    se_beta_corrected   = sd(boot_bc),
    se_gamma            = sd(boot_gamma),
    ci_beta_corrected   = quantile(boot_bc, probs = c(0.025, 0.975))
  )
}


# ── Bootstrap standard errors for adversarial debiasing ───────────────────────
# Resamples both labelled and unlabelled datasets on each replicate so that
# uncertainty in both the adversarial optimisation AND the final regression
# is propagated into the SE of beta_D.
#
# Same arguments as adversarial_debias, plus B and seed.
# Returns a list with:
#   boot_beta_D   numeric vector of length B
#   se_beta_D     bootstrap SE
#   ci_beta_D     named numeric[2]: 2.5% and 97.5% quantile CI

bootstrap_adversarial_se <- function(y_true_name, y_pred_name, x_name, W_names,
                                     labelled_data, unlabelled_data,
                                     alpha     = 1.0,
                                     use_feols = FALSE,
                                     cluster   = NULL,
                                     fe        = NULL,
                                     B         = 500,
                                     seed      = NULL) {

  if (!is.null(seed)) set.seed(seed)

  n_lab   <- nrow(labelled_data)
  n_unlab <- nrow(unlabelled_data)

  boot_beta_D <- numeric(B)

  for (b in seq_len(B)) {
    lab_b   <- labelled_data[  sample(n_lab,   replace = TRUE), , drop = FALSE]
    unlab_b <- unlabelled_data[sample(n_unlab, replace = TRUE), , drop = FALSE]

    res <- tryCatch(
      adversarial_debias(
        y_true_name     = y_true_name,
        y_pred_name     = y_pred_name,
        x_name          = x_name,
        W_names         = W_names,
        labelled_data   = lab_b,
        unlabelled_data = unlab_b,
        alpha           = alpha,
        use_feols       = use_feols,
        cluster         = cluster,
        fe              = fe
      ),
      error = function(e) NULL
    )

    boot_beta_D[b] <- if (!is.null(res)) res$beta_D else NA_real_
  }

  ok <- is.finite(boot_beta_D)
  boot_beta_D <- boot_beta_D[ok]

  list(
    boot_beta_D = boot_beta_D,
    se_beta_D   = sd(boot_beta_D),
    ci_beta_D   = quantile(boot_beta_D, probs = c(0.025, 0.975))
  )
}


# ── Approach 2: Adversarial debiasing (Gordon et al 2026) ─────────────
adversarial_loss <- function(par, y_true, W, x, alpha) {
  # par: coefficients of the PRIMARY prediction model (predicting y from W)
  b0 <- par[1]
  bW <- par[2:length(par)] 
  
  y_fitted <- b0 + as.matrix(W) %*% bW   # primary model predictions
  nu       <- y_true - y_fitted           # prediction errors
  
  mse       <- mean(nu^2)
  # Penalise covariance of errors with treatment x (equation 9 in paper)
  cov_term  <- mean((x - mean(x)) * nu)
  
  mse + alpha * cov_term^2
}

adversarial_debias <- function(y_true_name, y_pred_name, x_name, W_names,
                               labelled_data, unlabelled_data,
                               alpha     = 1.0,
                               use_feols = FALSE,
                               cluster   = NULL,
                               fe        = NULL) {
  lab   <- labelled_data
  unlab <- unlabelled_data

  # Validate W_names exist in both datasets and restrict to the intersection
  if (is.null(W_names)) W_names <- character(0)
  W_lab_names   <- intersect(W_names, names(lab))
  W_unlab_names <- intersect(W_names, names(unlab))

  if (length(W_lab_names) == 0) {
    stop(sprintf("adversarial_debias: none of the requested W_names are present in labelled_data. Available cols: %s",
                 paste(names(lab), collapse = ", ")))
  }
  if (length(W_unlab_names) == 0) {
    stop(sprintf("adversarial_debias: none of the requested W_names are present in unlabelled_data. Available cols: %s",
                 paste(names(unlab), collapse = ", ")))
  }

  # Use only the valid columns (keep ordering from W_names)
  W_lab_names   <- intersect(W_names, W_lab_names)
  W_unlab_names <- intersect(W_names, W_unlab_names)

  y_true <- lab[[y_true_name]]
  x_lab  <- lab[[x_name]]
  W_lab  <- lab[, W_lab_names, drop = FALSE]
  W_unlab <- unlab[, W_unlab_names, drop = FALSE]

  # ── Dummy-encode factor/character columns (model.matrix drops the intercept)
  # This converts e.g. soilgrid_grouped into numeric indicator columns so that
  # as.matrix(W) %*% bW works correctly in adversarial_loss.
  # Both matrices must share exactly the same columns: encode on the union of
  # factor levels, then align unlab to lab's column set (add 0-cols if missing,
  # drop any extra cols that appear in unlab but not in lab).
  .to_numeric_matrix <- function(df) {
    if (ncol(df) == 0) return(matrix(0, nrow = nrow(df), ncol = 0))
    model.matrix(~ . - 1, data = df)
  }
  W_lab_mat_raw   <- .to_numeric_matrix(W_lab)
  W_unlab_mat_raw <- .to_numeric_matrix(W_unlab)

  # Align columns: unlab gets exactly the same columns as lab (in same order)
  lab_cols <- colnames(W_lab_mat_raw)
  W_lab_mat <- W_lab_mat_raw   # already correct
  W_unlab_mat <- {
    m <- matrix(0, nrow = nrow(W_unlab_mat_raw), ncol = length(lab_cols),
                dimnames = list(NULL, lab_cols))
    shared <- intersect(lab_cols, colnames(W_unlab_mat_raw))
    m[, shared] <- W_unlab_mat_raw[, shared, drop = FALSE]
    m
  }

  # ── Starting values from plain OLS of y_true on W (encoded columns)
  # Use .lm.fit() directly on the matrix to avoid formula string-length limits
  # when soilgrid_grouped expands into many dummy columns.
  X_start   <- cbind(1, W_lab_mat)   # add intercept column
  ols_start <- as.numeric(.lm.fit(X_start, y_true)$coefficients)

  # ── Adversarial optimisation over PRIMARY model weights ────────────────────
  fit <- optim(
    par    = ols_start,
    fn     = adversarial_loss,
    y_true = y_true,
    W      = W_lab_mat,
    x      = x_lab,
    alpha  = alpha,
    method = "BFGS",
    control = list(maxit = 100)
  )

  if (fit$convergence != 0) warning("optim did not converge")

  # ── Generate debiased predictions for ALL data ─────────────────────────────
  b_hat <- fit$par

  # predictions on labelled data (for diagnostics)
  y_pred_lab_D <- b_hat[1] + W_lab_mat %*% b_hat[-1]

  # predictions on unlabelled data (used for final regression)
  y_pred_unlab_D <- b_hat[1] + W_unlab_mat %*% b_hat[-1]
  unlab$y_pred_D <- as.numeric(y_pred_unlab_D)

  # ── Estimate causal effect β using debiased predictions ────────────────────
  if (use_feols) {
    fe_str     <- if (!is.null(fe) && length(fe) > 0) paste(fe, collapse = " + ") else "0"
    causal_fml <- as.formula(paste("y_pred_D ~", x_name, "|", fe_str))
    causal_mod <- fixest::feols(causal_fml, data = unlab,
                                vcov = if (!is.null(cluster)) as.formula(paste0("~", cluster)) else "iid")
  } else {
    causal_mod <- lm(reformulate(x_name, response = "y_pred_D"), data = unlab)
  }

  # ── Check residual covariance with x on labelled data (diagnostics) ────────
  nu_lab_D   <- as.numeric(y_true - y_pred_lab_D)
  nu_lab_raw <- y_true - lab[[y_pred_name]]

  list(
    beta_D            = coef(causal_mod)[x_name],
    causal_model      = causal_mod,
    adversarial_par   = b_hat,
    convergence       = fit$convergence,
    # Diagnostics
    cov_nu_x_raw      = cov(nu_lab_raw, x_lab),   # should be nonzero
    cov_nu_x_debiased = cov(nu_lab_D,   x_lab),   # should be ~0
    predictions_unlab = unlab$y_pred_D
  )
}
