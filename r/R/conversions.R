# =============================================================================
# r/R/conversions.R — scale / model / direction conversions of EWAS sumstats
#             and phenotypic variance inference
# =============================================================================

#' Convert M-value-scale EWAS summary statistics to the beta-value scale
#'
#' Uses the first-order approximation dM/dbeta = 1 / (log(2) * beta * (1 -
#' beta)) evaluated at the mean methylation of each CpG. The approximation
#' degrades at extreme methylation, so very high or very low CpGs should be
#' filtered first.
#'
#' Best if the mean beta values are computed from the same samples as the EWAS
#' itself, since the conversion is a local linear approximation evaluated at
#' these values.
#'
#' @param ss EWAS summary statistics on m-value scale, with `CpGmarker`, `BETA`, `SE`.
#' @param mean_beta_df A `data.table` with columns `CpGmarker` and
#'   `beta_bar`.
#' @param ewas_direction Either `"M~Y"` or `"Y~M"`, already normalised.
#' @param verbose Logical.
#'
#' @return The summary statistics on the beta-value scale, restricted to
#'   CpGs present in `mean_beta_df`.
#' @noRd
convert_mvalue_ewas_to_beta_scale <- function(ss, mean_beta_df, ewas_direction, verbose = TRUE) {
  log_message(verbose, "INFO", "Converting M-value EWAS summary statistics to Beta-value scale.")
  
  # Check for extreme mean population beta values
  num_cpgs_low <- sum(mean_beta_df$beta_bar < 0.1, na.rm = TRUE)
  num_cpgs_high <- sum(mean_beta_df$beta_bar > 0.9, na.rm = TRUE)
  num_cpgs_total <- nrow(mean_beta_df)
  
  if (num_cpgs_low > 0 || num_cpgs_high > 0) {
    log_message(verbose, "WARN", "Mean beta reference contains %d CpGs with beta < 0.1 and %d CpGs with beta > 0.9 (out of %d total). Conversion from M-value to beta-value scale may be inaccurate for these extreme CpGs.", num_cpgs_low, num_cpgs_high, num_cpgs_total)
  }

  match_idx <- match(ss$CpGmarker, mean_beta_df$CpGmarker)
  missing_from_mean_beta_idx <- is.na(match_idx)

  if (any(missing_from_mean_beta_idx)) {
    log_message(verbose, "INFO", "Dropping %d CpGs from the EWAS summary stats because they are missing from the mean_beta reference.", sum(missing_from_mean_beta_idx))
    ss <- ss[!missing_from_mean_beta_idx, ]
    match_idx <- match_idx[!missing_from_mean_beta_idx]
  }

  beta_bar <- mean_beta_df$beta_bar[match_idx]
  cc <- log(2) * beta_bar * (1 - beta_bar)

  if (ewas_direction == "M~Y") {
    ss$BETA <- cc * ss$BETA
    ss$SE <- cc * ss$SE
  } else { # Y ~ M:
    ss$BETA <- ss$BETA / cc
    ss$SE <- ss$SE / cc
  }

  return(ss)
}

#' Rescale Cox summary statistics onto the linear Y ~ M scale
#'
#' Multiplies log hazard ratios by kappa(pi) = (1 - pi) * (-log(1 - pi)),
#' where pi is the proportion of cases at the end of follow-up. This is a
#' first-order approximation using the delta method.
#'     BETA' = kappa(pi) * BETA
#'     SE(BETA') = kappa(pi) * SE(BETA)
#'
#' @param ss EWAS summary statistics with `BETA`, `SE`, `N`.
#' @param cox_pi Proportion of cases at the end of the study, or `NULL`.
#' @param cox_E Number of events during follow-up, or `NULL`.
#' @param verbose Logical.
#'
#' @return A list with the rescaled `ss`, the `cox_pi` used, and `kappa_cox`.
#' @noRd
convert_cox_ewas_to_linear_scale <- function(ss, cox_pi, cox_E, verbose = TRUE) {
  log_message(verbose, "INFO", "Converting Cox regression summary statistics to linear Y ~ M model...")

  if (is.null(cox_pi) && is.null(cox_E)) {
    stop("[ERROR]: EWAS is Cox PH model, must provide either `cox_pi` or 'cox_E' and 'N'.", call. = FALSE)
  }

  if (is.null(cox_pi)) {
    cox_N <- stats::median(ss$N, na.rm = TRUE)
    cox_pi <- cox_E / cox_N
  }

  # Check for reasonable cox_pi values.
  if (!is.finite(cox_pi) || cox_pi <= 0 || cox_pi >= 1) {
    stop(sprintf("[ERROR] cox_pi must be strictly between 0 and 1; got %.6g.", cox_pi), call. = FALSE)
  }

  kappa_cox <- (1 - cox_pi) * (-log(1 - cox_pi))

  log_message(verbose, "INFO", "Cox conversion factor kappa(pi) = %.6g (pi = %.6g)", kappa_cox, cox_pi)

  ss$BETA <- ss$BETA * kappa_cox
  ss$SE <- ss$SE * kappa_cox

  return(list(ss = ss, cox_pi = cox_pi, kappa_cox = kappa_cox))
}

#' Convert linear Y ~ M summary statistics to linear M ~ Y
#'
#' Invert the regression per CpG using summary statistics:
#'    BETA' = BETA / (BETA^2 + (N - 2) * SE(BETA)^2)
#'    SE(BETA') = SE(BETA) / (BETA^2 + (N - 2) * SE(BETA)^2)
#' Derived from standard OLS formulas:
#'    BETA = Cov(Y,M) / Var(M)
#'    SE(BETA)² ≈ sigma_e² / ((N-2)·Var(M))
#'
#' @param ss EWAS summary statistics on beta-value, linear model scale, with `BETA`, `SE`, `N`.
#' @param verbose Logical.
#'
#' @return The summary statistics in the M ~ Y direction on the linear scale.
#' @noRd
convert_linear_y_on_m_ewas_to_m_on_y <- function(ss, verbose = TRUE) {
  log_message(verbose, "INFO", "Converting linear Y ~ M summary statistics to M ~ Y direction...")

  denom <- ss$BETA^2 + (ss$N - 2) * ss$SE^2

  ss$BETA <- ss$BETA / denom
  ss$SE <- ss$SE / denom

  return(ss)
}

#' Convert logistic Y ~ M summary statistics to linear M ~ Y
#'
#' Conversion uses the reference methylation variance:
#'    Let X = sqrt(1 + 4 * SIGMA * BETA² * Var(Y))
#'    BETA' = 2 * SIGMA * BETA / (1 + X)
#'    SE(BETA') = SE(BETA) * (2 * SIGMA / (1 + X)) / X
#'
#' Assumes small effects of Y on M and equal methylation variance in cases
#' and controls.
#'
#' @param ss EWAS summary statistics in the Y~M direction for a logistic regression
#'            using beta-value scale, with `CpGmarker`, `BETA`, `SE`.
#' @param Sigma Reference DNAm covariance matrix.
#' @param varY Covariate-adjusted phenotype variance
#' @param verbose Logical.
#'
#' @return The summary statistics in the M ~ Y direction on the linear scale.
#' @noRd
convert_logistic_y_on_m_ewas_to_linear_m_on_y <- function(ss, Sigma, varY, verbose = TRUE) {
  log_message(verbose, "INFO", "Converting logistic Y ~ M summary statistics to linear M ~ Y...")
  log_message(verbose, "NOTE", "Logistic to linear conversion assumes small effect sizes of Y on M and equal methylation variance between cases and controls.")

  # Convert logistic regression coefficients to linear regression coefficients
  Sigma_aligned = diag(Sigma)[ss$CpGmarker]
  aa <- sqrt(1 + 4 * ss$BETA^2 * Sigma_aligned * varY)
  a2b <- 2 * Sigma_aligned / (1 + aa)
  ss$BETA <- ss$BETA * a2b
  ss$SE <- ss$SE * (a2b / aa) # delta-method

  return(ss)
}

#' Infer Var(Y) from logistic Y ~ M summary statistics and a reference covariance
#' matrix.
#'
#' For each overlapping CpG, Var(Y) = 1/(N * Var(M) * SE(BETA)^2). The
#' median across CpGs is taken for stability.
#'
#' Note: If there are not enough CpGs overlapping both Sigma and the EWAS sumstats
#' (e.g. if the CpGs were already subset to the CpGs in a single clock), then the
#' Var(Y) estimate may be unstable/inaccurate/negative.
#'
#' @param ss EWAS summary statistics in the Y~M direction for a logistic regression
#'            using beta-value scale, with `CpGmarker`, `SE`, `N`.
#' @param Sigma Reference DNAm covariance matrix.
#' @param verbose Logical.
#'
#' @return A single numeric estimate of Var(Y).
#' @noRd
infer_phenotype_variance_from_logistic_y_on_m_ewas <- function(ss, Sigma, verbose = TRUE) {
  # Only use the cpgs in both the EWAS sumstats and the reference Sigma matrix
  cpg_subset <- intersect(ss$CpGmarker, rownames(Sigma))
  if (length(cpg_subset) == 0) {
    stop("[ERROR] Var(Y) cannot be inferred: no overlapping CpGs between the EWAS sumstats and DNAm_cov_matrix.", call. = FALSE)
  }
  ss_subset <- ss[match(cpg_subset, ss$CpGmarker), ]
  Sigma_subset <- Sigma[cpg_subset, cpg_subset]
  varY_hat <- stats::median(1 / (ss_subset$N * diag(Sigma_subset) * (ss_subset$SE^2)), na.rm = TRUE)

  log_message(verbose, "INFO", "Inferred covariate-adjusted Var(Y) from logistic EWAS sumstats and DNAm_cov_matrix: %.6f", varY_hat)
  return(varY_hat)
}

#' Infer Var(Y) from linear M ~ Y summary statistics and a reference covariance
#' matrix.
#'
#' For each overlapping CpG, Var(Y) = Var(M) / (BETA^2 + (N - 2) * SE(BETA)^2). The
#' median across CpGs is taken for stability.
#'
#' Note: If there are not enough CpGs overlapping both Sigma and the EWAS sumstats
#' (e.g. if the CpGs were already subset to the CpGs in a single clock), then the
#' Var(Y) estimate may be unstable/inaccurate/negative.
#'
#' @param ss EWAS summary statistics in the M ~ Y direction for a linear regression
#'            using beta-value scale.
#' @param Sigma Reference DNAm covariance matrix.
#' @param verbose Logical.
#'
#' @return A single numeric estimate of Var(Y).
#' @noRd
infer_phenotype_variance_from_linear_m_on_y_ewas <- function(ss, Sigma, verbose = TRUE) {
  if (!all(c("SE", "N") %in% names(ss))) {
    stop("[ERROR] Var(Y) cannot be inferred from EWAS sumstats because SE and/or N columns are missing.", call. = FALSE)
  }

  # Only use the cpgs in both the EWAS sumstats and the reference Sigma matrix
  cpg_subset <- intersect(ss$CpGmarker, rownames(Sigma))
  if (length(cpg_subset) == 0) {
    stop("[ERROR] Var(Y) cannot be inferred: no overlapping CpGs between the EWAS sumstats and DNAm_cov_matrix.", call. = FALSE)
  }
  ss_subset <- ss[match(cpg_subset, ss$CpGmarker), ]
  Sigma_subset <- Sigma[cpg_subset, cpg_subset]

  varY_hats <- diag(Sigma_subset) / (ss_subset$BETA^2 + (ss_subset$N - 2) * ss_subset$SE^2)
  varY_hat <- stats::median(varY_hats, na.rm = TRUE)

  log_message(verbose, "INFO", "Inferred covariate-adjusted Var(Y) from linear EWAS sumstats and DNAm_cov_matrix: %.6f", varY_hat)
  return(varY_hat)
}
