# =============================================================================
# R/src/estimate.R — CpG alignment and the core clock/phenotype
#             regression estimators
# =============================================================================

#' Drop intercept rows from a clock weight table
#'
#' Clock files frequently carry an "(Intercept)" row that is not a CpG and
#' must not be counted towards clock coverage.
#'
#' @param clock_weights Clock weights with a `CpGmarker` column.
#' @param verbose Logical.
#'
#' @return The clock weights with intercept rows removed.
#' @noRd
drop_intercept_rows <- function(clock_weights, verbose = TRUE) {
  n_intercept_rows <- sum(grepl("^\\(?intercept\\)?$", clock_weights$CpGmarker, ignore.case = TRUE))
  if (n_intercept_rows > 0) {
    log_message(verbose, "INFO", "Removed %d intercept row(s) from the clock weights.", n_intercept_rows)
    clock_weights <- clock_weights[!grepl("^\\(?intercept\\)?$", clock_weights$CpGmarker, ignore.case = TRUE), ]
  }
  clock_weights
}

#' Summarize CpG overlap across the three inputs
#'
#' @param ss EWAS summary statistics.
#' @param cw Clock weights, with the intercept row already removed.
#' @param Sigma Reference DNAm covariance matrix.
#' @param verbose Logical.
#'
#' @return A list of counts plus `cpg_intersect`, the CpGs common to all
#'   three.
#' @noRd
summarize_cpg_overlap <- function(ss, cw, Sigma, verbose = TRUE) {
  # Check that there are no intercept rows in the clock weights
  if (any(grepl("^\\(?intercept\\)?$", cw$CpGmarker, ignore.case = TRUE))) {
    stop("[ERROR] Clock weights contain intercept rows. Please remove them before calling summarize_cpg_overlap().", call. = FALSE)
  }

  N_CpGs_ewas <- length(ss$CpGmarker)
  N_CpGs_clock <- length(cw$CpGmarker)
  N_CpGs_ref <- dim(Sigma)[1]

  N_CpGs_ewas_clock <- length(intersect(ss$CpGmarker, cw$CpGmarker))
  N_CpGs_ewas_ref <- length(intersect(ss$CpGmarker, rownames(Sigma)))
  N_CpGs_clock_ref <- length(intersect(cw$CpGmarker, rownames(Sigma)))

  cpg_intersect <- intersect(intersect(ss$CpGmarker, cw$CpGmarker), rownames(Sigma))

  N_CpGs_final <- length(cpg_intersect)

  ClockCpG_coverage_ratio <- N_CpGs_final / N_CpGs_clock
  ClockCpG_missed_ratio <- 1 - N_CpGs_final / N_CpGs_clock

  log_message(verbose, "INFO", "Number of CpGs in EWAS: %s", format(N_CpGs_ewas, big.mark = ","))
  log_message(verbose, "INFO", "Number of CpGs for clock: %s", format(N_CpGs_clock, big.mark = ","))
  log_message(verbose, "INFO", "Number of CpGs in DNAm matrix: %s", format(N_CpGs_ref, big.mark = ","))

  log_message(verbose, "INFO", "Number of overlapping CpGs between EWAS and clock: %s", format(N_CpGs_ewas_clock, big.mark = ","))
  log_message(verbose, "INFO", "Number of overlapping CpGs between EWAS and DNAm matrix: %s", format(N_CpGs_ewas_ref, big.mark = ","))
  log_message(verbose, "INFO", "Number of overlapping CpGs between clock and DNAm matrix: %s", format(N_CpGs_clock_ref, big.mark = ","))

  log_message(verbose, "INFO", "Final number of overlapping CpGs: %s", format(N_CpGs_final, big.mark = ","))
  log_message(verbose, "INFO", "Percentage of missed clock CpGs: %s%%", round(ClockCpG_missed_ratio, digits = 2) * 100)

  if (N_CpGs_final == 0) {
    stop("\nERROR! Zero overlapping CpGs found between EWAS, Clock, and Covariance Matrix.", call. = FALSE)
  }

  return(list(
    N_CpGs_ewas = N_CpGs_ewas,
    N_CpGs_clock = N_CpGs_clock,
    N_CpGs_ref = N_CpGs_ref,
    N_CpGs_ewas_clock = N_CpGs_ewas_clock,
    N_CpGs_ewas_ref = N_CpGs_ewas_ref,
    N_CpGs_clock_ref = N_CpGs_clock_ref,
    cpg_intersect = cpg_intersect,
    N_CpGs_final = N_CpGs_final,
    ClockCpG_coverage_ratio = ClockCpG_coverage_ratio,
    ClockCpG_missed_ratio = ClockCpG_missed_ratio
  ))
}

#' Cap Var(Y) at the largest value compatible with the estimated clock
#' variance
#'
#' Var(Y) cannot exceed w'Sigma w / (w'beta)^2, which is the value attained
#' when the clock and the phenotype are perfectly correlated. Exceeding it
#' produces a negative residual variance and a missing standard error. As 
#' necessary, the estimate is shrunk just below this bound.
#'
#' @param varY Current Var(Y) estimate.
#' @param w_T_beta w'beta.
#' @param w_T_Sigma_w w'Sigma w.
#' @param verbose Logical.
#'
#' @return A list with the (possibly capped) `varY`, the diagnostic
#'   `varY_ratio` (original varY / maximum allowed varY), and a `note` describing any capping.
#' @noRd
cap_phenotype_variance <- function(varY, w_T_beta, w_T_Sigma_w, verbose = TRUE) {
  note_varY_capping <- NA_character_
  varY_ratio <- NA_real_

  # Check for valid inputs
  if (is.null(varY)) {
    stop("[ERROR] varY is NULL; cannot cap Var(Y) without a current estimate.", call. = FALSE)
  } else if (!is.finite(varY)) {
    stop("[ERROR] varY is not finite; cannot cap Var(Y) without a valid estimate.", call. = FALSE)
  } else if (varY <= 0) {
    stop("[ERROR] varY is not positive; cannot cap Var(Y) without a valid estimate.", call. = FALSE)
  } else if (is.null(w_T_beta) || is.null(w_T_Sigma_w)) {
    stop("[ERROR] w_T_beta or w_T_Sigma_w is NULL; cannot cap Var(Y) without clock statistics.", call. = FALSE)
  } else if (!is.finite(w_T_beta) || !is.finite(w_T_Sigma_w)) {
    stop("[ERROR] w_T_beta or w_T_Sigma_w is not finite; cannot cap Var(Y) without valid clock statistics.", call. = FALSE)
  } else if (w_T_Sigma_w <= 0) {
    stop("[ERROR] w_T_Sigma_w is not positive; cannot cap Var(Y) when the clock has no variance.", call. = FALSE)
  }

  if (w_T_beta == 0) {
    note_varY_capping <- "w'beta = 0; cannot cap Var(Y) when the clock has no correlation with the phenotype."
    return(list(varY = varY, varY_ratio = varY_ratio, note = note_varY_capping))
  }

  varY_max <- w_T_Sigma_w / w_T_beta^2
  varY_ratio <- varY / varY_max

  if (varY > varY_max) {
    # Cap varY to just below the max allowed
    msg_capping <- sprintf(
      "Estimated varY (%.6g) exceeds maximum (%.6g), ratio = varY/varY_max = %.6f. Capped to avoid negative error variance. If the ratio is much larger than 1: check that DNAm_cov_matrix matches the EWAS population and covariate adjustment.",
      varY, varY_max, varY_ratio
    )

    log_message(verbose, "INFO", "%s", msg_capping)
    varY <- varY_max * (1 - 1e-8)
    note_varY_capping <- msg_capping
  }

  return(list(varY = varY, varY_ratio = varY_ratio, note = note_varY_capping))
}

#' Estimate gamma: the effect of the clock on the phenotype, Y = C * gamma + e
#'
#' `gamma_hat` is derived directly from `w_T_beta`, `w_T_Sigma_w`, and `varY`
#' (summary stats only): `gamma_hat = w_T_beta * varY / w_T_Sigma_w`.
#' Its standard error uses the residual variance `sigma_e2 =
#' varY - gamma_hat^2 * w_T_Sigma_w` (also returned for use in 
#' [convert_gamma_to_logistic_scale()]).
#'
#' @param w_T_beta w'beta: the clock-weighted sum of per-CpG EWAS effect
#'   sizes. Equivalent to `lambda_hat`.
#' @param w_T_Sigma_w w'Sigma w: the estimated variance of the clock, Var(C).
#' @param varY Var(Y), the capped phenotype variance.
#' @param N EWAS sample size.
#'
#' @return A list with:
#'   \itemize{
#'     \item `estimate` - the gamma point estimate.
#'     \item `variance` - the estimated sampling variance of `estimate`.
#'     \item `se` - standard error (`NaN` if `variance < 0`).
#'     \item `z` - the Z statistic (`estimate / se`).
#'     \item `p` - two-sided p-value from `z`.
#'     \item `sigma_e2` - residual variance.
#'     \item `note` - `"var_gamma_hat < 0"` if `variance` came out negative,
#'       else `NA`.
#'   }
#' @noRd
estimate_gamma <- function(w_T_beta, w_T_Sigma_w, varY, N) {
  gamma_hat <- w_T_beta * varY / w_T_Sigma_w
  sigma_e2 <- varY - gamma_hat^2 * w_T_Sigma_w
  var_gamma_hat <- sigma_e2 / ((N - 2) * w_T_Sigma_w)

  note_gamma <- ifelse(var_gamma_hat < 0, "var_gamma_hat < 0", NA)

  se_gamma_hat <- sqrt_or_na(var_gamma_hat)
  z_gamma_hat <- gamma_hat / se_gamma_hat

  return(list(
    estimate = gamma_hat,
    variance = var_gamma_hat,
    se = se_gamma_hat,
    z = z_gamma_hat,
    p = two_sided_normal_pvalue(z_gamma_hat),
    sigma_e2 = sigma_e2,
    note = note_gamma
  ))
}

#' Estimate lambda: the effect of the phenotype on the clock, C = Y * lambda + u
#'
#' The estimate `lambda_hat` is just `w_T_beta`. Its standard error uses the
#' residual variance `sigma_u2 = w_T_Sigma_w - lambda_hat^2 * varY`.
#'
#' @param w_T_beta w'beta: the clock-weighted sum of per-CpG EWAS effect
#'   sizes. Used directly as `lambda_hat`.
#' @param w_T_Sigma_w w'Sigma w: the estimated variance of the clock, Var(C).
#' @param varY Var(Y), the capped phenotype variance.
#' @param N Sample size, used for the `N - 2` degrees-of-freedom adjustment.
#'
#' @return A list with:
#'   \itemize{
#'     \item `estimate` - the lambda point estimate.
#'     \item `variance` - the estimated sampling variance of `estimate`.
#'     \item `se` - standard error (`NaN` if `variance < 0`).
#'     \item `z` - the Z statistic (`estimate / se`).
#'     \item `p` - two-sided p-value from `z`.
#'     \item `note` - `"var_lambda_hat < 0"` if `variance` came out negative,
#'       else `NA`.
#'   }
#' @noRd
estimate_lambda <- function(w_T_beta, w_T_Sigma_w, varY, N) {
  lambda_hat <- w_T_beta
  sigma_u2 <- w_T_Sigma_w - lambda_hat^2 * varY
  var_lambda_hat <- sigma_u2 / ((N - 2) * varY)

  note_lambda <- ifelse(var_lambda_hat < 0, "var_lambda_hat < 0", NA)

  se_lambda_hat <- sqrt_or_na(var_lambda_hat)
  z_lambda_hat <- lambda_hat / se_lambda_hat

  return(list(
    estimate = lambda_hat,
    variance = var_lambda_hat,
    se = se_lambda_hat,
    z = z_lambda_hat,
    p = two_sided_normal_pvalue(z_lambda_hat),
    note = note_lambda
  ))
}

#' Convert the linear-scale gamma to a log-odds-ratio (logistic) scale
#'
#' The point estimate is gamma / sigma_e^2. Note: this conversion introduces
#' some downard bias to the point estimate as the true effect gets larger.
#' 
#' The standard error follows by the delta method. The Z statistic and p-value are
#' inherited from the linear-scale gamma (H0: gamma = 0 and H0: gamma_logistic = 0
#' are the same hypothesis). Note: This means Z != coef/SE. If Z was computed from
#' the delta-method SE, it would be attenuated towards 0 when the association is strong.
#' 
#' @param gamma Result of [estimate_gamma()].
#' @param varY Var(Y).
#' @param w_T_Sigma_w w'Sigma w: the estimated variance of the clock, Var(C).
#' @param varY_ratio varY / varY_max, used to detect non-identifiability.
#' @param verbose Logical.
#'
#' @return A list with:
#' \itemize{
#'   \item `estimate` - the logistic-scale point estimate.
#'   \item `variance` - the estimated sampling variance of `estimate`, computed using the delta method.
#'   \item `se` - standard error (`NaN` if `variance < 0`).
#'   \item `z` - the Z statistic (`estimate / se`).
#'   \item `p` - two-sided p-value from `z`.
#'   \item `note` - `"var_gamma_hat_logistic < 0"` if `variance` came out negative,
#'            a note about non-identifiability of `varY` is too large, and/or
#'            a note about the delta-method SE.
#' }
#' @noRd
convert_gamma_to_logistic_scale <- function(gamma, varY, w_T_Sigma_w, varY_ratio, verbose = TRUE) {
  log_message(verbose, "INFO", "Converting final gamma estimate from linear to logistic (log-odds) scale...")

  if (!is.finite(gamma$sigma_e2) || gamma$sigma_e2 <= 0 ||
      (is.finite(varY_ratio) && varY_ratio > 0.99)) {
    note_gamma_logistic <- sprintf(
      "sigma_e2 = %.6g (varY_ratio = %.6g); logistic-scale gamma not identifiable.",
      gamma$sigma_e2, varY_ratio
    )
    log_message(verbose, "WARN", "%s", note_gamma_logistic)

    return(list(
      estimate = NA_real_, variance = NA_real_, se = NA_real_,
      z = NA_real_, p = NA_real_, note = note_gamma_logistic
    ))
  }

  gamma_hat_logistic <- gamma$estimate / gamma$sigma_e2

  # Delta method for computing SE of logsitic-scale gamma
  # var(gamma_hat_logistic) = var(gamma_hat) * (d gamma_hat_logistic / d gamma_hat)^2
  dg <- (varY + gamma$estimate^2 * w_T_Sigma_w) / gamma$sigma_e2^2 # derivative of gamma_hat_logistic with respect to gamma_hat
  var_gamma_hat_logistic <- gamma$variance * dg^2
  se_gamma_hat_logistic <- sqrt_or_na(var_gamma_hat_logistic)

  note_gamma_logistic <- paste0("SE by delta-method; Z and P inherited from the linear-scale gamma")

  if (is.finite(var_gamma_hat_logistic) && var_gamma_hat_logistic < 0) {
    note_gamma_logistic <- paste0("var_gamma_hat_logistic < 0; ", note_gamma_logistic)
  }

  log_message(verbose, "NOTE", paste0(
    "For gamma in a logistic regression of Phenotype ~ Clock: ",
    "Z and P are inherited from the linear-scale gamma; ",
    "SE computed from delta-method, so Z != coef/SE for this row. ",
    "CI should be obtained by transforming the 'linear' gamma CI endpoints."
  ))

  return(list(
    estimate = gamma_hat_logistic,
    variance = var_gamma_hat_logistic,
    se = se_gamma_hat_logistic,
    z = gamma$z,
    p = gamma$p,
    note = note_gamma_logistic
  ))
}

#' Compute the clock-level log hazard ratio.
#'
#' Convert from linear scale gamma to Cox log-hazard-ratio by dividing by 
#' the cox conversion factor `kappa_cox`. This is simply undoing the
#' conversion applied to the summary stats during [convert_cox_to_linear_scale()].
#'    `gamma_hat_cox = gamma_hat / kappa_cox`
#'    `se_gamma_hat_cox = se_gamma_hat / kappa_cox`
#' 
#' Note: whether `varY` is inferred or user-supplied determines whether the
#' `kappa_cox`/`cox_pi` dependence cancels out or not. See the details below:
#' TODO: should move details elsewhere.
#' 
#' Tracing how `kappa_cox` propogates through the pipeline in `BETA` and `SE`:
#'    1. Converting cox sumstats to linear scale multiplies `BETA` and `SE` by `kappa_cox`.
#'    2. Converting linear Y~M to linear M~Y divides both by `denom = BETA^2 + (N-2)*SE(BETA)^2`, which scales as `kappa_cox^2`.
#'    => So the linear M~Y BETA values carry a factor of `1/kappa_cox`.
#'    => So `w_T_beta` carries a factor of `1/kappa_cox`
#' 
#' Tracing how `kappa_cox` propogates through the pipeline if `varY` is inferred:
#'    1. `varY = Var(M) / (BETA^2 + (N - 2) * SE(BETA)^2)`, so it scales as `kappa_cox^2` (same as `denom` above).
#' 
#' Tracing how `kappa_cox` propogates into `gamma_hat = w_T_beta * varY / w_T_Sigma_w`:
#'    1. If `varY` is inferred, `gamma_hat` carries a factor of `kappa_cox`
#'    2. If `varY` is user-supplied, `gamma_hat` carries a factor of `1/kappa_cox`
#' 
#' Tracing `kappa_cox` into the final log hazard ratio:
#'    1. If `varY` is inferred, `log_hr` is invariant to `kappa_cox` (and hence `cox_pi`).
#'    2. If `varY` is user-supplied, `log_hr` carries a factor of `1/kappa_cox^2` (so is sensitive to `cox_pi`).
#'
#' @param gamma Result of [estimate_gamma()], scaled as necessary for Horvath and PedBE clocks
#' @param kappa_cox Cox conversion factor from [convert_cox_to_linear_scale()].
#' @param varY_type Either `"Inferred"` or `"User_provided"` -- determines
#'   whether the `kappa_cox`/`cox_pi` dependence cancels (see above).
#' @param verbose Logical.
#'
#' @return A list with:
#' \itemize{
#'   \item `estimate` - the log-hazard-ratio scale point estimate.
#'   \item `se` - standard error in log-hazard-ratio scale.
#'   \item `z` - the Z statistic inherited from gamma.
#'   \item `p` - two-sided p-value from `z`.
#'   \item `note` - information about whether the estimate is invariant to `cox_pi` or not.
#' }
#' @noRd
compute_cox_log_hazard_ratio <- function(gamma, kappa_cox, varY_type, verbose = TRUE) {
  gamma_hat_cox <- gamma$estimate / kappa_cox
  se_gamma_hat_cox <- gamma$se / kappa_cox

  if (identical(varY_type, "Inferred")) {
    note_gamma_cox <- "Clock-level log hazard ratio; invariant to cox_pi. First-order approximation, attenuated toward the null."
  } else {
    note_gamma_cox <- "Clock-level log hazard ratio. WARNING: varY was user-supplied, so the cox_pi cancellation does NOT hold; this estimate requires varY = pi(1-pi) with a correct cox_pi."
  }

  log_message(verbose, "INFO", "Clock-level log HR = %.6g (SE %.6g); HR = %.4f", gamma_hat_cox, se_gamma_hat_cox, exp(gamma_hat_cox))

  return(list(
    estimate = gamma_hat_cox,
    se = se_gamma_hat_cox,
    z = gamma$z,
    p = gamma$p,
    note = note_gamma_cox
  ))
}
