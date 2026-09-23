# =============================================================================
# R/src/sardine.R: Exported, user-facing function
# =============================================================================

#' Summary-statistics-based association between epigenetic clock and phenotype
#'
#' Estimates the association between an epigenetic clock and a phenotype
#' using only EWAS summary statistics, the clock's CpG weights, and a
#' reference DNA methylation covariance matrix. No individual-level data is
#' required.
#'
#' Two coefficients are always returned: `lambda`, from the regression of
#' the clock on the phenotype, and `gamma`, from the regression of the
#' phenotype on the clock. Binary phenotypes additionally yield
#' `gamma_logistic`, a log-odds ratio. Cox-model inputs additionally
#' yield `log_HR`, a clock-level log hazard ratio. All rows share the same Z
#' statistic and p-value by construction.
#'
#' @param ewas_sumstats EWAS summary statistics: can be a dir/to/file or a
#'   data.table, data.frame; must contain columns CpG, BETA, SE, and N.
#' @param ewas_direction Defines the direction of the association (exposure 
#'   and outcome) in the EWAS model. Either M~Y (default) or Y~M.
#' @param ewas_model Whether the input sumstats are linear (default), logistic
#'   regression, or Cox (Cox proportional hazards model).
#' @param ewas_DNAm_scale Whether the EWAS DNA methylation data are
#'  beta values (default) or m values.
#' @param ewas_pheno_type Whether the EWAS phenotype is binary or continuous.
#'   If binary, will also output an estimate of the log-odds ratio of clock on
#'   phenotype.
#' @param clock Epigenetic clock weights: can be a dir/to/clock/file or a 
#'   data.frame/data.table; must contain columns CpG, CoefficientTraining.
#' @param DNAm_cov_matrix Variance-covariance matrix for DNAm. Must have the
#'   same covariate adjustments as the ewas summary stats.
#' @param mean_beta [optional] dir/to/file containing mean DNAm beta-value at each CpG
#'   site (expected to be in beta-scale). Needed for M-value to beta-value conversion.
#'   A data.frame or data.table with 2 columns: CpGmarker, beta_bar. Best if can match with
#'   the EWAS samples. Default: `NULL`; only needed for M-value to beta-value conversion.
#' @param varY Variance of the phenotype used in EWAS: raw phenotype if no
#'   covariates; residualized phenotype if adjusted for covariates in the EWAS
#'   model. Default: `NULL`; we can estimate varY using only EWAS sumstats and
#'   DNAm reference matrix.
#' @param scale Logical. For Horvath or PedBE clocks, we can scale the point estimate
#'   and SE; doesn't affect Pvalue.
#' @param age_adult Only used if `scale = TRUE`.
#' @param N [optional] EWAS total sample size. Default is to read this from the
#'   EWAS summary statistics.
#' @param cox_pi [`ewas_model` = "cox" only]. Proportion of cases at the end of the study.
#'   If not provided, this is computed as `cox_E/N`. Exactly one of `cox_pi` and
#'   `cox_E` must be provided if `ewas_model` = "cox".
#' @param cox_E [`ewas_model` = "cox" only]. Number of cases at the end of the study.
#'   Default `NULL`. Exactly one of `cox_pi` and `cox_E` must be provided
#'   if `ewas_model` = "cox".
#' @param verbose True/False. Whether to print progress and diagnostic messages.
#'
#' @section Cox models:
#' For a Cox Proportional Hazards model, `ewas_direction` must be `"Y~M"`, `ewas_model` must
#' be `"cox"`, `ewas_pheno_type` must be `"binary"`, and `N` should represent the
#' total sample size at the end of the study. Exactly one of `cox_pi` and `cox_E` must be provided.
#' `cox_pi = cox_E/N` is the proportion of cases at the end of the study. It is recommended that `varY`
#' is inferred rather than provided as user-supplied so that the `log_hr` estimate is invariant to `cox_pi`.
#' 
#' @section Logistic models:
#' for a Logistic model, `ewas_direction` must be `"Y~M"`, `ewas_model` must be `"logistic"`,
#' and `ewas_pheno_type` must be `"binary"`.
#'
#' @return An object of class `sardine_result`: a `data.frame` of coefficients (lambda, gamma,
#'   gamma_logistic, log_HR) with a `metadata` attribute carrying CpG overlap counts
#'   and Var(Y) diagnostics. Retrieve the metadata with `attr(result, "metadata")`.
#'   TODO: be more specific about the columns in the output and metadata
#'
#' @examples
#' \dontrun{
#' result <- sardine(
#'   ewas_sumstats   = "sumstats/ewas_sumstats_Age_std.txt",
#'   clock           = "clocks/Horvath1.csv",
#'   DNAm_cov_matrix = "sigma/Sigma_Age_std_N=10880.rds"
#' )
#' attr(result, "metadata")
#' }
#'
#' @export
sardine <- function(ewas_sumstats,
                    ewas_direction = "M~Y",
                    ewas_model = "linear",
                    ewas_DNAm_scale = "beta",
                    ewas_pheno_type = "continuous",
                    clock,
                    DNAm_cov_matrix,
                    mean_beta = NULL,
                    varY = NULL,
                    scale = FALSE,
                    age_adult = 20,
                    N = NULL,
                    cox_pi = NULL,
                    cox_E = NULL,
                    verbose = TRUE) {

  # Record start time
  start_time <- Sys.time()

  if (verbose) {
    cat("\n\n************************************************\n")
    cat("* SARDINE\n")
    cat("* (C) 2026 Ryan Moreno <rrmoreno@wisc.edu> and\n")
    cat("*     Yuchang Wu <ywu423@wisc.edu>, University of Wisconsin-Madison\n")
    cat("* GNU General Public License v3\n")
    cat("************************************************\n\n")

    cat("[INFO] Start time:", format(start_time, "%Y-%m-%d %H:%M:%S"), "\n")
  }

  # ============================================================
  # Check the validity of the parameters
  # ============================================================

  # ewas_direction: must be either Y~M or M~Y:
  if (!toupper(gsub(" ", "", ewas_direction)) %in% c("Y~M", "M~Y")) {
    stop("[ERROR] 'ewas_direction' must be either \"M~Y\" or \"Y~M\".", call. = FALSE)
  }
  ewas_direction <- toupper(gsub(" ", "", ewas_direction))

  # ewas_model must be either linear, logistic, or cox:
  if (!ewas_model %in% c("linear", "logistic", "cox")) {
    stop("[ERROR] 'ewas_model' must be one of \"linear\", \"logistic\" or \"cox\".", call. = FALSE)
  }

  # ewas_DNAm_scale must be either beta or m
  if (!tolower(ewas_DNAm_scale) %in% c("beta", "beta-value", "m-value", "mvalue", "m")) {
    stop("[ERROR] 'ewas_DNAm_scale' must denote either a beta-value or an M-value scale.", call. = FALSE)
  }

  # ewas_pheno_type must be either continuous or binary
  if (!ewas_pheno_type %in% c("continuous", "binary")) {
    stop("[ERROR] 'ewas_pheno_type' must be either \"continuous\" or \"binary\".", call. = FALSE)
  }

  # mean_beta must be provided if ewas_DNAm_scale is m-value
  if (tolower(ewas_DNAm_scale) %in% c("m-value", "mvalue", "m") && is.null(mean_beta)) {
    stop("[ERROR] 'mean_beta' must be provided when 'ewas_DNAm_scale' is an M-value scale.", call. = FALSE)
  }

  # if ewas_model is cox, ewas_direction must be Y~M, ewas_pheno_type must be binary,
  # and exactly one of cox_pi and cox_E must be set
  if (ewas_model == "cox") {
    if (ewas_direction != "Y~M") {
      stop("[ERROR] 'ewas_direction' must be 'Y~M' when 'ewas_model' is 'cox'.", call. = FALSE)
    }
    if (ewas_pheno_type != "binary") {
      stop("[ERROR] 'ewas_pheno_type' must be 'binary' when 'ewas_model' is 'cox'.", call. = FALSE)
    }
    if (is.null(cox_pi) == is.null(cox_E)) {
      stop("[ERROR] Exactly one of 'cox_pi' or 'cox_E' must be provided when 'ewas_model' is 'cox'.", call. = FALSE)
    }

    # Warn that it's recommended varY be inferred for Cox models
    if (!is.null(varY)) {
      if(verbose) cat("[WARN] ewas_model = 'cox' with user-supplied varY: `kappa_cox` no longer cancels, so `log_hr` is sensitive to `cox_pi`. Setting varY = NULL is recommended.\n")
    }
  }

  # if ewas_model is logistic, ewas_pheno_type must be binary
  if (ewas_model == "logistic" && ewas_pheno_type != "binary") {
    stop("[ERROR] 'ewas_pheno_type' must be 'binary' when 'ewas_model' is 'logistic'.", call. = FALSE)
  }

  # scale and age_adult should be set together
  if (is.null(scale) != is.null(age_adult)) {
    stop("[ERROR] 'scale' and 'age_adult' should be set together.", call. = FALSE)
  }

  # Alert the user if they seem to be using a Horvath or PedBE clock based
  # on the file name but did not scale
  if (is.character(clock)) {
    clock_name <- gsub("\\.csv", "", basename(clock))

    if (verbose & !scale & grepl("horvath|pedbe", tolower(clock_name))) {
      cat("[WARN] Your clock weight filename contains '", clock_name,
          "' yet you did not scale the point estimate and SE (Z-score and P-value are not affected). Make sure this was intentional.\n", sep = "")
    }
  }

  # Check numeric types
  if (!is.null(N) && !is.numeric(N)) {
    stop("[ERROR] 'N' must be numeric.", call. = FALSE)
  }
  if (!is.null(cox_E) && !is.numeric(cox_E)) {
    stop("[ERROR] 'cox_E' must be numeric.", call. = FALSE)
  }
  if (!is.null(cox_pi) && !is.numeric(cox_pi)) {
    stop("[ERROR] 'cox_pi' must be numeric.", call. = FALSE)
  }
  if (!is.null(age_adult) && !is.numeric(age_adult)) {
    stop("[ERROR] 'age_adult' must be numeric.", call. = FALSE)
  }
  if (!is.null(varY) && !is.numeric(varY)) {
    stop("[ERROR] 'varY' must be numeric.", call. = FALSE)
  }
  if (!is.null(varY) && varY < 0) {
    stop("[ERROR] 'varY' must be non-negative", call. = FALSE)
  }

  

  # ============================================================
  # Read and standardize input tables
  # ============================================================
  ss <- read_summary_statistics(ewas_sumstats, N = N, verbose = verbose)
  cw <- read_clock_weights(clock, verbose = verbose)
  Sigma <- read_dnam_covariance(DNAm_cov_matrix, verbose = verbose)


  # ============================================================
  # M-value to Beta-value transformation, as necessary
  # ============================================================
  if (tolower(ewas_DNAm_scale) %in% c("m-value", "mvalue", "m")) {
    ss <- convert_mvalue_ewas_to_beta_scale(ss, mean_beta, ewas_direction, verbose = verbose)
  }

  # ============================================================
  # Convert Cox to linear Y ~ M (Cox input only)
  # ewas_direction for cox must be "Y~M"
  # ============================================================
  kappa_cox <- NA_real_ # scaling factor for Cox model

  if (ewas_model == "cox") {
    cox_conversion <- convert_cox_ewas_to_linear_scale(
      ss, cox_pi, cox_E, verbose = verbose
    )
    ss <- cox_conversion$ss
    cox_pi <- cox_conversion$cox_pi
    kappa_cox <- cox_conversion$kappa_cox
  }

  # ============================================================
  # Estimate Var(Y) for logistic regression if not provided by the user.
  # Note: In the logistic case, it is more accurate to estimate Var(Y) before
  # converting the summary stats to linear M~Y.
  # ============================================================
  varY_type <- ifelse(!is.null(varY), "User_provided", NA)
  if (ewas_model == "logistic" && is.null(varY)) {
    # TODO: previous version always did this inference, even if varY was passed in from user. Is that necessary?
    varY <- infer_phenotype_variance_from_logistic_y_on_m_ewas(ss, Sigma, verbose = verbose)
    varY_type <- "Inferred"
  }

  # ============================================================
  # Convert Y ~ M to M ~ Y, as necessary
  # ============================================================
  if (ewas_direction == "Y~M") {
    if (ewas_model %in% c("linear", "cox")) {
      # note: Cox sumstats have already been converted to linear scale
      ss <- convert_linear_y_on_m_ewas_to_m_on_y(ss, verbose = verbose)
    } else if (ewas_model == "logistic") {
      ss <- convert_logistic_y_on_m_ewas_to_linear_m_on_y(ss, Sigma, varY, verbose = verbose)
    } else {
      stop("Error: 'ewas_model' must be either 'linear', 'logistic', or 'cox'.", call. = FALSE)
    }
  }

  # ============================================================
  # Estimate Var(Y) if not provided or already estimated
  # ============================================================
  if (is.null(varY)) {
    # note: all EWAS sumstats are now formatted as linear M~Y
    varY <- infer_phenotype_variance_from_linear_m_on_y_ewas(ss, Sigma, verbose=verbose)
    varY_type <- "Inferred"
  }

  # ============================================================
  # Align CpGs and compute clock/phenotype regression terms
  # ============================================================

  # Use the CpG intersection between EWAS sumstats, clock weights, and Sigma
  cw <- drop_intercept_rows(cw, verbose = verbose)
  overlap <- summarize_cpg_overlap(ss, cw, Sigma, verbose = verbose)
  cpg_intersect <- overlap$cpg_intersect

  # Align data to the same CpG order
  ss <- ss[match(cpg_intersect, ss$CpGmarker), ]
  cw <- cw[match(cpg_intersect, cw$CpGmarker), ]
  Sigma <- Sigma[cpg_intersect, cpg_intersect]

  # Find the sample size of the EWAS summary stats if not provided by the user
  if (is.null(N)) {
    # median N for the CpG sites we're interested in
    N <- stats::median(ss$N, na.rm = TRUE) 
  }

  # Compute clock-phenotype regression coefficient
  w <- cw$CoefficientTraining
  w_T_beta <- sum(w * ss$BETA)

  # Compute the variance of the clock sum on the reference population
  w_T_Sigma_w <- as.numeric(t(w) %*% Sigma %*% w)
  if (!is.finite(w_T_Sigma_w) || w_T_Sigma_w <= 0) {
    stop(sprintf("[ERROR] Estimated clock variance w'Sigma w = %.6g is not positive. Check that DNAm_cov_matrix is a valid covariance matrix on the same CpG scale as the EWAS.", w_T_Sigma_w), call. = FALSE)
  }

  # Cap the varY (otherwise, this can result in negative variance estimates for gamma and lambda)
  capping <- cap_phenotype_variance(varY, w_T_beta, w_T_Sigma_w, verbose = verbose)
  varY <- capping$varY
  varY_ratio <- capping$varY_ratio

  # ============================================================
  # Linear regression: Y = C * gamma + e
  # ============================================================
  gamma <- estimate_gamma(w_T_beta, w_T_Sigma_w, varY, N)

  # ============================================================
  # Logistic regression Y ~ C (only if phenotype is binary)
  # ============================================================
  if (ewas_pheno_type == "binary") {
    if (ewas_model == "linear") {
      if (verbose) cat("[NOTE] EWAS is a linear regression with a binary phenotype; we can still estimate log-odds ratio of Clock on this binary phenotype.\n")
    }
    gamma_logistic <- convert_gamma_to_logistic_scale(gamma, varY, w_T_Sigma_w, varY_ratio, verbose = verbose)

    # Add a note to clarify log odds for cumulative incidence versus log hazard ratio for Cox model
    if (ewas_model == "cox") {
      gamma_logistic$note <- append_note(
        gamma_logistic$note,
        "Log-odds ratio for cumulative incidence at end of follow-up; param = 'log_HR' is log of hazard ratio."
      )
    }
  } else {
    gamma_logistic <- list(estimate = NA_real_, variance = NA_real_, se = NA_real_,
                            z = NA_real_, p = NA_real_, note = NA_character_)
  }

  # ============================================================
  # Linear regression: C = Y * lambda + u
  # ============================================================
  lambda <- estimate_lambda(w_T_beta, w_T_Sigma_w, varY, N)

  # ============================================================
  # Re-scale effect size and SE, as necessary
  # This is necessary for Horvath1, Horvath2, and PedBE clocks.
  # TODO: describe why this is necessary
  # ============================================================
  if (scale) {
    scale_factor <- 1 + age_adult

    gamma$estimate <- gamma$estimate / scale_factor
    gamma$variance <- gamma$variance / (scale_factor^2)
    gamma$se <- gamma$se / scale_factor

    gamma_logistic$estimate <- gamma_logistic$estimate / scale_factor
    gamma_logistic$variance <- gamma_logistic$variance / (scale_factor^2)
    gamma_logistic$se <- gamma_logistic$se / scale_factor

    lambda$estimate <- lambda$estimate * scale_factor
    lambda$variance <- lambda$variance * (scale_factor^2)
    lambda$se <- lambda$se * scale_factor

    if (verbose) {
      cat("[INFO] Applied clock scaling factor =", scale_factor, "to point estimate and SE", "\n")
    }
  }

  # ============================================================
  # Clock-level log hazard ratio (Cox input only)
  # ============================================================
  if (ewas_model == "cox") {
    cox_log_hr <- compute_cox_log_hazard_ratio(gamma, kappa_cox, varY_type, verbose = verbose)
  } else {
    cox_log_hr <- list(estimate = NA_real_, se = NA_real_, z = NA_real_, p = NA_real_, note = NA_character_)
  }

  # ============================================================
  # Format the regression results into a sardine_result object
  # ============================================================
  res <- rbind(
    data.frame(
      model = "Clock ~ Phenotype (linear)",
      param = "lambda",
      coef = lambda$estimate,
      SE = lambda$se,
      Z = lambda$z,
      P = lambda$p,
      N = N,
      note = lambda$note,
      stringsAsFactors = FALSE
    ),
    data.frame(
      model = "Phenotype ~ Clock (linear)",
      param = "gamma",
      coef = gamma$estimate,
      SE = gamma$se,
      Z = gamma$z,
      P = gamma$p,
      N = N,
      note = gamma$note,
      stringsAsFactors = FALSE
    )
  )

  if (!is.na(gamma_logistic$estimate)) {
    res <- rbind(res, data.frame(
      model = "Phenotype ~ Clock (logistic)",
      param = "gamma_logistic",
      coef = gamma_logistic$estimate,
      SE = gamma_logistic$se,
      Z = gamma_logistic$z,
      P = gamma_logistic$p,
      N = N,
      note = gamma_logistic$note,
      stringsAsFactors = FALSE
    ))
  }

  # Note: Cox pathway's robustness ordering is log_HR > gamma >> lambda, gamma_logistic
  # they all have the same Z score and p-value though
  if (!is.na(cox_log_hr$estimate)) {
    res <- rbind(res, data.frame(
      model = "Phenotype ~ Clock (Cox)",
      param = "log_HR",
      coef = cox_log_hr$estimate,
      SE = cox_log_hr$se,
      Z = cox_log_hr$z,
      P = cox_log_hr$p,
      N = N,
      note = cox_log_hr$note,
      stringsAsFactors = FALSE
    ))

    # Add columns for hazard ratio and 95% CI, defined only for the log_HR row.
    # TODO: probably simpler to add this in the previous code chunk res <- ...
    is_hr <- res$param == "log_HR"
    z_crit <- stats::qnorm(0.975)
    res$HR <- ifelse(is_hr, exp(res$coef), NA_real_)
    res$HR_L95 <- ifelse(is_hr, exp(res$coef - z_crit * res$SE), NA_real_)
    res$HR_U95 <- ifelse(is_hr, exp(res$coef + z_crit * res$SE), NA_real_)
  }

  # Add a -log10p column to avoid underflow issues with very small p-values
  res$log10_P <- two_sided_neg_log10_pvalue(res$Z)

  # Add metadata (number of CpGs, varY diagnostics, etc.)
  res_metadata <- data.frame(
    N_CpGs_ewas = overlap$N_CpGs_ewas,
    N_CpGs_clock = overlap$N_CpGs_clock,
    N_CpGs_DNAmRef = overlap$N_CpGs_ref,

    N_CpGs_ewas_clock = overlap$N_CpGs_ewas_clock,
    N_CpGs_ewas_ref = overlap$N_CpGs_ewas_ref,
    N_CpGs_clock_ref = overlap$N_CpGs_clock_ref,

    N_CpGs_final = overlap$N_CpGs_final,
    ClockCpG_coverage_ratio = overlap$ClockCpG_coverage_ratio,
    ClockCpG_missed_ratio = overlap$ClockCpG_missed_ratio,

    cox_pi = ifelse(ewas_model == "cox", cox_pi, NA_real_),
    kappa_cox = ifelse(ewas_model == "cox", kappa_cox, NA_real_),

    varY = varY,
    varY_type = varY_type,
    varY_ratio = varY_ratio, # varY/varY_max
    varY_note = capping$note
  )
  attr(res, "metadata") <- res_metadata
  class(res) <- c("sardine_result", class(res))

  # Record end time and calculate elapsed time
  end_time <- Sys.time()
  elapsed_seconds <- as.numeric(difftime(end_time, start_time, units = "secs"))

  if (verbose) {
    cat("\n")
    cat("[INFO] Analysis complete!\n")
    cat("[INFO] End time:", format(end_time, "%Y-%m-%d %H:%M:%S"), "\n")
    cat("[INFO] Elapsed time:", format_elapsed_time(elapsed_seconds), "\n")
    cat("\n")
  }

  return(res)
}
