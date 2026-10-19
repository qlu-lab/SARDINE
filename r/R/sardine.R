# =============================================================================
# r/R/sardine.R: Exported, user-facing function
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
#' See <https://github.com/qlu-lab/SARDINE/blob/main/inputs.md> for details on
#' every input.
#'
#' @param ewas_sumstats EWAS summary statistics: can be a dir/to/file or a
#'   data.table, data.frame; must contain columns CpG, BETA, SE, and N.
#' @param ewas_direction Direction of the association (outcome ~ exposure) in
#'   the EWAS model: `"M~Y"` (methylation regressed on phenotype; default) or
#'   `"Y~M"` (phenotype regressed on methylation).
#' @param ewas_model EWAS regression model: `"linear"` (default), `"logistic"`,
#'   or `"cox"` (Cox proportional hazards).
#' @param ewas_DNAm_scale Scale of the EWAS DNA methylation data: `"beta"`
#'   (beta values; default) or `"m"` (M-values).
#' @param ewas_pheno_type EWAS phenotype type: `"continuous"` (default) or
#'   `"binary"`. If binary, will also output an estimate of the log-odds ratio
#'   of clock on phenotype.
#'
#'   For all four options above, matching ignores case and whitespace
#'   (e.g. `"y ~ m"` is accepted as `"Y~M"`).
#' @param clock Epigenetic clock weights: can be a dir/to/clock/file or a 
#'   data.frame/data.table; must contain columns CpG, CoefficientTraining.
#' @param DNAm_cov_matrix Variance-covariance matrix for DNAm. Must have the
#'   same covariate adjustments as the ewas summary stats. Can be a
#'   dir/to/file of a delimited table (e.g. TSV) whose first column holds the
#'   CpG IDs and whose header lists the same CpG IDs in the same order; a
#'   dir/to/file.rds containing a matrix; a matrix with CpG rownames and
#'   colnames; a data.frame of numeric columns with CpG rownames; or a
#'   data.frame/data.table laid out like the delimited file.
#' @param mean_beta Optional. dir/to/file containing mean DNAm beta-value at each CpG
#'   site (expected to be in beta-scale). Needed for M-value to beta-value conversion.
#'   A data.frame or data.table with 2 columns: CpGmarker, beta_bar. Best if computed
#'   from the same samples as the EWAS itself. Default: `NULL`; only needed for M-value
#'   to beta-value conversion.
#' @param varY Variance of the phenotype used in EWAS: raw phenotype if no
#'   covariates; residualized phenotype if adjusted for covariates in the EWAS
#'   model. Default: `NULL`; we can estimate varY using only EWAS sumstats and
#'   DNAm reference matrix.
#' @param scale Logical. For Horvath or PedBE clocks, we can scale the point estimate
#'   and SE; doesn't affect Pvalue.
#' @param age_adult Only used if `scale = TRUE`.
#' @param N Optional. EWAS total sample size. Default is to read this from the
#'   EWAS summary statistics.
#' @param cox_pi Only used when `ewas_model = "cox"`. Proportion of cases at the end of the study.
#'   If not provided, this is computed as `cox_E/N`. Exactly one of `cox_pi` and
#'   `cox_E` must be provided if `ewas_model` = "cox".
#' @param cox_E Only used when `ewas_model = "cox"`. Number of cases at the end of the study.
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
#' @return An object of class `sardine_result`: a `data.frame` with one row per
#'   estimated parameter (`lambda` and `gamma` always; `gamma_logistic` for binary
#'   phenotypes; `log_HR` for Cox EWAS) and columns:
#'   \describe{
#'     \item{model}{EWAS model used, e.g. `"Phenotype ~ Clock (Cox)"`.}
#'     \item{param}{Parameter being estimated, e.g. `"log_HR"`.}
#'     \item{coef}{Point estimate for the parameter.}
#'     \item{SE}{Estimated standard error for the point estimate.}
#'     \item{z}{Z statistic (estimate / SE). All rows share the z of the linear
#'       `gamma` test; for `gamma_logistic` the SE comes from the delta method,
#'       so `z != coef / SE` for that row.}
#'     \item{p}{Two-sided normal p-value for the estimate, testing the null
#'       hypothesis that the effect is 0.}
#'     \item{neglog10p}{-log10(p), computed directly from `z` so it stays finite
#'       when very small p-values underflow to 0.}
#'     \item{N}{EWAS sample size: the `N` argument if supplied, otherwise the
#'       median per-CpG N across the CpGs used.}
#'     \item{note}{Notes about the specific settings for this row (e.g. log-odds
#'       vs log_HR), or `NA`.}
#'   }
#'   Retrieve the metadata with `attr(result, "metadata")`: a one-row
#'   `data.frame` with columns:
#'   \describe{
#'     \item{sardine_version}{Version of SARDINE that produced the result.}
#'     \item{N_CpGs_ewas, N_CpGs_clock, N_CpGs_DNAmRef}{Number of CpGs in the
#'       EWAS, clock, and DNAm reference respectively.}
#'     \item{N_CpGs_ewas_clock, N_CpGs_ewas_ref, N_CpGs_clock_ref}{Number of CpGs
#'       overlapping between each pair of inputs.}
#'     \item{N_CpGs_final}{Number of CpGs in all three inputs; these are the CpGs
#'       actually used for the association analysis.}
#'     \item{ClockCpG_coverage_ratio}{Fraction of the clock's CpGs that ended up
#'       in the final set of CpGs used.}
#'     \item{ClockCpG_missed_ratio}{1 - `ClockCpG_coverage_ratio`.}
#'     \item{cox_pi}{Proportion of cases at the end of the study, either input
#'       directly or computed as `cox_E / median(N)`. `NA` if the EWAS model is
#'       not Cox.}
#'     \item{kappa_cox}{Cox conversion factor `(1 - cox_pi) * -log(1 - cox_pi)`
#'       used to convert between the linear-scale gamma and the Cox
#'       log hazard ratio. `NA` if the EWAS model is not Cox.}
#'     \item{varY}{Variance of the phenotype (after any capping), either input
#'       directly or estimated from the EWAS summary statistics and the DNAm
#'       reference matrix.}
#'     \item{varY_type}{`"user_provided"` or `"inferred"`.}
#'     \item{varY_ratio}{(varY prior to capping) / (theoretical maximum varY),
#'       where the maximum is the varY you'd get if the clock explained all of
#'       the phenotype's variance. Values above 1 trigger capping. `NA` if the
#'       clock has no correlation with the phenotype.}
#'     \item{varY_note}{Note about whether varY had to be capped at the
#'       theoretical maximum, or `NA` if no capping was needed.}
#'   }
#'
#' @examples
#' # Small synthetic inputs bundled with the package
#' inputs <- system.file("extdata", "fixtures", "inputs", package = "sardine")
#'
#' result <- sardine(
#'   ewas_sumstats   = file.path(inputs, "ewas_sumstats_baseline_synthetic.tsv"),
#'   clock           = file.path(inputs, "clock_weights_baseline_synthetic.tsv"),
#'   DNAm_cov_matrix = file.path(inputs, "sigma_baseline_synthetic.tsv"),
#'   verbose         = FALSE
#' )
#' result
#' attr(result, "metadata")
#'
#' @export
sardine <- function(ewas_sumstats,
                    ewas_direction = c("M~Y", "Y~M"),
                    ewas_model = c("linear", "logistic", "cox"),
                    ewas_DNAm_scale = c("beta", "m"),
                    ewas_pheno_type = c("continuous", "binary"),
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

  log_message(verbose, NULL, "%s", paste(
    "************************************************",
    sprintf("* SARDINE v%s", sardine_version()),
    "* (C) Ryan Moreno <rrmoreno@wisc.edu> and",
    "*     Yuchang Wu <ywu423@wisc.edu>, ",
    "*     University of Wisconsin-Madison, 2026",
    "* GNU General Public License v3",
    "************************************************",
    sep = "\n"
  ))

  log_message(verbose, "INFO", "Start time: %s", format(start_time, "%Y-%m-%d %H:%M:%S"))

  # ============================================================
  # Check the validity of the parameters
  # ============================================================

  # Resolve each option to its canonical value (case and whitespace are ignored)
  ewas_direction <- match_option(ewas_direction, c("M~Y", "Y~M"), "ewas_direction")
  ewas_model <- match_option(ewas_model, c("linear", "logistic", "cox"), "ewas_model")
  ewas_DNAm_scale <- match_option(ewas_DNAm_scale, c("beta", "m"), "ewas_DNAm_scale")
  ewas_pheno_type <- match_option(ewas_pheno_type, c("continuous", "binary"), "ewas_pheno_type")

  # mean_beta must be provided if ewas_DNAm_scale is m-value
  if (ewas_DNAm_scale == "m" && is.null(mean_beta)) {
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
      log_message(verbose, "WARN", "ewas_model = 'cox' with user-supplied varY: `kappa_cox` no longer cancels, so `log_hr` is sensitive to `cox_pi`. Leaving varY unset is recommended.")
    }
  }

  # if ewas_model is logistic, ewas_direction must be Y~M and ewas_pheno_type must be binary
  if (ewas_model == "logistic" && ewas_direction != "Y~M") {
    stop("[ERROR] 'ewas_direction' must be 'Y~M' when 'ewas_model' is 'logistic'.", call. = FALSE)
  }
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
      log_message(verbose, "WARN", "Your clock weight filename contains '%s' yet you did not scale the point estimate and SE (Z-score and P-value are not affected). Make sure this was intentional.", clock_name)
    }
  }

  # Check numeric types
  if (!is.null(N) && !is.numeric(N)) {
    stop("[ERROR] 'N' must be numeric.", call. = FALSE)
  }
  if (!is.null(cox_E) && !is.numeric(cox_E)) {
    stop("[ERROR] 'cox_E' must be numeric.", call. = FALSE)
  }
  if (!is.null(cox_E) && !(cox_E > 0)) {
    stop("[ERROR] 'cox_E' must be greater than 0.", call. = FALSE)
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
  if (!is.null(varY)) {
    log_message(verbose, "WARN", "User-supplied varY must be the variance of the phenotype after adjusting for the same covariates as the EWAS (the raw phenotype variance only if the EWAS had no covariates). Leave varY unset to estimate it from the EWAS summary statistics instead.")
  }

  # ============================================================
  # Read and standardize input tables
  # ============================================================
  ss <- read_summary_statistics(ewas_sumstats, N = N, verbose = verbose)
  cw <- read_clock_weights(clock, verbose = verbose)
  Sigma <- read_dnam_covariance(DNAm_cov_matrix, verbose = verbose)
  mean_beta_df <- if (!is.null(mean_beta)) read_mean_beta(mean_beta, verbose = verbose) else NULL


  # ============================================================
  # M-value to Beta-value transformation, as necessary
  # ============================================================
  if (ewas_DNAm_scale == "m") {
    ss <- convert_mvalue_ewas_to_beta_scale(ss, mean_beta_df, ewas_direction, verbose = verbose)
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
  varY_type <- ifelse(!is.null(varY), "user_provided", NA)
  if (ewas_model == "logistic" && is.null(varY)) {
    varY <- infer_phenotype_variance_from_logistic_y_on_m_ewas(ss, Sigma, verbose = verbose)
    varY_type <- "inferred"
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
    varY_type <- "inferred"
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
      log_message(verbose, "NOTE", "EWAS is a linear regression with a binary phenotype; we can still estimate log-odds ratio of Clock on this binary phenotype.")
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
  # This is necessary for Horvath1, Horvath2, and PedBE clocks
  # because they apply a scale factor to their clocks.
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

    log_message(verbose, "INFO", "Applied clock scaling factor = %s to point estimate and SE", scale_factor)
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
      z = lambda$z,
      p = lambda$p,
      N = N,
      note = lambda$note,
      stringsAsFactors = FALSE
    ),
    data.frame(
      model = "Phenotype ~ Clock (linear)",
      param = "gamma",
      coef = gamma$estimate,
      SE = gamma$se,
      z = gamma$z,
      p = gamma$p,
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
      z = gamma_logistic$z,
      p = gamma_logistic$p,
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
      z = cox_log_hr$z,
      p = cox_log_hr$p,
      N = N,
      note = cox_log_hr$note,
      stringsAsFactors = FALSE
    ))
  }

  # Add a -log10p column to avoid underflow issues with very small p-values
  res$neglog10p <- two_sided_neg_log10_pvalue(res$z)
  res <- res[, c("model", "param", "coef", "SE", "z", "p", "neglog10p", "N", "note")]


  # Add metadata (number of CpGs, varY diagnostics, etc.)
  res_metadata <- data.frame(
    sardine_version = sardine_version(),

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

  log_message(verbose, NULL, "")
  log_message(verbose, "INFO", "Analysis complete!")
  log_message(verbose, "INFO", "End time: %s", format(end_time, "%Y-%m-%d %H:%M:%S"))
  log_message(verbose, "INFO", "Elapsed time: %s", format_elapsed_time(elapsed_seconds))
  log_message(verbose, NULL, "")

  return(res)
}
