# library(tidyverse) # From original code, but wasn't used
library(data.table)

cat("\nSARDINE v.1.2.0 loading ...\n")

# utility functions: =====

standardize_column_name = function(dt, pattern, standard_name, col_desc) {
  # "standardize" EWAS sumstats or clock weights file column names
  # will change the column names inplace
  
  col_names = names(dt)
  matched_cols = grep(pattern, col_names, ignore.case = TRUE, value = TRUE)
  
  if (length(matched_cols) == 1) {
    data.table::setnames(dt, matched_cols, standard_name)
  } else if (length(matched_cols) > 1) {
    stop(paste0("\nERROR! Ambiguous ", col_desc, " columns found: ", 
                paste(matched_cols, collapse = ", "), ". Please rename manually."))
  } else {
    stop(paste0("\nERROR! Missing ", col_desc, " column. Looked for pattern: ", pattern))
  }
}

# append notes without producing the string "NA; message"
append_note = function(x, msg){
  if(is.na(msg)) return(x)
  ifelse(is.na(x), msg, paste0(x, "; ", msg))
}

# for testing codes:
# ewas = "ewas_sumstats_Age_std.txt"
# clock = "../clocks/Horvath1.csv"
# Sigma = "sigma/Sigma_Age_std_N=10880.rds"
# N = 10880
# varY = 1

# logs:
# 
# TBD: compare with R package
#
# - 20260822: when input cox sumstats, report estimate for log HR from individual-level h(t|C)
# - 20260817: fixed a potentially misleading result in the output log10P value for gamma_logistic.
#             now also exports Z score column and corresponding -log10Pvalue
#             now the outputs Z scores will all be identical
#             before this change: you're fine if you used the z-score for "lambda"!
# - 20260816: less biased conversion from logistic EWAS sumstats to linear EWAS sumstats
#             removed the K parameter as it's not useful when there're covariates and misleading
# - 20260813: fixed a bug which did not allow "cox"; and some minor improvements (N --> N - 2)
#             added parameter "K" for population prevalence which is relevant for binary phenotypes (TBD)
# - 20260807: address the case when the estimated varY is too large leading to negative NA SE
# - 20260731: allow EWAS to be Cox Proportional Hazard model
# - 20260726: allow EWAS DNAm in M-value 
# - 20260725: allow EWAS to be logistic Y ~ M; allows for logistic Y ~ Clock
# - 20260724: allow EWAS to be linear Y ~ M
# - 20260716: non-linear clock such as PedBE and Horvath: revise the entire derivations
#             linear part of Horvath and PedBE clocks: need to scale by 21 (default):
# 
#  - 20260624: 
#       clarified that should use EWAS N for the CY regression derivations
#       add a helper function to "standardize" input file column names
#  - 20260621: can now estimate Var(Y) using EWAS sumstats and DNAm matrix
#  - 20260620: version 0.0

# usages:
# EWAS sumstats file: must contain columns BETA, SE, and N
# clock file:         must contain CpGmarker column for CpG site and 
#                     CoefficientTraining column for corresponding clock weight
# DNAm_cov_matrix:    CpG covariance matrix; rownames and colnames must be CpG sites
# varY:               variance of the phenotype used in EWAS:
#                         raw phenotype if no covariates; 
#                         residulized phenotype if adjusted for covariates in EWAS model 
#                     default = NULL;
#                     we can infer varY using only EWAS sumstats and DNAm reference matrix
# scale:              For Horvarth or PedBE clocks, can scale the point estimate and SE; doesn't affect Pvalue
# age_adult:          only used if scale = TRUE


#' Summary Statistics Based Inference on Epigenetic Aging
#' 
#' @param ewas_sumstats EWAS summary statistics; can be a dir/to/file or a data.table, data.frame; must contain columns BETA, SE, and N
#' @param ewas_direction Who are the exposure and outcome in EWAS? Either M~Y (default) or Y~M.
#' @param ewas_model Input sumstats are linear (default) or logistic regression or cox?
#' @param ewas_DNAm_scale EWAS DNA methylation in beta (default) value or m value?
#' @param ewas_pheno_type Is EWAS phenotype binary or continuous? If binary, will also output estimate of log-odds ratio of clock on phenotype.
#' @param clock either a dir/to/clock/file or data.frame/data.table containing the clock CpGs and corresponding weights
#' @param DNAm_cov_matrix Variance-covariance matrix for DNAm. Must have same covariate adjustments.
#' @param mean_beta dir/to/file containing mean DNAm beta-value at each CpG sites. Needed for M-value to beta-value conversion.
#' a data.frame or data.table with 2 columns: CpGmarker, beta_bar. Best if can match with the EWAS samples.
#' @param N Allows user to manually set EWAS total sample size. Default is to read this in the EWAS summary statistics file.
#' @param cox_pi required only when the ewas_model is "cox" (Cox Proportional Hazard model), we allow user to provide the proportion of cases at the end of study. If not provide, `cox_pi` will be computed as `cox_E/N`
#' @param cox_E cox only; default NULL; the number of cases during the longitudinal study; If `cox_E` and total sample size `N` are supplied, `cox_pi` will be computed as `cox_E/N`
# note, for Cox PH model: must set 
#    ewas_direction = "Y~M"
#    ewas_model = "cox",
#    ewas_DNAm_scale: can be 'beta' or 'm'
#    ewas_pheno_type = "binary"
#    N: total sample size at the end of study (usually the number the controls at the beginning of the study)
#       can leave as NULL if the sumstats contains the total N column
#    cox_pi = cox_E/N, is the proportion of cases at the end of study;
#           if is NULL, must use cox_E/N to compute


f.aging = function(ewas_sumstats, 
                   ewas_direction = "M~Y", 
                   ewas_model = "linear",
                   ewas_DNAm_scale = "beta",
                   ewas_pheno_type = "continuous",
                   clock, 
                   DNAm_cov_matrix, 
                   mean_beta = NULL,
                   varY = NULL, 
                   scale = FALSE, age_adult = 20, 
                   N = NULL,
                   #model = "linear",
                   cox_pi = NULL,
                   cox_E = NULL,
                   verbose = TRUE){
  
  # Record start time
  start_time = Sys.time()
  
  if(verbose){
    cat("\n\n************************************************\n")
    cat("* SARDINE\n")
    cat("* Version 1.2.0\n")
    cat("* (C) 2026 Yuchang Wu | ywu423@wisc.edu\n")
    cat("* University of Wisconsin-Madison\n")
    cat("* GNU General Public License v3\n")
    cat("************************************************\n\n")
    
    cat("[INFO] Start time:", format(start_time, "%Y-%m-%d %H:%M:%S"), "\n")
  }
  
  
  # ========== check the validity of the flags ==========
  
  # ewas_direction: must be either Y~M or M~Y:
  stopifnot(toupper(gsub(" ", "", ewas_direction)) %in% c("Y~M", "M~Y"))
  ewas_direction = toupper(gsub(" ", "", ewas_direction))
  
  # ewas_model must be either linear or logistic:
  stopifnot(ewas_model %in% c("linear", "logistic", "cox"))
  
  # ewas_DNAm_scale must be either beta or m
  stopifnot(tolower(ewas_DNAm_scale) %in% c("beta", "beta-value", "m-value", "mvalue", "m"))
  
  # ewas_pheno_type must be either continuous or binary
  stopifnot(ewas_pheno_type %in% c("continuous", "binary"))
  
  # is user provided EWAS total sample size N, it must be numeric
  if(!is.null(N)){stopifnot(is.numeric(N) | is.integer(N))}
  
  varY_type = ifelse(is.numeric(varY), "User_provided", NA)
  
  
  # =========== read EWAS sumstats as a data.table ===========
  
  if(is.character(ewas_sumstats) && file.exists(ewas_sumstats)){
    ss = tryCatch(
      fread(ewas_sumstats),
      error = function(e) {
        stop("Failed to read ewas_sumstats as a data.table: ", conditionMessage(e), call. = F)
        }
      )
  }else if(is.data.frame(ewas_sumstats) || is.data.table(ewas_sumstats)){
    ss = as.data.table(ewas_sumstats)  # need to convert to data.table
  }else{
    stop("\nERROR! ewas_sumstats must be either EWAS sumstats file or a data.frame or data.table")
  }
  
  # Print summary information if verbose
  if (verbose) {
    cat("[INFO] EWAS sumstats:", nrow(ss), "rows and", ncol(ss), "columns\n")
    cat("[INFO] EWAS sumstats header:\n")
    print(head(ss, 6), digits = 3) 
  }
  
  # Standardize EWAS Column Names (Case-insensitive matching)
  standardize_column_name(ss, "^(cpgmarker|cpg|marker|probe|id)$", "CpGmarker", "CpG")
  standardize_column_name(ss, "^(beta|coef|coeff|estimate)$", "BETA", "BETA")
  standardize_column_name(ss, "^(se|std\\.error|std_error|standard_error)$", "SE", "SE")
  
  if(is.null(N)){
    standardize_column_name(ss, "^(n|samplesize|n_samples)$", "N", "N")
  }else{ # if user provided N, will use users instead
    N_user_supplied = N
    ss[, N := N_user_supplied]
  }
  
  

  # =========== read clock weights as a data.table ===========
  
  if(is.character(clock) && file.exists(clock)){
    ck = tryCatch(
      fread(clock),
      error = function(e) {
        stop("Failed to read clock as a data.table: ", conditionMessage(e), call. = F)
        }
      )
  }else if(is.data.frame(clock) || is.data.table(clock)){
    ck = as.data.table(clock)
  }else{
    stop("\nERROR! clock must be either epigenetic clock weight file or a data.frame or data.table")
  }
  
  # Print clock information if verbose
  if (verbose) {
    cat("[INFO] Clock file:", nrow(ck), "rows and", ncol(ck), "columns\n")
    cat("[INFO] Clock file header:\n")
    print(head(ck, 6)) 
  }
  
  # clock should have CpG column (CpGmarker) and weight column (CoefficientTraining):  
  standardize_column_name(ck, "^(cpgmarker|cpg|marker|probe|id|pred\\.var)$", "CpGmarker", "CpG")
  standardize_column_name(ck, "^(CoefficientTraining|beta|coef|coeff|estimate|weight)$", "CoefficientTraining", "clock weight")
  
  
  # =========== read DNAm covariance matrix ===========
  
  if(is.character(DNAm_cov_matrix) && file.exists(DNAm_cov_matrix)){
    # if exist, must be a .rds file for the DNAm covariance matrix
    Sigma = tryCatch(
      readRDS(DNAm_cov_matrix),
      error = function(e) {
        stop("Failed to read DNAm_cov_matrix as .rds: ", conditionMessage(e), call. = F)
        }
      )

    if (!is.matrix(Sigma)) {stop("\nERROR! DNAm_cov_matrix .rds must contain a matrix")}

  }else if(is.matrix(DNAm_cov_matrix)){
    Sigma = DNAm_cov_matrix
  }else{
    stop("\nERROR! DNAm_cov_matrix must be either a file path or a matrix")
  }
  
  # Covariance matrix must have row names
  if(is.null(rownames(Sigma)) || is.null(colnames(Sigma))){
    stop("\nERROR! DNAm_cov_matrix must have valid rownames and colnames corresponding to CpG sites.")
  }

  # rownames and colnames must be IDENTICAL, not merely both present.
  # diag() only propagates names when they match; otherwise diag(Sigma)[cpg_names]
  # silently returns all-NA and the logistic Y~M conversion fails without an error.
  if(!identical(rownames(Sigma), colnames(Sigma))){
    stop("\nERROR! DNAm_cov_matrix rownames and colnames must be identical and in the same order.")
  }
  
  # Print DNAm information if verbose
  if (verbose) {
    cat("[INFO] DNAm covariance matrix:", nrow(Sigma), "rows and", ncol(Sigma), "columns\n")
    cat("[INFO] DNAm covariance matrix header:\n")
    print(Sigma[c(1:min(6, nrow(Sigma))), c(1:min(6, ncol(Sigma)))], digits = 3) 
  }
  
  
  # =========== M-value to Beta-value transformation ===========
  #
  #  works for linear, logistic, or cox model
  
  if (tolower(ewas_DNAm_scale) %in% c("m-value", "mvalue", "m")) {
    if (is.null(mean_beta)) {
      stop("Error: 'mean_beta' must be provided when 'ewas_DNAm_scale' is 'm-value'.")
    }else{
      # Handle if user provides a file path or a pre-loaded data frame
      if (is.character(mean_beta)) {
        mean_beta_df = data.table::fread(mean_beta)  # Expects columns: CpGmarker, beta_bar
      } else {
        mean_beta_df = mean_beta
      }
    }
    
    if (verbose) message("Converting M-value EWAS summary statistics to Beta-value scale. It's highly recommended that you remove the CpGs with either very high (e.g. >90%) or very low (e.g. < 10%) methylations...")
    
    match_idx = match(ss$CpGmarker, mean_beta_df$CpGmarker)
    missing_idx = is.na(match_idx)
    
    if (any(missing_idx)) {
      if (verbose) message(sprintf("[INFO] Dropping %d CpGs missing from the mean_beta reference.", sum(missing_idx)))
      ss = ss[!missing_idx, ]
      match_idx = match_idx[!missing_idx]
    }
    
    beta_bar = mean_beta_df$beta_bar[match_idx]
    cc = log(2) * beta_bar * (1 - beta_bar)
    
    # Transform effect sizes and standard errors from M- to Beta-value scale
    if (ewas_direction == "M~Y") {
      ss$BETA = cc * ss$BETA
      ss$SE   = cc * ss$SE
    }else{  # Y ~ M:
      ss$BETA = ss$BETA / cc
      ss$SE   = ss$SE / cc
    }
  }
  
  # =========== convert Cox to linear Y ~ M: ==============
  #
  #  note, if ewas_model == "cox", ewas_direction must be "Y~M"
  
  # initialization:
  kappa_cox = NA_real_  # scaling factor for Cox model
  
  if(ewas_model == "cox"){
    if(verbose) cat("[INFO] Converting Cox regression summary statistics to Y ~ M model...\n")
    if(is.null(cox_pi) & is.null(cox_E)){
      stop("[ERROR]: EWAS is Cox PH model, must provide either `cox_pi` or 'cox_E' and 'N'.")
    }
    
    if(ewas_direction != "Y~M"){
      stop("[ERROR]: EWAS is Cox PH model, EWAS direction should be 'Y~M'!")
    }
    
    if(ewas_pheno_type != "binary"){
      stop("[ERROR]: EWAS is Cox PH model, EWAS phenotype must be binary!")
    }
    
    if(is.null(cox_pi)){
      # prioritiz the one provided by the user; if user did not provid, will use cox_E/N to compute
      #
      # note, using SE in the sumstats and Var(M) can back up cox_E
      #   but I haven't implemented it ... most study should report pi?
      cox_N = median(ss$N, na.rm = T)
      cox_pi = cox_E/cox_N
    }
    
    # kappa(pi) is undefined at pi = 0 (division by zero downstream) and at pi >= 1
    if(!is.finite(cox_pi) || cox_pi <= 0 || cox_pi >= 1){
      stop(sprintf("[ERROR] cox_pi must be strictly between 0 and 1; got %.6g.", cox_pi))
    }
    
    kappa_cox = (1 - cox_pi) * (-log(1 - cox_pi))
    
    if(verbose) cat(sprintf("[INFO] Cox conversion factor kappa(pi) = %.6g (pi = %.6g)\n",
                            kappa_cox, cox_pi))
    
    # convert:
    ss$BETA = ss$BETA * kappa_cox
    ss$SE = ss$SE * kappa_cox
    
  }
  
  
  # ===========  convert Y ~ M to M ~ Y: =========== 
  #
  # note, this should before estimating varY and after M-value -> Beta-value conversion
  
  if(ewas_direction == "Y~M"){
    
    if(ewas_model %in% c("linear", "cox")){ # regardless of whether Y is binary or not
      if(verbose) cat("[INFO] Converting Y ~ M summary statistics to M ~ Y direction...\n")
      
      denom = ss$BETA^2 + (ss$N - 2) * ss$SE^2
      
      # convert as if we ran M ~ Y:
      ss$BETA = ss$BETA / denom
      ss$SE   = ss$SE / denom
      
    }else if(ewas_model == "logistic"){# but if run logistic regression Y ~ M
      if(verbose) cat("[INFO] Converting logistic Y ~ M summary statistics to linear M ~ Y...\n")
      
      # Flag the mathematical assumption for the user
      if(verbose) cat("[NOTE] Logistic to linear conversion assumes small effect sizes of Y on M and equal methylation variance between cases and controls.")
      
      # Var(DNAm) from the reference panel
      # CpG intersection between Ref and EWAS to keep more CpGs for stability
      # Safely map Sigma variances to the EWAS data by CpG name
      sigma_vars = diag(Sigma)[ss$CpGmarker]
      alpha    = ss$BETA
      se_alpha = ss$SE
      
      # Stop if no intersecting CpGs are found
      if(all(is.na(sigma_vars))){
        stop("Error: No intersecting CpGs found between EWAS and the reference Sigma matrix.")
      }
      
      # CpGs matched; estimate varY
      varY_temp = median(1/(ss$N * sigma_vars * (se_alpha^2)), na.rm = TRUE)
      
      if(verbose) cat(sprintf("[INFO] Inferred covariate-adjusted Var(Y) from logistic SEs: %.6f\n", varY_temp))
      
      aa = sqrt(1 + 4 * alpha^2 * sigma_vars * varY_temp)
      a2b = 2 * sigma_vars/(1 + aa)
      ss$BETA = alpha * a2b
      ss$SE = se_alpha * (a2b / aa)  # delta-method
      
      # keep the infered var(Y)
      varY = varY_temp 
      varY_type = "Inferred"
      
    }else{
      stop("Error: 'ewas_model' must be either 'linear' or 'logistic'.")
    }
  }
  
  
  # =========== Var(Y) ===========
  
  # can infer varY from EWAS (M ~ Y) sumstats and DNAm reference matrix only
  #  note, we can already estimate it when input EWAS is logistic regressions
  if(is.null(varY)){
    if(all(c("SE", "N") %in% names(ss))){
      # can infer Var(Y) using EWAS sumstats (SE and N columns) and DNAm covariance matrix 
      #  using overlapping CpGs between sumstats and DNAm matrix
      cpg_varY = intersect(ss$CpGmarker, rownames(Sigma))
      ss_varY = ss[match(cpg_varY, ss$CpGmarker), ]
      Sigma_varY = Sigma[cpg_varY, cpg_varY]

      beta = ss_varY$BETA
      se = ss_varY$SE
      n = ss_varY$N
      sigma_diag = diag(Sigma_varY)
      
      varY_hat = sigma_diag/(beta^2 + (n - 2) * se^2)
      varY = median(varY_hat, na.rm = T)  # use median to improve stability
      
      # so here, it may numerically more stable if using a large number of CpGs in Sigma
      # if the Sigma is only for one clock, it may have too few CpGs 
      # leading to unstable/inaccurate/negative VarY
      
      varY_type = "Inferred"
    }
    else{
      stop("[ERROR] Var(Y) is not provided and cannot be inferred from EWAS sumstats either!\n")
    }
  }

  # get the correct number of CpGs included in clock:
  n_intercept_rows = sum(grepl("^\\(?intercept\\)?$", ck$CpGmarker, ignore.case = TRUE))
  if(n_intercept_rows > 0){
    if(verbose) cat("[INFO] Removed", n_intercept_rows, "intercept row(s) from the clock weights.\n")
    ck = ck[!grepl("^\\(?intercept\\)?$", ck$CpGmarker, ignore.case = TRUE), ]
  }
  
  # take CpG intersection between EWAS sumstats, clock weights, and Sigma
  N_CpGs_ewas = length(ss$CpGmarker)
  N_CpGs_clock = length(ck$CpGmarker)
  N_CpGs_ref = dim(Sigma)[1]
  
  N_CpGs_ewas_clock = length(intersect(ss$CpGmarker, ck$CpGmarker))
  N_CpGs_ewas_ref = length(intersect(ss$CpGmarker, rownames(Sigma)))
  N_CpGs_clock_ref = length(intersect(ck$CpGmarker, rownames(Sigma)))
  
  cpg_intersect = intersect(intersect(ss$CpGmarker, ck$CpGmarker), rownames(Sigma))
  
  N_CpGs_final = length(cpg_intersect)
  
  ClockCpG_coverage_ratio = N_CpGs_final / N_CpGs_clock
  ClockCpG_missed_ratio = 1 - N_CpGs_final/N_CpGs_clock
  
  if(verbose){
    cat("[INFO] Number of CpGs in EWAS:", format(N_CpGs_ewas, big.mark = ","), "\n")
    cat("[INFO] Number of CpGs for clock:", format(N_CpGs_clock, big.mark = ","), "\n")
    cat("[INFO] Number of CpGs in DNAm matrix:", format(N_CpGs_ref, big.mark = ","), "\n")
    
    cat("[INFO] Number of overlapping CpGs between EWAS and clock:", format(N_CpGs_ewas_clock, big.mark = ","), "\n")
    cat("[INFO] Number of overlapping CpGs between EWAS and DNAm matrix:", format(N_CpGs_ewas_ref, big.mark = ","), "\n")
    cat("[INFO] Number of overlapping CpGs between clock and DNAm matrix:", format(N_CpGs_clock_ref, big.mark = ","), "\n")
    
    cat("[INFO] Final number of overlapping CpGs:", format(N_CpGs_final, big.mark = ","), "\n")
    cat("[INFO] Percentage of missed clock CpGs: ", round(ClockCpG_missed_ratio, digits = 2)*100, "%", "\n", sep = "")
  }
  
  if(N_CpGs_final == 0){stop("\nERROR! Zero overlapping CpGs found between EWAS, Clock, and Covariance Matrix.")}
  
  
  # =========== estimate regression coefficients ===========
 
  # align ss, ck, and Sigma to the same CpG order
  ss = ss[match(cpg_intersect, ss$CpGmarker), ]
  ck = ck[match(cpg_intersect, ck$CpGmarker), ]
  Sigma = Sigma[cpg_intersect, cpg_intersect]

  # EWAS sumstats:
  beta_hat = ss$BETA  # EWAS effect size estimates
  
  if(is.null(N)){ # find total sample size in sumstats
    N = median(ss$N, na.rm = T)  # N used for the CpG sites we're interested in
  }else{  # user provide the total sample size
    N = N
  }
  
  # clock weight vector:
  w = ck$CoefficientTraining

  # core matrix:
  w_T_beta = sum(w * beta_hat, na.rm = T)
  w_T_Sigma_w = as.numeric(t(w) %*% Sigma %*% w)

  
  # Var(Clock) must be strictly positive
  if(!is.finite(w_T_Sigma_w) || w_T_Sigma_w <= 0){
    stop(sprintf("\nERROR! Estimated clock variance w'Sigma w = %.6g is not positive. Check that DNAm_cov_matrix is a valid covariance matrix on the same CpG scale as the EWAS.", w_T_Sigma_w))
  }
  

  ## if has negative var_gamma_hat or var_lambda_hat
  #   root cause of this is the estimated varY is too large
  #   solution: if so, shrink it to be a little smaller than max possible value
  
  note_varY_capping = NA
  varY_ratio = NA_real_   # varY / varY_max, always reported as a diagnostic
  
  if (!is.null(varY) && w_T_beta != 0) {
    varY_max <- w_T_Sigma_w / w_T_beta^2  # when clock and Y are perfectly correlated
    varY_ratio = varY / varY_max
    
    if (varY > varY_max) {
      msg_capping = sprintf("Estimated varY (%.6g) exceeds maximum (%.6g), ratio = varY/varY_max = %.6f. Capped to avoid negative error variance. If the ratio is much larger than 1: check that DNAm_cov_matrix matches the EWAS population and covariate adjustment.", varY, varY_max, varY_ratio)
      
      if (verbose) {
        cat(sprintf("[INFO] %s\n", msg_capping))
      }
      varY = varY_max * (1 - 1e-8) # shrink varY to be slightly less than max
      note_varY_capping = msg_capping
    }
  }
  
  
  ## linear regression: Y = C * gamma + e
  
  gamma_hat = w_T_beta * varY / w_T_Sigma_w
  sigma_e2 = varY - gamma_hat^2 * w_T_Sigma_w
  var_gamma_hat = sigma_e2/((N - 2) * w_T_Sigma_w)
  
  note_gamma = ifelse(var_gamma_hat < 0, "var_gamma_hat < 0", NA)
  
  se_gamma_hat = suppressWarnings(sqrt(var_gamma_hat))
  z_gamma_hat = gamma_hat / se_gamma_hat
  p_gamma_hat = 2 * pnorm(abs(z_gamma_hat), lower.tail = FALSE)
  
  
  # logistic regression Y ~ C =======
  #
  # if phenotype is binary and can convert to logistic scale 
  # (logit(P(Y = 1 | C)) = gamma_logistic * C + gamma_0):
  #
  # Initialize variables to avoid errors
  gamma_hat_logistic = NA_real_
  var_gamma_hat_logistic = NA_real_
  se_gamma_hat_logistic = NA_real_
  z_gamma_hat_logistic = NA_real_
  p_gamma_hat_logistic = NA_real_
  note_gamma_logistic = NA_character_
  
  if(ewas_pheno_type == "binary") {
    if(ewas_model == "linear") {
      if(verbose) cat("[NOTE] EWAS is a linear regression with a binary phenotype; we can still estimate log-odds ratio of Clock on this binary phenotype.\n")
    }
    
    if(verbose) cat("[INFO] Converting final gamma estimate from linear to logistic (log-odds) scale...\n")
    
    if (!is.finite(sigma_e2) || sigma_e2 <= 0 ||
        (is.finite(varY_ratio) && varY_ratio > 0.99)) {
      
      note_gamma_logistic = sprintf(
        "sigma_e2 = %.6g (varY_ratio = %.6g); logistic-scale gamma not identifiable.",
        sigma_e2, varY_ratio)
      if (verbose) cat("[WARN]", note_gamma_logistic, "\n")
      
    } else {  # note, this conversion introduces some downward bias under alternative
      gamma_hat_logistic = gamma_hat / sigma_e2
      
      dg = (varY + gamma_hat^2 * w_T_Sigma_w) / sigma_e2^2  
      var_gamma_hat_logistic = var_gamma_hat * dg^2
      se_gamma_hat_logistic = suppressWarnings(sqrt(var_gamma_hat_logistic))
      
      # H0: gamma_logistic = 0 and H0: gamma = 0 are the SAME hypothesis. 
      # The delta-method Wald statistic  is attenuated by 
      # sigma_e2/(varY + gamma^2 * w'Sigma w) and tends to 0
      # as the association tends to perfect (sigma_e2 -> 0). 
      # Use the linear scale statistic instead.
      
      z_gamma_hat_logistic = z_gamma_hat  
      p_gamma_hat_logistic = p_gamma_hat
      
      note_gamma_logistic = paste0("SE by delta-method; Z and P inherited from the linear-scale gamma")
      
      if(is.finite(var_gamma_hat_logistic) && var_gamma_hat_logistic < 0){
        note_gamma_logistic = paste0("var_gamma_hat_logistic < 0; ", note_gamma_logistic)
      }
      
    }
    if(verbose) cat(paste0("[NOTE] For gamma in a logistic regression of Phenotype ~ Clock: ",
                           "Z and P are inherited from the linear-scale gamma; ",
                           "SE computed from delta-method, so Z != coef/SE for this row. ",
                           "CI should be obtained by transforming the 'linear' gamma CI endpoints.\n"))
    
  }

  
  ## regression: C = Y * lambda + u
  
  lambda_hat = w_T_beta
  sigma_u2 = w_T_Sigma_w - lambda_hat^2 * varY
  var_lambda_hat = sigma_u2/((N - 2) * varY)
  
  note_lambda = ifelse(var_lambda_hat < 0, "var_lambda_hat < 0", NA)
  
  se_lambda_hat = suppressWarnings(sqrt(var_lambda_hat))
  z_lambda_hat = lambda_hat / se_lambda_hat
  p_lambda_hat = 2 * pnorm(abs(z_lambda_hat), lower.tail = FALSE)
  

  # if Horvath1, Horvath2, or PedBE clock, need to 
  # re-scale effect size and SE
  if(is.character(clock)){  # if can tell the clock is Horvath or PEDBE
    ck_nm = gsub("\\.csv", "", basename(clock))
    
    if(!scale & grepl("horvath|pedbe", tolower(ck_nm))){
      cat("[Note] Seems your clock weight filename contains '", ck_nm, 
          "' yet you did not scale the point estimate and SE (Z-score and P-value are not affected). Make sure this is intentional.\n", sep="")
    }
    
  }
   
  if(scale & (is.numeric(age_adult) | is.integer(age_adult))){
    scale_factor = 1 + age_adult
    
    gamma_hat = gamma_hat / scale_factor
    var_gamma_hat = var_gamma_hat / (scale_factor^2)
    se_gamma_hat = se_gamma_hat / scale_factor
    
    gamma_hat_logistic = gamma_hat_logistic / scale_factor
    var_gamma_hat_logistic = var_gamma_hat_logistic / (scale_factor^2)
    se_gamma_hat_logistic = se_gamma_hat_logistic / scale_factor
    
    lambda_hat = lambda_hat * scale_factor
    var_lambda_hat = var_lambda_hat * (scale_factor^2)
    se_lambda_hat = se_lambda_hat * scale_factor
    
    if (verbose) {
      cat("[INFO] Applied clock scaling factor =", scale_factor, "to point estimate and SE", "\n")
    }
  }
  
  
  # clock-level log hazard ratio (Cox input only): ============
  gamma_hat_cox    = NA_real_
  se_gamma_hat_cox = NA_real_
  note_gamma_cox   = NA_character_
  
  if(ewas_model == "cox"){
    gamma_hat_cox    = gamma_hat    / kappa_cox
    se_gamma_hat_cox = se_gamma_hat / kappa_cox
    
    if(identical(varY_type, "Inferred")){
      note_gamma_cox = "Clock-level log hazard ratio; invariant to cox_pi. First-order approximation, attenuated toward the null."
    }else{
      note_gamma_cox = "Clock-level log hazard ratio. WARNING: varY was user-supplied, so the cox_pi cancellation does NOT hold; this estimate requires varY = pi(1-pi) with a correct cox_pi."
      if(verbose) cat("[WARN] ewas_model = 'cox' with user-supplied varY: cox_pi no longer cancels. Setting varY = NULL is recommended.\n")
    }
    
    if(verbose){
      cat(sprintf("[INFO] Clock-level log HR = %.6g (SE %.6g); HR = %.4f\n",
                  gamma_hat_cox, se_gamma_hat_cox, exp(gamma_hat_cox)))
    }
  }
  
  # main results:
  res = rbind(
    data.frame(model = "Clock ~ Phenotype", 
               param = "lambda",
               coef = lambda_hat, 
               SE = se_lambda_hat,
               Z = z_lambda_hat,
               P = p_lambda_hat, 
               N = N,
               Note = note_lambda),
    data.frame(model = "Phenotype ~ Clock", 
               param = "gamma",
               coef = gamma_hat, 
               SE = se_gamma_hat,
               Z = z_gamma_hat,
               P = p_gamma_hat, 
               N = N,
               Note = note_gamma)
    )
  
  if(ewas_pheno_type == "binary" & !is.na(gamma_hat_logistic)){
    # if has gamma_logistic
    
    res = rbind(res,
                data.frame(
                  model = "Phenotype ~ Clock (logistic)", 
                  param = "gamma_logistic",
                  coef = gamma_hat_logistic, 
                  SE = se_gamma_hat_logistic,
                  Z = z_gamma_hat_logistic,
                  P = p_gamma_hat_logistic, 
                  N = N,
                  Note = note_gamma_logistic
                ))
    
  }
  
  # if input is cox sumstats, we estimate log HR:
  # note, Cox pathway's robustness ordering is log_HR > gamma ≫ lambda, gamma_logistic
  #   they all have a same Z score and p-value though
  if(ewas_model == "cox" & !is.na(gamma_hat_cox)){
    res = rbind(res,
                data.frame(
                  model = "Phenotype ~ Clock (Cox)",
                  param = "log_HR",
                  coef  = gamma_hat_cox,
                  SE    = se_gamma_hat_cox,
                  Z     = z_gamma_hat,  # z and pvalue from linear-scale gamma
                  P     = p_gamma_hat,
                  N     = N,
                  Note  = note_gamma_cox
                ))
    
    # HR and 95% CI: defined only for the log_HR row.
    is_hr  = res$param == "log_HR"
    z_crit = qnorm(0.975)
    res$HR     = ifelse(is_hr, exp(res$coef), NA_real_)
    res$HR_L95 = ifelse(is_hr, exp(res$coef - z_crit * res$SE), NA_real_)
    res$HR_U95 = ifelse(is_hr, exp(res$coef + z_crit * res$SE), NA_real_)
    
    idx = res$param == "gamma_logistic"
    if(any(idx)){
      res$Note[idx] = append_note(res$Note[idx],
                                  "Log-odds ratio for cumulative incidence at end of follow-up; param = 'log_HR' is log of hazard ratio.")
    }
  }
  
  # add -log10P to avoid pvalue = 0:
  # we're expecting same Z and Pvalues for all 2/3 rows in the output
  res$log10_P = -(pnorm(abs(res$Z), lower.tail = FALSE, log.p = TRUE) + log(2))/log(10)
  
  # add number of CpGs to final output:
  res$N_CpGs_ewas = N_CpGs_ewas
  res$N_CpGs_clock = N_CpGs_clock
  res$N_CpGs_DNAmRef = N_CpGs_ref
  res$N_CpGs_final = N_CpGs_final
  res$ClockCpG_coverage_ratio = ClockCpG_coverage_ratio
  
  res$varY = varY
  res$varY_type = varY_type
  res$varY_ratio = varY_ratio              # varY / varY_max; >1 means capping fired
  
  res$Note = append_note(res$Note, note_varY_capping)
  
  # put other useful info:
  res_metadata = data.frame(
    N_CpGs_ewas = N_CpGs_ewas,
    N_CpGs_clock = N_CpGs_clock,
    N_CpGs_DNAmRef = N_CpGs_ref,
    
    N_CpGs_ewas_clock = N_CpGs_ewas_clock,  # overlap between EWAS and clock
    N_CpGs_ewas_ref = N_CpGs_ewas_ref,      # overlap between EWAS and Sigma
    N_CpGs_clock_ref = N_CpGs_clock_ref,    # overlap between clock and Sigma
    
    N_CpGs_final = N_CpGs_final,
    ClockCpG_coverage_ratio = ClockCpG_coverage_ratio,
    ClockCpG_missed_ratio = ClockCpG_missed_ratio,  # ratio of Clock CpGs not used
    
    cox_pi     = ifelse(ewas_model == "cox", cox_pi,    NA_real_),
    kappa_cox  = ifelse(ewas_model == "cox", kappa_cox, NA_real_),  # for cox
    
    varY_ratio = varY_ratio  # varY/varY_max
  )
  

  # Record end time and calculate elapsed time
  end_time = Sys.time()
  elapsed_seconds = as.numeric(difftime(end_time, start_time, units = "secs"))
  
  # Format elapsed time: hours, minutes, seconds
  hours = floor(elapsed_seconds / 3600)
  minutes = floor((elapsed_seconds - hours * 3600) / 60)
  seconds = round(elapsed_seconds - hours * 3600 - minutes * 60)
  
  # Build time string
  time_parts = c()
  if (hours > 0) time_parts = c(time_parts, paste0(hours, " hour", ifelse(hours == 1, "", "s")))
  if (minutes > 0) time_parts = c(time_parts, paste0(minutes, " minute", ifelse(minutes == 1, "", "s")))
  if (seconds > 0 || length(time_parts) == 0) {
    time_parts = c(time_parts, paste0(seconds, " second", ifelse(seconds == 1, "", "s")))
  }
  elapsed_time_str = paste(time_parts, collapse = ", ")
  
  if(verbose){
    cat("\n")
    cat("[INFO] Analysis complete!\n")
    cat("[INFO] End time:", format(end_time, "%Y-%m-%d %H:%M:%S"), "\n")
    cat("[INFO] Elapsed time:", elapsed_time_str, "\n")
    cat("\n")
  }

  #attach the metadata as the attribute
  attr(res, "metadata") = res_metadata
  
  return(res)
}

# test:
# 
# t1 = as.data.frame(fread("20260403c/sumstats/ewas_sumstats_Adj_sex_age_age2_Age_std.txt"))
# t2 = as.data.frame(fread("data/clocks/Horvath1.csv"))
# t3 = readRDS("20260403c/sigma/Sigma_Insample_Adj_sex_age_age2_Age_std_N=8880.rds")
# 
# 
# f.aging(ewas_sumstats = "20260504/sumstats/ewas_sumstats_Adj_sex_age_age2_Age_std.txt",
#         clock = "data/clocks/Horvath1.csv",
#         DNAm_cov_matrix = "20260504/sigma/Sigma_Insample_Adj_sex_age_age2_Age_std_N=15677.rds",
#         verbose = TRUE)
# #               model  param      coef          SE        Z P     N Note  log10_P 
# # 1 Clock ~ Phenotype lambda 0.4947275 0.002008788 246.2816 0 15677   NA 13173.48           
# # 2 Phenotype ~ Clock  gamma 1.6062190 0.006521879 246.2816 0 15677   NA 13173.48
# 
# f.aging(ewas_sumstats = "20260504/sumstats/ewas_sumstats_Adj_sex_age_age2_Herbal_supplements.txt",
#         clock = "data/clocks/Horvath1.csv",
#         DNAm_cov_matrix = "20260504/sigma/Sigma_Insample_Adj_sex_age_age2_Herbal_supplements_N=15677.rds",
#         ewas_pheno_type = "binary",
#         verbose = TRUE)
#                          model          param        coef          SE        Z         P     N
# 1            Clock ~ Phenotype         lambda 0.005102679 0.004287306 1.190183 0.2339745 15677
# 2            Phenotype ~ Clock          gamma 0.017708526 0.014878827 1.190183 0.2339745 15677
# 3 Phenotype ~ Clock (logistic) gamma_logistic 0.081433541 0.068433395 1.190183 0.2339745 15677

# # test:
# t1 = f.aging(ewas_sumstats = "20260403c/sumstats/ewas_sumstats_Adj_sex_age_age2_Age_std.txt",
#              clock = "data/clocks/Horvath1.csv",
#              DNAm_cov_matrix = "20260403c/sigma/Sigma_Insample_Adj_sex_age_age2_Age_std_N=8880.rds")
# # note, insample DNAm file names are misleading; when age is phenotype, age is NOT adjusted in insample DNAm matrix
# 
# t2 = f.aging(ewas_sumstats = "20260403c/sumstats/ewas_sumstats_NoCovariates_Age_std.txt",
#              clock = "data/clocks/Horvath1.csv",
#              DNAm_cov_matrix = "20260403c/sigma/Sigma_ref_NoCovariates_EUR_N500.rds")
# 
# t3 = f.aging(ewas_sumstats = "20260403c/sumstats/ewas_sumstats_NoCovariates_Age_std.txt",
#              clock = "data/clocks/Horvath1.csv",
#              DNAm_cov_matrix = "20260403c/sigma/Sigma_ref_Adj_sex_only_EUR_N500.rds")
# 
# t4 = f.aging(ewas_sumstats = "20260403c/sumstats/ewas_sumstats_Adj_sex_age_age2_Age_std.txt",
#              clock = "data/clocks/Horvath1.csv",
#              DNAm_cov_matrix = "20260403c/sigma/Sigma_ref_Adj_sex_age_age2_EUR_N500.rds")
# 
# t5 = f.aging(ewas_sumstats = "20260403c/sumstats/ewas_sumstats_Adj_sex_age_age2_19CellTypes_Cigarettes_per_day_std.txt",
#              clock = "data/clocks/Horvath1.csv",
#              DNAm_cov_matrix = "20260403c/sigma/Sigma_ref_Adj_sex_age_age2_19CellTypes_EUR_N500.rds")


#               model  param      coef          SE P    N Note N_CpGs  varY varY_type
# 1 Phenotype ~ Clock  gamma 1.6074604 0.008896106 0 8880   NA    335  0.9688041  Inferred
# 2 Clock ~ Phenotype lambda 0.4890806 0.002706700 0 8880   NA    335  0.9688041  Inferred
# 
# Extract the hidden metadata data.frame
# meta_info = attr(tt, "metadata")

# 
# f.aging(ewas = "ewas_sumstats_BMI_std.txt", 
#         clock = "../clocks/Horvath1.csv",
#         DNAm_cov_matrix = "sigma/Sigma_BMI_std_N=9980.rds",
#         N = 9980)
# 




















