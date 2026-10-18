# SARDINE test fixture generator
#
# Builds a small synthetic dataset, then runs it through sardine.R to produce
# outputs. This was used to produce expected outputs from the original
# R code prior to restructuring. This should not be re-run unless the expected
# outputs to the test cases should change due to bug fixes, etc.

set.seed(20260917)
options(warn = 1) # print warnings immediately, so they land inside the right case's sink

suppressMessages(library(data.table))

this_script_path   = sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
root          = normalizePath(file.path(dirname(this_script_path), ".."))
fixtures_dir  = file.path(root, "fixtures")
inputs_dir    = file.path(fixtures_dir, "inputs")
expected_dir  = file.path(fixtures_dir, "expected")
dir.create(inputs_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(expected_dir, recursive = TRUE, showWarnings = FALSE)

package_files = c("utils.R", "io.R", "conversions.R", "estimate.R", "sardine.R", "print.R")
for (f in package_files) {
  script_path_i = file.path(root, "R", "src", f)
  stopifnot(file.exists(script_path_i))
  invisible(capture.output(source(script_path_i)))
}
stopifnot(exists("sardine"))

cat(sprintf("[generate_fixtures] R version: %s\n", R.version.string))
cat(sprintf("[generate_fixtures] data.table version: %s\n", as.character(packageVersion("data.table"))))

# ============================================================
# Synthetic CpGs: cg00000001 .. cg00000400, partitioned by which of
# {clock, DNAm reference, baseline EWAS sumstats} they fall in:
#
#   1-320    clock + DNAm reference + EWAS   (the final analysis set)
#   321-380  DNAm reference + EWAS, not clock
#   381-390  clock + EWAS, not DNAm reference
#   391-395  EWAS only
#   396-400  none of the three; used only in the separate
#            zero-overlap-error fixture (see zero_overlap_cpgs below)
# ============================================================

n_cpgs = 400
cpg_ids = sprintf("cg%08d", seq_len(n_cpgs))

sigma_cpgs = cpg_ids[1:380]              # DNAm reference covers 1-380
ewas_cpgs  = cpg_ids[1:395]              # baseline EWAS sumstats cover 1-395
clock_cpgs = cpg_ids[c(1:320, 381:390)]  # clock covers 1-320 and 381-390 (330 total)

zero_overlap_cpgs = cpg_ids[391:400]     # disjoint from clock_cpgs and sigma_cpgs

# ============================================================
# DNAm covariance matrix (random positive definite)
# ============================================================

A = matrix(rnorm(380 * 380), nrow = 380)        # standard normal noise
Sigma = crossprod(A) / 380 + diag(380) * 0.05   # A'A is PSD; add to diagonal for PD
dimnames(Sigma) = list(sigma_cpgs, sigma_cpgs) 
stopifnot(all(eigen(Sigma, only.values = TRUE)$values > 0)) # check that it is PD

# ============================================================
# Clock weights: 330 CpGs + an intercept row
# ============================================================

clock_weight_vals = rnorm(330, mean = 0, sd = 0.05)  # draw weights from normal distribution
clock_intercept = runif(1, min = 20, max = 60)       # arbitrary intercept, adult-age-ish range

clock_weight_dt = data.table(
  CpGmarker  = c("(Intercept)", clock_cpgs),
  coef = c(clock_intercept, clock_weight_vals)
)

# ============================================================
# EWAS summary statistics
# ============================================================

sample_size   = 5000
beta_vals = rnorm(395, mean = 0, sd = 0.01)
se_vals   = runif(395, min = 0.005, max = 0.02)
ewas_ss_dt = data.table(
  CpGmarker = ewas_cpgs,
  BETA      = beta_vals,
  SE        = se_vals,
  N         = sample_size
)

# mean beta-value reference, for the M-value -> beta-value conversion case
mean_beta_vals = runif(395, min = 0.2, max = 0.8)
mean_beta_dt = data.table(CpGmarker = ewas_cpgs, beta_bar = mean_beta_vals)

# sumstats with zero CpG overlap with clock or Sigma, for the error case
ewas_zero_overlap_ss_dt = data.table(
  CpGmarker = zero_overlap_cpgs,
  BETA      = rnorm(10, 0, 0.01),
  SE        = runif(10, 0.005, 0.02),
  N         = sample_size
)

# ============================================================
# Create varY below and over the varY_max = w'Sigma w / (w'beta)^2
# threshold. Used to check varY capping path.
# ============================================================

# subset to CpGs in all inputs
cpg_final = Reduce(intersect, list(ewas_ss_dt$CpGmarker, clock_cpgs, sigma_cpgs))
stopifnot(length(cpg_final) > 0)
ewas_ss_final         = ewas_ss_dt[match(cpg_final, ewas_ss_dt$CpGmarker)]
clock_weights_final   = clock_weight_dt$coef[match(cpg_final, clock_weight_dt$CpGmarker)]
Sigma_final           = Sigma[cpg_final, cpg_final]

# generate varY that are under and above the capping threshold
w_T_beta           = sum(clock_weights_final * ewas_ss_final$BETA, na.rm = TRUE)
w_T_Sigma_w        = as.numeric(t(clock_weights_final) %*% Sigma_final %*% clock_weights_final)
varY_threshold     = w_T_Sigma_w / w_T_beta^2
varY_user_uncapped = varY_threshold * 0.3
varY_user_capped   = varY_threshold * 5

cat(sprintf("[generate_fixtures] Num overlapped cpgs (baseline) = %d\n", length(cpg_final)))
cat(sprintf("[generate_fixtures] varY threshold (baseline)      = %.6g\n", varY_threshold))
cat(sprintf("[generate_fixtures] varY (uncapped case)           = %.6g\n", varY_user_uncapped))
cat(sprintf("[generate_fixtures] varY (capped case)             = %.6g\n", varY_user_capped))

# ============================================================
# Write input fixtures
# ============================================================

fwrite(ewas_ss_dt,      file.path(inputs_dir, "ewas_sumstats_baseline_synthetic.tsv"), sep = "\t")
fwrite(clock_weight_dt,   file.path(inputs_dir, "clock_weights_baseline_synthetic.tsv"), sep = "\t")
sigma_text = matrix(sprintf("%.17g", Sigma), nrow = nrow(Sigma), dimnames = dimnames(Sigma))
fwrite(data.table(CpGmarker = rownames(Sigma), sigma_text),
       file.path(inputs_dir, "sigma_baseline_synthetic.tsv"), sep = "\t", quote = FALSE)
fwrite(mean_beta_dt, file.path(inputs_dir, "mean_beta_ref_baseline_synthetic.tsv"), sep = "\t")
fwrite(ewas_zero_overlap_ss_dt, file.path(inputs_dir, "ewas_sumstats_zero_overlap_synthetic.tsv"), sep = "\t")

# Save clock weights with a Horvath-flavored filename since the original R script looks
# at the filename to warn about un-scaled Horvath/PedBE clocks.
fwrite(clock_weight_dt,   file.path(inputs_dir, "clock_baseline_horvath_named.tsv"), sep = "\t")

# ============================================================
# Test manifest: one row per test case, with the exact sardine() arguments.
# File paths are relative to fixtures/.
# ============================================================

manifest = data.table(
  case            = character(),
  description     = character(),
  expect_error    = logical(),
  ewas_sumstats   = character(),
  clock           = character(),
  DNAm_cov_matrix = character(),
  mean_beta       = character(),
  ewas_direction  = character(),
  ewas_model      = character(),
  ewas_DNAm_scale = character(),
  ewas_pheno_type = character(),
  varY            = numeric(),
  scale           = logical(),
  age_adult       = numeric(),
  N               = numeric(),
  cox_pi          = numeric(),
  cox_E           = numeric()
)

add_case = function(case, description, expect_error = FALSE,
                     ewas_sumstats = "inputs/ewas_sumstats_baseline_synthetic.tsv",
                     clock = "inputs/clock_weights_baseline_synthetic.tsv",
                     DNAm_cov_matrix = "inputs/sigma_baseline_synthetic.tsv",
                     mean_beta = NA_character_,
                     ewas_direction = "M~Y", ewas_model = "linear",
                     ewas_DNAm_scale = "beta", ewas_pheno_type = "continuous",
                     varY = NA_real_, scale = FALSE, age_adult = 20,
                     N = NA_real_, cox_pi = NA_real_, cox_E = NA_real_) {
  manifest <<- rbind(manifest, data.table(
    case = case, description = description, expect_error = expect_error,
    ewas_sumstats = ewas_sumstats, clock = clock, DNAm_cov_matrix = DNAm_cov_matrix,
    mean_beta = mean_beta, ewas_direction = ewas_direction, ewas_model = ewas_model,
    ewas_DNAm_scale = ewas_DNAm_scale, ewas_pheno_type = ewas_pheno_type,
    varY = varY, scale = scale, age_adult = age_adult, N = N,
    cox_pi = cox_pi, cox_E = cox_E
  ))
}

add_case("baseline_linear",
         "M~Y, linear, beta, continuous; varY inferred from EWAS + Sigma")

add_case("varY_user_supplied",
         "Same data as baseline_linear but with a user-supplied varY below varY capping threshold",
         varY = varY_user_uncapped)

add_case("reverse_direction",
         "Y~M, linear, beta, continuous; exercises the Y~M -> M~Y conversion",
         ewas_direction = "Y~M")

add_case("logistic_ewas",
         "Y~M, logistic, beta, binary; exercises the logistic->linear conversion and inferred varY",
         ewas_direction = "Y~M", ewas_model = "logistic", ewas_pheno_type = "binary")

add_case("cox_pi",
         "Y~M, cox, beta, binary; cox_pi supplied)",
         ewas_direction = "Y~M", ewas_model = "cox", ewas_pheno_type = "binary",
         cox_pi = 0.08)

add_case("cox_event_count",
         "Y~M, cox, beta, binary; cox_E supplied)",
         ewas_direction = "Y~M", ewas_model = "cox", ewas_pheno_type = "binary",
         cox_E = 400, N = 5000)

add_case("m_value_input",
         "M~Y, linear, m-value scale, continuous; mean_beta supplied for M->beta conversion",
         ewas_DNAm_scale = "m", mean_beta = "inputs/mean_beta_ref_baseline_synthetic.tsv")

add_case("horvath_scaling",
         "M~Y, linear, beta, continuous; scale=TRUE, age_adult=20, Horvath-named clock file",
         clock = "inputs/clock_baseline_horvath_named.tsv", scale = TRUE, age_adult = 20)

add_case("varY_capping_fires",
         "Same data as baseline_linear but varY is user-supplied above the capping threshold, so it gets capped",
         varY = varY_user_capped)

add_case("zero_overlap_error",
         "EWAS CpGs are entirely disjoint from clock and Sigma CpGs; raises the zero-overlap error",
         expect_error = TRUE,
         ewas_sumstats = "inputs/ewas_sumstats_zero_overlap_synthetic.tsv")

fwrite(manifest, file.path(fixtures_dir, "test_manifest.tsv"), sep = "\t")
cat(sprintf("[generate_fixtures] wrote test_manifest.tsv with %d cases\n", nrow(manifest)))

# ============================================================
# Run each case through the R script and save results as expected outputs
# ============================================================

build_args = function(row) {
  args = list(
    ewas_sumstats   = file.path(fixtures_dir, row$ewas_sumstats),
    clock           = file.path(fixtures_dir, row$clock),
    DNAm_cov_matrix = file.path(fixtures_dir, row$DNAm_cov_matrix),
    ewas_direction  = row$ewas_direction,
    ewas_model      = row$ewas_model,
    ewas_DNAm_scale = row$ewas_DNAm_scale,
    ewas_pheno_type = row$ewas_pheno_type,
    scale           = row$scale,
    age_adult       = row$age_adult,
    verbose         = TRUE
  )
  if (!is.na(row$mean_beta)) args$mean_beta = file.path(fixtures_dir, row$mean_beta)
  if (!is.na(row$varY))      args$varY      = row$varY
  if (!is.na(row$N))         args$N         = row$N
  if (!is.na(row$cox_pi))    args$cox_pi    = row$cox_pi
  if (!is.na(row$cox_E))     args$cox_E     = row$cox_E
  args
}

# ============================================================
# Runs sardine() with stdout and the message stream (which carries
# log_message()'s INFO/NOTE/WARN output) each captured separately.
# ============================================================
run_case = function(args) {
  stdout_con = textConnection("stdout_lines", open = "w", local = TRUE)
  stderr_con = textConnection("stderr_lines", open = "w", local = TRUE)
  sink(stdout_con, type = "output")
  sink(stderr_con, type = "message")

  out = tryCatch(
    list(ok = TRUE, value = do.call(sardine, args)),
    error = function(e) list(ok = FALSE, message = conditionMessage(e))
  )

  sink(type = "message")
  sink(type = "output")
  close(stdout_con)
  close(stderr_con)

  c(out, list(stdout = stdout_lines, stderr = stderr_lines))
}

# ============================================================
# Blanks out wall-clock-dependent text (start/end time, elapsed time) so the
# saved fixtures and later comparisons are stable across runs.
# ============================================================
normalize_output_lines = function(lines) {
  lines = gsub("(Start time:).*",   "\\1 <TIMESTAMP>", lines)
  lines = gsub("(End time:).*",     "\\1 <TIMESTAMP>", lines)
  lines = gsub("(Elapsed time:).*", "\\1 <ELAPSED>",    lines)
  lines
}

all_ok = TRUE

for (i in seq_len(nrow(manifest))) {
  row = manifest[i]
  args = build_args(row)

  out = run_case(args)

  writeLines(normalize_output_lines(out$stdout), file.path(expected_dir, paste0(row$case, "_stdout.txt")))
  writeLines(normalize_output_lines(out$stderr), file.path(expected_dir, paste0(row$case, "_stderr.txt")))

  if (row$expect_error) {
    if (out$ok) {
      cat(sprintf("[generate_fixtures] FAIL  %-22s expected an error but succeeded\n", row$case))
      all_ok = FALSE
    } else {
      writeLines(out$message, file.path(expected_dir, paste0(row$case, "_error.txt")))
      cat(sprintf("[generate_fixtures] OK    %-22s error captured: %s\n", row$case,
                  gsub("\n", " ", out$message)))
    }
  } else {
    if (!out$ok) {
      cat(sprintf("[generate_fixtures] FAIL  %-22s unexpected error: %s\n", row$case, out$message))
      all_ok = FALSE
    } else {
      res  = out$value
      meta = attr(res, "metadata")
      fwrite(res,  file.path(expected_dir, paste0(row$case, "_results.tsv")), sep = "\t")
      fwrite(meta, file.path(expected_dir, paste0(row$case, "_metadata.tsv")), sep = "\t")
      cat(sprintf("[generate_fixtures] OK    %-22s %d result row(s)\n", row$case, nrow(res)))
    }
  }
}

if (!all_ok) stop("[generate_fixtures] one or more fixture cases did not behave as expected — see FAIL lines above")

cat("[generate_fixtures] done.\n")
