# Runs every case in the shared fixture manifest and checks it against the
# expected outputs: the results and metadata (numbers within REL_TOL), the error
# message for cases that should fail, and the printed stdout/stderr line by line.
#
# The fixtures are copied into inst/extdata/fixtures from the repo-level
# fixtures/ by `pixi run -e dev sync`; edit them there, not here.

fixtures_dir <- system.file("extdata", "fixtures", package = "sardine")
expected_dir <- file.path(fixtures_dir, "expected")

# Tolerance allowed between actual and expected numbers. Loose enough to absorb
# floating-point rounding, which varY_capping_fires amplifies to ~6e-8 because
# capping leaves sigma_e2 = varY - gamma^2 * w'Sigma w as a near-cancellation
REL_TOL <- 1e-6

# The version is saved as this placeholder in the expected stdout and metadata,
# so a version bump doesn't require regenerating the expected outputs
VERSION <- as.character(utils::packageVersion("sardine"))
VERSION_PLACEHOLDER <- "<VERSION>"

read_tsv <- function(path) {
  data.table::fread(path, sep = "\t", na.strings = "")
}

read_expected_lines <- function(case, kind) {
  readLines(file.path(expected_dir, paste0(case, "_", kind, ".txt")), warn = FALSE)
}

build_args <- function(row) {
  args <- list(
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
  if (!is.na(row$mean_beta)) args$mean_beta <- file.path(fixtures_dir, row$mean_beta)
  if (!is.na(row$varY))      args$varY      <- row$varY
  if (!is.na(row$N))         args$N         <- row$N
  if (!is.na(row$cox_pi))    args$cox_pi    <- row$cox_pi
  if (!is.na(row$cox_E))     args$cox_E     <- row$cox_E
  args
}

# Runs sardine() and captures its result (or error), its stdout, and its
# message stream. Messages and warnings are recorded as the console would show
# them on stderr with options(warn = 1), which is how the expected files were made.
run_case <- function(args) {
  value <- NULL
  error <- NULL
  stderr_lines <- character(0)
  record <- function(text) stderr_lines <<- c(stderr_lines, strsplit(sub("\n$", "", text), "\n", fixed = TRUE)[[1]])

  stdout_lines <- utils::capture.output(
    withCallingHandlers(
      tryCatch(value <- do.call(sardine, args), error = function(e) error <<- conditionMessage(e)),
      message = function(m) {
        record(conditionMessage(m))
        invokeRestart("muffleMessage")
      },
      warning = function(w) {
        record(paste0("Warning: ", conditionMessage(w)))
        invokeRestart("muffleWarning")
      }
    )
  )
  list(value = value, error = error, stdout = stdout_lines, stderr = stderr_lines)
}

# Blanks out wall-clock-dependent text (start/end time, elapsed time) and the
# SARDINE version so the comparison is stable across runs and releases.
normalize_output_lines <- function(lines) {
  lines <- gsub(sprintf("* SARDINE v%s", VERSION), sprintf("* SARDINE v%s", VERSION_PLACEHOLDER), lines, fixed = TRUE)
  lines <- gsub("(Start time:).*",   "\\1 <TIMESTAMP>", lines)
  lines <- gsub("(End time:).*",     "\\1 <TIMESTAMP>", lines)
  lines <- gsub("(Elapsed time:).*", "\\1 <ELAPSED>",   lines)
  lines
}

# Compares column by column: numeric columns within REL_TOL, everything else exactly.
# fread reads an all-NA column as logical, so expected values are coerced to the
# actual column's type first.
expect_frame_matches <- function(actual, expected, label) {
  actual <- as.data.frame(actual)
  expected <- as.data.frame(expected)
  expect_identical(names(actual), names(expected), label = paste(label, "column names"))
  expect_identical(nrow(actual), nrow(expected), label = paste(label, "row count"))

  for (col in names(expected)) {
    a <- actual[[col]]
    if (is.numeric(a)) {
      expect_equal(a, as.numeric(expected[[col]]), tolerance = REL_TOL,
                   label = sprintf("%s column '%s'", label, col))
    } else {
      expect_identical(as.character(a), as.character(expected[[col]]),
                       label = sprintf("%s column '%s'", label, col))
    }
  }
}

# The printed output includes data.table's table previews, whose format changes
# between data.table versions (e.g. the column-type row), so it is only compared
# where the package's own dependencies are pinned, not on CRAN's machines.
expect_output_matches <- function(out, case) {
  skip_on_cran()
  expect_identical(normalize_output_lines(out$stdout), read_expected_lines(case, "stdout"),
                   label = paste(case, "stdout"))
  expect_identical(normalize_output_lines(out$stderr), read_expected_lines(case, "stderr"),
                   label = paste(case, "stderr"))
}

test_that("fixture manifest is available", {
  expect_true(nzchar(fixtures_dir))
  expect_true(file.exists(file.path(fixtures_dir, "test_manifest.tsv")))
})

manifest <- read_tsv(file.path(fixtures_dir, "test_manifest.tsv"))

for (i in seq_len(nrow(manifest))) {
  row <- as.list(manifest[i, ])

  test_that(sprintf("fixture case: %s", row$case), {
    out <- run_case(build_args(row))

    if (isTRUE(row$expect_error)) {
      expected_msg <- paste(read_expected_lines(row$case, "error"), collapse = "\n")
      expect_false(is.null(out$error), label = paste(row$case, "raised an error"))
      expect_identical(out$error, expected_msg, label = paste(row$case, "error message"))
    } else {
      expect_null(out$error, label = paste(row$case, "unexpected error"))
      meta <- attr(out$value, "metadata")

      # The expected metadata stores a placeholder in place of the version
      expect_identical(meta$sardine_version, VERSION)
      meta$sardine_version <- VERSION_PLACEHOLDER

      expect_frame_matches(out$value, read_tsv(file.path(expected_dir, paste0(row$case, "_results.tsv"))),
                           paste(row$case, "results"))
      expect_frame_matches(meta, read_tsv(file.path(expected_dir, paste0(row$case, "_metadata.tsv"))),
                           paste(row$case, "metadata"))
    }

    expect_output_matches(out, row$case)
  })
}
