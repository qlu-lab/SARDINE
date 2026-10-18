# SARDINE test fixture runner
#
# Runs every case in test_manifest.csv against the current R code, then diffs
# the output against the expected output in fixtures/expected/*. Relative
# error tolerance of 1e-6.
#
# Usage: Rscript fixtures/run_fixture_tests.R
# Exits with status 0 if every case passes, 1 otherwise.

suppressMessages(library(data.table))
options(warn = 1) # print warnings immediately, so they land inside the right case's sink

# Tolerance allowed between actual and expected outputs
REL_TOL = 1e-6

this_script_path  = sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
root         = normalizePath(file.path(dirname(this_script_path), ".."))
fixtures_dir = file.path(root, "fixtures")
expected_dir = file.path(fixtures_dir, "expected")

package_files = c("utils.R", "io.R", "conversions.R", "estimate.R", "sardine.R", "print.R")
for (f in package_files) {
  script_path_i = file.path(root, "R", "src", f)
  stopifnot(file.exists(script_path_i))
  invisible(capture.output(source(script_path_i)))
}
stopifnot(exists("sardine"))

manifest = fread(file.path(fixtures_dir, "test_manifest.tsv"), sep = "\t", na.strings = "")

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
# comparison against the saved fixtures is stable across runs.
# ============================================================
normalize_output_lines = function(lines) {
  lines = gsub("(Start time:).*",   "\\1 <TIMESTAMP>", lines)
  lines = gsub("(End time:).*",     "\\1 <TIMESTAMP>", lines)
  lines = gsub("(Elapsed time:).*", "\\1 <ELAPSED>",    lines)
  lines
}

read_expected_lines = function(case, kind) {
  readLines(file.path(expected_dir, paste0(case, "_", kind, ".txt")), warn = FALSE)
}

# ============================================================
# Computes the relative difference between two numbers
# ============================================================
rel_diff = function(a, b) {
  ifelse(is.na(a) & is.na(b), 0,
  ifelse(xor(is.na(a), is.na(b)), Inf,
  ifelse(a == b, 0, abs(a - b) / pmax(abs(a), abs(b), 1e-300))))
}

# ============================================================
# Compares actual and expected data frames
# ============================================================
compare_frames = function(actual, expected, case, kind) {
  problems = character(0)

  if (!identical(names(actual), names(expected))) {
    problems = c(problems, sprintf("[%s/%s] column names differ: actual=%s expected=%s",
                                    case, kind,
                                    paste(names(actual), collapse = ","),
                                    paste(names(expected), collapse = ",")))
    return(problems)
  }
  if (nrow(actual) != nrow(expected)) {
    problems = c(problems, sprintf("[%s/%s] row count differs: actual=%d expected=%d",
                                    case, kind, nrow(actual), nrow(expected)))
    return(problems)
  }

  for (col in names(expected)) {
    a = actual[[col]]
    e = expected[[col]]
    if (is.numeric(e)) {
      d = rel_diff(as.numeric(a), as.numeric(e))
      bad = which(d > REL_TOL)
      if (length(bad) > 0) {
        problems = c(problems, sprintf(
          "[%s/%s] column '%s' row(s) %s exceed rel tol %.1e: actual=%s expected=%s",
          case, kind, col, paste(bad, collapse = ","), REL_TOL,
          paste(signif(a[bad], 8), collapse = ";"),
          paste(signif(e[bad], 8), collapse = ";")))
      }
    } else {
      a_chr = as.character(a); a_chr[is.na(a_chr)] = ""
      e_chr = as.character(e); e_chr[is.na(e_chr)] = ""
      bad = which(a_chr != e_chr)
      if (length(bad) > 0) {
        problems = c(problems, sprintf(
          "[%s/%s] column '%s' row(s) %s differ: actual=%s expected=%s",
          case, kind, col, paste(bad, collapse = ","),
          paste(a_chr[bad], collapse = ";"), paste(e_chr[bad], collapse = ";")))
      }
    }
  }
  problems
}

# ============================================================
# Compares actual and expected captured output, line by line
# ============================================================
compare_lines = function(actual, expected, case, kind) {
  problems = character(0)

  if (length(actual) != length(expected)) {
    problems = c(problems, sprintf(
      "[%s/%s] line count differs: actual=%d expected=%d\n  actual:\n    %s\n  expected:\n    %s",
      case, kind, length(actual), length(expected),
      paste(actual, collapse = "\n    "), paste(expected, collapse = "\n    ")))
    return(problems)
  }

  bad = which(actual != expected)
  if (length(bad) > 0) {
    problems = c(problems, sprintf(
      "[%s/%s] line(s) %s differ:\n  actual:   %s\n  expected: %s",
      case, kind, paste(bad, collapse = ","),
      paste(actual[bad], collapse = " | "), paste(expected[bad], collapse = " | ")))
  }
  problems
}

# ============================================================
# Runs through all test cases
# ============================================================
n_pass = 0
n_fail = 0
all_problems = character(0)

for (i in seq_len(nrow(manifest))) {
  row  = manifest[i]
  case = row$case
  args = build_args(row)

  out = run_case(args)
  actual_stdout = normalize_output_lines(out$stdout)
  actual_stderr = normalize_output_lines(out$stderr)

  if (isTRUE(row$expect_error)) {
    expected_msg_path = file.path(expected_dir, paste0(case, "_error.txt"))
    expected_msg = paste(readLines(expected_msg_path, warn = FALSE), collapse = "\n")

    problems = character(0)
    if (out$ok) {
      problems = c(problems, sprintf("[%s] expected an error, but succeeded", case))
    } else if (!identical(out$message, expected_msg)) {
      problems = c(problems, sprintf(
        "[%s] error message changed:\n  actual:   %s\n  expected: %s",
        case, out$message, expected_msg))
    }

    problems = c(
      problems,
      compare_lines(actual_stdout, read_expected_lines(case, "stdout"), case, "stdout"),
      compare_lines(actual_stderr, read_expected_lines(case, "stderr"), case, "stderr")
    )

    if (length(problems) == 0) {
      n_pass = n_pass + 1
      cat(sprintf("[run_fixture_tests] PASS  %-22s (error message matches)\n", case))
    } else {
      n_fail = n_fail + 1
      all_problems = c(all_problems, problems)
      cat(sprintf("[run_fixture_tests] FAIL  %-22s (%d issue(s))\n", case, length(problems)))
    }
    next
  }

  if (!out$ok) {
    n_fail = n_fail + 1
    all_problems = c(all_problems, sprintf("[%s] unexpected error: %s", case, out$message))
    cat(sprintf("[run_fixture_tests] FAIL  %-22s unexpected error: %s\n", case, out$message))
    next
  }

  actual_res  = as.data.table(out$value)
  actual_meta = as.data.table(attr(out$value, "metadata"))

  expected_res  = fread(file.path(expected_dir, paste0(case, "_results.tsv")), sep = "\t", na.strings = "")
  expected_meta = fread(file.path(expected_dir, paste0(case, "_metadata.tsv")), sep = "\t", na.strings = "")

  problems = c(
    compare_frames(actual_res,  expected_res,  case, "results"),
    compare_frames(actual_meta, expected_meta, case, "metadata"),
    compare_lines(actual_stdout, read_expected_lines(case, "stdout"), case, "stdout"),
    compare_lines(actual_stderr, read_expected_lines(case, "stderr"), case, "stderr")
  )

  if (length(problems) == 0) {
    n_pass = n_pass + 1
    cat(sprintf("[run_fixture_tests] PASS  %-22s\n", case))
  } else {
    n_fail = n_fail + 1
    all_problems = c(all_problems, problems)
    cat(sprintf("[run_fixture_tests] FAIL  %-22s (%d issue(s))\n", case, length(problems)))
  }
}

cat(sprintf("\n[run_fixture_tests] %d passed, %d failed (of %d cases)\n", n_pass, n_fail, nrow(manifest)))

if (n_fail > 0) {
  cat("\n[run_fixture_tests] ---- details ----\n")
  cat(paste(all_problems, collapse = "\n\n"), "\n")
  quit(status = 1)
} else {
  quit(status = 0)
}
