# =============================================================================
# R/src/io.R — readers and column standardisation
# =============================================================================

#' Rename a column to a canonical name using a case-insensitive pattern
#'
#' EWAS summary statistics and clock weight files arrive with a wide variety
#' of column names. This resolves a column by pattern and renames it in place,
#' erroring when there are 0 or 2+ matches.
#'
#' @param dt A `data.table`, modified in place.
#' @param pattern A regular expression matched against the column names.
#' @param standard_name The name the matched column is renamed to.
#' @param col_optional Logical. If `TRUE`, no error is thrown if no match is found. Default FALSE.
#'
#' @return `dt`, with the column renamed (invisible - no console printing).
#' @noRd
standardize_column_name <- function(dt, pattern, standard_name, col_optional = FALSE) {
  col_names <- names(dt)
  matched_cols <- grep(pattern, col_names, ignore.case = TRUE, value = TRUE)

  if (length(matched_cols) == 1) {
    data.table::setnames(dt, matched_cols, standard_name)
  } else if (length(matched_cols) > 1) {
    stop(paste0(
      "\nERROR! Ambiguous ", standard_name, " columns found: ",
      paste(matched_cols, collapse = ", "), ". Please rename manually."
    ), call. = FALSE)
  } else {
    if (col_optional) {
      if (verbose) {
        cat("[INFO] Optional ", standard_name, " column not found. Looked for pattern: ", pattern, "\n", sep = "")
      }
    } else {
      stop(paste0(
        "\nERROR! Missing ", standard_name, " column. Looked for pattern: ", pattern
      ), call. = FALSE)
    }
  }

  invisible(dt)
}

#' Coerce a file path or in-memory table to a data.table
#'
#' @param input A file path, `data.frame`, or `data.table`.
#' @param table_name Name of table being read, for error messages.
#' @param read_fail_prefix Prefix used when `fread()` itself errors; the
#'   underlying `conditionMessage()` is appended after it.
#' @param invalid_type_message Full message used when `input` is none of the
#'   accepted types.
#'
#' @return A `data.table`.
#' @noRd
read_delimited_table <- function(input, table_name) {
  if (is.character(input) && length(input) == 1L) {
    if (file.exists(input)) {
      return(tryCatch(
        data.table::fread(input),
        error = function(e) {
          stop("\nERROR!Failed to read ", table_name, " as a data.table: ", conditionMessage(e), call. = FALSE)
        }
      ))
    }
    else {
      stop("\nERROR! File not found for ", table_name, ": ", input, call. = FALSE)
    }
  }

  if (is.data.frame(input) || data.table::is.data.table(input)) {
    return(data.table::as.data.table(input))
  }

  stop("\nERROR! ", table_name, " must be either a file, data.frame or data.table", call. = FALSE)
}

#' Read EWAS summary statistics and standardise column names
#'
#' @param ewas_sumstats A file path, `data.frame`, or `data.table`.
#' @param N [optional]. User-supplied total sample size. When given it overrides
#'   any N column present in the file.
#' @param verbose Logical.
#'
#' @return A `data.table` with columns `CpGmarker`, `BETA`, `SE` and `N`.
#' @noRd
read_summary_statistics <- function(ewas_sumstats, N = NULL, verbose = TRUE) {
  ss <- read_delimited_table(
    ewas_sumstats,
    table_name = "EWAS sumstats"
  )

  if (verbose) {
    cat("[INFO] EWAS sumstats:", nrow(ss), "rows and", ncol(ss), "columns\n")
    cat("[INFO] EWAS sumstats header:\n")
    print(head(ss, 6), digits = 3)
  }

  standardize_column_name(ss, "^(cpgmarker|cpg|marker|probe|id)$", "CpGmarker")
  standardize_column_name(ss, "^(beta|coef|coeff|estimate)$", "BETA")
  standardize_column_name(ss, "^(se|std\\.error|std_error|standard_error)$", "SE")
  standardize_column_name(ss, "^(n|samplesize|n_samples)$", "N")

  # If the user supplied a total sample size, override any N column in the file.
  if (!is.null(N)) {
    if (!is.null(ss$N)) {
      cat("[INFO] User-supplied N overrides the N column in the EWAS sumstats file.", call. = FALSE)
    }
    data.table::set(ss, j = "N", value = N)
  }

  # Remove any CpGs with NA summary stats
  if (any(is.na(ss$BETA))) {
    if (verbose) cat("[INFO] Dropping", sum(is.na(ss$BETA)), "CpGs with NA BETA from the EWAS summary statistics.\n")
    ss <- ss[!is.na(ss$BETA), ]
  }
  if (any(is.na(ss$SE))) {
    if (verbose) cat("[INFO] Dropping", sum(is.na(ss$SE)), "CpGs with NA SE from the EWAS summary statistics.\n")
    ss <- ss[!is.na(ss$SE), ]
  }

  return(ss)
}

#' Read epigenetic clock weights and standardise column names
#'
#' @param clock Path, `data.frame` or `data.table`.
#' @param verbose Logical.
#'
#' @return A `data.table` with columns `CpGmarker` and `CoefficientTraining`.
#' @noRd
read_clock_weights <- function(clock, verbose = TRUE) {
  ck <- read_delimited_table(
    clock,
    table_name = "clock weights"
  )

  if (verbose) {
    cat("[INFO] Clock weights:", nrow(ck), "rows and", ncol(ck), "columns\n")
    cat("[INFO] Clock weights header:\n")
    print(head(ck, 6))
  }

  standardize_column_name(
    ck, "^(cpgmarker|cpg|marker|probe|id|pred\\.var)$", "CpGmarker"
  )
  standardize_column_name(
    ck, "^(CoefficientTraining|beta|coef|coeff|estimate|weight)$", "CoefficientTraining"
  )

  return(ck)
}

#' Read the DNAm variance-covariance matrix and validate its row/col names
#'
#' @param DNAm_cov_matrix Either a path to an `.rds` file containing a
#'   matrix, or a matrix.
#' @param verbose Logical.
#'
#' @return A numeric matrix with identical row and column names.
#' @noRd
read_dnam_covariance <- function(DNAm_cov_matrix, verbose = TRUE) {
  if (is.character(DNAm_cov_matrix) && length(DNAm_cov_matrix) == 1L) {
    if (file.exists(DNAm_cov_matrix)) {
      Sigma <- tryCatch(
        readRDS(DNAm_cov_matrix),
        error = function(e) {
          stop("Failed to read DNAm_cov_matrix as .rds: ", conditionMessage(e), call. = FALSE)
        }
      )
      if (!is.matrix(Sigma)) {
        stop("\nERROR! DNAm_cov_matrix .rds must contain a matrix", call. = FALSE)
      }
  } else {
      stop("\nERROR! DNAm_cov_matrix file not found", call. = FALSE)
    }
  } else if (is.matrix(DNAm_cov_matrix)) {
    Sigma <- DNAm_cov_matrix
  } else {
    stop("\nERROR! DNAm_cov_matrix must be either a file path or a matrix", call. = FALSE)
  }

  if (is.null(rownames(Sigma)) || is.null(colnames(Sigma))) {
    stop("\nERROR! DNAm_cov_matrix must have valid rownames and colnames corresponding to CpG sites.",
         call. = FALSE)
  }

  # Row and column names must be identical.
  if (!identical(rownames(Sigma), colnames(Sigma))) {
    stop("\nERROR! DNAm_cov_matrix rownames and colnames must be identical and in the same order.",
         call. = FALSE)
  }

  if (verbose) {
    cat("[INFO] DNAm covariance matrix:", nrow(Sigma), "rows and", ncol(Sigma), "columns\n")
    cat("[INFO] DNAm covariance matrix header:\n")
    print(Sigma[seq_len(min(6, nrow(Sigma))), seq_len(min(6, ncol(Sigma)))], digits = 3)
  }

  return(Sigma)
}
