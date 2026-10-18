# =============================================================================
# r/R/io.R — readers and column standardisation
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
#' @param verbose Logical.
#'
#' @return `dt`, with the column renamed (invisible - no console printing).
#' @noRd
standardize_column_name <- function(dt, pattern, standard_name, col_optional = FALSE, verbose = TRUE) {
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
      log_message(verbose, "INFO", "Optional %s column not found. Looked for pattern: %s", standard_name, pattern)
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
#' @param N Optional. User-supplied total sample size. When given it overrides
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

  log_message(verbose, "INFO", "EWAS sumstats: %d rows and %d columns", nrow(ss), ncol(ss))
  log_table(verbose, "EWAS sumstats header", utils::head(ss, 6), digits = 3)

  standardize_column_name(ss, "^(cpgmarker|cpg|marker|probe|id)$", "CpGmarker", verbose = verbose)
  standardize_column_name(ss, "^(beta|coef|coeff|estimate)$", "BETA", verbose = verbose)
  standardize_column_name(ss, "^(se|std\\.error|std_error|standard_error)$", "SE", verbose = verbose)
  # The N column is only required when the user did not supply N
  standardize_column_name(ss, "^(n|samplesize|n_samples)$", "N", col_optional = !is.null(N), verbose = verbose)

  # If the user supplied a total sample size, override any N column in the file.
  if (!is.null(N)) {
    if (!is.null(ss$N)) {
      log_message(verbose, "WARN", "User-supplied N overrides the N column in the EWAS sumstats file.")
    }
    data.table::set(ss, j = "N", value = N)
  }

  # Remove any CpGs with NA summary stats
  if (any(is.na(ss$BETA))) {
    log_message(verbose, "INFO", "Dropping %d CpGs with NA BETA from the EWAS summary statistics.", sum(is.na(ss$BETA)))
    ss <- ss[!is.na(ss$BETA), ]
  }
  if (any(is.na(ss$SE))) {
    log_message(verbose, "INFO", "Dropping %d CpGs with NA SE from the EWAS summary statistics.", sum(is.na(ss$SE)))
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

  log_message(verbose, "INFO", "Clock weights: %d rows and %d columns", nrow(ck), ncol(ck))
  log_table(verbose, "Clock weights header", utils::head(ck, 6))

  standardize_column_name(
    ck, "^(cpgmarker|cpg|marker|probe|id|pred\\.var)$", "CpGmarker", verbose = verbose
  )
  standardize_column_name(
    ck, "^(CoefficientTraining|beta|coef|coeff|estimate|weight)$", "CoefficientTraining", verbose = verbose
  )

  return(ck)
}

#' Read the mean DNAm beta-value reference and standardise column names
#'
#' Used to convert M-value-scale EWAS summary statistics to the beta-value
#' scale (see `convert_mvalue_ewas_to_beta_scale()`).
#'
#' @param mean_beta Path, `data.frame` or `data.table` with columns
#'   `CpGmarker` and `beta_bar`.
#' @param verbose Logical.
#'
#' @return A `data.table` with columns `CpGmarker` and `beta_bar`.
#' @noRd
read_mean_beta <- function(mean_beta, verbose = TRUE) {
  mb <- read_delimited_table(
    mean_beta,
    table_name = "mean_beta"
  )

  log_message(verbose, "INFO", "mean_beta: %d rows and %d columns", nrow(mb), ncol(mb))
  log_table(verbose, "mean_beta header", utils::head(mb, 6), digits = 3)

  standardize_column_name(mb, "^(cpgmarker|cpg|marker|probe|id)$", "CpGmarker", verbose = verbose)
  standardize_column_name(mb, "^(beta_bar|mean_beta|beta_mean|mean_dnam_beta|avg_beta)$", "beta_bar", verbose = verbose)

  return(mb)
}

#' Convert a covariance table with a leading CpG ID column to a matrix
#'
#' @param dt A `data.table` whose first column holds CpG IDs and whose remaining
#'   columns hold the covariances, one column per CpG.
#'
#' @return A numeric matrix with rownames taken from the ID column.
#' @noRd
covariance_table_to_matrix <- function(dt) {
  if (ncol(dt) < 2L || !is.character(dt[[1]])) {
    stop("\nERROR! DNAm_cov_matrix table must have CpG IDs in its first column, ",
         "followed by one numeric column per CpG.", call. = FALSE)
  }

  values <- dt[, -1, with = FALSE]
  non_numeric <- names(values)[!vapply(values, is.numeric, logical(1))]
  if (length(non_numeric) > 0) {
    stop("\nERROR! DNAm_cov_matrix has non-numeric covariance columns: ",
         paste(utils::head(non_numeric, 5), collapse = ", "), call. = FALSE)
  }

  Sigma <- as.matrix(values)
  rownames(Sigma) <- dt[[1]]
  return(Sigma)
}

#' Read the DNAm variance-covariance matrix and validate its row/col names
#'
#' @param DNAm_cov_matrix One of:
#'   \itemize{
#'     \item a path to a delimited file (e.g. TSV) whose first column holds CpG
#'       IDs and whose header lists the same CpG IDs in the same order;
#'     \item a path to an `.rds` file containing a matrix;
#'     \item a matrix with CpG rownames and colnames;
#'     \item a `data.frame` of numeric columns with CpG rownames;
#'     \item a `data.frame`/`data.table` laid out like the delimited file.
#'   }
#' @param verbose Logical.
#'
#' @return A numeric matrix with identical row and column names.
#' @noRd
read_dnam_covariance <- function(DNAm_cov_matrix, verbose = TRUE) {
  is_path <- is.character(DNAm_cov_matrix) && length(DNAm_cov_matrix) == 1L

  if (is.matrix(DNAm_cov_matrix)) {
    Sigma <- DNAm_cov_matrix
  } else if (is_path && grepl("\\.rds$", DNAm_cov_matrix, ignore.case = TRUE)) {
    if (!file.exists(DNAm_cov_matrix)) {
      stop("\nERROR! File not found for DNAm_cov_matrix: ", DNAm_cov_matrix, call. = FALSE)
    }
    Sigma <- tryCatch(
      readRDS(DNAm_cov_matrix),
      error = function(e) {
        stop("\nERROR! Failed to read DNAm_cov_matrix as .rds: ", conditionMessage(e), call. = FALSE)
      }
    )
    if (!is.matrix(Sigma)) {
      stop("\nERROR! DNAm_cov_matrix .rds must contain a matrix", call. = FALSE)
    }
  } else if (is.data.frame(DNAm_cov_matrix) && all(vapply(DNAm_cov_matrix, is.numeric, logical(1)))) {
    # All-numeric data.frame: CpG IDs are expected in its rownames
    Sigma <- as.matrix(DNAm_cov_matrix)
  } else if (is_path || is.data.frame(DNAm_cov_matrix)) {
    Sigma <- covariance_table_to_matrix(read_delimited_table(DNAm_cov_matrix, table_name = "DNAm_cov_matrix"))
  } else {
    stop("\nERROR! DNAm_cov_matrix must be a delimited file, an .rds file, a matrix, ",
         "or a data.frame/data.table", call. = FALSE)
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

  log_message(verbose, "INFO", "DNAm covariance matrix: %d rows and %d columns", nrow(Sigma), ncol(Sigma))
  log_table(verbose, "DNAm covariance matrix header",
            Sigma[seq_len(min(6, nrow(Sigma))), seq_len(min(6, ncol(Sigma)))], digits = 3)

  return(Sigma)
}
