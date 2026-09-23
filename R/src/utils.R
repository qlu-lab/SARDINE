# =============================================================================
# R/src/utils.R — logging and small numeric helpers shared across the package
# =============================================================================

#' Emit a tagged log message, honouring the verbose flag
#'
#' Every diagnostic message in the package goes through this function so that
#' `verbose = FALSE` reliably silences all except for warnings.
#'
#' @param verbose Logical. When `FALSE` the message is suppressed.
#' @param level Character tag, one of "INFO", "NOTE" or "WARN".
#' @param fmt A `sprintf()` format string.
#' @param ... Values interpolated into `fmt`.
#'
#' @return Invisibly `NULL`
log_message <- function(verbose, level, fmt, ...) {
  if (toupper(level) == "WARN") {
    log_message(verbose, "INFO", fmt, ...) # Print the message as INFO too
    warning(sprintf("[%s] %s", level, sprintf(fmt, ...)), call. = FALSE)
    return(invisible(NULL))
  }

  if (!isTRUE(verbose)) {
    return(invisible(NULL))
  }
  message(sprintf("[%s] %s", level, sprintf(fmt, ...)))
  return(invisible(NULL))
}

#' Append a note to an existing note string without producing "NA; message"
#'
#' @param x Character vector of existing notes, possibly `NA`.
#' @param msg A single note to append, possibly `NA`.
#'
#' @return `x` with `msg` appended where both are non-missing.
#' @noRd
append_note <- function(x, msg) {
  if (is.na(msg)) {
    return(x)
  }
  ifelse(is.na(x), msg, paste0(x, "; ", msg))
}

#' Square root that returns NA rather than warning on negative input
#'
#' Note: Negative variance estimates are possible if Var(Y) is too
#' large relative to the clock variance.
#'
#' @param x Numeric vector.
#' @return `sqrt(x)`, with `NaN` where `x < 0`.
#' @noRd
sqrt_or_na <- function(x) {
  suppressWarnings(sqrt(x))
}

#' Two-sided normal p-value from a Z statistic
#' @noRd
two_sided_normal_pvalue <- function(z) {
  2 * stats::pnorm(abs(z), lower.tail = FALSE)
}

#' Two-sided -log10 p-value, computed in log space to avoid underflow to zero
#'
#' @noRd
two_sided_neg_log10_pvalue <- function(z) {
  -(stats::pnorm(abs(z), lower.tail = FALSE, log.p = TRUE) + log(2)) / log(10)
}

#' Format an elapsed duration in seconds as "1 hour, 2 minutes, 3 seconds"
#' @noRd
format_elapsed_time <- function(elapsed_seconds) {
  hours <- floor(elapsed_seconds / 3600)
  minutes <- floor((elapsed_seconds - hours * 3600) / 60)
  seconds <- round(elapsed_seconds - hours * 3600 - minutes * 60)

  parts <- c()
  if (hours > 0) parts <- c(parts, paste0(hours, " hour", ifelse(hours == 1, "", "s")))
  if (minutes > 0) parts <- c(parts, paste0(minutes, " minute", ifelse(minutes == 1, "", "s")))
  if (seconds > 0 || length(parts) == 0) {
    parts <- c(parts, paste0(seconds, " second", ifelse(seconds == 1, "", "s")))
  }
  paste(parts, collapse = ", ")
}
