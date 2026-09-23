# =============================================================================
# R/src/print.R — print method for the S3 result class
# =============================================================================

#' Print a SARDINE result
#'
#' @param x An object of class `sardine_result`.
#' @param ... Passed to `print.data.frame()`.
#'
#' @return `x`, invisibly.
#' @export
print.sardine_result <- function(x, ...) {
  metadata <- attr(x, "metadata")

  cat(sprintf("SARDINE result: %d CpGs used (%.1f%% of the clock)\n\n",
              metadata$N_CpGs_final, metadata$ClockCpG_coverage_ratio * 100))

  core_columns <- c("model", "param", "coef", "SE", "Z", "P", "log10_P", "N", "note")
  print(as.data.frame(x)[, intersect(core_columns, names(x)), drop = FALSE], ...)

  cat("\nUse attr(x, \"metadata\") for CpG overlap counts and Var(Y) diagnostics.\n")
  invisible(x)
}
