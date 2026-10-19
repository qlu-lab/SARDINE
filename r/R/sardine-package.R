# =============================================================================
# r/R/sardine-package.R — package-level documentation and settings

# Estimates the association between an epigenetic clock and a phenotype
# using only EWAS summary statistics, the clock's CpG weights, and a
# reference DNA methylation covariance matrix. No individual-level data is
# required.
# =============================================================================

#' @keywords internal
#' @importFrom R.utils gunzip
"_PACKAGE"

# R.utils is imported (above) only so that data.table::fread() can read
# .gz-compressed input files; it is not called directly.

# The package calls data.table functions with `data.table::` rather than
# importing them, so declare that `[` on a data.table should use data.table
# semantics (e.g. `dt[, -1, with = FALSE]`) inside this package.
.datatable.aware <- TRUE
