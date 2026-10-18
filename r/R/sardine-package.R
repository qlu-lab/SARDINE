# =============================================================================
# r/R/sardine-package.R — package-level documentation and settings
# =============================================================================

#' @keywords internal
"_PACKAGE"

# The package calls data.table functions with `data.table::` rather than
# importing them, so declare that `[` on a data.table should use data.table
# semantics (e.g. `dt[, -1, with = FALSE]`) inside this package.
.datatable.aware <- TRUE
