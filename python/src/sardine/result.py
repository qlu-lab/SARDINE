# =============================================================================
# python/src/sardine/result.py — the result class and how it prints
# (mirrors r/R/print.R)
# =============================================================================

import pandas as pd

from .utils import sprintf

CORE_COLUMNS = ["model", "param", "coef", "SE", "z", "p", "neglog10p", "N", "note"]


class SardineResult:
    """
    The result of sardine().

    Attributes:
        results: a DataFrame with one row per estimated parameter (lambda and
            gamma always; gamma_logistic for binary phenotypes; log_HR for Cox
            EWAS) and columns model, param, coef, SE, z, p, neglog10p, N and note.
        metadata: a one-row DataFrame with CpG overlap counts and Var(Y)
            diagnostics (mirrors the R package's attr(result, "metadata")).
    """

    def __init__(self, results: pd.DataFrame, metadata: pd.DataFrame):
        self.results = results
        self.metadata = metadata

    def __repr__(self):
        metadata = self.metadata.iloc[0]
        header = sprintf(
            "SARDINE result: %d CpGs used (%.1f%% of the clock)",
            metadata["N_CpGs_final"],
            metadata["ClockCpG_coverage_ratio"] * 100,
        )
        columns = [col for col in CORE_COLUMNS if col in self.results.columns]
        return (
            header
            + "\n\n"
            + self.results[columns].to_string()
            + "\n\n"
            + "Use result.metadata for CpG overlap counts and Var(Y) diagnostics."
        )
