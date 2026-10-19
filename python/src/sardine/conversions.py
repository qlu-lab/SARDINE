# =============================================================================
# python/src/sardine/conversions.py — scale / model / direction conversions of
# EWAS sumstats and phenotypic variance inference (mirrors r/R/conversions.R)
# =============================================================================

import numpy as np

from .utils import SardineError, intersect, log_message, match, sprintf


def _diag_by_name(Sigma, cpgs):
    """Diagonal of Sigma (the per-CpG variances) for the given CpG IDs, NaN where
    a CpG is not in Sigma (R's diag(Sigma)[cpgs])."""
    diag = np.diag(Sigma.to_numpy())
    idx = match(cpgs, Sigma.index)
    return np.where(idx >= 0, diag[np.maximum(idx, 0)], np.nan)


def convert_mvalue_ewas_to_beta_scale(ss, mean_beta_df, ewas_direction, verbose=True):
    """Convert M-value-scale EWAS summary statistics to the beta-value scale.

    Uses the first-order approximation dM/dbeta = 1 / (log(2) * beta * (1 -
    beta)) evaluated at the mean methylation of each CpG. The approximation
    degrades at extreme methylation, so very high or very low CpGs should be
    filtered first.

    Best if the mean beta values are computed from the same samples as the EWAS
    itself, since the conversion is a local linear approximation evaluated at
    these values.

    Returns the summary statistics on the beta-value scale, restricted to CpGs
    present in mean_beta_df.
    """
    log_message(verbose, "INFO", "Converting M-value EWAS summary statistics to Beta-value scale.")

    # Check for extreme mean population beta values
    num_cpgs_low = int((mean_beta_df["beta_bar"] < 0.1).sum())
    num_cpgs_high = int((mean_beta_df["beta_bar"] > 0.9).sum())
    num_cpgs_total = mean_beta_df.shape[0]

    if num_cpgs_low > 0 or num_cpgs_high > 0:
        log_message(verbose, "WARN", "Mean beta reference contains %d CpGs with beta < 0.1 and %d CpGs with beta > 0.9 (out of %d total). Conversion from M-value to beta-value scale may be inaccurate for these extreme CpGs.", num_cpgs_low, num_cpgs_high, num_cpgs_total)

    match_idx = match(ss["CpGmarker"], mean_beta_df["CpGmarker"])
    missing_from_mean_beta_idx = match_idx < 0

    if missing_from_mean_beta_idx.any():
        log_message(verbose, "INFO", "Dropping %d CpGs from the EWAS summary stats because they are missing from the mean_beta reference.", int(missing_from_mean_beta_idx.sum()))
        ss = ss[~missing_from_mean_beta_idx].reset_index(drop=True)
        match_idx = match_idx[~missing_from_mean_beta_idx]

    ss = ss.copy()
    beta_bar = mean_beta_df["beta_bar"].to_numpy(dtype=float)[match_idx]
    cc = np.log(2) * beta_bar * (1 - beta_bar)

    if ewas_direction == "M~Y":
        ss["BETA"] = cc * ss["BETA"]
        ss["SE"] = cc * ss["SE"]
    else:  # Y ~ M:
        ss["BETA"] = ss["BETA"] / cc
        ss["SE"] = ss["SE"] / cc

    return ss


def convert_cox_ewas_to_linear_scale(ss, cox_pi, cox_E, verbose=True):
    """Rescale Cox summary statistics onto the linear Y ~ M scale.

    Multiplies log hazard ratios by kappa(pi) = (1 - pi) * (-log(1 - pi)),
    where pi is the proportion of cases at the end of follow-up. This is a
    first-order approximation using the delta method.
        BETA' = kappa(pi) * BETA
        SE(BETA') = kappa(pi) * SE(BETA)

    Returns a dict with the rescaled ss, the cox_pi used, and kappa_cox.
    """
    log_message(verbose, "INFO", "Converting Cox regression summary statistics to linear Y ~ M model...")

    if cox_pi is None and cox_E is None:
        raise SardineError("[ERROR]: EWAS is Cox PH model, must provide either `cox_pi` or 'cox_E' and 'N'.")

    if cox_pi is None:
        cox_N = np.nanmedian(ss["N"].to_numpy(dtype=float))
        cox_pi = cox_E / cox_N

    # Check for reasonable cox_pi values.
    cox_pi = np.float64(cox_pi)
    if not np.isfinite(cox_pi) or cox_pi <= 0 or cox_pi >= 1:
        raise SardineError(sprintf("[ERROR] cox_pi must be strictly between 0 and 1; got %.6g.", cox_pi))

    kappa_cox = (1 - cox_pi) * (-np.log(1 - cox_pi))

    log_message(verbose, "INFO", "Cox conversion factor kappa(pi) = %.6g (pi = %.6g)", kappa_cox, cox_pi)

    ss = ss.copy()
    ss["BETA"] = ss["BETA"] * kappa_cox
    ss["SE"] = ss["SE"] * kappa_cox

    return {"ss": ss, "cox_pi": cox_pi, "kappa_cox": kappa_cox}


def convert_linear_y_on_m_ewas_to_m_on_y(ss, verbose=True):
    """Convert linear Y ~ M summary statistics to linear M ~ Y.

    Invert the regression per CpG using summary statistics:
       BETA' = BETA / (BETA^2 + (N - 2) * SE(BETA)^2)
       SE(BETA') = SE(BETA) / (BETA^2 + (N - 2) * SE(BETA)^2)
    Derived from standard OLS formulas:
       BETA = Cov(Y,M) / Var(M)
       SE(BETA)^2 ≈ sigma_e^2 / ((N-2)·Var(M))
    """
    log_message(verbose, "INFO", "Converting linear Y ~ M summary statistics to M ~ Y direction...")

    denom = ss["BETA"] ** 2 + (ss["N"] - 2) * ss["SE"] ** 2

    ss = ss.copy()
    ss["BETA"] = ss["BETA"] / denom
    ss["SE"] = ss["SE"] / denom

    return ss


def convert_logistic_y_on_m_ewas_to_linear_m_on_y(ss, Sigma, varY, verbose=True):
    """Convert logistic Y ~ M summary statistics to linear M ~ Y.

    Conversion uses the reference methylation variance:
       Let X = sqrt(1 + 4 * SIGMA * BETA^2 * Var(Y))
       BETA' = 2 * SIGMA * BETA / (1 + X)
       SE(BETA') = SE(BETA) * (2 * SIGMA / (1 + X)) / X

    Assumes small effects of Y on M and equal methylation variance in cases
    and controls.
    """
    log_message(verbose, "INFO", "Converting logistic Y ~ M summary statistics to linear M ~ Y...")
    log_message(verbose, "NOTE", "Logistic to linear conversion assumes small effect sizes of Y on M and equal methylation variance between cases and controls.")

    # Convert logistic regression coefficients to linear regression coefficients
    Sigma_aligned = _diag_by_name(Sigma, ss["CpGmarker"])
    aa = np.sqrt(1 + 4 * ss["BETA"].to_numpy() ** 2 * Sigma_aligned * varY)
    a2b = 2 * Sigma_aligned / (1 + aa)
    ss = ss.copy()
    ss["BETA"] = ss["BETA"] * a2b
    ss["SE"] = ss["SE"] * (a2b / aa)  # delta-method

    return ss


def infer_phenotype_variance_from_logistic_y_on_m_ewas(ss, Sigma, verbose=True):
    """Infer Var(Y) from logistic Y ~ M summary statistics and a reference
    covariance matrix.

    For each overlapping CpG, Var(Y) = 1/(N * Var(M) * SE(BETA)^2). The median
    across CpGs is taken for stability.

    Note: If there are not enough CpGs overlapping both Sigma and the EWAS
    sumstats (e.g. if the CpGs were already subset to the CpGs in a single
    clock), then the Var(Y) estimate may be unstable/inaccurate/negative.
    """
    # Only use the cpgs in both the EWAS sumstats and the reference Sigma matrix
    cpg_subset = intersect(ss["CpGmarker"], Sigma.index)
    if len(cpg_subset) == 0:
        raise SardineError("[ERROR] Var(Y) cannot be inferred: no overlapping CpGs between the EWAS sumstats and DNAm_cov_matrix.")
    ss_subset = ss.iloc[match(cpg_subset, ss["CpGmarker"])]
    diag_subset = _diag_by_name(Sigma, cpg_subset)
    with np.errstate(divide="ignore"):
        varY_hat = np.nanmedian(1 / (ss_subset["N"].to_numpy(dtype=float) * diag_subset * ss_subset["SE"].to_numpy() ** 2))

    log_message(verbose, "INFO", "Inferred covariate-adjusted Var(Y) from logistic EWAS sumstats and DNAm_cov_matrix: %.6f", varY_hat)
    return varY_hat


def infer_phenotype_variance_from_linear_m_on_y_ewas(ss, Sigma, verbose=True):
    """Infer Var(Y) from linear M ~ Y summary statistics and a reference
    covariance matrix.

    For each overlapping CpG, Var(Y) = Var(M) / (BETA^2 + (N - 2) * SE(BETA)^2).
    The median across CpGs is taken for stability.

    Note: If there are not enough CpGs overlapping both Sigma and the EWAS
    sumstats (e.g. if the CpGs were already subset to the CpGs in a single
    clock), then the Var(Y) estimate may be unstable/inaccurate/negative.
    """
    if not all(col in ss.columns for col in ["SE", "N"]):
        raise SardineError("[ERROR] Var(Y) cannot be inferred from EWAS sumstats because SE and/or N columns are missing.")

    # Only use the cpgs in both the EWAS sumstats and the reference Sigma matrix
    cpg_subset = intersect(ss["CpGmarker"], Sigma.index)
    if len(cpg_subset) == 0:
        raise SardineError("[ERROR] Var(Y) cannot be inferred: no overlapping CpGs between the EWAS sumstats and DNAm_cov_matrix.")
    ss_subset = ss.iloc[match(cpg_subset, ss["CpGmarker"])]
    diag_subset = _diag_by_name(Sigma, cpg_subset)

    beta = ss_subset["BETA"].to_numpy(dtype=float)
    se = ss_subset["SE"].to_numpy(dtype=float)
    n = ss_subset["N"].to_numpy(dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        varY_hats = diag_subset / (beta ** 2 + (n - 2) * se ** 2)
    varY_hat = np.nanmedian(varY_hats)

    log_message(verbose, "INFO", "Inferred covariate-adjusted Var(Y) from linear EWAS sumstats and DNAm_cov_matrix: %.6f", varY_hat)
    return varY_hat
