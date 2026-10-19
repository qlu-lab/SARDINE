# =============================================================================
# python/src/sardine/estimate.py — CpG alignment and the core clock/phenotype
# regression estimators (mirrors r/R/estimate.R)
# =============================================================================

import numpy as np

from .utils import (SardineError, intersect, log_message, sprintf, sqrt_or_na,
                    two_sided_normal_pvalue)

_INTERCEPT_PATTERN = r"^\(?intercept\)?$"


def _is_intercept(cpgs):
    return cpgs.astype(str).str.contains(_INTERCEPT_PATTERN, case=False, regex=True).to_numpy()


def drop_intercept_rows(clock_weights, verbose=True):
    """Drop intercept rows from a clock weight table.

    Clock files frequently carry an "(Intercept)" row that is not a CpG and
    must not be counted towards clock coverage.
    """
    is_intercept = _is_intercept(clock_weights["CpGmarker"])
    n_intercept_rows = int(is_intercept.sum())
    if n_intercept_rows > 0:
        log_message(verbose, "INFO", "Removed %d intercept row(s) from the clock weights.", n_intercept_rows)
        clock_weights = clock_weights[~is_intercept].reset_index(drop=True)
    return clock_weights


def summarize_cpg_overlap(ss, cw, Sigma, verbose=True):
    """Summarize CpG overlap across the three inputs.

    cw must already have its intercept row removed. Returns a dict of counts
    plus cpg_intersect, the CpGs common to all three.
    """
    # Check that there are no intercept rows in the clock weights
    if _is_intercept(cw["CpGmarker"]).any():
        raise SardineError("[ERROR] Clock weights contain intercept rows. Please remove them before calling summarize_cpg_overlap().")

    N_CpGs_ewas = len(ss["CpGmarker"])
    N_CpGs_clock = len(cw["CpGmarker"])
    N_CpGs_ref = Sigma.shape[0]

    N_CpGs_ewas_clock = len(intersect(ss["CpGmarker"], cw["CpGmarker"]))
    N_CpGs_ewas_ref = len(intersect(ss["CpGmarker"], Sigma.index))
    N_CpGs_clock_ref = len(intersect(cw["CpGmarker"], Sigma.index))

    cpg_intersect = intersect(intersect(ss["CpGmarker"], cw["CpGmarker"]), Sigma.index)

    N_CpGs_final = len(cpg_intersect)

    ClockCpG_coverage_ratio = N_CpGs_final / N_CpGs_clock
    ClockCpG_missed_ratio = 1 - ClockCpG_coverage_ratio

    log_message(verbose, "INFO", "Number of CpGs in EWAS: %s", f"{N_CpGs_ewas:,}")
    log_message(verbose, "INFO", "Number of CpGs for clock: %s", f"{N_CpGs_clock:,}")
    log_message(verbose, "INFO", "Number of CpGs in DNAm matrix: %s", f"{N_CpGs_ref:,}")

    log_message(verbose, "INFO", "Number of overlapping CpGs between EWAS and clock: %s", f"{N_CpGs_ewas_clock:,}")
    log_message(verbose, "INFO", "Number of overlapping CpGs between EWAS and DNAm matrix: %s", f"{N_CpGs_ewas_ref:,}")
    log_message(verbose, "INFO", "Number of overlapping CpGs between clock and DNAm matrix: %s", f"{N_CpGs_clock_ref:,}")

    log_message(verbose, "INFO", "Final number of overlapping CpGs: %s", f"{N_CpGs_final:,}")
    log_message(verbose, "INFO", "Percentage of missed clock CpGs: %s%%", round(ClockCpG_missed_ratio, 2) * 100)

    if N_CpGs_final == 0:
        raise SardineError("\nERROR! Zero overlapping CpGs found between EWAS, Clock, and Covariance Matrix.")

    return {
        "N_CpGs_ewas": N_CpGs_ewas,
        "N_CpGs_clock": N_CpGs_clock,
        "N_CpGs_ref": N_CpGs_ref,
        "N_CpGs_ewas_clock": N_CpGs_ewas_clock,
        "N_CpGs_ewas_ref": N_CpGs_ewas_ref,
        "N_CpGs_clock_ref": N_CpGs_clock_ref,
        "cpg_intersect": cpg_intersect,
        "N_CpGs_final": N_CpGs_final,
        "ClockCpG_coverage_ratio": ClockCpG_coverage_ratio,
        "ClockCpG_missed_ratio": ClockCpG_missed_ratio,
    }


def cap_phenotype_variance(varY, w_T_beta, w_T_Sigma_w, verbose=True):
    """Cap Var(Y) at the largest value compatible with the estimated clock variance.

    Var(Y) cannot exceed w'Sigma w / (w'beta)^2, which is the value attained
    when the clock and the phenotype are perfectly correlated. Exceeding it
    produces a negative residual variance and a missing standard error. As
    necessary, the estimate is shrunk just below this bound.

    Returns a dict with the (possibly capped) varY, the diagnostic varY_ratio
    (original varY / maximum allowed varY), and a note describing any capping.
    """
    note_varY_capping = None
    varY_ratio = np.nan

    # Check for valid inputs
    if varY is None:
        raise SardineError("[ERROR] varY is NULL; cannot cap Var(Y) without a current estimate.")
    elif not np.isfinite(varY):
        raise SardineError("[ERROR] varY is not finite; cannot cap Var(Y) without a valid estimate.")
    elif varY <= 0:
        raise SardineError("[ERROR] varY is not positive; cannot cap Var(Y) without a valid estimate.")
    elif w_T_beta is None or w_T_Sigma_w is None:
        raise SardineError("[ERROR] w_T_beta or w_T_Sigma_w is NULL; cannot cap Var(Y) without clock statistics.")
    elif not np.isfinite(w_T_beta) or not np.isfinite(w_T_Sigma_w):
        raise SardineError("[ERROR] w_T_beta or w_T_Sigma_w is not finite; cannot cap Var(Y) without valid clock statistics.")
    elif w_T_Sigma_w <= 0:
        raise SardineError("[ERROR] w_T_Sigma_w is not positive; cannot cap Var(Y) when the clock has no variance.")

    if w_T_beta == 0:
        note_varY_capping = "w'beta = 0; cannot cap Var(Y) when the clock has no correlation with the phenotype."
        return {"varY": varY, "varY_ratio": varY_ratio, "note": note_varY_capping}

    varY_max = w_T_Sigma_w / w_T_beta ** 2
    varY_ratio = varY / varY_max

    if varY > varY_max:
        # Cap varY to just below the max allowed
        msg_capping = sprintf(
            "Estimated varY (%.6g) exceeds maximum (%.6g), ratio = varY/varY_max = %.6f. Capped to avoid negative error variance. If the ratio is much larger than 1: check that DNAm_cov_matrix matches the EWAS population and covariate adjustment.",
            varY, varY_max, varY_ratio
        )

        log_message(verbose, "INFO", "%s", msg_capping)
        varY = varY_max * (1 - 1e-8)
        note_varY_capping = msg_capping

    return {"varY": varY, "varY_ratio": varY_ratio, "note": note_varY_capping}


def estimate_gamma(w_T_beta, w_T_Sigma_w, varY, N):
    """Estimate gamma: the effect of the clock on the phenotype, Y = C * gamma + e.

    gamma_hat is derived directly from w_T_beta, w_T_Sigma_w, and varY (summary
    stats only): gamma_hat = w_T_beta * varY / w_T_Sigma_w. Its standard error
    uses the residual variance sigma_e2 = varY - gamma_hat^2 * w_T_Sigma_w (also
    returned for use in convert_gamma_to_logistic_scale()).

    Returns a dict with estimate, variance, se (NaN if variance < 0), z, p,
    sigma_e2 and note ("var_gamma_hat < 0" if variance came out negative, else None).
    """
    gamma_hat = w_T_beta * varY / w_T_Sigma_w
    sigma_e2 = varY - gamma_hat ** 2 * w_T_Sigma_w
    var_gamma_hat = sigma_e2 / ((N - 2) * w_T_Sigma_w)

    note_gamma = "var_gamma_hat < 0" if var_gamma_hat < 0 else None

    se_gamma_hat = sqrt_or_na(var_gamma_hat)
    with np.errstate(divide="ignore", invalid="ignore"):
        z_gamma_hat = gamma_hat / se_gamma_hat

    return {
        "estimate": gamma_hat,
        "variance": var_gamma_hat,
        "se": se_gamma_hat,
        "z": z_gamma_hat,
        "p": two_sided_normal_pvalue(z_gamma_hat),
        "sigma_e2": sigma_e2,
        "note": note_gamma,
    }


def estimate_lambda(w_T_beta, w_T_Sigma_w, varY, N):
    """Estimate lambda: the effect of the phenotype on the clock, C = Y * lambda + u.

    The estimate lambda_hat is just w_T_beta. Its standard error uses the
    residual variance sigma_u2 = w_T_Sigma_w - lambda_hat^2 * varY.

    Returns a dict with estimate, variance, se (NaN if variance < 0), z, p and
    note ("var_lambda_hat < 0" if variance came out negative, else None).
    """
    lambda_hat = w_T_beta
    sigma_u2 = w_T_Sigma_w - lambda_hat ** 2 * varY
    var_lambda_hat = sigma_u2 / ((N - 2) * varY)

    note_lambda = "var_lambda_hat < 0" if var_lambda_hat < 0 else None

    se_lambda_hat = sqrt_or_na(var_lambda_hat)
    with np.errstate(divide="ignore", invalid="ignore"):
        z_lambda_hat = lambda_hat / se_lambda_hat

    return {
        "estimate": lambda_hat,
        "variance": var_lambda_hat,
        "se": se_lambda_hat,
        "z": z_lambda_hat,
        "p": two_sided_normal_pvalue(z_lambda_hat),
        "note": note_lambda,
    }


def convert_gamma_to_logistic_scale(gamma, varY, w_T_Sigma_w, varY_ratio, verbose=True):
    """Convert the linear-scale gamma to a log-odds-ratio (logistic) scale.

    The point estimate is gamma / sigma_e^2. Note: this conversion introduces
    some downward bias to the point estimate as the true effect gets larger.

    The standard error follows by the delta method. The Z statistic and p-value
    are inherited from the linear-scale gamma (H0: gamma = 0 and
    H0: gamma_logistic = 0 are the same hypothesis). Note: This means
    Z != coef/SE. If Z was computed from the delta-method SE, it would be
    attenuated towards 0 when the association is strong.

    Returns a dict with estimate, variance, se, z, p and note.
    """
    log_message(verbose, "INFO", "Converting final gamma estimate from linear to logistic (log-odds) scale...")

    if (not np.isfinite(gamma["sigma_e2"]) or gamma["sigma_e2"] <= 0
            or (np.isfinite(varY_ratio) and varY_ratio > 0.99)):
        note_gamma_logistic = sprintf(
            "sigma_e2 = %.6g (varY_ratio = %.6g); logistic-scale gamma not identifiable.",
            gamma["sigma_e2"], varY_ratio
        )
        log_message(verbose, "WARN", "%s", note_gamma_logistic)

        return {"estimate": np.nan, "variance": np.nan, "se": np.nan,
                "z": np.nan, "p": np.nan, "note": note_gamma_logistic}

    gamma_hat_logistic = gamma["estimate"] / gamma["sigma_e2"]

    # Delta method for computing SE of logistic-scale gamma
    # var(gamma_hat_logistic) = var(gamma_hat) * (d gamma_hat_logistic / d gamma_hat)^2
    dg = (varY + gamma["estimate"] ** 2 * w_T_Sigma_w) / gamma["sigma_e2"] ** 2  # derivative of gamma_hat_logistic with respect to gamma_hat
    var_gamma_hat_logistic = gamma["variance"] * dg ** 2
    se_gamma_hat_logistic = sqrt_or_na(var_gamma_hat_logistic)

    note_gamma_logistic = "SE by delta-method; Z and P inherited from the linear-scale gamma"

    if np.isfinite(var_gamma_hat_logistic) and var_gamma_hat_logistic < 0:
        note_gamma_logistic = "var_gamma_hat_logistic < 0; " + note_gamma_logistic

    log_message(verbose, "NOTE", (
        "For gamma in a logistic regression of Phenotype ~ Clock: "
        "Z and P are inherited from the linear-scale gamma; "
        "SE computed from delta-method, so Z != coef/SE for this row. "
        "CI should be obtained by transforming the 'linear' gamma CI endpoints."
    ))

    return {
        "estimate": gamma_hat_logistic,
        "variance": var_gamma_hat_logistic,
        "se": se_gamma_hat_logistic,
        "z": gamma["z"],
        "p": gamma["p"],
        "note": note_gamma_logistic,
    }


def compute_cox_log_hazard_ratio(gamma, kappa_cox, varY_type, verbose=True):
    """Compute the clock-level log hazard ratio.

    Convert from linear scale gamma to Cox log-hazard-ratio by dividing by the
    cox conversion factor kappa_cox. This is simply undoing the conversion
    applied to the summary stats during convert_cox_ewas_to_linear_scale().
       gamma_hat_cox = gamma_hat / kappa_cox
       se_gamma_hat_cox = se_gamma_hat / kappa_cox

    varY_type ("inferred" or "user_provided") determines whether the
    kappa_cox/cox_pi dependence cancels. Returns a dict with estimate, se, z, p
    and note.
    """
    gamma_hat_cox = gamma["estimate"] / kappa_cox
    se_gamma_hat_cox = gamma["se"] / kappa_cox

    if varY_type == "inferred":
        note_gamma_cox = "Clock-level log hazard ratio; invariant to cox_pi. First-order approximation, attenuated toward the null."
    else:
        note_gamma_cox = "Clock-level log hazard ratio. WARNING: varY was user-supplied, so the cox_pi cancellation does NOT hold; this estimate requires varY = pi(1-pi) with a correct cox_pi."

    log_message(verbose, "INFO", "Clock-level log HR = %.6g (SE %.6g); HR = %.4f", gamma_hat_cox, se_gamma_hat_cox, np.exp(gamma_hat_cox))

    return {
        "estimate": gamma_hat_cox,
        "se": se_gamma_hat_cox,
        "z": gamma["z"],
        "p": gamma["p"],
        "note": note_gamma_cox,
    }
