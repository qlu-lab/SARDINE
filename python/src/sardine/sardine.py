# =============================================================================
# python/src/sardine/sardine.py: Exported, user-facing function (mirrors r/R/sardine.R)
# =============================================================================

import numbers
import os
import re
from datetime import datetime
from importlib import resources

import numpy as np
import pandas as pd

from .conversions import (
    convert_cox_ewas_to_linear_scale,
    convert_linear_y_on_m_ewas_to_m_on_y,
    convert_logistic_y_on_m_ewas_to_linear_m_on_y,
    convert_mvalue_ewas_to_beta_scale,
    infer_phenotype_variance_from_linear_m_on_y_ewas,
    infer_phenotype_variance_from_logistic_y_on_m_ewas,
)
from .estimate import (
    cap_phenotype_variance,
    compute_cox_log_hazard_ratio,
    convert_gamma_to_logistic_scale,
    drop_intercept_rows,
    estimate_gamma,
    estimate_lambda,
    summarize_cpg_overlap,
)
from .io import (
    read_clock_weights,
    read_dnam_covariance,
    read_mean_beta,
    read_summary_statistics,
)
from .result import SardineResult
from .utils import (
    SardineError,
    append_note,
    format_elapsed_time,
    log_message,
    match,
    match_option,
    sardine_version,
    sprintf,
    two_sided_neg_log10_pvalue,
)


def _is_number(x):
    """
    True for a real number, dissallows bools to match R's is.numeric() behavior.
    """
    return isinstance(x, numbers.Real) and not isinstance(x, bool)


def sardine(
    ewas_sumstats,
    *,
    ewas_direction="M~Y",
    ewas_model="linear",
    ewas_DNAm_scale="beta",
    ewas_pheno_type="continuous",
    clock,
    DNAm_cov_matrix,
    mean_beta=None,
    varY=None,
    scale=False,
    age_adult=20,
    N=None,
    cox_pi=None,
    cox_E=None,
    verbose=True,
):
    """
    Summary-statistics-based association between epigenetic clock and phenotype.

    Estimates the association between an epigenetic clock and a phenotype using
    only EWAS summary statistics, the clock's CpG weights, and a reference DNA
    methylation covariance matrix. No individual-level data is required.

    Two coefficients are always returned: lambda, from the regression of the
    clock on the phenotype, and gamma, from the regression of the phenotype on
    the clock. Binary phenotypes additionally yield gamma_logistic, a log-odds
    ratio. Cox-model inputs additionally yield log_HR, a clock-level log hazard
    ratio. All rows share the same Z statistic and p-value by construction.

    All arguments after ewas_sumstats must be given by name. See
    https://github.com/qlu-lab/SARDINE/blob/main/inputs.md for details on every input.

    Args:
        ewas_sumstats: EWAS summary statistics: a path to a delimited file or a
            pandas DataFrame; must contain columns CpG, BETA, SE, and N.
        ewas_direction: Direction of the association (outcome ~ exposure) in the
            EWAS model: "M~Y" (methylation regressed on phenotype; default) or
            "Y~M" (phenotype regressed on methylation).
        ewas_model: EWAS regression model: "linear" (default), "logistic", or
            "cox" (Cox proportional hazards).
        ewas_DNAm_scale: Scale of the EWAS DNA methylation data: "beta" (beta
            values; default) or "m" (M-values).
        ewas_pheno_type: EWAS phenotype type: "continuous" (default) or
            "binary". If binary, will also output an estimate of the log-odds
            ratio of clock on phenotype. For all four options above, matching
            ignores case and whitespace (e.g. "y ~ m" is accepted as "Y~M").
        clock: Epigenetic clock weights: a path to a delimited file or a pandas
            DataFrame; must contain columns CpG, CoefficientTraining.
        DNAm_cov_matrix: Variance-covariance matrix for DNAm. Must have the same
            covariate adjustments as the ewas summary stats. Can be a path to a
            delimited table (e.g. TSV) whose first column holds the CpG IDs and
            whose header lists the same CpG IDs in the same order; a pandas
            DataFrame of numeric columns with CpG IDs as its index and column
            names; or a DataFrame laid out like the delimited file.
        mean_beta: Optional. A path to a delimited file or a DataFrame with
            columns CpGmarker and beta_bar, giving the mean DNAm beta-value at
            each CpG site. Only needed for M-value to beta-value conversion.
            Best if computed from the same samples as the EWAS itself.
        varY: Variance of the phenotype used in EWAS: raw phenotype if no
            covariates; residualized phenotype if adjusted for covariates in
            the EWAS model. Default None: estimate varY using only EWAS sumstats
            and the DNAm reference matrix.
        scale: For Horvath or PedBE clocks, we can scale the point estimate and
            SE; doesn't affect the p-value.
        age_adult: Only used if scale=True.
        N: Optional EWAS total sample size. Default is to read this from the
            EWAS summary statistics.
        cox_pi: Only used when ewas_model="cox". Proportion of cases at the end
            of the study. If not provided, this is computed as cox_E/N. Exactly
            one of cox_pi and cox_E must be provided if ewas_model="cox".
        cox_E: Only used when ewas_model="cox". Number of cases at the end of
            the study. Exactly one of cox_pi and cox_E must be provided if
            ewas_model="cox".
        verbose: Whether to print progress and diagnostic messages (to stderr).

    Returns:
        A SardineResult: result.results holds one row per estimated parameter,
        and result.metadata the CpG overlap counts and Var(Y) diagnostics.

    Raises:
        SardineError: if the inputs are invalid.

    Note on cox models:
        For a Cox Proportional Hazards model, ewas_direction must be "Y~M",
        ewas_model must be "cox", ewas_pheno_type must be "binary", and N should
        represent the total sample size at the end of the study. Exactly one of
        cox_pi and cox_E must be provided. cox_pi = cox_E/N is the proportion of
        cases at the end of the study. It is recommended that varY is inferred
        (left as None) rather than user-supplied, so that the log_HR estimate is
        invariant to cox_pi.

    Note on logistic models:
        For a logistic model, ewas_direction must be "Y~M", ewas_model must be
        "logistic", and ewas_pheno_type must be "binary".
    """
    # Record start time
    start_time = datetime.now()

    log_message(
        verbose,
        None,
        "%s",
        "\n".join(
            [
                "************************************************",
                sprintf("* SARDINE v%s", sardine_version()),
                "* (C) Ryan Moreno <rrmoreno@wisc.edu> and",
                "*     Yuchang Wu <ywu423@wisc.edu>, ",
                "*     University of Wisconsin-Madison, 2026",
                "* GNU General Public License v3",
                "************************************************",
            ]
        ),
    )

    log_message(
        verbose, "INFO", "Start time: %s", start_time.strftime("%Y-%m-%d %H:%M:%S")
    )

    # ============================================================
    # Check the validity of the parameters
    # ============================================================

    # Resolve each option to its canonical value (case and whitespace are ignored)
    ewas_direction = match_option(ewas_direction, ["M~Y", "Y~M"], "ewas_direction")
    ewas_model = match_option(ewas_model, ["linear", "logistic", "cox"], "ewas_model")
    ewas_DNAm_scale = match_option(ewas_DNAm_scale, ["beta", "m"], "ewas_DNAm_scale")
    ewas_pheno_type = match_option(
        ewas_pheno_type, ["continuous", "binary"], "ewas_pheno_type"
    )

    # mean_beta must be provided if ewas_DNAm_scale is m-value
    if ewas_DNAm_scale == "m" and mean_beta is None:
        raise SardineError(
            "[ERROR] 'mean_beta' must be provided when 'ewas_DNAm_scale' is an M-value scale."
        )

    # if ewas_model is cox, ewas_direction must be Y~M, ewas_pheno_type must be binary,
    # and exactly one of cox_pi and cox_E must be set
    if ewas_model == "cox":
        if ewas_direction != "Y~M":
            raise SardineError(
                "[ERROR] 'ewas_direction' must be 'Y~M' when 'ewas_model' is 'cox'."
            )
        if ewas_pheno_type != "binary":
            raise SardineError(
                "[ERROR] 'ewas_pheno_type' must be 'binary' when 'ewas_model' is 'cox'."
            )
        if (cox_pi is None) == (cox_E is None):
            raise SardineError(
                "[ERROR] Exactly one of 'cox_pi' or 'cox_E' must be provided when 'ewas_model' is 'cox'."
            )

        # Warn that it's recommended varY be inferred for Cox models
        if varY is not None:
            log_message(
                verbose,
                "WARN",
                "ewas_model = 'cox' with user-supplied varY: `kappa_cox` no longer cancels, so `log_hr` is sensitive to `cox_pi`. Leaving varY unset is recommended.",
            )

    # if ewas_model is logistic, ewas_direction must be Y~M and ewas_pheno_type must be binary
    if ewas_model == "logistic" and ewas_direction != "Y~M":
        raise SardineError(
            "[ERROR] 'ewas_direction' must be 'Y~M' when 'ewas_model' is 'logistic'."
        )
    if ewas_model == "logistic" and ewas_pheno_type != "binary":
        raise SardineError(
            "[ERROR] 'ewas_pheno_type' must be 'binary' when 'ewas_model' is 'logistic'."
        )

    # scale and age_adult should be set together
    if (scale is None) != (age_adult is None):
        raise SardineError("[ERROR] 'scale' and 'age_adult' should be set together.")

    # Alert the user if they seem to be using a Horvath or PedBE clock based
    # on the file name but did not scale
    if isinstance(clock, (str, os.PathLike)):
        clock_name = re.sub(r"\.csv", "", os.path.basename(os.fspath(clock)))

        if verbose and not scale and re.search("horvath|pedbe", clock_name.lower()):
            log_message(
                verbose,
                "WARN",
                "Your clock weight filename contains '%s' yet you did not scale the point estimate and SE (Z-score and P-value are not affected). Make sure this was intentional.",
                clock_name,
            )

    # Check numeric types
    if N is not None and not _is_number(N):
        raise SardineError("[ERROR] 'N' must be numeric.")
    if cox_E is not None and not _is_number(cox_E):
        raise SardineError("[ERROR] 'cox_E' must be numeric.")
    if cox_E is not None and not cox_E > 0:
        raise SardineError("[ERROR] 'cox_E' must be greater than 0.")
    if cox_pi is not None and not _is_number(cox_pi):
        raise SardineError("[ERROR] 'cox_pi' must be numeric.")
    if age_adult is not None and not _is_number(age_adult):
        raise SardineError("[ERROR] 'age_adult' must be numeric.")
    if varY is not None and not _is_number(varY):
        raise SardineError("[ERROR] 'varY' must be numeric.")
    if varY is not None and varY < 0:
        raise SardineError("[ERROR] 'varY' must be non-negative")
    if varY is not None:
        log_message(
            verbose,
            "WARN",
            "User-supplied varY must be the variance of the phenotype after adjusting for the same covariates as the EWAS (the raw phenotype variance only if the EWAS had no covariates). Leave varY unset to estimate it from the EWAS summary statistics instead.",
        )

    # Use numpy floats from here on, so that division by zero gives Inf/NaN as in R
    varY = None if varY is None else np.float64(varY)

    # ============================================================
    # Read and standardize input tables
    # ============================================================
    ss = read_summary_statistics(ewas_sumstats, N=N, verbose=verbose)
    cw = read_clock_weights(clock, verbose=verbose)
    Sigma = read_dnam_covariance(DNAm_cov_matrix, verbose=verbose)
    mean_beta_df = (
        read_mean_beta(mean_beta, verbose=verbose) if mean_beta is not None else None
    )

    # ============================================================
    # M-value to Beta-value transformation, as necessary
    # ============================================================
    if ewas_DNAm_scale == "m":
        ss = convert_mvalue_ewas_to_beta_scale(
            ss, mean_beta_df, ewas_direction, verbose=verbose
        )

    # ============================================================
    # Convert Cox to linear Y ~ M (Cox input only)
    # ewas_direction for cox must be "Y~M"
    # ============================================================
    kappa_cox = np.nan  # scaling factor for Cox model

    if ewas_model == "cox":
        cox_conversion = convert_cox_ewas_to_linear_scale(
            ss, cox_pi, cox_E, verbose=verbose
        )
        ss = cox_conversion["ss"]
        cox_pi = cox_conversion["cox_pi"]
        kappa_cox = cox_conversion["kappa_cox"]

    # ============================================================
    # Estimate Var(Y) for logistic regression if not provided by the user.
    # Note: In the logistic case, it is more accurate to estimate Var(Y) before
    # converting the summary stats to linear M~Y.
    # ============================================================
    varY_type = "user_provided" if varY is not None else None
    if ewas_model == "logistic" and varY is None:
        varY = infer_phenotype_variance_from_logistic_y_on_m_ewas(
            ss, Sigma, verbose=verbose
        )
        varY_type = "inferred"

    # ============================================================
    # Convert Y ~ M to M ~ Y, as necessary
    # ============================================================
    if ewas_direction == "Y~M":
        if ewas_model in ("linear", "cox"):
            # note: Cox sumstats have already been converted to linear scale
            ss = convert_linear_y_on_m_ewas_to_m_on_y(ss, verbose=verbose)
        elif ewas_model == "logistic":
            ss = convert_logistic_y_on_m_ewas_to_linear_m_on_y(
                ss, Sigma, varY, verbose=verbose
            )
        else:
            raise SardineError(
                "Error: 'ewas_model' must be either 'linear', 'logistic', or 'cox'."
            )

    # ============================================================
    # Estimate Var(Y) if not provided or already estimated
    # ============================================================
    if varY is None:
        # note: all EWAS sumstats are now formatted as linear M~Y
        varY = infer_phenotype_variance_from_linear_m_on_y_ewas(
            ss, Sigma, verbose=verbose
        )
        varY_type = "inferred"

    # ============================================================
    # Align CpGs and compute clock/phenotype regression terms
    # ============================================================

    # Use the CpG intersection between EWAS sumstats, clock weights, and Sigma
    cw = drop_intercept_rows(cw, verbose=verbose)
    overlap = summarize_cpg_overlap(ss, cw, Sigma, verbose=verbose)
    cpg_intersect = overlap["cpg_intersect"]

    # Align data to the same CpG order
    ss = ss.iloc[match(cpg_intersect, ss["CpGmarker"])].reset_index(drop=True)
    cw = cw.iloc[match(cpg_intersect, cw["CpGmarker"])].reset_index(drop=True)
    Sigma = Sigma.loc[cpg_intersect, cpg_intersect]

    # Find the sample size of the EWAS summary stats if not provided by the user
    if N is None:
        # median N for the CpG sites we're interested in
        N = np.nanmedian(ss["N"].to_numpy(dtype=float))
    N = np.float64(N)

    # Compute clock-phenotype regression coefficient
    w = cw["CoefficientTraining"].to_numpy(dtype=float)
    w_T_beta = np.sum(w * ss["BETA"].to_numpy(dtype=float))

    # Compute the variance of the clock sum on the reference population
    w_T_Sigma_w = np.float64(w @ Sigma.to_numpy() @ w)
    if not np.isfinite(w_T_Sigma_w) or w_T_Sigma_w <= 0:
        raise SardineError(
            sprintf(
                "[ERROR] Estimated clock variance w'Sigma w = %.6g is not positive. Check that DNAm_cov_matrix is a valid covariance matrix on the same CpG scale as the EWAS.",
                w_T_Sigma_w,
            )
        )

    # Cap the varY (otherwise, this can result in negative variance estimates for gamma and lambda)
    capping = cap_phenotype_variance(varY, w_T_beta, w_T_Sigma_w, verbose=verbose)
    varY = capping["varY"]
    varY_ratio = capping["varY_ratio"]

    # ============================================================
    # Linear regression: Y = C * gamma + e
    # ============================================================
    gamma = estimate_gamma(w_T_beta, w_T_Sigma_w, varY, N)

    # ============================================================
    # Logistic regression Y ~ C (only if phenotype is binary)
    # ============================================================
    if ewas_pheno_type == "binary":
        if ewas_model == "linear":
            log_message(
                verbose,
                "NOTE",
                "EWAS is a linear regression with a binary phenotype; we can still estimate log-odds ratio of Clock on this binary phenotype.",
            )
        gamma_logistic = convert_gamma_to_logistic_scale(
            gamma, varY, w_T_Sigma_w, varY_ratio, verbose=verbose
        )

        # Add a note to clarify log odds for cumulative incidence versus log hazard ratio for Cox model
        if ewas_model == "cox":
            gamma_logistic["note"] = append_note(
                gamma_logistic["note"],
                "Log-odds ratio for cumulative incidence at end of follow-up; param = 'log_HR' is log of hazard ratio.",
            )
    else:
        gamma_logistic = {
            "estimate": np.nan,
            "variance": np.nan,
            "se": np.nan,
            "z": np.nan,
            "p": np.nan,
            "note": None,
        }

    # ============================================================
    # Linear regression: C = Y * lambda + u
    # ============================================================
    lambda_ = estimate_lambda(w_T_beta, w_T_Sigma_w, varY, N)

    # ============================================================
    # Re-scale effect size and SE, as necessary
    # This is necessary for Horvath1, Horvath2, and PedBE clocks
    # because they apply a scale factor to their clocks.
    # ============================================================
    if scale:
        scale_factor = 1 + age_adult

        gamma["estimate"] = gamma["estimate"] / scale_factor
        gamma["variance"] = gamma["variance"] / (scale_factor**2)
        gamma["se"] = gamma["se"] / scale_factor

        gamma_logistic["estimate"] = gamma_logistic["estimate"] / scale_factor
        gamma_logistic["variance"] = gamma_logistic["variance"] / (scale_factor**2)
        gamma_logistic["se"] = gamma_logistic["se"] / scale_factor

        lambda_["estimate"] = lambda_["estimate"] * scale_factor
        lambda_["variance"] = lambda_["variance"] * (scale_factor**2)
        lambda_["se"] = lambda_["se"] * scale_factor

        log_message(
            verbose,
            "INFO",
            "Applied clock scaling factor = %s to point estimate and SE",
            scale_factor,
        )

    # ============================================================
    # Clock-level log hazard ratio (Cox input only)
    # ============================================================
    if ewas_model == "cox":
        cox_log_hr = compute_cox_log_hazard_ratio(
            gamma, kappa_cox, varY_type, verbose=verbose
        )
    else:
        cox_log_hr = {
            "estimate": np.nan,
            "se": np.nan,
            "z": np.nan,
            "p": np.nan,
            "note": None,
        }

    # ============================================================
    # Format the regression results into a SardineResult
    # ============================================================
    def result_row(model, param, estimate):
        return {
            "model": model,
            "param": param,
            "coef": estimate["estimate"],
            "SE": estimate["se"],
            "z": estimate["z"],
            "p": estimate["p"],
            "N": N,
            "note": estimate["note"],
        }

    rows = [
        result_row("Clock ~ Phenotype (linear)", "lambda", lambda_),
        result_row("Phenotype ~ Clock (linear)", "gamma", gamma),
    ]

    if not np.isnan(gamma_logistic["estimate"]):
        rows.append(
            result_row("Phenotype ~ Clock (logistic)", "gamma_logistic", gamma_logistic)
        )

    # Note: Cox pathway's robustness ordering is log_HR > gamma >> lambda, gamma_logistic
    # they all have the same Z score and p-value though
    if not np.isnan(cox_log_hr["estimate"]):
        rows.append(result_row("Phenotype ~ Clock (Cox)", "log_HR", cox_log_hr))

    res = pd.DataFrame(rows)

    # Add a -log10p column to avoid underflow issues with very small p-values
    res["neglog10p"] = two_sided_neg_log10_pvalue(res["z"].to_numpy(dtype=float))
    res = res[["model", "param", "coef", "SE", "z", "p", "neglog10p", "N", "note"]]

    # Add metadata (number of CpGs, varY diagnostics, etc.)
    res_metadata = pd.DataFrame(
        [
            {
                "sardine_version": sardine_version(),
                "N_CpGs_ewas": overlap["N_CpGs_ewas"],
                "N_CpGs_clock": overlap["N_CpGs_clock"],
                "N_CpGs_DNAmRef": overlap["N_CpGs_ref"],
                "N_CpGs_ewas_clock": overlap["N_CpGs_ewas_clock"],
                "N_CpGs_ewas_ref": overlap["N_CpGs_ewas_ref"],
                "N_CpGs_clock_ref": overlap["N_CpGs_clock_ref"],
                "N_CpGs_final": overlap["N_CpGs_final"],
                "ClockCpG_coverage_ratio": overlap["ClockCpG_coverage_ratio"],
                "ClockCpG_missed_ratio": overlap["ClockCpG_missed_ratio"],
                "cox_pi": cox_pi if ewas_model == "cox" else np.nan,
                "kappa_cox": kappa_cox if ewas_model == "cox" else np.nan,
                "varY": varY,
                "varY_type": varY_type,
                "varY_ratio": varY_ratio,  # varY/varY_max
                "varY_note": capping["note"],
            }
        ]
    )

    # Record end time and calculate elapsed time
    end_time = datetime.now()
    elapsed_seconds = (end_time - start_time).total_seconds()

    log_message(verbose, None, "")
    log_message(verbose, "INFO", "Analysis complete!")
    log_message(verbose, "INFO", "End time: %s", end_time.strftime("%Y-%m-%d %H:%M:%S"))
    log_message(
        verbose, "INFO", "Elapsed time: %s", format_elapsed_time(elapsed_seconds)
    )
    log_message(verbose, None, "")

    return SardineResult(res, res_metadata)


def example(verbose=False):
    """
    Run sardine() on the small synthetic dataset bundled with the package.

    The counterpart of example(sardine) in the R package; useful to check that the
    installation works.
    """
    inputs = resources.files("sardine") / "example_data"
    return sardine(
        str(inputs / "ewas_sumstats_baseline_synthetic.tsv"),
        clock=str(inputs / "clock_weights_baseline_synthetic.tsv"),
        DNAm_cov_matrix=str(inputs / "sigma_baseline_synthetic.tsv"),
        verbose=verbose,
    )
