# =============================================================================
# python/src/sardine/io.py — readers and column standardisation (mirrors r/R/io.R)
# =============================================================================

import gzip
import os
import re

import pandas as pd

from .utils import SardineError, log_message, log_table


def standardize_column_name(
    dt, pattern, standard_name, col_optional=False, verbose=True
):
    """
    Rename a column to a canonical name using a case-insensitive pattern.

    EWAS summary statistics and clock weight files arrive with a wide variety of column
    names. This resolves a column by pattern and renames it in place, erroring when
    there are 0 or 2+ matches.
    """
    matched_cols = [
        col for col in dt.columns if re.search(pattern, str(col), flags=re.IGNORECASE)
    ]

    if len(matched_cols) == 1:
        dt.rename(columns={matched_cols[0]: standard_name}, inplace=True)
    elif len(matched_cols) > 1:
        raise SardineError(
            "\nERROR! Ambiguous "
            + standard_name
            + " columns found: "
            + ", ".join(str(col) for col in matched_cols)
            + ". Please rename manually."
        )
    else:
        if col_optional:
            log_message(
                verbose,
                "INFO",
                "Optional %s column not found. Looked for pattern: %s",
                standard_name,
                pattern,
            )
        else:
            raise SardineError(
                "\nERROR! Missing "
                + standard_name
                + " column. Looked for pattern: "
                + pattern
            )

    return dt


def _open_text(path):
    return gzip.open(path, "rt") if path.lower().endswith(".gz") else open(path)


def _read_delimited_file(path):
    """Read a delimited file, detecting the delimiter from the header line the way
    data.table::fread() does for the common cases (tab, comma, semicolon, pipe,
    otherwise whitespace)."""
    with _open_text(path) as f:
        header = f.readline()
    counts = {sep: header.count(sep) for sep in ["\t", ",", ";", "|"]}
    sep = max(counts, key=counts.get)
    if counts[sep] == 0:
        sep = r"\s+"
    return pd.read_csv(path, sep=sep, float_precision="round_trip")


def read_delimited_table(input, table_name):
    """
    Coerce a file path or in-memory pandas DataFrame to a DataFrame.
    """
    if isinstance(input, (str, os.PathLike)):
        path = os.fspath(input)
        if os.path.exists(path):
            try:
                return _read_delimited_file(path)
            except Exception as e:
                raise SardineError(
                    "\nERROR!Failed to read " + table_name + " as a table: " + str(e)
                ) from e
        raise SardineError("\nERROR! File not found for " + table_name + ": " + path)

    if isinstance(input, pd.DataFrame):
        # Like R's as.data.table(), the row names (index) are dropped
        return input.reset_index(drop=True)

    raise SardineError(
        "\nERROR! " + table_name + " must be either a file or a pandas DataFrame"
    )


def read_summary_statistics(ewas_sumstats, N=None, verbose=True):
    """
    Read EWAS summary statistics and standardise column names.

    N is an optional user-supplied total sample size. When given it overrides any N
    column present in the file. Returns a DataFrame with columns CpGmarker, BETA, SE and
    N.
    """
    ss = read_delimited_table(ewas_sumstats, table_name="EWAS sumstats")

    log_message(
        verbose,
        "INFO",
        "EWAS sumstats: %d rows and %d columns",
        ss.shape[0],
        ss.shape[1],
    )
    log_table(verbose, "EWAS sumstats header", ss.head(6), digits=3)

    standardize_column_name(
        ss, r"^(cpgmarker|cpg|marker|probe|id)$", "CpGmarker", verbose=verbose
    )
    standardize_column_name(
        ss, r"^(beta|coef|coeff|estimate)$", "BETA", verbose=verbose
    )
    standardize_column_name(
        ss, r"^(se|std\.error|std_error|standard_error)$", "SE", verbose=verbose
    )
    # The N column is only required when the user did not supply N
    standardize_column_name(
        ss,
        r"^(n|samplesize|n_samples)$",
        "N",
        col_optional=N is not None,
        verbose=verbose,
    )

    # If the user supplied a total sample size, override any N column in the file.
    if N is not None:
        if "N" in ss.columns:
            log_message(
                verbose,
                "WARN",
                "User-supplied N overrides the N column in the EWAS sumstats file.",
            )
        ss["N"] = N

    # Remove any CpGs with NA summary stats
    if ss["BETA"].isna().any():
        log_message(
            verbose,
            "INFO",
            "Dropping %d CpGs with NA BETA from the EWAS summary statistics.",
            int(ss["BETA"].isna().sum()),
        )
        ss = ss[~ss["BETA"].isna()].reset_index(drop=True)
    if ss["SE"].isna().any():
        log_message(
            verbose,
            "INFO",
            "Dropping %d CpGs with NA SE from the EWAS summary statistics.",
            int(ss["SE"].isna().sum()),
        )
        ss = ss[~ss["SE"].isna()].reset_index(drop=True)

    return ss


def read_clock_weights(clock, verbose=True):
    """
    Read epigenetic clock weights and standardise column names.

    Returns a DataFrame with columns CpGmarker and CoefficientTraining.
    """
    ck = read_delimited_table(clock, table_name="clock weights")

    log_message(
        verbose,
        "INFO",
        "Clock weights: %d rows and %d columns",
        ck.shape[0],
        ck.shape[1],
    )
    log_table(verbose, "Clock weights header", ck.head(6))

    standardize_column_name(
        ck, r"^(cpgmarker|cpg|marker|probe|id|pred\.var)$", "CpGmarker", verbose=verbose
    )
    standardize_column_name(
        ck,
        r"^(CoefficientTraining|beta|coef|coeff|estimate|weight)$",
        "CoefficientTraining",
        verbose=verbose,
    )

    return ck


def read_mean_beta(mean_beta, verbose=True):
    """
    Read the mean DNAm beta-value reference and standardise column names.

    Used to convert M-value-scale EWAS summary statistics to the beta-value scale (see
    convert_mvalue_ewas_to_beta_scale()). Returns a DataFrame with columns CpGmarker and
    beta_bar.
    """
    mb = read_delimited_table(mean_beta, table_name="mean_beta")

    log_message(
        verbose, "INFO", "mean_beta: %d rows and %d columns", mb.shape[0], mb.shape[1]
    )
    log_table(verbose, "mean_beta header", mb.head(6), digits=3)

    standardize_column_name(
        mb, r"^(cpgmarker|cpg|marker|probe|id)$", "CpGmarker", verbose=verbose
    )
    standardize_column_name(
        mb,
        r"^(beta_bar|mean_beta|beta_mean|mean_dnam_beta|avg_beta)$",
        "beta_bar",
        verbose=verbose,
    )

    return mb


def covariance_table_to_matrix(dt):
    """
    Convert a covariance table with a leading CpG ID column to a matrix: a DataFrame of
    covariances whose index holds the CpG IDs.
    """
    if dt.shape[1] < 2 or not pd.api.types.is_string_dtype(dt.iloc[:, 0]):
        raise SardineError(
            "\nERROR! DNAm_cov_matrix table must have CpG IDs in its first column, "
            "followed by one numeric column per CpG."
        )

    values = dt.iloc[:, 1:]
    non_numeric = [
        str(col)
        for col in values.columns
        if not pd.api.types.is_numeric_dtype(values[col])
    ]
    if len(non_numeric) > 0:
        raise SardineError(
            "\nERROR! DNAm_cov_matrix has non-numeric covariance columns: "
            + ", ".join(non_numeric[:5])
        )

    # Built from one array so the DataFrame is a single block (fast for large matrices)
    return pd.DataFrame(
        values.to_numpy(dtype=float),
        index=dt.iloc[:, 0].tolist(),
        columns=list(values.columns),
    )


def read_dnam_covariance(DNAm_cov_matrix, verbose=True):
    """
    Read the DNAm variance-covariance matrix and validate its row/col names.

    DNAm_cov_matrix is one of:
      - a path to a delimited file (e.g. TSV) whose first column holds CpG IDs
        and whose header lists the same CpG IDs in the same order;
      - a pandas DataFrame of numeric columns with CpG IDs as its index and
        column names;
      - a pandas DataFrame laid out like the delimited file.

    Returns a DataFrame of covariances with identical index and column names.
    """
    is_path = isinstance(DNAm_cov_matrix, (str, os.PathLike))

    if is_path and re.search(
        r"\.rds$", os.fspath(DNAm_cov_matrix), flags=re.IGNORECASE
    ):
        raise SardineError(
            "\nERROR! DNAm_cov_matrix .rds files can only be read by the R package; "
            "save the matrix as a delimited file (e.g. TSV) instead."
        )
    elif isinstance(DNAm_cov_matrix, pd.DataFrame) and all(
        pd.api.types.is_numeric_dtype(DNAm_cov_matrix[col])
        for col in DNAm_cov_matrix.columns
    ):
        # All-numeric DataFrame: CpG IDs are expected in its index
        if isinstance(DNAm_cov_matrix.index, pd.RangeIndex):
            raise SardineError(
                "\nERROR! DNAm_cov_matrix must have valid rownames and colnames corresponding to CpG sites."
            )
        Sigma = DNAm_cov_matrix.astype(float)
    elif is_path or isinstance(DNAm_cov_matrix, pd.DataFrame):
        Sigma = covariance_table_to_matrix(
            read_delimited_table(DNAm_cov_matrix, table_name="DNAm_cov_matrix")
        )
    else:
        raise SardineError(
            "\nERROR! DNAm_cov_matrix must be a delimited file or a pandas DataFrame"
        )

    # Row and column names must be identical.
    if list(Sigma.index) != list(Sigma.columns):
        raise SardineError(
            "\nERROR! DNAm_cov_matrix rownames and colnames must be identical and in the same order."
        )

    log_message(
        verbose,
        "INFO",
        "DNAm covariance matrix: %d rows and %d columns",
        Sigma.shape[0],
        Sigma.shape[1],
    )
    log_table(verbose, "DNAm covariance matrix header", Sigma.iloc[:6, :6], digits=3)

    return Sigma
