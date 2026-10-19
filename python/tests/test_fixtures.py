# Runs every case in the shared fixture manifest and checks it against the
# expected outputs: the results and metadata (numbers within REL_TOL), the error
# message for cases that should fail, and the printed stdout/stderr line by line
# (except the table previews; see drop_table_previews()).
# Mirrors r/tests/testthat/test-fixtures.R.
#
# The fixtures are read straight from the repo-level fixtures/ directory, or from
# tests/fixtures/ inside an unpacked sdist (see _fixtures_dir.py).

import contextlib
import io
import math
import re
import warnings

import numpy as np
import pandas as pd
import pytest
from _fixtures_dir import FIXTURES_DIR

import sardine
from sardine import SardineError

EXPECTED_DIR = FIXTURES_DIR / "expected"

# Tolerance allowed between actual and expected numbers.
REL_TOL = 1e-6

# The version is saved as this placeholder in the expected stdout and metadata,
# so a version bump doesn't require regenerating the expected outputs
VERSION = sardine.__version__
VERSION_PLACEHOLDER = "<VERSION>"


def read_tsv(path):
    return pd.read_csv(path, sep="\t", float_precision="round_trip")


def read_expected_lines(case, kind):
    text = (EXPECTED_DIR / f"{case}_{kind}.txt").read_text()
    text = text.removesuffix("\n")
    return text.split("\n") if text else []


def build_args(row):
    args = {
        "ewas_sumstats": str(FIXTURES_DIR / row["ewas_sumstats"]),
        "clock": str(FIXTURES_DIR / row["clock"]),
        "DNAm_cov_matrix": str(FIXTURES_DIR / row["DNAm_cov_matrix"]),
        "ewas_direction": row["ewas_direction"],
        "ewas_model": row["ewas_model"],
        "ewas_DNAm_scale": row["ewas_DNAm_scale"],
        "ewas_pheno_type": row["ewas_pheno_type"],
        "scale": bool(row["scale"]),
        "age_adult": row["age_adult"],
        "verbose": True,
    }
    if not pd.isna(row["mean_beta"]):
        args["mean_beta"] = str(FIXTURES_DIR / row["mean_beta"])
    for name in ["varY", "N", "cox_pi", "cox_E"]:
        if not pd.isna(row[name]):
            args[name] = float(row[name])
    return args


def run_case(args):
    """
    Runs sardine() and captures its result (or error), its stdout, and its stderr.

    Warnings are recorded in stderr as "Warning: <message>", in order with the log,
    which is the form the expected files use.
    """
    stdout, stderr = io.StringIO(), io.StringIO()
    value, error = None, None

    def record_warning(message, category, filename, lineno, file=None, line=None):
        stderr.write(f"Warning: {message}\n")

    with warnings.catch_warnings():
        warnings.simplefilter("always")
        warnings.showwarning = record_warning
        with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
            try:
                value = sardine.sardine(**args)
            except SardineError as e:
                error = str(e)

    def lines(buffer):
        text = buffer.getvalue()
        text = text.removesuffix("\n")
        return text.split("\n") if text else []

    return {
        "value": value,
        "error": error,
        "stdout": lines(stdout),
        "stderr": lines(stderr),
    }


def normalize_output_lines(lines):
    """
    Blanks out wall-clock-dependent text (start/end time, elapsed time) and the SARDINE
    version so the comparison is stable across runs and releases.
    """
    lines = [
        line.replace(f"* SARDINE v{VERSION}", f"* SARDINE v{VERSION_PLACEHOLDER}")
        for line in lines
    ]
    lines = [re.sub(r"(Start time:).*", r"\1 <TIMESTAMP>", line) for line in lines]
    lines = [re.sub(r"(End time:).*", r"\1 <TIMESTAMP>", line) for line in lines]
    lines = [re.sub(r"(Elapsed time:).*", r"\1 <ELAPSED>", line) for line in lines]
    return lines


def drop_table_previews(lines):
    """
    Removes the rows printed under each "[INFO] ...

    header:" line. Their layout comes from each language's own table printing, so they
    can't be compared across R and Python. The header line itself is kept.
    """
    kept = []
    in_preview = False
    for line in lines:
        if in_preview and (line == "" or re.match(r"^(\[|Warning:)", line)):
            in_preview = False
        if not in_preview:
            kept.append(line)
        if re.match(r"^\[INFO\] .* header:$", line):
            in_preview = True
    return kept


def is_missing(x):
    return x is None or (isinstance(x, float) and math.isnan(x))


def assert_frame_matches(actual, expected, label):
    """
    Compares column by column: numeric columns within REL_TOL, everything else exactly.
    """
    assert list(actual.columns) == list(expected.columns), f"{label} column names"
    assert actual.shape[0] == expected.shape[0], f"{label} row count"

    for col in expected.columns:
        a, e = actual[col].tolist(), expected[col].tolist()
        if pd.api.types.is_numeric_dtype(actual[col]):
            a, e = np.asarray(a, dtype=float), np.asarray(e, dtype=float)
            both_nan = np.isnan(a) & np.isnan(e)
            close = np.abs(a - e) <= REL_TOL * np.maximum(np.abs(a), np.abs(e))
            assert np.all(both_nan | close), (
                f"{label} column '{col}': actual={a} expected={e}"
            )
        else:
            a = [None if is_missing(v) else str(v) for v in a]
            e = [None if is_missing(v) else str(v) for v in e]
            assert a == e, f"{label} column '{col}'"


def assert_output_matches(out, case):
    for kind in ["stdout", "stderr"]:
        actual = drop_table_previews(normalize_output_lines(out[kind]))
        expected = drop_table_previews(read_expected_lines(case, kind))
        assert actual == expected, f"{case} {kind}"


MANIFEST = read_tsv(FIXTURES_DIR / "test_manifest.tsv")


def test_fixture_manifest_is_available():
    assert (FIXTURES_DIR / "test_manifest.tsv").exists()


@pytest.mark.parametrize(
    "row", MANIFEST.to_dict("records"), ids=MANIFEST["case"].tolist()
)
def test_fixture_case(row):
    out = run_case(build_args(row))
    case = row["case"]

    if bool(row["expect_error"]):
        expected_msg = "\n".join(read_expected_lines(case, "error"))
        assert out["error"] is not None, f"{case} raised an error"
        assert out["error"] == expected_msg, f"{case} error message"
    else:
        assert out["error"] is None, f"{case} unexpected error: {out['error']}"
        metadata = out["value"].metadata.copy()

        # The expected metadata stores a placeholder in place of the version
        assert metadata["sardine_version"].iloc[0] == VERSION
        metadata["sardine_version"] = VERSION_PLACEHOLDER

        assert_frame_matches(
            out["value"].results,
            read_tsv(EXPECTED_DIR / f"{case}_results.tsv"),
            f"{case} results",
        )
        assert_frame_matches(
            metadata,
            read_tsv(EXPECTED_DIR / f"{case}_metadata.tsv"),
            f"{case} metadata",
        )

    assert_output_matches(out, case)
