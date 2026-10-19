# Python-specific input forms for the DNAm covariance matrix, which the shared
# fixtures (file inputs only) can't cover. Mirrors r/tests/testthat/test-inputs.R.


import pandas as pd
import pytest
from _fixtures_dir import FIXTURES_DIR
from sardine.io import read_dnam_covariance

from sardine import SardineError

SIGMA_TSV = FIXTURES_DIR / "inputs" / "sigma_baseline_synthetic.tsv"


def test_read_dnam_covariance_accepts_a_file_and_dataframes():
    from_file = read_dnam_covariance(str(SIGMA_TSV), verbose=False)
    assert list(from_file.index) == list(from_file.columns)

    # DataFrame with CpG IDs as its index and column names
    pd.testing.assert_frame_equal(
        read_dnam_covariance(from_file, verbose=False), from_file
    )

    # DataFrame laid out like the file, with CpG IDs in the first column
    with_id_col = from_file.reset_index(names="CpGmarker")
    pd.testing.assert_frame_equal(
        read_dnam_covariance(with_id_col, verbose=False), from_file
    )

    # pathlib.Path works as well as a string
    pd.testing.assert_frame_equal(
        read_dnam_covariance(SIGMA_TSV, verbose=False), from_file
    )


def test_read_dnam_covariance_rejects_mismatched_names_text_columns_and_other_types(
    tmp_path,
):
    mismatched = pd.DataFrame(
        [[1, 0.5], [0.5, 1]], index=["a", "b"], columns=["b", "a"]
    )
    with pytest.raises(SardineError, match="identical and in the same order"):
        read_dnam_covariance(mismatched, verbose=False)

    text_column = pd.DataFrame({"id": ["a", "b"], "a": ["x", "y"], "b": [1.0, 2.0]})
    with pytest.raises(SardineError, match="non-numeric covariance columns: a"):
        read_dnam_covariance(text_column, verbose=False)

    no_ids = pd.DataFrame([[1, 0.5], [0.5, 1]], columns=["a", "b"])
    with pytest.raises(SardineError, match="valid rownames and colnames"):
        read_dnam_covariance(no_ids, verbose=False)

    with pytest.raises(
        SardineError, match="must be a delimited file or a pandas DataFrame"
    ):
        read_dnam_covariance([1, 2, 3], verbose=False)

    with pytest.raises(SardineError, match="File not found"):
        read_dnam_covariance(str(tmp_path / "does_not_exist.tsv"), verbose=False)

    with pytest.raises(SardineError, match="can only be read by the R package"):
        read_dnam_covariance(str(tmp_path / "sigma.rds"), verbose=False)
