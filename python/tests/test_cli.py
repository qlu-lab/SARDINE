# The `sardine` command-line interface: runs fixture cases through the CLI and
# checks the TSV files it writes against the same expected outputs as test_fixtures.py.

import os
import subprocess
import sys
from pathlib import Path

import pytest
from _fixtures_dir import FIXTURES_DIR
from sardine.cli import main
from test_fixtures import VERSION_PLACEHOLDER, assert_frame_matches, read_tsv

import sardine

EXPECTED_DIR = FIXTURES_DIR / "expected"
INPUTS = FIXTURES_DIR / "inputs"


def run_cli(tmp_path, *extra_args):
    out = tmp_path / "out" / "result"
    exit_code = main(
        [
            "--ewas_sumstats",
            str(INPUTS / "ewas_sumstats_baseline_synthetic.tsv"),
            "--clock",
            str(INPUTS / "clock_weights_baseline_synthetic.tsv"),
            "--DNAm_cov_matrix",
            str(INPUTS / "sigma_baseline_synthetic.tsv"),
            "--out",
            str(out),
            *extra_args,
        ]
    )
    return exit_code, out


@pytest.mark.parametrize(
    "case, extra_args",
    [
        ("baseline_linear", []),
        (
            "cox_pi",
            [
                "--ewas_direction",
                "Y~M",
                "--ewas_model",
                "cox",
                "--ewas_pheno_type",
                "binary",
                "--cox_pi",
                "0.08",
            ],
        ),
        ("horvath_scaling", ["--scale", "--age_adult", "20"]),
    ],
)
def test_cli_writes_the_expected_results(tmp_path, case, extra_args):
    exit_code, out = run_cli(tmp_path, *extra_args)
    assert exit_code == 0

    results = read_tsv(f"{out}_results.tsv")
    metadata = read_tsv(f"{out}_metadata.tsv")
    assert metadata["sardine_version"].iloc[0] == sardine.__version__
    metadata["sardine_version"] = VERSION_PLACEHOLDER

    assert_frame_matches(
        results, read_tsv(EXPECTED_DIR / f"{case}_results.tsv"), f"{case} results"
    )
    assert_frame_matches(
        metadata, read_tsv(EXPECTED_DIR / f"{case}_metadata.tsv"), f"{case} metadata"
    )


def test_cli_reports_invalid_inputs_without_a_traceback(tmp_path, capsys):
    exit_code, out = run_cli(
        tmp_path,
        "--ewas_model",
        "cox",
        "--ewas_pheno_type",
        "binary",
        "--cox_pi",
        "0.08",
    )
    assert exit_code == 1
    assert (
        capsys.readouterr()
        .err.strip()
        .endswith("[ERROR] 'ewas_direction' must be 'Y~M' when 'ewas_model' is 'cox'.")
    )
    assert not (tmp_path / "out").exists()


def test_cli_quiet_prints_only_warnings(tmp_path):
    # Run as a real process: pytest captures warnings itself, so they would not reach stderr here
    env = {**os.environ, "PYTHONPATH": str(Path(sardine.__file__).resolve().parents[1])}
    completed = subprocess.run(
        [
            sys.executable,
            "-m",
            "sardine",
            "--ewas_sumstats",
            str(INPUTS / "ewas_sumstats_baseline_synthetic.tsv"),
            "--clock",
            str(INPUTS / "clock_weights_baseline_synthetic.tsv"),
            "--DNAm_cov_matrix",
            str(INPUTS / "sigma_baseline_synthetic.tsv"),
            "--out",
            str(tmp_path / "result"),
            "--varY",
            "6233.6",
            "--quiet",
        ],
        capture_output=True,
        text=True,
        env=env,
    )
    assert completed.returncode == 0
    assert completed.stdout == ""
    stderr_lines = [
        line for line in completed.stderr.splitlines() if not line.startswith("OMP:")
    ]
    assert stderr_lines == [
        "Warning: [WARN] User-supplied varY must be the variance of the phenotype after adjusting for "
        "the same covariates as the EWAS (the raw phenotype variance only if the EWAS had no covariates). "
        "Leave varY unset to estimate it from the EWAS summary statistics instead."
    ]


def test_cli_version(capsys):
    with pytest.raises(SystemExit) as exit_info:
        main(["--version"])
    assert exit_info.value.code == 0
    assert capsys.readouterr().out.strip() == f"sardine {sardine.__version__}"


def test_cli_example_runs_without_other_arguments(capsys):
    with pytest.raises(SystemExit) as exit_info:
        main(["--example"])
    assert exit_info.value.code == 0
    out = capsys.readouterr().out
    assert out.startswith("SARDINE result: 320 CpGs used (97.0% of the clock)")
    assert "lambda -0.006648" in out


def test_cli_requires_the_three_inputs_and_out(capsys):
    with pytest.raises(SystemExit) as exit_info:
        main(["--clock", "clock.tsv"])
    assert exit_info.value.code == 2
    assert "required" in capsys.readouterr().err


def test_output_paths_for_a_prefix_or_a_folder(tmp_path):
    from sardine.cli import output_paths

    # A prefix: <prefix>_results.tsv and <prefix>_metadata.tsv
    assert output_paths("results/my_trait") == (
        "results/my_trait_results.tsv",
        "results/my_trait_metadata.tsv",
    )

    # A folder, written with a trailing slash (whether or not it exists yet)
    assert output_paths("results/my_trait/") == (
        "results/my_trait/results.tsv",
        "results/my_trait/metadata.tsv",
    )

    # An existing folder, without the trailing slash
    existing = str(tmp_path)
    assert output_paths(existing) == (
        f"{existing}/results.tsv",
        f"{existing}/metadata.tsv",
    )


def test_cli_out_folder_writes_results_tsv_inside_it(tmp_path):
    out_dir = tmp_path / "new_folder"
    exit_code = main(
        [
            "--ewas_sumstats",
            str(INPUTS / "ewas_sumstats_baseline_synthetic.tsv"),
            "--clock",
            str(INPUTS / "clock_weights_baseline_synthetic.tsv"),
            "--DNAm_cov_matrix",
            str(INPUTS / "sigma_baseline_synthetic.tsv"),
            "--out",
            f"{out_dir}/",
            "--quiet",
        ]
    )
    assert exit_code == 0
    assert sorted(p.name for p in out_dir.iterdir()) == ["metadata.tsv", "results.tsv"]
