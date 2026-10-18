# SARDINE: Summary-statistics-based Association Regression for Deciphering INtegrated cpg Effects

This repository contains the source code for both the R package (`sardine`) and the Python package (`sardine-ewas`). Shared test fixtures ensure that the two packages behave the same, and both are released together under the same version number.

TODO: Python version not implemented yet, but placeholders are below in the README.

## Installation

TODO: This won't work until published publicly.

```r
# R
install.packages("sardine")
```

```sh
# Python
pip install sardine-ewas
```

To check that the installation works, run SARDINE on the small sample dataset bundled with the package:

```r
# R
library(sardine)
example(sardine)
```

```sh
# Python
# TODO
```

## Algorithm overview

SARDINE computes the association between epigenetic clocks and phenotypes using summary statistics from EWAS studies, clock weights, and a reference DNAm covariance matrix. In other words, SARDINE takes per-CpG effect estimates from EWAS summary statistics, combines them using clock weights and a covariance matrix, and outputs the effect estimate for the clock as a whole. See the accompanying manuscript [here](TODO: put link to preprint).

### Inputs

See [inputs.md](https://github.com/qlu-lab/SARDINE/blob/main/inputs.md) for full details on every input, including file formats and accepted column names.

| Argument | Required | Default | Accepted values | Description |
| --- | --- | --- | --- | --- |
| `ewas_sumstats` | yes | | file path or table ([format](https://github.com/qlu-lab/SARDINE/blob/main/inputs.md#ewas-summary-statistics-ewas_sumstats)) | EWAS summary statistics: one effect estimate, SE and N per CpG |
| `clock` | yes | | file path or table ([format](https://github.com/qlu-lab/SARDINE/blob/main/inputs.md#clock-weights-clock)) | Clock weights: one weight per CpG |
| `DNAm_cov_matrix` | yes | | file path, matrix or table ([format](https://github.com/qlu-lab/SARDINE/blob/main/inputs.md#dnam-covariance-matrix-dnam_cov_matrix)) | Reference DNAm variance-covariance matrix: use one of the reference matrices shared with SARDINE (TODO: link our data once we share it), or compute your own from a reference population |
| `ewas_direction` | no | `"M~Y"` | `"M~Y"`, `"Y~M"` | Direction of the EWAS regression (outcome ~ exposure) |
| `ewas_model` | no | `"linear"` | `"linear"`, `"logistic"`, `"cox"` | Regression model the EWAS used |
| `ewas_DNAm_scale` | no | `"beta"` | `"beta"`, `"m"` | Methylation scale of the EWAS: beta values or M-values |
| `ewas_pheno_type` | no | `"continuous"` | `"continuous"`, `"binary"` | Phenotype type; binary also produces `gamma_logistic` |
| `mean_beta` | only if `ewas_DNAm_scale = "m"` | none | file path or table ([format](https://github.com/qlu-lab/SARDINE/blob/main/inputs.md#mean-beta-values-mean_beta)) | Mean beta value per CpG, used to convert M-value sumstats to the beta scale |
| `varY` | no | inferred | non-negative number | Phenotype variance, with the same covariate adjustment as the EWAS; estimated from the sumstats and covariance matrix if omitted |
| `N` | no | median of the sumstats' `N` column | number | EWAS sample size |
| `cox_pi` | for Cox, exactly one of `cox_pi` / `cox_E` | none | number strictly between 0 and 1 | Proportion of cases at the end of the study |
| `cox_E` | for Cox, exactly one of `cox_pi` / `cox_E` | none | number $> 0$ | Number of cases at the end of the study |
| `scale` | no | `FALSE` | `TRUE` / `FALSE` | Rescale estimates for the Horvath and PedBE clocks |
| `age_adult` | no | `20` | number | Adult age used in that rescaling |
| `verbose` | no | `TRUE` | `TRUE` / `FALSE` | Print progress and diagnostic messages |

The EWAS model, direction and phenotype type must form one of these combinations:

| `ewas_model` | `ewas_direction` | `ewas_pheno_type` | Also required |
| --- | --- | --- | --- |
| `"linear"` | `"M~Y"` or `"Y~M"` | `"continuous"` or `"binary"` | |
| `"logistic"` | `"Y~M"` | `"binary"` | |
| `"cox"` | `"Y~M"` | `"binary"` (incidence by study end) | exactly one of `cox_pi` / `cox_E` |

- `ewas_DNAm_scale = "m"` works with any of these and also requires `mean_beta`.

### Flow chart

How `sardine()` routes inputs to final estimates.

![SARDINE workflow](https://raw.githubusercontent.com/qlu-lab/SARDINE/main/assets/sardine_workflow.png)

### Outputs

`sardine()` returns a table with one row per estimated parameter:

| `param` | `model` | Estimated when | Meaning |
| --- | --- | --- | --- |
| `lambda` | `Clock ~ Phenotype (linear)` | always | Linear regression of the clock on the phenotype |
| `gamma` | `Phenotype ~ Clock (linear)` | always | Linear regression of the phenotype on the clock |
| `gamma_logistic` | `Phenotype ~ Clock (logistic)` | binary phenotype | Log-odds ratio of the phenotype per unit of the clock, approximated from `gamma`. For a Cox EWAS, this is the log-odds ratio for cumulative incidence at the end of follow-up |
| `log_HR` | `Phenotype ~ Clock (Cox)` | Cox EWAS | Log hazard ratio per unit of the clock |

Its columns:

| Column | Description |
| --- | --- |
| `model` | EWAS model used (e.g. `Phenotype ~ Clock (Cox)`) |
| `param` | Parameter being estimated (e.g. `log_HR`) |
| `coef` | Point estimate for the parameter |
| `SE` | Estimated standard error for the point estimate |
| `z` | z-statistic (estimate / SE). All rows share the z of the linear `gamma` test |
| `p` | p-value for the estimate, testing the null hypothesis that the effect is 0 |
| `neglog10p` | -log10(p) (provided due to underflow with very small p-values) |
| `N` | EWAS sample size: the `N` argument if supplied, otherwise the median per-CpG N across the CpGs used |
| `note` | Any notes about the specific settings for this row (e.g. log-odds vs log_HR), or `NA` |

The result also carries a one-row metadata table, `attr(result, "metadata")` in R, with columns:

| Column | Description |
| --- | --- |
| `sardine_version` | Version of SARDINE that produced the result |
| `N_CpGs_ewas`, `N_CpGs_clock`, `N_CpGs_DNAmRef` | Number of CpGs in the EWAS, clock, and DNAm reference dataset respectively |
| `N_CpGs_ewas_clock`, `N_CpGs_ewas_ref`, `N_CpGs_clock_ref` | Number of CpGs overlapping between each pair of inputs |
| `N_CpGs_final` | Number of CpGs in all three inputs; these are the CpGs actually used for the association analysis |
| `ClockCpG_coverage_ratio` | Fraction of the clock's CpGs that ended up in the final set of CpGs used |
| `ClockCpG_missed_ratio` | 1 - `ClockCpG_coverage_ratio` |
| `cox_pi` | Proportion of cases at the end of the study, either input directly or computed from `cox_E`. `NA` if the EWAS model is not Cox |
| `kappa_cox` | Cox conversion factor used to convert from the linear scale gamma to the Cox log hazard ratio. `NA` if the EWAS model is not Cox |
| `varY` | Variance of the phenotype (after any capping), either input directly or estimated using the EWAS summary statistics and the DNAm reference matrix |
| `varY_type` | `user_provided` or `inferred` |
| `varY_ratio` | (varY prior to capping) / (theoretical max varY). The theoretical maximum is the varY you'd get if the clock explained all of the phenotype's variance |
| `varY_note` | Note about whether varY had to be capped at the maximum theoretical phenotypic variance |

### Additional notes

- See [this page](https://github.com/qlu-lab/SARDINE/blob/main/assets/cox_model_notes.md) for additional notes about using SARDINE with a Cox EWAS model.

## Development infrastructure and environment management

Set up once with `pixi install` (or let `pixi run` do it lazily on first use). Three environments are defined in `pixi.toml`:

- `dev` — everything: Python + R runtime deps, plus dev tooling for both (`pytest`, `hatchling`, `python-build`, `twine`; `devtools`, `usethis`, `testthat`, `roxygen2`, `covr`, `pkgdown`).
- `user-python` — only the Python package's runtime deps (`numpy`, `pandas`, `scipy`), plus `pip` for installing the built wheel. Approximates a real user's environment; use it to sanity-check that a built wheel/sdist actually installs and runs with nothing but its declared dependencies.
- `user-r` — only the R package's runtime dep (`data.table`). Same idea: sanity-check the built source tarball against a clean, dev-tool-free R install.

Run any command inside an environment with `pixi run -e <env> <command>`, or drop into an interactive shell with `pixi shell -e <env>`.

TODO: update the pixi.toml and these notes as packages are added.

### Running the tests

| Command | What it checks |
| --- | --- |
| `pixi run -e dev test-r` | Every R test: all shared fixture cases plus the R-only tests (e.g. covariance matrix given as an R object) |
| `pixi run -e dev check-r` | `R CMD check`: builds the package, runs `test-r` and the `example(sardine)` example, and checks the docs and package structure, as CRAN does |
| `pixi run -e dev test-py` | Every Python test, including the same shared fixture cases |
| `pixi run -e dev check-sync` | That `VERSION` and `fixtures/` match their copies in the packages |

`test-r`, `check-r` and `build-r` all run `sync` first, so they always use the current `VERSION` and `fixtures/`.

### Shared files

The repo root holds the files both packages share; `pixi run -e dev sync` copies them into the packages, and `pixi run -e dev check-sync` fails if any copy is out of date:

- `VERSION`: the single version number for both packages, copied into `r/DESCRIPTION` and `python/pyproject.toml`
- `fixtures/`: test inputs and expected outputs, copied into `r/inst/extdata/fixtures/` because `R CMD check` cannot see files outside `r/`. Edit fixtures in `fixtures/` only, never the copy in `r/inst/extdata/fixtures/`

Fixtures:

- Each fixture case is checked for its results, metadata, error message and printed stdout/stderr (`r/tests/testthat/test-fixtures.R`).
- Regenerate fixtures: `pixi run -e dev generate-fixtures`. Only do this when an expected output should change (e.g. a bug fix), and review the diff
- The version is stored as `<VERSION>` in the expected outputs so bumping the version never requires regenerating them, and time stamps are similarly stored as `<TIMESTAMP>`

### Python

- Run the test suite: `pixi run -e dev test-py`
- Build the package: `pixi run -e dev build-py` (produces a wheel + sdist under `python/dist/`)
- Sanity-check a clean install in the minimal `user-python` environment, which has only the runtime dependencies:
  1. `pixi run -e user-python python -m pip install python/dist/*.whl`
  2. Run the example in the same environment (TODO: add an `example-py` task once the Python package exists)

### R

The R package lives in `r/` (code in `r/R/`, tests in `r/tests/testthat/`). `NAMESPACE` and `man/` are generated by roxygen; don't edit them by hand.

- Document: `pixi run -e dev document-r`
- Run all R tests: `pixi run -e dev test-r` (see [Running the tests](#running-the-tests))
- Full package check: `pixi run -e dev check-r`
- Sanity-check a clean install in the minimal `user-r` environment, which has only the runtime dependencies:
  1. `pixi run -e dev build-r`: builds `dist/sardine_<VERSION>.tar.gz`
  2. `pixi run -e user-r install-r`: installs that tarball
  3. `pixi run -e user-r example-r`: runs `example(sardine)` on the bundled sample data, as a user would
- Check on CRAN's Windows development build (results arrive by email): `pixi run -e dev Rscript -e 'devtools::check_win_devel("r")'`

### Publishing an update

The R and Python packages are always released together, with the same version number, and only once both pass the same fixtures.

1. **Make and test the change in both packages.** Any change to behavior, inputs or outputs goes into both the R and Python code. If the change should alter an expected output, update `fixtures/` (manually or via `generate-fixtures` or a new case) and check that both packages match it.
2. **Bump `VERSION`** at the repo root, following [semantic versioning](https://semver.org): patch (`0.1.1`) for bug fixes, minor (`0.2.0`) for new features, major (`1.0.0`) for changes that break existing code. Use plain `X.Y.Z`, which is valid for both CRAN and PyPI. Never edit the version in `r/DESCRIPTION` or `python/pyproject.toml` directly.
3. **Add an entry to [CHANGELOG.md](https://github.com/qlu-lab/SARDINE/blob/main/CHANGELOG.md)** under the new version.
4. **Sync and regenerate:**
   - `pixi run -e dev sync`
   - `pixi run -e dev document-r`
5. **Run every check.** All must pass before releasing:
   - `pixi run -e dev check-sync`
   - `pixi run -e dev check-r`: runs all R tests and the example, plus CRAN's package checks; aim for 0 errors, 0 warnings, 0 notes
   - `devtools::check_win_devel("r")` for CRAN's other platforms (results arrive by email)
   - `pixi run -e dev test-py`
   - Install both built packages in the clean `user-r` and `user-python` environments and run the examples (see the R and Python sections above)
6. **Commit, tag and push:**
   - `git commit -am "Release vX.Y.Z"`
   - `git tag vX.Y.Z`
   - `git push && git push --tags`
7. **Publish both packages**, within the same day where possible:
   - R: from `r/`, run `devtools::submit_cran()`, then confirm the submission from the email CRAN sends the maintainer. CRAN review can take a few days, and they may ask for changes; if so, make them, bump the patch version, and repeat from step 2 for both packages.
   - Python: `pixi run -e dev build-py`, then `twine upload python/dist/*`
   - GitHub: create a release from the `vX.Y.Z` tag, with the CHANGELOG entry as its notes.
