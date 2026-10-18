# SARDINE: Summary-statistics-based Association Regression for Deciphering INtegrated cpg Effects

This repository contains the source code for both the R and python `sardine-ewas` packages. Test fixtures ensure that the behavior matches between the R and python packages.

## Algorithm overview

SARDINE computes the association between epigenetic clocks and phenotypes using summary statistics from EWAS studies, clock weights, and a reference DNAm covariance matrix. See the accompanying manuscript [here](TODO: put link to preprint).

### Inputs

Required inputs:

- EWAS summary statistics
- clock weights
- referance methylation variance-covariance matrix
  - TODO: link our data once we share it

Can handle various types of EWAS summary statistics, based on the following EWAS-generating models:

- linear model (linear regression of phenotype on methylation or methylation on phenotype)
  - can handle either direction of regression (Y~M and M~Y)
  - can handle binary or continuous phenotype
- logistic model (logistic regression of binary phenotype on methylation)
  - for binary phenotypes only
  - direction of regression must be Y~M
- Cox model (regression of time to binary event on methylation)
  - for binary phenotypes (incidence by study end) only
  - direction of regression must be Y~M
  
We then output the following regression estimates (in addition to metadata regarding Var(Y) and the coverage of CpGs):

- linear regression of Y~C (gamma)
- linear regression of C~Y (lambda)
- logistic regression of Y~C approximated from gamma (gamma_logistic)
  - this is for binary phenotypes only
  - this is the log-odds ratio
  - note: for a cox ewas model, this is the log-odds ratio for cumulative incidence at end of follow-up
- log-hazard-ratio (log_hr)
  - this is for a cox ewas model (thus binary phenotype) only

Optional inputs:

- Var(Y): used in the computation of gamma, gamma_logistic, and log_hr. Also needed to compute lambda given a logistic ewas model. Note: this should be the Var(Y) transformed for the same covariates as the EWAS
  - the user can input the phenotypic variance, or it will be estimated using the EWAS summary stats and the reference methylation variance-covariance matrix
- EWAS DNA methylation scale and mean beta: handles whether EWAS was computed on DNA methylation in beta-value or M-value scale
  - default is beta value
  - mean beta value for each CpG is used to convert M-scale EWAS to beta-scale ewas
- scale and age_adult: used to re-scale Horvath and PedBE clocks because they apply a scale factor to their clocks.
  - age_adult should be set iff scale is set to True
  - age_adult is the scale factor used
- cox_pi and cox_E: used to compute kappa_cox
  - only used when EWAS model is Cox
  - describe the proportion (or number) of cases at the end of the study
  - exactly one must be provided if the EWAS model is Cox

### Flow chart

How `sardine()` routes inputs to final estimates.

![SARDINE workflow](https://raw.githubusercontent.com/org/qlu-lab/main/assets/sardine_workflow.png)

### Outputs

Data frame of coefficient estimates (`lambda` and `gamma` are always estimated; `gamma_logistic` is estimated for binary phenotypes and `log_HR` is estimated for Cox EWAS), with columns:

- model: EWAS model used (`e.g. Phenotype ~ Clock (Cox)`)
- param: parameter being estimated (e.g. `log_HR`)
- coef: point estimate for the parameter
- SE: estimated standard error for the point estimate
- z: z-statistic (estimate / SE). All rows share the z of the linear `gamma` test
- p: p-value for the estimate, testing the null hypothesis that the effect is 0
- neglog10p: -log10(p) (provided due to underflow with very small p-values)
- N: EWAS sample size: the `N` argument if supplied, otherwise the median per-CpG N across the CpGs used.
- note: any notes about the specific settings for this row (e.g. log-odds vs log_HR), or `NA`.)

Metadata with columns:

- N_CpGs_ewas, N_CpGs_clock, N_CpGs_DNAmRef: number of CpGs in the EWAS, clock, and DNAm reference dataset respectively
- N_CpGs_ewas_clock, N_CpGs_ewas_ref, N_CpGs_clock_ref: number of CpGs overlapping between each pair of inputs
- N_CpGs_final: number of CpGs in all three inputs; these are the CpGs actually used for the association analysis
- ClockCpG_coverage_ratio: fraction of the clock's CpGs that ended up in the final set of CpGs used
- ClockCpG_missed_ratio: 1 - ClockCpG_coverage_ratio
- cox_pi (`NA` if EWAS model is not cox): proportion of cases at the end of the study. This was either input directly or computed from cox_E input
- kappa_cox (`NA` if EWAS model is not cox): Cox conversion factor used to convert from the linear scale gamma to the Cox log-hazard-ratio
- varY: variance of the phenotype (after any capping). This is either input directly or estimated using the EWAS summary statistics and the DNAm reference matrix
- varY_type: `User_provided` or `Inferred`
- varY_ratio: (varY prior to capping)/(theoretical max varY)
  - The theoretical maximum is the varY you'd get if the clock explained all of the phenotype's variance
- varY_note: Note about whether varY had to be capped at the maximum theoretical phenotypic variance

### Additional notes

- See [this page](https://raw.githubusercontent.com/org/qlu-lab/main/assets/cox_model_notes.png) for additional notes about using SARDINE with a Cox EWAS model.

## Development infrastructure and environment management

Set up once with `pixi install` (or let `pixi run` do it lazily on first use). Three environments are defined in `pixi.toml`:

- `dev` — everything: Python + R runtime deps, plus dev tooling for both (`pytest`, `hatchling`, `python-build`, `twine`; `devtools`, `usethis`, `testthat`, `roxygen2`, `covr`, `pkgdown`).
- `user-python` — only the Python package's runtime deps (`numpy`, `pandas`, `scipy`), plus `pip` for installing the built wheel. Approximates a real user's environment; use it to sanity-check that a built wheel/sdist actually installs and runs with nothing but its declared dependencies.
- `user-r` — only the R package's runtime dep (`data.table`). Same idea: sanity-check the built source tarball against a clean, dev-tool-free R install.

Run any command inside an environment with `pixi run -e <env> <command>`, or drop into an interactive shell with `pixi shell -e <env>`.

TODO: update the pixi.toml and these notes as packages are added.

### Python

- Run the test suite: `pixi run -e dev test-py`
- Build the package: `pixi run -e dev build-py` (produces a wheel + sdist under `python/dist/`)
- Sanity-check a clean install: `pixi run -e user-python python -m pip install python/dist/*.whl` then exercise the package from that environment

### R

- Regenerate/re-check golden fixtures: `pixi run -e dev generate-fixtures` / `pixi run -e dev run-test-fixtures`
- Document: `pixi run -e dev document-r`
- Run the test suite: `pixi run -e dev test-r`
- Full package check: `pixi run -e dev check-r`
- Sanity-check a clean install: build a source tarball (`pixi run -e dev Rscript -e 'devtools::build("r")'`), then `pixi run -e user-r Rscript -e 'install.packages("path/to/sardine_x.y.z.tar.gz", repos = NULL, type = "source")'`
- Perform CRAN's pre-submission checks: `devtools::check_win_devel()`