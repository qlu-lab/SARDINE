# SARDINE: Summary-statistics-based Association Regression for Deciphering INtegrated cpg Effects

This repository contains the source code for both the R and python `sardine-ewas` packages. Test fixtures ensure that the behavior matches between the R and python packages.

## Overview

Required inputs:

- EWAS summary statistics
- clock weights
- referance methylation variance-covariance matrix

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

- Var(Y): used in the computation of gamma, gamma_logistic, and log_hr. Also needed to compute lambda given a logistic ewas model.
  - the user can input the phenotypic variance, or it will be estimated using the reference methylation variance-covariance matrix
- EWAS DNA methylation scale and mean beta: handles whether EWAS was computed on DNA methylation in beta-value or M-value scale
  - default is beta value
  - mean beta value for each CpG is used to convert M-scale EWAS to beta-scale ewas
- scale and age_adult: used to re-scale Horvath and PedBE clocks
  - age_adult should be set iff scale is set to True
  - TODO: give a full description of this
- cox_pi and cox_E: used to compute kappa_cox
  - only used when EWAS model is Cox
  - describe the proportion (or number) of cases at the end of the study
  - exactly one must be provided if the EWAS model is Cox
  
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