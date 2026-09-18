# SARDINE: Summary-statistics-based Association Regression for Deciphering INtegrated cpg Effects

This repository contains the source code for both the R and python `sardine-ewas` packages. Test fixtures ensure that the behavior matches between the R and python packages.

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