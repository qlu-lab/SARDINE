# Changelog

All notable changes to SARDINE are listed here. The R (`sardine`) and Python (`sardine-ewas`) packages share one version number and are released together. Versions follow [semantic versioning](https://semver.org).

## 0.1.0 (unreleased)

First public release.

- `sardine()` estimates the clock-phenotype association (`lambda`, `gamma`, `gamma_logistic` for binary phenotypes, and `log_HR` for Cox EWAS) from EWAS summary statistics, clock weights and a reference DNAm covariance matrix.
- Supports linear, logistic and Cox EWAS models, either regression direction, and beta-value or M-value methylation.
- Available as an R package (`sardine`) and a Python package (`sardine-ewas`, imported as `sardine`).
- 