# SARDINE inputs

This page describes every input to SARDINE: what it is, which values it accepts, and the file formats for tabular inputs. It applies to both the R and Python packages, which accept the same values and file formats.

For a one-table overview of every input, its default and its accepted values, and for the valid combinations of EWAS model, direction and phenotype type, see the [Inputs section of the README](https://github.com/qlu-lab/SARDINE#inputs).

## Option values

`ewas_direction`, `ewas_model`, `ewas_DNAm_scale` and `ewas_pheno_type` each take one of a fixed set of values. Matching ignores case and whitespace.

- **`ewas_direction`**: the direction of the EWAS regression, written as outcome ~ exposure.
  - `"M~Y"`: methylation regressed on the phenotype.
  - `"Y~M"`: the phenotype regressed on methylation.
- **`ewas_model`**: the regression model the EWAS used.
  - `"linear"`: linear regression. Works with either direction and either phenotype type.
  - `"logistic"`: logistic regression of a binary phenotype on methylation.
  - `"cox"`: Cox proportional hazards regression of time to a binary event on methylation.
- **`ewas_DNAm_scale`**: the scale of the methylation values in the EWAS.
  - `"beta"`: beta values.
  - `"m"`: M-values. The summary statistics are converted to the beta scale, which requires `mean_beta`.
- **`ewas_pheno_type`**: the type of phenotype.
  - `"continuous"`.
  - `"binary"`: also produces `gamma_logistic`, a log-odds ratio of the phenotype on the clock.

Not every combination of these is valid; for example, a logistic or Cox EWAS must be `"Y~M"` with a binary phenotype. See the [valid combinations](https://github.com/qlu-lab/SARDINE#inputs) in the README.

## Tabular inputs

Every table can be given as a path to a delimited file (TSV, CSV or anything else `data.table::fread()` detects), or in R as a `data.frame` or `data.table`. Column names are matched case-insensitively against the names listed for each column.

### EWAS summary statistics (`ewas_sumstats`)

One row per CpG.

| Column | Accepted names | Notes |
| --- | --- | --- |
| CpG ID | `CpGmarker`, `CpG`, `marker`, `probe`, `ID` | |
| Effect estimate | `BETA`, `coef`, `coeff`, `estimate` | Rows with a missing value are dropped. |
| Standard error | `SE`, `std.error`, `std_error`, `standard_error` | Rows with a missing value are dropped. |
| Sample size | `N`, `samplesize`, `n_samples` | Optional when the `N` argument is given, which overrides it. |

Other columns are ignored.

`BETA` and `SE` are each CpG's own coefficient and standard error, exactly as the EWAS reported them, on the scale of the EWAS model:

| `ewas_model` | `BETA` for each CpG |
| --- | --- |
| `"linear"` | the linear regression slope |
| `"logistic"` | the log-odds ratio |
| `"cox"` | the log hazard ratio |

### Clock weights (`clock`)

One row per CpG.

| Column | Accepted names |
| --- | --- |
| CpG ID | `CpGmarker`, `CpG`, `marker`, `probe`, `ID`, `pred.var` |
| Weight | `CoefficientTraining`, `BETA`, `coef`, `coeff`, `estimate`, `weight` |

An intercept row (a CpG ID of `(Intercept)` or `intercept`, in any case) is dropped.

### DNAm covariance matrix (`DNAm_cov_matrix`)

The CpG-by-CpG covariance matrix of DNA methylation in a reference population. It should be computed with the same covariate adjustment as the EWAS. You can use the reference matrices shared with SARDINE (TODO: link once shared) or compute your own.

As a delimited file, the first column holds the CpG IDs and the header lists the same CpG IDs in the same order:

```
CpGmarker	cg00000001	cg00000002	cg00000003
cg00000001	0.0123	0.0011	-0.0004
cg00000002	0.0011	0.0098	0.0007
cg00000003	-0.0004	0.0007	0.0145
```

The name of the first column doesn't matter.

In R, these forms also work:

- a path to an `.rds` file containing a matrix;
- a matrix with CpG IDs as its row and column names;
- a `data.frame` of numeric columns with CpG IDs as its row names;
- a `data.frame` or `data.table` with CpG IDs in the first column, laid out like the tsv file.

Row and column names must be identical and in the same order.

### Mean beta values (`mean_beta`)

Only used when `ewas_DNAm_scale = "m"`. One row per CpG, giving the mean methylation beta value, ideally from the same samples as the EWAS.

| Column | Accepted names |
| --- | --- |
| CpG ID | `CpGmarker`, `CpG`, `marker`, `probe`, `ID` |
| Mean beta value | `beta_bar`, `mean_beta`, `beta_mean`, `mean_dnam_beta`, `avg_beta` |

## Which CpGs are used

SARDINE uses only the CpGs present in all three of the summary statistics, the clock weights and the covariance matrix. The output metadata reports statistics about how many CpGs were retained.

## Other numeric inputs

- **`varY`**: the variance of the phenotype. It should have the same covariate adjustment as the EWAS: the raw phenotype's variance if the EWAS had no covariates, otherwise the variance of the phenotype after regressing out those covariates. If omitted, the phenotypic variance is estimated from the summary statistics and the covariance matrix. Either way, it is capped at the largest value consistent with the inputs (the variance you'd get if the clock explained all of the phenotype's variance), and the metadata reports whether phenotypic variance capping occurred. For Cox models, leaving `varY` unset is recommended; see the [Cox model notes](assets/cox_model_notes.md).
- **`N`**: the EWAS sample size. Overrides the summary statistics' `N` column. For Cox models, this is the total sample size at the end of the study.
- **`cox_pi`**: the proportion of cases at the end of the study (Cox only).
- **`cox_E`**: the number of cases at the end of the study (Cox only). `cox_pi` is then computed as `cox_E / N`.
- **`scale`** and **`age_adult`**: for the Horvath and PedBE clocks, which predict a transformed age. With `scale = TRUE`, the clock-scale estimates are divided by `1 + age_adult` (and `lambda` is multiplied by it). This changes the estimates and standard errors, but not the z-statistics or p-values.
