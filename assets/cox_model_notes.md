# Cox model note

When the EWAS is a Cox model, whether `varY` is inferred or user-supplied determines whether the `kappa_cox`/`cox_pi` dependence cancels out or not. See the details below:

Tracing how `kappa_cox` propogates through the pipeline in `BETA` and `SE`:
   1. Converting cox sumstats to linear scale multiplies `BETA` and `SE` by `kappa_cox`.
   2. Converting linear Y~M to linear M~Y divides both by `denom = BETA^2 + (N-2)*SE(BETA)^2`, which scales as `kappa_cox^2`.
   => So the linear M~Y BETA values carry a factor of `1/kappa_cox`.
   => So `w_T_beta` carries a factor of `1/kappa_cox`

If `varY` is inferred, `kappa_cox` is part of the `varY` term:
   1. `varY = Var(M) / (BETA^2 + (N - 2) * SE(BETA)^2)`, so it scales as `kappa_cox^2` (same as `denom` above).

Tracing how `kappa_cox` propogates into `gamma_hat = w_T_beta * varY / w_T_Sigma_w`:
   1. If `varY` is inferred, `gamma_hat` carries a factor of `kappa_cox`
   2. If `varY` is user-supplied, `gamma_hat` carries a factor of `1/kappa_cox`

Tracing `kappa_cox` into the final log hazard ratio:
   1. If `varY` is inferred, `log_hr` is invariant to `kappa_cox` (and hence `cox_pi`).
   2. If `varY` is user-supplied, `log_hr` carries a factor of `1/kappa_cox^2` (so is sensitive to `cox_pi`)