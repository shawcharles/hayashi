# panel_fe_twoway_cluster_wagepan

This case is blocked because its raw two-way clustered covariance is materially indefinite. Fitted covariance inference is unavailable under the present full-matrix contract; rejection is not successful estimator inference.

## Model

```hayashi
xtset(df, nr, year)
let m = fe(lwage ~ union + married + d81 + d82 + d83 + d84 + d85 + d86 + d87, df, cluster=nr, cluster2=year)
```

Both references retain within-transformed OLS on the Wooldridge `wagepan` dataset and the Cameron-Gelbach-Miller decomposition:

```
V = V_entity + V_time - V_intersection
```

Entity, time and entity/time intersection score matrices retain the common small-sample correction `g/(g-1) * (n-1)/(n-k)`, where `g = min(G_entity, G_time)`. No diagonal clipping, eigenvalue projection or change to this covariance estimator is applied.

## Reference diagnostics

R and Python emit coefficients, the raw covariance as nested row/column names, `inference_available`, `inference_reason`, and `covariance_diagnostics`. The diagnostics preserve eigenvalues in ascending order, identify negative variance terms, and report numerical bounds. Standard errors are emitted only when covariance inference is admissible; the `standard_errors` field is absent otherwise. Non-finite raw entries are labelled as JSON strings rather than replaced by zero.

Admissibility requires finite entries, symmetry within `128*k*eps*max(abs(V))`, no negative diagonal variance, and no eigenvalue below `-128*k*eps*max(abs(eigenvalues))`. Here `eps` is binary64 machine epsilon and `k` is the covariance dimension. Eigenvalues are computed from the symmetric part solely for diagnosis; raw covariance entries remain unchanged. These declared scale-relative bounds address floating-point roundoff, not substantive negative variance.

On the reviewed 4,360-row snapshot, `d81` variance is approximately `-7.162277e-6`; four eigenvalues are negative, with minimum approximately `-6.862758e-5`. These values are materially below the roundoff bounds, so the references report unavailable inference and omit standard errors. A clipped `d81` standard error of zero would falsely imply exact precision for its nonzero coefficient.

## Qualification boundary

The declared coefficient and standard-error tolerances remain `1e-4`, but the blocked case cannot qualify either a fitted result or its uncertainty. Retaining raw diagnostic output does not promote the case to a numerical pass. Reconsidering inference would require an explicit, separately validated covariance policy rather than wider tolerances or an implicit repair.

The deterministic interpreter regression in `tests/panel_covariance_contract.rs` establishes finite inputs, a full-rank within design, adequate observations and valid two-way group counts before requiring covariance-specific rejection of a materially negative variance. That regression checks rejection behaviour; it does not replace numerical or inferential qualification on `wagepan`.
