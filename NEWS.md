# netmatchRI 0.1.0

## Initial release

* Provides dual-penalty, single-penalty, and covariate-only matching for
  observational network data.
* Provides matched-design diagnostics, comparison plots, and randomization-based
  inference under residual network dependence.
* Includes sensitivity grids, critical sensitivity curves, critical eta at
  rho = 1, the critical ratio, and simulated example network data.
* Distinguishes critical boundaries from cases with no crossing in the
  specified parameter range, and validates supplied network distances.
* Evaluates the design-based covariance bound by an exact finite-support sum rather
  than a quantile grid, using uniform marginals for full matching.
* Removes outcome-ties rejection without changing the strict rank comparison
  or the no-ties assumption underlying the null moments.
* Retains ATT, ATC, and ATE unit matching weights for diagnostics and
  matched-data summaries; sharp-null inference uses separate matched-set weights.
* Uses the open-source HiGHS solver by default, with optional Gurobi and GLPK
  backends.
