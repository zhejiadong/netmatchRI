# netmatchRI 0.1.0

* Stabilizes matched-design objects with full-length subclasses, estimand-specific
  normalized unit weights, matched indicators, original data and call, resolved
  network type, and backend-neutral solver information. Complete backend output
  is available with `include_solver = TRUE`.
* Adds `matched_data()`, factor-aware before/after absolute SMD diagnostics,
  disconnected within-set pair reporting, and structured matched-design
  summaries with sample counts, effective sample sizes, and set sizes.
* Strengthens treatment, covariate, network, and randomization-inference outcome
  validation.
* Renames the sensitivity-grid function to `sensitivity_grid()`.
* Initial release.
* Implements dual-penalty matching for observational network data.
* Provides randomization-based inference, sensitivity analysis, and critical
  sensitivity curves for network dependence and interference.
* Uses the open-source HiGHS solver by default, with optional Gurobi and GLPK
  backends.
