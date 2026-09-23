# Basic netmatchRI workflow with a small deterministic example.

library(netmatchRI)

sim <- simulate_netmatch_example(seed = 90141, n = 60, beta_z = 2)

m <- netmatch(
  data = sim$data,
  treat = "Z",
  covariates = c("X1", "X2", "X3"),
  network = sim$net_dist,
  method = "dual",
  kappa = 2,
  solver = "highs"
)

match_summary <- summary(m)
diagnostics <- diagnose_match(m)
ri_unadjusted <- RI_unadjusted(m, "Y")
ri_adjusted <- RI_adjusted(m, "Y", eta = 0.03, rho = 0.10)
ri_design <- RI_design(m, "Y")

sens <- sensitivity_grid(
  m,
  "Y",
  eta = seq(0, 0.5, by = 0.1),
  rho = seq(0, 0.3, by = 0.05)
)
critical <- critical_sensitivity(m, "Y", rho = seq(0, 1, by = 0.05))
print(critical)
critical$summary
plot_sensitivity(sens)
plot_sensitivity(critical)
