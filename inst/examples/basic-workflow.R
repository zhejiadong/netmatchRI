# Basic netmatchRI workflow with a small deterministic example.

library(netmatchRI)

sim <- simulate_netmatch_example(seed = 20260821, n = 32)

match <- netmatch(
  data = sim$data,
  treat = "Z",
  covariates = c("X1", "X2", "X3"),
  network = sim$net_dist,
  method = "dual",
  kappa = 2,
  solver = "highs"
)

match_summary <- summary(match)
diagnostics <- diagnose_match(match)
ri_naive <- RI_Naive(match, "Y")
ri_sensitivity <- RI_Sensitivity(match, "Y", eta = 0.03, rho = 0.10)
ri_design <- RI_Design(match, "Y")

grid <- sensitivity_grid(
  match,
  "Y",
  eta = c(0, 0.03),
  rho = c(0, 0.5, 1)
)
critical <- critical_sensitivity(match, "Y", rho = c(0, 0.5, 1))
