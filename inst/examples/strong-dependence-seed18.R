# 300-unit example using the default simulation setup.
# This mirrors the package illustration with i = 18.
# The implied seed is 90141.

library(netmatchRI)

sim <- simulate_netmatch_example(seed = 90141, beta_z = 0)
sim_dat <- sim$data
Adj <- sim$Adj
net_dist <- sim$net_dist
V <- sim$V

cat("Default simulation generated.\n")
cat("alpha1:", sim$alpha1, "\n")
cat("alpha2:", sim$alpha2, "\n")
cat("pin:", sim$pin, "\n")
cat("pout:", sim$pout, "\n")
cat("i_sim:", 18, "\n")
cat("seed:", sim$seed, "\n")
cat("treated:", sum(sim_dat$Z), "of", nrow(sim_dat), "\n")

# Fit the three package designs. Dual uses Gurobi if available, otherwise GLPK.
m_cov <- netmatch(sim_dat, "Z", c("X1", "X2", "X3"), net_dist,
                  method = "covariate", kappa = 2)
m_single <- netmatch(sim_dat, "Z", c("X1", "X2", "X3"), net_dist,
                     method = "single", kappa = 2)
m_dual <- netmatch(sim_dat, "Z", c("X1", "X2", "X3"), net_dist,
                   method = "dual", kappa = 2, solver = "auto")

print(m_cov)
print(m_single)
print(m_dual)

fit_naive <- netmatch_test(m_dual, "Y", method = "naive")
fit_decay <- netmatch_test(m_dual, "Y", method = "decay", eta = 0.03, rho = 0.10)
fit_design <- netmatch_test(m_dual, "Y", method = "design")

print(fit_naive)
print(fit_decay)
print(fit_design)

sens <- netmatch_sensitivity(
  m_dual,
  "Y",
  eta = c(0.01, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40),
  rho = c(0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60)
)
crit <- critical_sensitivity(m_dual, "Y")

if (requireNamespace("ggplot2", quietly = TRUE)) {
  print(plot_sensitivity(sens, type = "pvalue"))
  print(plot_sensitivity(crit, type = "critical"))
}

cat("Objects available: sim, sim_dat, Adj, net_dist, V, m_cov, m_single, m_dual, fit_naive, fit_decay, fit_design, sens, crit.\n")
