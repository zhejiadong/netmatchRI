\name{critical_sensitivity}
\alias{critical_sensitivity}
\title{Critical Sensitivity Curve and Critical Ratio}
\usage{
critical_sensitivity(match, outcome, rho = seq(0, 1, by = 0.01),
  alpha = 0.05, kappa = NULL, weight_type = c("ns", "ntc"))
}
\arguments{
\item{match}{A \code{netmatch} object. Critical-curve interpretation is intended for dual-penalty matching.}
\item{outcome}{Name of the numeric outcome column.}
\item{rho}{Sensitivity-decay values in \code{[0, 1]}. Default: \code{seq(0, 1, by = 0.01)}.}
\item{alpha}{Test level. Default: \code{0.05}.}
\item{kappa}{Network distance threshold. \code{NULL} uses \code{match$kappa}.}
\item{weight_type}{Matched-set weighting: \code{"ns"} or \code{"ntc"}. Default: \code{"ns"}.}
}
\value{
A \code{netmatch_critical_sensitivity} object. \code{summary} reports the
critical eta at \code{rho = 1} and the critical ratio; \code{curve} contains the
critical eta over the requested rho values; \code{detail} contains the supporting
counts and averages. The object also retains the analysis settings and matched
design.
}
\description{
Computes critical eta values where the p-value equals \code{alpha},
together with critical eta at
\code{rho = 1} and the critical ratio.
}
\details{
For matched sets at distance \eqn{d_{s,l}}, the sensitivity bound is
\eqn{M^*_{sl}(\eta,\rho)=\eta\rho^{d_{s,l}-1}M_{sl}}.
At a critical contour point, its value is the critical sensitivity bound
\eqn{M^{\mathrm{crit}}_{sl}(\rho)=M^*_{sl}(\eta^*(\rho),\rho)}.
Within the parameter range, eta values below the critical curve retain
rejection and values above it do not. On the curve, the p-value equals
\code{alpha}.

Let \eqn{S} be the number of matched sets and \eqn{N_\kappa} the number of
relevant matched-set pairs. The critical ratio is
\deqn{
R_{\mathrm{crit}}=
\frac{N_\kappa^{-1}\sum_{s<l:\,l\in\mathcal I_\kappa(s)}
w_sw_lM^{\mathrm{crit}}_{sl}(\rho)}
{S^{-1}\sum_{s=1}^S w_s^2\sigma_{ss}}.
}
The numerator averages one unordered matched-set pair contribution and excludes
the factor of two that appears in the full variance expansion. The ratio is
constant along the critical contour because the p-value threshold fixes the
total weighted contribution.

The summary's \code{critical_ratio} is \code{NA} when no critical
contour is attained within the range. If the baseline p-value is greater
than \code{alpha}, critical eta is \code{NA}. If the baseline p-value equals
\code{alpha} and the sensitivity denominator is positive, critical eta and
the ratio are zero.
}
\examples{
sim <- simulate_netmatch_example(seed = 90141, n = 60, beta_z = 2)
m <- netmatch(sim$data, "Z", c("X1", "X2", "X3"), sim$net_dist,
              method = "dual", kappa = 2)
crit <- critical_sensitivity(m, "Y", rho = seq(0, 1, by = 0.25))
crit
crit$summary
crit$curve
}
