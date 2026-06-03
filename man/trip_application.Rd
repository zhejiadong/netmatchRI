\name{trip_application}
\alias{trip_application}
\title{Run the TRIP Real-Data Application}
\usage{
trip_application(data_file,
  output_dir = file.path(getwd(), "trip-application-output"),
  methods = c("covariate", "single", "dual"), kappa = 2,
  eta_grid = c(0.01, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40),
  rho_grid = c(0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60),
  d0 = 2)
}
\arguments{
\item{data_file}{Path to \code{TRIP/data_clear.RData} or a compatible file.}
\item{output_dir}{Directory where application CSV outputs are written.}
\item{methods}{Matching methods to run.}
\item{kappa}{Network-distance threshold.}
\item{eta_grid}{Eta values for sensitivity analysis.}
\item{rho_grid}{Rho values for sensitivity analysis.}
\item{d0}{Maximum set distance for model-assisted RI.}
}
\value{A named list with matched objects and result tables.}
\description{
Runs the TRIP real-data illustration using a local data file containing
\code{ds_use} and \code{Adj}. Raw TRIP data are not bundled with the package.
}
\examples{
\dontrun{
trip_application("../TRIP/data_clear.RData")
}
}
