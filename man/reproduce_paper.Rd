\name{reproduce_paper}
\alias{reproduce_paper}
\title{Reproduce Paper Simulation Results}
\usage{
reproduce_paper(task = c("table1", "ri_grid", "matching_diagnostics",
  "power", "all"), n_rep = if (fast) 2 else 500, fast = FALSE,
  output_dir = file.path(getwd(), "netmatchRI-paper-output"), seed = 1)
}
\arguments{
\item{task}{Reproduction task.}
\item{n_rep}{Number of simulation replications.}
\item{fast}{If \code{TRUE}, run a small smoke-test setting.}
\item{output_dir}{Directory where CSV outputs are written.}
\item{seed}{Random seed.}
}
\value{A named list of raw result data frames.}
\description{
Runs paper-reproduction tasks such as Table 1, the eta-rho grid, matching
diagnostics, and power analysis.
}
\examples{
res <- reproduce_paper("table1", fast = TRUE, output_dir = tempdir())
paper_table1(res$table1)
}
