\name{paper_table1}
\alias{paper_table1}
\alias{paper_power_plot}
\alias{paper_grid_plot}
\title{Paper Reproduction Helpers}
\usage{
paper_table1(raw, alpha = 0.05)
paper_power_plot(raw, alpha = 0.05)
paper_grid_plot(raw, alpha = 0.05)
}
\arguments{
\item{raw}{Raw output from \code{reproduce_paper()}.}
\item{alpha}{Test level.}
}
\value{
\code{paper_table1()} returns a rejection-rate table. Plot helpers return a
plot object when \pkg{ggplot2} is installed, otherwise a base plot.
}
\description{
Summarize or plot raw paper reproduction outputs.
}
