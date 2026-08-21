\name{matched_data}
\alias{matched_data}
\title{Extract Data from a Matched Design}
\usage{
matched_data(object, drop_unmatched = TRUE)
}
\arguments{
\item{object}{A \code{netmatch} object.}
\item{drop_unmatched}{If \code{TRUE}, omit unmatched units. If \code{FALSE}, retain them with missing \code{subclass} and zero \code{weights}.}
}
\value{
A data frame containing the original columns plus full-length matched-set membership in \code{subclass} and unit matching weights in \code{weights}.
}
\description{
Returns the original unit-level data with matched-set membership and matching weights appended.
}
