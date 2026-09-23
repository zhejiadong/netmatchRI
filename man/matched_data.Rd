\name{matched_data}
\alias{matched_data}
\title{Extract Data from a Matched Design}
\usage{
matched_data(object, drop_unmatched = TRUE)
}
\arguments{
\item{object}{A \code{netmatch} object.}
\item{drop_unmatched}{If \code{TRUE}, omit unmatched units. Default: \code{TRUE}.}
}
\value{
The unit-level data with \code{subclass} and \code{weights}.
}
\description{
Returns the original unit-level data with matched-set membership and matching weights.
}
