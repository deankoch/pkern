% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bk.R
\name{anyNA.bk}
\alias{anyNA.bk}
\title{Check for presence of grid points with missing data (NAs)}
\usage{
\method{anyNA}{bk}(x)
}
\arguments{
\item{x}{a bk object}
}
\value{
logical
}
\description{
Returns a logical indicating if any of the grid points are NA
}
\examples{
g = bk_validate(list(gval=rnorm(4^2), gdim=4, gres=0.5))
anyNA(g)
g[1] = NA
anyNA(g)

}
