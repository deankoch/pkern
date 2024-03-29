% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bk.R
\name{as.logical.bk}
\alias{as.logical.bk}
\title{Coerce grid values to logical}
\usage{
\method{as.logical}{bk}(x, ...)
}
\arguments{
\item{x}{a bk object}

\item{...}{further arguments to as.logical}
}
\value{
a "bk" object with logical data values
}
\description{
Coerce grid values to logical
}
\examples{
g = bk_validate(list(gval=sample(c(0,1), 4^2, replace=TRUE), gdim=4, gres=0.5))
g[]
as.logical(g)[]

# "range" for logical is reported as integer
summary(as.logical(g))

}
