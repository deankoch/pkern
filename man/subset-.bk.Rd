% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bk.R
\name{[<-.bk}
\alias{[<-.bk}
\title{Single-bracket assign}
\usage{
\method{[}{bk}(x, i = NULL, j = NULL, ...) <- value
}
\arguments{
\item{x}{a bk object}

\item{i}{column-vectorized index}

\item{j}{index of layer (only for multi-layer x)}

\item{value}{the replacement values}

\item{...}{ignored}
}
\value{
the "bk" object with the specified subset replaced by value
}
\description{
Behavior depends on the class of i. For character vectors, this assigns to
the named list entries of x (as usual). For numeric indices, it assigns
vectorized grid data values. For multi-layer objects, specify the layer in j
and supply a matrix for replacement
}
\examples{
g = bk_validate(list(gval=rnorm(4^2), gdim=4, gres=0.5))
print(g)
g[1] = NA
print(g)

}
