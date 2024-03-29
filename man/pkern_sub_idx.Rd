% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pkern_index.R
\name{pkern_sub_idx}
\alias{pkern_sub_idx}
\title{Find column-vectorized index of a sub-grid}
\usage{
pkern_sub_idx(gdim, ij = NULL, idx = FALSE, nosort = FALSE)
}
\arguments{
\item{gdim}{integer vector, the number rows and columns in the full grid (in that order)}

\item{ij}{list containing vectors "i" and "j", the sub-grid row and column numbers}

\item{idx}{logical, indicates to return indices (default TRUE) versus logical vector}

\item{nosort}{logical, skips sorting the input vectors in \code{ij}}
}
\value{
integer or logical vector
}
\description{
Returns a logical vector indicating all grid points lying on the specified sub-grid.
A grid point is \code{TRUE} only if both its i and j grid lines are found in \code{ij}.
}
\details{
\code{ij} should be a list containing integer vectors named 'i', 'j', enumerating the
i and j grid lines of the desired sub-grid. If 'i' (or 'j') is missing, the function
automatically specifies all rows (or columns).

If \code{idx=TRUE}, the function computes the vectorized index of the sub-grid points
with respect to the full grid \code{gdim} (see \code{pkern_mat2vec}). Letting \code{i = c(i1, i2, ...im)}
and \code{j = c(j1, j2, ...in)} the function orders the sub-grid points as follows:

(i1, j1), (i2, j1), ... (im, j1), (i1, j2), (i2, j2), ..., (i1, j3), ... (in, jm)

This is the column-major vectorized order for the sub-grid (with y descending and
x ascending), provided the input grid line numbers in \code{ij} are in ascending order.
When \code{nosort=FALSE}, the function orders the input grid lines automatically.
}
\examples{

# example grid and a particular grid point
gdim = c(i=10, j=13)
ij_list = list(i=6, j=3)

# pkern_sub_idx returns a logical vector indexing the point (or the index itself)
is_pt = pkern_sub_idx(gdim, ij_list)
idx_sub = pkern_sub_idx(gdim, ij_list, idx=TRUE)
pkern_plot(is_pt, gdim, col_grid='white', zlab='index', breaks=c('other', idx_sub))

# equivalent call when ij_list is a single point
pkern_mat2vec(ij_list, gdim) == idx_sub

# if i or j are omitted from ij, the function returns the full row or column
is_col2 = pkern_sub_idx(gdim, ij_list['i'])
is_row3 = pkern_sub_idx(gdim, ij_list['j'])
pkern_plot(is_col2, gdim, col_grid='white', breaks=c('other', paste('row', ij_list['i'])))
pkern_plot(is_row3, gdim, col_grid='white', breaks=c('other', paste('col', ij_list['j'])))

# indices in column-vectorized order
pkern_sub_idx(gdim, list(i=2), idx=TRUE)
pkern_sub_idx(gdim, list(j=3), idx=TRUE)
pkern_sub_idx(gdim, idx=TRUE) # call without arguments returns all indices

# bigger sub-grid example
origin_sg = c(5, 2) # assign i,j of top left point
gdim_sg = c(3, 4) # sub-grid dimensions (make sure this fits in gdim!)
ij_list = stats::setNames(Map(\(d, o) o + seq(d) - 1, gdim_sg, origin_sg), c('i', 'j'))
is_sg = pkern_sub_idx(gdim, ij=ij_list)
pkern_plot(is_sg, gdim, col_grid='white', zlab='sub-grid')

# plot the index values: column major vectorization with y descending, x ascending
idx_sg = pkern_sub_idx(gdim, ij=ij_list, idx=TRUE)
vec_order = rep(NA, prod(gdim))
vec_order[is_sg] = as.character(idx_sg)
pkern_plot(vec_order, gdim, col_grid='black', zlab='vector idx')

# example with j indices supplied in reverse (descending) order
ij_list_xflip = modifyList(ij_list, list(j=rev(ij_list[['j']])))

# ordering in the vectors ij$i and ij$j doesn't matter if `nosort=FALSE` or `idx=FALSE`
identical(is_sg, pkern_sub_idx(gdim, ij=ij_list, nosort=TRUE))
all.equal(which(is_sg), pkern_sub_idx(gdim, ij=ij_list_xflip, idx=TRUE))

# when `nosort=TRUE` and `idx=TRUE` we get the same indices but in a different order
idx_sg_xflip = pkern_sub_idx(gdim, ij=ij_list_xflip, idx=TRUE, nosort=TRUE)
all.equal(sort(idx_sg), sort(idx_sg_xflip))
all.equal(idx_sg, idx_sg_xflip)
vec_order[is_sg] = as.character(idx_sg_xflip)
pkern_plot(vec_order, gdim, col_grid='black', zlab='vector index')

}
