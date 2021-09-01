#
# pkern_indexing.R
# Dean Koch, Sep 2021
# Miscellaneous functions for indexing vectors and matrices
#


#' Inverse column-vectorization indices
#'
#' Column vectorization (as in base::as.vector) builds a length(m*n) vector by stacking
#' the columns of an m X n matrix. The (i,j)th element of the matrix is mapped to the
#' (i + m*(j-1))th element of the vector. This function is the inverse of that mapping.
#' It finds the (i,j) associated with a given element after column-vectorization
#'
#' @param idx a vector of positive integers
#' @param m number of rows in the matrix
#'
#' @return a two column matrix of integers (row and column numbers) with length(idx) rows
#' @export
#'
#' @examples
#'
#' m = 5
#' n = 6
#' matrix.indexing = matrix(1:prod(m, n), m)
#' as.vector(matrix.indexing)
#' pkern_vec2mat(2, m) # finds the matrix indices for the second component
#' pkern_vec2mat(6:10, m) # the function is vectorized
pkern_vec2mat = function(idx, m)
{
  # compute column and row numbers
  cnum = ceiling( idx / m )
  rnum = idx - ( m * (cnum - 1) )

  # return as matrix
  rbind(i=rnum, j=cnum)
}


#' Inverse for base::which
#'
#' base::which takes a logical vector and returns the indices of TRUE elements.
#' This function does the opposite, creating a length-n logical vector that is
#' TRUE at the indices in idx, and FALSE otherwise
#
#' @param idx vector of positive integers no greater than n
#' @param n positive integer
#'
#' @return length-n vector
#' @export
#'
#' @examples
#' n = 25
#' foo = sample(c(TRUE, FALSE), size=n, replace=TRUE)
#' idx.foo = which(foo)
#' pkern_unwhich(idx.foo, n)
pkern_unwhich = function(idx, n)
{
  outvec = rep(FALSE, n)
  outvec[idx] = TRUE
  return(outvec)
}


#' Find column-vectorized index of a subgrid
#'
#' A grid with nx columns and ny rows has grid lines gx=1:nx and gy=1:ny.
#' A subgrid includes only a subset of these grid lines. This function
#' uses the grid line numbers to identifies points on the subgrid by their
#' column-vectorized index with respect to the full grid
#'
#' This is Useful when taking a subset of a dataframe or submatrix of a
#' covariance matrix for gridded data in column-vectorized order.
#'
#' gx or gy NA indicates to use all grid lines. For convenience, list(gx, gy)
#' can supplied in argument gx, but only when gy is NA
#'
#' @param dims c(nx, ny), the number of x and y grid lines in the full grid
#' @param gx vector of positive integers no greater than nx, the x grid lines of the subgrid
#' @param gy vector of positive integers no greater than ny, the y grid lines of the subgrid
#'
#' @return integer vector indexing the subgrid points with respect to the full grid
#' @export
#'
#' @examples
#' dims = c(5,6)
#' pkern_idx_sg(dims, c(1,3,5), c(2, 4))
pkern_idx_sg = function(dims, gx=NULL, gy=NULL)
{
  # handle input as list
  if( length(gx) == 2 )
  {
    if( !is.null(gy)  ) warning('argument gy ignored since gx had length 2')
    gy = gx[[2]]
    gx = gx[[1]]
  }

  # handle default gx, gy (select all grid lines on a dimension)
  if( is.null(gx) ) gx = 1:dims[1]
  if( is.null(gy) ) gy = 1:dims[2]

  # count desired subgrid dimensions and compute result
  ngx = length(gx)
  ngy = length(gy)

  # y coords cycle from highest to lowest in blocks, x coords increase blockwise
  return( rep(gy, ngx) + rep(dims[2] * ( gx - 1 ), each=ngy) )
}


#' Swap grid indices from row-vectorized to column-vectorized order
#'
#' "pkern" uses column-vectorization with y decreasing, to mirror standard mathematical
#' notation for matrices. This utility function is for converting to and from this ordering.
#'
#' dims can be a length-2 vector (as described below) or a RasterLayer object, from which
#' the dimensions are extracted. In the latter case
#'
#' @param dims c(nx, ny), the number of x and y grid lines in source grid
#' @param in.byrow logical indicating if source indexing is row-vectorized
#' @param out.byrow logical indicating if destination indexing is row-vectorized
#' @param flipx logical indicating to flip the grid left-to-right
#' @param flipy logical indicating to flip the grid top-to-bottom
#'
#' @return a reordering for the vectorized grid. ie a permutation of seq(prod(dims))
#' @export
#'
#' @examples
#' dims = c(12, 13)
#' permutation = pkern_r2c(dims, in.byrow=TRUE, out.byrow=FALSE)
#' matrix(permutation, nrow=dims[1], byrow=TRUE) |> as.vector()
pkern_r2c = function(dims, in.byrow=TRUE, out.byrow=FALSE, flipx=FALSE, flipy=FALSE)
{
  # handle raster input and unpack
  if( 'RasterLayer' %in% class(dims) ) dims = dim(dims)[2:1]
  nx = dims[1]
  ny = dims[2]

  # compute a matricized version of indicator index
  idx.mat = seq(prod(dims)) |> matrix(ny, byrow=in.byrow)

  # apply any requested flips
  if(flipx) idx.mat = idx.mat[rev(seq(nx)), ]
  if(flipy) idx.mat = idx.mat[, rev(seq(ny))]
  if(out.byrow) idx.mat = t(idx.mat)

  # vectorize and finish
  return(as.vector(idx.mat))
}

