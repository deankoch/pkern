#
# pkern_indexing.R
# Dean Koch, Oct 2021
# Miscellaneous functions for indexing vectors and matrices
#


#' Column-vectorization indices
#'
#' Returns the index in the vector associated with the supplied i, j values
#'
#' Column vectorization (as in base::as.vector) builds a length(mn) vector by stacking
#' the columns of an m X n matrix. The (i,j)th element of the matrix is mapped to the
#' `(i + m*(j-1))`th element of the vector.
#'
#' `ij` can be a matrix or a list of length-n vectors "i" and "j" (in that order), or a vector
#' representing a single point at the given row and column number. `ni` can be a vector
#' c(ni, nj), in which case its first element is selected.
#'
#' @param ij nx2 matrix, the row (y) and column (x) indices
#' @param ni number of rows in the matrix
#' @param simplified, if FALSE, the function returns an n x 1 matrix
#'
#' @return a vector of indices in the vectorized system
#' @export
#'
#' @examples
#'
#' # get column-vectorized ordering of points from matrix indices
#' ni = 10
#' g = expand.grid(i=seq(1, ni), j=seq(1:ni))
#' pkern_mat2vec(g, ni)
pkern_mat2vec = function(ij, ni, simplified=TRUE)
{
  # handle vector input to ni
  if( length(ni) > 1 ) ni = ni[1]

  # coerce list to matrix
  if( is.list(ij) )
  {
    if( !all( diff(sapply(ij, length)) == 0 ) ) stop('elements of list ij must have equal length')
    ij = do.call(cbind, ij)
  }

  # handle vector input (single point)
  if( is.vector(ij) ) ij = matrix(ij, 1)

  # coerce input to matrix
  ij = as.matrix(ij)

  # check for invalid input
  if( any(ij[,1] > ni) ) stop('ij contains "i" indices exceeding ni')

  # return the vectorized index
  idx = ij %*% c(1, ni) - ni
  if( !simplified ) return( idx )
  return( as.vector(idx) )
}


#' Inverse column-vectorization indices
#'
#' This function is for doing the inverse of `as.vector(amatrix)`. It returns the row and
#' column numbers (i,j) associated with a given vector element after column-vectorization.
#'
#' If a grid size vector `c(ni, nj)` is passed to `ni`, the function uses its first element.
#'
#' @param idx a vector of positive integers
#' @param ni number of rows in the matrix
#' @param out either 'matrix' or 'list'
#'
#' @return a two column matrix of integers (row and column numbers) with `length(idx)` rows
#' @export
#'
#' @examples
#'
#' # show how elements are ordered in `base::matrix`
#' ni = 5
#' nj = 6
#' matrix.indexing = matrix(1:prod(ni, nj), ni)
#' print(matrix.indexing)
#' as.vector(matrix.indexing)
#'
#' # doing the inverse
#' pkern_vec2mat(2, ni)
#' pkern_vec2mat(c(1,2,7), ni)
pkern_vec2mat = function(idx, ni, out='matrix')
{
  # handle vector input to ni
  if( length(ni) > 1 ) ni = ni[1]

  # compute column and row numbers
  cnum = ceiling( idx / ni )
  rnum = idx - ( ni * (cnum - 1) )

  # return as matrix
  if(out == 'matrix') return( cbind(i=rnum, j=cnum) )
  if(out == 'list') return( list(i=rnum, j=cnum) )
}


#' Find column-vectorized index of a subgrid
#'
#' Returns the column-vectorized index of a subgrid with respect to the full grid
#' of size `gdim`, based on the supplied grid line numbers `ij`. The returned vector
#' maps to points in the full grid (in column-vectorized order) lying at the
#' intersections of `ij` (a subset of `seq(n)`, where `n=prod(gdim)`).
#'
#' NA `ij[1]` indicates to use all rows and NA `ij[2]` indicates to use all columns.
#' If `j = c(j1, j2, ..., in)` and `i = c(i1, i2, ...im)` the function returns:
#'
#' (j1, i1), (j1, i2), ... (j1, im), (j2, i1), (j2, i2), ..., (j3, i1), ... (jm, in).
#'
#' By default, the function sorts the grid lines in `ij` into ascending order, so that
#' the output is in column-vectorized order. `nosort=TRUE` skips the sorting, allowing
#' alternative output orders.
#'
#' @param gdim c(ni, nj), the number rows and columns in the full grid
#' @param ij list containing vectors "i" and "j", the subgrid rows and columns
#' @param nosort logical, skips sorting the input vectors in `ij`
#'
#' @return integer vector, with length equal to the product of the lengths of "i" and "j"
#' @export
#'
#' @examples
#'
#' gdim = c(5,6)
#' pkern_idx_sg(gdim)
#' ij = list(i = c(2,4), j = c(1,3,5))
#' pkern_idx_sg(gdim, ij)
pkern_idx_sg = function(gdim, ij=NULL, nosort=FALSE)
{
  # check input and set defaults
  ijnm = c('i', 'j')
  if( is.null(ij) ) ij = stats::setNames(lapply(gdim, seq), ijnm)
  if( !any( ijnm %in% names(ij) ) ) ij = stats::setNames(ij, ijnm)
  if( !all( ijnm %in% names(gdim) ) ) gdim = stats::setNames(gdim, ijnm)

  # handle default i and j (select all grid lines along each dimension)
  if( is.null(ij[['i']]) ) ij[['i']] = seq( gdim['i'] )
  if( is.null(ij[['j']]) ) ij[['j']] = seq( gdim['j'] )

  # sort the inputs by default
  if( !nosort ) ij = lapply(ij, sort)

  # count desired subgrid dimensions
  nij = sapply(ij, length)

  # index mode
  return( rep(ij[['i']], nij['j']) + rep(gdim['i'] * (ij[['j']] - 1), each=nij['i']) )
}

#' Swap grid indices from row-vectorized to column-vectorized order
#'
#' Returns an indexing vector for swapping between various possible vectorization orderings.
#'
#' `pkern` uses column-vectorization, with y decreasing, to mirror mathematical notation
#' for matrices where rows are numbered from top to bottom. Packages `raster` and `graphics`
#' have different conventions, so this function is useful for reordering vectorizations when
#' switching between methods.
#'
#' `gdim` can be a length-2 vector or a RasterLayer object, from which the dimensions are
#' extracted.
#'
#' @param gdim c(ni, nj), the number of rows and columns in the grid
#' @param in.byrow logical indicating if source indexing is row-vectorized
#' @param out.byrow logical indicating if destination indexing is row-vectorized
#' @param flipx logical indicating to flip the grid left-to-right
#' @param flipy logical indicating to flip the grid top-to-bottom
#'
#' @return a reordering for the vectorized grid. ie a permutation of `seq(prod(gdim))`
#' @export
#'
#' @examples
#'
#' # with column vectorized input and no flips, the function does no reindexing
#' gdim = c(12, 13)
#' all( pkern_r2c(gdim, in.byrow=FALSE) == seq(prod(gdim)) )
#'
#' # flip a grid horizontally and vertically
#' matrix(pkern_r2c(gdim, in.byrow=FALSE, flipx=TRUE), gdim)
#' matrix(pkern_r2c(gdim, in.byrow=FALSE, flipy=TRUE), gdim)
#'
#' # default behaviour is to convert from row to column vectorization
#' roworder = pkern_r2c(gdim, in.byrow=FALSE, out.byrow=TRUE)
#' print(roworder)
#' print(matrix(roworder, gdim[1], byrow=TRUE))
pkern_r2c = function(gdim, in.byrow=TRUE, out.byrow=FALSE, flipx=FALSE, flipy=FALSE)
{
  # handle raster input and unpack
  if( 'RasterLayer' %in% class(gdim) ) gdim = dim(gdim)[1:2]

  # compute a matricized version of indicator index
  idx = seq(prod(gdim)) |> matrix(gdim[1], byrow=in.byrow)

  # apply any requested flips
  if( (!out.byrow & flipx) | (out.byrow & flipy) ) idx = idx[, rev(seq(gdim[2]))]
  if( (!out.byrow & flipy) | (out.byrow & flipx) ) idx = idx[rev(seq(gdim[1])), ]
  if(out.byrow) idx = t(idx)

  # vectorize and finish
  return(as.vector(idx))
}


#' Rotate a rectangular array by 45 degrees clockwise
#'
#' performs the rotation f(x,y) = (x + y, -x + y) about the center of a rectangular
#' array (45 degrees counterclockwise), with distances scaled by sqrt(2)/2 to snap
#' grid cells to those of a larger array.
#'
#' NAs are assigned to all output grid points not mapped to `z`.
#'
#' When `z` is an object containing the grid size (matrix, RasterLayer or list),
#' `gdim` is ignored (with a warning) and can be omitted. When `z` is a vector, it
#' should be in column-vectorized order and have length `prod(gdim)`.
#'
#' The function sets
#'
#' @param z either a numeric matrix, its column-vectorization, a RasterLayer, or a list
#' @param gdim integer vector c(ni, nj), the number of rows and columns in the grid
#'
#' @return a numeric matrix containing the rotated array
#' @export
#'
#' @examples
#'
#' # a wide example modified from `graphics::.filled.contour`
#' gdim = c(73,50)
#' y = seq(-pi, pi, length.out=gdim[1])
#' x = seq(-pi, pi, length.out=gdim[2])
#' d = sqrt(outer(y^2, x^2, '+'))
#' z = c( sin(d^2) * exp(-d) )
#' pkern_plot(z, gdim)
#' pkern_plot( pkern_r45(z, gdim) )
#'
#' # a tall example with raster
#' if( requireNamespace('raster') ) {
#' r = system.file('external/rlogo.grd', package='raster') |> raster::raster(band=1)
#' rg = pkern_fromRaster(r)
#' pkern_plot(rg)
#' pkern_plot( pkern_r45(rg) )
#' }
pkern_r45 = function(z, gdim=NULL)
{
  # handle raster input
  if( 'RasterLayer' %in% class(z) )
  {
    if( !is.null(gdim) ) warning('argument to gdim was ignored')
    gdim = pkern_fromRaster(z, 'gdim')
    z = pkern_fromRaster(z, 'values')
  }

  # handle list input
  if( is.list(z) )
  {
    if( !is.null(gdim) ) warning('argument to gdim was ignored')
    gdim = z[['gdim']]
    z = z[['values']]
  }

  # handle NULL gdim
  if( is.null(gdim) )
  {
    if( !is.matrix(z) ) stop('Failed to detect grid size. Try setting gdim')
    gdim = dim(z)
  }

  # map from original (padded) i, j to rotated coordinates in larger array
  ij = pkern_vec2mat(seq(prod(gdim)), gdim)
  ij.rot = ij %*% matrix(c(1,1,-1,1), 2)

  # translate to origin at (1,1), copy new grid size
  ij.rot = 1 + ij.rot - apply(ij.rot, 2, \(v) rep(min(v), length(v)))
  gdim.rot = apply(ij.rot, 2, max)

  # map points to their new positions in larger array
  idx.rot = pkern_mat2vec(ij.rot, gdim.rot)
  zrot = rep(NA, prod(gdim.rot))
  zrot[idx.rot] = z

  # return as matrix
  return( matrix(zrot, gdim.rot) )
}

#' Return coordinates of a grid of points in column-vectorized order
#'
#' Expands a set of "y" and "x" grid line locations in column vectorized order,
#' with "y" values decreasing (fastest, ie in cycles) and "x" values strictly
#' nondecreasing.
#'
#' This is similar to `base::expand.grid`, except input data are by default ordered
#' as described above. Sorting can be switched off for code optimization (`nosort=TRUE`)
#' but the coordinates in `g` must have "y" in decreasing order and "x" in increasing
#' order.
#'
#' `g` should either be a list of "y" and "x" coordinates, or an object from which they
#' can be extracted (RasterLayer, or list output from functions like `pkern_snap`,
#' `pkern_fromRaster`). If `g` is a vector of length 2, it is assumed to be the
#' dimensions (ni, nj) of the grid, and a sequence of integer grid line positions
#' is generated.
#'
#' @param g a list of grid line coordinates, or an object containing them (see details)
#' @param out character indicating return value type, either 'list' or 'matrix'
#' @param nosort logical indicating to assume input `g` is sorted correctly
#' @param quiet logical indicated to drop warnings
#'
#' @return a matrix or list of grid coordinates in column vectorized order
#' @export
#'
#' @examples
#' pkern_coords(g=c(3,2))
#' pkern_coords(g=list(y=1:5, x=2:3))
pkern_coords = function(g, out='matrix', nosort=FALSE, quiet=FALSE)
{
  # handle numeric vectors for g
  yxnm = c('y', 'x')
  if( is.numeric(g) )
  {
    err.msg = 'unrecognized input g'
    if( length(g) == 2 ) { g = stats::setNames(lapply(g, seq), yxnm) } else { stop(err.msg) }
  }

  # handle various input classes
  if( 'RasterLayer' %in% class(g) ) g = pkern_fromRaster(g, what='yx')
  if( !is.null( g[['yx']] ) ) g = g[['yx']]
  if( !is.null( g[['g']] ) ) g = g[['g']]

  # handle unnamed input
  if( !all(yxnm %in% names(g)) )
  {
    warning('Assuming first two elements of g are "y" and "x" values, in that order.')
    names(g)[1:2] = yxnm
  }

  # order "y" and "x" and find grid dimensions
  g = g[yxnm]
  gdim = sapply(g, length)
  ng = prod(gdim)

  # sort the coordinates unless otherwise specified
  if( !nosort )
  {
    g[['y']] = sort(g[['y']], decreasing=TRUE)
    g[['x']] = sort(g[['x']], decreasing=FALSE)
  }

  # sanity check warning if the first two elements not in correct order
  msg.order = 'y and x not in expected order (decreasing, increasing). try nosort=FALSE'
  if( !quiet & !all( sapply(g, \(gl) diff(gl[1:2]) ) * c(-1, 1) > 0 ) ) warning(msg.order)

  # return in requested class
  olist = Map(\(gl, i) gl[i], g, pkern_vec2mat(seq(ng), gdim, out='list'))
  if( out == 'list' ) return(olist)
  return( do.call(cbind, olist) )
}
