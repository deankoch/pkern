#
# pkern_indexing.R
# Dean Koch, Sep 2021
# Miscellaneous functions for indexing vectors and matrices
#


#' Column-vectorization indices
#'
#' Column vectorization (as in base::as.vector) builds a length(mn) vector by stacking
#' the columns of an m X n matrix. The (i,j)th element of the matrix is mapped to the
#' `(i + m*(j-1))`th element of the vector. This function returns the index in the vector
#' associated with the supplied i, j values
#'
#' If `ij` has named columns 'x' and 'y', these are interpreted as 'j' and 'i'. If column
#' names 'i' and 'j' also appear, these supercede the columns 'x' and 'y'. If no column names
#' are supplied, the function assumes they are in order 'j', 'i'
#'
#' @param ij n x 2 matrix of column (x) and row (y) indices
#' @param ny number of rows in the matrix
#'
#' @return a vector of indices in the vectorized system
#' @export
#'
#' @examples
#' ny = 10
#' g = expand.grid(i=seq(1, ny), j=seq(1:5))
#' pkern_mat2vec(g, ny)
pkern_mat2vec = function(ij, ny)
{
  # coerce list and other types to matrix
  if( is.list(ij) )
  {
    if( !all( diff(sapply(ij, length)) == 0 ) ) stop('list ij must have equal length entries')
    ij = do.call(cbind, ij)
  }

  # handle vector input (single point)
  if( is.vector(ij) ) ij = t(as.matrix(ij))

  # coerce input to matrix and extract any column names
  ij = as.matrix(ij)
  ij.nm = colnames(ij)

  # identify column of 'j' values (default is first column)
  j.col = ifelse('x' %in% ij.nm, which('x' == ij.nm), 1)
  j.col = ifelse('j' %in% ij.nm, which('j' == ij.nm), j.col)

  # identify column of 'i' values (default is second column)
  i.col = ifelse('y' %in% ij.nm, which('y' == ij.nm), 2)
  i.col = ifelse('i' %in% ij.nm, which('i' == ij.nm), i.col)

  # check for conflicts due to incomplete naming
  if(i.col == j.col) stop('ij should have named columns "x" and "y" or "i" and "j"')

  # check for invalid input ny
  if( any(ij[,i.col] > ny) ) stop('ij contains "i" indices exceeding ny')

  # return the vectorized index
  return( ij[,i.col] + ( ny * ( ij[,j.col] - 1 ) ) )

}

#' Inverse column-vectorization indices
#'
#' Column vectorization (as in base::as.vector) builds a length(mn) vector by stacking
#' the columns of an m X n matrix. The (i,j)th element of the matrix is mapped to the
#' `(i + ny*(j-1))`th element of the vector. This function is the inverse of that mapping.
#' It finds the (i,j) associated with a given element after column-vectorization
#'
#' @param idx a vector of positive integers
#' @param ny number of rows in the matrix
#'
#' @return a two column matrix of integers (row and column numbers) with length(idx) rows
#' @export
#'
#' @examples
#' m = 5
#' n = 6
#' matrix.indexing = matrix(1:prod(m, n), m)
#' as.vector(matrix.indexing)
#' pkern_vec2mat(2, m) # finds the matrix indices for the second component
#' pkern_vec2mat(6:10, m) # the function is vectorized
pkern_vec2mat = function(idx, ny)
{
  # compute column and row numbers
  cnum = ceiling( idx / ny )
  rnum = idx - ( ny * (cnum - 1) )

  # return as matrix
  cbind(i=rnum, j=cnum)
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
pkern_unwhich = function(idx, n, asinteger=FALSE)
{
  outvec = rep(FALSE, n)
  outvec[idx] = TRUE
  return(outvec)
}


#' Find column-vectorized index of a subgrid
#'
#' A grid with nx columns and ny rows has grid lines and i=1:ny, j=1:nx.
#' A subgrid includes only a subset of these grid lines. This function
#' uses the grid line numbers to identify points on the subgrid, returning
#' their column-vectorized index with respect to the full grid. This is often
#' helpful when extracting a subset of a dataframe or sub-matrix of a
#' covariance matrix for gridded data.
#'
#' The function returns indices in column-vectorized order, relative to the order of the
#' elements in its arguments `j` and `i`. This means that if `j = c(j1, j2, ..., in)`
#' and `i = c(i1, i2, ...im)` then the output may has the order
#'
#' (j1, i1), (j1, i2), ... (j1, im), (j2, i1), ... (jm, in).
#'
#' However the function does not order `i` or `j` prior to this mapping, so different
#' output orderings can be induced by changing the ordering of the `i` and `j`. For
#' example providing `i` in descending order results in a vertically flipped grid.
#'
#' NA `i` or `j` indicates to use all grid lines. For convenience, list(i,j)
#' can supplied in argument i. If
#'
#' @param dims c(nx, ny), the number of x and y grid lines in the full grid
#' @param j vector of positive integers no greater than nx, the x grid lines of the subgrid
#' @param i vector of positive integers no greater than ny, the y grid lines of the subgrid
#'
#' @return An integer vector indexing the subgrid points with respect to the full grid
#' @export
#'
#' @examples
#' dims = c(5,6)
#' pkern_idx_sg(dims, c(1,3,5), c(2, 4))
pkern_idx_sg = function(dims, i=NULL, j=NULL)
{
  # handle input as list
  if( length(i) == 2 )
  {
    if( !is.null(j)  ) warning('argument i ignored since j had length 2')
    j = i[[2]]
    i = i[[1]]
  }

  # handle default j, i (select all grid lines on a dimension)
  if( is.null(i) ) i = 1:dims[2]
  if( is.null(j) ) j = 1:dims[1]


  # count desired subgrid dimensions and compute result
  ni = length(i)
  nj = length(j)

  # y coords cycle from highest to lowest in blocks, x coords increase blockwise
  return( rep(i, nj) + rep(dims[2] * ( j - 1 ), each=ni) )
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
  if(flipx) idx.mat = idx.mat[rev(seq(ny)), ]
  if(flipy) idx.mat = idx.mat[, rev(seq(nx))]
  if(out.byrow) idx.mat = t(idx.mat)

  # vectorize and finish
  return(as.vector(idx.mat))
}


#' Snap a set of (possibly irregular) points to a larger grid
#'
#' A heuristics-based method of snapping points in `pts` to the grid `g`.
#'
#' If `pts` is a vector of x coordinates, the function first snaps them to `g`, the computes
#' the lag-1 differences in their their ordered indices in `g` and passes this vector to
#' base::kmeans for a k=2 clustering. The cluster having lower mean is interpreted as a
#' within-column group, whereas the higher mean group delineates jump points between columns.
#' This information is used to group `pts` into columns, which are then snapped to the
#' nearest grid line in `g` by least median distance.
#'
#' If `pts` is a vector, the function returns in a list: the grid line values ('gval'),
#' their index in the full grid ('gid'), and an indexing vector mapping elements (or rows)
#' of `pts` to the snapped grid line ('id'). If `pts` is a 2D set of points (matrix or sf
#' object), the the function is applied to each dimension separately with the results
#' returned as list elements 'x' and 'y'. In the 1D case, `g` should be a vector supplying
#' the grid line locations, and in the 2D case, a list of two vectors (the x and y grid
#' lines) or a RasterLayer.
#'
#' If `regular=TRUE`, the function regularizes the output grid lines (see `pkern_regular`)
#'
#' In the 1D case, `makeplot` can be passed a character string instead of `TRUE`, to set a
#' title for the plot.
#'
#' @param pts, numeric vector or n x 2 matrix (with columns for x and y) of point coordinates
#' @param g, RasterLayer, or list of x, y grid lines, or length-2 vector of dimensions (nx, ny)
#' @param regular, logical, indicating to select a regular set of grid lines
#' @param nstart, passed to base:kmeans
#' @param makeplot, logical, indicating to plot the ordered coordinates and assigned grid lines
#'
#' @return list
#' @export
#'
#' @examples
#' resx = 10
#' resy = 10
#' nx = 25
#' ny = 30
#' xy = list(resx*seq(nx), resy*seq(ny))
#' coords = expand.grid(xy) + rnorm(nx*ny)
#' gxy.snap = pkern_snap(pts=coords, g=xy, regular=F)
#' plot(coords)
#' abline(v = gxy.snap$x$gval)
#' abline(h = gxy.snap$y$gval)
pkern_snap = function(pts, g, regular=FALSE, nstart=25, makeplot=TRUE)
{
  # convert sf type input for pts to a matrix of coordinates, standardize column names
  if( any( c('sf', 'sfc', 'sfg') %in% class(pts) ) ) pts = sf::st_coordinates(pts)
  if( is.data.frame(pts) ) pts = as.matrix(pts)

  # convert RasterLayer type input for g to list of grid line positions
  g.crs = NA
  if( 'RasterLayer' %in% class(g) )
  {
    # copying crs in case makeraster=TRUE
    g.crs = crs(g)
    g = pkern_fromraster(g, 'xy')
  }

  # handle makeplot as title string
  ptitle = NULL
  if( is.character(makeplot) )
  {
    ptitle = makeplot
    makeplot = TRUE
  }

  # handle matrix input for pts
  if( is.matrix(pts) )
  {
    # a single matrix column or row is interpreted as a vector
    if( any(dim(pts) == 1) ) return( pkern_snap(as.vector(pts), g, regular=regular) )

    # order grid lines (x ascending, y descending)
    g[[1]] = sort(g[[1]], decreasing=FALSE)
    g[[2]] = sort(g[[2]], decreasing=TRUE)

    # set up titles and plot parameters for plots
    xtitle = ifelse(makeplot, 'x', FALSE)
    ytitle = ifelse(makeplot, 'y', FALSE)
    if( makeplot ) par(mfrow=c(1,2))

    # call this function recursively on each dimension
    result.x = pkern_snap(pts=as.vector(pts[,1]), g=g[[1]], regular, nstart, makeplot=xtitle)
    result.y = pkern_snap(pts=as.vector(pts[,2]), g=g[[2]], regular, nstart, makeplot=ytitle)
    if( makeplot ) par(mfrow=c(1,1))

    # return the results from the two dimensions in a list
    return( list(x=result.x, y=result.y) )
  }

  # unpack input
  n = length(pts)
  ng = length(g)

  # snap points to nearest grid line number and sort into ascending order
  g.snap = Rfast::Outer(pts, as.numeric(g), '-') |> abs() |> apply(2, which.min)
  g.order = order(g.snap)
  g.snap.order = g.snap[g.order]

  # k-means clustering (k=2) of adjacent grid line separation numbers
  g.kmeans = stats::kmeans(diff(g.snap.order), 2, nstart=nstart)

  # the cluster with the smaller mean should represent points on the same grid line
  idx.lowest = which.min(g.kmeans$centers)

  # separation distances not belonging to the lowest group imply a new grid line has started
  cluster.endpoints = c(which(g.kmeans$cluster != idx.lowest), n)
  cluster.n = length(cluster.endpoints)
  n.bycluster = c(cluster.endpoints[1], diff(cluster.endpoints))

  # grid line numbers of sorted input coordinates
  cluster.ids = do.call(c, lapply(seq(cluster.n), function(id) rep(id, n.bycluster[id])))

  # find the median position within each group and snap to nearest grid line number
  gid = sapply(split(g.snap.order, cluster.ids), function(x) round(stats::median(x)) )
  gval = g[gid]

  # replace the grid line index with a regularized version if requested
  if(regular) gid = pkern_regular(gid, ng)

  # make a diagnostic plot if requested
  if( makeplot )
  {
    plot(seq(n), pts, xlab='input order', ylab='coordinate', pch=16, main=ptitle)
    abline(h = g[gid])
  }

  # compute mapping to points in their original order
  id = cluster.ids[ match(seq(n), g.order) ]

  # return everything in a list
  return(list(gval=gval, gid=gid, id=id))
}


#' Snap an irregular subgrid to nearest regular version
#'
#' Grid line numbers `gid` should be a subset of `seq(ng)`. The function finds the
#' nearest regular subset, in terms of the sum of squared differences. ie it maps
#' the input `gid` to a new subset of `seq(ng)` in which subsequent elements are all
#' separated by an equal number of grid lines (`sep`). If `sep` is not supplied, it
#' is estimated as the median of the differences among the `gid`.
#'
#' Note that input grid lines will be merged when doing so results in a smaller sum of
#' squares.
#'
#' @param gid grid line numbers of the subgrid (sorted in ascending order)
#' @param ng positive integer no less than `max(gid)`, the number of grid lines in full grid
#' @param sep integer between 1 and `(ng-1)`, desired resolution, in number of grid lines
#'
#' @return vector of same length as `gid` containing the new (snapped) grid line numbers
#' @export
#'
#' @examples
#' ng = 100
#' gid = ( seq(1, ng, 5) + rnorm(20, 0, 1) ) |> ceiling() |> sort()
#' diff(gid)
#' gid.snap = pkern_regular(gid)
#' diff(gid.snap)
pkern_regular = function(gid, ng=max(gid), sep=NA)
{
  # default separation distance is the median among all pairs of adjacent grid line numbers
  if( is.na(sep) ) sep = round(stats::median(diff(gid)))

  # check that sep is valid
  sep = round(sep)
  if( ( sep > (ng-1) ) | (sep < 1) ) stop('invalid separation distance sep')

  # candidate grid lines for all possible origins
  g.test = lapply(seq(sep), function(x) seq(x, ng, by=sep))
  ng.test = sapply(g.test, length)

  # storage for the loop
  ss.test = rep(NA, sep)
  map.test = vector(mode='list', length=sep)

  # exhaustive search for best origin
  for(idx.test in seq(sep))
  {
    # map snapped points to candidate grid lines
    map.test[[idx.test]] = Rfast::Outer(gid, as.numeric(g.test[[idx.test]]), '-') |>
      abs() |> apply(2, which.min)

    # compute a sum of squares over grid line snapping distance
    ss.test[idx.test] = sum( ( g.test[[idx.test]][ map.test[[idx.test]] ] - gid )^2 )
  }

  # select alignment with least sum of squares
  idx.best = which.min(ss.test)
  map.test[[idx.best]]

  # return the grid line numbers after snapping
  return( g.test[[idx.best]][ map.test[[idx.test]] ] )
}


