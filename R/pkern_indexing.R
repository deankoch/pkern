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


#' Snap a set of (possibly irregular) points to a larger grid
#'
#' @param pts, input coordinates to snap
#' @param g, grid lines of larger grid
#' @param regular
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
pkern_snap = function(pts, g, regular=FALSE)
{

  # convert sf type input for pts to a matrix of coordinates, standardize column names
  if( any( c('sf', 'sfc', 'sfg') %in% class(pts) ) ) pts = sf::st_coordinates(pts)
  if( is.data.frame(pts) ) pts = as.matrix(pts)

  # convert RasterLayer type input for g to list of grid line positions
  if( 'RasterLayer' %in% class(g) ) g = pkern_fromRaster(g, 'xy')

  # handle matrix input for pts
  if( is.matrix(pts) )
  {
    # a single matrix column or row is interpreted as a vector
    if( any(dim(pts) == 1) ) return( pkern_snap(as.vector(pts), g, regular=regular) )

    # otherwise assume 2 sets of coordinates with matching arguments in xy
    result.x = pkern_snap(pts=as.vector(pts[,1]), g=g[[1]], regular=regular)
    result.y = pkern_snap(pts=as.vector(pts[,2]), g=g[[2]], regular=regular)
    return( list(x=result.x, y=result.y) )
  }

  # unpack input
  n = length(pts)
  ng = length(g)

  # snap the points to nearest grid line number and sort into ascending order
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

  # map (sorted) input data onto grid line numbers
  cluster.ids = do.call(c, lapply(seq(cluster.n), function(id) rep(id, n.bycluster[id])))

  # find the median position within each group and snap to nearest grid line number
  cluster.snap = sapply(split(g.snap.order, cluster.ids), function(x) round(stats::median(x)) )

  # if requested, find a regularized version
  if(regular)
  {
    # use median separation distance (in terms of grid line numbers)
    g.sep = round( stats::median(diff(cluster.snap)) )

    # candidate grid lines for all possible origins
    g.test = lapply(seq(g.sep), function(x) seq(x, ng, by=g.sep))
    ng.test = sapply(g.test, length)

    # storage for the loop
    ss.test = rep(NA, g.sep)
    map.test = vector(mode='list', length=g.sep)

    # exhaustive search for best origin
    for(idx.test in seq(g.sep))
    {
      # map snapped points to candidate grid lines
      map.test[[idx.test]] = Rfast::Outer(cluster.snap, as.numeric(g.test[[idx.test]]), '-') |>
        abs() |> apply(2, which.min)

      # compute a sum of squares over grid line snapping distance
      ss.test[idx.test] = sum( ( g.test[[idx.test]][ map.test[[idx.test]] ] - cluster.snap )^2 )
    }

    # select alignment with least sum of squares
    idx.best = which.min(ss.test)

    # trim outer grid lines which aren't mapped to anything
    g.start = min(map.test[[idx.best]])
    g.end = max(map.test[[idx.best]])

    # update the grid line index and points map
    cluster.snap = g.test[[idx.best]][g.start:g.end]
    cluster.ids = match(cluster.ids, ( map.test[[idx.best]] - g.start + 1 ))
  }

  # compute mapping to points in their original order and return list
  pts.map = cluster.ids[ match(seq(n), g.order) ]
  return(list(gval=g[cluster.snap], gid=cluster.snap, id=pts.map))
}






