#
# pkern_snap.R
# Dean Koch, Sep 2021
# Functions for snapping points to a grid
#



#' Plots a set of 1 or 2d points and their snapped versions
#'
#' The function plots the points in `pts` over the grid lines in `gyx`, with the option of
#' also drawing connecting lines between points and grid intersections to show snapping
#' distances.
#'
#' This is mainly intended for visualizing the output of `pkern_snap` and related functions.
#'
#' `pts` can be a vector of 1d positions (in which case the plot shows input order against
#' position), or a matrix whose columns are the "y" and "x" coordinates (in that order,
#' unless they are named), or a data frame (which is coerced to matrix), or a list of
#' coordinate vectors, or an `sf` points object.
#'
#' `gyx` can be a numeric vector in the 1d case; or list of two (in "y", "x" order); or a
#' list containing the grid line locations (as vector or list) in element "gyx", and,
#' optionally, "pmap" which provides a mapping for drawing connecting lines. For example
#' you can pass the return value of `pkern_snap` as `gxy`
#'
#' If supplied, `gyx$pmap` should be an integer vector (or list of two) with length(s)
#' matching those in `pts` - its elements index grid lines, so they should have values
#' drawn (possibly with repitition) from the sequences `seq_along(gyx$gyx$y)` and
#' `seq_along(gyx$gyx$x)`
#'
#' When `gyx$gyx` but not `gyx$pmap` is supplied, it is assumed that `gyx$gyx` and `pts` have
#' the same length and ordering, and lines are drawn connecting their respective elements.
#' When `gyx` is a list without named element "gyx", it is assumed to be a list of grid line
#' vectors (and no connecting lines are drawn).
#'
#' `ppars` is a list with any of the following (optional) named elements: "gcol", "mcol",
#' "dcol", and "pcol". They are the colors of, respectively, the grid lines, the snapped grid
#' lines connecting lines, and points (passed as parameter "col" to the appropriate `graphics`
#' call).
#'
#' @param gyx numeric vector, matrix, dataframe or list, the grid line coordinates (see details)
#' @param pts numeric vector, matrix, dataframe or list, the point coordinates (see details)
#' @param ppars list, (optional) plotting parameters "gcol", "mcol", "dcol", "pcol", and "yaxis"
#'
#' @return draws a plot and returns nothing
#' @export
#'
#' @examples
#' # 1d case
#' gyx = seq(1, 100, by=2)
#' pkern_snap_plot(gyx)
#' ng = length(gyx)
#' pts = gyx + stats::rnorm(ng)
#' pkern_snap_plot(gyx, pts)
#'
#' # plot again with with points mapped to grid line that generated them
#' pkern_snap_plot(gyx=list(gyx=g, pmap=seq(gyx)), pts)
#'
#' # plot again with a random mapping
#' pkern_snap_plot(list(gyx=gyx, pmap=sample(seq(ng), ng)), pts)
#'
#' # 2d case
#' gdim = c(30, 25)
#' sep = c(2, 3)
#' gyx = Map(\(d, s) seq(1, d, by=s), d=gdim, s=sep)
#' yx = expand.grid(gyx)
#' pts = yx + stats::rnorm(prod(dim(yx)))
#' pkern_snap_plot(gyx, pts)
#'
#' # plot again with with points mapped to grid line that generated them
#' pmap = Map(\(p, s) match(p, s), p=yx, s=gyx)
#' pkern_snap_plot(gyx=list(gyx=gyx, pmap=pmap), pts)
#'
#' # TODO: examples with pkern_snap
pkern_snap_plot = function(gyx, pts=NULL, ppars=list())
{
  # coerce various input types and set expected column order for 2d case
  sfnm = c('sf','sfc', 'sfg')
  if( any(sfnm %in% class(pts)) ) pts = sf::st_coordinates(pts)[,2:1]
  if( is.data.frame(pts) ) pts = as.matrix(pts)
  if( is.matrix(pts) ) pts = apply(pts, 2, identity, simplify=FALSE)
  yxnm = stats::setNames(nm=c('y', 'x'))

  # unpack gyx as list
  pmap = NULL
  is.1d = TRUE
  if( is.list(gyx) )
  {
    # unpack list element gl as list (connecting line mode)
    if( !is.null( gyx[['gyx']] ) )
    {
      # extract mapping and grid line vector(s)
      pmap = gyx[['pmap']]
      gyx = gyx[['gyx']]

      # set defaults as needed and catch some invalid input errors
      defaults_result = pkern_snap_plot_defaults(gyx, pmap, pts)
      pmap = defaults_result[['pmap']]
      pts = defaults_result[['pts']]
      is_1d = defaults_result[['is_1d']]

    } else { is.1d = FALSE }
  }

  # unpack plotting parameters and set defaults
  gcol = ifelse( is.null(ppars[['gcol']] ), NA, ppars[['gcol']])
  mcol = ifelse( is.null(ppars[['mcol']] ), 'grey40', ppars[['mcol']])
  dcol = ifelse( is.null(ppars[['dcol']] ), 'black', ppars[['dcol']])
  pcol = ifelse( is.null(ppars[['pcol']] ), 'grey50', ppars[['pcol']])

  # set defaults as needed and catch some invalid input errors
  defaults_result = pkern_snap_plot_defaults(gyx, pmap, pts)
  pmap = defaults_result[['pmap']]
  pts = defaults_result[['pts']]
  is_1d = defaults_result[['is_1d']]

  # handle 1D case
  if( is_1d )
  {
    # set limits for the plot
    pseq = seq_along(pts)
    ilim = range(pseq)
    plim = range(gyx, pts)

    # scatterplot of points and grid lines with connecting lines showing snap position
    plot(x=plim, y=ilim, xlab='x coordinate', ylab='input order', pch=NA, xaxt='n')
    axis(side=1, at=gyx)
    graphics::abline(v=gyx, col=gcol)
    if( !is.null(pmap) )
    {
      # draw sub-grid line ticks
      gx_occupied = gyx[unique(pmap)]
      axis(side=1, at=gx_occupied, col.ticks='black', labels=FALSE)

      # sub-grid lines
      graphics::abline(v=gx_occupied, col=mcol)
      lapply(pseq, \(x) graphics::lines(x=c(pts[x], gyx[pmap[x]]), y=rep(x,2), col=dcol) )
    }
    graphics::points(pts, pseq, col=pcol)

    # finish
    return( invisible() )
  }

  # set names as needed for 2d case
  if( !all( yxnm %in% names(gyx) ) ) names(gyx)[1:2] = yxnm
  if( !all( yxnm %in% names(pts) ) ) names(pts)[1:2] = yxnm
  if( !is.null(pmap) ) if( !all( yxnm %in% names(pmap) ) ) names(pmap)[1:2] = yxnm

  # set limits for the plot
  pseq = seq(length(pts[[1]]))
  plim = lapply(yxnm, \(d) range(gyx[[d]], pts[[d]]) )

  # create the plot, add grid line ticks
  plot(plim, pch=NA, xlab='x', ylab='y', asp=1, xaxt='n', yaxt='n')
  axis(side=1, at=gyx[[2]], col.ticks='darkgrey')
  axis(side=2, at=gyx[[1]], col.ticks='darkgrey')
  graphics::abline(h=gyx[['y']], v=gyx[['x']], col=gcol)

  # add occupied grid lines add snapping vectors
  if( !is.null(pmap) )
  {
    # add sub-grid line ticks
    gy_occupied = gyx[['y']][ unique(pmap[['y']]) ]
    gx_occupied = gyx[['x']][ unique(pmap[['x']]) ]
    axis(side=1, at=gx_occupied, col.ticks='black', labels=FALSE)
    axis(side=2, at=gy_occupied, col.ticks='black', labels=FALSE)

    # sub-grid lines
    graphics::abline(h=gy_occupied, v=gx_occupied, col=mcol)

    # add snapping lines
    yxvec = lapply(yxnm, \(d) lapply(pseq, \(i) c(pts[[d]][i], gyx[[d]][ pmap[[d]][i] ]) ) )
    Map(\(x, y) graphics::lines(x=x, y=y, col=dcol), x=yxvec[['x']], y=yxvec[['y']])
  }

  # draw the points and finish
  graphics::points(pts, col=pcol)
  return( invisible() )
}


#' Set defaults for pmap and pts in pkern_snap_plot
#'
#' Helper function for pkern_snap_plot
#'
#' Parses the input arguments `gyx`, `pmap` and `pts`, returning a modified
#' version of `pmap` and `pts` with defaults set as needed, along with the
#' dimension.
#'
#' @param gyx numeric vector of grid-line coordinates, or two such vectors (y, x) in a list
#' @param pmap integer vector of mapping indices, or two such vectors (y, x) in a list, or NULL
#' @param pts integer vector of point coordinates, or two such vectors (y, x) in a list, or NULL
#'
#' @return list containing `pmap`, `pts`, and a logical indicating if input was 1d
#' @export
#'
#' @examples
pkern_snap_plot_defaults = function(gyx, pmap, pts)
{
  # define some error messages
  err1 = 'mismatch in dimensions. Check that pmap and pts have the same class'
  err2 = 'mismatch in number of points and grid lines. Check pmap argument'

  # error 1: dimension mismatch check
  in_list = list(gyx=gyx, pmap=pmap, pts=pts)
  is_null = sapply(in_list, is.null)
  is_list = sapply(in_list, is.list)
  if( any(diff( is_list[!is_null] ) ) ) stop(err1)
  is_1d = !is_list[!is_null][1]

  # 1d case
  if( is_1d )
  {
    # unpack mapping info and assign defaults for 1d case
    if( is.null(pts) )
    {
      # treat grid-line intersections as points when pts is not supplied
      if( is.null(pmap) ) pmap = seq_along(gyx)
      pts = gyx[pmap]
    }

    # check for mismatched snapping vectors
    if( !is.null(pmap) & ( length(pmap) != length(pts) ) ) stop(err2)

    return(list(pmap=pmap, pts=pts, is_1d=TRUE))
  }

  # 2d case handled differently
  if( is.null(pts) )
  {
    # treat grid-line intersections as points when pts is not supplied
    if( is.null(pmap) ) pmap = Map(\(p, s) match(p, s), p=expand.grid(gyx), s=gyx)
    pts = Map(\(s, i) s[i], s=gyx, i=pmap)
  }

  # check for mismatched snapping vectors
  if( !is.null(pmap) & any( sapply(pmap, length) != sapply(pts, length) ) ) stop(err2)

  return(list(pmap=pmap, pts=pts, is_1d=FALSE))
}



#' Snap an irregular set of 2d points to the nearest regular subgrid
#'
#' Snaps points to a subgrid, either by calling `pkern_snap_1d` twice, or by using the
#' Hungarian algorithm to solve the grid assignment problem.
#'
#' `gyx` can be a RasterLayer, a list of y and x coordinates, or a list containing the
#' coordinates (in element `gyx$gyx`) as returned by `pkern_fromRaster`.
#'
#' `sep` is a positive integer (or vector of two) specifying the factor(s) by which to
#' multiply the outer grid resolution to get the subgrid. Equivalently, this is the
#' separation distance between adjacent subgrid lines, given in terms of the number of
#' outer grid lines - eg. if `sep = c(1,1)` then all outer grid lines are included in the
#' subgrid; And if `sep=c(2,1)`, then only every second y grid line is included. `pkern`
#' includes a helper function for choosing `sep` automatically. This is done by default
#' when `sep` is not specified.
#'
#' By default (`distinct=TRUE`) the function uses the Hungarian algorithm to find
#' the mapping which minimizes the total sum of squared snapping distances, under the
#' constraint that no more than one input point is assigned to each subgrid point.
#' When `distinct=FALSE` the function simply maps each point to the nearest x and y grid
#' lines separately (minimizing total Manhattan distance)
#'
#' The output list contains info about the configuration of the grid and subgrid, along
#' with several indexing vectors for mapping between the input points and the grid points -
#'
#'  `gdim`, the input grid dimensions (ny, nx)
#'  `gres`, the input grid resolution (copied from `gyx` or set to `c(1,1)`)
#'  `gyx`, the input grid line coordinates (copied from input `gyx`)
#'  `sep`, the selected subgrid spacing
#'  `pmap`, list of integer vectors, dimension-wise mapping from `pts` to `gyx`
#'  `gli`, list of integer vectors, dimension-wise grid line numbers forming the subgrid
#'  `pvec`, integer vector mapping from `pts` to full grid, in column-vectorized order
#'  `pdist`, numeric vector of Euclidean snapping distances
#'  `sg`, list with information about the subgrid:
#'    `gdim`, the subgrid dimensions (ny, nx)
#'    `gyx`, the subgrid line coordinates (subsets of `gyx` in parent list)
#'    `gres`, the subgrid resolution (equal to `sep * gres` in parent list)
#'    `pvec`, integer vector mapping from `pts` to subgrid, in column-vectorized order
#'    `ipvec`, integer vector, the inverse of the above mapping (with NAs for unmapped grid points)
#'
#' @param gyx list of two vectors, the y and x grid line coordinates (ascending order)
#' @param pts list of two vectors, the y and x point coordinates
#' @param sep integer vector, where `sep-1` is the number of grid lines between each subgrid line
#' @param distinct logical, indicating to snap no more than one input point to each grid point
#' @param quiet logical indicating to suppress console messages
#'
#' @return a large list containing info about grid and subgrid (see details)
#' @export
#'
#' @examples
#' # define a grid of coordinates and add jitter to a subgrid
#' gdim = c(y=100, x=100)
#' gres = c(y=2, x=2)
#' sep = c(y=5, x=8)
#' gyx = Map(\(p, r) seq(1, p, by=r), p=gdim, r=gres)
#' ng = prod(sapply(gyx, length))
#' pts = Map(\(p, r) seq(1, p, by=r), p=gdim, r=sep*gres) |>
#' expand.grid() + stats::rnorm(2*ng, 0 , 3)
#'
#' # plot grid and the point set generated from it
#' pkern_snap_plot(gyx, pts)
#'
#' # compute score for a range of y separation values, given sep_x=18
#' sep_test = seq_along(gyx[['y']])
#' g = list(gyx=gyx, gdim=gdim, gres=gres)
#' pts_list = pts |> as.data.frame() |> as.list()
#' test_score = sapply(sep_test, \(y) pkern_snap_score_2d(sep=c(y=y, x=8), g, pts_list))
#' plot(sep_test, test_score, pch=16, cex=0.5)
#' lines(sep_test, test_score)
#'
#' # estimate sep based on this score and the distance distribution among points
#' sep = pkern_estimate_sep(g, pts)
#' print(sep)
#'
#' # plot resulting snapping map
#' snap = pkern_snap(g, pts, sep)
#' pkern_snap_plot(snap, pts)
#'
#' # snap to grid lines, eliminate duplicates
#' snap.distinct = pkern_snap(g, pts, distinct=TRUE)
#' pkern_snap_plot(snap.distinct, pts)
#'
#' # another example with and more jitter
#' pts = expand.grid(gyx) + stats::rnorm(2*ng, sd=3)
#' pts = pts[sample(ng, 50),]
#' pkern_snap_plot(g, pts)
#' snap = pkern_snap(g, pts, sep=3)
#' pkern_snap_plot(snap, pts)
#'
#' # call again with smaller sep
#' snap.distinct = pkern_snap(g, pts, sep=2)
#' pkern_snap_plot(snap.distinct, pts)
#'
#' # call again with sep determined automatically
#' snap.distinct = pkern_snap(g, pts)
#' pkern_snap_plot(snap.distinct, pts)
pkern_snap = function(gyx, pts, sep=NULL, distinct=FALSE, quiet=FALSE)
{
  # initialize output list
  g = list(gres=c(1,1), sg=list())

  # default unit resolution and handle raster grid input
  spatnms = c('SpatRaster', 'RasterLayer', 'RasterStack')
  if( any( spatnms %in% class(gyx) ) )
  {
    # resolution and grid line locations
    g[['gres']] = pkern_fromRaster(gyx, 'gres')
    gyx = pkern_fromRaster(gyx, 'gyx')
  }

  # handle list output from `pkern_fromRaster` and similar
  if( !is.null( gyx[['gres']] ) ) g[['gres']] = gyx[['gres']]
  if( !is.null( gyx[['crs']] ) ) g[['crs']] = gyx[['crs']]
  if( !is.null( gyx[['gyx']] ) ) gyx = gyx[['gyx']]

  # grid dimensions
  g[['gdim']] = sapply(gyx, length)

  # expected object classes and names
  sfnm = c('sf','sfc', 'sfg')
  spnm = c('SpatialPoints', 'SpatialPointsDataFrame')
  yxnm = stats::setNames(nm=c('y', 'x'))

  # coerce various point input types and set expected column order
  if( any(sfnm %in% class(pts)) ) pts = sf::st_coordinates(pts)[,2:1]
  if( any(spnm %in% class(pts)) & requireNamespace('sp', quietly=TRUE)) pts = sp::coordinates(pts)[,2:1]
  if( is.data.frame(pts) ) pts = as.matrix(pts)
  if( is.matrix(pts) ) pts = apply(pts, 2, identity, simplify=FALSE)
  if( !all(yxnm %in% names(pts)) ) names(pts)[1:2] = yxnm

  # auto-detect sep or coerce as needed
  if( is.null(sep) ) sep = pkern_estimate_sep(modifyList(g, list(gyx=gyx)), pts)

  #if( is.null(sep) ) sep = Map(\(yx, p) pkern_estimate_sep(yx, p, distinct=F), gyx, pts)
  if( length(sep) == 1 ) sep = rep(sep, 2)
  if( !all( yxnm %in% names(sep) ) ) names(sep)[1:2] = yxnm

  # count grid lines and catch invalid requests
  msg.unbalanced = 'More input points than grid points. Try distinct=FALSE and/or decrease sep'
  npts = sapply(pts, length) |> unique()
  if( length(npts) > 1 ) stop('y and x coordinate vectors not of equal length')
  if( ( npts > prod(g[['gdim']]) ) & distinct ) stop(msg.unbalanced)

  # separately process x and y dimensions with 1d snapper,
  #snap.yx = Map(\(yx, p, s) pkern_snap_1d(yx, p, s, distinct=F), yx=gyx, p=pts, s=sep)
  #snap = Map(\(y, x) stats::setNames(list(y, x), yxnm), snap.yx[['y']], snap.yx[['x']])
  snap_yx = Map(\(yx, p, s) {

    off = pkern_est_offset(yx, p, s)
    pkern_snap_score_1d(yx, p, s, off, mapping=TRUE)

    }, yx=gyx, p=pts, s=sep)

  # reshape snapping output as list
  snap = Map(\(y, x) stats::setNames(list(y, x), yxnm), snap_yx[['y']], snap_yx[['x']])

  # to identify duplicates: convert y-x mapping to vectorized index of subgrid (y flipped)
  pmap = snap[['pmap']]
  pmap.vec = pkern_mat2vec(pmap, prod(g[['gdim']]))
  is.distinct = !any( diff(sort(pmap.vec)) == 0 )

  # handle requests for unique mapping when the dimension-wise snapping produced duplicates
  if( distinct & !is.distinct )
  {
    # define a padded extent for subgrid
    syx.range = Map(\(m, s, yx) c(max(1, min(m)-s), min(max(m)+s, length(yx))), pmap, sep, gyx)

    # define grid line coordinates and check that problem is solveable
    syx = Map(\(m, s, r) as.numeric(seq(r[1], r[2], s)), pmap, sep, syx.range)
    syx.gdim = sapply(syx, length)
    if( prod(syx.gdim) < npts ) stop(msg.unbalanced)

    # find coordinates of grid lines, then of individual (column vectorized) grid points
    yx = Map(\(gl, m) gl[m], gyx, syx)
    yx.coords = pkern_coords(yx, out='list', nosort=TRUE, quiet=TRUE)

    # find (x and y component) squared distances between subgrid points and input pts
    yx.sqdist.list = Map(\(p, yx) outer(p, yx, '-')^2, pts, yx.coords)

    # find optimal mapping by Hungarian algorithm
    result.hungarian = RcppHungarian::HungarianSolver(Reduce('+', yx.sqdist.list))
    pmap = result.hungarian[['pairs']][,2]

    # extract dimension-wise mapping to yx (a subset of the grid lines)
    idx.pmap = pkern_vec2mat(pmap, syx.gdim, out='list') |> stats::setNames(yxnm)

    # overwrite point mapping with new indices (to full set of grid lines)
    pmap = Map(\(m, sgl, gl) match(sgl[m], gl), m=idx.pmap, sgl=yx, gl=gyx)

    # recompute dimension-wise snapping distances
    snap[['dist']] = Map(\(gl, m, p) abs( p - gl[m] ), gl=yx, m=idx.pmap, p=pts)
  }

  # mapping of points to grid lines in original grid
  g[['gyx']] = gyx
  g[['pmap']] = pmap

  # list of all subgrid line numbers (indices of gyx, including empty subgrid lines)
  g[['gli']] = Map(\(m, s) seq(min(m), max(m), by=s), m=pmap, s=sep)

  # list mapping input points to y and x grid line numbers
  idx.pmap = Map(\(m,i) match(m,i), m=g[['pmap']], g[['gli']])

  # vector mapping input points to full grid data in column-vectorized order
  gmap.flipy = Map(\(gl, m) gl[m], g[['gli']], idx.pmap)
  gmap.flipy[['y']] = g[['gdim']][1] - gmap.flipy[['y']] + 1
  g[['pvec']] = pkern_mat2vec(gmap.flipy, g[['gdim']])

  # set subgrid spacing with respect to original grid, and subgrid resolution
  g[['sep']] = stats::setNames(unlist(sep), yxnm)
  g[['sg']][['gres']] = g[['sep']] * g[['gres']]

  # Euclidean snapping distance
  g[['sg']][['pdist']] = do.call(\(x,y) sqrt(x + y), snap[['d2snap']])

  # subgrid dimensions and grid line coordinates
  g[['sg']][['gdim']] = sapply(g[['gli']], length)
  g[['sg']][['gyx']] = Map(\(yx, m) yx[m], yx=gyx, m=g[['gli']])

  # vector mapping input points to subgrid data in column-vectorized order
  idx.pmap[['y']] = g[['sg']][['gdim']][1] - idx.pmap[['y']] + 1
  g[['sg']][['pvec']] = pkern_mat2vec(idx.pmap, g[['sg']][['gdim']])

  # vector mapping subgrid data to input points (inverse of pvec)
  g[['sg']][['ipvec']] = match(seq(prod(g[['sg']][['gdim']])), g[['sg']][['pvec']])

  # reorder output then return
  g[['sg']] = g[['sg']][c('gdim', 'gres', 'gyx', 'pvec', 'ipvec', 'pdist')]
  return(g[c('gdim', 'gres', 'gyx', 'sep', 'pmap', 'gli', 'pvec', 'sg')])
}


#' Compute squared snapping distances (and map) from 1d points to a 1d regular sub-grid
#'
#' This function constructs the sub-grid of `g` having spacing `sep` (ie `sep-1` grid
#' points lie between adjacent sub-grid points) and computes the sum of the squared
#' distances between elements of `pts` and their nearest (snapped) sub-grid point.
#'
#' When `mapping=TRUE` the function returns the individual squared snapping distances
#' for each point, along with information about the mapping from grid to points, in
#' a named list containing:
#'
#' gyx: numeric, the input grid-line coordinates (1d)
#' d2total: numeric, the total squared snapping distance
#' pmap: numeric vector, the element of `g` mapped to each point in `pts`
#' ngl: integer, the minimum dimension of the sub-grid containing all snapped points
#' noc: integer, the number of occupied sub-grid lines (can be less than ngl)
#'
#'
#' @param gl numeric vector, the grid line coordinate values in 1 dimension
#' @param pts numeric vector, the point coordinates in 1 dimension
#' @param sep integer vector of length 2, the sub-grid spacing indices
#' @param off positive integer, the sub-grid offset index, an element of `seq(sep)`
#' @param mapping logical, indicating to return extra information in a list
#'
#' @return numeric, the total squared snapping distance, or list (see DETAILS)
#' @export
#'
#' @examples
#' sep = 3
#' ng = 1e2
#' g = seq(ng)
#' sg = seq(sample(sep, 1), ng, by=sep)
#'
#' # sample part of the subgrid and add noise
#' npts = round( (3/4) * length(sg) )
#' pts = sort(sample(sg, npts)) + stats::rnorm(npts)
#' pkern_snap_plot(g, pts)
#'
#' # snap to subgrid and plot result
#' snap = pkern_snap_score_1d(gl=sg, pts, sep=1, mapping=TRUE)
#' pkern_snap_plot(snap, pts)
#'
#' # snap to outer grid and plot result
#' snap = pkern_snap_1d(g, pts, 1)
#' pkern_snap_plot(snap, pts)
#'
pkern_snap_score_1d = function(gl, pts, sep=1, off=NULL, mapping=FALSE)
{
  # coerce inputs to numeric
  pts = as.numeric(pts)
  gl = as.numeric(gl)

  # handle invalid sep cases
  if( sep > length(gl) ) stop('sep cannot be larger than length(gl)')
  if( sep < 1 ) stop('sep cannot be smaller than 1')

  # auto-detect offset if not supplied
  if( is.null(off) ) off = pkern_est_offset(gl, pts, sep)

  # find cross-distance matrices for pts vs subgrid lines
  dmat = abs(outer(pts, gl[seq.int(off, length(gl), by=sep)], '-'))^2

  # max.col(-d) equivalent to but faster than: apply(d, 1, which.min)
  map_dmin = max.col(-dmat)

  # extract squared snapping distances for each point, total score
  d2snap = sapply(seq_along(map_dmin), \(j) dmat[j, map_dmin[j]])
  d2total = sum(d2snap, na.rm=TRUE)
  if( !mapping ) return(d2total)

  # find the mapping to from `pts` to `gl` and compute some diagnostics
  pmap = seq.int(off, length(gl), by=sep)[ map_dmin ]
  gl_occupied = unique(pmap)
  n_gl = 1 + diff(range(gl_occupied)/sep)
  n_occupied = length(gl_occupied)

  # return in list
  return( list(gyx=gl, d2total=d2total, pmap=pmap, ngl=n_gl, noc=n_occupied, d2snap=d2snap) )
}


#' 2-d version of pkern_snap_score_1d (score function for optimizer in pkern_estimate_sep)
#'
#' Compute the total sum of squared snapping distances between the points in `pts` and
#' the grid lines of the sub-grid of `gl` with spacing `sep` (ie `sep-1` grid points lie
#' between adjacent sub-grid points) and offset `off` (ie the bottom left element of the
#' sub-grid is the `off`th x grid-line).
#'
#' This function is minimized in order to select an appropriate sub-grid spacing when
#' snapping points to a high-resolution grid.
#'
#' The tuning parameter `penalty>0` penalizes sparse grids with the aim of producing
#' simpler mappings and sub-grids of smaller dimension by adding a penalty term to the
#' squared distance sum equal to:
#'
#' `penalty` X (number of unoccupied sub-grid points) X (mean squared distance)
#'
#' set `penalty=0` to omit this term.
#'
#' `g_target` should be a list containing named elements 'gdim', the grid dimensions,
#'  and 'gyx', the grid line coordinates (as in the return value of pkern_fromRaster).
#'
#' @param sep integer vector of length 2, the sub-grid spacing indices
#' @param g_target list, the target grid configuration (see DETAILS)
#' @param pts list with named numeric vectors 'x', 'y', the point coordinates
#' @param penalty non-negative numeric, larger values penalize sparser grids
#'
#' @return numeric, the total squared snapping distance (plus penalty term)
#' @export
#'
#' @examples
#'
#' # define a grid of coordinates and add jitter to a subgrid
#' gdim = c(y=300, x=300)
#' gres = c(y=2, x=2)
#' sep = c(y=15, x=18)
#' gyx = Map(\(p, r) seq(1, p, by=r), p=gdim, r=gres)
#' ng = prod(sapply(gyx, length))
#' pts = Map(\(p, r) seq(1, p, by=r), p=gdim, r=sep*gres) |>
#'  expand.grid() + stats::rnorm(2*ng, 0 , 2)
#'  # plot grid and the point set generated from it
#'  pkern_snap_plot(gyx, pts)
#'
#'  # compute score for a range of x separation values
#'  sep_test = seq_along(gyx[['y']][-1])
#'  g = list(gyx=gyx, gdim=gdim, gres=gres)
#'  pts_list = pts |> as.data.frame() |> as.list()
#'  test_score = sapply(sep_test, \(y) pkern_snap_score_2d(sep=c(y=y, x=8), g, pts_list))
#'  plot(sep_test, test_score, pch=16, cex=0.5)
#'  lines(sep_test, test_score)
pkern_snap_score_2d = function(sep, g_target, pts, penalty=1)
{
  # handle invalid sep with worst (highest) possible score
  if( any(sep < 1) | !all( sep < g_target[['gdim']] ) ) return(Inf)

  # copy grid line locations,
  gx = g_target[['gyx']][['x']]
  gy = g_target[['gyx']][['y']]

  # compute squared snapping distances separately along x and y axes
  score_result = Map(\(g, p, s) pkern_snap_score_1d(g, p, s, mapping=TRUE),
                     g=g_target[['gyx']],
                     p=as.list(pts),
                     s=sep)

  # compute the squared distance sum total and total number of grid points
  d = sum( score_result[['x']][['d2snap']] + score_result[['y']][['d2snap']] )
  ngp = score_result[['x']][['ngl']] * score_result[['y']][['ngl']]

  # compute penalty and return score
  npts = length(pts[[1]])
  penalty_toadd = penalty * ( ngp - npts ) * ( d / npts )
  return(d + penalty_toadd)
}


#' Estimate separation distance between grid lines in 2-dimensional grids
#'
#' Estimates the subgrid resolution best matching an irregular set of points.
#' Helper function for `pkern_snap_1d`
#' TODO: flesh out documentation here and tidy inner loop in code
#'
#'
#' @return the estimated value of `sep`
#' @export
#'
#' @examples
#' sep = 25
#' ng = 10*sep
#' g = seq(ng)
#' sg = seq(sample(sep, 1), ng, by=sep)
#'
#' # sample part of the subgrid and add noise
#' npts = round( (3/4) * length(sg) )
#' pts = sort(sample(sg, npts)) + stats::rnorm(npts)
#' pkern_snap_plot(g, pts)
#'
#' # estimate separation and plot result
#' sep = pkern_estimate_sep(g, pts)
#' sep
#' snap = pkern_snap_1d(g, pts, sep)
#' pkern_snap_plot(snap, pts)
#'
#' # add duplication
#' ndupe = 10
#' pts = sapply(pts, \(p) p + stats::rnorm(ndupe)) |> c()
#'
#' # fit again
#' pkern_snap_plot(g, pts)
#' sep = pkern_estimate_sep(g, pts, distinct=FALSE)
#' sep
#' snap = pkern_snap_1d(g, pts, sep=sep, distinct=FALSE)
#' pkern_snap_plot(snap, pts)
#'
pkern_estimate_sep = function(g, pts, nmax=1e2, local_search=TRUE, penalty=1, control=NULL)
{
  # expected object classes and names
  sfnm = c('sf','sfc', 'sfg')
  spnm = c('SpatialPoints', 'SpatialPointsDataFrame')
  yxnm = stats::setNames(nm=c('y', 'x'))

  # coerce various point input types and set expected column order
  if( any(sfnm %in% class(pts)) ) pts = sf::st_coordinates(pts)[,2:1]
  if( any(spnm %in% class(pts)) & requireNamespace('sp', quietly=TRUE)) pts = sp::coordinates(pts)[,2:1]
  if( is.list(pts) ) pts = as.data.frame(pts)
  if( is.data.frame(pts) ) pts = as.matrix(pts)
  if( !all(yxnm %in% colnames(pts)) ) colnames(pts) = yxnm

  # take a sample as needed
  npts = nrow(pts)
  idx_sample = seq(npts)
  if(npts > nmax) idx_sample = sample.int(npts, size=nmax)
  pts_sample = pts[idx_sample,]
  npts_sample = length(idx_sample)

  # index for omitting 0's on diagonal of distance matrix
  idx_diag = seq(1, npts_sample^2, by=npts_sample+1)

  # compute all inter-point distances in 2d, then along each dimension
  dmat = stats::dist(pts_sample)
  ydmat = stats::dist(pts_sample[,'y'])
  xdmat = stats::dist(pts_sample[,'x'])

  # compute a quantile representing the inter-point distance of adjacent points
  p_adjacent = 4 * (npts - sqrt(npts)) / ( npts * (npts-1) )
  d_adjacent = quantile(dmat[-idx_diag], p_adjacent)
  yd_adjacent = quantile(ydmat[-idx_diag], p_adjacent/2)
  xd_adjacent = quantile(xdmat[-idx_diag], p_adjacent/2)
  # based on estimated probability of randomly drawing a point neighbour for square grid

  # handle 0-distance results
  if(yd_adjacent==0) yd_adjacent = sort(unique(ydmat[-idx_diag]))[2]
  if(xd_adjacent==0) xd_adjacent = sort(unique(xdmat[-idx_diag]))[2]

  # estimate y/x resolution ratio and inter-point distances along x, y dims
  yx_scale = yd_adjacent/xd_adjacent
  d_x = sqrt( (d_adjacent^2) / (1 + yx_scale^2) )
  d_y = yx_scale * d_x

  # set default resolution
  gres = g[['gres']]
  if( is.null(gres) ) gres = c(1,1)
  if( !all(yxnm %in% names(gres) ) ) names(gres) = yxnm

  # transform to number of grid cells and return as named vector
  sep_y = pmax(1, floor(d_x/gres['y']) )
  sep_x = pmax(1, floor(d_y/gres['x']) )
  sep = c(y=sep_y, x=sep_x)
  if( !local_search ) return(sep)

  # run further optimization
  if( is.null(control) ) control = list(maxit=5e2)
  optim_result = optim(sep, pkern_snap_score_2d,
                       g_target=g,
                       pts=as.list(as.data.frame(pts_sample)),
                       penalty=penalty,
                       control=control)

  # round to nearest valid integer
  sep_local = optim_result[['par']] |> round()
  idx_over = sep_local > g[['gdim']]
  idx_under = sep_local < c(1, 1)
  sep_local[idx_over] = g[['gdim']][idx_over]
  sep_local[idx_under] = 1
  return(sep_local)
}



#' Estimate the least-distance offset for 1D regular grid with respect to a set of points
#'
#' For a set of (possibly irregularly positioned) 1d points `pts` and a regular grid `g`,
#' the function attempts to find the integer-valued offset `off` (an element of `seq(sep)`)
#' that minimizes the snapping distance (sum) from `pts` to the subgrid of `g` with spacing
#' `sep`. The resulting subgrid points are located at `g[seq(off, length(g), by=sep)]`
#'
#' Note that `sep-1` specifies the number of grid lines (in `g`) skipped between adjacent
#' subgrid lines, ie it is an index, not a distance (which depends on the resolution of `g`).
#' Similarly, the output `off` is an offset index. However, in the algorithm, `off` is
#' treated like a continuous quantity lying in `c(1, sep)` for the purpose of sampling
#' test values, as described next:
#'
#' The algorithm is a simple iterative parameter sweep. In each iteration, the parameter space
#' (initialized to `c(1, sep)`) is sampled at `n_per_iter` uniformly spaced values, say off_i,
#' for i = 1, ..., `n_per_iter`, rounded to the nearest integer, and with duplicates omitted.
#' The best (least distance) option, off_j, is identified and a new reduced parameter space is
#' defined around it with bounds (off_{j-2}, off_{j+2}). The process is repeated until a
#' maximum `max_iter` iterations is reached, or until the set of test offsets (off_i) ceases
#' to change
#'
#' @param g numeric vector of grid point positions
#' @param pts numeric vector of point positions
#' @param sep integer specifying the spacing of grid lines
#' @param max_iter positive integer > 1, the maximum number of iterations
#' @param n_per_iter positive integer > 1, the maximum number of iterations
#'
#' @return positive integer, the offset index (an element of `seq(sep)`)
#' @export
#'
#' @examples
#' g = 1:100
#' off_test = 6
#' sep = 7
#' pts = seq(off_test, max(g), by=sep)
#' pkern_est_offset(g, pts, sep)
#'
#' # noisy example
#' pts_noisy = pts + rnorm(length(pts))
#' pkern_est_offset(g, pts_noisy, sep)
pkern_est_offset = function(g, pts, sep, max_iter=10L, n_per_iter=100L)
{
  # handle trivial requests (exhaustive search)
  if( !(sep > n_per_iter) )
  {
    # sample all offset values and return the best one
    test_out = sapply(seq(sep), \(off) pkern_snap_score_1d(g, pts, sep, off))
    return( which.min(test_out) )
  }

  # initial interval for optimizer
  min_z_new = 1
  max_z_new = sep

  # flags for exit condition
  i = 1
  is_finished = FALSE
  while(i < max_iter + 1 & !is_finished)
  {
    min_z = min_z_new
    max_z = max_z_new

    # test values are rounded and pruned to unique set
    off_test = seq(min_z, max_z, length.out=n_per_iter) |> round() |> unique()
    off_n = length(off_test)

    # evaluate scores
    test_out = sapply(off_test, \(off) pkern_snap_score_1d(g, pts, sep, off))
    idx_best = which.min(test_out)

    min_z_new = ifelse(idx_best < 3, min_z, off_test[idx_best-2])
    max_z_new = ifelse(idx_best > off_n - 2, max_z, off_test[idx_best+2])

    is_finished = (min_z == min_z_new) & (max_z == max_z_new)
    i = i + 1
  }

  return(off_test[idx_best])
}
