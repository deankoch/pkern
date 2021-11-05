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
#' call); and "yaxis" is a logical that applies only in the 1d case, indicating to plot point
#' locations along the y (instead of x) axis.
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
#' g = seq(1, 100, by=2)
#' ng = length(g)
#' pts = g + stats::rnorm(ng)
#' pkern_snap_plot(g, pts)
#'
#' # plot again with with default snapping (assumes g and pts have same order)
#' pkern_snap_plot(gyx=list(gyx=g), pts)
#'
#' # plot again with a random mapping
#' pkern_snap_plot(list(gyx=g, pmap=sample(seq(ng), ng)), pts)
#'
#' # 2d case
#' gdim = c(30, 25)
#' sep = c(2, 3)
#' gyx = Map(\(d, s) seq(1, d, by=s), d=gdim, s=sep)
#' yx = expand.grid(gyx)
#' pts = xy + stats::rnorm(prod(dim(xy)))
#' pkern_snap_plot(gyx, pts)
#'
#' # plot again with a mapping from pts to gyx
#' pkern_snap_plot(gyx=list(gyx=gyx), pts)
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
    # unpack list element gyx as list (connecting line mode)
    if( !is.null( gyx[['gyx']] ) )
    {
      # extract mapping and grid line vector(s)
      pmap = gyx[['pmap']]
      gyx = gyx[['gyx']]

      # error message for sanity checks below: pmap should be a mapping from pts to gyx
      err.msg = 'mismatch in number of points and grid lines. Check pmap argument'

      # handle 1d and 2d cases separately
      if( is.list(gyx) )
      {
        # unpack mapping info and/or assign defaults
        is.1d = FALSE
        if( is.null(pmap) ) pmap = Map(\(p, s) match(p, s), p=expand.grid(gyx), s=gyx)
        if( is.null(pts) ) pts = Map(\(s, i) s[i], s=gyx, i=pmap)
        if( any( sapply(pmap, length) != sapply(pts, length) ) ) stop(err.msg)

      } else {

        # unpack mapping info and assign defaults for 1d case
        if( is.null(pmap) ) pmap = seq_along(gyx)
        if( is.null(pts) ) pts = gyx[pmap]
        if( length(pmap) != length(pts) ) stop(err.msg)
      }

    } else { is.1d = FALSE }
  }

  # unpack plotting parameters and set defaults
  gcol = ifelse( is.null(ppars[['gcol']] ), NA, ppars[['gcol']])
  mcol = ifelse( is.null(ppars[['mcol']] ), 'grey40', ppars[['mcol']])
  dcol = ifelse( is.null(ppars[['dcol']] ), 'black', ppars[['dcol']])
  pcol = ifelse( is.null(ppars[['pcol']] ), 'grey50', ppars[['pcol']])
  yaxis = ifelse( is.null(ppars[['yaxis']] ), FALSE, ppars[['yaxis']])

  # handle 1D case
  if( is.1d )
  {
    # set limits for the plot
    pseq = seq_along(pts)
    ilim = range(pseq)
    plim = range(gyx, pts)

    # plot horizontal grid lines
    if( yaxis )
    {
      # scatterplot of points and grid lines with connecting lines showing snap position
      plot(x=ilim, y=plim, xlab='input order', ylab='y coordinate', pch=NA)
      graphics::abline(h=gyx, col=gcol)
      graphics::abline(h=gyx[unique(pmap)], col=mcol)
      if( !is.null(pmap) ) lapply(pseq, \(y) graphics::lines(x=rep(y,2), y=c(pts[y], gyx[pmap[y]]), col=dcol) )
      graphics::points(pseq, pts, col=pcol)

    } else {

      # same as above but with arguments reversed (vertical grid lines)
      plot(x=plim, y=ilim, xlab='x coordinate', ylab='input order', pch=NA)
      graphics::abline(v=gyx, col=gcol)
      graphics::abline(v=gyx[unique(pmap)], col=mcol)
      if( !is.null(pmap) ) lapply(pseq, \(x) graphics::lines(x=c(pts[x], gyx[pmap[x]]), y=rep(x,2), col=dcol) )
      graphics::points(pts, pseq, col=pcol)
    }

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

  # create the plot, add grid lines
  plot(plim, pch=NA, xlab='x', ylab='y', asp=1)
  graphics::abline(h=gyx[['y']], v=gyx[['x']], col=gcol)
  graphics::abline(h=gyx[['y']][ unique(pmap[['y']]) ], v=gyx[['x']][ unique(pmap[['x']]) ], col=mcol)

  # add snapping vectors and finish
  yxvec = lapply(yxnm, \(d) lapply(pseq, \(i) c(pts[[d]][i], gyx[[d]][ pmap[[d]][i] ]) ) )
  Map(\(x, y) graphics::lines(x=x, y=y, col=dcol), x=yxvec[['x']], y=yxvec[['y']])
  graphics::points(pts, col=pcol)
  return( invisible() )
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
#' nsep = c(y=5, x=8)
#' g = Map(\(p, s) seq(1, p, by=s), p=gdim, s=nsep)
#' ng = prod(sapply(g, length))
#' pts = expand.grid(g) + stats::rnorm(2*ng)
#'
#' # plot grid and the point set generated from it
#' pkern_snap_plot(g, pts)
#'
#' # snap to grid lines, allowing duplicates
#' snap = pkern_snap(g, pts, distinct=FALSE)
#' pkern_snap_plot(snap, pts)
#'
#' # snap to grid lines, eliminate duplicates
#' snap.distinct = pkern_snap(g, pts)
#' pkern_snap_plot(snap.distinct, pts)
#'
#' # another example with and more jitter and too many points
#' pts = expand.grid(g) + stats::rnorm(2*ng, sd=3)
#' pts = pts[sample(ng, 50),]
#' pkern_snap_plot(g, pts)
#' snap = pkern_snap(g, pts, sep=3, distinct=FALSE)
#' pkern_snap_plot(snap, pts)
#'
#' # uncomment this line to demonstrate error handling
#' #snap.distinct = pkern_snap(g, pts, sep=3)
#'
#' # call again with smaller sep
#' snap.distinct = pkern_snap(g, pts, sep=2)
#' pkern_snap_plot(snap.distinct, pts)
#'
#' # call again with sep determined automatically
#' snap.distinct = pkern_snap(g, pts)
#' pkern_snap_plot(snap.distinct, pts)
#'
pkern_snap = function(gyx, pts, sep=NULL, distinct=TRUE, quiet=FALSE)
{
  # initialize output list
  g = list(gres=c(1,1), sg=list())

  # default unit resolution and handle raster grid input
  if( 'RasterLayer' %in% class(gyx) )
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
  if( any(spnm %in% class(pts)) ) pts = sp::coordinates(pts)[,2:1]
  if( is.data.frame(pts) ) pts = as.matrix(pts)
  if( is.matrix(pts) ) pts = apply(pts, 2, identity, simplify=FALSE)
  if( !all(yxnm %in% names(pts)) ) names(pts)[1:2] = yxnm

  # auto-detect sep or coerce as needed
  if( is.null(sep) ) sep = Map(\(yx, p) pkern_estimate_sep(yx, p, distinct=F), gyx, pts)
  if( length(sep) == 1 ) sep = rep(sep, 2)
  if( !all( yxnm %in% names(sep) ) ) names(sep)[1:2] = yxnm

  # count grid lines and catch invalid requests
  msg.unbalanced = 'More input points than grid points. Try distinct=FALSE and/or decrease sep'
  npts = sapply(pts, length) |> unique()
  if( length(npts) > 1 ) stop('y and x coordinate vectors not of equal length')
  if( ( npts > prod(g[['gdim']]) ) & distinct ) stop(msg.unbalanced)

  # separately process x and y dimensions with 1d snapper, reshape the output as list
  snap.yx = Map(\(yx, p, s) pkern_snap_1d(yx, p, s, distinct=F), yx=gyx, p=pts, s=sep)
  snap = Map(\(y, x) stats::setNames(list(y, x), yxnm), snap.yx[['y']], snap.yx[['x']])

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
  g[['sg']][['pdist']] = do.call(\(x,y) sqrt(x^2 + y^2), snap[['dist']])

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


#' Snap an irregular set of 1d points to nearest regular subgrid
#'
#' This function finds a regular subgrid of the regular 1d grid `g` that is nearest
#' to the (noisy) point locations in `pts`, where the subgrid spacing `sep` specifies
#' that every `sep`th grid point in `g` belongs to the subgrid.
#'
#' The function returns a list containing:
#'  `gyx` vector of grid line coordinates
#'  `pmap` output mapping vector, of same length as `pts` but with values in `seq(g)`
#'  `dist` vector of squared snapping distances (corresponding to `pts`)
#'
#' The output mapping (`map`) is selected by least total cost (sum over all input points),
#' where the cost of a point is equal to the squared distance between the point and its
#' snapped (on-grid) location.
#'
#' To align the subgrid, the function calculates the cost of each possible offset and
#' selects the least cost option. eg. if `sep=2`, the least coordinate in the subgrid
#' can be set to `g[1]` (offset 0) or `g[2]` (offset 1). In general there are `sep`
#' possible offsets, which are all tested in a loop.
#'
#' When `distinct==TRUE`, the function calls an implementation of the Hungarian algorithm
#' (from R package `RcppHungarian`) to find an optimal assignment. Note that this can be
#' slow for large grids or large numbers of points.
#'
#' @param g vector of grid line locations in ascending order
#' @param pts vector of point locations in 1D
#' @param sep integer, where `sep-1` is the number of grid lines between subgrid lines
#' @param distinct logical, indicating to attempt a 1-to-1 mapping
#'
#' @return list with elements "g", "pmap", and "dist" (see below)
#' @export
#'
#' @examples
#' # define a grid and subgrid (every third grid line)
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
#' snap = pkern_snap_1d(sg, pts, sep, distinct=FALSE)
#' pkern_snap_plot(snap, pts)
#'
#' # snap to outer grid and plot result
#' snap = pkern_snap_1d(g, pts, sep, distinct=FALSE)
#' pkern_snap_plot(snap, pts)
#'
#' # repeat but avoid duplicate mappings (default `distinct=TRUE`)
#' snap = pkern_snap_1d(g, pts, sep)
#' pkern_snap_plot(snap, pts)
#'
pkern_snap_1d = function(g, pts, sep=1, distinct=TRUE)
{
  # coerce inputs to numeric and handle trivial case of single grid line
  pts = as.numeric(pts)
  g = as.numeric(g)
  ng = length(g)
  npts = length(pts)
  if( length(g) == 1 )  return(g=g, map=rep(1, npts), dist=g[[1]]-pts)

  # list of candidate subgrid lines, one for each possible origin (grid line numbers)
  sg.list = lapply(seq(sep), function(x) seq(x, length(g), by=sep))
  sg.n = sapply(sg.list, length)

  # find cross-distance matrices for pts vs subgrid lines
  dmat.list = lapply(sg.list, \(i) abs(Rfast::Outer(g[i], pts, '-'))^2 )

  # allow duplication in mapping or not?
  if( distinct )
  {
    # without duplication: Hungarian algorithm solves the unbalanced assignment problem
    result.hungarian = lapply(dmat.list, RcppHungarian::HungarianSolver)
    map.dmin = lapply(result.hungarian, \(m) m[['pairs']][, 2])

    # replace 0's (unassigned) with NA
    for(i in seq_along(map.dmin) ) { map.dmin[[i]][ map.dmin[[i]] == 0 ] = NA }

  } else {

    # with duplication: snap to nearest grid line
    map.dmin = lapply(dmat.list, \(d) apply(d, 1, which.min) )
  }

  # extract squared snapping distances for each point, total score for each candidate
  dsnap = Map(\(d, m) sapply(seq_along(m), \(j) d[j, m[j]]), d=dmat.list, m=map.dmin)
  dtotal = sapply(dsnap, \(d) sum(d, na.rm=TRUE))

  # find number of unassigned points and identify candidates having the fewest...
  norphan = sapply(map.dmin, \(m) sum(is.na(m)))
  is.fewest = norphan == min(norphan)

  # ...other candidates are scored so they never get selected
  if( distinct ) dtotal[ which(!is.fewest) ] = Inf

  # identify the eligible candidate subgrid with minimum total distance
  idx.best = which.min(dtotal)

  # find its mapping to `g`
  ptog = sg.list[[ idx.best ]][ map.dmin[[idx.best]] ]

  # return in list
  return( list(gyx=g, pmap=ptog, dist=sqrt(dsnap[[ idx.best ]])) )

}



#' Estimate separation distance between grid lines
#'
#' Estimates the subgrid resolution best matching an irregular set of points.
#' Helper function for `pkern_snap_1d`
#' TODO: flesh out documentation here and tidy inner loop in code
#'
#' @param g vector of grid line locations in ascending order
#' @param pts vector of point locations in 1D
#' @param distinct logical, indicating to look for 1-to-1 mappings
#' @param bias nonegative numeric, a tuning parameter
#' @param dlim length-2 positive numeric, lower/upper bounds for separation distance
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
pkern_estimate_sep = function(g, pts, distinct=TRUE, bias=10, dlim=NULL)
{
  # snap the points at finest detail level `(sep=1`) and sort the mapping
  snap.result = pkern_snap_1d(g, pts, sep=1, distinct=FALSE)
  pmap.order = order(snap.result[['pmap']])
  pmap = snap.result[['pmap']][pmap.order]
  dist = snap.result[['dist']][pmap.order]

  # find differences between adjacent mappings
  nmap = length(pmap)
  dmap = diff(pmap)

  # 2-means clustering
  if( nmap > 3 )
  {
    # identify the cluster representing among-grid-line differences
    gcl = stats::kmeans(dmap, 2, nstart=100)
    idx.cl = gcl[['cluster']] == which.max(gcl[['centers']])
    if(distinct) idx.cl = !idx.cl

    # find grid resolution and set default bounds for distance as needed
    gres = diff(g[1:2])
    if(is.null(dlim))
    {
      # guess distance limits based on input points - these won't be appropriate in every situation
      if(distinct) dlim = min(stats::dist(pts)) * c(1, 2)
      if(!distinct) dlim = range(dmap[idx.cl]) * gres
    }

    # find grid resolution and convert dlim to (integer) separation range to test
    seplim = pmax(round(dlim/gres), 1)
    if(distinct) seplim[2] = min(seplim[2], ceiling(length(g)/nmap))
    septest = seq(seplim[1], seplim[2])

    # identify the best separation
    idx.best = 1
    if(length(septest) > 1)
    {
      # express pointwise differences mod sep, compute score allowing some noise around sep
      pmapmod = lapply(septest, \(s) dmap[idx.cl] %% s )


     # septest |> head()
      #mapmod |> head()
     # Map(\(m, s) pmin(m, abs(m-s)), m=mapmod, s=septest) |> head()


      pmapmod2 = Map(\(m, s) pmin(m, abs(m-s)), m=pmapmod, s=septest)
      pmapcost = sapply(pmapmod2, stats::mad)




      #pmapcost = mapply(\(m, s)  stats::median( pmin(m, abs(m-s)) ), m=pmapmod, s=septest)

      #pmapcost = lapply(pmapmod, stats::median)


      # bias adjustment to favour sparser grids
      if( bias > 0 ) pmapcost = pmapcost * ( (septest-1)^(-bias^2) )


      #plot(septest, pmapcost, pch=NA)
      #lines(septest, pmapcost)
      #abline(v= septest[which.min(pmapcost)], col='red')



      idx.best = which.min(pmapcost)
    }

    # select the minimum cost option
    dsep = septest[idx.best]

  } else {

    # rounded median separation (used for nmap <= 3)
    dsep = stats::median(dmap) |> round() |> max(1)

  }

  return(dsep)
}



