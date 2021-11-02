#
# pkern_snap.R
# Dean Koch, Sep 2021
# Functions for snapping points to a grid
#



#' Plots a set of 1 or 2d points and their snapped versions
#'
#' The function plots the points in `pts` over the grid lines in `g`, with the option of
#' also drawing connecting lines between points and grid intersections, according to a
#' given mapping.
#'
#' This is mainly intended for visualizing the output of `pkern_snap` and related functions.
#'
#' `pts` can be a vector of 1d positions (in which case the plot shows input order against
#' position), or a matrix whose columns are the "y" and "x" coordinates (in that order,
#' unless they are named), or a data frame (which is coerced to matrix), or a list of
#' coordinate vectors.
#'
#' `g` can be a numeric vector, or list of two (in "y", "x" order), or a list containing the
#' grid line locations (as vector or list) in element "g", and, optionally,  "map" which
#' provides a mapping for drawing connecting lines. If supplied, "map" should be an integer
#' vector (or list of two) with length(s) matching (those of) `pts`, and with elements
#' indexing grid line vectors, ie from the sequence `1:length(g$g)`
#'
#' When `g$g` but not `g$map` is supplied, it is assumed that `g$g` and `pts` have
#' the same length and ordering, and lines are drawn connecting their respective elements.
#' When `g` is a list without named element "g", it is assumed to be a list of grid line
#' vectors (and no connecting lines are drawn).
#'
#' Different classes of `pts` and `g` (and "map") can be mixed but they must have the
#' same dimensionality; ie a 2d points set `pts` must be passed with a 2d grid line object
#' `g`, and if "map" is supplied, it must be a list of two mapping vectors.
#'
#' `ppars` is a list with any of the following (optional) named elements: "gcol", "mcol",
#' "dcol", and "pcol". They are the colors of, respectively, the grid lines, the snapped grid
#' lines connecting lines, and points (passed as parameter "col" to the appropriate `graphics`
#' call); and "yaxis" is a logical that applies only in the 1d case, indicating to plot point
#' locations along the y (instead of x) axis.
#'
#' @param g numeric vector, matrix, dataframe or list, the grid line coordinates (see details)
#' @param pts numeric vector, matrix, dataframe or list, the point coordinates (see details)
#' @param ppars list, (optional) plotting parameters "gcol", "mcol", "dcol", "pcol", and "yaxis"
#'
#' @return draws a plot and returns nothing
#' @export
#'
#' @examples
#' # vector g (1d case)
#' g = seq(1, 100, by=2)
#' ng = length(g)
#' pts = g + stats::rnorm(ng)
#' pkern_snap_plot(g, pts)
#'
#' # plot again with with default snapping (assumes g and pts have same order)
#' pkern_snap_plot(g=list(g=g), pts)
#'
#' # plot again with a random mapping
#' pkern_snap_plot(list(g=g, map=sample(seq(ng), ng)), pts)
#'
#' # matrix g (2d case)
#' gdim = c(30, 25)
#' sep = c(2, 3)
#' g = Map(\(d, s) seq(1, d, by=s), d=gdim, s=sep)
#' xy = expand.grid(g)
#' pts = xy + stats::rnorm(prod(dim(xy)))
#' pkern_snap_plot(g, pts)
#'
#' # plot again with a mapping from pts to g
#' pkern_snap_plot(g=list(g=g), pts)
#'
#' # TODO: examples with pkern_snap
pkern_snap_plot = function(g, pts=NULL, ppars=list())
{
  # coerce various input types and set expected column order for 2d case
  if( any(c('sf','sfc', 'sfg') %in% class(pts) ) ) pts = sf::st_coordinates(pts)[,2:1]
  if( is.data.frame(pts) ) pts = as.matrix(pts)
  if( is.matrix(pts) ) pts = apply(pts, 2, identity, simplify=FALSE)
  yxnm = stats::setNames(nm=c('y', 'x'))

  # unpack g as list
  map = NULL
  is.1d = TRUE
  if( is.list(g) )
  {
    # unpack list element g as list (connecting line mode)
    if( !is.null( g[['g']] ) )
    {
      # extract mapping and grid line vector(s)
      map = g[['map']]
      g = g[['g']]

      # error message for sanity checks below: map should be a mapping from pts to g
      err.msg = 'mismatch in number of points and grid lines. Check map argument'

      # handle 1d and 2d cases separately
      if( is.list(g) )
      {
        # unpack mapping info and/or assign defaults
        is.1d = FALSE
        if( is.null(map) ) map = Map(\(p, s) match(p, s), p=expand.grid(g), s=g)
        if( is.null(pts) ) pts = Map(\(s, i) s[i], s=g, i=map)
        if( any( sapply(map, length) != sapply(pts, length) ) ) stop(err.msg)

      } else {

        # unpack mapping info and assign defaults for 1d case
        if( is.null(map) ) map = seq_along(g)
        if( is.null(pts) ) pts = g[map]
        if( length(map) != length(pts) ) stop(err.msg)
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
    plim = range(g, pts)

    # plot horizontal grid lines
    if( yaxis )
    {
      # scatterplot of points and grid lines with connecting lines showing snap position
      plot(x=ilim, y=plim, xlab='input order', ylab='y coordinate', pch=NA)
      graphics::abline(h=g, col=gcol)
      graphics::abline(h=g[unique(map)], col=mcol)
      if( !is.null(map) ) lapply(pseq, \(y) graphics::lines(x=rep(y,2), y=c(pts[y], g[map[y]]), col=dcol) )
      graphics::points(pseq, pts, col=pcol)

    } else {

      # same as above but with arguments reversed (vertical grid lines)
      plot(x=plim, y=ilim, xlab='x coordinate', ylab='input order', pch=NA)
      graphics::abline(v=g, col=gcol)
      graphics::abline(v=g[unique(map)], col=mcol)
      if( !is.null(map) ) lapply(pseq, \(x) graphics::lines(x=c(pts[x], g[map[x]]), y=rep(x,2), col=dcol) )
      graphics::points(pts, pseq, col=pcol)
    }

    # finish
    return( invisible() )
  }

  # set names as needed for 2d case
  if( !all( yxnm %in% names(g) ) ) names(g)[1:2] = yxnm
  if( !all( yxnm %in% names(pts) ) ) names(pts)[1:2] = yxnm
  if( !is.null(map) ) if( !all( yxnm %in% names(map) ) ) names(map)[1:2] = yxnm

  # set limits for the plot
  pseq = seq(length(pts[[1]]))
  plim = lapply(yxnm, \(d) range(g[[d]], pts[[d]]) )

  # create the plot, add grid lines
  plot(plim, pch=NA, xlab='x', ylab='y', asp=1)
  graphics::abline(h=g[['y']], v=g[['x']], col=gcol)
  graphics::abline(h=g[['y']][ unique(map[['y']]) ], v=g[['x']][ unique(map[['x']]) ], col=mcol)

  # add snapping vectors and finish
  yxvec = lapply(yxnm, \(d) lapply(pseq, \(i) c(pts[[d]][i], g[[d]][ map[[d]][i] ]) ) )
  Map(\(x, y) graphics::lines(x=x, y=y, col=dcol), x=yxvec[['x']], y=yxvec[['y']])
  graphics::points(pts, col=pcol)
  return( invisible() )
}

# TODO: flesh out documentation and simplify output
#' Snap an irregular set of 2d points to nearest regular subgrid
#'
#' Snaps points to a subgrid, either by calling `pkern_snap_1d` twice, or by using the
#' hungarian algorithm on the 2d problem
#'
#' output list contains:
#'  `gdim`, the input grid dimensions (ny, nx)
#'  `sep`, the estimated spacing (number of grid lines between subgrid lines) for the subgrid
#'  `drange`, the range of snapping distances overall (in same units as `g`)
#'  `vec`, integer vector, the mapping from `pts` to `g` in vectorized order
#'  `g`, list of integer vectors, dimension-wise coordinates of the grid lines
#'  `map`, list of integer vectors, dimension-wise mapping from `pts` to `g`
#'  `dist`, list of numeric vectors, dimension-wise snapping distances
#'  `sg`, list with information about the subgrid:
#'    `gdim`, the subgrid dimensions (ny, nx)
#'    `g`, the subgrid lines (a subset of `g` in parent list)
#'    `map`, dimension-wise mapping from `pts` to `sg$g`
#'    `vec`, mapping from `pts` to `sg$g` in vectorized order
#'
#' @param g list of two vectors ("x" and "y"), the grid line coordinates (ascending order)
#' @param pts list of two vectors ("x" and "y"), the 2d point coordinates (any order)
#' @param sep integer vector, where `sep-1` is the number of grid lines between subgrid lines
#' @param distinct logical, indicating to attempt a 1-to-1 mapping
#' @param quiet logical indicating to suppress console messages
#' @param bias nonegative numeric passed to `pkern_estimate_sep` (when `sep` not supplied)
#'
#' @return a list (see details)
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
#' pkern_snap_plot(list(g=g), pts)
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
#' snap.distinct = pkern_snap(g, pts, sep=c(2,1))
#' pkern_snap_plot(snap.distinct, pts)
#'
#' # call again with sep determined automatically
#' snap.distinct = pkern_snap(g, pts)
#' pkern_snap_plot(snap.distinct, pts)
#'
pkern_snap = function(g, pts, sep=NULL, distinct=TRUE, quiet=FALSE, bias=10)
{
  # default unit resolution
  ds = c(1,1)

  # handle raster grids
  if( 'RasterLayer' %in% class(g) )
  {
    # resolution and grid line locations
    ds = pkern_fromRaster(g, 'ds')
    g = pkern_fromRaster(g, 'yx')
  }

  # handle list input
  if( !is.null( g[['ds']] ) ) ds = g[['ds']]
  if( !is.null( g[['yx']] ) ) g = g[['yx']]
  if( !is.null( g[['g']] ) ) g = g[['g']]

  # coerce various point input types and set expected column order
  yxnm = stats::setNames(nm=c('y', 'x'))
  if( any(c('sf','sfc', 'sfg') %in% class(pts) ) ) pts = sf::st_coordinates(pts)[,2:1]
  if( any(c('SpatialPoints', 'SpatialPointsDataFrame') %in% class(pts) ) ) pts = sp::coordinates(pts)[,2:1]
  if( is.data.frame(pts) ) pts = as.matrix(pts)
  if( is.matrix(pts) ) pts = apply(pts, 2, identity, simplify=FALSE)
  if( !all( yxnm %in% names(pts) ) ) names(pts)[1:2] = yxnm

  # auto-detect sep or coerce as needed
  if( is.null(sep) ) sep = Map(\(gl, p) pkern_estimate_sep(gl, p, distinct=FALSE, bias=bias), g, pts)
  if( length(sep) == 1 ) sep = rep(sep, 2)
  if( !all( yxnm %in% names(sep) ) ) names(sep)[1:2] = yxnm

  # count grid lines and catch invalid requests
  msg.unbalanced = 'More input points than grid points. Try distinct=FALSE and/or decrease sep'
  ng = sapply(g, length)
  npts = sapply(pts, length) |> unique()
  if( length(npts) > 1 ) stop('y and x coordinate vectors not of equal length')
  if( ( npts > prod(ng) ) & distinct ) stop(msg.unbalanced)

  # separately process x and y dimensions with 1d snapper, reshape the output list
  yxsnap = Map(\(yx, p, s) pkern_snap_1d(yx, p, s, distinct=F), yx=g, p=pts, s=sep)
  snap = Map(\(y, x) stats::setNames(list(y, x), yxnm), yxsnap[['y']], yxsnap[['x']])

  # to identify duplicates, convert y-x mapping to vectorized index of subgrid (y flipped)
  map.vec = pkern_mat2vec(snap[['map']], prod(sapply(snap[['g']], length)))
  is.distinct = !any( diff(sort(map.vec)) == 0 )

  # handle requests for unique mapping when the dimension-wise snapping produced duplicates
  if( distinct & !is.distinct )
  {
    syx.range = Map(\(m, s, gl) c(max(1, min(m)-s), min(max(m)+s, length(gl))), snap[['map']], sep, g)
    syx = Map(\(m, s, r) as.numeric(seq(r[1], r[2], s)), snap[['map']], sep, syx.range)
    syx.gdim = sapply(syx, length)
    if( prod(syx.gdim) < npts ) stop(msg.unbalanced)

    # find coordinates of grid lines, then of individual (column vectorized) grid points
    yx = Map(\(gl, m) gl[m], g, syx)
    yx.coords = pkern_coords(yx, out='list', nosort=TRUE, quiet=TRUE)

    # find (x and y component) squared distances between subgrid points and input pts
    yx.sqdist.list = Map(\(gl, p) Rfast::Outer(gl, p, '-')^2, yx.coords, pts)

    # find optimal mapping by Hungarian algorithm
    result.hungarian = RcppHungarian::HungarianSolver(Reduce('+', yx.sqdist.list))
    map = result.hungarian[['pairs']][, 2]

    # extract dimension-wise mapping and overwrite values in snap
    snap[['mapidx']] = pkern_vec2mat(map, syx.gdim, out='list') |> stats::setNames(yxnm)
    snap[['map']] = Map(\(m, sgl, gl) match(sgl[m], gl), m=snap[['mapidx']], sgl=yx, gl=g)

    # recompute dimension-wise snapping distances
    snap[['dist']] = Map(\(gl, m, p) abs( p - gl[m] ), gl=yx, m=snap[['mapidx']], p=pts)
  }

  # add some summary info to output list
  snap[['gdim']] = ng
  snap[['drange']] = range( Reduce('+', snap[['dist']]) )

  # copy separation distance
  snap[['sep']] = stats::setNames(sep, yxnm)
  snap[['ds']] = mapply(\(d, s) d*s, d=ds, s=sep)

  # trim to keep only the relevant subset of and replace g with mapping
  snap[['gidx']] = Map(\(m, s) seq(min(m), max(m), by=s), m=snap[['map']], s=sep)
  #snap[['map']] = Map(\(m, i) match(m, i), m=snap[['map']], i=snap[['gidx']])
  snap[['mapidx']] = Map(\(m, i) match(m, i), m=snap[['map']], i=snap[['gidx']])
  snap[['sgdim']] = sapply(snap[['gidx']], length)

  # add mapping to vectorized form of full grid
  gmap.flipy = Map(\(gl, m) gl[m], snap[['gidx']], snap[['mapidx']])
  gmap.flipy[['y']] = ng[1] - gmap.flipy[['y']] + 1
  snap[['vec']] = pkern_mat2vec(gmap.flipy, ng)

  # ... and of subgrid
  sgmap.flipy = snap[['mapidx']]
  sgmap.flipy[['y']] = snap[['sgdim']][1] - sgmap.flipy[['y']] + 1
  snap[['sgvec']] = pkern_mat2vec(sgmap.flipy, snap[['sgdim']])

  # reorder output in return
  return(snap[c('gdim', 'sep', 'ds', 'drange', 'vec', 'g', 'map', 'gidx', 'mapidx', 'dist', 'sgdim', 'sgvec')])
}


#' Snap an irregular set of 1d points to nearest regular subgrid
#'
#' This function finds a regular subgrid of the regular 1d grid `g` that is nearest
#' to the (noisy) point locations in `pts`, where the subgrid spacing `sep` specifies
#' that every `sep`th grid point in `g` belongs to the subgrid.
#'
#' The function returns a list containing:
#'  `g` vector of grid line coordinates
#'  `map` output mapping vector, of same length as `pts` but with values in `seq(g)`
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
#' @return list with elements "g", "map", and "dist" (see below)
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
  maptog = sg.list[[ idx.best ]][ map.dmin[[idx.best]] ]

  # return in list
  return( list(g=g, map=maptog, dist=sqrt(dsnap[[ idx.best ]])) )

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
pkern_estimate_sep = function(g, pts, distinct=TRUE, bias=0, dlim=NULL)
{
  # snap the points at finest detail level `(sep=1`) and sort the mapping
  snap.result = pkern_snap_1d(g, pts, sep=1, distinct=FALSE)
  map.order = order(snap.result[['map']])
  map = snap.result[['map']][map.order]
  dist = snap.result[['dist']][map.order]

  # find differences between adjacent mappings
  nmap = length(map)
  dmap = diff(map)

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
      mapmod = lapply(septest, \(s) dmap[idx.cl] %% s )


     # septest |> head()
      #mapmod |> head()
     # Map(\(m, s) pmin(m, abs(m-s)), m=mapmod, s=septest) |> head()


      mapmod2 = Map(\(m, s) pmin(m, abs(m-s)), m=mapmod, s=septest)
      mapcost = sapply(mapmod2, stats::mad)




      #mapcost = mapply(\(m, s)  stats::median( pmin(m, abs(m-s)) ), m=mapmod, s=septest)

      #mapcost = lapply(mapmod, stats::median)


      # bias adjustment to favour sparser grids
      if( bias > 0 ) mapcost = mapcost * ( (septest-1)^(-bias^2) )


      #plot(septest, mapcost, pch=NA)
      #lines(septest, mapcost)
      #abline(v= septest[which.min(mapcost)], col='red')



      idx.best = which.min(mapcost)
    }

    # select the minimum cost option
    dsep = septest[idx.best]

  } else {

    # rounded median separation (used for nmap <= 3)
    dsep = stats::median(dmap) |> round() |> max(1)

  }

  return(dsep)
}



