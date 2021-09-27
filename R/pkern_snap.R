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
#' `g` can be a numeric (vector or matrix) or a dataframe (similar to `pts`), or a list
#' containing the grid lines as element "g", and optionally "map" which provides a mapping
#' for drawing connecting lines. If supplied, "map" should be a vector (or list of vectors)
#' with length(s) matching (those of) `pts`, and mapping to elements of `g`.
#'
#' When `g` is a list and the mapping is not supplied, it is assumed that `g` and `pts` have
#' the same length and ordering, and lines are drawn connecting their respective elements.
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
#' pts = g + rnorm(ng)
#' pkern_snap_plot(g, pts)
#'
#' # plot again with with default snapping (assumes g and pts have same order)
#' pkern_snap_plot(list(g=g), pts)
#'
#' # plot again with a random mapping
#' pkern_snap_plot(list(g=g, map=sample(seq(ng), ng)), pts)
#'
#' # matrix g (2d case)
#' gdim = c(30, 25)
#' ds = c(2, 3)
#' yx = mapply(\(g, d) seq(1, g, by=d), g=gdim, d=ds)
#' g = expand.grid(yx)
#' pts = g + rnorm(prod(dim(g)))
#' pkern_snap_plot(g, pts)
#' pkern_snap_plot(g=list(g=g), pts)
#'
#' # TODO: examples with pkern_snap
pkern_snap_plot = function(g, pts=NULL, ppars=list())
{
  # coerce various input types and set expected column order for 2d case
  if( is.data.frame(g) ) g = as.matrix(g)
  if( is.data.frame(pts) ) pts = as.matrix(pts)
  if( is.matrix(pts) ) pts = apply(pts, 2, identity, simplify=FALSE)
  yxnm = setNames(nm=c('y', 'x'))

  # unpack g as list
  if( is.list(g) )
  {
    # unpack mapping info
    map = g[['map']]
    g = g[['g']]

    # convert list in element "g" as dataframe and matrix to list
    if( is.data.frame(g) ) g = as.matrix(g)
    if( is.matrix(g) ) g = apply(g, 2, identity, simplify=FALSE)

    # unpack list element g as list
    if( is.list(g) )
    {
      # unpack mapping info and/or assign defaults
      is.1d = FALSE
      if(is.null(map)) map = lapply(g, \(p) seq_along(p))
      if(is.null(pts)) pts = mapply(\(p, i) p[i], p=g, i=map, SIMPLIFY=FALSE)

    } else {

      # unpack mapping info and assign defaults
      is.1d = TRUE
      if(is.null(map)) map = seq_along(g)
      if(is.null(pts)) pts = g[map]
    }

  } else {

    # no mapping in this case
    map = NULL

    # 2D case: g is coordinates matrix
    if( is.matrix(g) )
    {
      is.1d = FALSE
      g = apply(g, 2, identity, simplify=FALSE)

    } else { is.1d = TRUE }
  }

  # unpack plotting parameters and set defaults
  gcol = ifelse( is.null(ppars[['gcol']] ), 'grey70', ppars[['gcol']])
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
      abline(h=g, col=gcol)
      abline(h=g[unique(map)], col=mcol)
      if( !is.null(map) ) lapply(pseq, \(y) lines(x=rep(y,2), y=c(pts[y], g[map[y]]), col=dcol) )
      points(pseq, pts, col=pcol)

    } else {

      # same as above but with arguments reversed (vertical grid lines)
      plot(x=plim, y=ilim, xlab='x coordinate', ylab='input order', pch=NA)
      abline(v=g, col=gcol)
      abline(v=g[unique(map)], col=mcol)
      if( !is.null(map) ) lapply(pseq, \(x) lines(x=c(pts[x], g[map[x]]), y=rep(x,2), col=dcol) )
      points(pts, pseq, col=pcol)
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
  plot(plim, pch=NA, xlab='x', ylab='y')
  abline(h=g[['y']], v=g[['x']], col=gcol)
  abline(h=g[['y']][ unique(map[['y']]) ], v=g[['x']][ unique(map[['x']]) ], col=mcol)

  # add snapping vectors and finish
  yxvec = lapply(yxnm, \(d) lapply(pseq, \(i) c(pts[[d]][i], g[[d]][ map[[d]][i] ]) ) )
  mapply(\(x, y) lines(x=x, y=y, col=dcol), x=yxvec[['x']], y=yxvec[['y']])
  points(pts, col=pcol)
  return( invisible() )
}


#' Snap an irregular set of 1d points to nearest regular subgrid
#'
#' This function finds a regular subgrid of the regular 1d grid `g` that is nearest
#' to the (noisy) point locations in `pts`, where the subgrid spacing `sep` specifies
#' that every `sep`th grid point in `g` belongs to the subgrid.
#'
#' The output mapping (say `map`) is selected by least sum of squared snapping distances
#' (SSSD), where the snapping distance for the `i`th point is `abs( pts[i] - g[map[i]] )`.
#'
#' Note that there are `sep` possible alignments for the subgrid - eg with `sep=2`, the
#' subgrid could begin at `g[1]`, or at `g[2]`. The function checks all options by brute
#' force, returning the one that yields least SSSD.
#'
#' When `distinct==TRUE`, the function first snaps by least distance, then reassigns
#' duplicate mappings (where a grid point in `g` is mapped to more than one point in `pts`)
#' using the following heuristic algorithm:
#'
#' 1) identify all duplicates in `map` and for each compute distances to nearest empty `g`
#' 2) select the least distance `dupe`, `empty` pair (or, in case of ties, the first)
#' 3) shift all mappings between `dupe` and `empty`, towards `empty` (see `pkern_shift`)
#' 4) if `map` still contains at least one duplicate and one empty grid point, go to (1)
#'
#' The function then returns the mapping with least SSSD among those with no duplicate
#' mappings remaining. Note that when `length(g)/sep < length(pts)`, all mappings will have
#' at least one duplicate (the algorithm will terminate when it runs out of empty grid
#' points), so the function then returns the best SSSD solution among those with the
#' least number of duplicates.
#'
#' @param g vector of grid line locations in ascending order
#' @param pts vector of point locations in 1D
#' @param sep integer
#' @param distinct logical, indicating to attempt a 1-to-1 mapping
#'
#' @return integer vector, mapping `pts` to `x`
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
#' pts = sort(sample(sg, npts)) + rnorm(npts)
#' pkern_snap_plot(g, pts)
#'
#' # snap to nearest subgrid (allowing duplicates) and plot result
#' map.list = pkern_snap_1d(g, pts, sep, distinct=FALSE)
#' pkern_snap_plot(map.list, pts)
#'
#' # repeat but avoid duplicate mappings (default `distinct=TRUE`)
#' map.list = pkern_snap_1d(g, pts, sep)
#' pkern_snap_plot(map.list, pts)
#'
pkern_snap_1d = function(g, pts, sep=1, distinct=TRUE)
{
  # coerce inputs to numeric and handle trivial case of single grid line
  pts = as.numeric(pts)
  g = as.numeric(g)
  ng = length(g)
  npts = length(pts)
  if( length(g) == 1 ) return( rep(1, npts) )

  # establish mapping to and from ascending order and sort `pts`
  idx.order = order(pts)
  idx.unorder = match(seq_along(pts), idx.order)
  pts = pts[idx.order]

  # list of candidate subgrid lines, one for each possible origin
  sg.list = lapply(seq(sep), function(x) seq(x, length(g), by=sep))
  sg.n = sapply(sg.list, length)

  # find cross-distance matrices for pts vs subgrid lines and minimum distance mapping
  dmat.list = lapply(sg.list, \(i) abs(Rfast::Outer(pts, g[i], '-'))^2 )
  map.dmin = lapply(dmat.list, \(d) apply(d, 2, which.min) )

  # for each mapping, identify empty and duplicate elements
  map.initial = mapply(\(m, sg) list(map=m, g=sg), m=map.dmin, sg=sg.list, SIMPLIFY=FALSE)
  table.initial = lapply(map.initial, \(map) pkern_shift(map, shift=FALSE))

  # identify candidates with duplicates and those that can be remedied
  is.dupe = sapply(table.initial, \(ini) length(ini[['dupe']]) > 0)
  must.dupe = npts > sg.n
  can.fix = is.dupe & (!must.dupe)

  # remap duplicate mappings to empty grid lines, if requested
  table.final = table.initial
  if( distinct & any(can.fix) )
  {
    # loop over fixable grids
    for( idx.sg in which(can.fix) )
    {
      can.shift = TRUE
      while( can.shift )
      {
        # unpack mapping table and coerce to numeric (class required by Rfast)
        map = table.final[[idx.sg]][['map']]
        empty = as.numeric( table.final[[idx.sg]][['empty']] )
        dupe = as.numeric( table.final[[idx.sg]][['dupe']] )

        # apply the shift when it is possible
        can.shift = ( length(empty) * length(dupe) ) > 0
        if( can.shift )
        {
          # identify nearest empty cell to shift towards
          empty.tofill = empty[ apply(abs(Rfast::Outer(dupe, empty, '-')), 2, which.min) ]
          dupe.tofix = dupe[ which.min( abs( dupe - empty.tofill ) ) ]

          # perform the shift and, overwriting mapping in table list
          table.final[[idx.sg]] = pkern_shift(map, g=sg.list[[idx.sg]], dupe.tofix)
        }
      }
    }
  }

  # compute sums of squared snap distance for each candidate
  dfinal = mapply(\(tab, d) sum(sapply(seq(length(tab[['map']])), \(j) d[tab[['map']][j],j])),
                  tab=table.final, d=dmat.list)

  # find minimum number of duplicates among the candidates
  ndupe = sapply(table.final, \(tab) length( tab[['dupe']] ))
  is.mindupe = ndupe == min(ndupe)

  # extract the candidate subgrid with minimum total distance
  if( distinct ) dfinal[ which(!is.mindupe) ] = Inf
  idx.best = which.min(dfinal)
  table.best = table.final[[ idx.best ]]

  # find its mapping to `g`, then reorder this mapping to match input `pts` order
  map.sorted = sg.list[[idx.best]][ table.best[['map']] ]
  table.best[['map']] = table.best[['map']][ idx.unorder ]

  # replace grid line numbers with grid line coordinates
  table.best[['g']] = g[ table.best[['g']] ]

  # reorder this mapping to match input `pts` order
  return( table.best )
}



#' Snap an irregular set of 2d points to nearest regular subgrid
#'
#' @param g
#' @param pts
#' @param sep
#'
#' @return
#' @export
#'
#' @examples
pkern_snap_2d = function(g, pts, sep=1)
{





}

#' Shift operation for discrete mapping (helper for `pkern_snap_1d`)
#'
#' This function identifies the unassigned integer in `g` (ie not appearing in `map`)
#' that is closest to the element `dupe` in `g[map]`. All mappings between this element
#' (call it `empty`) and `dupe` are then incremented towards `empty`.
#'
#' When `dupe` appears more than once in `g[map]`, the effect of the function is to
#' replace all but one of the duplicate mappings by shifting them - and any intervening
#' mappings - towards the nearest empty grid point number in `g`. This reassignment
#' results in exactly one point mapping to `empty`, and it eliminates `dupe` as a
#' duplicate.
#'
#' The function returns a list with elements:
#'
#'  "map", the new mapping vector (of same length as `map`)
#'  "g", grid line numbers for the outer grid (same as input `g`, if provided)
#'  "empty", the new indices of "g" values not appearing in "map"
#'  "dupe", the new duplicate "map" values (if any)
#'
#' Note that when `dupe` appears > 2 times in `map`, this operation generates a new,
#' smaller duplicate set in the output map with value `dupe +/- 1`. Since the net effect
#' is to reduce the number of duplicates by 1, the function can be run recursively to
#' eliminate all duplicates, provided that `length(map) < length(g)`
#'
#' For convenience, arguments `map`, `g`, and `dupe` can be passed together as a list
#' to the function's first argument `map`. If `g` is missing it is set to the integer
#' sequence from 1 to `max(map)`, and if `dupe` is missing it is set to the first
#' duplicate element in `map`.
#'
#' @param map integer vector, an indexing of `g`
#' @param g numeric vector, grid line locations in ascending order
#' @param dupe integer, the duplicate `x[map]` value to amend
#' @param shift logical, if FALSE, the function returns info about input mapping `map`
#'
#' @return list with elements "map", "g", "empty", "dupe"
#' @export
#'
#' @examples
#'
#' # create a mapping that contains at least one duplicate
#' g = seq(1, 20, by=2)
#' pts = sample(g, 11, replace=TRUE) |> sort()
#' map = match(pts, g)
#' pkern_snap_plot(list(g=g, map=map))
#'
#' # extract info about duplication
#' result0 = pkern_shift(map, g, shift=FALSE)
#' print(result0)
#'
#' # apply the shift operation twice
#' result1 = pkern_shift(result0)
#' result2 = pkern_shift(result1)
#'
#' # plot the reassignment done at each stage
#' par(mfrow=c(2,1))
#' pkern_snap_plot(list(g=g, map=result1$map), g[map])
#' pkern_snap_plot(list(g=g, map=result2$map), g[result1$map])
#' par(mfrow=c(1,1))
#'
#' # run the algorithm until there are no more empty slots
#' result = result0
#' while( ( length( result$empty ) > 0 ) & ( length( result$dupe ) > 0 ) ) {
#' result = pkern_shift(result)
#' }
#'
#' # plot the result
#' par(mfrow=c(3,1))
#' pkern_snap_plot(list(g=g, map=map))
#' pkern_snap_plot(result, g[map])
#' pkern_snap_plot(result, g[result$map])
pkern_shift = function(map, g=NULL, dupe=NULL, shift=TRUE)
{
  # unpack list input
  if( is.list(map) )
  {
    g = map[['g']]
    dupe = map[['dupe']]
    map = map[['map']]
  }

  # set defaults as needed
  if( is.null(g) ) g = seq( max(map) )
  if( is.null(dupe) ) dupe = g[ unique( map[ which( diff(map) == 0 ) ] ) ]

  # if dupe is a vector, use its first element
  if( length(dupe) > 1 ) dupe = dupe[1]

  # find empty grid points
  empty = g[ which( !( seq_along(g) %in% map ) ) ]

  # general sanity check
  if( !all( map %in% seq_along(g) ) ) stop('at least one element of map not found in g')

  # sanity checks for shift (note that dupe is ignored when shift=FALSE)
  if( shift )
  {
    if( length(empty) == 0 ) stop('all of g is already mapped in map')
    if( length(dupe) == 0 ) stop('no duplicates in map')
    if( !( dupe %in% g[map] ) ) stop('dupe was not found in map')
  }

  # overwrite `map` in ascending order but save inverse ordering
  idx.order = order(map)
  idx.unorder = match(seq_along(map), idx.order)
  map = map[idx.order]

  # apply shift operation if requested
  if( shift )
  {
    # size of the duplicate set
    n.dupe = sum(g[map]==dupe)

    # identify empty cell to shift towards, and shift direction
    empty.tofill = empty[ which.min( abs(dupe - empty) ) ]
    int.shift = (2 * as.integer(dupe < empty.tofill) ) - 1

    # identify set of points to modify in `map` and their indices
    toshift = dupe:(empty.tofill-int.shift)
    idx.toshift = which(g[map] %in% toshift)

    # select one of the duplicate mappings to remain unchanged
    idx.keep = match(dupe, g[map]) + as.integer( int.shift < 0 ) * (n.dupe - 1)
    idx.toshift = idx.toshift[ idx.toshift != idx.keep ]

    # shift the indices
    map[idx.toshift] = map[idx.toshift] + int.shift

  }

  # compute output vectors
  list.out = list(g = g,
                  map = map[idx.unorder],
                  dupe = unique( g[ map[ which( diff(map) == 0 ) ] ] ),
                  empty = g[ which( !( seq_along(g) %in% map ) ) ] )

  # finish
  return(list.out)
}



