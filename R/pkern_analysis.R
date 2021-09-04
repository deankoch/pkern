#
# pkern_analysis.R
# Dean Koch, Sep 2021
# development code and garbage heap
#


#' Cluster 1-dimensional coordinates onto grid lines of unknown position and number
#'
#' This function uses heuristics to find a best fitting irregular grid for the values
#' in `x`, with up to `nmax` grid lines and minimum separation distance `dmin`. It
#' does a brute force search of all `1:nmax` candidates for the number of grid lines,
#' using kmeans clustering (with `nstart` randomly assigned starting values) to find
#' the k grid line positions. Tuning parameter `cw` can be increased to favor simpler
#' grids.
#'
#' Different values of k are ranked according to two cost functions
#'
#' E(k) = sum of squared distances between points and their grid lines
#' C(k) = 1 + sum of squared differences between the expected and actual occupancy of grid line
#'
#' where if grid line i has ni points assigned, its actual occupancy is (n_i) / (n/k), and
#' the "expected" occupancy is that of a complete grid, or n/k. C(k) will be higher for
#' incomplete grids, particularly ones where high k causes grid lines to be assigned to
#' solitary points. The score function for k is the sum:
#'
#' S(k) = k + log( E(k) ) + `cw` * log( C(k) )
#'
#' This is evaluated for all k = 1,...`nmax` and the grid line set with smallest S(k) is
#' returned. If `nmax` is not supplied, it is set to the smaller of: n, or `4*sqrt(n)`,
#' or the number of unique coordinate positions. If `dmin` is not supplied, it is set
#' to `(range(x)/n)/4`
#'
#' @param x numeric vector, the positions of points to cluster
#' @param nmax positive integer, the maximum desired number of x grid lines
#' @param dmin numeric, lower limit on distance between grid lines
#' @param nstart positive integer, passed to stats:kmeans
#' @param cw non-negative numeric, the complexity exponent (see details)
#'
#' @return names list of: grid line positions "grid", and "input" their index in x
#' @export
#'
#' @examples
#' nx = 35
#' ny = 24
#' coords = expand.grid(seq(nx), seq(ny)) + rnorm(nx*ny, 0, 1/5)
#' plot(coords)
#' n = nrow(coords)
#' x = coords[,1] %>% as.vector
#' pkern_kmeans(x)
pkern_kmeans = function(x, nmax=NA, dmin=NA, nstart=25, cw=1)
{
  # compute some basic stats
  n = length(x)
  unx = unique(x)
  rx = range(unx)

  # default for minimum grid line separation distance
  if( is.na(dmin) ) dmin = (diff(rx) / n) / 4

  # default for maximum number of grid lines
  if( is.na(nmax) ) nmax = min(length(unx), round(4*sqrt(n)), n-1)

  # for allowable dimension size, find candidate grid line positions by kmeans clustering
  nx.test = seq(nmax)
  nx.cluster = lapply(nx.test, function(nx) stats::kmeans(x, nx, nstart=nstart) )

  # identify clusterings that do not violate the distance minimum
  nx.eligible = sapply(nx.cluster, function(nx) all( diff( sort( nx$centers ) ) > dmin ) )

  # catch dmin too high producing no eligible grid lines, then remove offenders
  if( !any(nx.eligible) ) stop('Found no eligible grids. Try lowering dmin')
  nx.test = nx.test[nx.eligible]
  nx.cluster = nx.cluster[nx.eligible]

  # compute a complexity score
  nx.cscore = sapply(nx.cluster, function(nx) sum( abs( 1-( length(nx$size) * nx$size / n ) )^2 ) )

  # compute net distance score for each set of grid lines
  nx.dscore = sapply(nx.cluster, function(nx) nx$tot.withinss)

  # compute exponential of overall score (lower is better) and identify best grid line set
  nx.score = nx.dscore * ( (1 + nx.cscore )^cw ) * seq_along(nx.test)
  idx.best = which.min( nx.score )
  nx = nx.test[ idx.best ]

  # find resulting grid line positions and map to input points
  gx = nx.cluster[[idx.best]]$centers
  idx = nx.cluster[[idx.best]]$cluster

  # return sorted grid line positions and the index mapping to input points
  return( list(grid=gx[order(gx)], input=match(idx, order(gx)) ) )
}


#' Ad-hoc method for snapping a set of points to the nearest grid
#'
#' This function finds a (possibly irregular) grid that most closely matches the matrix of
#' n points in `coords`. These points may form an incomplete subset of the grid. `coords`
#' should either be supported by `sf::st_coordinates`, or else coercible to a matrix with
#' column names 'x' (or 'X'), and 'y' (or 'Y'). Matrices with unnamed columns are assumed
#' to have x coordinates in the first column and y coordinates in the second.
#'
#' The algorithm works by using kmeans to identify candidate x and y grid lines, and
#' ranking them by a cost function that incorporates squared distances and a complexity
#' term. This is computed independently along each dimension, for a range of k. See
#' `?pkern_kmeans` for details on this process, and the tuning parameters `nmax`, `dmin`,
#' `nstart`, and `cw`. If a single number is supplied to any of these parameters, it is
#' duplicated for the x and y component analysis.
#'
#' The return list contains: the detected dimensions ("dims"), the x and y grid line
#' positions ("gxy"), and a dataframe containing the original and snapped coordinates
#' ("x", "xnew", "y", 'ynew") in column-vectorized order, along with their row and column
#' numbers ("i", "j"), the snapping distance ("dist"), and an indexing to recover the
#' original order of input coordinates ("input").
#'
#' @param coords n X 2 matrix of "x" and "y" coordinates or sf object
#' @param nmax length-2 integer vector, the maximum allowable number of x and y grid lines
#' @param dmin length-2 numeric vector, minimum x and y distances allowed between grid lines
#' @param nstart length-2 integer vector, passed to stats::kmeans
#' @param cw length-2 numeric vector, complexity weights for x and y, passed to pkern_kmeans
#'
#' @return named list of three elements "dims", "gxy", and "snap" (see details)
#' @export
#'
#' @examples
#' nx = 35
#' ny = 24
#' coords = expand.grid(seq(nx), seq(ny)) + rnorm(nx*ny, 0, 1/5)
#' test.grid = pkern_detect_grid(coords)
#' plot(coords)
#' abline(v = test.grid$gxy$x)
#' abline(h = test.grid$gxy$y)
pkern_detect_grid = function(coords, nmax=NA, dmin=NA, nstart=25, cw=1)
{
  # duplicate length-1 arguments to tuning parameters
  if( length(nmax) == 1 ) nmax = rep(nmax, 2)
  if( length(dmin) == 1 ) dmin = rep(dmin, 2)
  if( length(nstart) == 1 ) nstart = rep(nstart, 2)
  if( length(cw) == 1 ) cw = rep(cw, 2)

  # convert sf type input
  if( any( c('sf', 'sfc', 'sfg') %in% class(coords) ) ) coords = sf::st_coordinates(coords)

  # coerce to a two-column matrix and handle uppercase dimension names
  coords = as.matrix(coords)
  colnames(coords) = colnames(coords) |> tolower()

  # assume order x, y unless columns are named otherwise
  idx.col = match(c('x', 'y'), colnames(coords))
  if( is.na(idx.col[1]) ) idx.col[1] = 1
  if( is.na(idx.col[2]) ) idx.col[2] = ifelse(idx.col[1] == 1, 2, 1)

  # extract the input coordinate values and compute some basic stats
  x = coords[,idx.col[1]] %>% as.vector
  y = coords[,idx.col[2]] %>% as.vector

  # find best grid line positions
  gx.list = pkern_kmeans(x, nmax=nmax[1], dmin=dmin[1], nstart=nstart[1], cw=cw[1])
  gy.list = pkern_kmeans(y, nmax=nmax[2], dmin=dmin[2], nstart=nstart[2], cw=cw[2])

  # unpack output to get grid line assignments (ordered according to input)
  j = gx.list[['input']]
  i = gy.list[['input']]
  gx = gx.list[['grid']]
  gy = gy.list[['grid']]

  # build a dataframe of snapped coordinates and their i,j indices in the grid
  coords.snap = data.frame(input=seq_along(x), i=i, j=j, x=x, xnew=gx[j], y=y, ynew=gy[i])
  rownames(coords.snap) = NULL

  # add the snapping distance (error)
  coords.snap$dist = sqrt( (coords.snap$xnew - x)^2 + (coords.snap$ynew - y)^2 )

  # reorder to column-vectorized form, permuting the i indices to order of descending y
  coords.snap = coords.snap[order(coords.snap$i, coords.snap$j), ]
  coords.snap$i = max(coords.snap$i) - coords.snap$i + 1
  gy = rev(gy)

  # return list with all the info
  gxy=list(x=gx, y=gy)
  list(dims=sapply(gxy, length), gxy=gxy, snap=coords.snap)
}


#' Snap a set of points to a larger regular grid by least Manhattan distance
#'
#' @param coords n X 2 matrix of "x" and "y" coordinates or sf object
#' @param gxy list of two numeric vectors, the grid line positions "x" and "y"
#'
#' @return
#' @export
#'
#' @examples
pkern_snap_grid = function(coords, gxy)
{
  # The function returns a vector indexing grid points in grid.xy that best match the
  # point locations in gxy (shortest Manhattan distance). This is meant for snapping
  # data from point locations to a grid - eg if xy has data column 'foo', it can be
  # merged with the grid using...
  #
  # `grid.xy[ sk_snap_grid(xy, grid.xy), 'foo' ] = xy[,'foo']`
  #
  # ... or similar
  #
  # grid.xy should have coordinates columns 'x' and 'y', with coordinate pairs (rows)
  # ordered so that both 'x' and 'y' are increasing as you go down, and 'x' increases
  # the fastest. This is the ordering used by raster::coordinates and sf::st_coordinates.
  #
  # xy must also have 'x' and 'y' columns but its rows can be in any order.
  #
  # grid.xy must be a complete regular grid (ie no missing rows), as the function uses
  # the expected ordering for efficiency rather than matching numerical values of
  # 'x' and 'y'

  # convert RasterLayer and sf types for gxy to list of grid line positions
  if( 'sf' %in% class(gxy) ) gxy = sf::st_coordinates(gxy) |> apply(2, unique)
  if( 'RasterLayer' %in% class(gxy) ) gxy = raster::coordinates(gxy) |> apply(2, unique)

  # convert RasterLayer and sf types for coords to matrix of x, y values
  if( 'sf' %in% class(coords) ) coords = sf::st_coordinates(coords)
  if( 'RasterLayer' %in% class(coords) ) coords = raster::coordinates(coords)
  colnames(coords) = colnames(coords) |> tolower()

  # sort the grid line vectors (y descending, x ascending) and find dimensions
  gxy = list(x=sort(gxy[['x']], decreasing=FALSE), y=sort(gxy[['y']], decreasing=TRUE))
  dims = sapply(gxy, length)

  # define a grid structure for the input coordinates
  coords.grid = pkern_detect_grid(coords)
  coords.gxy = coords.grid$gxy
  coords.snap = coords.grid$snap

  # compute distance matrices and find least-distance matches to gxy grid lines
  idx.x = Rfast::Outer(coords.gxy[['x']], gxy[['x']], '-') |> abs() |> apply(2, which.min)
  idx.y = Rfast::Outer(coords.gxy[['y']], gxy[['y']], '-') |> abs() |> apply(2, which.min)
  gji = list(x=idx.x, y= idx.y)

  # update snapped grid line locations
  coords.gxy[['x']] = gxy[['x']][ idx.x ]
  coords.gxy[['y']] = gxy[['y']][ idx.y ]

  # define new snapped grid line locations
  xnew = coords.gxy[['x']][ coords.snap$j ]
  ynew = coords.gxy[['y']][ coords.snap$i ]

  # update output dataframe
  coords.snap$xnew = xnew
  coords.snap$ynew = ynew

  # define new distances (errors)
  coords.snap$dist = sqrt( (xnew - coords.snap$x)^2 + (ynew - coords.snap$y)^2 )

  # return list of outputs
  return( list(dims=coords.grid$dims, dims.outer=dims, gxy=coords.gxy, gji=gji, snap=coords.snap) )
}




