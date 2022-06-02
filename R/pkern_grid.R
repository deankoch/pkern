# pkern_grid.R
# Dean Koch, 2022
# Functions for building grids out of various inputs


#' Make a pkern grid list object
#'
#' Define a 2-dimensional spatial grid as a list of vectors.
#'
#' This function accepts (in its first argument) 'RasterLayer' and 'RasterStack' objects from
#' the `raster` package, 'SpatRaster' objects from `terra`, as well as any non-complex matrix,
#' or a list containing the vectorization of one, or a vector of grid dimensions.
#'
#' It returns a list with the following 3-5 elements:
#'
#' * gyx: `list(y, x)`, the coordinates of the y and x grid lines in vectors `y` and `x`
#' * gres: `c(y, x)`, the (numeric) y and x distances between grid lines
#' * gdim: `c(y, x)`, the (integer) number of y and x grid lines
#' * gval: vector, the data (if any) in column-major order with y descending, x ascending
#' * crs: character, the WKT string (if available) describing coordinate reference system
#'
#' where 'crs' is included only for geo-referenced inputs from `raster` and `terra`.
#'
#' The input `g` can itself be a list containing a subset of these elements (including at least
#' one of 'gdim' or 'gyx'), and the function will fill in missing entries with their defaults:
#' If 'gval' is missing, the function sets NAs in the data vector; If 'res' is missing, it is
#' computed from the first two grid lines in 'gyx'. If 'gyx' is missing, it is assigned the sequence
#' `1:n` (scaled by 'res', if available) for each `n` in 'gdim'; and if 'gdim' is missing, it
#' is set to equal the number of grid lines specified in (each vector of) 'gyx'.
#'
#' `g` can also be a vector of the form `c(y, x)` defining grid dimensions (as a shorthand for
#' the call `pkern_grid(g=list(gdim=gdim))`). Note that 1-dimensional grids are not supported,
#' ie. there must be at least 2 grid lines in both the x and y dimensions.
#'
#' @param g grid object such as a matrix or raster or vector of grid dimensions (see details)
#' @param vals logical indicating to include the data vector 'gval' in return list
#'
#' @return named list containing 'gyx', 'gres', 'gdim', and optionally 'gval' and 'crs'
#' @export
#'
#' @examples
#'
#' # simple grid construction from dimensions
#' gdim = c(12, 10)
#' g = pkern_grid(gdim)
#' str(g)
#' str(pkern_grid(gdim, vals=FALSE))
#'
#' # supply grid lines instead to get the same result
#' all.equal(g, pkern_grid(g=list(gyx=lapply(gdim, function(x) seq(x)-1L))) )
#'
#' # display coordinates and grid line indices
#' pkern_plot(g)
#' pkern_plot(g, ij=TRUE)
#'
#' # set a different resolution and notice the argument is ignored if conflicting gyx supplied
#' gres_new = c(3, 4)
#' pkern_plot(pkern_grid(g=list(gyx=lapply(gdim, seq), gres=gres_new)))
#' pkern_plot(pkern_grid(g=list(gdim=gdim, gres=gres_new)))
#'
#' # shorthand for square grids
#' all.equal(pkern_grid(g=2), pkern_grid(g=c(2,2)))
#'
#' # example with data
#' gdim = c(25, 25)
#' yx = as.list(expand.grid(lapply(gdim, seq)))
#' eg_vec = as.numeric( yx[[2]] %% yx[[1]] )
#' eg_mat = matrix(eg_vec, gdim)
#' g = pkern_grid(eg_mat)
#' pkern_plot(g, ij=T, zlab='j mod i')
#'
#' # y is in descending order
#' pkern_plot(g, xlab='x = j', ylab='y = 26 - i', zlab='j mod i')
#'
#' # data vectors should be in R's default matrix vectorization order
#' all.equal(eg_vec, as.vector(eg_mat))
#' all.equal(g, pkern_grid(list(gdim=gdim, gval=as.vector(eg_mat))))
#'
#' if( requireNamespace('raster') ) {
#'
#' # open example file as RasterLayer
#' r_path = system.file('external/rlogo.grd', package='raster')
#' r = raster::raster(r_path, band=1)
#'
#' # convert to pkern list
#' g = pkern_grid(r)
#' pkern_plot(g)
#' pkern_plot(g, ij=T)
#'
#' # open a RasterStack and notice only first band loaded
#' r_multi = raster::stack(r_path)
#' str(pkern_grid(r_multi))
#' str(g)
#'
#' # same with terra
#' if( requireNamespace('terra') ) {
#'
#' g = pkern_grid(terra::rast(r_path))
#' str(g)
#' pkern_plot(g, ij=T)
#'
#' }
#' }
pkern_grid = function(g, vals=TRUE)
{
  # names for dimensional components and entries of g
  nm_dim = c('y', 'x')
  nm_g = c('gyx', 'gres', 'gdim')

  # handle raster objects
  is_terra = any(c('SpatRaster') %in% class(g))
  is_raster = any(c('RasterLayer', 'RasterStack') %in% class(g))
  if( is_terra | is_raster )
  {
    # copy grid dimensions then do package-specific calls
    gdim = dim(g)[1:2]
    if(is_terra)
    {
      # terra class
      gcrs = terra::crs(g)
      gyx = list(y=terra::yFromRow(g, seq(gdim[1])), x=terra::xFromCol(g, seq(gdim[2])))
      gres = terra::res(g)[2:1] # order dy, dx
      if(vals) gval = terra::values(g)[ matrix(seq(prod(gdim)), gdim, byrow=TRUE) ]

    } else {

      # raster class
      gcrs = raster::wkt(g)
      gyx = list(y=raster::yFromRow(g, seq(gdim[1])), x=raster::xFromCol(g, seq(gdim[2])))
      gres = raster::res(g)[2:1] # order dy, dx
      if(vals) gval = raster::getValues(g)[ matrix(seq(prod(gdim)), gdim, byrow=TRUE) ]
    }

    # sort both sets of grid lines into ascending order
    gyx = lapply(gyx, sort)

    # build named list and return
    g_out = lapply(list(gyx, gres, gdim), function(r) stats::setNames(r, nm_dim))
    g_out = c(list(crs=gcrs), stats::setNames(g_out, nm_g))
    if( !vals ) return(g_out)
    return( c(list(gval=as.vector(gval)), g_out) )
  }

  # handle matrix objects
  if( is.matrix(g) )
  {
    # set unit resolution by default
    g_out = list(gres=c(y=1, x=1), gdim=dim(g))
    g_out[['gyx']] = lapply(g_out[['gdim']], function(d) as.numeric(seq(d)-1L))

    # build named list and return
    g_out = stats::setNames(lapply(g_out[nm_g], function(r) stats::setNames(r, nm_dim)), nm_g)
    if(!vals) return(g_out)
    return( c(list(gval=as.vector(g)), g_out) )
  }

  # handle list objects
  if( is.list(g) )
  {
    # shortcut to return the list unchanged when all required names are there
    if( all(nm_g %in% names(g)) )
    {
      # remove gval if not requested
      is_gval = names(g) == 'gval'
      if( any(is_gval) & !vals ) return(g[!is_gval])
      return(g)
    }

    # calculate both gdim and gres from gyx when gdim is missing
    if( !( 'gdim' %in% names(g) ) )
    {
      # require gyx in this case
      if( !( 'gyx' %in% names(g) ) ) stop('gdim and gyx both missing from input grid g')
      g[['gdim']] = sapply(g[['gyx']], length)
    }

    # set up default grid resolution
    g[['gdim']] = as.integer(g[['gdim']])
    if( is.null(g[['gres']]) )
    {
      g[['gres']] = c(1, 1)
      if( !is.null(g[['gyx']]) ) g[['gres']] = as.numeric(sapply(g[['gyx']], function(r) diff(r)[1]))

    } else {

      # resolution if it was supplied as a scalar
      if( length(g[['gres']]) == 1 ) g[['gres']] = stats::setNames(rep(g[['gres']], 2), nm_dim)
    }

    # set up default grid line positions
    if( is.null(g[['gyx']]) ) g[['gyx']] = Map(function(d, r) as.numeric(r*(seq(d)-1L)),
                                               d=g[['gdim']],
                                               r=g[['gres']])

    # set names and add crs and grid values if they're needed before returning
    g_out = lapply(g[nm_g], function(r) stats::setNames(r, nm_dim))
    if( !is.null(g[['crs']]) ) g_out = c(list(crs=g[['crs']]), g_out)
    if(!vals) return(g_out)
    if( is.null(g[['gval']]) ) g[['gval']] = rep(NA_real_, prod(g[['gdim']]))
    return(c(list(gval=g[['gval']]), g_out))
  }

  # handle numeric vectors
  if( is.vector(g) )
  {
    if( length(g) > 2 ) stop('numeric vector g must be of length 1 or 2')
    if( length(g) == 1 ) g = stats::setNames(rep(g, 2), nm_dim)
    return( pkern_grid(g=list(gdim=as.integer(g)), vals) )
  }

  # unrecognized objects are returned unchanged with a warning
  warning('input g was not recognized')
  return(NA)
}


#' Convert column-vectorized grid to SpatRaster
#'
#' @param g any object accepted or returned by `pkern_grid`
#' @param template character or RasterLayer/SpatRaster to set output type
#'
#' Converts a column-vectorized vector or matrix to a SpatRaster, or if terra is
#' unavailable, a RasterLayer.
#'
#' @return a RasterLayer or SpatRaster containing the data from `g` (or a sub-grid)
#' @export
#'
#' @examples
#'
#' if( requireNamespace('raster') ) {
#'
#' # open example file as RasterLayer
#' r_path = system.file('external/rlogo.grd', package='raster')
#' r = raster::raster(r_path, band=1)
#' g = pkern_grid(r)
#'
#' # convert back to RasterLayer and compare
#' r_from_g = pkern_export(g, 'raster')
#' print(r_from_g)
#' print(r)
#'
#' # layer name, band number, and various other metadata are lost
#' all.equal(r_from_g, r)
#'
#' # same with terra
#' if( requireNamespace('terra') ) {
#'
#' # SpatRaster is the default return class (when terra loaded)
#' r = terra::rast(r_path)
#' g = pkern_grid(r)
#' r_from_g = pkern_export(g)
#'
#' # notice only the first band is loaded by pkern_grid
#' print(r_from_g)
#' print(r)
#' }
#' }
#'
pkern_export = function(g, template='terra')
{
  # stop with an error message if raster/terra package is unavailable
  pkg_check = stats::setNames(nm=c('raster', 'terra'))
  pkg_msg = paste(pkg_check, collapse=' or ')
  msg_dep = paste(pkg_msg, 'package must be loaded first. Try `install.packages(terra)`')
  is_loaded = sapply(pkg_check, function(pkg) requireNamespace(pkg, quietly=TRUE))
  if( !any(is_loaded) ) stop(msg_dep)

  # load the input as pkern list
  g = pkern_grid(g)
  g[['crs']] = ifelse(is.null(g[['crs']]), '', g[['crs']])

  # extract grid cell boundaries as defined in raster/terra
  yx_bbox = Map(\(g, s) range(g) + (c(-1,1) * s/2), g=g[['gyx']], s=g[['gres']])

  # handle grid objects as templates
  if( !is.character(template) )
  {
    # check if template is a sub-grid or super-grid of `g`?
    # g_template = pkern_grid(template)
    # TODO: implement crop?

    # set template class name
    template = class(g)
    if( 'SpatRaster' %in% template ) template = 'terra'
    if( any(c('RasterLayer', 'RasterStack') %in% template ) ) template = 'raster'
    template = paste(template, collapse=', ')
  }

  # terra is preferred when available
  if( template == 'terra' )
  {
    g_ext = terra::ext(do.call(c, rev(yx_bbox)))
    r_out = terra::rast(extent=g_ext, resolution=rev(g[['gres']]), crs=g[['crs']])
    if( !is.null(g[['gval']]) ) r_out = terra::setValues(r_out, matrix(g[['gval']], g[['gdim']]))
    return(r_out)
  }

  # check for unknown class
  if( template == 'raster' )
  {
    # attempt to use raster if terra unavailable
    g_ext = raster::extent(do.call(c, rev(yx_bbox)))
    r_out = raster::raster(ext=g_ext, resolution=rev(g[['gres']]), crs=g[['crs']])
    if( !is.null(g[['gval']]) ) r_out = raster::setValues(r_out, matrix(g[['gval']], g[['gdim']]))
    return(r_out)
  }

  # error if we didn't get an expected class
  stop(paste('unrecognized template class', template))
}


#' Snap a set of points to a grid
#'
#' Maps the input points in `from` to the closest grid points in the extension of `g`
#' covering the bounding box of `from` (ie. the lattice of which `g` is a sub-grid).
#' In cases of duplicate mappings, the function returns the first matches only.
#'
#' `from` can be a geometry collection from packages `sf` or `sp`, or a matrix or list
#' of y and x coordinates (and, optionally, data values to copy to the snapped grid
#' points). When `from` is a matrix, its first two columns should by y and x, and the
#' (optional) third the data. When `from` is a list, the function expects (two or three)
#' vectors of equal length, ordered as above.
#'
#' When `from` is a geometry collection with a CRS string, points are first transformed
#' to the coordinate reference system of `g`. If one or both of `from` and `g` are missing
#' a CRS definition, the function assumes the same one is shared in both.
#'
#' `g` can be a raster geometry object (such as SpatRaster), in which case the function
#' behaves like `terra::rasterize`. It can also be a matrix (supplying dimensions) or a
#' list containing either `gdim` or`gres`, from which an appropriately spaced set of grid
#' lines is derived, centered under the bounding box of the points.
#'
#' `trim` controls the extent of the output grid. If `trim=FALSE` the function returns
#' the smallest regular grid containing both `g` (a sub-grid) and the snapped `from`
#' points. If `trim=TRUE` NA border columns and rows are omitted, so that the grid is
#' flush over the snapped points.
#'
#' @param from matrix, data frame, or points object from `sp` or `sf`, the source points
#' @param g any grid object accepted or returned by `pkern_grid`, the destination grid
#' @param trim logical, indicating to trim NA borders from the result.
#'
#' @return list of form returned by `pkern_grid`, defining a grid containing the snapped
#' points. These are assigned the corresponding data value in `from`, or if  `from` has no
#' data, an integer mapping to the points in `from`. Un-mapped grid points are set to NA.
#' @export
#'
#' @examples
#'
#' # functions to scale arbitrary inverval to (1, 2,... 100) and make color palettes
#' num_to_cent = function(x) 1L + floor(99*( x-min(x) ) / diff(range(x)))
#' my_pal = function(x) hcl.colors(x, 'Spectral', rev=T)
#' my_col = function(x) my_pal(1e2)[ num_to_cent(x) ]
#'
#' # create a grid object
#' gdim = c(40, 30)
#' g = pkern_grid(list(gdim=gdim, gres=1.1))
#'
#' # randomly position points within bounding box of g
#' n_pts = 10
#' from = lapply(g$gyx, function(yx) runif(n_pts, min(yx), max(yx)) )
#'
#' # translate away from g (no overlap is required)
#' from[['y']] = from[['y']] + 5
#' from[['x']] = from[['x']] + 15
#'
#' # add example data values and plot
#' from[['z']] = rnorm(length(from[['y']]))
#' pkern_plot(g, reset=FALSE)
#' graphics::points(from[c('x', 'y')], pch=16, col=my_col((from[['z']])))
#' graphics::points(from[c('x', 'y')])
#'
#' # snap points and plot
#' g_snap = pkern_snap(from, g)
#' pkern_plot(g_snap, col_grid='black', reset=FALSE)
#' graphics::points(from[c('x', 'y')], pch=16, col=my_col((from[['z']])))
#' graphics::points(from[c('x', 'y')])
#'
#' # find smallest subgrid enclosing all snapped grid points
#' g_snap = pkern_snap(from, g, trim=TRUE)
#' pkern_plot(g_snap, col_grid='black', reset=FALSE)
#' graphics::points(from[c('x', 'y')], pch=16, col=my_col((from[['z']])))
#' graphics::points(from[c('x', 'y')])
#'
#' # create a new grid of different resolution enclosing all input points
#' g_snap = pkern_snap(from, g=list(gres=c(0.5, 0.5)))
#' pkern_plot(g_snap, reset=FALSE)
#' graphics::points(from[c('x', 'y')], pch=16, col=my_col((from[['z']])))
#' graphics::points(from[c('x', 'y')])
#'
#' if( requireNamespace('sf') ) {
#'
#' # a different example, snapping misaligned subgrid
#' g_pts = pkern_grid(list(gdim=c(15, 8), gres=1.7), vals=FALSE)
#' g_pts[['gyx']][['y']] = g_pts[['gyx']][['y']] + 5
#' g_pts[['gyx']][['x']] = g_pts[['gyx']][['x']] + 5
#' n_pts = prod(g_pts$gdim)
#' from = pkern_coords(g_pts, out='list')
#'
#' # convert to sf
#' eg_sfc = sf::st_geometry(pkern_coords(g_pts, out='sf'))
#' pkern_plot(g, reset=FALSE)
#' plot(eg_sfc, add=TRUE)
#'
#' # generate example data and plot
#' eg_sf = sf::st_sf(data.frame(z=rnorm(n_pts)), geometry=eg_sfc)
#' pkern_plot(g, reset=FALSE)
#' plot(eg_sf, pch=16, add=TRUE, pal=my_pal)
#' plot(eg_sfc, add=TRUE)
#'
#' # snap points
#' g_snap = pkern_snap(from=eg_sf, g)
#' pkern_plot(g_snap, reset=FALSE, col_grid='black')
#' plot(eg_sf, pch=16, add=TRUE, pal=my_pal)
#' plot(eg_sfc, add=TRUE)
#'
#' # snapping points without data produces the mapping (non-NA values index "from")
#' g_snap = pkern_snap(from=eg_sfc, g)
#' pkern_plot(g_snap, ij=TRUE, reset=FALSE, col_grid='black')
#' plot(eg_sfc, add=TRUE)
#'
#' # with trim=TRUE
#' g_snap = pkern_snap(from=eg_sfc, g, trim=TRUE)
#' pkern_plot(g_snap, reset=FALSE, col_grid='black')
#' plot(eg_sfc, add=TRUE)
#'
#' # test with sp class
#' eg_sp = as(eg_sf,'Spatial')
#' g_snap = pkern_snap(from=eg_sp, g)
#' pkern_plot(g_snap, reset=FALSE, col_grid='black')
#' plot(eg_sf, pch=16, add=TRUE, pal=my_pal)
#' plot(eg_sfc, add=TRUE)
#'
#' }
#'
pkern_snap = function(from, g=NULL, trim=FALSE)
{
  # set up more appropriate grid lines later input g doesn't specify them
  gres_input = NULL
  if( is.matrix(g) ) g = list(gdim=dim(g))
  auto_gyx = (!is.list(g) & is.vector(g)) | is.list(g) & !('gyx' %in% names(g))

  # set up grid dimensions later when only gres specified (set dummy dimensions now)
  auto_gdim = ifelse(auto_gyx & is.list(g), ('gres' %in% names(g)) & !('gdim' %in% names(g)), FALSE)
  if(auto_gdim) g = modifyList(g, list(gdim=c(y=2L, x=2L)))

  # unpack second argument as pkern grid list but don't copy data values
  g = pkern_grid(g, vals=FALSE)
  to_crs = g[['crs']]
  vals = NULL

  # expected point object classes and names
  nm_yx = c('y', 'x')
  is_sf = any( c('sf','sfc', 'sfg') %in% class(from) )
  is_sp = any( c('SpatialPoints', 'SpatialPointsDataFrame') %in% class(from) )
  to_crs = NULL

  # get coordinate(s) from sf objects as matrix/vector
  if(is_sf)
  {
    # transform to destination coordinate system, if supplied (otherwise assume same as source)
    from_crs = sf::st_crs(from)[['wkt']]
    is_geo = !is.null(to_crs) & !is.na(from_crs)
    if( is_geo ) { from = sf::st_transform(from, to_crs) } else { to_crs = from_crs }

    # copy values if possible
    if( length(names(from)) > 1 ) vals = unlist(sf::st_drop_geometry(from)[,1])

    # reverse column order to get y, x coordinates as matrix
    from = apply(sf::st_coordinates(from)[,2:1], 2L, identity)
  }

  # the same for sp objects
  if(is_sp)
  {
    from_crs = sp::wkt(from)
    is_geo = !is.null(to_crs) & !is.null(from_crs)
    if( is_geo ) { from = sp::spTransform(from, to_crs) } else { to_crs = from_crs }

    # grab first column if there are any
    if( length(names(from)) > 0 ) vals = unlist(from[[names(from)]])
    from = apply(sp::coordinates(from)[,2:1], 2L, identity)
  }

  # convert (length-2) vector, matrix, and dataframe to list
  if( is.data.frame(from) ) from = as.matrix(from)
  if( is.list(from) ) from = do.call(cbind, from)
  if( is.vector(from) ) from = as.matrix(from, nrow=1)
  if( is.matrix(from) ) from = apply(from, 2, identity, simplify=FALSE)

  # copy data and find bounding box
  if( length(from) > 2 ) vals = from[[3]]
  from_yx = stats::setNames(from[1:2], nm_yx)
  from_n = length(from_yx[['x']])
  from_bds = lapply(list(min, max), function(f) sapply(from_yx, f))

  # automatically set grid lines
  if(auto_gyx)
  {
    # automatically set grid dimensions based on requested resolution
    if(auto_gdim)
    {
      # given gres, compute number of grid lines neeeded to span the bounding box
      g[['gdim']] = 1L + sapply(from_yx, function(yx) diff(range(yx))) %/% g[['gres']]

      # find offsets that center the grid lines under point bounding box
      auto_pad = ( sapply(from_yx, function(yx) diff(range(yx))) %% g[['gres']] ) / 2
      gyx = Map(function(yx, p, r, n) seq(yx + p, by=r, length.out=n),
                yx=from_bds[[1]], p=auto_pad, r=g[['gres']], n=g[['gdim']])

      # make grid object
      g = pkern_grid(list(gdim=g[['gdim']], gyx=gyx, gres=g[['gres']]))

    } else {

    # place grid lines to coincide with bounding box edge points, then recompute gres
    gyx = Map(function(yx, n) seq(min(yx), max(yx), length.out=n), yx=from_yx, n=g[['gdim']])
    g = pkern_grid(list(gdim=g[['gdim']], gyx=gyx, gres=NULL))

    }
  }

  # find bounding box of destination grid template and reshape to list of min and max
  g_bbox = lapply(pkern_coords(g, out='list', corner=TRUE), range)
  g_bds = lapply(list(min, max), function(f) sapply(g_bbox, f))

  # find the offsets between the two bounding boxes
  to_pad = Map(function(a, b) abs((a - b) %% g[['gres']]), a=from_bds, b=g_bds)
  to_min = from_bds[[1]] - to_pad[[1]] + as.integer(to_pad[[1]] > (g[['gres']]/2)) * g[['gres']]
  to_max = from_bds[[2]] - to_pad[[2]] + as.integer(to_pad[[2]] > (g[['gres']]/2)) * g[['gres']]

  # if not cropping, extend these grid lines to include all of g_bbox
  if( !trim )
  {
    to_min = pmin(to_min, sapply(g_bbox, min))
    to_max = pmax(to_max, sapply(g_bbox, max))
  }

  # compute new grid line locations and initialize the output grid list object
  to_yx = Map(function(a, b, r) seq(a, b, by=r), a=to_min, b=to_max, r=g[['gres']])
  g_out = pkern_grid(list(gyx=to_yx), vals=FALSE)
  if( !is.null(to_crs) ) g_out[['crs']] = to_crs

  # find cross-distance matrices for point coordinates and grid lines
  d_yx_all = Map(function(a, b) abs( outer(a, b, '-'))^2, a=from_yx, b=to_yx)

  # find minimum distance mappings
  ij_min = stats::setNames(lapply(d_yx_all, function(d) max.col(-d, ties.method='f')), c('i', 'j'))
  ij_min[['i']] = g_out[['gdim']]['y'] + 1L - ij_min[['i']]
  to_idx = pkern_mat2vec(ij_min, g_out[['gdim']])


  # handle multiple points mapping to a single grid-point
  is_dupe = duplicated(to_idx)
  if( any(is_dupe) ) warning( paste('omitting', sum(is_dupe), 'duplicate mapping(s)') )

  # match magic to get NAs at unassigned grid points
  to_all = seq( prod(g_out[['gdim']]) )
  from_idx = match(to_all, to_idx)
  if( is.null(vals) ) { gval = from_idx } else { gval = vals[from_idx] }
  return( c(list(gval=gval), g_out) )
}


# convenience function for downscaling
# automates model-fitting and kriging
pkern_down = function(g, down, X=NA)
{



}


