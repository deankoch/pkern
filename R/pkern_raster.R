#
# pkern_raster.R
# Dean Koch, Sep 2021
# Functions for converting gridded data to and from RasterLayer format
#

#' Check if terra or raster package is loaded and print error message otherwise
#'
#' @return logical vector, indicating if the package is loaded
#' @export
#'
#' @examples
#' pkern_checkRaster()
pkern_checkRaster = function()
{
  raster.loaded = requireNamespace('raster', quietly=TRUE)
  terra.loaded = requireNamespace('terra', quietly=TRUE)
  if( !(raster.loaded | terra.loaded) ) stop('terra/raster not loaded. Try library(terra)')
  return(invisible(c(raster=raster.loaded, terra=terra.loaded)) )
}


#' Load/vectorize a RasterLayer or SpatRaster (in column-vectorized order)
#'
#' The `values` function from the raster and terra packages returns data in row-vectorized
#' order. This function reorders them to column-vectorized order so they can be passed to
#' other pkern_* functions.
#'
#' In addition to the raster data, the function retrieves information on the spatial
#' configuration of the raster, in the y-x order expected by pkern_* functions
#'
#' The return value depends on argument `what`:
#'
#' "gdim": returns integer vector `c(ny, nx)`, the number of y and x grid lines
#' "gres": returns `c(dy, dx)`, the x and y distances between grid lines (ie the resolution)
#' "crs": returns the WKT string data on projections
#' "gyx": returns `list(y, x)`, a list of vectors containing the y and x grid line positions
#' "gval": returns a vector containing the raster data in column-vectorized order
#' "all": returns a named list containing the above five objects
#' "g": returns a named list containing everything but 'gval'
#'
#' @param r a RasterLayer, RasterStack, or SpatRaster to vectorize
#' @param what character, one of: 'gdim', 'gres', 'crs', 'gyx', 'gval', 'g', 'all'
#'
#' @return For default `what=='all'`, a named list containing 5 objects (see details)
#' @export
#'
#' @examples
#' if( requireNamespace('raster') ) {
#'
#' # open example file as RasterLayer, , then to SpatRaster
#' r.in = system.file('external/rlogo.grd', package='raster') |> raster::raster(band=1)
#'
#' # convert to list understood by pkern
#' pkern.in = pkern_fromRaster(r.in)
#' pkern_plot(pkern.in)
#' pkern.in$gval |> matrix(pkern.in$gdim) |> pkern_plot()
#'
#' # open a RasterStack and notice only first band loaded
#' r.in = system.file('external/rlogo.grd', package='raster') |> raster::stack()
#' str(pkern_fromRaster(r.in))
#'
#' # same with terra
#' if( requireNamespace('terra') ) {
#' str(pkern_fromRaster(terra::rast(r.in)))
#' }
pkern_fromRaster = function(r, what='all', quiet=FALSE)
{
  # stop with an error message if terra/raster packages unavailable
  pkern_checkRaster()

  # expected object classes and names
  raster_nm = c('RasterLayer', 'RasterStack')
  terra_nm = c('SpatRaster')

  # handle rasters
  if( any(raster_nm %in% class(r)) )
  {
    # extract dimensions and return if requested
    gdim = dim(r)[1:2] |> stats::setNames(c('y', 'x'))
    if( what == 'gdim' ) return(gdim)

    # notice return from raster::res is in reverse order (dx, dy)
    gres = raster::res(r)[2:1] |> stats::setNames(c('y', 'x'))
    if( what == 'gres' ) return(gres)

    # extract CRS string if available
    gcrs = raster::wkt(r)
    if( what == 'crs' ) return(gcrs)

    # extract grid line positions
    if( what != 'gval' )
    {
      gy = sort(raster::yFromRow(r, seq(gdim['y'])))
      gx = sort(raster::xFromCol(r, seq(gdim['x'])))
      gyx = list(y=gy, x=gx)

      # finish grid line mode
      if( what == 'gyx' ) return(gyx)
    }

    # extract data in column-vectorized order and finish
    gval = raster::getValues(r)[ pkern_r2c(gdim) ]
    if(what == 'gval') return(gval)
    return( list(gdim=gdim, gres=gres, crs=gcrs, gyx=gyx, gval=gval) )
  }

  # handle terra objects
  if( any(terra_nm %in% class(r)) )
  {
    # extract dimensions and return if requested
    gdim = dim(r)[1:2] |> stats::setNames(c('y', 'x'))
    if( what == 'gdim' ) return(gdim)

    # notice return from terra::res is in reverse order (dx, dy)
    gres = terra::res(r)[2:1] |> stats::setNames(c('y', 'x'))
    if( what == 'gres' ) return(gres)

    # extract CRS string if available
    gcrs = terra::crs(r)
    if( what == 'crs' ) return(gcrs)

    # extract grid line positions only when needed
    if( what != 'gval' )
    {
      gy = sort(terra::yFromRow(r, seq(gdim['y'])))
      gx = sort(terra::xFromCol(r, seq(gdim['x'])))
      gyx = list(y=gy, x=gx)

      # finish grid line mode
      if( what == 'gyx' ) return(gyx)
    }

    # return from no-values requests
    if( what == 'g' ) return( list(gdim=gdim, gres=gres, crs=gcrs, gyx=gyx) )

    # extract data in column-vectorized order and finish
    gval = terra::values(r)[ pkern_r2c(gdim) ]
    if(what == 'gval') return(gval)
    return( list(gdim=gdim, gres=gres, crs=gcrs, gyx=gyx, gval=gval) )

  }

  # unrecognized objects are returned unchanged with a warning
  if( !quiet ) warning('input r was not a raster or terra class object')
  return(r)
}

#' Convert column-vectorized grid to SpatRaster
#'
#' @param g either the column-vectorized data or a list containing it (in element "gval")
#' @param gdim c(ny, nx), the number rows and columns in the full grid
#' @param template optional RasterLayer to use as template (for setting crs, resolution etc)
#'
#' Converts a column-vectorized vector or matrix to a SpatRaster, or, if terra is unavailable,
#' a RasterLayer.
#'
#' `g` can either be a matrix (in which case `gdim` can be omitted) or its vectorized
#' equivalent, or a list containing elements returned by `pkern_fromRaster`. If the list
#' has an element named "gdim", it replaces any argument supplied to `gdim`. Other elements
#' of the list (eg 'crs') are passed on to the RasterLayer construction call as needed.
#'
#' When template is supplied, it supercedes any properties (eg 'crs') supplied in `rvec`.
#'
#' @return a RasterLayer containing the data from `rvec`
#' @export
#'
#' @examples
#' if( requireNamespace('raster') ) {
#' r.in = system.file('external/rlogo.grd', package='raster') |> raster::raster(band=1)
#' pkern.in = pkern_fromRaster(r.in)
#' r = pkern_toRaster(pkern.in)
#' print(r)
#' }
pkern_toRaster = function(g=NULL, gdim=NULL, template=NULL)
{
  # stop with an error message if raster/terra package is unavailable
  is_loaded = pkern_checkRaster()
  r_out = NULL

  # flags for different types of raster calls
  has.template = !is.null(template)
  has.properties = is.list(g)

  # assign dimensions from matrix input
  if( is.null(gdim) & is.matrix(g) ) gdim = rev(dim(g))

  # unpack list input
  gres = gyx = NULL
  rcrs = ''
  if( has.properties )
  {
    # unpack the data (overwrites argument to `gdim`)
    if( !is.null( g[['gdim']] ) ) gdim = g[['gdim']]
    if( !is.null( g[['crs']] ) ) rcrs = g[['crs']]
    gres = g[['gres']]
    gyx = g[['gyx']]

    # g should be a matrix or vector from here on
    g = g[['gval']]
  }

  # extract `gdim` from template as needed
  if( is.matrix(g) ) gdim = dim(g)
  if( is.null(gdim) )
  {
    # make sure we have a template then extract its dimensions in correct order
    if( !has.template ) stop('could not determine grid size. Try setting gdim')
    gdim = dim(template)[1:2]
  }

  # matricize data as needed
  has.values = !is.null(g)
  if( (!is.matrix(g)) & has.values ) g = matrix(g, gdim)

  # build from template if available
  if( has.template & has.values )
  {
    if( is_loaded['terra'] ) return( terra::setValues(terra::rast(template), values=g) )
    if( is_loaded['raster'] ) return( raster::raster(g, template=raster::raster(template)) )
    pkern_checkRaster()
  }

  # compute dimensional stuff when properties supplied
  if( has.properties )
  {
    # set defaults for gres and gyx
    if( is.null(gres) ) gres = c(1,1)
    if( is.null(gyx) ) gyx = lapply(setNames(gdim, c('y', 'x')), \(d) c(1, d))

    # extract coordinate "boundaries" as required by raster/terra
    yxb = Map(\(g,s) range(g) + (c(-1,1) * s/2), g=gyx, s=gres)

    # terra is preferred when available
    if( is_loaded['terra'] )
    {
      gext = terra::ext(do.call(c, rev(yxb)))
      r_out = terra::rast(extent=gext, resolution=rev(gres), crs=rcrs)
      if( has.values ) r_out = terra::setValues(r_out, g)

    } else if( is_loaded['raster'] ) {

    # attempt to use raster if terra unavailable
    gext = raster::extent(do.call(c, rev(yxb)))
    r_out = raster::raster(ext=gext, resolution=rev(gres), crs=rcrs)
    if( has.values ) r_out = raster::setValues(r_out, g)

    }
    return(r_out)
  }
}

# TODO: documentation for these

# unpack a points dataset into a list understood by other pkern functions
pkern_fromPoints = function(pts, what='all')
{
  # expected object classes and names
  sf_classes = c('sf','sfc', 'sfg')
  sp_classes = c('SpatialPoints', 'SpatialPointsDataFrame')
  yx_names = stats::setNames(nm=c('y', 'x'))

  # initialize outputs
  crs_out = NULL
  yx_out = NULL
  ext_out = NULL
  bbox_out = NULL
  is_spatial = FALSE

  # handle sf class objects
  if( any( sf_classes %in% class(pts) ) )
  {
    # copy crs wkt and convert matrix to list
    is_spatial = TRUE
    crs_out = sf::st_crs(pts)
    yx_out = sf::st_coordinates(pts)[,2:1] |> apply(2, identity, simplify=FALSE) |> setNames(yx_names)
  }

  # sp class objects are handled the same way, just with different function names
  if( any( sp_classes %in% class(pts) ) & requireNamespace('sp', quietly=TRUE) )
  {
    # copy crs wkt and convert matrix to list
    is_spatial = TRUE
    crs_out = sp::wkt(pts)
    yx_out = sp::coordinates(pts)[,2:1] |> apply(2, identity, simplify=FALSE) |> setNames(yx_names)
  }

  # convert matrices and dataframes to list (via matrix)
  if( is.data.frame(pts) & !is_spatial ) pts = as.matrix(pts)
  if( is.matrix(pts) )
  {
    yx_out = apply(pts, 2, identity, simplify=FALSE)
    if( !all(yx_names %in% names(yx_out) ) ) names(yx_out)[1:2] = yx_names
    yx_out = yx_out[1:2]
  }

  # read any user-supplied list-style inputs
  extra_out = list()
  if( is.list(pts) )
  {
    if( is.null(crs_out) ) crs_out = pts[['crs']]
    if( is.null(yx_out) ) yx_out = pts[['yx']]

    # copy any unrecognized list entries (to append to output)
    idx_overwrite = names(pts) %in% c('crs', 'yx', 'ext', 'bbox')
    if( sum(!idx_overwrite) > 0 ) extra_out = pts[!idx_overwrite]
  }

  # return coordinates or crs as requested
  if( what == 'yx') return(yx_out)
  if( what == 'crs') return(crs_out)

  # skip if coordinates not supplied
  if( !is.null(yx_out) )
  {
    # compute extent
    ext_out = lapply(yx_out, range)
    if( what == 'ext') return(ext_out)
    bbox_names = c('ymin', 'ymax', 'xmin', 'xmax')

    # build sf class bounding box polygon
    bbox_out = unlist(ext_out) |> setNames(bbox_names) |> sf::st_bbox(crs=crs_out) |> sf::st_as_sfc()
    if( what == 'bbox') return(bbox_out)
  }

  # warn about missing CRS as needed and finish
  if( is_spatial & is.null(crs_out) ) warning('no coordinate reference system found')
  return( c(list(yx=yx_out, crs=crs_out, ext=ext_out, bbox=bbox_out), extra_out) )
}

# defines a grid based number of grid lines and (optionally) an existing spatial extent
pkern_grid = function(gdim, gext=NULL)
{
  # defaults
  yx_names = stats::setNames(nm=c('y', 'x'))
  crs_ext = NULL
  gdim_ext = NULL

  # handle raster input to gext
  raster_classes = c('RasterLayer', 'RasterStack', 'SpatRaster')
  if( any(c(raster_classes) %in% class(gext)) )
  {
    crs_ext = pkern_fromRaster(gext, 'crs')
    gdim_ext = pkern_fromRaster(gext, 'gdim')
    gext = pkern_fromRaster(gext, 'gyx') |> lapply(range)
  }

  # handle list input to gext
  #if( is.list(gext) ) gext = as.data.frame(gext)

  # handle dataframe, matrix, sf, sp input to gext
  gext = pkern_fromPoints(gext)
  if( is.list(gext) )
  {
    crs_ext = gext[['crs']]
    gext = gext[['ext']]
  }

  # handle length-1 gdim argument
  if( length(gdim) == 1 ) gdim = rep(gdim, 2) |> stats::setNames(yx_names)

  # set default extent when no extent supplied (coordinates set to equal to grid line numbers)
  if( is.null(gext) ) gext = lapply(gdim, \(x) c(1, x) ) |> stats::setNames(yx_names)

  # compute new resolution and calculate equally spaced grid line positions
  gres = mapply(\(ext, d) diff(ext) / (d-1), ext=gext, d=gdim) |> stats::setNames(yx_names)
  gyx = Map(\(r, g) seq(r[1], r[2], by=diff(r)/g), r=gext, g=gdim-1) |> stats::setNames(yx_names)
  return( list(gyx=gyx, crs=crs_ext, gres=gres, gdim=gdim) )
}


