#
# pkern_raster.R
# Dean Koch, Sep 2021
# Functions for converting gridded data to and from RasterLayer format
#

#' Check if raster package is loaded and print error message otherwise
#'
#' @return nothing
#' @export
#'
#' @examples
#' pkern_checkRaster()
pkern_checkRaster = function()
{
  raster.loaded = requireNamespace('raster', quietly=TRUE)
  if( !raster.loaded ) stop('could not find raster package! Try library(raster)')
}


#' Vectorize a RasterLayer in column-vectorized order
#'
#' Calls to raster::values or raster::getValues produce the data in row-vectorized
#' order. This function reorders them to column-vectorized order so they can be
#' passed to other pkern_* functions
#'
#' The return value depends on argument `what`:
#'
#' "gdim": returns integer vector `c(ny, nx)`, the number of y and x grid lines
#' "gres": returns `c(dy, dx)`, the x and y distances between grid lines (ie the resolution)
#' "gyx": returns `list(y, x)`, a list of vectors containing the y and x grid line positions
#' "gval": returns a vector containing the raster data in column-vectorized order
#' "all": returns a named list containing the above four objects
#'
#' @param r a RasterLayer to vectorize
#' @param what character, one of: 'gdim', 'gres', 'crs', 'gyx', 'gval' or 'all'
#'
#' @return For default `what=='all'`, a named list containing 4 objects (see details)
#' @export
#'
#' @examples
#' if( requireNamespace('raster') ) {
#' r.in = system.file('external/rlogo.grd', package='raster') |> raster::raster(band=1)
#' pkern.in = pkern_fromRaster(r.in)
#' pkern_plot(pkern.in)
#' pkern.in$gval |> matrix(pkern.in$gdim) |> pkern_plot()
#' }
pkern_fromRaster = function(r, what='all')
{
  # stop with an error message if raster package is unavailable
  pkern_checkRaster()

  # extract dimensions and return if requested
  gdim = dim(r)[1:2] |> stats::setNames(c('y', 'x'))
  if( what == 'gdim' ) return(gdim)

  # notice return from raster::res is in reverse order (dx, dy)
  gres = raster::res(r)[2:1] |> stats::setNames(c('y', 'x'))
  if( what == 'gres' ) return(gres)

  # extract CRS string if available
  gcrs = raster::crs(r)
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

#' Convert column-vectorized grid to RasterLayer
#'
#' @param g either the column-vectorized data or a list containing it (in element "gval")
#' @param gdim c(ny, nx), the number rows and columns in the full grid
#' @param template optional RasterLayer to use as template (for setting crs, resolution etc)
#'
#' Produces a RasterLayer from a vector or matrix. This essentially a wrapper for various
#' `raster::raster` calls where the input data is in column-vectorized form.
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
#' print(pkern_toRaster(pkern.in))
#' print(r.in)
#' }
pkern_toRaster = function(g=NULL, gdim=NULL, template=NULL)
{
  # stop with an error message if raster package is unavailable
  pkern_checkRaster()

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
    # check for expected elements
    # nm.expected = c('gdim', 'gres', 'gyx', 'gval')
    # err.nm = paste('g must contain named elements', paste(nm.expected, collapse=', '))
    # if( !all( nm.expected %in% names(g) )) stop(err.nm)

    # unpack the data (overwrites argument to `gdim`)
    if( !is.null( g[['gdim']] ) ) gdim = g[['gdim']]
    if( !is.null( g[['crs']] ) ) rcrs = g[['crs']]
    gres = g[['gres']]
    gyx = g[['gyx']]

    # g should be a matrix or vector from here on
    g = g[['gval']]
  }

  # extract `gdim` from template when supplied
  if( is.null(gdim) )
  {
    # make sure we have a template then extract its dimensions in correct order
    if( !has.template ) stop('could not determine grid size. Try setting gdim')
    gdim = dim(template)[1:2]
  }

  # matricize data as needed
  has.values = !is.null(g)
  if( (!is.matrix(g)) & has.values ) g = matrix(g, gdim)

  # the raster call is simple with a template
  if( (has.template | !has.properties) & has.values ) return(raster::raster(g, template=template))

  # compute dimensional stuff when properties supplied
  if( has.properties )
  {
    # set defaults for gres and gyx
    if( is.null(gres) ) gres = c(1,1)
    if( is.null(gyx) ) gyx = Map(\(d,s) c(1,d)-(s/2), d=gdim, s=gres)

    # extract coordinate "boundaries" as required by raster
    yxb = Map(\(g,s) range(g) + (c(-1,1) * s/2), g=gyx, s=gres)

    # without a template the call has separate arguments for the four extent points
    if( has.values )
    {
      # build the raster and return
      return( raster::raster(g, yxb[[2]][1], yxb[[2]][2], yxb[[1]][1], yxb[[1]][2], crs=rcrs) )
    }
  }

  # without values we can use a simpler looking call by passing an extent object
  ext = raster::extent(do.call(c, yxb))
  return( raster::raster(crs=rcrs, ext=ext, resolution=gres) )
}


