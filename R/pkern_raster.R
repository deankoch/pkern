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
#' "ds": returns `c(dy, dx)`, the x and y distances between grid lines (ie the resolution)
#' "crs": returns the CRS code for the raster (if any)
#' "yx": returns `list(y, x)`, a list of vectors containing the y and x grid line positions
#' "values": returns a vector containing the raster data in column-vectorized order
#' "all": returns a named list containing the above three objects
#'
#' @param r a RasterLayer to vectorize
#' @param what character, one of: 'gdim', 'ds', 'u', 'crs', 'yx', 'values' or 'all'
#'
#' @return For default `what=='all'`, a named list containing 6 objects (see details)
#' @export
#'
#' @examples
#' if( requireNamespace('raster') ) {
#' r.in = system.file('external/rlogo.grd', package='raster') |> raster::raster(band=1)
#' pkern.in = pkern_fromRaster(r.in)
#' pkern_plot(pkern.in)
#' pkern.in$values |> matrix(pkern.in$gdim) |> pkern_plot()
#' }
pkern_fromRaster = function(r, what='all')
{
  # stop with an error message if raster package is unavailable
  pkern_checkRaster()

  # extract dimensions and return if requested
  gdim = dim(r)[1:2] |> stats::setNames(c('i', 'j'))
  if( what == 'gdim' ) return(gdim)

  # notice return from raster::res is in reverse order (dx, dy)
  ds = raster::res(r)[2:1] |> stats::setNames(c('y', 'x'))
  if( what == 'ds' ) return(ds)

  # extract CRS string if available
  crs.r = ifelse(requireNamespace('rgdal', quietly=TRUE), rgdal::CRSargs(raster::crs(r)), NA)
  if( what == 'crs' ) return(crs.r)

  # attempt to parse resolution units from CRS string
  u = NA
  regex.units = '+units\\='
  if( grepl(regex.units, crs.r) )
  {
    crs.split = strsplit(crs.r, regex.units)[[1]]
    u = substr(crs.split[ length(crs.split) ], 1, 1)
  }

  # return units if requested
  if( what == 'u' ) return(u)

  # extract grid line positions (as needed)
  if( what != 'values' )
  {
    yasc = sort( raster::yFromRow(r, seq(gdim[1]) ) )
    xasc = sort( raster::xFromCol(r, seq(gdim[2]) ) )
    yx = list(y=yasc, x=xasc)

    # finish grid line mode
    if( what == 'yx' ) return(yx)
  }

  # indexing vector to get column-vectorized version of values
  cvec.idx = pkern_r2c(gdim)

  # extract vectorized data and finish
  rvec = raster::getValues(r)[cvec.idx]
  if(what == 'values') return(rvec)
  return( list(gdim=gdim, ds=ds, crs=crs.r, u=u, yx=yx, values=rvec) )
}

#' Convert column-vectorized grid to RasterLayer
#'
#' @param rvec either the column-vectorized data or a list containing it (in element "values")
#' @param gdim c(ni, nj), the number rows and columns in the full grid
#' @param template optional RasterLayer to use as template (for setting crs, resolution etc)
#'
#' Produces a RasterLayer from a vector or matrix. This essentially a wrapper for various
#' `raster::raster` calls where the input data is in column-vectorized form.
#'
#' `rvec` can either be a matrix (in which case `gdim` can be omitted) or its vectorized
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
pkern_toRaster = function(rvec=NULL, gdim=NULL, template=NULL)
{
  # stop with an error message if raster package is unavailable
  pkern_checkRaster()

  # flags for different types of raster calls
  has.template = !is.null(template)
  has.properties = is.list(rvec)

  # assign dimensions from matrix input
  if( is.null(gdim) & is.matrix(rvec) ) gdim = rev(dim(rvec))

  # unpack list input
  if( is.list(rvec) )
  {
    # unpack the data (overwrites argument to `gdim`)
    if( !is.null( rvec[['gdim']] ) ) gdim = rvec[['gdim']]
    rcrs = ifelse( is.null(rvec[['crs']]), '', rvec[['crs']])
    ds = rvec[['ds']]
    yx = rvec[['yx']]

    # rvec should be a matrix or vector from here on
    rvec = rvec[['values']]
  }

  # extract `gdim` from template when supplied
  if( is.null(gdim) )
  {
    # make sure we have a template then extract its dimensions in correct order
    if( !has.template ) stop('could not determine grid size. Try setting gdim')
    gdim = dim(template)[1:2]
  }

  # matricize data as needed
  has.values = !is.null(rvec)
  if( (!is.matrix(rvec)) & has.values ) rvec = matrix(rvec, gdim)

  # the raster call is simple with a template
  if( (has.template | !has.properties) & has.values ) return(raster::raster(rvec, template=template))

  # compute dimensional stuff when properties supplied
  if( has.properties )
  {
    # set defaults for ds and yx
    if( is.null(ds) ) ds = c(1,1)
    if( is.null(yx) ) yx = Map(\(d,s) c(1, d)-(s/2), d=gdim, s=ds)

    # extract coordinate "boundaries" as required by raster
    yxb = Map(\(g,s) range(s*g) + (c(-1, 1) * s/2), g=yx, s=ds)

    # without a template the call has separate arguments for the four extent points
    if( has.values )
    {
      # build the raster and return
      return( raster::raster(rvec, yxb[[2]][1], yxb[[2]][2], yxb[[1]][1], yxb[[1]][2], crs=rcrs) )
    }
  }

  # without values we can use a simpler looking call by passing an extent object
  ext = raster::extent(do.call(c, yxb))
  return( raster::raster(crs=rcrs, ext=ext, resolution=ds) )
}


