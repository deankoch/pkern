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
#' "dims": returns `c(nx, ny)`, the length-2 vector of dimensions (in reverse order)
#' "ds": returns `c(dx, dy)`, the x and y distances between grid lines (ie the resolution)
#' "crs": returns the CRS code for the raster (if any)
#' "u": returns a character string representing the units for "ds" (extracted from "crs")
#' "xy": returns `list(x, y)`, a list of vectors containing the x and y grid line positions
#' "values": returns a vector containing the raster data in column-vectorized order
#' "all": returns a named list containing the above three objects
#'
#' @param r a RasterLayer to vectorize
#' @param what character, one of: 'dims', 'ds', 'u', 'crs', 'xy', 'values' or 'all'
#'
#' @return For default `what=='all'`, a named list of length six (see details)
#' @export
#'
#' @examples
#' require(raster)
#' r.in = system.file('external/rlogo.grd', package='raster') |> raster::raster(band=1)
#' pkern.in = pkern_fromRaster(r.in)
#' pkern.in$values |> matrix(nrow=pkern.in$dims[2]) |> pkern_plot()
#' pkern_toRaster(pkern.in, template=r.in)
#' pkern_toRaster(pkern_fromRaster(r.in, what='values'), template=r.in)
pkern_fromRaster = function(r, what='all')
{
  # stop with an error message if raster package is unavailable
  pkern_checkRaster()

  # extract dimensions and return if requested
  dims = dim(r)[2:1] |> stats::setNames(c('x', 'y'))
  if( what == 'dims' ) return(dims)

  # extract resolution and return if requested
  ds = raster::res(r) |> stats::setNames(c('x', 'y'))
  if( what == 'ds' ) return(ds)

  # extract CRS string if available
  crs.r = rgdal::CRSargs( raster::crs(r) )
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
  if( what != 'values' ) gxy = list(x=raster::xFromCol(r, seq(dims[1])),
                                    y=raster::yFromRow(r, rev(seq(dims[2]))) )

  # finish grid line mode
  if( what == 'xy' ) return(gxy)

  # indexing vector to get column-vectorized version of values
  cvec.idx = pkern_r2c(dims)

  # extract vectorized data and finish
  rvec = raster::values(r)[cvec.idx]
  if(what == 'values') return(rvec)
  return( list(dims=dims, ds=ds, crs=crs.r, u=u, xy=gxy, values=rvec) )
}


#' Convert column-vectorized grid to RasterLayer
#'
#' @param rvec either the column-vectorized data or a list containing it (in element "values")
#' @param dims c(nx, ny), the number of x and y grid lines in the full grid
#' @param template optional RasterLayer to use as template (for setting crs, resolution etc)
#'
#' Produces a RasterLayer from a vector or matrix. This essentially a wrapper for various
#' `raster::raster` calls where the input data is in column-vectorized form.
#'
#' `rvec` can either be a matrix (in which case `dims` can be omitted) or its vectorized
#' equivalent, or a list containing elements returned by `pkern_fromRaster`. If the list
#' as an element named "dims", it replaces any argument supplied to `dims`. Other elements
#' of the list (eg 'crs') are passed on to the RasterLayer construction call as needed.
#'
#' When template is supplied, it supercedes any properties (eg 'crs') supplied in `rvec`.
#'
#' @return a RasterLayer containing the data from `rvec`
#' @export
#'
#' @examples
#' require(raster)
#' r.in = system.file('external/rlogo.grd', package='raster') |> raster(band=1)
#' pkern.in = pkern_fromRaster(r.in)
#' pkern_toRaster(pkern.in)
#' r.in
pkern_toRaster = function(rvec=NULL, dims=NULL, template=NULL)
{
  # stop with an error message if raster package is unavailable
  pkern_checkRaster()

  # flags for different types of raster calls
  has.template = !is.null(template)
  has.properties = is.list(rvec)

  # assign dimensions from matrix input
  if( is.null(dims) & is.matrix(rvec) ) dims = rev(dim(rvec))

  # unpack list input
  if( is.list(rvec) )
  {
    # unpack the data (overwrites argument to dims)
    if( !is.null( rvec[['dims']] ) ) dims = rvec[['dims']]
    rcrs = ifelse( is.null(rvec[['crs']]), '', rvec[['crs']])
    ds = rvec[['ds']]
    xy = rvec[['xy']]

    # rvec should be a matrix or vector from here on
    rvec = rvec[['values']]
  }

  # extract `dims` from template when supplied
  if( is.null(dims) )
  {
    # make sure we have a template then extract its dimensions in correct order
    if( !has.template ) stop('could not determine grid size. Try setting dims')
    dims = raster::dim(template)[2:1]
  }

  # matricize data as needed
  has.values = !is.null(rvec)
  if( (!is.matrix(rvec)) & has.values ) rvec = matrix(rvec, dims[2])

  # the raster call is simple with a template
  if( (has.template | !has.properties) & has.values ) return(raster::raster(zmat, template=template))

  # compute dimensional stuff when properties supplied
  if( has.properties )
  {
    # set defaults for ds and xy
    if( is.null(ds) ) ds = c(1,1)
    if( is.null(xy) ) xy = mapply(\(d,s) c(1, d)-(s/2), d=dims, s=ds, SIMPLIFY=FALSE)

    # extract coordinate "boundaries" as required by raster
    xyb = mapply(\(g,s) range(s*g) + (c(-1, 1) * s/2), g=xy, s=ds, SIMPLIFY=FALSE)

    # without a template the call has separate arguments for the four extent points
    if( has.values )
    {
      # build the raster and return
      return( raster::raster(zmat, xyb[[1]][1], xyb[[1]][2], xyb[[2]][1], xyb[[2]][2], crs=rcrs) )
    }
  }

  # without values we can use a simpler looking call by passing an extent object
  ext = raster::extent(do.call(c, xyb))
  return( raster::raster(crs=rcrs, ext=ext, resolution=ds) )
}


#' Snap points to a subgrid in a RasterLayer
#'
#' @param pts, numeric vector or n x 2 matrix (with columns for x and y) of point coordinates
#' @param g, RasterLayer
#' @param regular, logical, indicating to select a regular set of grid lines
#'
#' @return a RasterLayer containing the data from `pts`
#' @export
#'
#' @examples
pkern_rasterize = function(pts, g, regular)
{



}

