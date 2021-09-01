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
#' order. This function reorders them to column-vectorized order, optionally returning
#' the dimensions in (ny, nx) format (the reverse of what you get with dim(raster))
#'
#' @param r a RasterLayer to vectorize
#' @param omode character indicating return value type, either 'list' or 'vector'
#'
#' @return either the data vector or a list containing it along with the dimensions
#' @export
#'
#' @examples
#' require(raster)
#' r.in = system.file('external/rlogo.grd', package='raster') |> raster(band=1)
#' pkern.in = pkern_fromRaster(r.in)
#' pkern.in$data |> matrix(ncol=pkern.in$dims[1]) |> image()
#' pkern_toRaster(pkern.in, template=r.in)
#' pkern_toRaster(pkern_fromRaster(r.in, omode='vector'), template=r.in)
pkern_fromRaster = function(r, omode='list')
{
  # stop with an error message if raster package is unavailable
  pkern_checkRaster()

  # extract dimensions and indexing vector to get column-vectorized version
  dims = dim(r)[2:1] |> stats::setNames(c('ny', 'nx'))
  cvec.idx = pkern_r2c(dims)

  # extract vectorized data in column-vectorized order
  rvec = raster::values(r)[cvec.idx]

  # output type depends on omode
  if( omode == 'list' ) { return( list(dims=dims, data=rvec) ) }
  if( omode == 'vector' ) { return( rvec ) }
}

#' Convert column-vectorized grid to RasterLayer
#'
#' @param rvec either the data vector or a list with elements "data" and "dims"
#' @param dims c(nx, ny), the number of x and y grid lines in the full grid
#' @param template optional RasterLayer to use as template (for setting crs, resolution etc)
#'
#' Produces a RasterLayer from a column-vectorized dataset. This is simply a
#' wrapper for the call `raster(matrix(rvec, ncol=dims[1]), template=template)`
#' with some error checking and the ability to pass both the data vector and its
#' dimensions together as a list (in "rvec").
#'
#' @return a RasterLayer containing the data from "rvec"
#' @export
#'
#' @examples
#' require(raster)
#' r.in = system.file('external/rlogo.grd', package='raster') |> raster(band=1)
#' pkern.in = pkern_fromRaster(r.in)
#' pkern_toRaster(pkern.in, template=r.in)
#' r.in
pkern_toRaster = function(rvec=NA, dims=NULL, template=NULL)
{

  # stop with an error message if raster package is unavailable
  pkern_checkRaster()

  # handle list input
  if( is.list(rvec) )
  {
    # expected list element names
    rvec.names = c('dims', 'data')

    # handle unnamed list input
    if( !all( rvec.names %in% names(rvec) ) )
    {
      # handle too many elements in unnamed list
      if( length(rvec) != 2 ) stop('list "rvec" should contain only two elements, "dims" and "data"')

      # identify the dimensions element as the shorter vector
      if( diff(sapply(rvec, length)) < 0 ) rvec.names = rev(rvec.names)
      names(rvec) = rvec.names
    }

    # warn the user if they supplied a list and an argument to dims
    if( !is.null(dims) ) warning('the first element of list "rvec" will supercede "dims"')

    # unpack the data
    dims = rvec$dims
    rvec = rvec$data
  }

  # handle NULL dims by extracting them from template
  if( is.null(dims) )
  {
    # make sure we have a template then extract its dimensions in correct order
    if( is.null(template) ) stop('one of "dims" or "template" must be supplied')
    dims = dim(template)[2:1]
  }

  # convert to raster, using template if provided
  return( raster::raster(matrix(rvec, ncol=dims[1]), template=template) )
}
