#
# pkern_raster.R
# Dean Koch, Sep 2021
# Functions for plotting objects created with pkern
#

#' Extract kernel parameters as plot-friendly strings
#'
#' Generates a list of useful strings based on the kernel names and parameter values
#' in `pars` for plot labels and such. `pars` should either specify kernel(s) by name
#' or else be a parameters list of the form accepted by `pkern_corr`.
#'
#' When `pars` is a character string, the function parenthesizes it with square brackets
#' and returns the resulting string. When `pars` is a vector of character strings, it
#' collapses the two parenthesized strings into a single vector separated by an "x".
#'
#' When `pars` is a list of kernel parameters for a single dimension, the function
#' returns a named list of two character strings: the kernel name ("k") formatted as
#' above, and the parameter(s) ("kp") (also square bracketed, in "name = value" format
#' where "value" is rounded to `nsig` significant digits).
#'
#' When `pars` is a list of two such kernel parameter lists (named "x" and "y"), the
#' function again returns a list with elements "k" and "kp", containing the strings from
#' each component collapsed with a " x " delimiter (indicating direct product). The third
#' element, "main", is a language object containing the bolded kernel name, and (if "nug"
#' is supplied) a fourth element supplies a string "nug" with the nugget effect in
#' parentheses.
#'
#' @param pars a kernel parameters list ("k", "kp"), or list of two of them ("x", "y")
#' @param nsig number of signficant figures to use when printing model parameter values
#'
#' @return a character string or list (depending on the contents of `pars`, see details)
#' @export
#'
#' @examples
#' kname = 'mat'
#' pkern_toString(kname)
#' pkern_toString(rep(kname, 2))
#' pkern_toString(list(k=kname))
#'
#' pars = pkern_corr(kname)
#' pkern_toString(pars)
#'
#' kname = c('mat', 'sph')
#' pars = pkern_corr(kname)
#' pkern_toString(pars)
#' pkern_toString(utils::modifyList(pars, list(nug=0.5)))
#'
pkern_toString = function(pars, nsig=3)
{
  # handle character input (kernel names)
  if( is.character(pars) )
  {
    # convert to expected list format for recursive call and handle unexpected input
    if( length(pars) == 1 ) return( paste0('[', pars, ']') )
    if( length(pars) == 2 ) return( paste(sapply(pars, pkern_toString), collapse=' x ') )
    stop('invalid argument to pars')
  }

  # single dimension case
  if( 'k' %in% names(pars) )
  {
    k = pkern_toString(pars[['k']])
    kp = '[]'
    if( 'kp' %in% names(pars) )
    {
      # single parameter case is always the range parameter rho
      kp.nm = names( pkern_corr(pars[['k']])[['kp']] )
      if( is.null(kp.nm) ) kp.nm = 'rho'

      # collapse name value pairs and parenthesize result
      kp.parts = mapply(\(nm, p) paste(nm, '=', format(p, digits=nsig)), nm=kp.nm, p=pars[['kp']])
      kp = paste0('[', paste(kp.parts, collapse=', '), ']')
    }

    # return the vector
    return( list(k=k, kp=kp) )
  }

  # handle 2-dimensional case by recursive call, concatenate results
  xy.list = lapply(pars[c('y', 'x')], pkern_toString)
  k = paste(sapply(xy.list, \(xy) xy[['k']]), collapse=' x ')
  kp = paste(sapply(xy.list, \(xy) xy[['kp']]), collapse=' x ')

  # make a title expression with bolding, adding nugget if supplied
  main = bquote(bold(.(k)))
  nug = pars[['nug']]
  if( !is.null(nug) ) nug = paste0('(nugget = ', format(nug, digits=nsig), ')')

  # return as named vector
  return( list(k=k, kp=kp, main=main, nug=nug) )
}


#' Plot gridded datasets
#'
#' A wrapper for graphics::filled.contour and graphics::image (work in progress). This
#' essentially does the same thing as `raster::plot`, but specialized for the pkern package
#'
#' `gdim` and `gyx` can be omitted when their values can be derived from `z` (eg. when `z`
#' is the output from `pkern_fromRaster` or when is a RasterLayer, but not when it is a vector).
#'
#' @param z numeric vector, matrix, RasterLayer, or list
#' @param gdim c(ni, nj), the number of rows and columns in the grid
#' @param gyx list of numeric vectors, the y and x grid line coordinates
#' @param ds positive numeric vector, distance scaling factors along y and x directions
#' @param ppars list of optional plotting parameters (see details)
#'
#' @return TODO
#' @export
#'
#' @examples
#' # TODO
pkern_plot = function(z, gdim=NULL, gyx=NULL, ds=1, ppars=list())
{
  # duplicate distance scaling factors as needed
  if( length(ds) == 1 ) ds = rep(ds, 2)

  # extract grid size from matrix/raster input
  if( any( c('matrix', 'RasterLayer') %in% class(z) ) )
  {
    # we use order (ny, nx)
    gdim = dim(z)[1:2]
    if( 'RasterLayer' %in% class(z) )
    {
      if( is.null(gyx) ) gyx = pkern_fromRaster(z, what='gyx')
      if( is.null(gdim) ) gdim = pkern_fromRaster(z, what='gdim')
      z = pkern_fromRaster(z, what='gval')
    }
  }

  # redirect list input
  if( is.list(z) )
  {
    # output from pkern_fromRaster handled by recursive call
    if( 'gval' %in% names(z) )
    {
      if( !is.null( z[['gdim']] ) ) gdim = z[['gdim']]
      if( !is.null( z[['gyx']] ) ) gyx = z[['gyx']]
      return(pkern_plot(z[['gval']], gdim, gyx, ds, ppars=ppars))
    }

    # otherwise assume list is kernel parameters and pass them on to `pkern_kplot`
    return(pkern_kplot(z, gdim=gdim, ds=ds, ppars=ppars))
  }

  # set up default grid line values when they are not supplied
  yxnm = stats::setNames(nm=c('y', 'x'))
  if( is.null(gyx) ) gyx = lapply(gdim, seq)
  if( length(gyx) != 2 ) stop('could not determine grid size for input z. Try setting gdim')
  if( is.null(names(gyx)) ) gyx = stats::setNames(gyx, yxnm)
  if( is.null(gdim) ) gdim = sapply(gyx, length)

  # unpack plotting parameters and/or set defaults
  asp = ifelse( is.null(ppars[['asp']]), 1, ppars[['asp']])
  leg = ifelse( is.null(ppars[['leg']]), 'z', ppars[['leg']])
  xlab = ifelse( is.null(ppars[['xlab']]), 'x', ppars[['xlab']])
  ylab = ifelse( is.null(ppars[['ylab']]), 'y', ppars[['ylab']])
  glcol = ifelse( is.null(ppars[['glcol']]), NA, ppars[['glcol']])
  pal = ifelse( is.null(ppars[['pal']]), 'Spectral', ppars[['pal']])
  bcol = ifelse( is.null(ppars[['bcol']]), 'black', ppars[['bcol']])
  fcol = ifelse( is.null(ppars[['fcol']]), 'black', ppars[['fcol']])
  discrete = ifelse( is.null(ppars[['discrete']]), FALSE, ppars[['discrete']])
  reset = ifelse( is.null(ppars[['reset']]), TRUE, ppars[['reset']])
  aggpars = ifelse( is.null(ppars[['aggpars']]), NA, ppars[['aggpars']])
  nmax = ifelse( is.null(ppars[['nmax']]), 1e5, ppars[['nmax']])
  smoothed = ifelse( is.null(ppars[['smoothed']]), FALSE, ppars[['smoothed']])

  # determine aggregation factors along x and y directions
  afact = ceiling( gdim/sqrt(nmax) )
  if( any( afact > 1 ) )
  {
    # aggregate the data as needed
    z.agg = pkern_agg(z, gdim, pars=aggpars, afact=afact, discrete=discrete)
    z = z.agg[['z']]
    ij = z.agg[['ij']]

    # update grid size and grid line locations
    gdim = sapply(ij, length)
    gyx = mapply(\(pos, idx) pos[idx], pos=gyx, idx=ij, SIMPLIFY=FALSE)
  }

  # convert vectorization order for compatibility with `graphics::image` and friends
  zmat = matrix(z[pkern_r2c(gdim, in.byrow=FALSE, out.byrow=TRUE, flipx=TRUE)], rev(gdim))

  # adjust gyx for supplied distance scaling and sort into ascending order
  gyx = lapply(gyx, sort)
  gyx = stats::setNames(mapply(\(d, g) g*d, g=gyx, d=ds, SIMPLIFY=FALSE), yxnm)

  # find grid line spacing halves on plot and positions of grid lines
  yxsp = sapply(gyx, \(g) diff(g[1:2]) / 2 )
  yxgl = lapply(yxnm, \(d) c( min(gyx[[d]]) - yxsp[[d]], gyx[[d]] + yxsp[[d]] ) )
  yxlim = lapply(yxgl, range)

  # set up a default continuous color palette
  zr = range(z, na.rm=TRUE, finite=TRUE)
  lvls = pretty(z, 1e2)

  # set up breakpoints in colorbar for discrete case
  if( discrete )
  {
    # TODO: replace this slow operation
    dlvls = sort(unique(z, na.rm=TRUE))

    # make a new z range with padding so that colorbar ticks don't wind up on the edge
    dpad = min(diff(dlvls))/2
    zr = zr + ( dpad * c(-1, 1) )

    # new levels are at midpoints between discrete values
    lvls = c(dlvls[1] - dpad, dlvls + c(diff(dlvls)/2, dpad))

  } else {

    # in continuous case we just pick a large number of breaks covering the range
    zr = range(z, na.rm=TRUE, finite=TRUE)
    lvls = pretty(z, 5e2)
  }

  # set up a color palette
  cols = grDevices::hcl.colors(length(lvls)-1, pal, rev=TRUE)

  # backing up current settings before we change anything
  old.par = graphics::par(no.readonly=TRUE)

  # below is based on `graphics::filled.contour`

  # new plot window layout will use 2 lines to the right of the plot for color bar
  w = ( 2 + old.par[['mar']][2] ) * ( 2.54 * graphics::par('csi') )
  graphics::layout(matrix(c(2, 1), ncol=2), widths = c(1, graphics::lcm(w)))

  # set margins and horizontal text for the colorbar labels
  graphics::par(las=1, mar=c(old.par[['mar']][1], 1, old.par[['mar']][c(3,2)]))

  # draw the color bar then draw a box around it and add legend title
  graphics::plot.new()
  graphics::plot.window(xlim=0:1, ylim=zr, xaxs='i', yaxs='i')
  graphics::rect(0, lvls[-length(lvls)], 1, lvls[-1], col=cols, border=NA)
  graphics::abline(v=1)
  graphics::axis(4, lwd=0, lwd.ticks=1)
  graphics::abline( v = 0 )
  graphics::abline( h = zr )
  graphics::mtext(side=2, line=0, adj=0, padj=-1, at=zr[2], leg)

  # prepare new plot window for data heatmap
  graphics::par(las=0, mar=c(old.par[['mar']][1:3], 1))

  # make the non-smoothed image plot
  if( !smoothed )
  {
    # draw the image using base graphics
    graphics::image(x=gyx[['x']], y=gyx[['y']], zmat, axes=FALSE, ann=FALSE, asp=asp, col=cols)

    # add grid lines
    graphics::abline( v = yxgl[['x']], col=glcol)
    graphics::abline( h = yxgl[['y']], col=glcol)
  }

  # make the smoothed (contour) plot
  if( smoothed )
  {
    # create plot window and coordinate system for bare bones version of graphics::filled.contour
    graphics::plot.new()
    graphics::plot.window(yxlim[['x']], yxlim[['y']], xaxs='i', yaxs='i', asp=asp)

    # call to draw contour plot
    graphics::.filled.contour(x=gyx[['x']], y=gyx[['y']], zmat, lvls, cols)
  }

  # draw x axis
  if( !is.na(xlab) )
  {
    xaxt = graphics::axTicks(1)
    graphics::axis(1, lwd=0, lwd.ticks=1, at=xaxt)
    graphics::title(xlab=xlab)
  }

  # draw y axis
  if( !is.na(ylab) )
  {
    yaxt = graphics::axTicks(2)
    graphics::axis(2, lwd=0, lwd.ticks=1, at=yaxt)
    graphics::title(ylab=ylab)
  }

  # add main title and subtitle
  tline = graphics::par('mgp')[1]-1
  graphics::mtext(ppars[['main']], side=3, line=tline, adj=0, font=2)
  graphics::mtext(ppars[['sub']], side=3, line=tline-1, adj=0)

  # draw inner frame around grid and outer frame (flush with axes)
  graphics::rect(min(yxgl[['x']]), min(yxgl[['y']]), max(yxgl[['x']]), max(yxgl[['y']]), border=bcol)
  if( !is.na(fcol) ) graphics::box(col=fcol)

  # restore old plot settings and finish
  if(reset) graphics::par(old.par)
  return( invisible() )
}


#' Make a heatmap of covariances/correlations around a grid's central point
#'
#' This function displays a separable kernel by building a grid of size `gdim`,
#' assigning to each cell the covariance (or correlation) with the central cell,
#' then passing the result to `graphics::image` for visualization. Covariance
#' parameters from both kernels are printed as a subtitle, with values rounded
#' to 3 decimal places.
#'
#' `pars` should be a list containing elements "y" and "x", which are lists of parameters
#' for the y and x component kernels (recognized by `pkern_corr`). Optionally, `pars` can
#' also include a (positive numeric) pointwise variance and nugget effect (numerics "v"
#' and "nug"). If `pars` contains an element named "ds", it supercedes any argument to `ds`.
#' If `v` and `nug` are supplied, the legend title is set to "covariance", otherwise it is
#' set to "correlation".
#'
#' When `gdim` is omitted, it is set to the smallest size such that the plot shows the
#' effective range of the kernel components in both the y and x directions. if either
#' entry of `gdim` is an even number, it incremented by 1 so that "central point"
#' is unambiguous.
#'
#' @param pars list of two parameter lists "y" and "x", and (optionally) "v", "nug"
#' @param gdim c(ni, nj), the number of y and x lags to include
#' @param ds c(dy, dx), the y and x distances between adjacent grid lines
#' @param ppars list of plotting arguments (passed to `pkern_plot`)
#'
#' @return returns nothing but prints a contour plot of the specified kernel
#' @export
#'
#' @examples
#' pars = pkern_bds(c('gxp', 'mat'))
#' pkern_kplot(pars)
#' gdim = c(30,50)
#' pkern_kplot(pars, gdim)
#' pkern_kplot(c(pars, list(v=2, nug=0.5)), gdim)
#' pkern_kplot(c(pars, list(v=2, nug=0.5)), gdim, ppars=list(smoothed=TRUE))
pkern_kplot = function(pars, gdim=NULL, ds=1, ppars=list())
{
  # duplicate distance scaling factors as needed
  if( 'ds' %in% names(pars) ) ds = pars[['ds']]
  if( length(ds) == 1 ) ds = rep(ds, 2)

  # determine if this is a correlation or covariance plot
  is.nug = !is.null(pars[['nug']])
  is.corr = is.null(pars[['v']]) & !is.nug

  # unpack v and nug, setting defaults as needed
  nug = ifelse( is.null(pars[['nug']]), 0, pars[['nug']])
  v = ifelse( is.null(pars[['v']]), 1, pars[['v']])
  ktype = ifelse(is.corr, 'correlation', 'covariance')

  # set up a plot title and subtitle
  titles = pkern_toString(pars)

  # when grid size not indicated, set it to show the effective range
  if( is.null(gdim) ) gdim = mapply( \(p, s) round(2*p[['kp']][1]/pkern_rho(0.05, p, s)),
                                     p=pars[c('y', 'x')],
                                     s=ds)

  # increment gdim as needed and find the row, column, and index of central cell
  gdim = gdim + 1 - (gdim %% 2)
  ijc = 1 + ( ( stats::setNames(gdim, c('i', 'j')) - 1 ) / 2 )
  idx.ijc = pkern_mat2vec(ijc, gdim)

  # compute the required component kernel values and their kronecker product
  ky = pkern_corrmat(pars[['y']], gdim[1], ds=ds[1], i=ijc[1])
  kx = pkern_corrmat(pars[['x']], gdim[2], ds=ds[2], j=ijc[2])
  z = nug + ( v * as.vector( kronecker(kx, ky) ) )

  # set up axis labels and titles
  gyx = Map( \(d, s, idx) s * ( seq(d) - idx ), d=gdim, s=ds, idx=ijc)
  ppars.default = list(main=bquote(.( titles[['main']] )~'kernel'~.( titles[['nug']] )),
                       sub=titles[['kp']],
                       leg=ktype,
                       xlab='dx',
                       ylab='dy',
                       glcol='black')

  # draw the plot and finish
  ppars.out = utils::modifyList(ppars.default, ppars)
  pkern_plot(z, gdim, gyx=gyx, ppars=ppars.out)
  return(invisible())
}


#' Plot 1-dimensional empirical semivariograms for a separable model
#'
#' Plots the output of `pkern_vario` or `pkern_xvario`, optionally includind the theoretical
#' semivariogram values for a model specified in `pars` (eg the return value of `pkern_vario_fit`).
#'
#' `pars` should be a list containing named elements: "y" and "x", lists of y and x kernel
#' parameters in form recognized by `pkern_corr`; "v", the pointwise variance in the absence
#' of a nugget effect and, optionally, "nug", the variance of the nugget effect.
#'
#' `plotpars` applies to all plots. It should be a list containing any of the following elements:
#'
#'    "title" character or vector of two, title(s) for the plot(s)
#'    "shade" logical, indicating to shade points according to the number of samples
#'    "dmax" positive numeric, upper limit for x in scatterplot
#'    "svlim" length-2 numeric vector, the y limits for the scatterplot
#'
#' @param vario list of semivariance data, the return value of `pkern_vario` or `pkern_xvario`
#' @param pars list of kernel parameters "x" and/or "y", variance "v", and optionally nugget "nug"
#' @param plotpars list, plotting options (see details)
#'
#' @return prints a plot and returns nothing
#' @export
#'
#' @examples
#' gdim = c(100, 50)
#' pars = pkern_corr('gau') |> utils::modifyList(list(v=1))
#' vec = pkern_sim(pars, gdim)
#'
#' # sample first 10 lags and plot
#' vario = pkern_vario(gdim, vec, sep=seq(10))
#'
#' pkern_vario_plot(vario[['x']], pars)
#' pkern_vario_plot(vario[['y']], pars)
#' pkern_vario_plot(vario, pars)
pkern_vario_plot = function(vario, pars=NULL, plotpars=NULL)
{
  # intialize empty list `plotpars` as needed and set some defaults
  if( is.null(plotpars) ) plotpars = list()
  if( is.null(plotpars[['shade']]) ) plotpars[['shade']] = TRUE

  # identify the semivariance directions provided in vario
  nm.yx = c(y='y', x='x')
  nm.dyx = c('d1', 'd2')
  nm.vario = c(nm.yx, nm.dyx)
  is.vario = stats::setNames(nm.vario %in% names(vario), nm.vario)

  # handle 2-dimensional case
  if( all( is.vario[nm.yx] ) )
  {
    # clean up pars if supplied
    if( !is.null(pars) )
    {
      # if a single kernel supplied, take its product
      if( all( c('k', 'kp') %in% names(pars) ) ) pars = list(y=pars, x=pars)

      # assume order "y" "x" unless otherwise indicated by names
      if( !all( nm.yx %in% names(pars) ) ) names(pars) = nm.yx
      if( is.null(pars[['v']]) ) pars[['v']] = pars[['y']][['v']]
      if( is.null(pars[['nug']]) ) pars[['nug']] = pars[['y']][['nug']]
    }

    # split parameters into component lists
    xpars = c( pars[['x']], list( v=pars[['v']], nug=pars[['nug']] ) )
    ypars = c( pars[['y']], list( v=pars[['v']], nug=pars[['nug']] ) )

    # set default plot titles
    title.def = FALSE
    if( is.null(plotpars[['title']]) )
    {
      title.def = TRUE
      plotpars[['title']] = nm.yx
    }

    # set default common x-axis limit for the scatterplots
    dmax.def = FALSE
    if( is.null(plotpars[['dmax']]) )
    {
      dmax.def = TRUE
      sepmax = sapply( vario[nm.yx], \(xy) max( xy[['sep']] ) )
      plotpars[['dmax']] = max(sepmax)
    }

    # same for y-axis
    svlim.def = FALSE
    if( is.null(plotpars[['svlim']]) )
    {
      svlim.def = TRUE
      idx.valid = lapply( vario[nm.yx], \(xy) xy[['sep']] < plotpars[['dmax']] )
      svmax = mapply(\(yx, i) max( yx[['sv']][i] ), yx=vario[nm.yx], i=idx.valid)
      plotpars[['svlim']] = c(0, max(svmax))
    }

    # update default plot settings for diagonals
    if( all( is.vario ) )
    {
      # update axis limits as needed
      if( svlim.def )
      {
        # x axis
        sepmax = sapply( vario[nm.dyx], \(xy) max( xy[['sep']] ) )
        if( dmax.def ) plotpars[['dmax']] = max(plotpars[['dmax']], sepmax)

        # y axis
        idx.valid = lapply( vario[nm.vario], \(xy) xy[['sep']] < plotpars[['dmax']] )
        svmax = mapply(\(xy, i) max( xy[['sv']][i] ), xy=vario[nm.vario], i=idx.valid)
        plotpars[['svlim']] = range(plotpars[['svlim']], svmax)
      }

      # assign plot titles as needed
      if( title.def )
      {
        dstr = paste0('(', round( atan( vario[['ds']][1] / vario[['ds']][2] )*180/pi, 1), ' degrees)')
        plotpars[['title']] = c(plotpars[['title']], d1=paste('y', dstr), d2=paste('x', dstr))
      }
    }

    # split `plotpars` into its components
    nm.toplot = stats::setNames(nm=names(is.vario)[is.vario])
    plotpars.list = lapply(nm.toplot, \(pp) {
      utils::modifyList(plotpars, list(title=plotpars[['title']][pp]))
    })

    # initialize either a 2 or 4 pane layout
    graphics::par(mfrow = c(sum(is.vario)/2, 2))

    # recursive calls to add "y" and "x" plots
    pkern_vario_plot(vario[['y']], pars=ypars, plotpars=plotpars.list[['y']])
    pkern_vario_plot(vario[['x']], pars=xpars, plotpars=plotpars.list[['x']])

    # recursive calls to add diagonal plots
    if( all( is.vario ) )
    {
      ds = vario[['ds']]
      pkern_vario_plot(c(vario[['d1']], list(ds=ds)), pars=pars, plotpars=plotpars.list[['d1']])
      pkern_vario_plot(c(vario[['d2']], list(ds=ds)), pars=pars, plotpars=plotpars.list[['d2']])
    }

    # reset plot panel layout before finishing, returning nothing
    graphics::par(mfrow=c(1,1))
    return(invisible())
  }

  # 1-dimension case:

  # unpack plotting limits
  if( is.null(plotpars[['dmax']]) ) plotpars[['dmax']] = max(vario[['sep']])
  xlim = c(0, plotpars[['dmax']])
  ylim = plotpars[['svlim']]
  main = plotpars[['title']]

  # catch empty variogram
  if( length(vario[['n']]) == 1 )
  {
    # initialize an empty plot
    plot(xlim, ylim, pch=NA, xlab='lag', ylab='semivariance', main=main)

  } else {

    # unpack arguments (exclude the zero point)
    n = vario[['n']][-1]
    sep = vario[['sep']][-1]
    sv = vario[['sv']][-1]

    # map point shading to the number of point pairs sampled if requested
    #colmap = 'black'
    #if( plotpars[['shade']] ) colmap = rev( grDevices::gray.colors( max(n) ) )[n]

    # make the scatterplot of sampled semivariances
    plot(sep, sv, xlab='lag', ylab='semivariance', xlim=xlim, ylim=ylim, main=main, pch=16)

    # add bubbles indicating sample sizes
    points(sep, sv, pch=1, cex=n/max(n))

    # # make the scatterplot of sampled semivariances
    # plot(sep, sv, xlab='lag', ylab='semivariance',
    #      xlim=xlim, ylim=ylim, main=main, pch=20, col=colmap)

  }

  # if kernel parameters are supplied, plot the theoretical values
  if( !is.null(unlist(pars)) )
  {
    # check for variance components in parameters list and set default nugget (0)
    nug = pars[['nug']]
    v = pars[['v']]
    if( is.null(v) ) stop('variance v not found in pars')
    if( is.null(nug) ) nug = 0

    #  handle x or y kernel component requests
    sep = seq(xlim[1], xlim[2], length.out=1e2)
    if( all( c('k', 'kp') %in% names(pars) ) ) fit = pkern_tvario(pars, v=v, nug=nug, d=sep)

    # handle requests along a diagonal (both y and x components included in pars)
    if( all( nm.yx %in% names(pars) ) )
    {
      # find unit distance along diagonals
      ds = vario[['ds']]
      udist = sqrt(sum(ds^2))

      # find the equivalent displacement along y and x directions
      sepy = ds[1] * sep / udist
      sepx = ds[2] * sep / udist

      # compute appropriately scaled kernel values
      fity = pkern_tvario(pars[['y']], v=sqrt(v), d=sepy)
      fitx = pkern_tvario(pars[['x']], v=sqrt(v), d=sepx)
      fit = nug + ( fitx * fity )
    }

    # add line plot showing semivariance at requested lags
    graphics::lines(sep, fit)
  }

  return(invisible())
}

