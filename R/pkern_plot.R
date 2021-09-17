#
# pkern_raster.R
# Dean Koch, Sep 2021
# Functions for plotting objects created with pkern
#


# wrapper for graphics::filled.contour and graphics::image (work in progress)
# compute coordinates and reshape z into a matrix oriented for `graphics::image` and friends
pkern_plot = function(z, dims=NULL, xy=NULL, ds=1, smoothed=FALSE, maxx=300, ppars=list())
{
  # extract grid size from matrix/raster input
  if( any( c('matrix', 'RasterLayer') %in% class(z) ) )
  {
    # we use order (nx, ny)
    dims = rev(dim(z)[1:2])
    if( 'RasterLayer' %in% class(z) )
    {
      if( is.null(xy) ) xy = pkern_fromRaster(z, what='xy')
      z = pkern_fromRaster(z, what='values')
    }
  }

  # redirect list input (interpreted as kernel parameters) to pkern_kplot
  if( is.list(z) ) return(pkern_kplot(z, dims=dims, ds=ds, smoothed=smoothed, maxx=maxx, ppars=ppars))

  # set up default grid line values when they are not supplied
  xynm = setNames(nm=c('x', 'y'))
  if( is.null(xy) ) xy = lapply(dims, seq)
  if( length(xy) != 2 ) stop('could not determine grid size for input z. Try setting dims')
  if( is.null(names(xy)) ) xy = setNames(xy, xynm)
  if( is.null(dims) ) dims = sapply(xy, length)

  # determine aggregation factors along x and y directions
  afact = ceiling( dims/maxx )
  if( any( afact > 1 ) )
  {
    # aggregate the data as needed
    z.agg = pkern_agg(z, dims, afact, tol=0.99)
    z = z.agg[['z']]
    ij = z.agg[['ij']]
    dims.sg = unname(rev(sapply(ij, length)))

    # update grid size and grid line locations
    dims = dims.sg
    xy = list(x=xy[['x']][ ij[['j']] ], y=xy[['y']][ ij[['i']] ])
  }

  # convert data to row-vectorized order for compatibility with `graphics::image` and friends
  zmat = matrix(z[pkern_r2c(dims, in.byrow=FALSE, out.byrow=TRUE, flipx=TRUE)], dims[1])

  # find grid line spacing halves on plot and positions of grid lines
  xysp = sapply(xy, \(g) diff(g[1:2]) / 2 )
  xygl = lapply(xynm, \(d) c(xy[[d]][1] - xysp[[d]], xy[[d]] + xysp[[d]]) )
  xylim = lapply(xygl, range)

  # unpack plotting parameters and/or set defaults
  asp = ifelse( is.null(ppars[['asp']]), 1, ppars[['asp']])
  leg = ifelse( is.null(ppars[['leg']]), 'z', ppars[['leg']])
  xlab = ifelse( is.null(ppars[['xlab']]), 'x', ppars[['xlab']])
  ylab = ifelse( is.null(ppars[['ylab']]), 'y', ppars[['ylab']])
  glcol = ifelse( is.null(ppars[['glcol']]), NA, ppars[['glcol']])
  pal = ifelse( is.null(ppars[['pal']]), 'Spectral', ppars[['pal']])
  bcol = ifelse( is.null(ppars[['bcol']]), 'black', ppars[['bcol']])
  fcol = ifelse( is.null(ppars[['fcol']]), 'black', ppars[['fcol']])

  # set up color palette
  zr = range(z, na.rm=TRUE, finite=TRUE)
  lvls = pretty(z, 1e2)
  cols = hcl.colors(length(lvls)-1, pal, rev = TRUE)

  # set up plot window ratios
  rd = sapply(xy, \(g) diff(range(g))) * rr
  rdd = diff(rd)

  # backing up current settings before we change anything
  old.par = par(no.readonly=TRUE)

  # new plot window layout will use 2 lines to the right of the plot for color bar
  w = ( 2 + old.par[['mar']][2] ) * ( 2.54 * par('csi') )
  layout(matrix(c(2, 1), ncol=2), widths = c(1, lcm(w)))

  # set margins and horizontal text for the colorbar labels
  par(las=1, mar=c(old.par[['mar']][1], 1, old.par[['mar']][c(3,2)]))

  # draw the color bar then draw a box around it and add legend title
  plot.new()
  plot.window(xlim=0:1, ylim=zr, xaxs='i', yaxs='i')
  rect(0, lvls[-length(lvls)], 1, lvls[-1], col=cols, border=NA)
  abline(v=1)
  axis(4, lwd=0, lwd.ticks=1)
  abline( v = 0 )
  abline( h = zr )
  mtext(side=2, line=0, adj=0, padj=-1, at=zr[2], leg)

  # prepare new plot window for data heatmap
  par(las=0, mar=c(old.par[['mar']][1:3], 1))

  # make the image plot
  if( !smoothed )
  {
    # draw the image using base graphics
    graphics::image(x=xy[['x']], y=xy[['y']], zmat, axes=FALSE, ann=FALSE, asp=asp, col=cols)

    # add grid lines
    abline( v = xygl[['x']], col=glcol)
    abline( h = xygl[['y']], col=glcol)
  }

  # make the smoothed (contour) plot
  if( smoothed )
  {
   # initialize plot window for barebones version of graphics::filled.contour
    plot.new()
    plot.window(xylim[['x']], xylim[['y']], xaxs='i', yaxs='i', asp=asp)

    # call to draw contour plot
    graphics::.filled.contour(x=xy[['x']], y=xy[['y']], zmat, lvls, cols)
  }

  # draw x axis
  if( !is.na(xlab) )
  {
    xaxt = axTicks(1)
    axis(1, lwd=0, lwd.ticks=1, at=xaxt)
    title(xlab=xlab)
  }

  # draw y axis
  if( !is.na(ylab) )
  {
    yaxt = axTicks(2)
    axis(2, lwd=0, lwd.ticks=1, at=yaxt)
    title(ylab=ylab)
  }

  # add main title and subtitle
  tline = par('mgp')[1]-1
  mtext(ppars[['main']], side=3, line=tline, adj=0, font=2)
  mtext(ppars[['sub']], side=3, line=tline-1, adj=0)

  # draw inner frame around grid and outer frame (flush with axes)
  rect(min(xygl[['x']]), min(xygl[['y']]), max(xygl[['x']]), max(xygl[['y']]), border=bcol)
  if( !is.na(fcol) ) box(col=fcol)

  # restore old plot settings and finish
  par(old.par)
  return( invisible() )
}


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
#' pkern_toString(modifyList(pars, list(nug=0.5)))
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
  xy.list = lapply(pars[c('x', 'y')], pkern_toString)
  k = paste(sapply(xy.list, \(xy) xy[['k']]), collapse=' x ')
  kp = paste(sapply(xy.list, \(xy) xy[['kp']]), collapse=' x ')

  # make a title expression with bolding, adding nugget if supplied
  main = bquote(bold(.(k)))
  nug = pars[['nug']]
  if( !is.null(nug) ) nug = paste0('(nugget = ', format(nug, digits=nsig), ')')

  # return as named vector
  return( list(k=k, kp=kp, main=main, nug=nug) )
}


#' Make a heatmap of covariances/correlations around a grid's central point
#'
#' This function displays a separable kernel by building a grid of size `dims`,
#' assigning to each cell the covariance (or correlation) with the central cell,
#' then passing the result to `graphics::image` for visualization.
#'
#' `pars` should be a list containing elements "x" and "y", which are lists of parameters
#' for the x and y component kernels (recognized by `pkern_corr`). Optionally, `pars` can
#' also include a (positive numeric) pointwise variance term in "v". Covariance parameters
#' from both kernels are printed as a subtitle, with values rounded to 3 decimal places.
#'
#' When `dims` is omitted, it is set to the smallest size such that the plot shows the
#' effective range of the kernel components in both the x and y directions. if either
#' entry of `dims` is an even number, it incremented by 1 in order so that "central"
#' is  unambiguous.
#'
#' If `v` and `nug` are not supplied, plot labels are changed to show "correlation".
#'
#' @param pars list of two parameter lists "x" and "y", and (optionally) "v", "nug"
#' @param dims c(nx, ny), the number of x and y lags to include
#' @param smoothed logical, whether to smooth the image (passed to `pkern_plot`)
#' @param maxx integer or integer vector, plotting resolution (passed to `pkern_plot`)
#' @param ppars list of plotting arguments (passed to `pkern_plot`)
#'
#' @return returns nothing but prints a contour plot of the specified kernel
#' @export
#'
#' @examples
#' pars = pkern_bds(c('gxp', 'mat'))
#' pkern_kplot(pars)
#' pkern_kplot(pars, dims=c(50,30))
#' pkern_kplot(c(pars, list(v=2, nug=0.5)), dims=c(50,30))
#' pkern_kplot(c(pars, list(v=2, nug=0.5)), dims=c(50,30), smoothed=TRUE)
pkern_kplot = function(pars, dims=NULL, ds=1, smoothed=FALSE, maxx=200, ppars=list())
{

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
  if( is.null(dims) )
  {
    nx = round(2*pars[['x']][['kp']][1]/pkern_rho(0.05, pars[['x']], ds[1]))
    ny = round(2*pars[['y']][['kp']][1]/pkern_rho(0.05, pars[['y']], ds[2]))
    dims = c(nx, ny)
  }

  # increment dims as needed and find the row, column, and index of central cell
  dims = dims + 1 - (dims %% 2)
  ij.central = setNames(rev( 1 + ( (dims - 1) / 2 ) ), c('i', 'j'))
  idx.central = pkern_mat2vec(ij.central, dims[2])

  # compute the required component kernel values and their kronecker product
  kx = pkern_corrmat(pars[['x']], dims[1], ds=ds[1], j=ij.central['j'])
  ky = pkern_corrmat(pars[['y']], dims[2], ds=ds[2], i=ij.central['i'])
  z = nug + ( v * as.vector( kronecker(kx, ky) ) )

  # set up axis labels and titles
  x = ds[1] * ( seq(dims[1]) - ij.central['j'] )
  y = ds[2] * ( seq(dims[2]) - ij.central['i'] )
  ppars.default = list(main=bquote(.( titles[['main']] )~'kernel'~.( titles[['nug']] )),
                       sub=titles[['kp']],
                       leg=ktype,
                       xlab='dx',
                       ylab='dy',
                       glcol='black')

  # draw the plot and finish
  pkern_plot(z, dims, xy=list(x,y),
             smoothed=smoothed,
             maxx=maxx,
             ppars=modifyList(ppars.default, ppars))

  return(invisible())
}


#' Plot 1-dimensional empirical semivariograms for a separable model
#'
#' Plots the output of `pkern_vario` or `pkern_xvario`, optionally includind the theoretical
#' semivariogram values for a model specified in `pars` (eg the return value of `pkern_vario_fit`).
#'
#' `pars` should be a list containing named elements: "x" and "y", lists of x and y kernel
#' parameters in form recognized by `pkern_corr`; "v", the pointwise variance in the absence
#' of a nugget effect and, optionally, "nug", the variance of the nugget effect.
#'
#' `plotpars` applies to all plots. It should be a list containing any of the following elements:
#'
#'    "title" character or vector of two, title(s) for the plot(s)
#'    "shade" logical, indicating to shade points according to the number of samples
#'    "seplim" length-2 numeric vector, the x limits for the scatterplot
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
pkern_vario_plot = function(vario, pars=NULL, plotpars=NULL)
{
  # intialize empty list `plotpars` as needed and set some defaults
  if( is.null(plotpars) ) plotpars = list()
  if( is.null(plotpars[['shade']]) ) plotpars[['shade']] = TRUE

  # identify the semivariance directions provided in vario
  nm.vario = c('x', 'y', 'd1', 'd2')
  is.vario = setNames(nm.vario %in% names(vario), nm.vario)

  # handle 2-dimensional case
  if( all( is.vario[c('x', 'y')] ) )
  {
    # assign names to pars if they are missing
    nm.pars = c('x', 'y', 'v')
    if( !is.null(pars) & !all( nm.pars %in% names(pars) ) ) names(pars) = nm.pars

    # extract componentwise parameters
    xpars = c( pars[['x']], list( v=pars[['v']], nug=pars[['nug']] ) )
    ypars = c( pars[['y']], list( v=pars[['v']], nug=pars[['nug']] ) )

    # set default plot titles
    title.def = FALSE
    if( is.null(plotpars[['title']]) )
    {
      title.def = TRUE
      plotpars[['title']] = c(x='x', y='y')
    }

    # set default (common) x-axis limits for the scatterplots
    seplim.def = FALSE
    if( is.null(plotpars[['seplim']]) )
    {
      seplim.def = TRUE
      sepmax = sapply( vario[c("x", "y")], \(xy) max( xy[['sep']] ) )
      plotpars[['seplim']] = c(0, max(sepmax))
    }

    # same for y-axis
    svlim.def = FALSE
    if( is.null(plotpars[['svlim']]) )
    {
      svlim.def = TRUE
      svmax = sapply( vario[c("x", "y")], \(xy) max( xy[['sv']] ) )
      plotpars[['svlim']] = c(0, max(svmax))
    }

    # update default plot settings for diagonals
    if( all( is.vario ) )
    {
      # update x-axis limits as needed
      if( svlim.def )
      {
        sepmax = sapply( vario[c("d1", "d2")], \(xy) max( xy[['sep']] ) )
        plotpars[['seplim']][2] = max(plotpars[['seplim']][2], sepmax)
      }

      # update y-axis limits as needed
      if( svlim.def )
      {
        svmax = sapply( vario[c("d1", "d2")], \(xy) max( xy[['sv']] ) )
        plotpars[['svlim']][2] = max(plotpars[['svlim']][2], svmax)
      }

      # assign plot titles as needed
      if( title.def )
      {
        dstr = paste0('(', round( atan( vario[['ds']][2] / vario[['ds']][1] )*180/pi, 1), ' degrees)')
        plotpars[['title']] = c(plotpars[['title']], d1=paste('x', dstr), d2=paste('y', dstr))
      }
    }

    # split `plotpars` into its components
    nm.toplot = setNames(nm=names(is.vario)[is.vario])
    plotpars.list = lapply(nm.toplot, \(pp) {
      modifyList(plotpars, list(title=plotpars[['title']][pp]))
    })

    # initialize either a 2 or 4 pane layout
    par( mfrow = c(sum(is.vario)/2, 2) )

    # recursive calls to add "x" and "y" plots
    pkern_vario_plot(vario[['x']], pars=xpars, plotpars=plotpars.list[['x']])
    pkern_vario_plot(vario[['y']], pars=ypars, plotpars=plotpars.list[['y']])

    # recursive calls to add diagonal plots
    if( all( is.vario ) )
    {
      ds = vario[['ds']]
      pkern_vario_plot(c(vario[['d1']], list(ds=ds)), pars=pars, plotpars=plotpars.list[['d1']])
      pkern_vario_plot(c(vario[['d2']], list(ds=ds)), pars=pars, plotpars=plotpars.list[['d2']])
    }

    # reset plot panel layout before finishing, returning nothing
    par(mfrow=c(1,1))
    return(invisible())
  }

  # 1-dimension case:

  # unpack plotting limits
  xlim = plotpars[['seplim']]
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
    colmap = 'black'
    if( plotpars[['shade']] ) colmap = rev( gray.colors( max(n) ) )[n]

    # make the scatterplot of sampled semivariances
    plot(sep, sv, xlab='lag', ylab='semivariance',
         xlim=xlim, ylim=ylim, main=main, pch=16, col=colmap)

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

    # handle requests along a diagonal (both x and y components included in pars)
    if( all( c('x', 'y') %in% names(pars) ) )
    {
      # find unit distance along diagonals
      ds = vario[['ds']]
      udist = sqrt(sum(ds^2))

      # find the equivalent displacement along x and y directions
      sepx = ds[1] * sep / udist
      sepy = ds[2] * sep / udist

      # compute appropriately scaled kernel values
      fitx = pkern_tvario(pars[['x']], v=sqrt(v), d=sepx)
      fity = pkern_tvario(pars[['y']], v=sqrt(v), d=sepy)
      fit = nug + ( fitx * fity )
    }

    # add line plot showing semivariance at requested lags
    lines(sep, fit)
  }

  return(invisible())
}

