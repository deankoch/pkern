#
# pkern_raster.R
# Dean Koch, Sep 2021
# Functions for plotting objects created with pkern
#


# wrapper for graphics::filled.contour and graphics::image (work in progress)

pkern_plot = function(z, dims=NULL, x=NULL, y=NULL, smoothed=FALSE,...)
{
  # extract dims from matrix input, vectorize
  if( is.matrix(z) )
  {
    dims = rev(dim(z))
    z = c(z)
  }

  # vector z is assumed to be in column vectorized order - convert to row-vectorized
  if( is.vector(z) ) z = z[pkern_r2c(dims, in.byrow=FALSE, out.byrow=TRUE, flipx=TRUE)]

  # default grid line labels
  if( is.null(x) ) x = seq(dims[1])
  if( is.null(y) ) y = seq(dims[2])

  # set up color palette
  zr = range(z)
  lvls = pretty(z, 1e2)
  cols = hcl.colors(length(lvls)-1, "YlOrRd", rev = TRUE)

  # backing up current settings then set up new plot window layout
  old.par = par(c('mar', 'las', 'mfrow'))
  w = ( 3 + old.par[['mar']][2] ) * ( 2.54 * par('csi') )
  layout(matrix(c(2, 1), ncol=2), widths = c(1, lcm(w)))
  par(las=1, mar = c(old.par[['mar']][1], 1, old.par[['mar']][c(3,2)]))

  # draw the color bar
  plot.new()
  plot.window(xlim=0:1, ylim=zr, xaxs='i', yaxs='i')
  rect(0, lvls[-length(lvls)], 1, lvls[-1], col=cols, border=NA)
  abline(v=1)
  axis(4, lwd=0, lwd.ticks=1)

  # add legend title
  ltitle = 'correlation'
  ltpos = zr[2]
  mtext(side=2, line=0, adj=0, padj=-1, at=ltpos, ltitle)

  # prepare new plot window for data heatmap
  par(mar = c(old.par[['mar']][1:3], 1))

  # make the plot, add subtitle, and finish
  zmat = matrix(z, dims[1])
  if( smoothed )
  {
    plot.new()
    plot.window(range(x), range(y), xaxs='i', yaxs='i')
    graphics::.filled.contour(x=x, y=y, zmat, lvls, cols)
  }

  if( !smoothed ) graphics::image(x=x, y=y, zmat, axes=FALSE, ann=FALSE, col=cols)

  # draw axes
  xaxt = pretty(x)
  yaxt = pretty(y)
  axis(1, lwd=0, lwd.ticks=1, at=xaxt)
  axis(2, lwd=0, lwd.ticks=1, at=yaxt)

  # restore old plot settings
  par(old.par)
}

#' Make a heatmap of covariances/correlations around a grid's central point
#'
#' This function displays a separable kernel by building a grid of size `dims`,
#' assigning to each cell the covariance (or correlation) with the central cell,
#' then passing the result to a contour plotter (`graphics::filled.contour`)
#' for visualization.
#'
#' `pars` should be a list containing elements "x" and "y", which are lists of parameters
#' for the x and y component kernels (recognized by `pkern_corr`). Optionally, `pars` can
#' also include a (positive numeric) pointwise variance term in "v". Covariance parameters
#' from both kernels are printed as a subtitle, with values rounded to 3 decimal places.
#'
#' if either entry of `dims` is an even number, it incremented by 1 in order to
#' make the notion of "central" unambiguous. If `v` is not supplied, it is set to
#' the default 1, and the legend label is changed to "correlation".
#'
#' Note that contour plots (which smooth their data) are not an appropriate tool for
#' visualizing the nugget effect (a discontinuity), so `pars$nug` is ignored by this
#' function.
#'
#' @param pars list of two parameter lists "x" and "y", and (optionally) "v", "nug"
#' @param dims c(nx, ny), the number of x and y lags to include
#' @param v numeric, the pointwise variance
#'
#' @return returns nothing but prints a contour plot of the specified kernel
#' @export
#'
#' @examples
#' pars = pkern_bds(c('exp', 'mat'))
#' pkern_kplot(pars)
#' pkern_kplot(pars, dims=c(10,5))
#' pkern_kplot(c(pars, list(v=2)), dims=c(10,5))
pkern_kplot = function(pars, dims=c(25, 25))
{
  # unpack v, setting defaults as needed, determine plot title
  ktype = ifelse( is.null(pars[['v']]), 'correlation', 'covariance')
  v = ifelse( is.null(pars[['v']]), 1, pars[['v']])

  # generate the kernel name string
  kname = sapply(pars[c('x', 'y')], \(x) paste0('[', x$k, ']') ) |> paste(collapse=' x ')

  # set the plot title
  ptitle = paste(ktype, 'about central grid point with kernel:', kname)

  # extract parameter names
  kpx.nm = names( pkern_corr(pars[['x']][['k']])[['kp']] )
  kpy.nm = names( pkern_corr(pars[['y']][['k']])[['kp']] )
  if( is.null(kpx.nm) ) kpx.nm = 'rho'
  if( is.null(kpy.nm) ) kpy.nm = 'rho'

  # create a subtitle
  kpx = mapply(\(nm, p) paste(nm, '=', round(p, 3)), nm=kpx.nm, p=pars[['x']][['kp']])
  kpy = mapply(\(nm, p) paste(nm, '=', round(p, 3)), nm=kpy.nm, p=pars[['y']][['kp']])
  xstr = paste0('[', paste(kpx, collapse=', '), ']')
  ystr = paste0('[', paste(kpy, collapse=', '), ']')
  stitle = paste('parameters', xstr, ' x ', ystr)

  # increment dims as needed and find the row, column, and index of central cell
  dims = dims + 1 - (dims %% 2)
  ij.central = setNames(rev( 1 + ( (dims - 1) / 2 ) ), c('i', 'j'))
  idx.central = pkern_mat2vec(ij.central, dims[2])

  # compute the required component kernel values and their kronecker product
  kx = pkern_corrmat(pars[['x']], dims[1], j=ij.central['j'])
  ky = pkern_corrmat(pars[['y']], dims[2], i=ij.central['i'])
  z = v * as.vector( kronecker(kx, ky) )

  # compute coordinates and reshape z into a matrix oriented for `graphics::image` and friends
  x = seq(dims[1]) - ij.central['j']
  y = seq(dims[2]) - ij.central['i']
  z.mat = matrix(z[pkern_r2c(dims, in.byrow=FALSE, out.byrow=TRUE)], dims[1])

  # make the plot, add subtitle, and finish
  graphics::filled.contour(x, y, z.mat, xlab='dx', ylab='dy', asp=1, nlevels=50)
  mtext(side=3, line=par('mgp')[1]-1, adj=1, ptitle)
  mtext(side=3, line=par('mgp')[1]-2, adj=1, cex=0.8, stitle)

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

