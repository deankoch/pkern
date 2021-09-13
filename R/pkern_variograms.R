#
# pkern_raster.R
# Dean Koch, Sep 2021
# Functions for making and fitting semivariograms
#


#' Theoretical semivariances along x or y direction
#'
#' Convenience function for getting theoretical semivariance values for a set of
#' spatial lags, assuming (intrinsic) stationarity.
#'
#' Kernel parameters `pars` are passed to the correlation function (`pkern_corr`)
#' to get correlation values at distances `d`. The result is converted to semivariance
#' assuming a pointwise variance of `v + nug` (an optional nugget effect) and returned
#' as a vector the same length as `d`.
#'
#' When `d` is NULL, the function instead returns the semivariance function as an
#' anonymous function of distance.
#'
#' @param pars list of kernel parameters, in form recognized by `pkern_corr`
#' @param v positive numeric, the pointwise variance absent a nugget effect
#' @param nug positive numeric, (optional) the variance of the nugget effect
#' @param d vector of nonegative numeric spatial lags to evaluate
#'
#' @return either a vector the same length as `d`, or an anonymous function of distance
#' @export
#'
#' @examples
#' pars.mat = pkern_corr('mat')
#' pkern_tvario(pars.mat, v=1, nug=0.5, 1:10)
#' sv = pkern_tvario(pars.mat, nug=0.5, v=1)
#' sv(1:10)
pkern_tvario = function(pars, v, nug=0, d=NULL)
{
  fn = function(d) { unname(v) * ( 1 + nug - pkern_corr(pars, d) ) }
  if( is.null(d) ) return(fn)
  return(fn(d))
}



#' Sample semivariogram for regular 2D gridded data along x direction
#'
#' The function identifies and computes the semivariance (in `vec`) of point pairs on a
#' regular grid of dimensions `dims` at the supplied x grid line lags, `sep`. Elements of
#' `sep` must be a subset of `1:(dims[1]-1)` (the default).
#'
#' Results are returned in a named list containing the lags ("sep"), the number of point
#' pairs ("n") and sample semivariances ("sv") at each lag, and optionally (if `simple`
#' is FALSE) a list of matrices containing the indices of point pairs sampled and
#' their absolute difference ("pairs").
#'
#' method "mean" is the classical method of Matheron (1962), "median" and "ch" are the
#' median and robust methods described in Cressie and Hawkins (1980).
#'
#' references:
#'
#' Cressie and Hawkins (1980): https://doi.org/10.1007/BF01035243
#'
#'
#' @param dims vector c(nx, ny) of positive integers, the number of grid lines in full grid
#' @param vec numeric vector of data, in column vectorized order (length `prod(dims)`)
#' @param sep positive integer vector, the grid line lags to sample
#' @param simple logical, if FALSE the function returns a list of point pair indices
#' @param method character, one of "mean", "median", "ch"
#'
#' @return a list with named elements "sep", "n", and "sv"
#' @export
#'
#' @examples
#'
#'
pkern_xvario = function(dims, vec, sep=NA, simple=TRUE, method='median')
{
  # set default sample lags
  if( is.na(sep) ) sep = seq( dims[1] - 1 )

  # identify missing values
  miss = is.na(vec)

  # sample along different grid line lags in a loop
  plist = vector(mode='list', length=length(sep))
  svx = rep(NA, length(sep))
  n = rep(0, length(sep))
  for( idx.dj in seq_along(sep) )
  {
    dj = sep[idx.dj]

    # total number of eligible point pairs
    idx.max = dims[2] * ( dims[1] - dj )

    # TODO: take a random sample of these
    idx.sample = seq( idx.max )

    # vectorization trick to convert this to samples of i,j indices
    ij.sample = pkern_vec2mat(idx.sample, dims[2])

    # translate sample points index into index of vec
    idx1 = ij.sample[,1] + ( dims[2] * (ij.sample[,2] - 1) )
    idx2 = idx1 + ( dims[2] * dj )

    # omit point pairs containing NAs
    idx.omit = miss[idx1] | miss[idx2]
    n[idx.dj] = sum(!idx.omit)
    plist[[idx.dj]] = cbind(idx1=idx1[!idx.omit], idx1=idx2[!idx.omit])

    # estimate semivariance
    if(n[idx.dj] > 0)
    {
      absdiff = apply(plist[[idx.dj]], 1, \(idx) abs(diff(vec[idx])) ) |> as.vector()
      plist[[idx.dj]] = cbind(plist[[idx.dj]], absdiff=absdiff)
      if(method=='mean') svx[idx.dj] = mean( absdiff^2 ) / 2
      if(method=='median') svx[idx.dj] = median( sqrt(absdiff) )^4 / (2 * 0.457)
      if(method=='ch')
      {
        ch.factor = 0.457 + ( 0.494 / n[idx.dj] ) + ( 0.494 / ( n[idx.dj]^2 ) )
        svx[idx.dj] = mean( sqrt(absdiff) )^4 / (2 * ch.factor)
      }
    }
  }

  if( simple ) return( list(sep=sep, n=n, sv=svx) )
  return( list(sep=sep, n=n, sv=svx, pairs=plist) )
}


#' Directional sample semivariograms for regular 2D gridded data
#'
#' A wrapper for calls to `pkern_xvario` to get the directional sample semivariograms
#' from regular lattice data along the x and y dimensions, and optionally along the
#' diagonals.
#'
#' The function makes the appropriate calls to `pkern_xvario` and returns the results
#' in list elements "x", "y". If `diagonal=TRUE`, the function also returns "d1" and
#' "d2", the x and y semivariograms after a 45 degree counterclockwise rotation of the
#' input grid (see `pkern_r45`), assuming equal x and y grid line spacings. For unequal
#' grid line spacings, the true rotation angle is `atan(ds[2]/ds[1])`
#'
#' Note that input `sep` is specified in terms of grid line indices (integers), not
#' distances (their product with the x or y resolution). However, non-integer separation
#' distances can be specified by setting the resolution of the input grid (dx, dy) in
#' argument `ds`. When `ds` is supplied, the output "sep" values are scaled accordingly
#' and `dmax` is interpreted as a distance, rather than a number of grid lines.
#'
#' When `ds` is not supplied, the function assumes that the x grid lines have the same
#' spacing distance (1) as the y grid lines. In this case, diagonally adjacent
#' grid points have separation distance `sqrt(2)`.
#'
#' @param vec numeric vector of data, in column vectorized order (length `prod(dims)`)
#' @param dims vector c(nx, ny) of positive integers, the number of grid lines
#' @param sep list c(dx, dy) of positive integer vectors, grid line lags to sample
#' @param simple logical, if FALSE the function includes point pair indices in return list
#' @param method character, one of "mean", "median", "ch", passed to `pkern_xvario`
#' @param diagonal logical, if TRUE, semivariogram results are returned also for diagonals
#' @param dmax numeric, maximum grid line separation distance (used to filter results)
#' @param ds length-2 vector of positive numeric, the spacing distance of x and y grid lines
#'
#' @return list of `pkern_xvario` output "x", "y" (and "d1", "d2"); sample variance "v"; and "ds"
#' @export
#'
#' @examples
pkern_vario = function(dims, vec, sep=NA, simple=TRUE, method='median', diagonal=TRUE, dmax=NA, ds=NA)
{
  # set default sample lags and compute sample variance
  if( anyNA(sep) ) sep = lapply(dims, \(d) seq(d-1))
  v = var(vec, na.rm=TRUE)

  # set default resolution (1,1)
  if( anyNA(ds) ) ds = c(1,1)

  # compute x semivariances, then repeat in row-vectorized order to get y semivariances
  xout = pkern_xvario(dims, vec, simple=simple, method=method)
  yout = pkern_xvario(rev(dims), vec[pkern_r2c(dims, FALSE, TRUE)], simple=simple, method=method)

  # scale sep by resolution to get distances, then bundle everything into an output list
  xout[['sep']] = ds[1] * xout[['sep']]
  yout[['sep']] = ds[2] * yout[['sep']]
  list.out = list(x=xout, y=yout)

  # add a zero-distance point so that the output sets are never empty
  list.out = lapply(list.out, \(d) lapply(d, \(z) c(0, z) ))

  # repeat with rotated coordinates, if requested
  if( diagonal )
  {
    # rotate the input array by 45 degrees clockwise
    r45.out = pkern_r45(vec, dims)
    vec45 = r45.out[['z']]
    dims45 = r45.out[['dims']]

    # recursive call to get component variograms
    vario45 = pkern_vario(dims45, vec45, sep=sep, simple=simple, method=method, diagonal=FALSE)

    # scale sep values to get distances, then bundle everything into output list
    dsdiag = sqrt( sum( ds^2 ) )
    vario45[['x']][['sep']] = dsdiag * vario45[['x']][['sep']]
    vario45[['y']][['sep']] = dsdiag * vario45[['y']][['sep']]
    list.out = c(list.out, list(d1=vario45[['y']], d2=vario45[['x']]))
  }

  # omit entries with NA semivariances
  list.out = lapply(list.out, \(d) lapply(d, \(z) z[ !is.na(d[['sv']]) ] ) )

  # truncate to requested distance maximum
  if( is.numeric(dmax) ) list.out = lapply(list.out, \(d) lapply(d, \(z) z[ !(d[['sep']] > dmax) ] ))

  # return results along with estimated variance
  return(c(list.out, list(v=v, ds=ds)))
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
      fitx = pkern_tvario(pars[['x']], v=sqrt(v), nug=nug, d=sepx)
      fity = pkern_tvario(pars[['y']], v=sqrt(v), nug=nug, d=sepy)
      fit = fitx * fity
    }

    # add line plot showing semivariance at requested lags
    lines(sep, fit)
  }

  return(invisible())
}



#' Fit a theoretical separable covariance model to a separated empirical semivariogram
#'
#' Pipe the results of a `pkern_vario` call to this function to estimate model
#' parameters for `xpars` and `ypars`, the x and y component kernels, by weighted
#' least squares (Cressie, 2015), where `stats::optim` ("L-BFGS-B") with default
#' settings is used to solve the minimization problem.
#'
#' `xpars` and `ypars` may be specified as character strings (the exponential
#' "exp" is the default), in which case they are assigned suggested initial values
#' and bounds. If parameter lists (recognized by `pkern_corr`), are supplied, the
#' function uses its "kp" entry for initial values, and if the list contains entries
#' "lower" and/or "upper", they are used as bounds.
#'
#' The pointwise variance is estimated simultaneously with the correlation kernel
#' parameters using the bounds and initial values in `v` and `nug`. If these are not
#' supplied, the function sets (for both) default default bounds of 0 to 4 times the
#' sample variance.
#'
#' The nugget effect can be disabled by setting its maximum to 0. Note that for
#' certain sampling configurations and kernels (particularly the Gaussian), this will
#' produce ill-conditioned covariance matrices, leading to inaccuracies or errors
#' in prediction functions like `pkern_cmean`.
#'
#' When `ninitial > 1`, the least squares optimization call is repeated for additional
#' sets of initial values, chosen uniformly at random within their upper and lower
#' bounds. Try increasing this to resolve issues of convergence to incorrect local optima.
#'
#' @param vario list, the return value from a call to `pkern_vario`
#' @param xpars character or list, recognizable (as `pars`) by `pkern_corr`
#' @param ypars character or list, recognizable (as `pars`) by `pkern_corr`
#' @param v length-3 numeric vector, c('min', 'ini', 'max'), bounds and initial value for variance
#' @param nug length-3 numeric vector, c('min', 'ini', 'max'), bounds and initial value for nugget
#' @param ninitial positive integer, the number of initial parameter sets tested (see details)
#'
#' @return
#' @export
#'
#' @examples
pkern_vario_fit = function(vario, xpars='mat', ypars=xpars, v=NULL, nug=NULL, ninitial=5)
{
  # handle kernel specification as string, convert to list
  if( is.character(xpars) ) xpars = pkern_corr(xpars)
  if( is.character(ypars) ) ypars = pkern_corr(ypars)

  # set defaults for variance and nugget
  if( is.null(v) ) v = c(min=0, ini=vario$v, max=4*vario$v)
  if( is.null(nug) ) nug = c(min=0, ini=v[2]/2, max=4*vario$v)

  # set default lower and upper bounds for kernel parameters as needed
  if( is.null(xpars$lower) ) xpars$lower = pkern_corr(xpars)[,'min']
  if( is.null(xpars$upper) ) xpars$upper = pkern_corr(xpars)[,'max']
  if( is.null(ypars$lower) ) ypars$lower = pkern_corr(ypars)[,'min']
  if( is.null(ypars$upper) ) ypars$upper = pkern_corr(ypars)[,'max']

  # vectorized bounds for all covariance parameters
  plower = c(v[1], nug[1], xpars$lower, ypars$lower)
  pinitial = c( v[2], nug[2], pkern_unpack(xpars, ypars) )
  pupper = c(v[3], nug[3], xpars$upper, ypars$upper)

  # build matrices to hold pertinent info from `vario` (omit zero point in first index)
  idx.vario = names(vario) %in% c('x', 'y', 'd1', 'd2')
  mvario = lapply(vario[idx.vario], \(vg) cbind(n=vg$n[-1], sep=vg$sep[-1], sv=vg$sv[-1]))

  # define anonymous objective function for optimizer
  fn = function(pvec, xp, yp, mv, ds)
  {
    # pvec, numeric vector of parameters
    # xp, x kernel parameter list (must element named "k")
    # yp, y kernel parameter list (must element named "k")
    # mv, list of matrices ("x", "y" and optionally "d1", "d2") with columns "n", "sep", "sv"
    # ds, numeric vector of grid line spacings in x and y directions

    # extract kernel parameters as lists (omit first element, the variance)
    pars = pkern_unpack(xp, yp, pvec[-(1:2)])

    # generate theoretical semivariance values along x and y for each lag
    xsv = pkern_tvario(pars=pars[['x']], v=pvec[1], nug=pvec[2], d=mv[['x']][,'sep'])
    ysv = pkern_tvario(pars=pars[['y']], v=pvec[1], nug=pvec[2], d=mv[['y']][,'sep'])

    # compute weighted sums of squares along both dimensions and return their sum
    #wss.x = sum( ( mv[['x']][,'n'] / (xsv)^2 ) * ( ( mv[['x']][,'sv'] - xsv )^2 ) )
    #wss.y = sum( ( mv[['y']][,'n'] / (ysv)^2 ) * ( ( mv[['y']][,'sv'] - ysv )^2 ) )
    wss.x = sum( ( mv[['x']][,'n'] ) * ( ( mv[['x']][,'sv'] - xsv )^2 ) )
    wss.y = sum( ( mv[['y']][,'n'] ) * ( ( mv[['y']][,'sv'] - ysv )^2 ) )

    # do the same for the 45 degree rotated version if it is supplied
    wss.d1 = wss.d2 = 0
    if( length(mv) == 4 )
    {
      # find unit distance along diagonals
      udist = sqrt(sum(ds^2))

      # generate theoretical semivariance values for each lag on first diagonal
      d1.xsep = ds[1] * mv[['d1']][,'sep'] / udist
      d1.ysep = ds[2] * mv[['d1']][,'sep'] / udist
      xsv.d1 = pkern_tvario(pars=pars[['x']], v=sqrt(pvec[1]), nug=pvec[2], d=d1.xsep)
      ysv.d1 = pkern_tvario(pars=pars[['y']], v=sqrt(pvec[1]), nug=pvec[2], d=d1.ysep)
      d1sv = xsv.d1 * ysv.d1

      # generate theoretical semivariance values for each lag on second diagonal
      d2.xsep = ds[1] * mv[['d2']][,'sep'] / udist
      d2.ysep = ds[2] * mv[['d2']][,'sep'] / udist
      xsv.d2 = pkern_tvario(pars=pars[['x']], v=sqrt(pvec[1]), nug=pvec[2], d=d2.xsep)
      ysv.d2 = pkern_tvario(pars=pars[['y']], v=sqrt(pvec[1]), nug=pvec[2], d=d2.ysep)
      d2sv = xsv.d2 * ysv.d2

      # compute weighted sums of squares along both dimensions and return their sum
      wss.d1 = sum( ( mv[['d1']][,'n'] ) * ( ( mv[['d1']][,'sv'] - d1sv )^2 ) )
      wss.d2 = sum( ( mv[['d2']][,'n'] ) * ( ( mv[['d2']][,'sv'] - d2sv )^2 ) )
      #wss.d1 = sum( ( mv[['d1']][,'n'] / (d1sv)^2 ) * ( ( mv[['d1']][,'sv'] - d1sv )^2 ) )
      #wss.d2 = sum( ( mv[['d2']][,'n'] / (d2sv)^2 ) * ( ( mv[['d2']][,'sv'] - d2sv )^2 ) )
    }

    # compute total of weighted sums of squares
    tss = wss.x + wss.y + wss.d1 + wss.d2
    if( is.na(tss) | is.infinite(tss) ) tss = 2^.Machine$double.digits
    return( tss )
  }

  # make a list of initial values to test
  np = length(pinitial)
  list.initial = c(list(pinitial), lapply(seq(ninitial), \(ini) plower + runif(np)*(pupper-plower)) )

  # run the optimizer for each one
  list.optim = lapply(list.initial, \(ini) stats::optim(par=ini,
                                                        f=fn,
                                                        method='L-BFGS-B',
                                                        lower=plower,
                                                        upper=pupper,
                                                        xp=xpars,
                                                        yp=ypars,
                                                        mv=mvario,
                                                        ds=vario[['ds']]))

  # select the best fit
  idx.best = which.min( sapply(list.optim, \(result) result$value) )
  result.optim = list.optim[[idx.best]]

  # unpack fitted parameters and return in a list
  vfit = setNames(result.optim$par[1], 'fitted')
  nugfit = setNames(result.optim$par[2], 'fitted')
  pfit = pkern_unpack(xpars, ypars, result.optim$par[-(1:2)])
  return( list(x=pfit[['x']], y=pfit[['y']], v=vfit, nug=nugfit, ds=vario[['ds']]) )
}
