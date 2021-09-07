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
#' assuming pointwise variance `v` and returned as a vector the same length as `d`.
#'
#' When `d` is NULL, the function instead returns the semivariance function as an
#' anonymous function of distance.
#'
#' @param pars list of kernel parameters, in form recognized by `pkern_corr`
#' @param v numeric, the pointwise variance
#' @param d vector of nonegative numeric spatial lags to evaluate
#'
#' @return either a vector the same length as `d`, or an anonymous function of distance
#' @export
#'
#' @examples
#' pars.mat = pkern_corr('mat')
#' pkern_tvario(pars.mat, v=1, 1:10)
#' sv = pkern_tvario(pars.mat, v=1)
#' sv(1:10)
pkern_tvario = function(pars, v, d=NULL)
{
  fn = function(d) { unname(v) * ( 1 - pkern_corr(pars, d) ) }
  if( is.null(d) ) return(fn)
  return(fn(d))
}



#' Sample semivariogram for regular 2D gridded data along x direction
#'
#' The function identifies and computes the semivariance (in `vec`) of point pairs on a
#' regular grid of dimensions `dims` at the supplied x grid line lags, `sep`. Elements of
#' `set` must be a subset of `1:(dims[1]-1)` (the default setting).
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
#' @return a list with named elements "sep", "vario", and "pairs"
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
#' A wrapper for two calls to `pkern_xvario` to get the directional sample semivariograms
#' along the x and y dimensions, on regular lattice data. Note that `sep` is given in terms
#' of grid line indices (integers), not distances (their product with the x or y resolution)
#'
#' The function returns its two (x and y component) results in a list along with
#' "v" the estimated (univariate) variance of `vec`
#'
#' @param vec numeric vector of data, in column vectorized order (length `prod(dims)`)
#' @param dims vector c(nx, ny) of positive integers, the number of grid lines
#' @param sep list c(dx, dy) of positive integer vectors, grid line lags to sample
#' @param nmin positive integer, the minimum number of point pairs needed to include an estimate
#' @param simple logical, if FALSE the function includes point pair indices in return list
#' @param method character, one of "mean", "median", "ch", passed to `pkern_xvario`
#'
#' @return list with elements "x" and "y" (returned from `pkern_xvario`) and "sig"
#' @export
#'
#' @examples
pkern_vario = function(dims, vec, sep=NA, nmin=2, simple=TRUE, method='median')
{
  # set default sample lags
  if( anyNA(sep) ) sep = lapply(dims, \(d) seq(d-1))

  # for y separations, switch to row vectorized order and relabel dimensions
  xout = pkern_xvario(dims, vec, simple=simple, method=method)
  yout = pkern_xvario(rev(dims), vec[pkern_r2c(dims, FALSE, TRUE)], simple=simple, method=method)

  # index of lags with more points than the threshold
  idx.x = xout$n > nmin - 1
  idx.y = yout$n > nmin - 1

  # return results in a list
  return(list(x=xout, y=yout, v=var(vec, na.rm=TRUE)))
}



#' Plot 1-dimensional empirical semivariograms for a separable model
#'
#' @param vario list of semivariance data, the return value of `pkern_vario` or `pkern_xvario`
#' @param pars list of kernel parameters including variance in named element "v"
#' @param nmin positive integer, the minimum number of point pairs needed to include in plot
#' @param ptitle character or vector of two, title(s) for the plot(s)
#' @param shade logical, indicating to shade points according to the number of samples
#'
#' @return prints a plot and returns nothing
#' @export
#'
#' @examples
pkern_vario_plot = function(vario, pars=NULL, nmin=1, ptitle=NULL, shade=TRUE)
{
  # handle 2-dimensional case
  if( all( c('x', 'y') %in% names(vario) ) )
  {
    # assign default plot titles as needed
    if( length(ptitle) != 2 ) ptitle = c(x='x', y='y')

    # assign names to pars if they are missing
    nm.pars = c('x', 'y', 'v')
    if( !is.null(pars) & !all( nm.pars %in% names(pars) ) ) names(pars) = nm.pars

    # make the plots
    par(mfrow=c(1,2))
    pkern_vario_plot(vario[['x']], list(x=pars[['x']], v=pars[['v']]), nmin, ptitle['x'], shade)
    pkern_vario_plot(vario[['y']], list(x=pars[['y']], v=pars[['v']]), nmin, ptitle['y'], shade)
    par(mfrow=c(1,1))
    return(invisible())
  }

  # handle 1-dimension case

  # unpack arguments
  idx = vario[['n']] > nmin - 1
  n = vario[['n']][idx]
  sep = vario[['sep']][idx]
  sv = vario[['sv']][idx]

  # map point shading to the number of point pairs sampled if requested
  colmap = 'black'
  if( shade ) colmap = rev( gray.colors( max(n) ) )[n]

  # make the scatterplot of sampled semivariances
  plot(sep, sv, xlab='lag', ylab='semivariance', main=ptitle, pch=16, col=colmap)

  # if kernel parameters are supplied, plot the theoretical values
  if( !is.null(unlist(pars)) )
  {
    # check for variance in parameters list
    idxv = names(pars) == 'v'
    if( !any(idxv) ) stop('variance v not found in pars')

    # extract kernel parameters and variance
    nm.pars = names(pars)
    v = pars[[ which(nm.pars %in% 'v')[1] ]]
    pars = pars[[ which(nm.pars != 'v')[1] ]]

    # evaluate semivariance at requested lags and add to plot as line
    fit = pkern_tvario(pars, v, d=sep)
    lines(sep, fit)
  }

  return(invisible())
}



#' Fit a theoretical separable covariance model to a separated empirical semivariogram
#'
#' Pipe the results of a `pkern_vario` call to this function to estimate model
#' parameters for `xpars` and `ypars`, the x and y component kernels, by weighted
#' least squares (Cressie, 2015), where `stats::optim` ("L-BFGS-B" method with default
#' settings is used to solve the minimization problem.
#'
#' `xpars` and `ypars` may be specified as character strings (the exponential
#' "exp" is the default), in which case they are assigned suggested initial values
#' and bounds. If parameter lists (recognized by `pkern_corr`), are supplied, the
#' function uses its "kp" entry for initial values, and if the list contains entries
#' "lower" and/or "upper", they are used as bounds.
#'
#' The variance is estimated simultaneously with the correlation kernel parameters
#' using bounds and initial value in `v`. If these are not supplied, the function
#' sets reasonable defaults.
#'
#' @param vario list, the return value from a call to `pkern_vario`
#' @param xpars character or list, recognizable (as `pars`) by `pkern_corr`
#' @param ypars character or list, recognizable (as `pars`) by `pkern_corr`
#' @param v length-3 numeric vector, c('min', 'ini', 'max'), bounds and initial value for variance
#' @param nmin positive integer, the minimum number of samples to include a lag
#'
#' @return
#' @export
#'
#' @examples
pkern_vario_fit = function(vario, xpars='exp', ypars=xpars, v=NA, nmin=2)
{
  # handle kernel specification as string, convert to list
  if( is.character(xpars) ) xpars = pkern_corr(xpars)
  if( is.character(ypars) ) ypars = pkern_corr(ypars)

  # set defaults for sigma
  if( anyNA(v) ) v = c(min=0, ini=1, max=2*vario$v)

  # vectorized initial values for covariance parameters
  pinitial = c( v[2], pkern_unpack(xpars, ypars) )

  # set default lower and upper bounds for kernel parameters as needed
  if( is.null(xpars$lower) ) xpars$lower = pkern_corr(xpars)[,'min']
  if( is.null(xpars$upper) ) xpars$upper = pkern_corr(xpars)[,'max']
  if( is.null(ypars$lower) ) ypars$lower = pkern_corr(ypars)[,'min']
  if( is.null(ypars$upper) ) ypars$upper = pkern_corr(ypars)[,'max']

  # vectorized bounds for all covariance parameters
  plower = c(v[1], xpars$lower, ypars$lower)
  pupper = c(v[3], xpars$upper, ypars$upper)

  # check that we have enough data in the empirical variograms
  is.ok = lapply(vario[c('x', 'y')], \(v) v$n > nmin - 1)
  if( any( sapply(is.ok, sum) < 2 ) ) stop('not enough lags sampled. Check nmin and vario')

  # build two matrices to hold pertinent info from `vario`
  mvario = mapply(\(v, idx) cbind(n=v$n[idx], sep=v$sep[idx], sv=v$sv[idx]),
                  v=vario[c('x', 'y')], idx=is.ok)

  # define anonymous function to compute scores for optimizer
  fn = function(pvec, xpars, ypars, mvario)
  {
    # extract kernel parameters as lists (omit first element, the variance)
    pars = pkern_unpack(xpars, ypars, pvec[-1])

    # generate theoretical semivariance values along x and y for each lag
    xsv = pkern_tvario(pars=pars[['xpars']], v=pvec[1], d=mvario[['x']][,'sep'])
    ysv = pkern_tvario(pars=pars[['ypars']], v=pvec[1], d=mvario[['y']][,'sep'])

    # compute weighted sums of squares along both dimensions and return their sum
    wss.x = sum( mvario[['x']][,'n'] * ( ( mvario[['x']][,'sv'] - xsv )^2 ) )
    wss.y = sum( mvario[['y']][,'n'] * ( ( mvario[['y']][,'sv'] - ysv )^2 ) )
    return( wss.x + wss.y )
  }

  # run optimizer
  result.optim = stats::optim(par=pinitial, f=fn, method='L-BFGS-B',
                              lower=plower, upper=pupper,
                              xpars=xpars, ypars=ypars, mvario=mvario)

  # unpack fitted parameters
  vfit = setNames(result.optim$par[1], 'fitted')
  pfit = pkern_unpack(xpars, ypars, result.optim$par[-1])

  # compute fitted values
  xfit = pkern_tvario(pars=pfit[['xpars']], v=vfit, d=mvario[['x']][,'sep'])
  yfit = pkern_tvario(pars=pfit[['ypars']], v=vfit, d=mvario[['y']][,'sep'])

  return( list(x=pfit[['xpars']], y=pfit[['ypars']], v=vfit) )
}
