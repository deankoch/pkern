#
# pkern_raster.R
# Dean Koch, Oct 2021
# Functions for making and fitting semivariograms
#


#' Theoretical semivariances along x or y direction
#'
#' Convenience function for getting theoretical semivariance values for a set of
#' spatial lags along one dimension, assuming (intrinsic) stationarity.
#'
#' Kernel parameters `pars` are passed to the correlation function (`pkern_corr`)
#' to get correlation values at distances `d`. The result is converted to semivariance
#' assuming a pointwise variance of `v + nug` (an optional nugget effect) and returned
#' as a vector the same length as `d`.
#'
#' When `d` is NULL, the function instead returns the semivariance function as an
#' anonymous function of distance.
#'
#' @param pars list of kernel parameters ("k" and "kp"), in form recognized by `pkern_corr`
#' @param v positive numeric, the pointwise variance absent a nugget effect
#' @param nug nonegative numeric, the variance of the nugget effect
#' @param d vector of nonegative numeric spatial lags to evaluate
#'
#' @return either a vector the same length as `d`, or an anonymous function of distance
#' @export
#'
#' @examples
#' d = seq(1e2)/10
#' pars = pkern_corr('mat')
#' pkern_tvario(pars, v=1, d=d)
#' pkern_tvario(pars, v=2, nug=0.5, d=d)
#'
#' sv = pkern_tvario(pars, v=2, nug=0.5)
#' plot(d, sv(d))
pkern_tvario = function(pars, v, nug=0, d=NULL)
{
  fn = function(d) { nug + ( unname(v) * ( 1 - pkern_corr(pars, d) ) ) }
  if( is.null(d) ) return(fn)
  return(fn(d))
}


#' Sample semivariogram for regular 2D gridded data along x direction
#'
#' The function computes the semivariance of point pairs located on the regular grid with
#' dimensions `gdim`, and data vector(s) `vec`, at the supplied x grid line lags, `lags`.
#' Elements of `lags` are must be a subset of `seq(gdim[2]-1)`.
#'
#' Results are returned in a named list containing the sampled lags (`lags`), the number of
#' point pairs sampled at each lag (`n`) and the sample semivariances (`sv`) at each lag;
#' If `simple=FALSE`, the function also returns `pairs`, a list of matrices containing the
#' indices of point pairs sampled and their absolute differences.
#'
#' `vec` can be a vector (in column vectorized order), or a list of such vectors for
#' repeated measurements at the same spatial locations. NAs are permitted but point
#' pairs containing NAs are excluded from computations.
#'
#' `nmax` specifies a maximum number of point pairs to process at each lag. When applicable,
#' a random sample of (`nmax`) eligible point pairs are taken at each lag. When `vec` is
#' a list of vectors, the function always samples at least one point pair from each of them.
#' Note that this means `nmax` is not respected when `vec` is a list and `nmax < length(vec)`.
#' If the function is too slow in this case, try taking a subset of `vec`.
#'
#' method "mean" is the classical method of Matheron (1962), "median" and "ch" are the
#' median and robust methods described in Cressie and Hawkins (1980).
#'
#' references:
#'
#' Cressie and Hawkins (1980): https://doi.org/10.1007/BF01035243
#'
#'
#' @param gdim vector c(ny, nx) of positive integers, the number of grid lines in full grid
#' @param vec numeric vector of data, in column vectorized order (length `prod(gdim)`)
#' @param lags positive integer vector, the grid line lags to sample
#' @param simple logical, if FALSE the function returns a list of point pair indices
#' @param method character, one of "mean", "median", "ch"
#' @param quiet logical, suppresses console output
#' @param nmax positive integer, the maximum number of point pairs to sample at each lag
#'
#' @return a list with named elements "lags", "n", "sv", and optionally "pairs"
#' @export
#'
#' @examples
#'
#' gdim = c(1e2, 50)
#' pars = pkern_corr('gau') |> utils::modifyList(list(v=1))
#' sim = pkern_sim(pars, gdim)
#'
#' # sample first 10 lags and plot
#' vario = pkern_xvario(gdim, sim, lags=seq(10))
#' pkern_vario_plot(vario)
#'
#' # sample all lags (the default)
#' vario = pkern_xvario(gdim, sim)
#' pkern_vario_plot(vario, pars)
#'
pkern_xvario = function(gdim, vec, simple=TRUE, method='median', quiet=FALSE, lags=NULL, nmax=1e4)
{
  # set default sample lags and express vector `vec` as length-1 list
  if( is.null(lags) ) lags = seq( gdim[2] - 1 )
  if( !is.list(vec) ) vec = list(vec)

  # prepare to sample along different grid line lags in a loop
  lags = sort(lags)
  n.lags = length(lags)
  plist = vector(mode='list', length=n.lags)
  svx = rep(NA, n.lags)
  n = rep(0, n.lags)
  miss = lapply(vec, is.na)

  # console information
  if( !quiet )
  {
    cat(paste('\nsampling semivariance at', n.lags, 'lags...\n'))
    pb = utils::txtProgressBar(max=n.lags, style=3)
  }

  # loop over separation distances
  for( idx.dj in seq_along(lags) )
  {
    if( !quiet ) utils::setTxtProgressBar(pb, idx.dj)
    dj = lags[idx.dj]

    # total number of eligible point pairs
    idx.max = gdim[1] * ( gdim[2] - dj )

    # vectorization trick to convert this to samples of i,j indices
    ij.sample = pkern_vec2mat(seq(idx.max), gdim, out='list')

    # translate sample points index into index of vec
    idx1 = ij.sample[['i']] + ( gdim[1] * (ij.sample[['j']] - 1) )
    idx2 = idx1 + ( gdim[1] * dj )

    # omit point pairs containing NAs, sample uniformly to satisfy nmax
    idx.valid = lapply(miss, \(m) !( m[idx1] | m[idx2] ) )
    n[idx.dj] = do.call(sum, idx.valid)
    if( n[idx.dj] > nmax )
    {
      # calculate number of point pairs to keep from each list element in vec
      n.per = max(1, round(  ( nmax / n[idx.dj] ) * sapply(idx.valid, sum) ) )

      # randomly sample n.per from each element of idx.valid to get nmax samples (or fewer)
      idx.valid = lapply(idx.valid, \(i) seq_along(i) %in% sample( which(i), min(n.per, sum(i)) ) )
      n[idx.dj] = do.call(sum, idx.valid)
    }

    # store the point pair indices
    plist[[idx.dj]] = lapply(idx.valid, \(valid) cbind(idx1=idx1[valid], idx1=idx2[valid]) )

    # estimate semivariance
    if(n[idx.dj] > 0)
    {
      # compute absolute differences of point pairs and append to plist
      dabs = Map(\(p, v) apply(p, 1, \(idx) abs(diff(v[idx])) ), p=plist[[idx.dj]], v=vec)
      if(!simple) plist[[idx.dj]] = Map(\(p, d) cbind(p, d=d), p=plist[[idx.dj]], d=dabs)

      # handle various methods
      dabs.vec = unlist(dabs, use.names=FALSE)
      if(method=='mean') svx[idx.dj] = ( mean(dabs.vec)^2 ) / 2
      if(method=='median') svx[idx.dj] = stats::median( sqrt(dabs.vec) )^4 / (2 * 0.457)
      if(method=='ch')
      {
        ch.factor = 0.457 + ( 0.494 / n[idx.dj] ) + ( 0.494 / ( n[idx.dj]^2 ) )
        svx[idx.dj] = mean( sqrt(dabs.vec) )^4 / (2 * ch.factor)
      }
    }
  }
  if( !quiet ) close(pb)

  if( simple ) return( list(lags=lags, n=n, sv=svx) )
  return( list(lags=lags, n=n, sv=svx, pairs=plist) )
}


#' Directional sample semivariograms for regular 2D gridded data
#'
#' A wrapper for calls to `pkern_xvario` to get the directional sample semivariograms
#' from regular lattice data along the y and x dimensions, and optionally along the
#' diagonals.
#'
#' The function makes the appropriate calls to `pkern_xvario` and returns the results
#' in list elements "y", "x". If `diagonal=TRUE`, the function also returns "d1" and
#' "d2", the y and x semivariograms after a 45 degree counterclockwise rotation of the
#' input grid (see `pkern_r45`) assuming equal y and x resolution. For unequal resolution
#' the rotation angle is `atan(gres[2]/gres[1])`.
#'
#' Note that input `lags` is specified in terms of grid line indices (integers), not
#' distances (their product with the x or y resolution). However, non-integer separation
#' distances can be specified by setting the resolution of the input grid (dy, dx) in
#' argument `gres`. When `gres` is supplied, the output "lags" values are scaled accordingly
#' and `dmax` is interpreted as a distance, rather than a number of grid lines.
#'
#' When `gres` is not supplied, the function assumes that the x grid lines have the same
#' spacing distance (1) as the y grid lines. In this case, diagonally adjacent
#' grid points have separation distance `sqrt(2)`.
#'
#' @param gdim vector c(ny, nx) of positive integers, the number of grid lines
#' @param vec numeric vector of data, in column vectorized order (length `prod(dims)`)
#' @param lags list c(dy, dx) of positive integer vectors, grid line lags to sample
#' @param method character, one of "mean", "median", "ch", passed to `pkern_xvario`
#' @param diagonal logical, if TRUE, semivariogram results are returned also for diagonals
#' @param dmax numeric, maximum grid line separation distance (used to filter results)
#' @param gres length-2 vector of positive numeric, the spacing distance of y and x grid lines
#' @param quiet logical, indicating to suppress console output
#' @param nmax positive integer, the maximum number of point pairs to process in each sample
#' @param vcalc logical, indicates to compute the sample variance (which can be relatively slow)
#'
#' @return list of `pkern_xvario` output "x", "y" (and "d1", "d2"); sample variance "v"; and "gres"
#' @export
#'
#' @examples
#' gdim = c(1e2, 2e2)
#' pars = pkern_corr('gau') |> utils::modifyList(list(v=1))
#' vec = pkern_sim(pars, gdim)
#'
#' # sample first 10 lags and plot
#' vario = pkern_vario(gdim, vec, lags=seq(10))
#' pkern_vario_plot(vario, pars)
#'
#' # example with unequal resolution and multiple replicates
#' gres = c(2,1)
#' nrep = 10
#' vec = lapply(seq(nrep), \(i) pkern_sim(pars, gdim, gres=gres, makeplot=FALSE))
#' vario = pkern_vario(gdim, vec, lags=seq(10), gres=gres)
#' pkern_vario_plot(vario, pars)
#'
pkern_vario = function(gdim, vec, lags=NA, method='median', diagonal=TRUE,
                       gres=NA, quiet=FALSE, nmax=1e4, vcalc=TRUE)
{
  # handle `gdim` as output from `pkern_snap`
  if( is.list(gdim) )
  {
    # if subgrid info is provided, it is assumed `vec` is a data vector from the subgrid
    if( !is.null(gdim[['sg']]) ) gdim = gdim[['sg']]
    if( is.na(gres) & !is.null(gdim[['gres']]) ) gres = gdim[['gres']]
    gdim = gdim[['gdim']]
  }

  # set default sample lags and compute sample variance
  if( anyNA(lags) ) lags = lapply(gdim, \(d) seq(d-1))
  if( !is.list(lags) ) lags = list(y=lags, x=lags)
  if( !is.list(vec) ) vec = list(vec)
  svar = ifelse(vcalc, stats::var(unlist(vec), na.rm=TRUE), NA)

  # set default resolution (1,1)
  yxnm = c('y', 'x')
  if( anyNA(gres) ) gres = c(1,1)
  if( length(gres) == 1 ) stop('gres must have length 2')
  if( !all( yxnm %in% names(gres) ) ) gres = setNames(gres, yxnm)

  # compute x semivariances
  if( !quiet ) cat('\nidentifying horizontal lags...')
  xout = pkern_xvario(gdim, vec, method=method, quiet=quiet, lags=lags[[2]], nmax=nmax)

  # repeat in row-vectorized order to get y semivariances
  if( !quiet ) cat('\nidentifying vertical lags...')
  vec.rv = lapply(vec, \(v) v[pkern_r2c(gdim, in.byrow=FALSE, out.byrow=TRUE)])
  yout = pkern_xvario(rev(gdim), vec.rv, method=method, quiet=quiet, lags=lags[[1]], nmax=nmax)

  # scale lags by resolution to get distances, then bundle everything into an output list
  yout[['lags']] = gres['y'] * yout[['lags']]
  xout[['lags']] = gres['x'] * xout[['lags']]
  list.out = list(y=yout, x=xout)

  # add a zero-distance point so that the output sets are never empty
  list.out = lapply(list.out, \(d) lapply(d, \(z) c(0, z) ))

  # repeat with rotated coordinates if requested
  if( diagonal )
  {
    # rotate the input array data by 45 degrees clockwise
    if( !quiet ) cat('\nrotating coordinate system by 45 degrees:\n')
    vec45 = lapply(vec, \(v) c(pkern_r45(v, gdim)))

    # find the dimensions rotated grid (based on first list entry of vec)
    gdim45 = dim( pkern_r45(vec[[1]], gdim) )

    # recursive call to get component variograms in new coordinate system
    vario45 = pkern_vario(gdim45, vec45, lags, method=method, quiet=quiet, diagonal=F, nmax=nmax, vcalc=FALSE)

    # scale lags values to get distances, then bundle everything into output list
    gresdiag = sqrt( sum( gres^2 ) ) / 2
    vario45[['y']][['lags']] = gresdiag * vario45[['y']][['lags']]
    vario45[['x']][['lags']] = gresdiag * vario45[['x']][['lags']]
    list.out = c(list.out, list(d1=vario45[['y']], d2=vario45[['x']]))
  }

  # omit entries with NA semivariances
  list.out = lapply(list.out, \(d) lapply(d, \(z) z[ !is.na(d[['sv']]) ] ) )

  # truncate to requested distance maximum
  #if( is.numeric(dmax) ) list.out = lapply(list.out, \(d) lapply(d, \(z) z[ !(d[['lags']] > dmax) ] ))

  # return results along with estimated variance
  return(c(list.out, list(v=svar, gres=gres)))
}


#' Fit a theoretical separable covariance model to a separated empirical semivariogram
#'
#' Pipe the results of a `pkern_vario` call to this function to estimate model
#' parameters for `xpars` and `ypars`, the x and y component kernels, by weighted
#' least squares using `stats::optim` ("L-BFGS-B") for numerical minimization.
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
#' When `add > 0`, the least squares optimization call is repeated for additional
#' sets of initial values, chosen uniformly at random within their upper and lower
#' bounds. The best result of these is selected an returned. Try increasing `add`
#' to resolve issues of convergence to incorrect local optima.
#'
#' `fit.method` controls the type of weights assigned to the sum of squares from each
#' spatial lag: '1' assigns unit weight to all lags (OLS); 'n' weights be sample size;
#' 'n/d' uses sample size divided by squared distance (the default in `gstat`); and
#' the default, 'n/v', uses sample size divided by squared theoretical semivariance,
#' (the robust estimator in Chapter 2.6.2 of Cressie, 1993).
#'
#' @param vario list, the return value from a call to `pkern_vario`
#' @param ypars character or list, recognizable (as `pars`) by `pkern_corr`
#' @param xpars character or list, recognizable (as `pars`) by `pkern_corr`
#' @param v length-3 numeric vector, c('min', 'ini', 'max'), bounds and initial value for variance
#' @param nug length-3 numeric vector, c('min', 'ini', 'max'), bounds and initial value for nugget
#' @param add positive integer, the number of initial parameter sets tested (see details)
#' @param dmax positive numeric, point pairs with separation distances exceeding `dmax` are ignored
#' @param fit.method character string, either '1', 'n', 'n/d', or 'n/v' (the default, see details)
#'
#' @return TODO
#' @export
#'
#' @examples
#' gdim = c(1e2, 2e2)
#' pars = pkern_corr('gau') |> utils::modifyList(list(v=1))
#' vec = pkern_sim(pars, gdim)
#'
#' # sample first 10 lags and plot
#' vario = pkern_vario(gdim, vec, lags=seq(10))
#' pkern_vario_plot(vario, pars)
#'
#' # fit the separable gaussian kernel
#' pars.fitted = pkern_vario_fit(vario)
#' pkern_vario_plot(vario, pars.fitted)
#' pkern_plot(pars.fitted)
#'
#' # fit a different separable kernel
#' pars.fitted2 = pkern_vario_fit(vario, ypars='sph', xpars='mat')
#' pkern_vario_plot(vario, pars.fitted2)
#' pkern_plot(pars.fitted2)
#'
pkern_vario_fit = function(vario, ypars='gau', xpars=ypars, v=NULL, nug=NULL, add=0, dmax=Inf,
                           fit.method='n/v')
{
  # set defaults for variance (first parameter) and nugget (second parameter)
  if( is.null(v) ) v = c(min=1e-9, ini=vario[['v']], max=4*vario[['v']])
  if( is.null(nug) ) nug = c(min=1e-3, ini=v[2]/2, max=4*vario[['v']])

  # set up defaults and bounds for kernel parameters (2-4 additional parameters)
  ypars = pkern_bds(ypars, vario[['gres']]['y'])
  xpars = pkern_bds(xpars, vario[['gres']]['x'])

  # vectorized bounds for all covariance parameters
  plower = c(v[1], nug[1], xpars[['lower']], ypars[['lower']])
  pinitial = c( v[2], nug[2], xpars[['kp']], ypars[['kp']] )
  pupper = c(v[3], nug[3], xpars[['upper']], ypars[['upper']])

  # build matrices to hold pertinent info from `vario` (omit zero point in first index)
  midx = names(vario) %in% c('x', 'y', 'd1', 'd2')
  mvario = lapply(vario[midx], \(vg) cbind(n=vg[['n']][-1], lags=vg[['lags']][-1], sv=vg[['sv']][-1]))

  # truncate at dmax
  mvario = Map(\(m, i) m[i,], mvario, sapply(mvario, \(m) m[, 'lags'] < dmax))
  msg.err = 'each direction must contain at least one lag. Try increasing dmax'
  if( any( sapply(mvario, length) == 0 ) ) stop(msg.err)

  # define anonymous objective function for optimizer
  fn = function(pvec, p, mv, gres, fit.method)
  {
    # pvec, numeric vector of parameters
    # p, list of "x", "y" kernel parameter lists associated with pvec
    # mv, list of matrices ("x", "y" and optionally "d1", "d2") with columns "n", "lags", "sv"
    # gres, numeric vector of grid line spacings in x and y directions

    # extract kernel parameters as lists (omit first elements, the variance components)
    p = pkern_unpack(p, pvec[-(1:2)])

    # generate theoretical semivariance values along x and y for each lag
    ysv = pkern_tvario(pars=p[['y']], v=pvec[1], nug=pvec[2], d=mv[['y']][,'lags'])
    xsv = pkern_tvario(pars=p[['x']], v=pvec[1], nug=pvec[2], d=mv[['x']][,'lags'])

    # set up weights for y direction
    wy = switch(fit.method,
                '1' = 1,
                'n' = mv[['y']][,'n'],
                'n/v'= mv[['y']][,'n'] / (ysv)^2,
                'n/d'=( mv[['y']][,'n'] / (mv[['y']][,'lags'])^2 ))

    # set up weights for x direction
    wx = switch(fit.method,
                '1' = 1,
                'n' = mv[['x']][,'n'],
                'n/v'= mv[['x']][,'n'] / (xsv)^2 ,
                'n/d'=( mv[['x']][,'n'] / (mv[['x']][,'lags'])^2 ))

    # compute weighted sums of squares along both dimensions and return their sum
    wss.y = sum( wy * ( ( mv[['y']][,'sv'] - ysv )^2 ) )
    wss.x = sum( wx * ( ( mv[['x']][,'sv'] - xsv )^2 ) )

    # do the same for the 45 degree rotated version if it is supplied
    wss.d1 = wss.d2 = 0
    if( length(mv) == 4 )
    {
      # find unit distance along diagonals
      udist = sqrt(sum(gres^2))

      # generate theoretical semivariance values for each lag on first diagonal
      d1.ylags = gres[1] * mv[['d1']][,'lags'] / udist
      d1.xlags = gres[2] * mv[['d1']][,'lags'] / udist
      ysv.d1 = pkern_tvario(pars=p[['y']], v=sqrt(pvec[1]), d=d1.ylags)
      xsv.d1 = pkern_tvario(pars=p[['x']], v=sqrt(pvec[1]), d=d1.xlags)

      # nugget gets added after the product
      d1sv = pvec[2] + ( xsv.d1 * ysv.d1 )

      # repeat for the other diagonal
      d2.ylags = gres[1] * mv[['d2']][,'lags'] / udist
      d2.xlags = gres[2] * mv[['d2']][,'lags'] / udist
      ysv.d2 = pkern_tvario(pars=p[['y']], v=sqrt(pvec[1]), d=d2.ylags)
      xsv.d2 = pkern_tvario(pars=p[['x']], v=sqrt(pvec[1]), d=d2.xlags)
      d2sv = pvec[2] + ( xsv.d2 * ysv.d2 )

      # set up weights for d1 direction
      wd1 = switch(fit.method,
                   '1' = 1,
                   'n' = mv[['d1']][,'n'],
                   'n/v'= mv[['d1']][,'n'] / (d1sv)^2,
                   'n/d'=( mv[['d1']][,'n'] / (mv[['d1']][,'lags'])^2 ))


      # set up weights for d2 direction
      wd2 = switch(fit.method,
                   '1' = 1,
                   'n' = mv[['d2']][,'n'],
                   'n/v'= mv[['d2']][,'n'] / (d2sv)^2 ,
                   'n/d'=( mv[['d2']][,'n'] / (mv[['d2']][,'lags'])^2 ))

      # compute weighted sums of squares along both dimensions and return their sum
      wss.d1 = sum( wd1 * ( ( mv[['d1']][,'sv'] - d1sv )^2 ) )
      wss.d2 = sum( wd2 * ( ( mv[['d2']][,'sv'] - d2sv )^2 ) )
    }

    # compute total of weighted sums of squares
    tss = wss.x + wss.y + wss.d1 + wss.d2
    if( is.na(tss) | is.infinite(tss) ) tss = 2^.Machine$double.digits
    return( tss )
  }

  # make a list of initial values to test
  np = length(pinitial)
  list.ini = list(pinitial)
  if(add>0) list.ini = c(list.ini, lapply(seq(add), \(ini) plower+stats::runif(np)*(pupper-plower)) )

  # run the optimizer for each one
  list.optim = lapply(list.ini, \(ini) stats::optim(par=ini,
                                                        f=fn,
                                                        method='L-BFGS-B',
                                                        lower=plower,
                                                        upper=pupper,
                                                        p=list(x=xpars, y=ypars),
                                                        mv=mvario,
                                                        gres=vario[['gres']],
                                                        fit.method=fit.method))

  # select the best fit
  idx.best = which.min( sapply(list.optim, \(result) result$value) )
  result.optim = list.optim[[idx.best]]

  # unpack fitted parameters and return in a list
  vfit = stats::setNames(result.optim$par[1], 'fitted')
  nugfit = stats::setNames(result.optim$par[2], 'fitted')
  pfit = pkern_unpack(list(y=ypars, x=xpars), result.optim$par[-(1:2)])
  return( list(y=pfit[['y']], x=pfit[['x']], v=vfit, nug=nugfit, gres=vario[['gres']]) )
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
#' vario = pkern_vario(gdim, vec, lags=seq(10))
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
      lagsmax = sapply( vario[nm.yx], \(xy) max( xy[['lags']] ) )
      plotpars[['dmax']] = max(lagsmax)
    }

    # same for y-axis
    svlim.def = FALSE
    if( is.null(plotpars[['svlim']]) )
    {
      svlim.def = TRUE
      idx.valid = lapply( vario[nm.yx], \(xy) xy[['lags']] < plotpars[['dmax']] )
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
        lagsmax = sapply( vario[nm.dyx], \(xy) max( xy[['lags']] ) )
        if( dmax.def ) plotpars[['dmax']] = max(plotpars[['dmax']], lagsmax)

        # y axis
        idx.valid = lapply( vario[nm.vario], \(xy) xy[['lags']] < plotpars[['dmax']] )
        svmax = mapply(\(xy, i) max( xy[['sv']][i] ), xy=vario[nm.vario], i=idx.valid)
        plotpars[['svlim']] = range(plotpars[['svlim']], svmax)
      }

      # assign plot titles as needed
      if( title.def )
      {
        dstr = paste0('(', round( atan( vario[['gres']][1] / vario[['gres']][2] )*180/pi, 1), ' degrees)')
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
      gres = vario[['gres']]
      pkern_vario_plot(c(vario[['d1']], list(gres=gres)), pars=pars, plotpars=plotpars.list[['d1']])
      pkern_vario_plot(c(vario[['d2']], list(gres=gres)), pars=pars, plotpars=plotpars.list[['d2']])
    }

    # reset plot panel layout before finishing, returning nothing
    graphics::par(mfrow=c(1,1))
    return(invisible())
  }

  # 1-dimension case:

  # unpack plotting limits
  if( is.null(plotpars[['dmax']]) ) plotpars[['dmax']] = max(vario[['lags']])
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
    lags = vario[['lags']][-1]
    sv = vario[['sv']][-1]

    # map point shading to the number of point pairs sampled if requested
    #colmap = 'black'
    #if( plotpars[['shade']] ) colmap = rev( grDevices::gray.colors( max(n) ) )[n]

    # make the scatterplot of sampled semivariances
    plot(lags, sv, xlab='lag', ylab='semivariance', xlim=xlim, ylim=ylim, main=main, pch=16)

    # add bubbles indicating sample sizes
    points(lags, sv, pch=1, cex=n/max(n))

    # # make the scatterplot of sampled semivariances
    # plot(lags, sv, xlab='lag', ylab='semivariance',
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
    lags = seq(xlim[1], xlim[2], length.out=1e2)
    if( all( c('k', 'kp') %in% names(pars) ) ) fit = pkern_tvario(pars, v=v, nug=nug, d=lags)

    # handle requests along a diagonal (both y and x components included in pars)
    if( all( nm.yx %in% names(pars) ) )
    {
      # find unit distance along diagonals
      gres = vario[['gres']]
      udist = sqrt(sum(gres^2))

      # find the equivalent displacement along y and x directions
      lagsy = gres[1] * lags / udist
      lagsx = gres[2] * lags / udist

      # compute appropriately scaled kernel values
      fity = pkern_tvario(pars[['y']], v=sqrt(v), d=lagsy)
      fitx = pkern_tvario(pars[['x']], v=sqrt(v), d=lagsx)
      fit = nug + ( fitx * fity )
    }

    # add line plot showing semivariance at requested lags
    graphics::lines(lags, fit)
  }

  return(invisible())
}


