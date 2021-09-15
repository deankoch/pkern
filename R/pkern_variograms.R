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
#' @param nug nonegative numeric, the variance of the nugget effect
#' @param d vector of nonegative numeric spatial lags to evaluate
#'
#' @return either a vector the same length as `d`, or an anonymous function of distance
#' @export
#'
#' @examples
#' d = 1:10
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
#' @param dims vector c(nx, ny) of positive integers, the number of grid lines
#' @param vec numeric vector of data, in column vectorized order (length `prod(dims)`)
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
  if( length(ds) == 1 ) stop('ds must have length 2')

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
  # set defaults for variance and nugget
  if( is.null(v) ) v = c(min=0, ini=vario$v, max=4*vario[['v']])
  if( is.null(nug) ) nug = c(min=1e-3, ini=v[2]/2, max=4*vario[['v']])

  # set up defaults and bounds for kernel parameters as needed
  xpars = pkern_bds(xpars, vario[['ds']][1])
  ypars = pkern_bds(ypars, vario[['ds']][2])

  # vectorized bounds for all covariance parameters
  plower = c(v[1], nug[1], xpars[['lower']], ypars[['lower']])
  pinitial = c( v[2], nug[2], xpars[['kp']], ypars[['kp']] )
  pupper = c(v[3], nug[3], xpars[['upper']], ypars[['upper']])

  # build matrices to hold pertinent info from `vario` (omit zero point in first index)
  midx = names(vario) %in% c('x', 'y', 'd1', 'd2')
  mvario = lapply(vario[midx], \(vg) cbind(n=vg[['n']][-1], sep=vg[['sep']][-1], sv=vg[['sv']][-1]))

  # define anonymous objective function for optimizer
  fn = function(pvec, p, mv, ds)
  {
    # pvec, numeric vector of parameters
    # p, list of "x", "y" kernel parameter lists associated with pvec
    # mv, list of matrices ("x", "y" and optionally "d1", "d2") with columns "n", "sep", "sv"
    # ds, numeric vector of grid line spacings in x and y directions

    # extract kernel parameters as lists (omit first element, the variance)
    p = pkern_unpack(p, pvec[-(1:2)])

    # generate theoretical semivariance values along x and y for each lag
    xsv = pkern_tvario(pars=p[['x']], v=pvec[1], nug=pvec[2], d=mv[['x']][,'sep'])
    ysv = pkern_tvario(pars=p[['y']], v=pvec[1], nug=pvec[2], d=mv[['y']][,'sep'])

    # compute weighted sums of squares along both dimensions and return their sum
    wss.x = sum( ( mv[['x']][,'n'] / (xsv)^2 ) * ( ( mv[['x']][,'sv'] - xsv )^2 ) )
    wss.y = sum( ( mv[['y']][,'n'] / (ysv)^2 ) * ( ( mv[['y']][,'sv'] - ysv )^2 ) )
    # wss.x = sum( ( mv[['x']][,'n'] ) * ( ( mv[['x']][,'sv'] - xsv )^2 ) )
    # wss.y = sum( ( mv[['y']][,'n'] ) * ( ( mv[['y']][,'sv'] - ysv )^2 ) )

    # do the same for the 45 degree rotated version if it is supplied
    wss.d1 = wss.d2 = 0
    if( length(mv) == 4 )
    {
      # find unit distance along diagonals
      udist = sqrt(sum(ds^2))

      # generate theoretical semivariance values for each lag on first diagonal
      d1.xsep = ds[1] * mv[['d1']][,'sep'] / udist
      d1.ysep = ds[2] * mv[['d1']][,'sep'] / udist
      xsv.d1 = pkern_tvario(pars=p[['x']], v=sqrt(pvec[1]), d=d1.xsep)
      ysv.d1 = pkern_tvario(pars=p[['y']], v=sqrt(pvec[1]), d=d1.ysep)

      # nugget gets added after the product
      d1sv = pvec[2] + ( xsv.d1 * ysv.d1 )

      # repeat for the other diagonal
      d2.xsep = ds[1] * mv[['d2']][,'sep'] / udist
      d2.ysep = ds[2] * mv[['d2']][,'sep'] / udist
      xsv.d2 = pkern_tvario(pars=p[['x']], v=sqrt(pvec[1]), d=d2.xsep)
      ysv.d2 = pkern_tvario(pars=p[['y']], v=sqrt(pvec[1]), d=d2.ysep)
      d2sv = pvec[2] + ( xsv.d2 * ysv.d2 )

      # compute weighted sums of squares along both dimensions and return their sum
      # wss.d1 = sum( ( mv[['d1']][,'n'] ) * ( ( mv[['d1']][,'sv'] - d1sv )^2 ) )
      # wss.d2 = sum( ( mv[['d2']][,'n'] ) * ( ( mv[['d2']][,'sv'] - d2sv )^2 ) )
      wss.d1 = sum( ( mv[['d1']][,'n'] / (d1sv)^2 ) * ( ( mv[['d1']][,'sv'] - d1sv )^2 ) )
      wss.d2 = sum( ( mv[['d2']][,'n'] / (d2sv)^2 ) * ( ( mv[['d2']][,'sv'] - d2sv )^2 ) )
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
                                                        p=list(x=xpars, y=ypars),
                                                        mv=mvario,
                                                        ds=vario[['ds']]))

  # select the best fit
  idx.best = which.min( sapply(list.optim, \(result) result$value) )
  result.optim = list.optim[[idx.best]]

  # unpack fitted parameters and return in a list
  vfit = setNames(result.optim$par[1], 'fitted')
  nugfit = setNames(result.optim$par[2], 'fitted')
  pfit = pkern_unpack(list(x=xpars, y=ypars), result.optim$par[-(1:2)])
  return( list(x=pfit[['x']], y=pfit[['y']], v=vfit, nug=nugfit, ds=vario[['ds']]) )
}
