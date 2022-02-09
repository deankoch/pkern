#
# pkern_raster.R
# Dean Koch, Oct 2021
# Functions for making and fitting semivariograms
#


#' Theoretical variogram function
#'
#' Returns either the variogram function for the (intrinsically stationary) process
#' specified in kernel parameter list `pars`; or, if distances `d` are supplied, the
#' vector of variogram values evaluated at `d`.
#'
#' `pars` and `d` are passed to the correlation function (`pkern_corr`) to get the
#' correlogram values `C(d)/C(0)` (where `C` is the covariance function), and the result
#' is converted to the variogram by the relation,
#'
#'  `2*gamma(d) = nug + 2 * C(0) * (1  - C(d)/C(0))`,
#'
#' where the nugget variance `nug` (by default 0) is specified in `pars$nug` and the
#' partial sill (`C(0)`, by default 1) in `pars$psill`.
#'
#' When `pars` is a 2d separable kernel (ie when it has elements "y" and "x" recognized
#' by `pkern_corr`), `d` should be a list or matrix containing the "y" and "x" distances
#' (vectors `x` and `y`). The correlogram `C(d)/C(0)` is then formulated as
#' `Cx(x) * Cy(y) ) / C(0)`, where `Cx` and `Cy` are the 1-dimensional correlograms for
#' the kernels `pars$x` and `pars$y`.
#'
#' When `d` is `NULL`, the function returns an anonymous function of distance (the
#' variogram). Otherwise `d` should be a numeric vector of distances, in which case the
#' function returns the variogram evaluated at those distances.
#'
#' @param pars list of kernel parameters (see details)
#' @param d vector of nonegative numeric spatial lags to evaluate
#'
#' @return either a vector the same length as `d`, or an anonymous function of distance
#' @export
#'
#' @examples
#' d = seq(1e2)/10
#' pars = pkern_corr('mat')
#' pkern_tvg(pars)
#' pkern_tvg(pars, d=d)
#'
#' vg = pkern_tvg(pars)
#' plot(d, vg(d))
pkern_tvg = function(pars, d=NULL)
{
  # check whether this is a 1 or 2d kernel definition
  yxnm = c('y', 'x')
  knm = c('k', 'kp')
  is.1d = all( knm %in% names(pars) )

  # extract partial sill and nugget, setting defaults when they are missing
  psill = ifelse(is.null(pars[['psill']]), 1, pars[['psill']]) |> unname()
  nug = ifelse(is.null(pars[['nug']]), 0, pars[['nug']]) |> unname()

  # check whether d (if supplied) is of the correct dimensionality
  if( !is.null(d) )
  {
    # in 1d case we expect a vector of distances
    if( is.1d )
    {
      if( !is.numeric(d) ) stop('d must be a numeric vector in 1d case')

      # coerce matrices to vector
      d = as.vector(d)

    } else {

      # in 2d case we coerce matrices and dataframes to lists of two vectors
      if( is.matrix(d) ) d = data.frame(d)
      if( is.data.frame(d) ) d = as.list(d)

      # assume (with a warning) order "y", "x" when d isn't named
      if( !all( yxnm %in% names(d) ) )
      {
        warning('assuming the first vector in d is "y" and second is "x"')
        d[1:2] = stats::setNames(d[1:2], yxnm)
      }
    }
  }

  # handle 1d case
  if( is.1d )
  {
    # variogram function
    fn = function(d) { 2 * ( nug + ( psill * ( 1 - pkern_corr(pars, d) ) ) ) }

    # return either the function or its evaluation
    if( is.null(d) ) return(fn)
    return( fn(d) )
  }

  # handle 2d case
  fn = function(d) {

    # taking product of correlation functions
    corrvals = do.call('*', Map(\(p, dyx) pkern_corr(p, dyx), p=pars[c('y', 'x')], dyx=d))
    return( 2 * ( nug + ( psill * ( 1 - corrvals ) ) ) )
  }

  # return either the function or its evaluation
  if( is.null(d) ) return(fn)
  return( fn(d) )

}


# TODO: optimization
#' Sample variogram for regular 2D gridded data along x direction
#'
#' The function estimates the variogram ("2-gamma") of point pairs located on the
#' regular grid with dimensions `gdim`, having (detrended) data vector(s) `vec`. Estimation
#' is done separately at each of the supplied x grid line lags, `lags`. Elements of `lags`
#' must be a subset of `seq(gdim[2]-1)`.
#'
#' Results are returned in a named list containing the sampled lags (`lags`), the number of
#' point pairs sampled at each lag (`n`) and the variogram estimate (`vg`) at each lag;
#' If `simple=FALSE`, the function also returns `pairs`, a list of matrices containing the
#' indices of point pairs sampled and their absolute differences.
#'
#' `vec` can be a vector (in column vectorized order), or a list of such vectors for
#' repeated measurements at the same spatial locations. NAs are permitted in `vec` (but
#' point pairs containing NAs are excluded from computations).
#'
#' `nmax` specifies a maximum number of point pairs to process at each lag. When applicable,
#' a random sample of (`nmax`) eligible point pairs are taken at each lag. When `vec` is
#' a list of vectors, the function always samples at least one point pair from each of them.
#' Note that this means `nmax` is not respected when `vec` is a list and `nmax < length(vec)`.
#' If the function is too slow in this case, try taking a subset of `vec`.
#'
#' fit.method "matheron" uses the method-of-moments formula of Matheron (1962), which takes the
#' average of the squared differences at each lag; "rmean" and "rmedian" are the robust methods
#' of Cressie and Hawkins (1980), which use mean and median in a 4th-root transformation with
#' a bias adjustmnet. See Cressie (1993, Ch. 2.4.3) for a review of all three methods.
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
#' @param fit.method character, one of "matheron", "rmean", "rmedian"
#' @param quiet logical, suppresses console output
#' @param nmax positive integer, the maximum number of point pairs to sample at each lag
#'
#' @return a list with named elements "lags", "n", "vg", and optionally "pairs"
#' @export
#'
#' @examples
#'
#' gdim = c(1e2, 50)
#' pars = pkern_corr('gau') |> utils::modifyList(list(psill=1))
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
pkern_xvario = function(gdim, vec, simple=TRUE, fit.method='rmedian', quiet=FALSE, lags=NULL, nmax=1e4)
{
  # set default sample lags and express vector `vec` as length-1 list
  if( is.null(lags) ) lags = seq( gdim[2] - 1 )
  if( !is.list(vec) ) vec = list(vec)

  # prepare to sample along different grid line lags in a loop
  lags = sort(lags)
  n.lags = length(lags)
  plist = vector(mode='list', length=n.lags)
  vgx = rep(NA, n.lags)
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

      # formulate the variogram statistic by various methods
      dabs.vec = unlist(dabs, use.names=FALSE)
      if(fit.method=='matheron') vgx[idx.dj] = ( mean(dabs.vec)^2 )
      if(fit.method=='rmedian') vgx[idx.dj] = (stats::median( sqrt(dabs.vec) )^4) / 0.457
      if(fit.method=='rmean')
      {
        ch.factor = 0.457 + ( 0.494 / n[idx.dj] ) #+ ( 0.494 / ( n[idx.dj]^2 ) )
        vgx[idx.dj] = mean( sqrt(dabs.vec) )^4 / ch.factor
      }
    }
  }
  if( !quiet ) close(pb)

  if( simple ) return( list(lags=lags, n=n, vg=vgx) )
  return( list(lags=lags, n=n, vg=vgx, pairs=plist) )
}


# TODO: refactoring "gdim" > "g" (list)
#' Directional sample variograms for regular 2D gridded data
#'
#' A wrapper for calls to `pkern_xvario` to get the directional sample variograms
#' from regular lattice data along the y and x dimensions, and optionally along the
#' diagonals.
#'
#' The function makes the appropriate calls to `pkern_xvario` and returns the results
#' in list elements "y", "x". If `diagonal=TRUE`, the function also returns "d1" and
#' "d2", the y and x semivariograms after a 45 degree counterclockwise rotation of the
#' input grid (see `pkern_r45`) assuming equal y and x resolution. For unequal resolution
#' the rotation angle is `atan(gres[2]/gres[1])`.
#'
#' For convenience, `gdim` can be a list containing both the 'gdim' and 'gres' arguments.
#' If the sub-list `gdim$sg` exists, the function copies these arguments from `gdim$sg$gdim`
#' and `gdim$sg$gres`. Otherwise it looks for them in `gdim$gdim` and `gdim$gres`.
#' For example, `gdim` can be the output of `pkern_snap` (to sample from snapped points) or
#' `pkern_fromRaster` (to sample from raster values).
#'
#' You can (and, with large grids, probably should) set `dmax` to specify a maximum
#' sampling distance. This both speeds computations, and omits pathological tail behavior
#' from variogram fitting downstream. `dmax` should be given in the same units as `gres`.
#' Output lags are also given in these units, as are the x-axis values drawn by
#' `pkern_vario_plot`.
#'
#' When `gres` is not supplied (directly, or in `gdim`), the function assumes
#' that the x and y grid lines have spacing 1.
#'
#' @param gdim vector c(ny, nx) of positive integers, the number of grid lines
#' @param vec numeric vector of data, in column vectorized order (length `prod(dims)`)
#' @param dmax numeric, the maximum lag distance (only shorter lags are sampled)
#' @param fit.method character, one of "matheron", "rmean", "rmedian", passed to `pkern_xvario`
#' @param diagonal logical, if TRUE, semivariogram results are returned also for diagonals
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
#' pars.1d = pkern_corr('gau')
#' pars.2d = list(y=pars.1d, x=pars.1d, psill=1)
#' vec = pkern_sim(pars.2d, gdim)
#'
#' # sample first 10 lags and plot
#' vario = pkern_vario(gdim, vec, lags=seq(10), quiet=TRUE)
#' pkern_vario_plot(vario, pars.2d)
#'
#' # example with unequal resolution and multiple replicates
#' gres = c(2,1)
#' nrep = 10
#' vec = lapply(seq(nrep), \(i) pkern_sim(pars.2d, gdim, gres=gres, makeplot=FALSE))
#' vario = pkern_vario(gdim, vec, lags=seq(10), gres=gres, quiet=TRUE)
#' pkern_vario_plot(vario, pars.2d)
#'
pkern_vario = function(gdim, vec, dmax=NA, fit.method='rmedian', diagonal=TRUE,
                       gres=NA, quiet=FALSE, nmax=1e4, vcalc=TRUE)
{
  # handle invalid fit.method
  if( !( fit.method %in% c('matheron', 'rmean', 'rmedian') ) ) stop('invalid fit method')

  # handle `gdim` as output from `pkern_snap`
  if( is.list(gdim) )
  {
    # if subgrid info is provided, it is assumed `vec` is a data vector from the subgrid
    if( !is.null(gdim[['sg']]) ) gdim = gdim[['sg']]
    if( is.na(gres) & !is.null(gdim[['gres']]) ) gres = gdim[['gres']]
    gdim = gdim[['gdim']]
  }

  # compute sample variance
  if( !is.list(vec) ) vec = list(vec)
  svar = ifelse(vcalc, stats::var(unlist(vec), na.rm=TRUE), NA)

  # set default resolution (1,1)
  yxnm = c('y', 'x')
  if( anyNA(gres) ) gres = c(1,1)
  if( length(gres) == 1 ) stop('gres must have length 2')
  if( !all( yxnm %in% names(gres) ) ) gres = stats::setNames(gres, yxnm)

  # set default sample lags
  #if( anyNA(lags) ) lags = lapply(gdim, \(d) seq(d-1))
  #if( !is.list(lags) ) lags = list(y=lags, x=lags)

  # set default dmax to diagonal span of grid (samples all lags)
  diagstep = sqrt( sum( gres^2 ) )
  if( is.na(dmax) ) dmax = min( gdim ) * diagstep

  # set sample lags based on dmax (always include at least one)
  lags = ceiling(dmax/gres) |> pmin(gdim-1) |> pmax(c(1,1)) |> as.list() |> lapply(seq)

  # compute x semivariances
  if( !quiet ) cat('\nidentifying horizontal lags...')
  xout = pkern_xvario(gdim, vec, fit.method=fit.method, quiet=quiet, lags=lags[['x']], nmax=nmax)

  # repeat in row-vectorized order to get y semivariances
  if( !quiet ) cat('\nidentifying vertical lags...')
  vec.rv = lapply(vec, \(v) v[pkern_r2c(gdim, in.byrow=FALSE, out.byrow=TRUE)])
  yout = pkern_xvario(rev(gdim), vec.rv, fit.method=fit.method, quiet=quiet, lags=lags[['y']], nmax=nmax)

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

    # find the dimensions of rotated grid - note that vario45 inflates distances by factor 2
    gdim45 = dim( pkern_r45(vec[[1]], gdim) )
    gres45 = c(y=diagstep/2, x=diagstep/2)

    # recursive call to get component variograms in new coordinate system
    vario45 = pkern_vario(gdim45, vec45, fit.method=fit.method, quiet=quiet, diagonal=F, gres=gres45, nmax=nmax, vcalc=F, dmax=dmax)
    list.out = c(list.out, list(d1=vario45[['y']], d2=vario45[['x']]))
  }

  # omit entries with NA semivariances
  list.out = lapply(list.out, \(d) lapply(d, \(z) z[ !is.na(d[['vg']]) ] ) )

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
#' Sill and nugget are estimated simultaneously with the correlation kernel parameters
#' using the bounds and initial values in `psill` and `nug`. If these are not
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
#' 'n/d' uses sample size divided by squared distance (the default here, as in `gstat`);
#' and 'n/v uses the robust estimator of Cressie (1993, Ch. 2.6.2), which divides sample
#' size by squared theoretical semivariance.
#'
#' when `xpars` is set to `NULL`, the function fits only the y parameters, and the x
#' kernel is made identical to the y kernel. For example with the Gaussian kernel, this
#' enforces isotropy.
#'
#'
#' @param vario list, the return value from a call to `pkern_vario`
#' @param ypars character or list, recognizable (as `pars`) by `pkern_corr`
#' @param xpars character or list, recognizable (as `pars`) by `pkern_corr`
#' @param psill length-3 numeric vector of bounds (lower, initial, upper) for the partial sill
#' @param nug length-3 numeric vector of bounds (lower, initial, upper) for the nugget effect
#' @param add positive integer, the number of initial parameter sets tested (see details)
#' @param dmax positive numeric, point pairs with separation distances exceeding `dmax` are ignored
#' @param fit.method character string, either '1', 'n', 'n/d', or 'n/v' (the default, see details)
#'
#' @return TODO
#' @export
#'
#' @examples
#' gdim = c(1e2, 2e2)
#' pars.1d = pkern_corr('gau')
#' pars.2d = list(y=pars.1d, x=pars.1d, psill=1)
#' vec = pkern_sim(pars.2d, gdim)
#'
#' # sample first 10 lags and plot
#' vario = pkern_vario(gdim, vec, lags=seq(10), nmax=1e3, quiet=TRUE)
#' pkern_vario_plot(vario, pars.2d)
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
pkern_vario_fit = function(vario, ypars='mat', xpars=ypars, psill=NULL, nug=NULL, add=0, dmax=Inf,
                           fit.method='n/d', maxit=1e4)
{

  # TODO: find better initials for nug and psill

  # set default bounds for nugget effect
  vario_l1med = median(sapply(vario[c('x', 'y', 'd1', 'd2')], \(x) x[['vg']][2] ))
  vario_mode = sapply(vario[c('x', 'y', 'd1', 'd2')], \(x) x[['vg']][-1] ) |> unlist()
  vario_v =  sqrt(vario[['v']])
  if( is.null(nug) ) nug = c(min=1e-6, ini=vario_l1med/2, max=vario_l1med)
  if( is.null(psill) ) psill = c(min=vario_l1med, ini=pmax(vario_l1med, vario_v), max=2*vario_v)

  # set defaults for partial sill (first parameter) and nugget (second parameter)
  #svmin = sapply(vario[c('y', 'x', 'd1', 'd2')], \(x) min(x[['vg']][-1])) |> min()

  # truncate max distance to max lag, compute number of cells this represents
  dmax = sapply(vario[c('y', 'x', 'd1', 'd2')], \(x) max(x[['lags']])) |> max() |> pmin(dmax)
  nmax = ceiling(dmax/vario[['gres']])

  # identical kernel parameter strings prompt initial constrained fit to get initial values
  ykp_ini = NULL
  if( is.character(ypars) & identical(ypars, xpars) )
  {
    ini_result = pkern_vario_fit(vario, ypars, NULL, psill, nug, add, dmax, fit.method)
    ykp_ini = ini_result[['y']][['kp']]
    nug['ini'] = ini_result[['nug']]
    psill['ini'] = ini_result[['psill']]
  }

  # set up defaults and bounds for y kernel parameters
  ypars = pkern_bds(ypars, gres=vario[['gres']]['y'], nmax=nmax['y'])
  if( !is.null(ykp_ini) ) ypars[['kp']] = ykp_ini

  # construct x kernel as needed
  if( !is.null(xpars) )
  {
    xpars = pkern_bds(xpars, gres=vario[['gres']]['x'], nmax=nmax['x'])
    if( !is.null(ykp_ini) ) xpars[['kp']] = ykp_ini
  }

  # vectorized bounds for all covariance parameters
  plower = c(psill[1], nug[1], ypars[['lower']], xpars[['lower']])
  pinitial = c( psill[2], nug[2], ypars[['kp']], xpars[['kp']] )
  pupper = c(psill[3], nug[3], ypars[['upper']], xpars[['upper']])

  # check for dmax set too low
  vario_valid = sapply(vario[c('x', 'y', 'd1', 'd2')], \(x) sum( x[['lags']] < dmax ) > 0 )
  vario_good = sapply(vario[c('x', 'y', 'd1', 'd2')], \(x) sum( x[['lags']] < dmax ) > 1 )
  if( !all(vario_valid) ) stop('Not enough lags sampled. Try increasing dmax')
  if( !all(vario_good) ) warning('Not enough lags sampled. Try increasing dmax')

  # define anonymous objective function for optimizer
  fn = function(pv, p, vario_in, gres, fit.method, dmax)
  {
    # pv, numeric vector of parameters
    # p, list of "x", "y" kernel parameter lists associated with pv
    # vario, list, the output of pkern_vario
    # gres, numeric vector of grid line spacing distancec in x and y directions

    # when x kernel parameters not supplied, the x kernel becomes a copy of the y kernel
    if( is.null(p[['x']]) )
    {
      # copy kernel structure and parameters
      p[['x']] = p[['y']]
      pv = c(pv, pv[-(1:2)])
    }

    # extract kernel parameters as lists
    pars_in = c(pkern_unpack(p, pv[-(1:2)]), list(psill=pv[1], nug=pv[2], gres=gres))
    pars_in_x = c(pars_in[['x']], modifyList(pars_in, list(x=NULL, y=NULL)))
    pars_in_y = c(pars_in[['y']], modifyList(pars_in, list(x=NULL, y=NULL)))

    # compute weighted sum of squares along x direction
    wss.x = pkern_vario_wss(vario_in, nm='x', pars=pars_in_x, fit.method=fit.method, dmax=dmax)

    # repeat for y direction
    wss.y = pkern_vario_wss(vario_in, 'y', pars=pars_in_y, fit.method=fit.method, dmax=dmax)

    # do the same for the 45 degree rotated version if it is supplied
    wss.d1 = wss.d2 = 0
    if( all( c('d1', 'd2') %in% vario_in ) )
    {
      wss.d1 = pkern_vario_wss(vario_in, nm='d1', pars=pars_in, fit.method=fit.method, dmax=dmax)
      wss.d2 = pkern_vario_wss(vario_in, nm='d2', pars=pars_in, fit.method=fit.method, dmax=dmax)
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
  list_optim = lapply(list.ini, \(ini) stats::optim(par=ini,
                                                    f=fn,
                                                    method='L-BFGS-B',
                                                    lower=plower,
                                                    upper=pupper,
                                                    p=list(y=ypars, x=xpars),
                                                    vario_in=vario,
                                                    gres=vario[['gres']],
                                                    fit.method=fit.method,
                                                    dmax=dmax,
                                                    control=list(maxit=maxit,
                                                                 parscale=pupper-plower)))

  # select the best fit
  idx_best = which.min( sapply(list_optim, \(result) result$value) )
  result_optim = list_optim[[idx_best]]

  # unpack fitted sill, nugget, copy the rest
  psill_optim = stats::setNames(result_optim[['par']][1], 'fitted')
  nug_optim = stats::setNames(result_optim[['par']][2], 'fitted')
  kp_optim =  result_optim[['par']][-(1:2)]

  # handle NULL x kernel parameters (copy from y)
  if( is.null(xpars) )
  {
    # copy kernel structure and parameters
    xpars = ypars
    kp_optim = rep(kp_optim, 2)
  }

  # unpack kernel parameters and finish
  pfit = pkern_unpack(list(y=ypars, x=xpars), kp_optim)
  return( list(y=pfit[['y']], x=pfit[['x']], psill=psill_optim, nug=nug_optim, gres=vario[['gres']]))
}

# TODO: write up documentation for this
# helper for pkern_vario_fit
# pars=c(p, list(psill=psill, nug=nug)
pkern_vario_wss = function(vario_in, nm, pars, fit.method, dmax=dmax)
{
  # check for valid direction name argument
  if( !( nm %in% c('x', 'y', 'd1', 'd2') ) ) stop('nm must be one of "x", "y", "d1", "d2"')

  # extract lags data from variogram list (omit lag-0 first entry)
  lags = vario_in[[nm]][['lags']][-1]
  idx_valid = lags < dmax
  if( sum(idx_valid) == 0 ) return(Inf)
  lags_valid = lags[idx_valid]

  # handle special case of diagonals - pars must be a list with both x and y component

  if( any( c('d1', 'd2') %in% nm ) )
  {
    # check for required input in  this case
    if( !( c('gres') %in% names(pars) ) ) stop('gres must be included in pars')

    # compute component-wise distances for each lag
    udist = sqrt( sum( pars[['gres']]^2 ) )
    lags_valid = lapply(pars[['gres']], \(gr) gr * lags_valid / udist)
  }

  # theoretical variogram values
  tvg = pkern_tvg(pars, d=lags_valid)

  # set up weights
  n_valid = vario_in[[nm]][['n']][-1][idx_valid]
  vg_valid = vario_in[[nm]][['vg']][-1][idx_valid]
  w = switch(fit.method,
             '1' = 1,
             'n' = n_valid,
             'n/v'= n_valid/(vg_valid^2),
             'n/d'= n_valid/(lags[idx_valid]^2) )

  # compute weighted sums of squares along both dimensions and return their sum
  return( sum( w * ( ( vg_valid - tvg )^2 ) ) )
}


#' Transform data forward and back
#'
#' Apply a transformation to the data vector(s) in `x` and return the back-transform function.
#'
#' `x` can be a numeric vector or a list of them. Data vector lists are handled in
#' `pken_transform` (and its returned `bfun` function) by applying the (back) transform to
#' each element using lapply, returning data as a list with the same structure as the input.
#'
#' `x` can also be the return value of pkern_transform (a list with elements "z", "nm", "bfun"),
#' in which case the transform is applied to `x$z`, and the back-transform function is composed
#' appropriately with the previous one. This allows calls to pkern_transform to be chained
#' together (see examples).
#'
#' Available functions (valid choices for `nm`) are:
#'
#'  "linear": applies `base::scale` with default options
#'  "log": uses log(x + `log_offset`)
#'  "nscore": uses the normal score quantiles returned (in `y`) by `base::qqnorm`
#'  "nscore_smooth": same as "nscore"
#'
#' "nscore" and "nscore_smooth" differ in their method of constructing the back-transform;
#' The former uses linear interpolation (`stats::approx`) and the latter uses spline
#' interpolation (`stats::spline`).In both cases, the interpolators are constrained to
#' predict within the range of `x`, by padding their training data with extreme values from
#' both the transformed `x` input and the input `z` values. This implementation of "nscore"
#' is adapted from `mstats::backtr` and `mstats::nscore_model`.
#'
#' `log_offset` is a small additive constant to ensure strict positivity of inputs to log().
#' It is ignored when `nm` is not "log".
#'
#' @param x, numeric data vector to transform (or list of them), or the output of pkern_transform
#' @param nm, the transformation name, one of 'linear', 'log', 'nscore', 'nscore_smooth'
#' @param log_offset, positive numeric, ignored if `nm` is not "log"
#'
#' @return list with: transformed data "z", transform name "nm", and back-transform function "bfun"
#' @export
#'
#' @examples
#'
#' # generate large spread of random data
#' random_data = runif(-1e6, 1e6, n=1e3)
#'
#' # log tranform:
#' transf_result = random_data |> pkern_transform('log')
#'
#' # the function returns in a list:
#' z = transf_result$z # the transformed data
#' transf_result$nm |> print() # the transform name
#' bfun = transf_result$bfun # the back-transform function
#'
#' # round trip accuracy is fine
#' (bfun(z) - random_data) |> abs() |> max()
#'
#' # transform calls can be composed
#' transf_result = random_data |> pkern_transform('linear') |> pkern_transform('log') |> pkern_transform('linear')
#' print(transf_result$nm)
#' (transf_result$bfun(transf_result$z) - random_data) |> abs() |> max()
#'
#' # both pkern_transform and its returned back-transform function are vectorized for lists
#' random_data = lapply(1:10, \(x) runif(-1e6, 1e6, n=1e3))
#' transf_result = random_data |> pkern_transform('linear') |> pkern_transform('log')
#' btrans_result = transf_result$bfun(transf_result$z)
#' str(btrans_result)
#' mapply(\(x, bx) max(abs(x-bx)), btrans_result, random_data) |> max()
#'
pkern_transform = function(x, nm='nscore_smooth', log_offset=1e-6)
{
  # process list arguments to x
  is_data_vector = FALSE
  model_x = identity
  nm_x = ''
  if( is.list(x) )
  {
    # check for data vector lists
    is_data_vector = all( sapply(x, is.numeric) )

    # handle list input to x that is not a data vector
    if( !is_data_vector )
    {
      # checking for valid list structure (expecting chained call)
      model_elements = c('z', 'nm', 'bfun')
      err_elements = paste('list input to x must contain:', paste(model_elements, collapse=', '))
      if( !all( c('z', 'bfun') %in% names(x) ) ) stop(err_elements)
      if( !is.function( x[['bfun']] ) ) stop('x$bfun must be a function')

      # copy inverse function for input list x then overwrite x with data values
      model_x = x[['bfun']]
      nm_x = x[['nm']]
      x = x[['z']]

      # set data vector flag for chained calls
      if( is.list(x) ) is_data_vector = all( sapply(x, is.numeric) )
    }
  }

  # handle data vector lists
  if( is_data_vector )
  {
    # recursive call with vectorized data
    x_recursive = list(z=unlist(x), nm=nm_x, bfun=model_x)
    transform_res = pkern_transform(x_recursive, nm, log_offset)

    # reshape results as list in 'z'
    transform_res[['z']] = utils::relist(transform_res[['z']], x)
    return(transform_res)
  }

  # drop names from x
  x = unname(x)

  # define valid model choices
  nm_str = c('linear', 'log', 'nscore', 'nscore_smooth')
  err_msg = paste('tpars$model must be one of:', paste(nm_str, collapse=', '))

  # check for valid model specification
  if( !(nm %in% nm_str) ) stop(err_msg)

  # scale and center
  if(nm == 'linear')
  {
    scale_out = scale(x)
    x_out = scale_out |> c()
    scale_centre = attr(scale_out, 'scaled:center')
    scale_mult = attr(scale_out, 'scaled:scale')

    # define the back-transform function
    model_back = \(z) {

      if( !is.list(z) ) z = list(z)
      z_out = lapply(z, \(z_in) scale_centre + ( scale_mult * z_in ) )
      if( length(z_out) == 1 ) z_out = unname( unlist(z_out) )
      return(z_out)
    }
  }

  # log transform
  if(nm == 'log')
  {
    # re-centre so that minimum is log_offset
    scale_min = log_offset - min(x)
    x_out = log(x + scale_min)

    # define the back-transform function
    model_back = \(z) {

      if( !is.list(z) ) z = list(z)
      z_out = lapply(z, \(z_in) exp(z_in) - scale_min )
      if( length(z_out) == 1 ) z_out = unname( unlist(z_out) )
      return(z_out)
    }
  }

  # normal score transforms
  if( nm %in% c('nscore', 'nscore_smooth') )
  {
    # compute normal scores and sorted x values with padding
    x_out = qqnorm(x, plot.it=FALSE)[['x']]
    x_pad = c(min(x), sort(x), max(x))

    # anonymous function to make y values matching x_pad
    y_pad = \(z) c(min(c(x_out, z)), sort(x_out), max(x_out, z))

    # use linear interpolation for back-transform
    if( nm == 'nscore' )
    {
      model_back = \(z) {

        if( !is.list(z) ) z = list(z)
        z_out = lapply(z, \(z_in) stats::approx(y_pad(z_in), x_pad, xout=z_in, ties=mean)[['y']] )
        if( length(z_out) == 1 ) z_out = unname( unlist(z_out) )
        return(z_out)
      }
    }

    # use spline interpolation for back-transform
    if( nm == 'nscore_smooth' )
    {
      model_back = \(z) {

        if( !is.list(z) ) z = list(z)
        z_out = lapply(z, \(z_in) stats::spline(y_pad(z_in), x_pad, xout=z_in, ties=mean)[['y']] )
        if( length(z_out) == 1 ) z_out = unname( unlist(z_out) )
        return(z_out)
      }
    }
  }

  # compose the new back-transform function with existing one
  model_out = \(z) { model_x(model_back(z)) }
  nm_out = paste0(nm, '(', nm_x, ')')

  # return the transformed values and the inverse function in a list
  return( list(z=x_out, nm=nm_out, bfun=model_out) )
}


# TODO: update docs
#' Plot 1-dimensional empirical semivariograms for a separable model
#'
#' Plots the output of `pkern_vario` or `pkern_xvario`, optionally including the theoretical
#' semivariogram values for a model specified in `pars` (eg the return value of `pkern_vario_fit`).
#' Note that the plot shows semivariance (gamma), whereas the "vg" elements returned by
#' `pkern_vario` are variogram values (2*gamma).
#'
#' `pars` should be a list containing named elements: "y" and "x", lists of y and x kernel
#' parameters in form recognized by `pkern_corr`; "psill", the partial sill and, optionally,
#' "nug", the nugget variance.
#'
#' `plotpars` applies to all plots. It should be a list containing any of the following elements:
#'
#'    "title" character or vector of two, title(s) for the plot(s)
#'    "shade" logical, indicating to shade points according to the number of samples
#'    "dmax" positive numeric, upper limit for x in scatterplot
#'    "vglim" length-2 numeric vector, the y limits for the scatterplot
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
#' pars.1d = pkern_corr('gau')
#' pars.2d = list(y=pars.1d, x=pars.1d, psill=1)
#' vec = pkern_sim(pars.2d, gdim)
#'
#' # sample first 10 lags and plot results
#' vario = pkern_vario(gdim, vec, lags=seq(10), quiet=TRUE)
#' pkern_vario_plot(vario, pars.2d)
#'
#' # plot componentwise
#' pkern_vario_plot(vario[['x']], c(pars.1d, list(psill=1)))
#' pkern_vario_plot(vario[['y']], c(pars.1d, list(psill=1)))
#'
pkern_vario_plot = function(vario, pars=NULL, plotpars=NULL)
{
  # intialize empty list `plotpars` as needed and set some defaults
  if( is.null(plotpars) ) plotpars = list()
  if( is.null(plotpars[['shade']]) ) plotpars[['shade']] = TRUE

  # identify the semivariance directions provided in vario
  nm.yx = c(y='y', x='x')
  nm.dyx = c(d1='d1', d1='d2')
  nm.vario = c(nm.yx, nm.dyx)
  is.vario = stats::setNames(nm.vario %in% names(vario), nm.vario)

  # handle 2-dimensional case
  if( all( is.vario[nm.yx] ) )
  {
    # clean up pars if supplied
    if( !is.null(pars) )
    {
      # if a 1d kernel supplied, assume we are assigning it to both y and x component
      if( all( c('k', 'kp') %in% names(pars) ) ) pars = utils::modifyList(pars, list(y=pars, x=pars))

      # assume order "y" "x" unless otherwise indicated by names
      if( !all( nm.yx %in% names(pars) ) ) names(pars) = nm.yx
    }

    # assign default partial sill and nugget if missing and split parameters into component lists
    if( !is.null(pars) & is.null(pars[['psill']]) ) pars[['psill']] = 1
    if( !is.null(pars) & is.null(pars[['nug']]) ) pars[['nug']] = 0
    xpars = c( pars[['x']], list( psill=pars[['psill']], nug=pars[['nug']] ) )
    ypars = c( pars[['y']], list( psill=pars[['psill']], nug=pars[['nug']] ) )

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
    vglim.def = FALSE
    if( is.null(plotpars[['vglim']]) )
    {
      vglim.def = TRUE
      idx.valid = lapply( vario[nm.yx], \(xy) xy[['lags']] < plotpars[['dmax']] )
      vgmax = mapply(\(yx, i) max( yx[['vg']][i] ), yx=vario[nm.yx], i=idx.valid)
      plotpars[['vglim']] = c(0, max(vgmax))
    }

    # update default plot settings for diagonals
    if( all( is.vario ) )
    {
      # update axis limits as needed
      if( vglim.def )
      {
        # x axis
        lagsmax = sapply( vario[nm.dyx], \(xy) max( xy[['lags']] ) )
        if( dmax.def ) plotpars[['dmax']] = max(plotpars[['dmax']], lagsmax)

        # y axis
        idx.valid = lapply( vario[nm.vario], \(xy) xy[['lags']] < plotpars[['dmax']] )
        vgmax = mapply(\(xy, i) max( xy[['vg']][i] ), xy=vario[nm.vario], i=idx.valid)
        plotpars[['vglim']] = range(plotpars[['vglim']], vgmax)
      }

      # assign plot titles as needed
      if( title.def )
      {
        degval = round( atan( vario[['gres']][1] / vario[['gres']][2] )*180/pi, 1)
        dstr = paste0('(', degval, ' degrees)')
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
  if( is.null(plotpars[['vglim']]) ) plotpars[['vglim']] = c(0, max(vario[['vg']]))
  xlim = c(0, plotpars[['dmax']])
  ylim = plotpars[['vglim']]
  main = plotpars[['title']]

  # catch empty variogram
  if( length( vario[['n']] ) == 1 )
  {
    # initialize an empty plot
    plot(xlim, ylim, pch=NA, xlab='lag', ylab='semivariance', main=main)

  } else {

    # unpack arguments (exclude the zero point)
    n = vario[['n']][-1]
    lags = vario[['lags']][-1]
    vg = vario[['vg']][-1]

    # map point shading to the number of point pairs sampled if requested
    #colmap = 'black'
    #if( plotpars[['shade']] ) colmap = rev( grDevices::gray.colors( max(n) ) )[n]

    # make the scatterplot of sampled semivariances
    plot(lags, vg/2, xlab='lag', ylab='semivariance', xlim=xlim, ylim=ylim/2, main=main, pch=16)

    # add bubbles indicating sample sizes
    #points(lags, vg, pch=1, cex=n/max(n))

    # # make the scatterplot of sampled semivariances
    # plot(lags, vg, xlab='lag', ylab='semivariance',
    #      xlim=xlim, ylim=ylim, main=main, pch=20, col=colmap)

  }

  # if kernel parameters are supplied, plot the theoretical values
  if( !is.null(unlist(pars)) )
  {
    # check for variance components in parameters list and set defaults
    if( is.null(pars[['psill']]) ) pars[['psill']] = 1
    if( is.null(pars[['nug']]) ) pars[['nug']] = 0

    #  handle x or y kernel component requests
    lags = seq(xlim[1], xlim[2], length.out=1e2)
    if( all( c('k', 'kp') %in% names(pars) ) ) fit = pkern_tvg(pars, d=lags)

    # handle requests along a diagonal (both y and x components included in pars)
    if( all( nm.yx %in% names(pars) ) )
    {
      # find unit distance along diagonals
      gres = vario[['gres']]
      udist = sqrt( sum( gres^2 ) )

      # find the true distance of lags
      dlags = lapply(gres, \(gr) gr * lags / udist)

      # compute appropriately scaled kernel values
      fit = pkern_tvg(pars, d=dlags)
      #fit = pkern_tvg(pars, d=lags)
    }

    # add line plot showing semivariance at requested lags
    graphics::lines(lags, fit/2)
  }

  return(invisible())
}


