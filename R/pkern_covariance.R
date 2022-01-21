#
# pkern_raster.R
# Dean Koch, Oct 2021
# Functions for handling separable covariance functions
#


#' Stationary 1D correlation kernels
#'
#' Computes stationary correlation function values for the n (nonegative) 1-dimensional
#' distances in `d`. Parameter list entry `pars$kp` supplies the kernel parameter(s).
#'
#' `pars$k` must be one of the following kernel names:
#'
#' (1 parameter)
#'
#' "exp": exponential (special case of "gex" with shape "p"=1)
#' "gau": gaussian/stable (special case of "gex" with shape "p"=2)
#' "sph": spherical (AKA stable/Gaussian for p=2)
#'
#' (2 parameters)
#'
#' "gex": gamma-exponential (with shape "p")
#' "mat": Whittle-Matern (Handcock and Wallis parameterization, with shape "kap")
#'
#' For the 1-parameter kernels, `pars$kp` is the range parameter value ("rho" ); For the
#' 2-parameter kernels, `pars$kp` is a vector whose first element is "rho", and second
#' element is the shape parameter ("p" or "kap"). Note that names in `pars$kp` are ignored
#' and only the order matters - the range parameter must come first.
#'
#' When `d` is NA (or contains at least one NA value), or if instead of a list `pars` is
#' a kernel name string (or a vector of two of them), then `d` is ignored and the function
#' passes `pars` to `pkern_bds` and returns the result. This sets suggested parameter values
#' and bounds wherever they are missing from `pars`.
#'
#' @param pars list with elements "k", the kernel name, and "kp" the parameter vector
#' @param d numeric vector of length n, the distances to evaluate
#'
#' @return length-n vector or a list of parameters and bounds (see details)
#' @export
#'
#' @examples
#'
#' # define a 1D kernel and grab values
#' pars = list(k='mat', kp=c(10, 2))
#' pkern_corr(pars, d=1:10)
#'
#' # get suggested bounds
#' pkern_corr(pars)
#' pkern_corr('mat')
#'
pkern_corr = function(pars, d=NA)
{
  # handle special requests for bounds and initial values
  if( is.character(pars) | anyNA(d) ) return( pkern_bds(pars) )

  # handle invalid pars
  if( !all( c('k', 'kp') %in% names(pars) ) ) stop('pars must be list with elements "k" and "kp"')

  # exponential
  if(pars[['k']] == 'exp')
  {
    # new parameter list for recursive call
    pars = list(k='gxp', kp=c(rho=pars[['kp']][1], p=1L))
    return( pkern_corr(pars, abs(d)) )
  }

  # gaussian/stable
  if(pars[['k']] == 'gau')
  {
    # new parameter list for recursive call
    pars = list(k='gxp', kp=c(rho=pars[['kp']][1], p=2L))
    return( pkern_corr(pars, d) )
  }

  # spherical
  if(pars[['k']] == 'sph')
  {
    # assign parameter and process truncation distance
    ds = d / pars[['kp']][1]
    cvals = rep(0, length(d))
    idx.nz = ds < 1

    # evaluate on non-truncated distances and return
    cvals[idx.nz] = 1 - ( (3/2) * ds[idx.nz] ) + ( (1/2) * ( ds[idx.nz]^3 ) )
    return( cvals )
  }

  # gamma-exponential
  if(pars[['k']] == 'gxp')
  {
    # return suggested bounds and initial value if requested
    if( length( pars[['kp']] ) != 2 ) stop(paste('pars$kp must be a vector of form c(rho, p)'))

    # assign parameters and evaluate
    return( exp( -( ( d / pars[['kp']][1] )^pars[['kp']][2] ) ) )
  }

  # Whittle-Matern
  if(pars[['k']] == 'mat')
  {
    # assign parameters and compute scaling constant
    kap = pars[['kp']][2]
    rho = pars[['kp']][1] / ( 2 * sqrt(kap) )
    sc = ( (2^( 1 - kap ) ) / gamma(kap) )

    # some preprocessing addressing numerical precision issues
    cvals = rep(0, length(d))
    ds = d / rho
    idx.z = ds == 0
    cvals[idx.z] = 1
    bk = besselK(ds, kap)
    idx.big = bk == 0

    # evaluate on well-behaved inputs and return
    idx.eval = !idx.big & !idx.z
    cvals[idx.eval] = sc * (ds[idx.eval]^kap) * bk[idx.eval]
    return(cvals)
  }

  # if we got this far, input `pars$k` didn't match anything
  stop(paste('kernel name', pars$k, 'not recognized'))
}



#' Solve for range parameter given a correlation value
#'
#' @param cval numeric between 0 and 1, the correlation
#' @param pars a 1D kernel parameter list (with elements "k", "kp") recognized by `pkern_corr`
#' @param d positive numeric, the distance at which to evaluate the kernel
#' @param upper positive numeric, an upper limit for the search
#'
#' @return the estimated "rho" value producing semivariance `cval` at distance `d`
#' @export
#'
#' @examples
#' pars = pkern_corr('exp')
#' cval = 0.5
#' pkern_rho(cval, pars, d=1)
#' pkern_rho(cval, pars, d=10)
#' pkern_corr(utils::modifyList(pars, list(kp=pkern_rho(0.5, pars, d=10))), d=10)
#'
pkern_rho = function(cval, pars, d=1, upper=1e3*d)
{
  # anonymous function for inverting kernel with respect to range parameter
  fn = function(rho) abs( pkern_corr(utils::modifyList(pars, list(kp=c(rho, pars[['kp']][-1]))), d) - cval)
  return(stats::optimize(fn, lower=0, upper=upper)$minimum)
}



#' Suggest bounds for parameters of y and x component kernels for a separable grid model
#'
#' `pars` should be a kernel parameter list, or a list of two of them (named "y", "x").
#' The function adds vectors "kp" (initial values), "lower" and "upper" (bounds), wherever
#' they are missing. As a shortcut, `pars` can also be a character string (or vector of two)
#' naming the desired kernel(s), in which case the function sets a suggested initial value
#' (see `pkern_corr`).
#'
#' Bounds for shape parameters are hard-coded, whereas the bounds for the range parameter
#' are determined based on the grid layout (`gres`, `gdim`) and the arguments `cmin`, `cmax`.
#' These are, respectively, the minimum correlation between adjacent grid points, and the
#' maximum correlation between points at a distance of `nmax` cells apart. The corresponding
#' range parameter values are then estimated using `pkern_rho`, which does a grid search over
#' the permissible shape parameters (if there are any).
#'
#' Initial values for the range and shape parameters, when missing from `pars`, are set to
#' the midpoints of their permissible ranges.
#'
#' By default `cmin` is set to 0.05, so that the effective range of the kernel is bounded
#' below by the shortest interpoint distance. The default for `cmax` was selected by trial
#' and error.
#'
#' @param pars kernel parameter list (recognized by `pkern_corr`) or list of two of them
#' @param gres vector c(dy, dx) or positive numeric, the distances between adjacent grid lines
#' @param nmax positive integer, the number of grid cells between "distant" points
#' @param cmin positive numeric, minimum allowable correlation between adjacent points
#' @param cmax positive numeric, maximum allowable correlation between "distant" points
#'
#'
#' @return kernel parameter list containing "kp", "lower", "upper", or list of two such lists
#' @export
#'
#' @examples
#'
#' # suggested bounds depend on grid resolution and nmax
#' pars = 'mat'
#' pkern_bds(pars)
#' pkern_bds(pars, gres=100)
#' pkern_bds(pars, gres=100, nmax=100)
#'
#' # specify separable 2D kernels with a vector of two names
#' pars = c('exp', 'mat')
#' pkern_bds(pars)
#'
#' # incomplete parameter definition lists are filled in
#' pars = list(k='mat', upper=c(10,10))
#' pkern_bds(pars)
#'
pkern_bds = function(pars, gres=NA, nmax=1, cmin=0.05, cmax=0.5)
{
  # handle character input (kernel names)
  if( is.character(pars) )
  {
    # convert to expected list format for recursive call and handle unexpected input
    if( length(pars) == 1 ) return( pkern_bds(list(k=pars), gres=gres, nmax=nmax, cmin=cmin, cmax=cmax) )
    if( length(pars) != 2 ) stop('expected vector of two kernel names in pars')
    pars = lapply(stats::setNames(pars, c('y', 'x')), \(p) list(k=p))
    return( pkern_bds(pars, gres=gres, nmax=nmax, cmin=cmin, cmax=cmax) )
  }

  # single dimension case
  if( 'k' %in% names(pars) )
  {
    # do nothing and return if upper, lower and initial are already assigned
    if( all( c('kp', 'lower', 'upper') %in% names(pars) ) ) return(pars)

    # check for invalid input to gres, set default as needed
    if( length(gres) > 1 ) stop('gres should have length 1')
    if( is.na(gres) ) gres = 1

    # handle shape parameter, if there is one
    if( pars[['k']] %in% c('mat', 'gxp') )
    {
      # Whittle-Matern (note "mat" becomes indistinguishable from "gau" as shp->Inf)
      if( pars[['k']] == 'mat' )
      {
        # these bounds are suggestions based on experimentation
        bds.shp = c(1, 25)
        shp.ini = 10
        nm.shp = 'kap'
      }

      # gamma-exponential
      if( pars[['k']] == 'gxp' )
      {
        # max of 2 to ensure positive definite covariance matrix
        bds.shp = c(0.5, 2)
        shp.ini = 1
        nm.shp = 'p'
      }

      # set up a grid of test values and find rho for each one, then take min/max
      shp.test = seq(bds.shp[1], bds.shp[2], length.out = 100)
      rhomin = min( sapply(shp.test, \(p) pkern_rho(sqrt(cmin), utils::modifyList(pars, list(kp=c(NA, p))), gres)) )
      rhomax = max( sapply(shp.test, \(p) pkern_rho(sqrt(cmax), utils::modifyList(pars, list(kp=c(NA, p))), gres * nmax)) )

      # append bounds as needed
      if( is.null( pars[['lower']] ) ) pars[['lower']] = c(rhomin, bds.shp[1])
      if( is.null( pars[['upper']] ) ) pars[['upper']] = c(rhomax, bds.shp[2])

      # compute and append initial values as needed
      if( is.null( pars[['kp']] ) )
      {
        #shp.ini = mean( c(pars[['lower']][2], pars[['upper']][2]) )
        #rho.ini = pkern_rho(sqrt(cini), utils::modifyList(pars, list(kp=c(NA, shp.ini))), gres)
        rho.ini = (rhomin + rhomax)/2
        pars[['kp']] = c(rho.ini, shp.ini)
      }

    } else {

      # no shape parameter case is simple 1d optimization
      nm.shp = NULL
      rhomin = pkern_rho(sqrt(cmin), pars, gres)
      rhomax = pkern_rho(sqrt(cmax), pars, nmax * gres)

      # append bounds and initial values as needed
      if( is.null( pars[['lower']] ) ) pars[['lower']] = rhomin
      if( is.null( pars[['upper']] ) ) pars[['upper']] = rhomax
      #if( is.null( pars[['kp']] ) ) pars[['kp']] = pkern_rho(sqrt(cini), pars, gres)
      if( is.null( pars[['kp']] ) ) pars[['kp']] = (rhomin + rhomax)/2
    }

    # name the elements of "kp", "lower", "upper"
    nm = c('rho', nm.shp)
    pars[['kp']] = stats::setNames(pars[['kp']], nm)
    pars[['lower']] = stats::setNames(pars[['lower']], nm)
    pars[['upper']] = stats::setNames(pars[['upper']], nm)

    # order the output
    idx.first = match(c('k', 'kp', 'lower', 'upper'), names(pars))
    pars = pars[ c(idx.first, seq( length(pars) )[-idx.first]) ]

    # finished with 1D case
    return(pars)
  }

  # nmax should be a vector of length 2
  if( length(nmax)==1 ) nmax = rep(nmax, 2)

  # else assume pars is a list containing two parameter lists "y" and "x"
  if( !all( c('y', 'x') %in% names(pars) ) ) stop('expected parameter lists "y" and "x" in pars')

  # recursive calls to build x and y component bounds and return results in list
  pars.new = Map(\(p, s, n) pkern_bds(p, s, n, cmin, cmax), p=pars[c('y', 'x')], s=gres, n=nmax)
  return( utils::modifyList(pars, stats::setNames(pars.new, c('y', 'x')) ) )
}


#' Vectorize x and y component parameters
#'
#' when `kp` is NULL, the function concatenates the two ("y" and "x") component kernel parameter
#' sets into a single vector. When this vector is passed back as `kp`, the function does the
#' inverse, copying parameter values from a vector to the corresponding list entry in `pars`
#' and returning them in a list
#'
#' @param pars list of "y" and "x" kernel parameter lists (each containing vector "kp")
#' @param kp vector of combined kernel parameters (or NA)
#'
#' @return numeric vector of parameters, or `pars` modified appropriately
#' @export
#'
#' @examples
#' pars = pkern_corr(c('mat', 'exp'))
#' pkern_unpack(pars)
#' pkern_unpack(pars, 1:3)
pkern_unpack = function(pars, kp=NULL)
{
  # vectorization in order y, x
  if( is.null(kp) ) return( c(pars[['y']][['kp']], pars[['x']][['kp']]) )

  # determine number of parameters from each kernel
  npy = length(pars[['y']][['kp']])
  npx = length(pars[['x']][['kp']])

  # unpack values from numeric vector kp
  kpy = kp[ seq(npy) ]
  kpx = kp[ npy + seq(npx) ]

  # overwrite x and y kernel parameter lists
  pars[['y']] = utils::modifyList(pars[['y']], list(kp=kpy))
  pars[['x']] = utils::modifyList(pars[['x']], list(kp=kpx))
  return(pars)
}


#' Construct 1D stationary correlation matrices for regularly spaced data
#'
#' An effient implementation that uses symmetry and Toeplitz structure arising
#' from assumption of stationarity.
#'
#' @param pars list of kernel parameters "k" and "kp" (see `pkern_corr`)
#' @param n positive integer, the number of points on the 1D line
#' @param gres positive numeric, the distance between adjacent grid lines
#' @param i vector, a subset of `seq(n)` indicating rows to return
#' @param j vector, a subset of `seq(n)` indicating columns to return
#'
#'
#' @return the n x n correlation matrix, or its subset as specified in `i`, `j`
#' @export
#'
#' @examples
#' pars = pkern_corr('exp')
#' pkern_corrmat(pars, n=10)
#' pkern_corrmat(pars, n=10, i=2:4, j=2:4)
#' pkern_corrmat(pars, n=3)
#' pkern_corrmat(pars, n=3, gres=2)
pkern_corrmat = function(pars, n, gres=1, i=seq(n), j=seq(n))
{
  # compute the set of distances over which we need to evaluate kernel
  du = gres * ( seq(n) - 1 )

  # compute kernel values for these distances
  dcorr = pkern_corr(pars, du)

  # build large vector to shift through in building the (Toeplitz) output matrix
  bigvec = c(dcorr[n:2], dcorr)

  # build and return the matrix
  return( Matrix::Matrix(sapply(j, function(x) bigvec[ (n-x) + i ])) )
}


