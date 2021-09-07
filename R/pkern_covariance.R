#
# pkern_raster.R
# Dean Koch, Sep 2021
# Functions for handling separable covariance functions
#


#' Stationary 1D correlation kernels
#'
#' Computes stationary correlation function values for (nonegative) 1-dimensional
#' distances `d`. Parameter list entry `pars$kp` supplies the kernel parameter(s),
#' where `pars$k` is one of the following kernel names
#'
#' "exp": exponential (special case of "gex" with p=1 and `d<-abs(d)`). `kp = rho`
#' "gau": gaussian/stable (special case of "gex" with p=2). `kp = rho`
#' "sph": spherical (AKA stable/Gaussian for p=2). `kp = rho`
#' "gex": gamma-exponential. `kp = c(rho, p)`
#' "mat": Whittle-Matern ((Handcock and Wallis parameterization). `kp = c(rho, kap)`
#'
#' @param d numeric vector of length n, the distances to evaluate
#' @param pars list with elements "k", the kernel name, and "p" the parameter vector
#'
#' @return the length-n vector of correlations
#' @export
#'
#' @examples
#' pkern_corr(d=1:10, pars=list(k='exp', kp=5))
#' pkern_corr(d=1:10, pars=list(k='mat', kp=c(5, 3)))
pkern_corr = function(d, pars)
{
  # if `d` is (or contains) NA, return suggested bounds for supplied pars
  #makebounds = anyNA(d)

  # handle invalid pars
  if( !all( c('k', 'kp') %in% names(pars) ) ) stop('pars must be list with elements "k" and "kp"')

  # exponential
  if(pars$k == 'exp')
  {
    # unpack names (or assume rho is first element of pars$kp)
    kp.idx = match('rho', names(pars$kp))
    if( is.na(kp.idx) ) kp.idx = 1
    rho = pars$kp[kp.idx]

    # new parameter list for recursive call
    pars = list(k='gxp', kp=c(rho=rho, p=1))
    return( pkern_corr(abs(d), pars) )
  }

  # gaussian/stable
  if(pars$k == 'gau')
  {
    # unpack names (or assume rho is first element of pars$kp)
    kp.idx = match('rho', names(pars$kp))
    if( is.na(kp.idx) ) kp.idx = 1
    rho = pars$kp[kp.idx]

    # new parameter list for recursive call
    pars = list(k='gxp', kp=c(rho=rho, p=2))
    return( pkern_corr(d, pars) )
  }

  # spherical
  if(pars$k == 'sph')
  {
    # unpack names (or assume rho is first element of pars$kp)
    kp.idx = match('rho', names(pars$kp))
    if( is.na(kp.idx) ) kp.idx = 1

    # assign parameter and identify truncations
    rho = pars$kp[kp.idx]
    ds = d/rho
    cvals = rep(0, length(d))
    idx.nz = ds < 1

    # evaluate on non-truncated distances and return
    cvals[idx.nz] = ( 1 - (3/2)*ds[idx.nz] + (1/2)*(ds[idx.nz]^3) )
    return( cvals )
  }

  # gamma-exponential
  if(pars$k == 'gxp')
  {
    # unpack names (or assume default order rho, p)
    kp.idx = match(c('rho', 'p'), names(pars$kp))
    if( anyNA(kp.idx) ) kp.idx = c(1,2)
    if( anyNA(pars$kp[kp.idx]) ) stop(paste('require', length(kp.idx), 'parameters in pars$kp'))

    # assign parameters and evaluate
    rho = pars$kp[kp.idx[1]]
    p = pars$kp[kp.idx[2]]
    return( exp( -( (d^p) / rho ) ) )
  }

  # Matern
  if(pars$k == 'mat')
  {
    # unpack names (or assume default order rho, p)
    kp.idx = match(c('rho', 'kap'), names(pars$kp))
    if( anyNA(kp.idx) ) kp.idx = c(1,2)
    if( anyNA(pars$kp[kp.idx]) ) stop(paste('require', length(kp.idx), 'parameters in pars$kp'))

    # assign parameters and compute scaling constant
    kap = pars$kp[kp.idx[2]]
    rho = pars$kp[kp.idx[1]] / ( 2 * sqrt(kap))
    sc = ( (2^( 1 - kap ) ) / gamma(kap) )

    # preprocessing to fix numerical precision issues
    cvals = rep(0, length(d))
    ds = d/rho
    idx.z = ds == 0
    cvals[idx.z] = 1
    bk = besselK(ds, kap)
    idx.big = bk == 0

    # evaluate on well-behaved inputs and return
    idx.eval = !idx.big & !idx.z
    cvals[idx.eval] = sc * (ds[idx.eval]^kap) * bk[idx.eval]
    return(cvals)
  }
}


#' Vectorize x and y component parameters
#'
#' @param x list of x component
#' @param y
#' @param p
#' @param sig
#'
#' @return
#' @export
#'
#' @examples
pkern_unpack = function(x, y, p=NULL)
{



}


#' Construct 1D stationary correlation matrices for regularly spaced data
#'
#' An effient implementation that uses symmetry and Toeplitz structure arising
#' from assumption of stationarity.
#'
#' @param n positive integer, the number of points on the 1D line
#' @param pars list or parameters (in the form accepted by `pkern_corr`)
#' @param ds positive numeric, the distance between adjacent grid lines
#' @param i vector, a subset of `seq(n)` indicating rows to return
#' @param j vector, a subset of `seq(n)` indicating columns to return
#'
#' @return the n x n correlation matrix, or its subset as specified in `i`, `j`
#' @export
#'
#' @examples
#' pars = list(k='exp', kp=5)
#' pkern_corrmat(10, pars)
#' pkern_corrmat(10, pars, c(1,2), c(3:4))
pkern_corrmat = function(n, pars, ds=1, i=seq(n), j=seq(n))
{
  # compute the set of distances over which we need to evaluate kernel
  du = ds * ( seq(n) - 1 )

  # compute kernel values for these distances
  dcorr = pkern_corr(du, pars)

  # build large vector to shift through in building the (Toeplitz) output matrix
  bigvec = c(dcorr[n:2], dcorr)

  # build and return the matrix
  return(sapply(i, function(x) bigvec[ (n-x) + j ]))
}



