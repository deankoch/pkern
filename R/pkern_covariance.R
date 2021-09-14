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
#' If any element of `d` is NA, the function ignores `d` and instead returns a list
#' containing the parameter names, along with suggested initial values and bounds to
#' use in optimization.
#'
#' If instead of a list, `pars` is one of the kernel name strings, the function returns
#' an example of the parameter list format for the kernel, with suggested initial values
#' (and `d` is ignored)
#'
#' @param pars list with elements "k", the kernel name, and "p" the parameter vector
#' @param d numeric vector of length n, the distances to evaluate
#'
#' @return length-n vector, or a matrix of parameter bounds, or a parameter list (see details)
#' @export
#'
#' @examples
#' pars.mat = pkern_corr('mat')
#' pars.mat
#' pkern_corr(pars.mat)
#' pkern_corr(pars.mat, d=1:10)
pkern_corr = function(pars, d=NA)
{
  # if `d` is (or contains) NA, return suggested bounds for supplied pars
  makebounds = anyNA(d)
  bds.rho = c(min=1e-2, ini=1, max=1e6)
  bds.p = c(min=1e-2, ini=1, max=2)
  bds.kap = c(min=1e-2, ini=1, max=40)
  bds.sph.rho = c(min=1e-2, ini=5, max=1e6)
  bds.exp.rho = c(min=1e-2, ini=2, max=1e6)

  # handle character input to pars as request for initial values
  if( is.character(pars) )
  {
    pars = list(k=pars, kp=NULL)
    kp.suggested = pkern_corr(pars)[,'ini']
    return( modifyList(pars, list(kp=kp.suggested)) )
  }

  # handle invalid pars
  if( !all( c('k', 'kp') %in% names(pars) ) ) stop('pars must be list with elements "k" and "kp"')

  # exponential
  if(pars[['k']] == 'exp')
  {
    # return suggested bounds and initial value if requested
    if(makebounds) return(rbind(rho=bds.exp.rho))

    # unpack names (or assume rho is first element of pars$kp)
    kp.idx = match('rho', names(pars$kp))
    if( is.na(kp.idx) ) kp.idx = 1
    rho = pars$kp[kp.idx]

    # new parameter list for recursive call
    pars = list(k='gxp', kp=c(rho=rho, p=1L))
    return( pkern_corr(pars, abs(d)) )
  }

  # gaussian/stable
  if(pars[['k']] == 'gau')
  {
    # return suggested bounds and initial value if requested
    if(makebounds) return(rbind(rho=bds.rho))

    # unpack names (or assume rho is first element of pars$kp)
    kp.idx = match('rho', names(pars$kp))
    if( is.na(kp.idx) ) kp.idx = 1
    rho = pars$kp[kp.idx]

    # new parameter list for recursive call
    pars = list(k='gxp', kp=c(rho=rho, p=2L))
    return( pkern_corr(pars, d) )
  }

  # spherical
  if(pars[['k']] == 'sph')
  {
    # return suggested bounds and initial value if requested
    if(makebounds) return(rbind(rho=bds.sph.rho))

    # unpack names (or assume rho is first element of pars$kp)
    kp.idx = match('rho', names(pars$kp))
    if( is.na(kp.idx) ) kp.idx = 1

    # assign parameter and identify truncations
    rho = pars$kp[kp.idx]
    ds = d / rho
    cvals = rep(0, length(d))
    idx.nz = ds < 1

    # evaluate on non-truncated distances and return
    cvals[idx.nz] = ( 1 - (3/2)*ds[idx.nz] + (1/2)*(ds[idx.nz]^3) )
    return( cvals )
  }

  # gamma-exponential
  if(pars[['k']] == 'gxp')
  {
    # return suggested bounds and initial value if requested
    if(makebounds) return(rbind(rho=bds.rho, p=bds.p))

    # unpack names (or assume default order rho, p)
    kp.idx = match(c('rho', 'p'), names(pars$kp))
    if( anyNA(kp.idx) ) kp.idx = c(1,2)
    if( anyNA(pars$kp[kp.idx]) ) stop(paste('require', length(kp.idx), 'parameters in pars$kp'))

    # assign parameters and evaluate
    rho = pars$kp[kp.idx[1]]
    p = pars$kp[kp.idx[2]]
    return( exp( -( ( d / rho )^p  ) ) )
  }

  # Matern
  if(pars[['k']] == 'mat')
  {
    # return suggested bounds and initial value if requested
    if(makebounds) return(rbind(rho=bds.rho, kap=bds.kap))

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
#' @param pars a kernel parameter list recognized by `pkern_corr`
#' @param d positive numeric, the distance at which to evaluate the kernel
#' @param upper positive numeric, an upper limit for the search
#'
#' @return
#' @export
#'
#' @examples
pkern_rho = function(cval, pars, d, upper=1e3*d)
{
  # anonymous function for inverting kernel with respect to range parameter
  fn = function(rho) abs( pkern_corr(modifyList( pars, list(kp=c(rho, pars[['kp']][-1]))), d) - cval)
  return(optimize(fn, lower=0, upper=upper)$minimum)
}


#' Suggest bounds for parameters of x and y component kernels for a grid model
#'
#' This function takes a kernel name or parameter list (`pars`) and adds vectors "lower"
#' and "upper", bounds for the kernel parameters, whenever they are missing. If the input
#' has no parameter values assigned (in `pars$kp`), the function sets them to the midpoint
#' (mean) of the respective bounds.
#'
#' Bounds for shape parameters are hard-coded, whereas the bounds for the range parameter
#' are determined based on the grid resolution (`ds`, the distance between grid lines in
#' x and y directions) and a user supplied correlation range; `cmin` is the minimum allowable
#' correlation between adjacent grid points, and `cmax` is the maximum. By default `cmin` is
#' set to 0.05, so that the effective range of the kernel is bounded
#' below by the shortest interpoint distance.
#'
#'
#' @param pars kernel parameter list (recognized by `pkern_corr`) or list of two of them
#' @param ds vector c(dx, dy) or positive numeric, the resolution of the sampled grid
#' @param cmin positive numeric, minimum allowable correlation between adjacent grid cells
#' @param cmax positive numeric, maximum allowable correlation between adjacent grid cells
#'
#'
#' @return kernel parameter list containing "kp", "lower", "upper", or list of two such lists
#' @export
#'
#' @examples
pkern_bds = function(pars, ds=NA, cmin=0.05, cmax=1-cmin)
{
  # handle character input (kernel names)
  if( is.character(pars) ) pars = list(k=pars, kp=NULL)

  # single dimension case
  if( all( c('k', 'kp') %in% names(pars) ) )
  {
    # do nothing and return if upper and lower are already assigned
    if( all( c('kp', 'lower', 'upper') %in% names(pars) ) ) return(pars)

    # check for invalid input to ds, set default as needed
    if( length(ds) > 1 ) stop('ds should have length 1')
    if( is.na(ds) ) ds = 1

    # handle shape parameter, if there is one
    if( pars[['k']] %in% c('mat', 'gxp') )
    {
      # set for numerical stability (note "mat" becomes indistinguishable from "gau" as shp->Inf)
      if( pars[['k']] == 'mat' ) bds.shp = c(0.5, 40)

      # max of 2 to ensure positive definite covariance matrix
      if( pars[['k']] == 'gxp' ) bds.shp = c(0.5, 2)

      # set up a grid of test values and find rho for each one, then take min/max
      shp.test = seq(bds.shp[1], bds.shp[2], length.out = 100)
      rhomin = min( sapply(shp.test, \(p) pkern_rho(cmin, modifyList(pars, list(kp=c(NA, p))), ds)) )
      rhomax = max( sapply(shp.test, \(p) pkern_rho(cmax, modifyList(pars, list(kp=c(NA, p))), ds)) )

      # append bounds and initial values as needed
      if( is.null( pars[['lower']] ) ) pars[['lower']] = c(rhomin, bds.shp[1])
      if( is.null( pars[['upper']] ) ) pars[['upper']] = c(rhomax, bds.shp[2])
      if( is.null( pars[['kp']] ) ) pars[['kp']] = ( pars[['lower']] + pars[['upper']] ) / 2

    } else {

      # no shape parameter case is simple 1d optimization
      rhomin = pkern_rho(cmin, pars, ds)
      rhomax = pkern_rho(cmax, pars, ds)

      # append bounds and initial values as needed
      if( is.null( pars[['lower']] ) ) pars[['lower']] = rhomin
      if( is.null( pars[['upper']] ) ) pars[['upper']] = rhomax
      if( is.null( pars[['kp']] ) ) pars[['kp']] = ( rhomin + rhomax ) / 2

    }

    # finished with 1D case
    return(pars)
  }

  # else assume pars is a list containing two parameter lists "x" and "y"
  if( !all( c('x', 'y') %in% names(pars) ) ) stop('expected parameter lists "x" and "y" in pars')

  # recursive calls to build x and y component bounds
  pars.out = modifyList(pars, list(x = pkern_bds(pars[['x']], ds[1], cmin=cmin, cmax=cmax),
                                   y = pkern_bds(pars[['y']], ds[2], cmin=cmin, cmax=cmax)) )

  # return results in list
  return( pars.out )
}




#' Vectorize x and y component parameters
#'
#' when `p` is NULL, the function concatenates the two (x and y) component kernel parameter
#' sets into a single vector. When this vector is passed back as `p`, the function does the
#' inverse, copying parameter values from a vector to the corresponding list entry in `xpars`
#' or `ypars`, and returning them in a list
#'
#' If `p` is set to "min", "ini", or "max", the function returns a vector of suggested
#' minima, inital values, or maxima for the parameters.
#'
#' @param xpars list of x component covariance paramaters (in the form accepted by `pkern_corr`)
#' @param ypars list of y component covariance paramaters (in the form accepted by `pkern_corr`)
#' @param p (optional) vector of combined kernel parameters, or one of "min", "ini", "max"
#'
#' @return numeric vector of parameters, or `list(xpars, ypars)` modified appropriately
#' @export
#'
#' @examples
#' pars.mat = pkern_corr('mat')
#' pars.exp = pkern_corr('exp')
#' p = pkern_unpack(pars.exp, pars.mat)
#' p
#' p[2:3] = 1
#' pkern_unpack(pars.exp, pars.mat, p)
#' pkern_unpack(pars.exp, pars.mat, 'min')
pkern_unpack = function(xpars, ypars, p=NULL)
{
  # vectorization in order x, y
  if( is.null(p) ) return( c(xpars$kp, ypars$kp) )

  # determine number of parameters from each kernel
  npx = length(xpars$kp)
  npy = length(ypars$kp)

  # differentiate special parameter set requests
  if( any( p %in% c('min', 'ini', 'max') ) )
  {
    # fetch the hard-coded parameter sets from another function
    kp.suggested = rbind( pkern_corr(xpars), pkern_corr(Nypars) )
    return( kp.suggested[, p] )

  } else {

    # unpack values from numeric vector p
    kpx = p[seq(npx)]
    kpy = p[npx + seq(npy)]
  }

  # generate x and y kernel parameter lists
  xpars = modifyList(xpars, list(kp=kpx))
  ypars = modifyList(ypars, list(kp=kpy))
  return( list(x=xpars, y=ypars) )
}


#' Construct 1D stationary correlation matrices for regularly spaced data
#'
#' An effient implementation that uses symmetry and Toeplitz structure arising
#' from assumption of stationarity.
#'
#' @param pars list or parameters (in the form accepted by `pkern_corr`)
#' @param n positive integer, the number of points on the 1D line
#' @param ds positive numeric, the distance between adjacent grid lines
#' @param i vector, a subset of `seq(n)` indicating rows to return
#' @param j vector, a subset of `seq(n)` indicating columns to return
#'
#' @return the n x n correlation matrix, or its subset as specified in `i`, `j`
#' @export
#'
#' @examples
#' pars = pkern_corr('exp')
#' pkern_corrmat(10, pars)
#' pkern_corrmat(10, pars, c(1,2), c(3:4))
pkern_corrmat = function(pars, n, ds=1, i=seq(n), j=seq(n))
{
  # compute the set of distances over which we need to evaluate kernel
  du = ds * ( seq(n) - 1 )

  # compute kernel values for these distances
  dcorr = pkern_corr(pars, du)

  # build large vector to shift through in building the (Toeplitz) output matrix
  bigvec = c(dcorr[n:2], dcorr)

  # build and return the matrix
  return(sapply(j, function(x) bigvec[ (n-x) + i ]))
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
#' pars = list(x=pkern_corr('mat'), y=pkern_corr('exp'))
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


