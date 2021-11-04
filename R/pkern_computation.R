#
# pkern_raster.R
# Dean Koch, Oct 2021
# Functions for computationally intensive linear algebra
#


#' efficient quadratic form computer: `t(x) %*% M %*% y`
#'
#' @param M numeric matrix
#' @param x numeric vector
#' @param y numeric vector
#'
#' @return TODO
#' @export
#'
#' @examples
#' # TODO
pkern_qf = function(M, x, y)
{
  # Efficient computation of the quadratic form `t(x) %*% M %*% Y`
  # where M is a numeric matrix with conforming dimensions
  #
  # based on the trick used in `emulator::quad.form`

  # coerce everything to matrices with appropriate dimensions
  ninner = dim(M)[1]
  nouter = dim(M)[2]
  if( is.vector(x) ) x = matrix(x, ninner)
  if( is.vector(y) ) y = matrix(y, nouter)

  # verify that dimensions conform (Rfast::Crossprod will crash R otherwise!)
  conforms.inner = ninner == dim(x)[1]
  conforms.right = nouter == dim(y)[1]
  if( !conforms.inner | !conforms.right ) stop('dimensions do not conform')

  # run computation
  return( as.vector(Matrix::crossprod(Matrix::crossprod(M, x), y)) )
}



#' Left-multiplication with a Kronecker product: evaluates (X (x) Y) * vector
#'
#' Evaluates a Kronecker product times a vector by the standard trick of representing
#' the vector as a matrix M and evaluating `Y %*% M %*% t(X)`.
#'
#' For computational efficiency, the function computes the result by a composition
#' of cross-products with the transposed input component matrices. The argument
#' `trans=TRUE` indicates that `t(X)` and `t(Y)` have been passed to `X` and `Y`,
#' so that the transpose step (which is slow for large matrices) can be skipped.
#'
#' @param X matrix or Matrix object, the x-component matrix
#' @param Y matrix or Matrix object, the y-component matrix
#' @param z numeric vector, the vector to left-multiply
#' @param trans logical, indicating the inputs `X` and `Y` are transposed
#'
#' @return TODO
#' @export
#'
#' @examples
#' # TODO
pkern_kprod = function(X, Y, z, trans=FALSE)
{
  # transposed input matrices are processed faster
  if(trans) return( pkern_qf(Matrix::Matrix(as.vector(z), nrow(Y), nrow(X)), Y, X) )

  # here we have to take transpose of X and Y
  return( pkern_qf(Matrix::Matrix(z, ncol(Y), ncol(X)), t(Y), t(X)) )
}


#' Likelihood function for sampled points in a subgrid
#'
#' `zobs` can either be a matrix (in which case `sgdim`, if supplied, must match `dim(zobs)`)
#' or a vector of length `prod(sgdim)`.
#'
#' @param zobs numeric vector or matrix of subgrid data (with NAs at unsampled points)
#' @param sgdim integer vector, the size of the subgrid (ny, nx), required if `zobs` is a vector
#' @param pars kernel parameter list, see `pkern_corr`
#' @param pc logical or list of matrices, precomputed objects (see details)
#'
#' @return
#' @export
#'
#' @examples
#' # TODO
pkern_LL = function(pars, zobs, sgdim, pc=FALSE)
{
  # extract covariance function parameters
  yx.nm = c('y', 'x')
  cpars = pars[yx.nm]
  pvar = ifelse( is.null(pars[['v']]), 1, sqrt(pars[['v']]))

  # set default subgrid resolution (distances between grid lines) as needed
  dsub = pars[['ds']]
  if( is.null(dsub) ) dsub = c(1,1)

  # set default nugget as needed
  pnug = ifelse( is.null(pars[['nug']]), 0, pars[['nug']])

  # subgrids with all points sampled are handled by more efficient means
  sobs = which( !is.na(zobs) )
  nobs = length(sobs)
  is.complete = nobs == length(zobs)

  # handle requests for precomputed objects
  if( is.logical(pc) )
  {
    # flag to continue and compute likelihood or else return precomputed objects
    continue = !pc

    # initialize precomputed object list with indices of sampled subgrid points
    pc = list( sobs = sobs )

    # build component marginal correlation matrices for the subgrid (at dsub resolution)
    pc[['vs']] = Map( \(p, n, d) { pvar * pkern_corrmat(p, n, d) },
                      p=cpars, n=sgdim, d=dsub) |> stats::setNames(nm=yx.nm)

    # when all of the subgrid is sampled: do eigendecompositions of component matrices
    if( length(sobs) == 0 ) { pc[['ed']] = lapply(vs, \(v) eigen(v, symmetric=TRUE)) } else {

      # missing data case: eigendecomposition on submatrix of full variance kronecker product
      pc[['ed']] = kronecker(pc[['vs']][['x']], pc[['vs']][['y']])[sobs, sobs] |>
        eigen(symmetric=TRUE)
    }

    # return precomputed object list if requested
    if( !continue ) return(pc)

  } else {

    # precomputed objects supplied: check for invalid input
    if( !is.list(pc) ) stop('pc must be a list')

    # check that we have same pattern of NAs in precomputed objects
    sobs.in = pc[['sobs']]
    if( !identical(sobs, sobs.in) ) stop('missing data in zobs must match that used to create pc')
  }

  # case of subgrid with no unsampled points
  if( is.complete )
  {
    # compute log-determinant as sum of log-eigenvalues (after adding nugget variance)
    ldet.yx = Map(\(ed) n * sum( log( ed[['values']] + pnug ) ), ed=pc[['ed']], n=sgdim)
    ldet = do.call(sum, ldet.yx)

    # TODO: finish this
    # compute quadratic form of observations with inverse covariance matrix
    #pkern_kprod(pc[['y']][['vectors']], pc[['ed']][['x']][['vectors']], zobs)
    lqf = NA
    warning('not yet implemented')

  } else {

    # missing data case

    # compute eigenvalues after adding nugget variance
    pwv = pc[['ed']][['values']] + pnug

    # numerical problems are indicated here by -Inf log-likelihood
    if( !all(pwv > 0) ) return(-Inf)

    # log-determinant is sum of log-eigenvalues
    ldet = sum( log(pwv) )

    # compute quadratic form of observations with inverse covariance matrix
    lqf = sum( ( ( t(pc[['ed']][['vectors']]) %*% zobs[sobs] ) / sqrt(pwv)  )^2 )
  }

  #  log-likelihood
  return( - ldet - lqf )

}


#' Fit a spatial covariance model by numerical maximum likelihood
#'
#' @param zobs numeric vector or matrix of subgrid data (with NAs at unsampled points)
#' @param gsnap list of vectors, output from `pkern_snap`
#' @param ypars character or list, recognizable (as `pars`) by `pkern_corr`
#' @param xpars character or list, recognizable (as `pars`) by `pkern_corr`
#' @param v length-3 numeric vector, c('min', 'ini', 'max'), bounds and initial value for variance
#' @param nug length-3 numeric vector, c('min', 'ini', 'max'), bounds and initial value for nugget
#' @param add positive integer, add additional initial parameter sets are tested (see details)
#'
#' @return
#' @export
#'
#' @examples
pkern_fit = function(zobs, gsnap, ypars='gau', xpars=ypars, v=NULL, nug=NULL, add=0, control=list())
{
  # set defaults for variance (first parameter) and nugget (second parameter)
  vini = var(zobs, na.rm=TRUE)
  if( length(v) == 0 ) v = vini
  if( length(v) != 3 ) v = c(min=1e-9, ini=v, max=4*vini)
  if( length(nug) == 0 ) nug = v[2]/2
  if( length(nug) != 3 ) nug = c(min=1e-3, ini=nug, max=v[3])

  # set up defaults and bounds for kernel parameters (2-4 additional parameters)
  ypars = pkern_bds(ypars, gsnap[['ds']][1])
  xpars = pkern_bds(xpars, gsnap[['ds']][2])

  # vectorized bounds for all covariance parameters
  plower = c(v[1], nug[1], xpars[['lower']], ypars[['lower']])
  pinitial = c( v[2], nug[2], xpars[['kp']], ypars[['kp']] )
  pupper = c(v[3], nug[3], xpars[['upper']], ypars[['upper']])

  # define anonymous objective function for optimizer
  fn = function(pvec, p, zobs, sgdim)
  {
    # pvec, numeric vector of parameters
    # p, list of "x", "y" kernel parameter lists associated with pvec
    # zobs, vectorized subgrid data (with NAs at unsampled points)
    # sgdim, integer vector, the size of the subgrid (ny, nx)

    # bundle kernel parameters into list
    ptest = pkern_unpack(p, pvec[-(1:2)])
    ptest[['v']] = pvec[1]
    ptest[['nug']] = pvec[2]

    # compute negative log-likelihood and handle invalid output
    nll = -pkern_LL(ptest, zobs, sgdim)
    if( is.na(nll) | is.infinite(nll) ) nll = 2^.Machine$double.digits
    return( nll )
  }

  # make a list of initial values to test
  np = length(pinitial)
  list.ini = list(pinitial)
  if(add>0) list.ini = c(list.ini, lapply(seq(add), \(ini) plower+stats::runif(np)*(pupper-plower)) )

  # run the optimizer for each one
  ptemp = list(x=xpars, y=ypars, ds=gsnap[['ds']])
  list.optim = lapply(list.ini, \(ini) stats::optim(par=ini,
                                                    f=fn,
                                                    method='L-BFGS-B',
                                                    lower=plower,
                                                    upper=pupper,
                                                    p=ptemp,
                                                    zobs=zobs,
                                                    sgdim=gsnap[['sgdim']],
                                                    control=control))


  # select the best fit
  idx.best = which.min( sapply(list.optim, \(result) result[['value']]) )
  result.optim = list.optim[[idx.best]]

  # unpack fitted parameters and return in a list
  vfit = stats::setNames(result.optim$par[1], 'fitted')
  nugfit = stats::setNames(result.optim$par[2], 'fitted')
  pfit = pkern_unpack(ptemp, result.optim$par[-(1:2)])
  return( list(y=pfit[['y']], x=pfit[['x']], v=vfit, nug=nugfit, ds=gsnap[['ds']]) )



}



#' Compute important indexing vectors relating to a subgrid within a larger grid
#'
#' Regular subgrids have computationally convenient properties that we can exploit
#' when predicting onto a supergrid
#'
#' TODO: flesh out documentation
#' (cross) covariance matrices
#'
#' The length of `obs` should equal the product of the lengths of the grid line vectors in
#' `gxy`, and its elements should be in column-vectorized order. Only the locations of the NAs
#' (and not the other contents) of `obs` matter.
#'
#' @param gdim integer vector c(ni, nj) the size of the outer grid
#' @param map list of integer vectors "i" and "j", row and column indices of the subgrid
#' @param pars (optional) list of "x" and "y" kernel parameters
#' @param zobs (optional) vector, containing NAs where data are missing in the subgrid
#' @param makev logical, indicates to compute variance matrix components for points off subgrid
#'
#' @return a large list of vectors and matrices
#' @export
#'
#' @examples
#' # TODO
pkern_precompute = function(gdim, map, pars=NULL, zobs=NULL, makev=FALSE)
{
  # ordering is i/y first, j/x second, by default assume all of subgrid is sampled
  yx.nm = c('y', 'x')
  cpars = pars[yx.nm]
  pvar = ifelse( is.null(pars[['v']]), 1, sqrt(pars[['v']]))
  nmap = sapply(map, length)
  if( is.null(zobs) ) zobs = seq( prod(nmap) )
  if( !all( c('i', 'j') %in% names(map) ) ) map = stats::setNames(map, c('i', 'j'))

  # find subgrid dimensions
  mapu = lapply(map, \(m) sort(unique(m)))
  dsg = sapply(mapu, \(s) sort(diff(s))[1])

  # find indices of grid lines on and off subgrid
  mapsg = Map(\(m, sep) seq(min(m), max(m), by=sep), m=mapu, sep=dsg)
  mapc = mapply(\(s, g) seq(g)[-s], g=gdim, s=mapsg)

  # index of subgrid points in full grid, and all points on one of its grid lines
  idx.s = pkern_idx_sg(gdim, mapsg)
  idx.sj = pkern_idx_sg(gdim, utils::modifyList(mapsg, list(i=NULL)))
  idx.si = pkern_idx_sg(gdim, utils::modifyList(mapsg, list(j=NULL)))

  # partition non-subgrid points into three sets depending on overlap with grid lines
  idx.o = seq( prod(gdim) )[ -c(idx.s, idx.sj, idx.si) ]
  idx.oj = idx.sj[ !( idx.sj %in% idx.s ) ]
  idx.oi = idx.si[ !( idx.si %in% idx.s ) ]

  # indices of sampled and unsampled points within the subgrid
  idx.sobs = which( !is.na(zobs) )
  idx.smiss = seq( length(zobs) )[-idx.sobs]

  # set default subgrid resolution (distances between grid lines) as needed
  dsub = pars[['ds']]
  if( is.null(dsub) ) dsub = c(1,1)

  # calculate full grid resolution (distances between grid lines)
  dfull = dsub / dsg

  # build component marginal correlation matrices for the subgrid (at dsub resolution)
  vs = Map( \(p, n, d) { pvar * pkern_corrmat(p, n, d) },
            p=cpars, n=nmap, d=dsub) |> stats::setNames(nm=yx.nm)

  # build component marginal correlation matrices for points off subgrid (slow!)
  vo = NULL
  if(makev) vo = Map(\(p, n, d, i, j) { pvar * pkern_corrmat(p, n, d, i, j) },
                     p=cpars, n=gdim, d=dfull, i=mapc, j=mapc) |> stats::setNames(nm=yx.nm)

  # build component cross-correlation matrices for points on vs off subgrid lines
  vso = Map( \(p, n, d, i, j) { pvar * pkern_corrmat(p, n, d, i, j) },
             p=cpars, n=gdim, d=dfull, i=map, j=mapc) |> stats::setNames(nm=yx.nm)

  # use faster method when all of the subgrid is sampled
  if( !anyNA(zobs) )
  {
    # eigendecompositions of component matrices
    ed = lapply(vs, \(v) eigen(v, symmetric=TRUE))
    vsc = NULL

  } else {

    # missing data case requires first evaluating the variance kronecker product
    vs.full = kronecker(vs[['x']], vs[['y']])

    # eigendecomposition of full covariance matrix of subgrid
    ed = eigen(vs.full[idx.sobs, idx.sobs], symmetric=TRUE)

    # we also need the (within subgrid) cross-covariance between sampled and unsampled points
    vsc = vs.full[idx.smiss, idx.sobs]

  }

  # return everything in list
  idx.list = list(s=idx.s, o=idx.o, oj=idx.oj, oi=idx.oi, sobs=idx.sobs)
  v.list = list(vs=vs, vo=vo, vso=vso, vsc=vsc, ed=ed)
  return(c(idx.list, v.list))
}


#' Compute conditional mean of points on a grid given a sample taken over a subgrid
#'
#' Given a separable spatial covariance model (`pars`), and a sample (`zobs`) on a
#' a regular subgrid (`gli`) of the grid of dimensions `gdim`, the function computes the
#' expected value of all n (`prod(gdim)`) points in the full grid.
#'
#' For this type of problem the full n-dimensional covariance matrix (V) can be represented
#' using kronecker products of dimension `gdim[1]^2` and `gdim[2]^2` (eg. as described in
#' Gilboa et al., 2015; and Koch et al., 2020), dramatically speeding up calculations and
#' reducing memory demands.
#'
#' `zobs` should be supplied in column-vectorized order, and its length must equal the
#' product of the lengths of the vectors in `gli` (ie the number of points in the subgrid).
#' Unsampled points in the subgrid should be indicated by NAs in `zobs`.
#'
#' The subgrid on which `zobs` is located must be specified in the indexing vectors `gli$y`
#' and `gli$x`, which gli to subsets of `seq(gdim[1])` and `seq(gdim[2])`.
#'
#' Nugget effects are supported by adding adding `pars$nug` to the diagonal of the eigenvalue
#' matrix in the eigendecomposition of V. Note that a small nugget effect will often solve
#' computational issues arising from numerically singular covariance matrices (a common problem
#' with large spatial covariance models).
#'
#' @param zobs numeric vector of observed data, possibly with NAs
#' @param gdim integer vector (ny, nx), the size of the full grid
#' @param pars list of kernel parameter lists "y" and "x", nugget "nug", and distance scaling "ds"
#' @param gli list of two integer vectors, indexing the y and x grid lines of the regular subgrid
#' @param pc either logical (default FALSE), or a list of precomputed objects (see details)
#'
#' @return TODO
#' @export
#'
#' @examples
#' # TODO
pkern_cmean = function(zobs, gdim, pars, gli=NULL, pc=FALSE)
{
  # when gidx is not supplied, it should be found in list input gdim
  if(is.null(gli))
  {
    err.gdim = 'either gli or gdim$gli must be supplied'
    err.nm = 'list gdim should have elements "gli" and "gdim"'
    if( !is.list(gdim) ) stop(err.gdim)
    if( !all(c('gdim', 'gli') %in% names(gdim)) ) stop(err.nm)
    gli = gdim[['gli']]
    gdim = gdim[['gdim']]
  }

  # the number of points in the sampled subgrid, and index of missing elements
  nsg = length(zobs)
  sobs = which(!is.na(zobs))
  nobs = length(sobs)

  # set default 0 nugget variance as needed
  pars[['nug']] = ifelse(is.null(pars[['nug']]), 0, pars[['nug']])

  # handle requests for precomputed objects
  if( is.logical(pc) )
  {
    # prepare matrices and indexing vectors
    pcout = pkern_precompute(gdim, gli, pars, zobs=zobs, makev=FALSE)
    if( pc ) { return(pcout) } else { return(pkern_cmean(zobs, gdim, pars, gli, pcout)) }

  } else {

    # precomputed objects supplied: check for invalid input
    if( !is.list(pc) ) stop('pc must be a list')

    # check that we have same pattern of NAs in precomputed objects
    sobs.in = pc[['sobs']]
    if( !identical(sobs, sobs.in) ) stop('missing data in zobs must match that used to create pc')
  }

  # initialize predictions vector and add the observed points
  zpred = vector(mode='numeric', length=prod(gdim))
  zpred[ pc[['s']] ] = zobs

  # omit any NA points from zobs
  zd = zobs[sobs]

  # copy eigenvector matrix
  edv = pc[['ed']][['vectors']]

  # pointwise variances are equal to the eigenvalues plus the nugget variance
  pwv = pars[['nug']] + pc[['ed']][['values']]

  # pad data vector with zeros (equivalent to subsetting cross-correlation matrices below)
  zd.mod = Matrix::Matrix(rep(0, nsg))

  # missing data case: not all points in subgrid are sampled
  if( nobs < nsg )
  {
    # ordinary case - reciprocal of eigenvalues produces the inverse
    zd.mod[sobs] = edv %*% ( (t(edv) %*% zd) / pwv )

    # multiply by cross covariance to get conditional mean of unsampled subgrid points
    zpred[ pc[['s']][-sobs] ] = pc[['vsc']] %*% zd.mod[sobs]


  } else {

    # kronecker product trick to get pointwise variances
    pwv = pars[['nug']] + kronecker(pc[['ed']][['x']][['values']], pc[['ed']][['y']][['values']])

    # kronecker product trick to get product with inverse covariance matrix
    edx = pc[['ed']][['x']][['vectors']]
    edy = pc[['ed']][['y']][['vectors']]
    zd.mod = pkern_kprod(edx, edy, pkern_kprod(edx, edy, zd, trans=T)/pwv, trans=F)
  }

  # compute predictions separately on the three subsets (third one is slowest by far)
  zpred[ pc[['oj']] ] = pkern_kprod(pc[['vs']][['x']], pc[['vso']][['y']], zd.mod, trans=T)
  zpred[ pc[['oi']] ] = pkern_kprod(pc[['vso']][['x']], pc[['vs']][['y']], zd.mod, trans=T)
  zpred[ pc[['o']] ] = pkern_kprod(pc[['vso']][['x']], pc[['vso']][['y']], zd.mod, trans=T)

  # finish
  return(zpred)
}





#
#' Simulate grid point values of a random field with separable covariance kernel
#'
#' Generates a simulation from the Gaussian random field with separable spatial
#' covariance kernel `pars` over the regular grid of size `dims`.
#'
#' Random vectors are generated by first drawing an iid normal vector of length
#' `n = prod(dims)` (using `base::stats::rnorm`), then taking its product with the covariance
#' matrix square root (via eigendecompositions). The result is returned (invisibly)
#' as a matrix.
#'
#' `pars` should be a list containing two kernel parameter definitions lists ("x" and
#' "y"), both compatible with `pkern_corr`, and, optionally, pointwise variance decomposed
#' as "v" and "nug" (the nugget variance). Defaults for "v" and "nug" are 1 and zero,
#' respectively.
#'
#' When `precompute==TRUE`, the function returns a list containing: "dims", the grid size;
#' "v", the list of ("x" and "y") component covariance matrices; and "ed" the list of their
#' eigendecompositions (as returned by `base::eigen`). This list can be passed
#' to subsequent `pkern_sim` calls in argument `precompute` to skip some computationally
#' expensive steps. In that case, note that argument `dims` is ignored (overwritten
#' by `precompute$dims`).
#'
#' Note also that the nugget effect is dealt with after precomputing - ie full covariance
#' matrix for the grid is equal to `kronecker(precompute$v$x, precompute$v$y)` plus the
#' nugget effect (added to diagonals).
#'
#' Note that eigendecomposition can produce negative eigenvalues (impossible in theory but
#' common in practice, due to numerical precision limits of the computer), a small nugget
#' effect is automatically added (with a warning) to make the component covariance matrices
#' numerically non-singular and allow simulation to proceed.
#'
#' @param pars a kernel parameters list ("k", "kp"), or list of two of them ("x", "y")
#' @param gdim integer vector c(ny, nx), the grid size
#' @param ds positive numeric or vector of two, the distance scaling factor(s)
#' @param precompute logical or list, for caching computationally expensive operations
#' @param makeplot logical, whether to plot the results
#'
#' @return a numeric matrix with dimensions `rev(dims)` or a list of precomputed objects
#' @export
#'
#' @examples
#' # basic usage with Gaussian kernel and default settings
#' pars = pkern_corr('gau')
#' sim = pkern_sim(pars)
#' utils::str(sim)
#'
#' # bigger example
#' gdim = c(ny=100, nx=50)
#' pkern_sim(pars, gdim=gdim)
#'
#' # store eigendecompositions to speed up repeated simulations
#' pre = pkern_sim(pars, gdim, precompute=TRUE, makeplot=FALSE)
#' utils::str(pre)
#' pkern_sim(pars, precompute=pre)
#'
#' # resolution can be adjusted via `ds`
#'  pkern_sim(pars, gdim, ds=c(1/4, 1/2))
#' # in this example nugget effect is added to the numerically singular matrix
pkern_sim = function(pars, gdim=c(25,25), ds=1, precompute=FALSE, makeplot=TRUE)
{
  # if only one component kernel is supplied, use it for both x and y
  if( all( c('k', 'kp') %in% names(pars) ) ) pars = list(y=pars, x=pars)

  # assign distance scaling with `pars` overwriting anything in argument `ds`
  if( 'ds' %in% names(pars) ) ds = pars[['ds']]
  if( length(ds) == 1 ) ds = rep(ds, 2)

  # set nugget effect to zero if it isn't supplied
  if( is.null( pars[['nug']] ) ) pars[['nug']] = 0

  # check if precomputed matrices and indices were supplied
  if( !is.logical(precompute) )
  {
    # check for invalid input
    if( !is.list(precompute) ) stop('"precompute" must be a list')

    # unpack x and y component correlation matrices and their eigendecompositions
    gdim = precompute[['gdim']]
    v = precompute[['v']]
    ed = precompute[['ed']]
    ds = precompute[['ds']]
    precompute = FALSE

  } else {

    # build x and y component correlation matrices, compute their eigendecompositions
    v = Map(\(p, d, s) pkern_corrmat(p, d, s), p=pars[c('y', 'x')], d=gdim, s=ds)
    ed = lapply(v, \(vmat) eigen(vmat, symmetric=TRUE))
  }

  # return the list of precomputed objects if requested
  if( precompute ) return( list(v=v, ed=ed, gdim=gdim, ds=ds) )

  # kronecker product trick to get pointwise variances
  pwv = pars[['nug']] + kronecker(ed[['x']][['values']], ed[['y']][['values']])

  # increase nugget to offset negative eigenvalues (numerical instability)
  if( any(pwv < 0) )
  {
    add.nug = -min(pwv)
    warning( paste('V is numerically singular. Nugget effect increased by', signif(add.nug, 3)) )
    pwv = pwv + add.nug
    pars[['nug']] = pars[['nug']] + add.nug
  }

  # generate independent standard normal data
  z = stats::rnorm( prod(gdim) )

  # equivalent to multiplying by full covariance matrix square root
  z.ortho = sqrt(pwv) * pkern_kprod(ed[['x']][['vectors']], ed[['y']][['vectors']], z, trans=TRUE)
  z.corr = pkern_kprod(ed[['x']][['vectors']], ed[['y']][['vectors']], z.ortho, trans=FALSE)

  # reshape as matrix and draw plot if requested
  z.corr.mat = matrix(z.corr, gdim[1])
  if( makeplot )
  {
    if( pars[['nug']] == 0 ) pars[['nug']] = NULL
    titles = pkern_toString(pars)
    main = bquote('simulated'~.(titles[['main']])~'random field'~.(titles[['nug']]))
    yx = Map(\(d, s) c(1, s*seq(d)[-d]), d=gdim, s=ds)
    pkern_plot(z.corr.mat, gdim, ds=ds, ppars=list(main=main, sub=titles[['kp']]))
  }

  # finish, returning data matrix
  return(invisible(z.corr.mat))
}


#
#' Reduce the resolution of a regular grid
#'
#' Returns a subgrid of `z` containing every `afact`th grid line. ie. only the rows
#' numbered `1`, `afact[1]`, `2*afact[1]`, .., and the columns numbered `1`, `afact[2]`,
#' `2*afact[2]`, ... are returned. This is mainly intended for reducing the size of
#' images before sending them to `graphics::image`.
#'
#' If `pars` is supplied, the function aggregates the data in `z` by convolving it
#' with the separable spatial kernel defined by `pars`. This can produce nicer looking
#' results than the default behaviour, and is also useful for revealing the locations of
#' sparse NAs in very large grids (as NAs are propegated in the convolution).
#'
#' @param z either a numeric matrix with dimensions `gdim`, or its column-vectorization
#' @param gdim integer vector c(ni, nj), the number of rows and columns in the grid
#' @param afact integer or vector of two integers, the aggregation factors along rows and columns
#' @param pars (optional) list of "y" and "x" kernel parameters or NA
#' @param discrete logical, indicating to round output to nearest integer
#'
#' @return numeric matrix of dimensions `floor(gdim/afact)`, the aggregated grid
#' @export
#'
#' @examples
#' # TODO
pkern_agg = function(z, gdim=NULL, afact=2, pars=NULL, discrete=FALSE)
{
  # duplicate aggregation factor and set grid size as needed
  if( length(afact) == 1 ) afact = rep(afact, 2)
  if( is.null(gdim) ) gdim = dim(z)
  if( is.null(gdim) ) stop('could not determine grid size. Check argument gdim')

  # set up dimensions of new grid
  ij = stats::setNames(mapply( \(d, b) seq(1, d, b), d=gdim, b=afact, SIMPLIFY=FALSE), c('i', 'j'))
  sgdim = sapply(ij, length)

  # set default kernel as needed
  if( is.null(pars) ) pars = 'gau'
  if( !anyNA(pars) )
  {
    # set default kernel parameters as needed
    yxnm = c('y', 'x')
    if( length(pars) == 1 | (!all(yxnm %in% names(pars))) ) pars = stats::setNames(rep(pars, 2), yxnm)
    pars = pkern_bds(pars)

    # adjust range parameters to assign (effectively) zero weight outside afact neighbourhood
    neighb.range = mapply( \(p, d) pkern_rho(1e-32, p, d), p=pars[yxnm], d=afact/2)
    pars[['y']][['kp']][1] = neighb.range[1]
    pars[['x']][['kp']][1] = neighb.range[2]

    # set up neighbourhood weights matrices
    Y = Matrix::Matrix( pkern_corrmat(pars[['y']], gdim[1], i=ij[['i']]) )
    X = Matrix::Matrix( pkern_corrmat(pars[['x']], gdim[2], j=ij[['j']]) )

    # compute normalization constants - set NAs wherever there is a division by zero
    cx = 1 / Matrix::colSums(X)
    cy = 1 / Matrix::rowSums(Y)
    cx[is.infinite(cx)] = NA
    cy[is.infinite(cy)] = NA

    # normalize convolution matrices - note that zero weight -> NA -> NA row
    Y = Matrix::Diagonal(sgdim[1], cy) %*% Y
    X = X %*% Matrix::Diagonal(sgdim[2], cx)

    # convolve with data
    zout = as.vector( Y %*% matrix(z, gdim) %*% X )

    # for discrete (integer-valued) data we want integer output
    if(discrete) zout = round(zout)

  } else {

    # when pars is NA we just do a regular sampling of the grid, a la `raster::plot`
    zout = z[pkern_idx_sg(gdim, ij)]
  }

  return( list(z=zout, ij=ij) )
}
