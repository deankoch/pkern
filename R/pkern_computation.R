#
# pkern_raster.R
# Dean Koch, Sep 2021
# Functions for computationally intensive stuff
#


# efficient quadratic form computer: `t(x) %*% M %*% y`
pkern_qf = function(M, x, y)
{
  # Efficient computation of the quadratic form `t(X) %*% M %*% Y`
  # where M is a numeric matrix with conforming dimensions
  #
  # based on the trick used in emulator::quad.form with Rfast::Crossprod
  # for efficiency with large matrices, where I see around a 2-4X speed speedup

  # coerce everything to matrices with appropriate dimensions
  M = as.matrix(M)
  ninner = dim(M)[1]
  nouter = dim(M)[2]
  if( is.vector(x) ) x = matrix(x, ninner)
  if( is.vector(y) ) y = matrix(y, nouter)

  # verify that dimensions conform (Rfast::Crossprod will crash R otherwise!)
  conforms.inner = ninner == dim(x)[1]
  conforms.right = nouter == dim(y)[1]
  if( !conforms.inner | !conforms.right ) stop('dimensions do not conform')

  # run computation
  return( as.vector(Rfast::Crossprod(Rfast::Crossprod(M, x), y)) )
}


# left-multiplication with a Kronecker product: evaluates (X (x) Y) * vector
pkern_kprod = function(X, Y, z, trans=FALSE)
{
  # handle transposed mode
  if(trans) return( pkern_qf(matrix(z, nrow(Y), nrow(X)), Y, X) )
  return( pkern_qf(matrix(z, ncol(Y), ncol(X)), Rfast::transpose(Y), Rfast::transpose(X)) )
}


# uses Rccp libraries in Rfast package and ...

# conditional mean of points on a grid given a sample taken over a subgrid
pkern_cmean = function(zobs, dims, pars, gxy, precompute=FALSE)
{
  # Fast conditional mean computer for problems where sampled locations form a subgrid
  #
  # Uses method of Gilboa et. al, 2015 where eigendecomposition of the covariance
  # matrix of sampled locations is used with a nugget effect
  #
  # ARGUMENTS:
  #
  # dims: vector of two integers (nx, ny), the dimensions of the full grid
  # pars: list containing kernel parameters "x" and "y", nugget size "nug", and resolution "ds"
  # zobs: numeric vector of observed data, possibly with NAs
  # gxy: list of two integer vectors, the x and y grid lines of the subgrid
  # precompute: either logical (default FALSE), or a list of precomputed objects (see details)
  #
  # DETAILS:
  #
  # zobs should be supplied in the usual column-vectorized order.
  #
  # The length of zobs must equal the product of the lengths of the vectors in gxy
  # (ie the number of grid line intersections), and the grid line numbers in gxy
  # should be ordered subsets of the sequences seq(dims[1]) and seq(dims[2])
  #
  # the function returns a numeric vector of length prod(dims), containing the
  # observed data from zobs and the mean values at all other grid points, conditional
  # on zobs (in column-vectorized order).
  #

  # the number of points in the sampled subgrid
  nsg = length(zobs)

  # check if precomputed matrices and indices were supplied
  if( !is.logical(precompute) )
  {
    # check for invalid input
    if( !is.list(precompute) ) stop('"precompute" must be a list')

    # unpack correlation matrices
    vx = precompute[['vx']]
    vy = precompute[['vy']]
    vx.cross = precompute[['vx.cross']]
    vy.cross = precompute[['vy.cross']]
    vs.cross = precompute[['vs.cross']]

    # eigendecomposition(s)
    ed = precompute[['ed']]

    # unpack indexing vectors
    idx.on = precompute[['idx.on']]
    idx.off = precompute[['idx.off']]
    idx.xoff = precompute[['idx.xoff']]
    idx.yoff = precompute[['idx.yoff']]
    idx.sg.obs = precompute[['idx.sg.obs']]

    # don't return the precomputed data in this case
    precompute = FALSE

    # TODO: verify that zobs agrees with the (possibly precomputed) index of NAs

  } else {

    # computationally expensive steps which are independent of `zobs`:

    ## compute important indices in vectorized output grid

    # sort x and y grid line indices for observed data and find their complements
    gx = sort( gxy[[1]] )
    gy = sort( gxy[[2]] )
    gxc = seq( dims[1] )[-gx]
    gyc = seq( dims[2] )[-gy]
    nx = length(gx)
    ny = length(gy)

    # index of subgrid points in full grid, and all points on one of its grid lines
    idx.on = pkern_idx_sg(dims, i=gy, j=gx)
    idx.on.x = pkern_idx_sg(dims, i=NULL, j=gx)
    idx.on.y = pkern_idx_sg(dims, i=gy, j=NULL)

    # partition non-subgrid points into three sets depending on overlap with grid lines
    idx.off = seq( prod(dims) )[ -c(idx.on, idx.on.x, idx.on.y) ]
    idx.xoff = idx.on.x[ !( idx.on.x %in% idx.on ) ]
    idx.yoff = idx.on.y[ !( idx.on.y %in% idx.on ) ]

    ## compute correlation matrices

    # find relative increase in resolution going from subgrid to full grid
    dsx = diff( gx[1:2] )
    dsy = diff( gy[1:2] )

    # regularize range parameters for computations over integer lattice
    pars[['x']][['kp']][1] = pars[['x']][['kp']][1] / pars[['ds']][1]
    pars[['y']][['kp']][1] = pars[['y']][['kp']][1] / pars[['ds']][2]

    # build component marginal correlation matrices for the subgrid
    vx = pkern_corrmat(pars[['x']], nx)
    vy = pkern_corrmat(pars[['y']], ny)

    # build component cross-correlation matrices for points on vs off subgrid lines
    vx.cross = pkern_corrmat(pars[['x']], dims[1], ds=1/dsx, j=gxc, i=gx)
    vy.cross = pkern_corrmat(pars[['y']], dims[2], ds=1/dsy, j=gyc, i=gy)

    ## compute eigendecomposition(s)

    # slower method when not all of the subgrid is sampled
    if( anyNA(zobs) )
    {
      # indices of sampled and unsampled within the subgrid
      idx.sg.obs = which( !is.na(zobs) )
      idx.sg.unobs = seq(nsg)[-idx.sg.obs]

      # missing data case requires first evaluating the variance kronecker product
      v = kronecker(vx, vy)
      ed = eigen(v[idx.sg.obs, idx.sg.obs], symmetric=TRUE)

      # we also need the (within subgrid) cross-covariance between sampled and unsampled points
      vs.cross = v[idx.sg.unobs, idx.sg.obs]

    } else {

      # with no missing data we can do eigendecompositions on component matrices
      idx.sg.obs = seq(nsg)
      ed = list(x=eigen(vx, symmetric=TRUE), y=eigen(vy, symmetric=TRUE))
      vs.cross = NULL
    }
  }

  # return the list of precomputed objects if requested
  if( precompute ) return(list(vx=vx,
                               vy=vy,
                               ed=ed,
                               vx.cross=vx.cross,
                               vy.cross=vy.cross,
                               vs.cross=vs.cross,
                               idx.on=idx.on,
                               idx.off=idx.off,
                               idx.xoff=idx.xoff,
                               idx.yoff=idx.yoff,
                               idx.sg.obs=idx.sg.obs))

  # initialize predictions vector and add the observed points
  zout = rep(as.numeric(NA), prod(dims))
  zout[idx.on] = zobs

  # missing data case
  if( length(idx.sg.obs) < nsg )
  {
    # omit the NA points from zobs
    zobs.crop = zobs[idx.sg.obs]

    # pointwise variances are equal to the eigenvalues plus the nugget variance
    pwv = pars[['nug']] + ed[['values']]

    # pad output vector with zeros (equivalent to subsetting cross-correlation matrices below)
    zobs.indep = rep(0, nsg)

    # transform `zobs` by product with inverse covariance to get mutual independence
    zobs.indep[idx.sg.obs] = ed[['vectors']] %*% ( ( t(ed[['vectors']]) %*% zobs.crop ) / pwv )

    # compute and assign conditional mean of unsampled subgrid points
    zout[ idx.on[-idx.sg.obs] ] = vs.cross %*% zobs.indep[idx.sg.obs]

  } else {

    # kronecker product trick to get pointwise variances
    pwv = pars[['nug']] + kronecker(ed[['x']][['values']], ed[['y']][['values']])

    # kronecker product trick to get product with covariance inverse
    zobs.ortho = pkern_kprod(ed[['x']][['vectors']], ed[['y']][['vectors']], zobs, trans=TRUE) / pwv
    zobs.indep = pkern_kprod(ed[['x']][['vectors']], ed[['y']][['vectors']], zobs.ortho, trans=FALSE)
  }

  # compute predictions separately on the three subsets
  zout[idx.xoff] = pkern_kprod(vx, vy.cross, zobs.indep, trans=TRUE)
  zout[idx.yoff] = pkern_kprod(vx.cross, vy, zobs.indep, trans=TRUE)
  zout[idx.off] = pkern_kprod(vx.cross, vy.cross, zobs.indep, trans=TRUE)

  # finish
  return(zout)
}


#
#' Simulate grid point values of a random field with separable covariance kernel
#'
#' Generates a simulation from the Gaussian random field with separable spatial
#' covariance kernel `pars` over the regular grid of size `dims`.
#'
#' Random vectors are generated by first drawing an iid normal vector of length
#' `n = prod(dims)` (using `base::rnorm`), then taking its product with the covariance
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
#' @param dims integer vector c(nx, ny), the grid size
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
#' str(sim)
#' pkern_sim(pars, dims=c(1e2, 2e2))
#'
#' # store eigendecompositions to speed up repeated simulations
#' pre = pkern_sim(pars, dims=c(5e2, 5e2), precompute=TRUE, makeplot=FALSE)
#' str(pre)
#' pkern_sim(pars, precompute=pre)
#'
#' # example of numerically singular matrices triggering nugget effect
#' pars = modifyList(pars, list(kp=1e2))
#' pkern_sim(pars)
pkern_sim = function(pars, dims=c(25,25), ds=1, precompute=FALSE, makeplot=TRUE)
{
  # if only one component kernel is supplied, use it for both x and y
  if( all( c('k', 'kp') %in% names(pars) ) ) pars = list(x=pars, y=pars)

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
    v = precompute[['v']]
    ed = precompute[['ed']]
    dims = precompute[['dims']]
    ds = precompute[['ds']]
    precompute = FALSE

  } else {

    # build x and y component correlation matrices, compute their eigendecompositions
    v = mapply(\(p, d, s) pkern_corrmat(p, d, s), p=pars[c('x', 'y')], d=dims, s=ds, SIMPLIFY=FALSE)
    ed = lapply(v, \(vmat) eigen(vmat, symmetric=TRUE))
  }

  # return the list of precomputed objects if requested
  if( precompute ) return( list(v=v, ed=ed, dims=dims, ds=ds) )

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
  z = rnorm( prod(dims) )

  # equivalent to multiplying by full covariance matrix square root
  z.ortho = sqrt(pwv) * pkern_kprod(ed[['x']][['vectors']], ed[['y']][['vectors']], z, trans=TRUE)
  z.corr = pkern_kprod(ed[['x']][['vectors']], ed[['y']][['vectors']], z.ortho, trans=FALSE)

  # reshape as matrix and draw plot if requested
  z.corr.mat = matrix(z.corr, dims[2])
  if( makeplot )
  {
    if( pars[['nug']] == 0 ) pars[['nug']] = NULL
    titles = pkern_toString(pars)
    main = bquote('simulated'~.(titles[['main']])~'random field'~.(titles[['nug']]))
    xy = mapply( \(d, s) c(1, s*seq(d)[-d]), d=dims, s=ds, SIMPLIFY=FALSE)
    pkern_plot(z.corr.mat, xy=xy, ppars=list(main=main, sub=titles[['kp']]))
  }

  # finish, returning data matrix
  return(invisible(z.corr.mat))
}


#
#' Reduce the resolution of a regular grid by separable convolution
#'
#' Aggregates the data in `z` by convolving it with the separable spatial kernel
#' defined by `pars`, producing a smaller output matrix containing every `afact`th
#' grid line. ie. only the x grid lines numbered 1, afact[1], 2*afact[1], ... etc.
#'
#' This is a heuristic algorithm intended for reducing the size of images before
#' sending them to `graphics::image`.
#'
#' @param z either a numeric matrix with dimensions `dims`, or its column-vectorization
#' @param dims integer vector c(nx, ny), the grid size
#' @param afact integer or vector of two integers, the aggregation factors
#' @param pars (optional) list of "x" and "y" kernel parameters
#' @param tol numeric between 0 and 1, tolerance for NAs
#' @param discrete logical, indicating to round output to nearest integer
#'
#' @return numeric matrix of dimensions `floor(dims/afact)`, the aggregated grid
#' @export
#'
#' @examples
pkern_agg = function(z, dims=NULL, afact=2, pars=NULL, tol=0.05, discrete=FALSE)
{
  # duplicate aggregation factor and set grid size as needed
  if( length(afact) == 1 ) afact = rep(afact, 2)
  if( is.null(dims) ) dims = dim(z)
  if( is.null(dims) ) stop('could not determine grid size. Check argument dims')

  # set default kernel parameters as needed
  xynm = c('x', 'y')
  if( is.null(pars) ) pars = 'gau'
  if( length(pars) == 1 | ( !all( xynm %in% names(pars) ) ) ) pars = setNames(rep(pars, 2), xynm)
  pars = pkern_bds(pars)

  # adjust range parameters to assign (effectively) zero weight outside afact neighbourhood
  pars[['x']][['kp']][1] = pkern_rho(1e-32, pars[['x']], d=afact[1]/2)
  pars[['y']][['kp']][1] = pkern_rho(1e-32, pars[['y']], d=afact[2]/2)

  # set up dimensions of new grid
  ij = setNames(mapply( \(d, b) seq(1, d, b), d=rev(dims), b=afact, SIMPLIFY=FALSE), c('i', 'j'))
  dims.sg = rev(sapply(ij, length))

  # set up neighbourhood weights matrices
  X = Matrix::Matrix( pkern_corrmat(pars[['x']], dims[1], j=ij[['j']]) )
  Y = Matrix::Matrix( pkern_corrmat(pars[['y']], dims[2], i=ij[['i']]) )

  # compute normalization constants - set NAs wherever there is a division by zero
  cx = 1 / Matrix::colSums(X)
  cy = 1 / Matrix::rowSums(Y)
  cx[is.infinite(cx)] = NA
  cy[is.infinite(cy)] = NA

  # normalize convolution matrices - note that zero weight -> NA -> NA row
  X = X %*% Matrix::Diagonal(dims.sg[1], cx)
  Y = Matrix::Diagonal(dims.sg[2], cy) %*% Y

  # # convolve with sparse indicator matrix of NAs
  # miss = Matrix::Matrix(0L, dims[2], dims[1])
  # miss[is.na(z)] = 1L
  # missout = Y %*% miss %*% X

  # convolve with data
  zout = as.vector( Y %*% matrix(z, dims[2]) %*% X )

  # for discrete (integer-valued) data we want integer output
  if(discrete) zout = round(zout)

  # assign NA to any point where over 5% of its weight comes from NA source points
  #zout[Matrix::which(missout > tol)] = NA
  return( list(z=zout, ij=ij) )
}
