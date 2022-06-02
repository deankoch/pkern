# pkern_var.R
# Dean Koch, 2022
# Functions for matrix algebra with correlation and covariance matrices

#' Stationary 1D correlation kernels
#'
#' Computes stationary correlation function values for the n (nonegative) 1-dimensional
#' distances in `d`. Parameter list entry `pars$kp` supplies the kernel parameter(s).
#'
#' `pars$k` must be one of the following kernel names:
#'
#' "exp": exponential (special case of "gex" with shape "p"=1)
#' "gau": gaussian/stable (special case of "gex" with shape "p"=2)
#' "sph": spherical (AKA stable/Gaussian for p=2)
#'
#' "gex": gamma-exponential (with shape "p")
#' "mat": Whittle-Matern (Handcock and Wallis parameterization, with shape "kap")
#'
#' where the first three have 1 single range parameter, and the last two have both a
#' range and shape parameter.
#'
#' For the 1-parameter kernels, `pars$kp` is the range parameter value ("rho" ); For the
#' 2-parameter kernels, `pars$kp` is a vector whose first element is "rho", and second
#' element is the shape parameter ("p" or "kap"). Note that names in `pars$kp` are ignored
#' and only the order matters - the range parameter always comes first.
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
#' pars = pkern_pars(10)[['x']]
#' pkern_corr(pars, d=1:10)
#'
pkern_corr = function(pars, d=NA)
{
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


#' Construct 1D stationary correlation matrices for regularly spaced data
#'
#' An effient implementation that uses symmetry and Toeplitz structure arising
#' from assumption of stationarity of the random field and regularity of the grid.
#'
#' `gres` scales the distance between adjacent points
#'
#' @param pars list of kernel parameters "k" and "kp" (see `pkern_corr`)
#' @param n positive integer, the number of points on the 1D line
#' @param gres positive numeric, the distance between adjacent grid lines
#' @param i vector, a subset of `seq(n)` indicating rows to return
#' @param j vector, a subset of `seq(n)` indicating columns to return
#'
#' @return the n x n correlation matrix, or its subset as specified in `i`, `j`
#' @export
#'
#' @examples
#' pars_x = pkern_pars(10, 'gxp')[['x']]
#' pkern_corr_mat(pars_x, n=10)
#' pkern_corr_mat(pars_x, n=10, i=2:4, j=2:4)
#' pkern_corr_mat(pars_x, n=3)
#' pkern_corr_mat(pars_x, n=3, gres=2)
pkern_corr_mat = function(pars, n, gres=1, i=seq(n), j=seq(n))
{
  # compute the set of distances over which we need to evaluate kernel
  du = gres * ( seq(n) - 1 )

  # compute kernel values for these distances
  dcorr = pkern_corr(pars, du)

  # build large vector to shift through in building the Toeplitz output matrix
  bigvec = c(dcorr[n:2], dcorr)

  # build and return the matrix
  return( sapply(j, function(x) bigvec[ (n-x) + i ]) )
}


#' Compute covariance matrix for a (sub)sample of gridded data
#'
#' Computes the covariance matrix `V` for observed grid data `g_obs` and covariance model
#' `pars` (default `method=='none'`) or the quadratic form `t(X) %*% V_inv %*% X` where
#' `X` is a covariates matrix and `V_inv` is the inverse of `V`.
#'
#' When `method=='eigen'` the function instead returns the eigen-decomposition of this
#' matrix; and when `method=='chol'` it returns the lower triangular Cholesky factor.
#'
#' NAs in the data vector `g_obs$gval` indicate to return the marginal covariance, ie
#' the sub-matrix of V where the rows and columns corresponding to the NA grid points
#' are omitted. If `g_obs$gval` is a matrix, it is assumed to be an array of data vectors
#' (in the columns), all having NA structure identical to the first column.
#'
#' The data vector(s) can be absent from `g_obs` (ie `is.null(g_obs$gval) == TRUE`),
#' in which case the function computes the full component correlation matrices (or their
#' factorization), for dimensions x and y, and returns them in a list (except when `X`
#' is supplied, since the output matrix is not generally separable in this case).
#'
#' When `scaled=TRUE` and the data vector is supplied in `g_obs`, the function returns
#' `V/pars$psill` (or the appropriate sub-matrix, or factorization). If `X` is also supplied,
#' then the product `t(X) %*% V_inv %*% X` is scaled by this factor (but `V_inv` is computed
#' without the scaling).
#'
#' @param g_obs list of form returned by `pkern_grid` (with entries 'gdim', 'gres', 'gval')
#' @param pars list of form returned by `pkern_pars` (with entries 'y', 'x', 'eps', 'psill')
#' @param scaled logical, if `TRUE` sets `pars$psill = 1` and `eps = pars$eps / pars$psill`
#' @param method character, the factorization to return, one of 'none', 'chol', 'eigen'
#' @param X numeric matrix, the `X` in `t(X) %*% V %*% X` (default is identity)
#'
#' @return either matrix `V`, or `t(X) %*% V_inv %*% X`, or a factorization ('chol' or 'eigen')
#' @export
#'
#' @examples
#' # define example grid with NAs and example predictors matrix
#' gdim = c(12, 13)
#' n = prod(gdim)
#' n_obs = floor(n/3)
#' idx_obs = sort(sample.int(n, n_obs))
#' g_obs = pkern_grid(gdim)
#' g_obs$gval[idx_obs] = rnorm(n_obs)
#' nX = 3
#' X_all = rnorm(nX * n) |> matrix(ncol=nX)
#'
#' # example kernel
#' psill = 0.3
#' pars = pkern_pars(g_obs) |> modifyList(list(psill=psill))
#'
#' # plot the full covariance matrix, its cholesky factor and eigen-decomposition
#' V = pkern_var(g_obs, pars)
#' V_chol = pkern_var(g_obs, pars, method='chol')
#' V_eigen = pkern_var(g_obs, pars, method='eigen')
#' pkern_plot(V)
#' pkern_plot(V_chol)
#' pkern_plot(V_eigen$vectors)
#'
#' # with no NAs (or no data at all) the function returns the correlation matrix components
#' g_nodata = modifyList(g_obs, list(gval=NULL))
#' str(pkern_var(g_nodata, pars))
#' str(pkern_var(g_nodata, pars, method='chol'))
#' str(pkern_var(g_nodata, pars, method='eigen'))
#'
#' # get the full covariance matrix with sep=FALSE...
#' V2 = pkern_var(g_nodata, pars, sep=FALSE)[idx_obs, idx_obs]
#' max(abs( V - V2 ))
#'
#' # ... or compute it yourself from the components
#' corr_components = pkern_var(g_nodata, pars)
#' corr_mat = kronecker(corr_components[['x']], corr_components[['y']])
#' nugget_effect = diag(pars$eps, n)
#' V3 = pars$psill * corr_mat + nugget_effect
#' max(abs(V - V3[idx_obs, idx_obs]))
#'
#' # test quadratic form with X
#' X = cbind(1, X_all)[idx_obs, ]
#' cprod = crossprod(X, chol2inv(chol(V))) %*% X
#' abs(max(pkern_var(g_obs, pars, X=X) - cprod ))
#'
#' # test products with inverse of quadratic form with X
#' z = rnorm(nX + 1)
#' cprod_inv = chol2inv(chol(cprod))
#' pars0 = pars |> modifyList(list(psill=1, eps=0))
#' cprod_inv_chol = pkern_var(g_obs, pars, X=X, method='chol')
#' pkern_var_mult(z, pars0, fac=cprod_inv_chol) - (cprod_inv %*% z)
#'
#' # `scaled` indicates to divide matrix by psill
#' print( pars[['eps']]/pars[['psill']] )
#' diag(pkern_var(g_obs, pars, scaled=TRUE)) # diagonal elements equal to 1 + eps/psill
#' ( pkern_var(g_obs, pars) - psill * pkern_var(g_obs, pars, scaled=TRUE) ) |> abs() |> max()
#' ( pkern_var(g_obs, pars, X=X, scaled=TRUE) - ( cprod/psill ) ) |> abs() |> max()
#'
#' # in cholesky factor this produces a scaling by square root of psill
#' max(abs( V_chol - sqrt(psill) * pkern_var(g_obs, pars, method='chol', scaled=TRUE) ))
#' # and in the eigendecomposition, a scaling of the eigenvalues
#' vals_scaled = pkern_var(g_obs, pars, method='eigen', scaled=TRUE)$values
#' max(abs( pkern_var(g_obs, pars, method='eigen')$values - psill*vals_scaled ))
#'
pkern_var = function(g_obs, pars=NULL, scaled=FALSE, method='none', X=NULL, fac=NULL, sep=TRUE)
{
  # default Gaussian kernel
  if(is.null(pars)) pars = pkern_pars(g_obs, 'gau')

  # check for invalid method
  nm_method = c('none', 'eigen', 'chol')
  if( is.null(method) ) method = 'chol'
  msg_method = paste('method must be one of: NA,', paste(nm_method, collapse=', '))
  if( !(method %in% nm_method) ) stop(msg_method)

  # identify NA grid-points and handle separable case (no missing data, or all missing)
  n = prod(g_obs[['gdim']])
  if(is.matrix(g_obs[['gval']])) g_obs[['gval']] = as.vector(g_obs[['gval']][,1])
  is_obs = !is.na(g_obs[['gval']])
  if( !any(is_obs) | all(is_obs) ) g_obs[['gval']] = NULL
  if( is.null(g_obs[['gval']]) ) is_obs = rep(TRUE, n)

  # predictor matrix case
  if( !is.null(X) )
  {
    print('DEBUGGING: X predictor case')
    # # separable case (no missing data)
    # if( all(is_obs) ) g_obs[['gval']] = NULL

    # call without X to get eigen-decomposition of variance
    X_method = ifelse(method=='none', 'eigen', method)
    if( is.null(fac) ) fac = pkern_var(g_obs, pars, scaled=TRUE, method=X_method)

    # check for invalid input
    msg_class = 'mu must be a matrix predictor columns'
    msg_mismatch = 'nrow(X) must equal the number of non-NA points in g_obs$gval'
    if( !is.matrix(X) ) stop(msg_class)
    if( nrow(X) != sum(is_obs) ) stop(msg_mismatch)

    # quadratic form of whitened data matrix (not the inverse)
    p_scale = ifelse(scaled, pars[['psill']], 1)
    X_quad = pkern_var_mult(X, pars, fac=fac, method=X_method, quad=TRUE) / p_scale

    # return requested decomposition
    if(method == 'none') return(X_quad)
    if(method == 'chol' ) return(t(chol(X_quad)))
    if(method == 'eigen' ) return(eigen(X_quad, symmetric=TRUE))
  }

  # unpack grid config and covariance parameters
  n_obs = sum(is_obs)
  gres = g_obs[['gres']]
  gdim = g_obs[['gdim']]
  eps = ifelse(scaled, pars[['eps']]/pars[['psill']], pars[['eps']])
  psill = ifelse(scaled, 1, pars[['psill']])

  # case of no data vector in arguments
  if( n == n_obs )
  {
    # return the full component correlation matrices (or their factorizations)
    cy = pkern_corr_mat(pars[['y']], gdim[['y']], gres[['y']])
    cx = pkern_corr_mat(pars[['x']], gdim[['x']], gres[['x']])
    if(method == 'chol') return( list( y=t(chol(cy)), x=t(chol(cx)) ) )
    if(method == 'eigen') return( list( y=eigen(cy), x=eigen(cx) ) )
    if(method == 'none')
    {
      if(sep) return( list(y=cy, x=cx) )
      return( eps * diag(1, n_obs) + psill * kronecker(cx, cy) )
    }
  }

  # build mapping from points to rows of component matrices
  yx_idx = which(is_obs) |> pkern_vec2mat(gdim['y'], out='list')

  # draw a selection of rows/columns from component correlation matrices
  cy_obs = pkern_corr_mat(pars[['y']], gdim[['y']], gres[['y']], i=yx_idx[['i']], j=yx_idx[['i']])
  cx_obs = pkern_corr_mat(pars[['x']], gdim[['x']], gres[['x']], i=yx_idx[['j']], j=yx_idx[['j']])

  # Hadamard product produces the correlation matrix for observed grid points
  if(method == 'none') return( ( psill * cy_obs * cx_obs ) + diag( rep(eps, n_obs) ) )

  # compute Cholesky factor (lower triangular) of correlation matrix
  if(method == 'chol') return( t( chol( psill * cy_obs * cx_obs + eps * diag(rep(1, n_obs)) ) ) )

  # method='eigen': compute eigen-decomposition of correlation matrix
  eigen_result = eigen(cy_obs*cx_obs, symmetric=TRUE)
  eigen_result[['values']] = ( psill * eigen_result[['values']] ) + eps
  return(eigen_result)
}

#' Multiply a vector by the inverse covariance matrix
#'
#' The function computes `W %*% z` where W=V^p, V is the covariance matrix for data `z`, and,
#' by default, `p`  is -1 (matrix inverse). In 'eigen' mode, the argument `p` can be set any
#' integer or fractional power. In 'chol' mode, argument `p` is ignored and W is the inverse
#' of V.
#'
#' Alternatively, `out='quad'` computes the quadratic form `t(z) %*% W %*% z`.
#'
#' `method` specifies the covariance matrix factorization to use for computations, either
#' 'chol' (Cholesky factorization), which is fastest, or 'eigen' (eigen-decomposition),
#' which supports matrix powers other than `p=-1`.
#'
#' Factorization is the slow part of the computation. It can be pre-computed
#' using `pkern_var(..., scaled=TRUE)` and passed to `pkern_var_mult` in argument `fac`.
#' This is the factorization of the covariance matrix after scaling by the partial sill
#' (see `?pkern_var`); it must either be the lower Cholesky factor (the transposed output
#' of chol), or a list of eigen-vectors and eigen-values (the output of eigen).
#'
#' Note that when `fac` is supplied, all entries in `pars` are ignored except for `psill`,
#' as they are baked into the eigen-decomposition already. `g_obs` can in this case be a
#' numeric vector, the vector of observed data (with NAs omitted).
#'
#' @param g_obs list of form returned by `pkern_grid` or numeric vector or matrix of non-NA data
#' @param pars list of form returned by `pkern_pars` (with entries 'y', 'x', 'eps', 'psill')
#' @param fac factorization of scaled covariance matrix of z (V divided by psill)
#' @param quad logical, if TRUE the function returns the quadratic form `t(z) %*% V_inv %*% z`
#' @param p numeric, the matrix power of V^p to multiply (ignored when `method=='chol'`)
#'
#' @return numeric matrix
#' @export
#'
#' @examples
#' # relative error comparing output x to reference y
#' rel_err = \(x, y) ifelse(y == 0, 0, abs( (x - y) / y ) )
#'
#' # define example grid and data
#' gdim = c(10, 15)
#' n = prod(gdim)
#' z_all = rnorm(n)
#' g_obs = modifyList(pkern_grid(gdim), list(gval = z_all))
#'
#' # define covariance parameters
#' pars = pkern_pars(g_obs, 'gau') |> modifyList(list(psill=2, eps=0.5))
#'
#' # COMPLETE CASE
#'
#' V = pkern_var(g_obs, pars, method='none', sep=FALSE)
#' V_inv = chol2inv(chol(V))
#' out_reference = V_inv %*% z_all
#' out_reference_quad = t(z_all) %*% out_reference
#' pkern_var_mult(g_obs, pars) |> rel_err(out_reference) |> max()
#' pkern_var_mult(g_obs, pars, quad=TRUE) |> rel_err(out_reference_quad)
#'
#' # pre-computed factorization on separable components of correlation matrix
#' fac_corr = pkern_var(modifyList(g_obs, list(gval=NULL)), pars, method='eigen')
#' pkern_var_mult(g_obs, pars, fac=fac_corr) |> rel_err(out_reference) |> max()
#' pkern_var_mult(g_obs, pars, fac=fac_corr, quad=TRUE) |> rel_err(out_reference_quad)
#'
#' # matrix powers
#' out_reference = V %*% z_all
#' pkern_var_mult(g_obs, pars, method='eigen', p=1) |> rel_err(out_reference) |> max()
#' pkern_var_mult(g_obs, pars, method='eigen', p=1, quad=TRUE) |> rel_err(t(z_all) %*% out_reference)
#'
#' # INCOMPLETE CASE
#'
#' n_sample = floor(n/10)
#' idx_sampled = sample.int(n, n_sample) |> sort()
#' z_obs = rep(NA, n)
#' z_obs[idx_sampled] = z_all[idx_sampled]
#' g_obs = modifyList(g_obs, list(gval = z_obs))
#' V = pkern_var(g_obs, pars, method='none')
#' pkern_plot(V)
#'
#' # correctness check
#' z = matrix(z_obs[idx_sampled], ncol=1)
#' V_inv = chol2inv(chol(V))
#' out_reference = (V_inv %*% z)
#' out_reference_quad = t(z) %*% out_reference
#'
#' # check error for two output types by Cholesky method
#' pkern_var_mult(g_obs, pars) |> rel_err(out_reference) |> max()
#' pkern_var_mult(g_obs, pars, quad=TRUE) |> rel_err(out_reference_quad)
#'
#' # check eigen-decomposition method
#' pkern_var_mult(g_obs, pars, method='eigen') |> rel_err(out_reference) |> max()
#' pkern_var_mult(g_obs, pars, quad=TRUE, method='eigen') |> rel_err(out_reference_quad)
#'
#' # supply data as a vector instead of list by pre-computing factorization
#' fac_chol = pkern_var(g_obs, pars, scaled=TRUE, method='chol')
#' fac_eigen = pkern_var(g_obs, pars, scaled=TRUE, method='eigen')
#' pkern_var_mult(z, pars, fac=fac_chol) |> rel_err(out_reference) |> max()
#' pkern_var_mult(g_obs, pars, fac=fac_eigen) |> rel_err(out_reference) |> max()
#' pkern_var_mult(z, pars, fac=fac_chol, quad=TRUE) |> rel_err(out_reference_quad)
#' pkern_var_mult(g_obs, pars, fac=fac_eigen, quad=TRUE) |> rel_err(out_reference_quad)
#'
#' # matrix powers in eigen mode
#' out_reference = V %*% z
#' pkern_var_mult(g_obs, pars, method='eigen', p=1) |> rel_err(out_reference) |> max()
#' pkern_var_mult(g_obs, pars, method='eigen', p=1, quad=TRUE) |> rel_err(t(z) %*% out_reference)
#' pkern_var_mult(g_obs, pars, method='eigen', p=2) |> rel_err(V %*% out_reference) |> max()
#'
#' # multiply g_obs twice by a square root of V
#' g_obs_sqrt = g_obs
#' g_obs_sqrt$gval[!is.na(g_obs$gval)] = pkern_var_mult(g_obs, pars, method='eigen', p=1/2)
#' pkern_var_mult(g_obs_sqrt, pars, method='eigen', p=1/2) |> rel_err(out_reference) |> max()
#'
pkern_var_mult = function(g_obs, pars, method=NULL, fac=NULL, quad=FALSE, p=-1)
{
  # detect method from class of factorization fac (matrix=chol, list=eigen)
  if(is.null(method)) method = c('chol', 'eigen')[ 1 + as.integer( is.list(fac) ) ]

  # unpack list g_obs as needed
  if( is.list(g_obs) )
  {
    # when there are missing data, omit from copy z
    is_obs = !is.na(g_obs[['gval']])
    z = matrix(g_obs[['gval']][is_obs], ncol=1)

    # complete and empty cases trigger separability option below
    is_all_obs = all(is_obs) | !any(is_obs)

  } else {

    # coerce g_obs to matrix
    #z = g_obs
    #if( !is.matrix(g_obs) ) z = matrix(g_obs, 1)
    z = matrix(g_obs, ncol=ifelse(is.vector(g_obs), 1, ncol(g_obs)))

    # check if the supplied factorization is a Kronecker product
    if( is.null(fac) ) stop('factorization fac must be supplied if g_obs is not a list')
    is_all_obs = is.list(fac) & all(c('y', 'x') %in% names(fac))
  }

  # separability property and eigen-decomposition used in no-missing case
  if( is_all_obs )
  {
    # eigen method is forced in this case
    if( is.null(fac) ) fac = pkern_var(g_obs[c('gres', 'gdim')], pars=pars, method='eigen')
    if( !all(c('y', 'x') %in% names(fac) ) ) stop('components "x" and "y" not found in fac')
    ny = length(fac[['y']][['values']])
    nx = length(fac[['x']][['values']])

    # eigenvalues of full correlation and covariance matrices
    ev_corr = kronecker(fac[['x']][['values']], fac[['y']][['values']])
    ev_p = ( pars[['eps']] + pars[['psill']] * as.numeric(ev_corr) )^p
    if( any( is.infinite(ev_p) ) ) stop('ill-conditioned covariance matrix (0 eigenvalue)')
    if( any( is.nan(ev_p) ) ) stop('ill-conditioned covariance matrix (eigenvalue < 0)')

    # left multiply data vector by square root of inverse correlation matrix
    g_evec = z |>
      apply(2, \(x) crossprod(fac[['y']][['vectors']], matrix(x, ny)), simplify=FALSE ) |>
      lapply(\(x) tcrossprod(x, t(fac[['x']][['vectors']]) ) )

    # quadratic form is the scaled vector multiplied by its transpose
    if(quad) return( crossprod(sqrt(ev_p) * sapply(g_evec, as.numeric)) )

    # scale by inverse covariance eigen-values and transform back
    g_prod = g_evec |>
      lapply(\(x) tcrossprod(fac[['y']][['vectors']], t(ev_p * x))) |>
      sapply(\(x) tcrossprod(x, fac[['x']][['vectors']]) )

    return( g_prod )
  }

  # check for expected class in fac, or compute it
  if( !is.null(fac) )
  {
    if( method == 'chol' & !is.matrix(fac) ) stop('Cholesky factor not found in fac')
    if( method == 'eigen' & !is.list(fac) ) stop('eigendecomposition (list) not found in fac')

  } else { fac = pkern_var(g_obs, pars=pars, scaled=TRUE, method=method) }

  # code below handles non-separable case
  if( is.null(z) ) stop('data vector not found g_obs')

  # computation via Cholesky factor (p=-1)
  if( method == 'chol' )
  {
    # solve for (t(C))(V_inv)(z) whose inner product is the quadratic form
    g_chol = forwardsolve(fac, z/pars[['psill']] )
    if(quad) return( pars[['psill']] * crossprod(g_chol) )
    g_prod = backsolve(t(fac), g_chol)
  }

  # computation via eigen-decomposition
  if( method == 'eigen' )
  {
    ev_sqrt = ( pars[['psill']] * fac[['values']] )^(p/2)
    if( any( is.infinite(ev_sqrt) ) ) stop('ill-conditioned covariance matrix (0 eigenvalue)')
    if( any( is.nan(ev_sqrt) ) ) stop('ill-conditioned covariance matrix (eigenvalue < 0)')
    g_evec = ev_sqrt * crossprod(fac[['vectors']], z)
    if(quad) return( crossprod(g_evec, g_evec) )

    # left-multiply by the remaining terms (transpose of the original transform)
    g_prod = tcrossprod(fac[['vectors']], t( ev_sqrt * g_evec ))
  }

  return( g_prod )
}


#' Efficiently compute yzx for symmetric Toeplitz matrices y and x
#'
#' This computes the product y %*% z (and, optionally, right-multiplies the result
#' by x) for symmetric Toeplitz matrices y (and x). It uses fast Fourier transforms to
#' reduce the memory footprint of computations, particularly when z is sparse.
#'
#' By default z is set to the identity matrix, so that pkern_toep_mult(y) returns the
#' matrix y. For vector input y, it returns the Toeplitz matrix generated by y, similar
#' to `base::toeplitz`.
#'
#' The function only requires the first row of y (and x). This is embedded in a
#' larger (zero-padded) vector representing a circulant matrix, whose action on a
#' zero-padded version of z is equivalent to element-wise product in Fourier space.
#' This allows the function to compute the desired matrix product without explicitly
#' forming y or x in memory.
#'
#' The function is optimized for grid data z that are sparse (many zeros). Before
#' computing any transformations it first scans for and removes columns and rows of
#' z which are all zero (replacing them afterwards).
#'
#' @param y numeric matrix or vector, the symmetric Toeplitz matrix y or its first row/column
#' @param z numeric matrix or vector with dimensionality matching `y` (and 'x')
#' @param x numeric matrix or vector, the symmetric Toeplitz matrix x or its first row/column
#'
#' @return numeric matrix, the product of yzx or yz (if x is NULL)
#' @export
#'
#' @examples
#' # define example 1D exponential variogram
#' n = 10
#' y = exp(1-seq(n))
#' y_mat = pkern_toep_mult(y)
#' ( y_mat - stats::toeplitz(y) ) |> abs() |> max()
#'
#' # multiply by random matrix and compare with default matrix multiply
#' z = rnorm(n^2) |> matrix(n)
#' result_default = y_mat %*% z
#' abs( result_default - pkern_toep_mult(y_mat, z) ) |> max()
#'
#' # save memory by passing only the first row of the Toeplitz matrix
#' abs( result_default - pkern_toep_mult(y, z) ) |> max()
#'
#' # sparsify z and repeat
#' idx_sparse = sample.int(n^2, n^2 - n)
#' z[idx_sparse] = 0
#' result_default = y_mat %*% z
#' abs( result_default - pkern_toep_mult(y, z) ) |> max()
#'
#' # right-multiply with another kernel
#' x = exp( 2 *( 1-seq(n) ) )
#' x_mat = pkern_toep_mult(x)
#' result_default = result_default %*% x_mat
#' abs( result_default - pkern_toep_mult(y, z, x) ) |> max()
#'
#' # z can also be supplied as vector of nonzero grid values
#' idx_obs = which(z != 0)
#' gdim = c(y=n, x=n)
#' abs( result_default - pkern_toep_mult(y, z=z[idx_obs], x, idx_obs, gdim) ) |> max()
#'
pkern_toep_mult = function(y, z=NULL, x=NULL, idx_obs=NULL, gdim=NULL)
{
  # convert (Toeplitz) matrix input y to vector
  if( is.matrix(y) ) y = y[1,] |> as.numeric()

  # find amount of zero padding needed to get next highest composite dimension
  n_y = length(y)
  n_pad = nextn(n_y) - n_y
  z_pad = rep(0, n_y + n_pad)

  # indexing and normalization constant for inverse fft
  z_idx = seq(n_y)
  norm_constant = 2*n_y + n_pad

  # set default for z (the identity)
  if( is.null(z) ) z = rep(1, n_y) |> diag()

  # reconstruct full grid with zeros for missing data
  if( !is.null(idx_obs) )
  {
    if( is.null(gdim) ) stop('gdim must be supplied with idx_obs')

    # initialize output vector with zeros then copy nonzero data
    z_nz = z
    z = matrix(0, gdim['y'], gdim['x'])
    z[idx_obs] = z_nz
  }

  # omit zero columns from computation
  if( !is.matrix(z) ) z = matrix(z, n_y)
  is_col_empty = colSums(abs(z)) == 0
  n_compute = ncol(z) - sum(is_col_empty)

  # initialize zero-padded z matrix and copy input data
  z_pad = matrix(numeric(1), norm_constant, n_compute)
  z_pad[z_idx,] = z[, !is_col_empty]

  # discrete FFT of circulant matrix containing y
  fy = stats::fft( c(y, rep(0, 1 + n_pad), y[n_y:2]) )

  # add zero padding, transform (by column), multiply, then transform back
  z[, !is_col_empty] = Re( stats::mvfft(fy * stats::mvfft(z_pad), inverse=TRUE)[z_idx,] )
  if( is.null(x) ) return(z/norm_constant)

  # do right-multiplication by transposing output z and passing back to pkern_toep_mult
  return( t( pkern_toep_mult(x, t(z/norm_constant)) ) )
}

