# pkern_model.R
# Dean Koch, 2022
# Functions for modeling with pkern


#' Compute likelihood function value for covariance parameters `pars` given data `g_obs`
#'
#' Returns the log-likelihood of model parameter list `pars`, given data grid `g_obs`.
#' This is equal to `-log( 2 * pi ) - ( 1/2 ) * ( log_det + quad_form )`, where `log_det`
#' is the log-determinant of covariance matrix V, and `quad_form` is z^T V^{-1} z, with
#' z equal to the observed data vector minus the model mean.
#'
#' If the model has a known trend (mean) vector, supply it in `X`. A spatially constant mean
#' can be supplied as (length-1) numeric. If `X` is a matrix, it is interpreted as a data
#' matrix of predictors for the non-NA values in `g_obs`; the function uses GLS to estimate
#' the mean. If `X=NA` it estimates a spatially constant mean.
#'
#' When `more=TRUE`, the function returns a list containing a count of the number of
#' observations, the likelihood function value, and its two major components; the
#' log-determinant `log_det`, and the quadratic form `quad_form`.
#'
#' `method` specifies how to factorize V, either using the Cholesky factor ('chol')
#' or eigen-decomposition ('eigen'). A pre-computed factorization `fac` can be supplied by
#' calling first `pkern_var` with the same `method`.
#'
#' @param pars list of form returned by `pkern_pars` (with entries 'y', 'x', 'eps', 'psill')
#' @param g_obs list of form returned by `pkern_grid` (with entries 'gdim', 'gres', 'gval')
#' @param X numeric, vector, matrix, or NA, a fixed mean value, or matrix of linear predictors
#' @param method character, the factorization to use: 'chol' (default) or 'eigen'
#' @param fac matrix or list, (optional) pre-computed covariance factorization
#' @param more logical, indicates to return list with likelihood components
#' @param quiet logical indicating to suppress console output
#'
#' @return numeric, the likelihood of `pars` given `g_obs`
#' @export
#'
#' @examples
#' # set up example grid, covariance parameters
#' gdim = c(25, 12)
#' n = prod(gdim)
#' g_obs = pkern_grid(gdim)
#' pars = modifyList(pkern_pars(g_obs, 'gau'), list(psill=0.7, eps=5e-2))
#'
#' # generate some covariates and complete data
#' n_betas = 3
#' betas = rnorm(n_betas)
#' X_all = cbind(1, matrix(rnorm(n*(n_betas-1)), n))
#' g_obs[['gval']] = as.vector( pkern_sim(g_obs) + (X_all %*% betas) )
#' z = g_obs[['gval']]
#'
#' # two methods for likelihood
#' LL_chol = pkern_LL(pars, g_obs, method='chol')
#' LL_eigen = pkern_LL(pars, g_obs, method='eigen')
#'
#' # compare to working directly with matrix inverse
#' V = pkern_var(g_obs, pars, method='none', sep=FALSE)
#' V_inv = chol2inv(chol(V))
#' g_trans = crossprod(V_inv, z)
#' quad_form = as.numeric( t(z) %*% g_trans )
#' log_det = as.numeric( determinant(V, logarithm=TRUE) )[1]
#' LL_naive = (-1/2) * ( n * log( 2 * pi ) + log_det + quad_form )
#'
#' # relative errors
#' abs( LL_naive - LL_chol ) / max(LL_naive, LL_chol)
#' abs( LL_naive - LL_eigen ) / max(LL_naive, LL_eigen)
#'
#'
#' # repeat with most data missing
#' n_obs = 50
#' idx_obs = sort(sample.int(n, n_obs))
#' z_obs = g_obs$gval[idx_obs]
#' g_miss = modifyList(g_obs, list(gval=rep(NA, n)))
#' g_miss[['gval']][idx_obs] = z_obs
#' LL_chol = pkern_LL(pars, g_miss, method='chol')
#' LL_eigen = pkern_LL(pars, g_miss, method='eigen')
#'
#' # working with matrix inverse
#' V = pkern_var(g_miss, pars, method='none')
#' z_mat = matrix(z_obs, ncol=1)
#' V_inv = chol2inv(chol(V))
#' g_trans = (V_inv %*% z_mat)
#' quad_form = t(z_mat) %*% g_trans
#' log_det = as.numeric( determinant(V, logarithm=TRUE) )[1]
#' LL_naive = (-1/2) * ( n_obs * log( 2 * pi ) + log_det + quad_form )
#' abs( LL_naive - LL_chol ) / max(LL_naive, LL_chol)
#' abs( LL_naive - LL_eigen ) / max(LL_naive, LL_eigen)
#'
#' # copy covariates (don't pass the intercept column in X)
#' X = X_all[idx_obs, -1]
#'
#' # use GLS to de-trend, with and without covariatea
#' g_detrend = g_detrend_X = g_miss
#' g_detrend[['gval']][idx_obs] = z_obs - pkern_GLS(g_miss, pars)
#' g_detrend_X[['gval']][idx_obs] = z_obs - pkern_GLS(g_miss, pars, X, out='z')
#'
#' # pass X (or NA) to pkern_LL to do this automatically
#' LL_detrend = pkern_LL(pars, g_detrend)
#' LL_detrend_X = pkern_LL(pars, g_detrend_X)
#' LL_detrend - pkern_LL(pars, g_miss, X=NA)
#' LL_detrend_X - pkern_LL(pars, g_miss, X=X)
#'
#' # equivalent sparse input specification
#' idx_grid = match(seq(n), idx_obs)
#' g_sparse = modifyList(g_obs, list(gval=matrix(z_obs, ncol=1), idx_grid=idx_grid))
#' LL_chol - pkern_LL(pars, g_sparse)
#' LL_eigen - pkern_LL(pars, g_sparse)
#' LL_detrend - pkern_LL(pars, g_sparse, X=NA)
#' LL_detrend_X - pkern_LL(pars, g_sparse, X=X)
#'
#'
#' # repeat with complete data
#'
#' # (don't pass the intercept column in X)
#' X = X_all[,-1]
#' LL_X_chol = pkern_LL(pars, g_obs, X=X)
#' LL_X_eigen = pkern_LL(pars, g_obs, method='eigen', X=X)
#' z_obs = g_obs$gval[!is.na(g_obs$gval)]
#' z_mat = matrix(z_obs, ncol=1)
#' V = pkern_var(g_obs, pars, sep=FALSE)
#' V_inv = chol2inv(chol(V))
#' X_tilde_inv = crossprod(crossprod(V_inv, X_all), X_all) |> chol() |> chol2inv()
#' g_trans = (V_inv %*% z_mat)
#' betas_gls = X_tilde_inv %*% crossprod(X_all, g_trans)
#' z_gls = z_mat - (X_all %*% betas_gls)
#' z_gls_trans = crossprod(V_inv, z_gls)
#' quad_form = as.numeric( t(z_gls) %*% z_gls_trans )
#' log_det = as.numeric( determinant(V, logarithm=TRUE) )[1]
#' LL_naive = (-1/2) * ( n * log( 2 * pi ) + log_det + quad_form )
#' abs( LL_naive - LL_X_chol ) / max(LL_naive, LL_X_chol)
#' abs( LL_naive - LL_X_eigen ) / max(LL_naive, LL_X_eigen)
#'
#' # return components of likelihood with more=TRUE
#' LL_result = pkern_LL(pars, g_obs, X=X, more=TRUE)
#' LL_result$LL - LL_X_chol
#' LL_result$q - quad_form
#' LL_result$d - log_det
#' LL_result$n - n
#'
#'
pkern_LL = function(pars, g_obs, X=0, method='chol', fac=NULL, quiet=TRUE, more=FALSE)
{
  # set flag for GLS and default one layer count
  if(is.data.frame(X)) X = as.matrix(X)
  use_GLS = is.matrix(X) | anyNA(X)

  # multi-layer support
  n_layer = 1
  is_multi = !is.null(g_obs[['idx_grid']])
  if(is_multi)
  {
    # coerce vector to 1-column matrix and identify non-NA points
    if( !is.matrix(g_obs[['gval']]) ) g_obs[['gval']] = matrix(g_obs[['gval']], ncol=1L)
    is_obs = !is.na(g_obs[['idx_grid']])

    # reorder the non-NA data matrix to grid-vectorized order
    reorder_z = g_obs[['idx_grid']][is_obs]
    n_layer = ncol(g_obs[['gval']])

    # matrix(.., ncol) avoids R simplifying to vector in 1 column case
    z = matrix(g_obs[['gval']][reorder_z,], ncol=n_layer)
    n_obs = nrow(z)

  } else {

    # single layer mode - copy non-NA data
    is_obs = as.vector(!is.na(g_obs[['gval']]))
    z = matrix(g_obs[['gval']][is_obs], ncol=1L)
    if( is.null(z) ) stop('No non-NA values found in g_obs')
    n_obs = length(z)
  }

  # complete data case triggers separability methods, which require eigen
  is_complete = all(is_obs)
  if(is_complete) method = 'eigen'

  # compute factorization (scaled=TRUE means partial sill is factored out)
  if( is.null(fac) ) fac = pkern_var(g_obs, pars=pars, scaled=TRUE, method=method)

  # GLS estimate of mean based on predictors in X
  if( use_GLS ) X = pkern_GLS(g_obs, pars, X=X, fac=fac, method=method, out='z')

  # matricize scalar and vector input to X
  if( !is.matrix(X) ) X = matrix(X, ncol=1L)
  if( nrow(X) == 1 ) X = matrix(rep(X, n_obs), ncol=ncol(X))

  # reorder X to match z
  if(is_multi) X = X[reorder_z,]

  # matrix of de-trended Gaussian random vectors to evaluate
  z_centered = matrix(z-X, ncol=n_layer)

  # Cholesky factor method is fastest
  if( method == 'chol' )
  {
    # check for bad fac input (it should be matrix t(C), where C is the output of chol)
    if( !is.matrix(fac) ) stop('Cholesky factor (matrix) not found in fac')

    # get quadratic form of de-trended variable(s)
    quad_form = apply(z_centered, 2, function(z_i) {
      as.numeric(pkern_var_mult(z_i, pars, fac=fac, method=method, quad=TRUE))
    })

    # determinant is the squared product of the diagonal elements in Cholesky factor t(C)
    log_det_corr = sum(log(diag(fac)))

    # partial sill was factored out with scaled=TRUE, so it is multiplied back in here
    log_det = n_obs * log(pars[['psill']]) + 2 * log_det_corr
  }

  # eigen-decomposition method
  if( method == 'eigen' )
  {
    # check for bad fac input (it should be the list output of eigen)
    if( !is.list(fac) ) stop('eigendecomposition (list) not found in fac')

    # get quadratic form of de-trended variable(s)
    quad_form = apply(z_centered, 2, function(z_i) {
      as.numeric(pkern_var_mult(z_i, pars, fac=fac, method=method, quad=TRUE))
    })

    # determinant is product of the eigenvalues
    if( !is_complete )
    {
      # partial sill was factored out earlier so we multiply it back in (on log scale)
      log_det = n_obs * log(pars[['psill']]) + sum(log(fac[['values']]))

    } else {

      # separable case: eigenvalues are a kronecker product plus diagonal nugget effect
      ev_scaled = kronecker(fac[['y']][['values']], fac[['x']][['values']])
      log_det = sum(log(pars[['eps']] + pars[['psill']] * ev_scaled))
    }
  }

  # compute log likelihood, print to console then return
  log_likelihood = (-1/2) * ( n_layer*( n_obs * log( 2 * pi ) + log_det ) + sum(quad_form) )
  if( !quiet ) cat( paste(round(log_likelihood, 5), '\n') )
  if(more) return(list(LL=log_likelihood, q=quad_form, d=log_det, n_obs=n_obs))
  return(log_likelihood)
}


#' Compute negative log-likelihood for parameter vector `p`
#'
#' Returns the negative log-likelihood of parameter vector `p` for the covariance
#' model `pars_fix`, given data grid `g_obs`.
#'
#' This is a wrapper for `-pkern_LL()` allowing parameters to be passed as a numeric
#' vector instead of a list (for use in optimization etc). parameters in `p` are copied
#' to `pars_fix` and passed to the likelihood computer.
#'
#' `p` is the vector of covariance parameters to test. It must be in the form (length,
#' order) accepted by `pkern_pars_update(pars_fix, p, na_omit=TRUE, iso=iso)`.
#'
#' @param p numeric vector of covariance parameters accepted by `pkern_pars_update`
#' @param g_obs list of form returned by `pkern_grid` (with entries 'gdim', 'gres', 'gval')
#' @param pars_fix list of form returned by `pkern_pars` (with entries 'y', 'x', 'eps', 'psill')
#' @param X numeric, vector, matrix, or NA, the mean or its linear predictors, passed to `pkern_LL`
#' @param iso logical, indicates to use identical kernels for x and y (`pars$x` is ignored)
#' @param quiet logical indicating to suppress console output
#'
#' @return numeric, the negative log-likelihood of `p` given `g_obs`
#' @export
#'
#' @examples
#' # set up example grid and data
#' g_obs = pkern_grid(10)
#' g_obs$gval = rnorm(10^2)
#'
#' # get some default parameters and vectorize them
#' pars_fix = pkern_pars(g_obs, 'gau')
#' p = pars_fix |> pkern_pars_update()
#' pkern_nLL(p, g_obs, pars_fix)
#' pkern_nLL(p, g_obs, pars_fix, eps_scaled=TRUE)
#'
pkern_nLL = function(p, g_obs, pars_fix, X=0, iso=FALSE, quiet=TRUE, log_scale=FALSE)
{
  # transform back from log scale
  if(log_scale)
  {
    # parameter vector and parameter list
    p = exp(p)
    pars_fix = pkern_pars_update(pars_fix, exp(pkern_pars_update(pars_fix)) )
  }

  # update covariance parameter list with new values then vectorize it
  pars_complete = pkern_pars_update(pars_fix, p, iso=iso, na_omit=TRUE)
  p_complete = pkern_pars_update(pars_complete)

  # print complete parameters vector then compute likelihood
  if(!quiet) cat( paste(paste(format(p_complete), collapse=', '), ' :: LL = ') )
  return(-pkern_LL(pars=pars_complete, g_obs=g_obs, X=X, quiet=quiet))
}


#' Generate random draw from multivariate normal distribution for grid points
#'
#' @param g any object accepted or returned by `pkern_grid`
#' @param pars list, covariance parameters in form returned by `pkern_pars`
#' @param fac (optional) list, precomputed factorization of component correlation matrices
#' @param add_frac numeric > 0, small constant for adjusting negative eigen-values (see details)
#' @param quiet logical, indicating to suppress warnings about negative eigen-values
#'
#' @return numeric vector, the vectorized grid data
#' @export
#'
#' @examples
#'
#' # example grid and covariance parameters
#' gdim = c(25,15)
#' g = pkern_grid(gdim)
#' pars = pkern_pars(g, 'mat')
#'
#' # this example has a large nugget effect
#' gval = pkern_sim(g, pars)
#' pkern_plot(matrix(gval, gdim))
#'
#' # plot with yx coordinates
#' g_sim = modifyList(g, list(gval=gval))
#' pkern_plot(g_sim)
#'
#' # repeat without nugget effect for smooth field
#' pars0 = modifyList(pars, list(eps=0))
#' gval = pkern_sim(g, pars0)
#' pkern_plot(matrix(gval, gdim))
#'
#' # plot with yx coordinates
#' g_sim = modifyList(g, list(gval=gval))
#' pkern_plot(g_sim)
#'
#' # parameters are automatically assigned if not supplied
#' pkern_pars(g) # default is gau x gau with ranges roughly equal to grid side lengths
#' pkern_plot(matrix(pkern_sim(g), gdim))
#'
pkern_sim = function(g_obs, pars=pkern_pars(g_obs), fac=NULL, add_frac=1e-16, quiet=FALSE)
{
  gdim = g_obs[['gdim']]
  n = prod(gdim)

  # eigen-decompositions of separable components of full grid correlation matrix
  g_empty = modifyList(g_obs, list(gval=NULL))
  if(is.null(fac)) fac = pkern_var(g_empty, pars, method='eigen')
  is_ev_negative = lapply(fac, function(eig) !( eig[['values']] > 0 ) )

  # warn of any non-positive eigenvalues and fix them
  is_dim_invalid = sapply(is_ev_negative, any)
  if( !quiet & any(is_dim_invalid) ) warning('component covariance had negative eigenvalue(s)')
  for( nm_dim in names(is_dim_invalid)[is_dim_invalid] )
  {
    nz_constant = add_frac * min( fac[[nm_dim]][['values']][ !is_ev_negative[[nm_dim]] ] )
    fac[[nm_dim]][['values']][ is_ev_negative[[nm_dim]] ] = nz_constant
  }

  return( as.vector(pkern_var_mult(rnorm(n), pars, fac=fac, method='eigen', p=1/2)) )
}

