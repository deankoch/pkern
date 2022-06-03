# pkern_estim.R
# Dean Koch, 2022
# Functions for parameter inference and spatial prediction

#' Generalized least squares (GLS) estimator
#'
#' Computes coefficients of the linear predictor E(Z) = Xb using the GLS equation
#'
#'  b = ( X^T V^{-1} X )^{-1} X^T V^{-1} z
#'
#' where V is the covariance matrix for data z, and X is a matrix of covariates. If
#' `X=NA`, X in the above is replaced by an intercept column (a vector of 1's) and this
#' expression is a length-1 numeric (the estimated spatially constant mean). Otherwise it
#' is a vector of coefficients whose first entry is the mean and subsequent entries are
#' coefficients to multiply the columns of the predictor matrix `X`
#'
#' If supplied, the columns of matrix `X` should be independent, and the number of rows should
#' match the number of observed grid points (non-`NA` values in `g_obs$gval`). An intercept
#' column is appended to this linear predictor matrix before evaluating the GLS expression
#' (without checking if it's there already). Do not include an intercept column in
#' argument `X`!
#'
#' When `out='z'`, the function multiplies the GLS coefficients by `X` (with appended
#' intercept column) to get the linear predictor for E(z)
#'
#' Note that the factorization returned in fac_X is for the model scaled by partial sill
#'
#' @param g_obs list of form returned by `pkern_grid` (with entries 'gdim', 'gres', 'gval')
#' @param pars list of form returned by `pkern_pars` (with entries 'y', 'x', 'eps', 'psill')
#' @param X matrix or NA, the linear predictors (in columns) excluding intercept
#' @param method character, the factorization to use: 'chol' (default) or 'eigen'
#' @param fac matrix or list, (optional) pre-computed covariance matrix factorization
#' @param out character, either 'c' (coefficients) or 'z' (linear predictor)
#'
#' @return numeric vector of length equal to one plus the number of columns in `X`
#' @export
#'
#' @examples
#' # set up example grid, and covariance parameters
#' gdim = c(45, 31)
#' n = prod(gdim)
#' g = pkern_grid(gdim)
#' pars = modifyList(pkern_pars(g, 'gau'), list(psill=12))
#'
#' # generate spatial noise
#' z = pkern_sim(g, pars, quiet=TRUE)
#' pkern_plot(modifyList(g, list(gval=z)))
#'
#' # generate some covariates and data
#' n_betas = 3
#' betas = rnorm(n_betas, 0, 10)
#' X_all = cbind(1, matrix(rnorm(n*(n_betas-1)), n))
#' g_obs = modifyList(g, list(gval=z + (X_all %*% betas)))
#'
#' # exclude intercept column in calls to pkern_GLS
#' X = X_all[,-1]
#'
#' # find the GLS coefficients
#' betas_est = pkern_GLS(g_obs, pars, X)
#' print(betas_est)
#' print(betas)
#'
#' # compute trend as product of betas with matrix X_all, or by setting out='z'
#' max(abs( pkern_GLS(g_obs, pars, X, out='z') - (X_all %*% betas_est) ))
#'
#' # missing data example
#' n_obs = 10
#' idx_rem = sort(sample.int(n, n-n_obs))
#' g_miss = g_obs
#' g_miss$gval[idx_rem] = NA
#' pkern_plot(g_miss)
#' betas_est = pkern_GLS(g_miss, pars, X)
#' print(betas_est)
#' print(betas)
#'
pkern_GLS = function(g_obs, pars, X=NA, fac=NULL, method='chol', out='b')
{
  # copy non-NA data
  is_obs = as.vector(!is.na(g_obs[['gval']]))
  if( is.null(g_obs[['gval']][is_obs]) ) stop('No non-NA values found in g_obs')
  z = g_obs[['gval']][is_obs]
  n_obs = sum(is_obs)
  n = length(g_obs[['gval']])

  # compute V factorization
  if( all(is_obs) ) method = 'eigen'
  if( is.null(fac) )
  {
    # check for sub-grid structure and substitute simpler equivalent problem if possible
    sg = pkern_sub_find(is_obs, g_obs[['gdim']])
    is_sg = !is.null(sg)
    if( is_sg )
    {
      # extract sub-grid layout and find separable covariance eigen-decomposition
      method = 'eigen'
      g_obs = list(gdim=sg[['gdim']], gres=g_obs[['gres']] * sg[['res_scale']])
      fac = pkern_var(g_obs, pars=pars, scaled=TRUE, method=method)

      # g_obs should have no missing (NA) data points now
      g_obs[['gval']] = z
      is_obs = rep(TRUE, n_obs)
      n = n_obs

    } else {

      # compute factorization (scaled=TRUE means partial sill is factored out)
      fac = pkern_var(g_obs, pars=pars, scaled=TRUE, method=method)
    }
  }

  # build matrix of covariate values from intercept column and (optionally) X
  if( !anyNA(X) )
  {
    # check for invalid input to X
    if( !is.matrix(X) ) stop('X must be a matrix of covariate values')
    msg_mismatch = 'nrow(X) did not match length(g_obs$gval) or the number of non-NA points'
    if( anyNA(X) ) { X = NULL } else { if( !( nrow(X) %in% c(n_obs, n) ) ) stop(msg_mismatch) }

    # handle case of observed data only
    if( nrow(X) == n_obs ) is_obs = rep(TRUE, n_obs)
    X = cbind(1, X)
    X_obs = matrix(X[is_obs,], nrow=n_obs)

  } else { X = X_obs = matrix(1, nrow=n_obs) }

  # return data matrix
  if(startsWith(out, 'x')) return(X_obs)

  # find the factorization of quadratic form with X (scaling by V inverse)
  fac_X = pkern_var(g_obs, pars, X=X_obs, method=method, scaled=TRUE, fac=fac)

  # compute GLS coefficients using whitened observation data
  z_trans = pkern_var_mult(z, pars, fac=fac, method=method)
  betas_gls = pkern_var_mult(t(crossprod(z_trans, X_obs)), pars, fac=fac_X, method=method)

  # return betas if requested
  if(startsWith(out, 'b')) return(as.numeric(betas_gls))

  # return linear predictor
  z_gls = as.numeric( tcrossprod(X, t(betas_gls)) )
  if(startsWith(out, 'z')) return(z_gls)

  # return everything in a list
  return(list(b = as.numeric(betas_gls), z = z_gls, x = X, fac_X = fac_X))
}


#' Fit a covariance model to the data by maximum likelihood
#'
#' An automated model fitting procedure
#'
#' documentation unfinished
#'
#' @param g_obs
#' @param pars
#' @param X
#'
#' @return
#' @export
#'
#' @examples
#'
#' # define a grid
#' gdim = c(50, 53)
#' g_empty = pkern_grid(gdim)
#' pars_src = pkern_pars(g_empty)
#' pars_src = modifyList(pars_src, list(eps=runif(1, 0, 1e1), psill=runif(1, 0, 1e2)))
#' pars_src[['y']][['kp']] = pars_src[['x']][['kp']] = runif(1, 1, 50)
#'
#' # generate example data and fit to it
#' gval = pkern_sim(g_empty, pars_src, quiet=F)
#' g_obs = modifyList(g_empty, list(gval=gval))
#' pkern_plot(g_obs)
#' fit_result = pkern_fit(g_obs, pars='gau')
#'
#' fit_result$pars |> pkern_pars_update()
#' pars_src |> pkern_pars_update()
#'
#' # check sequence of other psill values
#' pars_out = fit_result$pars
#' psill_test = ( 2^(seq(5) - 3) ) * pars_out[['psill']]
#' LL_test = sapply(psill_test, function(s) pkern_LL(modifyList(pars_out, list(psill=s)), g_obs) )
#' plot(psill_test, LL_test)
#' lines(psill_test, LL_test)
#' print(data.frame(psill=psill_test, likelihood=LL_test))
#'
#' # repeat with most data missing
#' n = prod(gdim)
#' n_obs = 200
#' gval = pkern_sim(g_empty, pars_src, quiet=TRUE)
#' g_obs = modifyList(g_empty, list(gval=gval))
#' idx_obs = sample.int(prod(gdim), n_obs)
#' g_miss = modifyList(g_obs, list(gval=rep(NA, n)))
#' g_miss$gval[idx_obs] = g_obs$gval[idx_obs]
#' pkern_plot(g_miss)
#'
#'
pkern_fit = function(g_obs, pars=NULL, X=NA, iso=TRUE, initial=NULL, quiet=FALSE)
{
  # unpack grid as list
  g = pkern_grid(g_obs)
  gdim = g[['gdim']]
  is_obs = !is.na(g[['gval']])
  if( !any(is_obs) ) stop('no data found in g_obs')

  n_all = prod(gdim)
  n_obs = sum(is_obs)

  # coerce to matrix of appropriate dimensions
  # if( is.matrix(X) )
  # {
  #   # pad with NAs as needed
  #   #idx_pad = match(seq(n_all), which(is_obs))
  #   # if( is.vector(X) ) if( length(X) == n_obs ) X = X[idx_pad]
  #   # if( is.matrix(X) ) if( nrow(X) == n_obs ) X = apply(X, 2, function(x) x[idx_pad])
  #   # X = matrix(X, nrow=n_all)
  #
  #   if( is.vector(X) ) if( length(X) == n_obs ) X = X[idx_pad]
  #   X = matrix(X, nrow=n_obs)
  #
  #   if( all(is.na(X)) ) X = NA
  # }

  # substitute equivalent sub-grid problem if possible
  sub_result = pkern_sub_find(g)
  if( !is.null(sub_result) )
  {
    # skip when full grid is non-NA
    if( any(sub_result[['gdim']] < gdim) )
    {
      # compute new grid configuration and copy to g
      gdim = sub_result[['gdim']]
      gres = g[['gres']] * sub_result[['res_scale']]
      gyx = Map(function(yx, idx) yx[idx], yx=g[['gyx']], idx=sub_result[['ij']])
      g = modifyList(g, list(gdim=gdim, gres=gres, gyx=gyx, gval=g[['gval']][is_obs]))
      is_obs = rep(TRUE, prod(gdim))
    }
  }

  # set covariance parameter defaults
  if( !is.list(pars) ) pars = pkern_pars(g, ifelse(is.null(pars), 'gau', pars))
  p_fixed = pkern_pars_update(pars, iso=TRUE)
  is_fitted = is.na(p_fixed)
  if( !any(is_fitted) ) is_fitted[] = TRUE

  # set initial value defaults
  nm_fitted = names(is_fitted)[is_fitted]
  nm_fixed = names(is_fitted)[!is_fitted]
  if( is.null(initial) ) initial = pkern_bds(pars, g, rows=nm_fitted)[['initial']]
  pars = pkern_pars_update(pars, p_fixed, iso=TRUE)

  # fit the model
  #v = var(g[['gval']], na.rm=TRUE)
  result_optim = pkern_optim(g_obs=g, pars_fix=pars, X=X, iso=TRUE, initial=initial, quiet=quiet)
  pars_out = result_optim[['pars']]

  # # check sequence of likely psill substitutions
  # psill_test = ( 2^(seq(5) - 3) ) * pars_out[['psill']]
  # pars_test = lapply(psill_test, function(s) modifyList(pars_out, list(psill=s)) )
  # LL_test = sapply(pars_test, function(p) pkern_LL(p, g, X=X) )
  # pars_out[['psill']] = psill_test[ which.max(LL_test) ]

  # de-trend the data by subtracting GLS estimate
  z_gls = 0
  if( anyNA(X) | is.matrix(X) ) z_gls = pkern_GLS(g, pars_out, X=X, out='z')
  g[['gval']][is_obs] = g[['gval']][is_obs] - z_gls

  # plot the semi-variogram for de-trended data
  vg_detrend = pkern_sample_vg(g)
  pkern_plot_semi(vg_detrend, pars_out)
  return(result_optim)


  if(0)
  {

  # scale the data
  z_std = scale(g[['gval']])
  z_centre = attr(z_std, 'scaled:center')
  z_scale = attr(z_std, 'scaled:scale')
  g = modifyList(g, list(gval=as.vector(z_std)))

  # variance parameters must also be scaled
  is_v_fitted = nm_fitted %in% c('eps', 'psill')
  is_v_fixed = nm_fixed %in% c('eps', 'psill')
  initial[is_v_fitted] = initial[is_v_fitted] / (2*z_scale^2)
  if( any(is_v_fixed) ) p_fixed[nm_fixed][is_v_fixed] = p_fixed[nm_fixed][is_v_fixed] / (2*z_scale^2)
  pars = pkern_pars_update(pars, p_fixed, iso=TRUE)

  # fit the model
  if( anyNA(X) ) X = NA
  result_optim = pkern_optim(g, pars_fix=pars, X=X, iso=TRUE, initial=initial)
  pars_out = result_optim[['pars']]

  # check sequence of likely psill values
  psill_test = ( 2^(seq(5) - 3) ) * pars_out[['psill']]
  pars_test = lapply(psill_test, function(s) modifyList(pars_out, list(psill=s)) )
  LL_test = sapply(pars_test, function(p) pkern_LL(p, g, X=X) )
  #print(z_scale)
  #print(data.frame(psill=psill_test, likelihood=LL_test))
  pars_out[['psill']] = psill_test[ which.max(LL_test) ]

  # de-trend the scaled data
  z_gls = pkern_GLS(g, pars_out, X=X, out='z')
  if(anyNA(X)) z_gls = z_gls[1]
  g_out = modifyList(g, list(gval = z_scale * (z_std-z_gls) ) )

  # transform scaled variance parameters back to scale of input data
  v_unscaled = list(eps = 2*z_scale^2 * pars_out[['eps']], psill = 2*z_scale^2 * pars_out[['psill']])
  pars_out = modifyList(pars_out, v_unscaled)

  # plot the semi-variogram on original scale
  vg_out = pkern_sample_vg(g_out)
  pkern_plot_semi(vg_out, pars_out)
  return(modifyList(result_optim, list(pars=pars_out)))
  }
}



#' Fit covariance parameters to data by maximum (profile) likelihood using optim
#'
#' documentation unfinished
#'
#' @param g_obs list of form returned by `pkern_grid` (with entries 'gdim', 'gres', 'gval')
#' @param pars_fix list of fixed kernel parameters, with NAs indicating parameters to fit
#' @param X numeric, vector, matrix, or NA, the mean or its linear predictors, passed to `pkern_LL`
#' @param iso logical, indicating to constrain the y and x kernel parameters to be the same
#' @param control list, passed to `stats::optim`
#' @param quiet logical, indicating to suppress console output
#' @param ... named arguments to pass to `pkern_bds`
#'
#' @return
#' @export
#'
#' @examples
#' # set up example grid and data
#' g_obs = pkern_grid(10)
#' g_obs$gval = rnorm(10^2)
#' pkern_optim(g_obs, quiet=TRUE)
#'
#' # repeat with one or more parameters fixed
#' pars_fix = pkern_pars_make('gau') # NA parameters list
#' pars_fix$psill = 1
#' pkern_optim(g_obs, pars_fix, quiet=TRUE)
#' pars_fix$y$kp = 1
#' pkern_optim(g_obs, pars_fix, quiet=TRUE)
#'
#' # iso mode constrains x parameters to equal y parameters
#' pkern_optim(g_obs, iso=T, quiet=TRUE)
#' pkern_optim(g_obs, pars_fix, iso=T, quiet=TRUE)
#'
pkern_optim = function(g_obs, pars_fix='gau', X=0, iso=F, control=list(), quiet=F,
                       log_scale=TRUE, method='L-BFGS-B', ...)
{
  # only L-BFGS-B accepts bounds, so log-scale is mandatory for other methods
  if( (method != 'L-BFGS-B' & !log_scale) )
  {
    warning('Setting log_scale=TRUE (required for methods other than L-BFGS-B)')
    log_scale = TRUE
  }

  # standardize input to pars_fix and set NAs for missing values
  pars_fix = pkern_pars_make(pars_fix)

  # extract parameter names and NA structure supplied in pars_fix
  pars_fix_vec = pkern_pars_update(pars_fix, iso=iso)
  nm_fit = names(pars_fix_vec)[ is.na(pars_fix_vec) ]

  # if all parameters are supplied in pars_fix, fit them all
  if( length(nm_fit) == 0 )
  {
    # set initial value automatically here?
    #initial = pkern_pars_update(pars_fix)
    pars_fix = pkern_pars_update(pars_fix, rep(NA, length(pars_fix_vec)), iso=iso)
    pars_fix_vec = pkern_pars_update(pars_fix, iso=iso)
    nm_fit = names(pars_fix_vec)
  }

  # get default initial values and bounds data frame
  pars_df = pkern_bds(pars_fix, g_obs, rows=nm_fit, ...)
  #pars_df = pkern_bds(pars_fix, g_obs, rows=nm_fit)

  # switch to log-scale if requested
  if(log_scale)
  {
    #bds_all_df = bds_all_df |> log()
    pars_df = log(pars_df)
    pars_fix_vec = log(pars_fix_vec)
  }

  # sacling constants for internal use by optimizer
  optim_scale = apply(pars_df, 1L, function(p) diff(range(p)))

  # TODO: check this
  # when iso=TRUE and pars_fix contains fixed y kernel parameter(s) they must be copied to x
  pars_fix = pkern_pars_update(pars_fix, pars_fix_vec, iso=iso)

  # evaluate objective at initial values and bounds
  bds_nLL = apply(pars_df, 2L, function(p) {
    pkern_nLL(p, g_obs, pars_fix, X, iso, quiet=TRUE, log_scale) })

  # 1d optimization
  if( nrow(pars_df) == 1 )
  {
    tol = control[['tol']]
    optimize_out = optimize(f = pkern_nLL,
                            interval = pars_df[c('lower', 'upper')],
                            tol = ifelse(is.null(tol), .Machine$double.eps^0.25, tol),
                            g_obs = g_obs,
                            pars_fix = pars_fix,
                            X = X,
                            iso = iso,
                            log_scale = log_scale,
                            quiet = quiet)

    optim_out = list(message='', par=optimize_out[['minimum']], value=optimize_out[['objective']])

  } else {

    # n-d optimization
    # set default control parameters and run the optimizer
    control = modifyList(list(maxit=1e3, parscale=optim_scale), control)
    optim_out = stats::optim(par = pars_df[['initial']],
                             lower = pars_df[['lower']],
                             upper = pars_df[['upper']],
                             f = pkern_nLL,
                             method = method,
                             g_obs = g_obs,
                             pars_fix = pars_fix,
                             X = X,
                             iso = iso,
                             quiet = quiet,
                             log_scale = log_scale,
                             control = control) |> suppressWarnings()
  }

  # unpack optimizer output
  obj_val = optim_out[['value']]
  pars_fitted_v = optim_out[['par']]
  if(!quiet) cat(paste('\n', optim_out[['message']]))

  # revert to initial value if optimize/optim result was not an improvement
  if(bds_nLL[['initial']] < obj_val) pars_fitted_v = pars_df[['initial']]
  # TODO: further checking if upper/lower bounds were also better

  # reshape as list and transform back if fitting on log-scale
  pars_fitted = pkern_pars_update(pars_fix, pars_fitted_v, iso=iso, na_omit=TRUE)
  pars_df['fitted'] = pars_fitted_v
  if(log_scale)
  {
    pars_df = exp(pars_df)
    pars_fitted = pkern_pars_update(pars_fitted, exp(pkern_pars_update(pars_fitted)))
  }

  # return parameters list and data frame in a list
  df_order = c('lower', 'initial', 'fitted', 'upper')
  return(list(pars=pars_fitted, df=pars_df[df_order], obj=obj_val))
}


#' Compute ordinary kriging predictor (or variance) for data on a grid
#'
#' Evaluates the ordinary kriging equations in section 3 of Cressie (1993) over the
#' grid defined in `g_obs`. These are the predicted values minimizing mean squared
#' prediction error under the covariance model specified by `pars`.
#'
#' Set `makev=TRUE` to return the pointwise kriging variance. This takes approximately
#' n_obs times longer to evaluate than `makev=FALSE`. A progress bar will be printed to
#' console unless `quiet=TRUE`.
#'
#' The covariance factorization `fac` can be pre-computed using `pkern_var(..., scaled=TRUE)`
#' to speed up repeated calls where only the observed data values change (ie same covariance
#' structure `pars`, and same NA structure in the data). Note that the kriging variance does
#' not change in this case and only needs to be computed once.
#'
#' @param g_obs list of form returned by `pkern_grid` (with entries 'gdim', 'gres', 'gval')
#' @param pars list of form returned by `pkern_pars` (with entries 'y', 'x', 'eps', 'psill')
#' @param X numeric, vector, matrix, or NA: the mean, or its linear predictors
#' @param out character, the return value, one of 'predictor', 'variance', or 'm'
#' @param fac (optional) pre-computed factorization of covariance matrix scaled by partial sill
#' @param quiet logical indicating to suppress console output
#'
#' @return numeric matrix, the predicted values (or their variance)
#' @export
#'
#' @examples
#' # make example grid and data
#' n = 25^2
#' n_obs = 10
#' g_obs = pkern_grid(sqrt(n))
#' idx_obs = sample.int(n, n_obs)
#' g_obs$gval[idx_obs] = rnorm(n_obs)
#' pars = pkern_pars('gau', g_obs)
#' g_pred = pkern_cmean(g_obs, pars)
#' g_var = pkern_cmean(g_obs, pars, makev=TRUE, quiet=TRUE)
#' #g_obs |> pkern_plot()
#' #g_obs |> modifyList(list(gval=g_pred)) |> pkern_plot()
#' #g_obs |> modifyList(list(gval=g_var)) |> pkern_plot()
#'
pkern_cmean = function(g_obs, pars, X=NA, fac=NULL, out='p', method='chol', quiet=FALSE)
{
  # check for expected objects in list in g_obs
  nm_expect = c('gdim', 'gres', 'gval')
  msg_expect = paste('g_obs must be a list with entries named', paste(nm_expect, collapse=', '))
  if( !is.list(g_obs) | !all( nm_expect %in% names(g_obs) ) ) stop(msg_expect)
  gdim = g_obs[['gdim']]
  gres = g_obs[['gres']]

  # copy non-NA data as needed
  is_obs = !is.na(g_obs[['gval']])
  idx_obs = which( is_obs )
  z = g_obs[['gval']][is_obs]
  if( is.null(z) ) stop('No non-NA values found in g_obs')
  n_obs = length(z)
  n = prod(gdim)

  # initialize parameters and check if we are estimating linear predictor for the model
  use_GLS = anyNA(X) | is.matrix(X)
  if( !is.matrix(X) ) mu_GLS = X
  m = rep(0, n)

  # construct first rows of the symmetric Toeplitz correlation matrices for y and x
  cy = cy_obs = sqrt(pars[['psill']]) * pkern_corr_mat(pars[['y']], gdim[['y']], gres[['y']], i=1)
  cx = cx_obs = sqrt(pars[['psill']]) * pkern_corr_mat(pars[['x']], gdim[['x']], gres[['x']], i=1)

  # set up factorization when not supplied. Variance mode forces eigen-decomposition
  if( startsWith(out, 'v') ) method = 'eigen'

  # check for sub-grid structure and substitute simpler equivalent problem if possible
  sg = pkern_sub_find(is_obs, g_obs[['gdim']])
  is_sg = !is.null(sg)
  if( is_sg )
  {
    # extract sub-grid layout and find separable covariance eigen-decomposition
    method = 'eigen'
    g_obs = list(gval=z, gdim=sg[['gdim']], gres=g_obs[['gres']] * sg[['res_scale']])
    if( is.null(fac) ) fac = pkern_var(g_obs, pars=pars, scaled=TRUE, method=method)
    cy_obs = cy[ sg[['ij']][['y']] ]
    cx_obs = cx[ sg[['ij']][['x']] ]

    # g_obs should have no missing (NA) data points now
    is_obs = rep(TRUE, n_obs)

  } else {

    # compute factorization (scaled=TRUE means partial sill is factored out)
    if( is.null(fac) ) fac = pkern_var(g_obs, pars=pars, scaled=TRUE, method=method)
  }

  # transform the observed data by left-multiplying with inverse covariance
  z_tilde = pkern_var_mult(g_obs[['gval']][is_obs], pars, fac=fac)

  # left-multiply by cross-covariance to get simple kriging predictor
  z_p_tilde = pkern_toep_mult(cy, z_tilde, cx, idx_obs, gdim) |> as.numeric()
  if( !use_GLS & startsWith(out, 'p') ) return(as.numeric(X) + z_p_tilde)

  # compute GLS coefficients and resulting adjustments to predictor as needed
  if(use_GLS)
  {
    # make a copy of the observed locations in X
    if( !anyNA(X) ) { X_obs = matrix(X[is_obs,], n_obs) } else {
      X = NULL
      X_obs = NA
    }

    # find betas, predictor, and the data matrix with an intercept column
    gls = pkern_GLS(g_obs, pars, X=X_obs, fac=fac, method=method, out='a')
    fac_X = gls[['fac_X']]

    # compute bias adjustment due to estimation of linear predictor
    X_p_tilde = gls[['x']] |>
      pkern_var_mult(pars, fac=fac) |>
      apply(2, \(x) pkern_toep_mult(cy, x, cx, idx_obs, gdim))
    X_p_tilde[idx_obs,] = gls[['x']]

    # compute trend and 'm'
    X_adj = cbind(rep(1, n), X) - X_p_tilde
    mu_GLS = tcrossprod(X_adj, t(gls[['b']])) |> as.numeric()
    if( startsWith(out, 'm') ) return( as.numeric(pkern_var_mult( t(X_adj), pars, fac=fac_X)) )
  }

  # universal kriging predictor
  if( startsWith(out, 'p') ) return( mu_GLS + z_p_tilde )

  # compute variance contribution from GLS (0 when not using GLS)
  v_gls = numeric(n)
  if(use_GLS)
  {
    # small loop over eigen-values in fac_X, adding each contribution to running total in v_gls
    for(idx in seq_along(fac_X[['values']]))
    {
      # eigen-values of inverse scaled by psill
      ev = 1 / ( pars[['psill']] * fac_X[['values']][idx] )
      v_gls[] = v_gls[] + ev * tcrossprod(fac_X[['vectors']][,idx], X_adj)^2
    }
  }

  # use a more efficient method when observed points form a sub-grid
  idx_ev = seq(n_obs)
  if(is_sg)
  {
    # check that the correct factorization (componentwise, in a list) was supplied
    if( !all( c('y', 'x') %in% names(fac) ) ) stop('supplied factorization was not separable')

    # compute eigenvalues for observed covariance matrix inverse
    ev_corr = kronecker(fac[['x']][['values']], fac[['y']][['values']])
    ev = 1 / ( (pars[['psill']] * ev_corr) + pars[['eps']] )

    # directly build and multiply the relatively small component covariance matrices
    c_cross_y = pkern_corr_mat(pars[['y']], gdim[['y']], gres[['y']], j=sg[['ij']][['y']])
    c_cross_x = pkern_corr_mat(pars[['x']], gdim[['x']], gres[['x']], j=sg[['ij']][['x']])
    add_y2 = ( c_cross_y %*% fac[['y']][['vectors']] )
    add_x2 = ( c_cross_x %*% fac[['x']][['vectors']] )

    # a different ordering for the loop below (largest eigenvalues first)
    idx_ev = order(ev)
  }

  # large loop over eigen-values of covariance matrix, iteratively adding to v_rem
  v_rem = numeric(n)
  if(!quiet) pb = utils::txtProgressBar(max=n_obs, style=3)
  for(idx in seq(n_obs))
  {
    # update progress bar then change to reordered index
    if(!quiet) utils::setTxtProgressBar(pb, idx)
    idx = idx_ev[idx]

    # non-separable case first
    if( !is_sg )
    {
      # eigen-values of inverse are scaled by psill, then slow multiplication with cross covariance
      ev = 1 / ( pars[['psill']] * fac[['values']][idx] )
      v_add = ev * pkern_toep_mult(cy, fac[['vectors']][, idx], cx, idx_obs, gdim)^2

    } else {

      # column indices in component correlation matrices corresponding to eigen-value idx
      idx_yx = pkern_vec2mat(idx, sg[['gdim']], out='list')

      # kronecker product of component columns to get large column vector
      add_yx = ( pars[['psill']] * kronecker(add_x2[, idx_yx[['j']]], add_y2[, idx_yx[['i']]]) )^2
      v_add = ev[idx] * add_yx
    }

    # add to running total in v
    v_rem[] = v_rem[] + as.numeric(v_add)
  }
  if(!quiet) close(pb)
  return(pars[['psill']] + pars[['eps']] + as.numeric(v_gls) - v_rem)
}











