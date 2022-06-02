# pkern_vario.R
# Dean Koch, 2022
# sample semi-variograms

#' Theoretical variogram function
#'
#' Computes the value of the semi-variogram function defined by covariance model `pars`
#' at the component (y and x) lags in list `d`.
#'
#' If `d` is a list, its 'y' and 'x' components should be equal length nonnegative numeric
#' vectors. The function returns the semi-variance in a vector of the same length.
#'
#' If `d` is a numeric vector, it is interpreted as a set of distances to at which to
#' evaluate the range of the (anisotropic) semi-variogram function. The function returns
#' a list with named numeric vectors 'min' and 'max'.
#'
#' @param pars list of the form returned by `pkern_pars` with entries 'y', 'x', 'eps', 'psill'
#' @param d numeric vector or list with vector entries 'y' and 'x', the distances to evaluate
#'
#' @return vector (for list `d`) or data.frame (for vector `d`)
#' @export
#'
#' @examples
#' gdim = c(10, 15)
#' d_max = sqrt(sum(gdim^2))
#' pars = pkern_pars(gdim, 'mat')
#'
#' d_test = seq(0, d_max, length.out=1e2)
#' vario = pkern_vario_fun(pars, d=list(y=0, x=d_test))
#' plot(d_test, vario)
#' vario_lims = pkern_vario_fun(pars, d=d_test)
#' lines(d_test, vario_lims[,1])
#' lines(d_test, vario_lims[,2])
#'
pkern_vario_fun = function(pars, d=NULL)
{
  # 1d vector case
  if( !is.list(d) )
  {
    # component distances to test along each axis
    d0 = rep(0, length(d))

    # compute three sets of component distances: varying x, y, and both equally
    cov_y = pkern_vario_fun( pars, list(y=d, x=d0) )
    cov_x = pkern_vario_fun( pars, list(y=d0, x=d) )
    cov_yx = pkern_vario_fun( pars, list(y=sqrt(2)*d, x=sqrt(2)*d))

    # compute range of theoretical semi-variogram at each test distance
    cov_min = pmin(cov_y, cov_x, cov_yx)
    cov_max = pmax(cov_y, cov_x, cov_yx)
    return( data.frame(min=cov_min, max=cov_max) )
  }

  # take product of correlation functions
  nm_yx = c('y', 'x')
  corrvals = do.call('*', Map(\(p, dyx) pkern_corr(p, dyx), p=pars[nm_yx], dyx=d[nm_yx]))

  # return the variogram
  return( pars[['eps']] + pars[['psill']] * ( 1 - corrvals ) )
}


#' Sub-grid point sampler for grid data
#'
#' Sample `n` locations from the non-NA points in the input grid `g`, optionally using
#' them as centers to place `n` sub-grids of the specified size and resolution.
#'
#' The function draws a sample of `n` locations (uniformly at random) from the non-NA
#' points in the input grid `g`, and returns their vector index. If there are fewer
#' than `n` locations available, they are all returned.
#'
#' When `lag_max > 1`, the function also identifies a regular sub-grid of the Moore
#' neighbourhood of (integer) radius `lag_max` around each of the sample points, with
#' `interval-1` grid lines separating each sub-grid line. `interval` is the factor by
#' which the resolution of the sub-grid is scaled to get the original grid resolution;
#' It must evenly divide `lag_max`.
#'
#' For a given `interval`, the grid `g` can be partitioned into `interval^2` distinct
#' non-overlapping sub-grids. When `over=FALSE` (the default), the function apportions
#' its `n` point samples as evenly as possible among these disjoint subsets. This ensures
#' that if `n` is less than or equal to `interval^2`, and there are no NAs, there can be
#' no repetition (or "over"lap) of points in the returned sub-grids.
#'
#' @param g any grid object accepted or returned by `pkern_grid`
#' @param n integer > 0, the maximum number of center points to sample
#' @param lag_max integer, Moore neighborhood radius (ie the maximum queen's distance)
#' @param interval integer > 0, the down-scaling factor for sub-grids of `g`
#' @param over logical, indicates to allow overlapping sub-grids (when they can be avoided)
#'
#' @return If `lag_max=0` (the default), the function returns the sample indices as a
#' length-`n` integer vector. If `lag_max` is positive, the function returns a list of
#' `n+1` indexing vectors; the first locates the `n` sub-grid center points, and the
#' following `n` vectors locate the points in each sub-grid (including the center points
#' itself)
#'
#' @export
#'
#' @examples
#' # define a grid
#' gdim = c(100, 100)
#' ng = prod(gdim)
#'
#' # get an ordinary random sample with default settings
#' idx_sample = pkern_sample_pt(g=gdim)
#' pkern_plot(seq(ng) %in% idx_sample, gdim)
#'
#' # reduce or increase number of center points from default 100
#' idx_sample = pkern_sample_pt(gdim, n=10)
#' pkern_plot(seq(ng) %in% idx_sample, gdim)
#'
#' # sampled from Moore neighbourhoods of radius 6
#' n = 10
#' idx_sample = pkern_sample_pt(gdim, n=n, lag_max=6L)
#' pkern_plot(seq(ng) %in% unlist(idx_sample), gdim, col_grid='white')
#'
#' # plot each list element a different color
#' group_sample = rep(0L, prod(gdim))
#' for(i in seq(1 + n)[-1]) group_sample[ idx_sample[[i]] ] = i-1L
#' pkern_plot(group_sample, gdim, breaks=c('not sampled', seq(n)), zlab='sub-grid')
#'
#' # When interval > 1, the function attempts to avoid overlap whenever possible
#' interval = 2
#' n = interval^2 # to get disjoint results n must be less than or equal to interval^2
#' lag_max = 10 * interval # vary to get larger/smaller subsets. max allowable: min(gdim)/2
#' idx_sample = pkern_sample_pt(gdim, n=n, interval=interval, lag_max=lag_max)
#' idx_overlap = rowSums( sapply(idx_sample[-1], \(i) seq(ng) %in% i) )
#' pkern_plot(as.integer(idx_overlap), gdim, zlab='times sampled')
#'
#' # plot each list element a different color
#' group_sample = rep(0L, prod(gdim))
#' for(i in seq(1 + interval^2)[-1]) group_sample[ idx_sample[[i]] ] = i-1L
#' pkern_plot(group_sample, gdim, breaks=c('not sampled', seq(interval^2)), zlab='sub-grid')
#'
#' # compare with over=TRUE (usually results in overlap - try running a few times)
#' idx_sample_compare = pkern_sample_pt(gdim, n=n, interval=interval, lag_max=lag_max, over=TRUE)
#' idx_overlap_compare = rowSums( sapply(idx_sample_compare[-1], \(i) seq(ng) %in% i) )
#' pkern_plot(as.integer(idx_overlap_compare), gdim, zlab='times sampled')
#'
#' # only non-NA points are eligible in initial sample of center points
#' g = pkern_grid(gdim)
#' g$gval = rep(NA, ng)
#' idx_obs = sample.int(ng, ng/1e2)
#' g$gval[idx_obs] = 'non-NA'
#' pkern_plot(g)
#'
#' # draw a sample of center points and indicate sub-grids in color
#' idx_sample = pkern_sample_pt(g, n=10, lag_max=6L, interval=2)
#' g$gval[unlist(idx_sample[-1])] = 'sub-grid'
#' g$gval[idx_sample[[1]]] = 'center point'
#' pkern_plot(g)
#'
pkern_sample_pt = function(g, n=1e2, lag_max=0, interval=1L, over=FALSE)
{
  # unpack the grid object
  g = pkern_grid(g)
  gdim = g[['gdim']]
  z = g[['gval']]

  # identify non-missing points in grid
  if( is.null(z) ) z = rep(0L, prod(gdim))
  is_obs = !is.na(z)

  # all-missing case treated as all-observed
  if( !any(is_obs) ) is_obs = rep(TRUE, length(is_obs))

  # check for valid parameters in request
  gdim_inner = gdim - 2 * lag_max
  if( !( all(gdim_inner > 0 ) ) ) stop('lag_max cannot exceed min(gdim)/2')

  # find subset for which Moore neighbourhood is completely inside grid
  ij_inner = lapply(setNames(gdim, c('i', 'j')), \(d) seq(1+lag_max, d-lag_max) )
  is_inner = pkern_sub_idx(gdim, ij=ij_inner)

  # compute their vector index and grid positions
  is_eligible = is_obs & is_inner
  idx_eligible = which(is_eligible)
  ij_obs = pkern_vec2mat(idx_eligible, gdim)

  # sample center points
  n = min(sum(is_eligible), n)
  idx_center_all = sample(idx_eligible, n)
  if(lag_max == 0) { idx_list = list(idx_center_all) } else {

    # pick a different sample that avoids overlap
    if(!over)
    {
      # the number of disjoint sub-grids and the index of points
      n_disjoint = interval^2
      idx_inner = which(is_inner)

      # identify all disjoint sub-grid origins within inner sub-grid,
      ij_origin = apply(expand.grid(seq(interval)-1, seq(interval)-1), 1, identity, simplify=FALSE)
      ij_all = lapply(ij_origin, \(ij) Map(\(d, r) seq(1L+r, d, by=interval), d=gdim_inner, r=ij))
      idx_disjoint = lapply(ij_all, \(ij) pkern_sub_idx(gdim_inner, unname(ij), idx=TRUE))

      # omit NA points
      idx_disjoint = lapply(idx_disjoint, \(ij) ij[ is_obs[idx_inner][ij] ] )

      # the number of non-NA points in each disjoint sub-grid
      n_obs_disjoint = sapply(idx_disjoint, length)

      # apportion samples evenly to remaining sample sites, building n_per in a loop
      n_per = integer(n_disjoint)
      n_remain = n
      while( (n_remain > 0) & any(n_obs_disjoint > 0) )
      {
        # find minimum number of non-NA points remaining in nonempty sub-grids
        is_nonempty = n_obs_disjoint > 0
        n_nonempty = sum(is_nonempty)

        # draw the same number of points from each one
        n_each = min(min(n_obs_disjoint[is_nonempty]), floor(n_remain/n_nonempty))

        # exit case: fewer points needed than available sub-grids
        if(n_each == 0)
        {
          # draw a single point from a subset of the sub-grids
          n_each = 1
          n_rem = sample(which(is_nonempty), n_remain)
          is_nonempty = seq(n_disjoint) %in% n_rem
        }

        # draw at most this many points from each non-empty sub-grid
        n_per[is_nonempty] = n_per[is_nonempty] + n_each
        n_obs_disjoint = n_obs_disjoint - n_each
        n_remain = n_remain - n_each * sum(is_nonempty)
      }

      # sample n_per[i] sub-grid origins from ith disjoint set
      idx_center_inner = unlist( Map(\(n, idx) sample(idx, n), n=n_per, idx=idx_disjoint) )

      # remap to original (not inner) grid
      idx_center_all = idx_inner[idx_center_inner]
    }

    # check for invalid down-scaling factor then make center box template
    if( ( lag_max %% interval ) != 0) stop('interval must divide lag_max')
    lag_max = min(c(lag_max, gdim))
    box_offset = seq(-lag_max, lag_max, by=as.integer(interval))

    # loop over center points, compute index of all points in Moore neighbourhood
    idx_list = lapply(idx_center_all, \(idx) {

      # find grid (i,j) indices for box template centered at center point
      ij_center = lapply(pkern_vec2mat(idx, gdim, out='list'), \(idx) idx + box_offset)

      # find vectorized index corresponding to this sample
      pkern_sub_idx(gdim, ij_center, idx=TRUE)
    })

    # set first entry to be the vector of all center points
    idx_list = c(list(idx_center_all), idx_list)
  }

  # filter length-1 entries and collapse length-1 results before returning
  idx_list = idx_list[sapply(idx_list, length) > 1]
  if( length(idx_list) == 1 ) idx_list = unlist(idx_list)
  return(idx_list)
}


#' Sample point pair absolute differences for use in semi-variogram estimation
#'
#' Computes the absolute differences for point pairs in `g`, along with their separation
#' distances. If no sample point index is supplied (in `idx`), the function samples points
#' at random using `pkern_sample_pt`.
#'
#' If `idx` has length n, there are n_pp(n) = (n^2 - n) / 2 possible point pairs. The inverse
#' n(n_pp), is used to determine the maximum number of sample points in `g` to use in order to
#' satisfy the argument to `n_pp` (a random sub-sample of `idx` is taken as needed).
#'
#' @param g any grid object accepted or returned by `pkern_grid`, containing non-NA data
#' @param idx optional integer vector indexing the points to sample
#' @param n_max integer maximum number of point pairs to sample
#' @param n_bin integer number of distance bins to assign (passed to `pkern_add_bins`)
#'
#' @return Results are returned in a data frame with each row representing one point pair.
#' Fields include 'dabs' and 'd', the absolute difference and distance mentioned earlier,
#' along with a number of indexing vectors for both point locations and relative separation.
#' 'bin' is an integer splitting distances into `n_bin` categories.
#' @export
#'
#' @examples
#' # make example grid and reference covariance model
#' gdim = c(22, 15)
#' n = prod(gdim)
#' g_obs = pkern_grid(gdim)
#' pars = pkern_pars(g_obs, 'mat')
#'
#' # generate sample data and sample semi-variogram
#' g_obs$gval = pkern_sim(g=g_obs, pars)
#' vg = pkern_sample_vg(g_obs)
#' str(vg)
#'
#' # pass to plotter and overlay the model that generated the data
#' pkern_plot_semi(vg, pars)
#'
#' # repeat with smaller sample size
#' pkern_plot_semi(pkern_sample_vg(g_obs, 1e2), pars)
#' pkern_plot_semi(pkern_sample_vg(g_obs, 1e3), pars)
#'
#' # use a set of specific points
#' n_sp = 10
#' ( n_sp^2 - n_sp ) / 2 # the number of point pairs
#' vg = pkern_sample_vg(g_obs, idx=sample.int(n, n_sp))
#' pkern_plot_semi(vg, pars)
#'
#' # repeat with all point pairs sampled (not recommended for big data sets)
#' vg = pkern_sample_vg(g_obs, n_pp=Inf)
#' pkern_plot_semi(vg, pars)
#' ( n^2 - n ) / 2 # the number of point pairs
#'
pkern_sample_vg = function(g, n_pp=1e4, idx=NULL, n_bin=25)
{
  # unpack expected inputs
  g = pkern_grid(g)
  gdim = g[['gdim']]
  gres = g[['gres']]
  z = g[['gval']]
  if( is.null(z) ) stop('g must have element "gval" (the data vector)')

  # the number of sample points required to get n point pairs, the inverse of f(n)=(n^2-n)/2
  n_obs = min(floor( (1 + sqrt(1 + 8*n_pp)) / 2 ), sum(!is.na(z)))

  # call point sampler if a subset of points was not specified
  if( is.null(idx) ) idx = pkern_sample_pt(g, n_obs)

  # verify that input sample point locations are non-NA and sub-sample as needed
  idx = idx[ !is.na(z[idx]) ]
  if( length(idx) > n_obs ) idx = sample(idx, n_obs)
  n_obs = length(idx)

  # compute grid indices (i, j) for the sample points
  ij_obs = pkern_vec2mat(idx, gdim)

  # anonymous function to compute lower triangular part of a 1d distance matrix
  vec_dist = \(x) c(stats::dist(x, method='manhattan', diag=T))

  # compute dimension-wise grid line (i,j) distances between all point pairs
  dmat = stats::setNames(apply(ij_obs, 2, vec_dist, simplify=F), c('di', 'dj'))

  # indexes the distance matrix entries vectorized by the `c` in vec_dist
  idx_lower = lower.tri(matrix(NA, n_obs, n_obs))

  # i,j indices corresponding to vectorized entries of idx_lower
  idx_i_lower = matrix(rep(seq(n_obs), n_obs), n_obs)[idx_lower]
  idx_j_lower = matrix(rep(seq(n_obs), each=n_obs), n_obs)[idx_lower]

  # i,j indices for each point pair, and their vectorized indices
  ij_pp = list(p1=ij_obs[idx_i_lower,], p2=ij_obs[idx_j_lower,])
  idx_pp = lapply(ij_pp, \(ij) pkern_mat2vec(ij, gdim))

  # compute absolute differences in data for each point pair, then separation distance
  dabs_all = abs( apply( sapply(idx_pp, \(i) z[i]), 1, diff ) )

  # compile everything in a data frame then append distance info
  vg = data.frame(dabs=dabs_all, idx_pp, ij_pp, dmat)
  vg[c('dy', 'dx')] = rep(gres, each=nrow(vg)) * vg[c('di', 'dj')]
  vg[['d']] = sqrt( rowSums( vg[c('dy', 'dx')]^2 ) )

  # split into distance bins before returning
  return( pkern_add_bins(vg, n_bin) )
}


#' Add bin labels to a variogram data frame
#'
#' Helper function for grouping the rows of input data frame `vg` into `n_bins` bins
#' according to the value of the (numeric) distance column `d`. This uses either `base::cut`
#' or, if `probs` is supplied, `stats::quantile`.
#'
#' By default, the function sets `probs` to a sequence of length `1+n_bin` evenly splitting
#' the interval [0,1] to ensure approximately equal sample sizes for each bin. Setting
#' `probs=NA` instead sets the bin endpoints such that the range of distances is split
#' evenly (note this may produce empty bins)
#'
#' The function is called by `pkern_sample_vg` and `pkern_plot_semi` (when column `bin` is
#' missing). It can also be used to recompute bins after an `rbind` of multiple variogram
#' data frames.
#'
#' @param vg data frame with numeric column 'd'
#' @param n_bin integer number of distance bins to assign
#' @param probs numeric vector of quantile probabilities to establish breakpoints (length `n_bin+1`)
#'
#' @return same as input `vg` but with integer column `bin` added/modified
#' @export
#'
#' @examples
#' distance_df = data.frame(d=runif(25))
#' pkern_add_bins(distance_df)
#'
#' # specify fewer bins and set up quantiles explicitly
#' pkern_add_bins(distance_df, n_bin = 5) # same as ...
#' pkern_add_bins(distance_df, n_bin = 5, probs=seq(0, 1, length.out=6))
#'
#' # break range of distances into evenly spaced bins (of varying sample sizes)
#' pkern_add_bins(distance_df, n_bin = 5, probs=NULL)
#'
pkern_add_bins = function(vg, n_bin=25, probs=NULL)
{
  # check for distance column
  d = vg[['d']]
  if( is.null(d) ) stop('missing distance column "d" in variogram data frame')

  # check for valid input
  if(n_bin > nrow(vg))
  {
    n_bin = nrow(vg)
    warning(paste('n_bin reduced to', n_bin))
  }

  # assign default quantile breaks
  if( is.null(NULL) ) probs = seq(0, 1, length.out=1L+n_bin)

  # n_bin=1 puts everything in bin 1
  if(n_bin == 1) return( modifyList(vg, list(bin=1L)) )

  # split range evenly into bins when `probs=NULL`
  if( anyNA(probs) ) return( modifyList(vg, list( bin=as.integer(cut(d, breaks=n_bin))) ) )

  # otherwise find bins of roughly equal content
  return( modifyList(vg, list(bin=findInterval(d, quantile(d, probs=probs)))) )
}
