# pkern_index.R
# Dean Koch, 2022
# Convenience functions for indexing grid points


#' Column-vectorization indices
#'
#' Maps matrix indices to vectorization indices. The function returns the vector
#' indices associated with the supplied matrix (i, j) indices.
#'
#' Column vectorization (as in base::as.vector) builds a length(mn) vector by stacking
#' the columns of an m X n matrix, with the leftmost column appearing first in the vector
#' and the rightmost column last. Matrix element i,j gets mapped to element `(i + m*(j-1))`
#' in the vector. This function returns that index.
#'
#' `ij` can be a matrix or a list of length-n vectors "i" and "j" (in that order), or a vector
#' representing a single point at the given row and column number. `ni` can be a vector of the
#' form `c(ni, nj)` (the return value of `dim` for example) in which case its first element
#' (the number of rows) is used.
#'
#' @param ij n x 2 matrix, the row and column indices
#' @param ni number of rows in the matrix
#' @param simplified, if FALSE, the function returns an n x 1 matrix
#'
#' @return a vector of indices in the vectorized system
#' @export
#'
#' @examples
#'
#' # get column-vectorized ordering of points from matrix indices
#' ni = 10
#' gyx = expand.grid(i=seq(ni), j=seq(ni))
#' pkern_mat2vec(gyx, ni)
pkern_mat2vec = function(ij, ni, simplified=TRUE)
{
  # handle vector input to ni
  if( length(ni) > 1 ) ni = ni[1]

  # coerce list to matrix
  if( is.list(ij) )
  {
    if( !all( diff(sapply(ij, length)) == 0 ) ) stop('elements of list ij must have equal length')
    ij = do.call(cbind, ij)
  }

  # handle vector input (single point)
  if( is.vector(ij) ) ij = matrix(ij, 1)

  # coerce input to matrix
  ij = as.matrix(ij)

  # check for invalid input
  if( any(ij[,1] > ni) ) stop('ij contains "i" indices exceeding ni')

  # return the vectorized index
  idx = ij %*% c(1, ni) - ni
  if( !simplified ) return( idx )
  return( as.vector(idx) )
}


#' Inverse column-vectorization indices
#'
#' This function is for doing the inverse of `as.vector(M)` where `M` is a matrix. It returns
#' the row and column numbers (i,j) associated with a given vector element after
#' column-vectorization.
#'
#' If a grid size vector `c(ni, nj)` is passed to `ni`, the function uses its first element.
#'
#' @param idx a vector of positive integers
#' @param ni number of rows in the matrix
#' @param out either 'matrix' or 'list'
#'
#' @return a two column matrix of integers (row and column numbers) with `length(idx)` rows
#' @export
#'
#' @examples
#'
#' # show how elements are ordered in `base::matrix`
#' ni = 5
#' nj = 6
#' matrix.indexing = matrix(1:prod(ni, nj), ni)
#' print(matrix.indexing)
#' as.vector(matrix.indexing)
#'
#' # doing the inverse
#' pkern_vec2mat(2, ni)
#' pkern_vec2mat(c(1,2,7), ni)
pkern_vec2mat = function(idx, ni, out='matrix')
{
  # handle vector input to ni
  if( length(ni) > 1 ) ni = ni[1]

  # compute column and row numbers
  cnum = ceiling( idx / ni ) |> as.integer()
  rnum = ( idx - ( ni * (cnum - 1) ) ) |> as.integer()

  # return as matrix
  if(out == 'matrix') return( cbind(i=rnum, j=cnum) )
  if(out == 'list') return( list(i=rnum, j=cnum) )
}


#' Find column-vectorized index of a sub-grid
#'
#' Returns the column-vectorized index of a sub-grid with respect to the column-major
#' vectorized order of the full grid of size `N = prod(gdim)`, based on the supplied
#' grid line numbers `ij`. The returned length-`N` logical vector identifies points
#' in the full grid lying at the intersections of grid lines in `ij`, a sub-grid with
#' dimensions `gdim_sg = sapply(ij, length)`.
#'
#' If `idx=TRUE`, the function returns an integer vector of length `prod(gdim_sg)`,
#' the positions of the sub-grid points in the vectorization. This is the same as piping
#' the default logical output to `which`; except if `nosort=FALSE`, in which case
#' the output order depends on the ordering in `ij`: Letting `j = c(j1, j2, ..., in)` and
#' `i = c(i1, i2, ...im)` the function returns:
#'
#' (j1, i1), (j1, i2), ... (j1, im), (j2, i1), (j2, i2), ..., (j3, i1), ... (jm, in)
#'
#' If `ij$i` is missing, the function returns results for all rows, and if `ij$j` it
#' returns results for all columns.
#'
#' @param gdim c(ni, nj), the number rows and columns in the full grid
#' @param ij list containing vectors "i" and "j", the sub-grid rows and columns
#' @param nosort logical, skips sorting the input vectors in `ij`
#' @param idx logical, indicates to return indices (default TRUE) versus logical vector
#'
#' @return integer vector (length `n_sg`) or logical vector (length `n`)
#' @export
#'
#' @examples
#'
#' # example grid and a particular grid point
#' gdim = c(i=10, j=13)
#' ij_list = list(i=6, j=3)
#'
#' # pkern_sub_idx returns a logical vector indexing the point (or the index itself)
#' is_pt = pkern_sub_idx(gdim, ij_list)
#' pkern_plot(is_pt, gdim, col_grid='white')
#' pkern_sub_idx(gdim, ij_list, idx=TRUE)
#' pkern_mat2vec(ij_list, gdim) # equivalent when ij_list is a single point
#'
#' # if i or j are omitted from ij, the function returns the full row or column
#' is_col2 = pkern_sub_idx(gdim, ij_list['i'])
#' is_row3 = pkern_sub_idx(gdim, ij_list['j'])
#' pkern_plot(is_col2, gdim, col_grid='white', breaks=c('other', paste('row', ij_list['i'])))
#' pkern_plot(is_row3, gdim, col_grid='white', breaks=c('other', paste('column', ij_list['j'])))
#' pkern_plot(is_col2 + is_row3, gdim, col_grid='white')
#'
#' # indices in column-vectorized order
#' pkern_sub_idx(gdim, list(i=2), idx=TRUE)
#' pkern_sub_idx(gdim, list(j=3), idx=TRUE)
#' pkern_sub_idx(gdim, idx=TRUE) # call without arguments returns all indices
#'
#' # bigger sub-grid example
#' origin_sg = c(5, 2) # assign i,j of top left point
#' gdim_sg = c(3, 4) # sub-grid dimensions (make sure this fits in gdim!)
#' ij_list = stats::setNames(Map(\(d, o) o + seq(d) - 1, gdim_sg, origin_sg), c('i', 'j'))
#' is_sg = pkern_sub_idx(gdim, ij=ij_list)
#' pkern_plot(is_sg, gdim, col_grid='white')
#'
#' # plot the index values: column major vectorization with y descending, x ascending
#' idx_sg = pkern_sub_idx(gdim, ij=ij_list, idx=TRUE)
#' vec_order = rep(NA, prod(gdim))
#' vec_order[is_sg] = as.character(idx_sg)
#' pkern_plot(vec_order, gdim, col_grid='black', zlab='vector index')
#'
#' # example with j indices supplied in reverse (descending) order
#' ij_list_xflip = modifyList(ij_list, list(j=rev(ij_list[['j']])))
#'
#' # ordering in the vectors ij$i and ij$j doesn't matter if `nosort=FALSE` or `idx=FALSE`
#' identical(is_sg, pkern_sub_idx(gdim, ij=ij_list, nosort=TRUE))
#' all.equal(which(is_sg), pkern_sub_idx(gdim, ij=ij_list_xflip, idx=TRUE))
#'
#' # when `nosort=TRUE` and `idx=TRUE` we get the same indices but in a different order
#' idx_sg_xflip = pkern_sub_idx(gdim, ij=ij_list_xflip, idx=TRUE, nosort=TRUE)
#' all.equal(sort(idx_sg), sort(idx_sg_xflip))
#' all.equal(idx_sg, idx_sg_xflip)
#' vec_order[is_sg] = as.character(idx_sg_xflip)
#' pkern_plot(vec_order, gdim, col_grid='black', zlab='vector index')
#'
pkern_sub_idx = function(gdim, ij=NULL, idx=FALSE, nosort=FALSE)
{
  # check input and set expected order
  ij_nm = c(y='i', x='j')
  if( !all( c('y', 'x') %in% names(gdim) ) ) gdim = stats::setNames(gdim, c('y', 'x'))
  if( is.null(ij) ) ij = stats::setNames(lapply(gdim, seq), ij_nm)
  if( is.null(names(ij)) & (length(ij) == 2) ) names(ij) = ij_nm
  if( sum( ij_nm %in% names(ij) ) != length(ij) ) stop('unexpected names in ij')

  # set default i to select all rows
  if( is.null(ij[['i']]) ) { ij[['i']] = seq( gdim['y'] ) } else {

    # otherwise coerce to integer and sort as needed
    ij[['i']] = as.integer(ij[['i']])
    if( !nosort ) ij[['i']] = sort(ij[['i']])
  }

  # set default j to select all columns
  if( is.null(ij[['j']]) ) { ij[['j']] = seq( gdim['x'] ) } else {

    # otherwise coerce to integer and sort as needed
    ij[['j']] = as.integer(ij[['j']])
    if( !nosort ) ij[['j']] = sort(ij[['j']])
  }

  # count desired sub-grid dimensions and compute the indexing vector
  n_ij = sapply(ij, length)
  idx_out = rep(ij[['i']], n_ij['j']) + rep(gdim['y'] * (ij[['j']] - 1L), each=n_ij['i'])
  if(idx) return(as.integer(unname(idx_out)))

  # compute the logical vector
  is_out = logical(prod(gdim))
  is_out[idx_out] = TRUE
  return(is_out)

}


#' Up-scale or down-scale a grid
#'
#' Changes the resolution of a grid by a factor of (positive integer) `up` or `down`.
#'
#' This constructs a new grid of based on `g`, where: if `up` is supplied, a sparser
#' (up-scaled) version of `g` is returned, with every `up`th grid line sampled; and
#' if `down` is supplied, a denser version of `g` is returned where every `down`th
#' grid line is copied over from `g` (and the rest initialized to NA).
#'
#' This effects the scaling of resolution (`g$gres`) by `up` (or `1/down`). It does
#' not do any interpolation to fill new grid-points created with `up` and it does not
#' change the bounding box of `g`.
#'
#' @param g any object accepted or returned by `pkern_grid`
#' @param up integer > 0, an up-scaling factor
#' @param down integer > 0, a down-scaling factor
#' @param vals logical, indicates to return data vector as well as grid info
#'
#' @return a grid list of the form returned by `pkern_grid`
#' @export
#'
#' @examples
#'
#' # example data
#' gdim = c(50, 53)
#' g = pkern_grid(gdim)
#' gval = pkern_sim(g, modifyList(pkern_pars(g), list(eps=0)), quiet=TRUE)
#' g_obs = modifyList(g, list(gval=gval))
#' pkern_plot(g_obs)
#'
#' # upscale
#' pkern_plot(pkern_rescale(g=g_obs, up=1)) # does nothing
#' pkern_plot(pkern_rescale(g=g_obs, up=2))
#' pkern_plot(pkern_rescale(g=g_obs, up=3))
#'
#' # downscale
#' pkern_plot(pkern_rescale(g=g_obs, down=1)) # does nothing
#' pkern_plot(pkern_rescale(g=g_obs, down=2))
#' pkern_plot(pkern_rescale(g=g_obs, down=3))
#'
#' # turn vals off to get grid configuration only
#' pkern_plot(g)
#' pkern_plot(pkern_rescale(g=g_obs, up=2, vals=FALSE))
#' pkern_plot(pkern_rescale(g=g_obs, down=2, vals=FALSE))
#'
pkern_rescale = function(g, up=NULL, down=NULL, vals=TRUE)
{
  # user has to pick one or the other
  is_up = !is.null(up)
  is_down = !is.null(down)
  if(is_up & is_down) stop('both up and down were specified')
  if( !(is_up | is_down) ) stop('either up or down must be specified')

  # unpack the grid object as list
  g = pkern_grid(g, vals)
  vals = vals & !is.null(g[['gval']])
  gdim = g[['gdim']]

  # up-scaling
  if(is_up)
  {
    # check for invalid up-scaling arguments
    msg_max = paste(gdim, collapse=', ')
    up = as.integer(up)
    if( any(up < 1) ) stop('upscaling factors cannot be less than 1')
    if( !any(up < gdim) ) stop( paste('upscaling factors must be less than', msg_max) )

    # set up dimensions of sub-grid
    ij_values = stats::setNames(Map(function(d, r) seq(1, d, r), d=gdim, r=up), c('i', 'j'))
    gdim_new = stats::setNames(sapply(ij_values, length), c('y', 'x'))

    # overwrite data vector with the subset then update other fields in g
    idx_sub = pkern_sub_idx(gdim, ij_values, idx=TRUE)
    if(vals) g[['gval']] = g[['gval']][idx_sub]
    g[['gres']] = g[['gres']] * up
    g[['gdim']] = gdim_new
    g[['gyx']] = Map(function(yx, idx) yx[idx], yx=g[['gyx']], idx=ij_values)

    return(g)
  }

  # check for invalid down-scaling arguments
  down = as.integer(down)
  if( any(down < 1) ) stop('downscaling factors cannot be less than 1')

  # set up dimensions of super-grid
  gdim_new = gdim + (down - 1L) * (gdim - 1L)
  gres_new = g[['gres']] / down
  g_new = pkern_grid(list(gdim=gdim_new, gres=gres_new), vals=FALSE)
  g_new[['gyx']] = Map(function(old, new) new + old[1], old=g[['gyx']], new=g_new[['gyx']])

  # copy data to the sub-grid of the output and return
  if(!vals) return(g_new)
  ij_sub = stats::setNames(Map(function(d, b) seq(1, d, b), d=gdim_new, b=down), c('i', 'j'))
  idx_sub = pkern_sub_idx(gdim_new, ij=ij_sub, idx=TRUE)
  if( !is.null(g[['gval']]) ) g_new[['gval']][idx_sub] = g[['gval']]
  return(g_new)
}


#' Check vectorized grid data for non-NA points that form a complete sub-grid
#'
#' If a gridded data-set `g_obs` has missing values (NAs), but the set of non-NA points
#' form a complete sub-grid, this function finds its grid lines, resolution scaling factor,
#' and dimensions. If no eligible sub-grids are found, the function returns NULL.
#'
#' A sub-grid is only eligible if it contains all of the non-NA points in `g_obs` and none
#' of the NAs; eg if a single point missing from the sub-grid, or a single non-NA point lies
#' outside the sub-grid, the function will fail to detect the sub-grid and return NULL. If no
#' points are NA, the function returns indices for the full grid.
#'
#' In the special case that `g_obs` is a logical vector, it is interpreted as indicating
#' the non-NA locations in a grid (and isn't itself checked for NAs). the argument `gdim`
#' must be supplied in this case. When grid dimensions can otherwise be derived from
#' `g_obs`, the function does so and overrides any argument to `gdim`.
#'
#' @param g_obs logical vector, or any other object accepted by `pkern_grid`
#' @param gdim integer vector, the grid dimensions (ny, nx)
#' @param g_out logical, indicates to return a grid list object
#'
#' @return NULL or list of information about the location and spacing of the sub-grid
#' within `g` (see details)
#' @export
#'
#' @examples
#'
#' # define a grid and example data
#' gdim = c(50, 53)
#' g = pkern_grid(gdim)
#' gval = pkern_sim(g, modifyList(pkern_pars(g), list(eps=0)), quiet=TRUE)
#' g_obs = modifyList(g, list(gval=gval))
#' pkern_plot(g_obs)
#'
#' # define a supergrid containing the original data and make sure we can find it
#' g_obs_big = pkern_rescale(g_obs, down=3)
#' str(pkern_sub_find(g_obs_big))
#'
#' # define a smaller sub-grid at random
#' spacing = sapply(floor(gdim/10), function(x) 1 + sample.int(x, 1))
#' gdim_sg = sapply(floor( (gdim - 1) / spacing), function(x) sample.int(x, 1))
#' ij_first = sapply(gdim - ( spacing * gdim_sg ), function(x) sample.int(x, 1))
#'
#' # find index of sub-grid lines and vectorized index of points
#' ij_sg = Map(\(idx, r, n) seq(idx, by=r, length.out=n), idx=ij_first, r=spacing, n=gdim_sg)
#' is_sg = pkern_sub_idx(gdim, ij_sg, idx=FALSE)
#'
#' # sub grids with side length 1 have no spacing defined along that dimension
#' spacing[gdim_sg==1] = NA
#'
#' # assign values to the sub-grid points
#' g = pkern_grid(gdim)
#' g$gval[is_sg] = TRUE
#' pkern_plot(g, zlab='sub-grid')
#'
#' # call the function on g and check for expected results
#' subgrid_result = pkern_sub_find(g)
#' all.equal(unname(subgrid_result$gdim), gdim_sg)
#' all.equal(unname(subgrid_result$ij), ij_sg)
#' all.equal(unname(subgrid_result$res_scale), spacing)
#'
#' # or call on the vector and supply gdim separately
#' identical(subgrid_result, pkern_sub_find(g$gval, g_obs$gdim))
#' identical(subgrid_result, pkern_sub_find(!is.na(g$gval), g$gdim))
#'
pkern_sub_find = function(g_obs, gdim=NULL)
{
  # handle vector input
  if( is.logical(g_obs) & is.vector(g_obs) )
  {
    # logical vectors interpreted as indicating non-NAs
    if( anyNA(g_obs) ) stop('g_obs vector of logical class cannot have NAs')
    gdim = stats::setNames(as.integer(gdim), c('y', 'x'))

  } else {

    # make sure gdim doesn't get ignored for non-logical vector input
    if( is.vector(g_obs) & !is.list(g_obs) ) g_obs = list(gval=g_obs, gdim=gdim)

    # open as pkern list object
    g_result = pkern_grid(g_obs)

    # process only the first column of multi-layer input
    if( is.matrix(g_result[['gval']]) ) g_result[['gval']] = as.vector(g_obs[['gval']][,1])
    g_obs = !is.na( g_result[['gval']] )
    gdim = g_result[['gdim']]

    # handle sparse indexing
    if( !is.null(g_result[['idx_obs']]) )
    {
      # decompress and replace NAs with FALSE
      g_obs = g_obs[ g_result[['idx_obs']] ]
      g_obs[is.na(g_obs)] = FALSE
    }
  }

  # checks for valid arguments
  n_obs = sum(g_obs)
  if(n_obs < 2) return(NULL)
  if( is.null(gdim) ) stop('full grid dimensions gdim must be supplied when g_obs is a vector')
  msg_len = paste('Expected', prod(gdim), 'but got', length(g_obs))
  if( prod(gdim) != length(g_obs)) stop(paste('input g_obs was the wrong length.', msg_len))

  # need this to get indices of first, second, and last elements in sub-grid
  idx_obs = which(g_obs)

  # find the dimensions of the smallest sub-grid enclosing all observed points
  ij_bbox = pkern_vec2mat(c(idx_obs[1], idx_obs[n_obs]), gdim)
  gdim_bbox = apply(ij_bbox, 2, function(x) diff(x) + 1L)
  if( !all(gdim_bbox > 1) ) return(NULL)

  # compute number of rows in sub-grid and do first existence check
  skip_i = diff(idx_obs[1:2]) - 1L
  ni = as.integer( ifelse(skip_i > -1L, ( gdim_bbox[['i']] + skip_i ) / (1L + skip_i), NA) )
  if( ni == 1 ) skip_i = NA
  if( (ni == 0) | ( ni %% 1L != 0 ) ) return(NULL)

  # compute number of columns in sub-grid and do second existence check
  nj = n_obs / ni
  if( nj %% 1L != 0 ) return(NULL)
  skip_j = as.integer( ifelse(nj > 1, (gdim_bbox[['j']] - nj) / (nj - 1), NA) )

  # sub-grid resolution scaling factor and dimensions
  nm_dim = c('y', 'x')
  res_ij = setNames(1L + c(skip_i, skip_j), nm_dim)
  res_scale = 1L / res_ij
  gdim_sub = setNames(as.integer(c(ni, nj)), nm_dim)

  # sub-grid line indices
  ij_obs = Map(\(idx, r, n) seq(idx, by=r, length.out=n), idx=ij_bbox[1,], r=res_ij, n=gdim_sub)
  ij_na = sapply(ij_obs, anyNA)
  ij_obs[ij_na] = as.list(ij_bbox[1,])[ ij_na ]

  # final existence check
  idx_sub = pkern_sub_idx(gdim, ij_obs)
  if( !all( g_obs[idx_sub] ) ) return(NULL)

  # return sub-grid info in a list
  return( list(ij=setNames(ij_obs, nm_dim), res_scale=res_ij, gdim=gdim) )
}


#' Return coordinates of a grid of points in column-vectorized order
#'
#' Expands sets of y and x grid line locations in the column-major order i that
#' pkern uses for grid data. This ordering sorts the data into blocks of length equal
#' to the vertical dimension of the grid, where: y descends within blocks within blocks
#' and x increases among blocks.
#'
#' This is similar to `base::expand.grid` but with the first dimension descending instead
#' of ascending. Optionally, setting `out='sf'` builds an `sf` simple features object
#' from the coordinates and populates it with data (if any) from `g`. This can be slow for
#' large grids.
#'
#' @param g any other object accepted/returned by `pkern_grid`
#' @param out character indicating return value type, either 'list', 'matrix', or 'sf'
#' @param corners logical, indicates to only return the corner points
#'
#' @return a matrix, list, or sf POINT collection, in column vectorized order
#' @export
#'
#' @examples
#'
#' # 5x5 grid example
#' str( pkern_coords(5) )
#' pkern_coords(5, out='list')
#' pkern_coords(5, corner=TRUE)
#' pkern_coords(5, corner=TRUE, out='list')
#'
#' # 3x2 example and sf output type
#' pkern_coords(c(3,2))
#' if( requireNamespace('sf') ) {
#' pkern_coords(c(13, 12), out='sf')
#' }
#'
pkern_coords = function(g, out='matrix', corner=FALSE)
{
  # unpack input and take subset of corner points if requested
  g = pkern_grid(g)
  if( is.null(g[['crs']]) ) g[['crs']] = ''
  if(corner)
  {
    g[['gdim']] = c(y=2L, x=2L)
    g[['gyx']] = lapply(g[['gyx']], range)
  }

  # sort grid lines
  g[['gyx']][['y']] = sort(g[['gyx']][['y']], decreasing=TRUE)
  g[['gyx']][['x']] = sort(g[['gyx']][['x']], decreasing=FALSE)

  # compute the coordinates
  ij = pkern_vec2mat(seq(prod(g[['gdim']])), g[['gdim']], out='list')
  out_list = stats::setNames(Map(function(gl, idx) gl[idx], g[['gyx']], idx=ij), c('y', 'x'))

  # return as list
  if( out == 'list' ) return(out_list)

  # return as matrix
  out_mat = do.call(cbind, out_list)
  if( out == 'matrix' ) return( out_mat )

  # return as sf
  if( !startsWith(out, 'sf') ) stop('Argument `out` must be either "list", "matrix", or "sf"')
  sf_loaded = requireNamespace('sf', quietly=TRUE)
  if( !sf_loaded ) stop('sf package not loaded. Try library(sf)')
  cat(paste('processing', prod(g[['gdim']]), 'grid points...\n'))
  sf_out = sf::st_as_sf(as.data.frame(out_mat), coords=c('x', 'y'), crs=g[['crs']])
  return(sf_out)
}

