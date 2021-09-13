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

# conditional mean of points on a grid given zobs, a sample taken over a subgrid
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

  nsg = length(zobs)


  # check if precomputed matrices and indices were supplied
  if( !is.logical(precompute) )
  {
    # check for invalid input
    if( !is.list(precompute) ) stop('"precompute" must be a list')

    # unpack correlation matrices and eigendecomposition of marginal variance
    vx = precompute[['vx']]
    vy = precompute[['vy']]
    vx.cross = precompute[['vx.cross']]
    vy.cross = precompute[['vy.cross']]
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
    vx = pkern_corrmat(pars[['x']], length( gxy[[1]] ) )
    vy = pkern_corrmat(pars[['y']], length( gxy[[2]] ) )

    # build component cross-correlation matrices for points on vs off subgrid lines
    vx.cross = pkern_corrmat(pars[['x']], dims[1], ds=1/dsx, j=gxc, i=gx)
    vy.cross = pkern_corrmat(pars[['y']], dims[2], ds=1/dsy, j=gyc, i=gy)

    ## compute eigendecomposition(s)

    # slower method when there is missing data
    if( anyNA(zobs) )
    {
      # missing data case requires first evaluating the variance kronecker product
      idx.sg.obs = which( !is.na(zobs) )
      ed = eigen(kronecker(vx, vy)[idx.sg.obs, idx.sg.obs], symmetric=TRUE)

    } else {

      # with no missing data we can do eigendecompositions on component matrices
      idx.sg.obs = seq(nsg)
      ed = list(x=eigen(vx, symmetric=TRUE), y=eigen(vy, symmetric=TRUE))
    }
  }

  # return the list of precomputed objects if requested
  if( precompute ) return(list(vx=vx,
                               vy=vy,
                               ed=ed,
                               vx.cross=vx.cross,
                               vy.cross=vy.cross,
                               idx.on=idx.on,
                               idx.off=idx.off,
                               idx.xoff=idx.xoff,
                               idx.yoff=idx.yoff,
                               idx.sg.obs=idx.sg.obs))

  # missing data case
  if( length(idx.sg.obs) < nsg )
  {
    # omit the NA values from zobs
    zobs.crop = zobs[idx.sg.obs]

    # pointwise variances are equal to the eigenvalues plus the nugget variance
    pwv = pars[['nug']] + ed[['values']]

    # pad output vector with zeros (equivalent to subsetting cross-correlation matrices below)
    zobs.indep = rep(0, nsg)

    # transform `zobs` by product with inverse covariance to get mutual independence
    zobs.indep[idx.sg.obs] = ed[['vectors']] %*% ( ( t(ed[['vectors']]) %*% zobs.crop ) / pwv )

  } else {

    # kronecker product trick to get pointwise variances
    pwv = pars[['nug']] + kronecker(ed[['x']][['values']], ed[['y']][['values']])

    # kronecker product trick to get product with covariance inverse
    zobs.ortho = pkern_kprod(ed[['x']][['vectors']], ed[['y']][['vectors']], zobs, trans=TRUE) / pwv
    zobs.indep = pkern_kprod(ed[['x']][['vectors']], ed[['y']][['vectors']], zobs.ortho, trans=FALSE)
  }

  ## TESTING EIGENVALUE-BASED METHOD ***********************************************

  # evx = eigen(vx, symmetric=TRUE)
  # evy = eigen(vy, symmetric=TRUE)
  # z1 = pkern_kprod(evx[['vectors']], evy[['vectors']], zobs, trans=TRUE)
  # z2 = z1 / ( nug + kronecker(evx[['values']], evy[['values']]) )
  # ztrans = pkern_kprod(evx[['vectors']], evy[['vectors']], z2, trans=FALSE)

  ## ***********************************************


  ## NON-KRONECKER VERSION FOR MISSING DATA ****************************************

  # idx.miss = is.na(zobs)
  # zcomplete = zobs[-which(idx.miss)]
  # evxy = eigen(kronecker(vx, vy)[!idx.miss, !idx.miss], symmetric=TRUE)
  #
  # xx = evxy[['vectors']] %*% ( ( t(evxy[['vectors']]) %*% zcomplete ) / ( nug + evxy[['values']] ) )
  # ztrans = rep(0, length(zobs))
  # ztrans[!idx.miss] = xx

  # works okay!

  ## ***********************************************

  # initialize predictions vector and add the observed points
  zout = rep(as.numeric(NA), prod(dims))
  zout[idx.on] = zobs

  # compute predictions separately on the three subsets
  zout[idx.xoff] = pkern_kprod(vx, vy.cross, zobs.indep, trans=TRUE)
  zout[idx.yoff] = pkern_kprod(vx.cross, vy, zobs.indep, trans=TRUE)
  zout[idx.off] = pkern_kprod(vx.cross, vy.cross, zobs.indep, trans=TRUE)

  # finish
  return(zout)
}
