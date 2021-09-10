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
pkern_cmean = function(dims, xpars, ypars, zobs, gxy, pre=NULL, nug=1e-9)
{
  # Fast conditional mean computer for problems where sampled locations form a subgrid
  #
  # Uses method of Gilboa et. al, 2015 where eigendecomposition of the covariance
  # matrix of sampled locations is used with a nugget effect
  #
  # ARGUMENTS:
  #
  # dims: vector of two integers (nx, ny), the dimensions of the full grid
  # xpars: list of kernel parameters for the x-dimension component
  # ypars: list of kernel parameters for the y-dimension component
  # zobs: numeric vector of observed data
  # gxy: list of two integer vectors, the x and y grid lines of the subgrid
  #
  # DETAILS:
  #
  # zobs should be supplied in the usual column-vectorized order.
  #
  # The length of zobs must equal the product of the lengths of the vectors in gxy
  # (ie the number of grid line intersections), and the grid line numbers in gxy
  # should be ordered subsets of the sequences seq(dims[1]) and seq(dims[2])
  #
  # see sk_kern for details on xpars, ypars
  #
  # RETURN:
  #
  # the function returns a numeric vector of length prod(dims), containing the
  # observed data from zobs and the mean values at all other grid points, conditional
  # on zobs (in column-vectorized order).
  #

  #### unpack input

  # unpack grid dimensions and create vectors indexing all its grid lines
  nx = dims[1]
  ny = dims[2]
  ngx = length(gxy[[1]])
  ngy = length(gxy[[2]])

  # unpack and sort x and y grid line indices for observed data and find their complements
  gx = sort(gxy[[1]])
  gy = sort(gxy[[2]])
  gxc = seq(nx)[-gx]
  gyc = seq(ny)[-gy]

  # compute resolution increase factors
  dsx = diff(gx[1:2])
  dsy = diff(gy[1:2])

  if( is.null(pre) )
  {

    #### compute subset indices

    # index of observed grid points in full grid
    idx.obs = pkern_idx_sg(dims, i=gy, j=gx)

    # initialize predictions vector and add the observed points
    zout = rep(as.numeric(NA), prod(dims))
    zout[idx.obs] = zobs

    # indexing on full grid of *all* points lying on grid lines of the subgrid
    idx.x = pkern_idx_sg(dims, i=NULL, j=gx)
    idx.y = pkern_idx_sg(dims, i=gy, j=NULL)

    # index of unobserved points sharing an x grid line with an observed one
    idx.unobs.x = idx.x[!(idx.x %in% idx.obs)]

    # index of unobserved points sharing a y grid line with an observed one
    idx.unobs.y = idx.y[!(idx.y %in% idx.obs)]

    # index of unobserved points sharing no grid lines with observations
    idx.unobs.no = seq(nx*ny)[-c(idx.obs, idx.x, idx.y)]


    #### matrix computations

    # build component marginal covariance matrices for observations
    vx = pkern_corrmat(xpars, ngx)
    vy = pkern_corrmat(ypars, ngy)

    # build components of cross covariance matrices between observations and the first two subsets
    vx.cross = pkern_corrmat(xpars, dims[1], ds=1/dsx, j=gxc, i=gx)
    vy.cross = pkern_corrmat(ypars, dims[2], ds=1/dsy, j=gyc, i=gy)

    # TODO: replace explicit inverse with faster methods
    vx.inv = chol2inv( Rfast::cholesky(vx) )
    vy.inv = chol2inv( Rfast::cholesky(vy) )


  } else {

    # unpack precomputed objects
    if( is.list(pre) )
    {
      vx = pre[[1]]
      vy = pre[[2]]
      vx.cross = pre[[3]]
      vy.cross = pre[[4]]
      vx.inv = pre[[5]]
      vy.inv = pre[[6]]
      idx.unobs.x = pre[[7]]
      idx.unobs.y = pre[[8]]
      idx.unobs.no = pre[[9]]

    } else {

      # any non-list input to `pre` interpreted as prompt to precompute expensive stuff
      pre = list(vx, vy, vx.cross, vy.cross, vx.inv, vy.inv, idx.unobs.x, idx.unobs.y, idx.unobs.no)
      return(pre)
    }
  }

  # left-multiply inverse marginal covariance to get transformed observations
  ztrans = pkern_kprod(vx.inv, vy.inv, zobs, trans=TRUE)

  ## TESTING QR BASED METHOD ***********************************************

  # ztrans.t = c( qr.solve(qr(vx), Rfast::transpose( qr.solve(qr(vy), matrix(zobs, ngy)) ) ) )
  # ztrans = ztrans.t[ pkern_r2c(rev(c(ngx, ngy)), FALSE, TRUE) ]

  # works okay!

  ## ***********************************************


  ## TESTING EIGENVALUE-BASED METHOD ***********************************************

  evx = eigen(vx, symmetric=TRUE)
  evy = eigen(vy, symmetric=TRUE)
  z1 = pkern_kprod(evx[['vectors']], evy[['vectors']], zobs, trans=TRUE)
  z2 = z1 / ( nug + kronecker(evx[['values']], evy[['values']]) )
  ztrans = pkern_kprod(evx[['vectors']], evy[['vectors']], z2, trans=FALSE)

  # works okay!

  ## ***********************************************



  # compute predictions separately on the three subsets
  zout[idx.unobs.x] = pkern_kprod(vx, vy.cross, ztrans, trans=TRUE)
  zout[idx.unobs.y] = pkern_kprod(vx.cross, vy, ztrans, trans=TRUE)
  zout[idx.unobs.no] = pkern_kprod(vx.cross, vy.cross, ztrans, trans=TRUE)

  # finish
  return(zout)
}
