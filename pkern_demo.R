#' ---
#' title: "pkern_demo"
#' author: "Dean Koch"
#' date: "`r format(Sys.Date())`"
#' output: github_document
#' ---
#'
#' **Mitacs UYRW project**
#'
#' **pkern_demo**: code snippets to use in a vignette
#'

library(devtools)
load_all()
document()

old.par = par(no.readonly=TRUE)


library(here)
library(raster)
library(sf)

# load some data - a raster DEM and temperature points geojson
path.testdata = 'D:/pkern/inst/testdata'

# smaller example
# dem = raster(file.path(path.testdata, 'dem_big_creek.tif'))
# tmin = st_read(file.path(path.testdata, 'tmin_big_creek.geojson'))
# plot(dem)
# plot(tmin, add=T, pch=16)
# par(old.par)

# larger example
dem = raster(file.path(path.testdata, 'dem_full.tif'))
tmin = st_read(file.path(path.testdata, 'tmin_full.geojson'))
# plot(dem)
# plot(tmin, add=T, pch=16)
# par(old.par)

# extract grid size
dims = pkern_fromraster(dem, 'dim')

# snap weather point to a grid
gxy.snap = pkern_snap(pts=tmin, g=dem, regular=TRUE, makeplot=FALSE)
nm.data = 'tmin'

# plot the larger example with snapped grid lines
# plot(dem)
# plot(tmin, add=T, pch=16)
# abline(v=gxy.snap$x$gval)
# abline(h=gxy.snap$y$gval)
# par(old.par)


# development of new raster function






# extract grid line mapping within the bounding box for the snapped points
xy = cbind(sapply(gxy.snap, \(d) d$id), id=seq(nrow(tmin)))
dims.bbox = Rfast::colMaxs(xy[,c('x', 'y')], value=TRUE)
res.inc = sapply(gxy.snap, \(xy) unique(diff(xy$gid)))

# plot a rasterized version of the tmin values after centering
tmin.mean = mean(tmin[[nm.data]])
vec.src = rep(NA, prod(dims.bbox))
vec.src[ pkern_mat2vec(xy, dims.bbox[2]) ] = tmin[[nm.data]] - tmin.mean
src = pkern_toraster(vec.src, dims.bbox)
plot(src, col=rainbow(100))

#

# another version with crs etc preserved from original raster
gext = sapply(gxy.snap, \(d) range(d$gval)) |> as.vector()
r.template = extent(gext) |> raster(ncols=dims.bbox[1], nrows=dims.bbox[2], crs=crs(dem))
r.src = pkern_toraster(vec.src, template=r.template)
plot(r.src, col=rainbow(100))



#

# development of sample variograms
vario = pkern_vario(dims.bbox, vec.src, dmax=17)
pkern_vario_plot(vario)

# fit a model
# xpars = pkern_corr('gau')
# xpars$lower = c(1e-1)
# xpars$upper = c(100)
# ypars = xpars
#
# xpars = 'mat'
# ypars = 'mat'

dims = pkern_fromraster(dem, 'dim')
pars = pkern_vario_fit(vario, 'mat')
pkern_vario_plot(vario, pars)
pars

# persian rug art (example from filled.contour docs)
# x <- y <- seq(-4*pi, 4*pi, length.out = 27)
# r <- sqrt(outer(x^2, y^2, "+"))
# filled.contour(cos(r^2)*exp(-r/(2*pi)), axes = FALSE)

pkern_kplot(pars, dims.bbox)


# compute conditional mean
gxy = lapply(gxy.snap, \(d) d$gid)
zobs = rep(0, prod(dims.bbox))
zobs[ pkern_mat2vec(xy, dims.bbox[2]) ] = tmin[[nm.data]] - tmin.mean
#pars$x$kp = c(0.5, 0.74)
#pars$y$kp = c(20, 3.12)
zpred = pkern_cmean(dims, xpars=pars[['x']], ypars=pars[['y']], zobs=zobs, gxy=gxy, nug=pars[['nug']])
rpred = pkern_toraster(zpred, template=dem)
#plot(crop(rpred, r.src), col=rainbow(100))

plot(rpred, col=rainbow(100))


library(microbenchmark)
microbenchmark(pkern_cmean(dims, xpars=pars[['x']], ypars=pars[['y']], zobs=zobs, gxy=gxy, nug=pars[['nug']]))


# notes:
#
# covariance matrices become numerically singular very quickly for large rho with the gau
# and mat models. It would be good to have a function that tests condition number for range
# of rho values, and setting bounds accordingly to avoid this problem
#
# semivariogram fitting function works pretty well if nmin is set high enough.

#

# development of eigenvalue decomposition based methods

# example problem
px = pkern_corr('gau', 10)
py = pkern_corr('gau', 10)
nx = 24
ny = 27
vx = pkern_corrmat(px, nx)
vy = pkern_corrmat(px, ny)
z = rnorm(nx * ny)

# naive evaluation
vx.inv = chol2inv(chol(vx))
vy.inv = chol2inv(chol(vy))
ans = c( kronecker(vx.inv, vy.inv) %*% c(z) )

# simplified version, still using inverses
ans2 = c( vy.inv %*% matrix(z, ny) %*% vx.inv )
max(abs(ans-ans2))

# evalue approach
evx = eigen(vx, symmetric=TRUE)
evy = eigen(vy, symmetric=TRUE)
z1 = pkern_kprod(evx[['vectors']], evy[['vectors']], z, trans=TRUE)
z2 = z1 / kronecker(evx[['values']], evy[['values']])
ans3 = pkern_kprod(evx[['vectors']], evy[['vectors']], z2, trans=FALSE)
max(abs(ans-ans3))






#

# development of Q-R decomposition based methods

# example problem
px = pars[['x']]
py = pars[['y']]
nx = 7
ny = 5
vx = pkern_corrmat(px, nx)
vy = pkern_corrmat(px, ny)
z = rnorm(nx * ny)

# naive evaluation
vx.inv = chol2inv(chol(vx))
vy.inv = chol2inv(chol(vy))
ans = c( kronecker(vx.inv, vy.inv) %*% c(z) )

# simplified version, still using inverses
ans2 = c( vy.inv %*% matrix(z, ny) %*% vx.inv )
max(abs(ans-ans2))

# qr based approach to the left multiplication
ans3 = c( vy.inv %*% matrix(z, ny) )
ans4 = qr.solve(qr(vy), matrix(z, ny))
max(abs(ans3 - ans4))

# and now the right
ans5.unsorted = c( qr.solve(qr(vx), Rfast::transpose( qr.solve(qr(vy), matrix(z, ny)) ) ) )
ans5 = ans5.unsorted[ pkern_r2c(rev(c(nx, ny)), FALSE, TRUE) ]
max(abs(ans5 - ans2))









