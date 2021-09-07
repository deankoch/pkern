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

old.par = par(no.readonly=TRUE)

library(here)
library(raster)
library(sf)

# load some data - a raster DEM and temperature points geojson
path.testdata = 'D:/pkern/inst/testdata'

# smaller example
dem = raster(file.path(path.testdata, 'dem_big_creek.tif'))
tmin = st_read(file.path(path.testdata, 'tmin_big_creek.geojson'))
plot(dem)
plot(tmin, add=T, pch=16)
par(old.par)

# larger example
dem = raster(file.path(path.testdata, 'dem_full.tif'))
tmin = st_read(file.path(path.testdata, 'tmin_full.geojson'))
# plot(dem)
# plot(tmin, add=T, pch=16)
# par(old.par)

# extract grid properties
dims = pkern_fromraster(dem, 'dim')
gxy = pkern_fromraster(dem, 'xy')

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

# plot a rasterized version of tmin
vec.src = rep(NA, prod(dims.bbox))
vec.src[ pkern_mat2vec(xy, dims.bbox[2]) ] = tmin[[nm.data]]
r.src = pkern_toraster(vec.src, dims.bbox)
plot(r.src)

# another version with crs etc preserved from original raster
gext = sapply(gxy.snap, \(d) range(d$gval)) |> as.vector()
r.template = extent(gext) |> raster(ncols=dims.bbox[1], nrows=dims.bbox[2], crs=crs(dem))
r.src = pkern_toraster(vec.src, template=r.template)
plot(r.src)

#

# development of sample variograms
dims = dims.bbox
vec = vec.src
vario = pkern_vario(dims, vec, nmin=25)

#

# pars.x = list(k='sph', kp=10)
# pars.y = list(k='sph', kp=10)
pars.x = list(k='mat', kp=c(10, 2))
pars.y = list(k='mat', kp=c(10, 2))
pvec = c(1, pars.x$kp, pars.y$kp)

# now try again using nlm
fn.vario = function(pvec, pars.x, pars.y, vario, makeplot=FALSE)
{
  # determine number of parameters from each kernel
  npx = length(pars.x$kp)
  npy = length(pars.y$kp)

  # extract parameters from input vector
  sig = pvec[1]
  kpx = pvec[1 + seq(npx)]
  kpy = pvec[1 + npx + seq(npy)]

  # ignore lags with no data
  inclx = vario$x$n > 0
  incly = vario$y$n > 0
  nx = vario$x$n[inclx]
  ny = vario$y$n[incly]

  # extract empirical semivariance data on dimension x
  sepx = vario$x$sep[inclx]
  svx = vario$x$semivariance[inclx]

  # and for y
  sepy = vario$y$sep[incly]
  svy = vario$y$semivariance[incly]

  # generate x and y kernel values along each lag
  semi.x = 2 * sig * (1 - pkern_corr(sepx, modifyList(pars.x, list(kp=kpx))) )
  semi.y = 2 * sig * (1 - pkern_corr(sepy, modifyList(pars.y, list(kp=kpy))) )

  # compute weighted sums of squares along both dimensions and return their sum
  wss.x = sum( nx * ( svx - semi.x )^2 )
  wss.y = sum( ny * ( svy - semi.y )^2 )

  if(makeplot)
  {
    old.par = par(no.readonly=TRUE)
    par(mfrow=c(1,2))

    plot(sepx, svx, xlab='lag', ylab='semivariance', main='x')
    lines(sepx, semi.x)
    plot(sepy, svy, xlab='lag', ylab='semivariance', main='y')
    lines(sepy, semi.y)

    par(old.par)
  }

  return( wss.x + wss.y )
}

result.nlm = nlm(fn.vario, p=pvec, pars.x=pars.x, pars.y=pars.y, vario=vario, iterlim=100)

# unpack results
npx = length(pars.x$kp)
npy = length(pars.y$kp)
sig = result.nlm$estimate[1]
kpx = result.nlm$estimate[1 + seq(npx)]
kpy = result.nlm$estimate[1 + npx + seq(npy)]

# use another method
# bds.lower = c(0, 0, 0)
# bds.upper = c(sqrt(var(vec, na.rm=TRUE)), 100, 100)
bds.lower = c(0, 1, 0.1, 1, 0.1)
bds.upper = c(sqrt(var(vec, na.rm=TRUE)), 1e3, 10, 1e3, 10)
result.optim = stats::optim(par=pvec,
                            f=fn.vario,
                            method='L-BFGS-B',
                            lower=bds.lower,
                            upper=bds.upper,
                            pars.x=pars.x,
                            pars.y=pars.y,
                            vario=vario)


fn.vario(result.optim$par, pars.x, pars.y, vario, makeplot=TRUE)



# extract separation distances and semivariances
idx.incl = vario$x$n > 0
sep = vario$x$sep[idx.incl]
sv = vario$x$semivariance[idx.incl]
n = vario$x$n[idx.incl]

# try fitting a covariogram to the y component
pars.y = list(k='sph', kp=10)
bds.sig = c(0, sqrt(var(vec, na.rm=TRUE)))

bds.kpy = list(rho=c(0, 90))

# test grids
n.each = 100
sig.test = seq(bds.sig[1], bds.sig[2], length.out=n.each)
kpy.test = sapply(bds.kpy, \(bds) seq(bds[1], bds[2], length.out=n.each))
wss.test = rep(NA, n.each^2)
for( idx.sig in seq(n.each) )
{
  sig = sig.test[idx.sig]

  for( idx.kpy in seq(n.each) )
  {
    idx.test = idx.kpy + n.each * (idx.sig - 1)

    kpy = kpy.test[idx.kpy,]

    pars = modifyList(pars.y, list(kp=kpy))

    wss.test[idx.test] = sum( n * ( ( sv - 2*sig*(1 - pkern_corr(sep, pars) ) )^2 ) )
  }
}

idx.best = pkern_vec2mat(which.min(wss.test), n.each)
idx.kpy = idx.best[1]
idx.sig = idx.best[2]
kpy = kpy.test[idx.kpy,]
pars = modifyList(pars.y, list(kp=kpy))

plot(sep, sv)
lines(sep, 2*sig*(1 - pkern_corr(sep, pars) ))

pars.x = list(k='gxp', kp=c(10,2))
pars.y = list(k='sph', kp=10)
pvec = c(1, pars.x$kp, pars.y$kp)

# now try again using nlm
fn.vario = function(pvec, pars.x, pars.y, vario)
{
  # determine number of parameters from each kernel
  npx = length(pars.x$kp)
  npy = length(pars.y$kp)

  # extract parameters from input vector
  sig = pvec[1]
  kpx = pvec[1 + seq(npx)]
  kpy = pvec[1 + npx + seq(npy)]

  # extract empirical semivariance data on dimension x
  inclx = vario$x$n > 0
  nx = vario$x$n[incl.x]
  sepx = vario$x$sep[incl.x]
  svx = vario$x$semivariance[inclx]

  # and for y
  incly = vario$y$n > 0
  ny = vario$y$n[incl.y]
  sepy = vario$y$sep[incl.y]
  svy = vario$y$semivariance[incly]

  # generate x and y kernel values along each lag
  kern.x = pkern_corr(sepx, modifyList(pars.x, list(kp=kpx)))
  kern.y = pkern_corr(sepy, modifyList(pars.y, list(kp=kpy)))

  # compute weighted sums of squares along both dimensions and return their sum
  wss.x = sum( nx * ( (svx  - 2 * sig * (1 - kern.x) )^2 ) )
  wss.y = sum( ny * ( (svy  - 2 * sig * (1 - kern.y) )^2 ) )
  return( wss.x + wss.y )
}

nlm(fn.vario, p=pvec, pars.x=pars.x, pars.y=pars.y, vario=vario)

# determine number of parameters from each kernel
npx = length(pars.x$kp)
npy = length(pars.y$kp)

# extract parameters from input vector
sig = pvec[1]
kpx = pvec[1 + seq(npx)]
kpy = pvec[1 + npx + seq(npy)]









# try fitting a correlation model to this data
pars.x = list(k='sph', kp=10)
pars.y = list(k='sph', kp=10)
bds.sig = c(0, sqrt(var(vec, na.rm=TRUE)))
bds.kpx = list(rho=c(0, 100))
bds.kpy = list(rho=c(0, 90))
n.each = 100
sig.test = seq(bds.sig[1], bds.sig[2], length.out=n.each)
kpx.test = lapply(bds.kpx, \(bds) seq(bds[1], bds[2], length.out=n.each))
kpy.test = lapply(bds.kpy, \(bds) seq(bds[1], bds[2], length.out=n.each))






2*(sig - pkern_corr(vario$x$sep, pars) )






pars= list(k='gex', kp=c(p=2, rho=100))
pkern_corr(1:10, pars)
pars= list(k='gau', kp=100)
pkern_corr(1:10, pars)
pars= list(k='sph', kp=100)
pkern_corr(1:10, pars)
pars= list(k='mat', kp=c(5, 3))
pkern_corr(1:10, pars)

pkern_xvario(dims.bbox, vec, simple=FALSE)

# sample variogram in y direction
pkern_xvario(rev(dims.bbox), vec[pkern_r2c(dims, FALSE, TRUE)])

install.packages('usdm')
library(usdm)

Variogram(r.src) |> str()




res.src = res(dem) * sapply(gxy.snap, \(g) unique(diff(g$gid)))
res(r.src) = res.src
load_all()

#
# # extract i,j indices as dataframe with NA values for missing grid points
# ij.all = expand.grid( setNames(sapply(rev(dims.bbox), seq), c('i', 'j')) )
# ij.all$x = ij.all$j
# ij.all$y = dims.bbox[2] - ij.all$i
# ij.all = merge(ij.all, as.data.frame(xy), all.x=TRUE)
#
# # make a plot to check if it looks okay
# idx.src = ij.all[['id']]
# vec.src = tmin[[nm.data]][ idx.src ]
# #vec.src[is.na(vec.src)] = mean(tmin[[nm.data]])
# pkern_toRaster(vec.src, dims.bbox) |> plot()







