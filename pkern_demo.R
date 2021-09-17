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


pars = pkern_corr(c('mat', 'mat'))
pkern_plot(pars, ppars=list(glcol='black'))
pkern_plot(pars, dims=c(1e3, 1e3), ds=c(0.3, 100))
pkern_sim(pars)


pkern_sim(pars, dims=c(1e2, 1e2), ds=ds)



pars$nug = 0

pkern_sim(pars, dims=c(500, 500))


pkern_plot(sim, maxx=500)
kname = c('mat', 'sph')
pars = pkern_corr(kname)
pkern_toString(pars)
pkern_plot(modifyList(pars, list(nug=0.5)), smoothed=F)

# nx / ny
rr = res(dem)[1]/res(dem)[2]
rr = 1

# extract grid size, grid line locations, vectorized data
dims.dem = pkern_fromRaster(dem, 'dim')
xy.dem = pkern_fromRaster(dem, 'xy')
z.dem = pkern_fromRaster(dem, 'values')

demt = t(dem)
pkern_plot(z=demt, rr=1, maxx=300)
pkern_plot(z=demt, rr=1, maxx=300)

pkern_plot(z=z.dem, dims=dims.dem, xy=xy.dem, maxx=300)

plot(demt)





library(microbenchmark)
microbenchmark( pkern_plot(z.dem, dims.dem, maxx=700), times=5 )
microbenchmark( plot(dem, maxpixels=700^2), times=5 )

# snap weather point to a grid
gxy.snap = pkern_snap(pts=tmin, g=dem, regular=TRUE, makeplot=FALSE)
nm.data = 'tmin'

# plot the larger example with snapped grid lines
# plot(dem)
# plot(tmin, add=T, pch=16, pal=\(n) rainbow(n, rev=TRUE))
# abline(v=gxy.snap$x$gval)
# abline(h=gxy.snap$y$gval)
# par(old.par)


# development of new raster function

# distance between grid lines in outer grid
res.src = res(dem)

# extract grid line mapping within the bounding box for the snapped points
xymap = cbind(sapply(gxy.snap, \(d) d$id), id=seq(nrow(tmin)))
res.inc = sapply(gxy.snap, \(xy) unique(diff(xy$gid)))
dims = Rfast::colMaxs(xymap[,c('x', 'y')], value=TRUE)
ds = res.inc * res.src

# plot a rasterized version of the tmin values after centering
tmin.mean = mean(tmin[[nm.data]])
vec.src = rep(NA, prod(dims))
vec.src[ pkern_mat2vec(xymap, dims[2]) ] = tmin[[nm.data]] - tmin.mean
ppars = list(asp=ds[2]/ds[1], leg='tmin', main='weather data after snapping', glcol='grey90', xlab=NA, ylab=NA)
pkern_plot(vec.src, dims, ppars=ppars)

# # another version with crs etc preserved from original raster
# gext = sapply(gxy.snap, \(d) range(d$gval)) |> as.vector()
# r.template = extent(gext) |> raster(ncols=dims[1], nrows=dims[2], crs=crs(dem))
# r.src = pkern_toraster(vec.src, template=r.template)
# plot(r.src, col=rainbow(100))



#

# development of sample variograms


#ds = setNames(res(r.src), c('dx', 'dy'))
gxy = lapply(gxy.snap, \(d) d$gid)

# vario = pkern_vario(dims.bbox, vec.src, ds=ds, diagonal=FALSE, dmax=8e4)
# pkern_vario_plot(vario)
# pars = pkern_vario_fit(vario, xpars='mat', ninitial=10)
# pkern_vario_plot(vario, pars)
# pars



# sample variogram
vario = pkern_vario(dims, vec.src, ds=ds)
pkern_vario_plot(vario)

# fit a model to the sample variograms
pars = pkern_vario_fit(vario)
pkern_vario_plot(vario, pars)
pkern_plot(pars)

# compute conditional mean and plot
zpred = pkern_cmean(vec.src, dims=dims.dem, pars=pars, gxy=gxy)
zpred[is.na(zpred)] = mean(zpred, na.rm=T) # temporary: zpred has NAs!
pkern_plot(zpred, dims.dem)

# repeat with a restriction on semivariance distance
vario = pkern_vario(dims, vec.src, ds=ds, dmax=10e4)
pars = pkern_vario_fit(vario, 'mat', ninitial=50)
pkern_vario_plot(vario, pars)
pkern_plot(pars)

# compute conditional mean and plot
zpred = pkern_cmean(vec.src, dims=dims.dem, pars=pars, gxy=gxy)
zpred[is.na(zpred)] = mean(zpred, na.rm=T) # temporary: zpred has NAs!
pkern_plot(zpred, dims.dem)





# TODO: implement sparse matrix multiplication to replace Crossprod

# TODO: recursive call within cmean to fill the NAs
#which(is.na(zpred))

# TODO: find a better way of dealing with large plot calls and test on:
#pkern_kplot(pars)
#pkern_plot(zpred, dims.dem)

# downscaling algorithm
#rpred = pkern_toraster(zpred, dims.dem)
#plot(rpred, col=rainbow(1e5))

library(Matrix)

px = list(k='gau', kp=10)
py = list(k='gau', kp=10)
newj = seq(1, dims.dem[1], by=100)
newi = seq(1, dims.dem[2], by=100)
X = pkern_corrmat(px, dims.dem[1], j=newj)
Y = pkern_corrmat(py, dims.dem[2], j=newi)
X = sweep(X, 2, colSums(X), '/')
Y = sweep(Y, 2, colSums(Y), '/')

z = zpred
z[is.na(z)] = 0
# zout = pkern_kprod(X, Y, z, trans=TRUE)
# zmout = matrix(zout, length(newi))
# zr = pkern_toraster( zmout )
#plot(zr, col=hcl.colors(1e2, 'Spectral'))










library(profvis)
profvis(pkern_kprod(X, Y, z, trans=TRUE))





z = pkern_sim(pars, dims=c(90, 100))
pkern_plot(z, smoothed=F, ppars=list(asp=1, glcol=NA))

# generate a matrix of data
z = pkern_sim(pkern_corr('mat'), dims=c(90, 100))
#plot(pkern_toraster(z), col=hcl.colors(100,"YlOrRd", rev = TRUE))
#pkern_plot(z)
pkern_plot(z, smoothed=F, ppars=list(asp=1, glcol=NA))




# persian rug art (example from filled.contour docs)
# x <- y <- seq(-4*pi, 4*pi, length.out = 27)
# r <- sqrt(outer(x^2, y^2, "+"))
# filled.contour(cos(r^2)*exp(-r/(2*pi)), axes = FALSE)
#pkern_kplot(pars, dims.bbox)





#############


###################################################
# development of resampling function

library(profvis)
library(microbenchmark)

# pkern_plot(zmout, smoothed=F, ppars=list(asp=1, glcol=NA, leg='elevation (m)'))



