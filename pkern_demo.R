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


# generate a matrix of data
z = pkern_sim(pkern_corr('mat'), dims=c(30, 40))
#plot(pkern_toraster(z), col=hcl.colors(100,"YlOrRd", rev = TRUE))
#pkern_plot(z)
pkern_plot(z, smoothed=F)

lvls = pretty(z, 10)
cols = hcl.colors(length(lvls)-1, "YlOrRd", rev = TRUE)
.filled.contour(x=x, y=y, zmat, lvls, cols)

graphics::image(x=x, y=y, zmat, axes=FALSE, ann=FALSE, asp=1)

methods(image)
getAnywhere(image.default)

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
# plot(tmin, add=T, pch=16, pal=\(n) rainbow(n, rev=TRUE))
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


ds = setNames(res(r.src), c('dx', 'dy'))
gxy = lapply(gxy.snap, \(d) d$gid)

# vario = pkern_vario(dims.bbox, vec.src, ds=ds, diagonal=FALSE, dmax=8e4)
# pkern_vario_plot(vario)
# pars = pkern_vario_fit(vario, xpars='mat', ninitial=10)
# pkern_vario_plot(vario, pars)
# pars



# sample variogram
vario = pkern_vario(dims.bbox, vec.src, ds=ds)
pkern_vario_plot(vario)

# fit a model to the sample variograms
pars = pkern_vario_fit(vario)
pkern_vario_plot(vario, pars)

# compute conditional mean
zpred = pkern_cmean(vec.src, dims=dims, pars=pars, gxy=gxy)
rpred = pkern_toraster(zpred, template=dem)
plot(crop(rpred, r.src), col=rainbow(100))



# repeat with a restriction on semivariance distance
vario = pkern_vario(dims.bbox, vec.src, ds=ds, dmax=10e4)
pars = pkern_vario_fit(vario, 'mat', ninitial=50)
pkern_vario_plot(vario, pars)
pars

# compute conditional mean
zpred = pkern_cmean(vec.src, dims=dims, pars=pars, gxy=gxy)
rpred = pkern_toraster(zpred, template=dem)
plot(crop(rpred, r.src), col=rainbow(1e2))
plot(rpred, col=rainbow(1e5))

#dims = dims.bbox





# persian rug art (example from filled.contour docs)
# x <- y <- seq(-4*pi, 4*pi, length.out = 27)
# r <- sqrt(outer(x^2, y^2, "+"))
# filled.contour(cos(r^2)*exp(-r/(2*pi)), axes = FALSE)
#pkern_kplot(pars, dims.bbox)





#############





