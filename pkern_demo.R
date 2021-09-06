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
plot(dem)
plot(tmin, add=T, pch=16)
abline(v=gxy.snap$x$gval)
abline(h=gxy.snap$y$gval)
par(old.par)


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

# development of






res.src = res(dem) * sapply(gxy.snap, \(g) unique(diff(g$gid)))
res(r.src) = res.src


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







