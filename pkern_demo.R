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
plot(dem)
plot(tmin, add=T, pch=16)
par(old.par)

# extract column-vectorized DEM data and grid properties
dem.pkern = pkern_fromRaster(dem)
dims = pkern_fromRaster(dem, 'dim')
gxy = pkern_fromRaster(dem, 'xy')

### testing

nstart = 100

# extract points coordinates
coords = st_coordinates(tmin)


# new functions
# pkern_sgrid()
#




gxy.snap = pkern_snap(pts=tmin, g=dem, regular=F)

plot(dem)
plot(tmin, add=T, pch=16)
abline(v = gxy.snap$x$gval)
abline(h = gxy.snap$y$gval)
par(old.par)

# example for docs
reso = c(10, 10)
nx = 25
ny = 30
xy = list(reso[1]*seq(nx), reso[2]*seq(ny))
coords = expand.grid(xy) + rnorm(nx*ny)
gxy.snap = pkern_snap(pts=coords, g=xy, regular=F)
plot(coords)
abline(v = gxy.snap$x$gval)
abline(h = gxy.snap$y$gval)



pkern_snap(coords[,'y'], gxy[['x']])



# order into row and column vectorized forms (with decreasing y, increasing x)
order.colvec = order(coords[,'y'], coords[,'x'], decreasing = c(TRUE, FALSE))
#order.rowvec = order(coords[,'x'], coords[,'y'], decreasing = c(TRUE, FALSE))

# snap coordinates to grid lines by least Manhattan distance
#gx.colvec = Rfast::Outer(coords[order.colvec,'x'], gxy[['x']], '-') |> abs() |> apply(2, which.min)
gy.colvec = Rfast::Outer(coords[order.colvec,'y'], gxy[['y']], '-') |> abs() |> apply(2, which.min)

# k-means clustering (k=2) of the (ordered) y grid line separations to find y resolution
kmeans.colvec = stats::kmeans(diff(gy.colvec), 2, nstart=nstart)
idx.lowest = which.min(kmeans.colvec$centers)

# cluster input coordinates as grid rows by binning into y grid lines
cluster.endpoints = c(which(kmeans.colvec$cluster != idx.lowest), n)
ncluster = length(cluster.endpoints)
n.bycluster = c(cluster.endpoints[1], diff(cluster.endpoints))
cluster.ids = do.call(c, lapply(seq(ncluster), function(id) rep(id, n.bycluster[id])))

# map to original coordinates
coords.id = cluster.ids[ match(seq(n), order.colvec) ]

# find regular snapped grid lines by nearest neighbour to median


sapply(split(gy.colvec, cluster.ids), stats::median)










cluster.endfactor = 1 + as.numeric( diff(kmeans.colvec$centers) > 0 )
cluster.endpoints = c(which(cluster.idx == cluster.endfactor), n)
ncol = length(cluster.breaks)
lapply( seq(cluster.endpoints), function(j) cluster.breaks[1]:cluster.breaks[2])





cumsum(idx.end)





idx.cluster = kmeans.colvec$cluster
idx.cluster = c(idx.cluster[1], idx.cluster)

diff(cluster.colvec)


center.colvec = as.vector(kmeans.colvec$centers)

idx.newcol = 1 + which(kmeans.colvec$cluster == 2)

idx.newcol




plot(seq(nrow(coords)), idx.y)



Rfast::Outer(gxy[['x']], coords[,'X'])

Rfast::rowMins( Rfast::Outer(gxy[['x']], coords[,'X']) )






# find their distance matrix and the range of interpoint distances
coords.dist = Rfast::Dist(coords)
coords.udist = coords.dist[upper.tri(coords.dist, diag=FALSE)]
ip.range = range(coords.udist)

# TODO: handle zeros here (coincident but separately recorded points)
ip.range[2]/ip.range[1]
exp(diff(log(ip.range)))


# given a target number of grid lines





idx.xorder = order(coords[,1])
plot(seq(nrow(coords)), coords[idx.xorder, 1])
plot(seq(nrow(coords)), sort(coords[,1]))



plot(seq(nrow(coords)), coords[idx.xorder, 2])




plot(seq_along(coords.udist), sort(coords.udist))






# need to snap at most one input point to grid point





upper.tri()






# snap weather grid to DEM

coords = tmin
gxy = dem

pkern_detect_grid(coords)

tmin.snap = pkern_snap_grid(tmin, dem)
tmin.snap
range(tmin.snap$snap$dist)
median(tmin.snap$snap$dist)

# plot the result
plot(dem)
abline(v = tmin.snap$gxy$x)
abline(h = tmin.snap$gxy$y)
plot(tmin, add=T, pch=16)


xx = rasterize(tmin, raster(extent(dem)), field='tmin')
yy = values(xx)

sum( is.na(yy) )

sum( !is.na(yy) )



# # test indexing
# j = tmin.snap$gji$x[tmin.snap$snap$j]
# i = tmin.snap$gji$y[tmin.snap$snap$i]
# zvec = dem.pkern$data
# zvec[ i + ( dem.pkern$dims[2] * (j-1) ) ] = NA
# pkern_toRaster(zvec, template=dem) %>% plot
# abline(v = tmin.snap$gxy$x)
# abline(h = tmin.snap$gxy$y)
# looks good

# example for docs
dev.off()
nx = 35
ny = 24
coords = expand.grid(seq(nx), seq(ny)) + rnorm(nx*ny, 0, 1/5)
plot(coords)
n = nrow(coords)
x = coords[,1] %>% as.vector
pkern_kmeans(x)

plot(coords)
#test.grid = pkern_detect_grid(coords, rnx=c(1,100), rny=c(1,100), dminx=0.0)
test.grid = pkern_detect_grid(coords)
test.grid$dims
abline(v = test.grid$gxy$x)
abline(h = test.grid$gxy$y)

nx = 5
ny = 5
coords = expand.grid(seq(nx), seq(ny)) + rnorm(nx*ny, 0, 1/8)
plot(coords)
test.grid = pkern_detect_grid(coords, cw=3, nstart=50)
test.grid$dims
abline(v = test.grid$gxy$x)
abline(h = test.grid$gxy$y)



plot(seq_along(coords[,1]), coords[,1])



abline(v = test.grid$gxy$x[5])
abline(h = test.grid$gxy$y[2])


plot(coords)
test.grid = pkern_detect_grid(coords, dmin=0.0)
test.grid$dims
abline(v = test.grid$gxy$x)
abline(h = test.grid$gxy$y)



