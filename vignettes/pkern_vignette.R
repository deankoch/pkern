#' ---
#' title: "pkern_vignette"
#' author: "Dean Koch"
#' date: "`r format(Sys.Date())`"
#' output: github_document
#' ---
#'
#' **Mitacs UYRW project**
#'
#' **pkern**: Fast kriging on gridded datasets
#'
#' Standard kriging methods don't scale well on problems with high resolution and large geographical
#' extent, such as when downscaling weather data. `pkern` gets around this issue using product kernels
#' to simplify kriging computations and substantially reduce their memory and CPU demands.
#' We use a combination of new methods and ones described in Gilboa et. al (2015), Koch et. al (2021).
#'
#' This vignette uses `pkern` to interpolate soil data from the "meuse" dataset (included
#' with `gstat`), showing how to use the main functions in the package in a
#' kriging workflow optimized for computational simplicity:
#'
#' * `pkern_fromRaster` defines a grid based on an existing raster (optional)
#' * `pkern_snap` snaps input point data to a regular subgrid
#' * `pkern_vario` and `pkern_vario_fit` fit a separable spatial covariance model to the data
#' * `pkern_cmean` computes conditional mean (kriging estimates) on full grid
#' * `pkern_variance` computes kriging variance (pointwise)
#'
#' ## Getting started
#'
#' `pkern` requires two additional packages: `Matrix` and `RcppHungarian`. The first helps reduce
#' memory usage, and the second allows snapping of input points to distinct grid points. Make sure these
#' can be loaded before using `pkern`:

# dependencies of the pkern package
library(Matrix)
library(RcppHungarian)

#' Install `pkern` by running:

#install.packages('pkern')
#library(pkern)

#' TEMPORARY: I'm still developing the package so for now the above lines above are commented
#' and I instead load all package functions by sourcing the following R files:

# source code for development version of the pkern package
library(here)
source(here('R/pkern_computation.R'))
source(here('R/pkern_covariance.R'))
source(here('R/pkern_indexing.R'))
source(here('R/pkern_plot.R'))
source(here('R/pkern_raster.R'))
source(here('R/pkern_snap.R'))
source(here('R/pkern_variograms.R'))

#' The `sf` and `raster` packages are also recommended for loading your spatial data
#' (they are used in the examples below), but they are not strictly required by `pkern`.

# some libraries for opening geospatial data
library(raster)
library(sf)

#' Finally, we need the `gstat` package to load the "meuse" dataset
library(gstat)

#'
#' ## Example data
#'
#' We will look at the "meuse" dataset (Pebesma, 2009) on heavy metal concentrations in the Meuse
#' river floodplain in the Netherlands. These data are included with the `gstat` package that we just
#' loaded (see this [tutorial pdf](https://cran.r-project.org/web/packages/gstat/vignettes/gstat.pdf)),
#' so they can now be accessed with the following call:

# load meuse dataframe into variable named `meuse`
data(meuse)
str(meuse)

#' Kriging assumes that a multivariate Gaussian process is generating the data. In our example the input data
#' will be zinc concentrations, which are strictly non-negative and skewed, so we first log-transform and
#' recenter the data to get something closer to a Gaussian random variable.

# histogram of the input data
hist(meuse$zinc, main='histogram of [zinc] in the meuse dataset')

#' This chunk does the transformation and appends it to `meuse` in the new column "logzinc". We save
#' the centering constant in variable `logzinc.mean` so we can invert the transformation later on.

# write new variable logzinc and copy the intercept
logzinc = scale(log(meuse$zinc), center=TRUE, scale=FALSE)
logzinc.mean = attr(logzinc, 'scaled:center')
meuse$logzinc = as.numeric(logzinc)

# histogram of the transformed data
hist(meuse$logzinc, main='histogram of zero-centred log[zinc] in the meuse dataset')

#' The `meuse` dataframe has columns "x" and "y" for the projected Easting and Northing coordinates,
#' so it can be converted to a spatial points object in R. Here's how this is done with package `sf`

# EPSG coordinate reference system code for meuse data
meuse.epsg = 28992

# convert meuse to an sf points object then plot log[zinc]
meuse.sf = st_as_sf(meuse, coords=c('x', 'y'), crs=meuse.epsg)
plot(meuse.sf['logzinc'], pch=16, main='zero-centred log[zinc]')

#' The river Meuse runs diagonally across this map (from bottom-left to top-right), and points
#' are sampled along its east bank. This vignette will demonstrate how to interpolate the zinc
#' data onto a high resolution grid covering the whole plot window.
#'
#' ## Prediction grid
#'
#' `pkern` needs a regular grid to predict over - in the examples below we will define one of
#' size 1000 x 1000 (dimensions `gdim`). The grid is positioned to cover the extent of the
#' input soil data coordinates, with coordinates of each grid line ("y" and "x") calculated
#' and stored in the list `gxy` below:

# define grid size
gdim = c(1e3, 1e3) |> setNames(c('y', 'x'))

# extract input points extent, switching from x-y (order returned by package sf) to y-x order
gext = apply(st_coordinates(meuse.sf), 2, range, simplify=FALSE) |> rev() |> setNames(c('y', 'x'))

# calculate resolution (distance between adjacent grid lines in y and x dimensions)
gres = Map(\(ext, d) diff(ext) / (d-1), ext=gext, d=gdim) |> setNames(c('y', 'x'))

# calculate equally spaced grid line positions
gyx = Map(\(r, g) seq(r[1], r[2], by=diff(r)/g), r=gext, g=gdim-1) |> setNames(c('y', 'x'))

#' We use `setNames` here to clarify the ordering of dimensions - all `pkern` objects store their
#' dimension-wise quantities in the order y, x. For example, the first entry of `gyx` is a vector of y
#' coordinates giving the (vertical) positions of the "y grid lines".

# print the first few entries of the first vector in gyx (y coordinates)
gyx[[1]] |> head()


#'
#' ## Ordering in pkern
#'
#' Why the order y, x? When I look at a raster image plot I see a matrix, and when we write out matrices
#' the vertical dimension, indexed by "i", is presented first in notation. Visually, "i" is the analogue of
#' "y" and therefore the usual x, y order seems backwards to me. Users should be aware that packages
#' `raster` and `sf` use the opposite convention: x, then y. This is why, for example we use `rev` after
#' `sf::st_coordinates` in the definition of `gres` above.
#'
#' Similarly, data values in `gfull` are stored in column-vectorized order, with y coordinates
#' decreasing and x coordinates increasing. This matches base R behaviour for vectorizing matrices (as in
#' `base::c` and `base::as.vector`). However, various other R packages (like `graphics` and `raster`) use
#' different vectorization orderings - if you are needing to convert between them, see `?pkern_r2c`.
#'
#' The average user shouldn't have to worry much about this, as most of the important functions in `pkern`
#' support `RasterLayer` and `sf` class arguments, and objects like `gres` and `gyx` are computed internally
#' in the correct order. In particular `pkern_fromRaster` and `pkern_toRaster` can be used
#' to switch to and from `RasterLayer` class, as we demonstrate next.
#'
#' ## Raster data
#'
#' I expect many users will already have an existing raster file (such as a DEM) that they want to predict
#' onto. These data are often loaded into R using the `raster` package. The following code shows how to
#' get the grid info from above (`gyx` etc) from a `RasterLayer`:

# make an example raster by supplying grid line coordinates (in x-y order!)
gfull.r = raster::rasterFromXYZ(expand.grid(rev(gyx)), crs=meuse.epsg)

# convert RasterLayer to `pkern`-style grid definition list and print structure
gfull = pkern_fromRaster(gfull.r)
str(gfull)

# convert back to RasterLayer
pkern_toRaster(gfull)


#'
#' ## PROJ4 warnings
#'
#' You'll notice the calls above issued a warning about datum strings. This is related to
#' R's [PROJ6/GDAL3 migration](https://cran.r-project.org/web/packages/rgdal/vignettes/PROJ6_GDAL3.html)
#' which deprecates the old PROJ4 strings in favour of modern well-known-text (WKT) and EPSG codes, creating
#' some dependency issues that are not completely resolved yet.
#'
#' Coordinate reference metadata is handled in `pkern_fromRaster(r)` by copying the WKT string from
#' `raster::wkt(r)` to the list element "crs". `pkern_toRaster` then passes this string to argument
#' "crs" of `raster:raster`, which appears to convert it to PROJ4 at some point, prompting a warning about
#' a dropped field. As PROJ4 is being phased out I imagine that updates to `raster` and `sp`
#' will soon make this warning disappear. For now (in this vignette) it can be ignored.
#'
#'
#' ## Snapping inputs to subgrid
#'
#' A key requirement of `pkern` is that sampled points lie on a regular subgrid of the full
#' prediction grid (`gfull`), where smaller subgrids (ie fewer grid lines) are to be preferred for
#' computational efficiency.  This means users should seek the smallest possible regular subgrid
#' while limiting positional errors to an acceptable level.
#'
#' `pkern_snap` finds this subgrid and snaps points to it. Argument `sep=c(y, x)` specifies the subgrid
#' spacing, in terms of the factor(s) by which to multiply the outer grid resolution to get the subgrid
#' resolution. For example, minimal spacing is `sep=c(1,1)`, which specifies to make the subgrid equal
#' to the full grid.
#'
#' In other words we simply snap to the nearest grid point, as in the following:

# snap zinc data to nearest grid point (`sep=1` is short for `sep=c(1,1)`)
gsnap.nearest = pkern_snap(gfull, pts=meuse.sf, sep=1, distinct=FALSE)

# print the range of snapping distances
gsnap.nearest$sg$pdist |> range()

# plot snapping pattern using helper function
pkern_snap_plot(gsnap.nearest, meuse.sf)

#' This minimizes total positional error (the largest snapping distance is about 2m) but results
#' in the most complicated subgrid possible (the full grid itself). This will lead to sluggish
#' computations and high memory demands.
#'
#' We can simplify by considering sparser subgrids: For example, the code below sets `sep=c(10,10)`,
#' forcing the snapping function to select a subgrid in which the spacing between adjacent subgrid
#' lines is equal to 10x the spacing in the outer grid - ie only every 10th grid line is used.

# snap zinc data to nearest subgrid point (10X coarser than full grid)
gsnap.nearest.10 = pkern_snap(gfull, pts=meuse.sf, sep=10, distinct=FALSE)

# print the sep values - the dimensions of the subgrid
print(gsnap.nearest.10$sg$gdim)

# and print snapping distance range
gsnap.nearest.10$sg$pdist |> range()

# plot snapping distance
pkern_snap_plot(gsnap.nearest.10, meuse.sf)


#' With wider spacing the positional error has increased - the largest snapping distance is now around 30m - but the
#' subgrid is now substantially simpler (only 100 grid lines are used along each dimension instead of 1000).
#' When `sep` is not supplied, `pkern_snap` will suggest a value automatically
#' (see `?pkern_estimate_sep`).

# snap with sep auto-detected, plot, and print snapping info
gsnap.nearest.auto = pkern_snap(gfull, pts=meuse.sf, distinct=FALSE)
pkern_snap_plot(gsnap.nearest.auto, meuse.sf)
gsnap.nearest.auto$sg$pdist |> range()
print(gsnap.nearest.auto$sg$gdim)

#' The function has selected a 35 X 29 subgrid with maximum positional error around 70m.
#' Users should of course verify that the suggested `sep` values are reasonable
#' by examining the snapping distances returned by `pkern_snap` and/or by visual inspection with
#' `pkern_snap_plot`.
#'
#' Looking closely at this plot, we can spot multiple input data points getting mapped to the
#' same grid point. One could do some kind of averaging or make a heirarchical model to deal with
#' co-located points, but there is a simpler way: `pkern_snap` includes the argument
#' `distinct` (TRUE by default), which specifies to map input points to distinct grid points.
#' The Hungarian algorithm (in package `RcppHungarian`) is used to solve the assignment
#' problem where total snapping distance is minimized under the distinctness constraint, as in
#' the following example:

# snap with no duplicates
gsnap = pkern_snap(gfull, pts=meuse.sf)
pkern_snap_plot(gsnap, meuse.sf)
gsnap$sg$pdist |> range()

#' Notice the maximum positional error has only increased marginally (from 72m to 76m), but now all
#' input points are mapped to a distinct grid point. The subgrid dimensions remain relatively small
#' (35 X 29), and this will make kriging computationally simple.

# print subgrid size
print(gsnap$sg$gdim)

#' ## Mapping subgrid points to input data
#'
#' The list (`gsnap`) returned by `pkern_snap` contains some useful indexing vectors that map input points
#' to the full grid as well as to/from the subgrid. For example the following code shows two ways to
#' represent the zinc data as a raster data vector containing the snapped values.

# indexing vectors mapping input points to and from the subgrid
idx.to = gsnap[['sg']][['pvec']]
idx.from = gsnap[['sg']][['ipvec']]

# initialize a subgrid data vector with NAs
zinc.sg = rep(NA, prod(gsnap[['sg']][['gdim']]))

# copy zinc data at snapped grid points
zinc.sg[idx.to] = meuse$logzinc

# plot the result
pkern_plot(zinc.sg, gdim=gsnap[['sg']][['gdim']], ppars=list(leg='log[zinc]'))

# you can also create the data vector in one line using inverse mapping
zinc.sg.compare = meuse$logzinc[idx.from]
identical(zinc.sg, zinc.sg.compare)

#' Note that the axis labels in the plot above refer to grid line numbers (not coordinates). To get
#' a plot in the original coordinate system, we can use the resolution and grid line position info
#' in `gsnap$sg`. The code below makes a copy of `gsnap` with grid data `zinc.sg`. This list is then
#' passed to the plotter function, which displays the grid with the correct coordinates
#' and aspect ratio

# make a copy of the `gsnap$sg`, adding values field, pipe to plotter
modifyList(gsnap[['sg']], list(gval=zinc.sg)) |> pkern_plot(ppars=list(leg='log[zinc]'))

#' The plotter function `pkern_plot` uses `graphics` to draw grid data, and is similar to
#' `raster::plot.raster` (but does not depend on it).

#'
#' ## Covariance models
#'
#' Before we can do kriging, we need a spatial covariance model. `pkern` requires that
#' we use a separable model (see Koch et al, 2020). In short, these models have
#' an x component and a y component, each of which is a covariance function of distance along
#' one dimension only.
#'
#' Component covariance functions are defined by starting from the 1d correlation function
#' (correlogram). In the example below we set up and visualize an example Gaussian correlogram using
#' the correlation function `pkern_corr`

# Set up a sequence of spatial lags
n = 25
d = seq(n)

# define a Gaussian correlation kernel with range 10
rho = 10
pars.1d = list(k='gau', kp=rho)

# plot correlation function
corr = pkern_corr(pars.1d, d)
plot(d, corr, pch=NA, main='1d Gaussian correlation function')
lines(d, corr)

# correlation matrix for the lag sequence above
pkern_corrmat(pars.1d, n) |> str()

#' The correlation matrices returned by `pkern_corrmat` are always given as (sparse) `Matrix`
#' class, which can be helpful on bigger problems (even though here the matrix is not sparse).
#'
#' A separable 2d correlogram is defined as a list with entries "y" and "x", each
#' a 1d kernel parameter list like `pars.1d`. These two component kernels need not have the same parameters
#' (or be from the same function family). For example, the code below doubles the range parameter
#' for the y component and plots the resulting 2d correlation pattern as a heatmap

# define the 1d component kernels
pars.1d.x = pars.1d
pars.1d.y = modifyList(pars.1d, list(kp=2*rho))

# define the 2d separable kernel and visualize correlation footprint
pars.2d = list(y=pars.1d.y, x=pars.1d.x)
pkern_plot(pars.2d)

#' This plot shows the correlation between a central point and all of the grid points
#' surrounding it, illustrating the spatial pattern of correlations under the model.
#' If you drew a transect from the middle of the heatmap to one of its edges and collected
#' correlation values, you would get the correlogram for that direction, similar to the
#' 1d correlation plot above.
#'
#' Notice that the heatmap is not radially symmetric, which means depending on the chosen
#' direction of your transect, the correlogram will look different. This is called
#' "anisotropy". All separable models are anisotropic, with the unique exception
#' of the Gaussian with equal x and y range parameters.
#'
#' To get a covariance function from a correlation function, you will need to specify the
#' partial sill parameter (`pars$psill`), and optionally a nugget variance (`pars$nug`)

# set a nugget effect and partial sill, plot semivariogram
pars.2d = modifyList(pars.2d, list(psill=2, nug=0.5))
pkern_plot(pars.2d)

#' The plot looks identical except that the color legend has a new range. The spatial
#' pattern of the covariance function is determined by the component correlograms, whereas the
#' nugget effect and partial sill essentially just introduce a scaling and non-zero intercept.
#'
#' ## Simulations
#'
#' The simulation function `pkern_sim` is another way to visualize a covariance model.
#' It generates random vectors from the specified Gaussian process which can be passed to the plotter
#' to get a sense of what type of spatial patterns are expected under the model:

# set random seed so we don't change the image every time I edit this script
set.seed(7)

# simulate data from the model defined above
pkern_sim(pars.2d)

# repeat with smaller nugget effect
set.seed(7)
pkern_sim(modifyList(pars.2d, list(nug=1e-2)))

#' The separable covariance model makes these computations pretty fast. The following
#' generates a simulation on the full (1 million pixel) grid in about 3 seconds on my desktop.

# repeat, this time setting desired grid dimensions
set.seed(7)
pkern_sim(modifyList(pars.2d, list(nug=1e-2)), gdim)

#' We will see later on that conditional mean (kriging) can be done even faster

#'
#' ## Likelihood
#'
#' In the next sections we will show how to fit the parameters of a covariance model to a
#' sample variogram by weighted least squares. This is not the most robust way of fitting covariance,
#' however it is easy to understand, and fast. For users preferring likelihood based-methods, we
#' also include the (log) likelihood function computer `pkern_LL`, which can be used in
#' numerical optimization.

# likelihood for the example model, given the meuse zinc data
pkern_LL(pars.2d, zinc.sg, gsnap[['sg']][['gdim']])


#' ## Sample variograms
#'
#' The function `pkern_vario` estimates sample variogram values
#' for the the input point dataset, after snapping to the grid. It produces four
#' sample variograms, corresponding to spatial lags along the x direction, the y
#' direction, and the two "diagonal" directions. These four directions are selected
#' for being computationally easy to work with when the sample data lie on a grid.

# compute sample variogram
vario = pkern_vario(gdim=gsnap, vec=zinc.sg, quiet=TRUE)

# a helper function can be used to plot the results
pkern_vario_plot(vario)

#' As an aside, these directional sample variograms reveal anisotropy in the random field. The
#' diagonal x direction (corresponding to the top-right to bottom-left direction) has a much lower
#' sill than the other directions, because we have not de-trended the data by accounting for
#' distance to the river Meuse - this will be illustrated in another vignette
#'
#' ## Fitting to sample variograms
#'
#' The function `pkern_vario_fit` can now be used to fit the variogram data to
#' a model by weighted least squares. The code below does so, and plots the resulting model
#' curves (expected variograms) over the scatterplots

# fit a separable (Gaussian) spatial covariance model to the variogram data
pars.vario = pkern_vario_fit(vario)

# plot the variograms
pkern_vario_plot(vario, pars.vario)

#' Here is what the fitted kernel looks like as a heatmap

# plot the fitted kernel
pkern_plot(pars.vario)

#' The fitted kernel is slightly anisotropic, with an effective range around 500m.
#'
#' By default `pkern_vario_fit` uses a Gaussian covariance model, the same one used in
#' the examples above, but any separable covariance family supported by `pkern_corr`
#' can be specified (along with parameter bounds) in arguments `ypars` and `xpars`.
#'
#' Users may also wish to set the `dmax` argument to ignore lags above a certain distance
#' (eg. 1km looks appropriate here), or the `fit.method` argument to change the method used
#' to weight sample variogram values. See `?pkern_vario_fit` for more details.
#'
#'
#' ## Conditional expectation
#'
#' The kriging function `pkern_cmean` takes the column vectorized subgrid data (in this case, `zinc.sg`),
#' the full grid dimensions (`gdim`), and the subgrid indexing vectors (`map`), and computes the
#' conditional mean at all grid points, given the sampled input data.
#'
#' As a shortcut, users can supply the return list of `pkern_snap` to `gdim` and the function will
#' extract `map` automatically:

# kriging to get conditional mean
zpred = pkern_cmean(zinc.sg, gsnap, pars.vario)

# plot with appropriate coordinate labels
modifyList(gsnap, list(gval=zpred)) |> pkern_plot(ppars=list(main='kriging prediction'))

#' Note that this computation happens quite fast compared to most kriging methods: The `pkern_cmean`
#' function call takes only around 1/6 of a second to complete on my desktop for the 1 million point
#' prediction problem.
#'
#' When conditional expectation will be computed repeatedly (such as with a long
#' time series sharing the same spatial covariance model), we can also pre-compute and recycle some
#' of the intermediate variance components to speed things up.
#'
#'
#' ## Variance components
#'
#' The code below demonstrates how to copy the pre-computed model objects, `pc` for later
#' use. On my computer, subsequent calls to `pkern_cmean` for different data vectors
#' only take about 20 milliseconds if we pass `pc` as an argument:

# pre-computed objects
pc = pkern_cmean(zinc.sg, gsnap, pars.vario, pc=TRUE)

#' The pre-computed model objects list `pc` contains grid and parameter info, along with various matrices
#' and factorizations derived from the full covariance matrix for the outer grid, which is usually too
#' large to work with directly.

# print the names of the covariance components
names(pc)

# they occupy very little space in memory
object.size(pc) |> print(units='Mb')

#' `pc` stores covariance and cross-covariance matrices for subsets of the grid, partitioned according
#' to their relationship with the subgrid ("vs", "vo", "vso", "vsc"). `pc` also includes indexing vectors
#' ("s", "o", "oj", "oi", "sobs") mapping the input point data to rows and columns of these matrices
#'
#'
#' ## Variance components example
#'
#' In the kriging problem one takes the covariance matrix for the input points, say `V`, and multiplies
#' its inverse by the input data vector. Rather than storing V or its inverse, `pkern` stores the
#' diagonalization of `V` (from `base::eigen`), absent the nugget effect:

# structure of the eigendecomposition for sampled data covariance V
str(pc$ed)

# notice there is one eigenvector per sampled point
meuse.sf |> nrow()

#' `V` can be recovered as a submatrix of the covariance matrix for the full subgrid (which
#' includes unsampled subgrid locations where we have NAs). This is the Kronecker product
#' of the component y and x covariance matrices, which are also stored in `pc`:

# structure of subgrid covariance components
str(pc$vs)

# notice their dimensions correspond to the subgrid dimensions
gsnap$sg$gdim |> print()

# compute the (nugget-free) sampled data covariance matrix as a submatrix of big Kronecker product
V = Matrix::kronecker(pc$vs$x, pc$vs$y)[pc$sobs, pc$sobs]

# verify the eigendecomposition works as expected
V.compare = pc$ed$vectors %*% diag(pc$ed$values) %*% t(pc$ed$vectors)
abs(V - V.compare) |> max() |> print()

#'
#' ## Where is the nugget effect?
#'
#' The variance components in `pc` define a model without a nugget. The nugget effect is
#' handled separately by specifying nugget variance in element "nug" of the kernel parameters
#' list `pars`.
#'
#' For example, the following code prints the nugget variance in the model fitted earlier
#' by `pkern_vario_fit` and computes the covariance matrix for the sampled points including this nugget
#' effect, by adding it to the diagonal matrix in the eigendecomposition.

# compute `V` for the full model with nugget
print(pars.vario$nug)
V.nugget = pc$ed$vectors %*% diag(pc$ed$values + pars.vario$nug) %*% t(pc$ed$vectors)

# notice this is quite different from the nugget-free matrix
abs(V - V.nugget ) |> max() |> print()

#' The reason for handling the nugget separately like this is that our full covariance matrix
#' is spatially separable only when the nugget effect is zero. Separability allows the use
#' of Kronecker products and various other computational shortcuts, so we omit the nugget
#' effect until it is needed - eg. by adding it diagonalization of `V` in functions like
#' `pkern_cmean` and `pkern_LL` (see Gilboa et al., 2015, for more on this)
#'
#' For more details on the objects in `pc`, see `?pkern_precompute`
#'
#'
#' ## Conditional variance
#'
#' The pointwise variance of the kriging predictions is computed by calling `pkern_cmean`
#' in a loop over a special set of inputs, one per input point. Expect this to
#' be slightly faster than running `pkern_cmean` n times, where n is the number of input points:

# compute pointwise variance
zv = pkern_variance(pc, quiet=TRUE)

# plot the result
modifyList(gsnap, list(gval=zv)) |> pkern_plot(ppars=list(main='kriging variance'))

#' This works by superimposing the n grid data vectors resulting from `pkern_cmean(vi)^2`, where
#' `vi` is the ith eigenvector of the sampled data covariance (scaled appropriately, and including
#' any nugget effect specified in `pc$pars`). If you're curious what these individual eigenvalue grids
#' look like, they can be viewed by using the `idx` argument of `pkern_variance`:

# pick an eigenvector index (1,2,..,50). They are ordered by decreasing eigenvalue magnitude
eigen.i = 45

# compute variance contribution related to this eigenvector
zv.i = pkern_variance(pc, quiet=TRUE, idx=eigen.i)

# plot the result
plot.title = paste('variance contribution from eigenvector', eigen.i)
modifyList(gsnap, list(gval=zv.i)) |> pkern_plot(ppars=list(main=plot.title))

#' The smallest eigenvalues (here, the 50th, 49th, etc) tend to have the most complex and
#' interesting appearance. I don't know of any analytical use for these plots, but they sure
#' look cool!
#'
#'
#' ## Back-transform
#'
#' To get predictions and variances on the original zinc concentration scale
#' we need to back-transform, then apply an adjustment for bias introduced by the
#' transform. This involves the pointwise variance, which we just computed
#' in the code above (object `zv`)
#'
#' In this case (logarithm transform) we can use the adjustment of Laurent (1963):

# back-transform and bias adjustments
zpred.adj = exp(zpred + logzinc.mean + (zv/2))
zv.adj =  exp( 2*(zpred + logzinc.mean) + zv ) * (exp(zv) - 1)

# plot kriging estimates in original coordinate system then plot
modifyList(gsnap, list(gval=zpred.adj)) |>
  pkern_plot(ppars=list(main='kriging predictions', leg='[zinc]'))

# repeat for pointwise variance
modifyList(gsnap, list(gval=zv.adj)) |>
  pkern_plot(ppars=list(main='kriging variance'))

#'
#' ## Summary
#'
#' This vignette is intended to get users up and running with `pkern` by demonstrating
#' a very simple example on a familiar dataset. In other vignettes we will look in more
#' detail at how `pkern_cmean` actually works, and how users can modify the workflow to
#' incorporate a trend model.


#'
#' ## Markdown
#'
#' This chunk below is used to create the markdown document you're reading from the
#' R script file (with same name, in this directory). It uses `rmarkdown` to create
#' the md file, then substitutes local image paths for github URLs so the document
#' will display the images properly on my repo page.

if(FALSE)
{
  # Restart session and run code chunk below to build the markdown file
  library(here)
  library(rmarkdown)

  # make the markdown document and delete the unwanted html
  path.input = here('vignettes/pkern_vignette.R')
  path.output = here('vignettes/pkern_vignette.md')
  path.garbage = here('vignettes/pkern_vignette.html')
  rmarkdown::render(path.input, clean=TRUE, output_file=path.output)
  unlink(path.garbage)

  # substitute local file paths for image files with URLs on github
  md.github = gsub('D:/pkern', 'https://github.com/deankoch/pkern/blob/main', readLines(path.output))
  writeLines(md.github, path.output)
}
