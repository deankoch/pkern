#' ---
#' title: "raster_vignette"
#' author: "Dean Koch"
#' date: "`r format(Sys.Date())`"
#' output: github_document
#' ---
#'
#' **Mitacs UYRW project**
#'
#' **pkern**: Fast kriging on gridded datasets
#'
#' This vignette shows how to downscale raster data using pkern
#'



#'
#' ## The mountain pine beetle
#'
#' This vignette looks at data on trees killed by mountain pine beetle in the Canadian province of
#' British Columbia (BC). If you're interested in learning more about pine beetles and their impact
#' check out the TRIA-Net program and their members' work. My research on this subject (see 1, 2, 3,
#' and my thesis) was completed with their support, and the `pkern` package itself is a revival
#' of code written during that project.
#'
#' The pine beetle is a native species in BC forests that feeds on and spends most of its life
#' inside mature pine trees. Its populations naturally undergo boom/bust cycles that can kill large
#' numbers of trees over a short period of time. Dead trees turn color and become visible in aerial
#' surveys as large clusters of "red tops" among otherwise healthy green crowns.
#'
#' BC conducts these surveys yearly and makes the data public. The period 2001-2018 of this historical
#' record can be downloaded using the `rasterbc` package. We will use it to look at the pine-rich
#' Merritt, BC, area in the year 2006, when beetle populations were rapidly climbing and outbreaks
#' widespread.
#'

#'
#' ## The rasterbc package
#'
#' This vignette uses data from the `rasterbc` package, along with the `sf` and `terra` libraries for
#' working with rasters and geometries.

#+ dependencies_show

# load pkern and geospatial libraries
library(sf)
library(terra)

library(rasterbc)

# pkern
library(devtools)
load_all()


#'
#' The first step with `rasterbc` is to specify a directory to put the GeoTIFFs it will download.
#' Change `rasterbc_storage_path` below to something appropriate for your local machine:
#'

# set storage directory for TIFF files
rasterbc_storage_path = 'D:/pkern/vignettes/rasterbc_example/data'
datadir_bc(rasterbc_storage_path, quiet=TRUE)

#'
#' BC is a very large province, so `rasterbc` downloads data in tiles covering small segments. You
#' can either specify the tile(s) you want, or supply a geo-referenced object (point or polygon) to
#' define your area of interest. See `rasterbc::ntspoly_bc` for a feature collection delineating all
#' available tiles.
#'

# identify an area of interest by its NTS code
SNRC_merritt = '092I'

# load a polygon to use as bounding box
bbox_merritt = ntspoly_bc[ntspoly_bc[['NTS_SNRC']] == SNRC_merritt,]

# display the polygon on a map of the BC area
title_bbox = paste('NTS tile number', bbox_merritt[['NTS_SNRC']], 'in south-central BC')
plot(st_geometry(ntspoly_bc), main=title_bbox)
plot(bbox_merritt, col='grey', add=TRUE)


#'
#' If you only have the place-name of the area of interest, I recommend using the `getbb()` function
#' from the `osmdata` package to look up coordinates using the Nominatim API. From there you can build
#' a polygon using `sf::st_bbox(...) |> sf::st_as_sfc(...)`.
#'
#' To list the available layer names, categories, and descriptions do `rasterbc::listdata_bc()`. This
#' vignette uses four layers: 'IBM_mid' from the 'fids' (Forest Insect and Disease Survey) collection;
#' along with 'vegTreed' 'needle' and 'pinusTotal' from the 'pine' collection.
#'
#' The first layer estimates the fraction of the land in each pixel that exhibits pine
#' beetle damage. This on its own is not a great estimator of the total number of killed pine, since
#' forest composition can be quite variable. For example areas with very little pine may be labelled
#' as 100% damaged, though the total number of trees killed is very small.
#'
#' So we scale by taking products with the last three layers, which together estimate the density of
#' pine within the surrounding vegetation. The result should be approximately proportional to the
#' number of pine beetle killed trees. We are missing a stems-per-unit-area factor but this will be
#' fine for demonstration purposes.
#'

# download and display midpoint of MPB damage (fraction) estimate in 2006
IBM_full_rast = rasterbc::opendata_bc(SNRC_merritt, 'fids', 'IBM_mid', 2006)

# download layers for pine forest density estimate (nearest available year)
treed_full_rast = rasterbc::opendata_bc(SNRC_merritt, 'pine', 'vegTreed', 2001)
needle_full_rast = rasterbc::opendata_bc(SNRC_merritt, 'pine', 'needle', 2001)
pinus_full_rast = rasterbc::opendata_bc(SNRC_merritt, 'pine', 'pinusTotal', 2001)

# estimate Pinus density as fraction of pixel area and plot
forest_full_rast = (treed_full_rast/100) * (needle_full_rast/100) * (pinus_full_rast/100)

# estimate pine beetle damage (fraction of area affected)
damage_full_rast = forest_full_rast * IBM_full_rast

#'
#' The `terra` package makes these kinds of raster-wide operations very easy. For most of this
#' vignette however we will use an unpacked (and reordered) version of the `terra::SpatRaster`
#' object, which you can get by passing the raster to `pkern_grid`. When you're done, you can
#' export a `pkern` list object back `terra::SpatRaster` using `pkern_export`.
#'

# create the list object
damage_full = pkern_grid(damage_full_rast)
str(damage_full)

# pkern has its own plotter with a slightly different layout
pkern_plot(damage_full)

#'
#' This heatmap shows the pattern of beetle-killed trees across the landscape in 2006. These
#' dead trees are closely linked to the state of the beetle population one year prior, so this
#' map also indicates the state of the pine beetle population in 2005, which was growing and
#' expanding. The area of this plot is about 120 X 150 km, so the extent of the damage is quite
#' large already.
#'
#'
#' find a complete sub-grid
#'
#' When your spatial dataset has many points, computation time and memory use
#' becomes a concern. `pkern` works most efficiently if there are many NAs (ie very few
#' observed points) or else none at all (ie the complete data case, all observed).
#' It will not be very efficient for the in-between situation, where large numbers of
#' `NA`s are mixed with large numbers of observed data. In that case we recommend the modeler
#' simplify their dataset by sub-sampling or splitting the problem into smaller pieces.
#'
#' In this example the source data have a large number of `NA`s around the boundary where a
#' mask was used to clip points lying outside of predefined area (a tile from Canada's
#' National Topographic System). The data are otherwise complete - inside the warped
#' rectangle there are no `NA`s - so we can simply crop the grid to this smaller extent.
#'
#' The function `pkern_sub` automates the process of trimming `NA` rows and columns from
#' the outer edge of the grid.
#'

# crop the damage layer to exclude outer NAs
damage = pkern_sub(damage_full)

# plot the damage raster and check for NAs
pkern_plot(damage)
anyNA(damage[['gval']])

# plot showing the selected sub-grid within outer grid
bbox_trim = st_as_sfc(st_bbox(pkern_coords(damage, out='sf', corner=TRUE)))
pkern_plot(damage_full, reset=FALSE)
plot(bbox_trim, col=adjustcolor('white', alpha.f=0.2), add=TRUE)

# another plot showing position in larger region
plot(st_geometry(ntspoly_bc), main=bbox_merritt[['NTS_SNRC']])
plot(bbox_merritt, col='grey', add=TRUE)
plot(bbox_trim, col='black', add=TRUE)


#'
#' As I understand it, finding a largest complete sub-grid in a grid with missing
#' data is an NP-hard problem. The algorithm implemented in `pkern_sub` is a heuristic solution
#' with no theoretical guarantees - it may fail to find a suitable result in your dataset,
#' and it will be slow when there are a lot of outer grid lines to trim. This is something
#' that may be improved upon in future releases of pkern.
#'
#' If you run into problems with automated cropping, you can always specify grid lines
#' to crop manually by number. In the code below, we first extract these numbers by calling
#' `pkern_sub` again with `idx=TRUE`, but you can use any selection of grid lines, provided
#' they result in a regular sub-grid.
#'

# find the row and column numbers to exclude
subgrid_ij = pkern_sub(damage_full, idx=TRUE)
str(subgrid_ij)

# passing this list back to `pkern_sub` produces the same output as before
damage_compare = pkern_sub(damage_full, subgrid_ij)
identical(damage, damage_compare)

#'
#' subgrid_ij has two elements, `keep` and `remove`, each a list of two integer vectors
#' `i` and `j`, the row and column numbers. Either (or both) of `keep` and `remove` can
#' be passed back to `pkern_sub` to crop any grid with same dimensions as `damage_full`.
#'
#'
#' ## Up-scale to simulate a coarser resolution dataset
#'
#' The terms up-scale and down-scale refer to changing the resolution of a gridded sampling
#' scheme by (respectively) omitting or introducing sample points. For example, down-scaling
#' by a factor of two means placing a new point midway between every two adjacent grid points,
#' producing a denser sample, at a finer resolution. Up-scaling does the inverse, dropping every
#' other grid point.
#'
#' In the next code chunk we will up-scale the pine beetle damage grid to get a coarser resolution
#' version to demonstrate interpolation. Most of the source dataset points are dropped in
#' this step, but we can use them as a test (or hold-out) set for gauging prediction error.
#'

# up-scale by a factor of 25
damage_coarse = pkern_rescale(damage, up=25)
pkern_plot(damage_coarse)


#' ## A warning about resolution
#'
#' The term "resolution" can, unfortunately, refer to both the spacing distance between
#' neighbouring sample points `s = c(dx, dy)` (pixel width, or sampling interval), and its
#' reciprocal, `1/s` (pixel density, or detail level).
#'
#' When reporting or measuring resolution for spatial data we typically use to pixel width.
#' This is how `pkern` keeps track of the resolution. For example, the source data in this
#' vignette are at "100m resolution", which you can verify by printing `damage$gres`
#'

# 100m resolution source
print(damage[['gres']])

# upscaled version at 2.5km resolution
print(damage_coarse[['gres']])

#'
#' However, when speaking about resolution with qualifiers like "high" and "low" - or
#' more generally in the context of image processing - we are usually referring
#' to pixel density. For example, a 10m grid has "higher" resolution than a 100m grid
#' because its points are spaced more densely, so the detail level is higher. To avoid
#' confusion I will try to use the unambiguous terms "fine" and "coarse" instead.
#'
#'
#' ## Fit a covariance model
#'
#' Our task below is to use the 2.5km resolution data in `damage_coarse` to estimate the
#' pine beetle damage level at 100m resolution. We already have the original survey data at
#' that fine level, so we can later compare it with our prediction to get a sense of error.
#'
#' A grid at 100m resolution will have about `25^2` or 625 times more points than
#' `damage_coarse`. We need to fill in these unobserved points with data interpolated from
#' the coarse dataset.
#'
#' To do this, we first fit a model that relates the separation distance between points to
#' their similarity - or more precisely, their covariance. This will allow us to compute an
#' expected value for all unknown points, called the "kriging predictor". To "krig" refers to
#' this entire workflow - model fitting then expected value to get an interpolated dataset.
#'
#' `pkern` krigs by means of `base::optim`, which is used to automatically find covariance
#' parameters that maximize the likelihood function for `damage_coarse`. Parameter bounds
#' are set automatically, based on the grid dimensions and the variance of the observed data.
#' All you need to do is pass the grid to `pkern_fit`. Readers interested in the details of
#' model fitting should check out the help pages for `pkern_LL`, `pkern_nLL`, and `pkern_optim`.
#'

# fit the model using `base::optim` on likelihood function
fit_result_iso = pkern_fit(damage_coarse)

# print fitted parameters
fit_result_iso[['df']]

#'
#' This data frame printed above shows the lower and upper bounds, along with fitted values,
#' for each of three covariance parameters. `eps` is the nugget effect, representing
#' (among other things) measurement error in the source data; `psill` is the partial sill,
#' which in combination with `eps` determines the variance; and the others, `y.rho` and
#' `x.rho` are range parameters, determining a distance along each axis at which
#' correlation begins to drop off significantly.
#'
#' Parameter bounds are set automatically by default, but can be set manually as needed.
#' Users should check this data frame to be sure the bounds are reasonable and that the
#' optimizer hasn't converged on a parameter bound (indicating an issue with the model fit
#' or bounds).
#'
#' The model fit can be visualized by plotting the covariance footprint of a single pixel,
#' using the function `pkern_plot_pars`
#'

# make a copy of the fitted parameters list and plot the estimated covariance footprint
fit_pars_iso = fit_result_iso[['pars']]
pkern_plot_pars(fit_pars_iso, damage_coarse)


# TODO: change pkern_fit to accept bds
# repeat with eps fixed around 0
# bds = pkern_bds(fit_pars_iso, damage_coarse)
# bds['eps',] = 1e-3 + ( 1e-16 * c(-1, 0, 1) )
#
# pars_fixed = pkern_pars_make()
# pars_fixed$eps = 1e-3
# fit_result_iso = pkern_fit(damage_coarse, pars_fixed)
# fit_result_iso$df









#'
#' The plot shows the degree of influence that any given observed data-point has on the
#' conditional mean of its neighbouring points. It hints at the overall "shape" we should
#' expect in the final interpolation product, as this will be built by super-imposing many
#' (scaled) copies of this same footprint image.
#'
#'
#' ## The semi-variogram diagnostic
#'
#' The semivariogram is another diagnostic for model fit. It shows the relationship between
#' covariance and distance, both theoretical (eg. as specified by `fit_pars_iso`) and as
#' sampled directly from the data. Use `pkern_sample_vg` to compute the latter, and
#' `pkern_plot_semi` to display either or both on the same graph.
#'

# compute and plot the semi-variogram along with the fitted (theoretical) curve
damage_vg = pkern_sample_vg(damage_coarse, n_pp=1e5, n_bin=100)
pkern_plot_semi(damage_vg, fit_pars_iso)

#'
#' Note that this plot will look different each time the vignette code is run, as there
#' is some sub-sampling happening in each call (as controlled by argument `n_pp`).
#'
#' Sample semi-variances are aggregated by distance and plotted as grey points. The
#' theoretical values, generated from the parameters `fit_pars_iso` are plotted as a blue
#' curve. Its y-intercept is the nugget effect `eps`; its long-distance limit is the sum
#' of `eps` and `psill` (the partial sill); and the point at which it levels off is the
#' range specified by `y.rho` and `x.rho`.
#'
#' Geo-statistics theory tells us the three features mentioned above (nugget, sill and range)
#' should roughly echo the pattern in the grey point cloud of sample values, provided we've
#' estimated them properly. In fact, traditionally, the theoretical (blue) curve is fitted
#' directly to the point cloud, either by eye or by some version of least squares. Thanks to
#' the computationally efficiency of `pkern` we were able to use the more elegant technique of
#' maximum likelihood; the semi-variogram shown here is only a diagnostic.
#'
#' Sample semi-variograms are fussy. There are many different ways of sub-sampling,
#' binning distances, and aggregating results, and different approaches can result in
#' dramatically different looking point clouds. This on its own is a good reason to prefer
#' the more consistent maximum likelihood approach implemented in `pkern_fit` over
#' semi-variograms based fitting. However, users are encouraged to understand the
#' semi-variogram and its relation to the model, and to explore the `pkern_sample_vg`
#' arguments `n_pp`, `n_bin`, `idx` and `fun` and their effects.
#'
#'
#' ## Anisotropy and kernel selection
#'
#' By default, `pkern_fit` fits an isotropic Gaussian covariance model. Isotropic means the
#' covariance footprint is constrained to be radially symmetric (ie circle shaped). This
#' constraint can be relaxed somewhat with `iso=FALSE`, which allows some stretching of
#' the covariance footprint along the x and y axes.
#'

# fit the model again allowing anisotropy
fit_result_aniso = pkern_fit(damage_coarse, iso=FALSE)
fit_result_aniso[['df']]

# make a copy of the fitted parameters list and plot the estimated covariance footprint
fit_pars_aniso = fit_result_aniso[['pars']]
pkern_plot_pars(fit_pars_aniso, damage_coarse)

# plot the semi-variogram
pkern_plot_semi(damage_vg, fit_pars_aniso)

#'
#' The difference in this case is very small, with the two range parameters only differing
#' by about 20 metres. For illustration, let's make the covariance model highly anisotropic
#' and see the effect on the semi-variogram and covariance footprint
#'

# double the range parameter for y but not x and plot result
fit_pars_demo = fit_pars_aniso
fit_pars_demo[['y']][['kp']] = 2 * fit_pars_demo[['y']][['kp']]
pkern_plot_pars(fit_pars_demo, damage_coarse)
pkern_plot_semi(damage_vg, fit_pars_demo)

#'
#' Notice that the fitted model in the semi-variogram graph is now a ribbon plot (instead
#' of a curve). This is because an anisotropic model can admit a range of correlations for
#' a given distance, depending the direction. `pkern_plot_semi` displays the full range
#' for any given distance.
#'
#' One can also specify alternative 1-d component models, to replace the Gaussian. These are
#' unconventional models for 2-dimensional problems, but they can provide quite a bit more
#' flexibility in terms of shape (see my paper).
#'
#' For example, a Matern X Matern model produces the following result

# fit the model again with a product-Matern kernel
fit_result = pkern_fit(damage_coarse, pars='mat', iso=FALSE)
fit_result[['df']]
fit_pars = fit_result[['pars']]

#'
#' Note that whenever a non-Gaussian (ie non-default) 1-d component is specified, the
#' resulting 2-d covariance model will be at least slightly anisotropic, even if `iso=TRUE`.
#' The Gaussian is the only separable covariance model in 2-d that is also isotropic.
#'
#' Notice that two additional parameters `y.kap` and `x.kap` are present in the results.
#' These are shape parameters which determine curvature in the semi-variogram.
#'

# plot covariance footprint
pkern_plot_pars(fit_pars, damage_coarse)

# plot the semi-variogram
pkern_plot_semi(damage_vg, fit_pars)


#'
#'
#' ## Interpolation
#'
#' Given a set of fitted covariance parameters, such as `fit_pars`, `pkern` makes it easy
#' to create a prediction at finer resolution. It does this by computing the conditional mean
#' given the observed data using an efficient implementation of the Kriging equation
#' `pkern_cmean`, designed specifically for grids.
#'
#' The observed data should be embedded in the desired output grid, with all non-observed
#' points set to NA. This can be done with `pkern_rescale`, which we used earlier make our
#' artificially coarse training data.
#'
#'

fit_pars = fit_result[['pars']]
fit_pars = fit_pars_aniso

#fit_pars$eps = 1e-16

g_observed = pkern_rescale(damage_coarse, down=25)
fit_cmean = pkern_cmean(g_obs=g_observed, pars=fit_pars)

# truncate to original range 0,1
fit_cmean[fit_cmean < 0] = 0
fit_cmean[fit_cmean > 1] = 1

modifyList(g_observed, list(gval=fit_cmean)) |> pkern_plot()

# find the subset of points that made it into the final interpolated product
ij = Map(function(obs, orig) which(orig %in% obs), obs=g_observed$gyx, orig=damage$gyx )
is_in_pred = pkern_sub_idx(damage$gdim, setNames(ij, c('i', 'j')))
is_unseen = is.na(g_observed$gval)

# compute error with pkern
err_pkern = fit_cmean[is_unseen] - damage[['gval']][is_in_pred][is_unseen]

# dataframe version for use with gstat
damage_coarse_df = as.data.frame(cbind(damage=damage_coarse$gval, pkern_coords(damage_coarse)))


#'
#' ## Comparison with Inverse Distance Weighting
#'

library(gstat)

mg = gstat(id='damage', formula = damage~1, locations=~x+y, data=damage_coarse_df, nmax=100, set=list(idp = 4))
r = pkern_export(g_observed)
z = interpolate(r, mg, debug.level=-1, index=1)
z_g = pkern_grid(z)
pkern_plot(z_g)

# compute error with idw
err_idw = z_g[['gval']][is_in_pred][is_unseen] - damage[['gval']][is_in_pred][is_unseen]

# mean squared error is nearly cut in half
mean(err_idw^2, na.rm=TRUE)
mean(err_pkern^2)


#'
#' ## Comparison with ordinary kriging in gstat
#'

v = variogram(damage~1, ~x+y, data=damage_coarse_df)
mv = fit.variogram(v, vgm(1, 'Gau', 6e3, 1.6e-3))
gOK = gstat(NULL, 'damage', damage~1, damage_coarse_df, locations=~x+y, model=mv)

# warning: this call took almost two hours to complete
#OK = interpolate(r, gOK, debug.level=-1, computeVar=FALSE)

OK_g = pkern_grid(OK$damage.pred)
pkern_plot(OK_g)

err_OK = OK_g[['gval']][is_in_pred][is_unseen] - damage[['gval']][is_in_pred][is_unseen]

# similar improvement with pkern, and way faster
mean(err_OK^2, na.rm=TRUE)
mean(err_pkern^2)



#'
#' ## Comparison with Voronoi polygons
#'
#'


pkern_plot(damage_coarse)

damage_coarse_sf = pkern_coords(damage_coarse, out='sf')
damage_coarse_terra = as(damage_coarse_sf, 'SpatVector')
#bbox_terra = as(bbox_trim, 'SpatVector')


damage_voronoi_terra = voronoi(damage_coarse_terra)

idx_relate = relate(damage_voronoi_terra, damage_coarse_terra, "contains") |> apply(1, which)
poly_vals = values(damage_coarse_terra)[['gval']][idx_relate]
damage_voronoi_terra[] = poly_vals

damage_voronoi_rast = rasterize(damage_voronoi_terra, damage_full_rast, field='gval')

damage_voronoi_full = pkern_grid(damage_voronoi_rast)
damage_voronoi = pkern_sub(damage_voronoi_full, subgrid_ij)
pkern_plot(damage_voronoi)


# find the subset of points that made it into the final interpolated product
ij = Map(function(obs, orig) which(orig %in% obs), obs=g_observed$gyx, orig=damage$gyx )
is_in_pred = pkern_sub_idx(damage$gdim, setNames(ij, c('i', 'j')))
is_unseen = is.na(g_observed$gval)

# compute error with voronoi
err_voronoi = damage_voronoi[['gval']][is_in_pred][is_unseen] - damage[['gval']][is_in_pred][is_unseen]

# compute error with pkern
err_pkern = fit_cmean[is_unseen] - damage[['gval']][is_in_pred][is_unseen]

# mean squared error is nearly cut in half!
mean(err_voronoi^2)
mean(err_pkern^2)

# examine only the errors at non-zero points
is_nz = damage[['gval']][is_in_pred][is_unseen] > 0
par(mfrow=c(1,2))
hist(err_voronoi[is_nz], breaks=100, ylim=c(0, 65e3))
hist(err_pkern[is_nz], breaks=100, ylim=c(0, 65e3))

# similar result
mean(err_voronoi[is_nz]^2)
mean(err_pkern[is_nz]^2)


if(0)
{
#'
#' kriging variance
#'
#'

kriging_var = pkern_cmean(g_obs=g_observed, pars=fit_pars, out='v')


# find mapping from points in `g_damage` (low resolution) to original resolution grid

g_damage
damage





#pkern_sub_idx()


# initialize a grid of NAs at original resolution and copy the up-scaled data
# g_observed = damage
# g_observed[['gval']][] = NA
# g_observed[['gval']][] = NA


g_observed = pkern_rescale(damage, up=10)
pkern_plot(g_observed)

fit_cmean = pkern_cmean(g_obs=g_observed, pars=fit_pars)

# truncate to original range 0,1
fit_cmean[fit_cmean < 0] = 0
fit_cmean[fit_cmean > 1] = 1

modifyList(g_observed, list(gval=fit_cmean)) |> pkern_plot()



pkern_plot(g_damage)

pkern_plot(damage)



#' we could also use `pkern_rescale(g_damage, down=25)` here, but the result might
#' aligned slightly differently from the original `damage` grid
#'


g_down = pkern_rescale(g_damage, down=10)
fit_cmean = pkern_cmean(g_obs=g_down, pars=fit_pars)

# truncate to original range 0,1
fit_cmean[fit_cmean < 0] = 0
fit_cmean[fit_cmean > 1] = 1

pkern_plot(g_damage)
modifyList(g_down, list(gval=fit_cmean)) |> pkern_plot()
pkern_plot(damage)




# # get matching foresta and elevation layers
# g_forest = pkern_rescale(forest, up=10)
# pkern_plot(g_forest)

#'
#'

# fit simple linear trend model
lm_forest = lm(y~x, list(y=g_damage[['gval']], x=g_forest[['gval']]))

# plot fitted values
fval_forest = c(cbind(1, g_forest[['gval']]) %*% lm_forest[['coefficients']])
pkern_plot(modifyList(g_damage, list(gval = fval_forest)))

# copy and plot residuals
res_forest = g_damage[['gval']] - fval_forest
g_res = pkern_grid(modifyList(g_damage, list(gval = res_forest)))


pkern_plot(g_res)


#'
#' fit a covariance model
#'
#'
fit_result = pkern_fit(g_damage, X=NA)
#fit_result = pkern_fit(g_res, X=0)
fit_pars = fit_result$pars


pkern_plot_pars(fit_pars, g_damage)

#'
#' interpolate residuals at higher resolution
#'

g_down = pkern_rescale(g_damage, down=10)
fit_cmean = pkern_cmean(g_obs=g_down, pars=fit_pars)

# truncate to original range 0,1
fit_cmean[fit_cmean < 0] = 0
fit_cmean[fit_cmean > 1] = 1

pkern_plot(g_damage)
modifyList(g_down, list(gval=fit_cmean)) |> pkern_plot()
pkern_plot(damage)

sum(damage$gval)

sum(g_damage$gval) * 100
sum(fit_cmean)



#'
#' add trend to get predictions
#'

fval_down = c(cbind(1, forest[['gval']]) %*% lm_forest[['coefficients']])
ij_sub = Map(\(xy, xy_sub) which(xy %in% xy_sub), forest$gyx,  g_down$gyx) |> setNames(c('i', 'j'))


# predictions
pred_down = fval_down[ pkern_sub_idx(forest$gdim, ij=ij_sub) ] + fit_cmean

# truncate to original range 0,1
pred_down[pred_down < 0] = 0
pred_down[pred_down > 1] = 1


g_pred = modifyList(g_down, list(gval=pred_down))
g_pred |> pkern_plot()

# compare with original and source points
pkern_plot(damage)
pkern_plot(g_damage)

}

#'
#' ## Markdown
#'
#' This chunk below is used to create the markdown document you're reading from the
#' R script file with same name, in this directory. It uses `rmarkdown` to create
#' the md file, then substitutes local image paths for github URLs so the document
#' will display the images properly on my repo page.

if(FALSE)
{
  # Restart session and run code chunk below to build the markdown file
  library(here)
  library(rmarkdown)

  # make the markdown document and delete the unwanted html
  path.input = here('vignettes/rasterbc_example/rasterbc_vignette.R')
  path.output = here('vignettes/rasterbc_example/rasterbc_vignette.md')
  path.garbage = here('vignettes/rasterbc_example/rasterbc_vignette.html')
  rmarkdown::render(path.input, clean=TRUE, output_file=path.output)
  unlink(path.garbage)

  # substitute local file paths for image files with URLs on github
  md.github = gsub('D:/pkern', 'https://github.com/deankoch/pkern/blob/main', readLines(path.output))
  writeLines(md.github, path.output)
}
