#' ---
#' title: "GHCND_data_preprocessing"
#' author: "Dean Koch"
#' date: "`r format(Sys.Date())`"
#' output: github_document
#' ---
#'
#' **Mitacs UYRW project**
#'
#' **pkern**: Fast kriging on gridded datasets
#'
#' This script prepares data for the "GHCND_data" shiny app
#'

#+ dependencies_hide, include=FALSE

# load pkern
library(devtools)
load_all()

#+ dependencies

# extra dependencies for loading and preparing spatial data
library(sf)
library(terra)
library(data.table)

# load helper functions and global variables
source('D:/pkern/vignettes/GHCN_vignette_helpers.R')

#'
#' ## Load and explore the weather data
#'
#' The weather station data are stored in a large CSV file. I start by opening this
#' table and adding some new attributes.
#'

#+ ghcn_open, include=FALSE

# open the CSV file and print attribute names
ghcn_dt = my_ghcn_dt()
names(ghcn_dt)

#'
#' Interpreting character strings as dates can be slow, so `my_ghcn_dt` computes
#' Julian day and year (integer attributes) for faster indexing. It also adds the
#' logarithm of elevation for use as a covariate later, and creates a new variable
#' 'LOG_PRCP', the logarithm of precipitation (in 10ths of a mm) plus 1.
#'
#' The function `my_ghcn_filter` is for drawing subsets of the dataset within a
#' specified data range. It returns one variable at a time, in the form of a matrix
#' with one row per date.
#'

# print instructions then extract data by variable
my_ghcn_filter(ghcn_dt)

# lookup table for printing titles
title_lookup = data.frame(vname=c('TMIN', 'TMAX', 'PRCP', 'LOG_PRCP'),
                          vname_print=c(paste('temperature', c('min', 'max'), '(C)'),
                                        rep('precipitation (mm)', 2)))

vname = 'PRCP'
data_mat = my_ghcn_filter(ghcn_dt, vname)
str(data_mat)

# Different dates form the "layers" here, which are stored in columns of a matrix
n_layers = ncol(data_mat)
n_station = nrow(data_mat)

# observed dates
times = colnames(data_mat) |> as.Date()
years = format(times, '%Y') |> as.integer()
jdates = format(times, '%j') |> as.integer()

# seasonal index
my_season = function(d, s=0) sin( 2 * pi * (d + s) / 365 )
seasonal_sin = jdates |> my_season()
seasonal_cos = jdates |> my_season(s=pi/2)

#'
#' Another helper function, `my_stations`, creates an `sf` POINT collection from the station
#' positions. Most of the stations are missing data on certain days (see below).
#'

# create the station points
stations_sf = my_stations(ghcn_dt)

# count the number of NAs at each station
n_NA = apply(data_mat, 1, \(x) sum(is.na(x)))
stations_sf[['n_missing']] = n_NA
is_complete = n_NA == 0




#'
#' `my_base_loader` loads some geographical reference points for plotting purposes
#'

# load reference data
geo_list = my_base_loader()
summary(geo_list)

#'
#' The next chunk plots the station locations and the number of missing dates (for precip)
#' at each one. The circled points are complete time series, with no missing data. There
#' are 35 such stations with complete records.
#'

# make a plot indicating completeness at each weather station
stations_sf['n_missing'] |> plot(pch=16, reset=FALSE, cex=0.5)
st_geometry(stations_sf)[is_complete] |> plot(add=TRUE, cex=2)
my_plot_annotations(geo_list=geo_list)

# copy the subset of complete stations
stations_complete_sf = stations_sf[is_complete,]

dev.off()

#'
#' The idea is to use these special stations to fit a simple covariance model by maximum
#' likelihood, then use the fitted parameters for interpolation from all stations.
#'
#'
#' ## Regression (simple)
#'
#' We initially de-trend the data using a linear regression on elevation and its logarithm.
#' This gives a good initial guess for the regression parameters, which are then improved
#' later on after we've estimated the covariance structure.
#'

# full covariates list
X_obs_all = stations_sf[c('elevation', 'log_elevation')] |>
  sf::st_drop_geometry() |>
  lapply(function(x) rep(x, n_layers)) |>
  c(list(year = rep(years, each=n_station),
         seasonal_sin = rep(seasonal_sin, each=n_station),
         seasonal_cos = rep(seasonal_cos, each=n_station)))


# subset with complete response data
y_obs = data_mat[is_complete,] |> as.vector()
X_obs = stations_complete_sf[c('elevation', 'log_elevation')] |>
  sf::st_drop_geometry() |>
  lapply(function(x) rep(x, n_layers)) |>
  c(list(year = rep(years, each=sum(is_complete)),
         seasonal_sin = rep(seasonal_sin, each=sum(is_complete)),
         seasonal_cos = rep(seasonal_cos, each=sum(is_complete))))

# simple linear regression
lm_elevation = lm(y~elevation+log_elevation+seasonal_sin+seasonal_cos, c(list(y=y_obs), X_obs))
summary(lm_elevation)

# matrices with fitted values and residuals for each layer in columns
pred_mat = predict(lm_elevation, newdata=X_obs_all) |> matrix(nrow=nrow(data_mat))
resid_mat = pred_mat - data_mat

# plot the fitted elevation relationship over the expected range in the UYRW area
elevation_test = seq(from=1e3, to=4e3, length.out=1e3)
test_X = data.frame(elevation=elevation_test, log_elevation=log(elevation_test), seasonal_sin=0, seasonal_cos=0)
y = predict(lm_elevation, newdata=test_X)
plot(y~elevation_test)

# plot the fitted seasonal relationship
jday = seq(365)
seasonal_sin_test = my_season(jday)
seasonal_cos_test = my_season(jday, s=pi/2)
test_X = data.frame(elevation=0, log_elevation=0, year=0, seasonal_sin=seasonal_sin_test, seasonal_cos=seasonal_cos_test)
y = predict(lm_elevation, newdata=test_X)
plot(y~jday)

# plot the mean of model residuals by date
# # center data at layer-wise means
# test_means = resid_mat |> colMeans(na.rm=TRUE)
# plot(y=test_means, x=as.Date(names(test_means)))
#
#


#'
#'
#' ## Covariance analysis
#'
#' `pkern` works with gridded data, so the first step in the analysis is to snap the
#' weather points to a grid. For the purposes of covariance estimation, this grid can be
#' fairly sparse. Later on we will increase the resolution.
#'

# # test on a single layer the old fashioned way
# idx = idx + 1
# stations_sf[['test']] = resid_mat[,idx]
# coarse_pkern2 = stations_sf['test'] |> pkern_snap(list(gres=gres_coarse))
# pkern_plot(coarse_pkern2)
#
# coarse_pkern2$gval = coarse_pkern2$gval - mean(coarse_pkern2$gval, na.rm=TRUE)
#
# fit_result = pkern_fit(coarse_pkern2, X=0)
# fit_result$pars |> pkern_plot_pars(g=coarse_pkern)
# fit_result$df
#
# summary( stations_sf[['test']] )
# sum( stations_sf[['test']] > 0, na.rm=TRUE)
# sum( !is.na(stations_sf[['test']]) )

# snap the complete weather points to a coarse grid
gres_coarse = c(y=500, x=500) # spacing in metres
coarse_pkern = stations_complete_sf['station_num'] |> pkern_snap(list(gres=gres_coarse))
is_obs = !is.na(coarse_pkern[['gval']])


#pkern_plot(coarse_pkern)


# this vector has a station key value at the indices of each observed point, and NAs otherwise
coarse_pkern[['idx_obs']] = coarse_pkern[['gval']]

# data matrices must be in column-vectorized order
idx_reorder = coarse_pkern[['gval']][is_obs]

# train on all layers
idx_train = seq(n_layers)
#idx_train = sample.int(n_layers, 1e3)
#idx_train = resid_mat[idx_reorder,] |> apply(2, \(x) sum(x>0) > 10) |> which() |> sample(size=2)
#
# #
# i = 2
# stations_complete_sf[['test']] = resid_mat[,idx_train[i]]
# coarse_pkern2 = stations_sf['test'] |> pkern_snap(list(gres=gres_coarse))
# pkern_plot(coarse_pkern2)


# supply multiple layers as columns in a matrix
coarse_pkern[['gval']] = resid_mat[idx_reorder, idx_train]

# center data at layer-wise means
coarse_means = coarse_pkern[['gval']] |> colMeans(na.rm=TRUE)
plot(y=coarse_means, x=as.Date(names(coarse_means)))
coarse_pkern[['gval']] = coarse_pkern[['gval']] |> sweep(2, coarse_means, '-')

# write list element idx_obs which tells pkern how to unpack gval
coarse_pkern[['idx_obs']][is_obs] = seq(sum(is_complete))

#
fit_result = pkern_fit(coarse_pkern, X=0, iso=F)



dev.off()
fit_result$pars |> pkern_plot_pars(g=coarse_pkern)
pars_interpolation



#
##
###
####


#'
#' ## Define the scale of interest
#'
#' Predictor values (when supplied) must be supplied for all grid points. Since we have a spatial
#' predictor that is already gridded - elevation - I will use a simplified version (an up-scaling)
#' of the DEM to establish the output grid.
#'
#' This output grid establishes the point locations where `pkern` will later compute expected
#' weather values by interpolation.
#'
#' The source DEM grid is very detailed, having around 100 million points. So I first
#' upscale by a factor of 25 to get a coarser grid with around 150 thousand points
#' (the 25 relates the square roots of these numbers).
#'

dem_src_rast = basemap_data_dir |> file.path('dem.tif') |> terra::rast()

# upscale the DEM grid for this demo. Later we will use it at full resolution
dem_agg_fac = 25
dem_rast = terra::aggregate(dem_src_rast, fact=dem_agg_fac)

# crop to points area
dem_rast = terra::crop(dem_rast, stations_sf)

#writeRaster(dem_rast, 'test.tif')
#
# plot the result
plot(dem_rast)

# convert to pkern object
dem_pkern = pkern_grid(dem_rast)

###

#
# testsf = st_crop(stations_sf, st_bbox(dem_rast))
# testsf$station_num |> length()
#
# plot(dem_rast, reset=FALSE)
# plot(testsf['station_num'], add=T)


# plot(stations_sf['station_num'], pch=1, cex=2)
# st_bbox(dem_rast) |> st_as_sfc() |> plot(add=T)
# plot(stations_sf['station_num'][stations_sf[['station_num']] == 230,], pch=16, cex=2, add=TRUE)


##
##

# from=stations_sf['station_num']
# g=dem_pkern
# crop_from=TRUE
# crop_g=FALSE
#

# snap all weather points to the grid
stations_pkern = stations_sf['station_num'] |> pkern_snap(dem_pkern, crop_from=TRUE)
pkern_plot(stations_pkern, reset=FALSE)
plot(stations_sf['station_num'], pch=1,  cex=2, add=T)


# plot(stations_sf['station_num'], pch=1, cex=2)
# st_bbox(dem_rast) |> st_as_sfc() |> plot(add=T)
#
# pkern_plot(stations_pkern, reset=FALSE)

# indexing on grid for weather stations
idx_map = stations_pkern[['gval']]
is_obs = !is.na(idx_map)
idx_reorder = idx_map[is_obs] # station numbers in grid vectorized order

# omit those falling outside of prediction region
#stations_cropped_sf = stations_sf[match(idx_reorder, stations_sf[['station_num']]),]
stations_cropped_sf = stations_sf[idx_reorder,]
pkern_plot(stations_pkern, reset=FALSE)
plot(stations_cropped_sf['station_num'], pch=1,  cex=2, add=T)

# compute constant-in-time predictors
x_dem = dem_pkern[['gval']]

# testing other parameters
pars_interpolation = fit_result$pars
# new_range = 3e4
# new_eps = 1
# pars_interpolation$psill = pars_interpolation$psill + pars_interpolation$eps - new_eps
# pars_interpolation$eps = new_eps
# pars_interpolation$x$kp = new_range
# pars_interpolation$y$kp = new_range

# # compute the variance
var_result = pkern_cmean(stations_pkern, pars_interpolation, X=0, out='v', quiet=F)
var_pkern = stations_pkern |> modifyList(list(gval=var_result))
pkern_plot(var_pkern)

# # find a region of acceptable variance
# pkern_plot(var_accept_pkern)
# is_acceptable = var_accept_pkern[['gval']] < 10^2
# var_accept_pkern[['gval']][!is_acceptable] = NA
# pkern_plot(var_accept_pkern)


#saveRDS(data_mat, 'test.rds')


# select an index (date) to start at
idx = 1960
idx = idx-1

my_pal = function(x) hcl.colors(x, 'Spectral', rev=T)
my_col = function(x) my_pal(1e2)[ num_to_cent(x) ]



bbox_dem = st_bbox(dem_rast) |> st_as_sfc()
n_dem = length(x_dem)

if(1)
{
  idx = idx + 1
  #dev.off()

  idx_date = as.Date(colnames(data_mat)[idx])
  yr = format(idx_date, '%Y') |> as.integer()
  jday = format(idx_date, '%j') |> as.integer()
  seasonal_sin_pred = my_season(jday)
  seasonal_cos_pred = my_season(jday, s=pi/2)

  # combine all predictors
  X_new = data.frame(elevation = x_dem,
                     log_elevation = log(x_dem),
                     year = rep(yr, n_dem),
                     seasonal_sin = rep(seasonal_sin_pred, n_dem),
                     seasonal_cos = rep(seasonal_cos_pred, n_dem))

  # compute their linear combination
  y_lm = predict(lm_elevation, newdata=X_new)


  ## log scale plots

  # plot the station data for this date as points
  #stations_sf[['obs_data']] = data_mat[,idx]
  #plot(stations_sf['obs_data'], pch=16, cex=2, reset=FALSE, pal=my_pal)

  #bbox_dem = st_bbox(dem_rast) |> st_as_sfc()
  #plot(bbox_dem, border='black', add=TRUE)

  #y_test = data_mat[,idx][idx_reorder]
  #stations_pkern[['gval']][is_obs] = y_test

  #dev.off()
  #pkern_plot(stations_pkern)


  # plot the residuals for the complete data regression
  #stations_sf[['resid_test']] = resid_mat[,idx]
  #plot(stations_sf['resid_test'], pch=16, cex=2, reset=FALSE)

  #bbox_dem = st_bbox(dem_rast) |> st_as_sfc()
  #plot(bbox_dem, border='black', add=TRUE)

  # compute all available linear model residuals
  y_residual = y_lm[is_obs] - data_mat[,idx][idx_reorder]
  stations_pkern[['gval']][is_obs] = y_residual
  #dev.off()
  #pkern_plot(stations_pkern, reset=FALSE)

  #plot(stations_sf['resid_test'], cex=2, add=TRUE)

  # pkern_plot(stations_pkern, reset=FALSE)
  # st_geometry(stations_sf)[is_complete] |> plot(add=TRUE, cex=2)
  # my_plot_annotations(geo_list=geo_list)

  # interpolate the anomaly
  fit_result_cmean = pkern_cmean(stations_pkern, pars_interpolation, X=0)
  zpred_pkern = stations_pkern |> modifyList(list(gval=fit_result_cmean))
  #pkern_plot(zpred_pkern)

  # kriging prediction
  krig = y_lm - fit_result_cmean

  # apply inverse transformations
  if(vname=='PRCP') krig[krig < 0] = 0
  if(vname=='LOG_PRCP') krig = (exp(krig + var_result/2) - 1)

  # change units from 10ths
  krig_plot = krig / 10


  # filter inapplicable pixels
  #krig_original[!is_acceptable] = NA

  # make a pkern object
  krig_plot_pkern = stations_pkern |> modifyList(list(gval=krig_plot))
  #pkern_plot(krig_plot_pkern, zlab='precip (mm)')

  # transform source point data
  krig_obs = data_mat[,idx]
  if(vname=='LOG_PRCP') krig_obs = ( exp(krig_obs) - 1 ) / 10
  krig_obs = krig_obs / 10

  # create points object from station observed
  stations_sf[['obs_data']] = krig_obs
  plot_sf = stations_sf['obs_data'][idx_reorder,]
  plot_sf = plot_sf[!is.na(plot_sf$obs_data),]
  zlim_min = min(c(plot_sf[['obs_data']], krig_original), na.rm=T)
  zlim_max = max(c(plot_sf[['obs_data']], krig_original), na.rm=T)
  zlim = c(zlim_min, zlim_max)

  vname_print = title_lookup[['vname_print']][title_lookup[['vname']] == vname]
  title_print = paste(as.character(colnames(data_mat)[idx]), 'daily', vname_print)
}


#zlim = range(plot_sf[['obs_data']], na.rm=T) + c(-5,5)
pkern_plot(krig_plot_pkern,
           main=title_print,
           zlab=vname,
           zlim=zlim,
           reset=FALSE)

mybreaks = seq(min(zlim), max(zlim), by=diff(zlim)/1e2)
#pkern_plot(krig_plot_pkern, zlab='precip (mm)', reset=FALSE)
plot(st_geometry(plot_sf), cex=2, add=TRUE)
plot(plot_sf, pch=16, cex=1, add=TRUE, pal=my_pal, breaks=mybreaks)

#my_plot_annotations(geo_list=geo_list, bw=TRUE)
plot(st_crop(geo_list[['watershed_boundary']], bbox_dem), col=NA, add=TRUE)
st_geometry(geo_list[['state_boundaries']]) |> st_crop(bbox_dem) |> plot(border='grey20', col=NA, add=TRUE)

