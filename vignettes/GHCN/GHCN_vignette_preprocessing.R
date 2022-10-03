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
source('D:/pkern/vignettes/GHCN/GHCN_vignette_helpers.R')

# path to the output file(s) written by this script
dest_rds_file = 'D:/pkern/vignettes/GHCN/GHCN_preprocessing_data.rds'


#'
#' `my_base_loader` loads some geographical reference points for plotting purposes
#'

# load reference data
geo_list = my_base_loader()
summary(geo_list)


#'
#' ## Load the weather data
#'
#' The weather station data are stored in a large CSV file. I start by opening this
#' table and adding some new attributes.
#'

#+ ghcn_open, include=FALSE

# open the CSV file and print attribute names
ghcn_dt = my_ghcn_dt()
names(ghcn_dt)

# print instructions for later
my_ghcn_filter(ghcn_dt)

#'
#' Interpreting character strings as dates can be slow, so `my_ghcn_dt` has computed
#' Julian day and year (integer attributes) for faster indexing. It also created a
#' new weather variable 'LOG_PRCP', the logarithm of precipitation (in 10ths of a mm)
#' plus a small constant (`log_const`) defined in the helper function file.
#'
#' Another helper function, `my_stations`, creates an `sf` POINT collection from the station
#' positions. Most of the stations are missing data on certain days (see below).
#'

# create the station points and count them
station_sf = my_stations(ghcn_dt)
n_station = nrow(station_sf)
cat(paste('\n', n_station, 'weather stations in loaded GHCN dataset\n'))


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
#' upscale by a factor of `dem_agg_fac=25` to get a coarser grid with around 150 thousand points
#' (`dem_agg_fac` relates the square roots of these numbers).
#'

# upscale the DEM grid for this demo. Later we will use it at full resolution
dem_src_rast = basemap_data_dir |> file.path('dem.tif') |> terra::rast()
dem_agg_fac = 100
dem_rast = terra::aggregate(dem_src_rast, fact=dem_agg_fac)

# crop to points area
dem_rast = terra::crop(dem_rast, station_sf)

# convert DEM to pkern grid object and plot
dem_pkern = pkern_grid(dem_rast)
pkern_plot(dem_pkern, main='DEM and weather data locations', zlab='elevation (m)')
n_grid = prod(dem_pkern[['gdim']])
cat(paste0('\noutput grid is ', paste(dem_pkern[['gdim']], collapse=' x '), ' (', n_grid, ' points)\n'))

# add GHCN points and their reported elevations
plot(station_sf['elevation'], pch=16, add=TRUE, pal=my_pal, breaks=my_breaks(dem_pkern[['gval']]))
plot(st_geometry(station_sf), add=TRUE)


#'
#' ## Snap weather points to the grid
#'
#' The snapping function `pkern_snap` will create a copy of the DEM with values set to
#' NA everywhere, except at grid cells having an overlying weather station, where the
#' associated station number is copied to the grid.
#'

# snap all weather points to the grid
station_pkern = station_sf['station_num'] |> pkern_snap(dem_pkern, crop_from=TRUE)

#'
#' Notice some "duplicate" stations were dropped. This happens when stations are too close
#' together and wind up snapped to the same grid point. Such ties are broken by keeping
#' the closest of the stations (to the target grid point), and omitting the others.
#'
#' If you want to keep ALL stations you will need to find a snapping scheme that prevents
#' multiple assignments to the same grid point, even when it is the closest.
#' This can be done using the Hungarian algorithm, but it is a bit too complicated to
#' implement here (or in pkern), so we will just accept a bit of information loss.
#'
#' The grid data values returned by `pkern_snap` define a mapping (`idx_grid`) that can
#' be used to snap any data associated with the input points (say data vector `z`) onto the
#' grid by simply doing `z[idx_grid]`. This reorders `z` (into grid column vectorized order)
#' and introduces NAs wherever there is a grid point with no mapped station.
#'
#' For example, this next chunk shows how to do this with the reported elevations in the
#' station dataset:
#'

# copy the mapping vector
idx_grid = station_pkern[['gval']]

# vector of elevation point data
station_sf[['elevation']] |> str()
station_sf[['elevation']] |> summary()

# square bracket indexing to create a new grid vector with elevations snapped
station_elevation_gval = station_sf[['elevation']][idx_grid]

# after snapping (notice all the NAs)
station_elevation_gval |> str()
station_elevation_gval |> summary()

# plot as grid object
station_pkern |> modifyList(list(gval=station_elevation_gval)) |>
  pkern_plot(main='snapped station elevation values', zlab='elevation (m)', reset=TRUE)


#'
#' ## Trend modeling
#'
#' The interpolation technique used here works by splitting the data into two components:
#' (1) a predictable trend plus (2) random deviations about that trend. `pkern` handles the
#' deviations; Our job is to model the trend.
#'
#' This example shows how to build a linear model for the trend using least squares (OLS and GLS).
#' This means we will assume the trend is a weighted sum of predictors, with the weights (but
#' not the predictors) being the same at every point. This simplifies the problem to that of
#' estimating a small number of weights, AKA regression coefficients or beta values. Once this
#' is done we can easily compute the predicted (or fitted) trend, and move on to interpolation.
#'
#' There are two obvious components to the trend in the weather data: seasonal (eg. winter vs
#' summer) and topographical (eg. mountain peaks get colder, wetter weather). We will represent
#' them with the following linear predictors:
#'
#'  * elevation and its square root
#'  * sine and cosine of Julian date (with 1 and 1/2 year periods)
#'
#' For convenience I have defined a pair of helper functions for building data matrices
#' for these predictors. The topographical stays constant in time, so it takes the argument
#' `n_layer` to duplicate rows when there are multiple layers (ie dates) in the analysis
#'
#' The second data matrix will have rows for every date in the analysis, so it accepts a `dates`
#' argument (either a Julian date integer vector or Date class vector). It also accepts an `n_points`
#' argument, which duplicates rows when there are multiple observations on each date.
#'
#' To demonstrate, open the temperature minimum data
#'

# load the data for the selected variable
vname = 'TMIN'
data_mat_src = my_ghcn_filter(ghcn_dt, vname)
times = as.Date(colnames(data_mat_src))

# count the dates loaded
n_layer = ncol(data_mat_src)
cat(paste0('\nloaded ', n_layer, ' dates ', '(', paste(range(times), collapse=' to '), ')\n'))

# compute the daily mean values
daily_mean = data_mat_src |> apply(2, \(x) mean(x, na.rm=TRUE))

#'
#' `my_ghcn_filter` returns a big matrix containing the station data for weather variable `vname`.
#' Different dates - which I'm calling "layers" - are stored in columns of a matrix, whereas rows map
#' to stations in `station_sf` (same order).
#'


# count the number of NAs at each station and order by completeness
is_na_mat = is.na(data_mat_src)
order_na = order(rowSums(is_na_mat))

# find dates at which the n most complete stations all have data
for(i in seq(n_station)[-1])
{
  idx = order_na[i]
  idx_prev = order_na[i-1]
  is_na_mat[idx,] = is_na_mat[idx,] | is_na_mat[idx_prev,]
}

# count the number of dates
n_date = rowSums(!is_na_mat)[order_na]

# given a desired number of dates, find a maximal set of stations having complete data
n_date_min = 1e3
n_station_train = which(n_date < n_date_min)[1] - 1

# index and count of dates and stations selected
n_date_train = n_date[n_station_train]
idx_date_train = which(!is_na_mat[order_na[n_station_train],])
station_num_train = order_na[seq(n_station_train)] |> sort()

# omit stations cropped in snapping procedure
station_num_train = station_num_train[station_num_train %in% idx_grid]
n_station_train = length(station_num_train)

# recompute the indexing vector in case stations were dropped
idx_grid_train = rep(NA, n_grid)
idx_grid_train[ match(station_num_train, idx_grid) ] = seq(n_station_train)

# console output
msg_station = paste(n_station_train, 'of', n_station, 'stations selected,')
msg_date = paste('having complete data on', n_date_train, 'of', length(times), 'dates')
cat(paste('\n', msg_station, msg_date))

# indicate selected station points in the earlier plot
plot(st_geometry(station_sf)[station_num_train], add=TRUE)



#'
#'
#' ## Model fitting
#'

# take subset of data and centre at daily mean
data_mat_train = data_mat_src[station_num_train, idx_date_train]
data_mat_train_c = sweep(data_mat_train, 2, daily_mean[idx_date_train], '-')
train_pkern = station_pkern |> modifyList( list(gval=data_mat_train_c, idx_grid=idx_grid_train) )

# build temporally constant predictors matrix
elevation_X = my_elevation_X(station_sf[['elevation']][station_num_train]) |> as.matrix()

# fit the model
fit_result_initial = pkern_optim(g_obs=train_pkern, X=elevation_X, iso=F)
pars_interpolate_initial = fit_result[['pars']]
print(fit_result_initial[['df']])
pars_interpolate_initial |> pkern_plot_pars(station_pkern)

# run GLS to recompute residuals
elevation_betas = pkern_GLS(train_pkern, pars_interpolate_initial, X=elevation_X)
elevation_lm_z = pkern_GLS(train_pkern, pars_interpolate_initial, X=elevation_X, out='z')

idx = 4
idx=idx+1

res_mat = sweep(data_mat_train_c, 1, elevation_lm_z, '-')
res_means = res_mat |> apply(2, mean)
res_mat_c = sweep(res_mat, 2, res_means, '-')

train_pkern2 = train_pkern |> modifyList( list(gval=res_mat_c) )
fit_result2 = pkern_optim(g_obs=train_pkern2, X=0, iso=F)







res_pkern = station_pkern |> modifyList( list(gval=as.vector(res_mat[,idx])[idx_grid_train]) )
res_pkern |> pkern_plot()


range(res_pkern$gval, na.rm=TRUE)
#pars_interpolate_initial$eps = 1
# pars_interpolate$y$kp = 2e4
# pars_interpolate$x$kp = 2e4
xx = pkern_cmean(res_pkern, pars_interpolate_initial, X=0)
xx[!is.na(res_pkern$gval)] = res_pkern$gval[!is.na(res_pkern$gval)]
gg2 = res_pkern |> modifyList( list(gval=xx) )
pkern_plot(gg2)


elevation_lm_z + daily_mean[idx_date_train]

gg = station_pkern |> modifyList( list(gval=as.vector(resid_mat_train[,idx])[idx_grid_train]) )







pkern_plot(gg)

range(gg$gval, na.rm=TRUE)
pars_interpolate$eps = 1
# pars_interpolate$y$kp = 2e4
# pars_interpolate$x$kp = 2e4
xx = pkern_cmean(gg, pars_interpolate, X=0)
xx[!is.na(gg$gval)] = gg$gval[!is.na(gg$gval)]
gg2 = gg |> modifyList( list(gval=xx) )
pkern_plot(gg2)







#'
#' We can take this matrix and stack its columns to get a single vector of length `n_layer * n_station`
#' using R's `c` function. The corresponding covariates matrix can be constructed by passing `n_layer`
#' to `my_elevation_X` and `n_station` to `season_X` and column-binding the result. From there it is
#' easy to fit a simple OLS model:
#'

# use helper functions to build full covariates matrix
elevation_X = my_elevation_X(station_sf[['elevation']], n_layer)
season_X = my_season_X(times, n_station)
all_X = cbind(elevation_X, season_X)
str(all_X)

# fit simple linear model
lm_all = lm(w~., cbind(w=c(data_mat_src), all_X))
summary(lm_all)

# copy coefficients
intercept_ols = lm_all[['coefficients']][1L]
elevation_ols_betas = lm_all[['coefficients']][ names(elevation_X) ]
season_ols_betas = lm_all[['coefficients']][ names(season_X) ]

#'
#' I've hidden some tedious code for plotting marginal effects below. This produces the following
#' plot, showing (a random sample of) the model residuals by date, along with marginal residuals
#' after separately accounting for elevation and seasonality.
#'

#+ hide=TRUE

# contribution to fitted values from elevation and seasonality
elevation_fitted = as.matrix(my_elevation_X(station_sf[['elevation']])) %*% elevation_ols_betas
season_fitted = as.matrix(my_season_X(times)) %*% season_ols_betas

# separate elevation effect from fitted values
elevation_residuals_df = data.frame(z=c(sweep(data_mat_src, 2, season_fitted, '-')), elevation=all_X[['elevation']])
elevation_effect_df = data.frame(fitted=elevation_fitted + intercept_ols, elevation=station_sf[['elevation']])

# separate seasonal effect from fitted values
season_residuals_df = data.frame(z=c(sweep(data_mat_src, 1, elevation_fitted, '-')), date=rep(times, each=n_station))
season_effect_df = data.frame(fitted=season_fitted + intercept_ols, date=times)

# fitted and residuals for complete model by date
all_fitted = intercept_ols + rep(season_fitted, each=n_station) + rep(elevation_fitted, n_layer)
all_residuals = data_mat_src - all_fitted
all_residuals_df = data.frame(z=c(all_residuals), date=rep(times, each=n_station))

# random sample to plot
idx_plot = sample.int(nrow(all_X), 1e4)

# plotting details
col_point = adjustcolor('black', alpha.f=0.05)
col_line = adjustcolor('blue', alpha.f=0.3)
residuals_title = 'OLS residuals by date'
seasonal_title = 'seasonal effect'
elevation_title = 'elevation effect'

# show all three plots in the same window
par(mfrow=c(1,3))

  # plot all residuals (divide by 10 to get from 10ths of a degree to degrees celsius)
  plot(c(z/10)~date, data=all_residuals_df[idx_plot,], ylab=vname, main=residuals_title, pch=16, col=col_point)

  # seasonality
  plot(c(z/10)~date, data=season_residuals_df[idx_plot,], ylab=vname, main=seasonal_title, pch=16, col=col_point)
  lines(c(fitted/10)~date, data=season_effect_df, lwd=3, col=col_line)
  lines(c(fitted/10)~date, data=season_effect_df, lwd=1, col=col_line)

  # elevation
  plot(c(z/10)~elevation, data=elevation_residuals_df[idx_plot,], ylab=vname, main=elevation_title, pch=16, col=col_point)
  lines(c(fitted/10)~elevation, data=elevation_effect_df[order(elevation_effect_df[['elevation']]),], lwd=3, col=col_line)
  lines(c(fitted/10)~elevation, data=elevation_effect_df[order(elevation_effect_df[['elevation']]),], lwd=1, col=col_line)

par(mfrow=c(1,1))

#'
#' The story of these effects is clear: temperatures are lower on the average at higher elevations, and there is a
#' strong periodicity aligning with seasons of the year.
#'
#' The plan is....
#'
#'
#' ## Data preparation
#'
#' find a complete subset
#'

# count the number of NAs at each station and order by completeness
is_na_mat = is.na(data_mat_src)
order_na = order(rowSums(is_na_mat))

# find dates at which the n most complete stations all have data
for(i in seq(n_station)[-1])
{
  idx = order_na[i]
  idx_prev = order_na[i-1]
  is_na_mat[idx,] = is_na_mat[idx,] | is_na_mat[idx_prev,]
}

# count the number of dates
n_date = rowSums(!is_na_mat)[order_na]

# given a desired number of dates, find a maximal set of stations having complete data
n_date_min = 1e3
n_station_train = which(n_date < n_date_min)[1] - 1

# index and count of dates and stations selected
n_date_train = n_date[n_station_train]
idx_date_train = which(!is_na_mat[order_na[n_station_train],])
station_num_train = order_na[seq(n_station_train)] |> sort()

# omit stations cropped in snapping procedure
station_num_train = station_num_train[station_num_train %in% idx_grid]
n_station_train = length(station_num_train)

# console output
msg_station = paste(n_station_train, 'of', n_station, 'stations selected,')
msg_date = paste('having complete data on', n_date_train, 'of', length(times), 'dates')
cat(paste('\n', msg_station, msg_date))



#'
#'
#' ## Model fitting
#'

# copy residuals for these stations and dates
resid_mat_train = all_residuals[station_num_train, idx_date_train]
idx_grid_train = rep(NA, n_grid)
idx_grid_train[ match(station_num_train, idx_grid) ] = seq(n_station_train)
train_pkern = station_pkern |> modifyList( list(gval=resid_mat_train, idx_grid=idx_grid_train) )

# exclude dates with unusual mean or variance
# daily_var = data_mat_src[,idx_date_train] |> apply(2, \(x) var(x, na.rm=TRUE))
# daily_mean = data_mat_src[,idx_date_train] |> apply(2, \(x) mean(x, na.rm=TRUE))
# is_var_typical = ( daily_var > quantile(daily_var, 0.25) ) & (daily_var < quantile(daily_var, 0.75) )
# is_mean_typical = ( daily_mean > quantile(daily_mean, 0.25) ) & (daily_mean < quantile(daily_mean, 0.75) )
#
# idx_sub_train = is_var_typical & is_mean_typical
# sum(idx_sub_train)
# train_sub_pkern = train_pkern |> modifyList( list(gval=resid_mat_train[,idx_sub_train]) )


# exclude outlier dates
#idx_sub_train = ( daily_mean > quantile(daily_mean, 0.05) ) & (daily_mean < quantile(daily_mean, 0.95) )
#idx_sub_train = ( daily_mean > quantile(daily_mean, 0.75) ) #& (daily_mean < quantile(daily_mean, 0.75) )
#idx_sub_train = (daily_mean < quantile(daily_mean, 0.75) )

resid_mat_train = data_mat_src[station_num_train, idx_date_train]
daily_resid_mean = resid_mat_train |> apply(2, \(x) mean(x, na.rm=TRUE))
train_sub_pkern = train_pkern |> modifyList( list(gval=sweep(resid_mat_train, 2, daily_resid_mean, '-')) )


plot(times[idx_date_train], daily_resid_mean)
lines(times[idx_date_train], daily_resid_mean)


# fit the model
X_obs = as.matrix(elevation_X[station_num_train,])
fit_result = pkern_optim(g_obs=train_sub_pkern, X=X_obs, iso=F)
#fit_result = pkern_optim(g_obs=train_sub_pkern, X=0, iso=F)
fit_result$df
fit_result$pars |> pkern_plot_pars(station_pkern)
pars_interpolate = fit_result$pars


idx = 4
idx=idx+1
gg = station_pkern |> modifyList( list(gval=as.vector(resid_mat_train[,idx])[idx_grid_train]) )
pkern_plot(gg)

range(gg$gval, na.rm=TRUE)
pars_interpolate$eps = 1
# pars_interpolate$y$kp = 2e4
# pars_interpolate$x$kp = 2e4
xx = pkern_cmean(gg, pars_interpolate, X=0)
xx[!is.na(gg$gval)] = gg$gval[!is.na(gg$gval)]
gg2 = gg |> modifyList( list(gval=xx) )
pkern_plot(gg2)


#

#
#



# run GLS to get coefficients and fitted values
pars_interpolation_first = pars_interpolation_second = fit_result$pars
elevation_betas = pkern_GLS(train_pkern, pars_interpolation_first, X=X_obs)
elevation_lm_z = pkern_GLS(train_pkern, pars_interpolation_first, X=X_obs, out='z')

# remove the GLS fitted trend and plot residuals
resid_mat_daily = resid_mat_train - elevation_lm_z
plot(times[idx_date_train], apply(resid_mat_daily, 2, mean))
abline(h=mean(elevation_lm_z))

# prefer dates that deviate from trend
abs_mean_daily = resid_mat_daily |> apply(2, mean) |> abs()
idx_train = quantile(abs_mean_daily, 0.5) > abs_mean_daily

# fit again with GLS betas fixed
train_sub_pkern = train_pkern |> modifyList(list(gval=resid_mat_daily[, idx_train]))
fit_result2 = pkern_optim(g_obs=train_sub_pkern, X=0, iso=F)
fit_result2$df
fit_result2$pars |> pkern_plot_pars(station_pkern)
pars_interpolation_second = fit_result2$pars

pars_interpolation_first |> pkern_pars_update()
pars_interpolation_second |> pkern_pars_update()



#'
#'
#'
#'
#'
#'
#'
#' This next chunk goes through each weather variable and estimates a covariance
#' function for the entire time series
#'
#'

# storage for outputs
stor_template = vector(mode='list', length=length(vnames)) |> setNames(vnames)
data_mat_all = lm_seasonal_all = lm_elevation_all = pars_interpolation_all = pars_interpolation = var_pkern_all = stor_template

vname = 'TMIN'
vname = 'TMAX'
vname = 'LOG_PRCP'
vname = 'PRCP'

# extract and analyse each variable in a loop
for(vname in vnames)
{
  # load the data for the selected variable
  data_mat_src = my_ghcn_filter(ghcn_dt, vname)

  # Different dates form the "layers" here, which are stored in columns of a matrix. Rows are stations
  times = colnames(data_mat_src) |> as.Date()
  n_layer = ncol(data_mat_src)
  n_station = nrow(data_mat_src)

  # build a temporal predictors matrix
  X_seasonal = data.frame(sin_j = my_season(format(times, '%j')),
                          cos_j = my_season(format(times, '%j'), s=pi/2),
                          sin_j2 = my_season(format(times, '%j'), f=2),
                          cos_j2 = my_season(format(times, '%j'), s=pi/2, f=2),
                          year = as.integer(format(times, '%Y')))

  # estimate temporal trend by OLS
  X_seasonal_all = X_seasonal[rep(seq_along(times), each=n_station),]
  lm_seasonal = lm(y~., data=data.frame(y=c(data_mat_src), X_seasonal_all))
  fitted_seasonal = predict(lm_seasonal, newdata=X_seasonal)

  # copy the daily means, plot them against OLS prediction line
  data_means = data_mat_src |> colMeans(na.rm=TRUE)
  plot(times, data_means/10, ylab=paste(vname, '(daily means)'), col=adjustcolor('black', alpha.f=0.2), pch=16)
  lines(times, fitted_seasonal/10, lwd=2)

  # count the number of NAs at each station and order by completeness
  is_na_mat = is.na(data_mat_src)
  order_na = order(rowSums(is_na_mat))

  # find dates at which the n most complete stations all have data
  for(i in seq(n_station)[-1])
  {
    idx = order_na[i]
    idx_prev = order_na[i-1]
    is_na_mat[idx,] = is_na_mat[idx,] | is_na_mat[idx_prev,]
  }

  # count the number of dates
  n_date = rowSums(!is_na_mat)[order_na]

  # given a desired number of dates, find a maximal set of stations having complete data
  n_date_min = 1e3
  n_station_train = which(n_date < n_date_min)[1] - 1

  # index and count of dates and stations selected
  n_date_train = n_date[n_station_train]
  idx_date_train = which(!is_na_mat[order_na[n_station_train],])
  station_num_train = order_na[seq(n_station_train)] |> sort()

  # snap the complete weather points to a coarse grid, possibly dropping some duplicate stations
  coarse_gdim = round(dem_pkern[['gdim']] / 4)
  coarse_pkern = station_sf[station_num_train, 'station_num'] |> pkern_snap(g=list(gdim=coarse_gdim))

  # mapped station_num, possibly a subset of candidate set (if there was cropping)
  idx_grid_coarse = coarse_pkern[['gval']]
  station_num_train = idx_grid_coarse[!is.na(idx_grid_coarse)] |> sort()
  n_station_train = length(station_num_train)

  # console output
  msg_station = paste(n_station_train, 'of', n_station, 'stations selected,')
  msg_date = paste('having complete data on', n_date_train, 'of', length(times), 'dates')
  cat(paste('\n', msg_station, msg_date))

  # copy the training subset and sweep out seasonal means
  data_mat_train = data_mat_src[station_num_train, idx_date_train]
  resid_mat_train = sweep(data_mat_train, 2, fitted_seasonal[idx_date_train], '-')

  # compute the daily residual averages
  daily_resid_mean = resid_mat_train |> apply(2, mean)

  # copy to sparse pkern object
  idx_grid_train = idx_grid_coarse
  idx_grid_train[ match(station_num_train, idx_grid_coarse) ] = seq(n_station_train)
  train_pkern = coarse_pkern |> modifyList( list(gval=resid_mat_train, idx_grid=idx_grid_train) )

  # copy covariates matrix for these stations
  X_obs = elevation_X[station_num_train,]

  #train_idx = sample.int(n_date_train, 1e2)

  # train_group_idx = daily_means |> cut(quantile(daily_means, seq(0, 1, length.out=25)), include.lowest=TRUE)
  # train_group = levels(train_group_idx)
  # test_out = list()
  # pb = utils::txtProgressBar(max=length(train_group), style=3)
  # for(i in seq_along(train_group))
  # {
  #   train_idx = train_group_idx == train_group[i]
  #   train_sub_pkern = train_pkern |> modifyList(list(gval=resid_mat_train[, train_idx]))
  #   fit_result = pkern_fit(g_obs=train_sub_pkern, pars=pars, X=X_obs, iso=F, quiet=T)
  #   test_out[[i]] = cbind(i=i, interval=train_group[i], data.frame(t(pkern_pars_update(fit_result$pars))))
  #   utils::setTxtProgressBar(pb, i)
  # }
  # close(pb)
  #
  # test_df = do.call(rbind, test_out)
  # plot(eps~i, data=test_df)
  # plot(psill~i, data=test_df)
  # plot(y.rho~i, data=test_df)
  # plot(x.rho~i, data=test_df)

  #train_idx = ( quantile(daily_vars, 1/3) < daily_vars ) & ( daily_vars < quantile(daily_vars, 2/3) )
  #train_idx = ( quantile(daily_means, 1/3) < daily_means ) & ( daily_means < quantile(daily_means, 2/3) )

  #train_idx = ( quantile(abs*daily_means, 1/3) < daily_means ) & ( daily_means < quantile(daily_means, 2/3) )

  #train_idx = abs(daily_resid_mean) < quantile(abs(daily_resid_mean), 0.55)

  #train_idx = seq(n_date_train)
  #sum(train_idx)

  # n_train = 1e2
  # train_idx = sample.int(n_date_train, n_train) |> sort()

  train_idx = seq(n_date_train)
  train_sub_pkern = train_pkern |> modifyList(list(gval=resid_mat_train[, train_idx]))

  # fit the model
  fit_result = pkern_optim(g_obs=train_sub_pkern, X=X_obs, iso=F)
  fit_result$df
  fit_result$pars |> pkern_plot_pars(coarse_pkern)

  # run GLS to get coefficients and fitted values
  pars_interpolation_first = pars_interpolation_second = fit_result$pars
  elevation_betas = pkern_GLS(train_pkern, pars_interpolation_first, X=X_obs)
  elevation_lm_z = pkern_GLS(train_pkern, pars_interpolation_first, X=X_obs, out='z')

  # remove the GLS fitted trend and plot residuals
  resid_mat_daily = resid_mat_train - elevation_lm_z
  plot(times[idx_date_train], apply(resid_mat_daily, 2, mean))
  abline(h=mean(elevation_lm_z))

  # # pre-compute variance factorization and evaluate likelihoods in a loop
  # fac = pkern_var(train_pkern, pars_interpolation, scaled=TRUE, method='chol')
  # obj_ref = sapply(seq(n_date_train), \(i) {
  #
  #   # copy single day layer and compute reference objective
  #   gval = as.vector(resid_mat_daily[,i])[train_pkern$idx_grid]
  #   g_i = modifyList(train_pkern, list(gval=gval, idx_grid=NULL))
  #   -pkern_LL(pars=pars_interpolation, g_obs=g_i, fac=fac)
  #
  # })


  # fit again with GLS betas fixed
  # idx_train = obj_ref > quantile(obj_ref, 0.95)

  # prefer dates that deviate from trend
  abs_mean_daily = resid_mat_daily |> apply(2, mean) |> abs()
  #ar_daily = resid_mat_daily |> apply(2, var)

  idx_train = quantile(abs_mean_daily, 0.5) > abs_mean_daily

  # fit again with GLS betas fixed
  train_sub_pkern = train_pkern |> modifyList(list(gval=resid_mat_daily[, idx_train]))
  fit_result2 = pkern_optim(g_obs=train_sub_pkern, X=0, iso=F)
  fit_result2$df
  fit_result2$pars |> pkern_plot_pars(coarse_pkern)
  pars_interpolation_second = fit_result2$pars

  pars_interpolation_first |> pkern_pars_update()
  pars_interpolation_second |> pkern_pars_update()

  # # fit dates individually in a loop to determine if optim is successful
  # tol = 0.05
  # pb = txtProgressBar(1, n_train, style=3)
  # obj_ref = sapply(seq_along(train_idx), \(i) {
  #
  #   # copy single day layer
  #   setTxtProgressBar(pb, i)
  #   idx = train_idx[i]
  #   gval = as.vector(resid_mat_daily[,idx])[train_pkern$idx_grid]
  #   g_i = modifyList(train_pkern, list(gval=gval, idx_grid=NULL))
  #   fit_result_i = pkern_optim(g_obs=g_i, pars=pars, X=0, iso=F, quiet=TRUE)
  #
  #   # check for convergence to an upper or lower parameter bound (mark as invalid)
  #   fit = fit_result_i[['df']][['fitted']]
  #   fit_width = fit_result_i[['df']]['upper'] - fit_result_i[['df']]['lower']
  #   fit_pad = ( fit_result_i[['df']][c('lower', 'upper')] - fit ) / cbind(-fit_width, fit_width)
  #   is_far_enough = all(fit_pad[['lower']] > tol) & all(fit_pad[['upper']] < (1-tol))
  #
  #   # return TRUE if valid fit
  #   is_far_enough & all(fit_pad > 0)
  # })
  # close(pb)

  # sum(obj_ref)
  #
  #
  # # train again on subset for which individual day could be fitted
  # train_sub_pkern = train_pkern |> modifyList(list(gval=resid_mat_train[, train_idx[obj_ref]]))
  #
  # # fit model
  # pars = pkern_pars_make('gau')
  # fit_result = pkern_optim(g_obs=train_sub_pkern, pars=pars, X=0, iso=F)
  # #fit_result = pkern_fit(g_obs=train_sub_pkern, pars=pars, X=X_obs, iso=F)
  # fit_result$df
  # fit_result$pars |> pkern_plot_pars(coarse_pkern)
  # pars_interpolation = fit_result$pars
  #
  # pars_interpolation_early |> pkern_pars_update()
  # pars_interpolation |> pkern_pars_update()
  #
  #



  # dev.off()
  # eee = seq(1e3, 4e3, length.out=1e3)
  # xxx = matrix(c(eee, sqrt(eee)), ncol=2)
  # yyy = betas[1] + (xxx %*% betas[-1])
  # plot(eee, yyy)

  cat('\n**********************************\n\n')
  cat(paste('computing kriging variance for: ', vname))

  # # compute the variance
  var_result = pkern_cmean(station_pkern, pars_interpolation_second, X=NA, out='v', quiet=F)
  var_pkern = station_pkern |> modifyList(list(gval=var_result))
  pkern_plot(var_pkern)

  # copy results to storage
  data_mat_all[[vname]] = data_mat_src
  lm_elevation_all[[vname]] = elevation_betas
  lm_seasonal_all[[vname]] = lm_seasonal[['coefficients']]
  pars_interpolation[[vname]] = pars_interpolation_second
  pars_interpolation_all[[vname]] = pars_interpolation_first
  var_pkern_all[[vname]] = var_pkern

  #
  #   #
  #   # fit another model to outlier dates
  #
  #
    # # pre-compute variance factorization and evaluate likelihoods in a loop
    # fac = pkern_var(train_pkern, pars_interpolation, scaled=TRUE, method='chol')
    # obj_ref = sapply(seq(n_date_train), \(i) {
    #
    #   # copy single day layer and compute reference objective
    #   gval = as.vector(resid_mat_daily[,i])[train_pkern$idx_grid]
    #   g_i = modifyList(train_pkern, list(gval=gval, idx_grid=NULL))
    #   -pkern_LL(pars=pars_interpolation, g_obs=g_i, fac=fac)
    #
    # })
  #
  #   # select the worst 5%
  #   is_outlier = quantile(obj_ref, 0.95) < obj_ref
  #   train_sub_pkern = train_pkern |> modifyList(list(gval=resid_mat_daily[, is_outlier]))
  #
  #   # re-fit model on these points
  #   pars = pkern_pars_make('gau')
  #   outlier_fit_result = pkern_optim(g_obs=train_sub_pkern, pars=pars, X=0, iso=F)
  #
  #   outlier_fit_result$df
  #
  #   fit_result$df
  #
  #
  #
  #   ### EXTRA
  #
  #   coarser_pkern = pkern_rescale(coarse_pkern, up=2)
  #
  #   # remove the GLS fitted trend
  #   resid_mat_daily = resid_mat_train - elevation_lm_z
  #   plot(times[idx_date_train], apply(resid_mat_daily, 2, mean))
  #
  #   # fit daily layers individually in a loop
  #   tol = 0.05
  #   aniso_max = 4
  #   pars = pkern_pars_make('gau')
  #   n_test = 8e2
  #   fit_result_list = vector(mode='list', length=n_test)
  #   pb = txtProgressBar(1, n_test, style=3)
  #   idx_test = seq(n_test)
  #   for(i in idx_test)
  #   {
  #     # copy single day layer and compute reference objective
  #     setTxtProgressBar(pb, i)
  #     g_i = modifyList(train_pkern, list(gval=resid_mat_daily[,i]))
  #     obj_ref = -pkern_LL(pars=pars_interpolation, g_obs=g_i)
  #
  #     # write NULL to storage entry i when a model fit throws an error
  #     fit_result_i = tryCatch({
  #
  #       # model fit on ith daily residuals
  #       optim_result = pkern_optim(g_i, pars=pars, X=0, iso=F, quiet=TRUE)
  #
  #       # run again with isotropy constraint if fitted anisotropy exceeded bound
  #       xy_rho = optim_result[['df']][c('y.rho', 'x.rho'),][['fitted']]
  #       is_aniso_valid = exp( diff( log(sort(xy_rho)) ) ) < aniso_max
  #       if( !is_aniso_valid )
  #       {
  #         optim_result = pkern_optim(g_i, pars=pars, X=0, iso=T, quiet=TRUE)
  #         optim_result[['df']]['x.rho',] = optim_result[['df']]['y.rho',]
  #       }
  #
  #       # tryCatch doesn't want return() here
  #       optim_result
  #
  #     }, error = function(e) list())
  #
  #     # obj = objective value local minimum, the best estimated negative log likelihood
  #     obj = fit_result_i[['obj']]
  #     if( is.null(obj) ) { fit_result_list[[i]] = NULL } else {
  #
  #       # check for convergence to an upper or lower parameter bound (mark as invalid)
  #       fit = fit_result_i[['df']][['fitted']]
  #       fit_width = fit_result_i[['df']]['upper'] - fit_result_i[['df']]['lower']
  #       fit_pad = ( fit_result_i[['df']][c('lower', 'upper')] - fit ) / cbind(-fit_width, fit_width)
  #       is_far_enough = all(fit_pad[['lower']] > tol) & all(fit_pad[['upper']] < (1-tol))
  #       is_valid = is_far_enough & all(fit_pad > 0)
  #
  #       # copy results
  #       time_i = times[idx_date_train][i]
  #       fit_result_list[[i]] = data.frame(i, time_i, obj_ref, obj, obj_ref > obj, t(is_valid), t(fit)) |>
  #         setNames(c('i', 'time', 'obj_ref', 'obj', 'is_better', 'valid_fit', 'eps', 'psill', 'y.rho', 'x.rho'))
  #
  #       title_i = NULL
  #       if(!is_valid) title_i = 'INVALID'
  #       fit_result_i[['pars']] |> pkern_plot_pars(coarser_pkern, main=title_i)
  #     }
  #   }
  #   close(pb)
  #
  #   df_result = do.call(rbind, fit_result_list)
  #   df_result['obs_var'] = resid_mat_daily[,idx_test] |> apply(2, var)
  #   df_result['obs_med'] = resid_mat_daily[,idx_test] |> apply(2, median)
  #   df_result['obs_mean'] = resid_mat_daily[,idx_test] |> apply(2, mean)
  #   #df_result['obs_kurtosis'] = resid_mat_daily[,idx_test] |> apply(2, kurtosis)
  #   df_result['major_range'] = pmax(df_result[['y.rho']], df_result[['x.rho']])
  #   df_result['aniso_factor'] = df_result[,c('y.rho', 'x.rho')] |> apply(1, \(yx) exp( diff( log(sort(yx)) ) ))
  #
  #   df_result |> head()
  #   is_valid = df_result[['valid_fit']] & df_result[['is_better']]
  #   df_pruned = df_result[is_valid, ]
  #
  #   # notice bias correction in variance (sample variance is usually an underestimate)
  #   plot(c(eps+psill)~obs_var, data=df_pruned)
  #   abline(0,1)
  #
  #   # bi-modal distribution in estimated range
  #   hist(c(df_pruned$x.rho + df_pruned$y.rho), breaks=50)
  #   abline(v=pars_interpolation$y$kp + pars_interpolation$x$kp, col='red')
  #
  #   # there is some temporal autocorrelation in anisotropy and range
  #   plot(aniso_factor~time, data=df_pruned, pch=NA)
  #   lines(aniso_factor~time, data=df_pruned, col=adjustcolor('black', alpha.f=0.4))
  #   lines(c(0.6 + major_range/max(major_range))~time, data=df_pruned, col=adjustcolor('blue', alpha.f=0.4))
  #
  #
  #   # identify points with unusually low likelihood
  #   plot(major_range ~ obj_ref, data=df_pruned)
  #
  #   # indicate scale of kurtosis
  #   #points(major_range ~ (obj_ref), data=df_pruned, cex=df_pruned[['obs_kurtosis']]/max(df_pruned[['obs_kurtosis']]), pch=16, col='violet')
  #
  #   # indicate scale of variance
  #   points(major_range ~ (obj_ref), data=df_pruned, cex=df_pruned[['obs_var']]/max(df_pruned[['obs_var']]), pch=16, col='red')
  #
  #   # indicate scale of median
  #   points(major_range ~ (obj_ref), data=df_pruned, cex=abs(df_pruned[['obs_med']])/max(abs(df_pruned[['obs_med']])), pch=16, col='blue')
  #
  #   plot(major_range ~ (obj_ref), data=df_pruned)
  #   idx_show = df_pruned$obj_ref > quantile(df_pruned$obj_ref, 0.95)
  #
  #   points(major_range ~ (obj_ref), data=df_pruned[idx_show,], pch=16)
  #   abline(h=median(df_pruned$major_range[idx_show]), col='blue')
  #   abline(h=max(pars_interpolation$y$kp, pars_interpolation$x$kp))
  #
  #   train_idx = which(is_valid)[idx_show]
  #   train_sub_pkern = train_pkern |> modifyList(list(gval=resid_mat_daily[, train_idx]))
  #
  #   # fit model
  #   pars = pkern_pars_make('gau')
  #   fit_result2 = pkern_optim(g_obs=train_sub_pkern, pars=pars, X=0, iso=F)
  #   #fit_result = pkern_fit(g_obs=train_sub_pkern, pars=pars, X=X_obs, iso=F)
  #
  #   abline(h=max(fit_result2$df$fitted[3:4]))
  #   fit_result2$pars
  #
  #   plot(c(x.rho + y.rho)~c(eps+psill), data=df_pruned)
  #   plot(c(x.rho + y.rho)~obs_var, data=df_pruned)
  #   abline(h=pars_interpolation$y$kp + pars_interpolation$x$kp)
  #
  #   #plot(c(x.rho + y.rho)~obs_kurtosis, data=df_pruned)
  #
  #
  #
  #   df_pruned$test = pmax(df_pruned$y.rho, df_pruned$x.rho)
  #   plot(test~obs_med, data=df_pruned)
  #
  #
  #
  #   plot(c(x.rho + y.rho)~obs_med, data=df_pruned)
  #   plot(c(x.rho + y.rho)~obs_mean, data=df_pruned)
  #
  #
  #
  #   plot(eps~time, data=df_pruned)
  #   plot(psill~time, data=df_pruned)
  #   plot(y.rho~time, data=df_pruned)
  #   #lines(eps~time, data=df_pruned)
  #   plot(y.rho~time, data=df_pruned)
  #   plot(x.rho~y.rho, data=df_pruned)
  #
  #
  #   plot(x.rho~c(obj-obj_ref), data=df_pruned)
}






bbox_dem = st_bbox(dem_rast) |> st_as_sfc()
uyrw_boundary = geo_list[['watershed_boundary']] |> st_crop(bbox_dem)
state_boundaries = st_geometry(geo_list[['state_boundaries']]) |> st_crop(bbox_dem)

# mask for pixels in UYRW
uyrw_mask_pkern = uyrw_boundary |> as('SpatVector') |> rasterize(dem_rast) |> pkern_grid()

# compute elevation range observed in UYRW
uyrw_elevation_range = dem_pkern[['gval']][ !is.na(uyrw_mask_pkern[['gval']]) ] |> range()



# inputs
stor = list(

  misc=list(pars_interpolation=pars_interpolation,
            pars_interpolation_all=pars_interpolation_all,
            title_lookup=title_lookup,
            my_season=my_season,
            lm_elevation=lm_elevation_all,
            lm_seasonal=lm_seasonal_all,
            idx_grid=idx_grid,
            uyrw_boundary=uyrw_boundary,
            bbox_dem=bbox_dem,
            station_sf=station_sf,
            state_boundaries=state_boundaries),

  pkern = list(dem_pkern=dem_pkern,
               station_pkern=station_pkern,
               var_pkern=var_pkern_all,
               data_mat=data_mat_all)

  )

saveRDS(stor, dest_rds_file)



## EXTRA

#
# snum_range = idx_grid |> range(na.rm=TRUE)
# my_breaks = seq(min(snum_range), max(snum_range), by=10)
#
# # verify this looks correct
# station_pkern |> pkern_plot(reset=FALSE, pal='Spectral')
# station_sf['station_num'] |> plot(pal=my_pal, add=TRUE, breaks=my_breaks)
#
# # plot the dem and overlay with station elevation points
# # unfilled points are stations lying outside the observed elevation range in the UYRW
# #my_breaks = seq(min(uyrw_elevation_range), max(uyrw_elevation_range), by=10)
# pkern_plot(dem_pkern, reset=FALSE, zlab='elevation (m)', main='GHCN stations')
# station_sf['elevation'] |> plot(pch=16, pal=my_pal, add=TRUE, breaks=my_breaks)
# station_sf |> st_geometry() |> plot(add=TRUE)
#
# # # to keep the model focused on URYW we will exclude these outlying elevations from training
# # stations_elev = station_sf[['elevation']]
# # station_sf['elevation_match'] = !(stations_elev < uyrw_elevation_range[1]) & !(uyrw_elevation_range[2] < stations_elev)
# # if(any(!station_sf[['elevation_match']])) cat(paste('excluding', sum(!station_sf[['elevation_match']]), 'sites'))
#
