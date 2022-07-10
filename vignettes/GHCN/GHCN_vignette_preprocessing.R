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

# variable names and ranges expected (ignore snow for now)
vnames = c('TMIN', 'TMAX', 'PRCP', 'LOG_PRCP')
reasonable_ranges = list(TMIN=c(-6e2, 4e2), TMAX=c(-6e2, 6e2), PRCP=c(0, 16e3))
reasonable_ranges[['LOG_PRCP']] = log(1 + reasonable_ranges[['PRCP']])

# path to the output file(s) written by this script
dest_rds_file = 'D:/pkern/vignettes/GHCN/GHCN_preprocessing_data.rds'


#'
#' `my_base_loader` loads some geographical reference points for plotting purposes
#'

# load reference data
geo_list = my_base_loader()
summary(geo_list)

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
#' Another helper function, `my_stations`, creates an `sf` POINT collection from the station
#' positions. Most of the stations are missing data on certain days (see below).
#'

# create the station points
stations_sf = my_stations(ghcn_dt)

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
dem_agg_fac = 50
dem_rast = terra::aggregate(dem_src_rast, fact=dem_agg_fac)

# crop to points area
dem_rast = terra::crop(dem_rast, stations_sf)

#writeRaster(dem_rast, 'test.tif')
#
# plot the result
plot(dem_rast)

# for later plots
bbox_dem = st_bbox(dem_rast) |> st_as_sfc()
uyrw_boundary = geo_list[['watershed_boundary']] |> st_crop(bbox_dem)
state_boundaries = st_geometry(geo_list[['state_boundaries']]) |> st_crop(bbox_dem)

# convert to pkern object
dem_pkern = pkern_grid(dem_rast)

# snap all weather points to the grid
stations_pkern = stations_sf['station_num'] |> pkern_snap(dem_pkern, crop_from=TRUE)
#pkern_plot(stations_pkern, reset=FALSE)

# indexing on grid for weather stations
idx_map = stations_pkern[['gval']]
is_obs = !is.na(idx_map)
idx_reorder = idx_map[is_obs] # station numbers in grid vectorized order

# omit those falling outside of prediction region
#stations_cropped_sf = stations_sf[match(idx_reorder, stations_sf[['station_num']]),]
stations_cropped_sf = stations_sf[idx_reorder,]
pkern_plot(stations_pkern, reset=FALSE)
plot(stations_cropped_sf['station_num'], pch=1,  cex=2, add=T)


#'
#' The function `my_ghcn_filter` is for drawing subsets of the dataset within a
#' specified data range. It returns one variable at a time, in the form of a matrix
#' with one row per date.
#'

# print instructions
my_ghcn_filter(ghcn_dt)

# storage for outputs
stor_template = vector(mode='list', length=length(vnames)) |> setNames(vnames)
data_mat_all = lm_elevation_all = pars_interpolation_all = var_pkern_all = stor_template

# extract and analyse each variable in a loop
for(vname in vnames)
{
  cat('\n**********************************\n\n')
  cat(paste('processing variable: ', vname))

  data_mat = my_ghcn_filter(ghcn_dt, vname)
  #str(data_mat)

  # overwrite unrealistic observations with NAs
  is_under = data_mat < reasonable_ranges[[vname]][1]
  is_over = data_mat > reasonable_ranges[[vname]][2]
  n_remove = sum(is_over | is_under, na.rm=TRUE)
  if(n_remove > 0) cat(paste('\nomitting', n_remove, 'unrealistic observation(s)\n'))
  data_mat[is_over | is_under] = NA

  # Different dates form the "layers" here, which are stored in columns of a matrix
  n_layers = ncol(data_mat)
  n_station = nrow(data_mat)

  # observed dates
  times = colnames(data_mat) |> as.Date()
  years = format(times, '%Y') |> as.integer()
  jdates = format(times, '%j') |> as.integer()

  # seasonal index
  seasonal_sin = jdates |> my_season()
  seasonal_cos = jdates |> my_season(s=pi/2)


  #'
  #' The next chunk plots the station locations and the number of missing dates (for precip)
  #' at each one. The circled points are complete-ish time series, with fewer than 5 missing
  #' dates.
  #'

  # count the number of NAs at each station
  n_NA = apply(data_mat, 1, \(x) sum(is.na(x)))
  stations_sf[['n_missing']] = n_NA
  is_complete = n_NA < 5

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

  # snap the complete weather points to a coarse grid
  gres_coarse = c(y=1e3, x=1e3) # spacing in metres
  #coarse_pkern = stations_sf['station_num'] |> pkern_snap(list(gres=gres_coarse))
  coarse_pkern = stations_complete_sf['station_num'] |> pkern_snap(list(gres=gres_coarse))
  is_obs = !is.na(coarse_pkern[['gval']])

  #831.295,    439.347, 114137.734, 103304.213  :: LL = -188151.3941
  #849.1719,    476.7852, 119429.8952, 102116.1375  :: LL = -377318.03353
  #pkern_plot(coarse_pkern)

  # this vector has a station key value at the indices of each observed point, and NAs otherwise
  coarse_pkern[['idx_obs']] = coarse_pkern[['gval']]

  # data matrices must be in column-vectorized order
  idx_reorder = coarse_pkern[['gval']][is_obs]

  # train on all layers
  #idx_train = seq(n_layers)

  idx_train = sample.int(n_layers, 1e3)

  # supply multiple layers as columns in a matrix
  coarse_pkern[['gval']] = resid_mat[idx_reorder, idx_train]

  # center data at layer-wise means
  coarse_means = coarse_pkern[['gval']] |> colMeans(na.rm=TRUE)
  plot(y=coarse_means, x=as.Date(names(coarse_means)))
  coarse_pkern[['gval']] = coarse_pkern[['gval']] |> sweep(2, coarse_means, '-')

  # set NAs to zero

  sum(is.na(coarse_pkern[['gval']]))
  coarse_pkern[['gval']][ is.na(coarse_pkern[['gval']]) ] = 0

  # write list element idx_obs which tells pkern how to unpack gval
  coarse_pkern[['idx_obs']][is_obs] = seq(sum(is_obs))
  #coarse_pkern[['idx_obs']][is_obs] = idx_reorder
  #coarse_pkern[['idx_obs']][is_obs] = seq(sum(is_complete))

  # fit the model
  fit_result = pkern_fit(coarse_pkern, X=0, iso=F)
  pars_interpolation = fit_result$pars

  dev.off()
  fit_result$pars |> pkern_plot_pars(g=coarse_pkern)


  ###
  ###

  cat('\n**********************************\n\n')
  cat(paste('computing kriging variance for: ', vname))

  # # compute the variance
  var_result = pkern_cmean(stations_pkern, pars_interpolation, X=0, out='v', quiet=F)
  var_pkern = stations_pkern |> modifyList(list(gval=var_result))
  #pkern_plot(var_pkern)

  # copy results to storage
  data_mat_all[[vname]] = data_mat
  lm_elevation_all[[vname]] = lm_elevation
  pars_interpolation_all[[vname]] = pars_interpolation
  var_pkern_all[[vname]] = var_pkern
}

# # find a region of acceptable variance
# pkern_plot(var_accept_pkern)
# is_acceptable = var_accept_pkern[['gval']] < 10^2
# var_accept_pkern[['gval']][!is_acceptable] = NA
# pkern_plot(var_accept_pkern)


###
### save preprocessed data for later

#saveRDS(data_mat, 'test.rds')



# lookup table for printing titles
title_lookup = data.frame(vname=c('TMIN', 'TMAX', 'PRCP', 'LOG_PRCP'),
                          vname_print=c(paste('temperature', c('min', 'max'), '(C)'),
                                        rep('precipitation (mm)', 2)))


# inputs
stor = list(

  misc=list(pars_interpolation=pars_interpolation_all,
            title_lookup=title_lookup,
            my_season=my_season,
            lm_elevation=lm_elevation_all,
            idx_map=idx_map,
            uyrw_boundary=uyrw_boundary,
            bbox_dem=bbox_dem,
            stations_sf=stations_sf,
            state_boundaries=state_boundaries),

  pkern = list(dem_pkern=dem_pkern,
               stations_pkern=stations_pkern,
               var_pkern=var_pkern_all,
               data_mat=data_mat_all)

  )

saveRDS(stor, dest_rds_file)

