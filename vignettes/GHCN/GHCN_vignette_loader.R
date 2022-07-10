#' ---
#' title: "GHCND_data_vignette"
#' author: "Dean Koch"
#' date: "`r format(Sys.Date())`"
#' output: github_document
#' ---
#'
#' **Mitacs UYRW project**
#'
#' **pkern**: Fast kriging on gridded datasets
#'
#' This vignette accompanies the "GHCND_data" example, and is required to generate the
#' data used there.
#'
#' ADD SOME INFO ABOUT DATA SOURCES HERE
#'
#' ## Outputs
#'
#' This script create spatial data files
#'
#' 1. ghcn_dt
#' 2. station_sf
#'
#'
#'
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

# print instructions then extract all precipitation data
my_ghcn_filter(ghcn_dt)
precip_mat = my_ghcn_filter(ghcn_dt, vname='LOG_PRCP')
str(precip_mat)

# Different dates form the "layers" here, which are stored in columns of a matrix
n_layers = ncol(precip_mat)
n_station = nrow(precip_mat)

# observed dates
times = colnames(precip_mat) |> as.Date()
years = format(times, '%Y') |> as.integer()
jdates = format(times, '%j') |> as.integer()

# seasonal index
my_season = function(d, p=1) -cos( p * 2 * pi * (d + 10) / 365 )
seasonal = jdates |> my_season(p=1)
seasonal_half = jdates |> my_season(p=2)
seasonal_quarter = jdates |> my_season(p=4)

#'
#' Another helper function, `my_stations`, creates an `sf` POINT collection from the station
#' positions. Most of the stations are missing data on certain days (see below).
#'

# create the station points
stations_sf = my_stations(ghcn_dt)

# count the number of NAs at each station
n_NA = apply(precip_mat, 1, \(x) sum(is.na(x)))
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
  c(list(year=rep(years, n_station),
         seasonal=rep(seasonal, n_station),
         seasonal_2=rep(seasonal_half, n_station),
         seasonal_4=rep(seasonal_quarter, n_station)))



# subset with complete response data
y_obs = precip_mat[is_complete,] |> as.vector()
X_obs = stations_complete_sf[c('elevation', 'log_elevation')] |>
  sf::st_drop_geometry() |>
  lapply(function(x) rep(x, n_layers)) |>
  c(list(year=rep(years, sum(is_complete)),
         seasonal=rep(seasonal, sum(is_complete)),
         seasonal_2=rep(seasonal_half, sum(is_complete)),
         seasonal_4=rep(seasonal_quarter, sum(is_complete))))

# simple linear regression
lm_elevation = lm(y~elevation+log_elevation+year+seasonal+seasonal_2+seasonal_4, c(list(y=y_obs), X_obs))
summary(lm_elevation)

# matrices with fitted values and residuals for each layer in columns
pred_mat = predict(lm_elevation, newdata=X_obs_all) |> matrix(nrow=nrow(precip_mat))
resid_mat = pred_mat - precip_mat

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

idx_train = seq(n_layers)
#idx_train = sample.int(n_layers, 1e2)
#idx_train = resid_mat[idx_reorder,] |> apply(2, \(x) sum(x>0) > 10) |> which() |> sample(size=2)
#
# #
# i = 2
# stations_complete_sf[['test']] = resid_mat[,idx_train[i]]
# coarse_pkern2 = stations_sf['test'] |> pkern_snap(list(gres=gres_coarse))
# pkern_plot(coarse_pkern2)


# supply multiple layers as columns in a matrix
coarse_pkern[['gval']] = resid_mat[idx_reorder, idx_train]

# scale by layer-wise means

coarse_pkern[['gval']] = coarse_pkern[['gval']] |> apply(2, \(x) x - mean(x, na.rm=TRUE))

# supply only the non-NA values by including idx
coarse_pkern[['idx_obs']][is_obs] = seq(sum(is_complete))

#
fit_result = pkern_fit(coarse_pkern, X=0, iso=FALSE)



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
#pkern_plot(var_pkern)


# select an index (date) to start at
idx = 1950
idx = idx-1

my_pal = function(x) hcl.colors(x, 'Spectral', rev=T)
my_col = function(x) my_pal(1e2)[ num_to_cent(x) ]

bbox_dem = st_bbox(dem_rast) |> st_as_sfc()

if(1)
{
  idx = idx + 1
  #dev.off()

  idx_date = as.Date(colnames(precip_mat)[idx])
  yr = format(idx_date, '%Y') |> as.integer()
  jday = format(idx_date, '%j') |> as.integer()
  s1 = my_season(jday)
  s2 = my_season(jday, 2)
  s4 = my_season(jday, 4)

  # combine all predictors
  X_new = data.frame(elevation = x_dem,
                     log_elevation = log(x_dem),
                     year = rep(yr, length(x_dem)),
                     seasonal = rep(s1, length(x_dem)),
                     seasonal_2 = rep(s2, length(x_dem)),
                     seasonal_4 = rep(s4, length(x_dem)))

  # compute their linear combination
  y_dem = predict(lm_elevation, newdata=X_new)


  ## log scale plots

  # plot the station data for this date as points
  #stations_sf[['precip_test']] = precip_mat[,idx]
  #plot(stations_sf['precip_test'], pch=16, cex=2, reset=FALSE, pal=my_pal)

  #bbox_dem = st_bbox(dem_rast) |> st_as_sfc()
  #plot(bbox_dem, border='black', add=TRUE)

  #y_test = precip_mat[,idx][idx_reorder]
  #stations_pkern[['gval']][is_obs] = y_test

  #dev.off()
  #pkern_plot(stations_pkern)


  # plot the residuals for the complete data regression
  #stations_sf[['resid_test']] = resid_mat[,idx]
  #plot(stations_sf['resid_test'], pch=16, cex=2, reset=FALSE)

  #bbox_dem = st_bbox(dem_rast) |> st_as_sfc()
  #plot(bbox_dem, border='black', add=TRUE)

  # compute all available linear model residuals
  y_residual = y_dem[is_obs] - precip_mat[,idx][idx_reorder]
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

  # kriging prediction on log scale
  krig = y_dem - fit_result_cmean

  # use variance to adjust for bias on original scale
  krig_original = (exp(krig + var_result/2) - 1) / 10
  krig_original_pkern = stations_pkern |> modifyList(list(gval=krig_original))
  #pkern_plot(krig_original_pkern, zlab='precip (mm)')

  # create original weather points for reference
  stations_sf[['precip_test']] = ( exp(precip_mat[,idx]) - 1 ) / 10
  plot_sf = stations_sf['precip_test'][idx_reorder,]
  plot_sf = plot_sf[!is.na(plot_sf$precip_test),]
  zlim_min = min(c(plot_sf[['precip_test']], krig_original), na.rm=T)
  zlim_max = max(c(plot_sf[['precip_test']], krig_original), na.rm=T)
  zlim = c(zlim_min, zlim_max)

  pkern_plot(krig_original_pkern,
             main=colnames(precip_mat)[idx],
             zlab='precip (mm)',
             zlim=zlim,
             reset=FALSE)

}

mybreaks = seq(zlim_min, zlim_max, by=diff(zlim)/1e2)
#pkern_plot(krig_original_pkern, zlab='precip (mm)', reset=FALSE)
plot(st_geometry(plot_sf), cex=2, add=TRUE)
plot(plot_sf, pch=16, cex=1, add=TRUE, pal=my_pal, breaks=mybreaks)

#my_plot_annotations(geo_list=geo_list, bw=TRUE)
plot(st_crop(geo_list[['watershed_boundary']], bbox_dem), col=NA, add=TRUE)
st_geometry(geo_list[['state_boundaries']]) |> st_crop(bbox_dem) |> plot(border='grey20', col=NA, add=TRUE)


#my_plot_annotations(add=TRUE)










# test on a single layer the old fashioned way
idx = 1000
idx = idx + 1
stations_sf[['test']] = precip_mat[,idx]
plot(stations_sf['test'], pch=16, cex=2)

fine_pkern[['gval']]
pkern_plot(coarse_pkern2)

# interpolate the anomaly
fit_result_cmean = pkern_cmean(coarse_pkern2, fit_result$pars, X=0)
zpred_pkern = coarse_pkern2 |> modifyList(list(gval=fit_result_cmean))
pkern_plot(zpred_pkern)


#
# g_obs = coarse_pkern
# pars = NULL
# X = 0
# iso = TRUE
# initial = NULL
# quiet = FALSE
#
# g_obs=g
# pars_fix=pars
# X=0
# control=list()
# log_scale=TRUE
# method='L-BFGS-B'
#
# p = pars_df[,2]
#
# pars=pars_complete
# method='chol'
# fac=NULL
# more=FALSE
#
#
# pkern_plot(coarse_pkern)






#'
#' extract data for one of the variables
#'


ghcn_fetch = function(vname=NULL, date_start=NULL, date_end=NULL, ghcn=ghcn_dt)
{
  if(is.null(vname))
  {
    year_range = range(ghcn[['year']])
    vnames = ghcn[['name']] |> unique()
    '\n select a vname from:' |> paste(paste(vnames, collapse=', ')) |> cat()
    '\n and (optionally) a date range in period:' |> paste(paste(year_range, collapse='-')) |> cat()
    return(invisible())
  }

  # filter the data table by variable
  vname_dt = ghcn[name==vname]

  # range for allowable start/end times
  year_range = range(vname_dt[['year']])
  date_min = year_range[1] |> paste(min(vname_dt[year==year_range[1]][['jday']])) |> as.Date('%Y %j')
  date_max = year_range[2] |> paste(max(vname_dt[year==year_range[2]][['jday']])) |> as.Date('%Y %j')
  date_allow = seq.Date(date_min, date_max, by='day')

  # set defaults
  if(is.null(date_start)) date_start = date_min
  if(is.null(date_end)) date_end = date_max
  date_request = date_allow[ date_allow %in% seq.Date(as.Date(date_start), as.Date(date_end), by='day')]

  # copy data to matrix
  out_mat = matrix(NA, nrow=length(date_request), ncol=length(station_ids))
  n_date = length(date_request)
  for( i in seq(n_date) )
  {
    dt_i = vname_dt[ collection_time %in% date_request[i] ]
    out_mat[i, dt_i[['station_num']]] = dt_i[['value']]
  }

  # one row per requested date
  rownames(out_mat) = as.character(date_request)
  return(out_mat)
}

# Example: precipitation (total in mm) on New Year's Day
vname = 'LOG_PRCP'




ghcn_mat = ghcn_fetch(vname)

# Note that the majority of stations are missing data on some or all of the allowed dates


ghcn_mat |> apply(2, \(x) sum(is.na(x)))

# identify stations with complete records



# filter to requested dates
#vname_dt[collection_time %in% date_request]









# filter the data table and match to existing unique points object
ghcnd_dt_filtered = ghcn_dt[collection_time==day][name==vname]
idx_filtered = ghcnd_dt_filtered[['station_id']] |> match(station_ids)
ghcnd_dt_filtered[['station_num']] = idx_filtered

# build and return the output sf object
sf::st_set_geometry(ghcnd_dt_filtered, sf::st_geometry(station_sf)[idx_filtered])



#'
#' For convenience I define a function that returns the data for a given
#' date and variable name. The function uses the objects `ghcn_dt` and `station_sf`
#' that we just loaded into the global R environment.
#'

# define a function to return a point collection for a given date/variable
my_ghcn_pts = function(day, vname='TMAX')
{
  # filter the data table and match to existing unique points object
  ghcnd_dt_filtered = ghcn_dt[collection_time==day][name==vname]
  idx_filtered = ghcnd_dt_filtered[['station_id']] |> match(station_ids)
  ghcnd_dt_filtered[['station_num']] = idx_filtered

  # build and return the output sf object
  sf::st_set_geometry(ghcnd_dt_filtered, sf::st_geometry(station_sf)[idx_filtered])
}

#'
#' Individual dates and variables can now be accessed with a call to my_ghcn_pts.
#' This returns an sf object that can be plotted directly
#'

# Example: precipitation (total in mm) on New Year's Day
day = '2022-01-01'
vname = 'LOG_PRCP'


ghcn_example = my_ghcn_pts(day, vname)
plot(ghcn_example['value'], pch=16)

# regression residuals from this example
ghcn_example[['z']] = ghcn_example[['value']] - predict(lm_by_vname[[vname]], newdata=ghcn_example)
ghcn_example = sf::st_set_geometry(ghcn_example, 'geometry')
plot(ghcn_example['z'], pch=16)






#'
#' Now create an sf points object from the coordinates of each station
#'

# load the DEM separately
dem_src_rast = basemap_data_dir |> file.path('dem.tif') |> terra::rast()
dem_bbox_sf = dem_src_rast |> sf::st_bbox() |> sf::st_as_sfc()















#'
#' Plots should include some clues about location and scale, so I load some reference data
#' about the area to use as base layers and annotations
#'

# convenience function for loading and transforming inputs
my_base_loader = function(data_dir, fname, epsg=epsg_dem)
{
  data_dir |>
    file.path(fname) |>
    sf::st_read() |>
    sf::st_transform(epsg)
}

# load the reference data and transform to the projection of the weather data
ynp_sf = basemap_data_dir |> my_base_loader('park_boundary.geojson')
uyrw_sf = basemap_data_dir |> my_base_loader('watershed_boundary.geojson')
states_sf = basemap_data_dir |> my_base_loader('state_boundaries.geojson')
sbasins_sf = basemap_data_dir |> my_base_loader('station_catchments.geojson')
channels_sf = basemap_data_dir |> my_base_loader('watershed_flowlines.geojson')
mainstem_sf = basemap_data_dir |> my_base_loader('watershed_mainstem.geojson')

# flow lines are very detailed (difficult to plot). Simplify them
channels_sf = channels_sf[channels_sf[['streamcalc']] > 2, ] |>
  st_collection_extract('LINESTRING') |>
  st_geometry()

#'
#' Now I can define another convenience function to add annotations to an existing
#' plot. This function looks in the global environment for the `sf` and `terra` objects
#' loaded above and adds them to the plot.
#'

my_plot_annotations = function(add=TRUE, bw=FALSE)
{
  # set colors for base layers
  ynp_col = adjustcolor('green', alpha=0.4)
  uyrw_col = adjustcolor('blue', alpha=0.4)
  states_col = adjustcolor('grey', alpha=0.4)
  #subbasin_col = adjustcolor('white', alpha=0.4)
  uyr_col = adjustcolor('blue', alpha=0.4)
  dem_bbox_col = adjustcolor('grey', alpha=0.4)
  station_col = adjustcolor('red', alpha=0.4)

  if(bw)
  {
    ynp_col = uyrw_col = states_col = uyr_col = dem_bbox_col = station_col = adjustcolor('black', alpha=0.2)
  }

  # first plot call initializes the plot with add=FALSE
  plot(st_geometry(dem_bbox_sf), col=dem_bbox_col, border=NA, add=add)

  # add state lines
  plot(st_geometry(states_sf), border=states_col, add=TRUE)

  # add YNP area, upper watershed
  plot(uyrw_sf, border=NA, col=uyrw_col, add=TRUE)
  plot(ynp_sf, border=NA, col=ynp_col, add=TRUE)
  #plot(st_geometry(sbasins_sf), border=subbasin_col, add=TRUE)

  # add channels
  plot(channels_sf, col=uyr_col, add=TRUE)

  # add point locations indicating the weather stations
  plot(st_geometry(station_sf), col=station_col, add=TRUE, pch=20, cex=0.75)
}

# display the base layers
my_plot_annotations(add=F)








#'
#' Regression on elevation
#'



# regression on four covariates: elevation and its logarithm, year, and time of year
vname_lm_covariates = c('elevation', 'log_elevation', 'season', 'year')
formula_lm = paste('value~', paste(vname_lm_covariates, collapse='+')) |> as.formula()
lm_by_vname = stats::setNames(nm=vnames) |> lapply(\(vn) lm(formula_lm, data=ghcn_dt[name==vn]))





#'
#' ## Define the scale of interest
#'
#' `pkern` works on grids of points. Since we have a spatial covariate that is already
#' gridded - elevation, in the DEM - I will use a simplified version of the DEM grid here.
#' This grid establishes the point locations at which `pkern` will later compute expected
#' weather values (by interpolation).
#'
#' The source DEM grid is very detailed, having around 100 million points. So I first
#' upscale by a factor of 25 to get a coarser grid with around 150 thousand points
#' (the 25 relates the square roots of these numbers).
#'

# upscale the DEM grid for this demo. Later we will use it at full resolution
dem_agg_fac = 25
dem_rast = terra::aggregate(dem_src_rast, fact=dem_agg_fac)

# plot the result
plot(dem_rast)

# convert pkern object
dem_pkern = pkern_grid(dem_rast)
dem_bbox_sf

#'
#' This is a 450 X 352 grid, with grid point spacing around 700m. It covers the same
#' extent as the original DEM data grid, but at a lower detail level.
#'
#' To snap the weather points to the grid, we have `pkern_snap`
#'

# first call snaps the station locations to the grid (copying station_num)
station_num_pkern = station_sf |> sf::st_crop(dem_bbox_sf) |> pkern_snap(dem_pkern)

#



zobs_pkern = ghcn_example['z'] |> sf::st_crop(dem_bbox_sf) |> pkern_snap(dem_pkern)
pkern_plot(zobs_pkern)

#'
#' ...so that the covariogram can be more easily estimated:
#'

X = ghcn_example[vname_lm_covariates] |> sf::st_drop_geometry() |> as.data.frame() |> as.matrix()

# fit a Gaussian covariogram
fit_result_SK = pkern_fit(g_obs=zobs_pkern, pars='gau', X=dem_pkern[['gval']], quiet=TRUE)

fit_result_SK = pkern_fit(g_obs=zobs_pkern, pars='gau', X=0, quiet=TRUE)


# plot the fitted kernel shape
pkern_plot_pars(fit_result_SK$pars, zobs_pkern)

#'
#' ## interpolation
#'

# linear predictor for this example
X_pred = data.frame(elevation=dem_pkern[['gval']], season=ghcn_example[['season']][1], year=ghcn_example[['year']][1])
lm_pred = predict(lm_by_vname[[vname]], newdata=X_pred)
lm_pred_pkern =  zobs_pkern |> modifyList(list(gval=lm_pred))
pkern_plot(lm_pred_pkern)

# interpolate the anomaly
fit_result_cmean = pkern_cmean(zobs_pkern, fit_result_SK$pars, X=0)
zpred_pkern = zobs_pkern |> modifyList(list(gval=fit_result_cmean))
pkern_plot(zpred_pkern)

# add them together to get prediction on log scale
krig_result = fit_result_cmean + lm_pred
krig_pkern = zobs_pkern |> modifyList(list(gval=krig_result))
pkern_plot(krig_pkern)

# compute variance
var_result = pkern_cmean(zobs_pkern, fit_result_SK$pars, X=0, out='v', quiet=F)
var_pkern = zobs_pkern |> modifyList(list(gval=var_result))
pkern_plot(var_pkern)

# use variance to adjust for bias on original scale
krig_original = exp(krig_result + var_result/2)
krig_original_pkern = zobs_pkern |> modifyList(list(gval=krig_original/10))
pkern_plot(krig_original_pkern, zlab='precip (mm)')



#'
#'
#' ## Extra
#'

# load the DEM separately
dem_src_rast = basemap_data_dir |> file.path('dem.tif') |> terra::rast()
dem_bbox_sf = dem_src_rast |> sf::st_bbox() |> sf::st_as_sfc()

# make a bounding box for everything
outer_bbox = sf::st_bbox(dem_src_rast, station_sf, ynp_sf, uyrw_sf) |> sf::st_as_sfc()


# create new attributes: Julian date, year, seasonality (cosine with minimum at Dec 21)
ghcn_dt[['jday']] = format(dates, '%j') |> as.integer()
ghcn_dt[['year']] = format(dates, '%Y') |> as.integer()
my_season = function(d) -cos( 2 * pi * (as.integer(format(as.Date(d), '%j')) + 10) / 365 )
ghcn_dt[['season']] = my_season(dates)


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
  path.input = here('vignettes/meuse_vignette.R')
  path.output = here('vignettes/meuse_vignette.md')
  path.garbage = here('vignettes/meuse_vignette.html')
  rmarkdown::render(path.input, clean=TRUE, output_file=path.output)
  unlink(path.garbage)

  # substitute local file paths for image files with URLs on github
  md.github = gsub('D:/pkern', 'https://github.com/deankoch/pkern/blob/main', readLines(path.output))
  writeLines(md.github, path.output)
}
