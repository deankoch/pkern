#' ---
#' title: "GHCND_vignette"
#' author: "Dean Koch"
#' date: "`r format(Sys.Date())`"
#' output: github_document
#' ---
#'
#' **Mitacs UYRW project**
#'
#' **pkern**: Fast kriging on gridded datasets
#'
#' This vignette shows how to use `pkern` to interpolate GHCND weather data on the Yellowstone flood.
#' I start with a CSV of the weather station data in the area and use pkern to create a covariance
#' model for fast interpolation. Then I prepare a shiny app for displaying the data on the web.
#'
#'

#+ dependencies_hide, include=FALSE

# load pkern
library(devtools)
load_all()

#+ dependencies

# load extra dependencies for the vignette
library(sf)
library(terra)
library(shiny)
library(data.table)

#'
#' ## Weather data
#'
#' The weather station data are stored in a large CSV file. I start by opening this
#' file, reshaping its contents, and creating an sf points object from the coordinates.
#'

#+ ghcn_open, include=FALSE

# open the file
gchn_source_path = 'D:/ghcnd_from_John/uyrw_Jan_2017_to_Jun_2022.csv'
ghcn_dt = fread(gchn_source_path)

# identify unique stations, variable names
station_ids = ghcn_dt[['station_id']] |> unique()
vnames = ghcn_dt[['name']] |> unique()

# create new attributes for Julian date and year
dates = ghcn_dt[['collection_time']] |> as.Date()
ghcn_dt[['jday']] = format(dates, '%j') |> as.integer()
ghcn_dt[['year']] = format(dates, '%Y') |> as.integer()

# seasonality represented by cosine with minimum centered at Dec 21
my_season = function(d) -cos( 2 * pi * (as.integer(format(as.Date(d), '%j')) + 10) / 365 )
ghcn_dt[['season']] = my_season(dates)

#


# EPSG code for WGS84 projection used in GHCN data and for the DEM loaded later
epsg_ghcnd = 4236
epsg_dem = 32612

# pull the latitude/longitude for each station and covert to sf points object
idx_first_instance = station_ids |> match(ghcn_dt[['station_id']])
station_sf = ghcn_dt[idx_first_instance,] |>
  sf::st_as_sf(coords=c('longitude', 'latitude')) |> # order matters here
  sf::st_set_crs(epsg_ghcnd) |> # input CRS
  sf::st_transform(epsg_dem) |> # transform to output CRS
  sf::st_geometry() |> # drop attributes
  cbind(data.frame(station_num = seq(station_ids))) |> # define a key for the stations
  sf::st_set_geometry('geometry')

#'
#' Regression on elevation
#'
#'

# create log precip variable
ghcn_add = ghcn_dt[name=='PRCP']
ghcn_add[['value']] = log(1 + ghcn_add[['value']])
ghcn_add[['name']] = 'LOG_PRCP'
ghcn_dt = rbind(ghcn_dt, ghcn_add)
vnames = c('LOG_PRCP', vnames)

# regression on three covariates: elevation, year, and time of year
vname_lm_covariates = c('elevation', 'season', 'year')
formula_lm = paste('value~', paste(vname_lm_covariates, collapse='+')) |> as.formula()
lm_by_vname = stats::setNames(nm=vnames) |> lapply(\(vn) lm(formula_lm, data=ghcn_dt[name==vn]))


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
#' These plots aren't much use without some clues about location and scale, so I'll
#' load some reference data about the area to use in annotating plots
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
basemap_data_dir = 'D:/UYRW_data/for_Patrick'
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

# load the DEM separately
dem_src_rast = basemap_data_dir |> file.path('dem.tif') |> terra::rast()
dem_bbox_sf = dem_src_rast |> sf::st_bbox() |> sf::st_as_sfc()

# make a bounding box for everything
outer_bbox = sf::st_bbox(dem_src_rast, station_sf, ynp_sf, uyrw_sf) |> sf::st_as_sfc()

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
#' ## other dates
#'
#' and improved plot
#'

day = '2022-06-14'
ghcn_example = my_ghcn_pts(day, vname)

# regression residuals from this example
ghcn_example[['z']] = ghcn_example[['value']] - predict(lm_by_vname[[vname]], newdata=ghcn_example)
ghcn_example = sf::st_set_geometry(ghcn_example, 'geometry')
zobs_pkern = ghcn_example['z'] |> sf::st_crop(dem_bbox_sf) |> pkern_snap(dem_pkern)

# linear predictor for this example
X_pred$season = ghcn_example[['season']][1]
X_pred$year = ghcn_example[['year']][1]
lm_pred = predict(lm_by_vname[[vname]], newdata=X_pred)
lm_pred_pkern = zobs_pkern |> modifyList(list(gval=lm_pred))

# interpolate the anomaly
#fit_result_SK$pars$eps=2.101979

fit_result_cmean = pkern_cmean(zobs_pkern, fit_result_SK$pars, X=0)
zpred_pkern = zobs_pkern |> modifyList(list(gval=fit_result_cmean))
pkern_plot(zpred_pkern)

# use variance to adjust for bias on original scale
range(fit_result_cmean + lm_pred, na.rm=T)
range(1+ghcn_example[['value']])


krig_original = exp(fit_result_cmean + lm_pred + var_result/2) - 1
krig_original_pkern = zobs_pkern |> modifyList(list(gval=krig_original/10))
pkern_plot(krig_original_pkern, zlab='precip (mm)')

# add YNP area, upper watershed
plot(uyrw_sf, border='black', add=TRUE)
plot(ynp_sf, border='black', add=TRUE)
plot(channels_sf, col='black', add=TRUE)
station_sf |> sf::st_geometry() |> sf::st_crop(dem_bbox_sf) |> plot(col='black', add=TRUE, pch=1, cex=0.75)





range(1+exp(ghcn_example[['value']]))/10



#'
#' ## temperature example
#'

day = '2022-01-01'
vname = 'TMAX'
ghcn_example_2 = my_ghcn_pts(day, vname)

# regression residuals from this example
ghcn_example_2[['z']] = ghcn_example_2[['value']] - predict(lm_by_vname[[vname]], newdata=ghcn_example_2)
ghcn_example_2 = sf::st_set_geometry(ghcn_example_2, 'geometry')
zobs_pkern = ghcn_example_2['z'] |> sf::st_crop(dem_bbox_sf) |> pkern_snap(dem_pkern)

# linear predictor for this example
X_pred$season = ghcn_example[['season']][1]
X_pred$year = ghcn_example[['year']][1]
lm_pred = predict(lm_by_vname[[vname]], newdata=X_pred)
lm_pred_pkern = zobs_pkern |> modifyList(list(gval=lm_pred))
pkern_plot(lm_pred_pkern)

# interpolate the anomaly
fit_result_cmean = pkern_cmean(zobs_pkern, fit_result_SK$pars, X=0)
zpred_pkern = zobs_pkern |> modifyList(list(gval=fit_result_cmean))
pkern_plot(zpred_pkern)

# add them together to get prediction
krig_pkern = zobs_pkern |> modifyList(list(gval=fit_result_cmean + lm_pred))
pkern_plot(krig_pkern)





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
