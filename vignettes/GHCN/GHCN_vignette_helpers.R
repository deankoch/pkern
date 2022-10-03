#' ---
#' title: "GHCND_data_helpers"
#' author: "Dean Koch"
#' date: "`r format(Sys.Date())`"
#' output: github_document
#' ---
#'
#' **Mitacs UYRW project**
#'
#' **pkern**: Fast kriging on gridded datasets
#'
#' Helper functions and global variables for GHCND_* vignette files
#'
#' requirements: pkern, sf, terra, data.table
#'
#' ## GLOBAL VARIABLES
#'

# EPSG codes for projections used in DEM and GHCN data
epsg_ghcnd = 4236
epsg_dem = 32612

# local paths for source data files
basemap_data_dir = 'D:/UYRW_data/for_Patrick'
gchn_source_path = 'D:/ghcnd_from_John/uyrw_Jan_2017_to_Jun_2022.csv' # weather data

# seasonal index and covariates matrix builder
my_season = function(d, s=0, f=1) sin( 2 * pi * f * (as.integer(d) + s) / 365 )
my_season_X = function(dates, n_points=1L)
{
  if(class(dates) == 'Date') dates = format(dates, '%j')
  jday = as.integer(dates)
  data.frame(sin_1 = rep(my_season(jday), each=n_points),
             cos_1 = rep(my_season(jday, s=pi/2), each=n_points),
             sin_2 = rep(my_season(jday, f=2), each=n_points),
             cos_2 = rep(my_season(jday, s=pi/2, f=2), each=n_points))
}

# topographical covariates matrix builder
my_elevation_X = function(elevation, n_layer=1L)
{
  data.frame(elevation = rep(elevation, n_layer), sqrt_elevation = rep(sqrt(elevation), n_layer))
}




# log adjustment constant
log_const = 1

# variable names expected (ignore snow for now) and look-up table for printing titles
vnames = c('TMIN', 'TMAX', 'PRCP', 'LOG_PRCP')
title_lookup = c(paste('temperature', c('min', 'max'), '(C)'), rep('precipitation (mm)', 2))
names(title_lookup) = vnames

# define plausible range for each variable (source data are in 10ths of degree, mm/10)
reasonable_ranges = list(TMIN=c(-6e2, 4e2), TMAX=c(-5e2, 5e2), PRCP=c(0, 16e3))
reasonable_ranges[['LOG_PRCP']] = log(log_const + reasonable_ranges[['PRCP']])

# color scale to use for plotting
my_pal = function(x) hcl.colors(x, 'Spectral', rev=TRUE)
my_breaks = function(x) seq(min(x, na.rm=TRUE), max(x, na.rm=TRUE), length.out=1e3)

#'
#' ## LOADING SOURCE DATA
#'
#' `my_ghcn_dt` loads the big CSV and append some attributes
#'
#' (log-elevation as a predictor may help later to avoid problems with extrapolation in
#' unobserved high alpine areas and log-precipitation produces a variable closer to the
#' Gaussian required in Kriging)
#'
my_ghcn_dt = function(source_path=gchn_source_path)
{
  # load the CSV
  paste0('\nloading CSV at ', gchn_source_path, '...') |> cat()
  ghcn_dt = data.table::fread(source_path)

  # create new attributes: Julian date, year
  cat('\nprocessing attributes...')
  dates = ghcn_dt[['collection_time']] |> as.Date()
  ghcn_dt[['jday']] = format(dates, '%j') |> as.integer()
  ghcn_dt[['year']] = format(dates, '%Y') |> as.integer()

  # create new log-elevation attribute
  ghcn_dt[['sqrt_elevation']] = sqrt(ghcn_dt[['elevation']])

  # create log precip variable
  ghcn_add = ghcn_dt[name=='PRCP']
  ghcn_add[['value']] = log(log_const + ghcn_add[['value']])
  ghcn_add[['name']] = 'LOG_PRCP'
  ghcn_dt = rbind(ghcn_dt, ghcn_add)

  # add key for unique weather stations
  station_ids = ghcn_dt[['station_id']] |> unique()
  ghcn_dt[['station_num']] = ghcn_dt[['station_id']] |> match(station_ids)
  cat(' done\n')
  return(ghcn_dt)
}

#'
#' `my_base_loader` load base layers for spatial reference in plot
#'
# convenience function for loading and transforming sf-compatible inputs
my_base_loader = function(data_dir=basemap_data_dir, epsg_out=epsg_dem, fname=NULL)
{
  # expected filename prefixes (all geojson)
  fname_expect = c('park_boundary', 'watershed_boundary', 'state_boundaries',
                   'station_catchments', 'watershed_flowlines', 'watershed_mainstem')

  # when filename not specified, load all and return in list
  if(is.null(fname))
  {
    fname_list = paste0(fname_expect, '.geojson')
    list_out = lapply(fname_list, \(f) my_base_loader(data_dir=data_dir, epsg_out=epsg_out, fname=f))
    return( setNames(list_out, fname_expect) )
  }

  # load and transform the geometry
  out_sf = data_dir |>
    file.path(fname) |>
    sf::st_read(quiet=TRUE) |>
    sf::st_transform(epsg_out)

  # flow lines are very detailed (difficult to plot). Simplify them
  if( fname=='watershed_flowlines.geojson' ) out_sf = out_sf[out_sf[['streamcalc']] > 2, ] |>
    st_collection_extract('LINESTRING') |>
    st_geometry()

  return(out_sf)

}

#'
#' ## PROCESSING
#'
#' `my_ghcn_fetch` filters the data.table object returned by `my_ghcn_dt`, and
#' creates a matrix of data points with one row per station in `station_ids`
#'
# extract data for one of the weather variables
my_ghcn_filter = function(ghcn_dt, vname=NULL, date_start=NULL, date_end=NULL, clean=TRUE)
{
  if(is.null(vname))
  {
    year_range = range(ghcn_dt[['year']])
    vnames = ghcn_dt[['name']] |> unique()
    '\n select a vname from:' |> paste(paste(vnames, collapse=', ')) |> cat()
    '\n and (optionally) a date range in period:' |> paste(paste(year_range, collapse='-')) |> cat()
    cat('\n')
    return(invisible())
  }

  # find the number of unique stations
  cat(paste('\nreading variable: ', vname))
  n_station = ghcn_dt[['station_id']] |> unique() |> length()

  # filter the data table by variable
  vname_dt = ghcn_dt[name==vname]

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
  out_mat = matrix(NA, nrow=n_station, ncol=length(date_request))
  n_date = length(date_request)
  for( i in seq(n_date) )
  {
    dt_i = vname_dt[ collection_time %in% date_request[i] ]
    out_mat[dt_i[['station_num']], i] = dt_i[['value']]
  }

  # one row per requested date
  rownames(out_mat) = as.character(seq(n_station))
  colnames(out_mat) = as.character(date_request)

  # clean up unrealistic datapoints
  if(clean)
  {
    # overwrite unrealistic observations with NAs
    is_under = out_mat < reasonable_ranges[[vname]][1]
    is_over = out_mat > reasonable_ranges[[vname]][2]
    n_remove = sum(is_over | is_under, na.rm=TRUE)
    if(n_remove > 0) cat(paste('\nomitting', n_remove, 'unrealistic observation(s)\n'))
    out_mat[is_over | is_under] = NA
  }

  return(out_mat)
}

#' `my_stations` creates an sf points object from the coordinates of each station
#'
# return an sf object with station locations
my_stations = function(ghcn_dt, epsg_in=epsg_ghcnd, epsg_out=epsg_dem)
{
  # identify unique stations
  station_ids = ghcn_dt[['station_id']] |> unique()

  # extract the latitude/longitude for each station and covert to sf points object with station key
  idx_first_instance = station_ids |> match(ghcn_dt[['station_id']])
  station_sf = ghcn_dt[idx_first_instance,] |>
    sf::st_as_sf(coords=c('longitude', 'latitude')) |> # order matters here
    sf::st_set_crs(epsg_in) |> # input CRS
    sf::st_transform(epsg_out) |> # transform to output CRS
    cbind(data.frame(station_num = seq(station_ids))) |> # define a key for the stations
    sf::st_set_geometry('geometry')

  # omit non-constant attributes from the station location points object
  nm_return = c('station_num', 'station_id', 'station_name', 'elevation', 'sqrt_elevation', 'geometry')
  station_sf = station_sf[, names(station_sf) %in% nm_return] |> sf::st_set_geometry('geometry')
  return(station_sf)
}

#'
#' ## PLOTTING
#'
# plots the geographical area with some reference points/lines (or adds them to a plot)
my_plot_annotations = function(add=TRUE, bw=FALSE, geo_list=my_base_loader(), simple=FALSE)
{
  # set colors for base layers
  ynp_col = adjustcolor('green', alpha=0.4)
  uyrw_col = adjustcolor('blue', alpha=0.4)
  states_col = adjustcolor('grey', alpha=0.4)
  uyr_col = adjustcolor('blue', alpha=0.4)

  # grey-scale only mode
  col_bw = adjustcolor('black', alpha=0.2)
  if(bw) ynp_col = uyrw_col = states_col = uyr_col = dem_bbox_col = station_col = col_bw

  # create a bounding box and plot if not adding to existing plot
  if(!add)
  {
    bbox_sf = sf::st_bbox(geo_list[['park_boundary']], geo_list[['watershed_boundary']]) |>
      st_as_sfc()

    plot(bbox_sf, border=NA)
  }

  # add state lines
  plot(st_geometry(geo_list[['state_boundaries']]), border=states_col, add=TRUE)

  # add YNP area, upper watershed
  if(simple)
  {
    plot(geo_list[['watershed_boundary']], border=uyrw_col, col=NA, add=TRUE)
    plot(geo_list[['park_boundary']], border=uyrw_col, col=NA, add=TRUE)
    return(invisible())

  } else {

    plot(geo_list[['watershed_boundary']], border=NA, col=uyrw_col, add=TRUE)
    plot(geo_list[['park_boundary']], border=ynp_col, col=NA, add=TRUE)
  }


  # add channels
  plot(geo_list[['watershed_flowlines']], col=uyr_col, add=TRUE)
}






#' returns the data for a given date and variable name. The function uses the objects `ghcn_dt` and `station_sf`
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
