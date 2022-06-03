# Meuse example

# replace this with pkern load before pushing
# ***
library(devtools)
load_all()

library(sp)
library(sf)
library(terra)

## ## ## ## ## ## ## ## ## ## ## ##  MEUSE EXAMPLE

# meuse example data is included in this package

# load the Meuse data into a convenient format
get_meuse = function(dfMaxLength = units::set_units(50, m))
{
  # EPSG code for the coordinate system
  epsg_meuse = 28992

  # open river location data
  utils::data(meuse.riv)
  crs_meuse = sf::st_crs(epsg_meuse)[['wkt']]

  # reshape the river (edge) point data as a more densely segmented polygon
  colnames(meuse.riv) = c('x', 'y')
  meuse_river_points = sf::st_as_sf(as.data.frame(meuse.riv), coords=c('x', 'y'), crs=crs_meuse)
  meuse_river_seg = sf::st_cast(sf::st_combine(meuse_river_points), 'LINESTRING')
  meuse_river_poly = sf::st_cast(st_segmentize(meuse_river_seg, dfMaxLength), 'POLYGON')

  # skeletonization trick to get a single linestring at center of the river
  meuse_river_voronoi = sf::st_cast(sf::st_voronoi(meuse_river_poly, bOnlyEdges=TRUE), 'POINT')
  meuse_river_skele = sf::st_intersection(meuse_river_voronoi, meuse_river_poly)
  n_skele = length(meuse_river_skele)

  # compute distance matrix
  dmat_skele = units::drop_units(sf::st_distance(meuse_river_skele))

  # re-order to start from northernmost point
  idx_first = which.max(st_coordinates(meuse_river_skele)[,2])
  idx_reorder = c(idx_first, integer(n_skele-1L))
  for(idx_skele in seq(n_skele-1L))
  {
    # find least distance match
    idx_tocheck = seq(n_skele) != idx_first
    idx_first = which(idx_tocheck)[ which.min(dmat_skele[idx_tocheck, idx_first]) ]
    idx_reorder[1L+idx_skele] = idx_first

    # modify distance matrix so the matching point is not selected again
    dmat_skele[idx_first, ] = Inf
  }

  # connect the points to get the spine
  meuse_river = sf::st_cast(sf::st_combine(meuse_river_skele[idx_reorder]), 'LINESTRING')

  # load soil points data
  utils::data(meuse)
  meuse_soils = sf::st_as_sf(meuse, coords=c('x', 'y'), crs=epsg_meuse)

  # add 'distance' (to river) and 'logzinc' columns
  meuse_soils[['distance']] = units::drop_units( sf::st_distance(meuse_soils, meuse_river))
  meuse_soils[['log_zinc']] = log(meuse_soils[['zinc']])

  # crop the river objects to buffered bounding box of soils data
  bbox_padded = st_buffer(sf::st_as_sfc(sf::st_bbox(meuse_soils)), units::set_units(500, m))
  meuse_river_poly = sf::st_crop(meuse_river_poly, bbox_padded)
  meuse_river = sf::st_crop(meuse_river, bbox_padded)

  # return three geometry objects in a list
  return( list(soils=meuse_soils, river_poly=meuse_river_poly, river_line=meuse_river) )
}

## load the data

# load data and plot using sf package
meuse = get_meuse()
plot(meuse[['river_poly']], col='lightblue', reset=FALSE)
plot(meuse[['river_line']], lwd=2, add=TRUE)
plot(meuse[['soils']]['zinc'], pch=16, add=TRUE)
n_meuse = nrow(meuse[['soils']])

## define a grid of a certain resolution and snap points to it

# desired resolution in units of metres
gres = c(y=25, x=25)

# snap points, copying values of dependent variable
g_meuse = pkern_snap(meuse[['soils']]['log_zinc'], g=list(gres=gres))
gdim = g_meuse[['gdim']]
is_obs = !is.na(g_meuse[['gval']])
n_obs = sum(is_obs)

# plot with source points indicated over their snapped grid location
pkern_plot(g_meuse, zlab='log(zinc) (ppm)', reset=FALSE)
plot(sf::st_geometry(meuse[['soils']]), add=TRUE)

# ordinary kriging: fit isotropic gaussian model by default
fit_result_OK = pkern_fit(g_obs=g_meuse, pars='gau', X=NA, quiet=TRUE)



pars_OK = fit_result_OK$pars

# plot the fitted kernel
pkern_plot_pars(pars_OK, g_meuse)

# plot the conditional mean of the spatial error
z_pred = pkern_cmean(g_obs=g_meuse, pars=pars_OK, X=NA)
g_meuse_pred = modifyList(g_meuse, list(gval=z_pred))
pkern_plot(g_meuse_pred)

## universal kriging: Adding a predictor yields a better fit

# get predictor values for entire grid
g_meuse_sf = pkern_coords(g_meuse, out='sf')
d2r_result = units::drop_units(st_distance(g_meuse_sf, meuse[['river_line']]))
#meuse_predictors = scale(cbind(d2r_result, sqrt(d2r_result)))
meuse_predictors = scale(cbind(d2r_result, sqrt(d2r_result)))

# subset of predictors at observed response locations
X = matrix(meuse_predictors[is_obs,], nrow=n_obs)

# snap distances to the grid then reshape into covariate matrix
# g_distance = pkern_snap(meuse[['soils']]['distance'], g=g_meuse)
# X_all = matrix(c(g_distance[['gval']], exp(-g_distance[['gval']])), ncol=2) |> scale()

# make a covariates matrix
#X = matrix(X_all[!is.na(g_meuse[['gval']]),], ncol=2)

# fit the model
fit_result_UK = pkern_fit(g_obs=g_meuse, pars='gau', X=X, quiet=TRUE)
pars_UK = fit_result_UK$pars
pkern_GLS(g_meuse, pars_UK, X=X, out='b')

# kriging over full grid

# GLS to get trend
z_gls = pkern_GLS(g_meuse, pars_UK, X=meuse_predictors, out='z')
g_meuse_gls = modifyList(g_meuse, list(gval=z_gls))
pkern_plot(g_meuse_gls)

# spatial error
g_meuse_detrend = modifyList(g_meuse, list(gval=g_meuse[['gval']]-z_gls))
z_spat = pkern_cmean(g_meuse_detrend, pars_UK, X=0)
g_meuse_spat = modifyList(g_meuse, list(gval=z_spat))
pkern_plot(g_meuse_spat)

# predictions
z_pred = z_gls + z_spat
g_meuse_pred = modifyList(g_meuse, list(gval=z_pred))
pkern_plot(g_meuse_pred)

# prediction variance
z_var = pkern_cmean(g_meuse_detrend, pars_UK, X=0, out='v', quiet=TRUE)
g_meuse_var = modifyList(g_meuse, list(gval=z_var))
pkern_plot(g_meuse_var)

# prediction bias adjustment from log scale
z_pred2 = exp(z_pred + z_var/2)
g_meuse_pred2 = modifyList(g_meuse, list(gval=z_pred2))
pkern_plot(g_meuse_pred2)

# mask outside observed range
pkern_plot(g_meuse_pred2, zlim=range(meuse[['soils']][['zinc']]))


##

# mask all distances outside observed range
distance_masked = ( meuse_predictors[,1] < min(X[,1]) ) | ( meuse_predictors[,1] > max(X[,1]) )
z_pred_masked = z_pred2
z_pred_masked[distance_masked] = NA
g_zinc_masked = modifyList(g_meuse, list(gval=z_pred_masked))
pkern_plot(g_zinc_masked)

# mask all distances outside observed range
z_var_masked = z_var
z_var_masked[distance_masked] = NA
g_zinc_var_masked = modifyList(g_meuse, list(gval=z_var_masked))
pkern_plot(g_zinc_var_masked)



##

# convert to SpatRaster
meuse_pred_rast = pkern_export(g_meuse_pred_zinc)
meuse_pred_var_rast = pkern_export(g_meuse_var)

river_rast = terra::rasterize(as(meuse[['river_poly']], 'SpatVector'), meuse_pred_rast)
terra::mask(meuse_pred_rast, river_rast, inverse=TRUE) |> pkern_grid() |> pkern_plot()
terra::mask(meuse_pred_var_rast, river_rast, inverse=TRUE) |> pkern_grid() |> pkern_plot()




#
