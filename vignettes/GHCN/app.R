library(shiny)

# load pkern
library(devtools)
load_all()

#+ dependencies

# extra dependencies for loading and preparing spatial data
library(sf)
library(data.table)

source('D:/pkern/vignettes/GHCN/GHCN_vignette_helpers.R')

# path to the output file written by GHCN_vignette_preprocessing.R
rds_data_file = 'D:/pkern/vignettes/GHCN/GHCN_preprocessing_data.rds'

# load the data
stor = readRDS(rds_data_file)
vnames = names(stor[['pkern']][['data_mat']])

# copy some of the small loaded objects
title_lookup = stor[['misc']][['title_lookup']]
idx_grid = stor[['misc']][['idx_grid']]
stations_sf = stor[['misc']][['stations_sf']]
stor[['misc']][['stations_sf']]

# layer with station mapping
#stations_pkern |> str()
# or: stations_pkern = dem_pkern |> modifyList(list(gval=idx_grid))

# grid elevation values
x_dem = stor[['pkern']][['dem_pkern']][['gval']]
n_dem = length(x_dem)
print(n_dem)

# derive mapping keys on grid
#is_obs = !is.na(idx_grid)
#idx_reorder = idx_grid[is_obs] # station numbers in grid vectorized order
stations_pkern = stor[['pkern']][['dem_pkern']] |> modifyList(list(gval=idx_grid))

# omit stations not appearing in grid
station_displayed = stations_sf[['station_num']] %in% idx_grid

# plotting parameters
my_pal = function(x) hcl.colors(x, 'Spectral', rev=T)
my_col = function(x) my_pal(1e2)[ num_to_cent(x) ]
dates_all_ini = as.Date(colnames(stor[['pkern']][['data_mat']][['TMIN']]))

# combine static predictors
X_new = data.frame(elevation = x_dem,
                   sqrt_elevation = sqrt(x_dem))

X_new_obs = data.frame(elevation = stations_sf[['elevation']],
                       sqrt_elevation = sqrt(stations_sf[['elevation']]))

# # adjust interpolation parameters
# for(v in vnames) {
#
#   stor$misc$pars_interpolation[[v]]
#   #pars_interpolation2 = pars_interpolation
#   stor$misc$pars_interpolation[[v]]$eps = 0.1
#   stor$misc$pars_interpolation$psill =  stor$misc$pars_interpolation$psill / 2
#   stor$misc$pars_interpolation[[v]]$x$kp = stor$misc$pars_interpolation[[v]]$x$kp / 2
#   stor$misc$pars_interpolation[[v]]$y$kp = stor$misc$pars_interpolation[[v]]$y$kp / 2
#
# }

my_lm = function(date_selected, vname, obs=FALSE)
{
  # compute seasonal covariates
  idx = match(date_selected, dates_all_ini)
  X_j = my_season_X(date_selected)
  # X_j = data.frame(sin_j = my_season(format(date_selected, '%j')),
  #                  cos_j = my_season(format(date_selected, '%j'), s=pi/2),
  #                  sin_j2 = my_season(format(date_selected, '%j'), f=2),
  #                  cos_j2 = my_season(format(date_selected, '%j'), s=pi/2, f=2),
  #                  year = as.integer(format(date_selected, '%Y')))

  # compute seasonal effect
  betas_j = stor[['misc']][['lm_seasonal']][[vname]]
  effect_j = as.vector( betas_j[1] + ( as.matrix(X_j) %*% betas_j[-1] ) )

  # compute elevation effect using original station metadata or grid data
  betas_elev = stor[['misc']][['lm_elevation']][[vname]]
  if(obs) { effect_elev = as.vector( betas_elev[1] + ( as.matrix(X_new_obs) %*% betas_elev[-1] ) ) } else {
    effect_elev = as.vector( betas_elev[1] + ( as.matrix(X_new) %*% betas_elev[-1] ) )
  }



  effect_j + effect_elev

  #as.vector( lm_selected[1] + ( as.matrix(X_new) %*% lm_selected[-1] ) )

  #predict(lm_selected, newdata=X_new)

  # seasonal and year terms
  # yr = format(date_selected, '%Y') |> as.integer()
  # jday = format(date_selected, '%j') |> as.integer()
  # seasonal_sin_pred = stor[['misc']][['my_season']](jday)
  # seasonal_cos_pred = stor[['misc']][['my_season']](jday, s=pi/2)

  # compile into dataframe
  # X_new_out = X_new |> cbind(data.frame(year = rep(yr, n_dem),
  #                                       seasonal_sin = rep(seasonal_sin_pred, n_dem),
  #                                       seasonal_cos = rep(seasonal_cos_pred, n_dem)))

  # compute linear combination
  # predict(lm_current, newdata=X_new_out)
}

my_exp = function(x, vname, stor) exp(x + stor[['pkern']][['var_pkern']][[vname]][['gval']]/2) - log_const


### testing:

# #select and copy some data
# idx = 1972
#
# idx = 1374
#
#
# idx= idx+1
# vname = 'TMAX'
# data_mat = stor$pkern$data_mat[[vname]]
# date_selected = as.Date(colnames(data_mat)[idx])
# date_selected
#
#
# pars_test = stor$misc$pars_interpolation[[vname]]
# vec_selected = data_mat[,idx] |> as.numeric()
# krig_lm = my_lm(date_selected, vname)
#
# # omit stations not appearing in grid
# station_displayed = stations_sf[['station_num']] %in% unique(idx_grid)
#
# # copy station data to points object and grid
# stations_sf[['test_value']] = vec_selected
#
#
# residuals_selected = krig_lm - vec_selected[idx_grid]
# residuals_selected_mean = residuals_selected |> mean(na.rm=T)
#
# stations_test = stations_pkern |> modifyList(list(gval=residuals_selected - residuals_selected_mean))
# residuals_interpolated = pkern_cmean(stations_test, pars_test, X=0)
#
# # plot to verify
# krig_pred = krig_lm - (residuals_interpolated - residuals_selected_mean)


# if(vname=='LOG_PRCP')
# {
#   krig_pred = my_exp(krig_pred, vname, stor=stor)
#   stations_sf[['test_value']] = exp(stations_sf[['test_value']]) - log_const
# }
#
# snum_range = c(krig_pred, stations_sf[['test_value']][station_displayed]) |> range(na.rm=TRUE)
# my_breaks = seq(min(snum_range), max(snum_range), length.out=1e2)
#
# stations_pkern |> modifyList(list(gval=krig_pred)) |> pkern_plot(zlim = snum_range)
# stations_sf['test_value'][station_displayed, ] |> plot(cex=2, pal=my_pal, breaks=my_breaks, add=TRUE)
# stations_sf['test_value'][station_displayed, ] |> plot(pch=16, pal=my_pal, breaks=my_breaks, add=TRUE)





# colnames(stor$pkern$data_mat$PRCP)[1976]
# stor$pkern$data_mat$PRCP[, 1976] |> range(na.rm=TRUE)
# stor$pkern$data_mat$LOG_PRCP[, 1976] |> exp() |> range(na.rm=TRUE)


#stor[['pkern']][['var_pkern']][['LOG_PRCP']] |> pkern_plot()

# shiny UI
ui = fluidPage(

  # centered title
  titlePanel( h1('GHCND weather data interpolation with pkern', align='center') ),
  br(),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      radioButtons('vname',
                   label = 'weather variable',
                   choices = vnames,
                   selected = 'TMIN'),

      radioButtons('model_output',
                   label = 'display',
                   choices = c('kriging predictor', 'kriging predictor early', 'kriging variance'),
                   selected = 'kriging predictor'),

      sliderInput('zlim_slider', 'display range', min=0, max=10, value = c(1,2), ticks=FALSE),

      #uiOutput('zlim_slider_output'),
      checkboxInput('zlim_auto', label='auto reset', value=TRUE),
      actionButton('zlim_reset', label='reset display range', value=TRUE),

    ),

    mainPanel(
      width = 9,
      plotOutput('krig_plot_display', height='70vh')
    )
  ),

  # centered date slider input
  fluidRow(
    column(width = 1, ''),
    column(width = 10, sliderInput('date',
                                  label = 'Date',
                                  min = dates_all_ini[1],
                                  max = dates_all_ini[length(dates_all_ini)],
                                  value = as.Date('2022-05-27'),
                                  timeFormat = '%Y-%m-%d',
                                  width = '100%')),
    column(width = 1, '')
  )
)

# extra (add in later)

# DEM (pkern)
#stor[['pkern']][['dem_pkern']] |> pkern_plot()

# precomputed variance (pkern)
#stor[['pkern']][['var_pkern']] |> pkern_plot()

# response data matrix
#stor[['pkern']][['data_mat']] |> str()

# input = list(vname='PRCP', date=as.Date('2020-10-01'))

server = function(input, output, session)
{
  observeEvent(input[['zlim_reset']], {

    krig_result = run_krig()
    zlim_outer = krig_result[['zlim_outer']]
    zlim_snapped = krig_result[['zlim_snapped']]
    step_size = 10^(ceiling(log(diff(zlim_outer), base=10))-3)

    updateSliderInput(session,
                      'zlim_slider',
                      label = 'display range',
                      min = zlim_outer[1],
                      max = zlim_outer[2],
                      step = step_size,
                      value = zlim_snapped)

  })


  # extract all observed data for selected variable
  copy_data = reactive({

    # copy data and date range
    vname = input[['vname']]
    data_mat = stor[['pkern']][['data_mat']][[vname]]
    dates_all = as.Date( colnames(data_mat) )
    data_range = data_mat[station_displayed,] |> range(na.rm=TRUE)

    # return in list
    return( list(vname=vname,
                 data_mat=data_mat,
                 dates_all=dates_all,
                 data_range=data_range) )
  })

  #set_pars = reactiveValues(ml=)

  # kriging for a selected date
  run_krig = reactive({


    # reactive inputs
    input[['model_output']]
    date_selected = input[['date']]
    data_in = copy_data()

    # find column for input date
    vname = data_in[['vname']]
    dates_all = data_in[['dates_all']]
    idx = match(date_selected, dates_all)

    # return everything in a list
    if( input[['model_output']] == 'kriging variance' )
    {
      # second division by 10 (since variance squares the constant)
      var_grid = stor[['pkern']][['var_pkern']][[vname]][['gval']] / 1e2
      zlim_snapped = range(var_grid, na.rm=TRUE)
      zlim_outer = zlim_snapped
      step_size = 10^(ceiling(log(diff(zlim_outer), base=10))-3)

      # update plot limits automatically when auto reset is on
      if( input[['zlim_auto']] )
      {
        # force output to wait until this input value is settled (see below)
        freezeReactiveValue(input, 'zlim_slider')

        # since zlim_slider is frozen, the plot output is delayed until this change propagates
        updateSliderInput(session,
                          'zlim_slider',
                          label = 'display range',
                          min = zlim_outer[1],
                          max = zlim_outer[2],
                          step = step_size,
                          value = zlim_snapped)
      }

      return( list(idx = idx,
                   date_print = as.character(date_selected),
                   pts = NA,
                   krig_lm = NA,
                   krig_resid = NA,
                   krig_pred = var_grid,
                   zlim_outer = zlim_outer,
                   zlim_snapped = zlim_snapped) )
    }

    # linear model predictions for this variable
    krig_lm = my_lm(date_selected, vname, obs=FALSE)
    krig_lm_obs = my_lm(date_selected, vname, obs=TRUE)

    # observed data vector
    # pts = data_in[['data_mat']][,idx]
    # data_obs = pts[idx_reorder]
    station_obs = data_in[['data_mat']][,idx]
    station_res = station_obs - krig_lm_obs
    station_res_grid = station_res[idx_grid]

    #station_obs_grid = station_obs[idx_grid]


    # compute residuals and copy to grid
    #station_residuals_grid = station_obs_grid - krig_lm
    #residuals_mean = mean(station_residuals_grid, na.rm=TRUE)
    #residuals_pkern = stations_pkern |> modifyList(list(gval=station_residuals_grid-residuals_mean))
    #residuals_pkern = stations_pkern |> modifyList(list(gval=station_residuals_grid))
    #krig_lm - station_obs_grid

    res_pkern = stations_pkern |> modifyList(list(gval=station_res_grid))

    #stations_pkern[['gval']][] = NA
    #stations_pkern[['gval']][is_obs] = krig_lm[is_obs] - data_obs

    # interpolate the anomaly
    pars_interpolation = stor[['misc']][['pars_interpolation']][[vname]]
    if( input[['model_output']] == 'kriging predictor early'  ) pars_interpolation = stor[['misc']][['pars_interpolation_all']][[vname]]
    res_interpolated = pkern_cmean(res_pkern, pars_interpolation, X=0)

    # kriging prediction
    #krig = krig_lm - (residuals_interpolated + residuals_mean)
    krig = krig_lm + res_interpolated

    # apply inverse transformations
    if(vname=='LOG_PRCP') krig = my_exp(krig, vname, stor=stor)
    if(vname=='PRCP') krig[krig < 0] = 0
    krig_range = range(krig, na.rm=TRUE)

    # observed range of all data and of current date
    range_all = data_in[['data_range']]
    range_obs = station_obs[station_displayed] |> range(na.rm=TRUE)
    if(vname=='LOG_PRCP')
    {
      range_all = exp(range_all) - log_const
      range_obs = exp(range_obs) - log_const
      station_obs = exp(station_obs) - log_const
    }

    # tight plot limits
    zlim_outer = range(c(range_all, krig), na.rm=TRUE) / 10
    zlim_snapped = range(c(range_obs, krig), na.rm=TRUE) / 10
    if( diff(zlim_outer) == 0 ) zlim_outer[2] = 1 + zlim_outer[2]
    if( diff(zlim_snapped) == 0 )  zlim_snapped[2] = 1 + zlim_snapped[2]

    # looser limits, rounded for easier labeling
    zlim_outer_pretty = c(floor(zlim_outer)[1], ceiling(zlim_outer[2]))
    step_size = 10^(ceiling(log(diff(zlim_outer_pretty), base=10))-3)

    #zlim_snapped = c(floor(zlim_snapped[1]), ceiling(zlim_snapped[2]))

    # update plot limits automatically when auto reset is on
    if( input[['zlim_auto']] )
    {
      # force output to wait until this input value is settled (see below)
      freezeReactiveValue(input, 'zlim_slider')

      # since zlim_slider is frozen, the plot output is delayed until this change propagates
      updateSliderInput(session,
                        'zlim_slider',
                        label = 'display range',
                        min = zlim_outer_pretty[1],
                        max = zlim_outer_pretty[2],
                        step = step_size,
                        value = zlim_snapped)
    }

    # return everything in a list
    return( list(idx = idx,
                 date_print = as.character(date_selected),
                 station_obs = station_obs / 10,
                 krig_lm = krig_lm,
                 krig_resid = res_interpolated,
                 krig_pred = krig / 10,
                 zlim_outer = zlim_outer_pretty,
                 zlim_snapped = zlim_snapped) )
  })

  # render the plot
  output[['krig_plot_display']] = renderPlot({

    # read user input and compute kriging estimates
    zlim_plot = input[['zlim_slider']]
    date_selected = input[['date']]
    vname = input[['vname']]
    krig_result = run_krig()

    print(date_selected)

    # add some padding to prevent clipping at plot limits
    zlim_padding = 0.05 * diff(zlim_plot)
    zlim_plot = zlim_plot + c(-1,1)*zlim_padding

    # make a pkern object, rescale to ordinary units
    krig_plot = krig_result[['krig_pred']]
    krig_plot_pkern = stations_pkern |> modifyList(list(gval=krig_plot))

    # plot text
    vname_print = title_lookup[[vname]]
    if( input[['model_output']] == 'kriging variance' ) vname_print = paste(vname_print, 'variance')
    title_print = paste(krig_result[['date_print']], 'daily', vname_print)

    # draw the plot
    pkern_plot(krig_plot_pkern,
               main=title_print,
               zlab=ifelse(vname=='LOG_PRCP', 'PRCP', vname),
               zlim=zlim_plot,
               pal='Spectral',
               reset=FALSE)

    # create points object from stations observed
    if( input[['model_output']] == 'kriging variance' )
    {
      plot(st_geometry(stations_sf), cex=2, add=TRUE)

    } else {

      stations_sf['station_obs'] = krig_result[['station_obs']]
      plot_sf = stations_sf['station_obs']

      #plot_sf = stations_sf['obs_data'][idx_reorder,]
      is_displayed = station_displayed & !is.na(stations_sf[['station_obs']])
      plot_sf = plot_sf[is_displayed,]

      mybreaks = seq(min(zlim_plot), max(zlim_plot), by=diff(zlim_plot)/1e2)
      plot(st_geometry(plot_sf), cex=2, add=TRUE)
      plot(plot_sf, pch=16, cex=1, add=TRUE, pal=my_pal, breaks=mybreaks)
    }

    #my_plot_annotations(geo_list=geo_list, bw=TRUE)
    stor[['misc']][['uyrw_boundary']] |> plot(col=NA, add=TRUE)
    stor[['misc']][['state_boundaries']] |> plot(border='grey20', col=NA, add=TRUE)
  })

}

shinyApp(ui=ui, server=server)
