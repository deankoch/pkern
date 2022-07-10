library(shiny)

# load pkern
library(devtools)
load_all()

#+ dependencies

# extra dependencies for loading and preparing spatial data
library(sf)
library(data.table)

# path to the output file written by GHCN_vignette_preprocessing.R
rds_data_file = 'D:/pkern/vignettes/GHCN/GHCN_preprocessing_data.rds'

# load the data
stor = readRDS(rds_data_file)
vnames = names(stor[['pkern']][['data_mat']])

# copy some of the small loaded objects
title_lookup = stor[['misc']][['title_lookup']]
idx_map = stor[['misc']][['idx_map']]
stations_sf = stor[['misc']][['stations_sf']]

# layer with station mapping
#stations_pkern |> str()
# or: stations_pkern = dem_pkern |> modifyList(list(gval=idx_map))

# grid elevation values
x_dem = stor[['pkern']][['dem_pkern']][['gval']]
n_dem = length(x_dem)

# derive mapping keys on grid
is_obs = !is.na(idx_map)
idx_reorder = idx_map[is_obs] # station numbers in grid vectorized order
stations_pkern = stor[['pkern']][['dem_pkern']] |> modifyList(list(gval=idx_map))

# plotting parameters
my_pal = function(x) hcl.colors(x, 'Spectral', rev=T)
my_col = function(x) my_pal(1e2)[ num_to_cent(x) ]
dates_all_ini = as.Date(colnames(stor[['pkern']][['data_mat']][['TMIN']]))

# shiny UI
ui = fluidPage(

  # centered title
  titlePanel( h1('GHCND weather data interpolation with pkern', align='center') ),
  br(),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      radioButtons('vname',
                   label = 'GHCND variable',
                   choices = vnames,
                   selected = 'TMIN'),

      checkboxInput('zlim_auto', label='automatic plot limits', value=TRUE),
      uiOutput('zlim_slider_output')


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
                                  value = dates_all_ini[1],
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

# input = list(vname='', date=as.Date('2020-10-01'))

server = function(input, output)
{

  # plot limit slider
  output[['zlim_slider_output']] = renderUI({

    zlim_selected = input[['zlim_slider_input']]
    if( is.null(zlim_selected) )  zlim_selected = zlim_inner()

    if( !input[['zlim_auto']] ) {

      sliderInput('zlim_slider_input',
                  label = 'display range',
                  min = zlim_outer()[1],
                  max = zlim_outer()[2],
                  value = zlim_selected)
    }
  })

  # define variable names and times
  vname = reactive({ input[['vname']] } )
  vname_print = reactive({ title_lookup[['vname_print']][title_lookup[['vname']] == vname()] })
  dates_all = reactive({ as.Date(colnames(stor[['pkern']][['data_mat']][[vname()]])) })
  idx = reactive({ match(input[['date']], dates_all()) })
  date_i = reactive({ dates_all()[idx()] })

  # extract observed data
  data_mat = reactive({ stor[['pkern']][['data_mat']][[vname()]] })

  # observed range all time
  zlim_outer = reactive({
    obs_range = range(data_mat(), na.rm=TRUE)
    if(vname()=='LOG_PRCP') obs_range = ( exp(obs_range) - 1 )
    c( floor(obs_range[1]/10), ceiling(obs_range[2]/10) )
  })

  # observed range on specific date
  zlim_inner = reactive({
    obs_range = range(data_mat()[, idx()], na.rm=TRUE)
    if(vname()=='LOG_PRCP') obs_range = ( exp(obs_range) - 1 )
    c( floor(obs_range[1]/10), ceiling(obs_range[2]/10) )
  })

  # plot limits
  zlim = reactive({
    zlim_min = min(c(zlim_inner(), krig_plot()), na.rm=T)
    zlim_max = max(c(zlim_inner(), krig_plot()), na.rm=T)
    zlim = c(zlim_min, zlim_max)
  })

  # extract and transform source point data
  krig_obs = reactive({
    ptval = stor[['pkern']][['data_mat']][[vname()]][,idx()]
    if(vname()=='LOG_PRCP') ptval = ( exp(ptval) - 1 )
    ptval / 10
  })

  # linear predictor computer
  krig_lm = reactive({

    yr = format(date_i(), '%Y') |> as.integer()
    jday = format(date_i(), '%j') |> as.integer()
    seasonal_sin_pred = stor[['misc']][['my_season']](jday)
    seasonal_cos_pred = stor[['misc']][['my_season']](jday, s=pi/2)

    # combine all predictors
    X_new = data.frame(elevation = x_dem,
                       log_elevation = log(x_dem),
                       year = rep(yr, n_dem),
                       seasonal_sin = rep(seasonal_sin_pred, n_dem),
                       seasonal_cos = rep(seasonal_cos_pred, n_dem))


    # compute their linear combination
    predict(stor[['misc']][['lm_elevation']][[vname()]], newdata=X_new)

  })

  # kriging computer
  krig_plot = reactive({

    # compute all available linear model residuals
    y_residual = krig_lm()[is_obs] - data_mat()[,idx()][idx_reorder]
    stations_pkern[['gval']][is_obs] = y_residual

    # interpolate the anomaly
    pars_interpolation = stor[['misc']][['pars_interpolation']][[vname()]]
    fit_result_cmean = pkern_cmean(stations_pkern, pars_interpolation, X=0)
    zpred_pkern = stations_pkern |> modifyList(list(gval=fit_result_cmean))

    # kriging prediction
    krig = krig_lm() - fit_result_cmean

    # apply inverse transformations
    if(vname()=='LOG_PRCP') krig = (exp(krig + stor[['pkern']][['var_pkern']][[vname()]][['gval']]/2) - 1)
    if(vname()=='PRCP') krig[krig < 0] = 0

    # change units from 10ths
    krig / 10
  })

  #
  output[['krig_plot_display']] = renderPlot({

      # make a pkern object
      krig_plot_pkern = stations_pkern |> modifyList(list(gval=krig_plot()))
      #pkern_plot(krig_plot_pkern, zlab='precip (mm)')

      # create points object from station observed
      stations_sf[['obs_data']] = krig_obs()
      plot_sf = stations_sf['obs_data'][idx_reorder,]
      plot_sf = plot_sf[!is.na(plot_sf$obs_data),]
      title_print = paste(as.character(date_i()), 'daily', vname_print())

      zlim_plot = range(c(krig_obs(), krig_plot()), na.rm=TRUE)
      if( !input[['zlim_auto']] ) zlim_plot = input[['zlim_slider_input']]
      if( is.null(zlim_plot) ) return(NULL)

      #zlim = range(plot_sf[['obs_data']], na.rm=T) + c(-5,5)
      pkern_plot(krig_plot_pkern,
                 main=title_print,
                 zlab=ifelse(vname()=='LOG_PRCP', 'PRCP', vname()),
                 zlim=zlim_plot,
                 reset=FALSE)

      mybreaks = seq(min(zlim_plot), max(zlim_plot), by=diff(zlim_plot)/1e2)
      #pkern_plot(krig_plot_pkern, zlab='precip (mm)', reset=FALSE)
      plot(st_geometry(plot_sf), cex=2, add=TRUE)
      plot(plot_sf, pch=16, cex=1, add=TRUE, pal=my_pal, breaks=mybreaks)

      #my_plot_annotations(geo_list=geo_list, bw=TRUE)
      stor[['misc']][['uyrw_boundary']] |> plot(col=NA, add=TRUE)
      stor[['misc']][['state_boundaries']] |> plot(border='grey20', col=NA, add=TRUE)

    })

}

shinyApp(ui=ui, server=server)
