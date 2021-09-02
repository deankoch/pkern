#' ---
#' title: "pkern_demo"
#' author: "Dean Koch"
#' date: "`r format(Sys.Date())`"
#' output: github_document
#' ---
#'
#' **Mitacs UYRW project**
#'
#' **pkern_demo**: code snippets to use in a vignette
#'

library(devtools)
load_all()


library(here)
library(raster)
library(sf)

# load some raster data
path.testdata = 'D:/pkern/inst/testdata'

# load the raster and geojson
dem = raster(file.path(path.testdata, 'dem_big_creek.tif'))
tmin = st_read(file.path(path.testdata, 'tmin_big_creek.geojson'))

#


