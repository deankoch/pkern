test_that('RasterLayer is loaded correctly', {

  # load raster package and print error on failure
  require(raster)
  pkern_checkRaster()

  # load an example raster and extract vectorized data
  r.in = system.file('external/rlogo.grd', package='raster') |> raster(band=1)
  pkern.in = pkern_fromRaster(r.in)
  pkern.in$data |> matrix(ncol=pkern.in$dims[1]) |> graphics::image()
  pkern_toRaster(pkern.in, template=r.in)

  r.out = pkern_toRaster(pkern_fromRaster(r.in, omode='vector'), template=r.in)
  expect_equal(dim(r.out), dim(r.in))
  expect_equal(raster::values(r.out), raster::values(r.in))

})
