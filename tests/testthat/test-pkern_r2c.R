test_that('default row to column vectorization functionality works', {

  dims = c(12, 13)
  x = pkern_r2c(dims, in.byrow=TRUE, out.byrow=FALSE) |> matrix(nrow=dims[1], byrow=TRUE) |> as.vector()
  expect_equal(x, seq(prod(dims)))

})
