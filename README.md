# pkern
R package for fast kriging on grids using product kernels

# Update June 2022

After some months of testing new ideas for improving this package I've finally committed a new version.
This is a major update since the last commit - code files are reorganized more sensibly, and almost every
function has been rewritten and renamed.

Nearly all of the package's functions now include extensive test code and examples, and I am working on
a new set of vignettes now. These will demonstrate new functionality for fast universal kriging prediction
and variance estimation.

# Features

* models Gaussian processes on 2-dimensional regular grids
* supports raster and geometry inputs from `sf` and `terra` 
* fast computation of the likelihood function, sample semo-variogram, and kriging predictor
* automated maximum likelihood model fitting
* user friendly helper functions for raster down-scaling and point interpolation


# Background

This is an R implementation of some methods I developed in [my thesis](https://doi.org/10.7939/r3-91zn-v276)
for speeding up geostatistical computations involving large covariance matrices. The central idea is to model
spatial dependence using a separable 2-dimensional covariance kernel, defined as the product of two (1-dimensional)
univariate covariance kernels. This introduces special symmetries and structure in the covariance matrix, which are
exploited in this package for fast and memory-efficient computations.

This package will accompany a paper on fast and efficient downscaling of weather grids so the focus is on a particular
application of kriging, but the underlying methods are applicable much more generally. See [[1](https://doi.org/10.7939/r3-g6qb-bq70)],
where I use product kernels to study directions of anisotropy in a nonstationary random fields, and
[[2](https://doi.org/10.1007/s11538-021-00899-z), [3](https://doi.org/10.1098/rsif.2020.0434)], where I apply it to
fit acovariance structure, and to speed up calculations of dispersal kernel convolutions.

