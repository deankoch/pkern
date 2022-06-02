# pkern
R package for modeling Gaussian processes on grids

# Update June 2022

After some months of testing new ideas for improving this package I've finally committed a new version.
This is a major update since the last commit - code files are reorganized more sensibly, and almost every
function has been rewritten and renamed.

Nearly all of the package's functions now include extensive test code and examples, and I will be working on
a new set of vignettes until around mid-summer. 

# Overview

`pkern` provides a computationally lean implementation of a 2-dimensional spatial correlation model for
gridded data. This can be useful when working with geo-referenced data, such as in earth sciences, where 
users often need a tool for interpolation or down-scaling, but want to steer clear of bias problems introduced
with Voronoi tiles (ie nearest neighbour) or inverse distance weighting.

More generally, `pkern` offers an fast and simple back-end for modeling with spatially correlated errors.
It works much faster than alternatives like `gstat`, at the expense of restricting the type of model
users can select. However I think for most users, the default model in `pkern` (a Gaussian process with a
Gaussian covariance kernel) is more than adequate in many situations.

`pkern` supports raster and geometry inputs from `sf` and `terra`, as well as simpler matrix and vector inputs.
These packages are suggested, but not required. `pkern` is written using only base dependencies (included by
default in R) like `graphics` and `stats`. 


# Technical Features

* models anisotropic Gaussian processes on 2-dimensional regular grids for a choice of covariance kernels
* fast computation of the likelihood function, generalized least squares, and kriging predictor/variance
* automated maximum likelihood model fitting with sample semi-variogram diagnostic
* user friendly helper functions for raster down-scaling and point interpolation

# Technical Background

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

