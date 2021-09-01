# pkern
R package for fast kriging on grids using product kernels
(work in progress)

This is an R implementation of some methods I developed in [my thesis](https://doi.org/10.7939/r3-91zn-v276)
for speeding up geostatistical computations involving large covariance matrices. The central idea is to model
spatial dependence using a separable 2-dimensional covariance kernel, defined as the product of two (1-dimensional)
univariate covariance kernels. This introduces special symmetries and structure in the covariance matrix, which are
exploited in this package for fast and memory-efficient computations.

This package will accompany a paper on fast and efficient downscaling of weather grids so the focus is on a particular
application of kriging, but the underlying methods are applicable much more generally. See [1(https://doi.org/10.7939/r3-g6qb-bq70)],
where I use product kernels to study directions of anisotropy in a nonstationary random fields, and
[2](https://doi.org/10.1007/s11538-021-00899-z), [3](https://doi.org/10.1098/rsif.2020.0434), where I apply it to fit a
covariance structure, and to speed up calculations of dispersal kernel convolutions.

