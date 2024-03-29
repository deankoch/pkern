% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pkern_var.R
\name{pkern_corr}
\alias{pkern_corr}
\title{Stationary 1D correlation kernels}
\usage{
pkern_corr(pars, d = NA)
}
\arguments{
\item{pars}{list with elements 'k', the kernel name, and 'kp' the parameter vector}

\item{d}{numeric vector of length n, the distances to evaluate}
}
\value{
length-n vector or a list of parameters and bounds (see details)
}
\description{
Computes stationary correlation function values for the n (non-negative) 1-dimensional
distances in \code{d}. Parameter list entry \code{pars$kp} supplies the kernel parameter(s).
}
\details{
\code{pars$k} must be one of the following kernel names:
\itemize{
\item 'exp': exponential (special case of 'gex' with shape p=1)
\item 'gau': gaussian/stable (special case of 'gex' with shape p=2)
\item 'sph': spherical (AKA stable/Gaussian for p=2)
\item 'gex': gamma-exponential (with shape p)
\item 'mat': Whittle-Matern (Handcock and Wallis parameterization, with shape kap)
}

where the first three kernels have only a range parameters, and the last two have both a
range and shape parameter.

For the 1-parameter kernels, \code{pars$kp} is the range parameter value ('rho'); For the
2-parameter kernels, \code{pars$kp} is a vector whose first element is 'rho', and second
element is the shape parameter ('p' or 'kap'). The names in \code{pars$kp} are ignored and
only the order matters - the range parameter always comes first.

Note that this function will not accept parameter lists \code{pars} of the form returned by
\code{pkern_pars(...)} etc, as these include a pair of 1d kernels (however the sub-lists
\code{pars$y} and \code{pars$x} are accepted).
}
\examples{

# define test distances, grid, and example kernel
n_test = 100
d_test = seq(n_test)-1
g_example = pkern_grid(n_test)
pars = pkern_pars(g_example, c('mat', 'gau'))
pars_x = pars[['x']]

# compute and plot the x component of the correlogram function
corr_x_example = pkern_corr(pars_x, d=d_test)
plot(d_test, corr_x_example, pch=NA)
lines(d_test, corr_x_example)

## show how this function gets used to build more complicated objects

# get the other component correlation, take product
pars_y = pars[['y']]
corr_y_example = pkern_corr(pars_y, d=d_test)
corr_example = corr_y_example * corr_x_example

# variogram
variogram_example = pkern_vario_fun(pars, d=list(y=d_test, x=d_test))
variogram_compare = 2 * pars$eps + pars$psill * (1 - corr_example)
max(abs( variogram_example - variogram_compare ))

# Toeplitz component matrices built entirely from these correlation vectors
variance_matrix_example = pkern_var(g_example, pars)
max(abs( variance_matrix_example[['y']][,1L] - corr_y_example ))
max(abs( variance_matrix_example[['x']][,1L] - corr_x_example ))


}
