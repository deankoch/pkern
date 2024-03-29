% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pkern_plot.R
\name{pkern_plot_semi}
\alias{pkern_plot_semi}
\title{Plot a semi-variogram}
\usage{
pkern_plot_semi(vg, pars = NULL, add = FALSE, fun = "classical", ...)
}
\arguments{
\item{vg}{data frame of sample absolute differences, with columns 'dabs', 'd' (and 'bin')}

\item{pars}{list of the form returned by \code{pkern_pars} with entries 'y', 'x', 'eps', 'psill'}

\item{add}{logical, indicates to draw on an existing plot rather than create a new one}

\item{fun}{character or function, the aggregation function (see details)}

\item{...}{further plotting parameters (see below)}
}
\value{
nothing
}
\description{
Plots a sample semi-variogram using the point pair difference data in \code{vg}.
Binned summary statistics are drawn as circles with size scaled to the sample sizes.
A covariance model (\code{pars}) is optionally drawn over the sample data as a ribbon plot.
}
\details{
If \code{vg} is a data frame, it should contain absolute differences (numeric 'dabs'),
inter-point distances (numeric 'd'), and, optionally, an assignment into distance bins
(integer 'bin') for a sample of point pairs. If 'bin' is missing, the function calls
\code{pkern_add_bins} to assign them automatically.

Function \code{fun} is the statistic to use for estimating the variogram (ie twice the
semi-variogram) from the distance-binned absolute differences in \code{vg}. If \code{fun} is a
function, it must accept sub-vectors of the numeric \code{vg$dabs} as its only argument,
returning a non-negative numeric scalar. \code{fun} can also be set to one of the names
'root_median', 'root_mean' (the default), or 'classical', as shorthand for the robust
fourth-root-based methods in section 2.4 of Cressie (1993), or the classical mean of
squares method of Matheron.

Optional list \code{pars} defines a theoretical semi-variogram to draw over the sample data.
When \code{add=TRUE}, the function overlays it on an existing plot (without changing the
legend, plot limits etc). Anisotropic models, which may assume a range of semi-variances
for any given distance, are drawn as a ribbon plot.

\code{add=TRUE} can only be used in combination with an earlier call to \code{pkern_plot_semi}
where \code{reset=FALSE} (which allows the function to change R's graphical parameters)

\code{vg} can be a grid object (anything understood by \code{pkern_grid}) rather than a
variogram data frame. When \code{add=FALSE}, the function uses it to set the distance limits
for an initial empty plot (the model semi-variance is then drawn if \code{pars} is supplied).
}
\section{Plotting parameters}{


The following style parameters are optional:

\describe{

\item{alpha_bin, alpha_model, alpha_bin_b, alpha_model_b}{numeric in \link{0,1}: respectively,
the transparency of the fill color for circles and ribbons (default 0.3), and their borders
(default 0.5 )}

\item{bty}{character: plot frame border type, passed to base::plot (default 'n' for no border)}

\item{col_bin, col_model}{character: respectively, the color to use for circles (default
'black') and ribbons (default 'blue')}

\item{cex_bin}{numeric > 0: scaling factor for circles (default 1.5)}

\item{d_max}{numeric > 0: x (distance) limit for plotting in input units}

\item{leg_main}{character: title for the sample bin legend (default 'model')}

\item{main}{character: a title}

\item{n_bin, n_test}{integer: respectively, the number of distance bins for the sample
(optional if \code{vg} has a 'bin' column, and ignored if \code{vg} is a grid object), and the
number of distances at which to evaluate the semi-variogram for model \code{pars} (default 5e3,
ignored if \code{pars} not supplied)}

\item{reset}{logical: indicates to reset graphical parameters to their original values
when finished (default TRUE)}

\item{unit_in, unit_out}{character: udunits2 strings specifying input distance units and
the desired units for distances in the plot (default for both is meters 'm')}

\item{xlab, ylab}{character: titles for the y and x axes. The default for y is
'semi-variance (gamma)', and for x 'distance (<unit_out>)'}

}
}

\examples{

# make example grid and reference covariance model
gdim = c(10, 15)
g_obs = pkern_grid(gdim)
pars = pkern_pars(g_obs, 'mat')

# plot a semivariance model
pkern_plot_semi(g_obs)
pkern_plot_semi(g_obs, pars)

# change annotations, sharpen ribbon border
pkern_plot_semi(g_obs, pars, main='title', xlab='x', ylab='y')
pkern_plot_semi(g_obs, pars, alpha_model_b=1, main='example title', xlab='x', ylab='y')

# input and output units are 'm' by default
pkern_plot_semi(g_obs, pars, unit_out='km')
pkern_plot_semi(g_obs, pars, unit_in='km', unit_out='km')

# generate sample data and sample semivariogram
g_obs$gval = pkern_sim(g_obs, pars)
vg = pkern_sample_vg(g_obs)
pkern_plot_semi(vg)

# different aggregation methods
pkern_plot_semi(vg, fun='root_median')
pkern_plot_semi(vg, fun='root_mean')
pkern_plot_semi(vg, fun='classical') # default
pkern_plot_semi(vg, fun=function(x) mean(x^2)) # same as classical

# plot again with reference model and adjust distance limits, number of bins
pkern_plot_semi(vg, pars)
pkern_plot_semi(vg, pars, d_max=10)
pkern_plot_semi(vg, pars, d_max=10, n_bin=1e2)

# add dashed line for sample variance (this tends to underestimate the sill)
pkern_plot_semi(vg, pars)
sample_var = var(g_obs[['gval']], na.rm=TRUE)
abline(h=sample_var, lty='dashed')

# initial call with reset=FALSE, then use add=TRUE to overlay the same model with a green fill
pkern_plot_semi(vg, pars, reset=FALSE)
pkern_plot_semi(vg, pars, add=TRUE, col_model='green', alpha_model_b=0)

# overlay several models with varying nugget effect
pars_vary = pars
for(i in seq(3))
{
  pars_vary$eps = 0.8 * pars_vary$eps
  pkern_plot_semi(vg, pars_vary, add=TRUE, alpha_model_b=0)
}
dev.off()

}
