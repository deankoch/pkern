% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pkern_estim.R
\name{pkern_fit}
\alias{pkern_fit}
\title{Fit a covariance model to the data by maximum likelihood}
\usage{
pkern_fit(
  g_obs,
  pars = NULL,
  X = NA,
  iso = TRUE,
  initial = NULL,
  quiet = FALSE,
  lower = NULL,
  upper = NULL,
  n_max = 1000
)
}
\arguments{
\item{g_obs}{todo}

\item{pars}{todo}

\item{X}{todo}
}
\description{
An automated model fitting procedure
}
\details{
documentation unfinished
}
\examples{

# define a grid
gdim = c(50, 53)
g_empty = pkern_grid(gdim)
pars_src = pkern_pars(g_empty)
pars_src = modifyList(pars_src, list(eps=runif(1, 0, 1e1), psill=runif(1, 0, 1e2)))
pars_src[['y']][['kp']] = pars_src[['x']][['kp']] = runif(1, 1, 50)

# generate example data and fit to it
gval = pkern_sim(g_empty, pars_src, quiet=F)
g_obs = modifyList(g_empty, list(gval=gval))
pkern_plot(g_obs)
fit_result = pkern_fit(g_obs, pars='gau')

fit_result$pars |> pkern_pars_update()
pars_src |> pkern_pars_update()

# check sequence of other psill values
pars_out = fit_result$pars
psill_test = ( 2^(seq(5) - 3) ) * pars_out[['psill']]
LL_test = sapply(psill_test, function(s) pkern_LL(modifyList(pars_out, list(psill=s)), g_obs) )
plot(psill_test, LL_test)
lines(psill_test, LL_test)
print(data.frame(psill=psill_test, likelihood=LL_test))

# repeat with most data missing
n = prod(gdim)
n_obs = 200
gval = pkern_sim(g_empty, pars_src, quiet=TRUE)
g_obs = modifyList(g_empty, list(gval=gval))
idx_obs = sample.int(prod(gdim), n_obs)
g_miss = modifyList(g_obs, list(gval=rep(NA, n)))
g_miss$gval[idx_obs] = g_obs$gval[idx_obs]
pkern_plot(g_miss)


}
