% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mJAM_SuSiE_helper.R
\name{susie_get_posterior_mean_v2}
\alias{susie_get_posterior_mean_v2}
\title{Get SuSiE posterior mean}
\usage{
susie_get_posterior_mean_v2(res, prior_tol = 1e-09)
}
\arguments{
\item{res}{A SuSiE fit object}

\item{prior_tol}{When the prior variance is estimated, compare the estimated value to prior_tol at the end of the computation, and exclude a single effect from PIP computation if the estimated prior variance is smaller than this tolerance value.}
}
\value{
A vector of posterior mean effects
}
\description{
Get SuSiE posterior mean
}
