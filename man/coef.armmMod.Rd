% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/armmMod_class.R
\name{coef.armmMod}
\alias{coef.armmMod}
\title{Coefficients from an \code{armmMod} object}
\usage{
\method{coef}{armmMod}(object, se_method = c("quantile", "hpdi"), ...)
}
\arguments{
\item{object}{an object of class \code{armmMod}, a result of a call to \code{armm}.}

\item{se_method}{A single string, for either using quantiles (\code{"quantile"})
or HPDI (\code{"hpdi"}) to compute the standard errors.
Defaults to \code{"quantile"}.}

\item{...}{Not used.}
}
\value{
A data frame of fixed and (if present) random effects estimates
with standard errors.
}
\description{
Coefficients from an \code{armmMod} object
}
