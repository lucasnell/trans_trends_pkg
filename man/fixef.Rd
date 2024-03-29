% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/armmMod_class.R
\docType{methods}
\name{fixef}
\alias{fixef}
\alias{fixef.armmMod}
\alias{fixed.effects}
\title{Extract fixed-effects estimates from an \code{armmMod} object}
\usage{
\method{fixef}{armmMod}(object, se_method = c("quantile", "hpdi"), ...)
}
\arguments{
\item{object}{an object of class \code{armmMod}, a result of a call to \code{armm}.}

\item{se_method}{A single string, for either using quantiles (\code{"quantile"})
or HPDI (\code{"hpdi"}) to compute the standard errors.
Defaults to \code{"quantile"}.}

\item{...}{Not used.}
}
\value{
A data frame of fixed-effects estimates and standard errors.
}
\description{
Extract fixed-effects estimates from an \code{armmMod} object
}
