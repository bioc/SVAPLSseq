% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/Methods.R
\docType{methods}
\name{sig.features}
\alias{sig.features}
\alias{sig.features}
\alias{sig.features.svplsTest}
\alias{sig.features,svplsTest-method}
\title{Accessor for the 'sig.features' slot of a 'svplsTest' object}
\usage{
sig.features(object)

\S4method{sig.features}{svplsTest}(object)
}
\arguments{
\item{object}{a \code{svplsTest} object}
}
\description{
Accessor for the 'sig.features' slot of a 'svplsTest' object
}
\examples{
data(sim.dat)
group = as.factor(c(rep(1, 10), rep(-1, 10)))
sv = svplsSurr(sim.dat, group, surr.select = "automatic")
surr = surr(sv)
fit = svplsTest(dat = sim.dat, group = group, surr = surr, normalization = "TMM", test = "t-test")
sig.features(fit)
}
