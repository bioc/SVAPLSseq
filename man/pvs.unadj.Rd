% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/Methods.R
\docType{methods}
\name{pvs.unadj}
\alias{pvs.unadj}
\alias{pvs.unadj}
\alias{pvs.unadj.svplsTest}
\alias{pvs.unadj,svplsTest-method}
\title{Accessor for the 'pvs.unadj' slot of a 'svplsTest' object}
\usage{
pvs.unadj(object)

\S4method{pvs.unadj}{svplsTest}(object)
}
\arguments{
\item{object}{a \code{svplsTest} object}
}
\description{
Accessor for the 'pvs.unadj' slot of a 'svplsTest' object
}
\examples{
data(sim.dat)
group = as.factor(c(rep(1, 10), rep(-1, 10)))
sv = svplsSurr(sim.dat, group, surr.select = "automatic")
surr = surr(sv)
fit = svplsTest(dat = sim.dat, group = group, surr = surr, normalization = "TMM", test = "t-test")
pvs.unadj(fit)
}
