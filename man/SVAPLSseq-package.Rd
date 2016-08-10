\name{SVAPLSseq-package}
\alias{SVAPLSseq-package}
\alias{SVAPLSseq}
\docType{package}
\title{
SVAPLSseq: An R package to adjust for the hidden factors of variability in 
differential gene expression studies based on RNAseq data.
}
\description{
The package \code{SVAPLSseq} contains functions that are intended for the 
identification and correction of the hidden variability owing
to a variety of unknown subject/sample specific and technical effects 
of residual heterogeneity in an RNAseq gene expression data. 
}
\details{
\tabular{ll}{
Package: \tab SVAPLSseq\cr
Type: \tab Package\cr
License: \tab GPL-3\cr
}
The package can be used to find the genes that are truly differentially expressed
between two groups of samples from an RNAseq data, after adjusting for different hidden factors of
expression heterogeneity. The function \code{svplsSurr} operates on the raw data matrix 
of gene level read counts and extracts the signatures of the underlying hidden variability in 
the form of a set of surrogate variables. The function \code{svplsTest} detects the truly positive 
genes after correcting for the hidden signals (surrogate variables) extracted by \code{svplsSurr}.
}

\seealso{
\code{\link{svplsSurr}}, \code{\link{svplsTest}}
}

\author{
Sutirtha Chakraborty.

Maintainer: Sutirtha Chakraborty <statistuta@gmail.com>
}
\references{
Boulesteix, A-L. and Strimmer, K. Partial least squares: a versatile tool for the analysis of high-dimensional 
genomic data. Briefings in Bioinformatics 2007; \bold{8}(1):32--44.
}

\keyword{models}

\examples{
##Loading the Simulated Data
data(sim.dat)

## Fitting a linear model with the surrogate variables and detecting the differentially expressed genes
group = as.factor(c(rep(1, 10), rep(-1, 10)))
sv <- svplsSurr(dat = sim.dat, group = group)
surr = surr(sv)
fit <- svplsTest(dat = sim.dat, group = group, surr = surr, test = "Wald")
head(sig.genes(fit))
}
