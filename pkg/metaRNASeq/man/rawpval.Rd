\name{rawpval}
\alias{rawpval}
\docType{data}
\title{
Simulated p-values
}
\description{
The p-values provided here result from the following procedure:
 1) simulation of two RNA-seq experiments with four replicates in each condition via the \code{\link{sim.function}}, 2) analysis of differentially expressed tags using the DESeq package.
}
\usage{data(rawpval)}
\format{
  The format is:
List of 2
 $ : num [1:14456] 8.31e-01 9.08e-01 5.86e-01 9.37e-09 2.87e-01 ...
 $ : num [1:14456] 4.42e-01 3.85e-01 4.89e-01 1.80e-71 6.99e-01 ...
}
\details{
It is possible to reproduce these p-values using the procedure described in the package vignette.
}
\references{
A. Rau, G. Marot and F. Jaffrezic (2012). Differential meta-analysis of RNA-seq data. (Manuscript in progress)
}
\examples{
data(rawpval)
## Maybe str(rawpval) ; plot(rawpval) ...
}
\keyword{datasets}
