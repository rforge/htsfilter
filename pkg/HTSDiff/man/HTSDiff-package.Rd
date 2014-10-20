\name{HTSDiff-package}
\alias{HTSDiff-package}
\docType{package}
\title{
Differential analysis of RNA-seq data with Poisson mixture models
}
\description{
This package implements a Poisson mixture model to identify differentially expressed 
genes from RNA-seq data.
}
\details{
\tabular{ll}{
Package: \tab HTSDiff\cr
Type: \tab Package\cr
Version: \tab 1.0.5\cr
Date: \tab 2014-10-20\cr
License: \tab GPL (>=3)\cr
LazyLoad: \tab yes\cr
}
}
\author{
Andrea Rau, Marie-Laure Martin-Magniette, Cathy Maugis-Rabusseau

Maintainer: Andrea Rau <\url{andrea.rau@jouy.inra.fr}>
}
\references{
S. Balzergue, G. Rigaill, V. Brunaud, E. Blondet, A. Rau, O. Rogier, J. Caius, C. Maugis-Rabusseau,
L. Soubigou-Taconnat, S. Aubourg, C. Lurin, E. Delannoy, and M.-L. Martin-Magniette. (2013) 
HTSDiff: a sensitive analysis of differential gene expression by RNA-seq (submitted).
}
\keyword{ models }
\examples{

set.seed(12345)

## Generate synthetic data: 2000 genes under H0
test <- syntheticData(H0number = 2000)

## Mixture model differential analysis
## DEtest <- HTSDiff(test, c(1,1,2,2))


}
