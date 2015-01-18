\name{HTSDiff-package}
\alias{HTSDiff-package}
\docType{package}
\title{
Differential Analysis of RNA-Seq Data with Poisson Mixture Models
}
\description{
A Poisson mixture model is implemented to identify differentially expressed 
genes from RNA-seq data.
}
\details{
\tabular{ll}{
Package: \tab HTSDiff\cr
Type: \tab Package\cr
Version: \tab 1.0.6\cr
Date: \tab 2015-01-18\cr
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
L. Soubigou-Taconnat, S. Aubourg, C. Lurin, E. Delannoy, and M.-L. Martin-Magniette. (2014) 
HTSDiff: A Model-Based Clustering Alternative to Test-Based Methods in Differential Gene Expression 
Analyses by RNA-Seq Benchmarked on Real and Synthetic Datasets (submitted).
}
\keyword{ models }
\examples{

set.seed(12345)

## Generate synthetic data: 2000 genes under H0
test <- syntheticData(H0number = 2000)

## Mixture model differential analysis
## DEtest <- HTSDiff(test, c(1,1,2,2))


}
