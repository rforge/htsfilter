\name{syntheticData}
\alias{syntheticData}
\title{
Generate synthetic data with a pre-specified number or proportion of genes under the null hypothesis
(non-differential expression)
}
\description{
This function makes use of the initialDataset included as data in this package to generate synthetic data
containing a pre-specified number or proportion of genes under the null hypothesis.
}
\usage{
syntheticData(H0number, plot = FALSE, plot.name = NA)
}
\arguments{
  \item{H0number}{
Number of genes (or alternatively, proportion of genes) to be generated under the null hypothesis.
}
  \item{plot}{
\code{TRUE} if the user wishes to visualize the synthetic data as a scatterplot of the log2 mean counts for each
tissue, with null hypothesis (H0), differentially expressed (DE), and non-differentially expressed (NDE) genes
highlighted.
}
  \item{plot.name}{
If a plot is to be generated (\code{plot = TRUE}) and the user wishes to save a PDF file to the current working
directory, the plot name to be used.
}
}
\details{
To generate synthetic data, we follow these two steps:

1) For genes validated as differentially expressed via qPCR, we retain the original counts for flower buds and leaves, respectively.

2) Among the remaining genes (identified as non-differentially expressed via qPCR or not validated via qPCR), we randomly sample
\code{H0number} (or \code{H0number} percent) genes. For these selected genes, we replace the read counts for the flower buds with those
for the leaves, and we replace the read counts for the leaves with those from the additional replicates performed for leaves. As such,
the selected genes have read counts arising from four replicates of leaves (rather than two replicates of flower buds and two replicates of
leaves, as is the case for the remaining genes); these genes thus represent true null hypothesis (non-differentially expressed) genes.

}
\value{
A \code{data.frame} of dimension 28,094 genes x 4 columns (with gene IDs as row names): read counts for two replicates of flower buds (BF1 and BF2) and
read counts for two replicates of leaves (F1 and F2). Note that the 251 genes validated as differentially expressed via qPCR have been marked in the data
by the addition of ".DE" to the end of their gene ID, the 81 validated as non-differentially expressed via
qPCR have been marked by an additional ".NDE" at the end of their gene ID, and the \code{H0number} (or \code{H0number} percent) genes
generated under the null hypothesis have been marked by an addition ".Hnull" at the end of their gene ID. 
}

\author{
Andrea Rau <\url{andrea.rau@jouy.inra.fr}> and Marie-Laure Martin-Magniette <\url{marie_laure.martin-magniette@agroparistech.fr}>
}

\examples{

set.seed(12345)

## Generate synthetic data: 2000 genes under H0
syn <- syntheticData(H0number = 2000)

## Generate synthetic data: 60% genes under H0
syn <- syntheticData(H0number= 0.60)

}
\keyword{ models }

