\name{HTSFilter-package}
\alias{HTSFilter-package}
\docType{package}
\title{
Filter replicated high-throughput sequencing expression data
}
\description{
This package implements a filtering procedure for replicated
high-throughput sequencing data based on a global Jaccard similarity index in 
order to identify genes with low, constant levels of expression across one or 
more experimental conditions.
}
\details{
\tabular{ll}{
Package: \tab HTSFilter\cr
Type: \tab Package\cr
Version: \tab 0.1.0\cr
Date: \tab 2012-10-16\cr
License: \tab Artistic-2.0\cr
}
}
\author{
Andrea Rau, Melina Gallopin, Gilles Celeux, and Florence Jaffrezic

Maintainer: Andrea Rau <andrea.rau@jouy.inra.fr>
}
\references{

R. Bourgon, R. Gentleman, and W. Huber. (2010) Independent filtering increases detection power for high-
throughput experiments. \emph{PNAS} \bold{107}(21):9546-9551.

P. Jaccard (1901). Etude comparative de la distribution orale dans une portion des Alpes et des Jura.
\emph{Bulletin de la Societe Vaudoise des Sciences Naturelles}, \bold{37}:547-549.

A. Rau, M. Gallopin, G. Celeux, F. Jaffrezic (2012). Independent data-based filtering 
for replicated high-throughput sequencing experiments. (Manuscript in progress).

}
\keyword{ package }
\seealso{
\code{\link{HTSFilter}}
}
\examples{

library(HTSFilter)
library(pasilla)
data("pasillaGenes")


## Matrix
dat <- counts(pasillaGenes)
conds <- pData(pasillaGenes)$condition
\dontrun{Jacc <- HTSFilter(dat, conds)}
\dontrun{dat <- Jacc$filteredData  ## Matrix}


## Count Data Set (DESeq)
cds <- pasillaGenes
\dontrun{Jacc <- HTSFilter(cds, conds = pData(cds)$condition)}
\dontrun{cds <- Jacc$filteredData  ## S4 class CountDataSet}


## DGEList (edgeR)
dge <- DGEList(counts=counts(pasillaGenes), group=pData(pasillaGenes)$condition)
\dontrun{Jacc <- HTSFilter(dge)}
\dontrun{dge <- Jacc$filteredData  ## S3 class DGEList}


## SeqExpressionSet (EDASeq)
ses <- newSeqExpressionSet(counts(pasillaGenes), phenoData = phenoData(pasillaGenes))
\dontrun{Jacc <- HTSFilter(ses, conds = pData(ses)$condition, s.len = 25)}
\dontrun{ses <- Jacc$filteredData  ## S4 class SeqExpressionSet}


}
