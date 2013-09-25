\name{initialDataset}
\docType{data}
\alias{initialDataset}
\title{RNA-seq data from Arabidopsis thaliana to be used to generate synthetic simulations}
\description{
This dataset represents RNA-seq data from Arabidopsis thaliana in two biological replicates in flower buds (BF1 and BF2) and leaves (F1 and F2), 
as well as an additional two biological replicates of leaves (Col0-1, Col0-2). In this package, these data are used as a basis to generate synthetic
data.
}
\usage{data(initialDataset)}
\format{
A \code{data.frame} containing 28,094 rows (genes) and 7 columns: 
ID (the Ensembl gene ID), BF1 and BF2 (the read counts for the two biological replicates in flower buds),
F1 and F2 (the read counts for tbe two biological replicates in leaves), and Col0-1 and Col0-2 (the
read counts for the two additional biological replicates in leaves).

Note that the 332 genes that have been validated as differentially expressed via qPCR have been marked in the data
by the addition of ".DE" to the end of their gene ID, and the 81 validated as non-differentially expressed via
qPCR have been marked by an additional ".NDE" at the end of their gene ID.
}
\keyword{datasets}

