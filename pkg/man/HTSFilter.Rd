\name{HTSFilter}
\docType{methods}
\alias{HTSFilter}
\alias{HTSFilter-methods}
\alias{HTSFilter,matrix-method}
\alias{HTSFilter,CountDataSet-method}
\alias{HTSFilter,DGEList-method}
\alias{HTSFilter,SeqExpressionSet-method}
\title{
Calculate data-based filtering threshold for replicated high-throughput sequencing data.
}
\description{
Calculate a data-based filtering threshold for replicated high-throughput
sequencing data through the pairwise Jaccard similarity index between pairs
of replicates within each experimental condition.
}

\usage{

\S4method{HTSFilter}{matrix}(x, conds, s.min=1, s.max=200, s.len=100, 
     loess.span=0.3, normalization=c("TMM", "DESeq", "none"), 
     plot=TRUE, plot.name=NA)

\S4method{HTSFilter}{CountDataSet}(x, conds=NA, s.min=1, s.max=200, s.len=100, 
     loess.span=0.3, normalization=c("TMM", "DESeq", "none"), 
     plot=TRUE, plot.name=NA)

\S4method{HTSFilter}{DGEList}(x, s.min=1, s.max=200, s.len=100, 
     loess.span=0.3, normalization=c("TMM", "DESeq", "none"), 
     plot=TRUE, plot.name=NA)

\S4method{HTSFilter}{SeqExpressionSet}(x, conds=NA, s.min=1, s.max=200, s.len=100, 
     loess.span=0.3, normalization=c("TMM", "DESeq", "none"), 
     plot=TRUE, plot.name=NA)


}

\arguments{
  \item{x}{
A numeric matrix representing the counts of dimension (\emph{g} x \emph{n}), for \emph{g} genes in \emph{n} samples;
a \code{\linkS4class{CountDataSet}} object; a \code{DGEList} object; or a \code{\linkS4class{SeqExpressionSet}} object.
}
  \item{conds}{
Vector of length \emph{n} identifying the experimental condition of each of the \emph{n} samples; required when sQuote(x)
is a numeric matrix, and optionally for \code{\linkS4class{CountDataSet}} and \code{\linkS4class{SeqExpressionSet}}.
}
  \item{s.min}{
Minimum value of filtering threshold to be considered, with default value equal to 1
}
  \item{s.max}{
Maximum value of filtering threshold to be considered, with default value equal to 200
}
  \item{s.len}{
Length of sequence of filtering thresholds to be considered (from \code{s.min} to \code{s.max}) 
for the calculation of the global similarity index
}
  \item{loess.span}{
Span of the loess curve to be fitted to the filtering thresholds and corresponding global similarity
indices, with default value equal to 0.3
}
  \item{normalization}{
Normalization method to be used to correct for differences in library sizes, with choices
\dQuote{TMM} (Trimmed Mean of M-values), \dQuote{DESeq} (normalization method proposed in the
DESeq package), and \dQuote{none} (to be used only if user is certain no normalization is required,
or if data have already been pre-normalized by an alternative method)
}
  \item{plot}{
If \dQuote{TRUE}, produce a plot of the calculated global similarity indices against the
filtering threshold with superimposed loess curve
}
  \item{plot.name}{
If \code{plot} = \dQuote{TRUE}, the name of the PDF file to be saved to the current working directory.
If \code{plot.name} = NA, the plot is drawn in the current window.
}
}


\details{

The Jaccard similarity index, which measures the overlap of two sets, is calculated as follows. 
Given two binary vectors, each of length \emph{n}, we define the following values:
\itemize{
  \item \emph{a} = the number of attributes with a value of 1 in both vectors
  \item \emph{b} = the number of attributes with a value of 1 in the first vector and 0 in the second
  \item \emph{c} = the number of attributes with a value of 0 in the first vector and 1 in the second
  \item \emph{d} = the number of attributes with a value of 0 in both vectors
}
We note that all attributes fall into one of these four quantities, so \eqn{a+b+c+d=n}. Given these
quantities, we may calculate the Jaccard similarity index between the two vectors as follows:
\deqn{J = \frac{a}{a+b+c}.}{J = a/(a+b+c).}

}

\value{
\item{filteredData }{An object of the same class as \code{x} containing the data that passed the filter.
When \code{x} is a \code{\linkS4class{CountDataSet}} object, the filtered data are contained in the \code{countData}
slot; when \code{x} is a \code{DGEList} object, the filtered data are contained in the \code{counts} slot; 
when \code{x} is a \code{\linkS4class{SeqExpressionSet}} object, the filtered data are contained in the \code{exprs} slot.}
%
\item{on }{A binary vector of length \emph{g}, where 1 indicates a gene with normalized expression
greater than the optimal filtering threshold \code{s.optimal} in at least one sample (irrespective of 
condition labels), and 0 indicates a gene with normalized expression less than or equal to the optimal 
filtering threshold in all samples}
%
\item{s }{The optimal filtering threshold as identified by the global similarity index}
%
\item{indexValues }{A matrix of dimension (\code{s.len} x 2) giving the tested filtering thersholds and the
corresponding global similarity indices}
%
\item{normFactor }{A vector of length \emph{n} giving the estimated library sizes estimated by the
normalization method specified in \code{normalization}}
%
\item{removedData }{A matrix containing the filtered data}
}

\references{

R. Bourgon, R. Gentleman, and W. Huber. (2010) Independent filtering increases detection power for high-
throughput experiments. \emph{PNAS} \bold{107}(21):9546-9551.

P. Jaccard (1901). Etude comparative de la distribution orale dans une portion des Alpes et des Jura.
\emph{Bulletin de la Societe Vaudoise des Sciences Naturelles}, \bold{37}:547-549.

A. Rau, M. Gallopin, G. Celeux, F. Jaffrezic (2012). Independent data-based filtering 
for replicated high-throughput sequencing experiments. (Manuscript in progress).

}
\author{
Andrea Rau, Melina Gallopin, Gilles Celeux, and Florence Jaffrezic
}
\note{
Filter should only be calculated on REPLICATED high-throughput sequencing data.
}

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
\keyword{ methods }

