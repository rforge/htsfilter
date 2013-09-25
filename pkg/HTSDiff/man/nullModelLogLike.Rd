\name{nullModelLogLike}
\alias{nullModelLogLike}
\title{
Calculation of loglikelihood for null model in differential analysis.
}
\description{
This function calculates the loglikelihood, BIC, and ICL for the null model, i.e.,
a single cluster of non-differentially expressed genes.
}
\usage{
nullModelLogLike(counts, conds, norm="DESeq")
}
\arguments{
  \item{counts}{
(\emph{n} x \emph{q}) matrix of observed counts for \emph{n} genes and \emph{q} samples
}
  \item{conds}{
Vector of length \emph{q} defining the condition (treatment group) for each variable (column) in \code{counts}
}
  \item{norm}{
The estimator to be used for the library size parameter: \dQuote{\code{TC}} for total count, \dQuote{\code{UQ}} 
for upper quantile, \dQuote{\code{Med}} for median, \dQuote{\code{DESeq}} for the normalization method in the 
DESeq package, and \dQuote{\code{TMM}} for the TMM normalization method.
}
}
\details{

This function implements the calculation of the loglikelihood, BIC, and ICL for the model
containing a single cluster of non-differentially expressed genes. Its main utility is to enable
model selection between this null model and a model containing at least one potential cluster of differentially
expressed genes.

}
\value{
\item{logLike }{Loglikelihood value}
\item{BIC }{Value of BIC}
\item{ICL }{Value of ICL criterion}
}

\references{
S. Balzergue, G. Rigaill, V. Brunaud, E. Blondet, A. Rau, O. Rogier, J. Caius, C. Maugis-Rabusseau,
L. Soubigou-Taconnat, S. Aubourg, C. Lurin, E. Delannoy, and M.-L. Martin-Magniette. (2013) 
HTSDiff: a sensitive analysis of differential gene expression by RNA-seq (submitted).
}
\author{
Andrea Rau <\url{andrea.rau@jouy.inra.fr}>
}

\examples{

set.seed(12345)

#### Generate synthetic data: 2000 genes under H0
## syn <- syntheticData(H0number = 2000)

#### Calculate criteria for null model
## nulltest <- nullModelLogLike(syn, conds=c(1,1,2,2),
##	norm="DESeq")

}
\keyword{ models }

