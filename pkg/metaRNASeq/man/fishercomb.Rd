\name{fishercomb}
\alias{fishercomb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
p-value combination using Fisher's method
}
\description{
Combines one sided p-values using Fisher's method.
}
\usage{
fishercomb(indpval, BHth = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{indpval}{
List of vectors of one sided p-values to be combined.
}
  \item{BHth}{
Benjamini Hochberg threshold. By default, the False Discovery Rate is controlled at 5\%.
}
}
\details{
%%  Mettre une formule?
}
\value{
\item{DEindices }{Indices of differentially expressed genes at the chosen Benjamini Hochberg threshold.}
\item{TestStatistic }{Vector with test statistics for differential expression in the meta-analysis.}
\item{rawpval }{Vector with raw p-values for differential expression in the meta-analysis.}
\item{adjpval }{Vector with adjusted p-values for differential expression in the meta-analysis.}
}
\references{
M. Brown  (1975). A method for combining non-independent, one-sided tests of significance. \emph{Biometrics} \bold{31}(4): 987–992.

A. Rau, G. Marot and F. Jaffrezic (2012). Differential meta-analysis of RNA-seq data. (Manuscript in progress)
}
\seealso{
\code{\link{metaRNASeq}}
}
\examples{
data(rawpval)
fishcomb <- fishercomb(rawpval, BHth = 0.05)
DE <- ifelse(fishcomb$adjpval<=0.05,1,0)
hist(fishcomb$rawpval,nclass=100)
##A more detailed example is given in the vignette of the package
vignette("metaRNASeq")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
\keyword{ models }% __ONLY ONE__ keyword per line
