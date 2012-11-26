\name{invnorm}
\alias{invnorm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
p-value combination using the inverse normal method
}
\description{
Combines one sided p-values using the inverse normal method.
}
\usage{
invnorm(pvalonesided, nrep, BHth = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pvalonesided}{
List of vectors of one sided p-values to be combined.
}
  \item{nrep}{
Vector of numbers of replicates used in each study to calculate the previous one-sided p-values. 
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
\references{ Hedges, L. and Olkin, I. (1985). Statistical Methods for Meta-Analysis. Academic
Press. 

A. Rau, G. Marot and F. Jaffrezic (2012). Differential meta-analysis of RNA-seq data. (Manuscript in progress)}
\note{
This function looks like \textit{directpvalcombi{metaMA}}. There is, however, one important difference to calculate p-values. 
%% A reprendre! For microarrays, it was advised to separate under- and over-expressed genes before the meta-analysis, which is not anymore possible with the use of Poisson like distributions in individual RNAseq studies. That is why the final p-values are one-sided p-values for RNA-seq data, while they were two sided p-values for microarrays.  
}
\seealso{
\code{\link{metaRNASeq}}
}
\examples{
data(rawpval)
invnormcomb <- invnorm(rawpval,nrep=c(8,8), BHth = 0.05)       #8 replicates had been simulated in each study
DE <- ifelse(invnormcomb$adjpval<=0.05,1,0)
hist(invnormcomb$rawpval,nclass=100)
##A more detailed example is given in the vignette of the package
vignette("metaRNASeq")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
\keyword{ models }% __ONLY ONE__ keyword per line
