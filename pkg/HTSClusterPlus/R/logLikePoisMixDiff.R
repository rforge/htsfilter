#' Log likelihood difference calculation for a Poisson mixture model
#' 
#' Function to calculate the
#' difference in log likelihoods for two different sets of parameters of a
#' Poisson mixture model.
#' 
#' The \code{logLikePoisMixDiff} function is used to calculate the difference
#' in log likelihood for two different sets of parameters in a Poisson mixture
#' model; it is used to determine convergence in the EM algorithm run by the
#' \code{\link{PoisMixClus_K}} function.
#' 
#' @param y (\emph{n} x \emph{q}) matrix of observed counts for \emph{n}
#' observations and \emph{q} variables
#' @param mean.new List of length \emph{K} containing the (\emph{n} x \emph{q})
#' matrices of conditional mean expression for all observations for one set of
#' parameters, as calculated by the \code{\link{PoisMixMean}} function, where
#' \emph{K} represents the number of clusters
#' @param mean.old List of length \emph{K} containing the (\emph{n} x \emph{q})
#' matrices of conditional mean expression for all observations for another set
#' of parameters, as calculated by the \code{\link{PoisMixMean}} function,
#' where \emph{K} represents the number of clusters
#' @param pi.new Vector of length \emph{K} containing one estimate for
#' \eqn{\hat{\ensuremath\boldsymbol{\pi}}}{\hat{\pi}}
#' @param pi.old Vector of length \emph{K} containing another estimate for
#' \eqn{\hat{\ensuremath\boldsymbol{\pi}}}{\hat{\pi}}
#' 
#' @return \item{ll }{The difference in log likelihoods for two different sets of parameters.}
#' @note In the \code{logLikePoisMixDiff} function, we make use of the
#' alternative mass function for a Poisson density proposed by Loader (2000) to
#' avoid computational difficulties. The \code{logLikePoisMixDiff} function
#' returns a default value of 100 if one or both of the log likelihoods
#' associated with the two parameter sets takes on a value of
#' \eqn{-\infty}{-\infty}.
#' @author Andrea Rau <\url{andrea.rau@@jouy.inra.fr}>
#' @seealso \code{\link{PoisMixClus_K}} for Poisson mixture model estimation and
#' model selection; \code{\link{PoisMixMean}} to calculate the per-cluster
#' conditional mean of each observation
#' @references
#' 
#' Loader, C. (2000) Fast and accurate computation of binomial probabilities.
#' Available at
#' \url{http://projects.scipy.org/scipy/raw-attachment/ticket/620/loader2000Fast.pdf}.
#' 
#' Rau, A., Maugis-Rabusseau, C., Martin-Magniette, M.-L., Celeux, G. (2015)
#' Co-expression analysis of high-throughput transcriptome sequencing data with
#' Poisson mixture models. Bioinformatics, doi: 10.1093/bioinformatics/btu845.
#' 
#' Rau, A., Celeux, G., Martin-Magniette, M.-L., Maugis-Rabusseau, C. (2011)
#' Clustering high-throughput sequencing data with Poisson mixture models.
#' Inria Research Report 7786. Available at
#' \url{http://hal.inria.fr/inria-00638082}.
#' @keywords methods
#' @example /inst/examples/logLikePoisMix.R
#' 
#' @export logLikePoisMixDiff

logLikePoisMixDiff <- function(y, mean.new, pi.new, mean.old, pi.old) {
  
  if(length(mean.old) != length(pi.old)) 
    stop(paste(sQuote("mean.old"), "must be a list of the same length as", sQuote("pi.old")))
  if(length(mean.new) != length(pi.new)) 
    stop(paste(sQuote("mean.new"), "must be a list of the same length as", sQuote("pi.new")))
  if(length(mean.old) != length(mean.new)) 
    stop(paste(sQuote("mean.new"), "must be a list of the same length as", sQuote("mean.old")))
  
  n <- dim(y)[1]; cols <- dim(y)[2]
  g <- length(pi.old)
  y <- matrix(y, nrow = n, ncol = cols)
  num.tmp <- matrix(unlist(lapply(mean.new, function(x) .myfxn(var1=y, var2=x)), use.names=F, 
                           recursive=F), nrow=n, ncol=g)
  den.tmp <- matrix(unlist(lapply(mean.old, function(x) .myfxn(var1=y, var2=x)), use.names=F, 
                           recursive=F), nrow=n, ncol=g)
  num <- rowSums(matrix(pi.new * matrix(exp(-num.tmp),byrow=T,nrow=g),byrow=T,ncol=g,nrow=n))
  den <- rowSums(matrix(pi.old * matrix(exp(-den.tmp),byrow=T,nrow=g),byrow=T,ncol=g,nrow=n))
  num <- ifelse(num == 0, NA, num)
  den <- ifelse(den == 0, NA, den)
  ll.tmp <- ifelse(is.nan(log(num) - log(den)) == TRUE, NA, log(num) - log(den))
  ll <- sum(ll.tmp, na.rm = TRUE)
  ## In case one or both of the log-likelihoods is infty or -infy, set difference to large number
  if(ll == -Inf | ll == Inf | is.nan(ll) == TRUE) {ll = 100}
  return(ll)
}
