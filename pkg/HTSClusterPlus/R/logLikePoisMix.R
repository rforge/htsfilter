#' Log likelihood calculation for a Poisson mixture model
#' 
#' Function to calculate the log likelihood for a Poisson mixture model.
#' 
#' The \code{logLikePoisMix} function
#' (taken largely from the \code{mylogLikePoisMix} function from the
#' \code{poisson.glm.mix} R package) calculates the log likelihood for a given
#' set of parameters in a Poisson mixture model and is used in the
#' \code{\link{PoisMixClus_K}} function for the calculation of the BIC and ICL.
#' 
#' @param y (\emph{n} x \emph{q}) matrix of observed counts for \emph{n}
#' observations and \emph{q} variables
#' @param mean List of length \emph{K} containing the (\emph{n} x \emph{q})
#' matrices of conditional mean expression for all observations, as calculated
#' by the \code{\link{PoisMixMean}} function, where \emph{K} represents the
#' number of clusters
#' @param pi Vector of length \emph{K} containing estimate for
#' \eqn{\hat{\ensuremath\boldsymbol{\pi}}}{\hat{\pi}}
#' 
#' @return 
#' \item{ll }{Log-likelihood}
#' \item{ll_obs }{Per-observation log-likelihood}
#' 
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
#' 
#' @example /inst/examples/logLikePoisMix.R
#' 
#' @export logLikePoisMix
#' @importFrom stats dpois
#' 
logLikePoisMix <- 
  function (y, mean, pi) 
  {
    if (length(mean) != length(pi)) 
      stop(paste(sQuote("mean"), "must be a list of the same length as", 
                 sQuote("pi")))
    g <- length(pi)
    n <- dim(y)[1]
    cols <- dim(y)[2]
    nas <- 0
    y <- matrix(y, nrow = n, ncol = cols)
    logLike <- rep(0, n)
    index <- 1:g
    epsilon <- exp(-720)
    thresh <- -745
    nn <- 0
    logpi <- log(pi)
    ef <- matrix(logpi, nrow = n, ncol = g, byrow = T)
    for (k in 1:g) {
      ef[, k] <- ef[, k] + rowSums(dpois(y, mean[[k]], log = T))
    }
    efmax <- apply(ef, 1, max)
    ef <- ef - efmax
    logLike <- efmax + log(rowSums(exp(ef)))
    return(list(ll_obs = logLike, ll = sum(logLike, na.rm = TRUE)))
  }

