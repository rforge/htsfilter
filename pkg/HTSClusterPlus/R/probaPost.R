#' Calculate the conditional probability of belonging to each cluster in a
#' Poisson mixture model
#' 
#' This function computes the conditional probabilities \eqn{t_{ik}} that an
#' observation \emph{i} arises from the \eqn{k^{\mathrm{th}}}{kth} component
#' for the current value of the mixture parameters.
#' 
#' @inheritParams PoisMixClus_K
#' @param lambda (\emph{d} x \code{K}) matrix containing the estimate of
#' \eqn{\hat{\ensuremath\boldsymbol{\lambda}}}{\hat{\lambda}}
#' @param pi Vector of length \emph{g} containing the estimate of
#' \eqn{\hat{\ensuremath\boldsymbol{\pi}}}{\hat{\pi}}
#' @param s Library size normalization factors
#' 
#' @return \item{t }{(\emph{n} x \code{K}) matrix made up of the conditional
#' probability of each observation belonging to each of the \code{K} clusters}

#' @note If all values of \eqn{t_{ik}} are 0 (or nearly zero), the observation
#' is assigned with probability one to belong to the cluster with the closest
#' mean (in terms of the Euclidean distance from the observation). To avoid
#' calculation difficulties, extreme values of \eqn{t_{ik}} are smoothed, such
#' that those smaller than 1e-10 or larger than 1-1e-10 are set equal to 1e-10
#' and 1-1e-10, respectively.
#' @author Andrea Rau <\url{andrea.rau@@jouy.inra.fr}>
#' @seealso \code{\link{PoisMixClus}} for Poisson mixture model estimation and
#' model selection; \code{\link{PoisMixMean}} to calculate the conditional
#' per-cluster mean of each observation
#' @references Rau, A., Maugis-Rabusseau, C., Martin-Magniette, M.-L., Celeux,
#' G. (2015) Co-expression analysis of high-throughput transcriptome sequencing
#' data with Poisson mixture models. Bioinformatics, doi:
#' 10.1093/bioinformatics/btu845.
#' 
#' Rau, A., Celeux, G., Martin-Magniette, M.-L., Maugis-Rabusseau, C. (2011).
#' Clustering high-throughput sequencing data with Poisson mixture models.
#' Inria Research Report 7786. Available at
#' \url{http://hal.inria.fr/inria-00638082}.
#' @keywords methods
#' @example /inst/examples/PoisMixClus.R
#' @importFrom stats dist
#' @export probaPost

probaPost <-
  function(y, K, conds, pi, s, lambda) {
    
    if(length(K) != 1 | K < 0 | round(K) != K) 
      stop(paste(sQuote("K"), "(the number of clusters) must be a nonnegative integer"))
    if(is.vector(conds) == FALSE | length(conds) != ncol(y))
      stop(paste(sQuote("conds"), "must be a vector the same length as the number of columns in", sQuote("y")))
    if(is.vector(pi) == FALSE | length(pi) != K) 
      stop(paste(sQuote("pi"), "must be a vector of length", sQuote("K")))
    if(length(s) != length(conds))
      stop(paste(sQuote("s"), "and", sQuote("conds"), 
                 "must be vectors the same length as the number of columns in", sQuote("y")))
    if(is.matrix(lambda) == FALSE | ncol(lambda) != K | nrow(lambda) != length(unique(conds)))
      stop(paste(sQuote("lambda"), "must be a (d x K) matrix"))
    
    
    n <- dim(y)[1]; cols <- dim(y)[2];
    t <- matrix(0, nrow = n, ncol = K)
    mean <- PoisMixMean(y, K, conds, s, lambda)
    t <- matrix(unlist(lapply(1:K, function(x) .myprobafxn(k=x, y=y, pi=pi, mean=mean)), use.names=F), 
                nrow=n, ncol=K)
    ## Fix problematic values of t (= 0 for all clusters)
    for(j in which(rowSums(t) == 0)) {
      mean.list <- matrix(unlist(lapply(mean, function(x) x[j,]), use.names=F), ncol=length(conds), byrow=T)
      distance <- as.matrix(dist(rbind(y[j,], mean.list)))[,1]
      distance <- distance[-1]
      ## If distances are exactly the same, arbitrarily pick the first as closest
      cluster <- which(distance == min(distance))[1]
      t[j,cluster] <- 1
    }
    ## Normalize t: I think this is an error here
    ##t <- apply(t, 2, function(x) x/rowSums(t))
    t <- t / rowSums(t)
    ## Smoothing prior to M-Step (0's set to 1e-10, 1's set to 1-1e-10)
    epsilon <- 1e-10;maxcut <- 1-epsilon; mincut <- epsilon
    t <- apply(t, 2, pmax, mincut); t <- apply(t, 2, pmin, maxcut);
    ## ADDED
    t <- t / rowSums(t)
    
    return(t)
    
  }

