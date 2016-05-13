#' Splitting EM parameter initialization for a Poisson mixture model.
#' 
#' This function implements the splitting small-EM initialization strategy
#' (\code{splitEMInit}) based on that described in Papastamoulis et al. (2014)
#' for a Poisson mixture model.
#' 
#' In practice, the user will not directly call the initialization functions
#' described here; they are indirectly called for a single number of clusters
#' through the \code{PoisMixClus_K} function (via \code{init.type}) or via the
#' \code{PoisMixClus} function for a sequence of cluster numbers (via
#' \code{Kmin.init} and \code{split.init}).
#' 
#' 
#' For the splitting small EM initialization strategy, we implement an approach
#' similar to that described in Papastamoulis et al. (2014), where the cluster
#' from the previous run (with \emph{K}-1 clusters) with the largest entropy is
#' chosen to be split into two new clusters, followed by a small EM run.
#' 
#' @inheritParams PoisMixClus_K
#' 
#' @return \item{pi.init }{Vector of length \code{K} containing the estimate
#' for \eqn{\hat{\ensuremath\boldsymbol{\pi}}}{\hat{\pi}} corresponding to the
#' highest log likelihood (or completed log likelihood) from the chosen
#' inialization strategy. }
#' 
#' \item{lambda.init }{(\emph{d} x \code{K}) matrix containing the estimate of
#' \eqn{\hat{\ensuremath\boldsymbol{\lambda}}}{\hat{\lambda}} corresponding to
#' the highest log likelihood (or completed log likelihood) from the chosen
#' initialization strategy, where \emph{d} is the number of conditions and
#' \code{K} is the number of clusters. }
#' 
#' @author Andrea Rau <\url{andrea.rau@@jouy.inra.fr}>
#' @seealso \code{\link{PoisMixClus_K}} for Poisson mixture model estimation for
#' a given number of clusters, \code{\link{PoisMixClus}} for Poisson
#' mixture model estimation and model selection for a sequence of cluster
#' numbers.
#' @references Anders, S. and Huber, W. (2010) Differential expression analysis
#' for sequence count data. \emph{Genome Biology}, \bold{11}(R106), 1-28.
#' 
#' Biernacki, C., Celeux, G., Govaert, G. (2003) Choosing starting values for
#' the EM algorithm for getting the highest likelhiood in multivariate Gaussian
#' mixture models. \emph{Computational Statistics and Data Analysis},
#' \bold{41}(1), 561-575.
#' 
#' MacQueen, J. B. (1967) Some methods for classification and analysis of
#' multivariate observations. In \emph{Proceedings of the 5th Berkeley
#' Symposium on Mathematical Statistics and Probability}, number 1, pages
#' 281-297. Berkeley, University of California Press.
#' 
#' Papastamoulis, P., Martin-Magniette, M.-L., and Maugis-Rabusseau, C. (2014).
#' On the estimation of mixtures of Poisson regression models with large number
#' of components. \emph{Computational Statistics and Data Analysis}: 3rd
#' special Issue on Advances in Mixture Models, DOI:
#' 10.1016/j.csda.2014.07.005.
#' 
#' Rau, A., Celeux, G., Martin-Magniette, M.-L., Maugis-Rabusseau, C. (2011).
#' Clustering high-throughput sequencing data with Poisson mixture models.
#' Inria Research Report 7786. Available at
#' \url{http://hal.inria.fr/inria-00638082}.
#' 
#' Rau, A., Maugis-Rabusseau, C., Martin-Magniette, M.-L., Celeux, G. (2015)
#' Co-expression analysis of high-throughput transcriptome sequencing data with
#' Poisson mixture models. Bioinformatics, doi: 10.1093/bioinformatics/btu845.
#' 
#' Robinson, M. D. and Oshlack, A. (2010) A scaling normalization method for
#' differential expression analysis of RNA-seq data. \emph{Genome Biology},
#' \bold{11}(R25).
#' @keywords models
#' @importFrom stats runif
#' @export
#' 
splitEMInit <- function(y, K, conds, norm, alg.type, fixed.lambda, equal.proportions, 
                        prev.labels, prev.probaPost, init.runs, init.iter, verbose) {
  
  ## K is the new number of clusters IN ADDITION TO FIXED LAMBDA
  ## NB: This function inspiried by init2.k() function from poisson.glm.mix package
  ## (written by Panos Papastamoulis)
  
  n <- dim(y)[1]
  unique.labels <- unique(prev.labels)
  ## Check whether any of the clusters has one or zero observations & remove from consideration
  tab <- table(prev.labels) 
  ## Fix to make sure that we only consider clusters with at least 2 observations
  unique.labels <- names(tab)
  if(length(which(tab < 2)) > 0) {
    unique.labels <- names(tab)[-which(tab < 2)]
  }
  
  if(class(fixed.lambda) == "list") {
    K <- K + length(fixed.lambda);	
  }	
  prev.K <- K - 1
  
  ## Calculate per-class entropy of previous model
  ## and choose cluster with largest entropy to split
  perEntropy <- rep(NA, length(unique.labels))
  names(perEntropy) <- unique.labels
  for(k in 1:length(unique.labels)) {
    perEntropy[k] <- -sum(log(prev.probaPost[which(prev.labels == as.numeric(unique.labels[k])),
                                             as.numeric(unique.labels[k])]))
  }
  cluster.choose <- as.numeric(unique.labels[which(perEntropy == max(perEntropy))])
  index1 <- which(prev.labels == cluster.choose)
  
  ## Random selection of observations within splitted component
  ## Repeat this init.runs times
  LL.all <- rep(NA, init.runs)
  init.all <- vector("list", init.runs)
  for(iter in 1:init.runs) {
    u.numbers <- runif(length(index1))
    t <- matrix(0, nrow = n, ncol = K)
    t[,1:prev.K] <- prev.probaPost
    t[index1,K] <- t[index1,cluster.choose] * u.numbers
    t[index1,cluster.choose] <- t[index1,cluster.choose] * (1-u.numbers)
    
    ## Smoothing t values
    epsilon <- 1e-10
    maxcut <- 1 - epsilon; mincut <- epsilon
    t <- apply(t, 2, pmax, mincut)
    t <- apply(t, 2, pmin, maxcut)
    t <- t/rowSums(t)
    
    ## Initialize small EM-algorithm with new z values
    initialize <- probaPostInit(y = y, K = K, norm=norm,
                                conds = conds, alg.type = alg.type,
                                fixed.lambda = fixed.lambda,
                                equal.proportions = equal.proportions, probaPost.init = t,
                                init.iter = init.iter, verbose = verbose)
    LL.all[iter] <- initialize$log.like
    init.all[[iter]] <- initialize
  }
  
  init.select <- which(LL.all == max(LL.all, na.rm = TRUE))[1]
  final.init <- init.all[[init.select]]
  lambda.init <- final.init$lambda
  pi.init <- final.init$pi
  
  return(list(lambda.init = lambda.init, pi.init = pi.init))
}



