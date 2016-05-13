#' Parameter initialization for a Poisson mixture model.
#' 
#' This function implements the Small EM initialization strategy
#' (\code{emInit}) described in Rau et al. (2011).
#' 
#' In practice, the user will not directly call the initialization functions
#' described here; they are indirectly called for a single number of clusters
#' through the \code{PoisMixClus_K} function (via \code{init.type}) or via the
#' \code{PoisMixClus} function for a sequence of cluster numbers (via
#' \code{Kmin.init} and \code{split.init}).
#' 
#' To initialize parameter values for the EM and CEM algorithms, for the
#' Small-EM strategy (Biernacki et al., 2003) we use the \code{emInit} function
#' as follows. For a given number of independent runs (given by
#' \code{init.runs}), the following procedure is used to obtain parameter
#' values: first, a K-means algorithm (MacQueen, 1967) is run to partition the
#' data into \code{g} clusters
#' (\eqn{\hat{\ensuremath\boldsymbol{z}}^{(0)}}{\hat{z}^(0)}). Second, initial
#' parameter values \eqn{\ensuremath\boldsymbol{\pi}^{(0)}}{\pi^(0)} and
#' \eqn{\ensuremath\boldsymbol{\lambda}^{(0)}}{\lambda^(0)} are calculated (see
#' Rau et al. (2011) for details). Third, a given number of iterations of an EM
#' algorithm are run (defined by \code{init.iter}), using
#' \eqn{\ensuremath\boldsymbol{\pi}^{(0)}}{\pi^(0)} and
#' \eqn{\ensuremath\boldsymbol{\lambda}^{(0)}}{\lambda^(0)} as initial values.
#' Finally, among the \code{init.runs} sets of parameter values, we use
#' \eqn{\hat{\ensuremath\boldsymbol{\lambda}}}{\hat{\lambda}} and
#' \eqn{\hat{\ensuremath\boldsymbol{\pi}}}{\hat{\pi}} corresponding to the
#' highest log likelihood or completed log likelihood to initialize the
#' subsequent full EM or CEM algorithms, respectively.
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
#' @export

emInit <- function(y, K, conds, norm, alg.type = "EM", 
                   init.runs, init.iter, fixed.lambda, equal.proportions, 
                   verbose) {
  
  if(alg.type != "EM" & alg.type != "CEM")
    stop(paste(sQuote("alg.type"), "must be one of", dQuote("EM"), "or", dQuote("CEM")))
  if(length(alg.type) > 1)
    stop(paste(sQuote("alg.type"), "must be one of", dQuote("EM"), "or", dQuote("CEM")))
  if(init.runs < 1 | length(init.runs) > 1 | round(init.runs) != init.runs) 
    stop(paste(sQuote("init.runs"), "must be a positive integer"))
  if(is.logical(verbose) == FALSE)
    stop(paste(sQuote("verbose"), "must be", dQuote("TRUE"), "or", dQuote("FALSE")))
  init.type1 <- "kmeans"
  lambda.init.all <- vector("list", init.runs)
  pi.init.all <- vector("list", init.runs)
  criterion.all <- rep(NA, init.runs)
  for(start in 1:init.runs) {
    em.init <- PoisMixClus_K(y = y, K=K, norm=norm, conds = conds, 
                           init.type = init.type1, alg.type = alg.type, iter = init.iter,
                           fixed.lambda = fixed.lambda, equal.proportions = equal.proportions,
                           wrapper=FALSE)
    lambda.init.all[[start]] <- em.init$lambda
    pi.init.all[[start]] <- em.init$pi
    criterion.all[start] <- em.init$log.like
    if(verbose == TRUE) print(paste("Initialization:", start))
  }
  ## If all criterion values are equal to NaN, then arbitrarily choose the first one
  if(sum(is.na(criterion.all)) == length(criterion.all)) {
    final.choice <- 1;
  }
  ## If two of the criterion values are exactly the same, pick only the first
  if(sum(is.na(criterion.all)) != length(criterion.all)) {
    final.choice <- which(criterion.all == min(criterion.all, na.rm = TRUE))[1]
  }
  lambda.init <- lambda.init.all[[final.choice]]
  pi.init <- pi.init.all[[final.choice]]
  return(list(pi.init = pi.init, lambda.init = lambda.init))
}

