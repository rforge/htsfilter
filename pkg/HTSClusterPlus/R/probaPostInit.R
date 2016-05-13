#' Parameter initialization for a Poisson mixture model.
#' 
#' This function implements a small-EM initialization strategy using
#' probabilities (\code{probaPostInit}) obtained from a previous run with one
#' fewer cluster following the splitting strategy.
#' 
#' In practice, the user will not directly call the initialization functions
#' described here; they are indirectly called for a single number of clusters
#' through the \code{PoisMixClus_K} function (via \code{init.type}) or via the
#' \code{PoisMixClus} function for a sequence of cluster numbers (via
#' \code{Kmin.init} and \code{split.init}).
#' 
#' @inheritParams PoisMixClus_K
#' @param probaPost.init (\emph{n} x \code{K-1}) matrix of posterior probabilities
#' from the model with \emph{K}-1 clusters
#' 
#' @return 
#' 
#' \item{lambda }{(\emph{d} x \code{K}) matrix containing the estimate of
#' \eqn{\hat{\ensuremath\boldsymbol{\lambda}}}{\hat{\lambda}} arising from the
#' splitting initialization and small EM run for a single split, where \emph{d}
#' is the number of conditions and \code{K} is the number of clusters. }
#' 
#' \item{pi }{Vector of length \code{K} containing the estimate for
#' \eqn{\hat{\ensuremath\boldsymbol{\pi}}}{\hat{\pi}} arising from the
#' splitting initialization and small EM run for a single split, where \code{K}
#' is the number of clusters.  }
#' 
#' \item{log.like }{Log likelihood arising from the splitting initialization
#' and small EM run for a single split. }
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
probaPostInit <- function(y, K, conds, norm,
                          alg.type = "EM", fixed.lambda, equal.proportions,
                          probaPost.init, init.iter, verbose) 
{
  
  ## fixed.lambda should be a list of length (number of fixed clusters)
  ## g gives the number of clusters IN ADDITION to the fixed clusters
  ## 	specified by fixed.lambda
  ## equal.proportions should be TRUE or FALSE
  
  if(is.vector(conds) == FALSE | length(conds) != ncol(y))
    stop(paste(sQuote("conds"), "must be a vector the same length as the number of columns in", sQuote("y")))
  if(alg.type != "EM" & alg.type != "CEM")
    stop(paste(sQuote("alg.type"), "must be one of", dQuote("EM"), "or", dQuote("CEM")))
  if(length(alg.type) > 1)
    stop(paste(sQuote("alg.type"), "must be one of", dQuote("EM"), "or", dQuote("CEM")))
  if(is.logical(verbose) == FALSE)
    stop(paste(sQuote("verbose"), "must be", dQuote("TRUE"), "or", dQuote("FALSE")))
  if(class(fixed.lambda) != "list" & is.na(fixed.lambda[1]) == FALSE)
    stop(paste(sQuote("fixed.lambda"), "must be", dQuote("NA") , "or a list."))
  if(is.matrix(probaPost.init) == FALSE)
    stop(paste(sQuote("z.init"), "must be a matrix of posterior probabilities."))
  
  conds.names <- unique(conds)
  
  d <- length(unique(conds))
  r <- as.vector(table(conds))
  diff <- 100 ## Convergence criterion
  if(length(rownames(y)) == 0) rn <- 1:nrow(y);
  if(length(rownames(y)) > 0) rn <- rownames(y);
  y <- as.matrix(y, nrow = nrow(y), ncol = ncol(y))
  rownames(y) <- rn;
  n <- dim(y)[1];cols <- dim(y)[2]
  
  w <- rowSums(y)
  ## Only calculate s values if they are not provided
  if(length(norm) == 1) {
    if(norm == "none") 	s <- rep(1, cols);
    if(norm == "TC") s <- colSums(y) / sum(as.numeric(y));
    if(norm == "UQ") s <- apply(y, 2, quantile, 0.75) / sum(apply(y, 2, quantile, 0.75));
    if(norm == "Med") s <- apply(y, 2, median) / sum(apply(y, 2, median));
    if(norm == "DESeq") {
      ## Code from DESeq, v1.8.3
      loggeomeans <- rowMeans(log(y))
      s <- apply(y, 2, function(x) exp(median((log(x)-loggeomeans)[is.finite(loggeomeans)])))
      s <- s / sum(s)
    }
    if(norm == "TMM") {
      f <- calcNormFactors(as.matrix(y), method = "TMM")
      s <- colSums(y)*f / sum(colSums(y)*f)
    } 
  }
  if(length(norm) == length(conds)) {
    s <- norm / sum(norm)
  }
  s.dot <- rep(NA, d) 
  for(j in 1:d) {
    s.dot[j] <- sum(s[which(conds == (unique(conds))[j])])
  }
  
  if(class(fixed.lambda) == "list") {
    K <- K + length(fixed.lambda);	
  }
  
  ## Inital values using probaPost.init
  pi <- pi.old <- rep(NA, K)
  lambda <- lambda.old <- matrix(NA, nrow = d, ncol = K)
  t <- probaPost.init
  
  for(index in 0:init.iter) {
    
    if(index > 0) {
      ############
      ## E-step ##
      ############
      t <- probaPost(y, K, conds, pi, s, lambda)
    }
    
    ############
    ## C-step ##
    ############
    if(alg.type == "CEM") {
      ## If two values of t_{ik} are map, 
      ## arbitrarily choose the first
      partition <- unlist(apply(t, 1, 
                                function(x) which(x == max(x, na.rm = TRUE))[1]))
      partition.mat <- matrix(0, nrow = n, ncol = K)
      for(i in 1:n) partition.mat[i,partition[i]] <- 1;
    }
    
    ############
    ## M-step ##
    ############
    if(alg.type == "CEM") {
      if(equal.proportions == FALSE) {
        for(k in 1:K) {
          pi[k] <- length(which(partition == k))/n
        }
      }
      if(equal.proportions == TRUE) {
        pi <- rep(1/K, K)
      }
      denom <- colSums(partition.mat * w)
      if(class(fixed.lambda) != "list") {
        for(j in 1:d) {
          denom.bis <- denom * s.dot[j]
          num <- colSums(partition.mat * 
                           rowSums(as.matrix(y[,which(conds == (unique(conds))[j])])))
          lambda[j,] <- num / denom.bis
        }
      }
      if(class(fixed.lambda) == "list") {
        for(ll in 1:length(fixed.lambda)) {
          lambda[,ll] <- fixed.lambda[[ll]]
        }
        for(j in 1:d) {
          denom.bis <- denom * s.dot[j]
          denom.bis <- denom.bis[-c(1:length(fixed.lambda))]
          num <- colSums(partition.mat * 
                           rowSums(as.matrix(y[,which(conds == (unique(conds))[j])])))
          num <- num[-c(1:length(fixed.lambda))]
          lambda[j,-c(1:length(fixed.lambda))] <- num / denom.bis
        }
      }
    }
    
    if(alg.type == "EM") {
      if(equal.proportions == FALSE) {
        pi <- colSums(t)/n
      }
      if(equal.proportions == TRUE) {
        pi <- rep(1/K, K)
      }
      denom <- colSums(t * w)
      if(class(fixed.lambda) != "list") {
        for(j in 1:d) {
          denom.bis <- denom * s.dot[j]
          num <- colSums(t * 
                           matrix(rep(rowSums(as.matrix(y[,which(conds == (unique(conds))[j])])),K), 
                                  ncol = K))
          lambda[j,] <- num / denom.bis
        }
      }
      if(class(fixed.lambda) == "list") {
        for(ll in 1:length(fixed.lambda)) {
          lambda[,ll] <- fixed.lambda[[ll]]
        }
        for(j in 1:d) {
          denom.bis <- denom * s.dot[j]
          denom.bis <- denom.bis[-c(1:length(fixed.lambda))]
          num <- colSums(t * 
                           matrix(rep(rowSums(as.matrix(y[,which(conds == (unique(conds))[j])])),K), 
                                  ncol = K))
          num <- num[-c(1:length(fixed.lambda))]
          lambda[j,-c(1:length(fixed.lambda))] <- num / denom.bis
        }
      }	
    }
    lambda.old <- lambda; pi.old <- pi;
  }
  
  #####################################
  ## Final estimates of lambda and p ##
  #####################################
  names(pi) <- paste("Cluster", 1:K)
  colnames(lambda) <- paste("Cluster", 1:K)
  rownames(lambda) <- conds.names
  lambda.final <- lambda
  pi.final <- pi
  
  ## Check to make sure one of the components is not degenerate
  if(min(pi) == 0 | is.nan(sum(lambda)) == TRUE) {
    probaPost <- NA; labels <- NA; BIC <- NA;	ICL <- NA
  }
  
  if(min(pi) > 0 | is.nan(sum(lambda)) == FALSE) {
    mean.calc <- PoisMixMean(y, K = K, conds, s, lambda)
    LL.tmp <- logLikePoisMix(y, mean.calc, pi)
    LL <- LL.tmp$ll
  }
  results <- list(lambda = lambda.final, pi = pi.final, log.like = LL)
  return(results)
}

