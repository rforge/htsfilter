#' Parameter initialization for a Poisson mixture model.
#' 
#' This function implements a K-means initialization
#' strategy (\code{kmeanInit}) for Poisson mixture models 
#' that may in turn be used to initialize the small EM
#' strategy. 
#' 
#' In practice, the user will not directly call the initialization functions
#' described here; they are indirectly called for a single number of clusters
#' through the \code{PoisMixClus_K} function (via \code{init.type}) or via the
#' \code{PoisMixClus} function for a sequence of cluster numbers (via
#' \code{Kmin.init} and \code{split.init}).
#' 
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
#' @importFrom stats kmeans
kmeanInit <- function(y, K, conds, norm, fixed.lambda,
                      equal.proportions) {
  
  n <- dim(y)[1];cols <- dim(y)[2];
  y <- as.matrix(y, nrow = n, ncol = cols)
  d <- length(unique(conds))
  r <- as.vector(table(conds))
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
  
  ## Use K-means to create initial partition
  g.init <- kmeans(y / w, K)
  partition <- g.init$cluster
  partition.mat <- matrix(0, nrow = n, ncol = K)
  for(i in 1:n) partition.mat[i,partition[i]] <- 1;
  
  ## Calculate lambda.init and p.init
  denom <- colSums(partition.mat * w)
  lambda.init <- matrix(NA, nrow = d, ncol = K)
  for(j in 1:d) {
    denom.bis <- denom * s.dot[j]
    num <- colSums(partition.mat *
                     rowSums(as.matrix(y[,which(conds == (unique(conds))[j])])))
    lambda.init[j,] <- num / denom.bis
  }
  if(class(fixed.lambda) == "list") {
    index <- lapply(fixed.lambda, function(x) {
      colSums((lambda.init - x)^2)})
    index.order <- lapply(index, order)
    fixed <- c()
    ## Pick which components correspond to fixed components
    ## Assign fixed lambdas to "closest" component (Euclidean distance),
    ## starting with the first fixed lambda
    for(ll in 1:length(fixed.lambda)) {
      while(index.order[[ll]][1] %in% fixed) {
        index.order[[ll]] <- index.order[[ll]][-1]
      }
      lambda.init[,index.order[[ll]][1]] <- fixed.lambda[[ll]]	
      fixed <- c(fixed, index.order[[ll]][1])
    }
    ## Rearrange lambda so that fixed values are in the first columns
    lambda.init <- cbind(lambda.init[,fixed], lambda.init[,-fixed])
  }
  if(equal.proportions == FALSE) {
    pi.init <- as.vector(table(g.init$cluster)/n)
    ## Rearrange pi so that it corresponds to lambda
    if(is.list(fixed.lambda) == TRUE) {
      pi.init <- c(pi.init[fixed], pi.init[-fixed])
    }
  }
  if(equal.proportions == TRUE) {
    pi.init <- rep(1/K, K)
  }
  return(list(pi.init = pi.init, lambda.init = lambda.init))
}
