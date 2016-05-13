#' Poisson mixture model estimation and selection for a series of cluster numbers
#' 
#' This function implements the EM and CEM algorithms for parameter
#' estimation in a Poisson mixture model for clustering high throughput
#' sequencing observations (e.g., genes) over a sequence of numbers of clusters.
#' Parameters are initialized using a Small-EM strategy as described in Rau et
#' al. (2011) or the splitting small-EM strategy described in Papastamoulis et
#' al. (2014), and model selection is performed using the BIC/ICL criteria or
#' the slope heuristics.
#' 
#' Output of \code{PoisMixClus} is an S3 object of class \code{PoisMixClus}.
#' 
#' In a Poisson mixture model, the data \eqn{\mathbf{y}}{y} are assumed to come
#' from \emph{K} distinct subpopulations (clusters), each of which is modeled
#' separately; the overall population is thus a mixture of these
#' subpopulations. In the case of a Poisson mixture model with \emph{K}
#' components, the model may be written as
#' \deqn{f(\mathbf{y};K,\ensuremath\boldsymbol{\Psi}_K) = \prod_{i=1}^n
#' \sum_{k=1}^K \pi_k \prod_{j=1}^{d}\prod_{l=1}^{r_j} P(y_{ijl} ;
#' \ensuremath\boldsymbol{\theta}_k)}{f(y;K,\psi_K) = \prod_{i=1}^n
#' \sum_{k=1}^K \pi_k \prod_{j=1}^{d}\prod_{l=1}^{r_j} P(y_{ijl} ; \theta_k)}
#' for \eqn{i = 1, \ldots, n} observations in \eqn{l = 1, \ldots, r_j}
#' replicates of \eqn{j = 1, \ldots, d} conditions (treatment groups), where
#' \eqn{P(\cdot)} is the standard Poisson density,
#' \eqn{\ensuremath\boldsymbol{\Psi}_K = (\pi_1,\ldots,\pi_{K-1},
#' \ensuremath\boldsymbol{\theta}^\prime)}{\psi_K = (\pi_1,\ldots,\pi_{K-1},
#' \theta^\prime)}, \eqn{\ensuremath\boldsymbol{\theta}^\prime}{\theta^\prime}
#' contains all of the parameters in
#' \eqn{\ensuremath\boldsymbol{\theta}_1,\ldots,\ensuremath\boldsymbol{\theta}_K}{\theta_1,\ldots,\theta_K}
#' assumed to be distinct, and \eqn{\ensuremath\boldsymbol{\pi} =
#' (\pi_1,\ldots,\pi_K)^\prime}{\pi = (\pi_1,\ldots,\pi_K)^\prime} are the
#' mixing proportions such that \eqn{\pi_k} is in (0,1) for all \emph{k} and
#' \eqn{\sum_k \pi_k = 1}. We consider the following parameterization for the mean
#' \eqn{\ensuremath\boldsymbol{\theta}_k = (\mu_{ijlk})}{\theta = (mu_{ijlk})}.
#' We consider \deqn{\mu_{ijlk} = w_i s_{jl} \lambda_{jk}} where \eqn{w_i}
#' corresponds to the expression level of observation \emph{i},
#' \eqn{\ensuremath\boldsymbol{\lambda}_k =
#' (\lambda_{1k},\ldots,\lambda_{dk})}{\lambda_k =
#' (\lambda_{1k},\ldots,\lambda_{dk})} corresponds to the clustering parameters
#' that define the profiles of the genes in cluster \emph{k} across all
#' variables, and \eqn{s_{jl}} is the normalized library size (a fixed
#' constant) for replicate \emph{l} of condition \emph{j}.
#' 
#' There are two approaches to estimating the parameters of a finite mixture
#' model and obtaining a clustering of the data: the estimation approach (via
#' the EM algorithm) and the clustering approach (via the CEM algorithm).
#' Parameter initialization is done using a Small-EM strategy as described in
#' Rau et al. (2011) via the \code{\link{emInit}} function. Model selection may
#' be performed using the BIC or ICL criteria, or the slope heuristics.
#' 
#' @inheritParams PoisMixClus_K
#' 
#' @param K Number of clusters (a single value or a sequence of values).
#' @param Kmin.init Type of initialization strategy to be used for the
#' minimum number of clusters in a sequence:
#' (\dQuote{\code{small-em}} for the Small-EM strategy described in Rau et al.
#' (2011), and \dQuote{\code{kmeans}} for a simple \emph{K}-means
#' initialization)
#' @param split.init If \code{TRUE}, the splitting initialization strategy of
#' Papastamoulis et al. (2014) will be used for cluster sizes larger than the minimum number, 
#' the initialization strategy specified in
#' \code{Kmin.init.type} is used otherwise.
#' @param parallel If \code{FALSE}, no parallelization. If \code{TRUE}, parallel 
#' execution using BiocParallel (see next argument \code{BPPARAM}). A note on running 
#' in parallel using BiocParallel: it may be advantageous to remove large, unneeded objects 
#' from the current R environment before calling the function, as it is possible that R's 
#' internal garbage collection will copy these files while running on worker nodes.
#' @param BPPARAM Optional parameter object passed internally to \code{bplapply} when 
#' \code{parallel=TRUE}. If not specified, the parameters last registered with \code{register} 
#' will be used.
#' @param ... Additional optional parameters to be passed to \code{\link{PoisMixClus_K}}.
#' 
#' @return 
#' \item{loglike.all }{Log likelihoods calculated for each of the fitted models}
#' \item{capushe }{Results of
#' capushe model selection, an object of class \code{"Capushe"}} 
#' \item{ICL.all
#' }{ICL values calculated for each of the fitted models}
#' \item{ICL.results }{Object of class
#' \code{PoisMixClus} giving the results from the model chosen via the ICL
#' criterion} 
#' \item{BIC.results }{Object of class \code{PoisMixClus} giving the
#' results from the model chosen via the BIC} 
#' \item{DDSE.results }{Object of
#' class \code{PoisMixClus} giving the results from the model chosen via the
#' DDSE slope heuristics criterion}
#'\item{Djump.results }{Object of class
#' \code{PoisMixClus} giving the results from the model chosen via the Djump
#' slope heuristics criterion} 
#' \item{all.results }{List of objects of class
#' \code{PoisMixClus} giving the results for all models}
#' \item{model.selection }{Type of criteria used
#' for model selection: \code{"DDSE"}, \code{"Djump"}, \code{"BIC"}, or
#' \code{"ICL"} }
#' 
#' @note Note that the \code{fixed.lambda} argument is primarily intended to be
#' used in the case when a single cluster is fixed to have equal clustering
#' parameters lambda across all conditions (i.e.,
#' \eqn{\lambda_{j1}=\lambda_{1}=1}); this is particularly useful when
#' identifying genes with non-differential expression across all conditions.. Alternatively, this
#' argument could be used to specify a cluster for which genes are only
#' expressed in a single condition (e.g., \eqn{\lambda_{11} = 1} and
#' \eqn{\lambda_{j1} = 0} for all \eqn{j > 1}). Other possibilities could be
#' considered, but note that the fixed values of lambda must satisfy the
#' constraint \eqn{\sum_j \lambda_{jk}s_{j.} = 1} for all \eqn{k} imposed in
#' the model; if this is not the case, a warning message will be printed.
#' 
#' @author Andrea Rau <\url{andrea.rau@@jouy.inra.fr}>
#' 
#' @seealso \code{\link{probaPost}} for the calculation of the conditional
#' probability of belonging to a cluster; \code{\link{PoisMixMean}} for the
#' calculation of the per-cluster conditional mean of each observation;
#' \code{\link{logLikePoisMixDiff}} for the calculation of the log likelihood
#' of a Poisson mixture model; \code{\link{emInit}} and \code{\link{kmeanInit}}
#' for the Small-EM parameter initialization strategy
#' 
#' @references Anders, S. and Huber, W. (2010) Differential expression analysis
#' for sequence count data. \emph{Genome Biology}, \bold{11}(R106), 1-28.
#' 
#' Papastamoulis, P., Martin-Magniette, M.-L., and Maugis-Rabusseau, C. (2014).
#' On the estimation of mixtures of Poisson regression models with large number
#' of components. \emph{Computational Statistics and Data Analysis}: 3rd
#' special Issue on Advances in Mixture Models, DOI:
#' 10.1016/j.csda.2014.07.005.
#' 
#' Rau, A., Maugis-Rabusseau, C., Martin-Magniette, M.-L., Celeux, G. (2015)
#' Co-expression analysis of high-throughput transcriptome sequencing data with
#' Poisson mixture models. Bioinformatics, doi: 10.1093/bioinformatics/btu845.
#' 
#' Rau, A., Celeux, G., Martin-Magniette, M.-L., Maugis-Rabusseau, C (2011).
#' Clustering high-throughput sequencing data with Poisson mixture models.
#' Inria Research Report 7786. Available at
#' \url{http://hal.inria.fr/inria-00638082}.
#' 
#' @keywords models cluster
#' @example /inst/examples/PoisMixClus.R
#' 
#' @importFrom BiocParallel bplapply
#' @importFrom BiocParallel bpparam
#' @importFrom capushe capushe
#' @importFrom capushe validation
#' @importFrom stats na.omit
#' @export 
#' 
PoisMixClus <- function(y, K, conds, norm="TMM", Kmin.init="small-em", split.init=FALSE, 
                        subset.index=NA, parallel=TRUE, BPPARAM = bpparam(), ...) 
{
  ## TODO: add check that K is a sequence of numbers
  
  if(is.matrix(y) == FALSE & is.data.frame(y) == FALSE) 
    stop(paste(sQuote("y"), "must be a matrix"))
  if(min(y) < 0 | sum(as.numeric(round(y))) != sum(as.numeric(y))) 
    stop(paste(sQuote("y"), "must be a matrix made up of nonnegative counts"))
  if(min(rowSums(y)) == 0)
    stop(paste("at least one observation in", sQuote("y"), 
               "contains all 0's and must be removed from the data"))
  if(is.vector(conds) == FALSE | length(conds) != ncol(y))
    stop(paste(sQuote("conds"), "must be a vector the same length as the number of columns in", sQuote("y")))
  if(Kmin.init != "small-em" & Kmin.init != "kmeans" & Kmin.init != "split.small-em") 
    stop(paste(sQuote("Kmin.init"), "must be one of", dQuote("small-em"), "or", 
               dQuote("kmeans"), "or", dQuote("split.small-em")))
  
  ## Parse ellipsis function
#  arg.user <- as.list(substitute(list(...)))[-1L]
  arg.user <- list(...)
  
  if(is.null(arg.user$init.runs)) arg.user$init.runs<-1;
  if(is.null(arg.user$init.iter)) arg.user$init.iter<-10;
  if(is.null(arg.user$alg.type)) arg.user$alg.type<-"EM";
  if(is.null(arg.user$cutoff)) arg.user$cutoff<-1e-05;
  if(is.null(arg.user$iter)) arg.user$iter<-1000;
  if(is.null(arg.user$fixed.lambda)) arg.user$fixed.lambda<-NA;
  if(is.null(arg.user$equal.proportions)) arg.user$equal.proportions<-FALSE;
  if(is.null(arg.user$verbose)) arg.user$verbose<-FALSE;
  if(is.null(arg.user$EM.verbose)) arg.user$EM.verbose<-FALSE;
  if(is.null(arg.user$interpretation)) arg.user$interpretation<-"sum";
  
  ## Grouping columns of y in order of condition (all replicates put together)
  o.ycols <- order(conds)
  y <- y[,o.ycols]
  conds <- conds[o.ycols]
  conds.names <- unique(conds)
  d <- length(unique(conds))
  r <- as.vector(table(conds))
  if(length(rownames(y)) == 0) rn <- 1:nrow(y);
  if(length(rownames(y)) > 0) rn <- rownames(y);
  y <- as.matrix(y, nrow = nrow(y), ncol = ncol(y))
  rownames(y) <- rn;
  
  if(length(norm) != 1 & length(norm) != length(conds)) 
    stop(paste(sQuote("norm"), "must be one of the following: none, TC, UQ, Med, DESeq, TMM, or a vector of the same length as", sQuote("conds")))
  ## If estimated from data, all genes should be used
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
  
  ## In case only a subset of data are to be used for analysis
  if(is.na(subset.index)[1] == FALSE) {
    y <- y[subset.index,]
    n <- dim(y)[1];cols <- dim(y)[2]
    w <- rowSums(y)
  }
  if(is.na(subset.index)[1] == TRUE) {
    n <- dim(y)[1];cols <- dim(y)[2]
    w <- rowSums(y)
  }
  
  all.results <- vector("list", length = length(K))
  names(all.results) <- paste("K=", K, sep = "")
  
  cat("Running K =", min(K), "...\n")
  all.results[[1]] <- PoisMixClus_K(y=y, K=min(K), conds=conds, norm=s, init.type=Kmin.init,  
                                          subset.index=NA, wrapper=TRUE, ...)
  
  index <- 2
  remainingK <- K[-which(K == min(K))]
  if(length(remainingK) > 0) {
    ## In the case where parallelization is NOT used
    if(!parallel) {
      for(k in remainingK) {
        cat("Running K =", k, "...\n")
        if(split.init == TRUE) {
          prev.labels <- all.results[[index-1]]$labels
          prev.probaPost <- all.results[[index-1]]$probaPost
          all.results[[index]] <- PoisMixClus_K(y=y, K=k, norm=s, conds = conds, 
                                                init.type = "split.small-em", subset.index=NA,
                                                wrapper=TRUE, prev.labels=prev.labels,
                                                prev.probaPost=prev.probaPost, ...)
        }
        if(split.init == FALSE) {
          all.results[[index]] <- PoisMixClus_K(y=y, K=k, conds=conds, norm=s, init.type=Kmin.init,  
                                                subset.index=NA, wrapper=TRUE, ...)
        }
        index <- index + 1
      }
    }  else if(parallel) {
    ## In the case where parallelization IS used
      if(split.init == TRUE) {
        warning("Splitting initialization is not compatible with parallelization.")
      }
      if(split.init == FALSE) {
        tmp <- bplapply(remainingK, function(ii) {
          cat("Running K =", ii, "...\n")
          res <- PoisMixClus_K(K=as.numeric(ii), y=y, conds=conds, norm=s, init.type="small-em", 
                                            subset.index=NA, wrapper=TRUE, prev.probaPost=NA, prev.labels=NA, 
                                            init.runs = arg.user$init.runs, init.iter = arg.user$init.iter, 
                                            alg.type = arg.user$alg.type, cutoff = arg.user$cutoff,
                                            iter = arg.user$iter, fixed.lambda = arg.user$fixed.lambda, 
                                            equal.proportions = arg.user$equal.proportions,
                                            verbose = arg.user$verbose,
                                            interpretation = arg.user$interpretation, 
                                            EM.verbose = arg.user$EM.verbose)
          return(res)}, BPPARAM=BPPARAM)
        Kmods <- paste("K=", unlist(lapply(tmp, function(x) ncol(x$lambda))), sep="")
        all.results[-1] <- tmp[na.omit(match(names(all.results), Kmods))]
      }
    }
  }

  logLike.all <- unlist(lapply(all.results, function(x) x$log.like))
  ICL.all <- unlist(lapply(all.results, function(x) x$ICL))
  ICL.choose <- which(ICL.all == max(ICL.all, na.rm = TRUE))
  select.results <- all.results[[ICL.choose]]
  select.results$model.selection <- "ICL"
  
  BIC.all <- unlist(lapply(all.results, function(x) x$BIC))
  BIC.choose <- which(BIC.all == max(BIC.all, na.rm = TRUE))
  select.results2 <- all.results[[BIC.choose]]
  select.results2$model.selection <- "BIC"
  
  # Apply capushe: only if at least 10 models are considered
  if(c(max(K) - min(K) + 1) <= 10) {
    message("Note: the slope heuristics approach for model selection may only be applied if more than 10 models are fit.")
    DDSE.results <- NA
    Djump.results <- NA
    capushe <- NA
    ResCapushe <- NA
  }
  if(c(max(K) - min(K) + 1) > 10) {
    message("Note: diagnostic plots for results corresponding to model selection via slope heuristics (Djump and DDSE) should be examined to ensure that sufficiently complex models have been considered.")
    Kchoice <- K
    np <- (Kchoice-1) + (length(unique(conds))-1)*(Kchoice)
    mat <- cbind(Kchoice, np/n, np/n, -logLike.all/n)
    ResCapushe <- suppressWarnings(capushe(mat, n))
    DDSE <- ResCapushe@DDSE@model
    Djump <- ResCapushe@Djump@model
    DDSE.results <- all.results[[paste("K=", DDSE, sep="")]]
    Djump.results <- all.results[[paste("K=", Djump, sep="")]]
    DDSE.results$model.selection <- "DDSE"
    Djump.results$model.selection <- "Djump"
  }
  
  RESULTS <- list(logLike.all = logLike.all, ICL.all = ICL.all,
                  capushe = ResCapushe, 
                  all.results = all.results,
                  DDSE.results = DDSE.results,
                  Djump.results = Djump.results,
                  ICL.results = select.results,
                  BIC.results = select.results2)
  class(RESULTS) <- "PoisMixClus"
  return(RESULTS)
}




