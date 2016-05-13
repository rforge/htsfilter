#' Title
#'
#' @param y (\emph{n} x \emph{q}) matrix of observed counts for \emph{n}
#' observations and \emph{q} variables
#' @param K Number of clusters (a single value). If \code{fixed.lambda}
#' contains a list of lambda values to be fixed, \code{K} corresponds to the
#' number of clusters in addition to those fixed.
#' @param conds Vector of length \emph{q} defining the condition (treatment
#' group) for each variable (column) in \code{y} 
#' @param norm The type of estimator to be used to normalize for differences in
#' library size: (\dQuote{\code{TC}} for total count, \dQuote{\code{UQ}} for
#' upper quantile, \dQuote{\code{Med}} for median, \dQuote{\code{DESeq}} for
#' the normalization method in the DESeq package, and \dQuote{\code{TMM}} for
#' the TMM normalization method (Robinson and Oshlack, 2010). Can also be a
#' vector (of length \emph{q}) containing pre-estimated library size estimates
#' for each sample. 
#' @param model Type of mixture model to use (\dQuote{\code{Poisson}} or \dQuote{\code{Normal}})
#' @param transformation Transformation type to be used: \dQuote{\code{voom}}, \dQuote{\code{logRPKM}}
#' (if \code{geneLength} is provided by user), \dQuote{\code{arcsin}}, \dQuote{\code{logit}},
#' \dQuote{\code{logMedianRef}}, \dQuote{\code{profile}}, \dQuote{\code{clrProfile}}, 
#' \dQuote{\code{log}}, \dQuote{\code{normlog}}, \dQuote{\code{logMeanRef}}, \dQuote{\code{logGMeanRef}},
#' \dQuote{\code{vst}}, \dQuote{\code{moderatedCPM}}, \dQuote{\code{none}}
#' @param subset.index Optional vector providing the indices of a subset of
#' genes that should be used for the co-expression analysis (i.e., row indices
#' of the data matrix \code{y}.
#' @param filterMean If \code{TRUE}, filter genes with normalized mean less than
#' the user-provided value in \dQuote{\code{filterCutoff}}
#' @param filterCutoff Value used to filter low mean normalized counts when 
#' \dQuote{\code{filterMean}} equals \code{TRUE}.
#' @param geneIDs Optional vector of gene ID's.
#' @param parallel If \code{FALSE}, no parallelization. If \code{TRUE}, parallel 
#' execution using BiocParallel (see next argument \code{BPPARAM}). A note on running 
#' in parallel using BiocParallel: it may be advantageous to remove large, unneeded objects 
#' from the current R environment before calling the function, as it is possible that R's 
#' internal garbage collection will copy these files while running on worker nodes.
#' @param BPPARAM Optional parameter object passed internally to \code{bplapply} when 
#' \code{parallel=TRUE}. If not specified, the parameters last registered with \code{register} 
#' will be used. 
#' @param ... Additional optional parameters to be passed to \code{\link{NormMixClus}} or 
#' \code{\link{PoisMixClus}}.
#'
#' @return
#' Object of class \code{HTSCluster} containing the full \code{results} (itself an object
#' of class \code{PoisMixClus} or \code{NormMixClus}), the \code{model} used, and if
#' applicable, the data \code{transformation} used.
#' @export
#'
#'
#' 
HTSCluster <- function(y, K, conds=NULL, norm="TMM", model="Normal", transformation="none", 
                       subset.index=NA, filterMean=FALSE, filterCutoff=50,
                       geneIDs=1:nrow(y), parallel=TRUE, BPPARAM=bpparam(), ...) {
  
  ## Parse ellipsis function
#  arg.user <- as.list(substitute(list(...)))[-1L]
  arg.user <- list()

  ## Optional parameters for PoisMixClus
  if(model == "Poisson") {
    if(is.null(arg.user$Kmin.init)) arg.user$Kmin.init<-"small-em";
    if(is.null(arg.user$split.init)) arg.user$split.init<-FALSE;
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
  }

  ## Optional parameters for transform_RNAseq
  if(is.null(arg.user$geneLength)) arg.user$geneLength<-NA;
  
  ## Optional parameters for NormMixClus
  if(model == "Normal") {
    if(is.null(arg.user$alg.type)) arg.user$alg.type<-"EM";
    if(is.null(arg.user$init.runs)) arg.user$init.runs<-50;
    if(is.null(arg.user$init.type)) arg.user$init.type<-"small-em";
    if(is.null(arg.user$init.iter)) arg.user$init.iter<-20;
    if(is.null(arg.user$iter)) arg.user$iter<-1000;
    if(is.null(arg.user$cutoff)) arg.user$cutoff<-0.001;
  }

  ########################
  ## POISSON MIXTURE MODEL
  ########################
  if(length(model) == 1 & model == "Poisson") {
    if(transformation != "none") stop("Poisson mixture model may only be applied on raw counts.")
    if(is.null(conds)) {
      message("Poisson mixture model fit assuming each sample is an independent condition.")
      conds <- 1:ncol(y)
    }
    tcounts <- transform_RNAseq(y, norm=norm, transformation="none", 
                                geneLength=arg.user$geneLength, filterMean=filterMean, 
                                filterCutoff=filterCutoff)

    run <- PoisMixClus(y=tcounts$tcounts, K=K, conds=conds, norm=norm, parallel=parallel, 
                       Kmin.init=arg.user$Kmin.init,
                       split.init=arg.user$split.init, subset.index=subset.index, 
                       BPPARAM=BPPARAM,
                       init.runs=arg.user$init.runs, init.iter=arg.user$init.iter, 
                       alg.type=arg.user$alg.type, cutoff=arg.user$cutoff, 
                       iter=arg.user$iter,
                       fixed.lambda=arg.user$fixed.lambda, 
                       equal.proportions=arg.user$equal.proportions,
                       verbose=arg.user$verbose, EM.verbose=arg.user$EM.verbose, 
                       interpretation=arg.user$interpretation)
  }
  
  ########################
  ## NORMAL MIXTURE MODEL
  ########################
  if(length(model) == 1 & model == "Normal") {
    
    tcounts <- transform_RNAseq(y, norm=norm, transformation=transformation, 
                                geneLength=arg.user$geneLength, filterMean=filterMean, 
                                filterCutoff=filterCutoff)
    run <- NormMixClus(y_profiles=tcounts$tcounts, K=K, subset.index=subset.index, parallel=parallel,
                       BPPARAM=BPPARAM, alg.type=arg.user$alg.type, init.runs=arg.user$init.runs,
                       init.type=arg.user$init.type, init.iter=arg.user$init.iter, 
                       iter=arg.user$iter, cutoff=arg.user$cutoff)
  }

  ####################################
  ## POISSON AND NORMAL MIXTURE MODELS? To be added later...
  ####################################
  # if("Normal" %in% model & "Poisson" %in% model) {
  #    
  # }
  
  
  ####################################
  ## RETURN RESULTS
  ####################################
  
  
  RESULTS <- list(results = run, model=model, transformation=transformation)
  class(RESULTS) <- "HTSCluster"
  return(RESULTS)
}