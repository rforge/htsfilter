HTSDiff <- function(counts, conds, DEclusters=4, norm="TMM", epsilon=0.8, EM.verbose=FALSE, ...)
{
  .x <-  as.list(substitute(list(...)))[-1L]
  counts <- as.matrix(counts)
  conds <- as.vector(conds)
  if(length(unique(conds)) != 2) {
    stop("The number of unique conditions must be 2.")
  }
  if(min(rowSums(counts)) == 0) {
    message(paste("Note that", length(min(rowSums(counts))),
                  "genes with 0 counts in all samples were removed prior to the analysis."))
    counts <- counts[-which(rowSums(counts) == 0),]
  }

  PMM.args <- c(list(y=counts, g=DEclusters, conds=conds, lib.size=TRUE, lib.type=norm,
				   EM.verbose = EM.verbose,
                   fixed.lambda=list(rep(1,length(unique(conds))))), .x)
  DE.PMM <- do.call(PoisMixClus, PMM.args)
  
  lambda <- DE.PMM$lambda
  pi <- DE.PMM$pi
  probaPost <- DE.PMM$probaPost
  colnames(probaPost) <- c("NDE", paste("DE", 1:DEclusters, sep="."))
  MAPlabels <- DE.PMM$labels
  ICL <- DE.PMM$ICL
  s <- DE.PMM$s
  
  ####################
  ## Identify DE genes
  ####################
  ## First identify NDE clusters using epsilon cutoff
  index<-which(abs(log2(lambda[1,])-log2(lambda[2,]))<=epsilon)        
  if(length(index)!=1) {
    probaNDE<-apply(matrix(probaPost[,index], nrow=nrow(counts)),1,sum)
    probaDE <- apply(matrix(probaPost[,-index], nrow=nrow(counts)), 1, sum)
  }
  if(length(index)==1) {
    probaNDE<-probaPost[,index];
    probaDE <- apply(matrix(probaPost[,-index], nrow=nrow(counts)), 1, sum)
  }
  ## Added September 11, 2014: all genes are NDE if all clusters have too small of a
  ## ratio between lambdas
  if(length(index) == DEclusters + 1) {
	probaNDE <- rep(1, length(probaNDE))
  }
  DE <- ifelse(probaNDE<=1e-8, TRUE, FALSE)
  
  ###############
  ## SAVE RESULTS
  ###############
  ## Fixed error in ID names: October 31, 2013
  if(is.null(rownames(counts)) == TRUE) {
    id <- 1:nrow(counts)
  }
  else {
    id <- rownames(counts)
  }
  ## Normalized baseMean, baseMeanA, and baseMeanB
  normCounts <- t(t(counts) / (s*length(conds)))
  baseMean <- rowMeans(normCounts)
  baseMeanA <- rowMeans(normCounts[,which(conds == unique(conds)[1])])
  baseMeanB <- rowMeans(normCounts[,which(conds == unique(conds)[2])])
  foldChange <- baseMeanB/baseMeanA
  log2FoldChange <- log2(foldChange)
  probaPostNDE <- probaPost[,1]
  
  res <- data.frame(id, baseMean, baseMeanA, baseMeanB, foldChange,
                    log2FoldChange, tauDE = probaDE, tauNDE = probaNDE, DE = DE)
  rownames(res) <- NULL
  
  return(list(res=res, PMM=DE.PMM, iterations=DE.PMM$iterations,
			  logLikeDiff=DE.PMM$logLikeDiff))
}

