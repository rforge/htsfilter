#' Title
#'
#' @param x blah
#' @param parallel blah
#' @param BPPARAM blah
#' @param plot blah
#' @param ... Additional parameters for corrplot
#'
#' @return Matrix of adjusted rand index values calculated between each pair of models.
#' @export
#' @importFrom HTSCluster highDimensionARI
#' @importFrom corrplot corrplot
#'
compareARI <- function(x, parallel=FALSE, BPPARAM=bpparam(), plot=TRUE, ...) {
  
  arg.user <- list(...)
  if(is.null(arg.user$digits)) arg.user$digits<-2;
  
  ## For class coseq
  if(class(x) == "coseq") {
    full_labels <- do.call("cbind", lapply(x$results$all.results, 
                                           function(y) apply(y$probaPost,1,which.max)))
  }
  
  ## For class NormMixClus and PoisMixClusWrapper
  if(class(x) == "NormMixClus" | class(x) == "PoisMixClusWrapper") {
    full_labels <- do.call("cbind", lapply(x$all.results, 
                                           function(y) apply(y$probaPost,1,which.max)))
  }
  
  if(class(x) == "data.frame" | class(x) == "matrix") {
    full_labels <- x
    if(length(colnames(full_labels)) == 0) {
      colnames(full_labels) <- paste("Model", 1:ncol(full_labels));
    }
  }

  ARI <- matrix(0, nrow=ncol(full_labels), ncol=ncol(full_labels))
  rownames(ARI) <- colnames(ARI) <- colnames(full_labels)

  index <- ARI
  index[upper.tri(index)] <- seq(1, (ncol(ARI) * nrow(ARI) - ncol(ARI))/2)
  
  if(!parallel) {
    tmp <- lapply(1:max(index, na.rm=TRUE), function(ii) {
      index2 <- which(index == ii, arr.ind=TRUE)
      ARItmp <- highDimensionARI(full_labels[,index2[1]], full_labels[,index2[2]])
      return(list(ARItmp=ARItmp, ii=ii))
    })  
  } else if(parallel) {
    tmp <- bplapply(1:max(index, na.rm=TRUE), function(ii) {
      index2 <- which(index == ii, arr.ind=TRUE)
      ARItmp <- highDimensionARI(full_labels[,index2[1]], full_labels[,index2[2]])
      return(list(ARItmp=ARItmp, ii=ii))
    }, BPPARAM=BPPARAM) 
  }
  
  for(i in 1:length(tmp)) {
    new_index <- which(index == tmp[[i]]$ii, arr.ind=TRUE)
    ARI[new_index] <- tmp[[i]]$ARItmp
  }

  diag(ARI) <- 1
  
  if(plot == TRUE) {
    corrplot(ARI, is.corr=FALSE, method="color", type="upper", p.mat=ARI, insig="p-value", 
             sig.level=-1, tl.pos="d", addgrid.col="white",
             tl.col="black", ...)
  }
  
  ARI <- round(ARI, digits = arg.user$digits)
  ARI[lower.tri(ARI)] <- ""
  ARI <- data.frame(ARI, check.names=FALSE)
  return(ARI)
}