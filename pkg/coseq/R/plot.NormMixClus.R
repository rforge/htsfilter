#' Visualize results from clustering using a Normal mixture model
#'
#' @param x An object of class \code{"NormMixClus"}
#' @param threshold Minimum threshold for conditional probabilities used in graphing
#' @param order If \code{TRUE}, orders boxplots by the median condiitonal probability value
#' @param graphs Type of graph to be included in plots. May be equal to \code{c("ICL",
#' "BIC")} for objects of class \code{"NormMixClus"}
#' @param ... Additional arguments
#' 
#' @author Cathy Maugis-Rabusseau
#'
#' @export
#' @importFrom graphics plot
#' @importFrom graphics points
#' @importFrom graphics boxplot
#' @importFrom graphics axis
#' @importFrom graphics barplot
plot.NormMixClus <- function(x, graphs=c("logLike", "ICL", "boxplots", "barplots"), 
                             threshold=0.8, order=FALSE, ...) {

  ## graphe de la logvraisemblance
  if("logLike" %in% graphs) {
    plot(x$nbClust.all, x$logLike.all, type="l",xlab="nbCluster",ylab="loglikelihood",main="")
    points(x$nbClust.all, x$logLike.all, pch=20)
  }
  
  ## graphe de ICL
  if("ICL" %in% graphs) {
    plot(x$nbClust.all, x$ICL.all,type="l",xlab="nbCluster",ylab="ICL",main="")
    points(x$nbClust.all, x$ICL.all,pch=20)
  }

  ## boxplot des probapost
  if("boxplots" %in% graphs) {
    probapost <- x$ICL.results$probaPost
    label <- apply(probapost,1,which.max)
    A <- boxplot(apply(probapost,1,max)~label, plot=F)
    if(order==TRUE){
      J <- sort.int(A$stat[3,],index.return=T,decreasing=T)$ix   #ordre par rapport à la mediane
    } else{
      J<-seq(1,ncol(probapost),1)
    }  
    boxplot(A$stat[,J],axes=F,outline=T)
    axis(side=1, at=seq(1,ncol(probapost)), labels=J, cex.lab=0.5)
    axis(side=2)
  }
  
  ## graphe barplot de la proportion de probapost >seuil et < seuil par classe
  if("barplots" %in% graphs) {
    A<-matrix(0, nrow=2, ncol=ncol(probapost))
    for (k in 1:ncol(probapost)){
      I<-which(label==k)
      A[,k]<-c(sum(probapost[I,k]>threshold), sum(probapost[I,k]<=threshold))
    }
    barplot(A[,J],names.arg=J)
  }

}
