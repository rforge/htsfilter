#' Visualize results from coseq clustering
#' 
#' Plot a coseq object.
#'
#' @param x An object of class \code{"coseq"}
#' @param y_profiles y (\emph{n} x \emph{q}) matrix of observed profiles for \emph{n}
#' observations and \emph{q} variables to be used for graphing results (optional for
#' \code{logLike}, \code{ICL}, \code{probapost_boxplots}, and \code{probapost_barplots},
#' and by default takes value \code{x$tcounts} if \code{NULL})
#' @param K If desired, the specific model to use for plotting. If \code{NULL},
#' the model chosen by ICL will be plotted
#' @param threshold Threshold used for maximum conditional probability; only observations
#' with maximum conditional probability greater than this threshold are visualized 
#' @param conds Condition labels, if desired
#' @param average_over_conds If \code{TRUE}, average values of \code{y_profiles} within
#' each condition identified by \code{conds} for the \code{profiles} and \code{boxplots}
#' plots
#' @param graphs Graphs to be produced, one (or more) of the following: 
#' \code{"logLike"} (log-likelihood plotted versus number of clusters),
#' \code{"ICL"} (ICL plotted versus number of clusters), 
#' \code{"profiles"} (line plots of profiles in each cluster), \code{"boxplots"} 
#' (boxplots of profiles in each cluster), \code{"probapost_boxplots"} (boxplots of
#' maximum conditional probabilities per cluster), \code{"probapost_barplots"} 
#' (number of observations with a maximum conditional probability greater than 
#' \code{threshold} per cluster), \code{"probapost_histogram"} (histogram of maximum
#' conditional probabilities over all clusters) ...
#' @param order If \code{TRUE}, order clusters in \code{probapost_boxplot} by median and
#' \code{probapost_barplot} by number of observations with maximum conditional probability
#' greater than \code{threshold}
#' @param ...  Additional optional plotting arguments
#' 
#' @author Andrea Rau, Cathy Maugis-Rabusseau
#'
#' @export
## TODO: COMPARISON OF DIFFERENT TRANSFORMATIONS?
plot.coseq <- function(x, y_profiles=NULL, K=NULL, threshold=0.8, conds=NULL,
                             average_over_conds=FALSE, 
                             graphs=c("logLike", "ICL", 
                                      "profiles", "boxplots", "probapost_boxplots",
                                      "probapost_barplots", "probapost_histogram"), 
                             order=FALSE, ...) {
  
  if(is.null(y_profiles) == TRUE) y_profiles <- x$tcounts
  
  plot(x$results, y_profiles=y_profiles, K=K, threshold=threshold, conds=conds,
       average_over_conds=average_over_conds, graphs=graphs, 
       order=order, ...)
  
}

# plot.coseq <- function(resarcsin, reslogit, reslogMedianRef, y_profiles) {
#   n=dim(y_profiles)[1]
#   p=dim(y_profiles)[2]
#   qarcsin<-(n*p*log(2)) + (0.5*sum(sum(log(y_profiles*(1-y_profiles)))))
#   qlogit<-(n*p*log(log(2))) + (sum(sum(log(y_profiles*(1-y_profiles)))))
#   qlogmedianref<-(n*p*log(log(2))) + sum(sum(log(y_profiles)))
#   
#   plot(resarcsin$nbClust,resarcsin$ICLvalue + (2*qarcsin),col="red",type="l",
#        xlab="nb cluster",ylab="ICL sur les pi")
#   points(reslogit$nbClust,reslogit$ICLvalue + (2*qlogit),col="blue",type="l")
#   points(reslogMedianRef$nbClust, reslogMedianRef$ICLvalue + (2*qlogmedianref),
#          col="magenta",type="l")
#   legend("topright",legend=c("arcsin","logit","logMedianRef"),
#          col=c("red","blue","magenta"),lty=1)
# }