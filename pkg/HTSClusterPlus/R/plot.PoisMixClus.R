#' Visualize results from clustering using a Poisson mixture model
#' 
#' A function to visualize the clustering results obtained from a Poisson
#' mixture model.
#' 
#' For objects of class \code{"PoisMixClus"}, the plotting function
#' provides the possibility for one or all of the following visualizations:
#' 
#' 1) ICL plot for all fitted models.
#' 
#' 2) BIC plot for all fitted models.
#' 
#' 5) Capushe diagnostic plots.
#' 
#' @param x An object of class \code{"PoisMixClus"}
#' @param file.name Optional file name if plots are to be saved in a PDF file.
#' @param graphs Type of graph to be included in plots. May be equal to \code{c("ICL",
#' "BIC")} for objects of class \code{"PoisMixClus"}
#' @param capushe.validation Optional number of clusters to use for capushe
#' validation (should be less than the maximum number of clusters specified in
#' the \code{"PoisMixClus"} object).
#' @param ... Additional arguments
#' @author Andrea Rau
#' @seealso \code{\link{PoisMixClus}}, \code{\link{PoisMixClus_K}}
#' @references
#' 
#' Rau, A., Maugis-Rabusseau, C., Martin-Magniette, M.-L., Celeux, G. (2015)
#' Co-expression analysis of high-throughput transcriptome sequencing data with
#' Poisson mixture models. Bioinformatics, doi: 10.1093/bioinformatics/btu845.
#' 
#' Andrea Rau, Gilles Celeux, Marie-Laure Martin-Magniette, and Cathy
#' Maugis-Rabusseau (2011).  Clustering high-throughput sequencing data with
#' Poisson mixture models. \emph{Technical report} RR-7786, Inria Saclay --
#' Ile-de-France.
#' @keywords methods models
#' @example /inst/examples/plot.PoisMixClus.R
#' 
#' @export
#' @importFrom capushe validation
#' @importFrom capushe capushe
#' @importFrom grDevices pdf
#' @importFrom grDevices dev.off
#' @importFrom graphics lines
#' @importFrom graphics points
#' @importFrom graphics par

plot.PoisMixClus <-
  function (x, file.name = FALSE, 
            graphs = c("capushe", "ICL", "BIC"), capushe.validation=NA, ...) 
  {	
    if (class(x) != "PoisMixClus") {
      stop(paste(sQuote("x"), sep = ""), " must be of class ", 
           paste(dQuote("PoisMixClus"), sep = ""), sep = "")
    }
    
    if(file.name != FALSE) pdf(paste(file.name));
    
    if("ICL" %in% graphs & "BIC" %in% graphs) {
      par(mfrow = c(1,2), mar = c(4,4,2,2))
      gpl <- unlist(lapply(strsplit(names(x$ICL.all), split="=", fixed=TRUE), function(xx) xx[2]))
      plot(gpl, x$ICL.all, xlab = "Number of clusters", ylab = "ICL", main="ICL", pch=19)
      lines(gpl, x$ICL.all, lwd=2)
      points(gpl[which(x$ICL.all == max(x$ICL.all))], max(x$ICL.all), col="red", pch="X", font=2, cex=2)
      plot(gpl, unlist(lapply(x$all.results, function(xx) xx$BIC)), xlab = "Number of clusters", ylab = "BIC",
           main="BIC", pch=19)
      lines(gpl, unlist(lapply(x$all.results, function(xx) xx$BIC)), lwd=2)          
    }
    
    if("ICL" %in% graphs & !"BIC" %in% graphs) {
      par(mar = c(4,4,2,2))
      gpl <- unlist(lapply(strsplit(names(x$ICL.all), split="=", fixed=TRUE), function(xx) xx[2]))
      plot(gpl, x$ICL.all, xlab = "Number of clusters", ylab = "ICL", main="ICL", pch=19)
      lines(gpl, x$ICL.all, lwd=2)
      points(gpl[which(x$ICL.all == max(x$ICL.all))], max(x$ICL.all), col="red", pch="X", font=2, cex=2)
    }
    
    if(!"ICL" %in% graphs & "BIC" %in% graphs) {
      par(mar = c(4,4,2,2))
      gpl <- unlist(lapply(strsplit(names(x$ICL.all), split="=", fixed=TRUE), function(xx) xx[2]))
      plot(gpl, unlist(lapply(x$all.results, function(xx) xx$BIC)), xlab = "Number of clusters", 
           ylab = "BIC",
           main="BIC", pch=19)
      lines(gpl, unlist(lapply(x$all.results, function(xx) xx$BIC)), lwd=2)          
    }
    
    if("capushe" %in% graphs) {
      plot(x$capushe, newwindow=FALSE)
      if(is.na(capushe.validation) == FALSE) {
        Kchoice <- as.numeric(unlist(lapply(strsplit(names(x$logLike), "="), function(x) x[2])))
        if(capushe.validation >= max(Kchoice)) 
          stop("Number of clusters for capushe validation should be less than largest number of clusters.");
        np <- (Kchoice-1) + (length(unique(x$all.result[[1]]$conds))-1)*(Kchoice)
        n <- nrow(x$all.result[[1]]$probaPost)
        mat <- cbind(Kchoice, np/n, np/n, -x$logLike.all)
        ResCapushe <- capushe(mat, n)
        validation(ResCapushe, mat[-c(which(Kchoice < capushe.validation)),], 
                   newwindow=FALSE)
      }
    }
    
    if(file.name != FALSE)  dev.off();
  }


