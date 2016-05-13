#' Summarize results from clustering using a Poisson mixture model
#' 
#' A function to summarize the clustering results obtained from a Poisson
#' mixture model.
#' 
#' The summary function for an object of class \code{"PoisMixClus"}
#' provides the number of clusters selected for the BIC, ICL, DDSE, and Djump
#' model selection approaches.
#' 
#' @param object An object of class \code{"PoisMixClus"} 
#' @param ... Additional arguments
#' 
#' @author Andrea Rau
#' @seealso \code{\link{PoisMixClus}}, \code{\link{PoisMixClus_K}}
#' @references Rau, A., Maugis-Rabusseau, C., Martin-Magniette, M.-L., Celeux,
#' G. (2015) Co-expression analysis of high-throughput transcriptome sequencing
#' data with Poisson mixture models. Bioinformatics, doi:
#' 10.1093/bioinformatics/btu845.
#' 
#' Rau, A., Celeux, G., Martin-Magniette, M.-L., Maugis-Rabusseau, C. (2011).
#' Clustering high-throughput sequencing data with Poisson mixture models.
#' Inria Research Report 7786. Available at
#' \url{http://hal.inria.fr/inria-00638082}.
#' @keywords methods
#' @example /inst/examples/PoisMixClus.R
#' @export
summary.PoisMixClus <-
  function (object, ...) 
  {
    x <- object
    if (class(x) != "PoisMixClus") {
      stop(paste(sQuote("x"), sep = ""), " must be of class ", 
           paste(dQuote("PoisMixClus"), sep = ""), sep = "")
    }
    cat("*************************************************\n")
    cat("Selected number of clusters via ICL = ", ncol(x$ICL.results$lambda), "\n", sep = "")
    cat("Selected number of clusters via BIC = ", ncol(x$BIC.results$lambda), "\n", sep = "")
    if(is.na(x$Djump.results) == FALSE) {
      cat("Selected number of clusters via Djump = ", ncol(x$Djump.results$lambda), "\n", sep = "")
    }
    if(is.na(x$DDSE.results) == FALSE) {
      cat("Selected number of clusters via DDSE = ", ncol(x$DDSE.results$lambda), "\n", sep = "")
    }
    if(is.na(x$Djump.results) == TRUE) {
      cat("Djump results not available \n", sep = "")
    }
    if(is.na(x$DDSE.results) == TRUE) {
      cat("DDSE results not available \n", sep = "")
    }
    cat("*************************************************\n")
  }




