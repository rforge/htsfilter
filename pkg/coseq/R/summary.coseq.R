#' Summarize results from clustering using a Poisson or Gaussian mixture model
#' 
#' A function to summarize the clustering results obtained from a Poisson or Gaussian
#' mixture model.
#' 
#' 
#' @param object An object of class \code{"coseq"}
#' @param ... Additional arguments
#' @author Andrea Rau
#' @seealso \code{\link{coseq}}
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
#' @example /inst/examples/coseq-package.R
#' @export
summary.coseq <-
  function (object,  ...) 
  {
    x <- object
    if (class(x) != "coseq") {
      stop(paste(sQuote("object"), sep = ""), " must be of class ", 
           paste(dQuote("coseq"), sep = ""), sep = "")
    }
    cat("*************************************************\n")
    cat("Model: ", x$model, "\n", sep = "")
    cat("Transformation: ", x$transformation, "\n", sep = "")

    if(class(x$results) == "NormMixClus") summary(object=x$results, y_profiles=x$tcounts, ...)
    if(class(x$results) == "HTSClusterWrapper") summary(x$results, ...)
  }




