#' Summarize results from clustering using a Poisson mixture model
#' 
#' A function to summarize the clustering results obtained from a Poisson
#' mixture model.
#' 
#' The summary function for an object of class \code{"PoisMixClus_K"} provides the
#' following summary of results:
#' 
#' 1) Number of clusters and model selection criterion used, if applicable.
#' 
#' 2) Number of observations across all clusters with a maximum conditional
#' probability greater than 90% (and corresponding percentage of total
#' observations) for the selected model.
#' 
#' 3) Number of observations per cluster with a maximum conditional probability
#' greater than 90% (and corresponding percentage of total observations per
#' cluster) for the selected model.
#' 
#' 4) \eqn{\ensuremath\boldsymbol{\lambda}}{\lambda} values for the selected
#' model.
#' 
#' 5) \eqn{\ensuremath\boldsymbol{\pi}}{\pi} values for the selected model.
#' 
#' @param object An object of class \code{"PoisMixClus_K"}
#' @param ... Additional arguments
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
`summary.PoisMixClus_K` <-
function (object, ...) 
{
	x <- object
    	if (class(x) != "PoisMixClus_K") {
        	stop(paste(sQuote("x"), sep = ""), " must be of class ", 
            paste(dQuote("PoisMixClus_K"), sep = ""), sep = "")
    	}

	probaPost <- x$probaPost
	labels <- x$labels
	lambda <- x$lambda
	pi <- x$pi
	g <- length(pi)

	map <- apply(probaPost, 1, max)
	length(which(map > 0.9))/length(map)

	cat("*************************************************\n")
	cat("Number of clusters = ", g, "\n", sep = "")
	if(is.null(x$model.selection) == FALSE) {
		cat("Model selection via ", x$model.selection, "\n", sep = "")
	}
	cat("*************************************************\n")
	tab <- table(labels)
	names(tab) <- paste("Cluster", names(tab))
	cat("Cluster sizes:\n"); print(tab); cat("\n")
	cat("Number of observations with MAP > 0.90 (% of total):\n")
	cat(length(which(map > 0.9)), " (", round(length(which(map > 0.9))/length(map)*100,2),
		"%)\n\n", sep = "")
	cat("Number of observations with MAP > 0.90 per cluster (% of total per cluster):\n"); 

	tab2 <- matrix(NA, nrow = 2, ncol = g)
	colnames(tab2) <- paste("Cluster", 1:g); rownames(tab2) <- rep("", 2)
	for(i in 1:g) {
		if(sum(labels == i) > 1) {
			map.clust <- apply(matrix(probaPost[labels == i,], ncol=g), 1, max)
			tab2[1,i] <- length(which(map.clust > 0.9))
			tab2[2,i] <- paste("(", round(100*length(which(map.clust > 0.9))/length(map.clust),2),
				"%)", sep = "")
		}
		if(sum(labels == i) == 1) {
			map.clust <- max(probaPost[labels == i,])
			tab2[1,i] <- length(which(map.clust > 0.9))
			tab2[2,i] <- paste("(", round(100*length(which(map.clust > 0.9))/length(map.clust),2),
				"%)", sep = "")
		}
		if(sum(labels == i) == 0) {
			tab2[1,i] <- "---"
			tab2[2,i] <- "---"
		}
	}
	print(tab2, quote = FALSE); cat("\n")

	cat("Lambda:\n"); print(round(lambda,2)); cat("\n")
	cat("Pi:\n"); print(round(pi,2)); cat("\n")
}
