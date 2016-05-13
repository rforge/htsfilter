#' Clustering high throughput sequencing (HTS) data
#' 
#' Mixture models are implemented to cluster genes from high-throughput
#' transcriptome sequencing (RNA-seq) data. Parameter estimation is performed
#' using either the EM algorithm, and model selection is performed using
#' either the slope heuristics or the integrated completed likelihood (ICL)
#' criterion. 
#' 
#' \tabular{ll}{ Package: \tab HTSCluster\cr Type: \tab Package\cr Version:
#' \tab 0.99.1\cr Date: \tab 2016-05-06\cr License: \tab GPL (>=3)\cr LazyLoad:
#' \tab yes\cr }
#' 
#' @name HTSClusterPlus-package
#' @aliases HTSClusterPlus-package HTSClusterPlus
#' @docType package
#' @author Andrea Rau, Gilles Celeux, Marie-Laure Martin-Magniette, Cathy
#' Maugis-Rabusseau
#' 
#' Maintainer: Andrea Rau <\url{andrea.rau@@jouy.inra.fr}>
#' @references Rau, A., Maugis-Rabusseau, C., Martin-Magniette, M.-L., Celeux,
#' G. (2015) Co-expression analysis of high-throughput transcriptome sequencing
#' data with Poisson mixture models. Bioinformatics, doi:
#' 10.1093/bioinformatics/btu845.
#' 
#' Rau, A., Celeux, G., Martin-Magniette, M.-L., Maugis-Rabusseau, C. (2011)
#' Clustering high-throughput sequencing data with Poisson mixture models.
#' Inria Research Report 7786. Available at
#' \url{http://hal.inria.fr/inria-00638082}.
#' @keywords models cluster
#' @example /inst/examples/HTSClusterPlus-package.R
#' @importFrom utils winMenuAddItem
NULL


