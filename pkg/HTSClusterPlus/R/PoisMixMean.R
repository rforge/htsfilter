#' Calculate the conditional per-cluster mean of each observation
#' 
#' This function is used to calculate the conditional per-cluster mean
#' expression for all observations. This value corresponds to
#' \eqn{\ensuremath\boldsymbol{\mu} = (\mu_{ijlk}) = (\hat{w}_i
#' \hat{\lambda}_{jk})}{\mu = (\mu_{ijlk}) = (\hat{w}_i \hat{\lambda}_{jk})}
#' for the PMM-I model and \eqn{\ensuremath\boldsymbol{\mu} = (\mu_{ijlk}) =
#' (\hat{w}_i s_{jl} \hat{\lambda}_{jk})}{\mu = (\mu_{ijlk}) = (\hat{w}_i
#' s_{jl}\hat{\lambda}_{jk})} for the PMM-II model.
#' 
#' @inheritParams PoisMixClus_K
#' @param lambda (\emph{d} x \code{K}) matrix containing the estimate of
#' \eqn{\hat{\ensuremath\boldsymbol{\lambda}}}{\hat{\lambda}} 
#' @param s Library size normalization factors
#' 
#' @return A list of length \code{K} containing the (\emph{n} x \emph{q})
#' matrices of mean expression for all observations, conditioned on each of the
#' \code{K} clusters
#' @author Andrea Rau <\url{andrea.rau@@jouy.inra.fr}>
#' @seealso \code{\link{PoisMixClus}} for Poisson mixture model estimation and
#' model selection
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
#' @example /inst/examples/PoisMixMean.R
#' @export PoisMixMean

PoisMixMean <-
function(y, K, conds, s, lambda) {
if(length(K) != 1)
	stop(paste(sQuote("K"), "(the number of clusters) must be a nonnegative integer"))
if(K < 0 | round(K) != K) 
	stop(paste(sQuote("K"), "(the number of clusters) must be a nonnegative integer"))
if(is.vector(conds) == FALSE | length(conds) != ncol(y))
	stop(paste(sQuote("conds"), "must be a vector the same length as the number of columns in", sQuote("y")))
if(length(s) != length(conds))
	stop(paste(sQuote("s"), "and", sQuote("conds"), 
	"must be vectors the same length as the number of columns in", sQuote("y")))
if(is.matrix(lambda) == FALSE | ncol(lambda) != K | nrow(lambda) != length(unique(conds)))
	stop(paste(sQuote("lambda"), "must be a (d x K) matrix"))
n <- dim(y)[1];cols <- dim(y)[2];
r <- as.vector(table(conds))
w <- rowSums(y)
mean.mat <- vector("list", K)
w.mat <- matrix(rep(w, times = cols), nrow = n, ncol = cols)
s.mat <- matrix(rep(s, each = n), nrow = n, ncol = cols)
mean.mat <- lapply(1:K, function(x)
	.myloopfxn(x, lambda=lambda, w.mat=w.mat, s.mat=s.mat, r=r, n=n, cols=cols))
return(mean.mat)
}

