#' Poisson mixture model estimation and model selection
#' 
#' This function implements the EM and CEM algorithms for parameter estimation
#' in a Poisson mixture model for clustering high throughput sequencing
#' observations (e.g., genes) for a single number of clusters. Parameters are initialized using a Small-EM
#' strategy as described in Rau et al. (2011) or the splitting small-EM
#' strategy described in Papastamoulis et al. (2014), and model selection is
#' performed using the BIC/ICL criteria or the slope heuristics. 
#' 
#' Output of \code{PoisMixClus_K} is an S3 object of class \code{PoisMixClus_K}.
#' 
#' In a Poisson mixture model, the data \eqn{\mathbf{y}}{y} are assumed to come
#' from \emph{K} distinct subpopulations (clusters), each of which is modeled
#' separately; the overall population is thus a mixture of these
#' subpopulations. In the case of a Poisson mixture model with \emph{K}
#' components, the model may be written as
#' 
#' \deqn{f(\mathbf{y};K,\ensuremath\boldsymbol{\Psi}_K) = \prod_{i=1}^n
#' \sum_{k=1}^K \pi_k \prod_{j=1}^{d}\prod_{l=1}^{r_j} P(y_{ijl} ;
#' \ensuremath\boldsymbol{\theta}_k)}{f(y;K,\psi_K) = \prod_{i=1}^n
#' \sum_{k=1}^K \pi_k \prod_{j=1}^{d}\prod_{l=1}^{r_j} P(y_{ijl} ; \theta_k)}
#' 
#' for \eqn{i = 1, \ldots, n} observations in \eqn{l = 1, \ldots, r_j}
#' replicates of \eqn{j = 1, \ldots, d} conditions (treatment groups), where
#' \eqn{P(\cdot)} is the standard Poisson density,
#' \eqn{\ensuremath\boldsymbol{\Psi}_K = (\pi_1,\ldots,\pi_{K-1},
#' \ensuremath\boldsymbol{\theta}^\prime)}{\psi_K = (\pi_1,\ldots,\pi_{K-1},
#' \theta^\prime)}, \eqn{\ensuremath\boldsymbol{\theta}^\prime}{\theta^\prime}
#' contains all of the parameters in
#' \eqn{\ensuremath\boldsymbol{\theta}_1,\ldots,\ensuremath\boldsymbol{\theta}_K}{\theta_1,\ldots,\theta_K}
#' assumed to be distinct, and \eqn{\ensuremath\boldsymbol{\pi} =
#' (\pi_1,\ldots,\pi_K)^\prime}{\pi = (\pi_1,\ldots,\pi_K)^\prime} are the
#' mixing proportions such that \eqn{\pi_k} is in (0,1) for all \emph{k} and
#' \eqn{\sum_k \pi_k = 1}.
#' 
#' We consider the following parameterization for the mean
#' \eqn{\ensuremath\boldsymbol{\theta}_k = (\mu_{ijlk})}{\theta = (mu_{ijlk})}.
#' We consider \deqn{\mu_{ijlk} = w_i s_{jl} \lambda_{jk}} where \eqn{w_i}
#' corresponds to the expression level of observation \emph{i},
#' \eqn{\ensuremath\boldsymbol{\lambda}_k =
#' (\lambda_{1k},\ldots,\lambda_{dk})}{\lambda_k =
#' (\lambda_{1k},\ldots,\lambda_{dk})} corresponds to the clustering parameters
#' that define the profiles of the genes in cluster \emph{k} across all
#' variables, and \eqn{s_{jl}} is the normalized library size (a fixed
#' constant) for replicate \emph{l} of condition \emph{j}.
#' 
#' There are two approaches to estimating the parameters of a finite mixture
#' model and obtaining a clustering of the data: the estimation approach (via
#' the EM algorithm) and the clustering approach (via the CEM algorithm).
#' Parameter initialization is done using a Small-EM strategy as described in
#' Rau et al. (2011) via the \code{\link{emInit}} function. Model selection may
#' be performed using the BIC or ICL criteria, or the slope heuristics.
#' 
#' @param y (\emph{n} x \emph{q}) matrix of observed counts for \emph{n}
#' observations and \emph{q} variables
#' @param K Number of clusters (a single value). If \code{fixed.lambda}
#' contains a list of lambda values to be fixed, \code{K} corresponds to the
#' number of clusters in addition to those fixed.
#' @param conds Vector of length \emph{q} defining the condition (treatment
#' group) for each variable (column) in \code{y}
#' @param norm The type of estimator to be used to normalize for differences in
#' library size: (\dQuote{\code{TC}} for total count, \dQuote{\code{UQ}} for
#' upper quantile, \dQuote{\code{Med}} for median, \dQuote{\code{DESeq}} for
#' the normalization method in the DESeq package, and \dQuote{\code{TMM}} for
#' the TMM normalization method (Robinson and Oshlack, 2010). Can also be a
#' vector (of length \emph{q}) containing pre-estimated library size estimates
#' for each sample. Note that if the user provides pre-calculated normalization
#' factors, the package will make use of \code{norm/sum(norm)} as normalization
#' factors.
#' @param init.type Type of initialization strategy to be used
#' (\dQuote{\code{small-em}} for the Small-EM strategy described in Rau et al.
#' (2011), and \dQuote{\code{kmeans}} for a simple \emph{K}-means
#' initialization)
#' @param init.runs Number of runs to be used for the Small-EM strategy
#' described in Rau et al. (2011), with a default value of 1
#' @param init.iter Number of iterations to be used within each run for the
#' Small-EM strategry, with a default value of 10
#' @param alg.type Algorithm to be used for parameter estimation
#' (\dQuote{\code{EM}} or \dQuote{\code{CEM}})
#' @param cutoff Cutoff to declare algorithm convergence (in terms of
#' differences in log likelihoods from one iteration to the next)
#' @param iter Maximum number of iterations to be run for the chosen algorithm
#' @param fixed.lambda If one (or more) clusters with fixed values of lambda is
#' desired, a list containing vectors of length \emph{d} (the number of
#' conditions).  specifying the fixed values of lambda for each fixed cluster.
#' @param equal.proportions If \code{TRUE}, the cluster proportions are set to
#' be equal for all clusters. Default is \code{FALSE} (unequal cluster
#' proportions).
#' @param prev.labels A vector of length \emph{n} of cluster labels obtained
#' from the previous run (K-1 clusters) to be used with the splitting small-EM
#' strategy described in described in Papastamoulis et al. (2014). For other
#' initialization strategies, this parameter takes the value NA
#' @param prev.probaPost An \emph{n} x (\emph{K}-1) matrix of the conditional
#' probabilities of each observation belonging to each of the \emph{K}-1
#' clusters from the previous run, to be used with the splitting small-EM
#' strategy of described in Papastamoulis et al. (2012). For other
#' initialization strategies, this parameter takes the value NA
#' @param verbose If \code{TRUE}, include verbose output
#' @param interpretation If \code{"sum"}, cluster behavior is interpreted with
#' respect to overall gene expression level (sums per gene), otherwise for
#' \code{"mean"}, cluster behavior is interpreted with respect to mean gene
#' expression (means per gene).
#' @param EM.verbose If \code{TRUE}, more informative output is printed about
#' the EM algorithm, including the number of iterations run and the difference
#' between log-likelihoods at the last and penultimate iterations.
#' @param subset.index Optional vector providing the indices of a subset of
#' genes that should be used for the co-expression analysis (i.e., row indices
#' of the data matrix \code{y}.
#' @param wrapper \code{TRUE} if the \code{PoisMixClus_K} function is run from
#' within the \code{PoisMixClus} main function, and \code{FALSE}
#' otherwise. This mainly helps to avoid recalculating parameters several times
#' that are used throughout the algorithm (e.g., library sizes, etc.)
#' 
#' 
#' @return 
#' \item{lambda }{(\emph{d} x \emph{K}) matrix containing the estimate
#' of \eqn{\hat{\ensuremath\boldsymbol{\lambda}}}{\hat{\lambda}}} 
#' \item{pi
#' }{Vector of length \emph{g} containing the estimate of
#' \eqn{\hat{\ensuremath\boldsymbol{\pi}}}{\hat{\pi}}} 
#' \item{labels }{Vector of
#' length \emph{n} containing the cluster assignments of the \emph{n}
#' observations} 
#' \item{probaPost }{Matrix containing the conditional
#' probabilities of belonging to each cluster for all observations}
#' \item{log.like }{Value of log likelihood} 
#' \item{BIC }{Value of BIC criterion} 
#' \item{ICL }{Value of ICL criterion} 
#' \item{alg.type }{Estimation algorithm used; matches the argument \code{alg.type} above)} 
#' \item{norm}{Library size normalization factors used} 
#' \item{conds }{Conditions specified by user} 
#' \item{iterations }{Number of iterations run}
#' \item{logLikeDiff }{Difference in log-likelihood between the last and
#' penultimate iterations of the algorithm} 
#' \item{subset.index }{If provided by the user, the indices of subset of genes 
#' used for co-expression analyses}
#' \item{K }{Number of clusters provided by the user}
#' 
#' @note Note that the \code{fixed.lambda} argument is primarily intended to be
#' used in the case when a single cluster is fixed to have equal clustering
#' parameters lambda across all conditions (i.e.,
#' \eqn{\lambda_{j1}=\lambda_{1}=1}); this may be useful when
#' identifying genes with non-differential expression across all conditions. 
#' Alternatively, this
#' argument could be used to specify a cluster for which genes are only
#' expressed in a single condition (e.g., \eqn{\lambda_{11} = 1} and
#' \eqn{\lambda_{j1} = 0} for all \eqn{j > 1}). Other possibilities could be
#' considered, but note that the fixed values of lambda must satisfy the
#' constraint \eqn{\sum_j \lambda_{jk}s_{j.} = 1} for all \eqn{k} imposed in
#' the model; if this is not the case, a warning message will be printed.
#' 
#' @author Andrea Rau <\url{andrea.rau@@jouy.inra.fr}>
#' @seealso \code{\link{probaPost}} for the calculation of the conditional
#' probability of belonging to a cluster; \code{\link{PoisMixMean}} for the
#' calculation of the per-cluster conditional mean of each observation;
#' \code{\link{logLikePoisMixDiff}} for the calculation of the log likelihood
#' of a Poisson mixture model; \code{\link{emInit}} and \code{\link{kmeanInit}}
#' for the Small-EM parameter initialization strategy
#' 
#' @references 
#' Anders, S. and Huber, W. (2010) Differential expression analysis
#' for sequence count data. \emph{Genome Biology}, \bold{11}(R106), 1-28.
#' 
#' Papastamoulis, P., Martin-Magniette, M.-L., and Maugis-Rabusseau, C. (2014).
#' On the estimation of mixtures of Poisson regression models with large number
#' of components. \emph{Computational Statistics and Data Analysis}: 3rd
#' special Issue on Advances in Mixture Models, DOI:
#' 10.1016/j.csda.2014.07.005.
#' 
#' Rau, A., Maugis-Rabusseau, C., Martin-Magniette, M.-L., Celeux, G. (2015)
#' Co-expression analysis of high-throughput transcriptome sequencing data with
#' Poisson mixture models. Bioinformatics, doi: 10.1093/bioinformatics/btu845.
#' 
#' Rau, A., Celeux, G., Martin-Magniette, M.-L., Maugis-Rabusseau, C (2011).
#' Clustering high-throughput sequencing data with Poisson mixture models.
#' Inria Research Report 7786. Available at
#' \url{http://hal.inria.fr/inria-00638082}.
#' @keywords models cluster
#' 
#' @example /inst/examples/PoisMixClus.R
#' 
#' @importFrom edgeR calcNormFactors
#' @importFrom stats quantile
#' @importFrom stats median
#' @export


PoisMixClus_K <- function(y, K, conds, norm="TMM",  
	init.type = "small-em", init.runs = 1, init.iter = 10, alg.type = "EM", cutoff = 10e-6, 
	iter = 1000, fixed.lambda = NA, equal.proportions = FALSE, prev.labels=NA, 
	prev.probaPost = NA, verbose = FALSE, interpretation="sum", EM.verbose=FALSE, wrapper=FALSE,
	subset.index = NA) {

	## Next loop only run if PoisMixClus_K is called directly (otherwise called by PoisMixClus)
	if(wrapper==FALSE) {
		if(is.matrix(y) == FALSE & is.data.frame(y) == FALSE) 
			stop(paste(sQuote("y"), "must be a matrix"))
		if(min(y) < 0 | sum(as.numeric(round(y))) != sum(as.numeric(y))) 
			stop(paste(sQuote("y"), "must be a matrix made up of nonnegative counts"))
		if(min(rowSums(y)) == 0)
			stop(paste("at least one observation in", sQuote("y"), 
			"contains all 0's and must be removed from the data"))
		if(is.vector(conds) == FALSE | length(conds) != ncol(y))
			stop(paste(sQuote("conds"), "must be a vector the same length as the number of columns in", sQuote("y")))
		if(length(init.type) > 1)
			stop(paste(sQuote("init.type"), "must be of length 1"))
		if(init.type != "small-em" & init.type != "kmeans" & init.type != "split.small-em") 
			stop(paste(sQuote("init.type"), "must be one of", dQuote("small-em"), "or", 
				dQuote("kmeans"), "or", dQuote("split.small-em")))
		if(alg.type != "EM" & alg.type != "CEM")
			stop(paste(sQuote("alg.type"), "must be one of", dQuote("EM"), "or", dQuote("CEM")))
		if(length(alg.type) > 1)
			stop(paste(sQuote("alg.type"), "must be one of", dQuote("EM"), "or", dQuote("CEM")))
		if(is.logical(verbose) == FALSE)
			stop(paste(sQuote("verbose"), "must be", dQuote("TRUE"), "or", dQuote("FALSE")))
		if(class(fixed.lambda) != "list" & is.na(fixed.lambda[1]) == FALSE)
			stop(paste(sQuote("fixed.lambda"), "must be", dQuote("NA") , "or a list."))
		if(is.vector(prev.labels) == FALSE & is.na(prev.labels[1]) == FALSE)
			stop(paste(sQuote("prev.labels"), "must be", dQuote("NA") , "or a vector of labels."))

		## Grouping columns of y in order of condition (all replicates put together)
		o.ycols <- order(conds)
		y <- y[,o.ycols]
		conds <- conds[o.ycols]
		conds.names <- unique(conds)
		d <- length(unique(conds))
		r <- as.vector(table(conds))
		if(length(rownames(y)) == 0) rn <- 1:nrow(y);
		if(length(rownames(y)) > 0) rn <- rownames(y);
		y <- as.matrix(y, nrow = nrow(y), ncol = ncol(y))
		rownames(y) <- rn;
		

		## Only calculate s values if they are not provided
		if(length(norm) != 1 & length(norm) != length(conds)) stop(paste(sQuote("norm"), "must be one of
		the following: none, TC, UQ, Med, DESeq, TMM, or a vector of the same length as", sQuote("conds")))
		## If estimated from data, all genes should be used
		if(length(norm) == 1) {
			if(norm == "none") 	s <- rep(1, cols);
			if(norm == "TC") s <- colSums(y) / sum(as.numeric(y));
			if(norm == "UQ") s <- apply(y, 2, quantile, 0.75) / sum(apply(y, 2, quantile, 0.75));
			if(norm == "Med") s <- apply(y, 2, median) / sum(apply(y, 2, median));
			if(norm == "DESeq") {
				## Code from DESeq, v1.8.3
				loggeomeans <- rowMeans(log(y))
				s <- apply(y, 2, function(x) exp(median((log(x)-loggeomeans)[is.finite(loggeomeans)])))
				s <- s / sum(s)
			}
			if(norm == "TMM") {
				f <- calcNormFactors(as.matrix(y), method = "TMM")
				s <- colSums(y)*f / sum(colSums(y)*f)
			} 
		}
		if(length(norm) == length(conds)) {
			s <- norm / sum(norm)
		}
		s.dot <- rep(NA, d) 
		for(j in 1:d) {
			s.dot[j] <- sum(s[which(conds == (unique(conds))[j])])
		}
		
		## In case only a subset of data are to be used for analysis
		if(is.na(subset.index)[1] == FALSE) {
			y <- y[subset.index,]
			n <- dim(y)[1];cols <- dim(y)[2]
			w <- rowSums(y)
		}
		if(is.na(subset.index)[1] == TRUE) {
			n <- dim(y)[1];cols <- dim(y)[2]
			w <- rowSums(y)
		}
	}
	
	if(wrapper==TRUE) {
		conds.names <- unique(conds)
		d <- length(unique(conds))
		r <- as.vector(table(conds))
		n <- dim(y)[1];cols <- dim(y)[2]
		w <- rowSums(y)
		s <- norm
		s.dot <- rep(NA, d) 
		for(j in 1:d) {
			s.dot[j] <- sum(s[which(conds == (unique(conds))[j])])
		}
	}

	if(class(fixed.lambda) == "list") {
		for(ll in 1:length(fixed.lambda)) {
			if(is.vector(fixed.lambda[[ll]]) == FALSE |
				length(fixed.lambda[[ll]]) != d)
				stop(paste(sQuote("fixed.lambda"), "must be", dQuote("NA") , 
					"or a list of vectors with length equal to the number of conditions."))
			if(length(which(fixed.lambda[[ll]] == 0)) > 0) {
				if(length(which(fixed.lambda[[ll]] == 1)) + 
					length(which(fixed.lambda[[ll]] == 0)) == length(fixed.lambda[[ll]])) {

					tmp <- 1/sum(s.dot[which(fixed.lambda[[ll]] == 1)])
					new <- fixed.lambda[[ll]]
					new[which(fixed.lambda[[ll]] == 1)] <- tmp
					message(cat("Fixed lambda\n", fixed.lambda[[ll]], "\n", "modified to\n",
						new, "\n", "to satisfy imposed parameter constraints.\n"))
					fixed.lambda[[ll]][which(fixed.lambda[[ll]] == 1)] <- tmp
				}
			}
			if(abs(as.numeric(fixed.lambda[[ll]] %*% s.dot)-1) > 10e-8)
				warning(paste("Check that constraint on lambda*s.dot is upheld for",sQuote("fixed.lambda")))
		}
		K <- K + length(fixed.lambda);	
	}
	diff <- 100 ## Convergence criterion
	index <- 0; go <- 1;

	## Inital values
	## init.type: "kmeans", "small-em", "split.small-em"

	if(init.type == "kmeans") {
		init.alg <- "kmeanInit";
		init.args <- list(y = y, K = K, conds = conds, norm=s, 
			fixed.lambda = fixed.lambda, equal.proportions = equal.proportions)
	}
	if(init.type == "small-em") {
		init.alg <- "emInit"
		init.args <- list(y = y, K = K, conds = conds, norm=s,
			alg.type = alg.type, init.run = init.runs,
			init.iter = init.iter, fixed.lambda = fixed.lambda, 
			equal.proportions = equal.proportions, verbose = verbose)
	}
	if(init.type == "split.small-em") {
		init.alg <- "splitEMInit"
		init.args <- list(y = y, K = K, conds = conds, norm=s,
			alg.type = alg.type,
			fixed.lambda = fixed.lambda,
			equal.proportions = equal.proportions, 
			prev.labels = prev.labels, prev.probaPost = prev.probaPost,
			init.iter = init.iter, init.runs = init.runs, verbose = verbose)
	}
	## Adding quote = TRUE to speed up do.call
	param.init <- do.call(init.alg, init.args, quote=TRUE)

	if(equal.proportions == FALSE) {
		pi <- pi.old <- param.init$pi.init
	}
	if(equal.proportions == TRUE) {
		pi <- pi.old <- rep(1/K, K)
	}
	lambda <- lambda.old <- param.init$lambda.init
	mean.calc <- mean.old <- PoisMixMean(y = y, K = K, conds = conds, 
		s = s, lambda = lambda)

	while(go == 1) {

		############
		## E-step ##
		############
		t <- probaPost(y, K, conds, pi, s, lambda)

		############
		## C-step ##
		############
		if(alg.type == "CEM") {
			## If two values of t_{ik} are map, 
			## arbitrarily choose the first
			partition <- unlist(apply(t, 1, 
				function(x) which(x == max(x, na.rm = TRUE))[1]))
			partition.mat <- matrix(0, nrow = n, ncol = K)
			for(i in 1:n) partition.mat[i,partition[i]] <- 1;
		}

		############
		## M-step ##
		############
		if(alg.type == "CEM") {
			if(equal.proportions == FALSE) {
				for(k in 1:K) {
					pi[k] <- length(which(partition == k))/n
				}
			}
			if(equal.proportions == TRUE) {
				pi <- rep(1/K, K)
			}
			denom <- colSums(partition.mat * w)
			if(class(fixed.lambda) != "list") {
				for(j in 1:d) {
					denom.bis <- denom * s.dot[j]
					num <- colSums(partition.mat * 
					rowSums(as.matrix(y[,which(conds == (unique(conds))[j])])))
					lambda[j,] <- num / denom.bis
				}
			}
			if(class(fixed.lambda) == "list") {
				for(ll in 1:length(fixed.lambda)) {
					lambda[,ll] <- fixed.lambda[[ll]]
				}
				for(j in 1:d) {
					denom.bis <- denom * s.dot[j]
					denom.bis <- denom.bis[-c(1:length(fixed.lambda))]
					num <- colSums(partition.mat * 
						rowSums(as.matrix(y[,which(conds == (unique(conds))[j])])))
					num <- num[-c(1:length(fixed.lambda))]
					lambda[j,-c(1:length(fixed.lambda))] <- num / denom.bis
				}
			}
		}

		if(alg.type == "EM") {
			if(equal.proportions == FALSE) {
				pi <- colSums(t)/n
			}
			if(equal.proportions == TRUE) {
				pi <- rep(1/K, K)
			}
			denom <- colSums(t * w)
			if(class(fixed.lambda) != "list") {
				denom.bis <- matrix(rep(denom, length(s.dot)) * rep(s.dot, each=K), byrow=T, ncol=K)
				num <- matrix(rowsum(matrix(y, ncol=n, nrow=length(conds), byrow=T), group=conds), 
					nrow=n, ncol=d, byrow=T)
				num <- matrix(unlist(lapply(lapply(1:d, function(x) t*num[,x]), colSums), recursive=FALSE,
					use.names=FALSE), nrow=d, ncol=K, byrow=T)
				lambda <- num/denom.bis
			}
			if(class(fixed.lambda) == "list") {
				for(ll in 1:length(fixed.lambda)) {
					lambda[,ll] <- fixed.lambda[[ll]]
				}
				## This loop could be improved for speed as above
				for(j in 1:d) {
					denom.bis <- denom * s.dot[j]
					denom.bis <- denom.bis[-c(1:length(fixed.lambda))]
					num <- colSums(t * 
						matrix(rep(rowSums(as.matrix(y[,which(conds == (unique(conds))[j])])),K), 
						ncol = K))
					num <- num[-c(1:length(fixed.lambda))]
					lambda[j,-c(1:length(fixed.lambda))] <- num / denom.bis
				}
			}	
		}

		#################
		## Convergence ##
		#################
		mean.calc <- PoisMixMean(y, K = K, conds, s, lambda)
		diff <- abs(logLikePoisMixDiff(y, mean.calc, pi, mean.old, pi.old))
		lambda.old <- lambda; pi.old <- pi; mean.old <- mean.calc;

		index <- index + 1
		if(verbose == TRUE) print(paste("Log-like diff:", diff))
		if(diff < cutoff) go <- 0;
		if(iter != FALSE & iter == index) go <- 0;
	}
	
	if(EM.verbose == TRUE) {
		cat("#####################################\n")
		cat("Number of EM iterations:", index, "\n")
		cat("Last log-likelihood difference:", diff, "\n")
		cat("#####################################\n")
	}

	#####################################
	## Final estimates of lambda and p ##
	#####################################
	names(pi) <- paste("Cluster", 1:K)
	colnames(lambda) <- paste("Cluster", 1:K)
	rownames(lambda) <- conds.names
	lambda.final <- lambda
	pi.final <- pi

	## Check to make sure one of the components is not degenerate
	if(min(pi) == 0 | is.nan(sum(lambda)) == TRUE) {
		probaPost <- NA
		labels <- NA
		BIC <- NA
		ICL <- NA
	}

	if(min(pi) > 0 | is.nan(sum(lambda)) == FALSE) {

		mean.calc <- PoisMixMean(y, K = K, conds, s, lambda)
		LL.tmp <- logLikePoisMix(y, mean.calc, pi)
		LL <- LL.tmp$ll
	
		######################
		## Determine labels ##
		######################
		t <- probaPost(y, K, conds, pi, s, lambda)
		## If two clusters have exactly identical map estimators,
		## arbitrarily choose the first one
		map <- unlist(apply(t, 1, function(x) which(x == max(x, 
			na.rm = TRUE))[1]))
		z <- matrix(0, nrow = n, ncol = K)
		for(i in 1:n) z[i,map[i]] <- 1;
		probaPost <- t
		labels <- map

		##############################
		## Calculate BIC, ICL       ##
		##############################
		if(equal.proportions == FALSE & class(fixed.lambda) != "list") {
#			np <- (K-1) + n + (d-1)*K 	# pi + w + lambda
			## CHANGE September 25, 2013: w is not considered as a parameter
			np <- (K-1) + (d-1)*K 	# pi + lambda
		}
		if(equal.proportions == TRUE & class(fixed.lambda) != "list") {
#			np <- n + (d-1)*K 	# w + lambda
			## CHANGE September 25, 2013: w is not considered as a parameter
			np <- (d-1)*K 		# lambda
		}
		if(equal.proportions == FALSE & class(fixed.lambda) == "list") {
#			np <- (K-1) + n + (d-1)*(K-length(fixed.lambda))	# pi + w + lambda not fixed
			## CHANGE September 25, 2013: w is not considered as a parameter
			np <- (K-1) + (d-1)*(K-length(fixed.lambda))		# pi + lambda not fixed
		}
		if(equal.proportions == TRUE & class(fixed.lambda) == "list") {
#			np <- n + (d-1)*(K-length(fixed.lambda))			# w + lambda not fixed
			## CHANGE September 25, 2013: w is not considered as a parameter
			np <- (d-1)*(K-length(fixed.lambda))				# lambda not fixed

		}
		BIC <- -LL + (np/2) * log(n)
#		entropy <- -2*sum(z*log(t), na.rm = TRUE)
		## CHANGE October 18, 2013: replace z with t in the entropy calculation for ICL
#		entropy <- -2*sum(t*log(t), na.rm = TRUE)
		## CHANGE July 25, 2014: typo in entropy (thanks, Melina Gallopin!)
		entropy <- -sum(t*log(t), na.rm = TRUE)
		ICL <- BIC + entropy
	}
	## Should cluster behavior (lambda) be interpretated wrt the gene means or sums?	
	if(interpretation == "mean") {
		s <- s * q
		w <- w / q
	}
	results <- list(lambda = lambda.final, pi = pi.final, labels = labels, 
		probaPost = probaPost, log.like = LL, BIC = -BIC, ICL = -ICL, 
		alg.type = alg.type, norm = s,
		conds = conds, iterations = index, logLikeDiff = diff, model.selection = NA,
		subset.index = subset.index, K=K)

	class(results) <- "PoisMixClus_K"
	return(results)
}

