#' Simulate data from a Poisson mixture model
#' 
#' This function simulates data from a Poisson mixture model, as described by
#' Rau et al. (2011). Data are simulated with varying expression level
#' (\eqn{w_i}) for 4 clusters. Clusters may be simulated with \dQuote{high} or
#' \dQuote{low} separation, and three different options are available for the
#' library size setting: \dQuote{equal}, \dQuote{A}, and \dQuote{B}, as
#' described by Rau et al. (2011).
#' 
#' 
#' @param n Number of observations
#' @param libsize The type of library size difference to be simulated
#' (\dQuote{\code{equal}}, \dQuote{\code{A}}, or \dQuote{\code{B}}, as
#' described by Rau et al. (2011))
#' @param separation Cluster separation (\dQuote{\code{high}} or
#' \dQuote{\code{low}}, as described by Rau et al. (2011))
#' @return \item{y }{(\emph{n} x \emph{q}) matrix of simulated counts for
#' \emph{n} observations and \emph{q} variables} \item{labels }{Vector of
#' length \emph{n} defining the true cluster labels of the simulated data}
#' \item{pi }{Vector of length 4 (the number of clusters) containing the true
#' value of \eqn{\ensuremath\boldsymbol{\pi}}{\pi}} \item{lambda }{(\emph{d} x
#' \emph{4}) matrix of \eqn{\ensuremath\boldsymbol{\lambda}}{\lambda} values
#' for \emph{d} conditions (3 in the case of \code{libsize =}
#' \dQuote{\code{equal}} or \dQuote{\code{A}}, and 2 otherwise) in 4 clusters
#' (see note below)} \item{w }{Row sums of \code{y} (estimate of
#' \eqn{\hat{w}})} \item{conditions }{Vector of length \emph{q} defining the
#' condition (treatment group) for each variable (column) in \code{y}}
#' @note If one or more observations are simulated such that all variables have
#' a value of 0, those rows are removed from the data matrix; as such, in some
#' cases the simulated data \code{y} may have less than \code{n} rows.
#' 
#' The PMM-I model includes the parameter constraint \eqn{\sum_k \lambda_{jk}
#' r_j = 1}, where \eqn{r_j} is the number of replicates in condition
#' (treatment group) \eqn{j}. Similarly, the parameter constraint in the PMM-II
#' model is \eqn{\sum_j \sum_l \lambda_{jk}s_{jl} = 1}, where \eqn{s_{jl}} is
#' the library size for replicate \emph{l} of condition \emph{j}. The value of
#' \code{lambda} corresponds to that used to generate the simulated data, where
#' the library sizes were set as described in Table 2 of Rau et al. (2011).
#' However, due to variability in the simulation process, the actually library
#' sizes of the data \code{y} are not exactly equal to these values; this means
#' that the value of \code{lambda} may not be directly compared to an estimated
#' value of \eqn{\hat{\ensuremath\boldsymbol{\lambda}}}{\hat{\lambda}} as
#' obtained from the \code{\link{PoisMixClus}} function.
#' @author Andrea Rau <\url{andrea.rau@@jouy.inra.fr}>
#' @references Rau, A., Celeux, G., Martin-Magniette, M.-L., Maugis-Rabusseau,
#' C. (2011). Clustering high-throughput sequencing data with Poisson mixture
#' models. Inria Research Report 7786. Available at
#' \url{http://hal.inria.fr/inria-00638082}.
#' @keywords datagen
#' @examples
#' ## Simulate data as shown in Rau et al. (2011)
#' ## Library size setting "A", high cluster separation
#' ## n = 200 observations
#' 
#' simulate <- PoisMixSim(n = 200, libsize = "A", separation = "high")
#' y <- simulate$y
#' conds <- simulate$conditions
#' 
#' @importFrom stats rexp
#' @importFrom stats runif
#' @importFrom stats rmultinom
#' @export PoisMixSim

PoisMixSim <-
function(n = 2000, libsize, separation) {

if(length(libsize) > 1)
	stop(paste(sQuote("libsize"), "must be equal to one of the following:", dQuote("A"), "or",
		dQuote("B"), "or", dQuote("equal"))) 
if(libsize != "A" & libsize != "B" & libsize != "equal")
	stop(paste(sQuote("libsize"), "must be equal to one of the following:", dQuote("A"), "or",
		dQuote("B"), "or", dQuote("equal"))) 
if(length(separation) > 1)
	stop(paste(sQuote("separation"), "must be equal to one of the following:", 
		dQuote("high"), "or", dQuote("low")))
if(separation != "high" & separation != "low")
	stop(paste(sQuote("separation"), "must be equal to one of the following:", 
		dQuote("high"), "or", dQuote("low")))
if(length(n) > 1)
	stop(paste(sQuote("n"), "must be a positive integer"))
if(n <= 0)
	stop(paste(sQuote("n"), "must be a positive integer"))
if(round(n) != n)
	stop(paste(sQuote("n"), "must be a positive integer"))

## libsize <- c("equal", "A", "B")
if(libsize == "A" | libsize == "equal") {
conds <- c(1, rep(2,4), rep(3, 3))
mean.expr <- 1640
s.norm <- c(0.156, 0.071, 0.248, 0.165, 0.014, 0.028, 0.206)
}

if(libsize == "B") {
conds <- c(rep(1, 4), rep(2,2))
mean.expr <- 1521
s.norm <- c(0.096, 0.084, 0.253, 0.205, 0.224, 0.138)
}

r <- table(conds)
cols <- length(conds)
d <- length(unique(conds))
g.true <- 4
w <- round(rexp(n, 1/mean.expr))

s.true <- ifelse(rep(libsize, cols) == "equal", rep(1/cols, cols), 
s.norm)
s.dot.true <- rep(NA, d)
for(j in 1:d) {
s.dot.true[j] <- sum(s.true[which(conds == (unique(conds))[j])])
} 
lambda.true <- matrix(NA, nrow = d, ncol = g.true)

##################################
## CHOOSING LAMBDA VALUES##
##################################

## High separation
if(separation == "high") {
tmp <- cbind(c(1,3,5), c(5,1,3), c(3,5,1), c(5,3,1))
}
## Low separation
if(separation == "low") {
tmp <- cbind(c(1,3,5), c(2,4,4), c(1,5,4), c(2,5,3))
}
if(libsize == "B") tmp <- tmp[1:d,];
## Choosing lambda values so that colSums(s.dot.true * lambda.true) = 1
for(k in 1:g.true) {
lambda.tmp <- tmp[,k]/sum(tmp[,k]);
lambda.true[,k] <- lambda.tmp/s.dot.true
}
pi.true <- c(.10, .20, .30, .40)

################################
## Simulating data##
################################

y <- matrix(NA, nrow = n, ncol = cols)
label.true <- rep(NA, n)
tmp <- runif(n); cp <- cumsum(pi.true);
for(i in 1:n) {
## Choose class label
lab <- 1
for(k in 2:g.true) {
if(tmp[i] < cp[k] & tmp[i] >= cp[k-1]) lab <- k;
}
label.true[i] <- lab
lambda.tmp <- rep(lambda.true[,label.true[i]], times = r)
y[i,] <- rmultinom(1, w[i], s.true*lambda.tmp)
}

## Remove rows with all zeros
if(min(rowSums(y) == 0)) {
y <- y[-which(rowSums(y) == 0),]
label.true <- label.true[-which(rowSums(y) == 0),]
w <- w[-which(rowSums(y) == 0),]
}

return(list(y = y, labels = label.true, pi = pi.true, lambda = lambda.true,
w = w, conditions = conds))
}

