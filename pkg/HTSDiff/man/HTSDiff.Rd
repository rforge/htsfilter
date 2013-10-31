\name{HTSDiff}
\alias{HTSDiff}
\title{
Differential analysis of RNA-seq data using a Poisson mixture model
}
\description{
This function implements a differential analysis of RNA-seq data using a Poisson mixture model,
where one cluster is fixed to represent genes with equal mean in each experimental condition (i.e.,
a cluster of non-differentially expressed genes).
}
\usage{
HTSDiff(counts, conds, DEclusters=4, norm="TMM", epsilon=0.8, ...)
}

\arguments{
  \item{counts}{
(\emph{n} x \emph{q}) matrix of observed counts for \emph{n} genes and \emph{q} samples, with row names corresponding
to gene IDs
}
  \item{conds}{
Vector of length \emph{q} defining the condition (treatment group) for each variable (column) in \code{counts}
}
  \item{DEclusters}{
Number of clusters to include to represent differentially expressed genes (default value of 4), in addition
to the cluster fixed to represent non-differentially expressed genes.
}
  \item{norm}{
The estimator to be used for the library size parameter: \dQuote{\code{TC}} for total count, \dQuote{\code{UQ}} 
for upper quantile, \dQuote{\code{Med}} for median, \dQuote{\code{DESeq}} for the normalization method in the 
DESeq package, and \dQuote{\code{TMM}} for the TMM normalization method (default).
}
  \item{epsilon}{
Cutoff used to identify whether the log2-ratio of cluster parameters between conditions is sufficiently large to
be declared as differentially expressed, with default value 0.8
}
  \item{...}{
Additional parameters to be passed to the HTSCluster package, if desired. These include notably the following: 1) \code{init.runs}, the number of 
small-EM algorithms to run in initialization of Poisson mixture model estimation, with default value of 10, and 2) \code{init.iter}, the 
number of iterations to run within each small-EM algorithm in initialization of Poisson mixture model estimation, with
default value of 10
}
}
\details{

In a Poisson mixture model, the data \eqn{\mathbf{y}}{y} are assumed to come from \emph{g} distinct subpopulations (clusters), each of which is modeled separately; the overall population is thus a mixture of these subpopulations. In the case of a Poisson mixture model with \emph{g} components, the model may be written as

\deqn{f(\mathbf{y};g,\ensuremath\boldsymbol{\Psi}_g) = \prod_{i=1}^n \sum_{k=1}^g \pi_k \prod_{j=1}^{d}\prod_{l=1}^{r_j} P(y_{ijl} ; \ensuremath\boldsymbol{\theta}_k)}{f(y;g,\psi_g) = \prod_{i=1}^n \sum_{k=1}^g \pi_k \prod_{j=1}^{d}\prod_{l=1}^{r_j} P(y_{ijl} ; \theta_k)}

for \eqn{i = 1, \ldots, n} observations in \eqn{l = 1, \ldots, r_j} replicates of \eqn{j = 1, \ldots, d} conditions (treatment groups), where \eqn{P(\cdot)} is the standard Poisson density, \eqn{\ensuremath\boldsymbol{\Psi}_g = (\pi_1,\ldots,\pi_{g-1}, \ensuremath\boldsymbol{\theta}^\prime)}{\psi_g = (\pi_1,\ldots,\pi_{g-1}, \theta^\prime)}, \eqn{\ensuremath\boldsymbol{\theta}^\prime}{\theta^\prime} contains all of the parameters in \eqn{\ensuremath\boldsymbol{\theta}_1,\ldots,\ensuremath\boldsymbol{\theta}_g}{\theta_1,\ldots,\theta_g} assumed to be distinct, and \eqn{\ensuremath\boldsymbol{\pi} = (\pi_1,\ldots,\pi_g)^\prime}{\pi = (\pi_1,\ldots,\pi_g)^\prime} are the mixing proportions such that \eqn{\pi_k} is in (0,1) for all \emph{k} and \eqn{\sum_k \pi_k = 1}. We consider
\deqn{\mu_{ijlk} = w_i s_{jl} \lambda_{jk}}
where \eqn{w_i} and \eqn{\ensuremath\boldsymbol{\lambda}_k}{\lambda_k} are as before and \eqn{s_{jl}} is the normalized library size (a fixed constant) for replicate \emph{l} of condition \emph{j}. See Rau et al. (2011) for more details on this model, including parameter estimation, algorithm initialization, and model selection.

In the case of differential analysis, we fix one of the clusters (typically the first, although this choice is arbitrary) to represent non-differentially expressed genes, i.e., \eqn{\lambda_{11} = ... = \lambda_{1d} = 1}. Typically we fix the number of remaining clusters (\code{DEclusters}) to be 4, although this choice may be modified by the user. In addition to
the fixed cluster, clusters for which the absolute value of \eqn{\log_2(\lambda_{1k} / \lambda{2k})} is less than \code{epsilon} (default value 0.8) are also considered to represent
non-differentially expressed genes.

Following clusering, a gene is declared differentially expressed if its conditional probability to be non-differentially expressed (i.e., to belong to a cluster of non-differentially
expressed genes) is less than 1e-8.

}
\value{
\item{res }{Results data frame containing the following information: \code{id} = gene IDs, \code{baseMean} = base mean (for normalized 
counts), \code{baseMeanA} = base mean for condition A (for normalized counts), \code{baseMeanB} = base mean for 
condition B (for normalized counts), \code{foldChange} = fold change between \code{baseMeanA} and \code{baseMeanB},
\code{log2FoldChange}, \code{tauDE} = conditional probability of differential expression, \code{tauNDE} = conditional probability
of non-differential expression, \code{DE} = TRUE if gene is declared differentially expressed using cutoff for conditional
probability of the non-differential cluster and FALSE otherwise}
\item{PMM}{Object of class \code{HTSCluster} containing parameter estimates and other information from the Poisson mixture
model estimation}
}

\references{
S. Balzergue, G. Rigaill, V. Brunaud, E. Blondet, A. Rau, O. Rogier, J. Caius, C. Maugis-Rabusseau,
L. Soubigou-Taconnat, S. Aubourg, C. Lurin, E. Delannoy, and M.-L. Martin-Magniette. (2013) 
HTSDiff: a sensitive analysis of differential gene expression by RNA-seq (submitted).
}
\author{
Andrea Rau <\url{andrea.rau@jouy.inra.fr}>
}

\examples{

set.seed(12345)

## Generate synthetic data: 2000 genes under H0
test <- syntheticData(H0number = 2000)

## Mixture model differential analysis
## DEtest <- HTSDiff(test, c(1,1,2,2))

}
\keyword{ models }

