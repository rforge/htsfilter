\name{sim.function}
\alias{sim.function}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sim.function(param, dispFuncs, nrep = 4, classes = NULL, inter.sd = 0.3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{param}{
%%     ~~Describe \code{param} here~~
}
  \item{dispFuncs}{
%%     ~~Describe \code{dispFuncs} here~~
}
  \item{nrep}{
%%     ~~Describe \code{nrep} here~~
}
  \item{classes}{
%%     ~~Describe \code{classes} here~~
}
  \item{inter.sd}{
%%     ~~Describe \code{inter.sd} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (param, dispFuncs, nrep = 4, classes = NULL, inter.sd = 0.3) 
{
    if (is.null(classes)) {
        classes = list(c(rep(1, nrep), rep(2, nrep)), c(rep(1, 
            nrep), rep(2, nrep)))
    }
    sim.means = param[, c("mucond1", "mucond2")]
    vecstudies <- unlist(lapply(1:length(classes), FUN = function(x) rep(x, 
        length(classes[[x]]))))
    vecstudies12 <- vecstudies
    vecstudies12[which(vecstudies \%in\% seq(1, 19, by = 2))] = 1
    vecstudies12[which(vecstudies \%in\% seq(2, 20, by = 2))] = 2
    varinter <- vector("list", length(classes))
    varinter <- lapply(varinter, function(x) matrix(rnorm(dim(param)[1] * 
        2, mean = 0, sd = inter.sd), ncol = 2))
    for (i in 1:length(varinter)) {
        varinter[[i]][which(param$DE == FALSE), ] <- cbind(varinter[[i]][which(param$DE == 
            FALSE), 1], varinter[[i]][which(param$DE == FALSE), 
            1])
    }
    calc.size1 <- function(mean) 1/(dispFuncs[[1]][1] + dispFuncs[[1]][2]/mean)
    calc.size2 <- function(mean) 1/(dispFuncs[[2]][1] + dispFuncs[[2]][2]/mean)
    matsim <- mapply(function(condition, study) {
        mean <- sim.means[, condition] * exp(varinter[[study]][, 
            condition])
        size <- get(paste("calc.size", study, sep = ""))(rowMeans(sim.means))
        size <- ifelse(mean > 10, size, 1e+11)
        rnbinom(dim(param)[1], mu = mean, size = size)
    }, unlist(classes), vecstudies12)
    colnames(matsim) = paste("study", vecstudies, "cond", unlist(classes), 
        sep = "")
    matsim
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
