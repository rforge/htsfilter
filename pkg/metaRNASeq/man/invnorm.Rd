\name{invnorm}
\alias{invnorm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
invnorm(pvalonesided, nrep, BHth = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pvalonesided}{
%%     ~~Describe \code{pvalonesided} here~~
}
  \item{nrep}{
%%     ~~Describe \code{nrep} here~~
}
  \item{BHth}{
%%     ~~Describe \code{BHth} here~~
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
function (pvalonesided, nrep, BHth = 0.05) 
{
    listres = vector("list", 4)
    nbstudies = length(pvalonesided)
    nbreptot = sum(nrep)
    weight = sqrt(nrep/nbreptot)
    qnormpval = do.call(cbind, lapply(pvalonesided, FUN = function(x) qnorm(1 - 
        x)))
    statc = apply(qnormpval, 1, FUN = function(x) sum(weight * 
        x, na.rm = TRUE))
    rpvalc = 1 - pnorm(statc)
    res = which(p.adjust(rpvalc, method = "BH") <= BHth)
    listres[[1]] = res
    listres[[2]] = statc
    listres[[3]] = rpvalc
    listres[[4]] = p.adjust(rpvalc, method = "BH")
    names(listres) = c("DEindices", "TestStatistic", "rawpval", 
        "adjpval")
    return(listres)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
