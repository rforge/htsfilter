\name{fishercomb}
\alias{fishercomb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fishercomb(indpval, BHth = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{indpval}{
%%     ~~Describe \code{indpval} here~~
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
function (indpval, BHth = 0.05) 
{
    listres = vector("list", 4)
    logpval = do.call(cbind, lapply(indpval, log))
    statc = apply(logpval, 1, FUN = function(x) -2 * sum(x, na.rm = TRUE))
    statc = apply(logpval, 1, FUN = function(x) -2 * sum(x))
    rpvalc = 1 - pchisq(statc, df = (2 * length(indpval)))
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
