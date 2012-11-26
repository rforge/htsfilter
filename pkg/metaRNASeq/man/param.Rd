\name{param}
\alias{param}
\docType{data}
\title{
Mean simulation parameters
}
\description{
Mean simulation parameters obtained from the analysis of a real dataset
}
\usage{data(param)}
\format{
  A data frame with 26408 observations on the following 3 variables.
  \describe{
    \item{\code{mucond1}}{a numeric vector with mean parameters for condition 1}
    \item{\code{mucond2}}{a numeric vector with mean parameters for condition 2}
    \item{\code{DE}}{a logical vector indicating which tags are differentially expressed (value 1)}
  }
}
\details{
Mean parameters provided in this package are empirical means obtained after normalization for library size differences of real data described in the following references.
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(param)
## maybe str(param) ; plot(param) ...
}
\keyword{datasets}
