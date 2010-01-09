\name{envelope-package}
\Rdversion{1.1}
\alias{envelope-package}
\alias{envelope}
\docType{package}
\title{
 Envelope normality test.
}
\description{
 Constructs total coverage confidence bands with a fixed level for a normality plot.
}
\details{
 \tabular{ll}{
  Package: \tab envelope\cr
  Type: \tab Package\cr
  Version: \tab 1.0\cr
  Date: \tab 2010-01-09\cr
  License: \tab GPL version 2 or newer\cr
  LazyLoad: \tab yes\cr
 }
 Available functions: \cr
 envl.plot: creates a normality plot with confidence bands of a dataset with the specified level of confidence.  \cr
 qqenvl: adds confidence bands to a plot created with \code{\link{qqnorm}}.  \cr
 envl.pval: approximates the p-value used to construct the envelope for the corresponding sample size and confidence level.  \cr
}
\author{
 Felipe Acosta
 Maintainer: Felipe Acosta <felip_aa@math.uprm.edu>
}
\references{
 W. Rolke (2003): An Extension of the Normal Probability Plot. Unpublished manuscript. \url{http://charma.uprm.edu/~rolke/publications.htm}
}
\keyword{ package}
\keyword{ htest }
\seealso{
 \code{\link{qqnorm}}
}
\examples{
 #95\% Confidence bands for a random sample of size 30.
 x=rnorm(30)
 envl.plot(x)

 #Add 90\% confidence bands to a plot generated with qqnorm()
 x=rnorm(20,10,5)
 qqnorm(x)
 qqenvl(x,conf=90)
}
