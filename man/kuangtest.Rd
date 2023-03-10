% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{kuangtest}
\alias{kuangtest}
\title{Runs tests for 0-1 series and series with a given median.}
\usage{
kuangtest(run, median = 0.200164, bigdata = FALSE)
}
\arguments{
\item{run}{0-1 series. If median is used, it can be any series to test its randomness.}

\item{median}{The defaulted value is 0.200164, it means nothing.}

\item{bigdata}{if the repetition of any element is over 16, it's better to set it as true. Here, central limit theorem will be applied to replace hyper-geometric distribution.}
}
\value{
The p-value of a one-sided test.
}
\description{
Runs tests for 0-1 series and series with a given median.
}
\examples{
library(KuangTest)
a=c(1,0,0,0,1,1,0,1,0,0,1,0,0,1,1,1)
kuangtest(a)   
b=c(0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1)
kuangtest(b)
c<- c(5.36, 5.29, 5.58, 5.65, 5.57, 5.53, 5.62, 5.29, 5.44, 5.34, 5.79, 5.10, 5.27, 5.39, 5.42, 5.47, 5.63, 5.34, 5.46, 5.30, 5.75, 5.68, 5.85)
kuangtest(c,median = median(c))
d<-c(1,0,0,0,1,1,0,1,0,0,1,0,0,1,1,1,0,0,0,1,1,0,1,0,0,1,0,0,1,1,1,0,0,0,1,1,0,1)
kuangtest(d,bigdata = TRUE)
}
\references{
Friedman J H, Rafsky L C. Multivariate generalizations of the Wald-Wolfowitz and Smirnov two-sample tests[J]. The Annals of Statistics, 1979: 697-717.
}
