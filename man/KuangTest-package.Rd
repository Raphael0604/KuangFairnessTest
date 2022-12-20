\name{KuangTest-package}
\alias{KuangTest-package}
\alias{KuangTest}
\docType{package}
\title{
\packageTitle{KuangTest}
}
\description{
\packageDescription{KuangTest}
}
\details{

Runs tests for 0-1 series and series with a given median. 
Runs tests for 1-n series with repetition, which is designed to judge the fairness of ranks.
When using kuangtest, we enter a 0-1 series, or a regular series with a median. When using KuangFairness Test, we enter a 1-n series with repetition.

}
\author{
\packageAuthor{KuangTest}

Maintainer: \packageMaintainer{KuangTest}
}
\references{
 Friedman J H, Rafsky L C. Multivariate generalizations of the Wald-Wolfowitz and Smirnov two-sample tests[J]. The Annals of Statistics, 1979: 697-717.


Bonferroni, C. E. "Il calcolo delle assicurazioni su gruppi di teste." In Studi in Onore del Professore Salvatore Ortu Carboni. Rome: Italy, pp. 13-60, 1935.


Benjamini, Yoav, and Yosef Hochberg. "Controlling the false discovery rate: a practical and powerful approach to multiple testing." Journal of the Royal statistical society: series B (Methodological) 57.1 (1995): 289-300.

}

\keyword{ package }
\seealso{
 tseries::runs.test
}
\examples{
library(KuangTest)
a=c(1,0,0,0,1,1,0,1,0,0,1,0,0,1,1,1)
kuangtest(a) 
b=c(0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1)
kuangtest(b)
c<- c(5.36, 5.29, 5.58, 5.65, 5.57, 5.53, 5.62, 5.29, 5.44, 5.34, 5.79)
kuangtest(c,median = median(c))
d<-c(1,0,0,0,1,1,0,1,0,0,1,0,0,1,1,1,0,0,0,1,1,0,1,0,0,1,0,0,1,1,1,0,0,0,1)
kuangtest(d,bigdata = TRUE)
test<-c(1,3,3,3,1,1,3,1,3,3,1,3,3,1,1,2,2,2,2,4,4,2,2,4,4,2,2,4,4,2)
Kfairnesstest(rank =test, max = 4,method = "BH")
Kfairnesstest(rank =test, max = 4,method = "B")
}
