\name{tsb}
\alias{tsb}
\title{Rao & Wu bootstrap}
\description{
  Implements Rao & Wu's (1988) bootstrap algorithm for two-stage
  cluster samples
}
\usage{
tsb(y, N, M.0, psu, psu.size, g, B, conf.level = 0.95)
}
\arguments{
  \item{y}{a dataframe whose columns are arguments used by g}
  \item{N}{the number of PSUs in the population}
  \item{M.0}{the number of SSUs in the population}
  \item{psu}{a vector giving the PSU index for each row in y}
  \item{psu.size}{a vector giving the number of SSUs contained in this
  PSU in the population}
  \item{g}{a function that computes a numeric value based on a row of y}
  \item{B}{the number of times to iterate the bootstrap algorithm}
  \item{conf.level}{the desired confidence level for a bootstrap
  confidence interval}
}
\value{
    \item{sample}{the set of \code{g} values computed during
      resampling}
    \item{mean}{the mean of boot sample}
    \item{SE}{the sample standard deviation of boot sample}
    \item{CI}{\code{conf.level/2} and \code{1-conf.level/2}
      quantiles of boot sample}
}
\references{
  Resampling Inference With Complex Survey Data
  Author(s): J. N. K. Rao and C. F. J. Wu
  Source: Journal of the American Statistical Association, Vol. 83,
  No. 401 (Mar., 1988), pp. 231-
  241
  Published by: American Statistical Association
  Stable URL: http://www.jstor.org/stable/2288945 .

}

\examples{
library(two.stage.boot)
data(tiny.eg)
tsb.out <- tsb(y=tiny.eg$y,N=tiny.eg$N,M.0=tiny.eg$M.0,
               psu=tiny.eg$psu,psu.size=tiny.eg$psu.size,
               g=tiny.eg$g,B=1000)
tsb.out$CI
}
\keyword{survey}
