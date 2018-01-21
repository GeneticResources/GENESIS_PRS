\name{numInterval}
\alias{numInterval}
\title{
Function to Calculate Number of SNPs with Effect Size in An Interval
}
\description{
This function calculates the number of SNPs with effect size falling into some interval according to the fitted mixture model.}
\usage{
numInterval(lower,upper,est,M=1070777)
}
\arguments{
  \item{\code{lower}}{the lower bound of the interval.}
  \item{\code{upper}}{the upper bound of the interval.}
  \item{\code{est}}{the parameter estimates by fitting either 2-component model, i.e., (pic, sigmasq, a); or 3-component model, i.e., (pic, p1, sigmasq1, sigmasq2, a).}
  \item{\code{M}}{total number of SNPs in the reference panel. By default, it is the total number of common SNPs in Hapmap3 reference panel, which is equal to 1070777. }
}
\value{
The number of SNPs with effect size falling into some interval (lower, upper) according to the fitted mixture model.
}
\author{
Yan Zhang <yandorazhang@gmail.com>
}
\examples{
numInterval(0.005,Inf,est=c(5.119735e-03,5.533230e-05,1.794370e-06))
numInterval(0.005,Inf,est=c(9.583307e-03,8.562964e-02,1.487684e-04,2.086576e-05,1.498790e-06)
}