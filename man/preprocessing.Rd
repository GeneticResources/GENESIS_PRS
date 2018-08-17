\name{preprocessing}
\alias{preprocessing}
\title{
Function to Calculate Number of SNPs with Effect Size Absolute Value Falling in an Interval
}
\description{
This function allows to preprocess the summary-level GWAS statistics. 
}
\usage{
preprocessing(summarydata, LDcutoff=0.1,LDwindow=1,filter=FALSE)
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
preprocessing(0.005,Inf,est=c(8.899809e-03, 9.476025e-02, 1.458650e-04, 2.227118e-05, 1.567643e-06))
}
