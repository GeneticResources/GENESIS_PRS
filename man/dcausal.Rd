\name{dcausal}
\alias{dcausal}
\title{
Mixture Normal Probability Density for the Causal SNPs.
}
\description{
This function allows to get the (mixture) normal probability density function for the causal SNPs.
}
\usage{
dcausal(x,est)
}
\arguments{
  \item{\code{x}}{number of quantiles}
  \item{\code{est}}{parameter estimates by fitting either 2-component model, i.e., (pic, sigmasq, a); or 3-component model, i.e., (pic, p1, sigmasq1, sigmasq2, a).}
}
\details{
This function allows to get the (mixture) normal probability density function for the causal SNPs. It is a single normal distribution for 2-component model; two-component mixture normal distribution for 3-component model. 
}
\value{
The density.
}
\references{
Zhang, Yan, et al. "Estimation of complex effect-size distributions using summary-level statistics from genome-wide association studies across 32 complex traits and implications for the future." bioRxiv (2017): 175406.
}
\author{
Yan Zhang <yandorazhang@gmail.com>
}
\examples{
dcausal(x,est=c(9.583307e-03,8.562964e-02,1.487684e-04,2.086576e-05,1.498790e-06))
}
