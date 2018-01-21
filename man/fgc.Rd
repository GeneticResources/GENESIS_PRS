\name{fgc}
\alias{fgc}
\title{
Future GC Factor Prediction Function
}
\description{
Function to predict future genomic control (GC) factor through simulations.
}
\usage{
fgc(est,n,nsim=100,M=1070777)
}
\arguments{
  \item{\code{est}}{parameter estimates by fitting either 2-component model, i.e., (pic, sigmasq, a); or 3-component model, i.e., (pic, p1, sigmasq1, sigmasq2, a).}
  \item{\code{n}}{specifided future GWAS sample size.}
  \item{\code{nsim}}{total number of simulations. By default, it is 100.}
  \item{\code{M}}{total number of SNPs in the reference panel. By default, it is the total number of common SNPs in Hapmap3 reference panel, which is equal to 1070777. }
}
\details{
This function allows to predict future genomic control (GC) factor through simulations. It will first generate joint effect size according to the fitted model; then get the marginal effect size by the LD structure across SNPs; then errors following a multivariate normal with mean zero and variance the LD correlation matrix are added to get the estimates of the marginal effect sizes. The predicted GC factor is then calculated. 
}
\value{
Mean predicted genomic factor value. 
}
\references{
Zhang, Yan, et al. "Estimation of complex effect-size distributions using summary-level statistics from genome-wide association studies across 32 complex traits and implications for the future." bioRxiv (2017): 175406.
}
\author{
Yan Zhang <yandorazhang@gmail.com>
}
\examples{
fgc(est=c(5.119735e-03,5.533230e-05,1.794370e-06), M=1070777,n=253288,nsim=1);
fgc(est=c(9.583307e-03,8.562964e-02,1.487684e-04,2.086576e-05,1.498790e-06),n=253288,nsim=1);
}