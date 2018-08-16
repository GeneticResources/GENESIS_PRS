\name{futuregc}
\alias{futuregc}
\title{
Function to Predict Future GC Factor
}
\description{
Function to predict future genomic control (GC) factor through simulations.
}
\usage{
futuregc(est,n,nsim=1,M=1070777)
}
\arguments{
  \item{\code{est}}{parameter estimates by fitting either 2-component model, i.e., (pic, sigmasq, a); or 3-component model, i.e., (pic, p1, sigmasq1, sigmasq2, a).}
  \item{\code{n}}{specifided future GWAS sample size.}
  \item{\code{nsim}}{total number of simulations. By default, it is 1.}
  \item{\code{M}}{total number of SNPs in the reference panel. By default, it is the total number of common SNPs in Hapmap3 reference panel, which is equal to 1070777. }
}
\details{
This function allows to predict future genomic control (GC) factor through simulations. It will first generate joint effect size according to the fitted model; then get the marginal effect size by the LD structure across SNPs; then errors following a multivariate normal with mean zero and variance of the LD correlation matrix are added to get the estimates of the marginal effect sizes. The predicted GC factor is then calculated based on the simulated estimates of marginal effect sizes.
}
\value{
Mean predicted genomic factor value. 
}
\references{
Zhang, Yan, et al. "Estimation of complex effect-size distributions using summary-level statistics from genome-wide association studies across 32 complex traits" Nature Genetics (2018).
}
\author{
Yan Zhang <yandorazhang@gmail.com>
}
\examples{
futuregc(est=c(8.899809e-03, 9.476025e-02, 1.458650e-04, 2.227118e-05, 1.567643e-06),n=253288,nsim=1);
}
