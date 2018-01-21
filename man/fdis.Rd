\name{fdis}
\alias{fdis}
\title{
Future Number of Discoveries Prediction Function
}
\description{
Function to predict number of discoveries in future GWAS and its confidence interval through bootstrap method. 
}
\usage{
fdis(est,v,n,alpha=0.025,M=1070777)
}
\arguments{
  \item{\code{est}}{parameter estimates by fitting either 2-component model, i.e., (pic, sigmasq, a); or 3-component model, i.e., (pic, p1, sigmasq1, sigmasq2, a).}
  \item{\code{v}}{covariance matrix of parameter estimates by fitting the 2- or 3-component model.}
  \item{\code{n}}{specifided future GWAS sample size.}
  \item{\code{alpha}}{significance level of the confidence interval. By default, it is 0.025.}
  \item{\code{M}}{total number of SNPs in the reference panel. By default, it is the total number of common SNPs in Hapmap3 reference panel, which is equal to 1070777. }
  \item{\code{nsim}}{total number of simulations. By default, it is 1000.}
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
v <- matrix(c(1.346867e-06, -1.036601e-05,  1.303674e-09, -2.444779e-09, -8.840004e-11, -8.087098e-06,  2.721732e-04,-1.428814e-07,  1.519391e-09,  7.685736e-10, -1.188552e-09, -1.270322e-07,  1.637225e-10,  8.813366e-12,-2.151989e-14, -2.584200e-09,7.023163e-09 , 3.743461e-12,  7.008698e-12,  7.664961e-14, -8.582952e-11,
8.923473e-10, -1.744992e-13,  6.291208e-14,  1.402886e-14),nrow=5,ncol=5)
fdis(est=c(9.583307e-03,8.562964e-02,1.487684e-04,2.086576e-05,1.498790e-06),v,n=253288)
}