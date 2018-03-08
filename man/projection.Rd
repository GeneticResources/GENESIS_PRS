\name{projection}
\alias{projection}
\title{
Future Projection Function
}
\description{
Function to make future projections about the number of expected discoveries and associated heritability explained in future studies based on analytic formulas according to the fitted model. 
}
\usage{
projection(est,v,n,gwsignificance=5e-8,tol=c(1e-12,1e-15),M=1070777,nsim=10000,alpha=0.025)
}
\arguments{
  \item{\code{est}}{parameter estimates by fitting either 2-component model, i.e., (pic, sigmasq, a); or 3-component model, i.e., (pic, p1, sigmasq1, sigmasq2, a).}
  \item{\code{v}}{covariance matrix of parameter estimates by fitting the 2- or 3-component model.}
  \item{\code{n}}{future GWAS sample size.}
  \item{\code{gwsignificance}}{genome-wide significance level, by default it is 5e-8.}
  \item{\code{tol}}{tolerance accuracy vector for intgrate() function. By default, it is c(1e-12,1e-15).}
  \item{\code{M}}{total number of SNPs in the reference panel. By default, it is the total number of common SNPs in Hapmap3 reference panel, which is equal to 1070777. }
  \item{\code{nsim}}{total number of simulations; by default, it is 10000.}
  \item{\code{alpha}}{significance level of the confidence interval; by default, it is 0.025, i.e., 95\% CI.}
}
\details{}
\value{
A list of future projections with components
\item{\code{Numdicoveries}}{predictions of future number of significant discoveries as well as its CI.}
\item{\code{GVpercentage}}{predictions of percentage of genetic variance explained in future GWAS as well as its CI.}
}
\references{
Zhang, Yan, et al. "Estimation of complex effect-size distributions using summary-level statistics from genome-wide association studies across 32 complex traits and implications for the future." bioRxiv (2017): 175406.
}
\author{
Yan Zhang <yandorazhang@gmail.com>
}
\examples{
projection(est=c(5.119735e-03,5.533230e-05,1.794370e-06),v, M=1070777, n=253288);
projection(est=c(9.583307e-03,8.562964e-02,1.487684e-04,2.086576e-05,1.498790e-06),v, n=253288);
}
