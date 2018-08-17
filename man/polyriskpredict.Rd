\name{polyriskpredict}
\alias{polyriskpredict}
\title{
Polygenic Risk Prediction Function
}
\description{
This function allows to get parameter estimates from fitting the mixture model.
}
\usage{
postmarginal(est,c0=10,summarydata,filter=FALSE)
}
\arguments{
  \item{\code{est}}{parameter estimates by fitting either 2-component model, i.e., (pic, sigmasq, a); or 3-component model, i.e., (pic, p1, sigmasq1, sigmasq2, a).}
    \item{\code{c0}}{an assumed maximum number of underlying susceptibility SNPs tagged by any individual GWAS marker. By default, c0 is set at 10.}
   \item{\code{summarydata}}{summay-level GWAS data, containing 3 columns: 
  SNP (SNP rsID), 
  Z (GWAS test z-statistic), 
  N (GWAS study sample size which can be different for different SNPs)}
  \item{\code{filter}}{logical; if TRUE, the input summary data will be filtered.}
}
\details{}
\value{
A list of future postmarginals with components
\item{\code{p_beta0}}{posterior probability for the marker to have a zero marginal effect.}
\item{\code{beta}}{posterior mean (best prediction) of marginal effect size.}
}
\references{
Zhang, Yan, et al. "Estimation of complex effect-size distributions using summary-level statistics from genome-wide association studies across 32 complex traits" Nature Genetics (2018).
}
\author{
Yan Zhang <yandorazhang@gmail.com>
}
\examples{
data(heightGWAS)
postmarginal(est=c(9.583307e-03,8.562964e-02,1.487684e-04,2.086576e-05,1.498790e-06),c0=10,df,filter=FALSE)
}