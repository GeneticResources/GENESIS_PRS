\name{postmarginal}
\alias{postmarginal}
\title{
Posterior inference of marginal effect size
}
\description{
This function allows to estimate posterior probability with a null marginal effect size and posterior mean of marginal effect size. 
}
\usage{
postmarginal(summarydata, est, c0=10, LDcutoff=0.1, LDwindow=1, filter=FALSE)
}
\arguments{
 \item{\code{summarydata}}{summay-level GWAS data, containing 3 columns: 
  SNP (SNP rsID), 
  Z (GWAS test z-statistic), 
  N (GWAS study sample size which can be different for different SNPs)}
  \item{\code{est}}{parameter estimates by fitting either 2-component model, i.e., (pic, sigmasq, a); or 3-component model, i.e., (pic, p1, sigmasq1, sigmasq2, a).}
    \item{\code{c0}}{an assumed maximum number of underlying susceptibility SNPs tagged by any individual GWAS marker. By default, c0 is set at 10.}
  \item{\code{LDcutoff}}{a number from (0.05, 0.1, 0.2); indicating LD score is calculated based on the particular cutoff. By default, it is 0.1.}
   \item{\code{LDwindow}}{a number from (0.5, 1, 2); indicating LD score is calculated based on the particular window size (MB). By default, it is 1 MB.}
  \item{\code{filter}}{logical; if TRUE, the input summary data will be filtered.}
}
\details{}
\value{
A list of future postmarginals with components
\item{\code{prob.beta0}}{posterior probability for the marker to have a zero marginal effect.}
\item{\code{posteriormean.marginal}}{posterior mean (best prediction) of marginal effect size.}
}
\references{
Zhang, Yan, et al. "Estimation of complex effect-size distributions using summary-level statistics from genome-wide association studies across 32 complex traits" Nature Genetics (2018).
}
\author{
Yan Zhang <yandorazhang@gmail.com>
}
\examples{
data(heightGWAS)
postmarginal(heightGWAS, est=c(9.583307e-03,8.562964e-02,1.487684e-04,2.086576e-05,1.498790e-06), c0=10, LDcutoff=0.1, LDwindow=1, filter=FALSE)
}
