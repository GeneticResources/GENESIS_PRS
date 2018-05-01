\name{genesis}
\alias{genesis}
\title{
GENetic Effect-Size distribution Inference from Summary-level data
}
\description{
Estimation of complex effect-size distributions using summary-level statistics from genome-wide association studies. 
}
\usage{
genesis(summarydata, filter=FALSE, modelcomponents=2, cores=10, LDcutoff=0.1,LDwindow=1,M=1070777,c0=10, BICgamma=0.5, print=TRUE, printfreq=10, starting=NA, startingpic=NA, tolerance=NA, qqplot=TRUE, qqplotCI=0.8, qqplotname="", nsim=100, summaryGWASdataLDsave=FALSE, qqplotdatasave=FALSE, qqaxis=10, siblingrisk=FALSE,herit_liability=FALSE,sample_prevalence=NA,population_prevalence=NA,stratification=TRUE)
}
\arguments{
  \item{\code{summarydata}}{summay-level GWAS data, containing 3 columns: 
  SNP (SNP rsID), 
  Z (GWAS test z-statistic), 
  N (GWAS study sample size which can be different for different SNPs)}
  \item{\code{filter}}{logical; if TRUE, the input summary data will be filtered.}
  \item{\code{modelcomponents}}{2 or 3, indicating fitting 2-component or 3-component model.}
  \item{\code{cores}}{number of CPU threads in parallel computing.}
  \item{\code{LDcutoff}}{a number from (0.05, 0.1, 0.2); indicating LD score is calculated based on the particular cutoff. By default, it is 0.1.}
   \item{\code{LDwindow}}{a number from (0.5, 1, 2); indicating LD score is calculated based on the particular window size (MB). By default, it is 1 MB.}
  \item{\code{M}}{the total number of SNPs in the reference panel used when calculating LD score; by default, it is 1070777.}
  \item{\code{c0}}{an assumed maximum number of underlying susceptibility SNPs tagged by any individual GWAS marker. By default, c0 is set at 10.}
  \item{\code{BICgamma}}{a tuning parameter in calculating BIC in (0,1). By default, BICgamma is set at 0.5.}
  \item{\code{print}}{logical; if TRUE, the EM algorithm details will be output.}
  \item{\code{printfreq}}{number indicating every printfreq steps, the EM results will be output.}
  \item{\code{starting}}{the starting values for the model. For 2-component model, the starting value of (pic, sigmasq, a); for 3-component model, (pic, p1, sigmasq1, sigmasq2, a).}
  \item{\code{startingpic}}{the starting value for pic when staring==NA.}
  \item{\code{tolerance}}{the accuracy of the tolerance. For 2-component model, it is a 6-dim vector with tolerance value for (pic,sigmasq,a,llk,maxEM,steps). For 3-component model, it is a 8-dim vector with tolerance value for (pic,p1,sigmasq1,sigmasq2,a,llk,maxEM,steps).}
  \item{\code{qqplot}}{logical; if TRUE, the QQ plot will be ploted.}
  \item{\code{qqplotCI}}{the threshold of confidence interval in the QQ plot.}
  \item{\code{qqplotname}}{the name of the QQ plot files.}
  \item{\code{nsim}}{total number of simulations based on which the QQ plot is obtained.}
  \item{\code{summaryGWASdataLDsave}}{logical; if TRUE, the summary GWAS data as well as the LD information will be saved.}
  \item{\code{qqplotdatasave}}{logical; if TRUE, the simulated data to generate the QQ plot will be saved.}
  \item{\code{siblingrisk}}{logical; if TRUE, the sibling risk will be calculated.}
  \item{\code{herit_liability}}{logical; if TRUE, the heritability in log-odds-ratio scale will be transferred to liability scale.}
  \item{\code{sample_prevelance}}{sample prevelance for the disease trait.}
  \item{\code{population_prevelance}}{population prevelance for the disease trait.}
  \ite{\code{stratification}}{logical; if TRUE, the population stratification effect is considered.}
}
\details{
If filter==TRUE, the input summary data will be filtered by 3 steps: (i) If sample size varies from SNP to SNP, remove SNPs with an effective sample size less than 0.67 times the 90th percentile of sample size. (ii) Remove SNPs within the major histocompatibility complex (MHR) region; filter SNPs to Hapmap3 SNPs. (iii) Remove SNPs with extremely large effect sizes (chi^2 > 80).
}
\value{
\item{\code{estimates}}{the list containing the parameter estimates of the fitted model.}
\item{\code{summaryGWASdataLD}}{the dataframe containing the filtered summary GWAS data to be input into the model.}
\item{\code{qqplotdatasave}}{the dataframe containing the simulated dataset to generate the QQ plot, i.e., the log10th of the observed p-values and expected p-values by the fitted model.}
}
\references{
Zhang, Yan, et al. "Estimation of complex effect-size distributions using summary-level statistics from genome-wide association studies across 32 complex traits and implications for the future." bioRxiv (2017): 175406.
}
\author{
Yan Zhang <yandorazhang@gmail.com>
}
\examples{
data(heightGWAS)
fit = genesis(heightGWAS,modelcomponents=2, cores=2, LDcutoff=0.1, LDwindow=1, c0=10,startingpic=0.01)
fit$estimates

fit3 = genesis(heightGWAS,modelcomponents=3, cores=20, LDcutoff=0.1, LDwindow=1, c0=10,startingpic=0.01)
fit3$estimates
}

