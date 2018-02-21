#' Illustration of projectionCI()
#'
#' This function allows to predict number of discoveries in future GWAS and its confidence interval through bootstrap method. 
#' @param est parameter estimates by fitting either 2-component model, i.e., (pic, sigmasq, a); or 3-component model, i.e., (pic, p1, sigmasq1, sigmasq2, a).
#' @param v covariance matrix of parameter estimates by fitting the 2- or 3-component model. 
#' @param n specifided future GWAS sample size.
#' @param alpha significance level of the confidence interval; by default, it is 0.025, i.e., 95% CI. 
#' @param M total number of SNPs in the reference panel; by default, it is the total number of common SNPs in Hapmap3 reference panel, which is equal to 1070777. 
#' @param nsim total number of simulations; by default, it is 1000. 
#' @keywords 
#' @export
#' @examples 
#' fdis(est=c(9.583307e-03,8.562964e-02,1.487684e-04,2.086576e-05,1.498790e-06),v,n=253288);

projectionCI <- function(est,v,n,alpha=0.025,M=1070777,nsim=10000){
  library(MASS)
  logest = log(est)
  logv = diag(1/est)%*%v%*%diag(1/est)
  estmat = exp(mvrnorm(nsim,mu=logest,Sigma=logv))
  
  pest = projection(est,n,M=M)$Numdicoveries;
  re = apply(estmat, 1, function(t) {projection(t,n,M=M)}$Numdicoveries)
  rere = apply(matrix(re,ncol=1),1,function(t) rbinom(1,size=M,prob=t/M))
  
  gvest = projection(est,n,M=M)$GVpercentage;
  regv = apply(estmat, 1, function(t) {projection(t,n,M=M)}$GVpercentage)

  return(list(Numdiscoveries = c(pest,quantile(rere,alpha),quantile(rere,1-alpha)),
              GVpercentage = c(gvest,quantile(regv,alpha),quantile(regv,1-alpha))))
}
