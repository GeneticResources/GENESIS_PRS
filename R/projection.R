#' Illustration of projection()
#'
#' This function allows to make future projections according to the fitted model.
#' @param est parameter estimates by fitting either 2-component model, i.e., (pic, sigmasq, a); or 3-component model, i.e., (pic, p1, sigmasq1, sigmasq2, a).
#' @param v covariance matrix of parameter estimates by fitting the 2- or 3-component model. 
#' @param n specifided future GWAS sample size.
#' @param gwsignificance genome-wide significance level, by default it is 5e-8. 
#' @param tol tolerance accuracy vector for intgrate() function.
#' @param M total number of SNPs in the reference panel; by default, it is the total number of common SNPs in Hapmap3 reference panel, which is equal to 1070777. 
#' @param nsim total number of simulations; by default, it is 10000.
#' @param alpha significance level of the confidence interval; by default, it is 0.025, i.e., 95% CI. 
#' @keywords 
#' @export
#' @examples projection(est=c(9.583307e-03,8.562964e-02,1.487684e-04,2.086576e-05,1.498790e-06),v,n=253288)

projection <- function(est,v,n,gwsignificance=5e-8,tol=c(1e-12,1e-15),M=1070777,nsim=10000,alpha=0.025){
  
  library(MASS)
  logest = log(est)
  logv = diag(1/est)%*%v%*%diag(1/est)
  estmat = exp(mvrnorm(nsim,mu=logest,Sigma=logv))
  
  tem = pp(est,n,gwsignificance,tol,M)
  tem1 = apply(estmat, 1, function(t) {pp(t,n,gwsignificance,tol,M)})
  
  pest = tem$Numdicoveries;
  gvest = tem(est,n,gwsignificance,tol,M)$GVpercentage;
  
  re = unlist(lapply(tem1,function(t)t[1]))
  rere = apply(matrix(re,ncol=1),1,function(t) rbinom(1,size=M,prob=t/M))
  regv = unlist(lapply(tem1,function(t)t[2]))
  
  return(list(Numdiscoveries = c(pest,quantile(rere,alpha),quantile(rere,1-alpha)),
              GVpercentage = c(gvest,quantile(regv,alpha),quantile(regv,1-alpha))))
}


