#' Illustration of projection()
#'
#' This function allows to make future projections according to the fitted model.
#' @param est parameter estimates by fitting either 2-component model, i.e., (pic, sigmasq, a); or 3-component model, i.e., (pic, p1, sigmasq1, sigmasq2, a).
#' @param n specifided future GWAS sample size.
#' @param alpha significance level. 
#' @param tol tolerance accuracy vector for intgrate() function.
#' @param M total number of SNPs in the reference panel; by default, it is the total number of common SNPs in Hapmap3 reference panel, which is equal to 1070777. 
#' @keywords 
#' @export
#' @examples projection(est=c(9.583307e-03,8.562964e-02,1.487684e-04,2.086576e-05,1.498790e-06), n=253288)

projection <- function(est,n,alpha=5e-8,tol=c(1e-12,1e-15),M=1070777){
  
  if(length(est)==3)components=2
  if(length(est)==5)components=3
  
  if(components==2){
    pic = est[1]
    sig = sqrt(est[2])
    den <- function(x){return(dnorm(x/sig)/sig )}
    herit0 <- pic*M*sig^2
  }
  
  if(components==3){
    pic = est[1]
    p1 = est[2]
    s1 = sqrt(est[3])
    s2 = sqrt(est[4])
    den <- function(x){return(p1 * dnorm(x/s1)/s1 + (1-p1)*dnorm(x/s2) /s2)}
    herit0 <- pic*M*(p1*est[3] + (1-p1)*est[4])
  }
  
  tem0 <- function(x){return(x^2*den(x))}
  ebeta2 <- integrate(tem0, -Inf, Inf,rel.tol=tol[1], abs.tol=tol[2])[[1]]
  
  c_alpha = abs(qnorm(alpha/2))
  pow <- function(x){return(1 - pnorm(c_alpha - sqrt(n)*x) + pnorm(-c_alpha - sqrt(n)*x) )}
  tem <- function(x){return(pow(x)*den(x))}
  Numdiscoveries = M*pic * integrate(tem, -Inf, Inf,rel.tol=tol[1], abs.tol=tol[2])[[1]]
  
  tem1 <- function(x){return(pow(x)*den(x)*x^2)}
  GVpercentage = M*pic * integrate(tem1, -Inf, Inf,rel.tol=tol[1], abs.tol=tol[2])[[1]]*100/herit0
  
  return(list(Numdicoveries=Numdiscoveries, GVpercentage=GVpercentage))
}
