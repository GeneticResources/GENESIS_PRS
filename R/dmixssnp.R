#' Illustration of dmixssnp()
#'
#' This function allows to get the (mixture) normal probability density function for the susceptibility SNP (sSNP).
#' @param x value of quantile. 
#' @param est parameter estimates by fitting either 2-component model, i.e., (pic, sigmasq, a); or 3-component model, i.e., (pic, p1, sigmasq1, sigmasq2, a).
#' @keywords 
#' @export
#' @examples dmixssnp(x,est=c(8.899809e-03, 9.476025e-02, 1.458650e-04, 2.227118e-05, 1.567643e-06))

dmixssnp <- function(x,est){
  
  if(length(est)==5){
    pic = est[1]
    p0 = est[2]
    s1 = sqrt(est[3])
    s2 = sqrt(est[4])
    den <- function(x){return((p0 * dnorm(x/s1)/s1 + (1-p0)*dnorm(x/s2) /s2))}
  }
  
  if(length(est)==3){
    pic = est[1]
    s1 = sqrt(est[2])
    den <- function(x){return(dnorm(x/s1)/s1)}
  }
  return(den(x))
}