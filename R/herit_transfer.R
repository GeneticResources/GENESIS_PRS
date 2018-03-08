#' Illustration of herit_transfer
#'
#' This function allows to make future projections according to the fitted model.
#' @param h2 heritability in log-odds-ratio scale. 
#' @param se standard error of heritability in log-odds-ratio scale.
#' @param population_prevalence population prevalence.
#' @param sample_prevalence sample prevalence. 
#' @keywords 
#' @export
#' @examples herit_transfer(.351, .075, 0.044, 17/(17+37))

herit_transfer <- function(h2,se,population_prevalence,sample_prevalence){
  
  P <- sample_prevalence
  K <- population_prevalence
  z <- dnorm(qnorm(1-K))
  
  hobs <- h2*P*(1-P)
  se.hobs <- se*P*(1-P)

  hliab <- hobs*(K^2)*((1-K)^2)/((z^2)*P*(1-P))
  se.hliab <- se.hobs*(K^2)*((1-K)^2)/((z^2)*P*(1-P))
  
  return(list(Hobserved=hobs, se_Hobserved=se.hobs,
              Hliability=hliab,se_Hliability=se.hliab))
}