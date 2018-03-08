#' Illustration of qqfigure()
#'
#' This function allows to predict future genomic control (GC) factor through simulations. 
#' @param qqplotdata The QQ plot object got from genesis() function.
#' @param seq_inx The QQ data to get the QQ plot will be thinned every seq_inx. 
#' @param qqaxis Numeric; the x- and y-axis limits is from 0 to qqaxis for the QQ plot. By default, it is 10. 
#' @keywords 
#' @export
#' @examples qqfigure(qqplotdata)

qqfigure <- function(qqplotdata,seq_inx=10,qqaxis=10){
  
  QQdata = data.frame(qqplotdata$data)
  obs_lambda = fit$qqplotdata$observedlambda
  m.lambda = fit$qqplotdata$meanEXPlambda
  l.lambda = fit$qqplotdata$lowEXPlambda
  h.lambda = fit$qqplotdata$highEXPlambda
  
  # ------------------------------------------------
  # QQ plot
  inx <- seq(1,nrow(QQdata),seq_inx)
  QQdata <- QQdata[inx,]
  plot(QQdata$mean_log_exp_pvalues, QQdata$log_obs_pvalues, type="l",xlab=expression(Expected~~-log[10](italic(P)~value)), xlim=c(0,qqaxis),ylim=c(0,qqaxis),ylab=expression(Observed~~-log[10](italic(P)~value)))
  polygon(c(QQdata$lower,rev(QQdata$upper)),c(QQdata$log_obs_pvalues,rev(QQdata$log_obs_pvalues)),col = "grey75", border = FALSE)
  # points(QQdata$mean_log_exp_pvalues, QQdata$log_obs_pvalues)
  abline(a=0,b=1)
  #add red lines on borders of polygon
  lines(y=QQdata$log_obs_pvalues,col="red", QQdata$upper,lty=2)
  lines(y=QQdata$log_obs_pvalues,col="red", QQdata$lower,lty=2)
  
  mylabel <- bquote(italic(lambda[obs])== .(formatC(obs_lambda,format="f", digits = 3)))
  text(x = 8.5, y = 2.5, labels = mylabel)
  
  mylabel <- bquote(italic(lambda[fit])== .(formatC(m.lambda,format="f", digits = 3)))
  text(x = 8.5, y = 1.5, labels = mylabel)
}
  