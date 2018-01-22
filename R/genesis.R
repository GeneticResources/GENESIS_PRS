#' Illustration of genesis()
#'
#' This function allows to get parameter estimates from fitting the mixture model.
#' @param summarydata susummay-level GWAS data, containing 3 columns: 
#' SNP (SNP rsID), 
#' Z (GWAS test z-statistic), 
#' N (GWAS study sample size which can be different for different SNPs)
#' @param modelcomponents 2 or 3, indicating fitting 2-component or 3-component model.
#' @param cores Number of CPU threads in parallel computing.
#' @param LDcutoff LDcutoff=1, LD score is calculated based on cutoff 0.1; cutoff=2, 0.2.
#' @param c0 An assumed maximum number of underlying susceptibility SNPs tagged by any individual GWAS marker. By default, c0 is set at 10.
#' @param BICgamma A tuning parameter in calculating BIC in (0,1). By default, BICgamma is set at 0.5.
#' @param print Logical; if TRUE, the EM algorithm details will be output.
#' @param printfreq  Number indicating every printfreq steps, the EM results will be output.
#' @param starting The starting values for the model. For 2-component model, the starting value of (pic, sigmasq, a); for 3-component model, (pic, p1, sigmasq1, sigmasq2, a). 
#' @param staringpic The starting value for pic when staring==NA. 
#' @param tolerance The accuracy of the tolerance. For 2-component model, it is a 6-dim vector with tolerance value for (pic,sigmasq,a,llk,maxEM,steps). For 3-component model, it is a 8-dim vector with tolerance value for (pic,p1,sigmasq1,sigmasq2,a,llk,maxEM,steps).
#' @param qqplot Logical; if TRUE, the QQ plot will be ploted.
#' @param qqplotCI The threshold of confidence interval in the QQ plot.
#' @param qqplotname The name of the QQ plot pdf files. 
#' @param nsim The total number of simulations based on which the QQ plot is obtained. 
#' @param summaryGWASLDdatasave Logical; if TRUE, the summary GWAS data as well as the LD information will be saved.
#' @param qqplotdatasave Logical; the simulated data to generate the QQ plot will be saved.
#' @param siblingrisk Logical; if TRUE, the sibling risk will also be calculated.
#' @param filter Logical; if TRUE, the input summary data will be filtered.
#' @keywords 
#' @export
#' 
#' @examples
#' genesis(summarydata,modelcomponents=2, cores=24, LDcutoff=1, c0=10, BICgamma=0.5,print=TRUE,printfreq=10, starting=NA,startingpic=0.01, tolerance=NA,qqplot=TRUE, qqplotCI=0.8, qqplotname=paste0(""),summaryGWASLDdatasave=FALSE,qqplotdatasave=T,siblingrisk=FALSE,filter=F)
                     
genesis <- function(summarydata, modelcomponents=2, cores=10, LDcutoff=1, c0=10, BICgamma=0.5, print=TRUE, printfreq=10, starting=NA, startingpic=NA, tolerance=NA, qqplot=TRUE, qqplotCI=0.8, qqplotname="", nsim=100, summaryGWASLDdatasave=FALSE, qqplotdatasave=FALSE, siblingrisk=FALSE, filter=FALSE){
  
  #----------------------------------------------------#----------------------------------------------------
  # I. Input values check 
  #----------------------------------------------------#----------------------------------------------------
  P = 5 # the total number of parameters in the 3-component model. 
  #----------------------------------------------------
  # summary GWAS data format check
  if(ncol(summarydata)!=3){
    stop("the summary GWAS data should have 3 columns with (SNP rsID, Z-statistic, Sample size)!")
  }
  
  if(any(!is.na(starting))){
    l = length(starting)
    if(modelcomponents==2){if(l!=3) stop("2-component model: the input starting values should be a 3-dim vector with value (pic,sigmasq,a)!")} 
    if(modelcomponents==3){if(l!=5) stop("3-component model: the input starting values should be a 5-dim vector with value (pic,p1,sigmasq1,sigmasq2,a)!")} 
  }
  
  #----------------------------------------------------
  # tolerance check 
  if(all(is.na(tolerance))){
    tolerance_pic=1e-6;
    tolerance_p1=1e-5;
    tolerance_sigsq1=1e-8;
    tolerance_sigsq2=1e-8;
    tolerance_a=1e-9;
    tolerance_llk=1e-6;
    maxEM=3e3;
    steps=2;
  }
  
  if(any(!is.na(tolerance))){
    l = length(tolerance)
    if(modelcomponents==2){
      if(l!=6) stop("2-component model: the tolerance values should be a 6-dim vector with tolerance value for (pic,sigmasq,a,llk,maxEM,steps)!")
      if(l==6) {tolerance_pic=tolerance[1]; tolerance_sigsq1=tolerance[2]; tolerance_a=tolerance[3]; tolerance_llk=tolerance[4]; maxEM=tolerance[5]; steps=tolerance[6];}
    } 
    if(modelcomponents==3){
      if(l!=8) stop("3-component model: the tolerance values should be a 8-dim vector with tolerance value for (pic,p1,sigmasq1,sigmasq2,a,llk,maxEM,steps)!")
      if(l==8) {tolerance_pic=tolerance[1]; tolerance_p1=tolerance[2]; tolerance_sigsq1=tolerance[3]; tolerance_sigsq2=tolerance[4]; tolerance_a=tolerance[5]; tolerance_llk=tolerance[6]; maxEM=tolerance[7]; steps=tolerance[8];}
    }   
  }
  
  #----------------------------------------------------
  # load the required R package
  library(doParallel)
  library(foreach)
  cl <- makeCluster(cores)  
  registerDoParallel(cl)  
  
  #----------------------------------------------------#----------------------------------------------------
  # II. preliminary summary GWAS data filtering. 
  #----------------------------------------------------#----------------------------------------------------
  colnames(summarydata) <- c("SNP","Z","N")
  
  if(filter==TRUE){
    #a. If sample size varies from SNP to SNP, remove SNPs with an effective sample size less than 0.67 times the 90th percentile of sample size.
    ikeep1 <- which(as.numeric(as.character(summarydata$N))>=0.67*quantile(as.numeric(as.character(summarydata$N)), 0.9))
    summarydata <- summarydata[ikeep1,]
    
    #b. Remove SNPs within the major histocompatibility complex (MHR) region; filter SNPs to Hapmap3 SNPs.
    data(w_hm3.noMHC.snplist)
    ikeep2 <- which(as.character(summarydata$SNP) %in% w_hm3.noMHC.snplist$SNP)
    summarydata <- summarydata[ikeep2,]
    
    #c. Remove SNPs with extremely large effect sizes (chi^2 > 80).
    ikeep3 <- which(as.numeric(as.character(summarydata$Z))^2 <=80)
    summarydata <- summarydata[ikeep3,]
  }

  #----------------------------------------------------#----------------------------------------------------
  # III. merge the summary GWAS data with the LD score data, and extract the variables needed for analysis.
  #----------------------------------------------------#----------------------------------------------------
  M <- 1070777 # the total number of SNPs in the reference panel used for calculating LD score (r^2>0.1), i.e., the Hapmap3 panel
  
  if(LDcutoff==1){data(LD_hm3_cutoff01); LDdata = LD_hm3_cutoff01}
  if(LDcutoff==2){data(LD_hm3_cutoff02); LDdata = LD_hm3_cutoff02}
  
  summarydata$SNP <- as.character(summarydata$SNP)
  summarydata$Z <- as.numeric(as.character(summarydata$Z))
  summarydata$N <- as.numeric(as.character(summarydata$N))
  
  # remove NA values.
  inx <- which(is.na(summarydata$Z))
  if(length(inx)>1) summarydata <- summarydata[-inx,]
  inx <- which(is.na(summarydata$N))
  if(length(inx)>1) summarydata <- summarydata[-inx,]
  
  # get the variable needed for our method.
  summarydata$betahat <- summarydata$Z/sqrt(summarydata$N)
  summarydata$varbetahat <- 1/summarydata$N
  df <- merge(summarydata, LDdata,by.x="SNP",by.y="SNPname",sort=F)

  betahat <- as.numeric(as.character(df$betahat))
  varbetahat <- as.numeric(as.character(df$varbetahat))
  ldscore <- as.numeric(as.character(df$LD.score.correct))
  Nstar <- as.numeric(as.character(df$Nstar))
  SNP <- df$SNP
  TaggingSNPs <- df$TaggingSNPs
  K <- length(betahat)
  n <- as.numeric(as.character(df$N))

  N_SNPs_summary <- length(betahat)
  #----------------------------------------------------#----------------------------------------------------
  # IV. starting value
  #----------------------------------------------------#----------------------------------------------------
  if(all(is.na(starting))){
    #----------------------------------------------------
    # LD score regression to get the starting values. 
    chisq <- (betahat/sqrt(varbetahat))^2
    tem <- ldscore*mean(1/varbetahat)/M
    ld.fit <- lm(chisq ~ tem)
    
    if(modelcomponents==2){
      starting <- rep(0,3)
      if(is.na(startingpic)){starting[1] <- 0.01}
      if(!is.na(startingpic)){starting[1] <- startingpic}
      starting[2] <- ld.fit$coefficients[2]/(M*starting[1])
      starting[3] <- (ld.fit$coefficients[1]-1)/M
      if(starting[2]<0){starting[2] <- 1e-5; starting[3] <- 2e-6}
      
      llk = loglikelihood(starting, betahat,varbetahat, ldscore,c0,Nstar,cores)
      if(is.na(llk)){
        if(is.na(startingpic)){starting[1] <- 5e-3}
        if(!is.na(startingpic)){starting[1] <- startingpic}
        starting[2] <- ld.fit$coefficients[2]/(M*starting[1])
        starting[3] <- (ld.fit$coefficients[1]-1)/M
        if(starting[2]<0){starting[2] <- 1e-5; starting[3] <- 2e-6}
      }
    }
    
    if(modelcomponents==3){
      #----------------------------------------------------
      # run 2-component model to get a rough starting value
      starting <- rep(0,3)
      if(is.na(startingpic)){starting[1] <- 0.01}
      if(!is.na(startingpic)){starting[1] <- startingpic}
      starting[2] <- ld.fit$coefficients[2]/(M*starting[1])
      starting[3] <- (ld.fit$coefficients[1]-1)/M
      if(starting[2]<0){starting[2] <- 1e-5; starting[3] <- 2e-6}
      llk = loglikelihood(starting, betahat,varbetahat, ldscore,c0,Nstar,cores)
      if(is.na(llk)){
        if(is.na(startingpic)){starting[1] <- 5e-3}
        if(!is.na(startingpic)){starting[1] <- startingpic}
        starting[2] <- ld.fit$coefficients[2]/(M*starting[1])
        starting[3] <- (ld.fit$coefficients[1]-1)/M
        if(starting[2]<0){starting[2] <- 1e-5; starting[3] <- 2e-6}
      }
      
      fit <- EM_func(starting, betahat,varbetahat,ldscore,Nstar,M, c0, max(tolerance_pic, 1e-6), max(tolerance_sigsq1,1e-7),max(tolerance_a,1e-7),max(tolerance_llk,1e-5),maxEM,steps,cores,print=F,printfreq)
      est <- fit[3:5];
      
      starting <- rep(0,5)
      starting[1] <- est[1]
      starting[2] <- 1.0/9
      starting[3] <- est[2]*5
      starting[4] <- starting[3]/10
      starting[5] <- est[3]
    }
  }
  zero.omit <- function(v){v[which(v!=0)]}
  
  #----------------------------------------------------#----------------------------------------------------
  # V. Data analysis, 2-component model. 
  #----------------------------------------------------#----------------------------------------------------
  if(modelcomponents == 2){
    
    time0 <- proc.time()[3]
    fit <- EM_func(starting, betahat,varbetahat,ldscore,Nstar,M, c0, tolerance_pic, tolerance_sigsq1,tolerance_a,tolerance_llk,maxEM,steps,cores,print,printfreq)
    est <- fit[3:5]; c0 = fit[7]
    runhour <- (proc.time()[3] - time0)/3600
    
    causalnum <- est[1]*M
    heritability <- est[1]*est[2]*M  
    
    #----------------------------------------------------
    # calculate variance
    #----------------------------------------------------
    time1 <- proc.time()[3]
    m_S <- SS(est,betahat, varbetahat, ldscore, c0, Nstar, cores); # score matrix K*3
    m_I <- I(est,betahat, varbetahat, ldscore, c0, Nstar, cores); # information matrix 3*3
    
    #----------------------------------------------------
    # get sum of scores in each neighbor of SNP, K*3 matrix
    m_Sbar <- matrix(0,K,length(est));
    inx_name <- apply(matrix(SNP,ncol=1), 1, function(t) as.numeric(strsplit(t, "rs")[[1]][2]))
    dictionary <- modification_loc(inx_name,K,max(inx_name)) 
    tem <- lapply(TaggingSNPs, function(t) {inx <- zero.omit(dictionary[as.vector(na.omit(as.numeric(unlist(strsplit(strsplit(t, ",")[[1]], "rs")))))]); colSums(matrix(m_S[inx,],ncol=length(est)))})
    m_Sbar = matrix(unlist(tem),ncol=ncol(m_S),byrow=T) + m_S

    #----------------------------------------------------
    inv_I <- solve(m_I,tol=1e-20);    
    J <- (t(m_S)%*%m_Sbar);
    var_est <- inv_I %*% J %*% inv_I; # variance matrix of parameter est
    sd_est <- sqrt(diag(var_est)); # standard error for each parameter estimate
    
    #----------------------------------------------------
    # calculate AIC, BIC - d_s
    #----------------------------------------------------
    llk = fit[2]
    ds = sum(diag(inv_I %*% J))
    aic = -2*llk + 2*ds
    bic = -2*llk + ds*log(mean(n)) + 2*BICgamma*log(P^ds)
   
    temtem <- matrix(c(M*est[2], M*est[1], 0), ncol=1)
    sd_heritability <- sqrt( t(temtem) %*% var_est %*% temtem)  # standard error of heritability
    sd_causalnum <- M*sd_est[1]
    runhour_var <- (proc.time()[3]-time1)/3600
    
    
    if(siblingrisk==T){
      risk <- sqrt(exp(heritability))
      sd_risk <- sqrt(risk*sd_heritability^2/2)
      
      estimates <- list("Number of sSNPs (sd)"=paste0(format(causalnum,digits=3)," (",format(sd_causalnum,digits=4), ")"),
                        "Total heritability (sd)"=paste0(format(heritability,digits=3)," (",format(sd_heritability,digits=4), ")"),
                        "Sibling risk (sd)"=paste0(format(risk,digits=3)," (",format(sd_risk,digits=4), ")"),
                        "parameter (pic, sigmasq, a) estimates" = est,
                        "S.D. of parameter estimates"=sd_est,
                        "Covariance matrix of parameter estimates"=var_est,
                        "log-likelihood of fitted model" = fit[2],
                        "AIC" = aic,
                        "BIC" = bic,
                        "ds" = ds, "c0" = c0, "Information matrix" = m_I, "J"=J,
                        "Total number of SNPs in the Hapmap3 reference panel"=M,
                        "Total number of SNPs in the GWAS study after quality control"=N_SNPs_summary
                        )
    }
    
    if(siblingrisk==F){
      estimates <- list("Number of sSNPs (sd)"=paste0(format(causalnum,digits=3)," (",format(sd_causalnum,digits=4), ")"),
                        "Total heritability (sd)"=paste0(format(heritability,digits=3)," (",format(sd_heritability,digits=4), ")"),
                        "parameter (pic, sigmasq, a) estimates" = est,
                        "S.D. of parameter estimates"=sd_est,
                        "Covariance matrix of parameter estimates" = var_est,
                        "log-likelihood of fitted model" = fit[2],
                        "AIC" = aic,
                        "BIC" = bic,
                        "ds" = ds,"c0" = c0, "Information matrix" = m_I,"J"=J,
                        "Total number of SNPs in the Hapmap3 reference panel"=M,
                        "Total number of SNPs in the GWAS study after quality control"=N_SNPs_summary
                        )
    }
    
    
    if(qqplot==T){
      a <- est[3]; if(a<0) a <- 0; 
      # ------------------------------------------------
      obs_z <- betahat/sqrt(varbetahat)
      obs_pvalues <- 2*pnorm(-abs(obs_z))
      obs_lambda <- median(obs_z^2)/0.456
      log_obs_pvalues <- -log10(obs_pvalues)
      log_obs_pvalues <- sort(log_obs_pvalues)
      # ------------------------------------------------
      SNPsum = SNP; data(list=paste0("error_iter1")); error.snplist = SNP; SNP = SNPsum
      inx <- which(!is.na(match(error.snplist, SNP)))
      te <- mixture_components_marginal(est, ldscore, c0, Nstar, cores)
      proportions <- te$proportions
      varcomponents <- te$varcomponents
      L <- ncol(proportions)
      # ------------------------------------------------
      log_exp_pvalues <- matrix(0,nsim,K);
      exp_z <- matrix(0,nsim,K)
      exp_lambda <- rep(0,nsim)
      temorder <- order(match(error.snplist, SNP))
      # ------------------------------------------------
      foreach(i=1:nsim)%do%{
        data(list=paste0("error_iter",i))
        set.seed(123*i)
        # ------------------------------------------------
        # generate causalnum causal SNPs, and hence then the true SNP effect size. 
        betamarginal <- rep(0,K)
        temtem <- for(k in 1:K){
          components <- sample(1:L, prob=proportions[k,],size=1, replace=T)
          mus <- rep(0, L)
          sds <- sqrt(varcomponents[k,])
          betamarginal[k] <- rnorm(n=1,mean=mus[components],sd=sds[components])
        }

        betahat <- betamarginal + error[temorder][1:K] /sqrt(n) + rnorm(1,mean=0,sd=sqrt(a))
        exp_z[i,] <- betahat*sqrt(n)
        log_exp_pvalues[i,] <- -log10(2*pnorm( -abs(exp_z[i,])) )
        exp_lambda[i] <- median(exp_z[i,]^2)/0.456
        log_exp_pvalues[i,] <- sort(log_exp_pvalues[i,])
      }
      
      mean_log_exp_pvalues <- apply(log_exp_pvalues, 2, mean) 
      lower <- apply(log_exp_pvalues, 2, function(t) quantile(t, (1-qqplotCI)/2)) 
      upper <- apply(log_exp_pvalues, 2, function(t) quantile(t, 1-(1-qqplotCI)/2))
      
      m.lambda <- mean(exp_lambda);
      l.lambda <- quantile(exp_lambda, (1-qqplotCI)/2)
      h.lambda <- quantile(exp_lambda, 1-(1-qqplotCI)/2)
      
      QQdata = data.frame(cbind(log_obs_pvalues,mean_log_exp_pvalues,lower,upper))
      # ------------------------------------------------
      # QQ plot
      pdf(file=paste0(qqplotname,"qq2com.pdf"))
      inx <- seq(1,nrow(QQdata),10)
      QQdata <- QQdata[inx,]
      plot(QQdata$mean_log_exp_pvalues, QQdata$log_obs_pvalues, type="l",xlab=expression(Expected~~-log[10](italic(P)~value)), xlim=c(0,10),ylim=c(0,10),ylab=expression(Observed~~-log[10](italic(P)~value)))
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
      dev.off()
      
      qqplotdata <- list(data=QQdata, observedlambda=obs_lambda,
                       meanEXPlambda=m.lambda, lowEXPlambda=l.lambda, highEXPlambda=h.lambda)
      
      if(qqplotdatasave==T & summaryGWASLDdatasave==T){result <- list(estimates=estimates,summaryGWASLDdata=df,qqplotdata=qqplotdata)}
      if(qqplotdatasave==T & summaryGWASLDdatasave==F){result <- list(estimates=estimates,qqplotdata=qqplotdata)}
      if(qqplotdatasave==F & summaryGWASLDdatasave==T){result <- list(estimates=estimates,summaryGWASLDdata=df)}
      if(qqplotdatasave==F & summaryGWASLDdatasave==F){result <- list(estimates=estimates)}
    }
    if(qqplot==F){
      result <- list(estimates=estimates)
    }
  }
  
  #----------------------------------------------------#----------------------------------------------------
  # 5. Data analysis, 3-component model. 
  #----------------------------------------------------#----------------------------------------------------
  if(modelcomponents == 3){
    #----------------------------------------------------
    # run 3-component model. 
    #----------------------------------------------------
    time0 <- proc.time()[3]
    fit <- EM_func3(starting,lower_pi=c(1e-7,1e-7),upper_pi=c(0.5,0.5),betahat,varbetahat,ldscore,Nstar,M,c0,tolerance_pic,tolerance_p1,tolerance_sigsq1,tolerance_sigsq2,tolerance_a,tolerance_llk,maxEM,steps,cores,print,printfreq)
    est <- fit[3:7]; c0 = fit[9]
    runhour <- (proc.time()[3] - time0)/3600
    llk = fit[2]
    
    causalnum <- est[1]*M
    largenum = M*est[1]*est[2]
    heritlarge = M*est[1]*est[2]*est[3]
    heritsmall = M*est[1]*(1-est[2])*est[4]
    heritability <- est[1]*(est[2]*est[3] + (1-est[2])*est[4])*M
    
    #----------------------------------------------------
    # calculate variance
    #----------------------------------------------------
    time1 <- proc.time()[3]
    m_S <- SS3(est,betahat, varbetahat, ldscore, c0, Nstar, cores); # score matrix K*3
    m_I <- I3(est,betahat, varbetahat, ldscore, c0, Nstar, cores); # information matrix 3*3
    
    #----------------------------------------------------
    # get sum of scores in each neighbor of SNP, K*3 matrix
    m_Sbar <- matrix(0,K,length(est));
    inx_name <- apply(matrix(SNP,ncol=1), 1, function(t) as.numeric(strsplit(t, "rs")[[1]][2]))
    dictionary <- modification_loc(inx_name,K,max(inx_name)) 
    tem <- lapply(TaggingSNPs, function(t) {inx <- zero.omit(dictionary[as.vector(na.omit(as.numeric(unlist(strsplit(strsplit(t, ",")[[1]], "rs")))))]); colSums(matrix(m_S[inx,],ncol=length(est)))})
    m_Sbar = matrix(unlist(tem),ncol=ncol(m_S),byrow=T) + m_S

    #----------------------------------------------------
    inv_I <- solve(m_I,tol=1e-20);    
    J <- (t(m_S)%*%m_Sbar);
    var_est <- inv_I %*% J %*% inv_I; # variance matrix of parameter est
    sd_est <- sqrt(diag(var_est)); # standard error for each parameter estimate
    
    #----------------------------------------------------
    # calculate AIC, BIC - d_s
    #----------------------------------------------------
    ds = sum(diag(inv_I %*% J))
    aic = -2*llk + 2*ds
    bic = -2*llk + ds*log(mean(n)) + 2*BICgamma*log(P^ds)
    
    temtem = matrix(c(M*(est[2]*est[3] + (1-est[2])*est[4]), 
                      est[1]*(est[3] - est[4])*M, 
                      est[1]*(est[2])*M,
                      est[1]*(1-est[2])*M,
                      0), ncol=1)
    sd_heritability = sqrt( t(temtem) %*% var_est %*% temtem) # standard error of heritability
    
    sd_causalnum = M*sd_est[1]
    tem_largenum = matrix(c(M*est[2], M*est[1], 0, 0, 0), ncol=1)
    sd_largenum = sqrt( t(tem_largenum) %*% var_est %*% tem_largenum)
    
    tem_heritlarge = matrix(c(M*est[2]*est[3], M*est[1]*est[3], M*est[1]*est[2], 0, 0),ncol=1)
    sd_heritlarge = sqrt( t(tem_heritlarge) %*% var_est %*% tem_heritlarge)
    
    tem_heritsmall = matrix(c(M*(1-est[2])*est[4],-M*est[1]*est[4], 0, M*est[1]*(1-est[2]), 0),ncol=1)
    sd_heritsmall = sqrt( t(tem_heritsmall) %*% var_est %*% tem_heritsmall)
    runhour_var <- (proc.time()[3]-time1)/3600
    
    if(siblingrisk==T){
      risk <- sqrt(exp(heritability))
      sd_risk <- sqrt(risk*sd_heritability^2/2)
      
      estimates <- list("Number of sSNPs (sd)"=paste0(format(causalnum,digits=3)," (",format(sd_causalnum,digits=4), ")"),
                        "Number of sSNPs in the cluster with larger variance component (sd)"=paste0(format(largenum,digits=3)," (",format(sd_largenum,digits=4), ")"),
                        "Total heritability (sd)"=paste0(format(heritability,digits=3)," (",format(sd_heritability,digits=4), ")"),
                        "Sibling risk (sd)"=paste0(format(risk,digits=3)," (",format(sd_risk,digits=4), ")"),
                        "Heritability explained by the cluster with larger variance component (sd)"=paste0(format(heritlarge,digits=3)," (",format(sd_heritlarge,digits=4), ")"),
                        "Heritability explained by the cluster with samller variance component"=paste0(format(heritsmall,digits=3)," (",format(sd_heritsmall,digits=4), ")"),
                        "parameter (pic, p1, sigmasq1, sigmasq2, a) estimates" = est, 
                        "S.D. of parameter estimates"=sd_est,
                        "Covariance matrix of parameter estimates" = var_est,
                        "log-likelihood of fitted model" = fit[2],
                        "AIC" = aic,
                        "BIC" = bic,
                        "ds" = ds,"c0" = c0, "Information matrix" = m_I,"J"=J,
                        "Total number of SNPs in the Hapmap3 reference panel"=M,
                        "Total number of SNPs in the GWAS study after quality control"=N_SNPs_summary
                        )
    }
    
    if(siblingrisk==F){
      estimates <- list("Number of sSNPs (sd)"=paste0(format(causalnum,digits=3)," (",format(sd_causalnum,digits=4), ")"),
                        "Number of sSNPs in the cluster with larger variance component (sd)"=paste0(format(largenum,digits=3)," (",format(sd_largenum,digits=4), ")"),
                        "Total heritability (sd)"=paste0(format(heritability,digits=3)," (",format(sd_heritability,digits=4), ")"),
                        "Heritability explained by the cluster with larger variance component (sd)"=paste0(format(heritlarge,digits=3)," (",format(sd_heritlarge,digits=4), ")"),
                        "Heritability explained by the cluster with samller variance component"=paste0(format(heritsmall,digits=3)," (",format(sd_heritsmall,digits=4), ")"),
                        "parameter (pic, p1, sigmasq1, sigmasq2, a) estimates" = est, 
                        "S.D. of parameter estimates"=sd_est,
                        "Covariance matrix of parameter estimates" = var_est,
                        "log-likelihood of fitted model" = fit[2],
                        "AIC" = aic,
                        "BIC" = bic,
                        "ds" = ds,"c0" = c0, "Information matrix" = m_I,"J"=J,
                        "Total number of SNPs in the Hapmap3 reference panel"=M,
                        "Total number of SNPs in the GWAS study after quality control"=N_SNPs_summary
                        )
    }
    
    if(qqplot==T){
      a <- est[5]; if(a<0) a <- 0; 
      # ------------------------------------------------
      obs_z <- betahat/sqrt(varbetahat)
      obs_pvalues <- 2*pnorm(-abs(obs_z))
      obs_lambda <- median(obs_z^2)/0.456
      log_obs_pvalues <- -log10(obs_pvalues)
      log_obs_pvalues <- sort(log_obs_pvalues)
      # ------------------------------------------------
      SNPsum = SNP; data(list=paste0("error_iter1")); error.snplist = SNP; SNP = SNPsum
      inx <- which(!is.na(match(error.snplist, SNP)))
      te <- mixture_3components_marginal(est, ldscore, c0, Nstar, cores)
      proportions <- te$proportions
      varcomponents <- te$varcomponents
      L <- ncol(proportions)
      # ------------------------------------------------
      log_exp_pvalues <- matrix(0,nsim,K);
      exp_z <- matrix(0,nsim,K)
      exp_lambda <- rep(0,nsim)
      temorder <- order(match(error.snplist, SNP))
      # ------------------------------------------------
      foreach(i=1:nsim)%do%{
        data(list=paste0("error_iter",i))
        set.seed(123*i)
        # ------------------------------------------------
        # generate causalnum causal SNPs, and hence then the true SNP effect size. 
        betamarginal <- rep(0,K)
        temtem <- for(k in 1:K){
          components <- sample(1:L, prob=proportions[k,],size=1, replace=T)
          mus <- rep(0, L)
          sds <- sqrt(varcomponents[k,])
          betamarginal[k] <- rnorm(n=1,mean=mus[components],sd=sds[components])
        }
        
        betahat <- betamarginal + error[temorder][1:K] /sqrt(n) + rnorm(1,mean=0,sd=sqrt(a))
        exp_z[i,] <- betahat*sqrt(n)
        log_exp_pvalues[i,] <- -log10(2*pnorm( -abs(exp_z[i,])) )
        exp_lambda[i] <- median(exp_z[i,]^2)/0.456
        log_exp_pvalues[i,] <- sort(log_exp_pvalues[i,])
      }
      
      mean_log_exp_pvalues <- apply(log_exp_pvalues, 2, mean) 
      lower <- apply(log_exp_pvalues, 2, function(t) quantile(t, (1-qqplotCI)/2)) 
      upper <- apply(log_exp_pvalues, 2, function(t) quantile(t, 1-(1-qqplotCI)/2))
      
      m.lambda <- mean(exp_lambda);
      l.lambda <- quantile(exp_lambda, (1-qqplotCI)/2)
      h.lambda <- quantile(exp_lambda, 1-(1-qqplotCI)/2)
      
      QQdata = data.frame(cbind(log_obs_pvalues,mean_log_exp_pvalues,lower,upper))
      # ------------------------------------------------
      # QQ plot
      pdf(file=paste0(qqplotname,"qq3com.pdf"))
      inx <- seq(1,nrow(QQdata),10)
      QQdata <- QQdata[inx,]
      plot(QQdata$mean_log_exp_pvalues, QQdata$log_obs_pvalues, type="l",xlab=expression(Expected~~-log[10](italic(P)~value)), xlim=c(0,10),ylim=c(0,10),ylab=expression(Observed~~-log[10](italic(P)~value)))
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
      dev.off()
      
      qqplotdata <- list(data=QQdata, observedlambda=obs_lambda,
                         meanEXPlambda=m.lambda, lowEXPlambda=l.lambda, highEXPlambda=h.lambda)
      
      if(qqplotdatasave==T & summaryGWASLDdatasave==T){result <- list(estimates=estimates,summaryGWASLDdata=df,qqplotdata=qqplotdata)}
      if(qqplotdatasave==T & summaryGWASLDdatasave==F){result <- list(estimates=estimates,qqplotdata=qqplotdata)}
      if(qqplotdatasave==F & summaryGWASLDdatasave==T){result <- list(estimates=estimates,summaryGWASLDdata=df)}
      if(qqplotdatasave==F & summaryGWASLDdatasave==F){result <- list(estimates=estimates)}
    }
    if(qqplot==F){
      result <- list(estimates=estimates)
    }
  }
  return(result)
}