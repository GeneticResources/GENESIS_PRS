#' Illustration of postmarginal()
#'
#' This function allows to estimate posterior probability with a null marginal effect size and posterior mean of marginal effect size. 
#' @param est parameter estimates by fitting either 2-component model, i.e., (pic, sigmasq, a); or 3-component model, i.e., (pic, p1, sigmasq1, sigmasq2, a).
#' @param c0 an assumed maximum number of underlying susceptibility SNPs tagged by any individual GWAS marker. By default, c0 is set at 10.
#' @param summarydata summay-level GWAS data, containing 3 columns: 
#' SNP (SNP rsID), 
#' Z (GWAS test z-statistic), 
#' N (GWAS study sample size which can be different for different SNPs)
#' @param filter logical; if TRUE, the input summary data will be filtered.
#' @keywords 
#' @export
#' @examples numInterval(0.005,Inf,est=c(9.583307e-03,8.562964e-02,1.487684e-04,2.086576e-05,1.498790e-06))

postmarginal <- function(est, c0, summarydata, filter=F){
 
  #----------------------------------------------------#----------------------------------------------------
  # I. preliminary summary GWAS data filtering. 
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
  # II. merge the summary GWAS data with the LD score data, and extract the variables needed for analysis.
  #----------------------------------------------------#----------------------------------------------------
  
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
  df <- merge(summarydata, dataLD,by.x="SNP",by.y="SNPname",sort=F)
  
  
  
  betahat <- as.numeric(as.character(df$betahat))
  varbetahat <- as.numeric(as.character(df$varbetahat))
  ldscore <- as.numeric(as.character(df$LD.score.correct))
  Nstar <- as.numeric(as.character(df$Nstar))
  SNP <- df$SNP
  TaggingSNPs <- df$TaggingSNPs
 
  
  zero.omit <- function(v){v[which(v!=0)]}
  SNPlist <- unlist(lapply(dataLD$SNPname, function(x) as.numeric(as.character(strsplit(x,"rs")[[1]][2]))))
  dictionary <- modification_loc(SNPlist,length(SNPlist),max(SNPlist))
  TaggingSNPinx <- lapply(dataLD$TaggingSNPs, function(t) {zero.omit(dictionary[as.vector(na.omit(as.numeric(unlist(strsplit(strsplit(t, ",")[[1]], "rs")))))])})
  pairwise_r <- lapply(dataLD$pairwise_r, function(t) {as.vector(as.numeric(unlist(strsplit(t, ",")[[1]])))})
  
  
  
  if(length(est)==3) components=2
  if(length(est)==5) components=3
  
  if(components==2){
    q0 = (1-est[1])^Nstar * 1/sqrt(2*pi*(varbetahat)) * exp(-0.5*betahat^2/(varbetahat))
    beta_exp = 0 # posterior expectation of marginal effect size (zero when Nk0=Nstar)
    # each SNP has a different N_k^*, therefore goes thru different number of loops.
    qsum = q0
    for (i in 1:c0){ # i is Nk1: first type of causal marker
      # print(i)
      qtemp = rep(0, length(Nstar))
      delta = i*est[2]*ldscore/Nstar 
      ind = Nstar>=i
      
      qtemp[ind] = choose(Nstar[ind],i) * 
        (1-est[1])^(Nstar[ind]-i)*(est[1]^i) / sqrt(2*pi*(delta[ind]+varbetahat[ind])) *
        exp(-0.5*betahat[ind]^2/(delta[ind]+varbetahat[ind]))
      
      qsum = qsum + qtemp # Summation of q's the denominator
      beta_exp = beta_exp + qtemp * (varbetahat)^(-1)*betahat/(1/delta+(varbetahat)^(-1))
    }
    
    p_beta0 = q0 / qsum
    beta_exp = beta_exp / qsum
    
  }
  
  if(components==3){
    q0 = (1-est[1])^Nstar * 1/sqrt(2*pi*(varbetahat)) * exp(-0.5*betahat^2/(varbetahat))
    beta_exp = 0 # posterior expectation of marginal effect size (zero when Nk0=Nstar)
    # each SNP has a different N_k^*, therefore goes thru different number of loops.
    qsum = q0
    for (i in 0:c0){ # i is Nk1: number of the first type of causal SNPs
      # print(i)
      for (j in 0:c0){ # j is Nk2: number of the second type of causal SNPs
        if (i==0 & j==0){
          next
        } else{
          print(paste(i,j))
          qtemp = rep(0, length(Nstar))
          delta = (i*est[3]+j*est[4])*ldscore/Nstar 
          ind = Nstar>=(i+j)
          
          qtemp[ind] = choose(Nstar[ind],i)*choose(Nstar[ind]-i,j) * 
            ((1-est[1])^(Nstar[ind]-i-j))*((est[2]*est[1])^i)*(((1-est[2])*est[1]))^j * 1/sqrt(2*pi*(delta[ind]+varbetahat[ind])) *
            exp(-0.5*betahat[ind]^2/(delta[ind]+varbetahat[ind]))
          
          qsum = qsum + qtemp # Summation of q's the denominator
          beta_exp = beta_exp + qtemp * (varbetahat)^(-1)*betahat/(1/delta+(varbetahat)^(-1))
        }
      }
    }
    p_beta0 = q0 / qsum
    beta_exp = beta_exp / qsum
  }
  posterior = list()
  posterior$p_beta0 = p_beta0
  posterior$beta = beta_exp
  return(posterior)
}
