Examples
===

## Load summary GWAS data - Height (allen2010)
The summary GWAS dataframe has 3 columns: (SNPname, z-statistic, sample size). 

```{r height}
library(GENESIS)
data(heightGWAS)
```

## Fit the 2-component model
### Fitting to the model

Note the startingpic value can be specifided at a list of values, and then the one with the largest log-likelihood is selected as the final model. 

```{r 2-component model}
fit2 <- genesis(df,modelcomponents=2, cores=2, LDcutoff=1, c0=10,startingpic=0.005,filter=F)
fit2$estimates

est <- fit2$estimates$`parameter (pic, sigmasq, a) estimates` # the model parameter estimates
v <- fit2$estimates$`Covariance matrix of parameter estimates` # the covariance matrix of model parameter estimtaes

# est <- c(5.119735e-03, 5.533230e-05, 1.794370e-06)
# v <- matrix(c(3.415864e-07, -2.695963e-09, -6.943308e-11, -2.695963e-09,  2.576368e-11,  5.008517e-13, -6.943308e-11,  5.008517e-13,  1.881822e-14),3,3)
```

### Get the density plot for the susceptibility SNPs 
```{r density plot}
x_seq <- seq(-0.02,0.02,length.out = 1000); 
y_seq <- apply(matrix(x_seq,ncol=1),1,function(t) dcausal(t,est))
plot(x_seq, y_seq,type="l",ylim=c(0,250),xlab="Joint effect size", ylab="Probability Density")
```

### Make future projections with specified sample size n
```{r future projections}
projection(est,n=253288)
```

### Predict number of discoveries in future GWAS and its confidence interval
```{r future projections}
fdis(est,v,n=253288)
```

### Calculate number of SNPs falling in an interval
```{r number of SNPs in an interval}
numInterval(0.005,Inf,est)
```

### Predict genomic control factor in future GWAS with specified sample size n
```{r prediction}
fgc(est,n=253288,nsim=10)
```



## Fit the 3-component model
### Fitting to the model

Note the startingpic value can be specifided at a list of values, and then the one with the largest log-likelihood is selected as the final model. 

```{r 3-component model}

# starting value of 3-component model comes from 2-component model estimates. 
starting <- rep(0,5)
starting[1] <- est[1]
starting[2] <- 1/9
starting[3] <- est[2]*5
starting[4] <- starting[3]/10
starting[5] <- est[3]

fit3 <- genesis(df,modelcomponents=3, cores=24, LDcutoff=1, c0=10,starting=starting,filter=F)
fit3$estimates

est <- fit3$estimates$`parameter (pic, p1, sigmasq1, sigmasq2, a) estimates` # the model parameter estimates
v <- fit3$estimates$`Covariance matrix of parameter estimates` # the covariance matrix of model parameter estimtaes

# est <- c(9.583307e-03,8.562964e-02,1.487684e-04,2.086576e-05,1.498790e-06)
# v <- matrix(c(1.346867e-06, -1.036601e-05,1.303674e-09, -2.444779e-09, -8.840004e-11, -8.087098e-06, 2.721732e-04,-1.428814e-07,  1.519391e-09,  7.685736e-10, -1.188552e-09, -1.270322e-07,  1.637225e-10, 8.813366e-12,-2.151989e-14, -2.584200e-09,7.023163e-09 , 3.743461e-12,  7.008698e-12,  7.664961e-14, -8.582952e-11,8.923473e-10, -1.744992e-13,  6.291208e-14,  1.402886e-14),nrow=5,ncol=5)
```

### Get the density plot for the susceptibility SNPs 
```{r density plot}
x_seq = seq(-0.02,0.02,length.out = 1000); 
y_seq = apply(matrix(x_seq,ncol=1),1,function(t) dcausal(t,est))
plot(x_seq, y_seq,type="l",ylim=c(0,250),xlab="Joint effect size", ylab="Probability Density")
```

### Make future projections with specified sample size n
```{r future projections}
projection(est,n=253288)
```

### Predict number of discoveries in future GWAS and its confidence interval
```{r future projections}
fdis(est,v,n=253288)
```

### Calculate number of SNPs falling in an interval
```{r number of SNPs in an interval}
numInterval(0.005,Inf,est)
```

### Predict genomic control factor in future GWAS with specified sample size n
```{r prediction}
fgc(est,n=253288,nsim=10)
```

