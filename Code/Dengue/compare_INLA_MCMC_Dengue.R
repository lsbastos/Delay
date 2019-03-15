# Reporting delay modelling
# Compare R-INLA and MCMC for the dengue data


# Reading dengue data and load required liraries
load("denguesinan.RData") 
require(INLA)
require(nimble)
require(coda)

# This function creates a data frame from the run-off triangle matrix to be used in INLA
make.df.trian <- function(M){
  Time <- nrow(M)
  Delay <- ncol(M)
  aux.df <- data.frame(Y = as.vector(as.matrix(M)), 
                       Time = rep(x = 1:Time, times = Delay),
                       Delay = rep(x = 0:(Delay-1), each=Time)
  )
  aux.df
}

# Calculating delay time
d$DelayDays <- difftime(as.Date(as.character(d$DT_DIGITA),format="%Y-%m-%d"), 
                        as.Date(as.character(d$DT_NOTIFIC),format="%Y-%m-%d"),
                        units = "days")

# Number of notifications greater than 6 months (>= 182 days)
sum(d$DelayDays > 183, na.rm = T) 

d <- na.exclude(d[d$DelayDays < 183, ])

sum(d$DelayDays > 183, na.rm = T) 

# Excluding 2010, 2013, 2014
sum(d$NU_ANO != 2010 & d$NU_ANO != 2013 & d$NU_ANO != 2014) 

d <- d[d$NU_ANO != 2010 & d$NU_ANO != 2013 & d$NU_ANO != 2014, ]
d$SEM_NOT <- factor(d$SEM_NOT)

# Delay in weeks
d$DelayWeeks <- floor( d$DelayDays / 7 )

Dmax <- max(10, as.numeric(quantile(d$DelayWeeks, probs = 0.90)) )

aux <- tapply(d$DelayWeeks >= 0 , INDEX = d$SEM_NOT, FUN = sum, na.rm = T)
delay.tbl <- data.frame(Notifications = aux[order(rownames(aux))])
for(k in 0:26){  
  aux <- tapply(d$DelayWeeks == k, INDEX = d$SEM_NOT, FUN = sum, na.rm = T)
  delay.tbl[paste("d",k, sep="")] <- aux[order(rownames(aux))]
}
delay.week <- paste("d",0:26, sep="")

# Data frame to calculate the actual observed data
delay.data <- delay.tbl[delay.week]


# 15th epidemic week of 2012 
Tactual <- 68


### The INLA model as per paper
half_normal_sd =function(sigma){
  return(
    paste("expression:
              sigma = ",sigma,";
              precision = exp(log_precision);
              logdens = -0.5*log(2*pi*sigma^2)-1.5*log_precision-1/(2*precision*sigma^2);
              log_jacobian = log_precision;
              return(logdens+log_jacobian);",sep='')
  )
} 

model <- Y ~ 1 +
  f(d, model = "rw1", constr=T, hyper = list("prec" = list(prior = half_normal_sd(1)))) +
  f(w, model = "rw2", hyper = list("prec" = list(prior = half_normal_sd(1))),cyclic=TRUE) +
  f(t, model = "rw1", constr= T, hyper = list("prec" = list(prior = half_normal_sd(0.1)))) +
  f(tt, model = "rw1", constr= T, hyper = list("prec" = list(prior = half_normal_sd(0.1))),replicate=d) 

  # Observed data
  delay.data.obs <- delay.data[1:Tactual,(0:Dmax)+1]
  
  # Time index of the unknown counts (Dmax+1,...,Tactual) 
  index.time <- (Tactual-Dmax+1):Tactual

  delay.data.obs.trian <- delay.data.obs
  
  # Creating the run-off triangle data frame
  delay.data.obs.trian[outer(1:Tactual, 0:Dmax, FUN = "+") > Tactual] <- NA

  # Creating a data frame for INLA
  delay.inla.trian <- make.df.trian(delay.data.obs.trian)
  # Find the missing values
  index.missing <- which(is.na(delay.inla.trian$Y))
    
  inla.data <- list(Y=delay.inla.trian$Y,t=delay.inla.trian$Time,tt=delay.inla.trian$Time,
                 d=delay.inla.trian$Delay+1)
  inla.data$w <- inla.data$t%%52
  
  # Fit the INLA model. Note this version takes very long (hours) since we ask for more values of the hyperparameters
  # to be sampled. This is to ensure full comparability with the MCMC.
  output <- inla(model, family = "nbinomial", data = inla.data,
                 control.predictor = list(link = 1, compute = T),
                 control.compute = list( config = T, waic=TRUE, dic=TRUE),
                 control.family = list( 
                   hyper = list("theta" = list(prior = "loggamma", param = c(1, 0.1)))
                 ), # theta prior changed to be more sensible
                 control.inla = list(int.strategy='grid', dz=0.7) # to get better joint posterior estimates
  )
  #save(output,file="/home/theo/Dropbox/Data_archive/Delay/INLAmodel.RData")
  
  
  # Produce samples of the random effects and hyperparameters and then use them to predict the counts.
  n.sim <- 10000 # Number of posterior samples
  delay.samples.list <- inla.posterior.sample(n = n.sim, output) # Posterior samples list
 
  # Sampling the missing triangle from inla output in vector format
  vector.samples <- lapply(X = delay.samples.list, FUN = function(x, idx = index.missing) rnbinom(n = idx, mu = exp(x$latent[idx]), size = x$hyperpar[1])) 
  
  # Creating a vectorized version of the triangle matrix
  delay.vec.trian <- inla.matrix2vector(as.matrix(delay.data.obs.trian[index.time,]))
  
  # Transforming back from the vector form to the matrix form
  matrix.samples <- lapply(vector.samples, FUN = function(xxx, data = delay.vec.trian){
    data[which(is.na(data))] <- xxx
    inla.vector2matrix(data, ncol = Dmax+1) } )
 
  # Samples of {N_t : t=Tactual-Dmax+1,...Tactual}
  predictions.samples <- sapply(matrix.samples, FUN = function(x) rowSums(x) )

  # Samples of n_{t,d}
  inla.preds <- matrix(nrow=ncol(predictions.samples),ncol=55) 
  for(i in 1:ncol(predictions.samples)){ inla.preds[i,] <- vector.samples[[i]] } # INLA samples of run-off triangle





## Now, fit the model using MCMC in Nimble
####################################################################
y <- delay.data.obs.trian
nT <- nrow(y)
nD <- ncol(y)
# second order CAR for seasonality (52 weeks)
nW <- 52
adj_w=c(51,52,2,3,52,1,3,4,as.numeric(rbind(1:48,2:49,4:51,5:52)),49,50,52,1,50,51,1,2)
num_w=rep(4,52)
weights_w=rep(c(-1,4,4,-1),52)
# Now put the data in different format
Consts <- list(nT=nT,nD=nD,nW=nW,adj_w=adj_w,num_w=num_w,weights_w=weights_w,LadjW=length(adj_w))
nimbleData <- list(y=as.matrix(y),w=c(1:52,1:16))
# nimble model 
nbModel <- nimbleCode({ 
  for(t in 1:nT){
    for(d in 1:nD){
      y[t,d] ~ dnegbin( p[t,d] , theta )
      log(lambda[t,d]) <- alpha[t] + gamma[t,d] + eta[w[t]] + beta[d]
      p[t,d] <- theta/(theta+lambda[t,d])
    }}
  alpha[1] <- 0 
  beta[1] <- 0
  for(i in 2:nD){
    beta[i] ~ dnorm(beta[i-1],tau_b)
  }
  for(i in 2:nT){
    alpha[i] ~ dnorm(alpha[i-1],tau_a)
  }
  for(d in 1:nD){
    gamma[1,d] <- 0 
    for(j in 2:nT){
      gamma[j,d] ~ dnorm(gamma[j-1,d],tau_g)
    }
  }
  eta[1:nW] ~ dcar_normal(adj=adj_w[1:LadjW],num=num_w[1:nW],weights=weights_w[1:LadjW],c=2,tau=tau_e,zero_mean=0)
  ## Priors
  tau_a <- 1/(sigA*sigA)
  tau_b <- 1/(sigB*sigB)
  tau_g <- 1/(sigG*sigG)
  tau_e <- 1/(sigG*sigE)
  theta ~ dexp(0.1)
  sigA ~ T(dnorm(0, 0.1), 0, )
  sigB ~ T(dnorm(0, 1), 0, )
  sigG ~ T(dnorm(0, 0.1), 0, )
  sigE ~ T(dnorm(0, 1), 0, )
})
########## 3 parallel chains
inits <- nimbleNBModel <- MCMCconfig <- nbMCMC <- Cnb <- CnbMCMC <- list()
seed <- c(29,22,26)
for(i in 1:3){
  inits[[i]] <- list(sigA=1,sigB=1,theta=2, sigG=1,sigE=1, alpha=c(NA,rnorm(nT-1)), beta=c(NA,rnorm(Consts$nD-1)),eta=rnorm(nW),
                gamma=rbind(rep(NA,nD),matrix(rnorm(nD*(nT-1)),ncol=nD,nrow=nT-1)), y=as.matrix(delay.data.obs) )
  nimbleNBModel[[i]] <- nimbleModel(code=nbModel, constants=Consts, data=nimbleData, inits = inits[[i]])
  MCMCconfig[[i]] <- configureMCMC(nimbleNBModel[[i]],monitors=c("theta","sigA","sigB","sigG","sigE","alpha","beta","gamma","eta","y")) # Default samplers
  # Compile
  nbMCMC[[i]] <- buildMCMC(MCMCconfig[[i]])
  Cnb[[i]] <- compileNimble(nimbleNBModel[[i]],resetFunctions = TRUE)
  CnbMCMC[[i]] <- compileNimble(nbMCMC[[i]], project = nimbleNBModel[[i]],resetFunctions = TRUE)
}
library(doParallel)
registerDoParallel(cores=3)
full_samples <- foreach(i=1:3)%dopar%{
	samples <- runMCMC(CnbMCMC[[i]],inits=inits[[i]],setSeed=seed[i],nchains = 1, nburnin=2e+06,niter = 2.5e+06,samplesAsCodaMCMC = TRUE, summary = FALSE, WAIC = FALSE, thin=10)
  return(samples)
}
#save(full_samples,file="/home/theo/Dropbox/Data_archive/Delay/MCMCresults.RData")
##########
results <- as.mcmc.list(full_samples)
plot(results,ask=T)
plot(results[,c("theta","sigA","sigB","sigG","sigE")],ask=T)
effectiveSize(results[,c("theta","sigA","sigB","sigE","sigB")])
libraryEDISON)
dum <- results[,1:884]
dum2 <- list()
for(i in 1:3){ dum2[[i]] <- t(as.matrix(dum[[i]]) ) }
Rhat <- psrf(dum2)

### Simulate from posterior pred of all n_{t,d}
sims.matrix <- do.call(rbind.data.frame, results)
n.sims <- nrow(sims.matrix)
y_pred_array <- array(dim=c(n.sims,Tactual,Dmax+1))  
theta <- sims.matrix[,"theta"]
for(i in 1:nT){
  print(paste(i," out of ",nT,sep=""))
  for(j in 1:nD){
    w_index <- nimbleData$w[i]
    a <- sims.matrix[,paste("alpha[",i,"]",sep="")]
    b <- sims.matrix[,paste("beta[",j,"]",sep="")]
    e <- sims.matrix[,paste("eta[",w_index,"]",sep="")]
    g <- sims.matrix[,paste("gamma[",i,", ",j,"]",sep="")]
    Mean <- exp(a + b + g + e)
    y_pred_array[,i,j] <- rnbinom(n.sims,size=theta,mu=Mean)
  }
}
# Aggregate to get samples of the totals
N_samples <- rowSums(y_pred_array,dims=2)




######### Compare predictions of totals with INLA preditions
week_index <- 59:68
delay_index <- list()
for(i in 1:10){
  delay_index[[i]] <- (12-i):11
}
preds <- list()
for(i in 1:length(week_index)){
  n <- length(delay_index[[i]])
  y_names <- paste0(rep("y[",n),rep(week_index[i],n),rep(", ",n),delay_index[[i]],rep("]",n),sep="")
  preds[[i]] <- as.matrix(sims.matrix[,y_names])
}
N <- lapply(preds,rowSums)

## QQ plots to compare INLA and nimble predictions for paper
x11(width=15,height=10)
par(mfrow=c(3,4),mar = c(2.2,2.2,1,1) + 0.1,cex=1.1,oma = c(2, 2, 0, 0))
for(i in 1:length(week_index)){
qqplot(predictions.samples[i,],N[[i]]+sum(y[week_index[i],],na.rm=T),pch=20,main=bquote(N[.(week_index[i])]),ylab="",xlab="")
abline(0,1)
mtext("MCMC", outer = TRUE, cex = 2, side=2)
mtext("R-INLA", outer = TRUE, cex = 2, side=1)
}
for(i in (length(week_index)-1):length(week_index)){
qqplot(inla.preds[,45+i],preds[[i]][,ncol(preds[[i]])],pch=20,ylab="",xlab="",main=bquote(n[.(week_index[i])][","][.(i)]))
abline(0,1)
}
dev.print(pdf,"comparisonsN2.pdf")



