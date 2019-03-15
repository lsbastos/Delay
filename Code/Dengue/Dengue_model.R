# Reporting delay modelling of Dengue

# Reading dengue data
load("denguesinan.RData") 
require(INLA)

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

# Function for plotting INLA random effects
plot.inla.re = function(outputRE, x = outputRE$ID){
  plot( x, y = outputRE$mean, type = "n", ylim = range(outputRE[,c(4,6)]), ylab="", xlab="" )
  polygon(x = c(x, rev(x)),
          y = c(outputRE$'0.025quant', rev(outputRE$'0.975quant')),
          border = "black", col = "gray")
  lines(x, outputRE$mean, lty=1, lwd=2)
  lines(x = range(x), y = rep(0,2), lty=2)
}

# Calculating delay time
d$DelayDays <- difftime(as.Date(as.character(d$DT_DIGITA),format="%Y-%m-%d"), 
                        as.Date(as.character(d$DT_NOTIFIC),format="%Y-%m-%d"),
                        units = "days")


d <- na.exclude(d[d$DelayDays < 183, ])

# Excluding 2010, 2013, 2014
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
half_normal_sd <- function(sigma){
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
  f(w, model = "rw2", constr=T, hyper = list("prec" = list(prior = half_normal_sd(1))),cyclic=TRUE) +
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
  
  output <- inla(model, family = "nbinomial", data = inla.data,
                 control.predictor = list(link = 1, compute = T),
                 control.compute = list( config = T, waic=TRUE, dic=TRUE),
                 control.family = list( 
                   hyper = list("theta" = list(prior = "loggamma", param = c(1, 0.1)))),
  )
  
  n.sim <- 10000 # Number of posterior samples
  delay.samples.list <- inla.posterior.sample(n = n.sim, output) # Posterior samples list
 
  # Sampling the missing triangule from inla output in vector format
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


####################
## Model Checking ##
####################
y.samples <- do.call('rbind',lapply(X = delay.samples.list, FUN = function(x) rnbinom(Tactual*(Dmax+1), mu = exp(x$latent), size = x$hyperpar[1]))) 

y.array <- array(y.samples,dim=c(n.sim,Tactual,Dmax+1))
  
total.samples <- rowSums(y.array,dims=2)

true.totals <- delay.tbl$Notifications[1:68]

# Fitted values plot
x11(width=9,height=6)
par(mar = c(4, 4, 1, 1),cex=1.2)
plot(apply(total.samples,2,mean),true.totals,xlab="Mean predicted totals",ylab='Observed Totals (N)',pch=20)
abline(a=0,b=1)
dev.print(pdf,"fitted.pdf")

# Prediction plot
plot(delay.tbl$Notifications[1:68], type = "n", axes=T, xlab="Time (weeks)", ylab="Number of dengue cases",ylim=c(0,18000))
lines(delay.tbl$Notifications[1:68], lwd=4)
lines(apply(total.samples,2,mean),col='blue',lwd=2)
lines(apply(total.samples,2,quantile,0.025),col='blue',lwd=2,lty=2)
lines(apply(total.samples,2,quantile,0.975),col='blue',lwd=2,lty=2)

## Posterior predictive checks

# Median, mean variance
x11(width=9,height=3)
par(mfrow=c(1,3),mar = c(4, 4, 1, 1),cex=1.2)
plot(density(apply(total.samples,1,median)),xlab='Sample Median',main='')
abline(v=median(true.totals))
title(paste('p =',mean(apply(total.samples,1,median)<=median(true.totals))))
plot(density(apply(total.samples,1,mean)),xlab='Sample Mean',main='')
abline(v=mean(true.totals)) 
title(paste('p =',mean(apply(total.samples,1,mean)<=mean(true.totals))))
plot(density(apply(total.samples,1,var)),xlab='Sample Variance',main='')
abline(v=var(true.totals))
title(paste('p =',mean(apply(total.samples,1,var)<=var(true.totals))))
dev.print(pdf,"MMV.pdf")

# Autocorrelation of log-totals
autocor=function(x,lag){
  l=length(x)
  cor(x[1:(l-lag)],x[(1+lag):l],use='pairwise.complete.obs')
}
x11(width=12,height=5)
par(mfrow=c(2,4),mar = c(4, 4, 1, 1),cex=1)
for(i in 1:8){
  plot(density(apply(log(total.samples),1,function(x)autocor(x,i))),xlab=paste('Sample Auto-correlation Lag-',i,sep=''),main='')
  abline(v=autocor(log(true.totals),i))
}
dev.print(pdf,"acf.pdf")


# Means of different delay levels
x11(width=8,height=6)
par(mfrow=c(4,3),mar=c(1,1,1,1))
for(i in 1:(Dmax+1)){
  plot(density(apply(y.array[,,i],1,mean)),xlab='Sample Mean',main='')
  abline(v=mean(delay.data.obs[1:Tactual,i]))
  title(paste('Delay',i-1,': p =',
              mean(apply(y.array[,,i],1,mean)<=mean(delay.data.obs[1:Tactual,i]))))
}

# Covariance of different delay levels
cov_matrix <- matrix(nrow=Dmax+1,ncol=Dmax+1)
for(j in 1:(Dmax+1)){
  for(i in 1:(Dmax+1)){
      cov_matrix[i,j] <- mean(apply(cbind(y.array[,,i],y.array[,,j]),1,
                                 function(x)cov(x[1:Tactual],x[(Tactual+1):(2*Tactual)]))<=
                             cov(delay.data.obs[1:Tactual,i],delay.data.obs[1:Tactual,j]))
  }
}
# Coverage of predictive intervals for covariances
mean(cov_matrix>=0.025 & cov_matrix<=0.975)

# Sorted values plot
x11(width=9,height=6)
par(mar = c(4, 4, 1, 1),cex=1.2)
sorted.samples <- t(apply(total.samples,1,sort))
par(mfrow=c(1,1))
plot( sort(true.totals), y = apply(sorted.samples,2,mean), type = "n", 
      ylim = c(min(apply(sorted.samples,2,quantile,0.025)),max(apply(sorted.samples,2,quantile,0.975))), 
      ylab="Predicted Totals (N)", xlab='Observed Totals (N)' )
polygon(x = c(sort(true.totals), rev(sort(true.totals))),
        y = c(apply(sorted.samples,2,quantile,0.025), rev(apply(sorted.samples,2,quantile,0.975))),
        border = "black", col = "gray")
lines(sort(true.totals), apply(sorted.samples,2,mean), lty=1, lwd=2)
abline(a=0,b=1,lwd=2)
title('Sorted Totals')
dev.print(pdf,"sorted.pdf")

## 15th epidemic week of 2012 figure for paper
x11(width=9,height=6)
par(mar = c(4, 4, 1, 1),cex=1.2)
plot(delay.tbl$Notifications[1:68], type = "n", axes=F, xlab="Time (in epidemic weeks)", ylab="Number of dengue cases",ylim=c(0,10000))
  polygon(
    x = c(index.time, rev(index.time)), 
    y = c(apply(predictions.samples,1,quantile,probs = 0.025), 
          rev(apply(predictions.samples,1,quantile,probs = 0.975))),
    col = "lightgray", border = NA, lty = 3
  )
  lines(index.time, rowMeans(predictions.samples), col="black", lwd=4, lty=3)
  lines(delay.tbl$Notifications[1:68], lwd=4)
  lines(index.time, rowSums(delay.data.obs.trian[index.time,],na.rm = T), lwd=4, col="red", lty=2)
  axis(2)
  axis(1, at = c(2, 54 , 106) ,labels = 2011:2013 )
  points(x=Tactual, y=-0, pch = 10, cex=2)

  legend("topleft", c("Eventually reported cases", "Currently reported cases", "Model predictions"), lty=c(1,2,3), col=c("black","red","black"), lwd=4)
dev.print(pdf,"15thWeek.pdf")
  

# Radnom effect plots
x11(width=9,height=3)
par(mfrow=c(1,3),mar = c(4, 4, 1, 1),cex=1)
plot.inla.re(output$summary.random$t)
title("Time")
plot.inla.re(output$summary.random$d)
title("Delay")
plot.inla.re(output$summary.random$w)
title("Seasonal Cycle")
dev.print(pdf,"RE.pdf")

# Time splines for each delay
par(mfrow=c(3,4))
for(d in 1:(Dmax+1)){
  plot.inla.re(output$summary.random$tt[((d-1)*Tactual+1):(d*Tactual),])
  title(paste('Delay',d-1))
}
round(output$summary.hyperpar,2)







####################################################################



## Produce stamp plot of 12 weeks of predictions
x11(width=15,height=10)
par(mfrow=c(3,4), mar = c(2.2,2.2,1,1) + 0.1,cex=1.1,oma = c(2, 2, 0, 0))

for(Tactual in 68:79){
  
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
  
  output <- inla(model, family = "nbinomial", data = inla.data,
                 control.predictor = list(link = 1, compute = T),
                 control.compute = list( config = T, waic=TRUE, dic=TRUE),
                 control.family = list( 
                   hyper = list("theta" = list(prior = "loggamma", param = c(1, 0.1)))),
  )

  delay.samples.list <- inla.posterior.sample(n = 1000, output)
  
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
  
  plot(54:86,delay.tbl$Notifications[54:86], type = "n", axes=F, xlab="", ylab="",ylim=c(0,10000),main=bquote(N[.(Tactual)]))
  polygon(
    x = c(index.time, rev(index.time)), 
    y = c(apply(predictions.samples,1,quantile,probs = 0.025), 
          rev(apply(predictions.samples,1,quantile,probs = 0.975))),
    col = "lightgray", border = NA, lty = 3
  )
  lines(54:86,delay.tbl$Notifications[54:86], lwd=4)
  lines(index.time, rowMeans(predictions.samples), col="black", lwd=4, lty=3)
  lines(index.time, rowSums(delay.data.obs.trian[index.time,],na.rm = T), lwd=4, col="red", lty=2)
  axis(2)
  axis(1, at = c(54, 62, 70, 78, 86) ,labels = c("Jan", "Mar", "May", "Jul", "Sep") )
  points(x=Tactual, y=-0, pch = 10, cex=2)
  
  #  legend("topleft", c("The actual number of cases", "Reported cases with delay", "Model predictions"), lty=c(1,2,2), col=1:3, lwd=4)
} 
mtext("Time (weeks)", outer = TRUE, cex = 2, side=1)
mtext("Number of dengue cases", outer = TRUE, cex = 2, side=2)
dev.print(pdf,"predictionsN.pdf")



