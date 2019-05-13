# Simulating from the prior distribution

# Reading dengue data
load("denguesinan.RData") 

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

# Observed data
delay.data.obs <- delay.data[1:Tactual,(0:Dmax)+1]

# Time index of the unknown counts (Dmax+1,...,Tactual) 
index.time <- (Tactual-Dmax+1):Tactual

delay.data.obs.trian <- delay.data.obs

# Creating the run-off triangle data frame
#delay.data.obs.trian[outer(1:Tactual, 0:Dmax, FUN = "+") > Tactual] <- NA

# Creating a data frame for INLA
delay.inla.trian <- make.df.trian(delay.data.obs.trian)

# Find the missing values
#index.missing <- which(is.na(delay.inla.trian$Y))

inla.data <- data.frame(Y=delay.inla.trian$Y,t=delay.inla.trian$Time,tt=delay.inla.trian$Time,
                        d=delay.inla.trian$Delay+1)
inla.data$w <- inla.data$t%%52


############################################
# Sampling

# Sampling the hiperparameters

# Half-normal samplers

rhn <- function(n, sd = 1){
  U = runif(n, 0, 1)
  X = sd * qnorm((U+1)/2)
  return(X)
}



sample.notifications = function(){

  sigma.alpha = rhn(1, .1)
  sigma.beta = rhn(1, 1)
  sigma.gamma = rhn(1, .1)
  sigma.eta = rhn(1, 1)
  phi = rgamma(1, 1, 0.1)
  
  alpha.t = rep(0,Tactual)
  for(t in 2:Tactual)
    alpha.t[t] = rnorm(n = 1, mean = alpha.t[t-1], sd = sigma.alpha)
  alpha.t = alpha.t - mean(alpha.t)
  
  
  beta.d = rep(0,Dmax+1)
  for(d in 2:(Dmax+1))
    beta.d[d] = rnorm(n = 1, mean = beta.d[d-1], sd = sigma.beta)
  beta.d = beta.d - mean(beta.d)
  
  gamma.td = rep(0, (Dmax+1)*Tactual )
  for(d in 2:(Dmax+1)){
    for(t in 2:Tactual)
      gamma.td[t + Tactual *(d-1)] = 
        rnorm(n = 1, 
              mean = gamma.td[t - 1 + Tactual *(d-1)], 
              sd = sigma.gamma)
  }
  gamma.td = gamma.td - mean(gamma.td)
  
  eta.w = rep(0, 52)
  for( w in 3:52)
    eta.w[w] = rnorm(1, 2*eta.w[w-1] - eta.w[w-2], sd=sigma.eta)
  eta.w = eta.w - mean(eta.w)
  
  mu = rnorm(1, mean = 0, sd = 10)
  
  loglambda = mu + 
    alpha.t[inla.data$t] + 
    beta.d[inla.data$d] + gamma.td + eta.w[inla.data$w+1]
  
  yy = rnbinom(n = length(loglambda), mu = exp(loglambda), size = phi)
  yy  
}


sample.tbl = data.frame(
  Time = inla.data$t,
  matrix(NA, 
         nrow = nrow(inla.data),
         ncol = 1000)
  )

for(k in 1:1000)
  sample.tbl[,k+1] = sample.notifications()

library(tidyverse)

aaa = sample.tbl %>% 
  bind_cols( Y = inla.data[,1]) %>% 
  group_by(Time) %>%
  summarise_all(sum, na.rm = T)

plot(aaa$Y, type = "l", ylim = c(0,8000), 
     ylab = "Notifications",
     xlab = "Time")
for(i in 1:1000)
  lines(aaa[,i+1], col="lightgray")
lines(aaa$Y)

save(list = c("aaa"), file = "sim.RData")
