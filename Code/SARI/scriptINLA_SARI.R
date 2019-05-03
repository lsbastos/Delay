# Script for SARI data
# Leo Bastos
# Fri Mar 15 08:44:24 2019

####################################################
# Required packages:
####################################################

library(tidyverse)
library(lubridate)

# Packages for spatial analysis
library(tmap)
library(rgdal)
library(spdep)
library(maptools)

# Package for modelling
library(INLA)

####################################################
# Reading data
####################################################


# Reading notification data for Parana state (PR)
sari <- read.csv("Data/clean_data_srag_epiweek_delay_table_PR.csv")

# Selecting the variables of interest
names(sari)

sari <- sari[,1:15]


episem2date <- function(epiyear, epiweek, weekday=0){
  # # Separate year and week:
  # epiyear <- stringi::stri_sub(epiyearweek, 1, 4)
  # epiweek <- as.integer(stringi::stri_sub(epiyearweek, -2, -1))
  
  # Obtain sunday of first epiweek of epiyear
  day.one <- as.Date(paste(epiyear, '01', '01', sep='-'))
  day.one.week <- as.numeric(format(day.one,"%w")) # 0: sunday
  
  # Check wether week day of Jan 1st was before or after a Wednesday
  # and set the start of the epiyear accordingly
  first.epiweek.day <- ifelse (day.one.week <=3, day.one - day.one.week ,
                               day.one + (7 - day.one.week) )
  first.epiweek.day <- as.Date(first.epiweek.day,
                               origin = '1970-01-01')
  return(first.epiweek.day + 7*(epiweek-1) + weekday)
}

# Epiweek to data
temp <- episem2date(epiyear = sari$epiyear, epiweek = sari$epiweek)

sari$EpiWeekStart <- temp






####################################################
# Maps and INLA stuff for spatial effects
####################################################


# Municipality level shapefiles
PR <- readOGR(layer = "41REGSAUD", dsn = "Data/PR_regsaud/")

# Health region spatial aggregation
PR.union <- unionSpatialPolygons(PR, PR$CO_REGSAUD)

tm_shape(PR) + tm_borders(lwd = 0.5) + tm_shape(PR.union) + tm_borders(lwd = 2)

# Creating the graph file for INLA 
graph.file <- "PRunion.graph"
PR.union %>% poly2nb() %>% nb2INLA(file = graph.file )


####################################################
# Preparing the data
####################################################


# Taking part of the dataset (not = notified cases)
# Making it as a long dataframe. 
# sari.not <- sari %>% filter( (epiyear >= 2012 & epiyear < 2014 & regionalsaude != 0 ) ) %>%
sari.not <- sari %>% filter( (epiyear >= 2016 & epiyear <= 2017 & regionalsaude != 0 ) ) %>%
  gather( key = Delay, value = Y, d0:d10) %>% 
  mutate(
    Delay = parse_number(Delay)
    )

# Maximum allowed delay (already set in gather function, d0:d10)
Dmax <- 10


# We will set today as the 26th epiweek  of 2013
#Today <- episem2date(2013,26)
Today <- episem2date(2017,14)

# Excluding the observations after "Today"
sari.obs <- sari.not %>% filter(EpiWeekStart <= Today)

# Setting as NA the occur-but-not-yet-reported cases

for( d in Dmax:1){
  sari.obs$Y[ (sari.obs$EpiWeekStart > (Today - 7*d)) & (sari.obs$Delay == d) ] <- NA
}

# Checking! it seems fine!
# sari.obs %>% filter( is.na(Y) == T) %>% group_by(Delay) %>% summarise(sum(is.na(Y)))

# Creating a Time variable. 
FirstDate <- min(sari.obs$EpiWeekStart)

# Creating the Time variable for INLA. First week is set as 1.
sari.obs <- sari.obs %>% mutate( Time = as.integer(difftime(EpiWeekStart, FirstDate, unit="weeks")) + 1 )





####################################################
# State level data - Needed for see the time series
####################################################


# Aggregated notified data
sari.state <- sari.not %>% filter(EpiWeekStart < Today + 1)  %>% group_by(EpiWeekStart) %>% summarise( Y = sum(Y)) 

# Time series
p0 <- ggplot(sari.state, aes(x = EpiWeekStart, y = Y, colour = "Actual number of cases")) + geom_line()

# Aggregated observed data at Today!
sari.state.obs <- sari.obs %>% filter( (EpiWeekStart >= Today - 7*Dmax)  ) %>%  
                  group_by(EpiWeekStart) %>% summarise( Y = sum(Y, na.rm = T) )


p0.s <- p0 + geom_line(data = sari.state.obs, aes(x = EpiWeekStart, y = Y, colour = "Reported cases with delay")) +
  ylab("Total number of SARI cases") + xlab("Time") + 
  scale_x_date(date_breaks = "3 months", date_labels = "%Y-%w") + #, limits = c(Today - 90, Today)) + ylim(c(0,100)) +
  scale_colour_manual(name = "", values = c("black", "red")) + 
  scale_linetype_manual(name = "", values = c("solid", "dashed")) + theme_bw( base_size = 14) + 
  theme( legend.position=c(0.8,.9) )
p0.s

# ggsave(plot = p0.s, filename = "sariPR.png", device = "png", dpi = 600)


####################################################
# Modelling using INLA
####################################################



# Find the indexes for the missing values 
index.missing <- which(is.na(sari.obs$Y))


#######################
# Building the models
#######################

# Half normal prior
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


# Adding an health region code form 1:R. (too simple here, need to generalise)
sari.obs <- sari.obs %>% mutate( GEOCOD = regionalsaude - 41000, 
                                 Delay.Region = paste(Delay,regionalsaude, sep="."),
                                 Time.Region = paste(Time,regionalsaude, sep=".") )

sari.obs <- mutate(sari.obs, 
                   Time2 = Time,
                   Delay2 = Delay + 1
                   )


# Model equation: space-delay effect
model.sd <- Y ~ 1 + 
  f(Time, model = "rw1", constr = T,
    hyper = list("prec" = list(prior = half_normal_sd(.1) ))
    ) + 
  f(Delay, model = "rw1", constr = T,
    hyper = list("prec" = list(prior = half_normal_sd(.1) ))
    ) +
  f(GEOCOD, model = "iid", constr = T,
    hyper = list("prec" = list(prior = half_normal_sd(.1) ))
  ) +
  f(Delay.Region, model = "iid", constr = T, 
    hyper = list("prec" = list(prior = half_normal_sd(.1))) 
  )

# Model equation: Time-delay effect
model.sd.b <- Y ~ 1 + 
  f(Time, model = "rw1", constr = T,
    hyper = list("prec" = list(prior = half_normal_sd(.1) ))
  ) + 
  f(Delay, model = "rw1", constr = T,
    hyper = list("prec" = list(prior = half_normal_sd(.1) ))
  ) +
  f(Time2, model = "rw1", constr= T, replicate = Delay2,
    hyper = list("prec" = list(prior = half_normal_sd(.1)))
  ) +
  f(GEOCOD, model = "iid", constr = T,
    hyper = list("prec" =list(prior = half_normal_sd(.1) ))
  ) +
  f(Delay.Region, model = "iid", constr = T, 
    hyper = list("prec" = list(prior = half_normal_sd(.1))) 
  )


# Model equation: Time-delay and spatio-delay effects
model.sd.b.bym <- Y ~ 1 + 
  f(Time, model = "rw1", constr = T,
    hyper = list("prec" = list(prior = half_normal_sd(.1) ))
  ) + 
  f(Delay, model = "rw1", constr = T,
    hyper = list("prec" = list(prior = half_normal_sd(1) ))
  ) +
  f(Time2, model = "rw1", constr= T, 
    hyper = list("prec" = list(prior = half_normal_sd(.1))), 
    replicate = Delay2
  ) +
  f(GEOCOD, model = "bym" , graph = graph.file, constr = T,
    hyper = list("prec.unstruct" = list(prior = half_normal_sd(.1) ),
                 "prec.spatial" = list(prior = half_normal_sd(.1) ))
  ) +
  f(Delay.Region, model = "iid", constr = T, 
    hyper = list("prec" = list(prior = half_normal_sd(.1) )) 
  )







#######################
# Running the model
#######################




# Running the Negative Binomial model in INLA (space-delay)
output <- inla(model.sd.b.bym, family = "nbinomial", data = sari.obs, num.threads = 4,
               control.predictor = list(link = 1, compute = T),
               control.compute = list( config = T),
               control.family = list( 
                 hyper = list("theta" = list(prior = "loggamma", param = c(1, 0.1)))
                 )
)





################
# Estimates
################

# Fixed effects 
output$summary.fixed

# Hyperparameters (negative binomial parameter, random effects precisions)
output$summary.hyperpar




# Algorithm to get samples for the predictive distribution for the number of cases

# Step 1: Sampling from the approximate posterior distribution using INLA
sari.samples.list <- inla.posterior.sample(n = 1000, output)


# Step 2: Sampling the missing triangle (in vector form) from the likelihood using INLA estimates
vector.samples <- lapply(X = sari.samples.list, 
                         FUN = function(x, idx = index.missing) rnbinom(n = idx, mu = exp(x$latent[idx]), size = x$hyperpar[1])) 


# Step 3: Calculate N_t for each triangle sample {N_t : t=Tactual-Dmax+1,...Tactual}

tibble.samples <- lapply( X = vector.samples,
                          FUN = function(x){
                            data.aux <- sari.obs
                            data.aux$Y[index.missing] <- x
                            data.aggregated <- data.aux %>% filter( EpiWeekStart >= Today - 70) %>% group_by(EpiWeekStart) %>% summarise( Y = sum(Y) )
                            data.aggregated
                          }
                        )

# Nowcasting
sari.state.pred <- bind_rows(tibble.samples)


# Epidemic threshold (given by the Ministry of Health influenza team)
# The threshold was calculated using MEM
Limiar <- 64.6

sari.state.pred <- sari.state.pred %>% group_by(EpiWeekStart) %>% 
  summarise( Median = median(Y), 
             LI = quantile(Y, probs = 0.025),
             LS = quantile(Y, probs = 0.975),
             Q75 = quantile(Y, probs = 0.75),
             Prob = mean(Y > Limiar))

sari.state.pred

# Time series
p0 <- ggplot(sari.state, 
             #aes(x = EpiWeekStart, y = Y, colour = "Actual number of cases", linetype = "Actual number of cases")) + 
             aes(x = EpiWeekStart, y = Y, colour = "Eventually reported cases", linetype = "Eventually reported cases")) + 
  geom_ribbon( data = sari.state.pred, aes( y = Median, ymin=LI, ymax=LS), fill = 'gray', color = 'gray', show.legend = F) +
  geom_line( size = 1)


p1.s <- p0 + 
  geom_line(data = sari.state.obs, aes(x = EpiWeekStart, y = Y, colour = "Currently reported cases", linetype = "Currently reported cases"), size = 1) +
  geom_line(data = sari.state.pred, aes(x = EpiWeekStart, y = Median, colour = "Model predictions", linetype = "Model predictions"), size = 1) +   
  geom_hline(yintercept = 64.6, color = 'blue', linetype = 'dashed', size = 1) +
  scale_colour_manual(name = "", values = c("red", "black", "black"), guide = guide_legend(reverse=F)) + 
  scale_linetype_manual(name = "", values = c("dashed", "solid", "dotted"), guide = guide_legend(reverse=F)) +
  ylab("Number of SARI cases") + xlab("Time (Week/year)") 

p1.sA <- p1.s +
  scale_x_date(date_breaks = "3 months", date_labels = "%W/%Y") +
  #scale_fill_manual(name = "", values = c("white","white","white")) + 
  theme_bw(base_size = 14) + theme( legend.position=c(0.75,.85), legend.key.width = unit(2.5,"line") ) +
  annotate("text",x=as.Date("2016-05-15"), y=70, size=4, label='Epidemic threshold')
p1.sA #+ guides(linetype = guide_legend(override.aes = list(fill = NA)))  

# ggsave(filename = "sariPR2a.png", plot = p1.sA, device = "png", dpi = 600)


# Zooming in
p1.sB <- p1.s %+% 
  subset(sari.state, EpiWeekStart > "2017-01-01") +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%W/%Y") + 
  theme_bw(base_size = 14) + theme( legend.position=c(0.3,.87), legend.key.width = unit(2.5,"line") ) +
  annotate("text",x=as.Date("2017-01-15"), y=63, size=4, label='Epidemic threshold')
  
p1.sB
# ggsave(filename = "sariPR2b.png", plot = p1.sB, device = "png", dpi = 600)



######################
# Time random effects
######################

# Setting dates for plotting the Time random effects
output$summary.random$Time <- output$summary.random$Time %>% mutate( Date = FirstDate + (ID - 1)*7 )

pTime <- ggplot( data = output$summary.random$Time, 
                 mapping = aes( x = Date, y = `0.5quant`, ymin = `0.025quant`, ymax = `0.975quant`)) 

pTime.s <- pTime + geom_line() + geom_ribbon(color = "gray", alpha = 0.5) + 
  theme_bw(base_size = 14) + xlab("Time (Week/year)") + ylab("Time random effects") +
  scale_x_date(date_breaks = "2 months", date_labels = "%W/%Y") 
  
pTime.s

# ggsave(filename = "pTime.png", plot = pTime.s, device = "png", dpi = 600)

######################
# Delay random effects
######################

# Setting dates for plotting the Time random effects

pDelay <- ggplot( data = output$summary.random$Delay, 
                  mapping = aes( x = ID, y = `0.5quant`, ymin = `0.025quant`, ymax = `0.975quant`)) 

pDelay.s <- pDelay + geom_line() + geom_ribbon(color = "gray", alpha = 0.5) + theme_bw(base_size = 14) +
  scale_x_continuous(name = "Delay (in weeks)", breaks = 0:10) + ylab("Delay random effects")

pDelay.s

# ggsave(filename = "pDelay.png", plot = pDelay.s, device = "png", dpi = 600)


######################
# Time-Delay effects
######################

# Setting dates for plotting the Time random effects


output$summary.random$Time2 <- output$summary.random$Time2 %>% 
  mutate( Date = FirstDate + (ID - 1)*7,
          Delay = rep(0:Dmax, each = 66))

pTimeDelay <- ggplot( data = output$summary.random$Time2, 
                  mapping = aes( x = Date, y = `0.5quant`, ymin = `0.025quant`, ymax = `0.975quant`))

pTimeDelay.s <- pTimeDelay + geom_line()  + theme_bw(base_size = 14) + geom_ribbon(color = "gray", alpha = 0.5) +
  scale_x_date(date_breaks = "5 months", date_labels = "%W/%Y") + ylab("Time-Delay random effects")

pTimeDelay.s <- pTimeDelay.s + facet_wrap( ~Delay )
pTimeDelay.s

# ggsave(filename = "pTimeDelay.png", plot = pTimeDelay.s, device = "png", dpi = 600)


#############################
# Space-Delay random effects
#############################


# Setting dates for plotting the Time random effects

temp <- output$summary.random$Delay.Region  %>% 
  separate(col = "ID", into = c("Delay", "Region"), sep = "\\.", remove = F) %>%
  mutate( Delay = as.numeric(Delay)) %>%
  left_join( y = output$summary.random$Delay[,c("ID", "mean")], by = c("Delay" = "ID") ) %>%
  mutate( Mean = mean.y + mean.x)

pSpaceDelay <- ggplot( data = temp, 
                       mapping = aes( x = Delay, y = Mean, color = Region)) 

pSpaceDelay.s <- pSpaceDelay + geom_line() + 
  scale_x_continuous(name = "Delay (in weeks)", breaks = 0:10) + 
  scale_y_continuous(name = "Posterior mean of delay + space-delay random effects") + 
  scale_color_discrete("Health region code") + 
  theme_bw(base_size = 14) + theme(legend.position = c(0.8, .65)) 


pSpaceDelay.s 

# ggsave(filename = "pSpaceDelay.png", plot =  pSpaceDelay.s, device = "png", dpi = 600)


######################
# Spatial random effects
######################

# Setting dates for plotting the Time random effects

pRegion <- ggplot( data = output$summary.random$GEOCOD[1:22,],
                  mapping = aes( x = factor(ID), y = `0.5quant`, ymin = `0.025quant`, ymax = `0.975quant`))

pRegion + geom_point() + geom_linerange()


# Mapping


# Reading Parana health region shapefile
#PR <- readOGR(layer = "41REGSAUD", dsn = "../Data/SARI/PR_regsaud/")



output$summary.random$GEOCOD <- output$summary.random$GEOCOD %>% mutate( GEOCOD = ifelse(ID > 22, ID-22, ID) + 41000)



PR@data <- PR@data %>% mutate(
  CO_REGSAUD = as.numeric(as.character(CO_REGSAUD))
) %>%
  left_join( y = output$summary.random$GEOCOD %>% 
                                    group_by(GEOCOD) %>% summarise(Mean = sum(mean)), 
                                  by = c("CO_REGSAUD" = "GEOCOD") )




PR.coords <- coordinates(PR)


PR.union <- unionSpatialPolygons(SpP = PR, IDs = PR$CO_REGSAUD )

PR.union$ID <- sapply(slot(PR.union, "polygons"), function(x) slot(x, "ID")) 

pMap.s <- tm_shape(PR) + 
  tm_fill(col = "Mean", style = "fixed", title = " ", midpoint = NA,
          breaks = c(-1.7, -1 ,-0.5, 0.5 , 1 , 3.6) 
          #, labels = c("< -1.0", "-1.0 to -0.5", "-0.5 to 0.5", "0.5 to 1.0", "> 1.0")
          ) +   
  tm_layout(frame=F, legend.position = c("right", "top"), legend.text.size = 1.25) + tm_shape(PR.union) + 
  tm_text( "ID") + 
  tm_borders() 


pMap.s

# tmap_save(filename = "pMap_label.png", tm = pMap.s, dpi = 600)



#####################################################
# Variance of spatial random effects accross regions


varRatio.samples <- sapply( X = sari.samples.list,
                            FUN = function(x){
                              spatialRE.idx <- which(startsWith(rownames(x$latent), prefix = "GEOCOD"))
                              auxBoth <- x$latent[spatialRE.idx][1:22]
                              auxUns <- x$latent[spatialRE.idx][23:44]
                              varSpatial <- var(auxBoth - auxUns) 
                              varBoth <- var(auxBoth)
                              ratio <- varSpatial / varBoth
                              c(varSpatial = varSpatial, varBoth = varBoth, Ratio = ratio )
                            } 
)

varRatio.tbl <- data.frame(t(varRatio.samples))

c( mean = mean(varRatio.tbl$Ratio), quantile(varRatio.tbl$Ratio, probs = c(0.5, 0.025, 0.975)) )


ppp <- ggplot( varRatio.tbl, aes(x = Ratio)) + geom_density() + 
  geom_density( fill = "red", alpha = .7) + theme_bw(base_size = 14) +
  xlab( expression(R == V(psi^{IAR}) / V(psi^{IAR} + psi^{ind}) )) +
  ylab("Ratio density") + xlim(c(0,1.2))

ppp

# ggsave(plot = ppp, filename = "densitityRatio.png", device = "png")



