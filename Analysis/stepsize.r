# Script to explore the numerical consquences of changes the Euler integration
# step size.

library(rstan)
library(doParallel)
library(plyr)
library(magrittr)
library(lubridate)
library(ggplot2)
library(gridExtra)

# Loading and processing mosquito and case data

tseries <- read.csv("Data/Vitoria.data.csv")

data <- ddply(tseries, .(tot.week), summarise,
              obs = sum(Cases, na.rm = T),
              q = sum(Mosquitoes, na.rm = T),
              tau = sum(Trap, na.rm = T),
              year = unique(Year),
              week = unique(Week))

# Population size
pop <- 327801

# Wrangling weather data and computing EIP
weather <- read.csv("Data/Vitoria.weather.csv") %>% 
  mutate(date = ymd(BRST), year = year(date), week = week(date))

weather <- subset(weather, date < date[1] + weeks(243))
weather$tot.week <- rep(c(1:243), each = 7)

covars <- ddply(weather, .(tot.week), summarise, 
                temp = mean(Mean.TemperatureC, na.rm = T))

rov <- 7 * exp(0.2 * covars$temp - 8)

#===============================================================================
# Loading mcmc samples

fit <- readRDS("Results/gamma_eip_dv0.rds")

samples <- rstan::extract(fit)[1:17]

extract_sample <- function(x, k){
  if(length(dim(x)) == 1) return(x[k])
  else return(x[k, ])
}

# Selecting a sample to run the simulation with
init <- list(lapply(samples, extract_sample, 1))

# Loading Stan code
model <- stan_model(file = "Code/gammaeip.stan")

#===============================================================================
# Simulating dynamics with smaller and smaller Euler steps

# Number of Euler steps per week
steps <- c(7, 14, 28, 56, 112, 224)

# Loop over each of the step sizes and simulate dynamics
runs <- foreach(k = 1:6, .packages = c("rstan", "magrittr"), .combine = "rbind") %do% {
  
  dat.stan <- list(T = 243,
                   steps = steps[k],
                   y =data$obs,
                   q = data$q,
                   tau = data$tau,
                   rov = rov,
                   control = matrix(1, nrow = 243, ncol = 3),
                   pop = pop)
  
  sim <- sampling(model,
                  data = dat.stan, 
                  init = init, 
                  iter = 1, 
                  chains = 1,
                  warmup = 0,
                  algorithm = "Fixed_param")
  
  cases <- diff(rstan::extract(sim, "state", permute = T)[[1]][1, , 11])
  
  return(cases)
}


#===============================================================================
# Compare number of cases per week for the different step sizes

plot(runs[6, ])
points(runs[5, ], pch = 20, col = "blue")

plot(runs[5, ])
points(runs[4, ], pch = 20, col = "blue")

plot(runs[4, ])
points(runs[3, ], pch = 20, col = "blue")

plot(runs[3, ])
points(runs[2, ], pch = 20, col = "blue")

plot(runs[2, ])
points(runs[1, ], pch = 20, col = "blue")
