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
              week = unique(Week)
)

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

fit <- readRDS("Results/oscillator.rds")

samples <- rstan::extract(fit)[1:16]

extract_sample <- function(x, k){
  if(length(dim(x)) == 1) return(x[k])
  else return(x[k, ])
}

init <- list(lapply(samples, extract_sample, 1))

model <- stan_model(file = "Code/seirs.stan")

#===============================================================================
# Simulating dynamics with smaller and smaller Euler steps

steps <- c(7, 14, 28, 56, 112, 224)

# Parallel for-loop over mcmc iterations
runs <- foreach(k = 1:6, .packages = c("rstan", "magrittr"), .combine = "rbind") %do% {
  
  dat.stan <- list(T = 243,
                   T_pred = 0,
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
  
  cases <- rstan::extract(sim, "system", permute = T)[[1]][1, , 8]
  
  return(cases)
}


#===============================================================================

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
