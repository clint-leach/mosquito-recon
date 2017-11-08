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

nmcmc <- length(samples[[1]])

extract_sample <- function(x, k){
  if(length(dim(x)) == 1) return(x[k])
  else return(x[k, ])
}

model <- stan_model(file = "Code/seirs.stan")

#===============================================================================
# In which week is adult control most effective?

# Setting up parallel
cl <- makeCluster(2, type = "SOCK")
registerDoParallel(cl)

# Parallel for-loop over mcmc iterations
reduction <- foreach(k = 1:nmcmc, .combine = "cbind", .packages = c("rstan", "magrittr")) %dopar% {
  
  init <- list(lapply(samples, extract_sample, k))
  
  cases <- matrix(NA, nrow = 52, ncol = 1)
  
  for(i in 1:52){
    
    control <- matrix(1, nrow = 243, ncol = 3)
    control[(data$week == i & data$year > 2008 & data$year < 2012), 1] <- 1.05
    
    dat.stan <- list(T = 243,
                     T_pred = 0,
                     y = data$obs,
                     q = data$q,
                     tau = data$tau,
                     rov = rov,
                     control = control,
                     pop = pop)
    
    sim <- sampling(model,
                    data = dat.stan, 
                    init = init, 
                    iter = 1, 
                    chains = 1,
                    warmup = 0,
                    algorithm = "Fixed_param")
    
    cases[i] <- rstan::extract(sim, "y_hat", permute = F)[, 1, ] %>% sum()
  }
  
  return(cases)
}

stopCluster(cl)

reduction <- sum(data$obs) - reduction

saveRDS(reduction, "Results/control.rds")

#===============================================================================
# Dynamics when control deployed in week 3

# Setting up parallel
cl <- makeCluster(3, type = "SOCK")
registerDoParallel(cl)

# Parallel for-loop over mcmc iterations
week3 <- foreach(k = 1:nmcmc, .packages = c("rstan", "magrittr")) %dopar% {
  
  init <- list(lapply(samples, extract_sample, k))
  
  control <- matrix(1.0, nrow = 243, ncol = 3)
  control[(data$week == 3 & data$year > 2008 & data$year < 2012), 1] <- 1.05
  
  dat.stan <- list(T = 243,
                   T_pred = 0,
                   y = data$obs,
                   q = data$q,
                   tau = data$tau,
                   rov = rov,
                   control = control,
                   pop = pop)
  
  sim <- sampling(model,
                  data = dat.stan, 
                  init = init, 
                  iter = 1, 
                  chains = 1,
                  warmup = 0,
                  algorithm = "Fixed_param")
  
  system <- rstan::extract(sim, "system", permute = T)[[1]][1, , ]
  
  return(system)
}

stopCluster(cl)

saveRDS(week3, "Results/optimal.rds")

#===============================================================================

peak_cases <- ddply(data, .(year), summarise, max = which.max(obs))

# Setting up parallel
cl <- makeCluster(3, type = "SOCK")
registerDoParallel(cl)

# Parallel for-loop over mcmc iterations
week12 <- foreach(k = 1:nmcmc, .packages = c("rstan", "magrittr")) %dopar% {
  
  init <- list(lapply(samples, extract_sample, k))
  
  control <- matrix(1.0, nrow = 243, ncol = 3)
  control[(data$week == 12 & data$year > 2008 & data$year < 2012), 1] <- 1.05
  
  dat.stan <- list(T = 243,
                   T_pred = 0,
                   y = data$obs,
                   q = data$q,
                   tau = data$tau,
                   rov = rov,
                   control = control,
                   pop = pop)
  
  sim <- sampling(model,
                  data = dat.stan, 
                  init = init, 
                  iter = 1, 
                  chains = 1,
                  warmup = 0,
                  algorithm = "Fixed_param")
  
  system <- rstan::extract(sim, "system", permute = T)[[1]][1, , ]
  
  return(system)
}

stopCluster(cl)

saveRDS(week12, "Results/peak.rds")

