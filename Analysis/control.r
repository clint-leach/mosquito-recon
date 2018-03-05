library(rstan)
library(doParallel)
library(plyr)
library(magrittr)
library(lubridate)

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

fit <- readRDS("Results/euler7.rds")
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
cl <- makeCluster(10, type = "SOCK")
registerDoParallel(cl)

# Parallel for-loop over mcmc iterations
reduction <- foreach(k = 1:nmcmc, .combine = "cbind", .packages = c("rstan", "magrittr")) %dopar% {
  
  init <- list(lapply(samples, extract_sample, k))
  
  cases <- matrix(NA, nrow = 52, ncol = 1)
  
  for(i in 1:52){
    
    control <- matrix(1, nrow = 243, ncol = 3)
    control[(data$week == i & data$year > 2008 & data$year < 2012), 3] <- 0.905
    
    first <- which(rowSums(control) != 3)[1]
    
    dat.stan <- list(T = 243,
                     steps = 7,
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
    
    cases[i] <- rstan::extract(sim, "y_meas", permute = F)[, 1, ] %>%
      extract(first:243) %>% 
      sum() %>% 
      divide_by(sum(data$obs[first:243]))
  }
  
  return(cases)
}

stopCluster(cl)

saveRDS(reduction, "Results/adult_control.rds")

#===============================================================================
# In which week is larval control most effective?

# Setting up parallel
cl <- makeCluster(10, type = "SOCK")
registerDoParallel(cl)

# Parallel for-loop over mcmc iterations
reduction <- foreach(k = 1:nmcmc, .combine = "cbind", .packages = c("rstan", "magrittr")) %dopar% {
  
  init <- list(lapply(samples, extract_sample, k))
  
  cases <- matrix(NA, nrow = 52, ncol = 1)
  
  for(i in 1:52){
    
    control <- matrix(1, nrow = 243, ncol = 3)
    control[(data$week == i & data$year > 2008 & data$year < 2012), 2] <- 0.95
    
    first <- which(rowSums(control) != 3)[1]
    
    dat.stan <- list(T = 243,
                     steps = 7,
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
    
    cases[i] <- rstan::extract(sim, "y_meas", permute = F)[, 1, ] %>%
      extract(first:243) %>% 
      sum() %>% 
      divide_by(sum(data$obs[first:243]))
  }
  
  return(cases)
}

stopCluster(cl)

saveRDS(reduction, "Results/larval_control.rds")

#===============================================================================
# Dynamics when control deployed at optimum for adult control

# Setting up parallel
cl <- makeCluster(10, type = "SOCK")
registerDoParallel(cl)

# Parallel for-loop over mcmc iterations
aopt <- foreach(k = 1:nmcmc, .packages = c("rstan", "magrittr")) %dopar% {
  
  init <- list(lapply(samples, extract_sample, k))
  
  control <- matrix(1.0, nrow = 243, ncol = 3)
  control[(data$week == 27 & data$year > 2008 & data$year < 2012), 3] <- 0.905
  
  dat.stan <- list(T = 243,
                   steps = 7,
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
  
  system <- rstan::extract(sim, "state", permute = T)[[1]][1, , ]
  
  return(system)
}

stopCluster(cl)

saveRDS(aopt, "Results/adult_optimal.rds")

#===============================================================================
# Dynamics when control deployed at optimum for larval control

# Setting up parallel
cl <- makeCluster(10, type = "SOCK")
registerDoParallel(cl)

# Parallel for-loop over mcmc iterations
lopt <- foreach(k = 1:nmcmc, .packages = c("rstan", "magrittr")) %dopar% {
  
  init <- list(lapply(samples, extract_sample, k))
  
  control <- matrix(1.0, nrow = 243, ncol = 3)
  control[(data$week == 4 & data$year > 2008 & data$year < 2012), 2] <- 0.95
  
  dat.stan <- list(T = 243,
                   steps = 7,
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
  
  system <- rstan::extract(sim, "state", permute = T)[[1]][1, , ]
  
  return(system)
}

stopCluster(cl)

saveRDS(lopt, "Results/larval_optimal.rds")

#===============================================================================
# At what case threshold is control most effective?

thresholds <- seq(20, 200, by = 10)

# Setting up parallel
cl <- makeCluster(10, type = "SOCK")
registerDoParallel(cl)

# Parallel for-loop over mcmc iterations
reduction <- foreach(k = 1:nmcmc, .combine = "cbind", .packages = c("rstan", "magrittr", "plyr")) %dopar% {
  
  init <- list(lapply(samples, extract_sample, k))
  
  cases <- matrix(NA, nrow = length(thresholds), ncol = 1)
  
  for(i in 1:length(thresholds)){
    
    data$above <- data$obs > thresholds[i]
    
    weeks <- ddply(subset(data, year > 2008 & year < 2012), .(year), summarise, 
                   week = tot.week[which(above)[1]])
    
    control <- matrix(1, nrow = 243, ncol = 3)
    control[weeks$week, 1] <- 1.05
    
    dat.stan <- list(T = 243,
                     steps = 7,
                     y =data$obs,
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
    
    cases[i] <- rstan::extract(sim, "y_meas", permute = F)[, 1, ] %>%
      extract(first:243) %>% 
      sum() %>% 
      divide_by(sum(data$obs[first:243]))  }
  
  return(cases)
}

stopCluster(cl)

saveRDS(reduction, "Results/case_control.rds")

#===============================================================================
# At what mosquito/trap threshold is control most effective?

thresholds <- seq(0.25, 1, by = 0.04)

data <- mutate(data, qt = q/tau)

# Setting up parallel
cl <- makeCluster(10, type = "SOCK")
registerDoParallel(cl)

# Parallel for-loop over mcmc iterations
reduction <- foreach(k = 1:nmcmc, .combine = "cbind", .packages = c("rstan", "magrittr", "plyr")) %dopar% {
  
  init <- list(lapply(samples, extract_sample, k))
  
  cases <- matrix(NA, nrow = length(thresholds), ncol = 1)
  
  for(i in 1:length(thresholds)){
    
    data$above <- data$qt > thresholds[i]
    
    weeks <- ddply(subset(data, year > 2008 & year < 2012), .(year), summarise, 
                   week = tot.week[which(above)[1]])
    
    control <- matrix(1, nrow = 243, ncol = 3)
    control[weeks$week, 1] <- 1.05
    
    dat.stan <- list(T = 243,
                     steps = 7,
                     y =data$obs,
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
    
    cases[i] <- rstan::extract(sim, "y_meas", permute = F)[, 1, ] %>%
      extract(first:243) %>% 
      sum() %>% 
      divide_by(sum(data$obs[first:243]))  }
  
  return(cases)
}

stopCluster(cl)

saveRDS(reduction, "Results/mosq_control.rds")

