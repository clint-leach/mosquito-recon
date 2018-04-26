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
              week = unique(Week),
              epiweek = 0,
              epiyear = 0
              )

# Generating epidemic weeks/years that count from peak to peak (week 16)
data$epiweek[16:243] <- c(1:52)
data$epiyear <- cumsum(data$epiweek == 1)

# Specifying the year windows over which we evaluate control
control_years = c(2009, 2010, 2011)

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
fitcases <- rstan::extract(fit, "system", permute = T)[[1]] %>% extract(, ,  8)

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
reduction <- foreach(k = 1:nmcmc, .combine = "rbind", .packages = c("rstan", "magrittr")) %dopar% {
  
  init <- list(lapply(samples, extract_sample, k))
  
  cases <- data.frame(rep = k, control = rep(1:52, each = 3), year = rep(control_years, 52), reduction = 0)
  
  for(i in 1:52){
    for(j in 1:3){
      
      control <- matrix(1, nrow = 243, ncol = 3)
      control[(data$epiweek == i & data$epiyear == j), 3] <- 0.95
      
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
      
      cases$reduction[(i - 1) * 3 + j] <- rstan::extract(sim, "state", permute = T)[[1]][1, , 8] %>% 
        diff() %>% 
        extract(data$year == control_years[j]) %>% 
        subtract(fitcases[k, data$year == control_years[j]]) %>% 
        sum()
    }
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
reduction <- foreach(k = 1:nmcmc, .combine = "rbind", .packages = c("rstan", "magrittr")) %dopar% {
  
  init <- list(lapply(samples, extract_sample, k))
  
  cases <- data.frame(rep = k, control = rep(1:52, each = 3), year = rep(control_years, 52), reduction = 0)
  
  for(i in 1:52){
    for(j in 1:3){
      
      control <- matrix(1, nrow = 243, ncol = 3)
      control[(data$epiweek == i & data$epiyear == j), 2] <- 0.95
      
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
      
      cases$reduction[(i - 1) * 3 + j] <- rstan::extract(sim, "state", permute = T)[[1]][1, , 8] %>% 
        diff() %>% 
        extract(data$year == control_years[j]) %>% 
        subtract(fitcases[k, data$year == control_years[j]]) %>% 
        sum()
    }
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
  control[(data$week == 23 & data$year > 2008 & data$year < 2012), 3] <- 0.905
  
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
