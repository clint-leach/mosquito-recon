# Script to simulate control interventions in fitted dengue transmission model

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

ends <- daply(data, .(year), function(df) max(df$tot.week))

# Population size
pop <- 327801

# Wrangling weather data and computing EIP
weather <- read.csv("Data/Vitoria.weather.csv") %>% 
  mutate(date = ymd(BRST), year = year(date), week = week(date))

weather <- subset(weather, date < date[1] + weeks(243))
weather$tot.week <- rep(c(1:243), each = 7)

covars <- ddply(weather, .(tot.week), summarise, 
                temp = mean(Mean.TemperatureC, na.rm = T))

rov <- 7 * exp(0.21 * covars$temp - 7.9)

#===============================================================================
# Loading mcmc samples

fit <- readRDS("Results/chain_phi12.rds")
samples <- rstan::extract(fit)[1:17]
fitcases <- rstan::extract(fit, "state", permute = T)[[1]] %>% extract(, ,  11)

nmcmc <- length(samples[[1]])

extract_sample <- function(x, k){
  if(length(dim(x)) == 1) return(x[k])
  else return(x[k, ])
}

model <- stan_model(file = "Code/gammaeip.stan")

#===============================================================================
# Adult control simulations

# Loops over each HMC iteration and simulates a 5% decresase in total mosquito 
# population size in each of the first 210 weeks of the time series.
# Summarizes the difference in number of cases observed in each year between the
# control and uncontrolled time series.

# Setting up parallel
cl <- makeCluster(10, type = "SOCK")
registerDoParallel(cl)

# Parallel for-loop over mcmc iterations
reduction <- foreach(k = 1:nmcmc, .combine = "rbind", .packages = c("rstan", "magrittr")) %:% 
  foreach(j = 1:210, .combine = "rbind") %dopar% {

    # Extract parameter values in kth iteration
    init <- list(lapply(samples, extract_sample, k))

    # Set-up control structure
    control <- matrix(1, nrow = 243, ncol = 3)
    control[j, 3] <- 0.95
    
    # Stan data object
    dat.stan <- list(T = 243,
                     steps = 7,
                     y = data$obs,
                     q = data$q,
                     tau = data$tau,
                     rov = rov,
                     control = control,
                     pop = pop)
    
    # Simulate trajectory with control
    sim <- sampling(model,
                    data = dat.stan, 
                    init = init, 
                    iter = 1, 
                    chains = 1,
                    warmup = 0,
                    algorithm = "Fixed_param")
    
    # Extract values of the state variables at the beginning of each year
    state <- rstan::extract(sim, "state", permute = T)[[1]][1, ends + 1, ]
    
    # Assemble results
    data.frame(rep = k,
               control = j,
               year = c(2009, 2010, 2011, 2012),
               S = state[1:4, 1],
               I = state[1:4, 3],
               reduction = state[, 11] %>%
                 diff() %>%
                 divide_by(diff(fitcases[k, ends + 1])))
}

stopCluster(cl)

saveRDS(reduction, "Results/gamma_control.rds")
