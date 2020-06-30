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

# Setting imports (if any)
firsts <- daply(data, .(year), function(df) min(df$tot.week))[2:5]
import <- rep(0, 243)
# import[firsts - 8] <- 10 / pop

model <- stan_model(file = "Code/gammaeip.stan", verbose = T, save_dso = TRUE, auto_write = TRUE)

#===============================================================================
# Loading mcmc samples

chain <- readRDS("Results/chain.rds")
samples <- rstan::extract(chain, permute = T)[1:18] 

fitcases <- rstan::extract(chain, "state", permute = T)[[1]] %>% 
  extract(, ,  11) %>% 
  apply(1, diff)

nmcmc <- length(samples[[1]])

extract_sample <- function(x, k){
  if(length(dim(x)) == 1) return(x[k])
  else return(x[k, ])
}


#===============================================================================
# Adult control simulations

# Loops over each HMC iteration and simulates a 5% decresase in total mosquito 
# population size in each of the first 210 weeks of the time series.
# Summarizes the difference in number of cases observed in each year between the
# control and uncontrolled time series.

# Setting up parallel
cl <- makeCluster(8, type = "SOCK")
registerDoParallel(cl)

# Parallel for-loop over mcmc iterations
reduction <- foreach(k = 1:nmcmc, .combine = "rbind", .packages = c("rstan", "magrittr", "lubridate")) %:% 
  foreach(j = 1:191, .combine = "rbind") %dopar% {

    # Extract parameter values in kth iteration
    inits <- list(lapply(samples, extract_sample, k))

    # Set-up control structure
    control <- matrix(1.0, nrow = 243, ncol = 3)
    control[j, 3] <- 0.95
    control[j + 1, 2] <- 1.05
    
    # Stan data object
    dat.stan <- list(T = 243,
                     steps = 7,
                     y = data$obs,
                     q = data$q,
                     tau = data$tau,
                     rov = rov,
                     control = control,
                     import = import,
                     pop = pop)
    
    # Simulate trajectory with control
    sim <- sampling(model,
                    data = dat.stan, 
                    init = inits, 
                    iter = 1, 
                    chains = 1,
                    warmup = 0,
                    algorithm = "Fixed_param")
    
    # Extract values of the state variables at the beginning of each year
    cases <- rstan::extract(sim, "state", permute = T)[[1]][1, , 11] %>% diff()

    # Assemble results
    data.frame(rep = k,
               control = j,
               dvmu = inits[[1]]$dv0,
               delta = inits[[1]]$delta_c,
               phi = inits[[1]]$phi_y,
               week = j:243,
               relweek = 0:(243 - j),
               year = year(ymd("2008-01-01") + weeks(j:243)),
               cases = cases[j:243],
               cases0 = fitcases[j:243, k])
}

stopCluster(cl)

saveRDS(reduction, "Results/control.rds")

# Effect of control in the number of cases in the following year
moving <- ddply(control, .(rep, control), summarise,
                tcases = sum(cases[relweek < 53]),
                tcases0 = sum(cases0[relweek < 53]),
                ratio = tcases / tcases0,
                diff = tcases0 - tcases,
                dvmu = mean(dvmu),
                delta = mean(delta),
                phi = mean(phi))

saveRDS(moving, "Results/moving_control.rds")
