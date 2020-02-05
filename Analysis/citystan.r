# Script to fit an SEIRS-SEI dengue transmission model to weekly trap counts
# and reports of "dengue-like illness" using Stan.

library(rstan)
library(plyr)
library(magrittr)
library(lubridate)

# Setting Stan options
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#===============================================================================
# Loading data and setting up model

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

# Weekly mean temperature
covars <- ddply(weather, .(tot.week), summarise, 
                temp = mean(Mean.TemperatureC, na.rm = T),
                humid = mean(Mean.Humidity, na.rm = T))

# Computing inverse of extrinsic incubation period
rov <- 7 * exp(0.21 * covars$temp - 7.9)

#===============================================================================
# Running stan

# Wrapping up all the data
dat.stan <- list(T = 243,
                 steps = 7,
                 y =data$obs,
                 q = data$q,
                 tau = data$tau,
                 rov = rov,
                 control = matrix(1.0, nrow = 243, ncol = 3),
                 pop = pop)

# Specifying initial conditions for 3 chains
inits = list(list(S0 = 0.4,
                  E0 = 80,
                  I0 = 40,
                  log_phi_q = -13,
                  phi_y = 1/12,
                  eta_inv_y = 0.1,
                  eta_inv_q = 0.1,
                  ro_c = 0.8,
                  gamma_c = 1,
                  delta_c = 0.8,
                  dvmu_c = 1, 
                  logNv = 1.5,
                  dv0 = 0,
                  rv0 = 0,
                  sigmadv = 0.2,
                  sigmarv = 0.2,
                  eps_dv = rep(0, 243),
                  eps_rv = rep(0, 243)),
             list(S0 = 0.3,
                  E0 = 100,
                  I0 = 80,
                  log_phi_q = -13.5,
                  phi_y = 1/6,
                  eta_inv_y = 1,
                  eta_inv_q = 1,
                  ro_c = 1,
                  gamma_c = 1.2,
                  delta_c = 1,
                  dvmu_c = 0.6,
                  logNv = 0.7,
                  dv0 = 0,
                  rv0 = 0,
                  sigmadv = 0.5,
                  sigmarv = 0.5,
                  eps_dv = rep(0, 243),
                  eps_rv = rep(0, 243)),
             list(S0 = 0.5,
                  E0 = 40,
                  I0 = 20,
                  log_phi_q = -12.5,
                  phi_y = 1/24,
                  eta_inv_y = 5,
                  eta_inv_q = 5,
                  ro_c = 0.9,
                  gamma_c = 0.5,
                  delta_c = 0.6,
                  dvmu_c = 1,
                  logNv = 1,
                  dv0 = 0,
                  rv0 = 0,
                  sigmadv = 0.01,
                  sigmarv = 0.01,
                  eps_dv = rep(0, 243),
                  eps_rv = rep(0, 243)))

# Running HMC
# Note non-default options for adapt_delta and max_treedepth, which 
# help ensure good mixing and avoid divergent transitions
# Also note that this will take on the order of 1 week to run

fit <- stan(file = "Code/gammaeip.stan", 
            data = dat.stan, 
            init = inits, 
            iter = 4000,
            chains = 3,
            control = list(adapt_delta = 0.99,
                           max_treedepth = 15))

saveRDS(fit, "Results/chain.rds")
