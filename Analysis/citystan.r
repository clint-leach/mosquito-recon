library(rstan)
library(plyr)
library(magrittr)
library(lubridate)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#===============================================================================
# Loading data and setting up model
  
tseries <- read.csv("Data/Vitoria.data.csv")
obs <- tapply(tseries$Cases, list(tseries$tot.week), sum, na.rm=TRUE)

# Population size
pop <- 327801

# Pull out mosquito counts and traps inspected
q <- tapply(tseries$Mosquitoes, list(tseries$tot.week), sum, na.rm = TRUE)
tau <- tapply(tseries$Trap, list(tseries$tot.week), sum, na.rm = TRUE)

# Computing eip forcing
weather <- read.csv("Data/Vitoria.weather.csv") %>% 
  mutate(date = ymd(BRST), year = year(date), week = week(date))

weather <- subset(weather, date < date[1] + weeks(243))
weather$tot.week <- rep(c(1:243), each = 7)

covars <- ddply(weather, .(tot.week), summarise,
                temp = mean(Mean.TemperatureC, na.rm = T),
                humidity = mean(Mean.Humidity, na.rm = T)
                )

covars <- covars[, -1]

rov <- 7 / exp(exp(1.9 - 0.04 * covars$temp) + 1 / (2 * 7))

#===============================================================================
# Running stan

dat.stan <- list(T = 243,
                 T_pred = 0,
                 D = dim(covars)[2],
                 y = obs,
                 q = q,
                 tau = tau,
                 rov = rov,
                 covars = scale(covars),
                 pop = pop)

inits = list(list(p0_raw = c(0.6, -1, 0.3, -15, 0.2),
                  log_phi_q = -13,
                  eta_inv_y = 5,
                  eta_inv_q = 5,
                  ro_c = 1,
                  gamma_c = 1,
                  delta_c = 1,
                  logNv0 = 0.7,
                  alpha0 = 1,
                  alpha1 = 1,
                  alpha2 = 1,
                  beta0 = 1,
                  beta = c(0.3, 0.1),
                  sigmad = 0.1,
                  z_d = rep(0, 243)),
             list(p0_raw = c(0.3, -1, 0.2, -15, 0.1),
                  log_phi_q = -13,
                  eta_inv_y = 10,
                  eta_inv_q = 10,
                  ro_c = 1,
                  gamma_c = 1,
                  delta_c = 1,
                  logNv0 = 1,
                  alpha0 = 0,
                  alpha1 = 2,
                  alpha2 = 3,
                  beta0 = 1,
                  beta = c(0.2, -0.1),
                  sigmad = 0.1,
                  z_d = rep(0, 243)),
             list(p0_raw = c(0.8, -1, 0.3, -15, 0.2),
                  log_phi_q = -13,
                  eta_inv_y = 2,
                  eta_inv_q = 2,
                  ro_c = 1,
                  gamma_c = 1,
                  delta_c = 1,
                  logNv0 = 0.5,
                  alpha0 = 1,
                  alpha1 = 0,
                  alpha2 = 2,
                  beta0 = 1,
                  beta = c(0.5, 0.05),
                  sigmad = 0.5,
                  z_d = rep(0, 243)))

fit <- stan(file = "Code/twosero.stan", 
            data = dat.stan, 
            init = inits, 
            iter = 1000, 
            chains = 3,
            control = list(adapt_delta = 0.9))
