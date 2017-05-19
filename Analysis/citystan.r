library(rstan)
library(plyr)
library(magrittr)

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
weather <- read.csv("Data/Vitoria.weather.csv")
temp <- weather$Mean.TemperatureC

week <- rep(c(1:243), each = 7)
meantemp <- tapply(temp[1:1701], week, mean, na.rm = T)

rov <- 7 / exp(exp(1.9 - 0.04 * meantemp) + 1 / (2 * 7))

#===============================================================================
# Running stan

dat.stan <- list(T = 243,
                 y = obs,
                 q = q,
                 tau = tau,
                 rov = rov,
                 sincos = matrix(c(sin(2 * pi * c(1:243) / 52),
                                   cos(2 * pi * c(1:243) / 52)),
                                 ncol = 2),
                 pop = pop)

inits = list(list(p0_raw = c(-0.4, -9, -9),
                  log_phi_q = -13,
                  eta_inv_y = 5,
                  eta_inv_q = 5,
                  ro_c = 1,
                  gamma_c = 1,
                  delta_c = 1,
                  logNv0 = 0.7,
                  alpha0 = 1,
                  alpha = c(0, 0),
                  beta0 = 0.39,
                  beta = c(0, 0),
                  sigmad = 0.1,
                  z_d = rep(0, 243)),
             list(p0_raw = c(-0.4, -9, -9),
                  log_phi_q = -13,
                  eta_inv_y = 10,
                  eta_inv_q = 10,
                  ro_c = 1,
                  gamma_c = 1,
                  delta_c = 1,
                  logNv0 = 0.7,
                  alpha0 = 0.5,
                  alpha = c(0.01, 0.01),
                  beta0 = 0.39,
                  beta = c(0.01, 0.01),
                  sigmad = 0.1,
                  z_d = rep(0, 243)),
             list(p0_raw = c(-0.4, -9, -9),
                  log_phi_q = -13,
                  eta_inv_y = 2,
                  eta_inv_q = 2,
                  ro_c = 1,
                  gamma_c = 1,
                  delta_c = 1,
                  logNv0 = 0.7,
                  alpha0 = 0,
                  alpha = c(0.5, 0.5),
                  beta0 = 0.1,
                  beta = c(0, 0),
                  sigmad = 0.5,
                  z_d = rep(0, 243)))

fit <- stan(file = "Code/seirs.stan", 
            data = dat.stan, 
            init = inits, 
            iter = 5000, 
            chains = 3,
            control = list(adapt_delta = 0.9))
