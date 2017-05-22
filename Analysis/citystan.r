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
weather <- read.csv("Data/Vitoria.weather.csv") %>% 
  mutate(date = ymd(BRST), year = year(date), week = week(date))

weather <- subset(weather, date < date[1] + weeks(243))
weather$tot.week <- rep(c(1:243), each = 7)

covars <- ddply(weather, .(tot.week), summarise,
                temp = mean(Mean.TemperatureC, na.rm = T),
                humidity = mean(Mean.Humidity, na.rm = T)
                )

rov <- 7 / exp(exp(1.9 - 0.04 * covars$temp) + 1 / (2 * 7))

#===============================================================================
# Running stan

dat.stan <- list(T = 243,
                 D = 2,
                 y = obs,
                 q = q,
                 tau = tau,
                 rov = rov,
                 sincos = matrix(c(sin(2 * pi * c(1:243) / 52), 
                                   cos(2 * pi * c(1:243) / 52)), 
                                 ncol = 2),
                 covars = scale(covars[, -1]),
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
                  alpha0 = 0,
                  alpha = c(4, 3),
                  beta0 = 0.1,
                  beta = c(0.2, 0.1),
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
                  alpha0 = 1,
                  alpha = c(1, 1),
                  beta0 = 0.1,
                  beta = c(0, 0),
                  sigmad = 0.5,
                  z_d = rep(0, 243)))

fit <- stan(file = "Code/seirs.stan", 
            data = dat.stan, 
            init = inits, 
            iter = 500, 
            chains = 3,
            control = list(adapt_delta = 0.9))
