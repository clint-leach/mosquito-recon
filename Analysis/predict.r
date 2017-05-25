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

# Wrangling weather data and computing EIP
weather <- read.csv("Data/Vitoria.weather.csv") %>% 
  mutate(date = ymd(BRST), year = year(date), week = week(date))

start <- weather$date[1]
T_pred <- 34
T_fit <- 243 - T_pred

weather <- subset(weather, date < date[1] + weeks(T_fit))
weather$tot.week <- rep(c(1:T_fit), each = 7)

covars <- ddply(weather, .(tot.week), summarise,
                temp = mean(Mean.TemperatureC, na.rm = T), 
                humidity = mean(Mean.Humidity, na.rm = T)
                )

annual <- ddply(weather, .(week), summarise,
                temp = mean(Mean.TemperatureC, na.rm = T),
                humidity = mean(Mean.Humidity, na.rm = T)
                )

predvars <- rbind(covars[, -1],
                  annual[week(start + weeks((T_fit + 1):243)), -1])

predvars <- cbind(predvars, matrix(c(sin(2 * pi * c(1:243) / 52), 
                                 cos(2 * pi * c(1:243) / 52),
                                 sin(4 * pi * c(1:243) / 52), 
                                 cos(4 * pi * c(1:243) / 52),
                                 sin(8 * pi * c(1:243) / 52), 
                                 cos(8 * pi * c(1:243) / 52),
                                 sin(16 * pi * c(1:243) / 52), 
                                 cos(16 * pi * c(1:243) / 52)), 
                               ncol = 8))

rov <- 7 / exp(exp(1.9 - 0.04 * predvars$temp) + 1 / (2 * 7))

#===============================================================================
# Running stan

dat.stan <- list(T = T_fit,
                 T_pred = T_pred,
                 D = dim(predvars)[2],
                 y = head(obs, T_fit),
                 q = head(q, T_fit),
                 tau = head(tau, T_fit),
                 rov = rov,
                 sincos = matrix(c(sin(2 * pi * c(1:243) / 52),
                                   cos(2 * pi * c(1:243) / 52)),
                                 ncol = 2),
                 covars = scale(predvars),
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
                  beta = rep(0, 10),
                  sigmad = 0.1,
                  z_d = rep(0, T_fit)))

fit <- stan(file = "Code/predict.stan", 
            data = dat.stan, 
            init = inits, 
            iter = 5000, 
            chains = 1,
            control = list(adapt_delta = 0.9))

