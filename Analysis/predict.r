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
T_pred <- 52
T_fit <- 243 - T_pred

weather <- subset(weather, date < date[1] + weeks(243))
weather$tot.week <- rep(c(1:243), each = 7)

covars <- ddply(weather, .(tot.week), summarise,
                temp = mean(Mean.TemperatureC, na.rm = T), 
                humidity = mean(Mean.Humidity, na.rm = T)
                )

annual <- ddply(weather, .(week), summarise,
                temp = mean(Mean.TemperatureC, na.rm = T),
                humidity = mean(Mean.Humidity, na.rm = T)
                )

covars <- covars[, -1]
# covars <- rbind(covars[, -1],
#                 annual[week(start + weeks((T_fit + 1):243)), -1])

rov <- 7 / exp(exp(1.9 - 0.04 * covars$temp + 1 / (2 * 7)))

#===============================================================================
# Running stan

dat.stan <- list(T = T_fit,
                 T_pred = T_pred,
                 D = dim(covars)[2],
                 y = head(obs, T_fit),
                 q = head(q, T_fit),
                 tau = head(tau, T_fit),
                 rov = rov,
                 covars = scale(covars),
                 week = week(weather$date[1] + weeks(c(0:242))),
                 pop = pop)

inits = list(list(S0 = 0.5,
                  E0 = 100,
                  I0 = 90,
                  log_phi_q = -13,
                  eta_inv_y = 5,
                  eta_inv_q = 5,
                  ro_c = 1,
                  gamma_c = 1,
                  delta_c = 1,
                  logNv = 0.7,
                  alpha1 = 0,
                  alpha2 = 0,
                  beta0 = 1,
                  beta = c(0.3, 0.1),
                  sigmad = 0.1,
                  z_d = rep(0, T_fit),
                  nu = rep(0, 53)),
             list(S0 = 0.4,
                  E0 = 40,
                  I0 = 50,
                  log_phi_q = -13,
                  eta_inv_y = 1,
                  eta_inv_q = 1,
                  ro_c = 1.2,
                  gamma_c = 1.2,
                  delta_c = 1.2,
                  logNv1 = 1,
                  alpha1 = 1,
                  alpha2 = 1,
                  beta0 = 0.8,
                  beta = c(0, 0),
                  sigmad = 0.5,
                  z_d = rep(0, T_fit),
                  nu = rep(0, 53)),
             list(S0 = 0.3,
                  E0 = 150,
                  I0 = 50,
                  log_phi_q = -13,
                  eta_inv_y = 10,
                  eta_inv_q = 10,
                  ro_c = 0.8,
                  gamma_c = 0.8,
                  delta_c = 0.8,
                  logNv = 0.5,
                  alpha1 = 0.5,
                  alpha2 = 0.5,
                  beta0 = 1.2,
                  beta = c(0.3, -0.1),
                  sigmad = 0.05,
                  z_d = rnorm(T_fit, 0, 1),
                  nu = rnorm(53, 0, 0.1)))

fit <- stan(file = "Code/seirs.stan", 
            data = dat.stan, 
            init = inits, 
            iter = 1000, 
            chains = 3,
            control = list(adapt_delta = 0.9,
                           max_treedepth = 12))

yhat <- rstan::extract(fit, "y_hat", permute = F) %>% apply(3, quantile, probs = c(0.1, 0.5, 0.9))
plot(yhat[2, ], type = "l")
lines(yhat[1, ])
lines(yhat[3, ])
points(obs, pch = 20)

dv <- rstan::extract(fit, "dv", permute = F) %>% apply(3, mean)
plot(dv, type = "l")

