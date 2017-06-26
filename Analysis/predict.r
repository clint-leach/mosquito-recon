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
                humidity = mean(Mean.Humidity, na.rm = T),
                dtr = mean(Max.TemperatureC - Min.TemperatureC, na.rm = T),
                wind = mean(Mean.Wind.SpeedKm.h, na.rm = T)
                )

rov <- 7 * exp(0.2 * covars$temp - 8)

covars <- scale(covars[, -1])

#===============================================================================
# Running stan

dat.stan <- list(T = T_fit,
                 T_pred = T_pred,
                 D = dim(covars)[2],
                 y = head(obs, T_fit),
                 q = head(q, T_fit),
                 tau = head(tau, T_fit),
                 rov = rov,
                 covars = covars,
                 week = week(weather$date[1] + weeks(c(0:242))),
                 pop = pop)

inits = list(list(S0 = 0.4,
                  E0 = 80,
                  I0 = 40,
                  log_phi_q = -13,
                  eta_inv_y = 0.1,
                  eta_inv_q = 0.1,
                  ro_c = 0.8,
                  gamma_c = 1,
                  delta_c = 0.8,
                  logNv = 0.7,
                  beta0 = 1,
                  beta = c(0, 0, 0, 0),
                  alpha = c(0, 0, 0, 0),
                  sigmad = 0.2,
                  eps_d = rep(0, T_fit)),
             list(S0 = 0.6,
                  E0 = 30,
                  I0 = 30,
                  log_phi_q = -13.5,
                  eta_inv_y = 5,
                  eta_inv_q = 5,
                  ro_c = 1.1,
                  gamma_c = 0.8,
                  delta_c = 0.6,
                  logNv = 0.4,
                  beta0 = 1.2,
                  beta = c(0, 0, 0, 0),
                  alpha = c(0, 0, 0, 0),
                  sigmad = 0.05,
                  eps_d = rep(0, T_fit)),
             list(S0 = 0.3,
                  E0 = 100,
                  I0 = 60,
                  log_phi_q = -12,
                  eta_inv_y = 1,
                  eta_inv_q = 1,
                  ro_c = 1,
                  gamma_c = 1,
                  delta_c = 1,
                  logNv = 1,
                  beta0 = 0.8,
                  beta = c(0, 0, 0, 0),
                  alpha = c(0, 0, 0, 0),
                  sigmad = 0.5,
                  eps_d = rep(0, T_fit)))

fit <- stan(file = "Code/seirs.stan", 
            data = dat.stan, 
            init = inits, 
            iter = 5000, 
            chains = 3,
            control = list(adapt_delta = 0.9,
                           max_treedepth = 12))

yhat <- rstan::extract(fit, "y_hat", permute = F) %>% apply(3, quantile, probs = c(0.1, 0.5, 0.9))
plot(yhat[3, ], type = "l")
lines(yhat[1, ])
lines(yhat[2, ])
points(obs, pch = 20)

dv <- rstan::extract(fit, "dv", permute = F) %>% apply(3, mean)
plot(dv, type = "l")

bv <- rstan::extract(fit, "bv", permute = F) %>% apply(3, mean)
plot(bv, type = "l")


