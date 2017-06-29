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
                 week = c(1:243),
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
                  ar_psi = 0,
                  ar_rv = 0,
                  alpha0 = 0,
                  beta0 = 0,
                  alpha = c(0, 0),
                  beta = c(0, 0),
                  sigmapsi = 0.2,
                  sigmarv = 0.2,
                  eps_psi = rep(0, T_fit),
                  rv = rep(0, T_fit)),
             list(S0 = 0.5,
                  E0 = 40,
                  I0 = 10,
                  log_phi_q = -12,
                  eta_inv_y = 1,
                  eta_inv_q = 1,
                  ro_c = 1,
                  gamma_c = 1.2,
                  delta_c = 1,
                  logNv = 0.6,
                  ar_psi = 0,
                  ar_rv = 0,
                  alpha0 = 0,
                  beta0 = 0,
                  alpha = c(0, 0),
                  beta = c(0, 0),
                  sigmapsi = 0.5,
                  sigmarv = 0.5,
                  eps_psi = rep(0, T_fit),
                  rv = rep(0, T_fit)),
             list(S0 = 0.3,
                  E0 = 100,
                  I0 = 80,
                  log_phi_q = -14,
                  eta_inv_y = 5,
                  eta_inv_q = 5,
                  ro_c = 1.2,
                  gamma_c = 0.8,
                  delta_c = 0.8,
                  logNv = 1,
                  ar_psi = 0,
                  ar_rv = 0,
                  alpha0 = 0,
                  beta0 = 0,
                  alpha = c(0, 0),
                  beta = c(0, 0),
                  sigmapsi = 0.05,
                  sigmarv = 0.05,
                  eps_psi = rep(0, T_fit),
                  rv = rep(0, T_fit)))

fit <- stan(file = "Code/seirs.stan", 
            data = dat.stan, 
            init = inits, 
            iter = 2000, 
            chains = 3,
            control = list(adapt_delta = 0.9,
                           max_treedepth = 12))

yhat <- rstan::extract(fit, "y_hat", permute = F) %>% apply(3, quantile, probs = c(0.1, 0.5, 0.9))
plot(yhat[3, ], type = "l")
lines(yhat[1, ])
lines(yhat[2, ])
points(obs, pch = 20)

qhat <- rstan::extract(fit, "q_hat", permute = F) %>% apply(3, quantile, probs = c(0.1, 0.5, 0.9))
plot(qhat[3, ], type = "l")
lines(qhat[1, ])
lines(qhat[2, ])
points(q, pch = 20)

dv <- rstan::extract(fit, "dv", permute = F) %>% apply(3, mean)
plot(dv, type = "l")

rv <- rstan::extract(fit, "rv", permute = F) %>% apply(3, mean)
plot(rv, type = "l")

psi <- rstan::extract(fit, "psi", permute = F) %>% apply(3, mean)
plot(psi, type = "l")

eps <- rstan::extract(fit, "eps_psi", permute = F) %>% apply(3, mean)
plot(eps, type = "l")


