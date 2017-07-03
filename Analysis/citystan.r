library(rstan)
library(plyr)
library(magrittr)
library(lubridate)

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

temp <- ddply(weather, .(tot.week), summarise, temp = mean(Mean.TemperatureC, na.rm = T))

rov <- 7 * exp(0.2 * temp$temp - 8)
#===============================================================================
# Running stan

dat.stan <- list(T = 243,
                 T_pred = 0,
                 y =data$ obs,
                 q = data$q,
                 tau = data$tau,
                 rov = rov,
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
                  seas_psi = 0,
                  seas_rv = 0,
                  sigmapsi = 0.2,
                  sigmarv = 0.2,
                  eps_psi = rep(0, 243),
                  eps_rv = rep(0, 243)),
             list(S0 = 0.3,
                  E0 = 100,
                  I0 = 80,
                  log_phi_q = -13.5,
                  eta_inv_y = 1,
                  eta_inv_q = 1,
                  ro_c = 1,
                  gamma_c = 1.2,
                  delta_c = 1,
                  logNv = 1,
                  ar_psi = 0.8,
                  ar_rv = 0.8,
                  seas_psi = 0.8,
                  seas_rv = 0.8,
                  sigmapsi = 0.5,
                  sigmarv = 0.5,
                  eps_psi = rep(0, 243),
                  eps_rv = rep(0, 243)),
             list(S0 = 0.5,
                  E0 = 40,
                  I0 = 20,
                  log_phi_q = -12.5,
                  eta_inv_y = 5,
                  eta_inv_q = 5,
                  ro_c = 1,
                  gamma_c = 1,
                  delta_c = 1,
                  logNv = 0.5,
                  ar_psi = 0.5,
                  ar_rv = 0.5,
                  seas_psi = 0,
                  seas_rv = 0,
                  sigmapsi = 0.01,
                  sigmarv = 0.01,
                  eps_psi = rep(0, 243),
                  eps_rv = rep(0, 243)))

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

qhat <- rstan::extract(fit, "q_hat", permute = F) %>% apply(3, quantile, probs = c(0.1, 0.5, 0.9))
plot(qhat[3, ], type = "l")
lines(qhat[1, ])
lines(qhat[2, ])
points(q, pch = 20)

rv <- rstan::extract(fit, "rv", permute = F) %>% apply(3, mean)
plot(rv, type = "l")

psi <- rstan::extract(fit, "psi", permute = F) %>% apply(3, mean)
plot(psi, type = "l")
