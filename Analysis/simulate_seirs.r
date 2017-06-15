library(rstan)
library(plyr)
library(magrittr)
library(lubridate)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#===============================================================================
# Loading data and setting up model

tseries <- read.csv("Data/Vitoria.data.csv")

# Population size
pop <- 327801

# Pull out mosquito counts and traps inspected
tau <- tapply(tseries$Trap, list(tseries$tot.week), sum, na.rm = TRUE)

# Computing eip forcing
weather <- read.csv("Data/Vitoria.weather.csv") %>% 
  mutate(date = ymd(BRST), year = year(date), week = week(date))

weather <- subset(weather, date < date[1] + weeks(243))
weather$tot.week <- rep(c(1:243), each = 7)

covars <- ddply(weather, .(tot.week), summarise,
                temp = mean(Mean.TemperatureC, na.rm = T), 
                humidity = mean(Mean.Humidity, na.rm = T))

covars <- covars[, -1]

rov <- 7 / exp(exp(1.9 - 0.04 * covars$temp) + 1 / (2 * 7))

#===============================================================================
# Simulating a new data set

dat.sim <- list(T = 243,
                 T_pred = 0,
                 D = dim(covars)[2],
                 y = rep(0, 243),
                 q = rep(0, 243),
                 tau = tau,
                 rov = rov,
                 covars = scale(covars),
                 pop = pop)

pars = list(list(p0_raw = c(0, -6, -6),
                  log_phi_q = -13,
                  eta_inv_y = 0.1,
                  eta_inv_q = 0.1,
                  ro_c = 0.7,
                  gamma_c = 1.2,
                  delta_c = 0.8,
                  logNv =0.7,
                  alpha1 = 0.4,
                  alpha2 = 0.4,
                  beta0 = 1,
                  beta = c(0.3, -0.05),
                  sigmad = 0.2,
                  z_d = rnorm(243, 0, 1)))

sim <- stan(file = "Code/seirs.stan", 
            data = dat.sim, 
            init = pars, 
            iter = 1, 
            chains = 1,
            algorithm = "Fixed_param",
            control = list(adapt_delta = 0.9))

ynew <- rstan::extract(sim, "y_hat", permute = F) %>% apply(3, mean)
plot(ynew, type = "l")

qnew <- rstan::extract(sim, "q_hat", permute = F) %>% apply(3, mean)
plot(qnew, type = "l")

dvtrue <- rstan::extract(sim, "dv", permute = F) %>% apply(3, mean)
#===============================================================================
# Trying to fit to the simulated data

T_pred <- 52
T_fit <- 243 - T_pred

dat.new <- list(T = T_fit,
                 T_pred = T_pred,
                 D = dim(covars)[2],
                 y = head(ynew, T_fit),
                 q = head(qnew, T_fit),
                 tau = head(tau, T_fit),
                 rov = rov,
                 covars = scale(covars),
                 pop = pop)

inits = list(list(p0_raw = c(-0.3, -8, -8),
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
                  z_d = rep(0, T_fit)),
             list(p0_raw = c(-0.4, -8, -8),
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
                  z_d = rep(0, T_fit)),
             list(p0_raw = c(-0.2, -8, -8),
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
                  z_d = rnorm(T_fit, 0, 1)))

new <- stan(file = "Code/seirs.stan", 
            data = dat.new, 
            init = inits, 
            iter = 500, 
            chains = 3,
            control = list(adapt_delta = 0.9))

yhat <- rstan::extract(new, "y_hat", permute = F) %>% apply(3, quantile, probs = c(0.1, 0.5, 0.9))
qhat <- rstan::extract(new, "q_hat", permute = F) %>% apply(3, mean)

dvhat <- rstan::extract(new, "dv", permute = F) %>% apply(3, mean)
dvpred <- rstan::extract(new, "d_pred", permute = F) %>% apply(3, mean)

plot(dvtrue, type = "l")
lines(dvhat, col = "blue")
lines(c(192:243), dvpred, col = "red")

plot(ynew, pch = 20, type = "b", col = "gray70")
lines(yhat[1, ])
lines(yhat[3, ])
lines(yhat[2, ])

plot(qnew, pch = 20)
lines(qhat)
