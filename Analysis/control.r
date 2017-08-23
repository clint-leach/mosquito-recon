library(rstan)
library(plyr)
library(magrittr)
library(lubridate)
library(ggplot2)
library(gridExtra)

# Loading and processing mosquito and case data

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

covars <- ddply(weather, .(tot.week), summarise, 
                temp = mean(Mean.TemperatureC, na.rm = T),
                humid = mean(Mean.Humidity, na.rm = T))

rov <- 7 * exp(0.2 * covars$temp - 8)

#===============================================================================
# Posterior means

fit <- readRDS("Results/oscillator.rds")

means = list(list(S0 = get_posterior_mean(fit, pars = "S0")[4],
                  E0 = get_posterior_mean(fit, pars = "E0")[4],
                  I0 = get_posterior_mean(fit, pars = "I0")[4],
                  log_phi_q = get_posterior_mean(fit, pars = "log_phi_q")[4],
                  eta_inv_y = get_posterior_mean(fit, pars = "eta_inv_y")[4],
                  eta_inv_q = get_posterior_mean(fit, pars = "eta_inv_q")[4],
                  ro_c = get_posterior_mean(fit, pars = "ro_c")[4],
                  gamma_c = get_posterior_mean(fit, pars = "gamma_c")[4],
                  delta_c = get_posterior_mean(fit, pars = "delta_c")[4],
                  logNv = get_posterior_mean(fit, pars = "logNv")[4],
                  psi0 = get_posterior_mean(fit, pars = "psi0")[4],
                  rv0 = get_posterior_mean(fit, pars = "rv0")[4],
                  sigmapsi = 1, #get_posterior_mean(fit, pars = "sigmapsi")[4],
                  sigmarv = 1, #get_posterior_mean(fit, pars = "sigmarv")[4],
                  eps_psi = get_posterior_mean(fit, pars = "eps_psi")[, 4],
                  eps_rv = get_posterior_mean(fit, pars = "eps_rv")[, 4]))

#===============================================================================
# Setting control series

data$psi <- rstan::extract(fit, "psi", permute = T)[[1]] %>% apply(2, mean)
data$qhat <- rstan::extract(fit, "q_hat", permute = T)[[1]] %>% apply(2, mean)

# Identifying local peaks
psipeaks <- which(diff(sign(diff(data$psi))) == -2) + 1
qpeaks <- which(diff(sign(diff(data$qhat))) == -2) + 1

cpsi <- ddply(data[psipeaks, ], .(year), summarise, peak = tot.week[order(psi, decreasing = T)[1:2]])
cpsi <- cpsi[1:8, ]

cq <- ddply(data[qpeaks, ], .(year), summarise, peak = tot.week[order(qhat, decreasing = T)[1]])
cq <- cq[1:8, ]

control_psi <- matrix(1, nrow = 243, ncol = 3)
control_psi[c(98, 150, 198), 1] <- 2

control_q <- matrix(1, nrow = 243, ncol = 3)
control_q[c(57, 107, 158), 1] <- 2
#===============================================================================
# Simulate model

dat.psi <- list(T = 243,
                 T_pred = 0,
                 y = data$obs,
                 q = data$q,
                 tau = data$tau,
                 rov = rov,
                 control = control_psi,
                 pop = pop)

sim_psi <- stan(file = "Code/seirs.stan", 
            data = dat.psi, 
            init = means, 
            iter = 500, 
            chains = 1,
            algorithm = "Fixed_param")

cases_psi <- rstan::extract(sim_psi, "y_hat", permute = F) %>% apply(1, sum)
cases_psi <- cases_psi - sum(data$obs)
hist(cases_psi)

#===============================================================================

dat.q <- list(T = 243,
                 T_pred = 0,
                 y = data$obs,
                 q = data$q,
                 tau = data$tau,
                 rov = rov,
                 control = control_q,
                 pop = pop)

sim_q <- stan(file = "Code/seirs.stan", 
                data = dat.q, 
                init = means, 
                iter = 500, 
                chains = 1,
                algorithm = "Fixed_param")

cases_q <- rstan::extract(sim_q, "y_hat", permute = F) %>% apply(1, sum)
cases_q <- cases_q - sum(data$obs)
hist(cases_q)

#===============================================================================
# In which week is adult control most effective?

reduction <- data.frame(week = c(1:52),
                        med = rep(0, 52),
                        min = rep(0, 52),
                        max = rep(0, 52))

for(i in 1:52){
  
  control <- matrix(1, nrow = 243, ncol = 3)
  control[(data$week == i & data$year > 2008 & data$year < 2012), 1] <- 1.05
  
  dat.stan <- list(T = 243,
                   T_pred = 0,
                   y = data$obs,
                   q = data$q,
                   tau = data$tau,
                   rov = rov,
                   control = control,
                   pop = pop)
  
  sim <- stan(file = "Code/seirs.stan", 
              data = dat.stan, 
              init = means, 
              iter = 200, 
              chains = 1,
              warmup = 0,
              algorithm = "Fixed_param")
  
  totals[i] <- rstan::extract(sim, "y_hat", permute = F)[, 1, ] %>% rowSums() %>% mean()
}

plot(totals - sum(data$obs))
abline(h = 0)

plot(scale(data$qhat), type = "l")
lines(scale(data$psi), col = "red")
abline(v = data$tot.week[data$week == which.min(totals)], col = "blue")

#===============================================================================
# How does the effect change as we increase control effort?

strength <- seq(1, 1.5, by = 0.01)
effect_psi <- vector(length = length(strength))
effect_q <- vector(length = length(strength))
effect_2 <- vector(length = length(strength))

for(i in 1:length(strength)){
  
  control_psi <- matrix(1, nrow = 243, ncol = 3)
  control_psi[c(98, 150, 198), 1] <- strength[i]
  
  control_q <- matrix(1, nrow = 243, ncol = 3)
  control_q[c(57, 107, 158), 1] <- strength[i]
  
  control_2 <- matrix(1, nrow = 243, ncol = 3)
  control_2[(data$week == 2 & data$year > 2008 & data$year < 2012), 1] <- strength[i]
  
  dat.psi <- list(T = 243,
                   T_pred = 0,
                   y = data$obs,
                   q = data$q,
                   tau = data$tau,
                   rov = rov,
                   control = control_psi,
                   pop = pop)
  
  sim_psi <- stan(file = "Code/seirs.stan", 
              data = dat.psi, 
              init = means, 
              iter = 200, 
              chains = 1,
              warmup = 0,
              algorithm = "Fixed_param")
  
  effect_psi[i] <- rstan::extract(sim_psi, "y_hat", permute = F)[, 1, ] %>% rowSums() %>% mean() - sum(data$obs)
  
  dat.q <- list(T = 243,
                  T_pred = 0,
                  y = data$obs,
                  q = data$q,
                  tau = data$tau,
                  rov = rov,
                  control = control_q,
                  pop = pop)
  
  sim_q <- stan(file = "Code/seirs.stan", 
              data = dat.q, 
              init = means, 
              iter = 200, 
              chains = 1,
              warmup = 0,
              algorithm = "Fixed_param")
  
  effect_q[i] <- rstan::extract(sim_q, "y_hat", permute = F)[, 1, ] %>% rowSums() %>% mean() - sum(data$obs)
  
  dat.2 <- list(T = 243,
                  T_pred = 0,
                  y = data$obs,
                  q = data$q,
                  tau = data$tau,
                  rov = rov,
                  control = control_2,
                  pop = pop)
  
  sim_2 <- stan(file = "Code/seirs.stan", 
                  data = dat.2, 
                  init = means, 
                  iter = 200, 
                  chains = 1,
                  warmup = 0,
                  algorithm = "Fixed_param")
  
  effect_2[i] <- rstan::extract(sim_2, "y_hat", permute = F)[, 1, ] %>% rowSums() %>% mean() - sum(data$obs)
  
}

plot(strength, effect_q, type = "l")
lines(strength, effect_psi, col = "blue")
lines(strength, effect_2, col = "red")
abline(h = 0)

#===============================================================================
# Control two weeks a year


for(i in 1:51){
  for(j in (i+1):52){
    
    control <- matrix(1, nrow = 243, ncol = 3)
    control[(data$week == i & data$year > 2008 & data$year < 2012), 1] <- 1.05
    control[(data$week == j & data$year > 2008 & data$year < 2012), 1] <- 1.05
    
    dat.stan <- list(T = 243,
                     T_pred = 0,
                     y = data$obs,
                     q = data$q,
                     tau = data$tau,
                     rov = rov,
                     control = control,
                     pop = pop)
    
    sim <- stan(file = "Code/seirs.stan", 
                data = dat.stan, 
                init = means, 
                iter = 100, 
                chains = 1,
                warmup = 0,
                algorithm = "Fixed_param")
    
    totals[i, j] <- rstan::extract(sim, "y_hat", permute = F)[, 1, ] %>% rowSums %>% mean()
  }
}

which.max(totals[])