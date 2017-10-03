library(rstan)
library(plyr)
library(magrittr)
library(lubridate)
library(doParallel)

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

weather <- subset(weather, date < date[1] + weeks(243))
weather$tot.week <- rep(c(1:243), each = 7)

temp <- ddply(weather, .(tot.week), summarise, temp = mean(Mean.TemperatureC, na.rm = T))

rov <- 7 * exp(0.2 * temp$temp - 8)
#===============================================================================
# Fitting to first four years

T_pred <- 0
T_fit <- 191

dat.stan <- list(T = T_fit,
                 T_pred = T_pred,
                 y = head(obs, T_fit),
                 q = head(q, T_fit),
                 tau = head(tau, T_fit),
                 rov = head(rov, T_fit),
                 control = matrix(1, nrow = 191, ncol = 3),
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
                  dv0 = 0,
                  rv0 = 0,
                  sigmadv = 0.2,
                  sigmarv = 0.2,
                  eps_dv = rep(0, T_fit),
                  eps_rv = rep(0, T_fit)),
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
                  dv0 = 0,
                  rv0 = 0,
                  sigmadv = 0.5,
                  sigmarv = 0.5,
                  eps_dv = rep(0, T_fit),
                  eps_rv = rep(0, T_fit)),
             list(S0 = 0.5,
                  E0 = 40,
                  I0 = 20,
                  log_phi_q = -12.5,
                  eta_inv_y = 5,
                  eta_inv_q = 5,
                  ro_c = 1,
                  gamma_c = 1,
                  delta_c = 1,
                  logNv = 0.7,
                  dv0 = 0,
                  rv0 = 0,
                  sigmadv = 0.01,
                  sigmarv = 0.01,
                  eps_dv = rep(0, T_fit),
                  eps_rv = rep(0, T_fit)))

fit <- stan(file = "Code/seirs.stan", 
            data = dat.stan, 
            init = inits, 
            iter = 2000, 
            chains = 3,
            control = list(adapt_delta = 0.99,
                           max_treedepth = 15))

saveRDS(fit, "Results/pred.rds")

#===============================================================================

fit <- readRDS("Results/pred.rds")

samples <- rstan::extract(fit)[1:23]

nmcmc <- length(samples[[1]])

#===============================================================================
# Setting up parallel

cl <- makeCluster(10, type = "SOCK")
registerDoParallel(cl)

#===============================================================================

extract_sample <- function(x, k){
  if(length(dim(x)) == 1) return(x[k])
  else return(x[k, ])
}

model <- stan_model(file = "Code/seirs.stan")

dat.stan <- list(T = T_fit,
                 T_pred = 52,
                 y = head(obs, T_fit),
                 q = head(q, T_fit),
                 tau = head(tau, T_fit),
                 rov = rov,
                 control = matrix(1, nrow = 243, ncol = 3),
                 pop = pop)

cases <- foreach(k = 1:nmcmc, .combine = "rbind", .packages = c("rstan", "magrittr")) %dopar% {
  
  init <- list(lapply(samples, extract_sample, k))
  
  sim <- sampling(model,
                  data = dat.stan, 
                  init = init, 
                  iter = 1, 
                  chains = 1,
                  warmup = 0,
                  algorithm = "Fixed_param")
  
  return(rstan::extract(sim, "y_hat", permute = F)[, 1, ])
}

foo <- colMeans(cases)
plot(obs)
lines(foo)

total <- apply(cases, 1, function(x) sum(x[191:243]))
hist(log10(total))
abline(v = log10(sum(obs[191:243])))
