# Code to generate plots for manuscript

library(rstan)
library(plyr)
library(magrittr)
library(ggplot2)
library(gridExtra)

# Loading and processing mosquito and case data
tseries <- read.csv("Data/Vitoria.data.csv")

post <- ddply(tseries, .(tot.week), summarise,
              qobs = sum(Mosquitoes, na.rm = T),
              yobs = sum(Cases, na.rm = T),
              tau = sum(Trap, na.rm = T),
              week = unique(Week),
              year = unique(Year))

# Loading and processing the weather data
weather <- read.csv("Data/Vitoria.weather.csv") %>% 
  mutate(date = ymd(BRST), year = year(date), week = week(date))

weather <- subset(weather, date < date[1] + weeks(243))
weather$tot.week <- rep(c(1:243), each = 7)

covars <- ddply(weather, .(tot.week), summarise,
                temp = mean(Mean.TemperatureC, na.rm = T),
                humidity = mean(Mean.Humidity, na.rm = T))

# Loading MCMC results
sim <- readRDS("Results/seasonal.rds")

#===============================================================================
# Figure 1: Observed and estimated time series

# Computing posterior median and 80% credible interval
yhat <- rstan::extract(sim, c("y_hat"), permute = F) %>% apply(3, quantile, c(0.1, 0.5, 0.9))
qhat <- rstan::extract(sim, c("q_hat"), permute = F) %>% apply(3, quantile, c(0.1, 0.5, 0.9))

# Assembling data frame for ggplot
post$yhat = yhat[2, ]
post$ymin = yhat[1, ]
post$ymax = yhat[3, ]
post$qhat = qhat[2, ]
post$qmin = qhat[1, ]
post$qmax = qhat[3, ]

# Plotting
ggplot(post, aes(tot.week, yobs)) + 
  geom_ribbon(aes(ymin = ymin, ymax = ymax), fill = "grey70") +
  geom_point(size = 0.5) + 
  geom_line(aes(tot.week, yhat)) + 
  theme_classic() +
  scale_y_continuous(expand = c(0, 0)) + 
  scale_x_continuous(expand = c(0, 1)) +
  ylab("case reports") + 
  xlab("week") +
  ggtitle("A") -> fig1.a

ggplot(post, aes(tot.week, qobs)) + 
  geom_ribbon(aes(ymin = qmin, ymax = qmax), fill = "grey70") +
  geom_point(size = 0.5) + 
  geom_line(aes(tot.week, qhat)) + 
  theme_classic() +
  scale_y_continuous(expand = c(0.05, 0)) + 
  scale_x_continuous(expand = c(0, 1)) +
  ylab("mosquitoes trapped") +
  xlab("week") +
  ggtitle("B") -> fig1.b

postscript("Manuscript/figures/fig1.eps",
     width = 5.2, height = 3,
     family = "ArialMT")

grid.arrange(fig1.a, fig1.b, ncol = 2)

dev.off()

#===============================================================================
# Figure 2: estimated demographic rates

b <- rstan::extract(sim, "bv", permute = T)[[1]] %>% apply(2, quantile, c(0.1, 0.5, 0.9))
d <- rstan::extract(sim, "dv", permute = T)[[1]] %>% apply(2, quantile, c(0.1, 0.5, 0.9))

post[, c("bmin", "b", "bmax")] <- t(b)
post[, c("dmin", "d", "dmax")] <- t(d)

ggplot(post, aes(tot.week, b)) + 
  geom_ribbon(aes(ymin = bmin, ymax = bmax), fill = "grey70") +
  geom_line() + 
  theme_classic() +
  scale_y_continuous(expand = c(0.05, 0)) + 
  scale_x_continuous(expand = c(0, 1)) +
  ylab("mosquito net emergence rate") +
  xlab("week") +
  ggtitle("A") -> fig2.a

ggplot(post, aes(tot.week, d)) + 
  geom_ribbon(aes(ymin = dmin, ymax = dmax), fill = "grey70") +
  geom_line() + 
  theme_classic() +
  scale_y_continuous(expand = c(0.05, 0)) + 
  scale_x_continuous(expand = c(0, 1)) +
  ylab("mosquito per-capita death rate") +
  xlab("week") +
  ggtitle("B") -> fig2.b

postscript("Manuscript/figures/fig2.eps",
           width = 5.2, height = 3,
           family = "ArialMT")

grid.arrange(fig2.a, fig2.b, ncol = 2)

dev.off()

betas <- rstan::extract(sim, c("beta0", "beta"), permute = T)
betas <- cbind(betas[[1]], betas[[2]])

X <- cbind(1, scale(covars[, -1]))
mud <- exp(X %*% t(betas)) %>% rowMeans()
res <- mud - d[2, ]

#===============================================================================
# Figure 3: Predicting 2012 outbreak

predsim <- readRDS("Results/covar_predict.rds")
post$fit <- c(rep(1, 243 - 34), rep(0, 34))

# Plotting the estimated and predicted cases
ypred <- rstan::extract(predsim, c("y_hat"), permute = F) %>% apply(3, quantile, c(0.1, 0.5, 0.9))

post[, c("predmin", "pred", "predmax")] <- t(ypred)

ggplot(post, aes(tot.week, yobs)) + 
  geom_ribbon(aes(ymin = predmin, ymax = predmax, fill = as.factor(fit))) +
  scale_fill_manual(values = c("0" = "grey90", "1" = "grey70"), guide = F) + 
  geom_point(size = 0.5) + 
  geom_line(aes(tot.week, pred)) + 
  theme_classic() +
  scale_y_continuous(expand = c(0, 0)) + 
  scale_x_continuous(expand = c(0, 1)) +
  ylab("case reports") + 
  xlab("week") +
  ggtitle("A")

# Plotting the estimated and predicted death rates
dpred <- rstan::extract(predsim, "d_pred", permute = F) %>% adply(3, quantile, .id = NULL, c(0.1, 0.5, 0.9))
names(dpred) <- c("predmin", "pred", "predmax")

dpred$week <- c((T_fit + 1):243)
dpred[, c("fitmin", "fit", "fitmax")] <- t(d[, (T_fit + 1):243])

ggplot(dpred, aes(week, pred)) + 
  geom_line(linetype = 2) + 
  geom_line(aes(week, fit)) + 
  theme_classic() +
  scale_y_continuous(expand = c(0, 0)) + 
  scale_x_continuous(expand = c(0, 1)) +
  ylab("mosquito per capita death rate") + 
  xlab("week") +
  ggtitle("B")

#===============================================================================
# Figure 3: observed and posterior autocorrelation functions

# Calculating the observed autocorrelation function
obsacf <- acf(post$qobs, lag.max = 55, plot = F)$acf

# Calculating the posterior acf function
acfs <- apply(mcmc, 1, function(x){acf(x, lag.max = 55, plot = F)$acf})

# Median and 80% credible interval
acfrange <- apply(acfs, 1, quantile, c(0.1, 0.5, 0.9))

# Assembling data frame for ggplot
acfdf <- data.frame("lag" = c(0:55), 
                    "obs" = obsacf, 
                    "med" = acfrange[2, ], 
                    "min" = acfrange[1, ], 
                    "max" = acfrange[3, ])

# Plot
postscript("Manuscript/figures/fig3.eps",
           width = 4, height = 3,
           family = "ArialMT")

ggplot(acfdf, aes(lag, obs)) + 
  geom_hline(yintercept = 0, color = "gray50") +
  geom_linerange(aes(lag, ymin = min, ymax = max)) +
  geom_point() +
  theme_classic() +
  xlab("lag (weeks)") +
  ylab("autocorrelation") + 
  scale_x_continuous(expand = c(0, 0.1))

dev.off()

