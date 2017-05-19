# Code to generate plots for manuscript

library(rstan)
library(plyr)
library(magrittr)
library(ggplot2)
library(gridExtra)

# Loading and processing data
tseries <- read.csv("Data/Vitoria.data.csv")

post <- ddply(tseries, .(tot.week), summarise,
              qobs = sum(Mosquitoes, na.rm = T),
              yobs = sum(Cases, na.rm = T),
              tau = sum(Trap, na.rm = T),
              week = unique(Week),
              year = unique(Year))

# Loading MCMC results
sim <- readRDS("Results/mcmc.rds")

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
  ggtitle("A")

ggplot(post, aes(tot.week, d)) + 
  geom_ribbon(aes(ymin = dmin, ymax = dmax), fill = "grey70") +
  geom_line() + 
  theme_classic() +
  scale_y_continuous(expand = c(0.05, 0)) + 
  scale_x_continuous(expand = c(0, 1)) +
  ylab("mosquito per-capita death rate") +
  xlab("week") +
  ggtitle("B")

#===============================================================================
# Figure 2: weekly posterior predictive p-value for replicate mosquito
# surveillance data

# Extracting full MCMC chain of replicate data
mcmc <- rstan::extract(sim, "q_hat", permute = T)[[1]]
nmcmc <- dim(mcmc)[1]

# Computing P(q^{rep}_t > q_t)
for(i in 1:243){
  post$p[i] <- sum(mcmc[, i] > post$qobs[i]) / nmcmc
}

# Plot
postscript("Manuscript/figures/fig2.eps",
           width = 4, height = 3,
           family = "ArialMT")

ggplot(post, aes(week, p, group = factor(year))) + 
  geom_point() + 
  geom_smooth(se = F, color = "gray40") +
  theme_classic() + 
  scale_y_continuous(expand = c(0, 0.01)) + 
  scale_x_continuous(expand = c(0, 1)) +
  xlab("week") + 
  ylab("p")

dev.off()
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

