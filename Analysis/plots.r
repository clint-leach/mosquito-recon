# Code to generate plots for manuscript

library(rstan)
library(plyr)
library(magrittr)
library(ggplot2)

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
# Figure 1

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
# Figure 2

mcmc <- rstan::extract(sim, "q_hat", permute = T)[[1]]
nmcmc <- dim(mcmc)[1]

for(i in 1:243){
  post$p[i] <- sum(mcmc[, i] > post$qobs[i]) / nmcmc
}

ggplot(post, aes(tot.week, p)) +
  geom_line() + 
  theme_classic() + 
  xlab("week") + 
  ylab("p")

#===============================================================================
# Figure 3