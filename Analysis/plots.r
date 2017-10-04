# Code to generate plots for manuscript

library(rstan)
library(plyr)
library(magrittr)
library(lubridate)
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

temp <- ddply(weather, .(tot.week), summarise, temp = mean(Mean.TemperatureC, na.rm = T))

post$rov <- 7 * exp(0.2 * temp$temp - 8)

# Loading MCMC results
sim <- readRDS("Results/oscillator.rds")

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
# Figure 2: Estimates of latent mosquito demographic rates

system <- rstan::extract(sim, "system", permute = T)[[1]]

dv <- system[, , 9] %>% + 0.39 %>% exp() %>% 
  adply(2, quantile, c(0.1, 0.5, 0.9), .id = NULL)
names(dv) <- c("dvmin", "dvmed", "dvmax")

rv <- system[, , 11] %>% apply(2, quantile, c(0.1, 0.5, 0.9), .id = NULL)
names(rv) <- c("rvmin", "rvmed", "rvmax")

post <- cbind(post, dv, rv)

#===============================================================================
# Figure 3: Effect of adult control implemented in different weeks

reduction <- readRDS("Results/control.rds") %>% 
  adply(1, quantile, probs = c(0.1, 0.5, 0.9)) %>% 
  mutate(week = as.numeric(X1))

names(reduction) <- c("X1", "min", "med", "max", "week")

ggplot(reduction, aes(week, med)) + 
  geom_ribbon(aes(ymin = min, ymax = max), fill = "grey70") + 
  geom_line() +
  geom_abline(slope = 0, color = "grey20") + 
  theme_classic() + 
  theme(text = element_text(size = 15)) + 
  scale_x_continuous(expand = c(0, 1)) +
  xlab("week of control") +
  ylab("cases prevented")
