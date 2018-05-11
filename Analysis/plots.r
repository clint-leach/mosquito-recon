# Code to generate plots for manuscript

library(rstan)
library(plyr)
library(reshape)
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
              year = unique(Year),
              epiweek = 0,
              epiyear = 0)

# Generating epidemic weeks/years that count from peak to peak (week 16)
post$epiweek[16:243] <- c(1:52)
post$epiyear <- cumsum(post$epiweek == 1)

# Population size
pop <- 327801

# Loading and processing the weather data
weather <- read.csv("Data/Vitoria.weather.csv") %>% 
  mutate(date = ymd(BRST), year = year(date), week = week(date))

weather <- subset(weather, date < date[1] + weeks(243))
weather$tot.week <- rep(c(1:243), each = 7)

temp <- ddply(weather, .(tot.week), summarise, temp = mean(Mean.TemperatureC, na.rm = T))

post$temp <- temp$temp
post$rov <- 7 * exp(0.2 * temp$temp - 8)

# Loading MCMC results
sim <- readRDS("Results/euler7.rds")

# Figure 1: Observed and estimated time series =================================

# Computing posterior median and 80% credible interval
yhat <- rstan::extract(sim, c("y_hat"), permute = F) %>% apply(3, quantile, c(0.1, 0.5, 0.9))
qhat <- rstan::extract(sim, c("q_hat"), permute = F) %>% apply(3, quantile, c(0.1, 0.5, 0.9))

# Assembling data frame for ggplot
post$yhat <- yhat[2, ]
post$ymin <- yhat[1, ]
post$ymax <- yhat[3, ]
post$qhat <- qhat[2, ]
post$qmin <- qhat[1, ]
post$qmax <- qhat[3, ]

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

# Figure 2: Estimates of latent mosquito mortality rate ========================

system <- rstan::extract(sim, "system", permute = T)[[1]]

dv <- system[, , 9] %>% + 0.39 %>% exp() %>% 
  adply(2, quantile, c(0.1, 0.5, 0.9), .id = NULL)
names(dv) <- c("dvmin", "dvmed", "dvmax")

post <- cbind(post, dv)

# Plotting

ggplot(post, aes(tot.week, dvmed)) + 
  geom_ribbon(aes(ymin = dvmin, ymax = dvmax), fill = "grey70") +
  geom_line() + 
  theme_classic() +
  scale_y_continuous(expand = c(0.05, 0), limits = c(0.6, 2.5)) + 
  scale_x_continuous(expand = c(0, 1)) +
  ylab("mosquito mortality rate") +
  xlab("week") +
  ggtitle("A") -> fig2.a

ggplot(post, aes(temp, dvmed)) +
  geom_point() + 
  theme_classic() + 
  scale_y_continuous(expand = c(0.05, 0), limits = c(0.6, 2.5)) + 
  scale_x_continuous(expand = c(0, 1)) +
  ylab("mosquito mortality rate") +
  xlab("weekly mean temperature") + 
  ggtitle("B") -> fig2.b

postscript("Manuscript/figures/fig2.eps",
           width = 5.2, height = 3,
           family = "ArialMT")

grid.arrange(fig2.a, fig2.b, ncol = 2)

dev.off()


# Figure 3: Effect of adult control implemented in different weeks ===================

adult <- readRDS("Results/adult_control.rds")

adult_best <- adult %>% 
  ddply(.(rep, year), summarise, 
        best = which.min(reduction),
        prevented =  1 - min(reduction))

delta <- 97 * rstan::extract(sim, "delta_c", permute = T)[[1]][1:1000]
adult_best$delta <- rep(delta, each = 3)

adult_best %>% 
  ggplot(aes(delta, prevented)) + 
  geom_point(size = 0.5) + 
  facet_grid(.~year) + 
  theme_classic() + 
  theme(strip.background = element_blank(), strip.text = element_blank()) +
  xlab("") +
  ylab("proportion of cases prevented") -> fig3.a

adult_best %>% 
  ggplot(aes(delta, best)) + 
  geom_point(size = 0.5) + 
  facet_grid(.~year) + 
  theme_classic() + 
  theme(strip.background = element_blank(), strip.text = element_blank()) +
  xlab("period of cross immunity (weeks)") +
  ylab("week of most effective control") -> fig3.b

# Plotting

postscript("Manuscript/figures/fig3.eps",
           width = 8 , height = 6,
           family = "ArialMT")

fig3.a + fig3.b + plot_layout(ncol = 1)

dev.off()

# Figure 4: Timing of control ==================================================

adult_bests <- ddply(adult_best, .(year), summarise, best = median(best))
adult_bests$epiyear <- c(1, 2, 3)
adult_bests$week <- 0

for(i in 1:3){
  adult_bests$week[i] <- post$tot.week[post$epiyear == adult_bests$epiyear[i] & post$epiweek == adult_bests$best[i]]
}

post <- mutate(post, 
               ymu = system[, , 8] %>% 
                 multiply_by(pop) %>% 
                 apply(2, median),
               qmu = system[, , 6] %>% 
                 apply(2, median),
               Smu = system[, , 1] %>% 
                 apply(2, median))

stacked <- melt(post, id.vars = c("tot.week", "epiweek", "year"), measure.vars = c("ymu", "qmu", "Smu"))

levels(stacked$variable) <- c("cases", "mosquitoes", "susceptibles")

postscript("Manuscript/figures/fig4.eps",
           width = 4, height = 6,
           family = "ArialMT")

subset(stacked, year < 2012) %>% 
  ggplot(aes(tot.week, value)) + 
  geom_line() + 
  geom_vline(xintercept = adult_bests$week, color = "grey50") +
  facet_grid(variable ~ ., scales = "free_y", switch = "y") + 
  theme_classic() +
  scale_color_brewer(palette = "Dark2", type = "qual", guide = F) + 
  scale_x_continuous(expand = c(0, 1), breaks = c(54, 106, 158), labels = c(2009, 2010, 2011)) +
  theme(strip.placement = "outside", strip.background = element_blank()) +
  ylab("") + 
  xlab("year")

dev.off()  

#===============================================================================
# Supplemental figures
#===============================================================================

# Figure S1: epsilons ==========================================================

eps_dv <- rstan::extract(sim, "eps_dv", permute = T)[[1]] %>% 
  adply(2, quantile, c(0.1, 0.5, 0.9)) %>% 
  mutate(week = as.numeric(X1))
names(eps_dv) <- c("X1", "min", "med", "max", "week")

eps_rv <- rstan::extract(sim, "eps_rv", permute = T)[[1]] %>% 
  adply(2, quantile, c(0.1, 0.5, 0.9)) %>% 
  mutate(week = as.numeric(X1))
names(eps_rv) <- c("X1", "min", "med", "max", "week")

postscript("Manuscript/figures/figS1.eps",
           width = 5.2, height = 3,
           family = "ArialMT")

ggplot(eps_dv, aes(week, med)) + 
  geom_ribbon(aes(ymin = min, ymax = max), fill = "grey70") +
  geom_line() + 
  theme_classic() + 
  scale_x_continuous(expand = c(0, 1)) +
  ylab(expression(epsilon[d])) +
  ggtitle("A") -> figs1.a

ggplot(eps_rv, aes(week, med)) + 
  geom_ribbon(aes(ymin = min, ymax = max), fill = "grey70") +
  geom_line() + 
  theme_classic() + 
  scale_x_continuous(expand = c(0, 1)) +
  ylab(expression(epsilon[r])) +
  ggtitle("B") -> figs2.b

grid.arrange(figs1.a, figs2.b, nrow = 1)

dev.off()

# Figure S2: sigmas ============================================================

postscript("Manuscript/figures/figS2.eps",
           width = 6, height = 4,
           family = "ArialMT")

par(mfrow = c(1, 2))
rstan::extract(sim, "sigmadv", permute = T)[[1]] %>% hist(main = "", xlab = expression(sigma[d]), freq = F)
rstan::extract(sim, "sigmarv", permute = T)[[1]] %>% hist(main = "", xlab = expression(sigma[r]), freq = F)

dev.off()

# Figure S3: epidemiological parameters ===================================================

epiparams <- rstan::extract(sim, c("ro_c", "gamma", "delta_c"), permute = T)

postscript("Manuscript/figures/figS3.eps",
           width = 6, height = 3,
           family = "ArialMT")

par(mfrow = c(1, 3))

hist(0.87 * epiparams[[1]], main = "", xlab = "latent period", freq = F)
abline(v = 0.87, lwd = 2)

hist(epiparams[[2]], main = "", xlab = "rate of infectious decay", freq = F)
abline(v = 3.5, lwd = 2)

hist(97 * epiparams[[3]], main = "", xlab = "period of cross-immunity", freq = F)
abline(v = 97, lwd = 2)

dev.off()

# Figure S4: initial conditions ===========================================================

ics <- rstan::extract(sim, c("S0", "E0","I0"), permute = T)

postscript("Manuscript/figures/figS4.eps",
           width = 6, height = 3,
           family = "ArialMT")

par(mfrow = c(1, 3))

hist(ics$S0, freq = F, main = "", xlab = expression(S[0]))
hist(ics$E0, freq = F, main = "", xlab = expression(E[0]))
hist(ics$I0, freq = F, main = "", xlab = expression(I[0]))

dev.off()


# Model checking ===============================================================

yrep <- rstan::extract(sim, "y_hat", permute = T)[[1]]
qrep <- rstan::extract(sim, "q_hat", permute = T)[[1]]

# Figure S5: ACF plots =========================================================

# Calculating the observed autocorrelation function
yacf <- acf(post$yobs, lag.max = 55, plot = F)$acf
qacf <- acf(post$qobs, lag.max = 55, plot = F)$acf

# Summarizing the posterior acf
yacfrange <- apply(yrep, 1, function(x){acf(x, lag.max = 55, plot = F)$acf}) %>% 
  apply(1, quantile, probs = c(0.1,0.5, 0.9))

qacfrange <- apply(qrep, 1, function(x){acf(x, lag.max = 55, plot = F)$acf}) %>% 
  apply(1, quantile, probs = c(0.1,0.5, 0.9))

# Assembling data frame for ggplot
acfdf <- data.frame("lag" = c(0:55), 
                    "yobs" = yacf,
                    "qobs" = qacf,
                    "ymed" = yacfrange[2, ], 
                    "ymin" = yacfrange[1, ], 
                    "ymax" = yacfrange[3, ],
                    "qmed" = qacfrange[2, ], 
                    "qmin" = qacfrange[1, ], 
                    "qmax" = qacfrange[3, ]
                    )

# Plots
ggplot(acfdf, aes(lag, yobs)) + 
  geom_hline(yintercept = 0, color = "gray50") +
  geom_linerange(aes(lag, ymin = ymin, ymax = ymax)) +
  geom_point() +
  theme_classic() +
  xlab("lag (weeks)") +
  ylab("case autocorrelation") + 
  scale_x_continuous(expand = c(0, 0.1)) + 
  ggtitle("A") -> figs5.a

ggplot(acfdf, aes(lag, qobs)) + 
  geom_hline(yintercept = 0, color = "gray50") +
  geom_linerange(aes(lag, ymin = qmin, ymax = qmax)) +
  geom_point() +
  theme_classic() +
  xlab("lag (weeks)") +
  ylab("mosquito autocorrelation") + 
  scale_x_continuous(expand = c(0, 0.1)) + 
  ggtitle("B") -> figs5.b

postscript("Manuscript/figures/figS5.eps",
           width = 5.2, height = 3,
           family = "ArialMT")

grid.arrange(figs5.a, figs5.b, nrow = 1)

dev.off()

# Figure S6: Totals ============================================================

postscript("Manuscript/figures/figS6.eps",
           width = 6, height = 4,
           family = "ArialMT")

par(mfrow = c(1, 2))

hist(rowSums(yrep), main = "", xlab = "total cases", freq = F, breaks = 20)
abline(v = sum(post$yobs), lwd = 2)

hist(rowSums(qrep), main = "", xlab = "total captured mosquitoes", freq = F, breaks = 20)
abline(v = sum(post$qobs), lwd = 2)

dev.off()

# Figure S7: Min and max ===========================================================

postscript("Manuscript/figures/figS7.eps",
           width = 6, height = 6,
           family = "ArialMT")

par(mfrow = c(2, 2))

hist(apply(yrep, 1, max), main = "", xlab = "maximum weekly cases", freq = F)
abline(v = max(post$yobs), lwd = 2)

hist(apply(qrep, 1, max), main = "", xlab = "maximum weekly trap count", freq = F)
abline(v = max(post$qobs), lwd = 2)

hist(apply(yrep, 1, min), main = "", xlab = "minimum weekly cases", freq = F)
abline(v = min(post$yobs), lwd = 2)

hist(apply(qrep, 1, min), main = "", xlab = "minimum weekly trap count", freq = F)
abline(v = min(post$qobs), lwd = 2)

dev.off()

# Figure S8: Larval control ====================================================

# Processing adult control experiment
larval_reduction <- readRDS("Results/larval_control.rds") %>% 
  adply(1, quantile, probs = c(0.1, 0.5, 0.9)) %>% 
  mutate(week = as.numeric(X1))

names(larval_reduction) <- c("X1", "min", "med", "max", "week")

postscript("Manuscript/figures/figS8.eps",
           width = 5.2, height = 3,
           family = "ArialMT")

ggplot(larval_reduction, aes(week, med)) + 
  geom_ribbon(aes(ymin = min, ymax = max), fill = "grey70") + 
  geom_line() +
  geom_abline(slope = 0, color = "grey20") + 
  theme_classic() + 
  scale_x_continuous(expand = c(0, 1)) +
  xlab("week of control") +
  ylab("cases prevented")

dev.off()


#===============================================================================
ggplot(post, aes(tot.week, rvmed)) + 
  geom_ribbon(aes(ymin = rvmin, ymax = rvmax), fill = "grey70") +
  geom_line() + 
  theme_classic() +
  scale_y_continuous(expand = c(0, 0)) + 
  scale_x_continuous(expand = c(0, 1)) +
  ylab("mosquito growth rate") + 
  xlab("week") +
  ggtitle("A") -> fig2.a