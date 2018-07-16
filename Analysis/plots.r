# Code to generate plots for manuscript

library(rstan)
library(plyr)
library(reshape)
library(magrittr)
library(lubridate)
library(ggplot2)
library(gridExtra)
library(patchwork)
library(truncnorm)

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

temp <- ddply(weather, .(tot.week), summarise, temp = mean(Mean.TemperatureC, na.rm = T),
              humid = mean(Mean.Humidity, na.rm = T))

post$temp <- temp$temp
post$rov <- 7 * exp(0.2 * temp$temp - 8)

# Loading MCMC results
sim <- readRDS("Results/gamma_eip_dv0.rds")

# Figure 1: Observed and estimated time series =================================

# Computing posterior median and 80% credible interval
yhat <- rstan::extract(sim, c("y_meas"), permute = F) %>% apply(3, quantile, c(0.1, 0.5, 0.9))
qhat <- rstan::extract(sim, c("q_meas"), permute = F) %>% apply(3, quantile, c(0.1, 0.5, 0.9))

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

dvmu <- 1.47 * rstan::extract(sim, "dvmu_c", permute = T)[[1]]
system <- rstan::extract(sim, "state", permute = T)[[1]]

dv <- system[, 2:244, 12] %>% exp() %>% t() %>% multiply_by_matrix(diag(dvmu)) %>% 
  adply(1, quantile, c(0.1, 0.5, 0.9), .id = NULL)
names(dv) <- c("dvmin", "dvmed", "dvmax")

post <- cbind(post, dv)

# Plotting

ggplot(post, aes(tot.week, dvmed)) + 
  geom_ribbon(aes(ymin = dvmin, ymax = dvmax), fill = "grey70") +
  geom_line() + 
  theme_classic() +
  scale_y_continuous(expand = c(0.05, 0), limits = c(0.25, 1.35)) + 
  scale_x_continuous(expand = c(0, 1)) +
  ylab("mosquito mortality rate") +
  xlab("week") +
  ggtitle("A") -> fig2.a

ggplot(post, aes(temp, dvmed)) +
  geom_point() + 
  theme_classic() + 
  scale_y_continuous(expand = c(0.05, 0), limits = c(0.25, 1.35)) + 
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

# First day of each year
firsts <- daply(post, .(year), function(df) min(df$tot.week))

# Extracting delta values
delta <- 97 * rstan::extract(sim, "delta_c", permute = T)[[1]]

# Adult control simulation results
adult <- readRDS("Results/gamma_control.rds")
adult <- mutate(adult, 
                prevented = 1 - reduction, 
                relweek = control - firsts[as.character(year)],
                delta = rep(delta, each = 840))

# Summarizing adult control simulations
adult_sum <- adult %>%
  ddply(.(year, relweek), summarise,
        med = median(prevented),
        high = quantile(prevented, 0.9),
        low = quantile(prevented, 0.1),
        absweek = unique(control))

# Plotting

postscript("Manuscript/figures/fig3.eps",
           width = 4, height = 5,
           family = "ArialMT")

subset(adult_sum, relweek > -53 & relweek < 53) %>% 
  ggplot(aes(relweek, med)) +
  geom_ribbon(aes(ymin = low, ymax = high), fill = "grey70") + 
  geom_line() + 
  facet_grid(year ~.) + 
  theme_classic() +
  theme(strip.background = element_blank()) + 
  scale_x_continuous(expand = c(0, 0.5)) + 
  xlab("relative week of control") + 
  ylab("proportion of cases prevented")

dev.off()

# Figure 4: Timing of control ==================================================

# Finding week when control is most effective for each rep
adult_best <- adult %>%
  subset(relweek > -53 & relweek < 53) %>% 
  ddply(.(rep, year), summarise,
        best = relweek[which.min(reduction)],
        ratio = max(prevented))

adult_best$delta <- rep(delta, each = 4)

# Finding median most effective week of control for each year
bests <- ddply(adult_best, .(year), summarise, best = median(best), prevented = median(ratio))
bests$week <- bests$best + firsts[2:5]

# Summarizing state variables
post <- mutate(post,
               ymu = system[, ,11] %>%
                 apply(1, diff) %>%
                 multiply_by(pop) %>%
                 apply(1, median),
               qmu = system[, 2:244, 9] %>%
                 apply(2, median),
               Smu = system[, 2:244, 1] %>%
                 apply(2, median))

# Stacking state variables
stacked <- melt(post, id.vars = c("tot.week", "epiweek", "year"), measure.vars = c("ymu", "Smu", "qmu"))
levels(stacked$variable) <- c("cases",  "susceptibles", "mosquitoes")

postscript("Manuscript/figures/fig4.eps",
           width = 4, height = 4,
           family = "ArialMT")

subset(stacked, year < 2012) %>%
  ggplot(aes(tot.week, value)) +
  geom_line() +
  geom_vline(xintercept = bests$week, color = "grey50") +
  facet_grid(variable ~ ., scales = "free_y", switch = "y") +
  theme_classic() +
  scale_color_brewer(palette = "Dark2", type = "qual", guide = F) +
  scale_x_continuous(expand = c(0, 5), breaks = c(1, 54, 106, 158, 210), labels = c(2008, 2009, 2010, 2011, 2012)) +
  theme(strip.placement = "outside", strip.background = element_blank()) +
  ylab("") +
  xlab("year")

dev.off()

# Figure 5: effect of delta on control =========================================

postscript("Manuscript/figures/fig5.eps",
           width = 6, height = 4,
           family = "ArialMT")

adult_best %>%
  ggplot(aes(delta, ratio)) +
  geom_bin2d(binwidth = c(5, 0.01)) +
  scale_fill_distiller(palette = "Oranges", direction = 1, guide = F) +
  # scale_fill_gradient(low = "#E6E6E6", high = "#000000", guide = FALSE) +
  facet_grid(.~year) +
  theme_classic() +
  theme(strip.background = element_blank()) +
  xlab("period of cross-immunity (weeks)") +
  ylab("maximum proportion of cases prevented")

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
  geom_hline(yintercept = 0, linetype = 2) + 
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

epiparams <- rstan::extract(sim, c("ro_c", "gamma_c", "delta_c", "dvmu_c"), permute = T)

postscript("Manuscript/figures/figS3.eps",
           width = 6, height = 6,
           family = "ArialMT")

par(mfrow = c(2, 2))

hist(epiparams[[1]], main = "", xlab = "scaled latent period", freq = F, breaks = 30, ylim = c(0, 1.2))
curve(dgamma(x, 8.3, 8.3), add = T, lwd = 2)

hist(epiparams[[2]], main = "", xlab = "scaled rate of infectious decay", freq = F, breaks = 30)
curve(dgamma(x, 100, 100), add = T, lwd = 2)

hist(epiparams[[3]], main = "", xlab = "scaled period of cross-immunity", freq = F, breaks = 30, ylim = c(0, 1.3))
curve(dgamma(x, 10, 10), add = T, lwd = 2)

hist(epiparams[[4]], main = "", xlab = "scaled mosquito mortality rate", freq = F, breaks = 30)
curve(dgamma(x, 100, 10), add = T, lwd = 2)

dev.off()

# Figure S4: initial conditions ===========================================================

ics <- rstan::extract(sim, c("S0", "E0","I0", "logNv", "dv0", "rv0"), permute = T)

postscript("Manuscript/figures/figS4.eps",
           width = 6, height = 5,
           family = "ArialMT")

par(mfrow = c(2, 3))

hist(ics$S0, freq = F, main = "", xlab = expression(S[0]))
curve(dbeta(x, 4, 6), add = T, lwd = 2)

hist(ics$E0, freq = F, main = "", xlab = expression(E[0]), ylim = c(0, 0.015))
curve(dgamma(x, 10, 0.1), add = T, lwd = 2)

hist(ics$I0, freq = F, main = "", xlab = expression(I[0]), ylim = c(0, 0.02))
curve(dgamma(x, 6, 0.1), add = T, lwd = 2)

hist(ics$logNv, freq = F, main = "", xlab = expression(log(V[N0])))
curve(dnorm(x, 0.7, 0.3), add = T, lwd = 2)

hist(ics$dv0, freq = F, main = "", xlab = expression(nu[0]))
curve(dnorm(x, 0, 0.5), add = T, lwd = 2)

hist(ics$rv0, freq = F, main = "", xlab = expression(r[0]))
curve(dnorm(x, 0, 0.5), add = T, lwd = 2)

dev.off()

# Figure S5: Measurement parameters ============================================

meas <- rstan::extract(sim, c("log_phi_q", "eta_inv_q", "eta_inv_y"), permute = T)

postscript("Manuscript/figures/figS5.eps",
           width = 6, height = 3,
           family = "ArialMT")

par(mfrow = c(1, 3))

hist(meas$log_phi_q, freq = F, main = "", xlab = expression(log(phi[q])))
curve(dnorm(x, -13, 0.5), add = T, lwd = 2)

hist(meas$eta_inv_q, freq = F, main = "", xlab = expression(eta[q]))
curve(dtruncnorm(x, a = 0, b = Inf, mean = 0, sd = 5), add = T)

hist(meas$eta_inv_y, freq = F, main = "", xlab = expression(eta[y]))
curve(dtruncnorm(x, a = 0, b = Inf, mean = 0, sd = 5), add = T)

dev.off()

# Model checking ===============================================================

yrep <- rstan::extract(sim, "y_meas", permute = T)[[1]]
qrep <- rstan::extract(sim, "q_meas", permute = T)[[1]]

# Figure S6: ACF plots =========================================================

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

postscript("Manuscript/figures/figS6.eps",
           width = 5.2, height = 3,
           family = "ArialMT")

grid.arrange(figs5.a, figs5.b, nrow = 1)

dev.off()

# Figure S7: Totals ============================================================

postscript("Manuscript/figures/figS7.eps",
           width = 6, height = 4,
           family = "ArialMT")

par(mfrow = c(1, 2))

hist(rowSums(yrep), main = "", xlab = "total cases", freq = F, breaks = 20)
abline(v = sum(post$yobs), lwd = 2)

hist(rowSums(qrep), main = "", xlab = "total captured mosquitoes", freq = F, breaks = 20)
abline(v = sum(post$qobs), lwd = 2)

dev.off()

# Figure S8: Min and max ===========================================================

postscript("Manuscript/figures/figS8.eps",
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
