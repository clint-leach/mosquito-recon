# Code to generate plots for manuscript
library(rstan)
library(plyr)
library(magrittr)
library(lubridate)
library(ggplot2)
library(gridExtra)
library(patchwork)
library(reshape2)

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
weather$control <- rep(c(1:243), each = 7)

temp <- ddply(weather, .(control), summarise, temp = mean(Mean.TemperatureC, na.rm = T),
              humid = mean(Mean.Humidity, na.rm = T))

post$temp <- temp$temp
post$rov <- 7 * exp(0.21 * temp$temp - 7.9)

# Loading MCMC results
sim <- readRDS("Results/chain.rds")

# Effect of control in the number of cases in the following year
moving <- readRDS("Results/moving_control.rds")

# Figure 1: Observed and estimated time series =================================

# Case reports
yhat <- rstan::extract(sim, c("y_meas"), permute = T)[[1]] %>% 
  reshape2::melt(varnames = c("rep", "week"), value.name = "yhat") %>% 
  ddply(.(week), summarise,
        ymed = median(yhat),
        ymin = quantile(yhat, 0.1), 
        ymax = quantile(yhat, 0.9))
  
# Trapped mosquitoes
qhat <- rstan::extract(sim, c("q_meas"), permute = T)[[1]] %>%
  reshape2::melt(varnames = c("rep", "week"), value.name = "qhat") %>% 
  ddply(.(week), summarise,
        qmed = median(qhat),
        qmin = quantile(qhat, 0.1), 
        qmax = quantile(qhat, 0.9))

# Mosquito mortality rate
system <- rstan::extract(sim, "state", permute = T)[[1]]

dvmu <- data.frame(rep = 1:dim(system)[1],
                   dvmu = 1.47 * rstan::extract(sim, "dvmu_c", permute = T)[[1]])

dv <- system[, 2:244, 12] %>% 
  reshape2::melt(varnames = c("rep", "control"), value.name = "dv") %>% 
  mutate(dvnat = exp(dv)) %>% 
  join(dvmu) %>% 
  mutate(dvnat = dvnat * dvmu) %>% 
  join(temp)

# Plotting
yhat %>% 
  ggplot(aes(week, ymin = ymin, ymax = ymax)) +
  geom_ribbon(fill = "grey70") + 
  geom_line(aes(week, ymed)) + 
  geom_point(aes(tot.week, yobs), data = post, inherit.aes = FALSE, size = 0.5) +
  theme_classic() +
  scale_y_continuous(expand = c(0, 0)) + 
  scale_x_continuous(expand = c(0, 1), breaks = c(1, 54, 106, 158, 210), labels = c(2008, 2009, 2010, 2011, 2012)) +
  ylab("case reports") + 
  xlab("year") +
  ggtitle("A") -> fig1a

qhat %>% 
  ggplot(aes(week, ymin = qmin, ymax = qmax)) +
  geom_ribbon(fill = "grey70") + 
  geom_line(aes(week, qmed)) + 
  geom_point(aes(tot.week, qobs), data = post, inherit.aes = FALSE, size = 0.5) +
  theme_classic() +
  scale_y_continuous(expand = c(0, 0)) + 
  scale_x_continuous(expand = c(0, 1), breaks = c(1, 54, 106, 158, 210), labels = c(2008, 2009, 2010, 2011, 2012)) +
  ylab("mosquitoes trapped") + 
  xlab("year") +
  ggtitle("B") -> fig1b

post %>% 
  ggplot(aes(tot.week, 1 / rov)) + 
  geom_line() + 
  theme_classic() +
  scale_y_continuous(expand = c(0.05, 0)) +
  scale_x_continuous(expand = c(0, 1), breaks = c(1, 54, 106, 158, 210), labels = c(2008, 2009, 2010, 2011, 2012)) +
  ylab("EIP (weeks)") +
  xlab("year") +
  ggtitle("C") -> fig1c

dv %>% 
  ggplot(aes(control, dvnat)) + 
  stat_summary(fun.ymin = function(x) quantile(x, 0.05),
               fun.ymax = function(x) quantile(x, 0.95),
               geom = "ribbon", fill = "gray70") +
  stat_summary(fun.y = "median", geom = "line") +
  theme_classic() +
  scale_y_continuous(expand = c(0.05, 0)) +
  scale_x_continuous(expand = c(0, 1), breaks = c(1, 54, 106, 158, 210), labels = c(2008, 2009, 2010, 2011, 2012)) +
  ylab("mosquito mortality rate") +
  xlab("year") +
  ggtitle("D") -> fig1d

postscript("Manuscript/figures/fig1.eps",
     width = 5, height = 7, paper = "special", horizontal = FALSE,
     family = "ArialMT")

fig1a + fig1b + fig1c + fig1d + plot_layout(ncol = 1)

dev.off()

# Figure 2: Estimates of latent mosquito mortality rate ========================

postscript("Manuscript/figures/fig2.eps",
           width = 4, height = 3, paper = "special", horizontal = FALSE,
           family = "ArialMT")

dv %>%
  ddply(.(control), summarise,
        dvnat = median(dvnat),
        temp = median(temp)) %>%
  ggplot(aes(temp, dvnat)) +
  geom_point() +
  theme_classic() +
  scale_y_continuous(expand = c(0.05, 0)) +
  scale_x_continuous(expand = c(0, 1)) +
  ylab("mosquito mortality rate") +
  xlab("weekly mean temperature")

dev.off()

# Figure 3: Effect of adult control implemented in different weeks ===================

# Adding week and year to moving window control results
moving <- moving %>% 
  dplyr::rename("dv0" = "dvmu") %>%
  mutate(week = week(ymd("2008-01-01") + weeks(control - 1)),
         year = year(ymd("2008-01-01") + weeks(control - 1))) %>% 
  join(dv)

# Computing the overall effect of control over 2008-2010
overall <- moving %>% 
  subset(year < 2011) %>% 
  ddply(.(rep, week), summarise,
        tcases = sum(tcases),
        tcases0 = sum(tcases0), 
        ratio = tcases / tcases0,
        diff = tcases0 - tcases,
        dv0 = mean(dv0),
        dvmu = mean(dvmu),
        delta = mean(delta),
        phi = mean(phi))

overall <- dplyr::bind_rows(overall, moving) %>% 
  mutate(year = factor(year, exclude = NULL, labels = c("2008","2009", "2010", "2011", "2012", "overall")))

# Computing the median overall best week for control
overall %>% 
  ddply(.(rep, year), summarise,
        argmin = week[which.min(ratio)]) %>% 
  ddply(.(year), summarise,
        argmin = median(argmin))

# Plotting

postscript("Manuscript/figures/fig3.eps",
           width = 4, height = 5, paper = "special", horizontal = FALSE,
           family = "ArialMT")

overall %>% 
  subset(year %in% c("2008", "2009", "2010", "overall") & week < 53) %>%
  ggplot(aes(week, ratio)) +
  stat_summary(fun.ymin = function(x) quantile(x, 0.1),
               fun.ymax = function(x) quantile(x, 0.9),
               geom = "ribbon", fill = "gray70") +
  stat_summary(fun.y = "median", geom = "line") +
  facet_grid(year ~.) +
  geom_hline(yintercept = 1.0, color = "gray50", linetype = 2) +
  theme_classic() + 
  theme(strip.background = element_blank()) + 
  scale_x_continuous(expand = c(0, 0)) + 
  xlab("week of control") +
  ylab("case ratio over following year")

dev.off()

# Figure 4: effect of phi on effectiveness of control ==========================

overall%>% 
  subset(week == 34 & year == "overall") %>% 
  ggplot(aes(dvmu, ratio)) + 
  geom_point(size = 0.5) + 
  theme_classic() + 
  theme(strip.background = element_blank()) + 
  labs(x = "mosquito mortality rate", y = "case ratio") -> fig4a

overall%>% 
  subset(week == 34 & year == "overall") %>% 
  ggplot(aes(phi, ratio)) + 
  geom_point(size = 0.5) + 
  theme_classic() + 
  theme(strip.background = element_blank()) + 
  labs(x = "case reporting probability", y = "case ratio") -> fig4b

postscript("Manuscript/figures/fig4.eps",
           width = 6, height = 3, paper = "special", horizontal = FALSE,
           family = "ArialMT")

fig4a + fig4b

dev.off()

# State var and control correlations ===========================================

# Cases
cases <- system[, , 11] %>% 
  aaply(1, diff) %>% 
  multiply_by(pop) %>% 
  reshape2::melt(varnames = c("rep", "control"), value.name = "cases") 

# Mosquitoes
mosq <- system[, 2:244, 9] %>% 
  reshape2::melt(varnames = c("rep", "control"), value.name = "mosq") 

# Probability of surviving EIP
infmosq <- system[, , 17] %>% 
  aaply(1, diff) %>% 
  melt(varnames = c("rep", "control"), value.name = "inf")

deadmosq <- system[, , 16] %>% 
  aaply(1, diff) %>% 
  reshape2::melt(varnames = c("rep", "control"), value.name = "dead")

surv <- join(infmosq, deadmosq) %>%
  mutate(total_exits = dead + inf,
         psurv = inf / total_exits) 

# Joining everything
moving <- join(moving, cases) %>% 
  join(mosq) %>% 
  join(surv)

postcorr <- moving %>% 
  subset(control < 191) %>% 
  ddply(.(rep), summarise,
        dv = cor(ratio, dvnat),
        mosq = cor(ratio, mosq), 
        cases = cor(ratio, cases), 
        psurv = cor(ratio, psurv),
        temp = cor(ratio, temp))

colwise(median)(postcorr)

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
  ylab(expression(epsilon[nu])) +
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
rstan::extract(sim, "sigmadv", permute = T)[[1]] %>% hist(main = "A", xlab = expression(sigma[d]), freq = F)
rstan::extract(sim, "sigmarv", permute = T)[[1]] %>% hist(main = "B", xlab = expression(sigma[r]), freq = F)

dev.off()

# Figure S3: epidemiological parameters ===================================================

epiparams <- rstan::extract(sim, c("ro_c", "gamma_c", "delta_c", "dvmu_c"), permute = T)

postscript("Manuscript/figures/figS3.eps",
           width = 6, height = 6,
           family = "ArialMT")

par(mfrow = c(2, 2))

hist(epiparams[[1]], main = "A", xlab = "scaled latent period", freq = F, breaks = 30, ylim = c(0, 1.2))
curve(dgamma(x, 8.3, 8.3), add = T, lwd = 2)

hist(epiparams[[2]], main = "B", xlab = "scaled rate of infectious decay", freq = F, breaks = 30)
curve(dgamma(x, 100, 100), add = T, lwd = 2)

hist(epiparams[[3]], main = "C", xlab = "scaled period of cross-immunity", freq = F, breaks = 30, ylim = c(0, 1.3))
curve(dgamma(x, 10, 10), add = T, lwd = 2)

hist(epiparams[[4]], main = "D", xlab = "scaled mosquito mortality rate", freq = F, breaks = 30)
curve(dgamma(x, 100, 10), add = T, lwd = 2)

dev.off()

# Figure S4: initial conditions ===========================================================

ics <- rstan::extract(sim, c("S0", "E0","I0", "logNv", "dv0", "rv0"), permute = T)

postscript("Manuscript/figures/figS4.eps",
           width = 6, height = 5,
           family = "ArialMT")

par(mfrow = c(2, 3))

hist(ics$S0, freq = F, main = "A", xlab = expression(S[0]))
curve(dbeta(x, 4, 6), add = T, lwd = 2)

hist(ics$E0, freq = F, main = "B", xlab = expression(E[0]), ylim = c(0, 0.015))
curve(dgamma(x, 10, 0.1), add = T, lwd = 2)

hist(ics$I0, freq = F, main = "C", xlab = expression(I[0]), ylim = c(0, 0.02))
curve(dgamma(x, 6, 0.1), add = T, lwd = 2)

hist(ics$logNv, freq = F, main = "D", xlab = expression(log(V[N0])))
curve(dnorm(x, 0.7, 0.3), add = T, lwd = 2)

hist(ics$dv0, freq = F, main = "E", xlab = expression(nu[0]))
curve(dnorm(x, 0, 0.5), add = T, lwd = 2)

hist(ics$rv0, freq = F, main = "F", xlab = expression(r[0]))
curve(dnorm(x, 0, 0.5), add = T, lwd = 2)

dev.off()

# Figure S5: Measurement parameters ============================================

library(truncnorm)

meas <- rstan::extract(sim, c("log_phi_q", "phi_y", "eta_inv_q", "eta_inv_y"), permute = T)

postscript("Manuscript/figures/figS5.eps",
           width = 6, height = 6,
           family = "ArialMT")

par(mfrow = c(2, 2))

hist(meas$log_phi_q, freq = F, main = "A", xlab = expression(log(phi[q])))
curve(dnorm(x, -13, 0.5), add = T, lwd = 2)

hist(meas$phi_y, freq = F, main = "B", xlab = expression(phi[y]))
curve(dbeta(x, 7, 77), add = T, lwd = 2)

hist(meas$eta_inv_q, freq = F, main = "C", xlab = expression(eta[q]))
curve(dtruncnorm(x, a = 0, b = Inf, mean = 0, sd = 5), add = T)

hist(meas$eta_inv_y, freq = F, main = "D", xlab = expression(eta[y]))
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

hist(rowSums(yrep), main = "A", xlab = "total cases", freq = F, breaks = 20)
abline(v = sum(post$yobs), lwd = 2)

hist(rowSums(qrep), main = "B", xlab = "total captured mosquitoes", freq = F, breaks = 20)
abline(v = sum(post$qobs), lwd = 2)

dev.off()

# Figure S8: Min and max ===========================================================

postscript("Manuscript/figures/figS8.eps",
           width = 6, height = 6,
           family = "ArialMT")

par(mfrow = c(2, 2))

hist(apply(yrep, 1, max), main = "A", xlab = "maximum weekly cases", freq = F)
abline(v = max(post$yobs), lwd = 2)

hist(apply(qrep, 1, max), main = "B", xlab = "maximum weekly trap count", freq = F)
abline(v = max(post$qobs), lwd = 2)

hist(apply(yrep, 1, min), main = "C", xlab = "minimum weekly cases", freq = F)
abline(v = min(post$yobs), lwd = 2)

hist(apply(qrep, 1, min), main = "D", xlab = "minimum weekly trap count", freq = F)
abline(v = min(post$qobs), lwd = 2)

dev.off()

# Figure S9: Probability of surviving EIP ======================================

postscript("Manuscript/figures/figS9.eps",
           width = 6, height = 5,
           family = "ArialMT")

moving %>% 
  ggplot(aes(control, psurv)) + 
  stat_summary(fun.ymin = function(x) quantile(x, 0.1),
               fun.ymax = function(x) quantile(x, 0.9),
               geom = "ribbon", fill = "gray70") +
  stat_summary(fun.y = "median", geom = "line") + 
  theme_classic() + 
  scale_x_continuous(expand = c(0, 1), breaks = c(1, 54, 106, 158, 210), labels = c(2008, 2009, 2010, 2011, 2012)) +
  xlab("year") + 
  ylab("probability of surviving EIP") -> figS9a

fig1d + ggtitle("") + figS9a + plot_layout(ncol = 1)
  
dev.off()

