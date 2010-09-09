# $Id$
# Copyright (C) 2010 Johannes Ranke
# Test script containing various snippets of R code for testing the mkin package
# Part of these are already integrated into the package vignette and/or the 
# unit tests under inst/unitTests

library(mkin)
source("trunk/R/DFOP.solution.R")
source("trunk/R/HS.solution.R")
source("trunk/R/mkinmod.R")
source("trunk/R/mkinfit.R")
SFO <- mkinmod(parent = list(type = "SFO"))
FOMC <- mkinmod(parent = list(type = "FOMC"))
SFORB <- mkinmod(parent = list(type = "SFORB"))
DFOP <- mkinmod(parent = list(type = "DFOP"))
HS <- mkinmod(parent = list(type = "HS"))

SFO_SFO <- mkinmod(parent = list(type = "SFO", to = "m1", sink = TRUE),
        m1 = list(type = "SFO", to = NULL, sink = TRUE))
FOMC_SFO <- mkinmod(parent = list(type="FOMC", to = "m1"), m1 = list(type="SFO"))
DFOP_SFO <- mkinmod(parent = list(type="DFOP", to = "m1"), m1 = list(type="SFO"))
HS_SFO <- mkinmod(parent = list(type="HS", to = "m1"), m1 = list(type="SFO"))
SFORB_SFO <- mkinmod(parent = list(type="SFORB", to = "m1"), m1 = list(type="SFO"))
ws <- mkinmod(water = list(type = "SFO", to = "sediment", sink = TRUE),
  sediment = list(type = "SFO"))
ws_back <- mkinmod(water = list(type = "SFO", to = "sediment", sink = TRUE),
  sediment = list(type = "SFO", to = "water"))

source("R/mkinfit.R")
summary(fit <- mkinfit(SFO, FOCUS_2006_A))
summary(fit <- mkinfit(FOMC, FOCUS_2006_C))

summary(fit <- mkinfit(SFORB, FOCUS_2006_C, eigen=FALSE))

summary(fit)


source("trunk/R/mkinfit.R")
source("trunk/R/mkinmod.R")
source("trunk/R/mkinplot.R")
summary(fit <- mkinfit(SFO, FOCUS_2006_C))
summary(fit <- mkinfit(FOMC, FOCUS_2006_C))
summary(fit <- mkinfit(DFOP, FOCUS_2006_A, parms.ini = c(0.1, 0.01, 0.5)))
summary(fit <- mkinfit(DFOP, FOCUS_2006_B, parms.ini = c(0.1, 0.01, 0.5)))
summary(fit <- mkinfit(DFOP, FOCUS_2006_C, parms.ini = c(0.1, 0.01, 0.5)))
summary(fit <- mkinfit(DFOP, FOCUS_2006_D, parms.ini = c(0.1, 0.01, 0.5)))
summary(fit <- mkinfit(HS, FOCUS_2006_C, parms.ini = c(0.1, 0.01, 10), plot=TRUE))
summary(fit <- mkinfit(HS, FOCUS_2006_D, parms.ini = c(0.1, 0.01, 10), plot=TRUE))
summary(fit <- mkinfit(SFORB, FOCUS_2006_C))
summary(fit <- mkinfit(SFO_SFO, FOCUS_2006_D))
summary(fit <- mkinfit(SFO_SFO, FOCUS_2006_D, fixed_parms = "k_parent_sink"))
summary(fit <- mkinfit(SFO_SFO, FOCUS_2006_E))
summary(fit <- mkinfit(SFO_SFO, FOCUS_2006_E, eigen=FALSE))
summary(fit <- mkinfit(FOMC_SFO, FOCUS_2006_D, plot=TRUE))
summary(fit <- mkinfit(DFOP, FOCUS_2006_D, parms.ini = c(0.7, 0.09, 0.1)))
summary(fit <- mkinfit(DFOP_SFO, FOCUS_2006_D, parms.ini = c(0.7, 0.09, 0.1, 0.1, 0.5), plot=TRUE))
summary(fit <- mkinfit(HS, FOCUS_2006_A, parms.ini = c(0.1, 0.01, 10), plot=TRUE))
summary(fit <- mkinfit(HS, FOCUS_2006_B, parms.ini = c(0.1, 0.01, 10), plot=TRUE))
summary(fit <- mkinfit(HS, FOCUS_2006_C, parms.ini = c(0.1, 0.01, 10), plot=TRUE))
summary(fit <- mkinfit(SFO, FOCUS_2006_D, plot=TRUE))
summary(fit <- mkinfit(HS, FOCUS_2006_D, parms.ini = c(0.1, 0.01, 10), plot=TRUE))
summary(fit <- mkinfit(HS_SFO, FOCUS_2006_D, parms.ini = c(0.09, 0.11, 14, 0.1, 0.5), 
  plot=TRUE))
summary(fit <- mkinfit(SFORB_SFO, FOCUS_2006_D))
summary(fit <- mkinfit(ws, FOCUS_2006_F, plot = TRUE))
summary(fit <- mkinfit(ws_back, FOCUS_2006_F, lower = 0, 
  parms.ini = c(0.05, 0.1, 0.1, 0.01), plot = TRUE))
summary(fit <- mkinfit(ws_back, FOCUS_2006_F, lower = 0, upper = c(100, 1, 1, 1, 1),
  method = "Pseudo", plot = TRUE))

plot(residual ~ time, data = fit$data, subset = variable == "water")
plot(residual ~ time, data = fit$data, subset = variable == "sediment")


# Debug the SFO model
mkinmod <- SFO
observed <- FOCUS_2006_A
parms.ini <- rep(0.1, 1)
state.ini <- c(100)
fixed_parms <- NULL
fixed_initials <- character(0)
lower = NULL
upper = Inf
plot = FALSE
quiet = FALSE
err = NULL
weight = "none"
scaleVar = FALSE

summary(fit <- mkinfit(SFORB_SFO, FOCUS_2006_D))
P <- c(state.ini.optim, parms.optim)
modCost(cost, P)
fit <- modFit(cost, c(state.ini.optim, parms.optim), lower = 0)
MCMC <- modMCMC(cost, c(state.ini.optim, parms.optim))

# Debug the DFOP model
source("trunk/R/mkinmod.R")
DFOP <- mkinmod(parent = list(type = "DFOP"))
mkinmod <- DFOP
observed <- FOCUS_2006_C
parms.ini <- c(0.45, 0.02, 0.85)
state.ini <- c(100)
lower = 0
upper = Inf
fixed_parms <- NULL
fixed_initials <- character(0)
plot = FALSE
quiet = FALSE
err = NULL
weight = "none"
scaleVar = FALSE
atol = 1e-6

P <- c(state.ini.optim, parms.optim)
summary(fit <- mkinfit(DFOP, FOCUS_2006_A, plot=TRUE, parms.ini = c(0.1, 0.01, 0.5)))
summary(fit <- mkinfit(DFOP, FOCUS_2006_B, plot=TRUE, parms.ini = c(0.1, 0.01, 0.5)))
summary(fit <- mkinfit(DFOP, FOCUS_2006_C, plot=TRUE, parms.ini = c(0.1, 0.01, 0.5)))

summary(fit <- mkinfit(SFORB_SFO, FOCUS_2006_D))

d <- FOCUS_2006_C
d2 <- FOCUS_2006_C
d2$name <- "m1"
d2$value <- cumsum(c(0, 9.5, 12, 6.8, 2, 1.2, 0.2, -0.1, -0.5))
d3 <- FOCUS_2006_C
d3$name <- "m2"
d3$value <- cumsum(c(0, 6.3, 7.2, 2.0, 1.0, -2, -3.5, -4, -4.5))
(observed <- rbind(d, d2, d3))

FOMC <- mkinmod(parent = list(type="FOMC"))
FOMC_SFO <- mkinmod(parent = list(type="FOMC", to = "m1"), m1 = list(type="SFO"))
FOMC_SFO2 <- mkinmod(parent = list(type="FOMC", to = c("m1", "m2")),
  m1 = list(type = "SFO"), m2 = list(type="SFO"))

summary(mkinfit(SFORB_SFO, observed, plot=TRUE)
summary(mkinfit(SFO_SFO, observed, plot=TRUE))
summary(mkinfit(FOMC, observed, plot=TRUE))
fit <- mkinfit(FOMC_SFO, observed, plot=TRUE)
fit <- mkinfit(FOMC_SFO, observed, parms.ini = c(1, 10, 0.1, 0.5), plot=TRUE)
fit <- mkinfit(FOMC_SFO2, observed, parms.ini = c(1, 10, 0.1, 0.1, 0.5, 0.5), plot = TRUE)
summary(fit)

schaefer07_complex_case <- data.frame(
  time = c(0, 1, 3, 7, 14, 30, 62, 100),
  parent = c(93.2, 89.4, 79.7, 61.1, 48.2, 15.9, 6.5, 6.0),
  A1 = c(NA, NA, 0.55, 6.87, 17.08, 21.68, 15.77, 13.63),
  B1 = c(NA, NA, NA, 0.55, 2.31, 15.76, 6.36, 3.74),
  C1 = c(NA, 0.55, 3.20, 5.46, 12.55, 10.45, 4.74, 4.33),
  A2 = c(NA, 0.55, 1.41, 0.55, 1.29, 1.95, 3.54, 3.86))

parameters = c("degradation rate", "DT50", "formation fraction")
schaefer07_complex_results <- data.frame(
  compound = c(rep("parent", 2), rep(paste("metabolite", c("A1", "B1", "C1", "A2")), each = 3)),
  parameter = parameters[c(1,2, rep(c(3,1,2), 4))],
  KinGUI = c(0.0496, 13.99, 0.3803, 0.0139, 49.96, 0.1866, 0.0175, 39.61,
     0.4331, 0.0638, 10.87, 0.4529, 0.0245, 28.24),
  ModelMaker = c(0.0506, 13.69, 0.3696, 0.0136, 50.89, 0.1818, 0.0172, 40.24,
     0.4486, 0.0700, 9.90, 0.4559, 0.0244, 28.45))
schaefer07_complex_results <- transform(schaefer07_complex_results, 
  deviation = abs(round(100 * ((KinGUI - ModelMaker)/ModelMaker), 1)))
save(schaefer07_complex_results, schaefer07_complex_case, 
  file = "schaefer07_complex_case.RData")

data <- mkin_wide_to_long(schaefer07_complex_case, time = "time")

schaefer07_complex_model <- mkinmod(
  parent = list(type = "SFO", to = c("A1", "B1", "C1"), sink = FALSE),
  A1 = list(type = "SFO", to = "A2"),
  B1 = list(type = "SFO"),
  C1 = list(type = "SFO"),
  A2 = list(type = "SFO"))


system.time(fit <- mkinfit(schaefer07_complex_model, data, plot=TRUE))
summary(fit, data=FALSE)


mkinmod <- SFORB_SFO
observed <- FOCUS_2006_D
parms.ini = rep(0.1, length(mkinmod$parms))
state.ini = c(100, rep(0, length(mkinmod$diffs) - 1))
fixed_parms <- NULL
fixed_initials = names(mkinmod$diffs)[-1]
lower = 0
upper = Inf
plot = FALSE
quiet = FALSE
err = NULL
weight = "none"
scaleVar = FALSE

P <- c(state.ini.optim, parms.optim)
fit <- modFit(cost, c(state.ini.optim, parms.optim), lower = 0)
