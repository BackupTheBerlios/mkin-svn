library(mkin)
source("R/mkinmod.R")
SFO <- mkinmod(parent = list(type = "SFO"))
FOMC <- mkinmod(parent = list(type = "FOMC"))
SFORB <- mkinmod(parent = list(type = "SFORB"))

SFO_SFO <- mkinmod(parent = list(type = "SFO", to = "m1", sink = TRUE),
        m1 = list(type = "SFO", to = NULL, sink = TRUE))
FOMC_SFO <- mkinmod(parent = list(type="FOMC", to = "m1"), m1 = list(type="SFO"))
SFORB_SFO <- mkinmod(parent = list(type="SFORB", to = "m1"), m1 = list(type="SFO"))
ws <- mkinmod(water = list(type = "SFO", to = "sediment", sink = TRUE),
  sediment = list(type = "SFO"))
ws_back <- mkinmod(water = list(type = "SFO", to = "sediment", sink = TRUE),
  sediment = list(type = "SFO", to = "water"))

source("R/mkinfit.R")
summary(fit <- mkinfit(SFO, FOCUS_2006_A))
fit <- mkinfit(SFO, FOCUS_2006_A)
fit <- mkinfit(FOMC, FOCUS_2006_C)

summary(fit <- mkinfit(SFORB, FOCUS_2006_C))

summary(fit)


source("R/mkinfit.R")
summary(fit <- mkinfit(SFO, FOCUS_2006_C))
summary(fit <- mkinfit(FOMC, FOCUS_2006_C))
summary(fit <- mkinfit(SFORB, FOCUS_2006_C))
summary(fit <- mkinfit(SFO_SFO, FOCUS_2006_D))
summary(fit <- mkinfit(SFO_SFO, FOCUS_2006_D, fixed_parms = "k_parent_sink"))
summary(fit <- mkinfit(SFO_SFO, FOCUS_2006_E))
summary(fit <- mkinfit(FOMC_SFO, FOCUS_2006_D))
summary(fit <- mkinfit(SFORB_SFO, FOCUS_2006_D))
summary(fit <- mkinfit(ws, FOCUS_2006_F, plot = TRUE))
summary(fit <- mkinfit(ws_back, FOCUS_2006_F, lower = 0, 
  parms.ini = c(0.05, 0.1, 0.1, 0.01), plot = TRUE))
summary(fit <- mkinfit(ws_back, FOCUS_2006_F, lower = 0, upper = c(100, 1, 1, 1, 1),
  method = "Pseudo", plot = TRUE))

plot(residual ~ time, data = fit$data, subset = variable == "water")
plot(residual ~ time, data = fit$data, subset = variable == "sediment")


mkinmod <- SFORB_SFO
observed <- FOCUS_2006_D
parms.ini <- rep(0.1, 5)
state.ini <- c(100, 0, 0)
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
