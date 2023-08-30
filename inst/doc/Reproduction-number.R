## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- include=FALSE-----------------------------------------------------------
library("EpiLPS")
library("EpiEstim")

## ---- out.width='100%', fig.align='center', fig.cap='', echo=FALSE------------
knitr::include_graphics('EpiLPS-Chart.png')

## ---- out.width='100%', fig.align='center', fig.cap='', echo=FALSE------------
knitr::include_graphics('EpiLPS-Arch.png')

## ----devepilps, eval=FALSE----------------------------------------------------
#  install.packages("devtools")
#  devtools::install_github("oswaldogressani/EpiLPS")

## ----si_interval, fig.align='center', fig.width=6.3---------------------------
si_spec <- Idist(mean = 2.7, sd = 1.3, dist = "gamma")
si <- si_spec$pvec
si
plot(si_spec, titlesize = 12)

## ----episim, fig.align='center', fig.width=7, fig.height=2.8------------------
set.seed(123)
datasim <- episim(si = si, Rpattern = 5, endepi = 40, dist = "negbin", overdisp = 15, plotsim = TRUE)
incidence <- datasim$y
incidence

## ----LPSfit1------------------------------------------------------------------
LPSfit <- estimR(incidence = incidence, si = si)
class(LPSfit)
knitr::kable(tail(LPSfit$RLPS[,1:8]))

## ----summaryLPSfit------------------------------------------------------------
summary(LPSfit)

## ----plotLPSfit, fig.align='center', fig.width=7, fig.height=4----------------
plot(LPSfit)

## ----estimMCMC, eval=FALSE----------------------------------------------------
#  LPSfitmcmc <- estimRmcmc(incidence = incidence, si = si, niter = 5000, burnin = 2000)

## ----estimMCMC2, echo = FALSE, results='hide'---------------------------------
LPSfitmcmc <- estimRmcmc(incidence = incidence, si = si, niter = 5000, burnin = 2000)

## ----estimMCMC3---------------------------------------------------------------
summary(LPSfitmcmc)

## ----addCori, fig.width=7, fig.height=4---------------------------------------
LPSfit2 <- estimR(incidence = incidence, si = si, CoriR = TRUE)
knitr::kable(tail(LPSfit2$RCori[,1:8]))
plot(LPSfit2, addfit = "Cori")

## ----traditionalplot, fig.width=7, fig.height=3.8-----------------------------
tt <- seq(8, 40, by = 1)
Rtrue <- sapply(tt, datasim$Rtrue)
plot(tt, Rtrue, type = "l", xlab = "Time", ylab = "R", ylim = c(0,4), lwd = 2)
lines(tt, LPSfit2$RLPS$R[-(1:7)], col = "red", lwd = 2)
lines(tt, LPSfit2$RCori$`Mean(R)`, col = "blue", lwd = 2)
legend("topright", col = c("black","red","blue"), 
       c("True R","EpiLPS","EpiEstim"), bty = "n", lty = c(1,1,1))

## ----customizeplot, fig.width=7, fig.height=5.3-------------------------------
gridExtra::grid.arrange(
  plot(LPSfit, col = "firebrick", legendpos = "top", cicol = "orange"),
  plot(LPSfit, col = "forestgreen", legendpos = "none", cicol = "green",
       theme = "light", title = "Reproduction number"),
  plot(LPSfit, col = "darkblue", legendpos = "none", cicol = "orchid", theme = "classic"),
  plot(LPSfit, col = "white", legendpos = "none", cicol = "gray",
       theme = "dark"), nrow = 2, ncol = 2)

## ----zikaload, fig.width=7, fig.height=3.8------------------------------------
# Loading the data
data("zika2015")
zika <- zika2015
plotIncidence <- epicurve(zika$incidence, dates = zika$dates, datelab = "14d", title = "Zika incidence", xtickangle = 70)
plotIncidence

## ----zikafit------------------------------------------------------------------
# Specification of the serial interval
si <- Idist(mean = 7, sd = 1.5)
siplot <- plot(si, titlesize = 11)
epifit <- estimR(zika$incidence, dates = zika$dates, si = si$pvec)
summary(epifit)

## ----zikaplot, fig.width=7, fig.height=5.3------------------------------------
# Plot the smoothed epidemic curve
plotsmooth <- epicurve(zika$incidence, dates = zika$dates, datelab = "14d", 
                       smooth = epifit, smoothcol = "orange", 
                       title = "Zika incidence (smoothed)",
                       xtickangle = 70)

# Plot of the estimated reproduction number
Rplot <- plot(epifit, datelab = "7d", xtickangle = 70, legendpos = "none", col = "forestgreen")

# Show all plots together
gridExtra::grid.arrange(plotIncidence, plotsmooth, siplot, Rplot, nrow = 2, ncol = 2)


## ----measles, fig.width=7.5, fig.height=5.7, warning=FALSE--------------------
data(Measles1861)
measlesDAT <- Measles1861
measles_incid <- measlesDAT$incidence
measles_si <- measlesDAT$si_distr
epifit_measles <- estimR(measles_incid, si = measles_si, CoriR = T)
epicurve_measles<- epicurve(measles_incid, datelab = "1d", title = "Measles, Hagelloch, Germany, 1861",
                            col = "lightblue3", smooth = epifit_measles, smoothcol = "dodgerblue4")
Rplot_measles <- plot(epifit_measles, timecut = 6, legendpos = "none")
Rplot_measles2 <- plot(epifit_measles, addfit = "Cori", timecut = 6, legendpos = "top")
gridExtra::grid.arrange(epicurve_measles, Rplot_measles, Rplot_measles2, nrow = 2, ncol = 2)

## ----flu1918, fig.width=7.5, fig.height=5.7, warning=FALSE--------------------
data(Flu1918)
fluDAT <- Flu1918
flu_incid <- fluDAT$incidence
flu_si <- fluDAT$si_distr[-1]
epifit_flu <- estimR(flu_incid, si = flu_si, CoriR = T)
epicurve_flu <- epicurve(flu_incid, datelab = "7d", title = "Influenza, Baltimore, 1918",
                            col = "orange", smooth = epifit_flu, smoothcol = "firebrick")
Rplot_flu <- plot(epifit_flu, legendpos = "none")
Rplot_flu2 <- plot(epifit_flu, addfit = "Cori", legendpos = "top")
siplot_flu <- plot(Idist(probs = flu_si), barcol = "indianred1")
gridExtra::grid.arrange(epicurve_flu, Rplot_flu, Rplot_flu2, 
                        siplot_flu, nrow = 2, ncol = 2)

## ----sars2003, fig.width=7.5, fig.height=5.7, warning=FALSE-------------------
data("SARS2003")
sarsDAT <- SARS2003
sars_incid <- sarsDAT$incidence
sars_si <- sarsDAT$si_distr[-1]
epifit_sars <- estimR(sars_incid, si = sars_si, CoriR = T)
epicurve_sars <- epicurve(sars_incid, datelab = "7d", title = "SARS, Hong Kong, 2003",
                         col = "ivory4", smooth = epifit_sars, smoothcol = "darkviolet")
Rplot_sars <- plot(epifit_sars, legendpos = "none")
Rplot_sars2 <- plot(epifit_sars, addfit = "Cori", legendpos = "top")
gridExtra::grid.arrange(epicurve_sars, Rplot_sars, Rplot_sars2, nrow = 2, ncol = 2)

## ----perf1, eval=FALSE--------------------------------------------------------
#  perfcheck <- perfRestim(si = "mers", scenario = 5, days = 60,  seed = 1905, overdisp = 50)
#  perfcheck$LPS

## ----perf1-bis, echo = FALSE, results='hide', eval = TRUE---------------------
perfcheck <- perfRestim(si = "mers", scenario = 5, days = 60,  seed = 1905, overdisp = 50)

## ----perf2, fig.width=7.8, fig.height=6.2, eval = TRUE------------------------
perfcheck$LPS
gridExtra::grid.arrange(perfcheck$inciplot, perfcheck$Rlpsplot,perfcheck$Repiestimplot, 
                        nrow = 2, ncol = 2)

