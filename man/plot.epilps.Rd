% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/S3_epilps_plot.R
\name{plot.epilps}
\alias{plot.epilps}
\title{Plot the EpiLPS fitted epidemic curve and reproduction number}
\usage{
\method{plot}{epilps}(x, plotout = c("rt", "epicurve"), dates = NULL,
     datelab = c("7d", "1m", "3m", "6m"), overlayEpiestim = FALSE, Rtitle = "",
     epititle = "", rtcol = "red", cicol = "gray", transparency = 0.5,
     epicol = "red", epiestimcol = "lightslateblue", incibars = FALSE, barwidth = 0.35,
     themetype = c("gray", "classic", "light", "dark"),  tcut = NULL, titlesize = 15,
     xtitlesize = 13, ytitlesize = 13, ...)
}
\arguments{
\item{x}{An object of class \code{epilps}.}

\item{plotout}{The type of plot, either "rt" for showing the reproduction
number or "epicurve" for showing the epidemic curve.}

\item{dates}{A vector of dates in format "YY-MM-DD".}

\item{datelab}{The spacing for ticks on the x-axis. Either 7 days, 1 month,
3 months or 6 months.}

\item{overlayEpiestim}{Should the EpiEstim fit be overlayed?}

\item{Rtitle}{The title for the plot of R.}

\item{epititle}{The title for the plot of the epidemic curve.}

\item{rtcol}{Color for the reproduction number curve fit.}

\item{cicol}{Color for shading the credible envelope.}

\item{transparency}{Controls the transparency of the credible envelope.}

\item{epicol}{The color for the epidemic curve.}

\item{epiestimcol}{The color for the EpiEstim Rt estimate.}

\item{incibars}{Should the bars of the incidence time series be shown?}

\item{barwidth}{The bar width associated to the incidence time series.}

\item{themetype}{Type of theme for the plot.}

\item{tcut}{Remove early estimates (starting day 8 in plot).}

\item{titlesize}{The size of the plot title. Default is 15.}

\item{xtitlesize}{The size of title and text on x axis. Default is 13.}

\item{ytitlesize}{The size of title and text on y axis. Default is 13.}

\item{...}{Further arguments to be passed to plot.}
}
\value{
A plot of the fitted time-varying reproduction number (default) or
 the epidemic curve.
}
\description{
This routine can be used to plot the estimated epidemic curve
 and reproduction number with EpiLPS.
}
\examples{
si <- c(0.05, 0.05, 0.1, 0.1, 0.1, 0.1, 0.1, 0.05, 0.05, 0.1, 0.1, 0.1)
epidemic <- episim(serial_interval = si, Rpattern = 2, endepi = 30)
epifit <- epilps(incidence = epidemic$y, K = 30, serial_interval = si)
# gridExtra::grid.arrange(plot(epifit, Rtitle = "Estimated R"),
#                         plot(epifit, plotout = "epicurve", epicol = "blue",
#                         epititle = "Epidemic curve"), nrow = 2)

}
\author{
Oswaldo Gressani \email{oswaldo_gressani@hotmail.fr}
}
