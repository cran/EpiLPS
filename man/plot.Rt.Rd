% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/S3_Rt_plot.R
\name{plot.Rt}
\alias{plot.Rt}
\title{Plot the estimated reproduction number}
\usage{
\method{plot}{Rt}(x, datelab = "7d", cilevel = 0.95, col = "black", cicol = "gray",
 xtickangle = 0, legendpos = "right", title = "Estimated R",
 addfit = c("none","Cori","WT"), theme = "gray", timecut = 0,...)
}
\arguments{
\item{x}{An object of class \code{Rt}.}

\item{datelab}{Spacing for the ticks on the x-axis.}

\item{cilevel}{Level of the credible interval.}

\item{col}{Color of the fitted \eqn{R_t} curve for LPS.}

\item{cicol}{Color for shading the credible envelope.}

\item{xtickangle}{Angle of the x-ticks. Default is 0 (horizontal).}

\item{legendpos}{Position of the legend.}

\item{title}{Title of the plot.}

\item{addfit}{Should an additional \eqn{R_t} fit be added?}

\item{theme}{Theme, either "gray", "classic", "light", "dark"}

\item{timecut}{Cut time points on plot.}

\item{...}{Further arguments to be passed to plot.}
}
\value{
A plot of the fitted time-varying reproduction number.
}
\description{
This routine can be used to plot the estimated reproduction
number based on an object of class \code{Rt}.
}
\examples{
si <- c(0.05, 0.05, 0.1, 0.1, 0.1, 0.1, 0.1, 0.05, 0.05, 0.1, 0.1, 0.1)
epidemic <- episim(si = si, Rpattern = 2, endepi = 30)
epifit <- estimR(incidence = epidemic$y, K = 30, si = si)
plot(epifit)

}
\author{
Oswaldo Gressani \email{oswaldo_gressani@hotmail.fr}
}
