#' Mortality data for Belgium in 2021
#'
#' @docType data
#'
#' @description A data frame with six columns including the calendar date of
#' death, the date on which cases have been reported and the number of cases/deaths.
#'
#' @usage data(cov19mort2021)
#'
#' @format
#' \describe{
#'  \item{\code{t}}{A numeric variable associated to a calendar date.}
#'  \item{\code{d}}{A numeric variable indicating the delay of reporting.}
#'  \item{\code{Date}}{The calendar date of the death event.}
#'  \item{\code{Rep.date}}{The calendar date for the reporting of the death event.}
#'  \item{\code{Cases}}{The number of cases/deaths.}
#'  \item{\code{Reported}}{Indicates whether cases are already reported or not.}
#' }
#'
#'@source \url{https://epistat.sciensano.be/covid/}
#'
"cov19mort2021"
