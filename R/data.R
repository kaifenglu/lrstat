#' A simulated time-to-event data set with 10 replications
#'
#' A simulated data set with stratification and delayed treatment effect:
#' \describe{
#'   \item{\code{iterationNumber}}{The iteration number}
#'   \item{\code{arrivalTime}}{The enrollment time for the subject}
#'   \item{\code{stratum}}{The stratum for the subject}
#'   \item{\code{treatmentGroup}}{The treatment group for the subject}
#'   \item{\code{timeUnderObservation}}{The time under observation since
#'   randomization}
#'   \item{\code{event}}{Whether the subject experienced the event}
#'   \item{\code{dropoutEvent}}{Whether the subject dropped out}
#' }
"rawdata"
