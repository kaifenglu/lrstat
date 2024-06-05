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

#' Acute myelogenous leukemia survival data from the survival package
#'
#' Survival in patients with acute myelogenous leukemia.
#' \describe{
#'   \item{\code{time}}{Survival or censoring time}
#'   \item{\code{status}}{censoring status}
#'   \item{\code{x}}{maintenance chemotherapy given or not}
#' }
"aml"

#' Stanford heart transplant data from the survival package
#'
#' Survival of patients on the waiting list for the Stanford heart
#' transplant program.
#' \describe{
#'   \item{\code{start, stop, event}}{entry and exit time and status for
#'   the time interval}
#'   \item{\code{age}}{age-48 years}
#'   \item{\code{year}}{year of acceptance (in years after Nov 1, 1967)}
#'   \item{\code{surgery}}{prior bypass surgery 1=yes, 0=no}
#'   \item{\code{transplant}}{received transplant 1=yes, 0=no}
#'   \item{\code{id}}{patient id}
#' }
"heart"
