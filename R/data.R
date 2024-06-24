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

#' Tobin's tobit data from the survival package
#'
#' Data from Tobin's original paper.
#' \describe{
#'   \item{\code{durable}}{Durable goods purchase}
#'   \item{\code{age}}{Age in years}
#'   \item{\code{quant}}{Liquidity ratio (x 1000)}
#' }
"tobin"

#' Simulated Concorde trial data from the rpsftm package
#'
#' Patients were randomly assigned to receive treatment immediately
#' or deferred, and those in the deferred arm could cross over and
#' receive treatment. The primary endpoint was time to disease progression.
#' \describe{
#'   \item{\code{id}}{Patient identification number}
#'   \item{\code{def}}{Indicator that the participant was assigned to
#'   the deferred treatment arm}
#'   \item{\code{imm}}{Indicator that the participant was assigned to
#'   the immediate treatment arm}
#'   \item{\code{censyrs}}{The censoring time, in years, corresponding to
#'   the close of study minus the time of entry for each patient}
#'   \item{\code{xo}}{Indicator that crossover occurred}
#'   \item{\code{xoyrs}}{The time, in years, from entry to switching, or
#'   0 for patients in the immediate arm}
#'   \item{\code{prog}}{Indicator of disease progression (1), or
#'   censoring (0)}
#'   \item{\code{progyrs}}{Time, in years, from entry to disease
#'   progression or censoring}
#'   \item{\code{entry}}{The time of entry into the study, measured in years
#'   from the date of randomisation}
#' }
"immdef"
