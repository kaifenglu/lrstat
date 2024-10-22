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

#' Simulated CONCORDE trial data from the rpsftm package
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

#' The randomized clinical trial SHIVA data in long format from the
#' ipcwswitch package
#'
#' The original SHIdat data set contains an anonymized excerpt of data from
#' the SHIVA01 trial. This was the first randomized clinical trial that
#' aimed at comparing molecularly targeted therapy based on tumor
#' profiling (MTA) versus conventional therapy (CT) for advanced cancer.
#' Patients were randomly assigned to receive the active or control
#' treatment and may switch to the other arm or subsequent anti-cancer
#' therapy upon disease progression.
#' The restructured data is in the long format.
#' \describe{
#'   \item{\code{id}}{The patient's identifier}
#'   \item{\code{tstart}}{The start of the time interval}
#'   \item{\code{tstop}}{The end of the time interval}
#'   \item{\code{event}}{Whether the patient died at the end of the interval}
#'   \item{\code{agerand}}{The patient's age (in years) at randomization}
#'   \item{\code{sex.f}}{The patients' gender, either Male or Female}
#'   \item{\code{tt_Lnum}}{The number of previous lines of treatment}
#'   \item{\code{rmh_alea.c}}{The Royal Marsden Hospital score segregated
#'   into two categories}
#'   \item{\code{pathway.f}}{The molecular pathway altered (the hormone
#'   receptors pathway, the PI3K/ AKT/mTOR pathway, and the RAF/MEK pathway)}
#'   \item{\code{bras.f}}{The patient's randomized arm, either MTA or CT}
#'   \item{\code{ps}}{The ECOG performance status}
#'   \item{\code{ttc}}{The presence of concomitant treatments}
#'   \item{\code{tran}}{The use of platelet transfusions}
#'   \item{\code{dpd}}{The relative day of a potential progression}
#'   \item{\code{dco}}{The relative day of treatment switching}
#'   \item{\code{ady}}{The relative day of the latest news}
#'   \item{\code{dcut}}{The relative day of administrative cutoff}
#'   \item{\code{pd}}{Whether the patient had disease progression}
#'   \item{\code{co}}{Whether the patient switched treatment}
#' }
"shilong"

#' The binary data from Cox and Snell (1989, pp. 10-11).
#'
#' The dataset consits of the number of ingots not ready for rolling
#' and the number of ingots ready for rolling for a number of
#' combinations of heating time and soaking time.
#'
#' \describe{
#'   \item{\code{Heat}}{The heating time}
#'   \item{\code{Soak}}{The soaking time}
#'   \item{\code{NotReady}}{Response indicator, with a value 1 for units
#'   not ready for rolling (event) and a value of 0 for units ready for
#'   rolling (nonevent)}
#'   \item{\code{Freq}}{The frequency of occurrence of each combination of
#'   \code{Heat}, \code{Soak}, and \code{NotReady}}
#' }
"ingots"

#' The repeated measures data from the "Six Cities" study of the health
#' effects of air pollution (Ware et al. 1984).
#'
#' The data analyzed are the 16 selected cases in Lipsitz et al. (1994).
#' The binary response is the wheezing status of 16 children at ages 9,
#' 10, 11, and 12 years. A value of 1 of wheezing status indicates the
#' occurrence of wheezing. The explanatory variables city of residence,
#' age, and maternal smoking status at the particular age.
#'
#' \describe{
#'   \item{\code{case}}{case id}
#'   \item{\code{city}}{city of residence}
#'   \item{\code{age}}{age of the child}
#'   \item{\code{smoke}}{maternal smoking status}
#'   \item{\code{wheeze}}{wheezing status}
#' }
"six"

#' Urinary tract infection data from the logistf package
#'
#' This data set deals with urinary tract infection in sexually active
#' college women, along with covariate information on age an
#' contraceptive use. The variables are all binary and coded
#' in 1 (condition is present) and 0 (condition is absent).
#'
#' \describe{
#'   \item{\code{case}}{urinary tract infection, the study outcome variable}
#'   \item{\code{age}}{>= 24 years}
#'   \item{\code{dia}}{use of diaphragm}
#'   \item{\code{oc}}{use of oral contraceptive}
#'   \item{\code{vic}}{use of condom}
#'   \item{\code{vicl}}{use of lubricated condom}
#'   \item{\code{vis}}{use of spermicide}
#' }
"sexagg"
