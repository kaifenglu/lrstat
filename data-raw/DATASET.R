library(lrstat)
library(dplyr)

lr1 <- lrsamplesize(
  beta = 0.2, kMax = 1,
  accrualTime = seq(0, 8),
  accrualIntensity = 26/9*seq(1, 9),
  piecewiseSurvivalTime = c(0, 6),
  stratumFraction = c(0.2, 0.8),
  lambda1 = c(0.0533, 0.0309, 1.5*0.0533, 1.5*0.0309),
  lambda2 = c(0.0533, 0.0533, 1.5*0.0533, 1.5*0.0533),
  gamma1 = -log(1-0.05)/12,
  gamma2 = -log(1-0.05)/12,
  accrualDuration = NA,
  followupTime = 18, fixedFollowup = FALSE)

sim1 <- lrsim(
  kMax = 1,
  criticalValues = 1.96,
  accrualTime = seq(0, 8),
  accrualIntensity = 26/9*seq(1, 9),
  piecewiseSurvivalTime = c(0, 6),
  stratumFraction = c(0.2, 0.8),
  lambda1 = c(0.0533, 0.0309, 1.5*0.0533, 1.5*0.0309),
  lambda2 = c(0.0533, 0.0533, 1.5*0.0533, 1.5*0.0533),
  gamma1 = -log(1-0.05)/12,
  gamma2 = -log(1-0.05)/12,
  accrualDuration = 22.9,
  plannedEvents = 376,
  maxNumberOfIterations = 1000,
  maxNumberOfRawDatasetsPerStage = 10,
  seed = 314159)

rawdata <- sim1$rawdata %>%
  select(iterationNumber, arrivalTime, stratum, treatmentGroup,
         timeUnderObservation, event, dropoutEvent)

aml <- survival::aml
heart <- survival::heart

# save to data/ folder
usethis::use_data(rawdata, overwrite = TRUE)
usethis::use_data(aml, overwrite = TRUE)
usethis::use_data(heart, overwrite = TRUE)
