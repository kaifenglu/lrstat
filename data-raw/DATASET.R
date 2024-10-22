library(dplyr)
library(lrstat)
library(ipcwswitch)

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

# save to data/ folder
usethis::use_data(rawdata, overwrite = TRUE)

# We start with data handling in the ipcwswitch package
# To obtain the times parameter, we can apply the timesTokeep function
# on the same dataframe in the wide format

# names of the repeated measurements and corresponding dates
vect.ps   <- c("myps.v2", "ps.v3", c(paste("ps1.v", seq(1,21), sep="")))
vect.ttc  <- c("myttc.v2", "ttc.v3", c(paste("ttc1.v", seq(1,21), sep="")))
vect.tran <- c("mytran.v1", paste("tran.v", seq(2,21), sep=""))
dates <- c("dexac.v2", "dexac.v3", c(paste("dexac1.v", seq(21), sep="")))
dates2 <- dates[!(dates %in% c("dexac.v2","dexac.v3"))]

# administrative cut off date
SHIdat <- ipcwswitch::SHIdat %>% mutate(cutoffdt = max(ddn))

# times to keep
kept.t <- timesTokeep(SHIdat, id = "id",
                      tstart = "dexac.v2", tstop = "ddn",
                      mes.cov = list(vect.ps, vect.ttc, vect.tran),
                      time.cov = list(dates, dates, dates2))

# Now, we can build the long format
shilong <- wideToLongTDC(
  SHIdat, id = "id",
  tstart = "dexac.v2", tstop = "ddn", event = "status",
  bas.cov = c("agerand", "sex.f","tt_Lnum", "rmh_alea.c", "pathway.f",
              "bras.f","debttCO","ddt.v1", "ddn", "progDate", "cutoffdt"),
  mes.cov = list(f1=vect.ps, f2=vect.ttc, f3=vect.tran),
  time.cov = list(dates, dates, dates2),
  times = kept.t[[1]])

# Put dates in numeric format with tstart at 0
tabi <- split(shilong, shilong$id)
L.tabi   <- length(tabi)
tablist <- lapply(1:L.tabi, function(i){
  refstart <- tabi[[i]]$tstart[1]
  tabi[[i]]$tstart  <- as.numeric(tabi[[i]]$tstart - refstart)
  tabi[[i]]$tstop <- as.numeric(tabi[[i]]$tstop - refstart)
  tabi[[i]]$dpd <- as.numeric(tabi[[i]]$progDate - refstart)
  tabi[[i]]$dco <- as.numeric(tabi[[i]]$debttCO - refstart)
  tabi[[i]]$ady <- as.numeric(tabi[[i]]$ddn - refstart)
  tabi[[i]]$dcut <- as.numeric(tabi[[i]]$cutoffdt - refstart)
  return(tabi[[i]])
})
shilong <- do.call( rbind, tablist )

# Eliminating patient not having initiated the treatment arm
shilong <- shilong[!is.na(shilong$ddt.v1),]

colnames(shilong)[16:18] <- c("ps", "ttc", "tran")
shilong <- shilong %>%
  select(-c(ddt.v1, debttCO, ddn, progDate, cutoffdt)) %>%
  mutate(event = as.integer(event), pd = !is.na(dpd), co = !is.na(dco))

usethis::use_data(shilong, overwrite = TRUE)


ingots <- readxl::read_xlsx("notes/ingots.xlsx")
usethis::use_data(ingots, overwrite = TRUE)


six <- readxl::read_excel("notes/sixcities.xlsx")
usethis::use_data(six, overwrite = TRUE)


aml = survival::aml
heart = survival::heart
tobin = survival::tobin
immdef = rpsftm::immdef

usethis::use_data(aml, overwrite = TRUE)
usethis::use_data(heart, overwrite = TRUE)
usethis::use_data(tobin, overwrite = TRUE)
usethis::use_data(immdef, overwrite = TRUE)

sexagg = logistf::sexagg
usethis::use_data(sexagg, overwrite = TRUE)
