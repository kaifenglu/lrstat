library(lrstat)
Rcpp::sourceCpp("src/lrsim2e.cpp")

# PFS: Arm A vs. C
lrpwr1 <- lrsamplesize(
  beta = 0.1,
  kMax = 2,
  informationRates = c(0.72, 1),
  alpha = 0.005,
  typeAlphaSpending = "sfOF",
  allocationRatio = 2,
  accrualTime = c(0, 8),
  accrualIntensity = c(10, 28)*3/5,  # Arm A + Arm C (2+1) out of (2+2+1)
  lambda1 = 0.60*log(2)/12,  # Arm A
  lambda2 = log(2)/12,  # Arm C
  accrualDuration = 30.143,
  followupTime = NA,
  typeOfComputation = "Schoenfeld")


# OS: Arm A vs. C
lrpwr2 <- lrsamplesize(
  beta = 0.25,
  kMax = 3,
  informationRates = c(0.55, 0.70, 1),
  alpha = 0.020,
  typeAlphaSpending = "sfOF",
  allocationRatio = 2,
  accrualTime = c(0, 8),
  accrualIntensity = c(10, 28)*3/5,  # Arm A + Arm C (2+1) out of (2+2+1)
  lambda1 = 0.65*log(2)/30,  # Arm A
  lambda2 = log(2)/30,  # Arm C
  accrualDuration = 30.143,
  followupTime = NA,
  typeOfComputation = "Schoenfeld")


# time to expect 186 and 259 PFS events
t1 = caltime(
  nevents = c(186, 259), 
  allocationRatioPlanned = 2,
  accrualTime = c(0, 8),
  accrualIntensity = c(10, 28)*3/5,  # Arm A + Arm C (2+1) out of (2+2+1)
  lambda1 = 0.60*log(2)/12,  # Arm A
  lambda2 = log(2)/12,  # Arm C
  accrualDuration = 30.143,
  followupTime = 100)


# time to expect 183 OS events
t2 = caltime(
  nevents = 183,
  allocationRatioPlanned = 2,
  accrualTime = c(0, 8),
  accrualIntensity = c(10, 28)*3/5,  # Arm A + Arm C (2+1) out of (2+2+1)
  lambda1 = 0.65*log(2)/30,  # Arm A
  lambda2 = log(2)/30,  # Arm C
  accrualDuration = 30.143,
  followupTime = 100)


t = c(t1, t2)

# expected number of events at each look by endpoint and comparision
lr13e1 <- lrstat(
  time = t,
  allocationRatioPlanned = 2,
  accrualTime = c(0, 8),
  accrualIntensity = c(10, 28)*3/5,  # Arm A + Arm C (2+1) out of (2+2+1)
  lambda1 = 0.60*log(2)/12,  # Arm A
  lambda2 = log(2)/12,  # Arm C
  accrualDuration = 30.143,
  followupTime = 100, predictEventOnly = 1)

lr23e1 <- lrstat(
  time = t,
  allocationRatioPlanned = 2,
  accrualTime = c(0, 8),
  accrualIntensity = c(10, 28)*3/5,  # Arm B + Arm C (2+1) out of (2+2+1)
  lambda1 = 0.70*log(2)/12,  # Arm B
  lambda2 = log(2)/12,  # Arm C
  accrualDuration = 30.143,
  followupTime = 100, predictEventOnly = 1)

lr12e1 <- lrstat(
  time = t,
  allocationRatioPlanned = 1,
  accrualTime = c(0, 8),
  accrualIntensity = c(10, 28)*4/5,  # Arm A + Arm B (2+2) out of (2+2+1)
  lambda1 = 0.60*log(2)/12,  # Arm A
  lambda2 = 0.70*log(2)/12,  # Arm B
  accrualDuration = 30.143,
  followupTime = 100, predictEventOnly = 1)


lr13e2 <- lrstat(
  time = t,
  allocationRatioPlanned = 2,
  accrualTime = c(0, 8),
  accrualIntensity = c(10, 28)*3/5,  # Arm A + Arm C (2+1) out of (2+2+1)
  lambda1 = 0.65*log(2)/30,  # Arm A
  lambda2 = log(2)/30,  # Arm C
  accrualDuration = 30.143,
  followupTime = 100, predictEventOnly = 1)

lr23e2 <- lrstat(
  time = t,
  allocationRatioPlanned = 2,
  accrualTime = c(0, 8),
  accrualIntensity = c(10, 28)*3/5,  # Arm B + Arm C (2+1) out of (2+2+1)
  lambda1 = 0.75*log(2)/30,  # Arm B
  lambda2 = log(2)/30,  # Arm C
  accrualDuration = 30.143,
  followupTime = 100, predictEventOnly = 1)

lr12e2 <- lrstat(
  time = t,
  allocationRatioPlanned = 1,
  accrualTime = c(0, 8),
  accrualIntensity = c(10, 28)*4/5,  # Arm A + Arm B (2+2) out of (2+2+1)
  lambda1 = 0.65*log(2)/30,  # Arm A
  lambda2 = 0.75*log(2)/30,  # Arm B
  accrualDuration = 30.143,
  followupTime = 100, predictEventOnly = 1)



# simulated log-rank test statistics
lr <- lrsim2e(
  kMax = 3,
  kMaxe1 = 2,
  allocation1 = 2,
  allocation2 = 2,
  allocation3 = 1,
  accrualTime = c(0, 8),
  accrualIntensity = c(10, 28),
  piecewiseSurvivalTime = 0,
  stratumFraction = 1,
  rho = 0,
  lambda1e1 = log(2)/12*0.60,
  lambda2e1 = log(2)/12*0.70,
  lambda3e1 = log(2)/12,
  lambda1e2 = log(2)/30*0.65,
  lambda2e2 = log(2)/30*0.75,
  lambda3e2 = log(2)/30,
  gamma1e1 = 0,
  gamma2e1 = 0,
  gamma3e1 = 0,
  gamma1e2 = 0,
  gamma2e2 = 0,
  gamma3e2 = 0,
  accrualDuration = 30.143,
  plannedEvents = c(186, 259, 183),
  maxNumberOfIterations = 1000,
  maxNumberOfRawDatasetsPerStage = 1,
  seed = 314159)

sum1 <- lr$sumdata
