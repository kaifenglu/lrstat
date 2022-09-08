library(lrstat)
library(dplyr)
library(tidyr)
library(tictoc)

tic("lrsim2e3a");

sim1 = lrsim2e3a(
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

toc()


df <- sim1$sumdata %>% 
  mutate(p13 = pnorm(logRankStatistic13), 
         p23 = pnorm(logRankStatistic23), 
         p12 = pnorm(logRankStatistic12),
         events13 = events1 + events3,
         events23 = events2 + events3,
         events12 = events1 + events2) %>% 
  select(iterationNumber, stageNumber, endpoint, 
         events13, events23, events12, 
         p13, p23, p12)


dfcomp <- tibble(comparison = c("13", "23", "12"),
                 order = c(1, 2, 3))

dfinfo <- df %>% arrange(iterationNumber, endpoint) %>% 
  pivot_longer(c("events13", "events23", "events12"), 
               names_to = "comparison_events", values_to = "events") %>% 
  mutate(comparison = substring(comparison_events, 7)) %>%
  left_join(dfcomp, by = c("comparison")) %>% 
  arrange(iterationNumber, order, endpoint, stageNumber)

dfp <- df %>% arrange(iterationNumber, endpoint) %>% 
  pivot_longer(c("p13", "p23", "p12"), 
               names_to = "comparison_p", values_to = "p") %>% 
  mutate(comparison = substring(comparison_p, 2)) %>%
  left_join(dfcomp, by = c("comparison")) %>% 
  arrange(iterationNumber, order, endpoint, stageNumber)


# strategy 1

tic("feqbon - strategy 1")

reject1 = fseqbon(
  w = c(0.2, 0.8, 0, 0, 0, 0), 
  G = matrix(c(0, 1, 0, 0, 0, 0,  
               0, 0, 0.2, 0.8, 0, 0, 
               0, 0, 0, 1, 0, 0,  
               0, 0, 0, 0, 0.2, 0.8,
               0, 0, 0, 0, 0, 1,  
               0, 0, 0, 0, 1, 0), 
             nrow=6, byrow=TRUE), 
  alpha = 0.025,
  asf = rep("sfHSD", 6),
  asfpar = rep(-4, 6),
  incid = matrix(c(1, 1, 0, 
                   1, 1, 1,
                   1, 1, 1,
                   1, 1, 1,
                   1, 1, 1,
                   1, 1, 1), 
                 nrow=6, ncol=3, byrow=TRUE),
  maxinfo = c(259, 183, 321, 195, 392, 230),
  info = matrix(dfinfo$events, 6000, 3, byrow=T),
  p = matrix(dfp$p, 6000, 3, byrow=T))

toc()

reject1a = matrix(rowSums(reject1, na.rm=TRUE), nrow=1000, ncol=6, byrow=TRUE)
(power1 = apply(reject1a, 2, mean))

# strategy 2

tic("feqbon - strategy 2")

eps = 1e-8
reject2 = fseqbon(
  w = c(0.1, 0.4, 0.1, 0.4, 0, 0), 
  G = matrix(c(0, 1, 0, 0, 0, 0,  
               0, 0, 0.2, 0.8, 0, 0, 
               0, 0, 0, 1, 0, 0,  
               0.2*(1-eps), 0.8*(1-eps), 0, 0, 0.2*eps, 0.8*eps,
               0, 0, 0, 0, 0, 1,  
               0, 0, 0, 0, 1, 0), 
             nrow=6, byrow=TRUE), 
  alpha = 0.025,
  asf = rep("sfHSD", 6),
  asfpar = rep(-4, 6),
  incid = matrix(c(1, 1, 0, 
                   1, 1, 1,
                   1, 1, 1,
                   1, 1, 1,
                   1, 1, 1,
                   1, 1, 1), 
                 nrow=6, ncol=3, byrow=TRUE),
  maxinfo = c(259, 183, 321, 195, 392, 230),
  info = matrix(dfinfo$events, 6000, 3, byrow=T),
  p = matrix(dfp$p, 6000, 3, byrow=T))

toc()

reject2a = matrix(rowSums(reject2, na.rm=TRUE), nrow=1000, ncol=6, byrow=TRUE)
(power2 = apply(reject2a, 2, mean))



# sample size calculation for PFS Arm A vs. Arm C
(lr1 <- lrsamplesize(
  beta = 0.1, kMax = 2, informationRates = c(186/259, 1), 
  alpha = 0.005, typeAlphaSpending = "sfHSD", parameterAlphaSpending = -4, 
  allocationRatioPlanned = 2, 
  accrualTime = c(0,8), accrualIntensity = c(10, 28)*3/5, 
  lambda1 = log(2)/12*0.60, lambda2 = log(2)/12, 
  accrualDuration = 30.143, followupTime = NA, 
  typeOfComputation = "Schoenfeld"))

# sample size calculation for OS Arm A vs. Arm C
lr2a <- lrstat(
  time = lr1$byStageResults$analysisTime, 
  allocationRatioPlanned = 2,
  accrualTime = c(0,8), accrualIntensity = c(10, 28)*3/5, 
  lambda1 = log(2)/30*0.65, lambda2 = log(2)/30,
  accrualDuration = 30.143, followupTime = 100,
  predictEventOnly = 1)

(lr2 <- lrsamplesize(
  beta = 0.25, kMax = 3, 
  informationRates = c(lr2a$nevents, 183)/183, 
  alpha = 0.020, typeAlphaSpending = "sfHSD", parameterAlphaSpending = -4, 
  allocationRatioPlanned = 2, 
  accrualTime = c(0,8), accrualIntensity = c(10, 28)*3/5, 
  lambda1 = log(2)/30*0.65, lambda2 = log(2)/30,
  accrualDuration = 30.143, followupTime = NA, 
  typeOfComputation = "Schoenfeld"))
