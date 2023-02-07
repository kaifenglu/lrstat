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
  lambda1e1 = log(2)/22,
  lambda2e1 = log(2)/17,
  lambda3e1 = log(2)/12,
  lambda1e2 = log(2)/50, 
  lambda2e2 = log(2)/40,
  lambda3e2 = log(2)/30,
  gamma1e1 = 0,
  gamma2e1 = 0,
  gamma3e1 = 0,
  gamma1e2 = 0,
  gamma2e2 = 0,
  gamma3e2 = 0,
  accrualDuration = 30.143,
  plannedEvents = c(185, 247, 195),
  maxNumberOfIterations = 1000,
  maxNumberOfRawDatasetsPerStage = 1,
  seed = 314159)

toc()


df <- sim1$sumdata %>% 
  mutate(events13 = events1 + events3,
         events23 = events2 + events3,
         events12 = events1 + events2,
         p13 = pnorm(logRankStatistic13), 
         p23 = pnorm(logRankStatistic23), 
         p12 = pnorm(logRankStatistic12)) %>% 
  select(iterationNumber, stageNumber, endpoint, 
         events13, events23, events12, 
         p13, p23, p12)


dfcomp <- tibble(comparison = c("13", "23", "12"),
                 order = c(1, 2, 3))

dfinfo <- df %>% 
  arrange(iterationNumber, endpoint) %>% 
  pivot_longer(
    c("events13", "events23", "events12"), 
    names_to = "c_events", 
    values_to = "events") %>% 
  mutate(comparison = substring(c_events, 7)) %>%
  left_join(dfcomp, by = c("comparison")) %>% 
  arrange(iterationNumber, order, endpoint, 
          stageNumber)

dfp <- df %>% 
  arrange(iterationNumber, endpoint) %>% 
  pivot_longer(
    c("p13", "p23", "p12"), 
    names_to = "c_p", 
    values_to = "p") %>% 
  mutate(comparison = substring(c_p, 2)) %>%
  left_join(dfcomp, by = c("comparison")) %>% 
  arrange(iterationNumber, order, endpoint, 
          stageNumber)


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
             nrow=6, ncol=6, byrow=TRUE), 
  alpha = 0.025,
  kMax = 3,
  typeAlphaSpending = rep("sfHSD", 6),
  parameterAlphaSpending = rep(-4, 6),
  incidenceMatrix = matrix(
    c(1,1,0, 1,1,1, 1,1,1, 1,1,1, 1,1,1, 1,1,1), 
    nrow=6, ncol=3, byrow=TRUE),
  maxInformation = c(247, 195, 341, 215, 411, 248),
  p = matrix(dfp$p, 6000, 3, byrow=TRUE),
  information = matrix(dfinfo$events, 6000, 3, byrow=TRUE))

toc()

reject1a = matrix(reject1>0, nrow=1000, ncol=6, byrow=TRUE)
(power1 = apply(reject1a, 2, mean))

# strategy 2

tic("feqbon - strategy 2")

reject2 = fseqbon(
  w = c(0.2, 0.8, 0, 0, 0, 0), 
  G = matrix(c(0, 0.8, 0.2, 0, 0, 0, 
               0.5, 0, 0, 0.5, 0, 0,
               0, 0.8, 0, 0, 0.2, 0, 
               0.72, 0, 0, 0, 0, 0.28,
               0, 1, 0, 0, 0, 0, 
               1, 0, 0, 0, 0, 0), 
             nrow=6, ncol=6, byrow=TRUE),
  alpha = 0.025,
  kMax = 3,
  typeAlphaSpending = rep("sfHSD", 6),
  parameterAlphaSpending = rep(-4, 6),
  incidenceMatrix = matrix(
    c(1,1,0, 1,1,1, 1,1,1, 1,1,1, 1,1,1, 1,1,1), 
    nrow=6, ncol=3, byrow=TRUE),
  maxInformation = c(247, 195, 341, 215, 411, 248),
  p = matrix(dfp$p, 6000, 3, byrow=TRUE),
  information = matrix(dfinfo$events, 6000, 3, byrow=TRUE))

toc()

reject2a = matrix(reject2>0, nrow=1000, ncol=6, byrow=TRUE)
(power2 = apply(reject2a, 2, mean))

