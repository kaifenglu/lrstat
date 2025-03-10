testthat::test_that("power with alpha-spending, weighted", {
  l = lrpower(kMax = 2, informationRates = c(0.8, 1),
              alpha = 0.025, typeAlphaSpending = "sfOF",
              allocationRatioPlanned = 1, accrualTime = seq(0, 8),
              accrualIntensity = 26/9*seq(1, 9),
              piecewiseSurvivalTime = c(0, 6),
              stratumFraction = c(0.2, 0.8),
              lambda1 = c(0.0533, 0.0309, 1.5*0.0533, 1.5*0.0309),
              lambda2 = c(0.0533, 0.0533, 1.5*0.0533, 1.5*0.0533),
              gamma1 = -log(1-0.05)/12,
              gamma2 = -log(1-0.05)/12, accrualDuration = 22,
              followupTime = 18, fixedFollowup = FALSE,
              rho1 = 0, rho2 = 1)
  testthat::expect_equal(round(l$overallResults$overallReject, 4), 0.9313)
})


testthat::test_that("power for stratified analysis", {
  p1 = c(0.28, 0.13, 0.25, 0.34)
  p2 = c(0.28, 0.72)
  p3 = c(0.43, 0.37, 0.2)
  stratumFraction = p1 %x% p2 %x% p3
  stratumFraction = stratumFraction/sum(stratumFraction)
  theta1 = c(1, 2.127, 0.528, 0.413)
  theta2 = c(1, 0.438)
  theta3 = c(1, 0.614, 0.159)
  lambda2 = 0.009211*exp(log(theta1) %x% log(theta2) %x% log(theta3))
  caltime(nevents = 66, accrualDuration = 24, accrualIntensity = 12,
          stratumFraction = stratumFraction,
          lambda1 = 0.4466*lambda2, lambda2 = lambda2,
          followupTime = 100)
  l = lrpower(kMax = 3,
              informationRates = (1:3)/3,
              alpha = 0.025, typeAlphaSpending = "sfOF",
              accrualIntensity = 12,
              stratumFraction = stratumFraction,
              lambda1 = 0.4466*lambda2,
              lambda2 = lambda2,
              accrualDuration = 24,
              followupTime = 30.92)
  testthat::expect_equal(round(l$overallResults$overallReject, 4), 0.8824)
})
