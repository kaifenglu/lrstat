library(survival)

testthat::test_that("rmdiff: unstratified test of rmst difference", {
  df1 <- rmdiff(veteran, treat="trt", time="time",
                event="status", milestone=90)

  # the following values are obtained from SAS PROC LIFETEST:
  #   proc lifetest data=veteran rmst(tau=90);
  #     time time*status(0);
  #     strata trt;
  #   run;
  rmst1 = 62.75618
  rmst2 = 56.45116
  stderr1 = 4.0321
  stderr2 = 3.9790
  rmstdiffchisq = 1.2388
  pvalue = 0.2657

  testthat::expect_equal(c(round(df1$rmst1, 5),
                           round(df1$rmst2, 5),
                           round(sqrt(df1$vrmst1), 4),
                           round(sqrt(df1$vrmst2), 4),
                           round(df1$rmstDiffZ^2, 4),
                           round(df1$rmstDiffPValue, 4)),
                         c(rmst1, rmst2, stderr1, stderr2,
                           rmstdiffchisq, pvalue))
})


testthat::test_that("rmdiff: stratified test of rmst difference", {
  df1 <- rmdiff(veteran, stratum="celltype", treat="trt",
                time="time", event="status", milestone=90)

  # Of note, the stratified results are different from SAS PROC LIFETEST:
  #   proc lifetest data=veteran rmst(tau=90);
  #     time time*status(0);
  #     strata celltype / group = trt;
  #   run;
  # This is because SAS adds up the rmst diffs and variances across strata,
  # while we use the number of subjects as the weight for each stratum.
  # Our approach yields a more meaningful rmst diff across strata.
  # To reproduce our results, we combine the stratum-specific information
  # from SAS PROC LIFETEST using the sample size weights across strata
  n = c(15, 20, 30, 18, 9, 18, 15, 12)
  rmst = c(68.55152, 65.95, 50.86667, 45.16667, 56.44444, 51.91667,
           84.8, 64.25)
  stderr = c(8.3317, 7.5975, 5.7978, 7.8103, 12.7986, 6.9569,
             5.0237, 8.0847)
  a = c(1, 3, 5, 7)
  b = c(2, 4, 6, 8)
  ns = n[a] + n[b]
  rmstDiffs = rmst[a] - rmst[b]
  vrmstDiffs = stderr[a]^2 + stderr[b]^2
  w = ns/sum(ns)
  rmstDiff = sum(w*rmstDiffs)
  sermstDiff = sqrt(sum(w*w*vrmstDiffs))
  rmstDiffZ = rmstDiff/sermstDiff

  testthat::expect_equal(c(round(df1$rmstDiff, 4),
                           round(df1$sermstDiff, 3),
                           round(df1$rmstDiffZ, 3)),
                         c(round(rmstDiff, 4),
                           round(sermstDiff, 3),
                           round(rmstDiffZ, 3)))
})
