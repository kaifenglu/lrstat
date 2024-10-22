library(dplyr, warn.conflicts = FALSE)

testthat::test_that("rmest: valid milestone time works", {
  df1 <- rmest(aml %>% filter(x == "Maintained"), time="time",
               event="status", milestone=161)

  # the following estimates can be obtained from SAS PROC LIFETEST
  #   proc lifetest data=aml rmst(tau=161);
  #     where x = "Maintained";
  #     time time*status(0);
  #   run;
  rmst = 52.64545
  stderr = 19.8286

  testthat::expect_equal(round(df1$rmst, 5), rmst)
  testthat::expect_equal(round(df1$stderr, 4), stderr)
})


testthat::test_that("rmest: bias correction for variance estimate", {
  df1 <- rmest(aml %>% filter(x == "Maintained"), time="time",
               event="status", milestone=161, biascorrection=TRUE)

  # the following estimates can be obtained from SAS PROC LIFETEST
  #   proc lifetest data=aml rmst(bc tau=161);
  #     where x = "Maintained";
  #     time time*status(0);
  #   run;
  stderr = 21.4173

  testthat::expect_equal(round(df1$stderr, 4), stderr)
})


testthat::test_that("rmest: milestone should be <= largest observed time", {
  testthat::expect_error(rmest(aml %>% filter(x == "Nonmaintained"),
                               time="time", event="status", milestone=161))
})
