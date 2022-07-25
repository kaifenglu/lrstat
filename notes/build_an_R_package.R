Rcpp::compileAttributes()
devtools::document()

devtools::test()

devtools::build()

system("R CMD check --as-cran ../lrstat_0.1.6.tar.gz")

