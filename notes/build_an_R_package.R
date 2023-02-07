Rcpp::compileAttributes()
devtools::document()
devtools::build()

system("R CMD check --as-cran ../lrstat_0.1.10.tar.gz")

