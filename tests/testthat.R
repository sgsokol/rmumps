library(testthat)
library(Matrix)
library(methods)
library(Rcpp)
if (getRversion() >= "3.4.0") {
  library(slam)
}
library(rmumps)
test_check("rmumps")
