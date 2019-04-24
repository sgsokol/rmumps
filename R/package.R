#' rmumps
#'
#' Rcpp wrapper for MUMPS library factorizing and solving linear systems with sparse matrices
#'
#' @docType package
#' @author Serguei Sokol
#' @import Rcpp RcppArmadillo
#' @importFrom Rcpp sourceCpp
#' @useDynLib rmumps
#' @name rmumps
NULL
# set useful constants
cnsts=c(
   "RMUMPS_PERM_AMD",
   "RMUMPS_PERM_AMF",
   "RMUMPS_PERM_SCOTCH",
   "RMUMPS_PERM_PORD",
   "RMUMPS_PERM_METIS",
   "RMUMPS_PERM_QAMD",
   "RMUMPS_PERM_AUTO"
)
