\name{rmumps-package}
\alias{rmumps-package}
\alias{rmumps}
\docType{package}
\title{
  Rcpp port of MUMPS library for LU or LDL^t factorization of sparse matrices
}
\description{
  Creates a MUMPS compatible object storing a sparse matrix.
  Gives a possibility to do separately symbolic analysis,
  factorization and system solving.
}
\details{
  Create a new Rmumps object with \code{A <- Rmumps$new(asparse)} then solve a
  linear system with one or many right hand sides \code{x <- solve(A, b)}
}
\author{
  Serguei Sokol, INRA

  Maintainer: Serguei Sokol (sokol at insa-toulouse.fr)
}
\references{
  MUMPS official site http://mumps.enseeiht.fr

  Sokol S (2020). _Rmumps: Rcpp port of MUMPS_. rmumps package version 5.2.1-8, <URL: http://CRAN.R-project.org/package=rmumps>.
}
\examples{
  \dontrun{
    A <- Rmumps$new(asparse)
    x <- solve(A, b)
  }
}
