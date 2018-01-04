\name{rmumps-package}
\alias{rmumps-package}
\alias{rmumps}
\docType{package}
\title{
  Rcpp port of MUMPS library for LU factorization of sparse matrices
}
\description{
  Creates a MUMPS compatible object storing a sparse matrix.
  Gives a possibility to do separately symbolic analysis,
  LU factorization and system solving.
}
\details{

  The DESCRIPTION file:
  Create a new Rmumps object with \code{A <- Rmumps$new(asparse)} then solve a
  linear system with one or many right hand sides \code{x <- solve(A, b)}
}
\author{
  Serguei Sokol, INRA

  Maintainer: Serguei Sokol (sokol at insa-toulouse.fr)
}
\references{
  MUMPS official site http://mumps.enseeiht.fr

  Sokol S (2018). _Rmumps: Rcpp port of MUMPS_. rmumps package version 5.1.2-3, <URL: http://CRAN.R-project.org/package=rmumps>.
}
\examples{
  \dontrun{
    A <- Rmumps$new(asparse)
    x <- solve(A, b)
  }
}
