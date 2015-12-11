\name{Rcpp_Rmumps-class}
\Rdversion{1.1}
\docType{class}
\alias{Rcpp_Rmumps-class}
\alias{Rmumps}
\title{
Rcpp Exported Class Wrapping MUMPS library
}
\description{
 This class can be used for storing sparse matrix
 and solving corresponding linear system with one or many right
 hand sides. There is a possibility to do separately symbolic analysis,
 LU factorization and system solving.
}
\examples{
 \dontrun{
  # prepare random sparse matrix
  library(Matrix)
  library(rmumps)
  n=2000
  a=Matrix(0, n, n)
  set.seed(7)
  ij=sample(1:(n*n), 15*n)
  a[ij]=runif(ij)
  diag(a)=0
  diag(a)=-rowSums(a)
  a[1,1]=a[1,1]-1
  am=Rmumps$new(a)
  b=as.double(a\%*\%(1:n)) # rhs for an exact solution vector 1:n
  # following time includes symbolic analysis, LU factorization and system solving
  system.time(x<-am$solve(b))
  bb=2*b
  # this second time should be much shorter
  # as symbolic analysis and LU factorization are already done
  system.time(xx<-am$solve(bb))
  # compare to Matrix corresponding times
  system.time(xm<-solve(a, b))
  system.time(xxm<-solve(a, bb))
  # compare to Matrix precision
  range(x-1:n)  # mumps
  range(xm-1:n) # Matrix

  # matrix inversion
  system.time(aminv <- am$inv())
  system.time(ainv <- solve(a)) # the same in Matrix

  # clean up by hand to avoid possible interference between gc() and
  # Rcpp object destructor after unloading this namespace
  rm(am)
 }
}
\references{
 MUMPS official site http://mumps.enseeiht.fr

 Sokol S (2015). _Rmumps: Rcpp port of MUMPS_. rmumps package version 5.0.1_1, <URL:
 http://CRAN.R-project.org/package=rmumps>.
}
\author{
 Serguei Sokol, INRA
}
\note{
 No need to call symbolic() and numeric() methods before a solve() call.

 If in constructor, parameter \code{copy} is set to FALSE, no rhs neither matrix copying is done.
 The solution is written "in place" thus overwriting rhs (watch out side effects)

 For a possible error messages (e.g. when factorizing a singular matrix), cf.
 MUMPS documentation on the official MUMPS cite).
}

\keyword{classes}
\section{Fields}{
  \describe{
    \item{\code{copy}:}{logical, copy or not rhs and matrix values}
    \item{\code{mrhs}:}{numeric matrix, multiple rhs (always overwritten with solution)}
    \item{\code{rhs}:}{numeric vector, single rhs (always overwritten with solution)}
  }
}
\section{Methods}{
  \describe{
    \item{\code{new(asp, copy=TRUE)}:}{constructor from dgTMatrix class (or from a matrix convertable to it)}
    \item{\code{new(i, j, x, copy=TRUE)}:}{constructor from triade rows, cols, vals}
    \item{\code{symbolic()}:}{do symbolic analysis (stored internally) }
    \item{\code{numeric()}:}{do LU factorization (stored internally)}
    \item{\code{solve(b)}:}{solve single rhs (if b is a vector) or
     multiple rhs if b is a matrix (can be dense or sparse). Return the solution(s).}
    \item{\code{inv()}:}{Return inverse of the matrix)}
    \item{\code{set_mat_data(x)}:}{updates matrix entries (same order as previously)}
  }
}
\keyword{classes}
\section{Fields}{
  \describe{
    \item{\code{copy}:}{Object of class \code{activeBindingFunction} ~~ }
    \item{\code{mrhs}:}{Object of class \code{activeBindingFunction} ~~ }
    \item{\code{rhs}:}{Object of class \code{activeBindingFunction} ~~ }
  }
}
