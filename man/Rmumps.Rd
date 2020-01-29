\name{Rcpp_Rmumps-class}
\Rdversion{1.1}
\docType{class}
\alias{Rcpp_Rmumps-class}
\alias{Rmumps}
\alias{solve.Rcpp_Rmumps}
\alias{solvet}
\alias{solvet.Rcpp_Rmumps}
\alias{determinant.Rcpp_Rmumps}
\alias{dim.Rcpp_Rmumps}
\alias{nrow.Rcpp_Rmumps}
\alias{ncol.Rcpp_Rmumps}
\alias{print.Rcpp_Rmumps}
\alias{show.Rcpp_Rmumps}
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
  system.time(x<-solve(am, b))
  bb=2*b
  # this second time should be much shorter
  # as symbolic analysis and LU factorization are already done
  system.time(xx<-solve(am, bb))
  # compare to Matrix corresponding times
  system.time(xm<-solve(a, b))
  system.time(xxm<-solve(a, bb))
  # compare to Matrix precision
  range(x-1:n)  # mumps
  range(xm-1:n) # Matrix

  # matrix inversion
  system.time(aminv <- solve(am))
  system.time(ainv <- solve(a)) # the same in Matrix
  
  # symmetric matrix
  asy=as(a+t(a), "symmetricMatrix")
  bs=as.double(asy\%*\%(1:n)) # rhs for 1:n solution
  au=asy
  # Here, we keep only diagonal and upper values of asy matrix.
  # It could be also diagonal and lower values.
  au[row(au)>col(au)]=0
  ams=Rmumps$new(au, sym=1)
  system.time(xs<-solve(ams, bs)) # rmumps
  system.time(xsm<-solve(asy, bs))# Matrix
  # compare to Matrix precision
  range(xs-1:n)  # mumps
  range(xsm-1:n) # Matrix

  # clean up by hand to avoid possible interference between gc() and
  # Rcpp object destructor after unloading this namespace
  rm(am, ams)
  gc()
 }
}
\references{
 MUMPS official site http://mumps.enseeiht.fr

 Sokol S (2020). _Rmumps: Rcpp port of MUMPS_. rmumps package version 5.2.1-X, <URL:
 http://CRAN.R-project.org/package=rmumps>.
}
\author{
 Serguei Sokol, INRA
}
\note{
 When creating a symmetric matrix (sym=1 or sym=2), the upper (or lower) mart of the input matrix must be zeroed.

 For meaning of entries in MUMPS vectors cntl, icntl, info, rinfo, infog and rinfog cf. original documentation of MUMPS project.

 No need to call symbolic() and numeric() methods before a solve() call.

 If in constructor, a parameter \code{copy} is set to FALSE, no rhs neither matrix copying is done.
 The solution is written "in place" thus overwriting rhs (watch out side effects)

 For a detailed error diagnostic (e.g. when factorizing a singular matrix), use method get_infos() and cf. MUMPS documentation on the official MUMPS site).
}

\keyword{classes}
\section{Fields}{
  \describe{
    \item{\code{sym}:}{integer (read only), 0=non symmetric matrix, 1=symmetric with pivots on diagonal or 2=general symmetric}
    \item{\code{copy}:}{logical, copy or not rhs and matrix values}
    \item{\code{mrhs}:}{numeric matrix, multiple rhs (always overwritten with solution)}
    \item{\code{rhs}:}{numeric vector, single rhs (always overwritten with solution)}
  }
}
\section{Methods}{
  \describe{
    \item{\code{new(asp, sym=0, copy=TRUE)}:}{constructor from Matrix::dgTMatrix class (or from convertible to it) and slam::simple_triplet_matrix class}
    \item{\code{new(i, j, x, n, copy=TRUE)}:}{constructor from triade rows, cols, vals}
    \item{\code{symbolic()}:}{do symbolic analysis (stored internally)}
    \item{\code{numeric()}:}{do LU or LDL^t factorization (stored internally)}
    \item{\code{solve(b)}:}{solve single rhs (if b is a vector) or
     multiple rhs if b is a matrix (can be dense or sparse). Return the solution(s).}
    \item{\code{solvet(b)}:}{same as solve() but solves with transposed matrix}
    \item{\code{det()}:}{Return determinant of the matrix}
    \item{\code{inv()}:}{Return inverse of the matrix)}
    \item{\code{set_mat_data(x)}:}{updates matrix entries (x must be in the same order as in previous calls}
    \item{\code{set_icntl(iv, ii)}:}{set ICNTL parameter vector}
    \item{\code{get_icntl()}:}{get ICNTL parameter vector}
    \item{\code{set_cntl(v, iv)}:}{set CNTL parameter vector}
    \item{\code{get_cntl()}:}{get CNTL parameter vector}
    \item{\code{get_infos()}:}{get a named list of information vectors: info, rinfo, infog and rinfog}
    \item{\code{dim()}:}{Return a dimension vector of the matrix}
    \item{\code{nrow()}:}{Return a row number of the matrix}
    \item{\code{ncol()}:}{Return a column number of the matrix}
    \item{\code{print()}:}{Print summary information on the matrix}
    \item{\code{show()}:}{Print summary information on the matrix}
    \item{\code{set_keep()}:}{Set KEEP array elements (undocumented feature of MUMPS)}
    \item{\code{get_keep()}:}{Get a copy of KEEP array elements (length=500)}
    \item{\code{set_permutation(perm)}:}{Set permutation type which can impact storage and factorization performances. Parameter perm can take one of the following predefined integer values RMUMPS_PERM_AMD, RMUMPS_PERM_AMF, RMUMPS_PERM_SCOTCH, RMUMPS_PERM_PORD, RMUMPS_PERM_METIS, RMUMPS_PERM_QAMD. This method should be called once and before symbolic analysis of the matrix. If it is called afterward, a new symbolic and numeric factorization will be performed when one of other methods (e.g. solve()) will request them. In other words, previous symbolic and numeric factorizations are canceled by this method.}
    \item{\code{get_permutation()}:}{get permutation type currently set in the object}
    \item{\code{mumps_version()}:}{Return a string with MUMPS version used in rmumps}
  }
}
