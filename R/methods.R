solve.Rcpp_Rmumps <- function (a, b, ...) if(missing(b)) a$inv() else a$solve(b)
determinant.Rcpp_Rmumps <- function (x, logarithm=TRUE, ...) {d=x$det(); m=li$rinfog[12]; li=x$get_infos(); z=list(); z$sign=sign(m); z$modulus=if (logarithm) li$infog[34]*log(2.)+log(abs(m)) else abs(d); z}
dim.Rcpp_Rmumps <- function (x) x$dim()
nrow.Rcpp_Rmumps <- function (x) x$nrow()
ncol.Rcpp_Rmumps <- function (x) x$ncol()
print.Rcpp_Rmumps <- function (x, ...) x$print()
show.Rcpp_Rmumps <- function (x) x$show()

solvet <- function(a, b, ...) UseMethod("solvet")
solvet.Rcpp_Rmumps <- function (a, b, ...) if(missing(b)) base::t(a$inv()) else a$solvet(b)
