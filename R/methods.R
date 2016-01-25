solve.Rcpp_Rmumps <- function(a,b) if(missing(b)) a$inv() else a$solve(b)
