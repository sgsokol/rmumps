#solve.Rcpp_Rmumps <- function(a,b) if(missing(b)) a$inv() else a$solve(b)
setMethod("solve",
   signature(a = "Rcpp_Rmumps", b="ANY"),
   function (a, b)
   {
      if(missing(b)) a$inv() else a$solve(b)
   }
)
