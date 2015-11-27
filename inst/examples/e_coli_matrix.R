rm(list=ls(all=TRUE)); unloadNamespace("rmumps"); gc()
require(Matrix)
#require(rmumps);
library(devtools); load_all("../..");
require(microbenchmark)

load("e_coli_ab0.RData")
a=as(lAb$A, "dgTMatrix")
b=lAb$b
bd=as.numeric(b)
n=ncol(a) # matrix size n x n

ad=as.matrix(a)

# define solution(s)
xsol=solve(a, b)

# init mumps object
am=Rmumps$new(a)
microbenchmark(am$solve(b), times=1)
#b2=as.matrix(b+1)
# new rhs but use LU already done
#microbenchmark(am$mrhs <- b2, times=1)
#microbenchmark(x2<-am$do_job(3), times=1)
microbenchmark(am$solve(as.numeric(b+1)), times=1)


# solve
microbenchmark(solve(ad, bd), solve(a, b), am$solves(b), times=1)# times > 1 is a pb as the solution is put in b (in place)
stopifnot(diff(range(am$solves(b)-xsol)) < 1.e-10)
