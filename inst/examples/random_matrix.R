rm(list=ls(all=TRUE)); unloadNamespace("rmumps"); gc()
require(Matrix)
#require(rmumps);
library(devtools); load_all();
require(microbenchmark)
set.seed(7)

n=5 # matrix size n x n
nsq=n*n
nz=round(0.5*nsq) # how many non zero
ax=-runif(nz)
ai=sample(1:n, nz, replace=T)
aj=sample(1:n, nz, replace=T)
# make sparse matrix from triplicate x, i, j
a=sparseMatrix(i=ai, j=aj, x=ax, dims=c(n, n), giveCsparse=FALSE)
ad=as.matrix(a)

# make a diagonal dominant
diag(a)=0
diag(a)=-rowSums(a)
a[1,1]=a[1,1]+1

# define solution(s)
# one rhs
xsol=rnorm(n)
b=as.numeric(a%*%xsol)

# many rhs
nrhs=round(n/2)
msol=matrix(rnorm(n*nrhs), nrow=n)
bm=as.matrix(a%*%msol)

# init mumps object
am=Rmumps$new(a)

# solve one rhs
microbenchmark(solve(ad, b), solve(a, b), am$solve(b), times=1)# times > 1 is a pb as the solution is put in b (in place)
stopifnot(diff(range(b-xsol)) < 1.e-10)

# solve many rhs
microbenchmark(solve(ad, bm), solve(a, bm), am$solvem(bm), times=1)
stopifnot(diff(range(bm-msol)) < 1.e-10)

# effect of pre-symbolic analysis
b=as.numeric(a%*%xsol)
am=Rmumps$new(a)
am$symbolic()
# solve one rhs
microbenchmark(solve(ad, b), solve(a, b), am$solve(b), times=1)# times > 1 is a pb as the solution is put in b (in place)
stopifnot(diff(range(b-xsol)) < 1.e-10)

# effect of pre-symbolic+numeric analysis
b=as.numeric(a%*%xsol)
am=Rmumps$new(a)
am$numeric()
# solve one rhs
microbenchmark(solve(ad, b), solve(a, b), am$solve(b), times=1)# times > 1 is a pb as the solution is put in b (in place)
stopifnot(diff(range(b-xsol)) < 1.e-10)
