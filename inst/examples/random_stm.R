# prepare random sparse matrix
library(slam)
library(rmumps)
n = 6000
a = simple_triplet_zero_matrix(n)
set.seed(7)
ij = sort(sample(1:(n*n), min(15*n, n*n)))
a$v = runif(ij)
a$i = as.integer((ij-1)%%n+1)
a$j = as.integer(floor((ij-1)/n)+1)
idi = a$i == a$j
a$v = a$v[!idi]; a$i = a$i[!idi]; a$j = a$j[!idi] #diag(a) = 0
a$v = c(a$v, -row_sums(a)); a$i=c(a$i, seq(n)); a$j=c(a$j, seq(n)) #diag(a) = -rowSums(a)
a[1,1] = as.vector(a[1,1]-1)
am = Rmumps$new(a)
b = as.double(matprod_simple_triplet_matrix(a, 1:n)) # rhs for an exact solution vector 1:n
# following time includes symbolic analysis, LU factorization and system solving
system.time(x<-am$solve(b))
bb = 2*b
# this second time should be much shorter
# as symbolic analysis and LU factorization are already done
system.time(xx<-am$solve(bb))
# see precision
range(x-1:n)  # mumps
# matrix inversion
system.time(aminv <- am$inv())

# symmetric matrix -> LDL^t decomposition
aa=as.simple_triplet_matrix(crossprod_simple_triplet_matrix(a))
aal=aa
aal$v[aal$j > aal$i]=0 # NB! upper part is zeroed, otherwise get a garbage in the solution vector
alm=Rmumps$new(aal, 1)
ba=as.double(matprod_simple_triplet_matrix(aa, 1:n)) # rhs for an exact solution vector 1:n
system.time(xl <- alm$solve(ba))
# the same result can be obtained with upper triangular matrix
aau=t(aal)
aum=Rmumps$new(aau, 1)
system.time(xu <- aum$solve(ba))

# clean up by hand to avoid a possible interference between gc() and
# Rcpp object destructor after unloading rmumps namespace
rm(am, alm, aum)
gc()
