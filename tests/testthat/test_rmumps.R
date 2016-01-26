context("random matrix")
n=10
a=Matrix(0, n, n)
set.seed(7)
ij=sample(1:(n*n), 5*n)
a[ij]=runif(ij)
diag(a)=0
diag(a)=-Matrix::rowSums(a)
a[1,1]=a[1,1]-1
am=Rmumps$new(a)
b=as.double(a%*%(1:n)) # rhs for an exact solution vector 1:n
x<-am$solve(b)
test_that("first solving", {
  expect_equal(x, 1:n)
})
bb=2*b
xx<-solve(am, bb)
# test basic usage
test_that("second solving", {
  expect_equal(xx, 2*(1:n))
})

# test matrix creation from ijv triade
a=as(a, "dgTMatrix")
ai=Rmumps$new(a@i, a@j, a@x, ncol(a))
xi=solve(ai, b)
test_that("testing a from i,j,v", {
  expect_equal(xi, 1:n)
})

# test error signaling on singular matrix
rm(am)
a=Matrix(diag(n)); a[1,1]=0; a[1,2]=1
am=Rmumps$new(a)
test_that("singular matrix", {
  expect_error(solve(am, b), "*rmumps: info\\[1\\]=-10*")
})
