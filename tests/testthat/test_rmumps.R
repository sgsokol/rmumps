context("random matrix")
n=10
a=Matrix::Matrix(0, n, n)
set.seed(7)
ij=sample(1:(n*n), 5*n)
a[ij]=runif(ij)
diag(a)=0
diag(a)=-Matrix::rowSums(a)
a[1,1]=a[1,1]-1
am=Rmumps$new(a)
b=as.double(a%*%(1:n)) # rhs for an exact solution vector 1:n
x=am$solve(b)
test_that("first solving", {
  expect_equal(x, 1:n)
})
bb=2*b
xx=solve(am, bb)
# test basic usage
test_that("second solving", {
  expect_equal(xx, 2*(1:n))
})
# test solving transposed system
bt=as.double(Matrix::crossprod(a, 1:n))
xt=solvet(am, bt)
test_that("solving transposed system", {
  expect_equal(xt, 1:n)
})

# test symmetric matrix solving
asy=a+Matrix::t(a)
bsy=as.double(asy%*%(1:n))
asl=asy
asl[row(asy)>col(asy)]=0.
ams=Rmumps$new(asl, sym=1)
xsy=solve(ams, bsy)
test_that("solving symmetric-1 system", {
  expect_equal(xsy, 1:n)
})
ams=Rmumps$new(asl, sym=2)
xsy=solve(ams, bsy)
test_that("solving symmetric-2 system", {
  expect_equal(xsy, 1:n)
})

# test matrix creation from ijv triade
a=as(a, "dgTMatrix")
ai=Rmumps$new(a@i, a@j, a@x, ncol(a))
ai$set_permutation(RMUMPS_PERM_SCOTCH)
xi=solve(ai, b)
test_that("testing a from i,j,v", {
  expect_equal(xi, 1:n)
})
rm(ai)

if (getRversion() >= "3.4.0") {
  # test matrix creation from slam::simple_triplet_matrix (which requires R-3.4.0+)
  asl=as.simple_triplet_matrix(a)
  ai=Rmumps$new(asl)
  xi=solve(ai, b)
  test_that("testing matrix from slam::simple_triplet_matrix", {
    expect_equal(xi, 1:n)
  })
  # test pointers
  code='
NumericVector solve_ptr(List a, NumericVector b) {
  IntegerVector ir=a["i"], jc=a["j"];
  NumericVector v=a["v"];
  int n=a["nrow"], nz=v.size();
  Rmumps rmu((int *) ir.begin(), (int *) jc.begin(), (double *) v.begin(), n, nz, 0);
  rmu.set_permutation(RMUMPS_PERM_SCOTCH);
  rmu.solveptr((double *) b.begin(), n, 1);
  return(b);
}
'
#cat("ls -R rmumps_path", list.files(path.package("rmumps"), recursive=TRUE), "", sep="\n")
  if (FALSE && Sys.info()[["sysname"]] != "Darwin" && .Platform$r_arch != "i386") {
    rso=paste0("rmumps", .Platform$dynlib.ext)
    rso_path=file.path(path.package("rmumps"), "libs", .Platform$r_arch, rso)
    cat("rso_path=", rso_path, "\n")
    if (!file.exists(rso_path))
      rso_path=file.path(path.package("rmumps"), "src", rso) # devtool context
#cat("rso_path='", rso_path, "'\n", sep="")
    Sys.setenv(PKG_LIBS=rso_path)
    cppFunction(code=code, depends="rmumps", verbose=TRUE)
    #sourceCpp(code=code, verbose=TRUE)
    xe=as.double(1:n)
    b0=slam::tcrossprod_simple_triplet_matrix(asl, t(xe))
    x=solve_ptr(asl, b0)
    test_that("testing solveptr() within Rcpp code", {
      expect_lt(diff(range(x-xe)), 1e-14)
    })
  }

  rm(asl)

  # test sparse rhs as slam::simple_triplet_matrix
  asl=slam::as.simple_triplet_matrix(a)
  eye=solve(ai, asl)
  test_that("testing solve() on slam::simple_triplet_matrix", {
    expect_equal(eye, diag(n), tol=1e-14)
  })
  rm(ai, asl)
}


# test error signaling on singular matrix
rm(am)
a=as(diag(n), "dgCMatrix")
a@p[2L]=0L
am=Rmumps$new(a)
test_that("singular matrix", {
  expect_error(solve(am, b), "rmumps: job=6, info\\[1\\]=-10*")
})

# test int size in KEEP
vkeep=am$get_keep()
sizeint=Rcpp::evalCpp("sizeof(int)")
test_that("int size", {
  expect_equal(vkeep[34], sizeint)
  expect_equal(vkeep[10], 8/sizeint)
})


rm(a, asy, am)
gc()

