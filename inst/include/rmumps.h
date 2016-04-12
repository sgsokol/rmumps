#ifndef __rmumps_h__
#define __rmumps_h__

#include <dmumps_c.h>
#include <Rcpp.h>
using namespace Rcpp;
class Rmumps {
public:
  std::vector<MUMPS_INT> irn;
  std::vector<MUMPS_INT> jcn;
  std::vector<MUMPS_INT> irhs_ptr;
  std::vector<MUMPS_INT> irhs_sparse;
  NumericVector rhs;
  NumericMatrix mrhs; // many rhs
  NumericVector rhs_sparse;
  NumericVector anz;
  bool copy;
  std::set<int> jobs;
  char buf[512];
  
  Rmumps(RObject mat);
  Rmumps(RObject mat, bool copy_);
  Rmumps(IntegerVector i, IntegerVector j, NumericVector x, int n);
  Rmumps(IntegerVector i, IntegerVector j, NumericVector x, int n, bool copy_);
  ~Rmumps();
  void clean();
  bool get_copy();
  void set_copy(bool copy_);
  void do_job(int job);
  void symbolic();
  void numeric();
  SEXP solve(RObject b);
  NumericVector solvev(NumericVector b);
  NumericMatrix solvem(NumericMatrix b);
  void solveptr(double* b, int nrow, int nrhs);
  NumericMatrix inv();
  NumericMatrix solves(S4 mat);
  NumericMatrix solvestm(List mat);
  void set_rhs(NumericVector b);
  void set_mrhs(NumericMatrix b);
  void set_mat_data(NumericVector x);
  NumericVector get_rhs();
  NumericMatrix get_mrhs();
  IntegerVector dim();
  int nrow();
  int ncol();
  void print();
  List triplet();
  std::string mumps_version();
private:
  DMUMPS_STRUC_C param;
  void new_mat(RObject mat, bool copy_);
  void new_ijv(IntegerVector i0, IntegerVector j0, NumericVector x, int n, bool copy_);
  void tri_init(MUMPS_INT *irn, MUMPS_INT *jcn, double *a);
};
#include "rmumps_RcppExports.h"

#endif // __rmumps_h__
