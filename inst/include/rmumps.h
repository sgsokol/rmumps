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
  int ncore;
  std::set<int> jobs;
  MUMPS_INT sym;
  char buf[512];
  
  Rmumps(RObject mat);
  Rmumps(RObject mat, int sym);
  Rmumps(RObject mat, int sym, bool copy_);
  Rmumps(IntegerVector i, IntegerVector j, NumericVector x, int n);
  Rmumps(IntegerVector i, IntegerVector j, NumericVector x, int n, int sym);
  Rmumps(IntegerVector i, IntegerVector j, NumericVector x, int n, int sym, bool copy_);
  ~Rmumps();
  void clean();
  bool get_copy();
  void set_copy(bool copy_);
  int get_ncore();
  void set_ncore(int ncore_);
  int get_sym();
  void set_sym(int sym);
  void do_job(int job);
  void symbolic();
  void numeric();
  SEXP solve(RObject b);
  SEXP solvet(RObject b);
  NumericVector solvev(NumericVector b);
  NumericMatrix solvem(NumericMatrix b);
  void solveptr(double* b, int nrow, int nrhs);
  NumericMatrix inv();
  NumericMatrix solves(S4 mat);
  NumericMatrix solvestm(List mat);
  void set_rhs(NumericVector b);
  void set_mrhs(NumericMatrix b);
  void set_mat_data(NumericVector x);
  void set_icntl(IntegerVector iv, IntegerVector ii);
  void set_cntl(NumericVector v, IntegerVector iv);
  IntegerVector get_icntl();
  NumericVector get_cntl();
  List get_infos();
  NumericVector get_rhs();
  NumericMatrix get_mrhs();
  IntegerVector dim();
  int nrow();
  int ncol();
  void print();
  List triplet();
  std::string mumps_version();
  double det();
private:
  DMUMPS_STRUC_C param;
  void new_mat(RObject mat, int sym, bool copy_);
  void new_ijv(IntegerVector i0, IntegerVector j0, NumericVector x, int n, int sym, bool copy_);
  void tri_init(MUMPS_INT *irn, MUMPS_INT *jcn, double *a, MUMPS_INT sym);
};
//#include "rmumps_RcppExports.h"

#endif // __rmumps_h__
