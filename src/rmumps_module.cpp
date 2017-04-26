// Port of sequential MUMPS to R.
// Only a subset of MUMPS is exposed to R.
// Author of port: Serguei Sokol, INRA, Toulouse, FRANCE
// Copyright 2015, INRA
// v 0.1 2015-09-02

#include <Rcpp.h>
using namespace Rcpp;

//#include <libseq/mpi.h> // useless in sequential mode
#define JOB_INIT -1
#define JOB_END -2
#define USE_COMM_WORLD -987654
#define ICNTL(I) icntl[(I)-1] /* macro s.t. indices match documentation */

#ifdef _OPENMP
#include <omp.h>
#endif

#include "../inst/include/rmumps.h"
Rmumps::~Rmumps() { if (ref == 0) clean(); }
void Rmumps::clean() {
  //Rprintf("clean() is called...\n");
  param.job=JOB_END;
  dmumps_c(&param);
  //Rprintf("clean(): done\n");
}
/* unfortunately, there is an unwanted extra destructor call
Rmumps Rmumps::shallow_copy(Rmumps a) {
  Rmumps x=a;
  x.ref++;
  return x;
}
* */
bool Rmumps::get_copy() {
  return copy;
}
void Rmumps::set_copy(bool copy_) {
  copy=copy_;
}
int Rmumps::get_ncore() {
  return ncore;
}
void Rmumps::set_ncore(int ncore_) {
#ifdef _OPENMP
  ncore=ncore_;
#else
  ncore=1;
#endif
}
int Rmumps::get_sym() {
  return param.sym;
}
void Rmumps::set_sym(int sym) {
  stop("Parameter 'sym' cannot be changed after matrix initialization");
}
void Rmumps::do_job(int job) {
  // downgrade job request if preliminaries are already done 
  switch (job) {
  case 4:
    if (jobs.count(1) == 1) {
       job=2;
    }
    break;
  case 5:
  case 6:
    if (jobs.count(2) == 1) {
       job=3;
    }
    break;
  }
  // do preliminary work if need
  switch (job) {
  case 2:
    if (jobs.count(1) != 1) {
       do_job(1);
    }
    break;
  case 3:
    if (jobs.count(1) != 1) {
       do_job(2);
    }
    break;
  }
  param.job=job;
  // set core number to use in openmp
//#ifdef _OPENMP
//#pragma omp parallel
  //{
  //omp_set_num_threads(ncore);
//#endif
  dmumps_c(&param);
//#ifdef _OPENMP
//  }
//#endif
  if (param.info[0] != 0) {
    //clean();
    stop("rmumps: info[1]=%d, info[2]=%d", param.info[0], param.info[1]);
  }
  /* combined jobs are split for the record */
  switch (job) {
  case 1:
  case 2:
    jobs.insert({job});
    break;
  case 3:
    jobs.insert({job});
    break;
  case 4:
    jobs.insert({1, 2});
    break;
  case 5:
    jobs.insert({2, 3});
    break;
  case 6:
    jobs.insert({1, 2, 3});
    break;
  }
}
void Rmumps::symbolic() {
  do_job(1);
}
void Rmumps::numeric() {
  do_job(4);
}
SEXP Rmumps::solvet(RObject b) {
  // solve A^t x = b
  param.ICNTL(9) = 2;
  SEXP res=solve(b);
  param.ICNTL(9) = 1;
  return res;
}
SEXP Rmumps::solve(RObject b) {
  /* wrapper for multiple type rhs */
  switch(b.sexp_type()) {
  case INTSXP:
  case REALSXP:
    if (b.hasAttribute("dim")) {
      // many dense rhs
      return solvem(as<NumericMatrix>(b));
    } else {
      // one dense rhs
      return solvev(as<NumericVector>(b));
    }
    break;
  case S4SXP:
    // sparse or dense many rhs
    if (b.inherits("dgeMatrix")) {
      // many dense rhs
      NumericVector bc(b.attr("x"));
      bc.attr("dim")=b.attr("Dim");
      return solvem(as<NumericMatrix>(bc));
    } else {
      // many sparse rhs
      return solves(as<S4>(b));
    }
    break;
  case VECSXP:
    // sparse or dense many rhs
    if (b.inherits("simple_triplet_matrix")) {
      // many sparse rhs
      return solvestm(as<List>(b));
    } else {
      sprintf(buf, "expected simple_triplet_matrix but got something else");
      stop(buf);
    }
    break;
  default:
    sprintf(buf, "unauthorized SEXP type of rhs (%d)", b.sexp_type());
    stop(buf);
    return R_NilValue;
  }
  /*return List::create(Named("is_null") = wrap(b.isNULL()),
                      Named("type") = wrap(b.sexp_type()),
                      Named("has_class") = wrap(b.isObject()),
                      Named("is_S4") = wrap(b.isS4()),
                      Named("attr_names") = b.attributeNames());
  */
  return R_NilValue;
}
NumericVector Rmumps::solvev(NumericVector b) {
  if (copy) {
    rhs=clone(b); // will hold the solution without altering b
  } else {
    rhs=b;
  }
  param.rhs=&*rhs.begin();
  param.nrhs=1;
  param.lrhs=0;
  param.ICNTL(20)=0; // rhs is dense
  do_job(6);
  return(rhs);
}
NumericMatrix Rmumps::solvem(NumericMatrix b) {
  if (copy) {
    mrhs=clone(b);
  } else {
    mrhs=b;
  }
  param.lrhs=b.nrow();
  param.ICNTL(20)=0; // rhs is dense
  param.rhs=&*mrhs.begin();
  param.nrhs=b.ncol();
#if 0 // #ifdef _OPENMP
  // unfortunately, MUMPS is not thread safe, there is a writing in working arrays during solving phase
  int ncol=b.ncol();
  int nthr=std::max(1, std::min(ncol , ncore)); // thread number to launch
  int ncol_th=b.ncol()/nthr; // col number per thread
  if (nthr == 1) {
    param.nrhs=ncol;
    do_job(6);
  } else {
    #pragma omp parallel for num_threads(nthr)
    for (int ith=0; ith < nthr; ith++) {
      Rmumps pth=*this;
      pth.ref++;
      /* Rcout << "jobs=" << std::endl;
      std::copy(jobs.begin(), jobs.end(), std::ostream_iterator<int>(Rcout, ", "));
      Rcout << std::endl;
      Rcout << "pth.jobs=" << std::endl;
      std::copy(pth.jobs.begin(), pth.jobs.end(), std::ostream_iterator<int>(Rcout, ", "));
      Rcout << std::endl;
      */
      pth.param.rhs=&*mrhs.begin()+ith*ncol_th*pth.param.lrhs;
      pth.param.nrhs=(ith == nthr-1 ? ncol-ith*ncol_th : ncol_th);
      /*
      Rcout << "param.rhs=" << param.rhs << ", nrhs=" << pth.param.nrhs << std::endl;
      Rcout << "pth.param.rhs=" << pth.param.rhs << std::endl;
      */
      pth.do_job(6);
    }
  }
#else
  do_job(6);
#endif
  return(mrhs);
}
void Rmumps::solveptr(double* b, int lrhs, int nrhs) {
  // in place solve
  param.rhs=b;
  param.nrhs=nrhs;
  param.lrhs=lrhs;
  param.ICNTL(20)=0; // rhs is dense
  do_job(6);
}
NumericMatrix Rmumps::inv() {
  MUMPS_INT n=param.n;
  MUMPS_INT nrhs=param.n;
  MUMPS_INT nz_rhs=param.n;
  rhs_sparse=NumericVector(n, 1.);
  irhs_ptr.resize(nrhs+1);
  irhs_sparse.resize(nz_rhs);
  for (int i=0; i < nz_rhs; i++) { // move to 1-based indices
    irhs_sparse[i]=i+1;
  }
  for (int i=0; i <= nrhs; i++) { // move to 1-based indices
    irhs_ptr[i]=i+1;
  }
  mrhs=NumericMatrix(n, nrhs);
  param.nz_rhs=nz_rhs;
  param.nrhs=nrhs;
  param.lrhs=n;
  param.irhs_ptr=&*irhs_ptr.begin();
  param.irhs_sparse=&*irhs_sparse.begin();
  param.rhs_sparse=&*rhs_sparse.begin();
  param.rhs=&*mrhs.begin(); // will hold the solution
  param.ICNTL(20)=1; // decide automaticaly about this sparsity exploit
  do_job(6);
  return(mrhs);
}
NumericMatrix Rmumps::solves(S4 mat) {
  // solve sparse rhs (may be multiple)
  // if dim(mat) == (0,0) then an inverse matrix is requested
  IntegerVector di(mat.slot("Dim"));
  if (!mat.inherits("dgCMatrix")) {
    Environment meth("package:methods");
    Function as_=meth["as"];
    mat=as_(mat, "dgCMatrix");
    //stop("matrix must be of dgTMatrix class (cf. pkg Matrix)");
  }
  //if (!mat.inherits("dgCMatrix")) {
  //  CharacterVector cstr(mat.slot("class"));
  //  sprintf(buf, "sparse rhs matrix must be of dgCMatrix class. Instead '%s' class is received (cf. pkg Matrix)", as<std::string>(cstr[0]).c_str());
  //  stop(buf);
  //}
  if (di[0] == 0 && di[1] == 0) {
     return(inv());
  }
  if (di[0] != param.n) stop("sparse rhs matrix must have the same number of rows than system matrix A");
  /* prepare a triplicate: sparse, pointer rhs_sparse */
  rhs_sparse=mat.slot("x");
  MUMPS_INT n=di[0];;
  MUMPS_INT nrhs=di[1];
  MUMPS_INT nz_rhs=rhs_sparse.size();
  IntegerVector mi(mat.slot("i"));
  IntegerVector mp(mat.slot("p"));
  irhs_ptr.resize(nrhs+1);
  irhs_sparse.resize(nz_rhs);
  for (int i=0; i < nz_rhs; i++) { // move to 1-based indices
    irhs_sparse[i]=mi[i]+1;
  }
  for (int i=0; i <= nrhs; i++) { // move to 1-based indices
    irhs_ptr[i]=mp[i]+1;
  }

  mrhs=NumericMatrix(n, nrhs);
  param.nz_rhs=nz_rhs;
  param.nrhs=nrhs;
  param.lrhs=n;
  param.irhs_ptr=&*irhs_ptr.begin();
  param.irhs_sparse=&*irhs_sparse.begin();
  param.rhs_sparse=&*rhs_sparse.begin();
  param.rhs=&*mrhs.begin(); // will hold the solution
  param.ICNTL(20)=1; // decide automaticaly about this sparsity exploit
  do_job(6);
  return(mrhs);
}
NumericMatrix Rmumps::solvestm(List mat) {
  // solve sparse rhs (may be multiple)
  // mat is expected to be of type slam::simple_triplet_matrix
  // if dim(mat) == (0,0) then an inverse matrix is requested
  int nrow=mat["nrow"];
  int ncol=mat["ncol"];
  if (!mat.inherits("simple_triplet_matrix")) {
    sprintf(buf, "solvestm() expects an rhs matrix of simple_triplet_matrix class");
    stop(buf);
  }
  if (nrow == 0 && ncol == 0) {
     return(inv());
  }
  if (nrow != param.n) stop("sparse rhs matrix must have the same number of rows than system matrix A");
  /* prepare a triplicate: sparse, pointer rhs_sparse */
  IntegerVector mi(as<IntegerVector>(mat["i"]));
  IntegerVector mj(as<IntegerVector>(mat["j"]));
  NumericVector mv(as<NumericVector>(mat["v"]));
//print(mi);
//print(mj);
//print(mv);

  MUMPS_INT n=nrow;
  MUMPS_INT nrhs=ncol;
  MUMPS_INT nz_rhs=mi.size();
  irhs_ptr.resize(nrhs+1);
  irhs_sparse.resize(nz_rhs);

  // bring to ccs format
  // first sort, as stm does not impose this
  IntegerVector iv1=mi+(mj-1)*nrow;
  IntegerVector iv1s=clone(iv1).sort();
  IntegerVector o=match(iv1s, iv1);
  mi=mi[o-1];
  mj=mj[o-1];
  rhs_sparse=mv[o-1];
//print(mi);
//print(mj);
//print(mv);
  for (int i=0; i < nz_rhs; i++) {
    irhs_sparse[i]=mi[i];
  }
  irhs_ptr[0]=1;
  int ip=0;
  for (int i=1; i <= nrhs; i++) {
    // count entries in this column
    int count=0;
    for (; ip < nz_rhs && mj[ip] == i; ip++, count++) {
//Rcout << "ip=" << ip << "; c=" << count << std::endl;
    }
//Rcout << "after ip=" << ip << "; c=" << count << std::endl;
    irhs_ptr[i]=irhs_ptr[i-1]+count;
  }
//Rcout << "nz_rhs=" << nz_rhs << std::endl;
//print(IntegerVector(irhs_ptr.begin(), irhs_ptr.end()));
  mrhs=NumericMatrix(n, nrhs);
  param.nz_rhs=nz_rhs;
  param.nrhs=nrhs;
  param.lrhs=n;
  param.irhs_ptr=&*irhs_ptr.begin();
  param.irhs_sparse=&*irhs_sparse.begin();
  param.rhs_sparse=&*rhs_sparse.begin();
  param.rhs=&*mrhs.begin(); // will hold the solution
  param.ICNTL(20)=1; // decide automaticaly about this sparsity exploit
  do_job(6);
  return(mrhs);
}
void Rmumps::set_rhs(NumericVector b) {
  // one dense rhs
  param.ICNTL(20)=0; // rhs is dense;
  if (copy) {
    rhs=clone(b);
  } else {
    rhs=b;
  }
  param.rhs=&*rhs.begin();
  param.nrhs=1;
  param.lrhs=0;
}
void Rmumps::set_mrhs(NumericMatrix b) {
  // many dense rhs
  param.ICNTL(20)=0; // rhs is dense;
  if (copy) {
    mrhs=clone(b);
  } else {
    mrhs=b;
  }
  param.rhs=&*mrhs.begin();
  param.nrhs=mrhs.ncol();
  param.lrhs=mrhs.nrow();
}
void Rmumps::set_mat_data(NumericVector x) {
  // for the same matrix pattern, set new entry values
  if (copy) {
    anz=clone(x);
  } else {
    anz=x;
  }
  param.a=&*anz.begin();
  // tell that previous numeric factorization is no more valid
  jobs.erase(2);
}
NumericVector Rmumps::get_rhs() {
  return rhs;
}
NumericMatrix Rmumps::get_mrhs() {
  return mrhs;
}
void Rmumps::set_icntl(IntegerVector iv, IntegerVector ii) {
  // set control vector ICNTL at positions in ii (1-based)) to the values in iv
  // only 1 <= ii <= 33 are effectively used
  if (iv.size() != ii.size()) {
    sprintf(buf, "set_icntl: length(iv) and length(ii) must be the same (got %d and %d respectively)", (int) iv.size(), (int) ii.size());
    stop(buf);
  }
  for (auto i=0; i < ii.size(); i++) {
    if (ii[i] > 33 || ii[i] < 1)
      continue;
    param.ICNTL(ii[i])=iv[i];
  }
}
void Rmumps::set_cntl(NumericVector v, IntegerVector iv) {
  // set control vector CNTL at positions in iv (1-based)) to the values in v
  // only 1 <= iv <= 5 are effectively used
  if (v.size() != iv.size()) {
    sprintf(buf, "set_cntl: length(v) and length(iv) must be the same (got %d and %d respectively)", (int) v.size(), (int) iv.size());
    stop(buf);
  }
  for (auto i=0; i < iv.size(); i++) {
    if (iv[i] > 5 || iv[i] < 1)
      continue;
    param.cntl[iv[i]-1]=v[i];
  }
}
IntegerVector Rmumps::get_icntl() {
  // return control vector ICNTL
  IntegerVector icntl(33);
  for (auto i=0; i < icntl.size(); i++) {
    icntl[i]=param.icntl[i];
  }
  return icntl;
}
NumericVector Rmumps::get_cntl() {
  // return control vector ICNTL
  NumericVector cntl(5);
  for (auto i=0; i < cntl.size(); i++) {
    cntl[i]=param.cntl[i];
  }
  return cntl;
}
List Rmumps::get_infos() {
  // return a named list with vectors info, rinfo, infog and rinfog
  NumericVector rinfo(3), rinfog(13);
  IntegerVector info(27), infog(34);
  for (auto i=0; i < rinfo.size(); i++)
    rinfo[i]=param.rinfo[i];
  for (auto i=0; i < info.size(); i++)
    info[i]=param.info[i];
  for (auto i=0; i < rinfog.size(); i++)
    rinfog[i]=param.rinfog[i];
  for (auto i=0; i < infog.size(); i++)
    infog[i]=param.infog[i];
  return List::create(
    _("info")=info,
    _("rinfo")=rinfo,
    _("infog")=infog,
    _("rinfog")=rinfog
  );
}
void Rmumps::set_keep(IntegerVector iv, IntegerVector ii) {
  // set control vector KEEP at positions in ii (1-based)) to the values in iv
  // only 1 <= ii <= 500 are effectively used
  if (iv.size() != ii.size()) {
    sprintf(buf, "set_keep: length(iv) and length(ii) must be the same (got %d and %d respectively)", (int) iv.size(), (int) ii.size());
    stop(buf);
  }
  for (auto i=0; i < ii.size(); i++) {
    if (ii[i] > 500 || ii[i] < 1)
      continue;
    param.keep[ii[i]-1]=iv[i];
  }
}

IntegerVector Rmumps::dim() {
  return IntegerVector::create(param.n, param.n);
}
int Rmumps::nrow() {
  return param.n;
}
int Rmumps::ncol() {
  return param.n;
}
void Rmumps::print() {
  Rcout << "A " << param.n << "x" << param.n << " Rmumps matrix" << std::endl;
  Rcout << "Decomposition(s) done: ";
  if (jobs.count(1) == 1) {
    Rcout << "symbolic";
  } else {
    Rcout << "none";
  }
  if (jobs.count(2) == 1) {
    Rcout << ", numeric";
  }
  Rcout << std::endl;
}
List Rmumps::triplet() {
  // returns a list with fields i, j, v of class slam::simple_triplet_matrix
  // representing the stored matrix. Objects i, j, v are cloned so the
  // returned object can be manipulated without risk for the original matrix
  List tr=List::create(
    _["i"]=IntegerVector(irn.begin(), irn.end()),
    _["j"]=IntegerVector(jcn.begin(), jcn.end()),
    _["v"]=NumericVector(anz.begin(), anz.end()),
    _["nrow"]=IntegerVector(1, param.n),
    _["ncol"]=IntegerVector(1, param.n),
    _["dimnames"]=R_NilValue
  );
  tr.attr("class")="simple_triplet_matrix";
  return tr;
}
std::string Rmumps::mumps_version() { return MUMPS_VERSION; }
double Rmumps::det() {
  if (jobs.count(2) != 1 || param.ICNTL(33) != 1) {
    param.ICNTL(33)=1;
    do_job(4);
  }
  return param.rinfog[12-1]*exp2(param.infog[34-1]);
}

/* constructors */
Rmumps::Rmumps(RObject mat, int sym, bool copy_) {
  new_mat(mat, sym, copy_);
}
Rmumps::Rmumps(RObject mat, int sym) {
  new_mat(mat, sym, true);
}
Rmumps::Rmumps(RObject mat) {
  new_mat(mat, 0, true); // default for sym=0 (unsymetric), copy_true
}
Rmumps::Rmumps(IntegerVector i0, IntegerVector j0, NumericVector x, int n) {
  new_ijv(i0, j0, x, n, 0, true);
}
Rmumps::Rmumps(IntegerVector i0, IntegerVector j0, NumericVector x, int n, int sym) {
  new_ijv(i0, j0, x, n, sym, true);
}
Rmumps::Rmumps(IntegerVector i0, IntegerVector j0, NumericVector x, int n, int sym, bool copy_) {
  new_ijv(i0, j0, x, n, sym, copy_);
}
/* end of constructors */

/* helpers */
void Rmumps::new_mat(RObject mat_, int sym, bool copy_) {
  MUMPS_INT n=-1;
  MUMPS_INT nz=-1;
  switch(mat_.sexp_type()) {
  case S4SXP: {
    S4 mat=as<S4>(mat_);
    if (!mat.inherits("dgTMatrix")) {
       Environment meth("package:methods");
       Function as_=meth["as"];
       mat=as_(mat, "dgTMatrix");
       //stop("matrix must be of dgTMatrix class (cf. pkg Matrix)");
    }
    IntegerVector di(mat.slot("Dim"));
    if (di[0] != di[1]) stop("matrix must be square");
    /* prepare a triplicate: irn, jcn, a */
    NumericVector x(mat.slot("x"));
    IntegerVector mi(mat.slot("i"));
    IntegerVector mj(mat.slot("j"));
    n=di[0];
    nz=x.size();
    irn.resize(nz);
    jcn.resize(nz);
    copy=copy_;
    if (copy) {
      anz=clone(x);
    } else {
      anz=x;
    }
    for (int i=0; i < nz; i++) { // move to 1-based indices
      irn[i]=mi[i]+1;
      jcn[i]=mj[i]+1;
    }
    break;
  }
  case VECSXP: {
    List mat=as<List>(mat_);
    if (!mat_.inherits("simple_triplet_matrix")) {
       stop("constructor input of type VECSXP must be of simple_triplet_matrix class (cf. pkg slam)");
    }
    n=mat["nrow"];
    if (n != as<int>(mat["ncol"])) stop("matrix must be square");
    /* prepare a triplicate: irn, jcn, a */
    NumericVector x(as<NumericVector>(mat["v"]));
    IntegerVector mi(as<IntegerVector>(mat["i"]));
    IntegerVector mj(as<IntegerVector>(mat["j"]));
    nz=x.size();
    irn.resize(nz);
    jcn.resize(nz);
    copy=copy_;
    if (copy) {
      anz=clone(x);
    } else {
      anz=x;
    }
    for (int i=0; i < nz; i++) { // already in 1-based indices
      irn[i]=mi[i];
      jcn[i]=mj[i];
    }
    break;
  }
  default:
    sprintf(buf, "constructor form a single object is expecting Matrix::dgTMatrix (i.e. S4SXP) or slam::simple_triplet_matrix (i.e. VECSXP) class as input. Got '%d' SEXP instead", mat_.sexp_type());
    stop(buf);
  }
  //Rf_PrintValue(wrap(irn));
  //Rf_PrintValue(wrap(jcn));
  //Rf_PrintValue(wrap(anz));
  tri_init(&*irn.begin(), &*jcn.begin(), &*anz.begin(), sym);
  param.n=n;
  param.nz=nz;
}
void Rmumps::new_ijv(IntegerVector i0, IntegerVector j0, NumericVector x, int n_, int sym, bool copy_) {
  MUMPS_INT nz=x.size();
  MUMPS_INT n=n_;
  irn.resize(nz);
  jcn.resize(nz);
  for (int i=0; i < nz; i++) { // move to 1-based indices
    irn[i]=i0[i]+1;
    jcn[i]=j0[i]+1;
  }
  copy=copy_;
  if (copy) {
    anz=clone(x);
  } else {
    anz=x;
  }
  //Rf_PrintValue(wrap(irn));
  //Rf_PrintValue(wrap(jcn));
  //Rf_PrintValue(wrap(anz));
  tri_init(&*irn.begin(), &*jcn.begin(), &*anz.begin(), sym);
  param.n=n;
  param.nz=nz;
}

void Rmumps::tri_init(MUMPS_INT *irn, MUMPS_INT *jcn, double *a, MUMPS_INT sym) {
  this->sym=sym;
  this->ncore=1;
  this->ref=0;
  /* Initialize a MUMPS instance. Use MPI_COMM_WORLD */
  param.job=JOB_INIT;
  param.keep[39]=0; // otherwise valgrind complaints
  param.par=1;
  param.sym=sym;
  param.comm_fortran=USE_COMM_WORLD;
  do_job(JOB_INIT);
  param.irn=irn;
  param.jcn=jcn;
  param.a=a;
  /* No outputs */
  param.ICNTL(1)=-1; param.ICNTL(2)=-1; param.ICNTL(3)=-1; param.ICNTL(4)=0;
  param.ICNTL(5)=0; param.ICNTL(18)=0; // matrix is centralized and assembled
  param.ICNTL(6)=0; // no column permutation
  param.ICNTL(7)=7; // automatic choice for symmetric permutation
  param.ICNTL(8)=77; // automatic choice for scaling
  param.ICNTL(9)=1; // solve a x = b (not a^t x = b)
  param.ICNTL(10)=0; // no iterative refinement
  param.ICNTL(11)=0; // no error analysis
  param.ICNTL(12)=0; // not used (as sym != 2)
  param.ICNTL(13)=0; // parallel factorization of the root node
  param.ICNTL(14)=50; // 50% working space increase during numeric phase
  // 15-17 are not used
  param.ICNTL(19)=0; // compleate factorization (no schur complement)
  param.ICNTL(20)=0; // dense rhs
  param.ICNTL(21)=0; // the solution is centralized (i.e. not distributed)
  param.ICNTL(22)=0; // in-core factorization
  param.ICNTL(23)=0; // decide him self for working space
  param.ICNTL(24)=0; // null pivot row will give an error
  param.ICNTL(25)=0; // normal solution step 
  param.ICNTL(26)=0; // not used as (19) == 0
  param.ICNTL(27)=-8; // automatic bloc size for mrhs
  param.ICNTL(28)=0; // automatic choice seq or parallel ordering
  param.ICNTL(29)=0; // automatic choice for tool of parallel ordering
  param.ICNTL(30)=0; // no selected entries of a^-1 are calculated
  param.ICNTL(31)=0; // matrix factors are kept (not discarded)
  param.ICNTL(32)=0; // standart factorization (without forward elim of rhs)
  param.ICNTL(33)=1; // compute the detreminant
}
//RCPP_EXPOSED_CLASS(Rmumps)
RCPP_MODULE(mod_Rmumps){
  using namespace Rcpp ;
  class_<Rmumps>("Rmumps")
  // expose the default constructor
  .constructor<SEXP>()
  .constructor<SEXP, int>()
  .constructor<SEXP, int, bool>()
  .constructor<IntegerVector, IntegerVector, NumericVector, int>()
  .constructor<IntegerVector, IntegerVector, NumericVector, int, int>()
  .constructor<IntegerVector, IntegerVector, NumericVector, int, int, bool>()
  //.finalizer(&cleanme) // crashes by double freeing of the same pointer
  
  .property("rhs", &Rmumps::get_rhs, &Rmumps::set_rhs)
  .property("mrhs", &Rmumps::get_mrhs, &Rmumps::set_mrhs)
  .field("copy", &Rmumps::copy, "copy or not input parameters")
  .field_readonly("sym", &Rmumps::sym)
  .field("ncore", &Rmumps::ncore, "how many cores to use within OpenMP regions")
  
  .method("symbolic", &Rmumps::symbolic , "Analyze sparsity pattern")
  .method("numeric", &Rmumps::numeric, "Factorize sparse matrix")
  .method("solve", &Rmumps::solve, "Solve sparse system with one or many, sparse or dense rhs")
  .method("solvet", &Rmumps::solvet, "Solve transpose of sparse system with one or many, sparse or dense rhs")
  .method("inv", &Rmumps::inv, "Calculate the inverse of a sparse matrix")
  .method("set_mat_data", &Rmumps::set_mat_data, "Update matrix entries keeping the non zero pattern untouched")
  .method("set_icntl", &Rmumps::set_icntl, "Set ICNTL parameter vector")
  .method("get_icntl", &Rmumps::get_icntl, "Get ICNTL parameter vector")
  .method("set_cntl", &Rmumps::set_cntl, "Set CNTL parameter vector")
  .method("get_cntl", &Rmumps::get_cntl, "Get CNTL parameter vector")
  .method("get_infos", &Rmumps::get_infos, "Get a named list of information vectors")
  .method("set_keep", &Rmumps::set_icntl, "Set KEEP parameter vector")
  .method("dim", &Rmumps::dim, "Return a vector with matrix dimensions")
  .method("nrow", &Rmumps::nrow, "Return an integer with matrix row number")
  .method("ncol", &Rmumps::ncol, "Return an integer with matrix column number")
  .method("print", &Rmumps::print, "Print the size of matrix and decompositions done")
  .method("show", &Rmumps::print, "Print the size of matrix and decompositions done")
  .method("triplet", &Rmumps::triplet, "Return an object of simple_triplet_matrix class with i, j, v fields representing the matrix")
  .method("det", &Rmumps::det, "Return determinant of the matrix")
  .method("mumps_version", &Rmumps::mumps_version, "Return determinant of the matrix")
  ;
}
