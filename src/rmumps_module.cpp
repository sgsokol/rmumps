// Sequential port of MUMPS to R.
// Only a subset of MUMPS is exposed to R.
// Author of port: Serguei Sokol, INRA, Toulouse, FRANCE
// Copyright 2015, INRA
// v 0.1 2015-09-02

#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::interfaces(r, cpp)]]

#include <dmumps_c.h>
//#include <libseq/mpi.h> // useless in sequential mode
#define JOB_INIT -1
#define JOB_END -2
#define USE_COMM_WORLD -987654
#define ICNTL(I) icntl[(I)-1] /* macro s.t. indices match documentation */

class Rmumps {
public:
  Rmumps(S4 mat, bool copy_);
  Rmumps(IntegerVector i, IntegerVector j, NumericVector x, int n, bool copy_);
  ~Rmumps() { clean(); };
  void clean() {
    //Rprintf("clean() is called...\n");
    param.job=JOB_END;
    dmumps_c(&param);
    //Rprintf("clean(): done\n");
  };
  bool get_copy() {
    return copy;
  }
  void set_copy(bool copy_) {
    copy=copy_;
  }
  void do_job(int job) { // later pass it to private scope
    /* do preliminary job if needed */
    switch (job) {
    case 2:
    case 5:
      if (jobs.count(1)==0) {
        // symbolic analysis was not yet done
        do_job(1);
      }
      break;
    case 3:
      if (jobs.count(2)==0) {
        // numeric factorization was not yet done
        do_job(2);
      }
      break;
    }
    /* donwgrade job request if preliminaries are already done */
    switch (job) {
    case 4:
      if (jobs.count(1)==1) {
         job=2;
      }
      break;
    case 5:
    case 6:
      if (jobs.count(2)==1) {
         job=3;
      }
      break;
    }
    param.job=job;
    dmumps_c(&param);
    if (param.info[0] != 0) {
      //clean();
      stop("rmumps: info[1]=%d, info[2]=%d", param.info[0], param.info[1]);
    }
    /* combined jobs are splited for the record */
    switch (job) {
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
  };
  /*void do_job(int job) {
    param.job=job;
    dmumps_c(&param);
    if (param.info[0] != 0) {
      //clean();
      stop("rmumps: info[1]=%d, info[2]=%d", param.info[0], param.info[1]);
    }
  }*/
  void symbolic() {
    do_job(1);
  };
  void numeric() {
    do_job(4);
  };
  SEXP solve(RObject b) {
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
  }
  NumericVector solvev(NumericVector b) {
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
  };
  NumericMatrix solvem(NumericMatrix b) {
    if (copy) {
      mrhs=clone(b);
    } else {
      mrhs=b;
    }
    param.rhs=&*mrhs.begin();
    param.nrhs=b.ncol();
    param.lrhs=b.nrow();
    param.ICNTL(20)=0; // rhs is dense
    do_job(6);
    return(mrhs);
  };
  NumericMatrix inv() {
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
  NumericMatrix solves(S4 mat) {
    // solve sparse rhs (may be multiple)
    // if dim(mat)==(0,0) then an inverse matrix is requested
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
  };
  void set_rhs(NumericVector b) {
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
  void set_mrhs(NumericMatrix b) {
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
  void set_mat_data(NumericVector x) {
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
  NumericVector get_rhs() {
    return rhs;
  }
  NumericMatrix get_mrhs() {
    return mrhs;
  }
  std::string mumps_version() { return MUMPS_VERSION; };
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
private:
  DMUMPS_STRUC_C param;
  void tri_init(MUMPS_INT *irn, MUMPS_INT *jcn, double *a);
  char buf[512];
};

/* constructors */
Rmumps::Rmumps(S4 mat, bool copy_=true) {
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
  MUMPS_INT n=di[0];
  MUMPS_INT nz=x.size();
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
  tri_init(&*irn.begin(), &*jcn.begin(), &*anz.begin());
  param.n=n;
  param.nz=nz;
}
Rmumps::Rmumps(IntegerVector i0, IntegerVector j0, NumericVector x, int n, bool copy_=true) {
  MUMPS_INT nz=x.size();
  std::vector<MUMPS_INT> irn(nz);
  std::vector<MUMPS_INT> jcn(nz);
  for (int i=0; i < nz; i++) { // move to 1-based indices
    irn[i]=i0[i]+1;
    jcn[i]=j0[i]+1;
  }
  tri_init(&*irn.begin(), &*jcn.begin(), &*x.begin());
  param.n=n;
  param.nz=nz;
  copy=copy_;
}

/* other methods */
void Rmumps::tri_init(MUMPS_INT *irn, MUMPS_INT *jcn, double *a) {
  /* Initialize a MUMPS instance. Use MPI_COMM_WORLD */
  param.job=JOB_INIT;
  param.par=1;
  param.sym=0;
  param.comm_fortran=USE_COMM_WORLD;
  do_job(JOB_INIT);
  param.irn=irn;
  param.jcn=jcn;
  param.a=a;
  /* No outputs */
  param.ICNTL(1)=-1; param.ICNTL(2)=-1; param.ICNTL(3)=-1; param.ICNTL(4)=0;
}
/*
void cleanme(Rmumps* obj) {
  obj->clean();
};
*/

RCPP_MODULE(mod_Rmumps){
  using namespace Rcpp ;
  class_<Rmumps>("Rmumps")
  // expose the default constructor
  .constructor<S4>()
  .constructor<S4, bool>()
  .constructor<IntegerVector, IntegerVector, NumericVector, int>()
  .constructor<IntegerVector, IntegerVector, NumericVector, int, bool>()
  //.finalizer(&cleanme) // makes crash by double freeing of the same pointer
  
  .property("rhs", &Rmumps::get_rhs, &Rmumps::set_rhs)
  .property("mrhs", &Rmumps::get_mrhs, &Rmumps::set_mrhs)
  .property("copy", &Rmumps::get_copy, &Rmumps::set_copy)
  
  .method("symbolic", &Rmumps::symbolic , "Analyze sparsity pattern")
  .method("numeric", &Rmumps::numeric, "Factorize sparse matrix")
  .method("solve", &Rmumps::solve, "Solve sparse system with one or many, sparse or dense rhs")
  .method("inv", &Rmumps::inv, "Calculate the inverse of a sparse matrix")
  .method("set_mat_data", &Rmumps::set_mat_data, "Epdate matrix entries keeping the non zero pattern untouched")
  //.method("solvev", &Rmumps::solvev, "Solve sparse system with one rhs")
  //.method("solvem", &Rmumps::solvem, "Solve sparse system with many rhs")
  //.method("solves", &Rmumps::solves, "Solve sparse system with many sparse rhs")
  //.method("do_job", &Rmumps::do_job, "Custom call to dmumps_c()") // comment it out it later
  ;
}
