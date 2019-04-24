#include <string>
#include "../inst/include/rmumps.h"

#define cnst_pair(c) {#c, c}

// [[Rcpp::export]]
int get_cnst(std::string s) {
  static std::map<std::string, int> dict={
    cnst_pair(RMUMPS_PERM_AMD),
    cnst_pair(RMUMPS_PERM_AMF),
    cnst_pair(RMUMPS_PERM_SCOTCH),
    cnst_pair(RMUMPS_PERM_PORD),
    cnst_pair(RMUMPS_PERM_METIS),
    cnst_pair(RMUMPS_PERM_QAMD),
    cnst_pair(RMUMPS_PERM_AUTO),
  };
  if (dict.count(s))
    return(dict[s]);
  else
    stop("get_cnst: constant '%s' is not in dictionary", s);
  return(NA_INTEGER);
}
