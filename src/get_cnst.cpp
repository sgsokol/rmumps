#include <string>
#include "../inst/include/rmumps.h"

#define cnst_pair(c) {#c, c}

// get value of a c++ macro constant
//
// Return an integer number defined as c++ macro if it is present in internal dictionary. The input parameter is the name of the constant to be returned. If the name is not found in the dictionary, an error is thrown. Constants currently present in the dictionary are following: RMUMPS_PERM_AMD, RMUMPS_PERM_AMF, RMUMPS_PERM_SCOTCH, RMUMPS_PERM_PORD, RMUMPS_PERM_METIS, RMUMPS_PERM_QAMD and RMUMPS_PERM_AUTO.
// 
// @param s character scalar, name of the constant
// @return integer scalar, value of the sonstant
// @aliases RMUMPS_PERM_AMD RMUMPS_PERM_AMF RMUMPS_PERM_SCOTCH RMUMPS_PERM_PORD RMUMPS_PERM_METIS RMUMPS_PERM_QAMD RMUMPS_PERM_AUTO cnsts
// @examples
// get_cnst("RMUMPS_PERM_AMD")
// [[Rcpp::export(".get_cnst")]]
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
    stop(".get_cnst: constant '%s' is not in dictionary", s);
  return(NA_INTEGER);
}
