#' Exported Constants
#'
#' @name RMUMPS_PERM
#' @aliases RMUMPS_PERM RMUMPS_PERM_AMD RMUMPS_PERM_AMF RMUMPS_PERM_SCOTCH RMUMPS_PERM_PORD RMUMPS_PERM_METIS RMUMPS_PERM_QAMD RMUMPS_PERM_AUTO
#' @description
#' Integer constants defining permutation types and exported from rmumps are following:
#' \itemize{
#' \item{\code{RMUMPS_PERM_AMD}}
#' \item{\code{RMUMPS_PERM_AMF}}
#' \item{\code{RMUMPS_PERM_SCOTCH}}
#' \item{\code{RMUMPS_PERM_PORD}}
#' \item{\code{RMUMPS_PERM_METIS}}
#' \item{\code{RMUMPS_PERM_QAMD}}
#' \item{\code{RMUMPS_PERM_AUTO}}
#' }
#' They are all regrouped in a named vector \code{RMUMPS_PERM} where names are items above and values are corresponding constants.
#' @examples
#' am=rmumps::Rmumps$new(slam::as.simple_triplet_matrix(diag(1:3)))
#' am$set_permutation(RMUMPS_PERM_SCOTCH)
#' am$solve(1:3)
NULL

# set useful constants
loadModule("mod_Rmumps", TRUE)
.onLoad <- function(libname, pkgname){
  env=parent.env(environment())
  for (cc in c(
   "RMUMPS_PERM_AMD",
   "RMUMPS_PERM_AMF",
   "RMUMPS_PERM_SCOTCH",
   "RMUMPS_PERM_PORD",
   "RMUMPS_PERM_METIS",
   "RMUMPS_PERM_QAMD",
   "RMUMPS_PERM_AUTO"
)) {
     var=.get_cnst(cc)
     env[[cc]]=var
     env[["RMUMPS_PERM"]]=c(env[["RMUMPS_PERM"]], structure(var, names=cc))
  }
}
