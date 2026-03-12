/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * timing.c
 *
 * This file contains routines that deal with timing Metis
 *
 * Started 7/24/97
 * George
 *
 * $Id: timing.c 13936 2013-03-30 03:59:09Z karypis $
 *
 */

#include "metislib.h"


/*************************************************************************
* This function clears the timers
**************************************************************************/
void InitTimers(ctrl_t *ctrl)
{
  gk_clearcputimer(ctrl->TotalTmr);
  gk_clearcputimer(ctrl->InitPartTmr);
  gk_clearcputimer(ctrl->MatchTmr);
  gk_clearcputimer(ctrl->ContractTmr);
  gk_clearcputimer(ctrl->CoarsenTmr);
  gk_clearcputimer(ctrl->UncoarsenTmr);
  gk_clearcputimer(ctrl->RefTmr);
  gk_clearcputimer(ctrl->ProjectTmr);
  gk_clearcputimer(ctrl->SplitTmr);
  gk_clearcputimer(ctrl->Aux1Tmr);
  gk_clearcputimer(ctrl->Aux2Tmr);
  gk_clearcputimer(ctrl->Aux3Tmr);
}



/*************************************************************************
* This function prints the various timers
**************************************************************************/
void PrintTimers(ctrl_t *ctrl)
{
  rcpp_warning("\nTiming Information -------------------------------------------------");
  rcpp_warning("\n Multilevel: \t\t %7.3"PRREAL"", gk_getcputimer(ctrl->TotalTmr));
  rcpp_warning("\n     Coarsening: \t\t %7.3"PRREAL"", gk_getcputimer(ctrl->CoarsenTmr));
  rcpp_warning("\n            Matching: \t\t\t %7.3"PRREAL"", gk_getcputimer(ctrl->MatchTmr));
  rcpp_warning("\n            Contract: \t\t\t %7.3"PRREAL"", gk_getcputimer(ctrl->ContractTmr));
  rcpp_warning("\n     Initial Partition: \t %7.3"PRREAL"", gk_getcputimer(ctrl->InitPartTmr));
  rcpp_warning("\n     Uncoarsening: \t\t %7.3"PRREAL"", gk_getcputimer(ctrl->UncoarsenTmr));
  rcpp_warning("\n          Refinement: \t\t\t %7.3"PRREAL"", gk_getcputimer(ctrl->RefTmr));
  rcpp_warning("\n          Projection: \t\t\t %7.3"PRREAL"", gk_getcputimer(ctrl->ProjectTmr));
  rcpp_warning("\n     Splitting: \t\t %7.3"PRREAL"", gk_getcputimer(ctrl->SplitTmr));
/*
  rcpp_warning("\n       Aux1Tmr: \t\t %7.3"PRREAL"", gk_getcputimer(ctrl->Aux1Tmr));
  rcpp_warning("\n       Aux2Tmr: \t\t %7.3"PRREAL"", gk_getcputimer(ctrl->Aux2Tmr));
  rcpp_warning("\n       Aux3Tmr: \t\t %7.3"PRREAL"", gk_getcputimer(ctrl->Aux3Tmr));
*/
  rcpp_warning("\n********************************************************************\n");
}



