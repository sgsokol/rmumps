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
  Rf_warning("\nTiming Information -------------------------------------------------");
  Rf_warning("\n Multilevel: \t\t %7.3"PRREAL"", gk_getcputimer(ctrl->TotalTmr));
  Rf_warning("\n     Coarsening: \t\t %7.3"PRREAL"", gk_getcputimer(ctrl->CoarsenTmr));
  Rf_warning("\n            Matching: \t\t\t %7.3"PRREAL"", gk_getcputimer(ctrl->MatchTmr));
  Rf_warning("\n            Contract: \t\t\t %7.3"PRREAL"", gk_getcputimer(ctrl->ContractTmr));
  Rf_warning("\n     Initial Partition: \t %7.3"PRREAL"", gk_getcputimer(ctrl->InitPartTmr));
  Rf_warning("\n     Uncoarsening: \t\t %7.3"PRREAL"", gk_getcputimer(ctrl->UncoarsenTmr));
  Rf_warning("\n          Refinement: \t\t\t %7.3"PRREAL"", gk_getcputimer(ctrl->RefTmr));
  Rf_warning("\n          Projection: \t\t\t %7.3"PRREAL"", gk_getcputimer(ctrl->ProjectTmr));
  Rf_warning("\n     Splitting: \t\t %7.3"PRREAL"", gk_getcputimer(ctrl->SplitTmr));
/*
  Rf_warning("\n       Aux1Tmr: \t\t %7.3"PRREAL"", gk_getcputimer(ctrl->Aux1Tmr));
  Rf_warning("\n       Aux2Tmr: \t\t %7.3"PRREAL"", gk_getcputimer(ctrl->Aux2Tmr));
  Rf_warning("\n       Aux3Tmr: \t\t %7.3"PRREAL"", gk_getcputimer(ctrl->Aux3Tmr));
*/
  Rf_warning("\n********************************************************************\n");
}



