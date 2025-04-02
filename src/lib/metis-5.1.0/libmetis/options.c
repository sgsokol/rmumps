/**
  \file
  \brief This file contains various routines for dealing with options and ctrl_t.

  \date   Started 5/12/2011
  \author George  
  \author Copyright 1997-2011, Regents of the University of Minnesota 
  \version\verbatim $Id: options.c 13901 2013-03-24 16:17:03Z karypis $ \endverbatim
  */

#include "metislib.h"


/*************************************************************************/
/*! This function creates and sets the run parameters (ctrl_t) */
/*************************************************************************/
ctrl_t *SetupCtrl(moptype_et optype, idx_t *options, idx_t ncon, idx_t nparts, 
            real_t *tpwgts, real_t *ubvec)
{
  idx_t i, j;
  ctrl_t *ctrl;

  ctrl = (ctrl_t *)gk_malloc(sizeof(ctrl_t), "SetupCtrl: ctrl");
  
  memset((void *)ctrl, 0, sizeof(ctrl_t));

  switch (optype) {
    case METIS_OP_PMETIS:
      ctrl->objtype = GETOPTION(options, METIS_OPTION_OBJTYPE, METIS_OBJTYPE_CUT);
      ctrl->rtype   = METIS_RTYPE_FM;
      ctrl->ncuts   = GETOPTION(options, METIS_OPTION_NCUTS,   1);
      ctrl->niter   = GETOPTION(options, METIS_OPTION_NITER,   10);

      if (ncon == 1) {
        ctrl->iptype    = GETOPTION(options, METIS_OPTION_IPTYPE,  METIS_IPTYPE_GROW);
        ctrl->ufactor   = GETOPTION(options, METIS_OPTION_UFACTOR, PMETIS_DEFAULT_UFACTOR);
        ctrl->CoarsenTo = 20;
      }
      else {
        ctrl->iptype    = GETOPTION(options, METIS_OPTION_IPTYPE,  METIS_IPTYPE_RANDOM);
        ctrl->ufactor   = GETOPTION(options, METIS_OPTION_UFACTOR, MCPMETIS_DEFAULT_UFACTOR);
        ctrl->CoarsenTo = 100;
      }

      break;


    case METIS_OP_KMETIS:
      ctrl->objtype = GETOPTION(options, METIS_OPTION_OBJTYPE, METIS_OBJTYPE_CUT);
      ctrl->iptype  = METIS_IPTYPE_METISRB;
      ctrl->rtype   = METIS_RTYPE_GREEDY;
      ctrl->ncuts   = GETOPTION(options, METIS_OPTION_NCUTS,   1);
      ctrl->niter   = GETOPTION(options, METIS_OPTION_NITER,   10);
      ctrl->ufactor = GETOPTION(options, METIS_OPTION_UFACTOR, KMETIS_DEFAULT_UFACTOR);
      ctrl->minconn = GETOPTION(options, METIS_OPTION_MINCONN, 0);
      ctrl->contig  = GETOPTION(options, METIS_OPTION_CONTIG,  0);
      break;


    case METIS_OP_OMETIS:
      ctrl->objtype  = GETOPTION(options, METIS_OPTION_OBJTYPE,  METIS_OBJTYPE_NODE);
      ctrl->rtype    = GETOPTION(options, METIS_OPTION_RTYPE,    METIS_RTYPE_SEP1SIDED);
      ctrl->iptype   = GETOPTION(options, METIS_OPTION_IPTYPE,   METIS_IPTYPE_EDGE);
      ctrl->nseps    = GETOPTION(options, METIS_OPTION_NSEPS,    1);
      ctrl->niter    = GETOPTION(options, METIS_OPTION_NITER,    10);
      ctrl->ufactor  = GETOPTION(options, METIS_OPTION_UFACTOR,  OMETIS_DEFAULT_UFACTOR);
      ctrl->compress = GETOPTION(options, METIS_OPTION_COMPRESS, 1);
      ctrl->ccorder  = GETOPTION(options, METIS_OPTION_CCORDER,  0);
      ctrl->pfactor  = 0.1*GETOPTION(options, METIS_OPTION_PFACTOR,  0);

      ctrl->CoarsenTo = 100;
      break;

    default:
      gk_errexit(SIGERR, "Unknown optype of %d\n", optype);
  }

  /* common options */
  ctrl->ctype   = GETOPTION(options, METIS_OPTION_CTYPE, METIS_CTYPE_SHEM);
  ctrl->no2hop  = GETOPTION(options, METIS_OPTION_NO2HOP, 0);
  ctrl->seed    = GETOPTION(options, METIS_OPTION_SEED, -1);
  ctrl->dbglvl  = GETOPTION(options, METIS_OPTION_DBGLVL, 0);
  ctrl->numflag = GETOPTION(options, METIS_OPTION_NUMBERING, 0);

  /* set non-option information */
  ctrl->optype  = optype;
  ctrl->ncon    = ncon;
  ctrl->nparts  = nparts;
  ctrl->maxvwgt = ismalloc(ncon, 0, "SetupCtrl: maxvwgt");

  /* setup the target partition weights */
  if (ctrl->optype != METIS_OP_OMETIS) {
    ctrl->tpwgts = rmalloc(nparts*ncon, "SetupCtrl: ctrl->tpwgts");
    if (tpwgts) {
      rcopy(nparts*ncon, tpwgts, ctrl->tpwgts);
    }
    else {
      for (i=0; i<nparts; i++) {
        for (j=0; j<ncon; j++)
          ctrl->tpwgts[i*ncon+j] = 1.0/nparts;
      }
    }
  }
  else {  /* METIS_OP_OMETIS */
    /* this is required to allow the pijbm to be defined properly for
       the edge-based refinement during initial partitioning */
    ctrl->tpwgts = rsmalloc(2, .5,  "SetupCtrl: ctrl->tpwgts");
  }


  /* setup the ubfactors */
  ctrl->ubfactors = rsmalloc(ctrl->ncon, I2RUBFACTOR(ctrl->ufactor), "SetupCtrl: ubfactors");
  if (ubvec)
    rcopy(ctrl->ncon, ubvec, ctrl->ubfactors);
  for (i=0; i<ctrl->ncon; i++)
    ctrl->ubfactors[i] += 0.0000499;

  /* Allocate memory for balance multipliers. 
     Note that for PMETIS/OMETIS routines the memory allocated is more 
     than required as balance multipliers for 2 parts is sufficient. */
  ctrl->pijbm = rmalloc(nparts*ncon, "SetupCtrl: ctrl->pijbm");

  InitRandom(ctrl->seed);

  IFSET(ctrl->dbglvl, METIS_DBG_INFO, PrintCtrl(ctrl));

  if (!CheckParams(ctrl)) {
    FreeCtrl(&ctrl);
    return NULL;
  }
  else {
    return ctrl;
  }
}


/*************************************************************************/
/*! Computes the per-partition/constraint balance multipliers */
/*************************************************************************/
void SetupKWayBalMultipliers(ctrl_t *ctrl, graph_t *graph)
{
  idx_t i, j;

  for (i=0; i<ctrl->nparts; i++) {
    for (j=0; j<graph->ncon; j++)
      ctrl->pijbm[i*graph->ncon+j] = graph->invtvwgt[j]/ctrl->tpwgts[i*graph->ncon+j];
  }
}


/*************************************************************************/
/*! Computes the per-partition/constraint balance multipliers */
/*************************************************************************/
void Setup2WayBalMultipliers(ctrl_t *ctrl, graph_t *graph, real_t *tpwgts)
{
  idx_t i, j;

  for (i=0; i<2; i++) {
    for (j=0; j<graph->ncon; j++)
      ctrl->pijbm[i*graph->ncon+j] = graph->invtvwgt[j]/tpwgts[i*graph->ncon+j];
  }
}


/*************************************************************************/
/*! This function prints the various control fields */
/*************************************************************************/
void PrintCtrl(ctrl_t *ctrl)
{
  idx_t i, j, modnum;

  Rf_warning(" Runtime parameters:\n");

  Rf_warning("   Objective type: ");
  switch (ctrl->objtype) {
    case METIS_OBJTYPE_CUT:
      Rf_warning("METIS_OBJTYPE_CUT\n");
      break;
    case METIS_OBJTYPE_VOL:
      Rf_warning("METIS_OBJTYPE_VOL\n");
      break;
    case METIS_OBJTYPE_NODE:
      Rf_warning("METIS_OBJTYPE_NODE\n");
      break;
    default:
      Rf_warning("Unknown!\n");
  }

  Rf_warning("   Coarsening type: ");
  switch (ctrl->ctype) {
    case METIS_CTYPE_RM:
      Rf_warning("METIS_CTYPE_RM\n");
      break;
    case METIS_CTYPE_SHEM:
      Rf_warning("METIS_CTYPE_SHEM\n");
      break;
    default:
      Rf_warning("Unknown!\n");
  }

  Rf_warning("   Initial partitioning type: ");
  switch (ctrl->iptype) {
    case METIS_IPTYPE_GROW:
      Rf_warning("METIS_IPTYPE_GROW\n");
      break;
    case METIS_IPTYPE_RANDOM:
      Rf_warning("METIS_IPTYPE_RANDOM\n");
      break;
    case METIS_IPTYPE_EDGE:
      Rf_warning("METIS_IPTYPE_EDGE\n");
      break;
    case METIS_IPTYPE_NODE:
      Rf_warning("METIS_IPTYPE_NODE\n");
      break;
    case METIS_IPTYPE_METISRB:
      Rf_warning("METIS_IPTYPE_METISRB\n");
      break;
    default:
      Rf_warning("Unknown!\n");
  }

  Rf_warning("   Refinement type: ");
  switch (ctrl->rtype) {
    case METIS_RTYPE_FM:
      Rf_warning("METIS_RTYPE_FM\n");
      break;
    case METIS_RTYPE_GREEDY:
      Rf_warning("METIS_RTYPE_GREEDY\n");
      break;
    case METIS_RTYPE_SEP2SIDED:
      Rf_warning("METIS_RTYPE_SEP2SIDED\n");
      break;
    case METIS_RTYPE_SEP1SIDED:
      Rf_warning("METIS_RTYPE_SEP1SIDED\n");
      break;
    default:
      Rf_warning("Unknown!\n");
  }

  Rf_warning("   Perform a 2-hop matching: %s\n", (ctrl->no2hop ? "Yes" : "No"));

  Rf_warning("   Number of balancing constraints: %"PRIDX"\n", ctrl->ncon);
  Rf_warning("   Number of refinement iterations: %"PRIDX"\n", ctrl->niter);
  Rf_warning("   Random number seed: %"PRIDX"\n", ctrl->seed);

  if (ctrl->optype == METIS_OP_OMETIS) {
    Rf_warning("   Number of separators: %"PRIDX"\n", ctrl->nseps);
    Rf_warning("   Compress graph prior to ordering: %s\n", (ctrl->compress ? "Yes" : "No"));
    Rf_warning("   Detect & order connected components separately: %s\n", (ctrl->ccorder ? "Yes" : "No"));
    Rf_warning("   Prunning factor for high degree vertices: %"PRREAL"\n", ctrl->pfactor);
  }
  else {
    Rf_warning("   Number of partitions: %"PRIDX"\n", ctrl->nparts);
    Rf_warning("   Number of cuts: %"PRIDX"\n", ctrl->ncuts);
    Rf_warning("   User-supplied ufactor: %"PRIDX"\n", ctrl->ufactor);

    if (ctrl->optype == METIS_OP_KMETIS) {
      Rf_warning("   Minimize connectivity: %s\n", (ctrl->minconn ? "Yes" : "No"));
      Rf_warning("   Create contigous partitions: %s\n", (ctrl->contig ? "Yes" : "No"));
    }

    modnum = (ctrl->ncon==1 ? 5 : (ctrl->ncon==2 ? 3 : (ctrl->ncon==3 ? 2 : 1)));
    Rf_warning("   Target partition weights: ");
    for (i=0; i<ctrl->nparts; i++) {
      if (i%modnum == 0)
        Rf_warning("\n     ");
      Rf_warning("%4"PRIDX"=[", i);
      for (j=0; j<ctrl->ncon; j++) 
        Rf_warning("%s%.2e", (j==0 ? "" : " "), (double)ctrl->tpwgts[i*ctrl->ncon+j]);
      Rf_warning("]");
    }
    Rf_warning("\n");
  }

  Rf_warning("   Allowed maximum load imbalance: ");
  for (i=0; i<ctrl->ncon; i++) 
    Rf_warning("%.3"PRREAL" ", ctrl->ubfactors[i]);
  Rf_warning("\n");

  Rf_warning("\n");
}


/*************************************************************************/
/*! This function checks the validity of user-supplied parameters */
/*************************************************************************/
int CheckParams(ctrl_t *ctrl)
{
  idx_t i, j;
  real_t sum;
  mdbglvl_et  dbglvl=METIS_DBG_INFO;

  switch (ctrl->optype) {
    case METIS_OP_PMETIS:
      if (ctrl->objtype != METIS_OBJTYPE_CUT) {
        IFSET(dbglvl, METIS_DBG_INFO, Rf_warning("Input Error: Incorrect objective type.\n"));
        return 0;
      }
      if (ctrl->ctype != METIS_CTYPE_RM && ctrl->ctype != METIS_CTYPE_SHEM) {
        IFSET(dbglvl, METIS_DBG_INFO, Rf_warning("Input Error: Incorrect coarsening scheme.\n"));
        return 0;
      }
      if (ctrl->iptype != METIS_IPTYPE_GROW && ctrl->iptype != METIS_IPTYPE_RANDOM) {
        IFSET(dbglvl, METIS_DBG_INFO, Rf_warning("Input Error: Incorrect initial partitioning scheme.\n"));
        return 0;
      }
      if (ctrl->rtype != METIS_RTYPE_FM) {
        IFSET(dbglvl, METIS_DBG_INFO, Rf_warning("Input Error: Incorrect refinement scheme.\n"));
        return 0;
      }
      if (ctrl->ncuts <= 0) {
        IFSET(dbglvl, METIS_DBG_INFO, Rf_warning("Input Error: Incorrect ncuts.\n"));
        return 0;
      }
      if (ctrl->niter <= 0) {
        IFSET(dbglvl, METIS_DBG_INFO, Rf_warning("Input Error: Incorrect niter.\n"));
        return 0;
      }
      if (ctrl->ufactor <= 0) {
        IFSET(dbglvl, METIS_DBG_INFO, Rf_warning("Input Error: Incorrect ufactor.\n"));
        return 0;
      }
      if (ctrl->numflag != 0 && ctrl->numflag != 1) {
        IFSET(dbglvl, METIS_DBG_INFO, Rf_warning("Input Error: Incorrect numflag.\n"));
        return 0;
      }
      if (ctrl->nparts <= 0) {
        IFSET(dbglvl, METIS_DBG_INFO, Rf_warning("Input Error: Incorrect nparts.\n"));
        return 0;
      }
      if (ctrl->ncon <= 0) {
        IFSET(dbglvl, METIS_DBG_INFO, Rf_warning("Input Error: Incorrect ncon.\n"));
        return 0;
      }

      for (i=0; i<ctrl->ncon; i++) {
        sum = rsum(ctrl->nparts, ctrl->tpwgts+i, ctrl->ncon);
        if (sum < 0.99 || sum > 1.01) {
          IFSET(dbglvl, METIS_DBG_INFO, 
              Rf_warning("Input Error: Incorrect sum of %"PRREAL" for tpwgts for constraint %"PRIDX".\n", sum, i));
          return 0;
        }
      }
      for (i=0; i<ctrl->ncon; i++) {
        for (j=0; j<ctrl->nparts; j++) {
          if (ctrl->tpwgts[j*ctrl->ncon+i] <= 0.0) {
            IFSET(dbglvl, METIS_DBG_INFO, 
                Rf_warning("Input Error: Incorrect tpwgts for partition %"PRIDX" and constraint %"PRIDX".\n", j, i));
            return 0;
          }
        }
      }

      for (i=0; i<ctrl->ncon; i++) {
        if (ctrl->ubfactors[i] <= 1.0) {
          IFSET(dbglvl, METIS_DBG_INFO, 
              Rf_warning("Input Error: Incorrect ubfactor for constraint %"PRIDX".\n", i));
          return 0;
        }
      }

      break;

    case METIS_OP_KMETIS:
      if (ctrl->objtype != METIS_OBJTYPE_CUT && ctrl->objtype != METIS_OBJTYPE_VOL) {
        IFSET(dbglvl, METIS_DBG_INFO, Rf_warning("Input Error: Incorrect objective type.\n"));
        return 0;
      }
      if (ctrl->ctype != METIS_CTYPE_RM && ctrl->ctype != METIS_CTYPE_SHEM) {
        IFSET(dbglvl, METIS_DBG_INFO, Rf_warning("Input Error: Incorrect coarsening scheme.\n"));
        return 0;
      }
      if (ctrl->iptype != METIS_IPTYPE_METISRB) {
        IFSET(dbglvl, METIS_DBG_INFO, Rf_warning("Input Error: Incorrect initial partitioning scheme.\n"));
        return 0;
      }
      if (ctrl->rtype != METIS_RTYPE_GREEDY) {
        IFSET(dbglvl, METIS_DBG_INFO, Rf_warning("Input Error: Incorrect refinement scheme.\n"));
        return 0;
      }
      if (ctrl->ncuts <= 0) {
        IFSET(dbglvl, METIS_DBG_INFO, Rf_warning("Input Error: Incorrect ncuts.\n"));
        return 0;
      }
      if (ctrl->niter <= 0) {
        IFSET(dbglvl, METIS_DBG_INFO, Rf_warning("Input Error: Incorrect niter.\n"));
        return 0;
      }
      if (ctrl->ufactor <= 0) {
        IFSET(dbglvl, METIS_DBG_INFO, Rf_warning("Input Error: Incorrect ufactor.\n"));
        return 0;
      }
      if (ctrl->numflag != 0 && ctrl->numflag != 1) {
        IFSET(dbglvl, METIS_DBG_INFO, Rf_warning("Input Error: Incorrect numflag.\n"));
        return 0;
      }
      if (ctrl->nparts <= 0) {
        IFSET(dbglvl, METIS_DBG_INFO, Rf_warning("Input Error: Incorrect nparts.\n"));
        return 0;
      }
      if (ctrl->ncon <= 0) {
        IFSET(dbglvl, METIS_DBG_INFO, Rf_warning("Input Error: Incorrect ncon.\n"));
        return 0;
      }
      if (ctrl->contig != 0 && ctrl->contig != 1) {
        IFSET(dbglvl, METIS_DBG_INFO, Rf_warning("Input Error: Incorrect contig.\n"));
        return 0;
      }
      if (ctrl->minconn != 0 && ctrl->minconn != 1) {
        IFSET(dbglvl, METIS_DBG_INFO, Rf_warning("Input Error: Incorrect minconn.\n"));
        return 0;
      }

      for (i=0; i<ctrl->ncon; i++) {
        sum = rsum(ctrl->nparts, ctrl->tpwgts+i, ctrl->ncon);
        if (sum < 0.99 || sum > 1.01) {
          IFSET(dbglvl, METIS_DBG_INFO, 
              Rf_warning("Input Error: Incorrect sum of %"PRREAL" for tpwgts for constraint %"PRIDX".\n", sum, i));
          return 0;
        }
      }
      for (i=0; i<ctrl->ncon; i++) {
        for (j=0; j<ctrl->nparts; j++) {
          if (ctrl->tpwgts[j*ctrl->ncon+i] <= 0.0) {
            IFSET(dbglvl, METIS_DBG_INFO, 
                Rf_warning("Input Error: Incorrect tpwgts for partition %"PRIDX" and constraint %"PRIDX".\n", j, i));
            return 0;
          }
        }
      }

      for (i=0; i<ctrl->ncon; i++) {
        if (ctrl->ubfactors[i] <= 1.0) {
          IFSET(dbglvl, METIS_DBG_INFO, 
              Rf_warning("Input Error: Incorrect ubfactor for constraint %"PRIDX".\n", i));
          return 0;
        }
      }

      break;



    case METIS_OP_OMETIS:
      if (ctrl->objtype != METIS_OBJTYPE_NODE) {
        IFSET(dbglvl, METIS_DBG_INFO, Rf_warning("Input Error: Incorrect objective type.\n"));
        return 0;
      }
      if (ctrl->ctype != METIS_CTYPE_RM && ctrl->ctype != METIS_CTYPE_SHEM) {
        IFSET(dbglvl, METIS_DBG_INFO, Rf_warning("Input Error: Incorrect coarsening scheme.\n"));
        return 0;
      }
      if (ctrl->iptype != METIS_IPTYPE_EDGE && ctrl->iptype != METIS_IPTYPE_NODE) {
        IFSET(dbglvl, METIS_DBG_INFO, Rf_warning("Input Error: Incorrect initial partitioning scheme.\n"));
        return 0;
      }
      if (ctrl->rtype != METIS_RTYPE_SEP1SIDED && ctrl->rtype != METIS_RTYPE_SEP2SIDED) {
        IFSET(dbglvl, METIS_DBG_INFO, Rf_warning("Input Error: Incorrect refinement scheme.\n"));
        return 0;
      }
      if (ctrl->nseps <= 0) {
        IFSET(dbglvl, METIS_DBG_INFO, Rf_warning("Input Error: Incorrect nseps.\n"));
        return 0;
      }
      if (ctrl->niter <= 0) {
        IFSET(dbglvl, METIS_DBG_INFO, Rf_warning("Input Error: Incorrect niter.\n"));
        return 0;
      }
      if (ctrl->ufactor <= 0) {
        IFSET(dbglvl, METIS_DBG_INFO, Rf_warning("Input Error: Incorrect ufactor.\n"));
        return 0;
      }
      if (ctrl->numflag != 0 && ctrl->numflag != 1) {
        IFSET(dbglvl, METIS_DBG_INFO, Rf_warning("Input Error: Incorrect numflag.\n"));
        return 0;
      }
      if (ctrl->nparts != 3) {
        IFSET(dbglvl, METIS_DBG_INFO, Rf_warning("Input Error: Incorrect nparts.\n"));
        return 0;
      }
      if (ctrl->ncon != 1) {
        IFSET(dbglvl, METIS_DBG_INFO, Rf_warning("Input Error: Incorrect ncon.\n"));
        return 0;
      }
      if (ctrl->compress != 0 && ctrl->compress != 1) {
        IFSET(dbglvl, METIS_DBG_INFO, Rf_warning("Input Error: Incorrect compress.\n"));
        return 0;
      }
      if (ctrl->ccorder != 0 && ctrl->ccorder != 1) {
        IFSET(dbglvl, METIS_DBG_INFO, Rf_warning("Input Error: Incorrect ccorder.\n"));
        return 0;
      }
      if (ctrl->pfactor < 0.0 ) {
        IFSET(dbglvl, METIS_DBG_INFO, Rf_warning("Input Error: Incorrect pfactor.\n"));
        return 0;
      }

      for (i=0; i<ctrl->ncon; i++) {
        if (ctrl->ubfactors[i] <= 1.0) {
          IFSET(dbglvl, METIS_DBG_INFO, 
              Rf_warning("Input Error: Incorrect ubfactor for constraint %"PRIDX".\n", i));
          return 0;
        }
      }

      break;

    default:
      IFSET(dbglvl, METIS_DBG_INFO, Rf_warning("Input Error: Incorrect optype\n"));
      return 0;
  }

  return 1;
}

  
/*************************************************************************/
/*! This function frees the memory associated with a ctrl_t */
/*************************************************************************/
void FreeCtrl(ctrl_t **r_ctrl)
{
  ctrl_t *ctrl = *r_ctrl;

  FreeWorkSpace(ctrl);

  gk_free((void **)&ctrl->tpwgts, &ctrl->pijbm, 
          &ctrl->ubfactors, &ctrl->maxvwgt, &ctrl, LTERM);

  *r_ctrl = NULL;
}


