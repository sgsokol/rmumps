/* Copyright 2004,2007,2012,2014,2016 IPB, Universite de Bordeaux, INRIA & CNRS
**
** This file is part of the Scotch software package for static mapping,
** graph partitioning and sparse matrix ordering.
**
** This software is governed by the CeCILL-C license under French law
** and abiding by the rules of distribution of free software. You can
** use, modify and/or redistribute the software under the terms of the
** CeCILL-C license as circulated by CEA, CNRS and INRIA at the following
** URL: "http://www.cecill.info".
** 
** As a counterpart to the access to the source code and rights to copy,
** modify and redistribute granted by the license, users are provided
** only with a limited warranty and the software's author, the holder of
** the economic rights, and the successive licensors have only limited
** liability.
** 
** In this respect, the user's attention is drawn to the risks associated
** with loading, using, modifying and/or developing or reproducing the
** software by the user in light of its specific status of free software,
** that may mean that it is complicated to manipulate, and that also
** therefore means that it is reserved for developers and experienced
** professionals having in-depth computer knowledge. Users are therefore
** encouraged to load and test the software's suitability as regards
** their requirements in conditions enabling the security of their
** systems and/or data to be ensured and, more generally, to use and
** operate it in the same conditions as regards security.
** 
** The fact that you are presently reading this means that you have had
** knowledge of the CeCILL-C license and that you accept its terms.
*/
/************************************************************/
/**                                                        **/
/**   NAME       : library_random.c                        **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module is the API for the random   **/
/**                generator control routine.              **/
/**                                                        **/
/**   DATES      : # Version 4.0  : from : 15 jan 2005     **/
/**                                 to     15 jun 2005     **/
/**                # Version 6.0  : from : 08 oct 2012     **/
/**                                 to     19 mar 2016     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define LIBRARY

#include "module.h"
#include "common.h"
#include "scotch.h"

/************************************/
/*                                  */
/* These routines are the C API for */
/* the random handling routines.    */
/*                                  */
/************************************/

/*+ This routine loads a random state.
*** It returns:
*** - 0  : if state successfully loaded.
*** - 1  : state cannot be loaded.
*** - 2  : on error.
+*/

int
SCOTCH_randomLoad (
FILE *                      stream)
{
  return (intRandLoad (stream));
}

/*+ This routine saves the random state.
*** It returns:
*** - 0  : if state successfully saved.
*** - 1  : state cannot be saved.
*** - 2  : on error.
+*/

int
SCOTCH_randomSave (
FILE *                      stream)
{
  return (intRandSave (stream));
}

/*+ This routine sets the process number that
*** is used to generate a different seed across
*** all processes.
*** It returns:
*** - void  : in all cases.
+*/

void
SCOTCH_randomProc (
int                         procnum)
{
  intRandProc (procnum);
}

/*+ This routine resets the random generator
*** to simulate a start from scratch.
*** It returns:
*** - void  : in all cases.
+*/

void
SCOTCH_randomReset (void)
{
  intRandReset ();
}

/*+ This routine sets the value of the
*** random seed.
*** It returns:
*** - void  : in all cases.
+*/

void
SCOTCH_randomSeed (
INT                         seedval)
{
  intRandSeed (seedval);
}
