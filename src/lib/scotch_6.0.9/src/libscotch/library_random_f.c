/* Copyright 2004,2007,2012,2014,2018 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : library_random_f.c                      **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This file contains the Fortran API for  **/
/**                the random generator handling routines  **/
/**                of the libSCOTCH library.               **/
/**                                                        **/
/**   DATES      : # Version 4.0  : from : 21 nov 2005     **/
/**                                 to     23 nov 2005     **/
/**                # Version 6.0  : from : 08 oct 2012     **/
/**                                 to     25 apr 2018     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define LIBRARY

#include "module.h"
#include "common.h"
#include "scotch.h"

/**************************************/
/*                                    */
/* These routines are the Fortran API */
/* for the random handling routines.  */
/*                                    */
/**************************************/

/*
**
*/

SCOTCH_FORTRAN (                      \
RANDOMPROC, randomproc, (             \
const int * const           procnum), \
(procnum))
{
  SCOTCH_randomProc (*procnum);
}

/*
**
*/

SCOTCH_FORTRAN (              \
RANDOMRESET, randomreset, (void), \
())
{
  SCOTCH_randomReset ();
}

/*
**
*/

SCOTCH_FORTRAN (                    \
RANDOMSEED, randomseed, (           \
const SCOTCH_Num * const  seedptr), \
(seedptr))
{
  SCOTCH_randomSeed (*seedptr);
}
