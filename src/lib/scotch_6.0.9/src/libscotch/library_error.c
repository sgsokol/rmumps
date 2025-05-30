/* Copyright 2004,2007,2008,2011 ENSEIRB, INRIA & CNRS
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
/**   NAME       : library_error.c                         **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module provides error handling     **/
/**                routines to process errors generated by **/
/**                the routines of the libSCOTCH library.  **/
/**                                                        **/
/**   DATES      : # Version 3.3  : from : 02 oct 1998     **/
/**                                 to     02 oct 1998     **/
/**                # Version 3.4  : from : 01 nov 2001     **/
/**                                 to     01 nov 2001     **/
/**                # Version 5.0  : from : 06 mar 2008     **/
/**                                 to     24 may 2008     **/
/**                # Version 5.1  : from : 27 sep 2008     **/
/**                                 to     17 jul 2011     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define LIBRARY_ERROR

#include "module.h"
#include "common.h"
#include "scotch.h"

/********************************/
/*                              */
/* The error handling routines. */
/*                              */
/********************************/

static char                 _SCOTCHerrorProgName[32] = "";

/* This routine sets the program name for
** error reporting.
** It returns:
** - VOID  : in all cases.
*/

void
SCOTCH_errorProg (
const char * const          progstr)              /*+ Program name +*/
{
  int                 charnbr;
  const char *        nsrcptr;
  char *              ndstptr;

  nsrcptr = progstr;
  ndstptr = _SCOTCHerrorProgName;
  charnbr = strlen (progstr);
  if (charnbr > 31) {
    _SCOTCHerrorProgName[0] =
    _SCOTCHerrorProgName[1] =
    _SCOTCHerrorProgName[2] = '.';
    ndstptr += 3;
    nsrcptr += charnbr - 28;
    charnbr  = 28;
  }
  strncpy (ndstptr, nsrcptr, charnbr);
  _SCOTCHerrorProgName[31] = '\0';
}

/* This routine prints an error message with
** a variable number of arguments, as printf ()
** does, and exits.
** It returns:
** - void  : in all cases.
*/

void
SCOTCH_errorPrint (
const char * const          errstr,               /*+ printf-like variable argument list */
...)
{
  va_list             errlist;                    /* The argument list of the call */
#ifdef SCOTCH_PTSCOTCH
  int                 proclocnum;
#endif /* SCOTCH_PTSCOTCH */

  /*fprintf  (stderr, "%s", _SCOTCHerrorProgName);*/
#ifdef SCOTCH_PTSCOTCH
  if ((MPI_Initialized (&proclocnum) == MPI_SUCCESS) &&
      (proclocnum != 0)                              &&
      (MPI_Comm_rank (MPI_COMM_WORLD, &proclocnum) == MPI_SUCCESS))
    /*fprintf (stderr, "(%d): ", proclocnum);*/
  else
    /*fprintf (stderr, ": ");*/;
#else /* SCOTCH_PTSCOTCH */
  if (_SCOTCHerrorProgName[0] != '\0')
    /*fprintf  (stderr, ": ");*/;
#endif /* SCOTCH_PTSCOTCH */
  /*fprintf  (stderr, "ERROR: ");*/
  va_start (errlist, errstr);
  /*vfprintf (stderr, errstr, errlist);*/             /* Print arguments */
  Rf_error(errstr, errlist);
  va_end   (errlist);
  /*fprintf  (stderr, "\n");
  fflush   (stderr);*/                              /* In case it has been set to buffered mode */
}

/* This routine prints a warning message with
** a variable number of arguments, as printf ()
** does.
** It returns:
** - VOID  : in all cases.
*/

void
SCOTCH_errorPrintW (
const char * const          errstr,               /*+ printf-like variable argument list */
...)
{
  va_list             errlist;                    /* The argument list of the call */
#ifdef SCOTCH_PTSCOTCH
  int                 proclocnum;
#endif /* SCOTCH_PTSCOTCH */

  /*fprintf  (stderr, "%s", _SCOTCHerrorProgName);*/;
#ifdef SCOTCH_PTSCOTCH
  /*if ((MPI_Initialized (&proclocnum) == MPI_SUCCESS) &&
      (proclocnum != 0)                              &&
      (MPI_Comm_rank (MPI_COMM_WORLD, &proclocnum) == MPI_SUCCESS))
    fprintf (stderr, "(%d): ", proclocnum);
  else
    fprintf (stderr, ": ");*/
#else /* SCOTCH_PTSCOTCH */
  /*if (_SCOTCHerrorProgName[0] != '\0')
    fprintf  (stderr, ": ");*/
#endif /* SCOTCH_PTSCOTCH */
  /*fprintf  (stderr, "WARNING: ");*/
  va_start (errlist, errstr);
  /*vfprintf (stderr, errstr, errlist);*/             /* Print arguments */
  Rf_warning(errstr, errlist);
  va_end   (errlist);
  /*fprintf  (stderr, "\n");
  fflush   (stderr);*/                              /* In case it has been set to buffered mode */
}
