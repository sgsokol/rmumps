/*
 *
 *  This file is part of MUMPS 5.1.1, released
 *  on Mon Mar 20 14:34:33 UTC 2017
 *
 *
 *  Copyright 1991-2017 CERFACS, CNRS, ENS Lyon, INP Toulouse, Inria,
 *  University of Bordeaux.
 *
 *  This version of MUMPS is provided to you free of charge. It is
 *  released under the CeCILL-C license:
 *  http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html
 *
 */
/* Utility to automatically get the sizes of Fortran types */
#include "mumps_size.h"
void  MUMPS_CALL MUMPS_SIZE_C(char *a, char *b, MUMPS_INT *diff)
{
    *diff = (MUMPS_INT) (b - a);
}
