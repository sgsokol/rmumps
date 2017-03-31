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
#include "mumps_scotch_int.h"
#if defined(scotch) || defined(ptscotch)
#  include <stdio.h>
#  include "scotch.h"
#endif
void MUMPS_CALL
MUMPS_SCOTCH_INTSIZE(MUMPS_INT *scotch_intsize)
  {
#    if defined(scotch) || defined(ptscotch)
         *scotch_intsize=8*sizeof(SCOTCH_Num);
#    else
       *scotch_intsize=-99999;
#    endif
}
