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
#ifndef MUMPS_METIS_INT_H
#define MUMPS_METIS_INT_H
#include "mumps_common.h" /* includes mumps_compat.h and mumps_c_types.h */
#define MUMPS_METIS_IDXSIZE \
  F_SYMBOL(metis_idxsize,METIS_IDXSIZE)
void MUMPS_CALL
MUMPS_METIS_IDXSIZE(MUMPS_INT *metis_idx_size);
#endif
