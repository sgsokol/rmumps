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
#include "mpi.h"
LIBSEQ_INT LIBSEQ_CALL MPI_Init(LIBSEQ_INT *pargc, char ***pargv)
{
  return 0;
}

LIBSEQ_INT LIBSEQ_CALL MPI_Comm_rank( MPI_Comm comm, LIBSEQ_INT *rank)
{
  *rank=0;
  return 0;
}
LIBSEQ_INT LIBSEQ_CALL MPI_Finalize(void)
{
   return 0;
}

