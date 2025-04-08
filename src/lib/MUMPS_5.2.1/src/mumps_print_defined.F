C
C  This file is part of MUMPS 5.2.1, released
C  on Fri Jun 14 14:46:05 UTC 2019
C
C
C  Copyright 1991-2019 CERFACS, CNRS, ENS Lyon, INP Toulouse, Inria,
C  Mumps Technologies, University of Bordeaux.
C
C  This version of MUMPS is provided to you free of charge. It is
C  released under the CeCILL-C license:
C  http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html
C
      SUBROUTINE MUMPS_PRINT_IF_DEFINED(MPG)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: MPG
      IF (MPG.LE.0) RETURN
c~       write(MPG,*) "================================================="
c~ #if defined(ZERO_TRIANGLE)
c~       write(MPG,*) "MUMPS compiled with option -DZERO_TRIANGLE"
c~ #endif
c~ #if defined(GEMMT_AVAILABLE)
c~       write(MPG, *) "MUMPS compiled with option -DGEMMT_AVAILABLE"
c~ #endif
c~ #if defined(DETERMINISTIC_PARALLEL_GRAPH)
c~       write(MPG,*) "MUMPS compiled with option"
c~      &     ," -DDETERMINISTIC_PARALLEL_GRAPH"
c~ #endif
c~ #if defined(metis)
c~       write(MPG,*) "MUMPS compiled with option -Dmetis"
c~ #endif
c~ #if defined(metis4)
c~       write(MPG,*) "MUMPS compiled with option -Dmetis4"
c~ #endif
c~ #if defined(MUMPS_F2003)
c~       write(MPG,*) "MUMPS compiled with option -DMUMPS_F2003"
c~ #endif
c~ #if defined(OLD_OOC_NOPANEL)
c~       write(MPG,*) "MUMPS compiled with option -DOLD_OOC_NOPANEL"
c~ #endif
c~ #if defined(parmetis)
c~       write(MPG,*) "MUMPS compiled with option -Dparmetis"
c~ #endif
c~ #if defined(parmetis3)
c~       write(MPG,*) "MUMPS compiled with option -Dparmetis3"
c~ #endif
c~ #if defined(ptscotch)
c~       write(MPG,*) "MUMPS compiled with option -Dptscotch"
c~ #endif
c~ #if defined(scotch)
c~       write(MPG,*) "MUMPS compiled with option -Dscotch"
c~ #endif
c~ #if defined(MUMPS_USE_BLAS2)
c~       write(MPG,*) "MUMPS compiled with option -DMUMPS_USE_BLAS2"
c~ #endif
c~       write(MPG,*) "This MUMPS version includes code for SAVE_RESTORE"
c~ #if defined(NODYNAMICCB)
c~       write(MPG,*) "MUMPS compiled with option -DNODYNAMICCB"
c~ #endif
c~ #if defined(BLR_MT)
c~       write(MPG,*) "MUMPS compiled with option -DBLR_MT"
c~ #endif
c~       write(MPG,*) "================================================="
      RETURN
      END SUBROUTINE MUMPS_PRINT_IF_DEFINED
