C
C  This file is part of MUMPS 5.1.1, released
C  on Mon Mar 20 14:34:33 UTC 2017
C
C
C  Copyright 1991-2017 CERFACS, CNRS, ENS Lyon, INP Toulouse, Inria,
C  University of Bordeaux.
C
C  This version of MUMPS is provided to you free of charge. It is
C  released under the CeCILL-C license:
C  http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html
C
      SUBROUTINE MUMPS_SET_VERSION( VERSION_STR )
      IMPLICIT NONE
      CHARACTER(LEN=*) :: VERSION_STR
      CHARACTER(LEN=*) :: V;
      PARAMETER (V = "5.1.1" )
      IF ( len(V) .GT. 30 ) THEN
         WRITE(*,*) "Version string too long ( >30 characters )"
         CALL MUMPS_ABORT()
      END IF
      VERSION_STR = V
      RETURN
      END SUBROUTINE MUMPS_SET_VERSION
