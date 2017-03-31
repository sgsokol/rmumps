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
      MODULE DMUMPS_LR_TYPE
      IMPLICIT NONE
      TYPE LRB_TYPE
        DOUBLE PRECISION,POINTER,DIMENSION(:,:) :: Q,R
        INTEGER :: LRFORM,K,M,N,KSVD
        LOGICAL :: ISLR
      END TYPE LRB_TYPE
      CONTAINS
      SUBROUTINE DEALLOC_LRB(LRB_OUT,KEEP8,IS_FACTOR)
        TYPE(LRB_TYPE), INTENT(INOUT) :: LRB_OUT
        LOGICAL, INTENT(IN) :: IS_FACTOR
        INTEGER(8) :: KEEP8(150)
        INTEGER :: MEM
        MEM = 0
        IF (LRB_OUT%ISLR) THEN
           IF(associated(LRB_OUT%Q)) MEM = MEM + size(LRB_OUT%Q)
           IF(associated(LRB_OUT%R)) MEM = MEM + size(LRB_OUT%R)
        ELSE
           IF(associated(LRB_OUT%Q)) MEM = MEM + size(LRB_OUT%Q)
        ENDIF
        KEEP8(70) = KEEP8(70) + int(MEM,8)
        IF (.NOT.IS_FACTOR) THEN
          KEEP8(71) = KEEP8(71) + int(MEM,8)
        ENDIF
        IF (LRB_OUT%ISLR) THEN
          IF (associated(LRB_OUT%Q)) THEN
            DEALLOCATE (LRB_OUT%Q)
            NULLIFY(LRB_OUT%Q)
          ENDIF
          IF (associated(LRB_OUT%R)) THEN
            DEALLOCATE (LRB_OUT%R)
            NULLIFY(LRB_OUT%R)
          ENDIF
        ELSE
          IF (associated(LRB_OUT%Q)) THEN
            DEALLOCATE (LRB_OUT%Q)
            NULLIFY(LRB_OUT%Q)
          ENDIF
        ENDIF
      END SUBROUTINE DEALLOC_LRB
      SUBROUTINE DEALLOC_BLR_PANEL(BLR_PANEL, NB_BLR, KEEP8, IS_FACTOR)
        INTEGER, INTENT(IN)           :: NB_BLR
        TYPE(LRB_TYPE), INTENT(INOUT) :: BLR_PANEL(:) 
        INTEGER(8) :: KEEP8(150)
        LOGICAL, INTENT(IN) :: IS_FACTOR
        INTEGER :: I
        IF (NB_BLR.GT.0) THEN
          IF (BLR_PANEL(1)%M.NE.0) THEN
            DO I=1, NB_BLR
              CALL DEALLOC_LRB(BLR_PANEL(I), KEEP8, IS_FACTOR)
            ENDDO
          ENDIF
        ENDIF
      END SUBROUTINE DEALLOC_BLR_PANEL
      END MODULE DMUMPS_LR_TYPE
