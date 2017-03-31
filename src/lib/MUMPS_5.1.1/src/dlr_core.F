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
      MODULE DMUMPS_LR_CORE
      USE MUMPS_LR_COMMON
      USE DMUMPS_LR_TYPE
      USE DMUMPS_LR_STATS
!$    USE OMP_LIB
      IMPLICIT NONE
      CONTAINS
      SUBROUTINE INIT_LRB(LRB_OUT,K,KSVD,M,N,ISLR)
        TYPE(LRB_TYPE), INTENT(OUT) :: LRB_OUT
        INTEGER,INTENT(IN) :: K,KSVD,M,N
        LOGICAL,INTENT(IN) :: ISLR
C This routine simply initializes a LR block but does NOT allocate it
        LRB_OUT%M = M
        LRB_OUT%N = N
        LRB_OUT%K = K
        LRB_OUT%KSVD = KSVD
        LRB_OUT%ISLR = ISLR
        NULLIFY(LRB_OUT%Q)
        NULLIFY(LRB_OUT%R)
        IF (ISLR) THEN
          LRB_OUT%LRFORM = 1
        ELSE 
          LRB_OUT%LRFORM = 0
        ENDIF
      END SUBROUTINE INIT_LRB
      SUBROUTINE IS_FRONT_BLR_CANDIDATE(INODE, NFRONT, NASS, K486, K489,
     &                    K490, K491, K492, N, LRGROUPS, LRSTATUS)
        INTEGER,INTENT(IN) :: INODE, NFRONT, NASS, K486, K489, K490,
     &                        K491, K492
        INTEGER,INTENT(IN) :: N, LRGROUPS(N)
        INTEGER,INTENT(OUT):: LRSTATUS 
C
C     Local variables
        LOGICAL :: COMPRESS_PANEL, COMPRESS_CB
        COMPRESS_PANEL = .FALSE.
        IF ((K486.GT.0).and.( 
     &        ((K492.LT.0).and.INODE.EQ.abs(K492))
     &        .or.
     &        ( (K492.GT.0).and.(K491.LE.NFRONT)
     &        .and.(K490.LE.NASS)))) THEN
          COMPRESS_PANEL = .TRUE.
C         Compression for NASS =1 is useless
          IF (NASS.EQ.1) COMPRESS_PANEL =.FALSE. 
          IF (LRGROUPS (INODE) .LT. 0) COMPRESS_PANEL = .FALSE.
        ENDIF
        COMPRESS_CB = .FALSE.
        IF ((K492.GT.0).AND.(K489.EQ.1).AND.(NFRONT-NASS.GT.K491)) THEN
          COMPRESS_CB = .TRUE.
        ENDIF
        IF (COMPRESS_PANEL.OR.COMPRESS_CB) THEN
          IF (COMPRESS_CB.AND.(.NOT.COMPRESS_PANEL)) THEN
            LRSTATUS = 1 
          ELSE IF (COMPRESS_PANEL.AND.(.NOT.COMPRESS_CB)) THEN
            LRSTATUS = 2 
          ELSE
            LRSTATUS = 3 
          ENDIF
        ELSE 
         LRSTATUS = 0
        ENDIF
      END SUBROUTINE IS_FRONT_BLR_CANDIDATE
      SUBROUTINE ALLOC_LRB(LRB_OUT,K,KSVD,M,N,ISLR,IFLAG,IERROR,KEEP8)
        TYPE(LRB_TYPE), INTENT(OUT) :: LRB_OUT
        INTEGER,INTENT(IN) :: K,KSVD,M,N
        INTEGER,INTENT(OUT) :: IFLAG, IERROR
        LOGICAL,INTENT(IN) :: ISLR
        INTEGER(8) :: KEEP8(150)
        INTEGER :: MEM, allocok
        DOUBLE PRECISION :: ZERO
        PARAMETER (ZERO = 0.0D0)
        IF (ISLR) THEN
          IF (K.EQ.0) THEN
            nullify(LRB_OUT%Q)
            nullify(LRB_OUT%R)
          ELSE
            allocate(LRB_OUT%Q(M,K),LRB_OUT%R(K,N),stat=allocok)
            IF (allocok > 0) THEN
              IFLAG  = -13
              IERROR = K*(M+N)
              write(*,*) 'Allocation problem in BLR routine ALLOC_LRB:',
     &             ' not enough memory? memory requested = ' , IERROR
              RETURN
            ENDIF
          ENDIF
        ELSE
          allocate(LRB_OUT%Q(M,N),stat=allocok)
          IF (allocok > 0) THEN
            IFLAG  = -13
            IERROR = M*N
            write(*,*) 'Allocation problem in BLR routine ALLOC_LRB:',
     &           ' not enough memory? memory requested = ' , IERROR
            RETURN
          ENDIF
          nullify(LRB_OUT%R)
        ENDIF
        LRB_OUT%M = M
        LRB_OUT%N = N
        LRB_OUT%K = K
        LRB_OUT%KSVD = KSVD
        LRB_OUT%ISLR = ISLR
        IF (ISLR) THEN
          LRB_OUT%LRFORM = 1
        ELSE 
          LRB_OUT%LRFORM = 0
        ENDIF
        IF (ISLR) THEN
          MEM = M*K + N*K
        ELSE
          MEM = M*N
        ENDIF
        KEEP8(70) = KEEP8(70) - int(MEM,8)
        KEEP8(68) = min(KEEP8(70), KEEP8(68))
        KEEP8(71) = KEEP8(71) - int(MEM,8)
        KEEP8(69) = min(KEEP8(71), KEEP8(69))
      END SUBROUTINE ALLOC_LRB
      SUBROUTINE REGROUPING2(CUT, NPARTSASS, NASS,
     &                   NPARTSCB, NCB, IBCKSZ, ONLYCB, K472)
        INTEGER, INTENT(IN) :: IBCKSZ, NASS, NCB
        INTEGER, INTENT(INOUT) :: NPARTSCB, NPARTSASS
        INTEGER, POINTER, DIMENSION(:) :: CUT
        INTEGER, POINTER, DIMENSION(:) :: NEW_CUT
        INTEGER :: I, INEW, MINSIZE, NEW_NPARTSASS
        LOGICAL :: ONLYCB, TRACE
        INTEGER, INTENT(IN) :: K472
        INTEGER :: IBCKSZ2
        ALLOCATE(NEW_CUT(max(NPARTSASS,1)+NPARTSCB+1))
        CALL COMPUTE_BLR_VCS(K472, IBCKSZ2, IBCKSZ, NASS)
        MINSIZE = int(IBCKSZ2 / 2)
        NEW_NPARTSASS = max(NPARTSASS,1)
        IF (.NOT. ONLYCB) THEN
           NEW_CUT(1) = 1
           INEW = 2
           I = 2
           DO WHILE (I .LE. NPARTSASS + 1)
              NEW_CUT(INEW) = CUT(I)
              TRACE = .FALSE.
              IF (NEW_CUT(INEW) - NEW_CUT(INEW-1) .GT. MINSIZE) THEN
                 INEW = INEW + 1
                 TRACE = .TRUE.
              ENDIF
              I = I + 1
           END DO
           IF (TRACE) THEN
              INEW = INEW - 1 
           ELSE
              IF (INEW .NE. 2) THEN
                 NEW_CUT(INEW-1) = NEW_CUT(INEW)
                 INEW = INEW - 1
              ENDIF
           ENDIF
           NEW_NPARTSASS = INEW - 1
        ENDIF
        IF (ONLYCB) THEN
           DO I=1,max(NPARTSASS,1)+1
              NEW_CUT(I) = CUT(I)
           ENDDO
        ENDIF
        IF (NCB .EQ. 0) GO TO 50
        INEW = NEW_NPARTSASS+2
        I = max(NPARTSASS,1) + 2
        DO WHILE (I .LE. max(NPARTSASS,1) + NPARTSCB + 1)
              NEW_CUT(INEW) = CUT(I)
              TRACE = .FALSE.
              IF (NEW_CUT(INEW) - NEW_CUT(INEW-1) .GT. MINSIZE) THEN
                 INEW = INEW + 1
                 TRACE = .TRUE.
              ENDIF
              I = I + 1
        END DO
        IF (TRACE) THEN
           INEW = INEW - 1 
        ELSE
           IF (INEW .NE.  NEW_NPARTSASS+2) THEN
           NEW_CUT(INEW-1) = NEW_CUT(INEW)
              INEW = INEW - 1
           ENDIF
        ENDIF
        NPARTSCB = INEW - 1 - NEW_NPARTSASS
 50     CONTINUE       
        NPARTSASS = NEW_NPARTSASS
        DEALLOCATE(CUT)
        ALLOCATE(CUT(NPARTSASS+NPARTSCB+1))
        DO I=1,NPARTSASS+NPARTSCB+1
           CUT(I) = NEW_CUT(I)
        ENDDO
        DEALLOCATE(NEW_CUT)
      END SUBROUTINE REGROUPING2
      SUBROUTINE DMUMPS_LRGEMM_SCALING(LRB, SCALED, A, LA, POSELTD, 
     &          LD_DIAG, IW2, POSELTT, NFRONT, BLOCK, MAXI_CLUSTER) 
C This routine does the scaling (for the symmetric case) before 
C computing the LR product (done in DMUMPS_LRGEMM3)        
        TYPE(LRB_TYPE),INTENT(IN) :: LRB
        INTEGER(8), intent(in)  :: LA
        DOUBLE PRECISION, intent(inout)  :: A(LA)
        DOUBLE PRECISION, intent(inout), DIMENSION(:,:)  :: SCALED
        INTEGER,INTENT(IN) :: LD_DIAG, NFRONT, IW2(*)
        INTEGER(8), INTENT(IN) :: POSELTD, POSELTT
        INTEGER, INTENT(IN) :: MAXI_CLUSTER
        DOUBLE PRECISION, intent(inout)  :: BLOCK(MAXI_CLUSTER)
        INTEGER :: J, NROWS
        DOUBLE PRECISION :: PIV1, PIV2, OFFDIAG
        IF (LRB%LRFORM.EQ.1) THEN
            NROWS = LRB%K
        ELSE ! Full Rank Block
            NROWS = LRB%M
        ENDIF
        J = 1
        DO WHILE (J <= LRB%N)
            IF (IW2(J) > 0) THEN
                SCALED(1:NROWS,J) = A(POSELTD+LD_DIAG*(J-1)+J-1) 
     &           * SCALED(1:NROWS,J)
                J = J+1
            ELSE !2x2 pivot
                PIV1    = A(POSELTD+LD_DIAG*(J-1)+J-1)
                PIV2    = A(POSELTD+LD_DIAG*J+J)
                OFFDIAG = A(POSELTD+LD_DIAG*(J-1)+J)
                BLOCK(1:NROWS)    = SCALED(1:NROWS,J)
                SCALED(1:NROWS,J) = PIV1 * SCALED(1:NROWS,J)
     &            + OFFDIAG * SCALED(1:NROWS,J+1)
                SCALED(1:NROWS,J+1) = OFFDIAG * BLOCK(1:NROWS)
     &            + PIV2 * SCALED(1:NROWS,J+1)
                 J=J+2
            ENDIF
        END DO
      END SUBROUTINE DMUMPS_LRGEMM_SCALING
      SUBROUTINE DMUMPS_LRGEMM3(TRANSB1, TRANSB2, ALPHA,
     &           LRB1, LRB2, BETA, A, LA, POSELTT, NFRONT, SYM, NIV,
     &           IFLAG, IERROR,
     &           COMPRESS_MID_PRODUCT, TOLEPS, KPERCENT, RANK, BUILDQ,
     &           POSELTD, LD_DIAG, IW2, BLOCK, MAXI_CLUSTER)
        TYPE(LRB_TYPE),INTENT(IN) :: LRB1,LRB2
        INTEGER(8), intent(in)  :: LA
        DOUBLE PRECISION, intent(inout)  :: A(LA)
        INTEGER,INTENT(IN) :: NFRONT, SYM, NIV
        INTEGER,INTENT(OUT) :: IFLAG, IERROR
        INTEGER(8), INTENT(IN) :: POSELTT
        INTEGER(8), INTENT(IN), OPTIONAL :: POSELTD
        INTEGER,INTENT(IN), OPTIONAL :: LD_DIAG, IW2(*)
        INTEGER, INTENT(IN), OPTIONAL :: MAXI_CLUSTER
        CHARACTER(len=1),INTENT(IN) :: TRANSB1, TRANSB2
        INTEGER,intent(in) :: COMPRESS_MID_PRODUCT, KPERCENT
        DOUBLE PRECISION, intent(in) :: TOLEPS
        DOUBLE PRECISION :: ALPHA,BETA
        DOUBLE PRECISION, intent(inout), OPTIONAL  :: BLOCK(:)
        DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: XY_YZ
        DOUBLE PRECISION, ALLOCATABLE, TARGET, DIMENSION(:,:) :: XQ, R_Y
        DOUBLE PRECISION, POINTER, DIMENSION(:,:) :: X, Y, Y1, Y2, Z
        CHARACTER(len=1) :: SIDE, TRANSX, TRANSY, TRANSZ
        INTEGER :: M_X, K_XY, K_YZ, N_Z, LDX, LDY, LDY1, LDY2, LDZ, K_Y
        INTEGER :: I, J, RANK, MAXRANK, INFO, LWORK
        LOGICAL :: BUILDQ
        DOUBLE PRECISION,    ALLOCATABLE :: RWORK_RRQR(:)
        DOUBLE PRECISION, ALLOCATABLE :: WORK_RRQR(:), TAU_RRQR(:), 
     &                          Y_RRQR(:,:)
        INTEGER, ALLOCATABLE :: JPVT_RRQR(:)
        INTEGER :: T1, T2, CR
        INTEGER :: allocok, MREQ
        DOUBLE PRECISION :: LOC_UPDT_TIME_OUT
        DOUBLE PRECISION, EXTERNAL ::dnrm2
        DOUBLE PRECISION :: ONE, MONE, ZERO
        PARAMETER (ONE = 1.0D0, MONE=-1.0D0)
        PARAMETER (ZERO=0.0D0)
        IF (LRB2%M.EQ.0) THEN
          write(*,*) "Internal error in DMUMPS_LRGEMM3, LRB2%M=0"
          CALL MUMPS_ABORT() 
        ENDIF
        IF ((SYM.NE.0).AND.((TRANSB1.NE.'N').OR.(TRANSB2.NE.'T'))) THEN
            WRITE(*,*) "SYM > 0 and (", TRANSB1, ",", TRANSB2,
     &                ") parameters found. Symmetric LRGEMM is only ",
     &                 "compatible with (N,T) parameters"
            CALL MUMPS_ABORT()
        ENDIF
        RANK = 0
        BUILDQ = .FALSE.
        IF ((LRB1%LRFORM==1).AND.(LRB2%LRFORM==1)) THEN 
            IF ((LRB1%K.EQ.0).OR.(LRB2%K.EQ.0)) GOTO 700
            allocate(Y(LRB1%K,LRB2%K),stat=allocok)
            IF (allocok > 0) THEN
              MREQ = LRB1%K*LRB2%K
              GOTO 860
            ENDIF
            IF (TRANSB1 == 'N') THEN
                X    => LRB1%Q
                LDX  =  LRB1%M
                M_X  =  LRB1%M 
                K_Y  =  LRB1%N
                IF (SYM .EQ. 0) THEN
                    Y1  => LRB1%R
                ELSE
                    allocate(Y1(LRB1%K,LRB1%N),stat=allocok)
                    IF (allocok > 0) THEN
                      MREQ = LRB1%K*LRB1%N
                      GOTO 860
                    ENDIF
                    DO J=1,LRB1%N
                        DO I=1,LRB1%K
                            Y1(I,J) = LRB1%R(I,J)
                        ENDDO
                    ENDDO
                    CALL DMUMPS_LRGEMM_SCALING(LRB1, Y1, A, LA, POSELTD,
     &                     LD_DIAG, IW2, POSELTT, NFRONT, BLOCK, 
     &                     MAXI_CLUSTER) 
                ENDIF
                LDY1 =  LRB1%K
            ELSE  !TRANSB1 == 'T'
                M_X  =  LRB1%N
                X    => LRB1%R
                LDX  =  LRB1%K
                K_Y  =  LRB1%M
                Y1   => LRB1%Q
                LDY1 =  LRB1%M
            ENDIF
            IF (TRANSB2 == 'N') THEN
                Z    => LRB2%R
                LDZ  =  LRB2%K
                N_Z  =  LRB2%N
                Y2   => LRB2%Q
                LDY2 =  LRB2%M
            ELSE  !TRANSB2 == 'T'
                N_Z  =  LRB2%M
                Z    => LRB2%Q
                LDZ  =  LRB2%M
                Y2   => LRB2%R
                LDY2 =  LRB2%K
            ENDIF
            TRANSZ = TRANSB2
            CALL dgemm(TRANSB1 , TRANSB2 , LRB1%K , LRB2%K, K_Y, ONE,
     &            Y1(1,1), LDY1, Y2(1,1), LDY2, ZERO, Y(1,1), LRB1%K )
            BUILDQ = .FALSE.
            IF (COMPRESS_MID_PRODUCT.GE.1) THEN 
                LWORK = MAX(LRB2%K**2, M_X**2)
                allocate(Y_RRQR(LRB1%K,LRB2%K),
     &               WORK_RRQR(LWORK), RWORK_RRQR(2*LRB2%K), 
     &               TAU_RRQR(MIN(LRB1%K,LRB2%K)),
     &               JPVT_RRQR(LRB2%K),stat=allocok)
                IF (allocok > 0) THEN
                  MREQ = LRB1%K*LRB2%K + LWORK + 2*LRB2%K +
     &                   MIN(LRB1%K,LRB2%K) + LRB2%K
                  GOTO 860
                ENDIF
                DO J=1,LRB2%K
                    DO I=1,LRB1%K
                        Y_RRQR(I,J) = Y(I,J)
                    ENDDO
                ENDDO
                MAXRANK = MIN(LRB1%K, LRB2%K)-1
                MAXRANK = max (1, int((MAXRANK*KPERCENT/100)))
                JPVT_RRQR = 0
                CALL DMUMPS_TRUNCATED_RRQR(LRB1%K, LRB2%K, Y_RRQR(1,1),
     &               LRB1%K, JPVT_RRQR(1), TAU_RRQR(1), WORK_RRQR(1),
     &               LRB2%K, RWORK_RRQR(1), TOLEPS, RANK, MAXRANK, INFO)
                IF ((RANK.GT.MAXRANK).OR.(RANK.EQ.0)) THEN 
                    deallocate(Y_RRQR, WORK_RRQR, RWORK_RRQR, TAU_RRQR, 
     &                     JPVT_RRQR)   
                    BUILDQ = .FALSE.
                ELSE
                    BUILDQ = .TRUE.
                ENDIF
                IF (BUILDQ) THEN ! Successfully compressed middle block
                  allocate(XQ(M_X,RANK), R_Y(RANK,LRB2%K),stat=allocok)
                  IF (allocok > 0) THEN
                    MREQ = M_X*RANK + RANK*LRB2%K
                    GOTO 860
                  ENDIF
                    DO J=1, LRB2%K
                       R_Y(1:MIN(RANK,J),JPVT_RRQR(J)) =
     &                   Y_RRQR(1:MIN(RANK,J),J)
                       IF(J.LT.RANK) R_Y(MIN(RANK,J)+1:
     &                   RANK,JPVT_RRQR(J))= ZERO
                    END DO
                    CALL dorgqr 
     &                  (LRB1%K, RANK, RANK, Y_RRQR(1,1),
     &                  LRB1%K, TAU_RRQR(1),  
     &                  WORK_RRQR(1), LWORK, INFO )
                    CALL dgemm(TRANSB1, 'N', M_X, RANK, LRB1%K, ONE,
     &                    X(1,1), LDX, Y_RRQR(1,1), LRB1%K, ZERO, 
     &                    XQ(1,1), M_X)
                    deallocate(Y_RRQR, WORK_RRQR, RWORK_RRQR, TAU_RRQR, 
     &                         JPVT_RRQR)   
                    nullify(X)
                    X      => XQ
                    LDX    =  M_X
                    K_XY   =  RANK
                    TRANSX =  'N'
                    deallocate(Y)
                    nullify(Y)
                    Y      => R_Y
                    LDY    =  RANK
                    K_YZ   =  LRB2%K
                    TRANSY =  'N'
                    SIDE   = 'R'
                ENDIF
            ENDIF
            IF (.NOT.BUILDQ) THEN 
                LDY    = LRB1%K
                K_XY   = LRB1%K
                K_YZ   = LRB2%K
                TRANSX = TRANSB1
                TRANSY = 'N'
                IF (LRB1%K .GE. LRB2%K) THEN
                    SIDE = 'L'
                ELSE ! LRB1%K < LRB2%K
                    SIDE = 'R'
                ENDIF
            ENDIF
        ENDIF
        IF ((LRB1%LRFORM==1).AND.(LRB2%LRFORM==0)) THEN 
            IF (LRB1%K.EQ.0) GOTO 700
            SIDE   =  'R'
            K_XY   =  LRB1%K
            TRANSX =  TRANSB1 
            TRANSY =  TRANSB1 
            Z      => LRB2%Q 
            LDZ    =  LRB2%M
            TRANSZ =  TRANSB2 
            IF (TRANSB1 == 'N') THEN
                X   => LRB1%Q
                LDX =  LRB1%M
                M_X =  LRB1%M
                LDY =  LRB1%K
                IF (SYM .EQ. 0) THEN
                    Y   => LRB1%R
                ELSE
                    allocate(Y(LRB1%K,LRB1%N),stat=allocok)
                    IF (allocok > 0) THEN
                      MREQ = LRB1%K*LRB1%N
                      GOTO 860
                    ENDIF
                    DO J=1,LRB1%N
                        DO I=1,LRB1%K
                            Y(I,J) = LRB1%R(I,J)
                        ENDDO
                    ENDDO
                    CALL DMUMPS_LRGEMM_SCALING(LRB1, Y, A, LA, POSELTD,
     &                     LD_DIAG, IW2, POSELTT, NFRONT, BLOCK, 
     &                     MAXI_CLUSTER) 
                ENDIF
            ELSE ! TRANSB1 == 'T'
                X   => LRB1%R
                LDX =  LRB1%K
                M_X =  LRB1%N
                Y   => LRB1%Q
                LDY =  LRB1%M
            ENDIF
            IF (TRANSB2 == 'N') THEN
                K_YZ = LRB2%M
                N_Z  = LRB2%N
            ELSE ! TRANSB2 == 'T'
                K_YZ = LRB2%N
                N_Z  = LRB2%M
            ENDIF
        ENDIF
        IF ((LRB1%LRFORM==0).AND.(LRB2%LRFORM==1)) THEN 
            IF (LRB2%K.EQ.0) GOTO 700
            SIDE   =  'L'
            K_YZ   =  LRB2%K
            X      => LRB1%Q 
            LDX    =  LRB1%M
            TRANSX =  TRANSB1 
            TRANSY =  TRANSB2 
            TRANSZ =  TRANSB2 
            IF (TRANSB1 == 'N') THEN
                M_X  = LRB1%M
                K_XY = LRB1%N
            ELSE ! TRANSB1 == 'T'
                M_X  = LRB1%N
                K_XY = LRB1%M
            ENDIF
            IF (TRANSB2 == 'N') THEN
                Y   => LRB2%Q
                LDY =  LRB2%M
                Z   => LRB2%R
                LDZ =  LRB2%K
                N_Z =  LRB2%N
            ELSE ! TRANSB2 == 'T'
                IF (SYM .EQ. 0) THEN
                    Y  => LRB2%R
                ELSE ! Symmetric case: column scaling of R2 is done
                    allocate(Y(LRB2%K,LRB2%N),stat=allocok)
                    IF (allocok > 0) THEN
                      MREQ = LRB2%K*LRB2%N
                      GOTO 860
                    ENDIF
                    DO J=1,LRB2%N
                        DO I=1,LRB2%K
                            Y(I,J) = LRB2%R(I,J)
                        ENDDO
                    ENDDO
                    CALL DMUMPS_LRGEMM_SCALING(LRB2, Y, A, LA, POSELTD,
     &                     LD_DIAG, IW2, POSELTT, NFRONT, BLOCK, 
     &                     MAXI_CLUSTER) 
                ENDIF
                LDY =  LRB2%K
                Z   => LRB2%Q
                LDZ =  LRB2%M
                N_Z =  LRB2%M
            ENDIF
        ENDIF
        IF ((LRB1%LRFORM==0).AND.(LRB2%LRFORM==0)) THEN 
            IF (SYM .EQ. 0) THEN
                X => LRB1%Q
            ELSE
                allocate(X(LRB1%M,LRB1%N),stat=allocok)
                IF (allocok > 0) THEN
                  MREQ = LRB1%M*LRB1%N
                  GOTO 860
                ENDIF
                DO J=1,LRB1%N
                    DO I=1,LRB1%M
                        X(I,J) = LRB1%Q(I,J)
                    ENDDO
                ENDDO
                CALL DMUMPS_LRGEMM_SCALING(LRB1, X, A, LA, POSELTD,
     &                   LD_DIAG, IW2, POSELTT, NFRONT, BLOCK, 
     &                   MAXI_CLUSTER) 
            ENDIF
            SIDE   =  'N'
            LDX    =  LRB1%M
            TRANSX =  TRANSB1
            Z      => LRB2%Q
            LDZ    =  LRB2%M
            TRANSZ =  TRANSB2
            IF (TRANSB1 == 'N') THEN
                M_X  = LRB1%M
                K_XY = LRB1%N  
            ELSE ! TRANSB1 == 'T'
                M_X  = LRB1%N
                K_XY = LRB1%M  
            ENDIF
            IF (TRANSB2 == 'N') THEN
                N_Z =  LRB2%N
            ELSE ! TRANSB2 == 'T'
                N_Z =  LRB2%M
            ENDIF
        ENDIF
        IF (SIDE == 'L') THEN ! LEFT: XY_YZ = X*Y; A = XY_YZ*Z
            allocate(XY_YZ(M_X,K_YZ),stat=allocok)
            IF (allocok > 0) THEN
              MREQ = M_X*K_YZ
              GOTO 860
            ENDIF
            CALL dgemm(TRANSX , TRANSY , M_X , K_YZ, K_XY, ONE,
     &             X(1,1), LDX, Y(1,1), LDY, ZERO, XY_YZ(1,1), M_X)
            CALL SYSTEM_CLOCK(T1)
            CALL dgemm('N', TRANSZ, M_X, N_Z, K_YZ, ALPHA,
     &             XY_YZ(1,1), M_X, Z(1,1), LDZ, BETA, A(POSELTT),
     &             NFRONT)
            CALL SYSTEM_CLOCK(T2,CR)
            LOC_UPDT_TIME_OUT = dble(T2-T1)/dble(CR)
            CALL UPDATE_UPDT_TIME_OUT(LOC_UPDT_TIME_OUT)
            deallocate(XY_YZ)
        ELSEIF (SIDE == 'R') THEN ! RIGHT: XY_YZ = Y*Z; A = X*XY_YZ
            allocate(XY_YZ(K_XY,N_Z),stat=allocok)
            IF (allocok > 0) THEN
              MREQ = K_XY*N_Z
              GOTO 860
            ENDIF
            CALL dgemm(TRANSY , TRANSZ , K_XY , N_Z, K_YZ, ONE,
     &             Y(1,1), LDY, Z(1,1), LDZ, ZERO, XY_YZ(1,1), K_XY)
            CALL SYSTEM_CLOCK(T1)
            CALL dgemm(TRANSX, 'N', M_X, N_Z, K_XY, ALPHA,
     &             X(1,1), LDX, XY_YZ(1,1), K_XY, BETA, A(POSELTT),
     &             NFRONT)
            CALL SYSTEM_CLOCK(T2,CR)
            LOC_UPDT_TIME_OUT = dble(T2-T1)/dble(CR)
            CALL UPDATE_UPDT_TIME_OUT(LOC_UPDT_TIME_OUT)
            deallocate(XY_YZ)
        ELSE ! SIDE == 'N' : NONE; A = X*Z
            CALL dgemm(TRANSX, TRANSZ, M_X, N_Z, K_XY, ALPHA,
     &             X(1,1), LDX, Z(1,1), LDZ, BETA, A(POSELTT),
     &             NFRONT)
        ENDIF
        GOTO 870
  860 CONTINUE        
C       Alloc NOT ok!!        
        write(*,*) 'Allocation problem in BLR routine DMUMPS_LRGEMM3: ',
     &    'not enough memory? memory requested = ' , MREQ        
        IFLAG  = - 13
        IERROR = MREQ
        RETURN
  870 CONTINUE       
C       Alloc ok!!        
        IF ((LRB1%LRFORM==0).AND.(LRB2%LRFORM==0)) THEN 
            IF (SYM .NE. 0) deallocate(X)
        ELSEIF ((LRB1%LRFORM==0).AND.(LRB2%LRFORM==1)) THEN 
            IF (SYM .NE. 0) deallocate(Y)
        ELSEIF ((LRB1%LRFORM==1).AND.(LRB2%LRFORM==0)) THEN 
            IF (SYM .NE. 0) deallocate(Y)
        ELSE ! 1 AND 1
            IF ((TRANSB1=='N').AND.(SYM .NE. 0)) deallocate(Y1)
            IF ((COMPRESS_MID_PRODUCT.GE.1).AND.BUILDQ) THEN
                deallocate(XQ)
                deallocate(R_Y)
            ELSE
                deallocate(Y)
            ENDIF
        ENDIF
 700    CONTINUE       
      END SUBROUTINE DMUMPS_LRGEMM3
      SUBROUTINE MAX_CLUSTER(CUT,CUT_SIZE,MAXI_CLUSTER)
        INTEGER, INTENT(IN) :: CUT_SIZE
        INTEGER, intent(out) :: MAXI_CLUSTER
        INTEGER, POINTER, DIMENSION(:) :: CUT
        INTEGER :: I
        MAXI_CLUSTER = 0
        DO I = 1, CUT_SIZE
          IF (CUT(I+1) - CUT(I) .GE. MAXI_CLUSTER) THEN
            MAXI_CLUSTER = CUT(I+1) - CUT(I)
          END IF
        END DO
      END SUBROUTINE MAX_CLUSTER
      END MODULE DMUMPS_LR_CORE
      SUBROUTINE DMUMPS_TRUNCATED_RRQR( M, N, A, LDA, JPVT, TAU, WORK,
     &     LDW, RWORK, TOLEPS, RANK, MAXRANK, INFO)
C     This routine computes a Rank-Revealing QR factorization of a dense
C     matrix A. The factorization is truncated when the absolute value of
C     a diagonal coefficient of the R factor becomes smaller than a
C     prescribed threshold TOLEPS. The resulting partial Q and R factors
C     provide a rank-k approximation of the input matrix A with accuracy
C     TOLEPS.
C     
C     This routine is obtained by merging the LAPACK
C     (http://www.netlib.org/lapack/) CGEQP3 and CLAQPS routines and by
C     applying a minor modification to the outer factorization loop in
C     order to stop computations as soon as possible when the required
C     accuracy is reached.
C
C     The authors of the LAPACK library are:
C     - Univ. of Tennessee 
C     - Univ. of California Berkeley 
C     - Univ. of Colorado Denver 
C     - NAG Ltd. 
      IMPLICIT NONE
      INTEGER            ::  INFO, LDA, LDW, M, N, RANK, MAXRANK
      DOUBLE PRECISION               ::  TOLEPS
      INTEGER            ::  JPVT(*)
      DOUBLE PRECISION, TARGET       ::  RWORK(*)
      DOUBLE PRECISION            ::  A(LDA,*), TAU(*)
      DOUBLE PRECISION, TARGET    ::  WORK(LDW,*)
      INTEGER, PARAMETER ::  INB=1, INBMIN=2
      INTEGER            :: J, JB, MINMN, NB
      INTEGER            :: OFFSET, ITEMP
      INTEGER            :: LSTICC, PVT, K, RK, LDF
      DOUBLE PRECISION               :: TEMP, TEMP2, TOL3Z
      DOUBLE PRECISION, POINTER      :: VN1(:), VN2(:)
      DOUBLE PRECISION, POINTER   :: AUXV(:,:), F(:,:)
      DOUBLE PRECISION            :: AKK
      DOUBLE PRECISION, PARAMETER    :: RZERO=0.0D+0, RONE=1.0D+0
      DOUBLE PRECISION :: ZERO
      DOUBLE PRECISION :: ONE
      PARAMETER          ( ONE = 1.0D+0 )
      PARAMETER          ( ZERO = 0.0D+0 ) 
      DOUBLE PRECISION               :: dlamch
      INTEGER            :: ilaenv, idamax
      EXTERNAL           :: idamax, dlamch
      EXTERNAL           dgeqrf, dormqr, xerbla
      EXTERNAL           ilaenv
      EXTERNAL           dgemm, dgemv, dlarfg, dswap
      DOUBLE PRECISION, EXTERNAL :: dnrm2
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.EQ.0 ) THEN
         IF( LDW.LT.N ) THEN
            INFO = -8
         END IF
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CGEQP3', -INFO )
         RETURN
      END IF
      MINMN = MIN(M,N)
      IF( MINMN.EQ.0 ) THEN
         RETURN
      END IF
      NB = ILAENV( INB, 'CGEQRF', ' ', M, N, -1, -1 )
      VN1  => RWORK(1:N)
      VN2  => RWORK(N+1:2*N)
      AUXV => WORK(1:LDW,1:1)
      F    => WORK(1:LDW,2:NB+1)
      LDF  =  LDW
*     Initialize partial column norms. The first N elements of work
*     store the exact column norms.
      DO J = 1, N
         VN1( J ) = dnrm2( M, A( 1, J ), 1 )
         VN2( J ) = VN1( J )
         JPVT(J) = J
      END DO
      OFFSET = 0
      TOL3Z  = SQRT(dlamch('Epsilon'))
      DO 
         JB     = MIN(NB,MINMN-OFFSET)
         LSTICC = 0
         K      = 0
         DO 
            IF(K.EQ.JB) EXIT
            K   = K+1
            RK  = OFFSET+K
            PVT = ( RK-1 ) + IDAMAX( N-RK+1, VN1( RK ), 1 )
            IF(VN1(PVT).LT.TOLEPS) THEN
               RANK = RK-1
               RETURN
            END IF
            IF (RK.GT.MAXRANK) THEN
               RANK = RK
               INFO = RK
               RETURN
            END IF
            IF( PVT.NE.RK ) THEN
               CALL dswap( M, A( 1, PVT ), 1, A( 1, RK ), 1 )
               CALL dswap( K-1, F( PVT-OFFSET, 1 ), LDF,
     &              F( K, 1 ), LDF )
               ITEMP     = JPVT(PVT)
               JPVT(PVT) = JPVT(RK)
               JPVT(RK)  = ITEMP
               VN1(PVT)  = VN1(RK)
               VN2(PVT)  = VN2(RK)
            END IF
*     Apply previous Householder reflectors to column K:
*     A(RK:M,RK) := A(RK:M,RK) - A(RK:M,OFFSET+1:RK-1)*F(K,1:K-1)**H.
            IF( K.GT.1 ) THEN
               CALL dgemv( 'No transpose', M-RK+1, K-1, -ONE,
     &              A(RK,OFFSET+1), LDA, F(K,1), LDF,
     &              ONE, A(RK,RK), 1 )
            END IF
*     Generate elementary reflector H(k).
            IF( RK.LT.M ) THEN
               CALL dlarfg( M-RK+1, A(RK,RK), A(RK+1,RK), 1, TAU(RK) )
            ELSE
               CALL dlarfg( 1, A(RK,RK), A(RK,RK), 1, TAU(RK) )
            END IF
            AKK      = A(RK,RK)
            A(RK,RK) = ONE
*     Compute Kth column of F:
*     F(K+1:N,K) := tau(K)*A(RK:M,K+1:N)**H*A(RK:M,K).
            IF( RK.LT.N ) THEN
               CALL dgemv( 'Transpose', M-RK+1, N-RK, TAU(RK),
     &              A(RK,RK+1), LDA, A(RK,RK), 1, ZERO,
     &              F( K+1, K ), 1 )
            END IF
*     Padding F(1:K,K) with zeros.
            DO J = 1, K
               F( J, K ) = ZERO
            END DO
*     Incremental updating of F:
*     F(1:N,K) := F(1:N-OFFSET,K) - tau(RK)*F(1:N,1:K-1)*A(RK:M,OFFSET+1:RK-1)**H*A(RK:M,RK).
            IF( K.GT.1 ) THEN
               CALL dgemv( 'Transpose', M-RK+1, K-1, -TAU(RK),
     &              A(RK,OFFSET+1), LDA, A(RK,RK), 1, ZERO,
     &              AUXV(1,1), 1 )
               CALL dgemv( 'No transpose', N-OFFSET, K-1, ONE,
     &              F(1,1), LDF, AUXV(1,1), 1, ONE, F(1,K), 1 )
            END IF
*     Update the current row of A:
*     A(RK,RK+1:N) := A(RK,RK+1:N) - A(RK,OFFSET+1:RK)*F(K+1:N,1:K)**H.
            IF( RK.LT.N ) THEN
               CALL dgemv( 'No Transpose', N-RK, K, -ONE, F( K+1, 1 ), 
     &              LDF,
     &              A( RK, OFFSET+1 ), LDA, ONE, A( RK, RK+1 ), LDA )
            END IF
*     Update partial column norms.
*     
            IF( RK.LT.MINMN ) THEN
               DO J = RK + 1, N
                  IF( VN1( J ).NE.RZERO ) THEN
*     
*     NOTE: The following 4 lines follow from the analysis in
*     Lapack Working Note 176.
*
                     TEMP = ABS( A( RK, J ) ) / VN1( J )
                     TEMP = MAX( RZERO, ( RONE+TEMP )*( RONE-TEMP ) )
                     TEMP2 = TEMP*( VN1( J ) / VN2( J ) )**2
                     IF( TEMP2 .LE. TOL3Z ) THEN
                        VN2( J ) = dble( LSTICC )
                        LSTICC = J
                     ELSE
                        VN1( J ) = VN1( J )*SQRT( TEMP )
                     END IF
                  END IF
               END DO
            END IF
            A( RK, RK ) = AKK
            IF (LSTICC.NE.0) EXIT
         END DO
*     Apply the block reflector to the rest of the matrix:
*     A(RK+1:M,RK+1:N) := A(RK+1:M,RK+1:N) -
*     A(RK+1:M,OFFSET+1:RK)*F(K+1:N-OFFSET,1:K)**H.
         IF( RK.LT.MIN(N,M) ) THEN
            CALL dgemm( 'No transpose', 'Transpose', M-RK,
     &           N-RK, K, -ONE, A(RK+1,OFFSET+1), LDA,
     &           F(K+1,1), LDF, ONE, A(RK+1,RK+1), LDA )
         END IF
*     Recomputation of difficult columns.
         DO WHILE( LSTICC.GT.0 ) 
            ITEMP = NINT( VN2( LSTICC ) )
            VN1( LSTICC ) = dnrm2( M-RK, A( RK+1, LSTICC ), 1 )
*     
*     NOTE: The computation of VN1( LSTICC ) relies on the fact that 
*     SNRM2 does not fail on vectors with norm below the value of
*     SQRT(DLAMCH('S')) 
*     
            VN2( LSTICC ) = VN1( LSTICC )
            LSTICC = ITEMP
         END DO
         IF(RK.GE.MINMN) EXIT
         OFFSET = RK
      END DO
      RANK = RK
      END SUBROUTINE DMUMPS_TRUNCATED_RRQR
