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
      RECURSIVE SUBROUTINE DMUMPS_TRAITER_MESSAGE_SOLVE
     &     ( BUFR, LBUFR, LBUFR_BYTES,
     &     MSGTAG, MSGSOU, MYID, SLAVEF, COMM,
     &     N, NRHS, IPOOL, LPOOL, III, LEAF,
     &     NBFIN, NSTK_S, IW, LIW, A, LA, PTRIST,
     &     PTRFAC, IWCB, LIWCB,
     &     WCB, LWCB, POSWCB,
     &     PLEFTWCB, POSIWCB,
     &     PTRICB,
     &     INFO, KEEP, KEEP8, DKEEP, STEP, PROCNODE_STEPS, 
     &     RHSCOMP, LRHSCOMP, POSINRHSCOMP_FWD
     &     , FROM_PP )
      USE DMUMPS_OOC 
      USE DMUMPS_BUF 
      IMPLICIT NONE
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER MSGTAG, MSGSOU, MYID, SLAVEF, COMM
      INTEGER LIW
      INTEGER(8), INTENT(IN) :: LA, LWCB
      INTEGER N, NRHS, LPOOL, III, LEAF, NBFIN, LRHSCOMP
      INTEGER LIWCB, POSIWCB
      INTEGER(8) :: POSWCB, PLEFTWCB
      INTEGER INFO( 40 ), KEEP( 500)
      INTEGER(8) KEEP8(150)
      DOUBLE PRECISION, INTENT(INOUT) :: DKEEP(230)
      INTEGER BUFR( LBUFR )
      INTEGER IPOOL( LPOOL ),  NSTK_S( N )
      INTEGER IWCB( LIWCB )
      INTEGER IW( LIW )
      INTEGER PTRICB(KEEP(28)),PTRIST(KEEP(28))
      INTEGER(8) :: PTRFAC(KEEP(28))
      INTEGER STEP(N)
      INTEGER PROCNODE_STEPS(KEEP(28))
      DOUBLE PRECISION WCB( LWCB ), A( LA )
#if defined(RHSCOMP_BYROWS)
      DOUBLE PRECISION RHSCOMP( NRHS, LRHSCOMP )
#else
      DOUBLE PRECISION RHSCOMP( LRHSCOMP, NRHS )
#endif
      INTEGER, intent(in) :: POSINRHSCOMP_FWD(N)
      LOGICAL, intent(in) :: FROM_PP
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER(8) :: PTRX, PTRY, IFR8
      INTEGER IERR, K, JJ, JBDEB, JBFIN, NRHS_B
      INTEGER FINODE, FPERE, LONG, NCB, POSITION, NCV, NPIV
      INTEGER PDEST, I, IPOSINRHSCOMP
      INTEGER J1
      INTEGER(8) :: APOS
      LOGICAL DUMMY
      LOGICAL FLAG
!$    LOGICAL :: OMP_FLAG
      EXTERNAL MUMPS_PROCNODE
      INTEGER  MUMPS_PROCNODE
      DOUBLE PRECISION :: TIME_TMP
      DOUBLE PRECISION ALPHA, ONE
      PARAMETER (ONE = 1.0D0, ALPHA=-1.0D0)
      INCLUDE 'mumps_headers.h'
      IF ( MSGTAG .EQ. RACINE_SOLVE ) THEN
         NBFIN = NBFIN - 1
         IF ( NBFIN .eq. 0 ) GOTO 270
      ELSE  IF (MSGTAG .EQ. ContVec ) THEN
#if defined(ALLTIMERS)
         TIME_TMP = MPI_WTIME()
#endif
         POSITION = 0
         CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &        FINODE, 1, MPI_INTEGER, COMM, IERR )
         CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &        FPERE, 1, MPI_INTEGER, COMM, IERR )
         CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &        NCB, 1, MPI_INTEGER, COMM, IERR )
         CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &        JBDEB, 1, MPI_INTEGER, COMM, IERR )
         CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &        JBFIN, 1, MPI_INTEGER, COMM, IERR )
         CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &        LONG, 1, MPI_INTEGER, COMM, IERR )
         NRHS_B = JBFIN-JBDEB+1
          IF ( NCB .eq. 0 ) THEN
             PTRICB(STEP(FINODE)) = -1
             NSTK_S(STEP(FPERE)) = NSTK_S(STEP(FPERE)) - 1
             IF ( NSTK_S(STEP(FPERE)) .EQ. 0 ) THEN
                   IPOOL( LEAF ) = FPERE
                LEAF = LEAF + 1
                IF ( LEAF > LPOOL ) THEN
                   WRITE(*,*) 'Internal error 41r2 : Pool is too small.'
                   CALL MUMPS_ABORT()
                END IF
             END IF
          ELSE
             IF ( PTRICB(STEP(FINODE)) .EQ. 0 ) THEN
                PTRICB(STEP(FINODE)) = NCB + 1
             END IF
             IF ( ( POSIWCB - LONG ) .LT. 0 ) THEN
                INFO( 1 ) = -14
                INFO( 2 ) = LONG
                GOTO 260
             END IF
             IF ( POSWCB - PLEFTWCB + 1_8 .LT.
     &            int(LONG,8) * int(NRHS_B,8)) THEN
                INFO( 1 ) = -11
                CALL MUMPS_SET_IERROR(PLEFTWCB-POSWCB-1_8+
     &              int(LONG,8) * int(NRHS_B,8),
     &              INFO(2))
                GOTO 260
             END IF
             IF (LONG .GT. 0) THEN
                CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &               IWCB( 1 ),
     &               LONG, MPI_INTEGER, COMM, IERR )
                DO K = 1, NRHS_B
                   CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                  WCB( PLEFTWCB ),
     &                  LONG, MPI_DOUBLE_PRECISION, COMM, IERR )
                   DO I = 1, LONG
                   IPOSINRHSCOMP= abs(POSINRHSCOMP_FWD(IWCB(I)))
#if defined(RHSCOMP_BYROWS)
                   RHSCOMP(JBDEB+K-1,IPOSINRHSCOMP) = 
     &             RHSCOMP(JBDEB+K-1,IPOSINRHSCOMP) + 
     &                 WCB(PLEFTWCB+I-1)
#else
                   RHSCOMP(IPOSINRHSCOMP,JBDEB+K-1) = 
     &             RHSCOMP(IPOSINRHSCOMP,JBDEB+K-1) + 
     &                 WCB(PLEFTWCB+I-1)
#endif
                   ENDDO
                END DO
                PTRICB(STEP(FINODE)) = PTRICB(STEP(FINODE)) - LONG
             ENDIF
             IF ( PTRICB(STEP(FINODE)) == 1 ) THEN
                NSTK_S(STEP(FPERE)) = NSTK_S(STEP(FPERE)) - 1
             END IF
             IF ( NSTK_S(STEP(FPERE)) .EQ. 0 ) THEN
                   IPOOL( LEAF ) = FPERE
                LEAF = LEAF + 1
                IF ( LEAF > LPOOL ) THEN
                   WRITE(*,*) 'Internal error 41r2 : Pool is too small.'
                   CALL MUMPS_ABORT()
                END IF
             ENDIF
          END IF
#if defined(ALLTIMERS)
          DKEEP(128)=DKEEP(128)+MPI_WTIME()-TIME_TMP
#endif
       ELSEIF ( MSGTAG .EQ. Master2Slave ) THEN
#if defined(ALLTIMERS)
          TIME_TMP=MPI_WTIME()
#endif
          POSITION = 0
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &         FINODE, 1, MPI_INTEGER, COMM, IERR )
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &         FPERE, 1, MPI_INTEGER, COMM, IERR )
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &         NCV, 1, MPI_INTEGER, COMM, IERR )
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &         NPIV, 1, MPI_INTEGER, COMM, IERR )
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &         JBDEB, 1, MPI_INTEGER, COMM, IERR )
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &         JBFIN, 1, MPI_INTEGER, COMM, IERR )
          NRHS_B = JBFIN-JBDEB+1
          PTRY = PLEFTWCB
          PTRX = PLEFTWCB + int(NCV,8) * int(NRHS_B,8)
          PLEFTWCB = PLEFTWCB + (NPIV + NCV) * NRHS_B
          IF ( POSWCB - PLEFTWCB + 1 .LT. 0 ) THEN
             INFO(1) = -11
             CALL MUMPS_SET_IERROR(-POSWCB+PLEFTWCB-1_8,INFO(2))
             GO TO 260
          END IF
          DO K=1, NRHS_B
             CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &            WCB( PTRY + (K-1) * NCV ), NCV,
     &            MPI_DOUBLE_PRECISION, COMM, IERR )
          ENDDO
          IF ( NPIV .GT. 0 ) THEN
             DO K=1, NRHS_B
                CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &               WCB( PTRX + (K-1)*NPIV ), NPIV,
     &               MPI_DOUBLE_PRECISION, COMM, IERR )
             END DO
          END IF
#if defined(ALLTIMERS)
          DKEEP(129)=DKEEP(129)+MPI_WTIME()-TIME_TMP
#endif
          IF (KEEP(201).GT.0) THEN
             CALL DMUMPS_SOLVE_GET_OOC_NODE(
     &            FINODE,PTRFAC,KEEP,A,LA,STEP,
     &            KEEP8,N,DUMMY,IERR)
             IF(IERR.LT.0)THEN
                INFO(1)=IERR
                INFO(2)=0
                GOTO 260
             ENDIF
          ENDIF
          APOS = PTRFAC(STEP(FINODE))
#if defined(ALLTIMERS)
          TIME_TMP=MPI_WTIME()
#endif
          IF (KEEP(201).EQ.1) THEN
#if defined(MUMPS_USE_BLAS2)
             IF ( NRHS_B == 1 ) THEN
                CALL dgemv( 'N', NCV, NPIV, ALPHA, A(APOS), NCV,
     &               WCB( PTRX ), 1, ONE,
     &               WCB( PTRY ), 1 )
             ELSE
#endif
                CALL dgemm( 'N', 'N', NCV, NRHS_B, NPIV, ALPHA,
     &               A(APOS), NCV,
     &               WCB( PTRX), NPIV, ONE,
     &               WCB( PTRY), NCV )
#if defined(MUMPS_USE_BLAS2)
             ENDIF
#endif
          ELSE                  
#if defined(MUMPS_USE_BLAS2)
             IF ( NRHS_B == 1 ) THEN
                CALL dgemv( 'T', NPIV, NCV, ALPHA, A(APOS), NPIV,
     &               WCB( PTRX ), 1, ONE,
     &               WCB( PTRY ), 1 )
             ELSE
#endif
                CALL dgemm( 'T', 'N', NCV, NRHS_B, NPIV, ALPHA,
     &               A(APOS), NPIV,
     &               WCB( PTRX), NPIV, ONE,
     &               WCB( PTRY), NCV )
#if defined(MUMPS_USE_BLAS2)
             ENDIF
#endif
          ENDIF
#if defined(ALLTIMERS)
          DKEEP(130)=DKEEP(130)+MPI_WTIME()-TIME_TMP
#endif
          IF (KEEP(201).GT.0) THEN
             CALL DMUMPS_FREE_FACTORS_FOR_SOLVE(FINODE,PTRFAC,
     &            KEEP(28),A,LA,.TRUE.,IERR)
             IF(IERR.LT.0)THEN
                INFO(1)=IERR
                INFO(2)=0
                GOTO 260
             ENDIF
          ENDIF
          PLEFTWCB = PLEFTWCB - NPIV * NRHS_B
          PDEST = MUMPS_PROCNODE( PROCNODE_STEPS(STEP(FPERE)),
     &                            SLAVEF )
          IF ( PDEST .EQ. MYID ) THEN
#if defined(ALLTIMERS)
             TIME_TMP=MPI_WTIME()
#endif
             IF ( PTRICB(STEP(FINODE)) .EQ. 0 ) THEN
                NCB = IW( PTRIST(STEP(FINODE)) + 2 + KEEP(IXSZ) )
                PTRICB(STEP(FINODE)) = NCB + 1
             END IF
             IF (KEEP(350).EQ.0) THEN
              DO I = 1, NCV
                JJ=IW(PTRIST(STEP(FINODE))+3+I+ KEEP(IXSZ) )
                IPOSINRHSCOMP= abs(POSINRHSCOMP_FWD(JJ))
                DO K=1, NRHS_B
#if defined(RHSCOMP_BYROWS)
                   RHSCOMP(JBDEB+K-1,IPOSINRHSCOMP)= 
     &              RHSCOMP(JBDEB+K-1,IPOSINRHSCOMP) +
     &              WCB(PTRY+I-1+(K-1)*NCV)
#else
                   RHSCOMP(IPOSINRHSCOMP,JBDEB+K-1)= 
     &              RHSCOMP(IPOSINRHSCOMP,JBDEB+K-1) +
     &              WCB(PTRY+I-1+(K-1)*NCV)
#endif
                ENDDO
              END DO
             ELSE IF (KEEP(350).EQ.1.OR.KEEP(350).EQ.2) THEN
                J1 = PTRIST(STEP(FINODE))+3+KEEP(IXSZ) 
!$              OMP_FLAG = ( JBFIN-JBDEB+1.GE.KEEP(362) .AND.
!$   &                     (NCV*(JBFIN-JBDEB+1) .GE. KEEP(363) ) )
!$OMP PARALLEL DO PRIVATE(I,JJ,IFR8,IPOSINRHSCOMP) IF(OMP_FLAG)
                DO K=1, NRHS_B
                  IFR8 = PTRY+int(K-1,8)*int(NCV,8)
                  DO I = 1,NCV
                  JJ = IW(J1+I)
                  IPOSINRHSCOMP= abs(POSINRHSCOMP_FWD(JJ))
#if defined(RHSCOMP_BYROWS)
                  RHSCOMP(JBDEB+K-1,IPOSINRHSCOMP)= 
     &              RHSCOMP(JBDEB+K-1,IPOSINRHSCOMP)
#else
                  RHSCOMP(IPOSINRHSCOMP,JBDEB+K-1)= 
     &              RHSCOMP(IPOSINRHSCOMP,JBDEB+K-1)
#endif
     &              + WCB(IFR8+int(I-1,8))
                  ENDDO
                ENDDO
!$OMP END PARALLEL DO
              ELSE
                  WRITE(*,*) "UNKNOWN VERSION OF KEEP(350)"
                  CALL MUMPS_ABORT()
              ENDIF
#if defined(ALLTIMERS)
             DKEEP(131)=DKEEP(131)+MPI_WTIME()-TIME_TMP
#endif
             PTRICB(STEP(FINODE)) =
     &            PTRICB(STEP(FINODE)) - NCV
             IF ( PTRICB( STEP( FINODE ) ) == 1 ) THEN
                NSTK_S(STEP(FPERE)) = NSTK_S(STEP(FPERE)) - 1
             END IF
             IF ( NSTK_S(STEP(FPERE)) .EQ. 0 ) THEN
                   IPOOL( LEAF ) = FPERE
                LEAF = LEAF + 1
                IF ( LEAF > LPOOL ) THEN
                   WRITE(*,*) 'INTERNAL Error 41r: Pool is too small.'
                   CALL MUMPS_ABORT()
                END IF
             ENDIF
          ELSE
 210         CONTINUE
#if defined(ALLTIMERS)
             TIME_TMP=MPI_WTIME()
#endif
             CALL DMUMPS_BUF_SEND_VCB( NRHS_B, FINODE, FPERE,
     &            IW(PTRIST(STEP( FINODE )) + 2 + KEEP(IXSZ) ), NCV,NCV,
     &            IW(PTRIST(STEP(FINODE))+4+ KEEP(IXSZ) ),
     &            WCB( PTRY ), JBDEB, JBFIN,
     &            RHSCOMP, 1, 1, -9999, -9999,
     &            KEEP, PDEST, ContVec, COMM, IERR )
             IF (IERR.EQ.0) THEN
#if defined(ALLTIMERS)
               DKEEP(132)=DKEEP(132)+MPI_WTIME()-TIME_TMP
#endif
             ELSE
#if defined(ALLTIMERS)
               DKEEP(134)=DKEEP(134)+MPI_WTIME()-TIME_TMP
#endif
             ENDIF
             IF ( IERR .EQ. -1 ) THEN
                CALL DMUMPS_SOLVE_RECV_AND_TREAT( .FALSE., FLAG,
     &               BUFR, LBUFR, LBUFR_BYTES,
     &               MYID, SLAVEF, COMM,
     &               N, NRHS, IPOOL, LPOOL, III, LEAF,
     &               NBFIN, NSTK_S, IW, LIW, A, LA, PTRIST, PTRFAC,
     &               IWCB, LIWCB,
     &               WCB, LWCB, POSWCB, PLEFTWCB, POSIWCB,
     &               PTRICB, INFO, KEEP,KEEP8, DKEEP, STEP,
     &               PROCNODE_STEPS, 
     &               RHSCOMP, LRHSCOMP, POSINRHSCOMP_FWD
     &               , FROM_PP )
                IF ( INFO( 1 )  .LT. 0 )  GOTO 270
                GOTO 210
             ELSE IF ( IERR .EQ. -2 ) THEN
                INFO( 1 ) = -17
                INFO( 2 ) = ( NCV + 4 ) * KEEP( 34 ) +
     &               NCV * KEEP( 35 )
                GOTO 260
             ELSE IF ( IERR .EQ. -3 ) THEN
                INFO( 1 ) = -20
                INFO( 2 ) = ( NCV + 4 ) * KEEP( 34 ) +
     &               NCV * KEEP( 35 )
             END IF
          END IF
          PLEFTWCB = PLEFTWCB - NCV * NRHS_B
       ELSEIF ( MSGTAG .EQ. TERREUR ) THEN
          INFO(1) = -001
          INFO(2) = MSGSOU
          GOTO 270
       ELSE IF ( (MSGTAG.EQ.UPDATE_LOAD).OR.
     &         (MSGTAG.EQ.TAG_DUMMY) ) THEN
          GO TO 270
       ELSE
          INFO(1)=-100
          INFO(2)=MSGTAG
          GO TO 260
       ENDIF
       GO TO 270
 260   CONTINUE
       CALL DMUMPS_BDC_ERROR( MYID, SLAVEF, COMM, KEEP )
 270   CONTINUE
       RETURN
       END SUBROUTINE DMUMPS_TRAITER_MESSAGE_SOLVE
      SUBROUTINE DMUMPS_SOLVE_NODE( INODE,
     &     BUFR, LBUFR, LBUFR_BYTES,
     &     MSGTAG, MSGSOU, MYID, SLAVEF, COMM,
     &     N, IPOOL, LPOOL, III, LEAF,
     &     NBFIN, NSTK_S,
     &     IWCB, LIWCB,
     &     WCB, LWCB, A, LA, IW, LIW,
     &     NRHS, POSWCB,
     &     PLEFTWCB, POSIWCB,
     &     PTRICB, PTRIST, PTRFAC, PROCNODE_STEPS,
     &     FILS, STEP, FRERE, DAD,
     &     MYROOT,
     &     INFO, KEEP,KEEP8, DKEEP, RHS_ROOT, MTYPE,
     &     RHSCOMP, LRHSCOMP, POSINRHSCOMP_FWD,
     &     
     &     ISTEP_TO_INIV2, TAB_POS_IN_PERE
     &    , RHS_BOUNDS, LRHS_BOUNDS, DO_NBSPARSE 
     &     , FROM_PP )
      USE DMUMPS_OOC
      USE DMUMPS_BUF
      IMPLICIT NONE
      INTEGER MTYPE
      INTEGER INODE, LBUFR, LBUFR_BYTES
      INTEGER MSGTAG, MSGSOU, MYID, SLAVEF, COMM
      INTEGER LIWCB, LIW, POSIWCB
      INTEGER(8) :: POSWCB, PLEFTWCB, LWCB
      INTEGER(8) :: LA
      INTEGER N, LPOOL, III, LEAF, NBFIN
      INTEGER MYROOT
      INTEGER INFO( 40 ), KEEP( 500)
      INTEGER(8) KEEP8(150)
      DOUBLE PRECISION, INTENT(INOUT) :: DKEEP(230)
      INTEGER BUFR( LBUFR )
      INTEGER IPOOL( LPOOL ), NSTK_S(KEEP(28))
      INTEGER IWCB( LIWCB ), IW( LIW )
      INTEGER NRHS
      DOUBLE PRECISION WCB( LWCB ), A( LA )
      DOUBLE PRECISION RHS_ROOT( * )
      INTEGER PTRICB(KEEP(28)), PTRIST(KEEP(28))
      INTEGER(8) :: PTRFAC(KEEP(28))
      INTEGER PROCNODE_STEPS(KEEP(28))
      INTEGER FILS( N ), STEP( N ), FRERE(KEEP(28)), DAD(KEEP(28))
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     &     TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
      INTEGER POSINRHSCOMP_FWD(N), LRHSCOMP
#if defined(RHSCOMP_BYROWS)
      DOUBLE PRECISION RHSCOMP(NRHS, LRHSCOMP)
#else
      DOUBLE PRECISION RHSCOMP(LRHSCOMP, NRHS)
#endif
      DOUBLE PRECISION VALPIV, A11, A22, A12, DETPIV
      LOGICAL, intent(in) :: DO_NBSPARSE
      INTEGER, intent(in) :: LRHS_BOUNDS
      INTEGER, intent(in) :: RHS_BOUNDS(LRHS_BOUNDS)
      LOGICAL, intent(in) :: FROM_PP
      EXTERNAL dgemv, dtrsv, dgemm, dtrsm, MUMPS_PROCNODE
      INTEGER MUMPS_PROCNODE
      DOUBLE PRECISION ALPHA,ONE,ZERO
      PARAMETER (ZERO=0.0D0, ONE = 1.0D0, ALPHA=-1.0D0)
      DOUBLE PRECISION TIME_TMP
      INTEGER JBDEB, JBFIN, NRHS_B 
      INTEGER(8) :: APOS, APOS1, APOS2, APOSOFF, IFR8, IFR_ini8
      INTEGER I, J, K, IPOS, NSLAVES, J1, J2, J3, FPERE, NPIV, NCB,
     &     IERR,
     &     LIELL, JJ,
     &     NELIM
      INTEGER(8) :: PCB_COURANT, PPIV_COURANT, PPIV_PANEL, PCB_PANEL
      INTEGER IPOSINRHSCOMP, IPOSINRHSCOMP_TMP
      INTEGER Effective_CB_Size, NUPDATE, ISLAVE, PDEST, FirstIndex
      LOGICAL FLAG
!$    LOGICAL :: OMP_FLAG
      INCLUDE 'mumps_headers.h'
      INTEGER(8) :: POSWCB1, POSWCB2
      INTEGER(8) :: APOSDEB
      INTEGER TempNROW, TempNCOL, PANEL_SIZE, LIWFAC, 
     &     JFIN, NBJ, NUPDATE_PANEL,
     &     NBK, NBK_ini, TYPEF
      INTEGER LD_WCBPIV         
      INTEGER LD_WCBCB          
      INTEGER LDAJ, LDAJ_ini, LDAJ_FIRST_PANEL
      INTEGER TMP_NBPANELS,
     &     I_PIVRPTR, I_PIVR, IPANEL
      LOGICAL MUST_BE_PERMUTED
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER DUMMY( 1 )
      DUMMY(1)=1
      IF (DO_NBSPARSE) THEN
       JBDEB= RHS_BOUNDS(2*STEP(INODE)-1)
       JBFIN= RHS_BOUNDS(2*STEP(INODE))
       NRHS_B = JBFIN-JBDEB+1
      ELSE
       JBDEB = 1
       JBFIN = NRHS
       NRHS_B = NRHS
      ENDIF
      IF (DO_NBSPARSE) THEN
       if (JBDEB.GT.JBFIN) then
         write(6,*) " Internal error 1 in nbsparse :", 
     &    JBDEB, JBFIN
         CALL MUMPS_ABORT()
       endif
       IF (JBDEB.LT.1 .OR. JBDEB.GT.NRHS .or. 
     &    JBFIN.LT.1 .OR. JBFIN.GT.NRHS ) THEN
         write(6,*) " Internal error 2 in nbsparse :", 
     &    JBDEB, JBFIN
         CALL MUMPS_ABORT()
       endif
      ENDIF
      IF ( INODE .eq. KEEP( 38 ) .OR. INODE .eq.KEEP( 20 ) ) THEN
         LIELL = IW( PTRIST( STEP(INODE)) + 3 + KEEP(IXSZ))
         NPIV  = LIELL
         NELIM = 0
         NSLAVES = 0
         IPOS = PTRIST( STEP(INODE)) + 5 + KEEP(IXSZ)
      ELSE
        IPOS = PTRIST(STEP(INODE)) + 2 + KEEP(IXSZ)
        LIELL = IW(IPOS-2)+IW(IPOS+1)
        NELIM = IW(IPOS-1)
        NSLAVES = IW( PTRIST(STEP(INODE)) + 5 + KEEP(IXSZ) )
        IPOS = IPOS + 1
        NPIV = IW(IPOS)
        IPOS = IPOS + 1
        IF (KEEP(201).GT.0) THEN
           CALL DMUMPS_SOLVE_GET_OOC_NODE(
     &          INODE,PTRFAC,KEEP,A,LA,STEP,
     &          KEEP8,N,MUST_BE_PERMUTED,IERR)
           IF(IERR.LT.0)THEN
              INFO(1)=IERR
              INFO(2)=0
              GOTO 260
           ENDIF
           IF (KEEP(201).EQ.1 .AND. KEEP(50).NE.1) THEN
           CALL DMUMPS_OOC_PP_CHECK_PERM_FREED(
     &                 IW(IPOS+1+2*LIELL+1+NSLAVES),
     &                 MUST_BE_PERMUTED )
           ENDIF
        ENDIF                     
        NSLAVES = IW( PTRIST(STEP(INODE)) + 5 + KEEP(IXSZ))
        IPOS = IPOS + 1 + NSLAVES
      END IF
      IF ( MTYPE .EQ. 1 .OR. KEEP(50) .NE. 0 ) THEN
         J1 = IPOS + 1
         J2 = IPOS + LIELL
         J3 = IPOS + NPIV
      ELSE
         J1 = IPOS + LIELL + 1
         J2 = IPOS + 2 * LIELL
         J3 = IPOS + LIELL + NPIV
      END IF
      NCB = LIELL-NPIV
      IF ( INODE .eq. KEEP( 38 ) .OR. INODE .eq. KEEP(20) ) THEN
        IFR8 = 0_8
        IPOSINRHSCOMP_TMP = POSINRHSCOMP_FWD(IW(J1)) 
        IF (KEEP(350).EQ.0) THEN
         DO JJ = J1, J3
            IFR8 = IFR8 + 1_8
            DO K=JBDEB,JBFIN
               RHS_ROOT(IFR8+int(NPIV,8)*int(K-1,8)) =
#if defined(RHSCOMP_BYROWS)
     &                      RHSCOMP(K,IPOSINRHSCOMP_TMP)
#else
     &                      RHSCOMP(IPOSINRHSCOMP_TMP,K) 
#endif
            END DO
            IPOSINRHSCOMP_TMP = IPOSINRHSCOMP_TMP + 1
         END DO
        ELSE IF (KEEP(350).EQ.1.OR.KEEP(350).EQ.2) THEN
          IFR_ini8 = IFR8
!$        OMP_FLAG = ( JBFIN-JBDEB+1.GE.KEEP(362) .AND.
!$   &                     (J3-J1+1)*(JBFIN-JBDEB+1) .GE. KEEP(363) )
!$OMP PARALLEL DO PRIVATE(IFR8,JJ) IF(OMP_FLAG)
          DO K=JBDEB,JBFIN
           IFR8 = IFR_ini8 + int(K-1,8)*int(NPIV,8)
           DO JJ = J1, J3
             RHS_ROOT(IFR8+int(JJ-J1+1,8)) = 
#if defined(RHSCOMP_BYROWS)
     &               RHSCOMP(K,IPOSINRHSCOMP_TMP+JJ-J1)
#else
     &               RHSCOMP(IPOSINRHSCOMP_TMP+JJ-J1,K)
#endif
           ENDDO
          ENDDO
!$OMP END PARALLEL DO
        ELSE
             WRITE(*,*) "UNKNOWN VERSION OF KEEP(350)"
             CALL MUMPS_ABORT()
        ENDIF
         IF ( NPIV .LT. LIELL ) THEN
            WRITE(*,*) ' Internal error in SOLVE_NODE for Root node'
            CALL MUMPS_ABORT()
         END IF
         MYROOT = MYROOT - 1
         IF ( MYROOT .EQ. 0 ) THEN
            NBFIN = NBFIN - 1
            IF (SLAVEF .GT. 1) THEN
               CALL DMUMPS_MCAST2(DUMMY, 1, MPI_INTEGER, MYID,
     &              COMM, RACINE_SOLVE, SLAVEF, KEEP)
            ENDIF
         END IF
         GO TO 270
      END IF
      APOS = PTRFAC(STEP(INODE))
      IF (KEEP(201).EQ.1) THEN  
        IF (MTYPE.EQ.1) THEN
            IF ((MTYPE.EQ.1).AND.NSLAVES.NE.0) THEN
              TempNROW= NPIV+NELIM
              TempNCOL= NPIV
              LDAJ_FIRST_PANEL=TempNROW
            ELSE
              TempNROW= LIELL
              TempNCOL= NPIV
              LDAJ_FIRST_PANEL=TempNROW
            ENDIF
            TYPEF=TYPEF_L
        ELSE 
            TempNCOL= LIELL
            TempNROW= NPIV
            LDAJ_FIRST_PANEL=TempNCOL
            TYPEF= TYPEF_U
        ENDIF
        LIWFAC =  IW(PTRIST(STEP(INODE))+XXI)
        PANEL_SIZE = DMUMPS_OOC_PANEL_SIZE( LDAJ_FIRST_PANEL )
      ENDIF                     
      PPIV_COURANT = PLEFTWCB
      PLEFTWCB = PLEFTWCB + int(LIELL,8) * int(NRHS_B,8)
      IF ( POSWCB - PLEFTWCB + 1_8 .LT. 0 ) THEN
         INFO(1) = -11
         CALL MUMPS_SET_IERROR(PLEFTWCB-POSWCB-1_8, INFO(2))
         GO TO 260
      END IF
      IF (KEEP(201).EQ.1) THEN  
         LD_WCBPIV = LIELL
         LD_WCBCB  = LIELL
         PCB_COURANT = PPIV_COURANT + NPIV
         DO K=JBDEB, JBFIN
            IFR8 = PPIV_COURANT+int(K-JBDEB,8)*int(LD_WCBPIV,8)-1_8
            IPOSINRHSCOMP = POSINRHSCOMP_FWD(IW(J1)) 
            DO JJ = J1, J3
               IFR8 = IFR8 + 1_8
#if defined(RHSCOMP_BYROWS)
               WCB(IFR8) = RHSCOMP(K,IPOSINRHSCOMP) 
#else
               WCB(IFR8) = RHSCOMP(IPOSINRHSCOMP,K) 
#endif
               IPOSINRHSCOMP = IPOSINRHSCOMP + 1
            ENDDO
            IF (NCB.GT.0) THEN
               DO JJ = J3+1, J2
                  J = IW(JJ)
                  IFR8 = IFR8 + 1_8
                  IPOSINRHSCOMP = abs(POSINRHSCOMP_FWD(J))
#if defined(RHSCOMP_BYROWS)
                  WCB(IFR8) = RHSCOMP(K,IPOSINRHSCOMP)
                  RHSCOMP (K,IPOSINRHSCOMP) = ZERO
#else
                  WCB(IFR8) = RHSCOMP(IPOSINRHSCOMP,K) 
                  RHSCOMP (IPOSINRHSCOMP,K) = ZERO
#endif
               ENDDO
            ENDIF
         ENDDO
      ELSE                      
         LD_WCBPIV = NPIV
         LD_WCBCB  = NCB
         PCB_COURANT = PPIV_COURANT + NPIV*NRHS_B
         IFR8 = PPIV_COURANT - 1_8
#if defined(ALLTIMERS)
         TIME_TMP=MPI_WTIME()
#endif
         IFR_ini8 = IFR8
         IPOSINRHSCOMP = POSINRHSCOMP_FWD(IW(J1)) 
         IF (KEEP(350).EQ.0) THEN
!$        OMP_FLAG = NRHS_B.GT.4 .AND. .FALSE.
!$OMP PARALLEL DO PRIVATE(J,IFR8,K) IF(OMP_FLAG)
          DO 130 JJ = J1, J3
            J = IW(JJ)
            IFR8 = IFR_ini8 + int(JJ-J1+1,8)
            DO K=JBDEB, JBFIN
#if defined(RHSCOMP_BYROWS)
               WCB(IFR8+(K-JBDEB)*NPIV) = 
     &                RHSCOMP(K,IPOSINRHSCOMP+JJ-J1)
#else
               WCB(IFR8+(K-JBDEB)*NPIV) = 
     &                RHSCOMP(IPOSINRHSCOMP+JJ-J1,K)
#endif
            END DO
 130      CONTINUE
!$OMP END PARALLEL DO
         ELSE IF (KEEP(350).EQ.1.OR.KEEP(350).EQ.2) THEN
!$        OMP_FLAG = ( JBFIN-JBDEB+1.GE.KEEP(362) .AND.
!$   &                     (J3-J1+1)*(JBFIN-JBDEB+1) .GE. KEEP(363) )
!$OMP PARALLEL DO PRIVATE(JJ,IFR8) IF(OMP_FLAG)
          DO K=JBDEB, JBFIN
           IFR8 = IFR_ini8 + (K-JBDEB)*NPIV
           DO  JJ = J1, J3
#if defined(RHSCOMP_BYROWS)
            WCB(IFR8+int(JJ-J1+1,8)) = 
     &           RHSCOMP(K,IPOSINRHSCOMP+JJ-J1)
#else
            WCB(IFR8+int(JJ-J1+1,8)) = 
     &           RHSCOMP(IPOSINRHSCOMP+JJ-J1,K)
#endif
           ENDDO
          ENDDO
!$OMP END PARALLEL DO
         ELSE
             WRITE(*,*) "UNKNOWN VERSION OF KEEP(350)"
             CALL MUMPS_ABORT()
         ENDIF
#if defined(ALLTIMERS)
         DKEEP(123)=DKEEP(123)+MPI_WTIME()-TIME_TMP
#endif
         IFR8 = PCB_COURANT - 1_8
         IF (NPIV .LT. LIELL) THEN
            IFR_ini8 = IFR8
#if defined(ALLTIMERS)
            TIME_TMP=MPI_WTIME()
#endif
            IF (KEEP(350).EQ.0) THEN
!$OMP PARALLEL DO PRIVATE(J,IFR8,K,IPOSINRHSCOMP) IF(OMP_FLAG)
            DO 140 JJ = J3 + 1, J2
               J = IW(JJ)
               IFR8 = IFR_ini8 + (JJ-J3)
               IPOSINRHSCOMP = abs(POSINRHSCOMP_FWD(J))
               DO K=JBDEB, JBFIN
#if defined(RHSCOMP_BYROWS)
                  WCB(IFR8+(K-JBDEB)*NCB) = RHSCOMP(K,IPOSINRHSCOMP)
#else
                  WCB(IFR8+(K-JBDEB)*NCB) = RHSCOMP(IPOSINRHSCOMP,K)
#endif
#if defined(RHSCOMP_BYROWS)
                  RHSCOMP(K,IPOSINRHSCOMP)=ZERO
#else
                  RHSCOMP(IPOSINRHSCOMP,K)=ZERO
#endif
               ENDDO
 140        CONTINUE
!$OMP END PARALLEL DO
           ELSE IF (KEEP(350).EQ.1.OR.KEEP(350).EQ.2) THEN
!$            OMP_FLAG = ( JBFIN-JBDEB+1.GE.KEEP(362) .AND.
!$   &                     (J2-J3)*(JBFIN-JBDEB+1) .GE. KEEP(363) )
!$OMP PARALLEL DO PRIVATE (IFR8, JJ, J, IPOSINRHSCOMP) IF (OMP_FLAG)
              DO K=JBDEB, JBFIN
                IFR8 = IFR_ini8+(K-JBDEB)*NCB
                DO JJ = J3 + 1, J2
                  J = IW(JJ)
                  IPOSINRHSCOMP = abs(POSINRHSCOMP_FWD(J))
#if defined(RHSCOMP_BYROWS)
                  WCB(IFR8+int(JJ-J3,8)) = RHSCOMP(K,IPOSINRHSCOMP)
#else
                  WCB(IFR8+int(JJ-J3,8)) = RHSCOMP(IPOSINRHSCOMP,K)
#endif
#if defined(RHSCOMP_BYROWS)
                  RHSCOMP(K,IPOSINRHSCOMP)=ZERO
#else
                  RHSCOMP(IPOSINRHSCOMP,K)=ZERO
#endif
                ENDDO
              ENDDO
!$OMP END PARALLEL DO
           ELSE
             WRITE(*,*) "UNKNOWN VERSION OF KEEP(350)"
             CALL MUMPS_ABORT()
           ENDIF
#if defined(ALLTIMERS)
         DKEEP(124)=DKEEP(124)+MPI_WTIME()-TIME_TMP
#endif
         ENDIF
      ENDIF                     
#if defined(ALLTIMERS)
      TIME_TMP=MPI_WTIME()  
#endif
      IF ( NPIV .NE. 0 ) THEN
         IF (KEEP(201).EQ.1) THEN 
        APOSDEB = APOS
        J = 1
        IPANEL = 0
  10    CONTINUE
          IPANEL = IPANEL + 1
          JFIN    = min(J+PANEL_SIZE-1, NPIV)
          IF (IW(IPOS+ LIELL + JFIN) < 0) THEN
            JFIN=JFIN+1
          ENDIF
          NBJ     = JFIN-J+1
          LDAJ    = LDAJ_FIRST_PANEL-J+1 
          IF ( (KEEP(50).NE.1).AND. MUST_BE_PERMUTED ) THEN
           CALL DMUMPS_GET_OOC_PERM_PTR(TYPEF, TMP_NBPANELS,
     &            I_PIVRPTR, I_PIVR, IPOS+1+2*LIELL, IW, LIW)
               IF (NPIV.EQ.(IW(I_PIVRPTR+IPANEL-1)-1)) THEN
                  MUST_BE_PERMUTED=.FALSE. 
               ELSE
                  CALL DMUMPS_PERMUTE_PANEL(
     &                 IW( I_PIVR+ IW(I_PIVRPTR+IPANEL-1)-
     &                 IW(I_PIVRPTR)),
     &                 NPIV-IW(I_PIVRPTR+IPANEL-1)+1, 
     &                 IW(I_PIVRPTR+IPANEL-1)-1,
     &                 A(APOSDEB),
     &                 LDAJ, NBJ, J-1 ) 
               ENDIF
            ENDIF 
            NUPDATE_PANEL = LDAJ - NBJ
            PPIV_PANEL = PPIV_COURANT+int(J-1,8)
            PCB_PANEL  = PPIV_PANEL+int(NBJ,8)
            APOS1 = APOSDEB+int(NBJ,8)
            IF  (MTYPE.EQ.1) THEN
#if defined(MUMPS_USE_BLAS2)
               IF ( NRHS_B == 1 ) THEN
                  CALL dtrsv( 'L', 'N', 'U', NBJ, A(APOSDEB), LDAJ, 
     &                 WCB(PPIV_PANEL), 1 )
                  IF (NUPDATE_PANEL.GT.0) THEN
                     CALL dgemv('N', NUPDATE_PANEL,NBJ,ALPHA, A(APOS1),
     &                    LDAJ,  WCB(PPIV_PANEL), 1, ONE,
     &                    WCB(PCB_PANEL), 1)
                  ENDIF
               ELSE
#endif
                  CALL dtrsm( 'L','L','N','U', NBJ, NRHS_B, ONE,
     &                 A(APOSDEB), LDAJ, WCB(PPIV_PANEL),
     &                 LIELL )
                  IF (NUPDATE_PANEL.GT.0) THEN
                     CALL dgemm('N', 'N', NUPDATE_PANEL, NRHS_B, NBJ, 
     &                    ALPHA,
     &                    A(APOS1), LDAJ, WCB(PPIV_PANEL), LIELL, ONE,
     &                    WCB(PCB_PANEL), LIELL)
                  ENDIF
#if defined(MUMPS_USE_BLAS2)
               ENDIF
#endif
            ELSE
#if defined(MUMPS_USE_BLAS2)
               IF (NRHS_B == 1) THEN
                  CALL dtrsv( 'L', 'N', 'N', NBJ, A(APOSDEB), LDAJ,
     &                 WCB(PPIV_PANEL), 1 )
                  IF (NUPDATE_PANEL.GT.0) THEN
                     CALL dgemv('N',NUPDATE_PANEL, NBJ, ALPHA, A(APOS1),
     &                    LDAJ, WCB(PPIV_PANEL), 1,
     &                    ONE, WCB(PCB_PANEL), 1 )
                  ENDIF
               ELSE
#endif
                  CALL dtrsm('L','L','N','N',NBJ, NRHS_B, ONE,
     &                 A(APOSDEB), LDAJ, WCB(PPIV_PANEL),
     &                 LIELL)
                  IF (NUPDATE_PANEL.GT.0) THEN
                     CALL dgemm('N', 'N', NUPDATE_PANEL, NRHS_B, NBJ, 
     &                    ALPHA,
     &                    A(APOS1), LDAJ, WCB(PPIV_PANEL), LIELL, ONE,
     &             WCB(PCB_PANEL), LIELL)
                  ENDIF
#if defined(MUMPS_USE_BLAS2)
               ENDIF
#endif
            ENDIF
            APOSDEB = APOSDEB+int(LDAJ,8)*int(NBJ,8)
            J=JFIN+1
            IF ( J .LE. NPIV ) GOTO 10
         ELSE                   
            IF (KEEP(50).NE.0) THEN
#if defined(MUMPS_USE_BLAS2)
               IF ( NRHS_B == 1 ) THEN
                  CALL dtrsv( 'U', 'T', 'U', NPIV, A(APOS), NPIV,
     &                   WCB(PPIV_COURANT), 1 )
               ELSE
#endif
                  CALL dtrsm( 'L','U','T','U', NPIV, NRHS_B, ONE,
     &                   A(APOS), NPIV, WCB(PPIV_COURANT),
     &                   NPIV )
#if defined(MUMPS_USE_BLAS2)
               ENDIF
#endif
            ELSE
               IF ( MTYPE .eq. 1 ) THEN
#if defined(MUMPS_USE_BLAS2)
                  IF ( NRHS_B == 1)  THEN
                     CALL dtrsv( 'U', 'T', 'U', NPIV, A(APOS), LIELL, 
     &                    WCB(PPIV_COURANT), 1 )
                  ELSE
#endif
                     CALL dtrsm( 'L','U','T','U', NPIV, NRHS_B, ONE,
     &                    A(APOS), LIELL, WCB(PPIV_COURANT),
     &                    NPIV )
#if defined(MUMPS_USE_BLAS2)
                  ENDIF
#endif
               ELSE
#if defined(MUMPS_USE_BLAS2)
                  IF (NRHS_B == 1) THEN
                     CALL dtrsv( 'L', 'N', 'N', NPIV, A(APOS), LIELL,
     &                    WCB(PPIV_COURANT), 1 )
                  ELSE
#endif
                     CALL dtrsm('L','L','N','N',NPIV, NRHS_B, ONE,
     &                    A(APOS), LIELL, WCB(PPIV_COURANT),
     &                    NPIV)
#if defined(MUMPS_USE_BLAS2)
                  ENDIF
#endif
               END IF
            END IF              
         END IF                 
      END IF                    
      NCB   = LIELL - NPIV
      IF ( MTYPE .EQ. 1 ) THEN
         IF ( KEEP(50) .eq. 0 ) THEN
            APOS1 = APOS  + int(NPIV,8) * int(LIELL,8)
         ELSE
            APOS1 = APOS + int(NPIV,8) * int(NPIV,8)
         END IF
         IF ( NSLAVES .EQ. 0 .OR. NPIV .eq. 0 ) THEN
            NUPDATE = NCB
         ELSE
            NUPDATE = NELIM
         END IF
      ELSE
         APOS1 = APOS + int(NPIV,8)
         NUPDATE = NCB
      END IF
      IF (KEEP(201).NE.1) THEN  
         IF ( NPIV .NE. 0 .AND. NUPDATE.NE.0 ) THEN
            IF ( MTYPE .eq. 1 ) THEN
#if defined(MUMPS_USE_BLAS2)
               IF ( NRHS_B == 1 ) THEN
                  CALL dgemv('T', NPIV, NUPDATE, ALPHA, A(APOS1),
     &            NPIV,  WCB(PPIV_COURANT), 1, ONE,
     &            WCB(PCB_COURANT), 1)
               ELSE
#endif
                  CALL dgemm('T', 'N', NUPDATE, NRHS_B, NPIV, ALPHA,
     &            A(APOS1), NPIV, WCB(PPIV_COURANT), NPIV, ONE,
     &            WCB(PCB_COURANT), NCB)
#if defined(MUMPS_USE_BLAS2)
               END IF
#endif
            ELSE                
#if defined(MUMPS_USE_BLAS2)
               IF ( NRHS_B == 1 ) THEN
                  CALL dgemv('N',NUPDATE, NPIV, ALPHA, A(APOS1),
     &                 LIELL, WCB(PPIV_COURANT), 1,
     &                 ONE, WCB(PCB_COURANT), 1 )
               ELSE
#endif
                  CALL dgemm('N', 'N', NUPDATE, NRHS_B, NPIV, ALPHA,
     &                 A(APOS1), LIELL, WCB(PPIV_COURANT), NPIV, ONE,
     &                 WCB(PCB_COURANT), NCB)
#if defined(MUMPS_USE_BLAS2)
               END IF
#endif
            END IF
         END IF
      END IF                    
#if defined(ALLTIMERS)
      DKEEP(125)=DKEEP(125)+MPI_WTIME()-TIME_TMP
      TIME_TMP=MPI_WTIME()
#endif
      IPOSINRHSCOMP =  POSINRHSCOMP_FWD(IW(J1)) 
      IF ( KEEP(50) .eq. 0 ) THEN
         IF (KEEP(350).EQ.0) THEN
           DO K=JBDEB,JBFIN
             IFR8 =  PPIV_COURANT + int(K-JBDEB,8)*int(LD_WCBPIV,8)
#if defined(RHSCOMP_BYROWS)
             RHSCOMP(K,IPOSINRHSCOMP:IPOSINRHSCOMP+NPIV-1) =
#else
             RHSCOMP(IPOSINRHSCOMP:IPOSINRHSCOMP+NPIV-1,K) =
#endif
     &            WCB(IFR8:IFR8+int(NPIV-1,8))
           ENDDO
         ELSE IF (KEEP(350).EQ.1.OR.KEEP(350).EQ.2) THEN
#if defined(RHSCOMP_BYROWS)
!$         OMP_FLAG = (NPIV.GE.KEEP(362).AND.NRHS_B*NPIV.GE.KEEP(363))
!$OMP      PARALLEL DO PRIVATE(IFR8,K) IF (OMP_FLAG)
           DO I=1,NPIV
             IFR8 =  PPIV_COURANT + I-1
             DO K=JBDEB,JBFIN
               RHSCOMP(K,IPOSINRHSCOMP+I-1) =
     &              WCB(IFR8+(K-JBDEB)*LD_WCBPIV)
             ENDDO
           ENDDO
!$OMP      END PARALLEL DO
#else
!$         OMP_FLAG=(NRHS_B.GE.KEEP(362).AND.NRHS_B*NPIV.GE.KEEP(363))
!$OMP PARALLEL DO PRIVATE(IFR8) IF (OMP_FLAG)
           DO K=JBDEB,JBFIN
             IFR8 =  PPIV_COURANT + (K-JBDEB)*LD_WCBPIV
             RHSCOMP(IPOSINRHSCOMP:IPOSINRHSCOMP+NPIV-1, K) =
     &            WCB(IFR8:IFR8+int(NPIV-1,8))
           ENDDO
!$OMP END PARALLEL DO
#endif
         ELSE
           WRITE(*,*) "UNKNOWN VERSION OF KEEP(350)"
           CALL MUMPS_ABORT()
         ENDIF
      ELSE
         IFR8 = PPIV_COURANT - 1_8
         IF (KEEP(201).EQ.1) THEN 
            LDAJ = TempNROW  
         ELSE                
            LDAJ = NPIV 
         ENDIF
         APOS1 = APOS
         JJ    = J1
         IF (KEEP(201).EQ.1) THEN
            NBK   = 0           
         ENDIF
        IF (KEEP(350).EQ.0) THEN
         DO 
            IF(JJ .GT. J3) EXIT
            IFR8 = IFR8 + 1_8
            IF(IW(JJ+LIELL) .GT. 0) THEN
              VALPIV  = ONE/A( APOS1 )
              DO K=JBDEB, JBFIN
#if defined(RHSCOMP_BYROWS)
                  RHSCOMP(K, IPOSINRHSCOMP+JJ-J1) = 
     &                 WCB( IFR8+(K-JBDEB)*LD_WCBPIV ) * VALPIV
#else
                  RHSCOMP(IPOSINRHSCOMP+JJ-J1 , K ) = 
     &                 WCB( IFR8+(K-JBDEB)*LD_WCBPIV ) * VALPIV
#endif
              END DO
              IF (KEEP(201).EQ.1) THEN
                NBK = NBK+1
                IF (NBK.EQ.PANEL_SIZE) THEN
                  NBK = 0
                  LDAJ = LDAJ - PANEL_SIZE
                ENDIF
              ENDIF
              APOS1 = APOS1 + int(LDAJ + 1,8)
              JJ = JJ+1
            ELSE
              IF (KEEP(201).EQ.1) THEN
                NBK = NBK+1
              ENDIF
              APOS2 = APOS1+int(LDAJ+1,8)
              IF (KEEP(201).EQ.1) THEN
                APOSOFF = APOS1+int(LDAJ,8)
              ELSE
                APOSOFF=APOS1+1_8
              ENDIF
              A11 = A(APOS1)
              A22 = A(APOS2)
              A12 = A(APOSOFF)
              DETPIV = A11*A22 - A12**2
              A22 = A11/DETPIV
              A11 = A(APOS2)/DETPIV
              A12 = -A12/DETPIV
              DO K=JBDEB, JBFIN
                  POSWCB1 = IFR8+int(K-JBDEB,8)*int(LD_WCBPIV,8)
                  POSWCB2 = POSWCB1+1_8
#if defined(RHSCOMP_BYROWS)
                  RHSCOMP(K,IPOSINRHSCOMP+JJ-J1) =
     &               WCB(POSWCB1)*A11
     &               + WCB(POSWCB2)*A12
                  RHSCOMP(K,IPOSINRHSCOMP+JJ-J1+1) = 
     &                 WCB(POSWCB1)*A12
     &                 + WCB(POSWCB2)*A22
#else
                  RHSCOMP(IPOSINRHSCOMP+JJ-J1,K) =
     &               WCB(POSWCB1)*A11
     &               + WCB(POSWCB2)*A12
                  RHSCOMP(IPOSINRHSCOMP+JJ-J1+1,K) = 
     &                 WCB(POSWCB1)*A12
     &                 + WCB(POSWCB2)*A22
#endif
              END DO
              IF (KEEP(201).EQ.1) THEN
                 NBK = NBK+1
                 IF (NBK.GE.PANEL_SIZE) THEN
                    LDAJ = LDAJ - NBK
                    NBK = 0
                 ENDIF
              ENDIF
              APOS1 = APOS2 + int(LDAJ + 1,8)
              JJ = JJ+2
              IFR8 = IFR8+1_8
            ENDIF
         ENDDO
        ELSE IF (KEEP(350).EQ.1.OR.KEEP(350).EQ.2) THEN
         IFR_ini8 = PPIV_COURANT - 1_8
         LDAJ_ini = LDAJ   
         IF (KEEP(201).EQ.1) NBK_ini = NBK  
!$           OMP_FLAG = ( JBFIN-JBDEB+1.GE.KEEP(362) .AND.
!$   &                  ((J3-J1+1)*(JBFIN-JBDEB+1) .GE. KEEP(363)))
!$OMP  PARALLEL DO PRIVATE(JJ,IFR8,NBK,APOS1,APOS2,APOSOFF,VALPIV,
!$OMP&      POSWCB1, POSWCB2,A11,A22,A12,DETPIV,LDAJ)  IF(OMP_FLAG)
         DO K = JBDEB, JBFIN
            IFR8  = IFR_ini8 + int(K-JBDEB,8)*int(LD_WCBPIV,8)
            NBK   = NBK_ini
            APOS1 = APOS
            LDAJ  = LDAJ_ini
            JJ    = J1
          DO 
            IF (JJ .GT. J3) EXIT
            IFR8  = IFR8 + 1_8
            IF (IW(JJ+LIELL) .GT. 0) THEN
               VALPIV  = ONE/A( APOS1 )
#if defined(RHSCOMP_BYROWS)
                RHSCOMP(K, IPOSINRHSCOMP+JJ-J1) = 
     &                 WCB( IFR8 ) * VALPIV
#else
                 RHSCOMP(IPOSINRHSCOMP+JJ-J1 , K ) = 
     &                 WCB( IFR8 ) * VALPIV
#endif
              IF (KEEP(201).EQ.1) THEN
                NBK = NBK+1
                IF (NBK.EQ.PANEL_SIZE) THEN
                  NBK = 0
                  LDAJ = LDAJ - PANEL_SIZE
                ENDIF
              ENDIF
              APOS1 = APOS1 + int(LDAJ + 1,8)
              JJ = JJ+1
            ELSE
              IF (KEEP(201).EQ.1) THEN
                NBK = NBK+1
              ENDIF
              APOS2 = APOS1+int(LDAJ+1,8)
              IF (KEEP(201).EQ.1) THEN
                APOSOFF = APOS1+int(LDAJ,8)
              ELSE
                APOSOFF=APOS1+1_8
              ENDIF
              A11 = A(APOS1)
              A22 = A(APOS2)
              A12 = A(APOSOFF)
              DETPIV = A11*A22 - A12**2
              A22 = A11/DETPIV
              A11 = A(APOS2)/DETPIV
              A12 = -A12/DETPIV
              POSWCB1 = IFR8
              POSWCB2 = POSWCB1+1_8
#if defined(RHSCOMP_BYROWS)
              RHSCOMP(K,IPOSINRHSCOMP+JJ-J1) =
     &               WCB(POSWCB1)*A11
     &               + WCB(POSWCB2)*A12
              RHSCOMP(K,IPOSINRHSCOMP+JJ-J1+1) = 
     &                 WCB(POSWCB1)*A12
     &                 + WCB(POSWCB2)*A22
#else
              RHSCOMP(IPOSINRHSCOMP+JJ-J1,K) =
     &               WCB(POSWCB1)*A11
     &               + WCB(POSWCB2)*A12
              RHSCOMP(IPOSINRHSCOMP+JJ-J1+1,K) = 
     &                 WCB(POSWCB1)*A12
     &                 + WCB(POSWCB2)*A22
#endif
              IF (KEEP(201).EQ.1) THEN
                 NBK = NBK+1
                 IF (NBK.GE.PANEL_SIZE) THEN
                    LDAJ = LDAJ - NBK
                    NBK = 0
                 ENDIF
              ENDIF
              APOS1 = APOS2 + int(LDAJ + 1,8)
              JJ = JJ+2
              IFR8 = IFR8+1_8
            ENDIF  
           ENDDO   
         ENDDO     
!$OMP END PARALLEL DO
         ELSE
             WRITE(*,*) "UNKNOWN VERSION OF KEEP(350)"
             CALL MUMPS_ABORT()
         ENDIF
      END IF
#if defined(ALLTIMERS)
      DKEEP(126)=DKEEP(126)+MPI_WTIME()-TIME_TMP
#endif
      IF (KEEP(201).GT.0) THEN
         CALL DMUMPS_FREE_FACTORS_FOR_SOLVE(INODE,PTRFAC,KEEP(28),
     &        A,LA,.TRUE.,IERR)
         IF(IERR.LT.0)THEN
            INFO(1)=IERR
            INFO(2)=0
            GOTO 260
         ENDIF
      END IF
      FPERE = DAD(STEP(INODE))
      IF ( FPERE .EQ. 0 ) THEN
         MYROOT = MYROOT - 1
         PLEFTWCB = PLEFTWCB - LIELL *NRHS_B
         IF ( MYROOT .EQ. 0 ) THEN
            NBFIN = NBFIN - 1
            IF (SLAVEF .GT. 1) THEN
               CALL DMUMPS_MCAST2(DUMMY, 1, MPI_INTEGER, MYID,
     &             COMM, RACINE_SOLVE, SLAVEF, KEEP)
            ENDIF
         END IF
         GO TO 270
      ENDIF
      IF ( NUPDATE .NE. 0 .OR. NCB.eq.0 ) THEN
         IF (MUMPS_PROCNODE(PROCNODE_STEPS(STEP(FPERE)),
     &        SLAVEF) .EQ. MYID) THEN
            IF ( NCB .ne. 0 ) THEN
             PTRICB(STEP(INODE)) = NCB + 1
             IF (KEEP(350).EQ.0) THEN
!$             OMP_FLAG = .FALSE. 
#if defined(ALLTIMERS)
               TIME_TMP=MPI_WTIME()
#endif
!$OMP PARALLEL DO PRIVATE(K,IPOSINRHSCOMP_TMP) IF(OMP_FLAG)
               DO 190 I = 1, NUPDATE
                  IPOSINRHSCOMP_TMP = 
     &                 abs(POSINRHSCOMP_FWD(IW(J3 + I)))
                  DO K=JBDEB, JBFIN
#if defined(RHSCOMP_BYROWS)
                     RHSCOMP( K, IPOSINRHSCOMP_TMP ) = 
     &                    RHSCOMP( K, IPOSINRHSCOMP_TMP )
#else
                     RHSCOMP( IPOSINRHSCOMP_TMP, K ) = 
     &                    RHSCOMP( IPOSINRHSCOMP_TMP, K )
#endif
     &             + WCB(PCB_COURANT + I-1 +(K-JBDEB)*LD_WCBCB)
                  ENDDO
 190           CONTINUE
!$OMP END PARALLEL DO
             ELSE IF (KEEP(350).EQ.1.OR.KEEP(350).EQ.2) THEN
!$              OMP_FLAG = ( JBFIN-JBDEB+1.GE.KEEP(362) .AND.
!$   &                     (NUPDATE*(JBFIN-JBDEB+1) .GE. KEEP(363)) )
!$OMP PARALLEL DO PRIVATE(I,IFR8,IPOSINRHSCOMP_TMP) IF(OMP_FLAG)
                DO K = JBDEB, JBFIN
                  IFR8 = PCB_COURANT + int(K-JBDEB,8)*int(LD_WCBCB,8)
                  DO I = 1, NUPDATE
                    IPOSINRHSCOMP_TMP = 
     &                 abs(POSINRHSCOMP_FWD(IW(J3 + I)))
#if defined(RHSCOMP_BYROWS)
                     RHSCOMP( K, IPOSINRHSCOMP_TMP ) = 
     &                    RHSCOMP( K, IPOSINRHSCOMP_TMP )
#else
                     RHSCOMP( IPOSINRHSCOMP_TMP, K ) = 
     &                    RHSCOMP( IPOSINRHSCOMP_TMP, K )
#endif
     &             + WCB(IFR8 + int(I-1,8))
                  ENDDO
                ENDDO
!$OMP END PARALLEL DO
         ELSE
             WRITE(*,*) "UNKNOWN VERSION OF KEEP(350)"
             CALL MUMPS_ABORT()
         ENDIF
#if defined(ALLTIMERS)
               DKEEP(127)=DKEEP(127)+MPI_WTIME()-TIME_TMP
#endif
               PTRICB(STEP( INODE )) = PTRICB(STEP( INODE )) - NUPDATE
               IF ( PTRICB(STEP(INODE)) == 1 ) THEN
                  NSTK_S(STEP(FPERE)) = NSTK_S(STEP(FPERE)) - 1
                  IF (NSTK_S(STEP(FPERE)) .EQ. 0) THEN
                        IPOOL( LEAF ) = FPERE
                     LEAF = LEAF + 1
                  ENDIF
               END IF
            ELSE
               PTRICB(STEP( INODE )) = -1
               NSTK_S(STEP(FPERE)) = NSTK_S(STEP(FPERE)) - 1
               IF (NSTK_S(STEP(FPERE)) .EQ. 0) THEN 
                     IPOOL( LEAF ) = FPERE 
                  LEAF = LEAF + 1
               ENDIF            
            ENDIF
         ELSE
 210        CONTINUE
#if defined(ALLTIMERS)
            TIME_TMP=MPI_WTIME()
#endif
            CALL DMUMPS_BUF_SEND_VCB( NRHS_B, INODE, FPERE, 
     &           NCB, LD_WCBCB,
     &           NUPDATE,
     &           IW( J3 + 1 ), WCB( PCB_COURANT ), JBDEB, JBFIN,
     &           RHSCOMP, 1, 1, -9999, -9999,
     &           KEEP,
     &           MUMPS_PROCNODE(PROCNODE_STEPS(STEP(FPERE)), SLAVEF),
     &           ContVec,
     &           COMM, IERR )
#if defined(ALLTIMERS)
            IF (IERR.EQ.0) THEN
              DKEEP(132)=DKEEP(132)+MPI_WTIME()-TIME_TMP
            ELSE
              DKEEP(134)=DKEEP(134)+MPI_WTIME()-TIME_TMP
            ENDIF
#endif
            IF ( IERR .EQ. -1 ) THEN
               CALL DMUMPS_SOLVE_RECV_AND_TREAT( .FALSE., FLAG,
     &              BUFR, LBUFR, LBUFR_BYTES,
     &              MYID, SLAVEF, COMM,
     &              N, NRHS, IPOOL, LPOOL, III, LEAF,
     &              NBFIN, NSTK_S, IW, LIW, A, LA, PTRIST, PTRFAC,
     &              IWCB, LIWCB,
     &              WCB, LWCB, POSWCB, PLEFTWCB, POSIWCB,
     &              PTRICB, INFO, KEEP,KEEP8, DKEEP, STEP,
     &              PROCNODE_STEPS,
     &              RHSCOMP, LRHSCOMP, POSINRHSCOMP_FWD
     &              , FROM_PP )
               IF ( INFO( 1 )  .LT. 0 )  GOTO 270
               GOTO 210
            ELSE IF ( IERR .EQ. -2 ) THEN
               INFO( 1 ) = -17
               INFO( 2 ) = NUPDATE * KEEP( 35 ) +
     &              ( NUPDATE + 3 ) * KEEP( 34 )
               GOTO 260
            ELSE IF ( IERR .EQ. -3 ) THEN
               INFO( 1 ) = -20
               INFO( 2 ) = NUPDATE * KEEP( 35 ) +
     &              ( NUPDATE + 3 ) * KEEP( 34 )
               GOTO 260
            END IF
         ENDIF
      END IF
      IF ( NSLAVES .NE. 0 .AND. MTYPE .eq. 1
     &     .and. NPIV .NE. 0 ) THEN
         DO ISLAVE = 1, NSLAVES
            PDEST = IW( PTRIST(STEP(INODE)) + 5 + ISLAVE +KEEP(IXSZ))
            CALL MUMPS_BLOC2_GET_SLAVE_INFO( 
     &           KEEP,KEEP8, INODE, STEP, N, SLAVEF,
     &           ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &           ISLAVE, NCB - NELIM, 
     &           NSLAVES, 
     &           Effective_CB_Size, FirstIndex )
 222        CONTINUE
#if defined(ALLTIMERS)
                 TIME_TMP=MPI_WTIME()
#endif
                 CALL DMUMPS_BUF_SEND_MASTER2SLAVE( NRHS_B,
     &           INODE, FPERE,
     &           Effective_CB_Size, LD_WCBCB, LD_WCBPIV, NPIV,
     &           JBDEB, JBFIN, 
     &           WCB( PCB_COURANT + NELIM + FirstIndex - 1 ),
     &           WCB( PPIV_COURANT ),
     &           PDEST, COMM, KEEP, IERR )
#if defined(ALLTIMERS)
                 IF (IERR.EQ.0) THEN
                   DKEEP(133)=DKEEP(133)+MPI_WTIME()-TIME_TMP
                 ELSE
                   DKEEP(134)=DKEEP(134)+MPI_WTIME()-TIME_TMP
                 ENDIF
#endif
            IF ( IERR .EQ. -1 ) THEN
               CALL DMUMPS_SOLVE_RECV_AND_TREAT( .FALSE., FLAG,
     &              BUFR, LBUFR, LBUFR_BYTES,
     &              MYID, SLAVEF, COMM,
     &              N, NRHS, IPOOL, LPOOL, III, LEAF,
     &              NBFIN, NSTK_S, IW, LIW, A, LA, PTRIST,PTRFAC,
     &              IWCB, LIWCB,
     &              WCB, LWCB, POSWCB, PLEFTWCB, POSIWCB,
     &              PTRICB, INFO, KEEP,KEEP8, DKEEP, STEP,
     &              PROCNODE_STEPS,
     &              RHSCOMP, LRHSCOMP, POSINRHSCOMP_FWD
     &              , FROM_PP )
               IF ( INFO( 1 )  .LT. 0 )  GOTO 270
               GOTO 222
            ELSE IF ( IERR .EQ. -2 ) THEN
               INFO( 1 ) = -17
               INFO( 2 ) = (NPIV+Effective_CB_Size)*NRHS_B*KEEP(35) +
     &               6 * KEEP( 34 )
               GOTO 260
            ELSE IF ( IERR .EQ. -3 ) THEN
               INFO( 1 ) = -20
               INFO( 2 ) = (NPIV+Effective_CB_Size)*NRHS_B*KEEP(35) +
     &              6 * KEEP( 34 )
               GOTO 260
            END IF
         END DO
      END IF
      PLEFTWCB = PLEFTWCB - LIELL*NRHS_B
 270  CONTINUE
      RETURN
 260  CONTINUE
      CALL DMUMPS_BDC_ERROR( MYID, SLAVEF, COMM, KEEP )
      RETURN
      END SUBROUTINE DMUMPS_SOLVE_NODE
      RECURSIVE SUBROUTINE DMUMPS_SOLVE_RECV_AND_TREAT( BLOQ, FLAG,
     &           BUFR, LBUFR, LBUFR_BYTES,
     &           MYID, SLAVEF, COMM,
     &           N, NRHS, IPOOL, LPOOL, III, LEAF,
     &           NBFIN, NSTK_S, IW, LIW, A, LA, PTRIST,PTRFAC,
     &           IWCB, LIWCB,
     &           WCB, LWCB, POSWCB,
     &           PLEFTWCB, POSIWCB,
     &           PTRICB, INFO, KEEP,KEEP8, DKEEP, STEP, PROCNODE_STEPS,
     &           RHSCOMP, LRHSCOMP, POSINRHSCOMP_FWD
     &            , FROM_PP)
      IMPLICIT NONE
      LOGICAL BLOQ
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER MYID, SLAVEF, COMM
      INTEGER N, NRHS, LPOOL, III, LEAF, NBFIN
      INTEGER LIWCB, POSIWCB
      INTEGER(8) :: POSWCB, PLEFTWCB
      INTEGER LIW
      INTEGER(8), INTENT(IN) :: LA, LWCB
      INTEGER INFO( 40 ), KEEP( 500)
      INTEGER(8) KEEP8(150)
      DOUBLE PRECISION, INTENT(INOUT) :: DKEEP(230)
      INTEGER BUFR( LBUFR ), IPOOL(LPOOL)
      INTEGER NSTK_S( KEEP(28) )
      INTEGER IWCB( LIWCB )
      INTEGER IW( LIW )
      DOUBLE PRECISION WCB( LWCB ), A( LA )
      INTEGER PTRICB(KEEP(28)), PTRIST(KEEP(28))
      INTEGER(8) :: PTRFAC(KEEP(28))
      INTEGER STEP(N)
      INTEGER PROCNODE_STEPS(KEEP(28))
      LOGICAL FLAG
      INTEGER LRHSCOMP, POSINRHSCOMP_FWD(N)
#if defined(RHSCOMP_BYROWS)
      DOUBLE PRECISION RHSCOMP(NRHS,LRHSCOMP)
#else
      DOUBLE PRECISION RHSCOMP(LRHSCOMP,NRHS)
#endif
      LOGICAL, intent(in) :: FROM_PP
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER :: IERR
      INTEGER :: STATUS(MPI_STATUS_SIZE)
      INTEGER MSGSOU, MSGTAG, MSGLEN
      DOUBLE PRECISION :: TIME_TMP
#if defined(ALLTIMERS)
      TIME_TMP=MPI_WTIME()
#endif
      FLAG = .FALSE.
      IF ( BLOQ ) THEN
        CALL MPI_PROBE( MPI_ANY_SOURCE, MPI_ANY_TAG,
     &                   COMM, STATUS, IERR )
        FLAG = .TRUE.
      ELSE
        CALL MPI_IPROBE( MPI_ANY_SOURCE, MPI_ANY_TAG, COMM,
     &                   FLAG, STATUS, IERR )
      END IF
#if defined(ALLTIMERS)
      DKEEP(134)=DKEEP(134)+MPI_WTIME()-TIME_TMP
#endif
      IF ( FLAG ) THEN
         KEEP(266) = KEEP(266) -1
         MSGSOU = STATUS( MPI_SOURCE )
         MSGTAG = STATUS( MPI_TAG )
         CALL MPI_GET_COUNT( STATUS, MPI_PACKED, MSGLEN, IERR )
         IF ( MSGLEN .GT. LBUFR_BYTES ) THEN
           INFO(1) = -20
           INFO(2) = MSGLEN
           CALL DMUMPS_BDC_ERROR( MYID, SLAVEF, COMM, KEEP )
         ELSE
#if defined(ALLTIMERS)
           TIME_TMP=MPI_WTIME()
#endif
           CALL MPI_RECV( BUFR, LBUFR_BYTES, MPI_PACKED,
     &                  MSGSOU, MSGTAG, COMM, STATUS, IERR )
#if defined(ALLTIMERS)
           DKEEP(134)=DKEEP(134)+MPI_WTIME()-TIME_TMP
#endif
           CALL DMUMPS_TRAITER_MESSAGE_SOLVE( BUFR, LBUFR, LBUFR_BYTES,
     &          MSGTAG, MSGSOU, MYID, SLAVEF, COMM,
     &          N, NRHS, IPOOL, LPOOL, III, LEAF,
     &          NBFIN, NSTK_S, IW, LIW, A, LA, PTRIST, PTRFAC,
     &          IWCB, LIWCB,
     &          WCB, LWCB, POSWCB,
     &          PLEFTWCB, POSIWCB,
     &          PTRICB, INFO, KEEP,KEEP8, DKEEP, STEP,
     &          PROCNODE_STEPS, 
     &          RHSCOMP, LRHSCOMP, POSINRHSCOMP_FWD
     &          , FROM_PP )
         END IF
      END IF
      RETURN
      END SUBROUTINE DMUMPS_SOLVE_RECV_AND_TREAT
