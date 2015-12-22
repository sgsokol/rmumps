C
C  This file is part of MUMPS 5.0.1, released
C  on Thu Jul 23 17:08:29 UTC 2015
C
C
C  Copyright 1991-2015 CERFACS, CNRS, ENS Lyon, INP Toulouse, Inria,
C  University of Bordeaux.
C
C  This version of MUMPS is provided to you free of charge. It is
C  released under the CeCILL-C license:
C  http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html
C
      MODULE DMUMPS_FAC_PAR_M
      CONTAINS
      SUBROUTINE DMUMPS_FAC_PAR(N,IW,LIW,A,LA,
     &             NSTK_STEPS, NBPROCFILS,ND,FILS,STEP,
     &             FRERE, DAD, CAND,
     &             ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &             MAXFRT, NTOTPV, NMAXNPIV, PTRIST, PTRAST,
     &             PIMASTER, PAMASTER, PTRARW, PTRAIW,
     &             ITLOC, RHS_MUMPS, IPOOL, LPOOL,
     &             RINFO, POSFAC ,IWPOS, LRLU, IPTRLU,
     &             LRLUS, LEAF, NBROOT, NBRTOT,
     &             UU, ICNTL, PTLUST, PTRFAC, NSTEPS, INFO,
     &             KEEP,KEEP8,
     &             PROCNODE_STEPS,SLAVEF,MYID, COMM_NODES,
     &             MYID_NODES,
     &             BUFR,LBUFR,LBUFR_BYTES,INTARR,DBLARR,root,
     &             PERM, NELT, FRTPTR, FRTELT, LPTRAR,
     &             COMM_LOAD, ASS_IRECV, SEUIL, SEUIL_LDLT_NIV2,
     &             MEM_DISTRIB, NE,
     &             DKEEP,PIVNUL_LIST,LPN_LIST
     &             )
      USE DMUMPS_LOAD
      USE DMUMPS_OOC
      USE DMUMPS_FAC_ASM_MASTER_M
      USE DMUMPS_FAC_ASM_MASTER_ELT_M
      USE DMUMPS_FAC1_LDLT_M
      USE DMUMPS_FAC2_LDLT_M
      USE DMUMPS_FACTO_NIV1_M
      USE DMUMPS_FACTO_NIV2_M
      IMPLICIT NONE
      INCLUDE 'dmumps_root.h'
      TYPE (DMUMPS_ROOT_STRUC) :: root
      INTEGER N,NTOTPV,MAXFRT,LIW, LPTRAR, NMAXNPIV,
     &        NSTEPS, INFO(40)
      INTEGER(8) :: LA
      DOUBLE PRECISION, TARGET :: A(LA)
      INTEGER SLAVEF, COMM_NODES, MYID, MYID_NODES
      INTEGER, DIMENSION(0: SLAVEF - 1) :: MEM_DISTRIB
      INTEGER KEEP(500), ICNTL(40)
      INTEGER(8) KEEP8(150)
      INTEGER LPOOL
      INTEGER PROCNODE_STEPS(KEEP(28))
      INTEGER ITLOC(N+KEEP(253))
      DOUBLE PRECISION :: RHS_MUMPS(KEEP(255))
      INTEGER IW(LIW), NSTK_STEPS(KEEP(28)), NBPROCFILS(KEEP(28))
      INTEGER PTRARW(LPTRAR), PTRAIW(LPTRAR), ND(KEEP(28))
      INTEGER FILS(N),PTRIST(KEEP(28))
      INTEGER STEP(N), FRERE(KEEP(28)), DAD(KEEP(28))
      INTEGER PIMASTER(KEEP(28))
      INTEGER PTLUST(KEEP(28)), PERM(N)
      INTEGER CAND(SLAVEF+1,max(1,KEEP(56)))
      INTEGER   ISTEP_TO_INIV2(KEEP(71)),
     &          TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
      INTEGER IPOOL(LPOOL)
      INTEGER NE(KEEP(28))
      DOUBLE PRECISION RINFO(40)
      INTEGER(8) :: PAMASTER(KEEP(28)), PTRAST(KEEP(28))
      INTEGER(8) :: PTRFAC(KEEP(28))
      INTEGER(8) :: POSFAC, LRLU, LRLUS, IPTRLU
      INTEGER IWPOS, LEAF, NBROOT, NBRTOT
      INTEGER COMM_LOAD, ASS_IRECV
      DOUBLE PRECISION UU, SEUIL, SEUIL_LDLT_NIV2
      INTEGER NELT
      INTEGER FRTPTR( N+1 ), FRTELT( NELT )
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER BUFR( LBUFR )
      INTEGER        INTARR( max(1,KEEP(14)) )
      DOUBLE PRECISION DBLARR( max(1,KEEP(13)) )
      LOGICAL IS_ISOLATED_NODE
      INTEGER LPN_LIST
      INTEGER PIVNUL_LIST(LPN_LIST)
      DOUBLE PRECISION DKEEP(130)
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER :: STATUS(MPI_STATUS_SIZE)
      INTEGER :: IERR
      DOUBLE PRECISION, PARAMETER :: DZERO = 0.0D0, DONE = 1.0D0
      INTEGER INODE
      INTEGER IWPOSCB
      INTEGER FPERE, TYPEF
      INTEGER MP, LP, DUMMY(1)
      INTEGER NBFIN, NBROOT_TRAITEES
      INTEGER NFRONT, IOLDPS, NASS, HF, XSIZE
      INTEGER(8) NFRONT8
      INTEGER(8) :: POSELT
      INTEGER IPOSROOT, IPOSROOTROWINDICES
      INTEGER GLOBK109
      INTEGER(8) :: LBUFRX
      DOUBLE PRECISION, POINTER, DIMENSION(:) :: BUFRX
      LOGICAL :: IS_BUFRX_ALLOCATED
      DOUBLE PRECISION FLOP1
      INTEGER TYPE
      LOGICAL SON_LEVEL2, SET_IRECV, BLOCKING,
     &        MESSAGE_RECEIVED
      LOGICAL AVOID_DELAYED
      LOGICAL LAST_CALL
      INTEGER MASTER_ROOT
      INTEGER LOCAL_M, LOCAL_N
      INTEGER LRHS_CNTR_MASTER_ROOT, FWD_LOCAL_N_RHS
      LOGICAL ROOT_OWNER
      EXTERNAL MUMPS_TYPENODE, MUMPS_PROCNODE
      INTEGER MUMPS_TYPENODE, MUMPS_PROCNODE
      LOGICAL MUMPS_INSSARBR,MUMPS_ROOTSSARBR
      EXTERNAL MUMPS_INSSARBR,MUMPS_ROOTSSARBR
      LOGICAL DMUMPS_POOL_EMPTY
      EXTERNAL DMUMPS_POOL_EMPTY, DMUMPS_EXTRACT_POOL
      LOGICAL STACK_RIGHT_AUTHORIZED
      INTEGER numroc
      EXTERNAL numroc
      INTEGER MAXFRW, NPVW, NOFFW, NELVAW, COMP,
     &        JOBASS, ETATASS
      INTEGER(8) :: LAFAC
      INTEGER LIWFAC, STRAT, TYPEFile, NextPiv2beWritten,
     &        IDUMMY
      INTEGER(8) :: ITMP8
      TYPE(IO_BLOCK) :: MonBloc
      INCLUDE 'mumps_headers.h'
      DOUBLE PRECISION    OPASSW, OPELIW
      ASS_IRECV = MPI_REQUEST_NULL
      ITLOC(1:N+KEEP(253)) =0
      PTRIST(1:KEEP(28))=0
      PTLUST(1:KEEP(28))=0
      PTRAST(1:KEEP(28))=0_8
      PTRFAC(1:KEEP(28))=-99999_8
      MP = ICNTL(2)
      LP = ICNTL(1)
      MAXFRW = 0
      NPVW   = 0
      NOFFW  = 0
      NELVAW = 0
      COMP  = 0
      OPASSW = DZERO
      OPELIW = DZERO
      IWPOSCB = LIW
      STACK_RIGHT_AUTHORIZED = .TRUE.
      CALL DMUMPS_ALLOC_CB( .FALSE., 0_8,
     &     .FALSE., .FALSE., MYID_NODES, N, KEEP, KEEP8, DKEEP,
     &     IW, LIW, A, LA, LRLU, IPTRLU, IWPOS, IWPOSCB,
     &     PTRIST, PTRAST, STEP, PIMASTER,
     &     PAMASTER, KEEP(IXSZ), 0_8, -444, -444, .true.,
     &     COMP, LRLUS,
     &     INFO(1), INFO(2)
     &     )
      JOBASS  = 0
      ETATASS = 0
      NBFIN = NBRTOT
      NBROOT_TRAITEES = 0
      NBPROCFILS(1:KEEP(28)) = 0
#if ! defined(NO_XXNBPR)
      KEEP(121)=0
#endif
      IF ( KEEP(38).NE.0 ) THEN
        IF (root%yes) THEN
            CALL DMUMPS_ROOT_ALLOC_STATIC(
     &        root, KEEP(38), N, IW, LIW,
     &        A, LA,
     &        FILS, MYID_NODES, PTRAIW, PTRARW,
     &        INTARR, DBLARR,
     &        LRLU, IPTRLU,
     &        IWPOS, IWPOSCB, PTRIST, PTRAST,
     &        STEP, PIMASTER, PAMASTER, ITLOC, RHS_MUMPS,
     &        COMP, LRLUS, INFO(1), KEEP,KEEP8, DKEEP, INFO(2) )
        ENDIF
        IF ( INFO(1) .LT. 0 ) GOTO 635
      END IF
 20   CONTINUE
      NIV1_FLAG=0
      SET_IRECV = .TRUE.
      BLOCKING = .FALSE.
      MESSAGE_RECEIVED = .FALSE.
      CALL DMUMPS_TRY_RECVTREAT(
     &      COMM_LOAD, ASS_IRECV, BLOCKING, SET_IRECV,
     &      MESSAGE_RECEIVED,
     &      MPI_ANY_SOURCE, MPI_ANY_TAG,
     &      STATUS, BUFR, LBUFR,
     &      LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     &      IWPOS, IWPOSCB, IPTRLU,
     &      LRLU, LRLUS, N, IW, LIW, A, LA,
     &      PTRIST, PTLUST, PTRFAC,
     &      PTRAST, STEP, PIMASTER, PAMASTER, NSTK_STEPS,
     &      COMP, INFO(1), INFO(2), COMM_NODES, NBPROCFILS,
     &      IPOOL, LPOOL, LEAF, NBFIN, MYID_NODES, SLAVEF,
     &      root, OPASSW, OPELIW, ITLOC, RHS_MUMPS, FILS,
     &      PTRARW, PTRAIW,
     &      INTARR, DBLARR, ICNTL, KEEP,KEEP8,DKEEP, ND, FRERE,
     &      LPTRAR, NELT, FRTPTR, FRTELT,
     &      ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &      STACK_RIGHT_AUTHORIZED
     &       )
      CALL DMUMPS_LOAD_RECV_MSGS(COMM_LOAD, KEEP)
      IF (MESSAGE_RECEIVED) THEN
          IF ( INFO(1) .LT. 0 ) GO TO 640
          IF ( NBFIN .eq. 0 ) GOTO 640
      ELSE
          IF ( .NOT. DMUMPS_POOL_EMPTY( IPOOL, LPOOL) )THEN
            CALL DMUMPS_EXTRACT_POOL( N, IPOOL, LPOOL,
     &      PROCNODE_STEPS,
     &      SLAVEF, STEP, INODE, KEEP,KEEP8, MYID_NODES, ND,
     &      (.NOT. STACK_RIGHT_AUTHORIZED) )
            STACK_RIGHT_AUTHORIZED = .TRUE.
            IF (KEEP(47) .GE. 3) THEN
              CALL DMUMPS_LOAD_POOL_UPD_NEW_POOL(
     &              IPOOL, LPOOL,
     &              PROCNODE_STEPS, KEEP,KEEP8, SLAVEF, COMM_LOAD,
     &              MYID_NODES, STEP, N, ND, FILS )
            ENDIF
            IF (KEEP(47).EQ.4) THEN
               IF(INODE.GT.0.AND.INODE.LE.N)THEN
                  IF((NE(STEP(INODE)).EQ.0).AND.
     &                 (FRERE(STEP(INODE)).EQ.0))THEN
                     IS_ISOLATED_NODE=.TRUE.
                  ELSE
                     IS_ISOLATED_NODE=.FALSE.
                  ENDIF
               ENDIF
               CALL DMUMPS_LOAD_SBTR_UPD_NEW_POOL(
     &              IS_ISOLATED_NODE,INODE,IPOOL,LPOOL,
     &              MYID_NODES,SLAVEF,COMM_LOAD,KEEP,KEEP8)
            ENDIF
            IF ((( KEEP(80) == 2 .OR. KEEP(80)==3 ) .AND.
     &           ( KEEP(47) == 4 )).OR.
     &           (KEEP(80) == 1 .AND. KEEP(47) .GE. 1)) THEN
               CALL DMUMPS_UPPER_PREDICT(INODE,STEP,KEEP(28),
     &         PROCNODE_STEPS,FRERE,ND,COMM_LOAD,SLAVEF,
     &         MYID_NODES,KEEP,KEEP8,N)
            END IF
            GOTO 30
          ENDIF
      ENDIF
      GO TO 20
 30   CONTINUE
      IF ( INODE .LT. 0 ) THEN
        INODE = -INODE
        FPERE = DAD(STEP(INODE))
        GOTO 130
      ELSE IF (INODE.GT.N) THEN
       INODE = INODE - N
       IF (INODE.EQ.KEEP(38)) THEN
         NBROOT_TRAITEES = NBROOT_TRAITEES + 1
         IF ( NBROOT_TRAITEES .EQ. NBROOT ) THEN
            NBFIN = NBFIN - NBROOT
            IF (SLAVEF.GT.1) THEN
                DUMMY(1) = NBROOT
                CALL DMUMPS_MCAST2(DUMMY, 1, MPI_INTEGER, MYID_NODES,
     &          COMM_NODES, RACINE, SLAVEF)
            END IF
         ENDIF
         IF (NBFIN.EQ.0) GOTO 640
         GOTO 20
       ENDIF
       TYPE = MUMPS_TYPENODE(PROCNODE_STEPS(STEP(INODE)),SLAVEF)
       IF (TYPE.EQ.1) GOTO 100
       FPERE = DAD(STEP(INODE))
       AVOID_DELAYED = ( (FPERE .eq. KEEP(20) .OR. FPERE .eq. KEEP(38))
     &                   .AND. KEEP(60).ne.0 )
       IF ( KEEP(50) .eq. 0 ) THEN
         CALL  DMUMPS_FACTO_NIV2( COMM_LOAD, ASS_IRECV,
     &        N, INODE, FPERE, IW, LIW, A, LA, UU, NOFFW, NPVW,
     &        COMM_NODES, MYID_NODES, BUFR, LBUFR,LBUFR_BYTES,
     &        NBFIN,LEAF, INFO(1), INFO(2), IPOOL,LPOOL,
     &        SLAVEF, POSFAC, IWPOS, IWPOSCB, IPTRLU, LRLU,
     &        LRLUS, COMP, PTRIST, PTRAST, PTLUST, PTRFAC,
     &        STEP, PIMASTER, PAMASTER,
     &        NSTK_STEPS,NBPROCFILS,PROCNODE_STEPS,
     &        root, OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &        FILS, PTRARW, PTRAIW,
     &        INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE,
     &        LPTRAR, NELT, FRTPTR, FRTELT, SEUIL,
     &        ISTEP_TO_INIV2, TAB_POS_IN_PERE, AVOID_DELAYED,
     &        DKEEP(1),PIVNUL_LIST(1),LPN_LIST
     &             )
          IF ( INFO(1) .LT. 0 ) GOTO 640
       ELSE
         CALL  DMUMPS_FAC2_LDLT( COMM_LOAD, ASS_IRECV,
     &             N, INODE, FPERE, IW, LIW, A, LA, UU, NOFFW, NPVW,
     &             COMM_NODES, MYID_NODES, BUFR, LBUFR,LBUFR_BYTES,
     &             NBFIN,LEAF, INFO(1), INFO(2), IPOOL,LPOOL,
     &             SLAVEF, POSFAC, IWPOS, IWPOSCB, IPTRLU, LRLU,
     &             LRLUS, COMP, PTRIST, PTRAST, PTLUST, PTRFAC,
     &             STEP, PIMASTER, PAMASTER,
     &             NSTK_STEPS,NBPROCFILS,PROCNODE_STEPS,
     &             root, OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &             FILS, PTRARW, PTRAIW,
     &             INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE,
     &             LPTRAR, NELT, FRTPTR, FRTELT, SEUIL_LDLT_NIV2,
     &             ISTEP_TO_INIV2, TAB_POS_IN_PERE, AVOID_DELAYED,
     &              DKEEP(1),PIVNUL_LIST(1),LPN_LIST
     &             )
        IF ( INFO(1) .LT. 0 ) GOTO 640
        IF ( IW( PTLUST(STEP(INODE)) + KEEP(IXSZ) + 5 ) .GT. 1 ) THEN
             GOTO 20
        END IF
       END IF
       GOTO 130
      ENDIF
      IF (INODE.EQ.KEEP(38)) THEN
         CALL  DMUMPS_LAST_RTNELIND( COMM_LOAD, ASS_IRECV,
     &    root, FRERE,
     &    INODE,
     &    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     &    IWPOS, IWPOSCB, IPTRLU,
     &    LRLU, LRLUS, N, IW, LIW, A, LA, PTRIST,
     &    PTLUST, PTRFAC,
     &    PTRAST, STEP, PIMASTER, PAMASTER, NSTK_STEPS, COMP,
     &    INFO(1), INFO(2), COMM_NODES,
     &    NBPROCFILS,
     &    IPOOL, LPOOL, LEAF,
     &    NBFIN, MYID_NODES, SLAVEF,
     &
     &    OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &    FILS, PTRARW, PTRAIW,
     &    INTARR, DBLARR,ICNTL,KEEP,KEEP8,DKEEP, ND,
     &    LPTRAR, NELT, FRTPTR, FRTELT,
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE
     &      )
         IF ( INFO(1) .LT. 0 ) GOTO 640
        GOTO 20
      ENDIF
      TYPE = MUMPS_TYPENODE(PROCNODE_STEPS(STEP(INODE)),SLAVEF)
      IF (TYPE.EQ.1) THEN
        IF (KEEP(55).NE.0) THEN
         CALL DMUMPS_FAC_ASM_NIV1_ELT( COMM_LOAD, ASS_IRECV,
     &        NELT, FRTPTR, FRTELT,
     &        N,INODE,IW,LIW,A,LA,
     &        INFO(1),ND,
     &        FILS,FRERE,DAD,MAXFRW,root,OPASSW, OPELIW,
     &     PTRIST,PTLUST,PTRFAC,PTRAST,STEP, PIMASTER,PAMASTER,
     &        PTRARW,PTRAIW,
     &        ITLOC, RHS_MUMPS, NSTEPS, SON_LEVEL2,
     &        COMP, LRLU, IPTRLU,
     &        IWPOS,IWPOSCB, POSFAC, LRLUS,
     &        ICNTL, KEEP,KEEP8,DKEEP, INTARR, DBLARR,
     &    NSTK_STEPS,NBPROCFILS, PROCNODE_STEPS, SLAVEF,
     &    COMM_NODES, MYID_NODES,
     &    BUFR, LBUFR, LBUFR_BYTES, NBFIN, IPOOL, LPOOL, LEAF,
     &    PERM, ISTEP_TO_INIV2, TAB_POS_IN_PERE
     &    )
        ELSE
         JOBASS = 0
         CALL DMUMPS_FAC_ASM_NIV1(COMM_LOAD, ASS_IRECV,
     &        N,INODE,IW,LIW,A,LA,
     &        INFO(1),ND,
     &        FILS,FRERE,DAD,MAXFRW,root,OPASSW, OPELIW,
     &      PTRIST,PTLUST,PTRFAC,PTRAST,STEP, PIMASTER,PAMASTER,
     &        PTRARW,PTRAIW,
     &        ITLOC, RHS_MUMPS, NSTEPS, SON_LEVEL2,
     &        COMP, LRLU, IPTRLU,
     &        IWPOS,IWPOSCB, POSFAC, LRLUS,
     &        ICNTL, KEEP,KEEP8,DKEEP, INTARR, DBLARR,
     &    NSTK_STEPS,NBPROCFILS, PROCNODE_STEPS, SLAVEF,
     &    COMM_NODES, MYID_NODES,
     &    BUFR, LBUFR, LBUFR_BYTES, NBFIN, IPOOL, LPOOL, LEAF,
     &    PERM,
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE, JOBASS,ETATASS
     &    )
        ENDIF
       IF ( INFO(1) .LT. 0 ) GOTO 640
#if ! defined(NO_XXNBPR)
        CALL CHECK_EQUAL(NBPROCFILS(STEP(INODE)),
     &       IW(PTLUST(STEP(INODE))+XXNBPR) )
        IF ((IW(PTLUST(STEP(INODE))+XXNBPR).GT.0).OR.(SON_LEVEL2)) THEN
#else
        IF ((NBPROCFILS(STEP(INODE)).GT.0).OR.(SON_LEVEL2)) THEN
#endif
          GOTO 20
        ENDIF
      ELSE
        IF ( KEEP(55) .eq. 0 ) THEN
          CALL DMUMPS_FAC_ASM_NIV2(COMM_LOAD, ASS_IRECV,
     &    N, INODE, IW, LIW, A, LA,
     &    INFO(1),
     &    ND, FILS, FRERE, DAD, CAND,
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &    MAXFRW,
     &    root, OPASSW, OPELIW, PTRIST, PTLUST, PTRFAC,
     &    PTRAST, STEP, PIMASTER, PAMASTER, PTRARW, NSTK_STEPS,
     &    PTRAIW, ITLOC, RHS_MUMPS, NSTEPS,
     &    COMP, LRLU, IPTRLU, IWPOS, IWPOSCB, POSFAC, LRLUS,
     &    ICNTL, KEEP,KEEP8,DKEEP,INTARR,DBLARR,
     &    NBPROCFILS, PROCNODE_STEPS, SLAVEF, COMM_NODES,
     &    MYID_NODES,
     &    BUFR, LBUFR, LBUFR_BYTES,
     &    NBFIN, LEAF, IPOOL, LPOOL, PERM,
     &    MEM_DISTRIB(0)
     &    )
        ELSE
          CALL DMUMPS_FAC_ASM_NIV2_ELT( COMM_LOAD, ASS_IRECV,
     &    NELT, FRTPTR, FRTELT,
     &    N, INODE, IW, LIW, A, LA, INFO(1),
     &    ND, FILS, FRERE, DAD, CAND,
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &    MAXFRW,
     &    root, OPASSW, OPELIW, PTRIST, PTLUST, PTRFAC,
     &    PTRAST, STEP, PIMASTER, PAMASTER, PTRARW, NSTK_STEPS,
     &    PTRAIW, ITLOC, RHS_MUMPS, NSTEPS,
     &    COMP, LRLU, IPTRLU, IWPOS, IWPOSCB, POSFAC, LRLUS,
     &    ICNTL, KEEP,KEEP8,DKEEP,INTARR,DBLARR,
     &    NBPROCFILS, PROCNODE_STEPS, SLAVEF, COMM_NODES,
     &    MYID_NODES,
     &    BUFR, LBUFR, LBUFR_BYTES,
     &    NBFIN, LEAF, IPOOL, LPOOL, PERM,
     &    MEM_DISTRIB(0)
     &     )
        END IF
        IF (INFO(1).LT.0) GOTO 640
        GOTO 20
      ENDIF
 100  CONTINUE
       FPERE = DAD(STEP(INODE))
      IF ( INODE .eq. KEEP(20) ) THEN
        POSELT = PTRAST(STEP(INODE))
        IF (PTRFAC(STEP(INODE)).NE.POSELT) THEN
          WRITE(*,*) "ERROR 2 in DMUMPS_FAC_PAR", POSELT
          CALL MUMPS_ABORT()
        ENDIF
        CALL DMUMPS_CHANGE_HEADER
     &       ( IW(PTLUST(STEP(INODE))+KEEP(IXSZ)), KEEP(253) )
        GOTO 200
      END IF
      POSELT = PTRAST(STEP(INODE))
      IOLDPS = PTLUST(STEP(INODE))
      XSIZE =  KEEP(IXSZ)
      HF = 6 + IW(IOLDPS+5+XSIZE)+XSIZE
      NFRONT = IW(IOLDPS+XSIZE) 
      NASS   = iabs(IW(IOLDPS+2+XSIZE))
      AVOID_DELAYED = ( (FPERE .eq. KEEP(20) .OR. FPERE .eq. KEEP(38))
     &                   .AND. KEEP(60).ne.0 )
        IF (KEEP(50).EQ.0) THEN
           CALL DMUMPS_FACTO_NIV1 ( 
     &               N, INODE, IW, LIW, A, LA,
     &               IOLDPS, POSELT,
     &               INFO(1), UU, NOFFW, NPVW,
     &               KEEP,KEEP8,
     &               STEP, PROCNODE_STEPS, MYID_NODES, SLAVEF,
     &               SEUIL, AVOID_DELAYED, ETATASS,
     &              DKEEP(1),PIVNUL_LIST(1),LPN_LIST, IWPOS 
     &           )
      ELSE  
            IW( IOLDPS+4+KEEP(IXSZ) ) = 1
              CALL DMUMPS_FAC1_LDLT( N, INODE,
     &           IW, LIW, A, LA,
     &           IOLDPS, POSELT,
     &           INFO(1), UU, NOFFW, NPVW,
     &           KEEP,KEEP8, MYID_NODES, SEUIL, AVOID_DELAYED,
     &           ETATASS,
     &           DKEEP(1),PIVNUL_LIST(1),LPN_LIST, IWPOS
     &           )
            IW( IOLDPS+4+KEEP(IXSZ) ) = STEP(INODE)
          ENDIF 
          JOBASS = ETATASS          
          IF (JOBASS.EQ.1) THEN
              CALL DMUMPS_FAC_ASM_NIV1(COMM_LOAD, ASS_IRECV,
     &        N,INODE,IW,LIW,A,LA,
     &        INFO(1),ND,
     &        FILS,FRERE,DAD,MAXFRW,root,OPASSW, OPELIW,
     &        PTRIST,PTLUST,PTRFAC,PTRAST,STEP,PIMASTER,PAMASTER,
     &        PTRARW,PTRAIW,
     &        ITLOC, RHS_MUMPS, NSTEPS, SON_LEVEL2,
     &        COMP, LRLU, IPTRLU,
     &        IWPOS,IWPOSCB, POSFAC, LRLUS,
     &        ICNTL, KEEP,KEEP8,DKEEP, INTARR, DBLARR,
     &        NSTK_STEPS,NBPROCFILS, PROCNODE_STEPS, SLAVEF,
     &        COMM_NODES, MYID_NODES,
     &        BUFR, LBUFR, LBUFR_BYTES, NBFIN, IPOOL, LPOOL, LEAF,
     &        PERM,
     &        ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &        JOBASS,ETATASS
     &         )
          ENDIF
      IF (INFO(1).LT.0) GOTO 635
 130  CONTINUE
      TYPE  = MUMPS_TYPENODE(PROCNODE_STEPS(STEP(INODE)),SLAVEF)
      IF ( FPERE .NE. 0 ) THEN
        TYPEF = MUMPS_TYPENODE(PROCNODE_STEPS(STEP(FPERE)),SLAVEF)
      ELSE
        TYPEF = -9999
      END IF
      CALL DMUMPS_FAC_STACK( COMM_LOAD, ASS_IRECV,
     &       N,INODE,TYPE,TYPEF,LA,IW,LIW,A,
     &       INFO(1),INFO(2),OPELIW,NELVAW,NMAXNPIV,
     &       PTRIST,PTLUST,PTRFAC,
     &       PTRAST, STEP, PIMASTER, PAMASTER,
     &       NE, POSFAC,LRLU,
     &       LRLUS,IPTRLU,ICNTL,KEEP,KEEP8,DKEEP,COMP,IWPOS,IWPOSCB,
     &       PROCNODE_STEPS,SLAVEF,FPERE,COMM_NODES,MYID_NODES,
     &       IPOOL, LPOOL, LEAF,
     &       NSTK_STEPS, NBPROCFILS, BUFR, LBUFR, LBUFR_BYTES, NBFIN,
     &       root, OPASSW, ITLOC, RHS_MUMPS, FILS, PTRARW, PTRAIW,
     &       INTARR, DBLARR,
     &       ND, FRERE, LPTRAR, NELT, FRTPTR, FRTELT,
     &       ISTEP_TO_INIV2, TAB_POS_IN_PERE
     &       )
      IF (INFO(1).LT.0) GOTO 640
 200  CONTINUE
      IF ( INODE .eq. KEEP(38) ) THEN
        WRITE(*,*) 'Error .. in DMUMPS_FAC_PAR: ',
     &             ' INODE == KEEP(38)'
        Stop
      END IF
      IF ( FPERE.EQ.0 ) THEN
        NBROOT_TRAITEES = NBROOT_TRAITEES + 1
        IF ( NBROOT_TRAITEES .EQ. NBROOT ) THEN
           IF (KEEP(201).EQ.1) THEN 
              CALL DMUMPS_OOC_FORCE_WRT_BUF_PANEL(IERR)
           ELSE IF ( KEEP(201).EQ.2) THEN 
              CALL DMUMPS_FORCE_WRITE_BUF(IERR)
           ENDIF
            NBFIN = NBFIN - NBROOT
            IF ( NBFIN .LT. 0 ) THEN
              WRITE(*,*) ' ERROR 1 in DMUMPS_FAC_PAR: ',
     &                   ' NBFIN=', NBFIN
              CALL MUMPS_ABORT()
            END IF
            IF ( NBROOT .LT. 0 ) THEN
              WRITE(*,*) ' ERROR 1 in DMUMPS_FAC_PAR: ',
     &                   ' NBROOT=', NBROOT
              CALL MUMPS_ABORT()
            END IF
            IF (SLAVEF.GT.1) THEN
                DUMMY(1) = NBROOT
                CALL DMUMPS_MCAST2( DUMMY(1), 1, MPI_INTEGER,
     &          MYID_NODES, COMM_NODES, RACINE, SLAVEF)
            END IF
        ENDIF
        IF (NBFIN.EQ.0)THEN
           GOTO 640
        ENDIF
      ELSEIF ( FPERE.NE.KEEP(38) .AND.
     &         MUMPS_PROCNODE(PROCNODE_STEPS(STEP(FPERE)),SLAVEF).EQ.
     &         MYID_NODES ) THEN
        NSTK_STEPS(STEP(FPERE)) = NSTK_STEPS(STEP(FPERE))-1
        IF ( NSTK_STEPS( STEP( FPERE )).EQ.0) THEN
          IF (KEEP(234).NE.0 .AND.
     &      MUMPS_INSSARBR(PROCNODE_STEPS(STEP(INODE)),SLAVEF))
     &      THEN
            STACK_RIGHT_AUTHORIZED = .FALSE.
          ENDIF
          CALL DMUMPS_INSERT_POOL_N(N, IPOOL, LPOOL,
     &         PROCNODE_STEPS, SLAVEF, KEEP(28), KEEP(76),
     &         KEEP(80), KEEP(47), STEP, FPERE )
          IF (KEEP(47) .GE. 3) THEN
             CALL DMUMPS_LOAD_POOL_UPD_NEW_POOL(
     &            IPOOL, LPOOL,
     &            PROCNODE_STEPS, KEEP,KEEP8, SLAVEF, COMM_LOAD,
     &            MYID_NODES, STEP, N, ND, FILS )
          ENDIF
          CALL MUMPS_ESTIM_FLOPS( FPERE, N, PROCNODE_STEPS,SLAVEF,
     &           ND, FILS, FRERE, STEP, PIMASTER, KEEP(28),
     &           KEEP(50), KEEP(253), FLOP1,
     &           IW, LIW, KEEP(IXSZ) )
          IF (FPERE.NE.KEEP(20))
     &    CALL DMUMPS_LOAD_UPDATE(1,.FALSE.,FLOP1,KEEP,KEEP8)
        ENDIF
      ENDIF
      GO TO 20
 635  CONTINUE
      CALL DMUMPS_BDC_ERROR( MYID_NODES, SLAVEF, COMM_NODES )
 640  CONTINUE
        CALL DMUMPS_MPI_CANCEL( INFO(1),
     &       ASS_IRECV, BUFR, LBUFR,
     &       LBUFR_BYTES,
     &       COMM_NODES,
     &       MYID_NODES, SLAVEF)
       CALL DMUMPS_CLEAN_PENDING( INFO(1),
     &      BUFR, LBUFR,
     &      LBUFR_BYTES,
     &      COMM_NODES, COMM_LOAD, SLAVEF, MP)
      CALL MPI_BARRIER( COMM_NODES, IERR )
       IF ( INFO(1) .GE. 0 ) THEN
          IF( KEEP(38) .NE. 0 .OR. KEEP(20).NE.0) THEN
            MASTER_ROOT = MUMPS_PROCNODE(
     &                  PROCNODE_STEPS(STEP(max(KEEP(38),KEEP(20)))),
     &                  SLAVEF)
            ROOT_OWNER  = (MASTER_ROOT .EQ. MYID_NODES)
            IF ( KEEP(38) .NE. 0 ) THEN
               IF (KEEP(60).EQ.0) THEN
                 IOLDPS  = PTLUST(STEP(KEEP(38)))
                 LOCAL_M = IW(IOLDPS+2+KEEP(IXSZ))
                 LOCAL_N = IW(IOLDPS+1+KEEP(IXSZ))
               ELSE
                 IOLDPS  = -999
                 LOCAL_M = root%SCHUR_MLOC
                 LOCAL_N = root%SCHUR_NLOC
               ENDIF
               ITMP8   = int(LOCAL_M,8)*int(LOCAL_N,8)
               LBUFRX = min(int(root%MBLOCK,8)*int(root%NBLOCK,8),
     &            int(root%TOT_ROOT_SIZE,8)*int(root%TOT_ROOT_SIZE,8) )
               IF ( LRLU .GT. LBUFRX ) THEN
                   BUFRX => A(POSFAC:POSFAC+LRLU-1_8)
                   LBUFRX=LRLU
                   IS_BUFRX_ALLOCATED = .FALSE.
               ELSE
                   ALLOCATE( BUFRX( LBUFRX ), stat = IERR )
                   IF (IERR.gt.0) THEN
                         INFO(1) = -9
                         CALL MUMPS_SET_IERROR(LBUFRX, INFO(2) )
                         IF (LP > 0 )
     &                   write(LP,*) ' Error allocating, real array ',
     &                   'of size before DMUMPS_FACTO_ROOT',  LBUFRX
                         CALL MUMPS_ABORT()
                   ENDIF
                   IS_BUFRX_ALLOCATED = .FALSE.
               ENDIF
               CALL DMUMPS_FACTO_ROOT( MYID_NODES, MASTER_ROOT,
     &               root, N, KEEP(38),
     &               COMM_NODES, IW, LIW, IWPOS + 1,
     &               A, LA, PTRAST, PTLUST, PTRFAC, STEP,
     &               INFO(1), KEEP(50), KEEP(19),
     &               BUFRX(1), LBUFRX, KEEP,KEEP8, DKEEP,
     &               OPELIW )
               IF (IS_BUFRX_ALLOCATED) DEALLOCATE ( BUFRX )
               NULLIFY(BUFRX)
                IF ( MYID_NODES .eq. 
     &               MUMPS_PROCNODE(PROCNODE_STEPS(STEP(KEEP(38))),
     &                              SLAVEF)
     &             ) THEN
                   IF ( INFO(1) .EQ. -10 .OR. INFO(1) .EQ. -40 ) THEN
                      NPVW = NPVW + INFO(2)
                   ELSE
                      NPVW = NPVW + root%TOT_ROOT_SIZE
                      NMAXNPIV = max(NMAXNPIV,root%TOT_ROOT_SIZE)
                   END IF
                END IF
                IF (root%yes.AND.KEEP(60).EQ.0) THEN
                  IF (KEEP(252).EQ.0) THEN
                  IF (KEEP(201).EQ.1) THEN 
                    CALL MUMPS_GETI8(LAFAC, IW(IOLDPS+XXR))
                    LIWFAC    = IW(IOLDPS+XXI)
                    TYPEFile  = TYPEF_L
                    NextPiv2beWritten = 1 
                    MonBloc%INODE    = KEEP(38)   
                    MonBloc%MASTER   = .TRUE.
                    MonBloc%Typenode = 3
                    MonBloc%NROW     = LOCAL_M
                    MonBloc%NCOL     = LOCAL_N
                    MonBloc%NFS      = MonBloc%NCOL
                    MonBloc%Last     = .TRUE.   
                    MonBloc%LastPiv  =  MonBloc%NCOL
                    MonBloc%LastPanelWritten_L=-9999 
                    MonBloc%LastPanelWritten_U=-9999 
                    NULLIFY(MonBloc%INDICES)
                    STRAT        = STRAT_WRITE_MAX
                    MonBloc%Last = .TRUE.
                    LAST_CALL = .TRUE.
                    CALL DMUMPS_OOC_IO_LU_PANEL
     &                                 ( STRAT, TYPEFile,
     &                                  A(PTRFAC(STEP(KEEP(38)))),
     &                                  LAFAC, MonBloc,
     &                                  NextPiv2beWritten, IDUMMY,
     &                                  IW(IOLDPS), LIWFAC,
     &                                  MYID, KEEP8(31), IERR,LAST_CALL)
                  ELSE IF (KEEP(201).EQ.2) THEN
                    KEEP8(31)=KEEP8(31)+ ITMP8
                    CALL DMUMPS_NEW_FACTOR(KEEP(38),PTRFAC,
     &              KEEP,KEEP8,A,LA, ITMP8, IERR)
                    IF(IERR.LT.0)THEN
                      WRITE(*,*)MYID,
     &                ': Internal error in DMUMPS_NEW_FACTOR'
                      CALL MUMPS_ABORT()
                    ENDIF
                  ENDIF 
                  ENDIF 
                  IF (KEEP(201).NE.0 .OR. KEEP(252).NE.0) THEN
                     LRLUS = LRLUS + ITMP8 
                     IF (KEEP(252).NE.0) THEN
                       CALL DMUMPS_LOAD_MEM_UPDATE(.FALSE.,.FALSE.,
     &                 LA-LRLUS
     &                 ,0_8,-ITMP8,
     &                 KEEP,KEEP8,LRLUS)
                     ELSE         
                       CALL DMUMPS_LOAD_MEM_UPDATE(.FALSE.,.FALSE.,
     &                 LA-LRLUS
     &                 ,ITMP8,    
     &                 0_8,
     &                 KEEP,KEEP8,LRLUS)
                     ENDIF
                     IF (PTRFAC(STEP(KEEP(38))).EQ.POSFAC-ITMP8) THEN
                       POSFAC = POSFAC  - ITMP8
                       LRLU   = LRLU    + ITMP8
                     ENDIF
                  ELSE
                       CALL DMUMPS_LOAD_MEM_UPDATE(.FALSE.,.FALSE.,
     &                 LA-LRLUS
     &                 ,ITMP8,    
     &                 0_8,
     &                 KEEP,KEEP8,LRLUS)
                  ENDIF
                ENDIF  
                IF (root%yes. AND. KEEP(252) .NE. 0 .AND.
     &              (KEEP(60).EQ.0 .OR. KEEP(221).EQ.1)) THEN
                  IF (MYID_NODES .EQ. MASTER_ROOT) THEN
                    LRHS_CNTR_MASTER_ROOT = root%TOT_ROOT_SIZE*KEEP(253)
                  ELSE
                    LRHS_CNTR_MASTER_ROOT = 1
                  ENDIF
                  ALLOCATE(root%RHS_CNTR_MASTER_ROOT( 
     &                     LRHS_CNTR_MASTER_ROOT), stat=IERR )
                  IF (IERR.gt.0) THEN
                    INFO(1) = -13
                    CALL MUMPS_SET_IERROR(LRHS_CNTR_MASTER_ROOT,INFO(2))
                    IF (LP > 0 )
     &              write(LP,*) ' Error allocating, real array ',
     &              'of size before DMUMPS_FACTO_ROOT',
     &              LRHS_CNTR_MASTER_ROOT
                    CALL MUMPS_ABORT()
                  ENDIF
                  FWD_LOCAL_N_RHS = numroc(KEEP(253), root%NBLOCK,
     &            root%MYCOL, 0, root%NPCOL)
                  FWD_LOCAL_N_RHS = max(1,FWD_LOCAL_N_RHS)
                  CALL DMUMPS_GATHER_ROOT( MYID_NODES,
     &            root%TOT_ROOT_SIZE, KEEP(253),
     &            root%RHS_CNTR_MASTER_ROOT(1), LOCAL_M,
     &            FWD_LOCAL_N_RHS, root%MBLOCK, root%NBLOCK,
     &            root%RHS_ROOT(1,1), MASTER_ROOT,
     &            root%NPROW, root%NPCOL, COMM_NODES )
     &             
                ENDIF
            ELSE
                IF (KEEP(19).NE.0) THEN
                  CALL MPI_REDUCE(KEEP(109), GLOBK109, 1,
     &                 MPI_INTEGER, MPI_SUM,
     &                 MASTER_ROOT,
     &                 COMM_NODES, IERR)
                ENDIF
                IF (ROOT_OWNER) THEN
                   IPOSROOT = PTLUST(STEP(KEEP(20)))
                   NFRONT   = IW(IPOSROOT+KEEP(IXSZ)+3)   
                   NFRONT8  = int(NFRONT,8)
                   IPOSROOTROWINDICES=IPOSROOT+6+KEEP(IXSZ)+ 
     &                             IW(IPOSROOT+5+KEEP(IXSZ)) 
                   NPVW = NPVW + NFRONT 
                   NMAXNPIV = max(NMAXNPIV,NFRONT)
                END IF
               IF (ROOT_OWNER.AND.KEEP(60).NE.0) THEN 
                  IF ( PTRFAC(STEP(KEEP(20))) .EQ. POSFAC -
     &                 NFRONT8*NFRONT8 ) THEN
                    POSFAC = POSFAC - NFRONT8*NFRONT8
                    LRLUS  = LRLUS  + NFRONT8*NFRONT8
                    LRLU   = LRLUS  + NFRONT8*NFRONT8
                    CALL DMUMPS_LOAD_MEM_UPDATE(.FALSE.,.FALSE.,
     &              LA-LRLUS,0_8,-NFRONT8*NFRONT8,KEEP,KEEP8,LRLUS)
                  ENDIF
               ENDIF
            END IF
          END IF
       END IF
       IF ( KEEP(38) .NE. 0 ) THEN
         IF (MYID_NODES.EQ. 
     &        MUMPS_PROCNODE(PROCNODE_STEPS(STEP(KEEP(38))),SLAVEF)
     &      ) THEN
           MAXFRW = max ( MAXFRW, root%TOT_ROOT_SIZE)
         END IF
       END IF
       MAXFRT       = MAXFRW
       NTOTPV       = NPVW
       INFO(12)     = NOFFW
       RINFO(2)     = dble(OPASSW)
       RINFO(3)     = dble(OPELIW)
       INFO(13)     = NELVAW
       INFO(14)     = COMP
      RETURN
      END SUBROUTINE DMUMPS_FAC_PAR
      SUBROUTINE DMUMPS_CHANGE_HEADER( HEADER, KEEP253 )
        INTEGER HEADER( 6 ), KEEP253
        INTEGER NFRONT, NASS
        NFRONT = HEADER(1)
        IF ( HEADER(2) .ne. 0 ) THEN
          WRITE(*,*) ' *** CHG_HEADER ERROR 1 :',HEADER(2)
          CALL MUMPS_ABORT()
        END IF
        NASS   = abs( HEADER( 3 ) )
        IF ( NASS .NE. abs( HEADER( 4 ) ) ) THEN
          WRITE(*,*) ' *** CHG_HEADER ERROR 2 :',HEADER(3:4)
          CALL MUMPS_ABORT()
        END IF
        IF ( NASS+KEEP253 .NE. NFRONT ) THEN
          WRITE(*,*) ' *** CHG_HEADER ERROR 3 : not root',
     &    NASS, KEEP253, NFRONT
          CALL MUMPS_ABORT()
        END IF
        HEADER( 1 ) = KEEP253 
        HEADER( 2 ) = 0
        HEADER( 3 ) = NFRONT 
        HEADER( 4 ) = NFRONT-KEEP253    
        RETURN
      END SUBROUTINE DMUMPS_CHANGE_HEADER
      END MODULE DMUMPS_FAC_PAR_M