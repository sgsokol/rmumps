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
      MODULE DMUMPS_FAC2_LDLT_M
      CONTAINS
      SUBROUTINE DMUMPS_FAC2_LDLT( COMM_LOAD, ASS_IRECV, 
     &           N, INODE, FPERE, IW, LIW, A, LA,
     &           UU, NOFFW,
     &           NPVW,
     &             COMM, MYID, BUFR, LBUFR,LBUFR_BYTES,NBFIN,LEAF,
     &             IFLAG, IERROR, IPOOL,LPOOL,
     &             SLAVEF, POSFAC, IWPOS, IWPOSCB, IPTRLU, LRLU,
     &             LRLUS, COMP,
     &             PTRIST, PTRAST, PTLUST_S, PTRFAC, STEP,
     &             PIMASTER, PAMASTER,
     &             NSTK_S,NBPROCFILS,PROCNODE_STEPS, root,
     &             OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &             FILS, PTRARW, PTRAIW,
     &             INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE,
     &             LPTRAR, NELT, FRTPTR, FRTELT, SEUIL,
     &             ISTEP_TO_INIV2, TAB_POS_IN_PERE, AVOID_DELAYED,
     &             DKEEP,PIVNUL_LIST,LPN_LIST
     &     , LRGROUPS
     &             )
      USE DMUMPS_FAC_FRONT_AUX_M
      USE DMUMPS_FAC_FRONT_TYPE2_AUX_M
      USE DMUMPS_OOC
      USE DMUMPS_FAC_LR
      USE DMUMPS_LR_TYPE
      USE DMUMPS_LR_STATS
      USE DMUMPS_ANA_LR
!$    USE OMP_LIB
      IMPLICIT NONE
      INCLUDE 'dmumps_root.h'
      INTEGER COMM_LOAD, ASS_IRECV
      INTEGER N, INODE, FPERE, LIW, NOFFW, NPVW
      INTEGER(8) :: LA
      INTEGER, TARGET :: IW( LIW )
      DOUBLE PRECISION A( LA )
      DOUBLE PRECISION UU, SEUIL
      TYPE (DMUMPS_ROOT_STRUC) :: root
      INTEGER COMM, MYID, LBUFR, LBUFR_BYTES
      INTEGER LPTRAR, NELT
      INTEGER ICNTL(40), KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER NBFIN, SLAVEF,
     &        IFLAG, IERROR, LEAF, LPOOL
      INTEGER(8) :: POSFAC, IPTRLU, LRLU, LRLUS
      INTEGER IWPOS, IWPOSCB, COMP 
      INTEGER NB_BLOC_FAC
      INTEGER FRTPTR( N + 1 ), FRTELT( NELT )
      INTEGER BUFR( LBUFR ), IPOOL(LPOOL),
     &        ITLOC(N+KEEP(253)), FILS(N),
     &        ND( KEEP(28) ), FRERE( KEEP(28) )
      INTEGER(8), INTENT(IN) :: PTRARW(LPTRAR), PTRAIW(LPTRAR)
      DOUBLE PRECISION :: RHS_MUMPS(KEEP(255))
      INTEGER(8) :: PTRAST(KEEP(28))
      INTEGER(8) :: PTRFAC(KEEP(28))
      INTEGER(8) :: PAMASTER(KEEP(28))
      INTEGER PTRIST(KEEP(28)), PTLUST_S(KEEP(28)),
     &        STEP(N), PIMASTER(KEEP(28)),
     &        NSTK_S(KEEP(28)), NBPROCFILS(KEEP(28)),
     &        PROCNODE_STEPS(KEEP(28))
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     &        TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
      DOUBLE PRECISION OPASSW, OPELIW
      DOUBLE PRECISION DBLARR(KEEP8(26))
      INTEGER INTARR(KEEP8(27))
      LOGICAL AVOID_DELAYED
      INTEGER LPN_LIST
      INTEGER PIVNUL_LIST(LPN_LIST)
      DOUBLE PRECISION DKEEP(230)
      INTEGER :: LRGROUPS(N)
      INTEGER(8) :: POSELT
      INTEGER IOLDPS, allocok, K263
      INTEGER INOPV, IFINB, NFRONT, NPIV, IEND_BLOCK
      INTEGER NASS, LDAFS, IBEG_BLOCK
      INTEGER :: IBEG_BLOCK_FOR_IPIV
      LOGICAL LASTBL, LR_ACTIVATED
      INTEGER NBLR_ORIG, IBEG_BLR, IEND_BLR, CURRENT_BLR
      INTEGER Inextpiv
      LOGICAL RESET_TO_ONE
      INTEGER K109_SAVE
      INTEGER XSIZE, NBKJIB_ORIG
      DOUBLE PRECISION UUTEMP
      INCLUDE 'mumps_headers.h'
      INTEGER , ALLOCATABLE, DIMENSION ( : ) :: IPIV
      DOUBLE PRECISION , ALLOCATABLE, DIMENSION ( : )    :: DIAG_ORIG
      INTEGER    :: SIZEDIAG_ORIG
      INTEGER(8) :: LAFAC
      INTEGER LIWFAC, STRAT, TYPEFile, NextPiv2beWritten,
     &        IDUMMY, NELIM
      TYPE(IO_BLOCK) :: MonBloc 
      LOGICAL LAST_CALL
      INTEGER PP_FIRST2SWAP_L, IFLAG_OOC
      INTEGER PP_LastPIVRPTRFilled 
      INTEGER HF, NPARTSASS, NPARTSCB, NB_BLR
      INTEGER MAXI_CLUSTER, LWORK
      INTEGER T1, T2, COUNT_RATE, T1P, T2P, CRP
      INTEGER TTOT1, TTOT2, COUNT_RATETOT
      INTEGER TTOT1FR, TTOT2FR, COUNT_RATETOTFR
      DOUBLE PRECISION :: LOC_UPDT_TIME, 
     &        LOC_PROMOTING_TIME, LOC_DEMOTING_TIME,
     &        LOC_CB_DEMOTING_TIME, LOC_FRPANELS_TIME,
     &        LOC_TRSM_TIME,
     &        LOC_FRFRONTS_TIME, LOC_LR_MODULE_TIME,
     &        LOC_FAC_I_TIME, LOC_FAC_MQ_TIME, LOC_FAC_SQ_TIME
      INTEGER, POINTER, DIMENSION(:) :: BEGS_BLR
      TYPE(LRB_TYPE), POINTER, DIMENSION(:) :: BLR_L, BLR_SEND
      TYPE(LRB_TYPE), DIMENSION(1), TARGET  :: BLR_DUMMY
      DOUBLE PRECISION, ALLOCATABLE :: WORK(:), TAU(:)
      INTEGER, ALLOCATABLE :: JPVT(:)
      DOUBLE PRECISION, ALLOCATABLE :: RWORK(:)
      DOUBLE PRECISION, ALLOCATABLE :: BLOCK(:,:)
      INTEGER :: OMP_NUM
      INTEGER PIVOT_OPTION
      EXTERNAL DMUMPS_BDC_ERROR
      LOGICAL STATICMODE
      DOUBLE PRECISION SEUIL_LOC
      DOUBLE PRECISION GW_FACTCUMUL
      INTEGER PIVSIZ,IWPOSPIV
      DOUBLE PRECISION ONE
      PARAMETER (ONE = 1.0D0)
      NULLIFY(BLR_L) 
      LOC_UPDT_TIME = 0.D0
      LOC_PROMOTING_TIME = 0.D0
      LOC_DEMOTING_TIME = 0.D0
      LOC_CB_DEMOTING_TIME = 0.D0
      LOC_FRPANELS_TIME = 0.0D0
      LOC_FRFRONTS_TIME = 0.0D0
      LOC_TRSM_TIME = 0.D0
      LOC_LR_MODULE_TIME = 0.D0
      LOC_FAC_I_TIME = 0.D0
      LOC_FAC_MQ_TIME = 0.D0
      LOC_FAC_SQ_TIME = 0.D0
      IF (KEEP(206).GE.1) THEN
        Inextpiv = 1   
      ELSE 
        Inextpiv = 0   
      ENDIF
      INOPV = 0
      IF(KEEP(97) .EQ. 0) THEN
         STATICMODE = .FALSE.
      ELSE
         STATICMODE = .TRUE.
      ENDIF
      IF (AVOID_DELAYED) THEN
        STATICMODE = .TRUE.
        UUTEMP=UU
        SEUIL_LOC = max(SEUIL,epsilon(SEUIL))
      ELSE
        SEUIL_LOC=SEUIL
        UUTEMP=UU
      ENDIF
      PIVOT_OPTION = MIN(2,KEEP(468))
      IF (UUTEMP == 0.0D0 .AND. KEEP(201).NE.1) THEN
      ENDIF
      RESET_TO_ONE = ((KEEP(110).GT.0).AND.(DKEEP(2).LE.0.0D0))
      IF (RESET_TO_ONE) THEN
        K109_SAVE = KEEP(109)
      ENDIF
      IBEG_BLOCK  = 1
      NB_BLOC_FAC = 0
      XSIZE  = KEEP(IXSZ)
      IOLDPS = PTLUST_S(STEP( INODE ))
      POSELT = PTRAST(STEP( INODE ))
      NFRONT = IW(IOLDPS+XSIZE)
      NASS   = iabs(IW(IOLDPS+2+XSIZE))
      LDAFS  = NASS
      IW(IOLDPS+3+XSIZE) =  -99999
      LR_ACTIVATED= .FALSE. 
      NULLIFY(BEGS_BLR)
      LR_ACTIVATED   = (IW(IOLDPS+XXLR).GT.0)
      IF (NASS.LT.KEEP(4)) THEN
        NBKJIB_ORIG = NASS
      ELSE IF (NASS .GT. KEEP(3)) THEN
        NBKJIB_ORIG = min( KEEP(6), NASS )
      ELSE
        NBKJIB_ORIG = min( KEEP(5), NASS )
      ENDIF
      IF (.not.LR_ACTIVATED) THEN
          NBLR_ORIG     = KEEP(420)
      ELSE
          NBLR_ORIG  = -9999 
      ENDIF
      IF (LR_ACTIVATED) THEN
         K263 = 1   
      ELSE
         K263 = KEEP(263)
         IF (K263 .NE. 0 .AND. NASS/NBLR_ORIG < 4) THEN
           K263 = 0
         ENDIF
      ENDIF
      IEND_BLOCK  = 0
      IEND_BLR    = 0
      CURRENT_BLR = 0
      ALLOCATE( IPIV( NASS ), stat = allocok )
      IF ( allocok .GT. 0 ) THEN
        WRITE(*,*) MYID, ' : DMUMPS_FAC2_LDLT failed to allocate ',
     &  NASS, ' integers'
        IFLAG = -13
        IERROR=NASS
        GO TO 490
      END IF
      IF (KEEP(219).GE.3) THEN
       SIZEDIAG_ORIG = NASS
      ELSE
       SIZEDIAG_ORIG = 1
      ENDIF
      ALLOCATE ( DIAG_ORIG(SIZEDIAG_ORIG), stat = allocok )
      IF ( allocok .GT. 0 ) THEN
          WRITE(*,*) MYID,
     &      ' : FAC_NIV2 failed to allocate ',
     &      NASS, ' REAL/COMPLEX entries'
          IFLAG=-13
          IERROR=NASS
          GO TO 490
      END IF
      IF (KEEP(201).EQ.1) THEN 
        IDUMMY    = -9876
        CALL MUMPS_GETI8(LAFAC,IW(IOLDPS+XXR))
        LIWFAC    = IW(IOLDPS+XXI)
        TYPEFile  = TYPEF_L
        NextPiv2beWritten = 1 
        PP_FIRST2SWAP_L = NextPiv2beWritten 
        MonBloc%LastPanelWritten_L = 0 
        MonBloc%INODE    = INODE
        MonBloc%MASTER   = .TRUE.
        MonBloc%Typenode = 2
        MonBloc%NROW     = NASS
        MonBloc%NCOL     = NASS
        MonBloc%NFS      = NASS
        MonBloc%Last     = .FALSE.   
        MonBloc%LastPiv  = -66666    
        MonBloc%INDICES =>
     &  IW(IOLDPS+6+NFRONT+XSIZE+IW(IOLDPS+5+XSIZE)
     &    :IOLDPS+5+2*NFRONT+XSIZE+IW(IOLDPS+5+XSIZE))
      ENDIF
      IF (LR_ACTIVATED) THEN
             CNT_NODES = CNT_NODES + 1 
             CALL SYSTEM_CLOCK(TTOT1)
      ELSE IF (KEEP(486).GT.0) THEN
             CALL SYSTEM_CLOCK(TTOT1FR)
      ENDIF
      HF = 6 + IW(IOLDPS+5+XSIZE)+XSIZE
      IF (KEEP(201).EQ.1) THEN  
       IF (PIVOT_OPTION.LT.2) PIVOT_OPTION=2
      ENDIF
      IF (LR_ACTIVATED) THEN
         CALL GET_CUT(IW(IOLDPS+HF:IOLDPS+HF+NFRONT-1), NASS,
     &        0, LRGROUPS, NPARTSCB, 
     &        NPARTSASS, BEGS_BLR)
         CALL REGROUPING2(BEGS_BLR, NPARTSASS, NASS, NPARTSCB,
     &        0, KEEP(488), .FALSE., KEEP(472))     
         NB_BLR = NPARTSASS + NPARTSCB
         call MAX_CLUSTER(BEGS_BLR,NB_BLR,MAXI_CLUSTER)
         LWORK = MAXI_CLUSTER*MAXI_CLUSTER
         OMP_NUM = 1
#if defined(BLR_MT)
!$       OMP_NUM = OMP_GET_MAX_THREADS()
#endif
         ALLOCATE(BLOCK(MAXI_CLUSTER, OMP_NUM*MAXI_CLUSTER),
     &             RWORK(2*MAXI_CLUSTER*OMP_NUM), 
     &             TAU(MAXI_CLUSTER*OMP_NUM),
     &             JPVT(MAXI_CLUSTER*OMP_NUM), 
     &             WORK(LWORK*OMP_NUM),stat=allocok)
         IF (allocok > 0) THEN
           IFLAG  = -13
           IERROR = OMP_NUM*(LWORK + MAXI_CLUSTER*(MAXI_CLUSTER+4))
           GOTO 480
         ENDIF
      ENDIF
      LASTBL = .FALSE.
      DO WHILE (IEND_BLR < NASS ) 
        CURRENT_BLR = CURRENT_BLR + 1
        IBEG_BLR = IW(IOLDPS+1+KEEP(IXSZ)) + 1 
        IF (.NOT. LR_ACTIVATED)THEN
          IEND_BLR = min(IEND_BLR + NBLR_ORIG, NASS)
        ELSE
          IEND_BLR = min(BEGS_BLR(CURRENT_BLR+1)-1, NASS)
          BEGS_BLR( CURRENT_BLR ) = IBEG_BLR
          IF ( IEND_BLR - IBEG_BLR + 1 .GT. MAXI_CLUSTER ) THEN
            MAXI_CLUSTER = IEND_BLR - IBEG_BLR + 1
            LWORK = MAXI_CLUSTER*MAXI_CLUSTER
            DEALLOCATE(BLOCK, WORK, RWORK, TAU, JPVT)
            ALLOCATE(BLOCK(MAXI_CLUSTER, OMP_NUM*MAXI_CLUSTER),
     &             RWORK(2*MAXI_CLUSTER*OMP_NUM), 
     &             TAU(MAXI_CLUSTER*OMP_NUM),
     &             JPVT(MAXI_CLUSTER*OMP_NUM), 
     &             WORK(LWORK*OMP_NUM),stat=allocok)
            IF (allocok > 0) THEN
              IFLAG  = -13
              IERROR = OMP_NUM*(LWORK + MAXI_CLUSTER*(MAXI_CLUSTER+4))
              GOTO 480
            ENDIF
          ENDIF
        ENDIF
        IF (LR_ACTIVATED) THEN
          CALL SYSTEM_CLOCK(T1)
        ENDIF
        DO WHILE (IEND_BLOCK < IEND_BLR ) 
          IBEG_BLOCK = IW(IOLDPS+1+KEEP(IXSZ)) + 1
          IEND_BLOCK = min(IEND_BLOCK + NBKJIB_ORIG, IEND_BLR)
  50      CONTINUE  
            IF (K263.EQ.0) THEN
              IBEG_BLOCK_FOR_IPIV = IBEG_BLOCK
            ELSE
              IBEG_BLOCK_FOR_IPIV = IBEG_BLR
            ENDIF
            IF (LR_ACTIVATED) THEN
                CALL SYSTEM_CLOCK(T1P)
            ENDIF
            CALL DMUMPS_FAC_I_LDLT_NIV2(
     &                DIAG_ORIG, SIZEDIAG_ORIG, GW_FACTCUMUL,
     &                NFRONT,NASS,IBEG_BLOCK_FOR_IPIV,
     &                IBEG_BLOCK, IEND_BLOCK,
     &                NASS, IPIV,
     &                N,INODE,IW,LIW,A,LA,NOFFW,INOPV,
     &                IFLAG,IOLDPS,POSELT,UU, SEUIL_LOC,
     &                KEEP,KEEP8,PIVSIZ,
     &           DKEEP(1),PIVNUL_LIST(1),LPN_LIST,
     &           PP_FIRST2SWAP_L, MonBloc%LastPanelWritten_L,
     &           PP_LastPIVRPTRFilled, 
     &           PIVOT_OPTION,
     &           Inextpiv, IEND_BLR)
            IF (LR_ACTIVATED) THEN
                CALL SYSTEM_CLOCK(T2P,CRP)
                LOC_FAC_I_TIME = LOC_FAC_I_TIME + 
     &                           dble(T2P-T1P)/dble(CRP)
            ENDIF
            IF (IFLAG.LT.0) GOTO 490
            IF(KEEP(109).GT. 0) THEN
              IF(PIVNUL_LIST(KEEP(109)).EQ.-1) THEN
                IWPOSPIV = IOLDPS+IW(IOLDPS+1+XSIZE)+6
     &              +IW(IOLDPS+5+XSIZE)
                PIVNUL_LIST(KEEP(109)) = IW(IWPOSPIV+XSIZE)
              ENDIF
            ENDIF
            IF (INOPV.EQ. 1) THEN
              IF (STATICMODE) THEN
                INOPV = -1
                GOTO 50
              ENDIF
             LASTBL = .TRUE.
            ELSE IF (INOPV .LE. 0) THEN 
              NPVW = NPVW + PIVSIZ
              IF (LR_ACTIVATED) THEN
                  CALL SYSTEM_CLOCK(T1P)
              ENDIF
              CALL DMUMPS_FAC_MQ_LDLT_NIV2(IEND_BLOCK,
     &             NASS, IW(IOLDPS+1+XSIZE), INODE,A,LA,
     &             LDAFS, POSELT,IFINB,
     &             PIVSIZ,
     &             KEEP(219),
     &             PIVOT_OPTION, IEND_BLR)
              IF (LR_ACTIVATED) THEN
                  CALL SYSTEM_CLOCK(T2P,CRP)
                  LOC_FAC_MQ_TIME = LOC_FAC_MQ_TIME + 
     &                              dble(T2P-T1P)/dble(CRP)
              ENDIF
              IF(PIVSIZ .EQ. 2) THEN
                IWPOSPIV = IOLDPS+XSIZE+IW(IOLDPS+1+XSIZE)+6+
     &                     IW(IOLDPS+5+XSIZE)
                IW(IWPOSPIV+NFRONT) = -IW(IWPOSPIV+NFRONT)
              ENDIF
              IW(IOLDPS+1+XSIZE) = IW(IOLDPS+1+XSIZE) + PIVSIZ
            IF (IFINB.EQ.0) THEN
              GOTO 50 
            ELSE IF (IFINB .EQ. -1) THEN
              LASTBL = .TRUE.
            ENDIF
          ENDIF
          NPIV = IW(IOLDPS+1+XSIZE)
          IF (  KEEP(201).EQ.1.AND.(PIVOT_OPTION.GE.3)
     &        .AND.
     &          ( .NOT. LR_ACTIVATED .OR.
     &            ( (KEEP(485).EQ.0) .AND. (PIVOT_OPTION.GT.2) )
     &          ) 
     &       ) THEN  
            IF (.NOT.RESET_TO_ONE.OR.K109_SAVE.EQ.KEEP(109)) THEN
              MonBloc%Last   = .FALSE.
              MonBloc%LastPiv= NPIV
              LAST_CALL=.FALSE.
              CALL DMUMPS_OOC_IO_LU_PANEL(
     &        STRAT_TRY_WRITE,
     &        TYPEFile, A(POSELT),
     &        LAFAC, MonBloc, NextPiv2beWritten, IDUMMY, IW(IOLDPS),
     &        LIWFAC, MYID, KEEP8(31), IFLAG_OOC,LAST_CALL )
              IF (IFLAG_OOC .LT. 0 ) IFLAG = IFLAG_OOC
              IF (IFLAG .LT. 0 ) RETURN
            ENDIF
          ENDIF
          IF (K263.eq.0) THEN
            NELIM = IEND_BLR-NPIV
            CALL DMUMPS_SEND_FACTORED_BLK( COMM_LOAD, ASS_IRECV,
     &             N, INODE, FPERE, IW, LIW, 
     &             IOLDPS, POSELT, A, LA, LDAFS,
     &             IBEG_BLOCK, NPIV, IPIV, NASS,LASTBL, NB_BLOC_FAC,
     &             COMM, MYID, BUFR, LBUFR, LBUFR_BYTES,NBFIN,LEAF,
     &             IFLAG, IERROR, IPOOL,LPOOL,
     &             SLAVEF, POSFAC, IWPOS, IWPOSCB, IPTRLU, LRLU,
     &             LRLUS, COMP,
     &             PTRIST, PTRAST, PTLUST_S, PTRFAC, STEP,
     &             PIMASTER, PAMASTER,
     &             NSTK_S,NBPROCFILS,PROCNODE_STEPS, root,
     &             OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &             FILS, PTRARW, PTRAIW,
     &             INTARR,DBLARR,ICNTL,KEEP,KEEP8,DKEEP,ND,FRERE,
     &             LPTRAR, NELT, FRTPTR, FRTELT, 
     &             ISTEP_TO_INIV2, TAB_POS_IN_PERE
     &             , NELIM, .FALSE. 
     &             , NPARTSASS, CURRENT_BLR, BLR_DUMMY, LRGROUPS
     & )
            IF ( IFLAG .LT. 0 ) GOTO 500
            IF (RESET_TO_ONE.AND.K109_SAVE.LT.KEEP(109)) THEN
              CALL DMUMPS_RESET_TO_ONE( 
     &        IW(IOLDPS+KEEP(IXSZ)+IW(IOLDPS+5+KEEP(IXSZ))+6),
     &        NPIV, IBEG_BLOCK,
     &        K109_SAVE, KEEP(109), PIVNUL_LIST, LPN_LIST,
     &        A, POSELT, LA, LDAFS)
            ENDIF
            IF (KEEP(201).EQ.1.AND.(PIVOT_OPTION.GE.3)) THEN
              MonBloc%Last  = .FALSE.
              MonBloc%LastPiv= NPIV
              LAST_CALL=.FALSE.
              CALL DMUMPS_OOC_IO_LU_PANEL(
     &        STRAT_TRY_WRITE,
     &        TYPEFile, A(POSELT),
     &        LAFAC, MonBloc, NextPiv2beWritten, IDUMMY, IW(IOLDPS),
     &        LIWFAC, MYID, KEEP8(31), IFLAG_OOC,LAST_CALL )
              IF (IFLAG_OOC .LT. 0 ) THEN
                IFLAG = IFLAG_OOC
                RETURN
              ENDIF
            ENDIF
          ENDIF
          IF ( IEND_BLR .GT. IEND_BLOCK ) THEN
              IF (LR_ACTIVATED) THEN
                  CALL SYSTEM_CLOCK(T1P)
              ENDIF
              CALL DMUMPS_FAC_SQ_LDLT(IBEG_BLOCK,IEND_BLOCK,NPIV,
     &             NASS,NASS,IEND_BLR,INODE,A,LA,
     &             LDAFS, POSELT,
     &             KEEP,KEEP8,
     &             PIVOT_OPTION, .FALSE.)
              IF (LR_ACTIVATED) THEN
                  CALL SYSTEM_CLOCK(T2P,CRP)
                  LOC_FAC_SQ_TIME = LOC_FAC_SQ_TIME +
     &                              dble(T2P-T1P)/dble(CRP)
              ENDIF
          ENDIF
        END DO 
        NPIV   = IW(IOLDPS+1+XSIZE)
        IF (LR_ACTIVATED) THEN
          CALL SYSTEM_CLOCK(T2,COUNT_RATE)
          LOC_FRPANELS_TIME = LOC_FRPANELS_TIME +
     &                        dble(T2-T1)/dble(COUNT_RATE)
          CALL UPDATE_FLOP_STATS_PANEL(NFRONT - IBEG_BLR + 1,
     &                                 NPIV   - IBEG_BLR + 1, 2, 1)
        ENDIF
        IF (LR_ACTIVATED) THEN
          NELIM = IEND_BLOCK - NPIV
          IF (IEND_BLR.NE.IEND_BLOCK) THEN
            WRITE(*,*) "Internal error 1 in DMUMPS_FAC2_LDLT",
     &      IEND_BLR, IEND_BLOCK
            CALL MUMPS_ABORT()
          ENDIF
          IF (NELIM .EQ. IEND_BLR - IBEG_BLR + 1) THEN
            GOTO 101
          ENDIF
          ALLOCATE(BLR_L(NB_BLR-CURRENT_BLR)) 
          CALL SYSTEM_CLOCK(T1)
#if defined(BLR_MT)          
!$OMP PARALLEL 
#endif
          CALL DMUMPS_COMPRESS_PANEL(A, LA, POSELT, IFLAG, IERROR, NASS,
     &         BEGS_BLR, NB_BLR, DKEEP(8), KEEP(473), BLR_L, 
     &         CURRENT_BLR, 'V', WORK, TAU, JPVT, LWORK, RWORK,
     &         BLOCK, MAXI_CLUSTER, NELIM, 
     &         .FALSE., 0, 0,
     &         2, KEEP(483), KEEP(470), KEEP8
     &        )
          IF (IFLAG.LT.0) GOTO 400
#if defined(BLR_MT)          
!$OMP BARRIER
!$OMP MASTER
#endif          
          CALL SYSTEM_CLOCK(T2,COUNT_RATE)
          LOC_DEMOTING_TIME = LOC_DEMOTING_TIME + 
     &         DBLE(T2-T1)/DBLE(COUNT_RATE)
          CALL STATS_STORE_BLR_PANEL_MRY(BLR_L,
     &               NB_BLR-CURRENT_BLR-NPARTSCB,
     &               NPARTSCB, 'V', 2)
          CALL SYSTEM_CLOCK(T1)
#if defined(BLR_MT)          
!$OMP END MASTER
#endif          
 400      CONTINUE         
#if defined(BLR_MT)          
!$OMP END PARALLEL
#endif          
          IF (IFLAG.LT.0) GOTO 490
        ENDIF
 101    CONTINUE       
        IF (K263.NE.0) THEN
          NELIM = IEND_BLR-NPIV
          BLR_SEND=>BLR_DUMMY
          IF (associated(BLR_L)) THEN
            BLR_SEND=>BLR_L
          ENDIF
          CALL DMUMPS_SEND_FACTORED_BLK( COMM_LOAD, ASS_IRECV,
     &             N, INODE, FPERE, IW, LIW, 
     &             IOLDPS, POSELT, A, LA, LDAFS,
     &             IBEG_BLR, NPIV, IPIV, NASS,LASTBL, NB_BLOC_FAC,
     &
     &             COMM, MYID, BUFR, LBUFR, LBUFR_BYTES,NBFIN,LEAF,
     &             IFLAG, IERROR, IPOOL,LPOOL,
     &             SLAVEF, POSFAC, IWPOS, IWPOSCB, IPTRLU, LRLU,
     &             LRLUS, COMP,
     &             PTRIST, PTRAST, PTLUST_S, PTRFAC, STEP,
     &             PIMASTER, PAMASTER,
     &             NSTK_S,NBPROCFILS,PROCNODE_STEPS, root,
     &             OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &             FILS, PTRARW, PTRAIW,
     &             INTARR,DBLARR,ICNTL,KEEP,KEEP8,DKEEP,ND,FRERE,
     &             LPTRAR, NELT, FRTPTR, FRTELT, 
     &             ISTEP_TO_INIV2, TAB_POS_IN_PERE
     &             , NELIM, LR_ACTIVATED
     &             , NPARTSASS, CURRENT_BLR , BLR_SEND , LRGROUPS
     &             )
          IF ( IFLAG .LT. 0 ) GOTO 500
          IF (RESET_TO_ONE.AND.K109_SAVE.LT.KEEP(109)) THEN
              CALL DMUMPS_RESET_TO_ONE( 
     &        IW(IOLDPS+KEEP(IXSZ)+IW(IOLDPS+5+KEEP(IXSZ))+6),
     &        NPIV, IBEG_BLR,
     &        K109_SAVE, KEEP(109), PIVNUL_LIST, LPN_LIST,
     &        A, POSELT, LA, LDAFS)
          ENDIF
          IF (KEEP(201).EQ.1.AND.(PIVOT_OPTION.GE.3)) THEN
              MonBloc%Last  = .FALSE.
              MonBloc%LastPiv= NPIV
              LAST_CALL=.FALSE.
              CALL DMUMPS_OOC_IO_LU_PANEL(
     &        STRAT_TRY_WRITE,
     &        TYPEFile, A(POSELT),
     &        LAFAC, MonBloc, NextPiv2beWritten, IDUMMY, IW(IOLDPS),
     &        LIWFAC, MYID, KEEP8(31), IFLAG_OOC,LAST_CALL )
              IF (IFLAG_OOC .LT. 0 ) THEN
                IFLAG = IFLAG_OOC
                RETURN
              ENDIF
          ENDIF
        ENDIF
        IF (.NOT. LR_ACTIVATED) THEN
          CALL DMUMPS_FAC_SQ_LDLT(IBEG_BLR,IEND_BLR,NPIV,
     &             NASS,NASS,NASS,INODE,A,LA,
     &             LDAFS, POSELT,
     &             KEEP,KEEP8,PIVOT_OPTION, .TRUE.)
        ELSE
          NELIM = IEND_BLOCK - NPIV
          IF (IEND_BLR.NE.IEND_BLOCK) CALL MUMPS_ABORT()
#if defined(BLR_MT)          
!$OMP PARALLEL 
#endif
          IF (NELIM .EQ. IEND_BLR - IBEG_BLR + 1) GOTO 450
            CALL DMUMPS_BLR_UPDATE_TRAILING_LDLT(A, LA, POSELT, 
     &        IFLAG, IERROR, NASS,
     &        BEGS_BLR, NB_BLR, CURRENT_BLR, BLR_L, NELIM,
     &        IW(HF+IOLDPS+NFRONT+IBEG_BLR-1), BLOCK,
     &        MAXI_CLUSTER, NPIV,
     &        2, 
     &        KEEP(481), DKEEP(8), KEEP(477) 
     &        )
            IF (IFLAG.LT.0) GOTO 450
 450      CONTINUE
#if defined(BLR_MT)          
!$OMP END PARALLEL
#endif          
          IF (IFLAG.LT.0) GOTO 490
          IF (NELIM .EQ. IEND_BLR - IBEG_BLR + 1) GOTO 100
          CALL SYSTEM_CLOCK(T2,COUNT_RATE)
          LOC_UPDT_TIME = LOC_UPDT_TIME + 
     &               DBLE(T2-T1)/DBLE(COUNT_RATE)
          IF (PIVOT_OPTION.LE.2) THEN
             CALL SYSTEM_CLOCK(T1)
             CALL DMUMPS_DECOMPRESS_PANEL(A, LA, POSELT, NASS,
     &              .TRUE.,  
     &       BEGS_BLR(CURRENT_BLR),
     &       BEGS_BLR(CURRENT_BLR+1), NB_BLR, BLR_L, CURRENT_BLR, 'V',
     &       NASS, KEEP(470))
             CALL SYSTEM_CLOCK(T2,COUNT_RATE)
            LOC_PROMOTING_TIME = LOC_PROMOTING_TIME + 
     &                dble(T2-T1)/dble(COUNT_RATE)
          ELSE
            IF (KEEP(485).NE.0) THEN
            CALL SYSTEM_CLOCK(T1)
            CALL DMUMPS_DECOMPRESS_PANEL(A, LA, POSELT, NASS,
     &         .FALSE.,   
     &         BEGS_BLR(CURRENT_BLR),
     &         BEGS_BLR(CURRENT_BLR+1), NB_BLR, BLR_L, CURRENT_BLR, 'V',
     &         NASS, KEEP(470))
            CALL SYSTEM_CLOCK(T2,COUNT_RATE)
            LOC_PROMOTING_TIME = LOC_PROMOTING_TIME + 
     &                DBLE(T2-T1)/DBLE(COUNT_RATE)
            END IF
          ENDIF
            CALL DEALLOC_BLR_PANEL (BLR_L, NB_BLR-CURRENT_BLR, KEEP8,
     &                             .TRUE.)
            DEALLOCATE(BLR_L)
          NULLIFY(BLR_L)
        ENDIF
        IF (KEEP(201).EQ.1.AND.(PIVOT_OPTION.GE.3)) THEN
          MonBloc%Last   = .FALSE.
          MonBloc%LastPiv= NPIV
          LAST_CALL=.FALSE.
          CALL DMUMPS_OOC_IO_LU_PANEL(
     &        STRAT_TRY_WRITE,
     &        TYPEFile, A(POSELT),
     &        LAFAC, MonBloc, NextPiv2beWritten, IDUMMY, IW(IOLDPS),
     &        LIWFAC, MYID, KEEP8(31), IFLAG_OOC,LAST_CALL )
          IF (IFLAG_OOC < 0 ) THEN
              IFLAG = IFLAG_OOC
              GOTO 490
          ENDIF
        ENDIF
  100   CONTINUE
      END DO 
      IF (KEEP(201).EQ.1) THEN 
          STRAT        = STRAT_WRITE_MAX   
          MonBloc%Last = .TRUE.
          MonBloc%LastPiv  = IW(IOLDPS+1+XSIZE)
          LAST_CALL    = .TRUE.
          CALL DMUMPS_OOC_IO_LU_PANEL
     &          ( STRAT, TYPEFile, 
     &           A(POSELT), LAFAC, MonBloc,
     &           NextPiv2beWritten, IDUMMY,
     &           IW(IOLDPS), LIWFAC, 
     &           MYID, KEEP8(31), IFLAG_OOC, LAST_CALL )
          IF (IFLAG_OOC .LT. 0 ) THEN
             IFLAG = IFLAG_OOC
             RETURN
          ENDIF
          CALL DMUMPS_OOC_PP_TRYRELEASE_SPACE (IWPOS, 
     &      IOLDPS, IW, LIW, MonBloc , NFRONT, KEEP)
      ENDIF
      GOTO 500
 480  CONTINUE
        write(*,*) 'Allocation problem in BLR routine
     &     DMUMPS_FAC_FRONT_LDLT_TYPE2: ',
     &    'not enough memory? memory requested = ' , IERROR
 490  CONTINUE
      CALL DMUMPS_BDC_ERROR( MYID, SLAVEF, COMM, KEEP )
 500  CONTINUE
      IF(allocated(IPIV)) DEALLOCATE( IPIV )
      IF (allocated(DIAG_ORIG)) DEALLOCATE(DIAG_ORIG)
      IF (LR_ACTIVATED) THEN
         CALL STATS_COMPUTE_MRY_FRONT_TYPE2(NASS, NFRONT, 1, INODE,
     &                                      NELIM)
         CALL STATS_COMPUTE_FLOP_FRONT_TYPE2(NFRONT, NASS, KEEP(50),
     &                      INODE, NELIM)
         CALL SYSTEM_CLOCK(TTOT2,COUNT_RATETOT)
         LOC_LR_MODULE_TIME = DBLE(TTOT2-TTOT1)/DBLE(COUNT_RATETOT)
         IF (allocated(RWORK))  DEALLOCATE(RWORK)
         IF (allocated(WORK))  deallocate(WORK)
         IF (allocated(TAU))   deallocate(TAU)
         IF (allocated(JPVT))  deallocate(JPVT)
         IF (allocated(BLOCK)) deallocate(BLOCK)
         IF (associated(BEGS_BLR)) THEN
           DEALLOCATE(BEGS_BLR)
           NULLIFY(BEGS_BLR)
         ENDIF
      ELSE
          CALL SYSTEM_CLOCK(TTOT2FR,COUNT_RATETOTFR)
          LOC_FRFRONTS_TIME =
     &       DBLE(TTOT2FR-TTOT1FR)/DBLE(COUNT_RATETOTFR)
          CALL UPDATE_FLOP_STATS_FRFRONTS(NFRONT, NPIV, NASS, KEEP(50), 
     &                                       2)
      ENDIF
      CALL UPDATE_ALL_TIMES(INODE,LOC_UPDT_TIME,LOC_PROMOTING_TIME,
     &                        LOC_DEMOTING_TIME, LOC_CB_DEMOTING_TIME,
     &                        LOC_FRPANELS_TIME, LOC_FRFRONTS_TIME,
     &                        LOC_TRSM_TIME, LOC_LR_MODULE_TIME,
     &                        LOC_FAC_I_TIME, LOC_FAC_MQ_TIME,
     &                        LOC_FAC_SQ_TIME)
      RETURN
      END SUBROUTINE DMUMPS_FAC2_LDLT
      SUBROUTINE DMUMPS_RESET_TO_ONE(FRONT_INDEX_LIST, NPIV,
     & IBEG_BLOCK, K109_SAVE, K109, PIVNUL_LIST, LPN_LIST,
     & A, POSELT, LA, LDAFS)
      INTEGER, INTENT(IN) :: NPIV, IBEG_BLOCK
      INTEGER, INTENT(IN) :: FRONT_INDEX_LIST(NPIV)
      INTEGER, INTENT(IN) :: K109
      INTEGER, INTENT(INOUT) :: K109_SAVE
      INTEGER, INTENT(IN) :: LPN_LIST
      INTEGER, INTENT(IN) :: PIVNUL_LIST(LPN_LIST)
      INTEGER(8), INTENT(IN) :: POSELT, LA
      INTEGER, INTENT(IN) :: LDAFS
      DOUBLE PRECISION, INTENT(INOUT) :: A(LA)
      LOGICAL :: TO_UPDATE
      INTEGER :: I, JJ, K
      DOUBLE PRECISION ONE
      PARAMETER (ONE = 1.0D0)
      DO K = K109_SAVE+1, K109
        TO_UPDATE = .FALSE. 
        I = PIVNUL_LIST(K)  
        DO JJ=IBEG_BLOCK, NPIV
          IF (FRONT_INDEX_LIST(JJ) .EQ.I) THEN
            TO_UPDATE=.TRUE. 
            EXIT
          ENDIF
        ENDDO
        IF (TO_UPDATE) THEN
          A(POSELT+int(JJ,8)+int(LDAFS,8)*int(JJ-1,8))= ONE
          TO_UPDATE=.FALSE. 
        ELSE
          write(*,*) ' Internal error related ', 
     &               'to null pivot row detection'
          CALL MUMPS_ABORT()
        ENDIF
      ENDDO
      K109_SAVE = K109
      RETURN
      END SUBROUTINE DMUMPS_RESET_TO_ONE
      END MODULE DMUMPS_FAC2_LDLT_M
