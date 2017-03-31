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
      MODULE DMUMPS_FAC2_LU_M
      CONTAINS
      SUBROUTINE DMUMPS_FAC2_LU( COMM_LOAD, ASS_IRECV, 
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
     &               , LRGROUPS
     &             )
!$    USE OMP_LIB
      USE DMUMPS_FAC_FRONT_AUX_M
      USE DMUMPS_FAC_FRONT_TYPE2_AUX_M
      USE DMUMPS_OOC
      USE DMUMPS_BUF, ONLY : DMUMPS_BUF_TEST
      USE DMUMPS_FAC_LR
      USE DMUMPS_LR_CORE
      USE DMUMPS_LR_TYPE
      USE DMUMPS_LR_STATS
      USE DMUMPS_ANA_LR
!$    USE OMP_LIB
      IMPLICIT NONE
      INCLUDE 'dmumps_root.h'
      INTEGER COMM_LOAD, ASS_IRECV
      INTEGER N, INODE, FPERE, LIW, NOFFW, NPVW
      INTEGER(8) :: LA
      INTEGER IW( LIW )
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
      INTEGER INOPV, IFINB, NFRONT, NPIV, IBEG_BLOCK, IEND_BLOCK
      INTEGER :: IBEG_BLOCK_FOR_IPIV
      INTEGER NASS, NBKJIB_ORIG, XSIZE
      INTEGER NBLR_ORIG, IBEG_BLR, IEND_BLR
      INTEGER Inextpiv
      LOGICAL LASTBL 
      INTEGER(8) :: POSELT
      INTEGER IOLDPS, allocok, K263
      INTEGER idummy 
      DOUBLE PRECISION    UUTEMP
      LOGICAL STATICMODE
      DOUBLE PRECISION SEUIL_LOC
      INTEGER , ALLOCATABLE, DIMENSION ( : ) :: IPIV
      INTEGER(8) :: LAFAC
      INTEGER LIWFAC, STRAT, LNextPiv2beWritten, 
     &        UNextPiv2beWritten, IFLAG_OOC,
     &        PP_FIRST2SWAP_L, PP_FIRST2SWAP_U,
     &        PP_LastPIVRPTRFilled_L,
     &        PP_LastPIVRPTRFilled_U
      TYPE(IO_BLOCK) :: MonBloc 
      LOGICAL LAST_CALL
      INTEGER CURRENT_BLR, NELIM
      LOGICAL LR_ACTIVATED
      INTEGER PIVOT_OPTION, LAST_COL
      INTEGER HF, NPARTSASS, NPARTSCB, NB_BLR, END_I
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
      TYPE(LRB_TYPE), POINTER, DIMENSION(:) :: BLR_L, BLR_U, BLR_SEND
      TYPE(LRB_TYPE), DIMENSION(1), TARGET  :: BLR_DUMMY
      DOUBLE PRECISION, ALLOCATABLE :: WORK(:), TAU(:)
      INTEGER, ALLOCATABLE :: JPVT(:)
      DOUBLE PRECISION, ALLOCATABLE :: RWORK(:)
      DOUBLE PRECISION, ALLOCATABLE :: BLOCK(:,:)
      INTEGER :: OMP_NUM
      INTEGER :: NOMP
      INCLUDE 'mumps_headers.h'
      NULLIFY(BLR_L,BLR_U) 
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
      NOMP=1
!$    NOMP=OMP_GET_MAX_THREADS()
      INOPV = 0
      SEUIL_LOC = SEUIL
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
         UUTEMP=UU
      ENDIF
      PIVOT_OPTION = KEEP(468)
      IF (UUTEMP == 0.0D0 .AND. KEEP(201).NE.1) THEN
      ENDIF
      idummy  = 0
      IOLDPS = PTLUST_S(STEP( INODE ))
      POSELT = PTRAST(STEP( INODE ))
      XSIZE  = KEEP(IXSZ)
      NFRONT = IW(IOLDPS+XSIZE)
      NASS   = iabs(IW(IOLDPS+2+XSIZE))
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
        WRITE(*,*) MYID,' : DMUMPS_FAC2_LU :failed to allocate ',
     &  NASS, ' integers'
        IFLAG  = -13
        IERROR =NASS
        GO TO 490
      END IF
      IF (KEEP(201).EQ.1) THEN 
          CALL MUMPS_GETI8(LAFAC,IW(IOLDPS+XXR))
          LIWFAC    = IW(IOLDPS+XXI)
          LNextPiv2beWritten = 1 
          UNextPiv2beWritten = 1 
          PP_FIRST2SWAP_L = LNextPiv2beWritten 
          PP_FIRST2SWAP_U = UNextPiv2beWritten 
          MonBloc%LastPanelWritten_L = 0 
          MonBloc%LastPanelWritten_U = 0        
          MonBloc%INODE    = INODE
          MonBloc%MASTER   = .TRUE.
          MonBloc%Typenode = 2
          MonBloc%NROW     = NASS
          MonBloc%NCOL     = NFRONT
          MonBloc%NFS      = NASS
          MonBloc%Last     = .FALSE.   
          MonBloc%LastPiv  = -68877    
          NULLIFY(MonBloc%INDICES)
      ENDIF
      IF (LR_ACTIVATED) THEN
             PIVOT_OPTION = 4
             CNT_NODES = CNT_NODES + 1 
             CALL SYSTEM_CLOCK(TTOT1)
      ELSE IF (KEEP(486).GT.0) THEN
             CALL SYSTEM_CLOCK(TTOT1FR)
      ENDIF
      HF = 6 + IW(IOLDPS+5+XSIZE)+XSIZE
      IF (KEEP(201).EQ.1) THEN  
       IF (PIVOT_OPTION.LT.3) PIVOT_OPTION=3
      ENDIF
      IF (LR_ACTIVATED) THEN
         CALL GET_CUT(IW(IOLDPS+HF:IOLDPS+HF+NFRONT-1), NASS,
     &        NFRONT-NASS, LRGROUPS, NPARTSCB, 
     &        NPARTSASS, BEGS_BLR)
         CALL REGROUPING2(BEGS_BLR, NPARTSASS, NASS, NPARTSCB,
     &        NFRONT-NASS, KEEP(488), .FALSE., KEEP(472))     
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
              GOTO 490
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
            CALL DMUMPS_FAC_I(NFRONT,NASS,NASS,
     &      IBEG_BLOCK_FOR_IPIV,IEND_BLOCK,N,INODE,
     &      IW,LIW,A,LA,INOPV,NOFFW,
     &      IFLAG,IOLDPS,POSELT,UU,SEUIL_LOC,KEEP,KEEP8,
     &      DKEEP(1),PIVNUL_LIST(1),LPN_LIST,
     &      PP_FIRST2SWAP_L,  MonBloc%LastPanelWritten_L,
     &      PP_LastPIVRPTRFilled_L,
     &      PP_FIRST2SWAP_U,  MonBloc%LastPanelWritten_U,
     &      PP_LastPIVRPTRFilled_U,
     &      PIVOT_OPTION, IEND_BLR, 
     &      Inextpiv,
     &      IPIV 
     &      )
            IF (LR_ACTIVATED) THEN
                CALL SYSTEM_CLOCK(T2P,CRP)
                LOC_FAC_I_TIME = LOC_FAC_I_TIME + 
     &                           dble(T2P-T1P)/dble(CRP)
            ENDIF
            IF (IFLAG.LT.0) GOTO 490   
          IF (INOPV.EQ.1) THEN
              IF (STATICMODE) THEN
                INOPV = -1
                GOTO 50
              ENDIF
              LASTBL = .TRUE.
          ELSE IF (INOPV .LE. 0) THEN 
            IF (LR_ACTIVATED) THEN
                CALL SYSTEM_CLOCK(T1P)
            ENDIF
            IF (PIVOT_OPTION.GE.3) THEN
              LAST_COL = NFRONT
            ELSEIF (PIVOT_OPTION.EQ.2) THEN
              LAST_COL = NASS
            ELSE
              LAST_COL = IEND_BLR
            ENDIF
            CALL DMUMPS_FAC_MQ(IBEG_BLOCK, IEND_BLOCK,
     &             NFRONT, NASS, IW(IOLDPS+1+XSIZE),
     &             LAST_COL, A, LA, POSELT, IFINB)
            IF (LR_ACTIVATED) THEN
                CALL SYSTEM_CLOCK(T2P,CRP)
                LOC_FAC_MQ_TIME = LOC_FAC_MQ_TIME + 
     &                            dble(T2P-T1P)/dble(CRP)
            ENDIF
            IW(IOLDPS+1+XSIZE) = IW(IOLDPS+1+XSIZE) + 1
            NPVW = NPVW + 1
            IF (IFINB.EQ.0) THEN
              GOTO 50 
            ELSE IF (IFINB .EQ. -1) THEN
              LASTBL = .TRUE.
            ENDIF
          ENDIF
          NPIV = IW(IOLDPS+1+XSIZE)
          IF (K263.EQ.0) THEN
            NELIM = IEND_BLR - NPIV
            CALL DMUMPS_SEND_FACTORED_BLK( COMM_LOAD, ASS_IRECV, 
     &             N, INODE, FPERE, IW, LIW, 
     &             IOLDPS, POSELT, A, LA, NFRONT,
     &             IBEG_BLOCK, NPIV, IPIV, NASS,LASTBL, idummy, 
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
     &             , NPARTSASS, CURRENT_BLR
     &             , BLR_DUMMY, LRGROUPS
     & )
          END IF
          IF ( IFLAG .LT. 0 ) GOTO 500
          IF ( KEEP(201).EQ.1.AND.(PIVOT_OPTION.GE.3) 
     &        .AND.
     &          ( .NOT. LR_ACTIVATED .OR.
     &            (KEEP(485).EQ.0) 
     &          ) 
     &       ) THEN  
            MonBloc%LastPiv= IW(IOLDPS+1+XSIZE)
            STRAT          = STRAT_TRY_WRITE
            LAST_CALL      = .FALSE.
            CALL DMUMPS_OOC_IO_LU_PANEL
     &          ( STRAT, TYPEF_U,
     &           A(POSELT), LAFAC, MonBloc,
     &           LNextPiv2beWritten, UNextPiv2beWritten,
     &           IW(IOLDPS), LIWFAC, 
     &           MYID, KEEP8(31), IFLAG_OOC,LAST_CALL )
            IF (IFLAG_OOC < 0 ) THEN
              IFLAG=IFLAG_OOC
              GOTO 490
            ENDIF
          ENDIF
          CALL DMUMPS_BUF_TEST()
          NPIV       =  IW(IOLDPS+1+XSIZE)
          IF ( IEND_BLR .GT. IEND_BLOCK ) THEN
              IF (LR_ACTIVATED) THEN
                  CALL SYSTEM_CLOCK(T1P)
              ENDIF
              IF (PIVOT_OPTION.GE.3) THEN
                LAST_COL = NFRONT
              ELSEIF (PIVOT_OPTION.EQ.2) THEN
                LAST_COL = NASS
              ELSE
                LAST_COL = IEND_BLR
              ENDIF
              CALL DMUMPS_FAC_SQ(IBEG_BLOCK, IEND_BLOCK,
     &            NPIV, NFRONT, IEND_BLR, LAST_COL,
     &            A, LA, POSELT, .FALSE., .TRUE.,
     &            .FALSE. ) 
              IF (LR_ACTIVATED) THEN
                  CALL SYSTEM_CLOCK(T2P,CRP)
                  LOC_FAC_SQ_TIME = LOC_FAC_SQ_TIME +
     &                              dble(T2P-T1P)/dble(CRP)
              ENDIF
          ENDIF
          CALL DMUMPS_BUF_TEST()
        END DO 
        NPIV   = IW(IOLDPS+1+XSIZE)
        IF (LR_ACTIVATED) THEN
          CALL SYSTEM_CLOCK(T2,COUNT_RATE)
          LOC_FRPANELS_TIME = LOC_FRPANELS_TIME +
     &                        dble(T2-T1)/dble(COUNT_RATE)
          CALL UPDATE_FLOP_STATS_PANEL(NFRONT - IBEG_BLR + 1, 
     &                                 NPIV   - IBEG_BLR + 1, 2, 0)
        ENDIF
        IF (LR_ACTIVATED) THEN
          NELIM = IEND_BLR - NPIV
          IF (NELIM .EQ. IEND_BLR - IBEG_BLR + 1) THEN
            GOTO 101
          ENDIF
            END_I=NB_BLR
          ALLOCATE(BLR_U(NB_BLR-CURRENT_BLR)) 
          CALL SYSTEM_CLOCK(T1)
#if defined(BLR_MT)          
!$OMP PARALLEL 
#endif
          CALL DMUMPS_COMPRESS_PANEL(A, LA, POSELT, IFLAG, IERROR, 
     &       NFRONT,
     &       BEGS_BLR, NB_BLR, DKEEP(8), KEEP(473), BLR_U, CURRENT_BLR,
     &       'H', WORK, TAU, JPVT, LWORK, RWORK,
     &       BLOCK, MAXI_CLUSTER, NELIM, 
     &       .FALSE., 0, 0, 2, KEEP(483), KEEP(470), KEEP8,
     &       END_I_IN=END_I
     &        )
          IF (IFLAG.LT.0) GOTO 300
#if defined(BLR_MT)          
!$OMP BARRIER
!$OMP MASTER
#endif          
          CALL SYSTEM_CLOCK(T2,COUNT_RATE)
          LOC_DEMOTING_TIME = LOC_DEMOTING_TIME + 
     &              DBLE(T2-T1)/DBLE(COUNT_RATE)
          CALL STATS_STORE_BLR_PANEL_MRY(BLR_U,
     &               NB_BLR-CURRENT_BLR-NPARTSCB,
     &               NPARTSCB, 'H', 2)
            CALL SYSTEM_CLOCK(T1)
#if defined(BLR_MT)          
!$OMP END MASTER
#endif          
 300      CONTINUE         
#if defined(BLR_MT)          
!$OMP END PARALLEL
#endif          
        ENDIF
 101    CONTINUE       
        IF (K263.NE.0) THEN
          NELIM = IEND_BLR - NPIV
          BLR_SEND=>BLR_DUMMY
          IF (associated(BLR_U)) THEN
            BLR_SEND=>BLR_U
          ENDIF
          CALL DMUMPS_SEND_FACTORED_BLK( COMM_LOAD, ASS_IRECV, 
     &             N, INODE, FPERE, IW, LIW, 
     &             IOLDPS, POSELT, A, LA, NFRONT,
     &             IBEG_BLR, NPIV, IPIV, NASS,LASTBL, idummy, 
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
     &             , NPARTSASS, CURRENT_BLR
     &             , BLR_SEND, LRGROUPS
     &             )
        ENDIF
        IF (.NOT. LR_ACTIVATED) THEN
          IF (PIVOT_OPTION.GE.3) THEN
            LAST_COL = NFRONT
          ELSE
            LAST_COL = NASS
          ENDIF
          IF (IEND_BLR .LT. NASS) THEN
              CALL DMUMPS_FAC_SQ(IBEG_BLR, IEND_BLR,
     &            NPIV, NFRONT, NASS, LAST_COL,
     &            A, LA, POSELT, (PIVOT_OPTION.LT.2), .TRUE.
     &            , (KEEP(377) .EQ. 1) 
     &            )
          ENDIF
        ELSE
         NELIM = IEND_BLR - NPIV
         IF (NELIM .EQ. IEND_BLR - IBEG_BLR + 1) THEN
            GOTO 100
         ENDIF
           CALL SYSTEM_CLOCK(T1)
           IF (IEND_BLR.LT.NFRONT) THEN
              CALL DMUMPS_FAC_SQ(IBEG_BLR, IEND_BLR,
     &            NPIV, NFRONT, NASS, 
     &            -77777, 
     &            A, LA, POSELT, .FALSE., .FALSE.,
     &            .FALSE. ) 
           ENDIF
           CALL SYSTEM_CLOCK(T2,COUNT_RATE)
           LOC_TRSM_TIME = LOC_TRSM_TIME + 
     &              DBLE(T2-T1)/DBLE(COUNT_RATE)
         CALL SYSTEM_CLOCK(T1)
         ALLOCATE(BLR_L(NPARTSASS-CURRENT_BLR)) 
#if defined(BLR_MT)          
!$OMP PARALLEL 
#endif
         CALL DMUMPS_COMPRESS_PANEL(A, LA, POSELT, IFLAG, IERROR, 
     &        NFRONT,
     &        BEGS_BLR, NPARTSASS, DKEEP(8), KEEP(473), BLR_L,
     &        CURRENT_BLR, 'V', WORK, TAU, JPVT, LWORK, RWORK,
     &        BLOCK, MAXI_CLUSTER, NELIM,
     &        .FALSE., 0, 0,
     &        2, KEEP(483), KEEP(470), KEEP8 
     &        )
         IF (IFLAG.LT.0) GOTO 400
#if defined(BLR_MT)          
!$OMP BARRIER
!$OMP MASTER
#endif          
          CALL SYSTEM_CLOCK(T2,COUNT_RATE)
          LOC_DEMOTING_TIME = LOC_DEMOTING_TIME + 
     &              DBLE(T2-T1)/DBLE(COUNT_RATE)
          CALL SYSTEM_CLOCK(T1)
#if defined(BLR_MT)          
!$OMP END MASTER
#endif          
            CALL DMUMPS_BLR_UPDATE_TRAILING(A, LA, POSELT, 
     &        IFLAG, IERROR, NFRONT,
     &        BEGS_BLR, BEGS_BLR, CURRENT_BLR, BLR_L, NPARTSASS,
     &        BLR_U, NB_BLR, NELIM, .FALSE., 0,
     &        2, 0, KEEP(470), 
     &        KEEP(481), DKEEP(8), KEEP(477) 
     &        )
 400      CONTINUE
#if defined(BLR_MT)          
!$OMP END PARALLEL
#endif          
          IF (IFLAG.LT.0) GOTO 490
          CALL SYSTEM_CLOCK(T2,COUNT_RATE)
          LOC_UPDT_TIME = LOC_UPDT_TIME + 
     &               DBLE(T2-T1)/DBLE(COUNT_RATE)
          CALL STATS_STORE_BLR_PANEL_MRY(BLR_L,
     &               NB_BLR-CURRENT_BLR-NPARTSCB,
     &               0, 'V', 2)
           IF (KEEP(485).NE.0) THEN
             CALL SYSTEM_CLOCK(T1)
             CALL DMUMPS_DECOMPRESS_PANEL(A, LA, POSELT, NFRONT,
     &             .FALSE., 
     &        BEGS_BLR(CURRENT_BLR), BEGS_BLR(CURRENT_BLR+1), 
     &        NPARTSASS, BLR_L, CURRENT_BLR, 'V', NFRONT, KEEP(470))
             CALL SYSTEM_CLOCK(T2,COUNT_RATE)
             LOC_PROMOTING_TIME = LOC_PROMOTING_TIME + 
     &               DBLE(T2-T1)/DBLE(COUNT_RATE)
           ENDIF 
           IF (KEEP(485).NE.0) THEN
             CALL SYSTEM_CLOCK(T1)
             CALL DMUMPS_DECOMPRESS_PANEL(A, LA, POSELT, NFRONT,
     &        .FALSE.,   
     &         BEGS_BLR(CURRENT_BLR),
     &         BEGS_BLR(CURRENT_BLR+1), NB_BLR, BLR_U, CURRENT_BLR, 'H',
     &         NFRONT, KEEP(470),
     &         END_I_IN=END_I
     &         )
             CALL SYSTEM_CLOCK(T2,COUNT_RATE)
             LOC_PROMOTING_TIME = LOC_PROMOTING_TIME + 
     &               DBLE(T2-T1)/DBLE(COUNT_RATE)
           ENDIF
            CALL DEALLOC_BLR_PANEL (BLR_U, NB_BLR-CURRENT_BLR, KEEP8,
     &                             .TRUE.)
            CALL DEALLOC_BLR_PANEL (BLR_L, NPARTSASS-CURRENT_BLR, KEEP8,
     &                             .TRUE.)
            DEALLOCATE(BLR_U,BLR_L)
          NULLIFY(BLR_L)
          NULLIFY(BLR_U)
        ENDIF
        IF (KEEP(201).EQ.1) THEN 
          STRAT            = STRAT_TRY_WRITE
          MonBloc%LastPiv  = NPIV
          LAST_CALL= .FALSE.
          CALL DMUMPS_OOC_IO_LU_PANEL
     &          ( STRAT, TYPEF_BOTH_LU, 
     &           A(POSELT), LAFAC, MonBloc,
     &           LNextPiv2beWritten, UNextPiv2beWritten,
     &           IW(IOLDPS), LIWFAC, 
     &           MYID, KEEP8(31), IFLAG_OOC,LAST_CALL )
             IF (IFLAG_OOC < 0 ) THEN
                IFLAG=IFLAG_OOC
                GOTO 490
             ENDIF
        ENDIF
 100    CONTINUE
      END DO 
      IF (LR_ACTIVATED) THEN
         CALL STATS_COMPUTE_MRY_FRONT_TYPE2(NASS, NFRONT, 0, INODE,
     &                                      NELIM)
         CALL SYSTEM_CLOCK(TTOT2,COUNT_RATETOT)
         CALL STATS_COMPUTE_FLOP_FRONT_TYPE2(NFRONT, NASS, KEEP(50),
     &                      STEP_STATS(INODE), NELIM )
         LOC_LR_MODULE_TIME = DBLE(TTOT2-TTOT1)/DBLE(COUNT_RATETOT)
        IF (allocated(RWORK))  DEALLOCATE(RWORK)
        IF (allocated(work)) DEALLOCATE(WORK)
        IF (allocated(TAU)) DEALLOCATE(TAU)
        IF (allocated(JPVT)) DEALLOCATE(JPVT)
        IF (allocated(BLOCK)) DEALLOCATE(BLOCK)
         IF (associated(BEGS_BLR)) THEN
           DEALLOCATE(BEGS_BLR)
           NULLIFY(BEGS_BLR)
         ENDIF
      ELSE 
          CALL SYSTEM_CLOCK(TTOT2FR,COUNT_RATETOTFR)
          LOC_FRFRONTS_TIME =
     &       DBLE(TTOT2FR-TTOT1FR)/DBLE(COUNT_RATETOTFR)
          CALL UPDATE_FLOP_STATS_FRFRONTS(NFRONT, NPIV, NASS, 0, 2)
      ENDIF
      CALL UPDATE_ALL_TIMES(INODE,LOC_UPDT_TIME,LOC_PROMOTING_TIME,
     &                        LOC_DEMOTING_TIME, LOC_CB_DEMOTING_TIME,
     &                        LOC_FRPANELS_TIME, LOC_FRFRONTS_TIME,
     &                        LOC_TRSM_TIME, LOC_LR_MODULE_TIME,
     &                        LOC_FAC_I_TIME, LOC_FAC_MQ_TIME,
     &                        LOC_FAC_SQ_TIME)
      IF (KEEP(201).EQ.1) THEN 
          STRAT        = STRAT_WRITE_MAX   
          MonBloc%Last = .TRUE.
          MonBloc%LastPiv  = IW(IOLDPS+1+XSIZE)
          LAST_CALL    = .TRUE.
          CALL DMUMPS_OOC_IO_LU_PANEL
     &          ( STRAT, TYPEF_BOTH_LU,
     &           A(POSELT), LAFAC, MonBloc,
     &           LNextPiv2beWritten, UNextPiv2beWritten,
     &           IW(IOLDPS), LIWFAC, 
     &           MYID, KEEP8(31), IFLAG_OOC, LAST_CALL )
          IF (IFLAG_OOC < 0 ) THEN
            IFLAG=IFLAG_OOC
            GOTO 490
          ENDIF
          CALL DMUMPS_OOC_PP_TRYRELEASE_SPACE (IWPOS, 
     &      IOLDPS, IW, LIW, MonBloc , NFRONT, KEEP)
      ENDIF
      GOTO 500
 480  CONTINUE
      write(*,*) 'Allocation problem in BLR routine
     &     DMUMPS_FAC_FRONT_LU_TYPE2: ',
     &     'not enough memory? memory requested = ' , IERROR
 490  CONTINUE
      CALL DMUMPS_BDC_ERROR( MYID, SLAVEF, COMM, KEEP )
 500  CONTINUE
      DEALLOCATE( IPIV )
      RETURN
      END SUBROUTINE DMUMPS_FAC2_LU
      END MODULE DMUMPS_FAC2_LU_M
