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
      MODULE DMUMPS_FAC1_LDLT_M
      CONTAINS
      SUBROUTINE DMUMPS_FAC1_LDLT( N, INODE, IW, LIW, A, LA,
     &                           IOLDPS, POSELT, IFLAG, IERROR,
     &                           UU, NNEG, NPVW,
     &                           KEEP,KEEP8,
     &                           MYID, SEUIL, AVOID_DELAYED, ETATASS,
     &     DKEEP,PIVNUL_LIST,LPN_LIST, IWPOS
     &     , LRGROUPS
     &     )
      USE DMUMPS_FAC_FRONT_AUX_M
      USE DMUMPS_OOC
      USE DMUMPS_FAC_LR  
      USE DMUMPS_LR_TYPE
      USE DMUMPS_LR_STATS
      USE DMUMPS_ANA_LR
!$    USE OMP_LIB
      IMPLICIT NONE
      INTEGER(8) :: LA, POSELT
      INTEGER N, INODE, LIW, IFLAG, IERROR, NNEG, NPVW
      INTEGER MYID, IOLDPS
      INTEGER KEEP( 500 )
      INTEGER(8) KEEP8(150)
      DOUBLE PRECISION UU, SEUIL
      DOUBLE PRECISION A( LA )
      INTEGER, TARGET :: IW( LIW )
      LOGICAL AVOID_DELAYED
      INTEGER ETATASS, IWPOS
      INTEGER LPN_LIST
      INTEGER PIVNUL_LIST(LPN_LIST)
      DOUBLE PRECISION DKEEP(230)
      INTEGER :: LRGROUPS(N)
      INTEGER INOPV, IFINB, NFRONT, NPIV, IBEG_BLOCK, IEND_BLOCK
      INTEGER NASS, NBKJIB_ORIG, XSIZE
      INTEGER :: LDA
      DOUBLE PRECISION UUTEMP
      LOGICAL STATICMODE
      DOUBLE PRECISION SEUIL_LOC
      LOGICAL IS_MAXFROMM_AVAIL
      INTEGER PIVOT_OPTION
      DOUBLE PRECISION MAXFROMM
      INTEGER(8) :: LAFAC
      INTEGER LIWFAC, STRAT, NextPiv2beWritten, IFLAG_OOC,
     &        IDUMMY, PP_FIRST2SWAP_L, PP_LastPIVRPTRFilled
      TYPE(IO_BLOCK) :: MonBloc 
      LOGICAL LAST_CALL
      LOGICAL LASTBL
      INTEGER CURRENT_BLR
      LOGICAL LR_ACTIVATED 
      LOGICAL COMPRESS_CB, COMPRESS_PANEL
      INTEGER HF, NPARTSASS, NPARTSCB, NB_BLR
      INTEGER MAXI_CLUSTER, LWORK, NELIM
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
      TYPE(LRB_TYPE), POINTER, DIMENSION(:) :: BLR_L
      DOUBLE PRECISION, ALLOCATABLE :: WORK(:), TAU(:)
      INTEGER, ALLOCATABLE :: JPVT(:)
      DOUBLE PRECISION,ALLOCATABLE :: RWORK(:)
      DOUBLE PRECISION, ALLOCATABLE :: BLOCK(:,:)
      INTEGER :: allocok
      INTEGER :: OMP_NUM
      DOUBLE PRECISION :: ONE, MONE, ZERO
      PARAMETER (ONE = 1.0D0, MONE=-1.0D0)
      PARAMETER (ZERO=0.0D0)
      INCLUDE 'mumps_headers.h'
      INTEGER NBLR_ORIG, IBEG_BLR, IEND_BLR
      INTEGER Inextpiv
      INTEGER PIVSIZ,IWPOSP2
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
      IS_MAXFROMM_AVAIL = .FALSE.
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
        UUTEMP=UU
        SEUIL_LOC = SEUIL
      ENDIF
      PIVOT_OPTION = KEEP(468)
      LAFAC  = -9999_8  
      XSIZE  = KEEP(IXSZ)
      NFRONT = IW(IOLDPS+XSIZE)
      LDA    = NFRONT
      NASS   = iabs(IW(IOLDPS+2+XSIZE))
      IW(IOLDPS+3+XSIZE) =  -99999
      LR_ACTIVATED= .FALSE.        
      COMPRESS_PANEL = .FALSE.
      COMPRESS_CB = .FALSE.
      NULLIFY(BEGS_BLR)
      COMPRESS_PANEL = (IW(IOLDPS+XXLR).GE.2)
      COMPRESS_CB    = ((IW(IOLDPS+XXLR).EQ.1).OR.
     &                  (IW(IOLDPS+XXLR).EQ.3))
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
      IEND_BLOCK  = 0
      IEND_BLR    = 0
      CURRENT_BLR = 0
      LASTBL      = .FALSE.
      IF (KEEP(201).EQ.1) THEN 
          IDUMMY    = -8765
          CALL MUMPS_GETI8(LAFAC,IW(IOLDPS+XXR))
          LIWFAC    = IW(IOLDPS+XXI)
          NextPiv2beWritten = 1 
          PP_FIRST2SWAP_L = NextPiv2beWritten 
          MonBloc%LastPanelWritten_L = 0 
          PP_LastPIVRPTRFilled       = 0
          MonBloc%INODE    = INODE
          MonBloc%MASTER   = .TRUE.
          MonBloc%Typenode = 1
          MonBloc%NROW     = NFRONT
          MonBloc%NCOL     = NFRONT
          MonBloc%NFS      = NASS
          MonBloc%Last     = .FALSE.   
          MonBloc%LastPiv  = -77777    
          MonBloc%INDICES  =>
     &              IW(IOLDPS+6+NFRONT+XSIZE:
     &                 IOLDPS+5+NFRONT+XSIZE+NFRONT)
      ENDIF
      IF (LR_ACTIVATED) THEN
             PIVOT_OPTION = 3
             CNT_NODES = CNT_NODES + 1 
             CALL INIT_STATS_FRONT(NFRONT, STEP_STATS(INODE), NASS,
     &                       NFRONT-NASS)
             CALL SYSTEM_CLOCK(TTOT1)
      ELSE IF (KEEP(486).GT.0) THEN
            CALL INIT_STATS_FRONT(-NFRONT, STEP_STATS(INODE), NASS, 
     &           NFRONT-NASS)
             CALL SYSTEM_CLOCK(TTOT1FR)
      ENDIF
      IF (KEEP(201).EQ.1) THEN   
       IF (PIVOT_OPTION.LT.3) PIVOT_OPTION=3
      ENDIF
      HF = 6 + IW(IOLDPS+5+XSIZE)+XSIZE
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
           GOTO 490
         ENDIF
      ENDIF
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
            IF (LR_ACTIVATED) THEN
                CALL SYSTEM_CLOCK(T1P)
            ENDIF
            CALL DMUMPS_FAC_I_LDLT(NFRONT,NASS,INODE,
     &                IBEG_BLOCK, IEND_BLOCK,
     &                IW,LIW,A,LA,
     &                INOPV, NNEG, IFLAG,IOLDPS,POSELT,UUTEMP,
     &                SEUIL_LOC,KEEP,KEEP8,PIVSIZ,
     &      DKEEP(1),PIVNUL_LIST(1),LPN_LIST, XSIZE,
     &      PP_FIRST2SWAP_L, MonBloc%LastPanelWritten_L,
     &      PP_LastPIVRPTRFilled, MAXFROMM, IS_MAXFROMM_AVAIL,
     &      PIVOT_OPTION, IEND_BLR, Inextpiv)
            IF (LR_ACTIVATED) THEN
                CALL SYSTEM_CLOCK(T2P,CRP)
                LOC_FAC_I_TIME = LOC_FAC_I_TIME + 
     &                           dble(T2P-T1P)/dble(CRP)
            ENDIF
            IF (IFLAG.LT.0) GOTO 500
            IF(KEEP(109).GT. 0) THEN
              IF(PIVNUL_LIST(KEEP(109)).EQ.-1) THEN
                IWPOSP2 = IOLDPS+IW(IOLDPS+1+XSIZE)+6+XSIZE
     &                  +IW(IOLDPS+5+XSIZE)
                PIVNUL_LIST(KEEP(109)) = IW(IWPOSP2)
              ENDIF
            ENDIF
          IF (INOPV.EQ.1) THEN
            IF(STATICMODE) THEN
              INOPV = -1
              GOTO 50 
            ENDIF
            LASTBL = .TRUE.
          ELSE IF ( INOPV.LE.0 ) THEN 
            NPVW = NPVW + PIVSIZ
            IF (LR_ACTIVATED) THEN
                CALL SYSTEM_CLOCK(T1P)
            ENDIF
            CALL DMUMPS_FAC_MQ_LDLT(IEND_BLOCK,
     &             NFRONT, NASS, IW(IOLDPS+1+XSIZE),
     &             INODE,A,LA,
     &             LDA, 
     &             POSELT,IFINB,
     &             PIVSIZ, MAXFROMM,
     &             IS_MAXFROMM_AVAIL, (UUTEMP.NE.0.0D0),
     &             KEEP(253), 
     &             PIVOT_OPTION, IEND_BLR
     &             )
            IF (LR_ACTIVATED) THEN
                CALL SYSTEM_CLOCK(T2P,CRP)
                LOC_FAC_MQ_TIME = LOC_FAC_MQ_TIME + 
     &                            dble(T2P-T1P)/dble(CRP)
            ENDIF
            IF(PIVSIZ .EQ. 2) THEN
              IWPOSP2 = IOLDPS+IW(IOLDPS+1+XSIZE)+6
              IW(IWPOSP2+NFRONT+XSIZE) =
     &                              -IW(IWPOSP2+NFRONT+XSIZE)
            ENDIF
            IW(IOLDPS+1+XSIZE) = IW(IOLDPS+1+XSIZE) + PIVSIZ
            IF (IFINB.EQ.0) THEN
              GOTO 50 
            ELSE IF (IFINB.EQ.-1) THEN
              LASTBL = .TRUE.
            ENDIF
          ENDIF
          IF ( KEEP(201).EQ.1.AND.(PIVOT_OPTION.GE.3)
     &        .AND.
     &          ( .NOT. LR_ACTIVATED .OR. (.NOT. COMPRESS_PANEL) .OR.
     &            (KEEP(485).EQ.0) 
     &          ) 
     &       ) THEN
            MonBloc%Last = LASTBL
            MonBloc%LastPiv= IW(IOLDPS+1+XSIZE)
            LAST_CALL=.FALSE.
            CALL DMUMPS_OOC_IO_LU_PANEL(
     &        STRAT_TRY_WRITE,
     &        TYPEF_L, A(POSELT),
     &        LAFAC, MonBloc, NextPiv2beWritten, IDUMMY,
     &        IW(IOLDPS), LIWFAC,
     &        MYID, KEEP8(31), IFLAG_OOC,LAST_CALL )
            IF (IFLAG_OOC < 0 ) THEN
              IFLAG=IFLAG_OOC
              GOTO 500
            ENDIF
          ENDIF
          NPIV       =  IW(IOLDPS+1+XSIZE)
          IF ( IEND_BLR .GT. IEND_BLOCK ) THEN
              IF (LR_ACTIVATED) THEN
                  CALL SYSTEM_CLOCK(T1P)
              ENDIF
              CALL DMUMPS_FAC_SQ_LDLT(IBEG_BLOCK,IEND_BLOCK,
     &             NPIV, NFRONT,NASS,IEND_BLR,INODE,A,LA,
     &             LDA, POSELT,
     &             KEEP,KEEP8,
     &             PIVOT_OPTION, .FALSE.)
              IF (LR_ACTIVATED) THEN
                  CALL SYSTEM_CLOCK(T2P,CRP)
                  LOC_FAC_SQ_TIME = LOC_FAC_SQ_TIME + 
     &                              dble(T2P-T1P)/dble(CRP)
              ENDIF
          ENDIF
        END DO 
        NPIV   =  IW(IOLDPS+1+XSIZE)
        IF (.NOT. LR_ACTIVATED
     &      .OR. (.NOT. COMPRESS_PANEL)
     &     ) THEN
          CALL DMUMPS_FAC_SQ_LDLT(IBEG_BLR,IEND_BLR,NPIV,
     &             NFRONT,NASS,NASS,INODE,A,LA,
     &             LDA, POSELT,
     &             KEEP,KEEP8, PIVOT_OPTION, .TRUE.)
        ELSE
          CALL SYSTEM_CLOCK(T2,COUNT_RATE)
          LOC_FRPANELS_TIME = LOC_FRPANELS_TIME +
     &                        dble(T2-T1)/dble(COUNT_RATE)
          CALL UPDATE_FLOP_STATS_PANEL(NFRONT - IBEG_BLR + 1,
     &                                 NPIV   - IBEG_BLR + 1, 1, 1)
          NELIM = IEND_BLOCK - NPIV
          IF (NELIM .EQ. IEND_BLR - IBEG_BLR + 1) THEN
            GOTO 100
          ENDIF
            ALLOCATE(BLR_L(NB_BLR-CURRENT_BLR)) 
          CALL SYSTEM_CLOCK(T1)
#if defined(BLR_MT)          
!$OMP PARALLEL 
#endif
          CALL DMUMPS_COMPRESS_PANEL(A, LA, POSELT, IFLAG, IERROR, 
     &        NFRONT,
     &        BEGS_BLR, NB_BLR, DKEEP(8), KEEP(473), BLR_L, CURRENT_BLR,
     &        'V', WORK, TAU, JPVT, LWORK, RWORK,
     &        BLOCK, MAXI_CLUSTER, NELIM, 
     &        .FALSE., 0, 0,
     &        1, KEEP(483), KEEP(470), KEEP8,
     &        K480=KEEP(480)
     &        )
          IF (IFLAG.LT.0) GOTO 400
#if defined(BLR_MT)          
!$OMP BARRIER
!$OMP MASTER
#endif          
          CALL SYSTEM_CLOCK(T2,COUNT_RATE)
          LOC_DEMOTING_TIME = LOC_DEMOTING_TIME + 
     &              dble(T2-T1)/dble(COUNT_RATE)
          CALL SYSTEM_CLOCK(T1)
#if defined(BLR_MT)          
!$OMP END MASTER
#endif          
            CALL DMUMPS_BLR_UPDATE_TRAILING_LDLT(A, LA, POSELT, 
     &        IFLAG, IERROR, NFRONT,
     &        BEGS_BLR, NB_BLR, CURRENT_BLR, BLR_L, NELIM,
     &        IW(HF+IOLDPS+NFRONT+IBEG_BLR-1), BLOCK,
     &        MAXI_CLUSTER, NPIV,
     &        1, 
     &        KEEP(481), DKEEP(8), KEEP(477) 
     &        )
            IF (IFLAG.LT.0) GOTO 400
 400      CONTINUE         
#if defined(BLR_MT)          
!$OMP END PARALLEL
#endif          
          IF (IFLAG.LT.0) GOTO 500
          CALL SYSTEM_CLOCK(T2,COUNT_RATE)
          LOC_UPDT_TIME = LOC_UPDT_TIME + 
     &               dble(T2-T1)/dble(COUNT_RATE)
          CALL STATS_STORE_BLR_PANEL_MRY(BLR_L,
     &               NB_BLR-CURRENT_BLR-NPARTSCB,
     &               NPARTSCB, 'V',1)
            IF (KEEP(485).NE.0) THEN
              CALL SYSTEM_CLOCK(T1)
              CALL DMUMPS_DECOMPRESS_PANEL(A, LA, POSELT, NFRONT,
     &          .FALSE.,   
     &          BEGS_BLR(CURRENT_BLR),
     &          BEGS_BLR(CURRENT_BLR+1), NB_BLR, BLR_L, CURRENT_BLR,'V',
     &          NFRONT, KEEP(470))
              CALL SYSTEM_CLOCK(T2,COUNT_RATE)
              LOC_PROMOTING_TIME = LOC_PROMOTING_TIME + 
     &               dble(T2-T1)/dble(COUNT_RATE)
            ENDIF
          CALL DEALLOC_BLR_PANEL (BLR_L, NB_BLR-CURRENT_BLR, KEEP8,
     &                             .TRUE.)
          DEALLOCATE(BLR_L)
          NULLIFY(BLR_L)
        ENDIF
        IF (KEEP(201).EQ.1.AND.(PIVOT_OPTION.GE.3)) THEN 
             MonBloc%Last = LASTBL
             MonBloc%LastPiv= NPIV
             LAST_CALL=.FALSE.
             CALL DMUMPS_OOC_IO_LU_PANEL(
     &          STRAT_TRY_WRITE,
     &          TYPEF_L, A(POSELT),
     &          LAFAC, MonBloc, NextPiv2beWritten, IDUMMY, IW(IOLDPS),
     &          LIWFAC, MYID, KEEP8(31), IFLAG_OOC,LAST_CALL )
             IF (IFLAG_OOC < 0 ) THEN
                IFLAG=IFLAG_OOC
                GOTO 500
             ENDIF
        ENDIF
 100    CONTINUE
      END DO 
      IF (LR_ACTIVATED)  THEN
        IF (COMPRESS_CB) THEN
             CALL DMUMPS_FAKE_COMPRESS_CB(A, LA, POSELT, NFRONT,
     &         BEGS_BLR, NPARTSCB+NPARTSASS, 
     &         BEGS_BLR, NPARTSCB+NPARTSASS, NPARTSASS,
     &         DKEEP(8), NASS, NFRONT-NASS,
     &         KEEP(50), WORK, TAU, JPVT, LWORK, RWORK, BLOCK,
     &         MAXI_CLUSTER, STEP_STATS(INODE), 1, 
     &         .FALSE., 0, KEEP(484)) 
        END IF
         CALL SYSTEM_CLOCK(TTOT2,COUNT_RATETOT)
         CALL STATS_COMPUTE_MRY_FRONT_TYPE1(NASS, NFRONT-NASS,
     &                      KEEP(50), INODE, NASS-NPIV )
         CALL STATS_COMPUTE_FLOP_FRONT_TYPE1(NFRONT, NASS, NPIV, 
     &                      KEEP(50), INODE)
         LOC_LR_MODULE_TIME = dble(TTOT2-TTOT1)/dble(COUNT_RATETOT)
         IF (allocated(RWORK))  DEALLOCATE(RWORK)
         IF (allocated(work))  deallocate(WORK)
         IF (allocated(TAU))   deallocate(TAU)
         IF (allocated(JPVT))  deallocate(JPVT)
         IF (allocated(BLOCK)) deallocate(BLOCK)
         IF (associated(BEGS_BLR)) THEN
           DEALLOCATE(BEGS_BLR)
           NULLIFY(BEGS_BLR)
         ENDIF
      ENDIF
      IF (COMPRESS_PANEL)  THEN
        IF ( PIVOT_OPTION.NE.3
     &   )  THEN
          WRITE(*,*) "Internal error 1 in DMUMPS_FACTO_NIV1"
          CALL MUMPS_ABORT()
        ENDIF
      ELSE
        CALL DMUMPS_FAC_T_LDLT(NFRONT,NASS,IW,LIW,A,LA,
     &         LDA, IOLDPS,POSELT, KEEP,KEEP8,
     &         (PIVOT_OPTION.NE.3), ETATASS,
     &         TYPEF_L, LAFAC, MonBloc, NextPiv2beWritten,
     &         LIWFAC, MYID, IFLAG, IOLDPS+6+XSIZE+NFRONT)
      ENDIF
      IF (.NOT.LR_ACTIVATED) THEN
          CALL SYSTEM_CLOCK(TTOT2FR,COUNT_RATETOTFR)
          LOC_FRFRONTS_TIME =
     &       DBLE(TTOT2FR-TTOT1FR)/DBLE(COUNT_RATETOTFR)
          CALL UPDATE_FLOP_STATS_FRFRONTS(NFRONT, NPIV, NASS, 1, 1)
      ELSE
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
     &          ( STRAT, TYPEF_L, 
     &           A(POSELT), LAFAC, MonBloc,
     &           NextPiv2beWritten, IDUMMY,
     &           IW(IOLDPS), LIWFAC, 
     &           MYID, KEEP8(31), IFLAG_OOC,LAST_CALL )
          IF (IFLAG_OOC < 0 ) THEN
            IFLAG=IFLAG_OOC
            GOTO 500
          ENDIF
          CALL DMUMPS_OOC_PP_TRYRELEASE_SPACE (IWPOS, 
     &      IOLDPS, IW, LIW, MonBloc , NFRONT, KEEP)
      ENDIF
      GOTO 500
 490  CONTINUE
      write(*,*) 'Allocation problem in BLR routine
     &     DMUMPS_FAC_FRONT_LDLT_TYPE1: ',
     &     'not enough memory? memory requested = ' , IERROR
 500  CONTINUE
      RETURN
      END SUBROUTINE DMUMPS_FAC1_LDLT
      END MODULE DMUMPS_FAC1_LDLT_M
