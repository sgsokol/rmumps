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
      RECURSIVE SUBROUTINE DMUMPS_PROCESS_BLOCFACTO(
     &   COMM_LOAD, ASS_IRECV,
     &   BUFR, LBUFR,
     &   LBUFR_BYTES, PROCNODE_STEPS, MSGSOU,
     &   SLAVEF, IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, N, IW, LIW,
     &   A, LA, PTRIST, PTRAST, NSTK_S, NBPROCFILS,
     &   COMP, STEP, PIMASTER, PAMASTER, POSFAC,
     &   MYID, COMM, IFLAG, IERROR, NBFIN,
     &
     &    PTLUST_S, PTRFAC, root, OPASSW, OPELIW,
     &    ITLOC, RHS_MUMPS, FILS,  
     &    PTRARW, PTRAIW, INTARR, DBLARR,
     &    ICNTL, KEEP,KEEP8, DKEEP, 
     &    IPOOL, LPOOL, LEAF, ND, FRERE_STEPS,
     &    LPTRAR, NELT, FRTPTR, FRTELT, 
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE
     &               , LRGROUPS
     &    )
      USE DMUMPS_OOC
      USE DMUMPS_LOAD
      USE DMUMPS_LR_STATS
      USE DMUMPS_LR_CORE
      USE DMUMPS_LR_TYPE
      USE DMUMPS_FAC_LR, ONLY : DMUMPS_DECOMPRESS_PANEL,
     &                          DMUMPS_COMPRESS_PANEL,
     &                          DMUMPS_BLR_UPDATE_TRAILING, 
     &                          DMUMPS_FAKE_COMPRESS_CB
      USE DMUMPS_ANA_LR, ONLY : GET_CUT
!$    USE OMP_LIB
      IMPLICIT NONE
      INCLUDE 'dmumps_root.h'
      INCLUDE 'mumps_headers.h'
      TYPE (DMUMPS_ROOT_STRUC) :: root
      INTEGER ICNTL( 40 ), KEEP( 500 )
      INTEGER(8) KEEP8(150)
      DOUBLE PRECISION DKEEP(230)
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER COMM_LOAD, ASS_IRECV
      INTEGER BUFR( LBUFR )
      INTEGER N, SLAVEF, IWPOS, IWPOSCB, LIW
      INTEGER(8) :: IPTRLU, LRLU, LRLUS, LA
      INTEGER(8) :: POSFAC
      INTEGER COMP
      INTEGER IFLAG, IERROR, NBFIN, MSGSOU
      INTEGER PROCNODE_STEPS(KEEP(28)), PTRIST(KEEP(28)),
     &        NSTK_S(KEEP(28))
      INTEGER(8) :: PAMASTER(KEEP(28))
      INTEGER(8) :: PTRAST(KEEP(28))
      INTEGER(8) :: PTRFAC(KEEP(28))
      INTEGER NBPROCFILS( KEEP(28) ), STEP(N), 
     & PIMASTER(KEEP(28))
      INTEGER IW( LIW )
      DOUBLE PRECISION A( LA )
      INTEGER, intent(in) :: LRGROUPS(N)
      INTEGER COMM, MYID
      INTEGER NELT, LPTRAR
      INTEGER FRTPTR( N+1 ), FRTELT( NELT )
      INTEGER PTLUST_S(KEEP(28)),
     &        ITLOC(N+KEEP(253)), FILS(N), ND(KEEP(28))
      DOUBLE PRECISION :: RHS_MUMPS(KEEP(255))
      INTEGER(8), INTENT(IN) :: PTRAIW( LPTRAR ), PTRARW( LPTRAR )
      INTEGER FRERE_STEPS(KEEP(28))
      DOUBLE PRECISION OPASSW, OPELIW
      DOUBLE PRECISION FLOP1
      INTEGER INTARR( KEEP8(27) )
      DOUBLE PRECISION DBLARR( KEEP8(26) )
      INTEGER LEAF, LPOOL 
      INTEGER IPOOL( LPOOL )
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     &        TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER :: STATUS(MPI_STATUS_SIZE)
      LOGICAL :: I_HAVE_SET_K117
      INTEGER INODE, POSITION, NPIV, IERR, LP
      INTEGER NCOL
      INTEGER(8) :: POSBLOCFACTO
      INTEGER :: LD_BLOCFACTO 
      INTEGER(8) :: LAELL
      INTEGER(8) :: POSELT
      INTEGER IOLDPS, LCONT1, NASS1, NROW1, NCOL1, NPIV1
      INTEGER NSLAV1, HS, ISW
      INTEGER (8) :: LPOS, LPOS1, LPOS2, IPOS, KPOS
      INTEGER ICT11
      INTEGER I, IPIV, FPERE
      LOGICAL LASTBL
      LOGICAL BLOCKING, SET_IRECV, MESSAGE_RECEIVED
      DOUBLE PRECISION ONE,ALPHA
      PARAMETER (ONE = 1.0D0, ALPHA=-1.0D0)
      INTEGER(8) :: LAFAC
      INTEGER LIWFAC, STRAT, NextPivDummy
      TYPE(IO_BLOCK) :: MonBloc
      LOGICAL LAST_CALL
      INTEGER LRELAY_INFO
      INTEGER :: SEND_LR_INT, NELIM, NPARTSASS_MASTER, 
     &           CURRENT_BLR_PANEL, 
     &           CURRENT_BLR, 
     &           NB_BLR_L, NB_BLR_U
      TYPE (LRB_TYPE), DIMENSION(:), POINTER :: BLR_U, BLR_L
      LOGICAL :: SEND_LR
      INTEGER, POINTER, DIMENSION(:) :: BEGS_BLR_L, BEGS_BLR_U
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: WORK, TAU
      INTEGER, ALLOCATABLE, DIMENSION(:) :: JPVT
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: RWORK
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: BLOCK
      INTEGER :: OMP_NUM
      INTEGER NPARTSASS, NPARTSCB, MAXI_CLUSTER, LWORK,
     &        MAXI_CLUSTER_L, MAXI_CLUSTER_U
      INTEGER T1, T2, COUNT_RATE
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DYN_BLOCFACTO
      INTEGER, DIMENSION(:), ALLOCATABLE :: DYN_PIVINFO
      LOGICAL :: DYNAMIC_ALLOC
      INTEGER MUMPS_PROCNODE
      EXTERNAL MUMPS_PROCNODE
      I_HAVE_SET_K117 = .FALSE.
      DYNAMIC_ALLOC = .FALSE.
      FPERE    = -1
      POSITION = 0
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, INODE, 1,
     &                 MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, NPIV, 1,
     &                 MPI_INTEGER, COMM, IERR )
      LASTBL = (NPIV.LE.0)
      IF (LASTBL) THEN 
         NPIV = -NPIV
         CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, FPERE, 1,
     &                 MPI_INTEGER, COMM, IERR )
      ENDIF
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, NCOL, 1,
     &                 MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, NELIM, 1,
     &                 MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, 
     &                 NPARTSASS_MASTER , 1,
     &                 MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, CURRENT_BLR_PANEL,
     &                 1, MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, SEND_LR_INT, 1,
     &                 MPI_INTEGER, COMM, IERR )
      IF ( SEND_LR_INT .EQ. 1) THEN
        SEND_LR = .TRUE.
      ELSE
        SEND_LR = .FALSE.
      ENDIF
      IF ( SEND_LR ) THEN
        LAELL = int(NPIV,8) * int(NPIV+NELIM,8)
      ELSE
        LAELL = int(NPIV,8) * int(NCOL,8)
      ENDIF
      IF ( LRLU .LT. LAELL .OR. IWPOS + NPIV - 1 .GT. IWPOSCB ) THEN
        IF ( LRLUS .LT. LAELL ) THEN
          IFLAG = -9
          CALL MUMPS_SET_IERROR(LAELL - LRLUS, IERROR)
          IF (ICNTL(1).GT.0 .AND. ICNTL(4).GE.1) THEN
            LP=ICNTL(1)
            WRITE(LP,*)
     &" FAILURE, WORKSPACE TOO SMALL DURING DMUMPS_PROCESS_BLOCFACTO"
          ENDIF
          GOTO 700
        END IF
        CALL DMUMPS_COMPRE_NEW(N, KEEP(28), IW, LIW, A, LA,
     &      LRLU, IPTRLU,
     &      IWPOS, IWPOSCB, PTRIST, PTRAST,
     &      STEP, PIMASTER, PAMASTER, KEEP(216),LRLUS,
     &      KEEP(IXSZ),COMP,DKEEP(97),MYID)
        IF ( LRLU .NE. LRLUS ) THEN
        WRITE(*,*) 'PB compress DMUMPS_PROCESS_BLOCFACTO, LRLU,LRLUS='
     &       ,LRLU,LRLUS
             IFLAG = -9
             CALL MUMPS_SET_IERROR( LAELL-LRLUS, IERROR )
             GOTO 700
        END IF
        IF ( IWPOS + NPIV - 1 .GT. IWPOSCB ) THEN
          IF (ICNTL(1).GT.0 .AND. ICNTL(4).GE.1) THEN
            LP=ICNTL(1)
            WRITE(LP,*)
     &" FAILURE IN INTEGER ALLOCATION DURING DMUMPS_PROCESS_BLOCFACTO"
          ENDIF
          IFLAG = -8
          IERROR = IWPOS + NPIV - 1 - IWPOSCB
          GOTO 700
        END IF
      END IF
      LRLU  = LRLU - LAELL
      LRLUS = LRLUS - LAELL
      KEEP8(67) = min(LRLUS, KEEP8(67))
      KEEP8(70) = KEEP8(70) - LAELL
      KEEP8(68) = min(KEEP8(70), KEEP8(68))
      KEEP8(71) = KEEP8(71) - LAELL
      KEEP8(69) = min(KEEP8(71), KEEP8(69))
      POSBLOCFACTO = POSFAC
      POSFAC = POSFAC + LAELL
      CALL DMUMPS_LOAD_MEM_UPDATE(.FALSE., .FALSE.,
     &               LA-LRLUS,0_8,LAELL,KEEP,KEEP8,LRLUS)
      IF ((NPIV .EQ. 0) 
     &     ) THEN
        IPIV=1
      ELSE
        IPIV = IWPOS
        IWPOS = IWPOS + NPIV
        IF (NPIV .GT. 0) THEN
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                 IW( IPIV ), NPIV,
     &                 MPI_INTEGER, COMM, IERR )
        ENDIF
        IF ( SEND_LR ) THEN
            CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                 A(POSBLOCFACTO), NPIV*(NPIV+NELIM),
     &                 MPI_DOUBLE_PRECISION,
     &                 COMM, IERR )
            LD_BLOCFACTO = NPIV+NELIM
            CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                 NB_BLR_U, 1, MPI_INTEGER,
     &                 COMM, IERR )
            ALLOCATE(BLR_U(max(NB_BLR_U,1)))
            ALLOCATE(BEGS_BLR_U(NB_BLR_U+2))
            CALL DMUMPS_MPI_UNPACK_LR(BUFR, LBUFR, LBUFR_BYTES, 
     &                             POSITION, NPIV, NELIM, 'H',
     &                             BLR_U(1:), NB_BLR_U, KEEP(470), 
     &                             BEGS_BLR_U(1:),
     &                             KEEP8, COMM, IERR, IFLAG, IERROR)
            IF (IFLAG.LT.0) GOTO 700
        ELSE
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                 A(POSBLOCFACTO), NPIV*NCOL,
     &                 MPI_DOUBLE_PRECISION,
     &                 COMM, IERR )
          LD_BLOCFACTO = NCOL
        ENDIF
      ENDIF 
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, 
     &                 LRELAY_INFO, 1,
     &                 MPI_INTEGER, COMM, IERR )
      IF (PTRIST(STEP( INODE )) .EQ. 0) THEN
          CALL DMUMPS_TREAT_DESCBAND( INODE, COMM_LOAD,
     &    ASS_IRECV,
     &    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     &    IWPOS, IWPOSCB, IPTRLU,
     &    LRLU, LRLUS, N, IW, LIW, A, LA, PTRIST,
     &    PTLUST_S, PTRFAC,
     &    PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP,
     &    IFLAG, IERROR, COMM,
     &    NBPROCFILS, IPOOL, LPOOL, LEAF, NBFIN, MYID, SLAVEF,
     &
     &    root, OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &    FILS, PTRARW, PTRAIW,
     &    INTARR, DBLARR, ICNTL,KEEP,KEEP8,DKEEP,ND, FRERE_STEPS, 
     &    LPTRAR, NELT, FRTPTR, FRTELT, 
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE, .TRUE.
     &               , LRGROUPS
     &    )
          IF ( IFLAG .LT. 0 ) GOTO 600
      ENDIF
      IF ( IW( PTRIST(STEP(INODE)) + 3 +KEEP(IXSZ)) .EQ. 0 ) THEN
#if ! defined(NO_XXNBPR)
       CALL CHECK_EQUAL(NBPROCFILS(STEP(INODE)),
     &                  IW(PTRIST(STEP(INODE))+XXNBPR))
       DO WHILE ( IW(PTRIST(STEP(INODE)) + XXNBPR) .NE. 0)
#else
       DO WHILE ( NBPROCFILS( STEP(INODE)) .NE. 0 )
#endif
        BLOCKING = .TRUE.
        SET_IRECV = .FALSE.
        MESSAGE_RECEIVED = .FALSE.
        CALL DMUMPS_TRY_RECVTREAT( COMM_LOAD,
     &    ASS_IRECV, BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     &    MPI_ANY_SOURCE, CONTRIB_TYPE2,
     &    STATUS,
     &    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     &    IWPOS, IWPOSCB, IPTRLU,
     &    LRLU, LRLUS, N, IW, LIW, A, LA, PTRIST,
     &    PTLUST_S, PTRFAC,
     &    PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP,
     &    IFLAG, IERROR, COMM,
     &    NBPROCFILS, IPOOL, LPOOL, LEAF, NBFIN, MYID, SLAVEF,
     &
     &    root, OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &    FILS, PTRARW, PTRAIW,
     &    INTARR, DBLARR, ICNTL,KEEP,KEEP8,DKEEP,ND, FRERE_STEPS,
     &    LPTRAR, NELT, FRTPTR, FRTELT, 
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE, .TRUE. 
     &               , LRGROUPS
     &      )
        IF ( IFLAG .LT. 0 ) GOTO 600
       END  DO
      ENDIF
        SET_IRECV = .TRUE.
        BLOCKING  = .FALSE.
        MESSAGE_RECEIVED = .TRUE.
        CALL DMUMPS_TRY_RECVTREAT( COMM_LOAD, ASS_IRECV,
     &    BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     &    MPI_ANY_SOURCE, MPI_ANY_TAG, 
     &    STATUS,
     &    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     &    IWPOS, IWPOSCB, IPTRLU,
     &    LRLU, LRLUS, N, IW, LIW, A, LA, PTRIST,
     &    PTLUST_S, PTRFAC,
     &    PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP,
     &    IFLAG, IERROR, COMM,
     &    NBPROCFILS, IPOOL, LPOOL, LEAF, NBFIN, MYID, SLAVEF,
     &
     &    root, OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &    FILS, PTRARW, PTRAIW,
     &    INTARR, DBLARR,ICNTL,KEEP,KEEP8,DKEEP,ND, FRERE_STEPS,
     &    LPTRAR, NELT, FRTPTR, FRTELT,
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE, .TRUE.  
     &               , LRGROUPS
     &    )
      IOLDPS = PTRIST(STEP(INODE))
      POSELT = PTRAST(STEP(INODE))
      LCONT1 = IW( IOLDPS +KEEP(IXSZ))
      NASS1  = IW( IOLDPS + 1 +KEEP(IXSZ))
      IF ( NASS1 < 0 ) THEN
        NASS1 = -NASS1
        IW( IOLDPS + 1 + KEEP(IXSZ)) = NASS1
        IF (KEEP(55) .EQ. 0) THEN 
          CALL DMUMPS_ASM_SLAVE_ARROWHEADS(INODE, N, IW, LIW,
     &       IOLDPS, A, LA, POSELT, KEEP, KEEP8, ITLOC, FILS, PTRAIW,
     &       PTRARW, INTARR, DBLARR, KEEP8(27), KEEP8(26), RHS_MUMPS)
        ELSE
          CALL DMUMPS_ASM_SLAVE_ELEMENTS(INODE, N, NELT, IW, LIW,
     &       IOLDPS, A, LA, POSELT, KEEP, KEEP8, ITLOC, FILS, PTRAIW,
     &       PTRARW, INTARR, DBLARR, KEEP8(27), KEEP8(26),
     &       FRTPTR, FRTELT, RHS_MUMPS)
        ENDIF
      ENDIF
      NROW1  = IW( IOLDPS + 2 +KEEP(IXSZ))
      NPIV1  = IW( IOLDPS + 3 +KEEP(IXSZ))
      NSLAV1 = IW( IOLDPS + 5 +KEEP(IXSZ))
      HS     = 6 + NSLAV1 + KEEP(IXSZ)
      NCOL1  = LCONT1 + NPIV1
      IF (NPIV.GT.0) THEN
        ICT11 = IOLDPS+HS+NROW1+NPIV1 - 1
        IF (DYNAMIC_ALLOC) THEN
          DO I = 1, NPIV
            IF (DYN_PIVINFO(I).EQ.I) CYCLE
            ISW = IW(ICT11+I)
            IW(ICT11+I) = IW(ICT11+DYN_PIVINFO(I))
            IW(ICT11+DYN_PIVINFO(I)) = ISW
            IPOS = POSELT + int(NPIV1 + I - 1,8)
            KPOS = POSELT + int(NPIV1 + DYN_PIVINFO(I) - 1,8)
            CALL dswap(NROW1, A(IPOS), NCOL1, A(KPOS), NCOL1)
          ENDDO
        ELSE
          DO I = 1, NPIV
            IF (IW(IPIV+I-1).EQ.I) CYCLE
            ISW = IW(ICT11+I)
            IW(ICT11+I) = IW(ICT11+IW(IPIV+I-1))
            IW(ICT11+IW(IPIV+I-1)) = ISW
            IPOS = POSELT + int(NPIV1 + I - 1,8)
            KPOS = POSELT + int(NPIV1 + IW(IPIV+I-1) - 1,8)
            CALL dswap(NROW1, A(IPOS), NCOL1, A(KPOS), NCOL1)
          ENDDO
        ENDIF
        LPOS2 = POSELT + int(NPIV1,8)
        LPOS  = LPOS2 + int(NPIV,8)
        CALL SYSTEM_CLOCK(T1)
        IF (DYNAMIC_ALLOC) THEN
          CALL dtrsm('L','L','N','N',NPIV, NROW1, ONE, 
     &           DYN_BLOCFACTO, LD_BLOCFACTO, A(LPOS2), NCOL1)
        ELSE
          CALL dtrsm('L','L','N','N',NPIV, NROW1, ONE, 
     &           A(POSBLOCFACTO), LD_BLOCFACTO, A(LPOS2), NCOL1)
        ENDIF
        CALL SYSTEM_CLOCK(T2,COUNT_RATE)
        ACC_TRSM_TIME = ACC_TRSM_TIME + 
     &            DBLE(T2-T1)/DBLE(COUNT_RATE)
      ENDIF 
      IF ( (NPIV .GT. 0)
     &   ) THEN
        IF (SEND_LR) THEN
             CALL GET_CUT(IW(IOLDPS+HS:IOLDPS+HS+NROW1-1), 0,
     &                    NROW1, LRGROUPS, NPARTSCB, 
     &                    NPARTSASS, BEGS_BLR_L)
              CALL REGROUPING2(BEGS_BLR_L, NPARTSASS, 0, NPARTSCB,
     &                        NROW1-0, KEEP(488), .TRUE., KEEP(472))
             NB_BLR_L =  NPARTSCB
        call MAX_CLUSTER(BEGS_BLR_L,NB_BLR_L+1,MAXI_CLUSTER_L)
        call MAX_CLUSTER(BEGS_BLR_U,NB_BLR_U+1,MAXI_CLUSTER_U)
        MAXI_CLUSTER=max(MAXI_CLUSTER_U,MAXI_CLUSTER_L)
        LWORK = MAXI_CLUSTER*MAXI_CLUSTER
        OMP_NUM = 1
#if defined(BLR_MT)
!$      OMP_NUM = OMP_GET_MAX_THREADS()
#endif
         ALLOCATE(BLOCK(MAXI_CLUSTER, OMP_NUM*MAXI_CLUSTER),
     &            RWORK(2*MAXI_CLUSTER*OMP_NUM), 
     &            TAU(MAXI_CLUSTER*OMP_NUM),
     &            JPVT(MAXI_CLUSTER*OMP_NUM), 
     &            WORK(LWORK*OMP_NUM))
        CURRENT_BLR=1 
        ALLOCATE(BLR_L(NB_BLR_L))
        CALL SYSTEM_CLOCK(T1)
#if defined(BLR_MT)          
!$OMP PARALLEL
#endif
        CALL DMUMPS_COMPRESS_PANEL
     &        (A, LA, POSELT, IFLAG, IERROR, NCOL1,
     &        BEGS_BLR_L, NB_BLR_L+1, DKEEP(8), KEEP(473), BLR_L, 
     &        CURRENT_BLR, 'V', WORK, TAU, JPVT, LWORK, RWORK,
     &        BLOCK, MAXI_CLUSTER, NELIM, 
     &        .TRUE.,  
     &        NPIV, NPIV1,
     &        2, KEEP(483), KEEP(470), KEEP8
     &        )
        IF (IFLAG.LT.0) GOTO 300
#if defined(BLR_MT)          
!$OMP BARRIER
!$OMP MASTER
#endif          
          CALL SYSTEM_CLOCK(T2,COUNT_RATE)
          ACC_DEMOTING_TIME = ACC_DEMOTING_TIME + 
     &              DBLE(T2-T1)/DBLE(COUNT_RATE)
          CALL SYSTEM_CLOCK(T1)
#if defined(BLR_MT)          
!$OMP END MASTER
#endif          
 300      CONTINUE
#if defined(BLR_MT)          
!$OMP END PARALLEL
#endif          
          IF (IFLAG.LT.0) GOTO 700
        ENDIF
      ENDIF
      IF ( (KEEP(201).eq.1) .AND.
     &    ( .NOT. SEND_LR .OR. (NPIV .EQ.0) .OR.
     &    (KEEP(485).EQ.0) ) 
     &   ) THEN
        MonBloc%INODE = INODE
        MonBloc%MASTER = .FALSE.
        MonBloc%Typenode = 2
        MonBloc%NROW = NROW1
        MonBloc%NCOL = NCOL1  
        MonBloc%NFS  = NASS1
        MonBloc%LastPiv = NPIV1 + NPIV 
        MonBloc%LastPanelWritten_L=-9999 
        MonBloc%LastPanelWritten_U=-9999 
        NULLIFY(MonBloc%INDICES)
        MonBloc%Last = LASTBL
        STRAT = STRAT_TRY_WRITE 
        NextPivDummy      = -8888 
        LIWFAC = IW(IOLDPS+XXI)
        CALL MUMPS_GETI8(LAFAC, IW(IOLDPS+XXR))
        LAST_CALL = .FALSE.
        CALL DMUMPS_OOC_IO_LU_PANEL( STRAT, TYPEF_L, A(POSELT),
     &       LAFAC, MonBloc, NextPivDummy, NextPivDummy,
     &       IW(IOLDPS), LIWFAC, MYID, KEEP8(31), IFLAG,LAST_CALL)
      ENDIF
      IF ( (NPIV .GT. 0)
     &   ) THEN
        IF (SEND_LR) THEN
        IF (NELIM.GT.0) THEN
          IF (DYNAMIC_ALLOC) THEN
            LPOS1 = int(NPIV+1,8)
            CALL dgemm('N','N', NELIM,NROW1,NPIV,
     &             ALPHA,DYN_BLOCFACTO(LPOS1),LD_BLOCFACTO,
     &             A(LPOS2),NCOL1,ONE, A(LPOS),NCOL1)
          ELSE
            LPOS1 = POSBLOCFACTO+int(NPIV,8)
            CALL dgemm('N','N', NELIM,NROW1,NPIV,
     &             ALPHA,A(LPOS1),LD_BLOCFACTO,
     &             A(LPOS2),NCOL1,ONE, A(LPOS),NCOL1)
          ENDIF
        ENDIF
#if defined(BLR_MT)          
!$OMP PARALLEL
#endif          
         CALL DMUMPS_BLR_UPDATE_TRAILING(A, LA, POSELT, 
     &        IFLAG, IERROR, NCOL1,
     &        BEGS_BLR_L, BEGS_BLR_U, CURRENT_BLR, 
     &        BLR_L, NB_BLR_L+1, 
     &        BLR_U, NB_BLR_U+1, 
     &        0,      
     &        .TRUE., 
     &        NPIV1,  
     &        2, 0, KEEP(470), 
     &        KEEP(481), DKEEP(8), KEEP(477) 
     &        )
 400      CONTINUE          
#if defined(BLR_MT)          
!$OMP END PARALLEL
#endif          
          IF (IFLAG.LT.0) GOTO 700
          CALL SYSTEM_CLOCK(T2,COUNT_RATE)
          ACC_UPDT_TIME = ACC_UPDT_TIME + 
     &               DBLE(T2-T1)/DBLE(COUNT_RATE)
          CALL STATS_STORE_BLR_PANEL_MRY(BLR_L,
     &               0, NPARTSCB, 'V', 2)
            IF (KEEP(485).NE.0) THEN
              CALL SYSTEM_CLOCK(T1)
              CALL DMUMPS_DECOMPRESS_PANEL(A, LA, POSELT, NCOL1,
     &        .FALSE.,   
     &        NPIV1+1,   
     &        1,         
     &        NB_BLR_L+1, BLR_L, CURRENT_BLR, 'V', NCOL1, KEEP(470))
              CALL SYSTEM_CLOCK(T2,COUNT_RATE)
              ACC_PROMOTING_TIME = ACC_PROMOTING_TIME + 
     &                  DBLE(T2-T1)/DBLE(COUNT_RATE)
             IF (KEEP(201).eq.1) THEN
              MonBloc%INODE = INODE
              MonBloc%MASTER = .FALSE.
              MonBloc%Typenode = 2
              MonBloc%NROW = NROW1
              MonBloc%NCOL = NCOL1  
              MonBloc%NFS  = NASS1
              MonBloc%LastPiv = NPIV1 + NPIV 
              MonBloc%LastPanelWritten_L=-9999 
              MonBloc%LastPanelWritten_U=-9999 
              NULLIFY(MonBloc%INDICES)
              MonBloc%Last = LASTBL
              STRAT = STRAT_TRY_WRITE 
              NextPivDummy      = -8888 
              LIWFAC = IW(IOLDPS+XXI)
              CALL MUMPS_GETI8(LAFAC, IW(IOLDPS+XXR))
              LAST_CALL = .FALSE.
              CALL DMUMPS_OOC_IO_LU_PANEL( STRAT, TYPEF_L, A(POSELT),
     &          LAFAC, MonBloc, NextPivDummy, NextPivDummy,
     &          IW(IOLDPS), LIWFAC, MYID, KEEP8(31), IFLAG,LAST_CALL)
             ENDIF
            ENDIF
        ELSE
          IF (DYNAMIC_ALLOC) THEN
            LPOS1 = int(NPIV+1,8)
            CALL dgemm('N','N', NCOL-NPIV,NROW1,NPIV,
     &             ALPHA,DYN_BLOCFACTO(LPOS1),NCOL,
     &             A(LPOS2),NCOL1,ONE, A(LPOS),NCOL1)
          ELSE
            LPOS1 = POSBLOCFACTO+int(NPIV,8)
            CALL dgemm('N','N', NCOL-NPIV,NROW1,NPIV,
     &             ALPHA,A(LPOS1),NCOL,
     &             A(LPOS2),NCOL1,ONE, A(LPOS),NCOL1)
          ENDIF
        ENDIF
      ENDIF
      IW(IOLDPS+KEEP(IXSZ) ) = IW(IOLDPS+KEEP(IXSZ) ) - NPIV
      IW(IOLDPS + 3+KEEP(IXSZ) ) = IW(IOLDPS+3+KEEP(IXSZ) ) + NPIV
      IF (LASTBL) THEN
        IW(IOLDPS+1+KEEP(IXSZ) ) = IW(IOLDPS + 3+KEEP(IXSZ) )
      ENDIF
      IF ( .not. LASTBL .AND. 
     &  (IW(IOLDPS+1+KEEP(IXSZ)) .EQ. IW(IOLDPS + 3+KEEP(IXSZ))) ) THEN
        write(*,*) 'Internal ERROR 1 **** IN BLACFACTO '
        CALL MUMPS_ABORT()
      ENDIF
      IF (SEND_LR) THEN
        IF ((NPIV.GT.0)
     &     ) THEN
          CALL DEALLOC_BLR_PANEL( BLR_U, NB_BLR_U, KEEP8, .FALSE.)
          DEALLOCATE(BLR_U)
          CALL DEALLOC_BLR_PANEL( BLR_L, NB_BLR_L, KEEP8, .TRUE.)
          DEALLOCATE(BLR_L)
      ENDIF 
      ENDIF 
      IF (DYNAMIC_ALLOC) THEN
        DEALLOCATE(DYN_BLOCFACTO)
        DEALLOCATE(DYN_PIVINFO)
      ELSE
        LRLU  = LRLU + LAELL
        LRLUS = LRLUS + LAELL
        KEEP8(70) = KEEP8(70) + LAELL
        KEEP8(71) = KEEP8(71) + LAELL
        POSFAC = POSFAC - LAELL
        CALL DMUMPS_LOAD_MEM_UPDATE(.FALSE.,.FALSE.,
     &             LA-LRLUS,0_8,-LAELL,KEEP,KEEP8,LRLUS)
        IWPOS = IWPOS - NPIV
      ENDIF
      FLOP1 = dble( NPIV1*NROW1 ) +
     &        dble(NROW1*NPIV1)*dble(2*NCOL1-NPIV1-1)
     &   -
     &        dble((NPIV1+NPIV)*NROW1 ) -
     &        dble(NROW1*(NPIV1+NPIV))*dble(2*NCOL1-NPIV1-NPIV-1)
      CALL DMUMPS_LOAD_UPDATE( 1, .FALSE., FLOP1, KEEP,KEEP8 )
      IF (LASTBL) THEN
       IF (SEND_LR) THEN
         CALL STATS_COMPUTE_FLOP_SLAVE_TYPE2(NROW1, NCOL1, NASS1,
     &              KEEP(50), INODE)
       ELSE
         CALL UPD_FLOP_FRFRONT_SLAVE(NROW1, NCOL1, NASS1,
     &              KEEP(50), INODE)
       ENDIF
       IF (SEND_LR) THEN
         IF (KEEP(489) .EQ. 1) THEN
           CALL DMUMPS_FAKE_COMPRESS_CB(A, LA, POSELT, NCOL1,
     &          BEGS_BLR_L, NB_BLR_L+1,
     &          BEGS_BLR_U, NB_BLR_U+1, 1,
     &          DKEEP(8), NASS1, NROW1,
     &          KEEP(50), WORK, TAU, JPVT, LWORK, RWORK,
     &          BLOCK, MAXI_CLUSTER, STEP_STATS(INODE), 2, 
     &          .TRUE., NPIV1, KEEP(484)) 
         ENDIF
       ENDIF
         CALL DMUMPS_END_FACTO_SLAVE(
     &    COMM_LOAD, ASS_IRECV, 
     &    N, INODE, FPERE, 
     &    root, 
     &    MYID, COMM,
     &    
     &    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     &    IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, IW, LIW, A, LA,
     &    PTRIST, PTLUST_S, PTRFAC, PTRAST, STEP, PIMASTER,
     &    PAMASTER,
     &    NSTK_S, COMP, IFLAG, IERROR, NBPROCFILS,
     &    IPOOL, LPOOL, LEAF, NBFIN, SLAVEF,
     &    OPASSW, OPELIW, ITLOC, RHS_MUMPS, FILS, PTRARW, PTRAIW,
     &    INTARR, DBLARR,ICNTL,KEEP,KEEP8,DKEEP,ND, FRERE_STEPS,
     &    LPTRAR, NELT, FRTPTR, FRTELT, 
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE
     &               , LRGROUPS
     &     )
      ENDIF 
      IF (SEND_LR) THEN
        IF (allocated(RWORK))  DEALLOCATE(RWORK)
        IF (allocated(work)) DEALLOCATE(WORK)
        IF (allocated(TAU)) DEALLOCATE(TAU)
        IF (allocated(JPVT)) DEALLOCATE(JPVT)
        IF (allocated(BLOCK)) DEALLOCATE(BLOCK)
        IF ((NPIV.GT.0)
     &     ) THEN
          IF (associated(BEGS_BLR_L)) DEALLOCATE(BEGS_BLR_L)
          IF (associated(BEGS_BLR_U)) DEALLOCATE(BEGS_BLR_U)
        ENDIF
      ENDIF
 600  CONTINUE
      RETURN
 700  CONTINUE
      CALL DMUMPS_BDC_ERROR( MYID, SLAVEF, COMM, KEEP )
      RETURN
      END SUBROUTINE DMUMPS_PROCESS_BLOCFACTO
      SUBROUTINE DMUMPS_MPI_UNPACK_LR(
     &           BUFR, LBUFR, LBUFR_BYTES, POSITION,
     &                             NPIV, NELIM, DIR, 
     &                             BLR_U, NB_BLOCK_U, K470,
     &                             BEGS_BLR_U, KEEP8,
     &                             COMM, IERR, IFLAG, IERROR)
      USE DMUMPS_LR_CORE, ONLY : LRB_TYPE, ALLOC_LRB
      USE DMUMPS_LR_TYPE
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: LBUFR
      INTEGER, INTENT(IN) :: LBUFR_BYTES
      INTEGER, INTENT(IN) :: BUFR(LBUFR)
      INTEGER, INTENT(INOUT) :: POSITION
      INTEGER, INTENT(IN)    :: NB_BLOCK_U, NELIM, NPIV, K470
      CHARACTER(len=1) :: DIR
      INTEGER, INTENT(IN) :: COMM
      INTEGER, INTENT(OUT) :: IERR, IFLAG, IERROR
      TYPE (LRB_TYPE), INTENT(OUT), 
     &          DIMENSION(max(NB_BLOCK_U,1)):: BLR_U
      INTEGER, INTENT(OUT), DIMENSION(NB_BLOCK_U+2)  :: BEGS_BLR_U 
      INTEGER(8) :: KEEP8(150)
      LOGICAL :: ISLR
      INTEGER :: ISLR_INT, I
      INTEGER :: LRFORM, K, M, N, KSVD
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      IERR = 0
      IF (size(BLR_U) .NE. 
     &    MAX(NB_BLOCK_U,1) ) THEN
        WRITE(*,*) "Internal error 1 in DMUMPS_MPI_UNPACK",
     &             NB_BLOCK_U,size(BLR_U)
        CALL MUMPS_ABORT()
      ENDIF
      BEGS_BLR_U(1) = 1
      BEGS_BLR_U(2) = NPIV+NELIM+1 
      DO I = 1, NB_BLOCK_U
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                 ISLR_INT, 1, MPI_INTEGER, COMM, IERR )
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                 LRFORM, 1,
     &                 MPI_INTEGER, COMM, IERR )
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                 K, 1,
     &                 MPI_INTEGER, COMM, IERR )
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                 M, 1,
     &                 MPI_INTEGER, COMM, IERR )
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                 N, 1,
     &                 MPI_INTEGER, COMM, IERR )
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                 KSVD, 1,
     &                 MPI_INTEGER, COMM, IERR )
        IF (DIR.EQ.'H') THEN
          IF (K470.EQ.1) THEN
            BEGS_BLR_U(I+2) = BEGS_BLR_U(I+1) + M
          ELSE
            BEGS_BLR_U(I+2) = BEGS_BLR_U(I+1) + N
          ENDIF
        ELSE
         BEGS_BLR_U(I+2) = BEGS_BLR_U(I+1) + M
        ENDIF
        IF (ISLR_INT .eq. 1) THEN
          ISLR = .TRUE.
        ELSE
          ISLR = .FALSE.
        ENDIF
        CALL ALLOC_LRB( BLR_U(I), K, KSVD, M, N, ISLR, 
     &             IFLAG, IERROR, KEEP8 )
        IF (IFLAG.LT.0) RETURN
        IF (LRFORM .NE. BLR_U(I)%LRFORM) THEN
          WRITE(*,*) "Internal error 2 in ALLOC_LRB",
     &    LRFORM, BLR_U(I)%LRFORM
        ENDIF
        IF (ISLR) THEN
          IF (K .GT. 0) THEN
            CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                     BLR_U(I)%Q(1,1), M*K, MPI_DOUBLE_PRECISION,
     &                     COMM, IERR )
            CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                     BLR_U(I)%R(1,1), N*K, MPI_DOUBLE_PRECISION,
     &                     COMM, IERR)
          ENDIF
        ELSE
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                     BLR_U(I)%Q(1,1), M*N, MPI_DOUBLE_PRECISION,
     &                     COMM, IERR)
        ENDIF
      ENDDO
      RETURN
      END SUBROUTINE DMUMPS_MPI_UNPACK_LR
