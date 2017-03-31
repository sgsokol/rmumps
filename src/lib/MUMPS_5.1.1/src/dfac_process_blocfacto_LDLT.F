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
      RECURSIVE SUBROUTINE DMUMPS_PROCESS_SYM_BLOCFACTO( 
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
     &    ICNTL,KEEP,KEEP8,DKEEP,IPOOL, LPOOL, LEAF, ND, FRERE_STEPS,
     &    LPTRAR, NELT, FRTPTR, FRTELT, 
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE  
     &               , LRGROUPS
     &    )
      USE DMUMPS_BUF
      USE DMUMPS_LOAD
      USE DMUMPS_OOC
      USE DMUMPS_LR_CORE
      USE DMUMPS_LR_TYPE
      USE DMUMPS_LR_STATS
      USE DMUMPS_FAC_LR
      USE DMUMPS_ANA_LR
      USE DMUMPS_LR_DATA_M
!$    USE OMP_LIB
      IMPLICIT NONE
      INCLUDE 'dmumps_root.h'
      INCLUDE 'mumps_headers.h'
      TYPE (DMUMPS_ROOT_STRUC) :: root
      INTEGER ICNTL( 40 ), KEEP( 500 )
      INTEGER(8) KEEP8(150)
      DOUBLE PRECISION    DKEEP(230)
      INTEGER COMM_LOAD, ASS_IRECV
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER BUFR( LBUFR )
      INTEGER N, SLAVEF, IWPOS, IWPOSCB, LIW
      INTEGER(8) IPTRLU, LRLU, LRLUS, LA, POSFAC
      INTEGER COMP
      INTEGER IFLAG, IERROR, NBFIN, MSGSOU
      INTEGER PROCNODE_STEPS(KEEP(28)), PTRIST(KEEP(28)),
     &        NSTK_S(KEEP(28))
      INTEGER(8) PTRAST(KEEP(28)), PTRFAC(KEEP(28)), PAMASTER(KEEP(28))
      INTEGER NBPROCFILS( KEEP(28) ), STEP(N), 
     & PIMASTER(KEEP(28))
      INTEGER IW( LIW )
      DOUBLE PRECISION A( LA )
      INTEGER, intent(in) :: LRGROUPS(N)
      INTEGER LPTRAR, NELT
      INTEGER FRTPTR( N+1 ), FRTELT( NELT )
      INTEGER COMM, MYID
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
      INTEGER PIVI
      INTEGER (8) POSPV1,POSPV2,OFFDAG,LPOS1
      INTEGER J2
      DOUBLE PRECISION MULT1,MULT2, A11, DETPIV, A22, A12
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER :: STATUS(MPI_STATUS_SIZE)
      INTEGER LP
      INTEGER INODE, POSITION, NPIV, IERR
      INTEGER NCOL, LD_BLOCFACTO
      INTEGER(8) LAELL, POSBLOCFACTO
      INTEGER(8) POSELT
      INTEGER IOLDPS, LCONT1, NASS1, NROW1, NCOL1, NPIV1
      INTEGER NSLAV1, HS, ISW, DEST
      INTEGER ICT11
      INTEGER(8) LPOS, LPOS2, DPOS, UPOS
      INTEGER (8) IPOS, KPOS
      INTEGER I, IPIV, FPERE, NSLAVES_TOT,
     &        NSLAVES_FOLLOW, NB_BLOC_FAC
      INTEGER IPOSK, JPOSK, NPIVSENT, Block, IROW, BLSIZE
      INTEGER allocok, TO_UPDATE_CPT_END
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: UIP21K
      INTEGER, DIMENSION(:), ALLOCATABLE :: LIST_SLAVES_FOLLOW
      LOGICAL LASTBL
      INTEGER SRC_DESCBAND
      LOGICAL BLOCKING, SET_IRECV, MESSAGE_RECEIVED
      DOUBLE PRECISION ONE,ALPHA
      PARAMETER (ONE = 1.0D0, ALPHA=-1.0D0)
      INTEGER(8) :: LAFAC
      INTEGER LIWFAC, STRAT, NextPivDummy
      LOGICAL LAST_CALL
      TYPE(IO_BLOCK) :: MonBloc
      INTEGER LRELAY_INFO
      LOGICAL COUNTER_WAS_HUGE
      INTEGER TO_UPDATE_CPT_RECUR
      LOGICAL :: SEND_LR 
      INTEGER :: XSIZE, CURRENT_BLR, NSLAVES_PREC, INFO_TMP(2)
      INTEGER :: SEND_LR_INT, NELIM, NB_BLR_LM, NB_BLR_LS,  
     &           MAXI_CLUSTER_LM, MAXI_CLUSTER_LS, MAXI_CLUSTER, 
     &           NPARTSASS, NPARTSCB, NPARTSCB_COL, NPARTSASS_COL, 
     &           NB_BLR_COL, MAXI_CLUSTER_COL
       INTEGER :: NPARTSASS_MASTER, IPANEL, NB_ACCESSES_INIT
      TYPE (LRB_TYPE), DIMENSION(:), ALLOCATABLE :: BLR_LM 
      TYPE (LRB_TYPE), DIMENSION(:), POINTER     :: BLR_LS
      INTEGER, POINTER, DIMENSION(:) :: BEGS_BLR_LM, BEGS_BLR_LS, 
     &                                  BEGS_BLR_COL
      LOGICAL KEEP_BEGS_BLR_LS, KEEP_BEGS_BLR_COL, KEEP_BLR_LS
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: WORK, TAU
      INTEGER, ALLOCATABLE, DIMENSION(:) :: JPVT
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: BLOCKLR
      INTEGER T1, T2, COUNT_RATE, LWORK
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: RWORK
      INTEGER :: OMP_NUM, MY_NUM
      INTEGER MUMPS_PROCNODE
      EXTERNAL MUMPS_PROCNODE
      LP = ICNTL(1)
      IF (ICNTL(4) .LE. 0) LP = -1
      POSITION = 0
      TO_UPDATE_CPT_END = -654321
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, INODE, 1,
     &                 MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, NPIV, 1,
     &                 MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, FPERE, 1,
     &                 MPI_INTEGER, COMM, IERR )
      LASTBL = (NPIV.LE.0)
      IF (LASTBL) THEN 
         NPIV = -NPIV
         CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, NSLAVES_TOT, 1,
     &                 MPI_INTEGER, COMM, IERR )
         CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, NB_BLOC_FAC, 1,
     &                 MPI_INTEGER, COMM, IERR )
      ENDIF
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, NCOL, 1,
     &                 MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, NELIM, 1,
     &                 MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, 
     &                 NPARTSASS_MASTER, 1,
     &                 MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, IPANEL,
     &                 1, MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, SEND_LR_INT, 1,
     &                 MPI_INTEGER, COMM, IERR )
      IF ( SEND_LR_INT .EQ. 1) THEN
        SEND_LR = .TRUE.
      ELSE
        SEND_LR = .FALSE.
      ENDIF
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, NSLAVES_TOT, 1,
     &                 MPI_INTEGER, COMM, IERR )
      XSIZE  = KEEP(IXSZ)
      KEEP_BEGS_BLR_LS  =.FALSE. 
      KEEP_BEGS_BLR_COL =.FALSE.
      KEEP_BLR_LS       =.FALSE.
      IF ( SEND_LR ) THEN
        LAELL = int(NPIV,8) * int(NPIV+NELIM,8)
       ELSE
         LAELL = int(NPIV,8) * int(NCOL,8)
      ENDIF
      IF ( NPIV.GT.0 ) THEN
       IF ( LRLU .LT. LAELL .OR. IWPOS + NPIV - 1 .GT. IWPOSCB ) THEN
        IF ( LRLUS .LT. LAELL ) THEN
          IFLAG = -9
          CALL MUMPS_SET_IERROR(LAELL-LRLUS, IERROR)
          IF (LP > 0 ) WRITE(LP,*) MYID,
     &": FAILURE IN DMUMPS_PROCESS_SYM_BLOCFACTO,
     & REAL WORKSPACE TOO SMALL"
          GOTO 700
        END IF
        CALL DMUMPS_COMPRE_NEW(N, KEEP(28), IW, LIW, A, LA,
     &       LRLU, IPTRLU,
     &       IWPOS, IWPOSCB, PTRIST, PTRAST,
     &       STEP, PIMASTER, PAMASTER, KEEP(216),LRLUS,
     &       KEEP(IXSZ),COMP,DKEEP(97),MYID)
        IF ( LRLU .NE. LRLUS ) THEN
             WRITE(*,*) 'PB compress DMUMPS_PROCESS_SYM_BLOCFACTO,",
     &       " LRLU,LRLUS='
     &       ,LRLU,LRLUS
             IFLAG = -9
             CALL MUMPS_SET_IERROR(LAELL-LRLUS,IERROR)
             GOTO 700
        END IF
        IF ( IWPOS + NPIV - 1 .GT. IWPOSCB ) THEN
          IF (LP > 0 ) WRITE(LP,*) MYID,
     &": FAILURE IN DMUMPS_PROCESS_SYM_BLOCFACTO,
     & INTEGER WORKSPACE TOO SMALL"
          IFLAG = -8
          IERROR = IWPOS + NPIV - 1 - IWPOSCB
          GOTO 700
        END IF
       END IF
       LRLU  = LRLU - LAELL
       LRLUS = LRLUS - LAELL
       KEEP8(70) = KEEP8(70) - LAELL
       KEEP8(71) = KEEP8(71) - LAELL
      ENDIF
      KEEP8(67) = min(LRLUS, KEEP8(67))
      KEEP8(68) = min(KEEP8(70), KEEP8(68))
      KEEP8(69) = min(KEEP8(71), KEEP8(69))
      POSBLOCFACTO = POSFAC
      POSFAC = POSFAC + LAELL
      CALL DMUMPS_LOAD_MEM_UPDATE(.FALSE.,.FALSE.,
     &                           LA-LRLUS,0_8,LAELL,KEEP,KEEP8,LRLUS)
      IF ( NPIV.EQ.0 ) THEN
        IPIV = 1 
        LD_BLOCFACTO = NPIV+NELIM
      ELSE
        IPIV = IWPOS
        IWPOS = IWPOS + NPIV
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                 IW( IPIV ), NPIV,
     &                 MPI_INTEGER, COMM, IERR )
      IF ( SEND_LR ) THEN
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                 A(POSBLOCFACTO), NPIV*(NPIV+NELIM),
     &                 MPI_DOUBLE_PRECISION,
     &                 COMM, IERR )
          LD_BLOCFACTO = NPIV+NELIM
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                 NB_BLR_LM, 1, MPI_INTEGER,
     &                 COMM, IERR )
          ALLOCATE(BLR_LM(max(NB_BLR_LM,1))) 
          ALLOCATE(BEGS_BLR_LM(NB_BLR_LM+2))
          CALL DMUMPS_MPI_UNPACK_LR(
     &          BUFR, LBUFR, LBUFR_BYTES, POSITION, NPIV, NELIM, 
     &          'V', BLR_LM, NB_BLR_LM, KEEP(470), 
     &          BEGS_BLR_LM(1), KEEP8, COMM, IERR, IFLAG, IERROR)
          IF (IFLAG.LT.0) GOTO 700
      ELSE
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &              A(POSBLOCFACTO), NPIV*NCOL, MPI_DOUBLE_PRECISION,
     &              COMM, IERR )
        LD_BLOCFACTO = NCOL
      ENDIF
      ENDIF 
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, 
     &                 LRELAY_INFO, 1,
     &                 MPI_INTEGER, COMM, IERR )
      IF (PTRIST(STEP( INODE )) .EQ. 0) THEN
        SRC_DESCBAND =
     &      MUMPS_PROCNODE( PROCNODE_STEPS(STEP(INODE)), SLAVEF )
          CALL DMUMPS_TREAT_DESCBAND( INODE, COMM_LOAD, ASS_IRECV,
     &      BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     &      IWPOS, IWPOSCB, IPTRLU,
     &      LRLU, LRLUS, N, IW, LIW, A, LA, PTRIST,
     &      PTLUST_S, PTRFAC,
     &      PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP,
     &      IFLAG, IERROR, COMM,
     &      NBPROCFILS,
     &      IPOOL, LPOOL, LEAF,
     &      NBFIN, MYID, SLAVEF,
     &
     &      root, OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &      FILS, PTRARW, PTRAIW,
     &      INTARR, DBLARR,ICNTL,KEEP,KEEP8,DKEEP,ND, FRERE_STEPS,
     &      LPTRAR, NELT, FRTPTR, FRTELT, 
     &      ISTEP_TO_INIV2, TAB_POS_IN_PERE, .TRUE. 
     &               , LRGROUPS
     &        )
          IF ( IFLAG .LT. 0 ) GOTO 600
      ENDIF
      IF ( IW( PTRIST(STEP(INODE)) + 3 + KEEP(IXSZ)) .EQ. 0 ) THEN
#if ! defined(NO_XXNBPR)
       CALL CHECK_EQUAL(NBPROCFILS(STEP(INODE)),
     &                  IW(PTRIST(STEP(INODE))+XXNBPR))
       DO WHILE ( IW(PTRIST(STEP(INODE)) + XXNBPR) .NE. 0)
#else
       DO WHILE ( NBPROCFILS( STEP(INODE)) .NE. 0 )
#endif
        BLOCKING = .TRUE.
        SET_IRECV=.FALSE.
        MESSAGE_RECEIVED = .FALSE.
        CALL DMUMPS_TRY_RECVTREAT( COMM_LOAD, ASS_IRECV,
     &    BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
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
     &    INTARR, DBLARR,ICNTL,KEEP,KEEP8,DKEEP,ND, FRERE_STEPS,
     &    LPTRAR, NELT, FRTPTR, FRTELT, 
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE, .TRUE. 
     &               , LRGROUPS
     &     )
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
     &       )
      IOLDPS = PTRIST(STEP(INODE))
      POSELT = PTRAST(STEP(INODE))
      LCONT1 = IW( IOLDPS + KEEP(IXSZ))
      NASS1  = IW( IOLDPS + 1 + KEEP(IXSZ))
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
      NROW1  = IW( IOLDPS + 2 + KEEP(IXSZ))
      NPIV1  = IW( IOLDPS + 3 + KEEP(IXSZ))
      NSLAV1 = IW( IOLDPS + 5 + KEEP(IXSZ))
      NSLAVES_FOLLOW = NSLAV1 - XTRA_SLAVES_SYM
      HS     = 6 + NSLAV1 + KEEP(IXSZ)
      NCOL1  = LCONT1 + NPIV1
      IF ( LASTBL ) THEN
        TO_UPDATE_CPT_END = ( NSLAVES_TOT - NSLAVES_FOLLOW - 1 ) * 
     &                       NB_BLOC_FAC
      END IF
      IF (NPIV.GT.0) THEN
        IF ( NPIV1 + NCOL .NE. NASS1 ) THEN
          WRITE(*,*) 'SymBLFC Error: NPIV1 + NCOL .NE. NASS1 :',
     &               NPIV1,NCOL,NASS1
          CALL MUMPS_ABORT()
        END IF
        ICT11 = IOLDPS+HS+NROW1+NPIV1 - 1
        DO I = 1, NPIV
          PIVI = abs(IW(IPIV+I-1))
          IF (PIVI.EQ.I) CYCLE
          ISW = IW(ICT11+I)
          IW(ICT11+I) = IW(ICT11+PIVI)
          IW(ICT11+PIVI) = ISW
          IPOS = POSELT + int(NPIV1 + I - 1,8)
          KPOS = POSELT + int(NPIV1 + PIVI - 1,8)
          CALL dswap(NROW1, A(IPOS), NCOL1, A(KPOS), NCOL1)
        ENDDO
        IF (.NOT.SEND_LR) THEN
        ALLOCATE( UIP21K( NPIV * NROW1 ), stat = allocok )
        IF ( allocok .GT. 0 ) THEN
            IF (LP > 0 ) WRITE(LP,*) MYID,
     &": ALLOCATION FAILURE FOR UIP21K IN DMUMPS_PROCESS_SYM_BLOCFACTO"
          IFLAG = -13
          IERROR = NPIV * NROW1
          GOTO 700
        END IF
        ELSE
         ALLOCATE( UIP21K( 1 ), stat = allocok )
         IF ( allocok .GT. 0 ) THEN
            IF (LP > 0 ) WRITE(LP,*) MYID,
     &": ALLOCATION FAILURE FOR UIP21K IN DMUMPS_PROCESS_SYM_BLOCFACTO"
          IFLAG = -13
          IERROR = NPIV * 1
          GOTO 700
        END IF
        ENDIF
        IF ( NSLAVES_FOLLOW .NE. 0 .and. NPIV .NE. 0 ) THEN
          ALLOCATE( LIST_SLAVES_FOLLOW ( NSLAVES_FOLLOW ),
     &            stat = allocok )
          IF ( allocok .GT. 0 ) THEN
            IF (LP > 0 ) WRITE(LP,*) MYID,
     &": ALLOCATION FAILURE FOR LIST_SLAVES_FOLLOW
     & IN DMUMPS_PROCESS_SYM_BLOCFACTO"
            IFLAG = -13
            IERROR = NSLAVES_FOLLOW
            GOTO 700
          END IF
          LIST_SLAVES_FOLLOW(1:NSLAVES_FOLLOW)=
     &    IW(IOLDPS+6+XTRA_SLAVES_SYM+KEEP(IXSZ):
     &     IOLDPS+5+XTRA_SLAVES_SYM+KEEP(IXSZ)+NSLAVES_FOLLOW)
        END IF
          CALL SYSTEM_CLOCK(T1)
          CALL dtrsm( 'L', 'U', 'T', 'U', NPIV, NROW1, ONE,
     &               A( POSBLOCFACTO ), LD_BLOCFACTO,
     &               A(POSELT+int(NPIV1,8)), NCOL1 )
          CALL SYSTEM_CLOCK(T2,COUNT_RATE)
          ACC_TRSM_TIME = ACC_TRSM_TIME + 
     &            DBLE(T2-T1)/DBLE(COUNT_RATE)
        IF (.NOT.SEND_LR) THEN
         LPOS = POSELT + int(NPIV1,8)
         UPOS = 1_8
         DO I = 1, NROW1
          UIP21K( UPOS: UPOS + int(NPIV-1,8) ) = 
     &                       A(LPOS: LPOS+int(NPIV-1,8))
          LPOS = LPOS + int(NCOL1,8)
          UPOS = UPOS + int(NPIV,8)
         END DO
        ENDIF
        LPOS = POSELT + int(NPIV1,8)
        DPOS = POSBLOCFACTO
        I = 1
        DO
          IF(I .GT. NPIV) EXIT
          IF(IW(IPIV+I-1) .GT. 0) THEN
          A11 = ONE/A(DPOS)
            CALL dscal( NROW1, A11, A(LPOS), NCOL1 )
            LPOS = LPOS + 1_8
            DPOS = DPOS + int(LD_BLOCFACTO + 1,8)
            I = I+1
          ELSE
            POSPV1 = DPOS
            POSPV2 = DPOS+ int(LD_BLOCFACTO + 1,8)
            OFFDAG = POSPV1+1_8
            A11 = A(POSPV1)
            A22 = A(POSPV2)
            A12 = A(OFFDAG)
            DETPIV = A11*A22 - A12**2
            A22 = A11/DETPIV
            A11 = A(POSPV2)/DETPIV
            A12 = -A12/DETPIV
            LPOS1 = LPOS
            DO J2 = 1,NROW1
               MULT1 = A11*A(LPOS1)+A12*A(LPOS1+1_8)
               MULT2 = A12*A(LPOS1)+A22*A(LPOS1+1_8)
               A(LPOS1) = MULT1
               A(LPOS1+1_8) = MULT2
               LPOS1 = LPOS1 + int(NCOL1,8)
            ENDDO
            LPOS = LPOS + 2_8
            DPOS = POSPV2 + int(LD_BLOCFACTO + 1,8)
            I = I+2
          ENDIF
        ENDDO
      ENDIF
      IF (SEND_LR) THEN
        NSLAVES_PREC = NSLAVES_TOT - NSLAVES_FOLLOW -1
      ENDIF
      IF (NPIV.GT.0) THEN
       IF (NROW1.LE.0) CALL MUMPS_ABORT()  
       IF (SEND_LR) THEN
        IF (NPIV1.NE.0) THEN
           CALL DMUMPS_BLR_RETRIEVE_BEGS_BLR_L (IW(IOLDPS+XXF), 
     &                  BEGS_BLR_LS)
           KEEP_BEGS_BLR_LS = .TRUE.  
           NB_BLR_LS = size(BEGS_BLR_LS) - 2
           NPARTSCB  = NB_BLR_LS
        ELSE
             CALL GET_CUT(IW(IOLDPS+HS:IOLDPS+HS+NROW1-1), 0,
     &                    NROW1, LRGROUPS, NPARTSCB, 
     &                    NPARTSASS, BEGS_BLR_LS)
              CALL REGROUPING2(BEGS_BLR_LS, NPARTSASS, 0, NPARTSCB,
     &                        NROW1-0, KEEP(488), .TRUE., KEEP(472))
             NB_BLR_LS = NPARTSCB
        ENDIF
        call MAX_CLUSTER(BEGS_BLR_LM,NB_BLR_LM+1,MAXI_CLUSTER_LM)
        call MAX_CLUSTER(BEGS_BLR_LS,NB_BLR_LS+1,MAXI_CLUSTER_LS)
        MAXI_CLUSTER=max(MAXI_CLUSTER_LS,MAXI_CLUSTER_LM,NPIV)
        IF (KEEP(489).EQ.1) THEN
         IF (NPIV1.EQ.0) THEN
          CALL GET_CUT(IW(IOLDPS+HS+NROW1:IOLDPS+HS+NROW1+NCOL1-1), 
     &                    NASS1,
     &                    NCOL1-NASS1, LRGROUPS, NPARTSCB_COL, 
     &                    NPARTSASS_COL, BEGS_BLR_COL)
          CALL REGROUPING2(BEGS_BLR_COL, NPARTSASS_COL, NASS1, 
     &                     NPARTSCB_COL,
     &                     NCOL1-NASS1, KEEP(488), .FALSE., KEEP(472))
          NB_BLR_COL = NPARTSCB_COL + NPARTSASS_COL
         ELSE
            CALL DMUMPS_BLR_RETRIEVE_BEGS_BLR_C (IW(IOLDPS+XXF), 
     &                  BEGS_BLR_COL,  NPARTSASS_MASTER)
            KEEP_BEGS_BLR_COL = .TRUE.  
            NB_BLR_COL   = size(BEGS_BLR_COL) - 1
            NPARTSCB_COL = NB_BLR_COL - NPARTSASS_MASTER 
         ENDIF
         CALL MAX_CLUSTER(BEGS_BLR_COL,NB_BLR_COL,MAXI_CLUSTER_COL)
         MAXI_CLUSTER = max(MAXI_CLUSTER,MAXI_CLUSTER_COL)
        ELSE
         NULLIFY(BEGS_BLR_COL)
        ENDIF
        IF (NPIV1.EQ.0)  THEN
          INFO_TMP(1) = IFLAG
          INFO_TMP(2) = IERROR
          NB_ACCESSES_INIT=0
            IF (NSLAVES_PREC.GT.0) THEN
              NB_ACCESSES_INIT=NSLAVES_PREC+1
            ENDIF
          CALL DMUMPS_BLR_INIT_FRONT(IW(IOLDPS+XXF), 
     &              .TRUE., .TRUE., .TRUE., NPARTSASS_MASTER, 
     &              BEGS_BLR_LS, BEGS_BLR_COL, NB_ACCESSES_INIT, 
     &              INFO_TMP)
         IFLAG  = INFO_TMP(1) 
         IERROR = INFO_TMP(2) 
         IF (IFLAG.LT.0) GOTO 700
        ENDIF
        LWORK = MAXI_CLUSTER*MAXI_CLUSTER
        OMP_NUM = 1
#if defined(BLR_MT)
!$      OMP_NUM = OMP_GET_MAX_THREADS()
#endif
        ALLOCATE(BLOCKLR(MAXI_CLUSTER, OMP_NUM*MAXI_CLUSTER),
     &            RWORK(2*MAXI_CLUSTER*OMP_NUM), 
     &            TAU(MAXI_CLUSTER*OMP_NUM),
     &            JPVT(MAXI_CLUSTER*OMP_NUM), 
     &            WORK(LWORK*OMP_NUM),
     &            stat=allocok)
        IF (allocok > 0 ) THEN
           IFLAG  = -13
           IERROR = OMP_NUM*(LWORK + MAXI_CLUSTER*(MAXI_CLUSTER+4))
           write(6,*) 'ERROR 2 allocate temporary BLR blocks during', 
     &         ' DMUMPS_PROCESS_SYM_BLOCFACTO', IERROR
           GOTO 700
        ENDIF
          CURRENT_BLR = 1
          ALLOCATE(BLR_LS(NB_BLR_LS)) 
          CALL SYSTEM_CLOCK(T1)
          MY_NUM=0
#if defined(BLR_MT)          
!$OMP PARALLEL PRIVATE(MY_NUM)
!$        MY_NUM = OMP_GET_THREAD_NUM()
#endif
          CALL DMUMPS_COMPRESS_PANEL(A, LA, POSELT, IFLAG, IERROR, 
     &         NCOL1,
     &         BEGS_BLR_LS, NB_BLR_LS+1, DKEEP(8), KEEP(473), BLR_LS, 
     &         CURRENT_BLR, 'V', WORK, TAU, JPVT, LWORK, RWORK,
     &         BLOCKLR, MAXI_CLUSTER, NELIM, 
     &         .TRUE.,          
     &         NPIV, NPIV1,
     &         2, KEEP(483), KEEP(470), KEEP8 
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
     &    ( .NOT. SEND_LR .OR. (NPIV.EQ.0) .OR.
     &    (KEEP(485).EQ.0) ) 
     &   ) THEN
        MonBloc%INODE = INODE
        MonBloc%MASTER = .FALSE.
        MonBloc%Typenode = 2
        MonBloc%NROW = NROW1  
        MonBloc%NCOL = NCOL1  
        MonBloc%NFS  = NASS1
        MonBloc%LastPiv = NPIV1 + NPIV 
        MonBloc%LastPanelWritten_L = -9999 
        MonBloc%LastPanelWritten_U = -9999 
        NULLIFY(MonBloc%INDICES)
        MonBloc%Last = LASTBL
        STRAT = STRAT_TRY_WRITE 
        NextPivDummy      = -8888 
        LIWFAC = IW(IOLDPS+XXI)
        CALL MUMPS_GETI8(LAFAC, IW(IOLDPS+XXR))
        LAST_CALL=.FALSE.
        CALL DMUMPS_OOC_IO_LU_PANEL( STRAT, TYPEF_L, A(POSELT),
     &       LAFAC, MonBloc, NextPivDummy, NextPivDummy,
     &       IW(IOLDPS), LIWFAC, MYID, KEEP8(31), IFLAG,LAST_CALL)
      ENDIF
      IF (NPIV.GT.0) THEN
       IF (SEND_LR) THEN
          IF (NELIM.GT.0) THEN
            LPOS2 = POSELT + int(NPIV1,8)
            UPOS = POSBLOCFACTO+int(NPIV,8)
            LPOS  = LPOS2 + int(NPIV,8)
            CALL dgemm('N','N', NELIM,NROW1,NPIV,ALPHA,
     &           A(UPOS),LD_BLOCFACTO,
     &           A(LPOS2),NCOL1,ONE, A(LPOS),NCOL1)
          ENDIF
#if defined(BLR_MT)          
!$OMP PARALLEL PRIVATE(MY_NUM)
!$        MY_NUM = OMP_GET_THREAD_NUM()
#endif          
          CALL DMUMPS_SLAVE_BLR_UPD_TRAIL_LDLT(A, LA, POSELT, 
     &        IFLAG, IERROR, NCOL1, NROW1,
     &        POSBLOCFACTO,
     &        LD_BLOCFACTO, 
     &        BEGS_BLR_LM, NB_BLR_LM+1, BLR_LM, NPIV1, 
     &        BEGS_BLR_LS, NB_BLR_LS+1, BLR_LS, 0, 
     &        CURRENT_BLR, CURRENT_BLR,   
     &        IW(IPIV), 
     &        BLOCKLR(1:MAXI_CLUSTER,MY_NUM*MAXI_CLUSTER+1),
     &        MAXI_CLUSTER,
     &        KEEP(481), DKEEP(8), KEEP(477) 
     &        )
#if defined(BLR_MT)          
!$OMP END PARALLEL
#endif          
          CALL SYSTEM_CLOCK(T2,COUNT_RATE)
          ACC_UPDT_TIME = ACC_UPDT_TIME + 
     &               DBLE(T2-T1)/DBLE(COUNT_RATE)
          CALL STATS_STORE_BLR_PANEL_MRY(BLR_LS,
     &               0, NPARTSCB, 'V', 2)
            IF (KEEP(485).NE.0) THEN
              CALL SYSTEM_CLOCK(T1)
              CALL DMUMPS_DECOMPRESS_PANEL(A, LA, POSELT, NCOL1,
     &          .FALSE.,   
     &          NPIV1+1,   
     &          1,         
     &          NB_BLR_LS+1, BLR_LS, 
     &          CURRENT_BLR, 'V', NCOL1, KEEP(470))
              CALL SYSTEM_CLOCK(T2,COUNT_RATE)
              ACC_PROMOTING_TIME = ACC_PROMOTING_TIME + 
     &               DBLE(T2-T1)/DBLE(COUNT_RATE)
              IF (KEEP(201).eq.1) THEN
               MonBloc%INODE = INODE
               MonBloc%MASTER = .FALSE.
               MonBloc%Typenode = 2
               MonBloc%NROW = NROW1  
               MonBloc%NCOL = NCOL1  
               MonBloc%NFS  = NASS1
               MonBloc%LastPiv = NPIV1 + NPIV 
               MonBloc%LastPanelWritten_L = -9999 
               MonBloc%LastPanelWritten_U = -9999 
               NULLIFY(MonBloc%INDICES)
               MonBloc%Last = LASTBL
               STRAT = STRAT_TRY_WRITE 
               NextPivDummy      = -8888 
               LIWFAC = IW(IOLDPS+XXI)
               CALL MUMPS_GETI8(LAFAC, IW(IOLDPS+XXR))
               LAST_CALL=.FALSE.
               CALL DMUMPS_OOC_IO_LU_PANEL( STRAT, TYPEF_L, A(POSELT),
     &           LAFAC, MonBloc, NextPivDummy, NextPivDummy,
     &           IW(IOLDPS), LIWFAC, MYID, KEEP8(31), IFLAG,LAST_CALL)
              ENDIF
            ENDIF
          CALL DEALLOC_BLR_PANEL (BLR_LM, NB_BLR_LM, KEEP8, .FALSE.)
          DEALLOCATE(BLR_LM)
          IF (NSLAVES_PREC.GT.0) THEN
            CALL DMUMPS_BLR_SAVE_PANEL_LORU(
     &          IW(IOLDPS+XXF),
     &          0,   
     &          IPANEL,BLR_LS)
            KEEP_BLR_LS = .TRUE.
          ENDIF
       ELSE 
        LPOS2 = POSELT + int(NPIV1,8)
        UPOS = POSBLOCFACTO+int(NPIV,8)
        LPOS  = LPOS2 + int(NPIV,8)
        CALL dgemm('N','N', NCOL-NPIV,NROW1,NPIV,ALPHA,A(UPOS),NCOL,
     &           A(LPOS2),NCOL1,ONE, A(LPOS),NCOL1)
        DPOS = POSELT + int(NCOL1 - NROW1,8)
        IF ( NROW1 .GT. KEEP(7) ) THEN
          BLSIZE = KEEP(8)
        ELSE
          BLSIZE = NROW1
        ENDIF
        IF ( NROW1 .GT. 0 ) THEN
          DO IROW = 1, NROW1, BLSIZE
            Block = min( BLSIZE, NROW1 - IROW + 1 )
            DPOS  = POSELT + int(NCOL1 - NROW1,8)
     &            + int( IROW - 1, 8 ) * int( NCOL1 + 1, 8 )
            LPOS2 = POSELT + int(NPIV1,8)
     &            + int( IROW - 1, 8 ) * int( NCOL1, 8 )
            UPOS  = int( IROW - 1, 8 ) * int(NPIV, 8) + 1_8
            DO I = 1, Block
              CALL dgemv( 'T', NPIV, Block-I+1, ALPHA,
     &                A( LPOS2 + int(I - 1,8) * int(NCOL1,8) ), NCOL1,
     &                UIP21K( UPOS + int(NPIV,8) * int( I - 1, 8 ) ),
     &                1, ONE, A(DPOS+int(NCOL1+1,8)*int(I-1,8)),NCOL1 )
            END DO
           IF ( NROW1-IROW+1-Block .ne. 0 )
     &     CALL dgemm( 'T', 'N', Block, NROW1-IROW+1-Block, NPIV, ALPHA,
     &             UIP21K( UPOS ), NPIV,
     &             A( LPOS2 + int(Block,8) * int(NCOL1,8) ), NCOL1, ONE,
     &             A( DPOS + int(Block,8) * int(NCOL1,8) ), NCOL1 )
          ENDDO
        ENDIF
        ENDIF
        FLOP1 = dble(NROW1) * dble(NPIV) *
     &           dble( 2 * NCOL  - NPIV + NROW1 +1 )
        FLOP1 = -FLOP1
        CALL DMUMPS_LOAD_UPDATE( 1, .FALSE., FLOP1, KEEP,KEEP8 )
      ENDIF 
      IW(IOLDPS+KEEP(IXSZ)) = IW(IOLDPS+KEEP(IXSZ)) - NPIV
      IW(IOLDPS + 3+KEEP(IXSZ)) = IW(IOLDPS+3+KEEP(IXSZ)) + NPIV
      IF (LASTBL) IW(IOLDPS+1+KEEP(IXSZ)) = IW(IOLDPS + 3+KEEP(IXSZ))
      IF ( .NOT. SEND_LR ) THEN
      LRLU  = LRLU + LAELL
      LRLUS = LRLUS + LAELL
      KEEP8(70) = KEEP8(70) + LAELL
      KEEP8(71) = KEEP8(71) + LAELL
      POSFAC = POSFAC - LAELL
      IWPOS = IWPOS - NPIV
      CALL DMUMPS_LOAD_MEM_UPDATE(.FALSE.,.FALSE.,
     &                           LA-LRLUS,0_8,-LAELL,KEEP,KEEP8,LRLUS)
      ENDIF
      IF ( NSLAVES_FOLLOW .NE. 0 .and. NPIV .NE. 0 ) THEN
         IPOSK = NPIV1 + 1
         JPOSK = NCOL1 - NROW1 + 1
           NPIVSENT = NPIV
           IERR = -1
           DO WHILE ( IERR .eq. -1 )
            CALL DMUMPS_BUF_SEND_BLFAC_SLAVE(
     &                    INODE, NPIVSENT, FPERE,
     &                    IPOSK, JPOSK,
     &                    UIP21K, NROW1,
     &                    NSLAVES_FOLLOW,
     &                    LIST_SLAVES_FOLLOW(1),
     &                    COMM, KEEP,
     &             SEND_LR, BLR_LS, IPANEL, 
     &             A, LA, POSBLOCFACTO, LD_BLOCFACTO,
     &             IW(IPIV), MAXI_CLUSTER,
     &                    IERR )
            IF (IERR .EQ. -1 ) THEN
              IOLDPS = PTRIST(STEP(INODE))
              IF ( IW(IOLDPS+6+KEEP(IXSZ)) .eq.
     &              huge(IW(IOLDPS+6+KEEP(IXSZ))) ) THEN
                    COUNTER_WAS_HUGE=.TRUE.
                    IW(IOLDPS+6+KEEP(IXSZ)) = 1
              ELSE
                    COUNTER_WAS_HUGE=.FALSE.
              ENDIF
              TO_UPDATE_CPT_RECUR =
     &                      ( NSLAVES_TOT - NSLAVES_FOLLOW - 1 ) *
     &                       (2*NASS1/KEEP(6))
              IW(IOLDPS+6+KEEP(IXSZ)) =
     &             IW(IOLDPS+6+KEEP(IXSZ)) - TO_UPDATE_CPT_RECUR - 10
              BLOCKING = .FALSE.
              SET_IRECV= .TRUE.
              MESSAGE_RECEIVED = .FALSE.
              CALL DMUMPS_TRY_RECVTREAT( COMM_LOAD, ASS_IRECV,
     &         BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     &         MPI_ANY_SOURCE, MPI_ANY_TAG,
     &         STATUS, 
     &         BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     &         IWPOS, IWPOSCB, IPTRLU,
     &         LRLU, LRLUS, N, IW, LIW, A, LA, PTRIST,
     &         PTLUST_S, PTRFAC,
     &         PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP,
     &         IFLAG, IERROR, COMM,
     &         NBPROCFILS, IPOOL, LPOOL, LEAF, NBFIN, MYID, SLAVEF,
     &         root, OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &         FILS, PTRARW, PTRAIW,
     &         INTARR, DBLARR,ICNTL,KEEP,KEEP8,DKEEP,ND, FRERE_STEPS,
     &         LPTRAR, NELT, FRTPTR, FRTELT, 
     &         ISTEP_TO_INIV2, TAB_POS_IN_PERE, .TRUE. 
     &               , LRGROUPS
     &           )
              IOLDPS = PTRIST(STEP(INODE))
              IW(IOLDPS+6+KEEP(IXSZ)) =
     &             IW(IOLDPS+6+KEEP(IXSZ)) + TO_UPDATE_CPT_RECUR + 10
              IF ( COUNTER_WAS_HUGE .AND.
     &             IW(IOLDPS+6+KEEP(IXSZ)).EQ.1 ) THEN
                IW(IOLDPS+6+KEEP(IXSZ)) = huge(IW(IOLDPS+6+KEEP(IXSZ)))
              ENDIF
              IF ( IFLAG .LT. 0 ) GOTO 600
            END IF
           END DO
#if defined(IBC_TEST)
           WRITE(*,*) MYID,":Send2slave worked"
#endif
           IF ( IERR .eq. -2 ) THEN
              IF (LP > 0 ) WRITE(LP,*) MYID,
     &": FAILURE, SEND BUFFER TOO SMALL DURING
     & DMUMPS_PROCESS_SYM_BLOCFACTO"
             WRITE(LP,*) "NPIV=", NPIV, "NROW1=",NROW1
             IFLAG = -17
             IERROR = 5 * KEEP(34) + NPIV * NROW1 * KEEP(35)
             GOTO 700
           END IF
           IF ( IERR .eq. -3 ) THEN
              IF (LP > 0 ) WRITE(LP,*) MYID,
     &": FAILURE, RECV BUFFER TOO SMALL DURING
     & DMUMPS_PROCESS_SYM_BLOCFACTO"
             IFLAG = -20
             IERROR = 5 * KEEP(34) + NPIV * NROW1 * KEEP(35)
             GOTO 700
           END IF
           DEALLOCATE(LIST_SLAVES_FOLLOW)
      END IF
      IF ( NPIV.GT. 0 .AND. SEND_LR ) THEN
        IF (NSLAVES_PREC.GT.0) THEN
          IOLDPS = PTRIST(STEP(INODE))
          CALL DMUMPS_BLR_DEC_AND_TRYFREE_L(IW(IOLDPS+XXF),IPANEL,
     &                       KEEP8, .TRUE.)
        ENDIF 
      LRLU  = LRLU + LAELL
      LRLUS = LRLUS + LAELL
      KEEP8(70) = KEEP8(70) + LAELL
      KEEP8(71) = KEEP8(71) + LAELL
      POSFAC = POSFAC - LAELL
      IWPOS = IWPOS - NPIV
      CALL DMUMPS_LOAD_MEM_UPDATE(.FALSE.,.FALSE.,
     &                           LA-LRLUS,0_8,-LAELL,KEEP,KEEP8,LRLUS)
      ENDIF
      IF ( NPIV .NE. 0 )  THEN
        IF (allocated(UIP21K)) DEALLOCATE( UIP21K )
      ENDIF
      IOLDPS = PTRIST(STEP(INODE))
      IF (LASTBL) THEN
         IF (SEND_LR) THEN
           CALL STATS_COMPUTE_FLOP_SLAVE_TYPE2(NROW1, NCOL1, NASS1,
     &             KEEP(50), INODE)
         ELSE
           CALL UPD_FLOP_FRFRONT_SLAVE(NROW1, NCOL1, NASS1,
     &             KEEP(50), INODE)
         ENDIF
         IF ( IW(IOLDPS+6+KEEP(IXSZ)).EQ.
     &     huge(IW(IOLDPS+6+KEEP(IXSZ))) ) THEN
           IW(IOLDPS+6+KEEP(IXSZ)) =  1
         ENDIF
         IW(IOLDPS+6+KEEP(IXSZ)) = IW(IOLDPS+6+KEEP(IXSZ))
     &                           - TO_UPDATE_CPT_END 
     &                           - 1 
         IF ( IW(IOLDPS+6+KEEP(IXSZ) ) .eq. 0
     &       .and. KEEP(50) .ne. 0 .and. NSLAVES_FOLLOW .eq. 0
     &       .and. NSLAVES_TOT.NE.1 ) THEN
          DEST = MUMPS_PROCNODE( PROCNODE_STEPS(STEP(INODE)), SLAVEF )
          CALL DMUMPS_BUF_SEND_1INT( INODE, DEST, END_NIV2_LDLT,
     &                              COMM, KEEP, IERR )
          IF ( IERR .LT. 0 ) THEN
            write(*,*) ' Internal error in PROCESS_SYM_BLOCFACTO.'
            IFLAG = -99
            GOTO 700
          END IF
        ENDIF
      END IF
        IF (IW(IOLDPS+6+KEEP(IXSZ)) .eq. 0 ) THEN 
         IF (SEND_LR) THEN
          IF (KEEP(489) .EQ. 1) THEN
          CALL DMUMPS_FAKE_COMPRESS_CB(A, LA, POSELT, NCOL1,
     &          BEGS_BLR_LS, NB_BLR_LS+1,
     &          BEGS_BLR_COL, NB_BLR_COL, NPARTSASS_MASTER,
     &          DKEEP(8), NASS1, NROW1,
     &          KEEP(50), WORK, TAU, JPVT, LWORK, RWORK,
     &          BLOCKLR, MAXI_CLUSTER, STEP_STATS(INODE), 2, 
     &          .TRUE., 0, KEEP(484)) 
          ENDIF
         ENDIF
          CALL DMUMPS_END_FACTO_SLAVE( COMM_LOAD, ASS_IRECV, 
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
     &    INTARR,DBLARR,ICNTL,KEEP,KEEP8,DKEEP,ND,FRERE_STEPS,
     &    LPTRAR, NELT, FRTPTR, FRTELT, 
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE
     &               , LRGROUPS
     &      )
        ENDIF
        IF (SEND_LR) THEN
          IF (allocated(RWORK))  DEALLOCATE(RWORK)
          IF (allocated(work)) DEALLOCATE(WORK)
          IF (allocated(TAU)) DEALLOCATE(TAU)
          IF (allocated(JPVT)) DEALLOCATE(JPVT)
          IF (allocated(BLOCKLR)) DEALLOCATE(BLOCKLR)
          IF (NPIV.GT.0) THEN
            IF (.NOT.KEEP_BEGS_BLR_LS) THEN
              IF (associated(BEGS_BLR_LS)) DEALLOCATE(BEGS_BLR_LS)
            ENDIF
            IF (.NOT.KEEP_BLR_LS) THEN
              CALL DEALLOC_BLR_PANEL (BLR_LS, NB_BLR_LS, KEEP8, .TRUE.)
              IF (associated(BLR_LS)) DEALLOCATE(BLR_LS)
            ENDIF
            IF (associated(BEGS_BLR_LM)) DEALLOCATE(BEGS_BLR_LM)
            IF (.NOT.KEEP_BEGS_BLR_COL) THEN
              IF (KEEP(489).EQ.1) THEN
                IF (associated(BEGS_BLR_COL)) THEN 
                  DEALLOCATE( BEGS_BLR_COL)
                ENDIF
              ENDIF
            ENDIF
          ENDIF
        ENDIF
 600  CONTINUE
#if defined(IBC_TEST)
      write(6,*) MYID,' :Exiting DMUMPS_PROCESS_SYM_BLOCFACTO for
     &INODE=', INODE
#endif
      RETURN
 700  CONTINUE
      CALL DMUMPS_BDC_ERROR( MYID, SLAVEF, COMM, KEEP )
      RETURN
      END SUBROUTINE DMUMPS_PROCESS_SYM_BLOCFACTO
