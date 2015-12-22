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
      RECURSIVE SUBROUTINE DMUMPS_TRAITER_MESSAGE(
     &    COMM_LOAD, ASS_IRECV,
     &    MSGSOU, MSGTAG, MSGLEN,
     &    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     &    IWPOS, IWPOSCB, IPTRLU,
     &    LRLU, LRLUS, N, IW, LIW, A, LA, PTRIST,
     &    PTLUST, PTRFAC,
     &    PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP,
     &    IFLAG, IERROR, COMM,
     &    NBPROCFILS,
     &    IPOOL, LPOOL, LEAF,
     &    NBFIN, MYID, SLAVEF,
     &
     &    root, OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &    FILS, PTRARW, PTRAIW,
     &    INTARR, DBLARR, ICNTL, KEEP,KEEP8,DKEEP, ND, FRERE,
     &    LPTRAR, NELT, FRTPTR, FRTELT,
     &
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE
     &    )
      USE DMUMPS_LOAD
      IMPLICIT NONE
      INCLUDE 'dmumps_root.h'
      INCLUDE 'mumps_headers.h'
      TYPE (DMUMPS_ROOT_STRUC) :: root
      INTEGER MSGSOU, MSGTAG, MSGLEN
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER BUFR( LBUFR )
      INTEGER KEEP(500), ICNTL( 40 )
      INTEGER(8) KEEP8(150)
      DOUBLE PRECISION    DKEEP(130)
      INTEGER(8) :: POSFAC, IPTRLU, LRLU, LRLUS, LA
      INTEGER IWPOS, IWPOSCB
      INTEGER N, LIW
      INTEGER IW( LIW )
      DOUBLE PRECISION A( LA )
      INTEGER(8) :: PTRFAC(KEEP(28))
      INTEGER(8) :: PTRAST(KEEP(28))
      INTEGER(8) :: PAMASTER(KEEP(28))
      INTEGER PTRIST(KEEP(28)), PTLUST(KEEP(28))
      INTEGER STEP(N), PIMASTER(KEEP(28))
      INTEGER COMP
      INTEGER NSTK_S(KEEP(28)), PROCNODE_STEPS( KEEP(28) )
      INTEGER NBPROCFILS( KEEP(28) )
      INTEGER IFLAG, IERROR, COMM
      INTEGER LPOOL, LEAF
      INTEGER IPOOL( LPOOL )
      INTEGER COMM_LOAD, ASS_IRECV
      INTEGER MYID, SLAVEF, NBFIN
      DOUBLE PRECISION OPASSW, OPELIW
      INTEGER NELT, LPTRAR
      INTEGER FRTPTR( N+1), FRTELT( NELT )
      INTEGER ITLOC( N+KEEP(253) ), FILS( N )
      DOUBLE PRECISION :: RHS_MUMPS(KEEP(255))
      INTEGER PTRARW( LPTRAR ), PTRAIW( LPTRAR )
      INTEGER ND( KEEP(28) ), FRERE( KEEP(28) )
      INTEGER ISTEP_TO_INIV2(KEEP(71)),
     &        TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
      INTEGER INTARR( max(1,KEEP(14)) )
      DOUBLE PRECISION DBLARR( max(1,KEEP(13)) )
      INTEGER INIV2, ISHIFT, IBEG
      INTEGER ISHIFT_HDR
      INTEGER MUMPS_PROCNODE, MUMPS_TYPENODE
      EXTERNAL MUMPS_PROCNODE, MUMPS_TYPENODE
      LOGICAL FLAG
      INTEGER MP, LP
      INTEGER TMP( 2 )
      INTEGER NBRECU, POSITION, INODE, ISON, IROOT
      INTEGER NSLAVES_PERE, NFRONT_PERE, NASS_PERE,
     &     LMAP, FPERE, NELIM,
     &     HDMAPLIG,NFS4FATHER,
     &     TOT_ROOT_SIZE, TOT_CONT_TO_RECV
      DOUBLE PRECISION FLOP1
      CHARACTER(LEN=35) :: SUBNAME
      INCLUDE 'mumps_tags.h'
      INCLUDE 'mpif.h'
      INTEGER :: IERR
      INTEGER :: STATUS(MPI_STATUS_SIZE)
      MP = ICNTL(2)
      LP = ICNTL(1)
      SUBNAME="??????"
      CALL DMUMPS_LOAD_RECV_MSGS(COMM_LOAD, KEEP)
      IF ( MSGTAG .EQ. RACINE ) THEN
          POSITION = 0
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, NBRECU,
     &     1, MPI_INTEGER, COMM, IERR)
          NBRECU = BUFR( 1 )
          NBFIN =  NBFIN - NBRECU
      ELSEIF ( MSGTAG .EQ. NOEUD ) THEN
          CALL DMUMPS_PROCESS_NODE( MYID, KEEP, KEEP8, DKEEP,
     &    BUFR, LBUFR, LBUFR_BYTES,
     &    IWPOS, IWPOSCB, IPTRLU,
     &    LRLU, LRLUS, N, IW, LIW, A, LA, PTRIST, PTRAST,
     &    STEP, PIMASTER, PAMASTER,
     &    NSTK_S, COMP, FPERE, FLAG, IFLAG, IERROR, COMM,
     &    ITLOC, RHS_MUMPS )
          SUBNAME="DMUMPS_PROCESS_NODE"
          IF ( IFLAG .LT. 0 ) GO TO 500
          IF ( FLAG ) THEN
            CALL DMUMPS_INSERT_POOL_N(N, IPOOL, LPOOL,
     &           PROCNODE_STEPS, SLAVEF, KEEP(28), KEEP(76),
     &           KEEP(80), KEEP(47), STEP, FPERE )
            IF (KEEP(47) .GE. 3) THEN
               CALL DMUMPS_LOAD_POOL_UPD_NEW_POOL(
     &              IPOOL, LPOOL,
     &              PROCNODE_STEPS, KEEP,KEEP8, SLAVEF, COMM_LOAD,
     &              MYID, STEP, N, ND, FILS )
            ENDIF
            CALL MUMPS_ESTIM_FLOPS( FPERE, N,
     &           PROCNODE_STEPS,SLAVEF,
     &           ND, FILS, FRERE, STEP, PIMASTER,
     &           KEEP(28), KEEP(50), KEEP(253), FLOP1,
     &           IW, LIW, KEEP(IXSZ) )
            IF (FPERE.NE.KEEP(20))
     &        CALL DMUMPS_LOAD_UPDATE(1,.FALSE.,FLOP1,KEEP,KEEP8)
          ENDIF
      ELSEIF ( MSGTAG .EQ. END_NIV2_LDLT ) THEN
          INODE = BUFR( 1 )
          CALL DMUMPS_INSERT_POOL_N(N, IPOOL, LPOOL,
     &         PROCNODE_STEPS, SLAVEF, KEEP(28), KEEP(76),
     &         KEEP(80), KEEP(47),
     &         STEP, -INODE )
          IF (KEEP(47) .GE. 3) THEN
             CALL DMUMPS_LOAD_POOL_UPD_NEW_POOL(
     &            IPOOL, LPOOL,
     &            PROCNODE_STEPS, KEEP,KEEP8, SLAVEF, COMM_LOAD,
     &            MYID, STEP, N, ND, FILS )
          ENDIF
      ELSEIF ( MSGTAG .EQ. TERREUR ) THEN
          IFLAG  = -001
          IERROR = MSGSOU
          GOTO 100
      ELSEIF ( MSGTAG .EQ. MAITRE_DESC_BANDE ) THEN
        CALL DMUMPS_PROCESS_DESC_BANDE( MYID,BUFR, LBUFR,
     &    LBUFR_BYTES, IWPOS,
     &    IWPOSCB,
     &    IPTRLU, LRLU, LRLUS, NBPROCFILS,
     &    N, IW, LIW, A, LA,
     &    PTRIST, PTRAST, STEP, PIMASTER, PAMASTER, COMP,
     &    KEEP, KEEP8, DKEEP, ITLOC, RHS_MUMPS, ISTEP_TO_INIV2, 
#if ! defined (NO_FDM_DESCBAND)
     &    -1,
#endif
     &    IFLAG, IERROR )
          SUBNAME="DMUMPS_PROCESS_DESC_BANDE"
        IF ( IFLAG .LT. 0 ) GO to 500
      ELSEIF ( MSGTAG .EQ. MAITRE2           ) THEN
        CALL DMUMPS_PROCESS_MASTER2(MYID,BUFR, LBUFR, LBUFR_BYTES,
     &    PROCNODE_STEPS, SLAVEF, IWPOS, IWPOSCB,
     &    IPTRLU, LRLU, LRLUS, N, IW, LIW, A, LA,
     &    PTRIST, PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP,
     &    IFLAG, IERROR, COMM, COMM_LOAD,
     &    IPOOL, LPOOL, LEAF,
     &    KEEP, KEEP8, DKEEP, ND, FILS, FRERE, ITLOC, RHS_MUMPS,
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE )
          SUBNAME="DMUMPS_PROCESS_MASTER2"
        IF ( IFLAG .LT. 0 ) GO to 500
      ELSEIF ( MSGTAG .EQ. BLOC_FACTO  .OR.
     &         MSGTAG .EQ. BLOC_FACTO_RELAY ) THEN
        CALL DMUMPS_PROCESS_BLOCFACTO( COMM_LOAD, ASS_IRECV,
     &   BUFR,  LBUFR, LBUFR_BYTES,
     &   PROCNODE_STEPS, MSGSOU,
     &   SLAVEF, IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, N, IW, LIW,
     &   A, LA, PTRIST, PTRAST, NSTK_S, NBPROCFILS,
     &   COMP, STEP, PIMASTER, PAMASTER, POSFAC,
     &   MYID, COMM , IFLAG, IERROR, NBFIN,
     &
     &    PTLUST, PTRFAC, root, OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &    FILS, PTRARW, PTRAIW, INTARR, DBLARR,
     &    ICNTL, KEEP,KEEP8,DKEEP, IPOOL, LPOOL, LEAF, ND, FRERE,
     &    LPTRAR, NELT, FRTPTR, FRTELT,
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE  
     &     )
      ELSEIF ( MSGTAG .EQ. BLOC_FACTO_SYM_SLAVE    ) THEN
        CALL DMUMPS_PROCESS_BLFAC_SLAVE( COMM_LOAD, ASS_IRECV,
     &   BUFR, LBUFR,
     &   LBUFR_BYTES, PROCNODE_STEPS, MSGSOU,
     &   SLAVEF, IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, N, IW, LIW,
     &   A, LA, PTRIST, PTRAST, NSTK_S, NBPROCFILS,
     &   COMP, STEP, PIMASTER, PAMASTER, POSFAC,
     &   MYID, COMM, IFLAG, IERROR, NBFIN,
     &
     &    PTLUST, PTRFAC, root, OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &    FILS, PTRARW, PTRAIW, INTARR, DBLARR,
     &    ICNTL, KEEP,KEEP8,DKEEP, IPOOL, LPOOL, LEAF, ND, FRERE,
     &    LPTRAR, NELT, FRTPTR, FRTELT,
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE  
     &     )
      ELSEIF ( MSGTAG .EQ. BLOC_FACTO_SYM    ) THEN
        CALL DMUMPS_PROCESS_SYM_BLOCFACTO( COMM_LOAD, ASS_IRECV,
     &   BUFR, LBUFR,
     &   LBUFR_BYTES, PROCNODE_STEPS, MSGSOU,
     &   SLAVEF, IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, N, IW, LIW,
     &   A, LA, PTRIST, PTRAST, NSTK_S, NBPROCFILS,
     &   COMP, STEP, PIMASTER, PAMASTER, POSFAC,
     &   MYID, COMM, IFLAG, IERROR, NBFIN,
     &
     &    PTLUST, PTRFAC, root, OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &    FILS, PTRARW, PTRAIW, INTARR, DBLARR,
     &    ICNTL,KEEP,KEEP8,DKEEP,IPOOL, LPOOL, LEAF, ND, FRERE,
     &    LPTRAR, NELT, FRTPTR, FRTELT,
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE
     & )
      ELSEIF ( MSGTAG .EQ. CONTRIB_TYPE2    ) THEN
        CALL DMUMPS_PROCESS_CONTRIB_TYPE2( COMM_LOAD, ASS_IRECV,
     &       MSGLEN, BUFR, LBUFR,
     &       LBUFR_BYTES, PROCNODE_STEPS,
     &       SLAVEF, IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, POSFAC,
     &       N, IW, LIW, A, LA, PTRIST,
     &       PTLUST, PTRFAC, PTRAST,
     &       STEP, PIMASTER, PAMASTER, NBPROCFILS, COMP, root,
     &       OPASSW, OPELIW, ITLOC, RHS_MUMPS, NSTK_S,
     &       FILS, PTRARW, PTRAIW, INTARR, DBLARR, NBFIN, MYID, COMM,
     &       ICNTL,KEEP,KEEP8,DKEEP,IFLAG, IERROR, IPOOL, LPOOL, LEAF,
     &       ND, FRERE, LPTRAR, NELT, FRTPTR, FRTELT,
     &       ISTEP_TO_INIV2, TAB_POS_IN_PERE
     &       )
        IF ( IFLAG .LT. 0 ) GO TO 100
      ELSEIF ( MSGTAG .EQ. MAPLIG            ) THEN
         HDMAPLIG = 7
         INODE        = BUFR( 1 )
         ISON         = BUFR( 2 )
         NSLAVES_PERE = BUFR( 3 )
         NFRONT_PERE  = BUFR( 4 )
         NASS_PERE    = BUFR( 5 )
         LMAP         = BUFR( 6 )
         NFS4FATHER   = BUFR( 7 )
         IF ( NSLAVES_PERE.NE.0 ) THEN
            INIV2 = ISTEP_TO_INIV2 ( STEP(INODE) )
            ISHIFT = NSLAVES_PERE+1
            TAB_POS_IN_PERE(1:NSLAVES_PERE+1, INIV2) =
     &           BUFR(HDMAPLIG+1:HDMAPLIG+1+NSLAVES_PERE)
            TAB_POS_IN_PERE(SLAVEF+2, INIV2) = NSLAVES_PERE
         ELSE
            ISHIFT = 0
         ENDIF
         IBEG = HDMAPLIG+1+ISHIFT
         CALL DMUMPS_MAPLIG( COMM_LOAD, ASS_IRECV,
     &    BUFR, LBUFR, LBUFR_BYTES,
     &    INODE, ISON, NSLAVES_PERE,
     &    BUFR(IBEG),
     &    NFRONT_PERE, NASS_PERE, NFS4FATHER,LMAP,
     &    BUFR(IBEG+NSLAVES_PERE),
     &    PROCNODE_STEPS, SLAVEF, POSFAC, IWPOS, IWPOSCB,
     &    IPTRLU, LRLU, LRLUS, N, IW, LIW, A, LA,
     &    PTRIST, PTLUST, PTRFAC, PTRAST, STEP, PIMASTER, PAMASTER,
     &    NSTK_S, COMP,
     &    IFLAG, IERROR, MYID, COMM, NBPROCFILS,
     &    IPOOL, LPOOL, LEAF, NBFIN, ICNTL, KEEP,KEEP8,DKEEP, root,
     &    OPASSW, OPELIW,
     &    ITLOC, RHS_MUMPS, FILS, PTRARW, PTRAIW, INTARR, DBLARR,
     &    ND, FRERE, LPTRAR, NELT, FRTPTR, FRTELT,
     &
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE
     &    )
         IF ( IFLAG .LT. 0 ) GO TO 100
      ELSE IF ( MSGTAG .EQ. ROOT_CONT_STATIC ) THEN
        CALL DMUMPS_PROCESS_CONTRIB_TYPE3(
     &        BUFR, LBUFR, LBUFR_BYTES,
     &        root, N, IW, LIW, A, LA, NBPROCFILS,
     &        LRLU, IPTRLU, IWPOS, IWPOSCB,
     &        PTRIST, PTLUST, PTRFAC, PTRAST,
     &        STEP, PIMASTER, PAMASTER,
     &        COMP, LRLUS, IPOOL, LPOOL, LEAF,
     &        FILS, MYID, PTRAIW, PTRARW, INTARR, DBLARR,
     &        KEEP, KEEP8, DKEEP, IFLAG, IERROR, COMM, COMM_LOAD,
     &        ITLOC, RHS_MUMPS,
     &        ND, PROCNODE_STEPS, SLAVEF)
        SUBNAME="DMUMPS_PROCESS_CONTRIB_TYPE3"
        IF ( IFLAG .LT. 0 ) GO TO 500
      ELSE IF ( MSGTAG .EQ. ROOT_NON_ELIM_CB ) THEN
        IROOT  = KEEP( 38 )
        MSGSOU = MUMPS_PROCNODE( PROCNODE_STEPS(STEP(IROOT)), 
     &           SLAVEF )
        IF ( PTLUST( STEP(IROOT)) .EQ. 0 ) THEN
          CALL MPI_RECV( TMP, 2 * KEEP(34), MPI_PACKED,
     &                   MSGSOU, ROOT_2SLAVE,
     &                   COMM, STATUS, IERR )
          CALL DMUMPS_PROCESS_ROOT2SLAVE( TMP( 1 ), TMP( 2 ),
     &    root,
     &    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     &    IWPOS, IWPOSCB, IPTRLU,
     &    LRLU, LRLUS, N, IW, LIW, A, LA, PTRIST,
     &    PTLUST, PTRFAC,
     &    PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP,
     &    IFLAG, IERROR, COMM, COMM_LOAD,
     &    NBPROCFILS,
     &    IPOOL, LPOOL, LEAF,
     &    NBFIN, MYID, SLAVEF,
     &
     &    OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &    FILS, PTRARW, PTRAIW,
     &    INTARR, DBLARR, ICNTL, KEEP,KEEP8, DKEEP,ND )
          SUBNAME="DMUMPS_PROCESS_ROOT2SLAVE"
          IF ( IFLAG .LT. 0 ) GOTO 500
        END IF
        CALL DMUMPS_PROCESS_CONTRIB_TYPE3(
     &       BUFR, LBUFR, LBUFR_BYTES,
     &       root, N, IW, LIW, A, LA, NBPROCFILS,
     &       LRLU, IPTRLU, IWPOS, IWPOSCB,
     &       PTRIST, PTLUST, PTRFAC, PTRAST, STEP, PIMASTER, PAMASTER,
     &       COMP, LRLUS, IPOOL, LPOOL, LEAF,
     &       FILS, MYID, PTRAIW, PTRARW, INTARR, DBLARR,
     &       KEEP, KEEP8, DKEEP, IFLAG, IERROR, COMM, COMM_LOAD,
     &       ITLOC, RHS_MUMPS,
     &       ND, PROCNODE_STEPS, SLAVEF )
          SUBNAME="DMUMPS_PROCESS_CONTRIB_TYPE3"
        IF ( IFLAG .LT. 0 ) GO TO 500
      ELSE IF ( MSGTAG .EQ. ROOT_2SON ) THEN
         ISON  = BUFR( 1 )
         NELIM = BUFR( 2 )
         CALL DMUMPS_PROCESS_ROOT2SON( COMM_LOAD, ASS_IRECV,
     &    ISON, NELIM, root,
     &    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     &    IWPOS, IWPOSCB, IPTRLU,
     &    LRLU, LRLUS, N, IW, LIW, A, LA, PTRIST,
     &    PTLUST, PTRFAC,
     &    PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP,
     &    IFLAG, IERROR, COMM,
     &    NBPROCFILS,
     &    IPOOL, LPOOL, LEAF,
     &    NBFIN, MYID, SLAVEF,
     &
     &    OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &    FILS, PTRARW, PTRAIW,
     &    INTARR,DBLARR,ICNTL,KEEP,KEEP8,DKEEP,ND, FRERE,
     &    LPTRAR, NELT, FRTPTR, FRTELT,
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE 
     &     )
          IF ( IFLAG .LT. 0 ) GO TO 100
         IF (MYID.NE.MUMPS_PROCNODE(PROCNODE_STEPS(STEP(ISON)), 
     &          SLAVEF)) THEN
          IF (KEEP(50).EQ.0) THEN
            ISHIFT_HDR = 6
          ELSE
            ISHIFT_HDR = 8
          ENDIF
          IF (IW(PTRIST(STEP(ISON))+ISHIFT_HDR+KEEP(IXSZ)).EQ.
     &                                 S_REC_CONTSTATIC) THEN
             IW(PTRIST(STEP(ISON))+ISHIFT_HDR+KEEP(IXSZ)) =
     &                                        S_ROOT2SON_CALLED
          ELSE
             CALL DMUMPS_FREE_BAND( N, ISON, PTRIST, PTRAST,
     &       IW, LIW, A, LA, LRLU, LRLUS, IWPOSCB,
     &       IPTRLU, STEP, MYID, KEEP,
     &       MUMPS_TYPENODE(PROCNODE_STEPS(STEP(ISON)),SLAVEF)
     &       )
          ENDIF
         ENDIF
      ELSE IF ( MSGTAG .EQ. ROOT_2SLAVE ) THEN
          TOT_ROOT_SIZE    = BUFR( 1 )
          TOT_CONT_TO_RECV = BUFR( 2 )
          CALL DMUMPS_PROCESS_ROOT2SLAVE( TOT_ROOT_SIZE,
     &    TOT_CONT_TO_RECV, root,
     &    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     &    IWPOS, IWPOSCB, IPTRLU,
     &    LRLU, LRLUS, N, IW, LIW, A, LA, PTRIST,
     &    PTLUST, PTRFAC,
     &    PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP,
     &    IFLAG, IERROR, COMM, COMM_LOAD,
     &    NBPROCFILS,
     &    IPOOL, LPOOL, LEAF,
     &    NBFIN, MYID, SLAVEF,
     &
     &    OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &    FILS, PTRARW, PTRAIW,
     &    INTARR, DBLARR, ICNTL, KEEP,KEEP8, DKEEP, ND )
          IF ( IFLAG .LT. 0 ) GO TO 100
      ELSE IF ( MSGTAG .EQ. ROOT_NELIM_INDICES ) THEN
         ISON         = BUFR( 1 )
         NELIM        = BUFR( 2 )
         NSLAVES_PERE = BUFR( 3 )
         CALL DMUMPS_PROCESS_RTNELIND( root,
     &    ISON, NELIM, NSLAVES_PERE, BUFR(4), BUFR(4+BUFR(2)),
     &    BUFR(4+2*BUFR(2)),
     &
     &    PROCNODE_STEPS,
     &    IWPOS, IWPOSCB, IPTRLU,
     &    LRLU, LRLUS, N, IW, LIW, A, LA, PTRIST,
     &    PTLUST, PTRFAC,
     &    PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S,
     &    ITLOC, RHS_MUMPS, COMP,
     &    IFLAG, IERROR,
     &    IPOOL, LPOOL, LEAF, MYID, SLAVEF,
     &    KEEP, KEEP8, DKEEP,
     &    COMM, COMM_LOAD, FILS, ND)
          SUBNAME="DMUMPS_PROCESS_RTNELIND"
         IF ( IFLAG .LT. 0 ) GO TO 500
      ELSE IF ( MSGTAG .EQ. UPDATE_LOAD ) THEN
         WRITE(*,*) "Internal error 3 in DMUMPS_TRAITER_MESSAGE"
         CALL MUMPS_ABORT()
      ELSE IF ( MSGTAG .EQ. TAG_DUMMY   ) THEN
      ELSE
         IF ( LP > 0 )
     &     WRITE(LP,*) MYID,
     &': Internal error, routine DMUMPS_TRAITER_MESSAGE.',MSGTAG
         IFLAG = -100
         IERROR= MSGTAG
         GOTO 500
      ENDIF
 100  CONTINUE
      RETURN
 500  CONTINUE
      IF ( ICNTL(1) .GT. 0 .AND. ICNTL(4).GE.1 ) THEN
        LP=ICNTL(1)
        IF (IFLAG.EQ.-9) THEN
         WRITE(LP,*) 'FAILURE, WORKSPACE TOO SMALL DURING ',SUBNAME
        ENDIF
        IF (IFLAG.EQ.-8) THEN
         WRITE(LP,*) 'FAILURE IN INTEGER ALLOCATION DURING ',SUBNAME
        ENDIF
        IF (IFLAG.EQ.-13) THEN
         WRITE(LP,*) 'FAILURE IN DYNAMIC ALLOCATION DURING ',SUBNAME
        ENDIF
      ENDIF
      CALL DMUMPS_BDC_ERROR( MYID, SLAVEF, COMM )
      RETURN
      END SUBROUTINE DMUMPS_TRAITER_MESSAGE
      RECURSIVE SUBROUTINE DMUMPS_RECV_AND_TREAT(
     &    COMM_LOAD, ASS_IRECV,
     &    STATUS,
     &    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     &    IWPOS, IWPOSCB, IPTRLU,
     &    LRLU, LRLUS, N, IW, LIW, A, LA, PTRIST,
     &    PTLUST, PTRFAC,
     &    PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP,
     &    IFLAG, IERROR, COMM,
     &    NBPROCFILS,
     &    IPOOL, LPOOL, LEAF,
     &    NBFIN, MYID, SLAVEF,
     &
     &    root, OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &    FILS, PTRARW, PTRAIW,
     &    INTARR, DBLARR, ICNTL, KEEP,KEEP8,DKEEP, ND, FRERE,
     &    LPTRAR, NELT, FRTPTR, FRTELT ,
     &
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE
     &    )
      IMPLICIT NONE
      INCLUDE 'dmumps_root.h'
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      TYPE (DMUMPS_ROOT_STRUC) :: root
      INTEGER :: STATUS(MPI_STATUS_SIZE)
      INTEGER KEEP(500), ICNTL(40)
      INTEGER(8) KEEP8(150)
      DOUBLE PRECISION       DKEEP(130)
      INTEGER COMM_LOAD, ASS_IRECV
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER BUFR( LBUFR )
      INTEGER(8) :: POSFAC, LA, IPTRLU, LRLU, LRLUS
      INTEGER IWPOS, IWPOSCB
      INTEGER N, LIW
      INTEGER IW( LIW )
      DOUBLE PRECISION A( LA )
      INTEGER(8) :: PTRFAC(KEEP(28))
      INTEGER(8) :: PTRAST(KEEP(28))
      INTEGER(8) :: PAMASTER(KEEP(28))
      INTEGER PTRIST( KEEP(28) ),
     &        PTLUST( KEEP(28) )
      INTEGER STEP(N), PIMASTER(KEEP(28))
      INTEGER COMP
      INTEGER NSTK_S(KEEP(28)), PROCNODE_STEPS( KEEP(28) )
      INTEGER NBPROCFILS( KEEP(28) )
      INTEGER IFLAG, IERROR, COMM
      INTEGER LPOOL, LEAF
      INTEGER IPOOL( LPOOL )
      INTEGER MYID, SLAVEF, NBFIN
      DOUBLE PRECISION OPASSW, OPELIW
      INTEGER NELT, LPTRAR
      INTEGER FRTPTR( N+1 ), FRTELT( NELT )
      INTEGER ITLOC( N+KEEP(253) ), FILS( N )
      DOUBLE PRECISION :: RHS_MUMPS(KEEP(255))
      INTEGER PTRARW( LPTRAR ), PTRAIW( LPTRAR )
      INTEGER ND( KEEP(28) ), FRERE( KEEP(28) )
      INTEGER ISTEP_TO_INIV2(KEEP(71)),
     &        TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
      INTEGER INTARR( max(1,KEEP(14)) )
      DOUBLE PRECISION DBLARR( max(1,KEEP(13)) )
      INTEGER MSGSOU, MSGTAG, MSGLEN, IERR
      MSGSOU = STATUS( MPI_SOURCE )
      MSGTAG = STATUS( MPI_TAG )
      CALL MPI_GET_COUNT( STATUS, MPI_PACKED, MSGLEN, IERR )
      IF ( MSGLEN .GT. LBUFR_BYTES ) THEN
        IFLAG  = -20
        IERROR = MSGLEN
         WRITE(*,*) ' RECEPTION BUF TOO SMALL, Msgtag/len=',
     &                MSGTAG,MSGLEN
        CALL DMUMPS_BDC_ERROR( MYID, SLAVEF, COMM )
        RETURN
       ENDIF
       CALL MPI_RECV( BUFR, LBUFR_BYTES, MPI_PACKED, MSGSOU,
     &                 MSGTAG,
     &                 COMM, STATUS, IERR )
       CALL DMUMPS_TRAITER_MESSAGE(
     &      COMM_LOAD, ASS_IRECV,
     &      MSGSOU, MSGTAG, MSGLEN, BUFR, LBUFR,
     &      LBUFR_BYTES,
     &      PROCNODE_STEPS, POSFAC,
     &      IWPOS, IWPOSCB, IPTRLU,
     &      LRLU, LRLUS, N, IW, LIW, A, LA, PTRIST,
     &      PTLUST, PTRFAC,
     &      PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP, IFLAG,
     &      IERROR, COMM,
     &      NBPROCFILS,
     &      IPOOL, LPOOL, LEAF, NBFIN, MYID, SLAVEF,
     &
     &      root, OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &      FILS, PTRARW, PTRAIW,
     &      INTARR, DBLARR, ICNTL, KEEP,KEEP8,DKEEP, ND, FRERE,
     &      LPTRAR, NELT, FRTPTR, FRTELT,
     &
     &      ISTEP_TO_INIV2, TAB_POS_IN_PERE
     &      )
      RETURN
      END SUBROUTINE DMUMPS_RECV_AND_TREAT
      RECURSIVE SUBROUTINE DMUMPS_TRY_RECVTREAT(
     &    COMM_LOAD, ASS_IRECV, BLOCKING, SET_IRECV,
     &    MESSAGE_RECEIVED, MSGSOU, MSGTAG,
     &    STATUS,
     &    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     &    IWPOS, IWPOSCB, IPTRLU,
     &    LRLU, LRLUS, N, IW, LIW, A, LA, PTRIST,
     &    PTLUST, PTRFAC,
     &    PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP,
     &    IFLAG, IERROR, COMM, NBPROCFILS,
     &    IPOOL, LPOOL, LEAF, NBFIN, MYID, SLAVEF,
     &
     &    root, OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &    FILS, PTRARW, PTRAIW,
     &    INTARR, DBLARR, ICNTL, KEEP,KEEP8,DKEEP, ND, FRERE,
     &    LPTRAR, NELT, FRTPTR, FRTELT,
     &
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &    STACK_RIGHT_AUTHORIZED
     & )
      USE DMUMPS_LOAD
      IMPLICIT NONE
      INCLUDE 'dmumps_root.h'
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      TYPE (DMUMPS_ROOT_STRUC) :: root
      INTEGER :: STATUS(MPI_STATUS_SIZE)
      LOGICAL, INTENT (IN)  :: BLOCKING
      LOGICAL, INTENT (IN)  :: SET_IRECV
      LOGICAL, INTENT (INOUT) :: MESSAGE_RECEIVED
      INTEGER, INTENT (IN) :: MSGSOU, MSGTAG
      INTEGER KEEP(500), ICNTL(40)
      INTEGER(8) KEEP8(150)
      DOUBLE PRECISION       DKEEP(130)
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER COMM_LOAD, ASS_IRECV
      INTEGER BUFR( LBUFR )
      INTEGER(8) :: LA, POSFAC, IPTRLU, LRLU, LRLUS
      INTEGER IWPOS, IWPOSCB
      INTEGER N, LIW
      INTEGER IW( LIW )
      DOUBLE PRECISION A( LA )
      INTEGER(8) :: PTRAST(KEEP(28))
      INTEGER(8) :: PTRFAC(KEEP(28))
      INTEGER(8) :: PAMASTER(KEEP(28))
      INTEGER PTRIST( KEEP(28) ),
     &        PTLUST(KEEP(28))
      INTEGER STEP(N),
     & PIMASTER(KEEP(28))
      INTEGER COMP
      INTEGER NSTK_S(KEEP(28)), PROCNODE_STEPS( KEEP(28) )
      INTEGER NBPROCFILS( KEEP(28) )
      INTEGER IFLAG, IERROR, COMM
      INTEGER LPOOL, LEAF
      INTEGER IPOOL( LPOOL )
      INTEGER MYID, SLAVEF, NBFIN
      DOUBLE PRECISION OPASSW, OPELIW
      INTEGER NELT, LPTRAR
      INTEGER FRTPTR( N+1 ), FRTELT( NELT )
      INTEGER ITLOC( N + KEEP(253) ), FILS( N )
      DOUBLE PRECISION :: RHS_MUMPS(KEEP(255))
      INTEGER PTRARW( LPTRAR ), PTRAIW( LPTRAR )
      INTEGER ND( KEEP(28) ), FRERE( KEEP(28) )
      INTEGER ISTEP_TO_INIV2(KEEP(71)),
     &        TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
      INTEGER INTARR( max(1,KEEP(14)) )
      DOUBLE PRECISION DBLARR( max(1,KEEP(13)) )
      LOGICAL, intent(in) :: STACK_RIGHT_AUTHORIZED
       LOGICAL FLAG, RIGHT_MESS, FLAGbis
       INTEGER LP, MSGSOU_LOC, MSGTAG_LOC, MSGLEN_LOC
       INTEGER IERR
       INTEGER :: STATUS_BIS(MPI_STATUS_SIZE)
       INTEGER, SAVE :: RECURS = 0
      CALL DMUMPS_LOAD_RECV_MSGS(COMM_LOAD, KEEP)
      IF ( .NOT. STACK_RIGHT_AUTHORIZED ) THEN
          RETURN
      ENDIF
      RECURS = RECURS + 1
      LP = ICNTL(1)
      IF (ICNTL(4).LT.1) LP=-1
      IF ( MESSAGE_RECEIVED ) THEN
        MSGSOU_LOC = MPI_ANY_SOURCE
        MSGTAG_LOC = MPI_ANY_TAG
        GOTO 250
      ENDIF
      IF ( ASS_IRECV .NE. MPI_REQUEST_NULL) THEN
        IF (KEEP(117).NE.0) THEN
         WRITE(*,*) "Problem of active IRECV with KEEP(117)=",KEEP(117)
         CALL MUMPS_ABORT()
        ENDIF
        RIGHT_MESS = .TRUE.
        IF (BLOCKING) THEN
          CALL MPI_WAIT(ASS_IRECV,
     &                STATUS, IERR)
          FLAG = .TRUE.
          IF ( ( (MSGSOU.NE.MPI_ANY_SOURCE) .OR.
     &      (MSGTAG.NE.MPI_ANY_TAG) )  ) THEN
            IF ( MSGSOU.NE.MPI_ANY_SOURCE) THEN
              RIGHT_MESS = MSGSOU.EQ.STATUS(MPI_SOURCE)
            ENDIF
            IF ( MSGTAG.NE.MPI_ANY_TAG) THEN
              RIGHT_MESS =
     &        ( (MSGTAG.EQ.STATUS(MPI_TAG)).AND.RIGHT_MESS )
            ENDIF
            IF (.NOT.RIGHT_MESS) THEN
              CALL MPI_PROBE(MSGSOU,MSGTAG,
     &           COMM, STATUS_BIS, IERR)
            ENDIF
          ENDIF
        ELSE 
          CALL MPI_TEST(ASS_IRECV,
     &             FLAG, STATUS, IERR)
        ENDIF
        IF (IERR.LT.0) THEN
          IFLAG = -20
          IF (LP.GT.0)
     &    write(LP,*) ' Error return from MPI_TEST ',
     &     IFLAG, ' in DMUMPS_TRY_RECVTREAT'
          CALL DMUMPS_BDC_ERROR( MYID, SLAVEF, COMM )
          RETURN
        ENDIF
        IF ( FLAG ) THEN
          MESSAGE_RECEIVED = .TRUE.
          MSGSOU_LOC = STATUS( MPI_SOURCE )
          MSGTAG_LOC = STATUS( MPI_TAG )
          CALL MPI_GET_COUNT( STATUS, MPI_PACKED, MSGLEN_LOC, IERR )
           IF (.NOT.RIGHT_MESS) RECURS = RECURS + 10
          CALL DMUMPS_TRAITER_MESSAGE( COMM_LOAD, ASS_IRECV,
     &      MSGSOU_LOC, MSGTAG_LOC, MSGLEN_LOC, BUFR, LBUFR,
     &      LBUFR_BYTES,
     &      PROCNODE_STEPS, POSFAC,
     &      IWPOS, IWPOSCB, IPTRLU,
     &      LRLU, LRLUS, N, IW, LIW, A, LA,
     &      PTRIST, PTLUST, PTRFAC,
     &      PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP, IFLAG,
     &      IERROR, COMM,
     &      NBPROCFILS,
     &      IPOOL, LPOOL, LEAF, NBFIN, MYID, SLAVEF,
     &
     &      root, OPASSW, OPELIW, ITLOC, RHS_MUMPS, FILS,
     &      PTRARW, PTRAIW,
     &      INTARR, DBLARR, ICNTL, KEEP,KEEP8,DKEEP, ND, FRERE,
     &      LPTRAR, NELT, FRTPTR, FRTELT,
     &      ISTEP_TO_INIV2, TAB_POS_IN_PERE
     &      )
          IF (.NOT.RIGHT_MESS) RECURS = RECURS - 10
          IF ( IFLAG .LT. 0 ) RETURN
          IF (.NOT.RIGHT_MESS) THEN
            IF (ASS_IRECV .NE. MPI_REQUEST_NULL) THEN
                CALL MUMPS_ABORT()
            ENDIF
            CALL MPI_IPROBE(MSGSOU,MSGTAG,
     &           COMM, FLAGbis, STATUS, IERR)
            IF (FLAGbis) THEN
               MSGSOU_LOC = STATUS( MPI_SOURCE )
               MSGTAG_LOC = STATUS( MPI_TAG )
               CALL DMUMPS_RECV_AND_TREAT( COMM_LOAD, ASS_IRECV,
     &            STATUS, BUFR, LBUFR,
     &            LBUFR_BYTES,
     &            PROCNODE_STEPS, POSFAC,
     &            IWPOS, IWPOSCB, IPTRLU,
     &            LRLU, LRLUS, N, IW, LIW, A, LA,
     &            PTRIST, PTLUST, PTRFAC,
     &            PTRAST, STEP, PIMASTER, PAMASTER,
     &            NSTK_S, COMP, IFLAG,
     &            IERROR, COMM,
     &            NBPROCFILS,
     &            IPOOL, LPOOL, LEAF,
     &            NBFIN, MYID, SLAVEF,
     &
     &            root, OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &            FILS, PTRARW, PTRAIW,
     &            INTARR, DBLARR, ICNTL,
     &            KEEP,KEEP8, DKEEP,ND, FRERE,
     &            LPTRAR, NELT, FRTPTR, FRTELT,
     &            ISTEP_TO_INIV2, TAB_POS_IN_PERE
     &            )
                  IF ( IFLAG .LT. 0 ) RETURN
            ENDIF
          ENDIF
       ENDIF
      ELSE
        IF (BLOCKING) THEN
           CALL MPI_PROBE(MSGSOU,MSGTAG,
     &           COMM, STATUS, IERR)
           FLAG = .TRUE.
        ELSE
           CALL MPI_IPROBE( MPI_ANY_SOURCE, MPI_ANY_TAG,
     &           COMM, FLAG, STATUS, IERR)
        ENDIF
        IF (FLAG) THEN
          MSGSOU_LOC = STATUS( MPI_SOURCE )
          MSGTAG_LOC = STATUS( MPI_TAG )
          MESSAGE_RECEIVED = .TRUE.
          CALL DMUMPS_RECV_AND_TREAT( COMM_LOAD, ASS_IRECV,
     &      STATUS, BUFR, LBUFR,
     &      LBUFR_BYTES,
     &      PROCNODE_STEPS, POSFAC,
     &      IWPOS, IWPOSCB, IPTRLU,
     &      LRLU, LRLUS, N, IW, LIW, A, LA,
     &      PTRIST, PTLUST, PTRFAC,
     &      PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP, IFLAG,
     &      IERROR, COMM,
     &      NBPROCFILS,
     &      IPOOL, LPOOL, LEAF, NBFIN, MYID, SLAVEF,
     &
     &      root, OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &      FILS, PTRARW, PTRAIW,
     &      INTARR, DBLARR, ICNTL, KEEP,KEEP8,DKEEP, ND, FRERE,
     &      LPTRAR, NELT, FRTPTR, FRTELT,
     &      ISTEP_TO_INIV2, TAB_POS_IN_PERE
     &   )
          IF ( IFLAG .LT. 0 ) RETURN
        ENDIF
      ENDIF
 250  CONTINUE
      RECURS  = RECURS - 1
      IF ( NBFIN .EQ. 0 ) RETURN
      IF ( RECURS .GT. 3 ) RETURN
      IF ( KEEP(36).EQ.1 .AND. SET_IRECV  .AND.
     &      (ASS_IRECV.EQ.MPI_REQUEST_NULL) .AND.
     &    MESSAGE_RECEIVED ) THEN
       CALL MPI_IRECV ( BUFR(1),
     &      LBUFR_BYTES, MPI_PACKED, MPI_ANY_SOURCE,
     &      MPI_ANY_TAG, COMM,
     &      ASS_IRECV, IERR )
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_TRY_RECVTREAT
      SUBROUTINE DMUMPS_MPI_CANCEL( INFO1,
     &    ASS_IRECV,
     &    BUFR, LBUFR, LBUFR_BYTES,
     &    COMM,
     &    MYID, SLAVEF)
      USE DMUMPS_COMM_BUFFER
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER ASS_IRECV
      INTEGER BUFR( LBUFR )
      INTEGER COMM
      INTEGER MYID, SLAVEF, INFO1, DEST
      INTEGER :: STATUS(MPI_STATUS_SIZE)
      LOGICAL NO_ACTIVE_IRECV
      INTEGER IERR, DUMMY
      INTRINSIC mod
      IF (SLAVEF .EQ. 1) RETURN
      IF (ASS_IRECV.EQ.MPI_REQUEST_NULL) THEN
        NO_ACTIVE_IRECV=.TRUE.
      ELSE
        CALL MPI_TEST(ASS_IRECV, NO_ACTIVE_IRECV,
     &                STATUS, IERR)
      ENDIF
      CALL MPI_BARRIER(COMM,IERR)
      DUMMY = 1
      DEST = mod(MYID+1, SLAVEF)
      CALL DMUMPS_BUF_SEND_1INT
     &    (DUMMY, DEST, TAG_DUMMY, COMM, IERR)
      IF (NO_ACTIVE_IRECV) THEN
        CALL MPI_RECV( BUFR, LBUFR,
     &             MPI_INTEGER, MPI_ANY_SOURCE,
     &             TAG_DUMMY, COMM, STATUS, IERR )
      ELSE
        CALL MPI_WAIT(ASS_IRECV,
     &                STATUS, IERR)
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_MPI_CANCEL
      SUBROUTINE DMUMPS_CLEAN_PENDING(
     &    INFO1, BUFR, LBUFR, LBUFR_BYTES,
     &    COMM_NODES, COMM_LOAD, SLAVEF, MP )
      USE DMUMPS_COMM_BUFFER
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER BUFR( LBUFR )
      INTEGER COMM_NODES, COMM_LOAD, SLAVEF, INFO1, MP
      INTEGER :: STATUS(MPI_STATUS_SIZE)
      LOGICAL FLAG, BUFFERS_EMPTY, BUFFERS_EMPTY_ON_ALL_PROCS
      INTEGER MSGSOU_LOC, MSGTAG_LOC, COMM_EFF
      INTEGER IERR
      INTEGER IBUF_EMPTY, IBUF_EMPTY_ON_ALL_PROCS
      IF (SLAVEF.EQ.1) RETURN
      BUFFERS_EMPTY_ON_ALL_PROCS = .FALSE.
 10   CONTINUE
      FLAG = .TRUE.
      DO WHILE ( FLAG )
        COMM_EFF = COMM_NODES
        CALL MPI_IPROBE(MPI_ANY_SOURCE,MPI_ANY_TAG,
     &       COMM_NODES, FLAG, STATUS, IERR)
        IF ( .NOT. FLAG ) THEN
          COMM_EFF = COMM_LOAD
          CALL MPI_IPROBE( MPI_ANY_SOURCE, MPI_ANY_TAG,
     &         COMM_LOAD, FLAG, STATUS, IERR)
        END IF
        IF (FLAG) THEN
            MSGSOU_LOC = STATUS( MPI_SOURCE )
            MSGTAG_LOC = STATUS( MPI_TAG )
               CALL MPI_RECV( BUFR, LBUFR_BYTES,
     &             MPI_PACKED, MSGSOU_LOC,
     &             MSGTAG_LOC, COMM_EFF, STATUS, IERR )
           ENDIF
         END DO
        IF (BUFFERS_EMPTY_ON_ALL_PROCS) THEN
        RETURN
        ENDIF
        CALL DMUMPS_BUF_ALL_EMPTY(BUFFERS_EMPTY)
        IF ( BUFFERS_EMPTY ) THEN
          IBUF_EMPTY = 0
        ELSE
          IBUF_EMPTY = 1
        ENDIF
        CALL MPI_ALLREDUCE(IBUF_EMPTY,
     &                     IBUF_EMPTY_ON_ALL_PROCS,
     &                     1, MPI_INTEGER, MPI_MAX,
     &                     COMM_NODES, IERR)
        IF ( IBUF_EMPTY_ON_ALL_PROCS == 0) THEN
          BUFFERS_EMPTY_ON_ALL_PROCS = .TRUE.
        ELSE
          BUFFERS_EMPTY_ON_ALL_PROCS = .FALSE.
        ENDIF
        GOTO 10
      END SUBROUTINE DMUMPS_CLEAN_PENDING