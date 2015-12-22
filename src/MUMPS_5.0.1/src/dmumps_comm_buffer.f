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
        MODULE DMUMPS_COMM_BUFFER
        PRIVATE
        PUBLIC :: DMUMPS_BUF_TRY_FREE_CB, DMUMPS_BUF_INIT,
     &   DMUMPS_BUF_INI_MYID,
     &   DMUMPS_BUF_ALLOC_CB ,       DMUMPS_BUF_DEALL_CB ,
     &   DMUMPS_BUF_ALLOC_SMALL_BUF, DMUMPS_BUF_DEALL_SMALL_BUF,
     &   DMUMPS_BUF_ALLOC_LOAD_BUFFER,DMUMPS_BUF_DEALL_LOAD_BUFFER,
     &   DMUMPS_BUF_SEND_CB,     DMUMPS_BUF_SEND_VCB,
     &   DMUMPS_BUF_SEND_1INT,       DMUMPS_BUF_SEND_DESC_BANDE,
     &   DMUMPS_BUF_SEND_MAPLIG, DMUMPS_BUF_SEND_MAITRE2,
     &   DMUMPS_BUF_SEND_CONTRIB_TYPE2,
     &   DMUMPS_BUF_SEND_BLOCFACTO, DMUMPS_BUF_SEND_BLFAC_SLAVE,
     &   DMUMPS_BUF_SEND_MASTER2SLAVE,
     &   DMUMPS_BUF_SEND_CONTRIB_TYPE3, DMUMPS_BUF_SEND_RTNELIND,
     &   DMUMPS_BUF_SEND_ROOT2SLAVE, DMUMPS_BUF_SEND_ROOT2SON,
     &   DMUMPS_BUF_SEND_BACKVEC,DMUMPS_BUF_SEND_UPDATE_LOAD, 
     &   DMUMPS_BUF_DIST_IRECV_SIZE,
     &   DMUMPS_BUF_BCAST_ARRAY, DMUMPS_BUF_ALL_EMPTY,
     &   DMUMPS_BUF_BROADCAST, DMUMPS_BUF_SEND_NOT_MSTR,
     &   DMUMPS_BUF_SEND_FILS ,DMUMPS_BUF_DEALL_MAX_ARRAY
     &   ,DMUMPS_BUF_MAX_ARRAY_MINSIZE
     &   ,DMUMPS_BUF_TEST
        INTEGER NEXT, REQ, CONTENT, OVHSIZE
        PARAMETER( NEXT = 0, REQ = 1, CONTENT = 2, OVHSIZE = 2 )
        INTEGER, SAVE :: SIZEofINT, SIZEofREAL, BUF_MYID
        TYPE DMUMPS_COMM_BUFFER_TYPE
          INTEGER LBUF, HEAD, TAIL,LBUF_INT, ILASTMSG
          INTEGER, DIMENSION(:),POINTER :: CONTENT
        END TYPE DMUMPS_COMM_BUFFER_TYPE
        TYPE ( DMUMPS_COMM_BUFFER_TYPE ), SAVE :: BUF_CB
        TYPE ( DMUMPS_COMM_BUFFER_TYPE ), SAVE :: BUF_SMALL
        TYPE ( DMUMPS_COMM_BUFFER_TYPE ), SAVE :: BUF_LOAD
        INTEGER, SAVE :: SIZE_RBUF_BYTES
        INTEGER BUF_LMAX_ARRAY
        DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: BUF_MAX_ARRAY
        PUBLIC :: BUF_LMAX_ARRAY, BUF_MAX_ARRAY
      CONTAINS
        SUBROUTINE DMUMPS_BUF_TRY_FREE_CB()
        CALL DMUMPS_BUF_TRY_FREE(BUF_CB)
        RETURN
        END SUBROUTINE DMUMPS_BUF_TRY_FREE_CB
        SUBROUTINE DMUMPS_BUF_TRY_FREE(B)
        IMPLICIT NONE
        TYPE ( DMUMPS_COMM_BUFFER_TYPE ) :: B
        INCLUDE 'mpif.h'
        LOGICAL :: FLAG
        INTEGER :: IERR
        INTEGER :: STATUS(MPI_STATUS_SIZE)
        IF ( B%HEAD .NE. B%TAIL ) THEN
 10       CONTINUE
          CALL MPI_TEST( B%CONTENT( B%HEAD + REQ ), FLAG, STATUS, IERR )
          IF ( FLAG ) THEN
            B%HEAD = B%CONTENT( B%HEAD + NEXT )
            IF ( B%HEAD .EQ. 0 ) B%HEAD = B%TAIL
            IF ( B%HEAD .NE. B%TAIL ) GOTO 10
          END IF
        END IF
        IF ( B%HEAD .EQ. B%TAIL ) THEN
          B%HEAD = 1
          B%TAIL = 1
          B%ILASTMSG = 1
        END iF
        RETURN
        END SUBROUTINE DMUMPS_BUF_TRY_FREE
        SUBROUTINE DMUMPS_BUF_INI_MYID( MYID )
        IMPLICIT NONE
        INTEGER MYID
        BUF_MYID  = MYID
        RETURN
        END SUBROUTINE DMUMPS_BUF_INI_MYID
        SUBROUTINE DMUMPS_BUF_INIT( IntSize, RealSize )
        IMPLICIT NONE
        INTEGER IntSize, RealSize
        SIZEofINT = IntSize
        SIZEofREAL = RealSize
        NULLIFY(BUF_CB  %CONTENT)
        NULLIFY(BUF_SMALL%CONTENT)
        NULLIFY(BUF_LOAD%CONTENT)
        BUF_CB%LBUF     = 0
        BUF_CB%LBUF_INT = 0
        BUF_CB%HEAD     = 1
        BUF_CB%TAIL     = 1
        BUF_CB%ILASTMSG = 1
        BUF_SMALL%LBUF     = 0
        BUF_SMALL%LBUF_INT = 0
        BUF_SMALL%HEAD     = 1
        BUF_SMALL%TAIL     = 1
        BUF_SMALL%ILASTMSG = 1
        BUF_LOAD%LBUF     = 0
        BUF_LOAD%LBUF_INT = 0
        BUF_LOAD%HEAD     = 1
        BUF_LOAD%TAIL     = 1
        BUF_LOAD%ILASTMSG = 1
        RETURN
        END SUBROUTINE DMUMPS_BUF_INIT
        SUBROUTINE DMUMPS_BUF_ALLOC_CB( SIZE, IERR )
        IMPLICIT NONE
        INTEGER SIZE, IERR
        CALL BUF_ALLOC( BUF_CB, SIZE, IERR )
        RETURN
        END SUBROUTINE DMUMPS_BUF_ALLOC_CB
        SUBROUTINE DMUMPS_BUF_ALLOC_SMALL_BUF( SIZE, IERR )
        IMPLICIT NONE
        INTEGER SIZE, IERR
        CALL BUF_ALLOC( BUF_SMALL, SIZE, IERR )
        RETURN
        END SUBROUTINE DMUMPS_BUF_ALLOC_SMALL_BUF
        SUBROUTINE DMUMPS_BUF_ALLOC_LOAD_BUFFER( SIZE, IERR )
        IMPLICIT NONE
        INTEGER SIZE, IERR
        CALL BUF_ALLOC( BUF_LOAD, SIZE, IERR )        
        RETURN
        END SUBROUTINE DMUMPS_BUF_ALLOC_LOAD_BUFFER
        SUBROUTINE DMUMPS_BUF_DEALL_LOAD_BUFFER( IERR )
        IMPLICIT NONE
        INTEGER IERR
        CALL BUF_DEALL( BUF_LOAD, IERR )
        RETURN
        END SUBROUTINE DMUMPS_BUF_DEALL_LOAD_BUFFER
        SUBROUTINE DMUMPS_BUF_DEALL_MAX_ARRAY()
        IMPLICIT NONE
        IF (allocated( BUF_MAX_ARRAY)) DEALLOCATE( BUF_MAX_ARRAY )
        RETURN
        END SUBROUTINE DMUMPS_BUF_DEALL_MAX_ARRAY
        SUBROUTINE DMUMPS_BUF_MAX_ARRAY_MINSIZE(NFS4FATHER,IERR)
        IMPLICIT NONE
        INTEGER IERR, NFS4FATHER
        IERR = 0
        IF (allocated( BUF_MAX_ARRAY)) THEN
          IF (BUF_LMAX_ARRAY .GE. NFS4FATHER) RETURN
          DEALLOCATE( BUF_MAX_ARRAY )
        ENDIF
        ALLOCATE(BUF_MAX_ARRAY(NFS4FATHER),stat=IERR)
        BUF_LMAX_ARRAY=NFS4FATHER
        RETURN
        END SUBROUTINE DMUMPS_BUF_MAX_ARRAY_MINSIZE
        SUBROUTINE DMUMPS_BUF_DEALL_CB( IERR )
        IMPLICIT NONE
        INTEGER IERR
        CALL BUF_DEALL( BUF_CB, IERR )
        RETURN
        END SUBROUTINE DMUMPS_BUF_DEALL_CB
        SUBROUTINE DMUMPS_BUF_DEALL_SMALL_BUF( IERR )
        IMPLICIT NONE
        INTEGER IERR
        CALL BUF_DEALL( BUF_SMALL, IERR )
        RETURN
        END SUBROUTINE DMUMPS_BUF_DEALL_SMALL_BUF
        SUBROUTINE BUF_ALLOC( BUF, SIZE, IERR )
        IMPLICIT NONE
        TYPE ( DMUMPS_COMM_BUFFER_TYPE ) :: BUF
        INTEGER SIZE, IERR
        IERR         = 0
        BUF%LBUF     = SIZE
        BUF%LBUF_INT = ( SIZE + SIZEofINT - 1 ) / SIZEofINT
        IF ( associated ( BUF%CONTENT ) ) DEALLOCATE( BUF%CONTENT )
        ALLOCATE( BUF%CONTENT( BUF%LBUF_INT ), stat = IERR )
        IF (IERR .NE. 0) THEN
          NULLIFY( BUF%CONTENT )
          IERR         = -1
          BUF%LBUF     =  0
          BUF%LBUF_INT =  0
        END IF
        BUF%HEAD     = 1
        BUF%TAIL     = 1
        BUF%ILASTMSG = 1
        RETURN
        END SUBROUTINE BUF_ALLOC
        SUBROUTINE BUF_DEALL( BUF, IERR )
        IMPLICIT NONE
        TYPE ( DMUMPS_COMM_BUFFER_TYPE ) :: BUF
        INCLUDE 'mpif.h'
        INTEGER :: IERR
        INTEGER :: STATUS(MPI_STATUS_SIZE)
        LOGICAL :: FLAG
        IF ( .NOT. associated ( BUF%CONTENT ) ) THEN
          BUF%HEAD     = 1
          BUF%LBUF     = 0
          BUF%LBUF_INT = 0
          BUF%TAIL     = 1
          BUF%ILASTMSG = 1
          RETURN
        END IF
        DO WHILE ( BUF%HEAD.NE.0 .AND. BUF%HEAD .NE. BUF%TAIL )
          CALL MPI_TEST(BUF%CONTENT( BUF%HEAD + REQ ), FLAG,
     &                  STATUS, IERR)
          IF ( .not. FLAG ) THEN
            WRITE(*,*) '** Warning: trying to cancel a request.'
            WRITE(*,*) '** This might be problematic'
            CALL MPI_CANCEL( BUF%CONTENT( BUF%HEAD + REQ ), IERR )
            CALL MPI_REQUEST_FREE( BUF%CONTENT( BUF%HEAD + REQ ), IERR )
          END IF
          BUF%HEAD = BUF%CONTENT( BUF%HEAD + NEXT )
        END DO
        DEALLOCATE( BUF%CONTENT )
        NULLIFY( BUF%CONTENT )
        BUF%LBUF     = 0
        BUF%LBUF_INT = 0
        BUF%HEAD     = 1
        BUF%TAIL     = 1
        BUF%ILASTMSG = 1
        RETURN
        END SUBROUTINE BUF_DEALL
        SUBROUTINE DMUMPS_BUF_SEND_CB( NBROWS_ALREADY_SENT,
     &                                INODE, FPERE, NFRONT, LCONT,
     &                                NASS, NPIV,
     &                                IWROW, IWCOL, A, COMPRESSCB,
     &                                DEST, TAG, COMM, IERR )
        IMPLICIT NONE
        INTEGER DEST, TAG, COMM, IERR
        INTEGER NBROWS_ALREADY_SENT
        INTEGER INODE, FPERE, NFRONT, LCONT, NASS, NPIV 
        INTEGER IWROW( LCONT ), IWCOL( LCONT )
        DOUBLE PRECISION A( * )
        LOGICAL COMPRESSCB
        INCLUDE 'mpif.h'
        INTEGER NBROWS_PACKET
        INTEGER POSITION, IREQ, IPOS, I, J1
        INTEGER SIZE1, SIZE2, SIZE_PACK, SIZE_AV, SIZE_AV_REALS
        INTEGER IZERO, IONE
        INTEGER SIZECB
        INTEGER LCONT_SENT
        INTEGER DEST2(1)
        PARAMETER( IZERO = 0, IONE = 1 )
        LOGICAL RECV_BUF_SMALLER_THAN_SEND
        DOUBLE PRECISION TMP
        DEST2(1) = DEST
        IERR = 0
        IF (NBROWS_ALREADY_SENT .EQ. 0) THEN
          CALL MPI_PACK_SIZE( 11 + LCONT + LCONT, MPI_INTEGER,
     &                        COMM, SIZE1,  IERR)
        ELSE
          CALL MPI_PACK_SIZE( 5, MPI_INTEGER, COMM, SIZE1, IERR)
        ENDIF
        CALL DMUMPS_BUF_SIZE_AVAILABLE( BUF_CB, SIZE_AV )
        IF ( SIZE_AV .LT. SIZE_RBUF_BYTES ) THEN
          RECV_BUF_SMALLER_THAN_SEND = .FALSE.
        ELSE
          SIZE_AV = SIZE_RBUF_BYTES
          RECV_BUF_SMALLER_THAN_SEND = .TRUE.
        ENDIF
        SIZE_AV_REALS = ( SIZE_AV - SIZE1 ) / SIZEofREAL
        IF (SIZE_AV_REALS < 0 ) THEN
          NBROWS_PACKET = 0
        ELSE
          IF (COMPRESSCB) THEN
            TMP=2.0D0*dble(NBROWS_ALREADY_SENT)+1.0D0
            NBROWS_PACKET = int(
     &                      ( sqrt( TMP * TMP
     &                        + 8.0D0 * dble(SIZE_AV_REALS)) - TMP )
     &                        / 2.0D0 )
          ELSE
            IF (LCONT.EQ.0) THEN
              NBROWS_PACKET = 0
            ELSE
              NBROWS_PACKET = SIZE_AV_REALS / LCONT
            ENDIF
          ENDIF
        ENDIF
 10     CONTINUE
        NBROWS_PACKET = max(0,
     &            min(NBROWS_PACKET, LCONT - NBROWS_ALREADY_SENT))
        IF (NBROWS_PACKET .EQ. 0 .AND. LCONT .NE. 0) THEN
          IF (RECV_BUF_SMALLER_THAN_SEND) THEN
            IERR = -3
            GOTO 100
          ELSE
            IERR = -1
            GOTO 100
          ENDIF
        ENDIF
        IF (COMPRESSCB) THEN
          SIZECB = (NBROWS_ALREADY_SENT*NBROWS_PACKET)+(NBROWS_PACKET
     &             *(NBROWS_PACKET+1))/2
        ELSE
          SIZECB = NBROWS_PACKET * LCONT
        ENDIF
        CALL MPI_PACK_SIZE( SIZECB, MPI_DOUBLE_PRECISION,
     &                    COMM, SIZE2,  IERR )
        SIZE_PACK = SIZE1 + SIZE2
        IF (SIZE_PACK .GT. SIZE_AV ) THEN
          NBROWS_PACKET = NBROWS_PACKET - 1
          IF (NBROWS_PACKET > 0) THEN
             GOTO 10
          ELSE
             IF (RECV_BUF_SMALLER_THAN_SEND) THEN
               IERR=-3
               GOTO 100
             ELSE
               IERR = -1
               GOTO 100
             ENDIF
          ENDIF
        ENDIF
        IF (NBROWS_PACKET + NBROWS_ALREADY_SENT.NE.LCONT .AND.
     &     SIZE_PACK  .LT. SIZE_RBUF_BYTES / 4
     &    .AND. 
     &    .NOT. RECV_BUF_SMALLER_THAN_SEND)
     &    THEN
            IERR = -1
            GOTO 100
        ENDIF
        CALL BUF_LOOK( BUF_CB, IPOS, IREQ, SIZE_PACK, IERR, 
     &                 IONE , DEST2
     &               )
        IF (IERR .EQ. -1 .OR. IERR .EQ. -2) THEN
          NBROWS_PACKET = NBROWS_PACKET - 1
          IF ( NBROWS_PACKET > 0 )  GOTO 10
        ENDIF
        IF ( IERR .LT. 0 ) GOTO 100
        POSITION = 0
        CALL MPI_PACK( INODE, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                        POSITION, COMM, IERR )
        CALL MPI_PACK( FPERE, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                        POSITION, COMM, IERR )
        IF (COMPRESSCB) THEN
          LCONT_SENT=-LCONT
        ELSE
          LCONT_SENT=LCONT
        ENDIF
        CALL MPI_PACK( LCONT_SENT, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                        POSITION, COMM, IERR )
        CALL MPI_PACK( NBROWS_ALREADY_SENT, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                        POSITION, COMM, IERR )
        CALL MPI_PACK( NBROWS_PACKET, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                        POSITION, COMM, IERR )
        IF (NBROWS_ALREADY_SENT == 0) THEN
          CALL MPI_PACK( LCONT, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                        POSITION, COMM, IERR )
          CALL MPI_PACK( NASS-NPIV, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                        POSITION, COMM, IERR )
          CALL MPI_PACK( LCONT , 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                        POSITION, COMM, IERR )
          CALL MPI_PACK( IZERO, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                        POSITION, COMM, IERR )
          CALL MPI_PACK( IONE,  1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                        POSITION, COMM, IERR )
          CALL MPI_PACK( IZERO, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                        POSITION, COMM, IERR )
          CALL MPI_PACK( IWROW, LCONT, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                        POSITION, COMM, IERR )
          CALL MPI_PACK( IWCOL, LCONT, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                        POSITION, COMM, IERR )
        ENDIF
        IF ( LCONT .NE. 0 ) THEN
          J1 = 1 + NBROWS_ALREADY_SENT * NFRONT
          IF (COMPRESSCB) THEN
           DO I = NBROWS_ALREADY_SENT+1,
     &            NBROWS_ALREADY_SENT+NBROWS_PACKET
            CALL MPI_PACK( A( J1 ), I, MPI_DOUBLE_PRECISION,
     &                        BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                        POSITION, COMM, IERR )
             J1 = J1 + NFRONT
           END DO
          ELSE
           DO I = NBROWS_ALREADY_SENT+1,
     &            NBROWS_ALREADY_SENT+NBROWS_PACKET
            CALL MPI_PACK( A( J1 ), LCONT, MPI_DOUBLE_PRECISION,
     &                        BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                        POSITION, COMM, IERR )
             J1 = J1 + NFRONT
           END DO
          ENDIF
        END IF
        CALL MPI_ISEND( BUF_CB%CONTENT( IPOS ), POSITION, MPI_PACKED,
     &                DEST, TAG, COMM, BUF_CB%CONTENT( IREQ ), IERR )
        IF ( SIZE_PACK .LT. POSITION ) THEN
          WRITE(*,*) 'Error Try_send_cb: SIZE, POSITION=',SIZE_PACK,
     &               POSITION
          CALL MUMPS_ABORT()
        END IF
        IF ( SIZE_PACK .NE. POSITION )
     &    CALL BUF_ADJUST( BUF_CB, POSITION )
        NBROWS_ALREADY_SENT = NBROWS_ALREADY_SENT + NBROWS_PACKET
        IF (NBROWS_ALREADY_SENT .NE. LCONT ) THEN
          IERR = -1
          RETURN
        ENDIF
 100    CONTINUE
        RETURN
        END SUBROUTINE DMUMPS_BUF_SEND_CB
        SUBROUTINE DMUMPS_BUF_SEND_MASTER2SLAVE( NRHS, INODE, IFATH,
     &             EFF_CB_SIZE, LD_CB, LD_PIV, NPIV, 
     &             JBDEB, JBFIN,
     &             CB, SOL,
     &             DEST, COMM, IERR )
        IMPLICIT NONE
        INTEGER NRHS, INODE, IFATH, EFF_CB_SIZE, LD_CB, LD_PIV, NPIV 
        INTEGER DEST, COMM, IERR, JBDEB, JBFIN
        DOUBLE PRECISION CB( LD_CB*(NRHS-1)+EFF_CB_SIZE )
        DOUBLE PRECISION SOL( max(1, LD_PIV*(NRHS-1)+NPIV) )
        INCLUDE 'mpif.h'
        INCLUDE 'mumps_tags.h'
        INTEGER SIZE, SIZE1, SIZE2, K
        INTEGER POSITION, IREQ, IPOS
        INTEGER IONE
        INTEGER DEST2(1)
        PARAMETER ( IONE=1 )
        DEST2(1) = DEST
        IERR = 0
        CALL MPI_PACK_SIZE( 6, MPI_INTEGER, COMM, SIZE1, IERR )
        CALL MPI_PACK_SIZE( NRHS * (EFF_CB_SIZE + NPIV),
     &                      MPI_DOUBLE_PRECISION, COMM,
     &                      SIZE2, IERR )
        SIZE = SIZE1 + SIZE2
        CALL BUF_LOOK( BUF_CB, IPOS, IREQ, SIZE, IERR, 
     &                 IONE , DEST2
     &               )
        IF ( IERR .LT. 0 ) THEN
           RETURN
        ENDIF
        POSITION = 0
        CALL MPI_PACK( INODE, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE,
     &                        POSITION, COMM, IERR )
        CALL MPI_PACK( IFATH, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE,
     &                        POSITION, COMM, IERR )
        CALL MPI_PACK( EFF_CB_SIZE  , 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE,
     &                        POSITION, COMM, IERR )
        CALL MPI_PACK( NPIV , 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE,
     &                        POSITION, COMM, IERR )
        CALL MPI_PACK( JBDEB , 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE,
     &                        POSITION, COMM, IERR )
        CALL MPI_PACK( JBFIN , 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE,
     &                        POSITION, COMM, IERR )
        DO K = 1, NRHS
               CALL MPI_PACK( CB ( 1 + LD_CB * (K-1) ),
     &                        EFF_CB_SIZE, MPI_DOUBLE_PRECISION,
     &                        BUF_CB%CONTENT( IPOS ), SIZE,
     &                        POSITION, COMM, IERR )
        END DO
        IF ( NPIV .GT. 0 ) THEN
          DO K=1, NRHS
          CALL MPI_PACK( SOL(1+LD_PIV*(K-1)),
     &                         NPIV, MPI_DOUBLE_PRECISION,
     &                         BUF_CB%CONTENT( IPOS ), SIZE,
     &                         POSITION, COMM, IERR )
          ENDDO
        END IF
        CALL MPI_ISEND( BUF_CB%CONTENT( IPOS ), POSITION, MPI_PACKED,
     &                  DEST, Master2Slave, COMM,
     &                  BUF_CB%CONTENT( IREQ ), IERR )
        IF ( SIZE .LT. POSITION ) THEN
          WRITE(*,*) 'Try_send_master2slave: SIZE, POSITION = ',
     &               SIZE, POSITION
          CALL MUMPS_ABORT()
        END IF
        IF ( SIZE .NE. POSITION ) CALL BUF_ADJUST( BUF_CB, POSITION )
        RETURN
        END SUBROUTINE DMUMPS_BUF_SEND_MASTER2SLAVE
        SUBROUTINE DMUMPS_BUF_SEND_VCB( NRHS, NODE1, NODE2, NCB, LDW,
     &             LONG,
     &             IW, W, JBDEB, JBFIN,
     &             DEST, TAG, COMM, IERR )
        IMPLICIT NONE
        INTEGER LDW, DEST, TAG, COMM, IERR
        INTEGER NRHS, NODE1, NODE2, NCB, LONG, JBDEB, JBFIN
        INTEGER IW( max( 1, LONG ) )
        DOUBLE PRECISION W( max( 1, LDW * NRHS ) )
        INCLUDE 'mpif.h'
        INTEGER POSITION, IREQ, IPOS
        INTEGER SIZE1, SIZE2, SIZE, K
        INTEGER IONE
        INTEGER DEST2(1)
        PARAMETER ( IONE=1 )
        DEST2(1)=DEST
        IERR = 0
        IF ( NODE2 .EQ. 0 ) THEN
         CALL MPI_PACK_SIZE( 4+LONG, MPI_INTEGER, COMM, SIZE1, IERR )
        ELSE
         CALL MPI_PACK_SIZE( 6+LONG, MPI_INTEGER, COMM, SIZE1, IERR )
        END IF
        SIZE2 = 0
        IF ( LONG .GT. 0 ) THEN
          CALL MPI_PACK_SIZE( NRHS*LONG, MPI_DOUBLE_PRECISION,
     &                        COMM, SIZE2, IERR )
        END IF
        SIZE = SIZE1 + SIZE2
        CALL BUF_LOOK( BUF_CB, IPOS, IREQ, SIZE, IERR, 
     &                 IONE , DEST2
     &               )
        IF ( IERR .LT. 0 ) THEN
           RETURN
        ENDIF
        POSITION = 0
        CALL MPI_PACK( NODE1, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE,
     &                        POSITION, COMM, IERR )
        IF ( NODE2 .NE. 0 ) THEN
          CALL MPI_PACK( NODE2, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE,
     &                        POSITION, COMM, IERR )
          CALL MPI_PACK( NCB, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE,
     &                        POSITION, COMM, IERR )
        ENDIF
        CALL MPI_PACK( JBDEB, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE,
     &                        POSITION, COMM, IERR )
        CALL MPI_PACK( JBFIN, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE,
     &                        POSITION, COMM, IERR )
        CALL MPI_PACK( LONG,  1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE,
     &                        POSITION, COMM, IERR )
        IF ( LONG .GT. 0 ) THEN
          CALL MPI_PACK( IW, LONG, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE,
     &                        POSITION, COMM, IERR )
          DO K=1, NRHS
          CALL MPI_PACK( W(1+(K-1)*LDW), LONG, MPI_DOUBLE_PRECISION,
     &                        BUF_CB%CONTENT( IPOS ), SIZE,
     &                        POSITION, COMM, IERR )
          END DO
        END IF
        CALL MPI_ISEND( BUF_CB%CONTENT( IPOS ), POSITION, MPI_PACKED,
     &                  DEST, TAG, COMM, BUF_CB%CONTENT( IREQ ), IERR )
        IF ( SIZE .NE. POSITION ) CALL BUF_ADJUST( BUF_CB, POSITION )
        RETURN
        END SUBROUTINE DMUMPS_BUF_SEND_VCB
        SUBROUTINE DMUMPS_BUF_SEND_1INT( I, DEST, TAG, COMM, IERR )
        IMPLICIT NONE
        INTEGER I
        INTEGER DEST, TAG, COMM, IERR
        INCLUDE 'mpif.h'
        INTEGER IPOS, IREQ, MSG_SIZE, POSITION
        INTEGER IONE
        INTEGER DEST2(1)
        PARAMETER ( IONE=1 )
        DEST2(1)=DEST
        IERR = 0
        CALL MPI_PACK_SIZE( 1, MPI_INTEGER,
     &                      COMM, MSG_SIZE, IERR )
        CALL BUF_LOOK( BUF_SMALL, IPOS, IREQ, MSG_SIZE, IERR, 
     &                 IONE , DEST2
     &               )
        IF ( IERR .LT. 0 ) THEN 
         write(6,*) ' Internal error in DMUMPS_BUF_SEND_1INT',
     &       ' Buf size (bytes)= ',BUF_SMALL%LBUF
         RETURN
        ENDIF
        POSITION=0
        CALL MPI_PACK( I, 1,
     &                 MPI_INTEGER, BUF_SMALL%CONTENT( IPOS ),
     &                 MSG_SIZE,
     &                 POSITION, COMM, IERR )
        CALL MPI_ISEND( BUF_SMALL%CONTENT(IPOS), MSG_SIZE,
     &                  MPI_PACKED, DEST, TAG, COMM,
     &                  BUF_SMALL%CONTENT( IREQ ), IERR )
        RETURN
        END SUBROUTINE DMUMPS_BUF_SEND_1INT
        SUBROUTINE DMUMPS_BUF_ALL_EMPTY(FLAG)
        LOGICAL FLAG
        LOGICAL FLAG1, FLAG2, FLAG3
        CALL DMUMPS_BUF_EMPTY( BUF_SMALL, FLAG1 )
        CALL DMUMPS_BUF_EMPTY( BUF_CB, FLAG2 )
        CALL DMUMPS_BUF_EMPTY( BUF_LOAD, FLAG3 )
        FLAG = FLAG1 .AND. FLAG2 .AND. FLAG3
        RETURN
        END SUBROUTINE DMUMPS_BUF_ALL_EMPTY
        SUBROUTINE DMUMPS_BUF_EMPTY( B, FLAG )
        TYPE ( DMUMPS_COMM_BUFFER_TYPE ) :: B
        LOGICAL :: FLAG
        INTEGER SIZE_AVAIL
        CALL DMUMPS_BUF_SIZE_AVAILABLE(B, SIZE_AVAIL)
        FLAG = ( B%HEAD == B%TAIL )
        RETURN
        END SUBROUTINE DMUMPS_BUF_EMPTY
        SUBROUTINE DMUMPS_BUF_SIZE_AVAILABLE( B, SIZE_AV )
        IMPLICIT NONE
        TYPE ( DMUMPS_COMM_BUFFER_TYPE ) :: B
        INTEGER SIZE_AV
        INCLUDE 'mpif.h'
        INTEGER :: IERR
        INTEGER :: STATUS(MPI_STATUS_SIZE)
        LOGICAL :: FLAG
        IF ( B%HEAD .NE. B%TAIL ) THEN
 10       CONTINUE
          CALL MPI_TEST( B%CONTENT( B%HEAD + REQ ), FLAG, STATUS, IERR )
          IF ( FLAG ) THEN
            B%HEAD = B%CONTENT( B%HEAD + NEXT )
            IF ( B%HEAD .EQ. 0 ) B%HEAD = B%TAIL
            IF ( B%HEAD .NE. B%TAIL ) GOTO 10
          END IF
        END IF
        IF ( B%HEAD .EQ. B%TAIL ) THEN
          B%HEAD = 1
          B%TAIL = 1
          B%ILASTMSG = 1
        END IF
        IF ( B%HEAD .LE. B%TAIL ) THEN
           SIZE_AV = max( B%LBUF_INT - B%TAIL, B%HEAD - 2 )
        ELSE
           SIZE_AV = B%HEAD - B%TAIL - 1
        END IF
        SIZE_AV = min(SIZE_AV - OVHSIZE, SIZE_AV)
        SIZE_AV = SIZE_AV * SIZEofINT
        RETURN
        END SUBROUTINE DMUMPS_BUF_SIZE_AVAILABLE
        SUBROUTINE DMUMPS_BUF_TEST()
        INTEGER :: IPOS, IREQ, IERR
        INTEGER, PARAMETER :: IONE=1
        INTEGER :: MSG_SIZE
        INTEGER :: DEST2(1)
        DEST2=-10
        MSG_SIZE=1
        CALL BUF_LOOK( BUF_CB, IPOS, IREQ, MSG_SIZE, IERR, 
     &                 IONE , DEST2,.TRUE.)
        RETURN
        END SUBROUTINE DMUMPS_BUF_TEST
        SUBROUTINE BUF_LOOK( B, IPOS, IREQ, MSG_SIZE, IERR, 
     &    NDEST , PDEST, TEST_ONLY
     &         )
        IMPLICIT NONE
        TYPE ( DMUMPS_COMM_BUFFER_TYPE ) :: B
        INTEGER, INTENT(IN)        :: MSG_SIZE
        INTEGER, INTENT(OUT)       :: IPOS, IREQ, IERR
        LOGICAL, INTENT(IN), OPTIONAL :: TEST_ONLY
        INTEGER NDEST
        INTEGER, INTENT(IN)        :: PDEST(max(1,NDEST))
        INCLUDE 'mpif.h'
        INTEGER :: MSG_SIZE_INT
        INTEGER :: IBUF
        LOGICAL :: FLAG
        INTEGER :: STATUS(MPI_STATUS_SIZE)
        IERR = 0
        IF ( B%HEAD .NE. B%TAIL ) THEN
 10       CONTINUE
          CALL MPI_TEST( B%CONTENT( B%HEAD + REQ ), FLAG, STATUS, IERR )
          IF ( FLAG ) THEN
            B%HEAD = B%CONTENT( B%HEAD + NEXT )
            IF ( B%HEAD .EQ. 0 ) B%HEAD = B%TAIL
            IF ( B%HEAD .NE. B%TAIL ) GOTO 10
          END IF
        END IF
        IF ( B%HEAD .EQ. B%TAIL ) THEN
          B%HEAD = 1
          B%TAIL = 1
          B%ILASTMSG = 1
        END iF
        MSG_SIZE_INT = ( MSG_SIZE + ( SIZEofINT - 1 ) ) / SIZEofINT
        MSG_SIZE_INT = MSG_SIZE_INT + OVHSIZE
        IF (present(TEST_ONLY)) RETURN
        FLAG = (     ( B%HEAD .LE. B%TAIL )
     &               .AND. (
     &                 ( MSG_SIZE_INT .LE. B%LBUF_INT - B%TAIL )
     &                 .OR. ( MSG_SIZE_INT .LE. B%HEAD - 2 ) ) )
     &         .OR.
     &               ( ( B%HEAD .GT. B%TAIL )
     &               .AND. ( MSG_SIZE_INT .LE. B%HEAD - B%TAIL - 1 ) )
        IF ( .NOT. FLAG
     &    ) THEN
          IERR = -1
          IF ( MSG_SIZE_INT .GT. B%LBUF_INT - 1 ) THEN
            IERR = -2
          ENDIF
          IPOS = -1
          IREQ = -1
          RETURN
        END IF
        IF ( B%HEAD .LE. B%TAIL ) THEN
          IF ( MSG_SIZE_INT .LE. B%LBUF_INT - B%TAIL + 1 ) THEN
            IBUF = B%TAIL
          ELSE IF ( MSG_SIZE_INT .LE. B%HEAD - 1 ) THEN
            IBUF = 1
          END IF
        ELSE
          IBUF = B%TAIL
        END IF
        B%CONTENT( B%ILASTMSG + NEXT ) = IBUF
        B%ILASTMSG = IBUF
        B%TAIL = IBUF + MSG_SIZE_INT
        B%CONTENT( IBUF + NEXT ) = 0
        IPOS = IBUF + CONTENT
        IREQ = IBUF + REQ
        RETURN
        END SUBROUTINE BUF_LOOK
        SUBROUTINE BUF_ADJUST( BUF, SIZE )
        IMPLICIT NONE
        TYPE ( DMUMPS_COMM_BUFFER_TYPE ) :: BUF
        INTEGER SIZE
        INTEGER SIZE_INT
        SIZE_INT = ( SIZE + SIZEofINT - 1 ) / SIZEofINT
        SIZE_INT = SIZE_INT + OVHSIZE
        BUF%TAIL = BUF%ILASTMSG + SIZE_INT
        RETURN
        END SUBROUTINE BUF_ADJUST
      SUBROUTINE DMUMPS_BUF_SEND_DESC_BANDE(
     &             INODE, NBPROCFILS, NLIG, ILIG, NCOL, ICOL,
     &             NASS, NSLAVES, LIST_SLAVES,
     &             DEST, IBC_SOURCE, NFRONT, COMM, IERR
     &)
      IMPLICIT NONE
        INTEGER COMM, IERR, NFRONT
        INTEGER INODE
        INTEGER NLIG, NCOL, NASS, NSLAVES
        INTEGER NBPROCFILS, DEST
        INTEGER ILIG( NLIG )
        INTEGER ICOL( NCOL )
        INTEGER, INTENT(IN) :: IBC_SOURCE
        INTEGER LIST_SLAVES( NSLAVES )
        INCLUDE 'mpif.h'
        INCLUDE 'mumps_tags.h'
        INTEGER SIZE_INT, SIZE_BYTES, POSITION, IPOS, IREQ
        INTEGER IONE
        INTEGER DEST2(1)
        PARAMETER ( IONE=1 )
        DEST2(1) = DEST
        IERR = 0
        SIZE_INT = ( 7 + NLIG + NCOL + NSLAVES + 1 )
        SIZE_BYTES = SIZE_INT * SIZEofINT
        IF (SIZE_INT.GT.SIZE_RBUF_BYTES ) THEN
         IERR = -2
         RETURN
        END IF
        CALL BUF_LOOK( BUF_CB, IPOS, IREQ, SIZE_BYTES, IERR, 
     &                 IONE , DEST2
     &               )
        IF ( IERR .LT. 0 ) THEN
           RETURN
        ENDIF
        POSITION = IPOS
        BUF_CB%CONTENT( POSITION ) = SIZE_INT
        POSITION = POSITION + 1
        BUF_CB%CONTENT( POSITION ) = INODE
        POSITION = POSITION + 1
        BUF_CB%CONTENT( POSITION ) = NBPROCFILS
        POSITION = POSITION + 1
        BUF_CB%CONTENT( POSITION ) = NLIG
        POSITION = POSITION + 1
        BUF_CB%CONTENT( POSITION ) = NCOL
        POSITION = POSITION + 1
        BUF_CB%CONTENT( POSITION ) = NASS
        POSITION = POSITION + 1
        BUF_CB%CONTENT( POSITION ) = NFRONT
        POSITION = POSITION + 1
        BUF_CB%CONTENT( POSITION ) = NSLAVES
        POSITION = POSITION + 1
        IF (NSLAVES.GT.0) THEN
         BUF_CB%CONTENT( POSITION: POSITION + NSLAVES - 1 ) = 
     &   LIST_SLAVES( 1: NSLAVES )
         POSITION = POSITION + NSLAVES
        ENDIF
        BUF_CB%CONTENT( POSITION:POSITION + NLIG - 1 ) = ILIG
        POSITION = POSITION + NLIG
        BUF_CB%CONTENT( POSITION:POSITION + NCOL - 1 ) = ICOL
        POSITION = POSITION + NCOL
        POSITION = POSITION - IPOS
        IF ( POSITION * SIZEofINT .NE. SIZE_BYTES ) THEN
          WRITE(*,*) 'Error in DMUMPS_BUF_SEND_DESC_BANDE :',
     &               ' wrong estimated size'
          CALL MUMPS_ABORT()
        END IF
        CALL MPI_ISEND( BUF_CB%CONTENT( IPOS ), SIZE_BYTES,
     &                  MPI_PACKED,
     &                  DEST, MAITRE_DESC_BANDE, COMM,
     &                  BUF_CB%CONTENT( IREQ ), IERR )
        RETURN
        END SUBROUTINE DMUMPS_BUF_SEND_DESC_BANDE
        SUBROUTINE DMUMPS_BUF_SEND_MAITRE2( NBROWS_ALREADY_SENT,
     &  IPERE, ISON, NROW,
     &  IROW, NCOL, ICOL, VAL, LDA, NELIM, TYPE_SON,
     &  NSLAVES, SLAVES, DEST, COMM, IERR, 
     & 
     &  SLAVEF, KEEP,KEEP8, INIV2, TAB_POS_IN_PERE )
        IMPLICIT NONE
        INTEGER NBROWS_ALREADY_SENT
        INTEGER LDA, NELIM, TYPE_SON
        INTEGER IPERE, ISON, NROW, NCOL, NSLAVES
        INTEGER IROW( NROW )
        INTEGER ICOL( NCOL )
        INTEGER SLAVES( NSLAVES )
        DOUBLE PRECISION VAL(LDA, *)
        INTEGER IPOS, IREQ, DEST, COMM, IERR
        INTEGER SLAVEF, KEEP(500), INIV2
        INTEGER(8) KEEP8(150)
        INTEGER TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
        INCLUDE 'mpif.h'
        INCLUDE 'mumps_tags.h'
        INTEGER SIZE1, SIZE2, SIZE3, SIZE_PACK, POSITION, I
        INTEGER NBROWS_PACKET, NCOL_SEND
        INTEGER SIZE_AV
        LOGICAL RECV_BUF_SMALLER_THAN_SEND
        INTEGER IONE
        INTEGER DEST2(1)
        PARAMETER ( IONE=1 )
        DEST2(1) = DEST
        IERR = 0
        IF ( NELIM .NE. NROW ) THEN
          WRITE(*,*) 'Error in TRY_SEND_MAITRE2:',NELIM, NROW
          CALL MUMPS_ABORT()
        END IF
        IF (NBROWS_ALREADY_SENT .EQ. 0) THEN
          CALL MPI_PACK_SIZE( NROW+NCOL+7+NSLAVES, MPI_INTEGER,
     &                      COMM, SIZE1, IERR )
          IF ( TYPE_SON .eq. 2 ) THEN
          CALL MPI_PACK_SIZE( NSLAVES+1, MPI_INTEGER,
     &                          COMM, SIZE3, IERR )
          ELSE
            SIZE3 = 0
          ENDIF
          SIZE1=SIZE1+SIZE3
        ELSE
          CALL MPI_PACK_SIZE(7, MPI_INTEGER,COMM,SIZE1,IERR)
        ENDIF
        IF ( KEEP(50).ne.0  .AND. TYPE_SON .eq. 2 ) THEN
          NCOL_SEND = NROW
        ELSE
          NCOL_SEND = NCOL
        ENDIF
        CALL DMUMPS_BUF_SIZE_AVAILABLE( BUF_CB, SIZE_AV )
        IF (SIZE_AV .LT. SIZE_RBUF_BYTES) THEN
          RECV_BUF_SMALLER_THAN_SEND = .FALSE.
        ELSE
          RECV_BUF_SMALLER_THAN_SEND = .TRUE.
          SIZE_AV = SIZE_RBUF_BYTES
        ENDIF
        IF (NROW .GT. 0 ) THEN 
         NBROWS_PACKET = (SIZE_AV - SIZE1) / NCOL_SEND / SIZEofREAL
         NBROWS_PACKET = min(NBROWS_PACKET, NROW - NBROWS_ALREADY_SENT)
         NBROWS_PACKET = max(NBROWS_PACKET, 0)
        ELSE
          NBROWS_PACKET =0
        ENDIF
        IF (NBROWS_PACKET .EQ. 0 .AND. NROW .NE. 0) THEN
          IF (RECV_BUF_SMALLER_THAN_SEND) THEN
              IERR=-3
              GOTO 100
          ELSE
              IERR=-1
              GOTO 100
          ENDIF
        ENDIF
 10     CONTINUE
        CALL MPI_PACK_SIZE( NBROWS_PACKET * NCOL_SEND,
     &           MPI_DOUBLE_PRECISION,
     &           COMM, SIZE2, IERR )
        SIZE_PACK = SIZE1 + SIZE2
        IF (SIZE_PACK .GT. SIZE_AV) THEN
          NBROWS_PACKET = NBROWS_PACKET - 1
          IF ( NBROWS_PACKET .GT. 0 ) THEN
            GOTO 10
          ELSE
            IF (RECV_BUF_SMALLER_THAN_SEND) THEN
                IERR = -3
                GOTO 100
            ELSE
                IERR = -1
                GOTO 100
            ENDIF
          ENDIF
        ENDIF
       IF (NBROWS_PACKET + NBROWS_ALREADY_SENT.NE.NROW .AND.
     &   SIZE_PACK - SIZE1  .LT. ( SIZE_RBUF_BYTES - SIZE1 ) / 2
     &  .AND. 
     &   .NOT. RECV_BUF_SMALLER_THAN_SEND)
     &   THEN
           IERR = -1
           GOTO 100
       ENDIF
        CALL BUF_LOOK( BUF_CB, IPOS, IREQ, SIZE_PACK, IERR, 
     &                 IONE , DEST2
     &               )
        IF ( IERR .LT. 0 ) THEN
          GOTO 100
        ENDIF
        POSITION = 0
        CALL MPI_PACK( IPERE, 1, MPI_INTEGER,
     &                 BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                 POSITION, COMM, IERR )
        CALL MPI_PACK( ISON,  1, MPI_INTEGER,
     &                 BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                 POSITION, COMM, IERR )
        CALL MPI_PACK( NSLAVES, 1, MPI_INTEGER,
     &                 BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                 POSITION, COMM, IERR )
        CALL MPI_PACK( NROW, 1, MPI_INTEGER,
     &                 BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                 POSITION, COMM, IERR )
        CALL MPI_PACK( NCOL, 1, MPI_INTEGER,
     &                 BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                 POSITION, COMM, IERR )
        CALL MPI_PACK( NBROWS_ALREADY_SENT, 1, MPI_INTEGER,
     &                 BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                 POSITION, COMM, IERR )
        CALL MPI_PACK( NBROWS_PACKET, 1, MPI_INTEGER,
     &                 BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                 POSITION, COMM, IERR )
        IF (NBROWS_ALREADY_SENT .EQ. 0) THEN
          IF (NSLAVES.GT.0) THEN
            CALL MPI_PACK( SLAVES, NSLAVES, MPI_INTEGER,
     &                BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                POSITION, COMM, IERR )
          ENDIF
          CALL MPI_PACK( IROW, NROW, MPI_INTEGER,
     &                 BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                 POSITION, COMM, IERR )
          CALL MPI_PACK( ICOL, NCOL, MPI_INTEGER,
     &                 BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                 POSITION, COMM, IERR )
          IF ( TYPE_SON .eq. 2 ) THEN
            CALL MPI_PACK( TAB_POS_IN_PERE(1,INIV2), NSLAVES+1, 
     &                 MPI_INTEGER,
     &                 BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                 POSITION, COMM, IERR )
          ENDIF
        ENDIF
        IF (NBROWS_PACKET.GE.1) THEN
          DO I=NBROWS_ALREADY_SENT+1,
     &                   NBROWS_ALREADY_SENT+NBROWS_PACKET
            CALL MPI_PACK( VAL(1,I), NCOL_SEND, 
     &               MPI_DOUBLE_PRECISION,
     &               BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &               POSITION, COMM, IERR )
          ENDDO
        ENDIF
        CALL MPI_ISEND( BUF_CB%CONTENT( IPOS ), POSITION, MPI_PACKED,
     &                  DEST, MAITRE2, COMM,
     &                  BUF_CB%CONTENT( IREQ ), IERR )
        IF ( SIZE_PACK .LT. POSITION ) THEN
          write(*,*) 'Try_send_maitre2, SIZE,POSITION=',
     &                SIZE_PACK,POSITION
          CALL MUMPS_ABORT()
        END IF
        IF ( SIZE_PACK .NE. POSITION )
     &    CALL BUF_ADJUST( BUF_CB, POSITION )
        NBROWS_ALREADY_SENT = NBROWS_ALREADY_SENT + NBROWS_PACKET
        IF ( NBROWS_ALREADY_SENT .NE. NROW ) THEN
          IERR = -1
        ENDIF
 100    CONTINUE
        RETURN
        END SUBROUTINE DMUMPS_BUF_SEND_MAITRE2
        SUBROUTINE DMUMPS_BUF_SEND_CONTRIB_TYPE2(NBROWS_ALREADY_SENT,
     &  DESC_IN_LU,
     &  IPERE, NFRONT_PERE, NASS_PERE, NFS4FATHER,
     &  NSLAVES_PERE,
     &  ISON, NBROW, LMAP, MAPROW, PERM, IW_CBSON, A_CBSON,
     &  ISLAVE, PDEST, PDEST_MASTER, COMM, IERR, 
     &  
     & KEEP,KEEP8, STEP, N, SLAVEF,
     & ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     & COMPRESSCB, KEEP253_LOC )
        IMPLICIT NONE
        INTEGER NBROWS_ALREADY_SENT
        INTEGER, INTENT (in) :: KEEP253_LOC
        INTEGER IPERE, ISON, NBROW 
        INTEGER PDEST, ISLAVE, COMM, IERR
        INTEGER PDEST_MASTER, NASS_PERE, NSLAVES_PERE,
     &       NFRONT_PERE, LMAP
        INTEGER MAPROW( LMAP ), PERM( max(1, NBROW ))
        INTEGER IW_CBSON( * )
        DOUBLE PRECISION A_CBSON( * )
        LOGICAL DESC_IN_LU, COMPRESSCB
       INTEGER   KEEP(500), N , SLAVEF
       INTEGER(8) KEEP8(150)
       INTEGER   STEP(N), 
     &          ISTEP_TO_INIV2(KEEP(71)), 
     &          TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER NFS4FATHER,SIZE3,PS1,NCA,LROW1
      INTEGER(8) :: ASIZE
      LOGICAL COMPUTE_MAX
      INTEGER NBROWS_PACKET 
      INTEGER MAX_ROW_LENGTH
      INTEGER LROW, NELIM
      INTEGER(8) :: SIZFR, ITMP8
      INTEGER NPIV, NFRONT, HS
      INTEGER SIZE_PACK, SIZE1, SIZE2, POSITION,I
      INTEGER SIZE_INTEGERS, B, SIZE_REALS, TMPSIZE, ONEorTWO, SIZE_AV
      INTEGER NBINT, L
      INTEGER(8) :: APOS, SHIFTCB_SON, LDA_SON8
      INTEGER IPOS_IN_SLAVE
      INTEGER STATE_SON
      INTEGER INDICE_PERE, NROW, IPOS, IREQ, NOSLA
      INTEGER IONE, J, THIS_ROW_LENGTH
      INTEGER SIZE_DESC_BANDE, DESC_BANDE_BYTES
      LOGICAL RECV_BUF_SMALLER_THAN_SEND
      LOGICAL NOT_ENOUGH_SPACE
      INTEGER PDEST2(1)
      PARAMETER ( IONE=1 )
      INCLUDE 'mumps_headers.h'
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO = 0.0D0)
      COMPUTE_MAX = (KEEP(219) .NE. 0) .AND.
     &              (KEEP(50) .EQ. 2) .AND.
     &              (PDEST.EQ.PDEST_MASTER)
      IF (NBROWS_ALREADY_SENT == 0) THEN 
        IF (COMPUTE_MAX) THEN
          CALL DMUMPS_BUF_MAX_ARRAY_MINSIZE(NFS4FATHER,IERR)
          IF (IERR .NE. 0) THEN
            IERR         = -4
            RETURN
          ENDIF
        ENDIF
      ENDIF
      PDEST2(1) = PDEST
      IERR   = 0
      LROW   = IW_CBSON( 1 + KEEP(IXSZ))
      NELIM  = IW_CBSON( 2 + KEEP(IXSZ))
      NPIV   = IW_CBSON( 4 + KEEP(IXSZ))
      IF ( NPIV .LT. 0 ) THEN
          NPIV = 0
      END IF
      NROW   = IW_CBSON( 3 + KEEP(IXSZ))
      NFRONT = LROW + NPIV
      HS     = 6 + IW_CBSON( 6 + KEEP(IXSZ)) + KEEP(IXSZ)
      CALL MUMPS_GETI8( SIZFR, IW_CBSON( 1 + XXR ) )
      STATE_SON = IW_CBSON(1+XXS)
      IF (STATE_SON .EQ. S_NOLCBCONTIG) THEN
               LDA_SON8    = int(LROW,8)
               SHIFTCB_SON = int(NPIV,8)*int(NROW,8)
      ELSE IF (STATE_SON .EQ. S_NOLCLEANED) THEN
               LDA_SON8    = int(LROW,8)
               SHIFTCB_SON = 0_8
      ELSE
               LDA_SON8     = int(NFRONT,8)
               SHIFTCB_SON = int(NPIV,8)
      ENDIF
      CALL DMUMPS_BUF_SIZE_AVAILABLE( BUF_CB, SIZE_AV )
      IF (PDEST .EQ. PDEST_MASTER) THEN
        SIZE_DESC_BANDE=0 
      ELSE
        SIZE_DESC_BANDE=(7+SLAVEF+KEEP(127)*2)
        SIZE_DESC_BANDE=SIZE_DESC_BANDE+int(dble(KEEP(12))*
     &                  dble(SIZE_DESC_BANDE)/100.0D0)
        SIZE_DESC_BANDE=max(SIZE_DESC_BANDE,
     &     7+NSLAVES_PERE+NFRONT_PERE+NFRONT_PERE-NASS_PERE)
      ENDIF
      DESC_BANDE_BYTES=SIZE_DESC_BANDE*SIZEofINT
      IF ( SIZE_AV .LT. SIZE_RBUF_BYTES-DESC_BANDE_BYTES ) THEN
        RECV_BUF_SMALLER_THAN_SEND = .FALSE.
      ELSE
        RECV_BUF_SMALLER_THAN_SEND = .TRUE.
        SIZE_AV = SIZE_RBUF_BYTES-DESC_BANDE_BYTES
      ENDIF
      SIZE1=0
      IF (NBROWS_ALREADY_SENT==0) THEN
          IF(COMPUTE_MAX) THEN
               CALL MPI_PACK_SIZE(1, MPI_INTEGER,
     &            COMM, PS1, IERR )
               IF(NFS4FATHER .GT. 0) THEN
                CALL MPI_PACK_SIZE( NFS4FATHER, MPI_DOUBLE_PRECISION,
     &             COMM, SIZE1, IERR )
               ENDIF
               SIZE1 = SIZE1+PS1
          ENDIF
      ENDIF
      IF (KEEP(50) .EQ. 0) THEN
        ONEorTWO = 1
      ELSE
        ONEorTWO = 2
      ENDIF
      IF (PDEST .EQ.PDEST_MASTER) THEN
        L = 0
      ELSE IF (KEEP(50) .EQ. 0) THEN
        L = LROW
      ELSE
        L = LROW + PERM(1) - LMAP + NBROWS_ALREADY_SENT - 1
        ONEorTWO=ONEorTWO+1
      ENDIF
      NBINT = 6 + L
      CALL MPI_PACK_SIZE( NBINT, MPI_INTEGER,
     &                    COMM, TMPSIZE, IERR )
      SIZE1 = SIZE1 + TMPSIZE
      SIZE_AV = SIZE_AV - SIZE1
      NOT_ENOUGH_SPACE=.FALSE.
      IF (SIZE_AV .LT.0 ) THEN
        NBROWS_PACKET = 0
        NOT_ENOUGH_SPACE=.TRUE.
      ELSE
        IF ( KEEP(50) .EQ. 0 ) THEN
          NBROWS_PACKET =
     &       SIZE_AV / ( ONEorTWO*SIZEofINT+LROW*SIZEofREAL)
        ELSE
          B = 2 * ONEorTWO + 
     &      ( 1 + 2 *  LROW + 2 * PERM(1) + 2 * NBROWS_ALREADY_SENT )
     &      * SIZEofREAL / SIZEofINT
          NBROWS_PACKET=int((dble(-B)+sqrt((dble(B)*dble(B))+
     &        dble(4)*dble(2)*dble(SIZE_AV)/dble(SIZEofINT) *
     &        dble(SIZEofREAL/SIZEofINT)))*
     &        dble(SIZEofINT) / dble(2) / dble(SIZEofREAL))
        ENDIF
      ENDIF
 10   CONTINUE
      NBROWS_PACKET = max( 0,
     &           min( NBROWS_PACKET, NBROW - NBROWS_ALREADY_SENT))
      NOT_ENOUGH_SPACE = NOT_ENOUGH_SPACE .OR.
     &                   (NBROWS_PACKET .EQ.0.AND. NBROW.NE.0)
      IF (NOT_ENOUGH_SPACE) THEN
        IF (RECV_BUF_SMALLER_THAN_SEND) THEN
          IERR = -3
          GOTO 100
        ELSE
          IERR = -1
          GOTO 100
        ENDIF
      ENDIF
      IF (KEEP(50).EQ.0) THEN
        MAX_ROW_LENGTH = -99999
        SIZE_REALS = NBROWS_PACKET * LROW
      ELSE
        SIZE_REALS = (  LROW + PERM(1) + NBROWS_ALREADY_SENT ) *
     &  NBROWS_PACKET + ( NBROWS_PACKET * ( NBROWS_PACKET + 1) ) / 2
        MAX_ROW_LENGTH = LROW+PERM(1)-LMAP+NBROWS_ALREADY_SENT
     &                 + NBROWS_PACKET-1
      ENDIF
      SIZE_INTEGERS = ONEorTWO* NBROWS_PACKET
      CALL MPI_PACK_SIZE( SIZE_REALS, MPI_DOUBLE_PRECISION,
     &                    COMM, SIZE2,  IERR)
      CALL MPI_PACK_SIZE( SIZE_INTEGERS, MPI_INTEGER,
     &                    COMM, SIZE3,  IERR)
      IF (SIZE2 + SIZE3 .GT. SIZE_AV ) THEN
         NBROWS_PACKET = NBROWS_PACKET -1
         IF (NBROWS_PACKET .GT. 0 ) THEN
           GOTO 10
         ELSE
           IF (RECV_BUF_SMALLER_THAN_SEND) THEN
             IERR = -3
             GOTO 100
           ELSE
             IERR = -1
             GOTO 100
           ENDIF
         ENDIF
      ENDIF
        SIZE_PACK = SIZE1 + SIZE2 + SIZE3
        IF (NBROWS_PACKET + NBROWS_ALREADY_SENT.NE.NBROW .AND.
     &       SIZE_PACK  .LT. SIZE_RBUF_BYTES / 4 .AND.
     &    .NOT. RECV_BUF_SMALLER_THAN_SEND)
     &    THEN
            IERR = -1
            GOTO 100
        ENDIF
        CALL BUF_LOOK( BUF_CB, IPOS, IREQ, SIZE_PACK, IERR, 
     &                 IONE , PDEST2
     &               )
        IF (IERR .EQ. -1 .OR. IERR.EQ. -2) THEN
          NBROWS_PACKET = NBROWS_PACKET - 1
          IF (NBROWS_PACKET > 0 ) GOTO 10
        ENDIF
        IF ( IERR .LT. 0 ) GOTO 100
          IF (SIZE_PACK.GT.SIZE_RBUF_BYTES ) THEN
             IERR = -3
             GOTO 100
          ENDIF
        POSITION = 0
        CALL MPI_PACK( IPERE, 1, MPI_INTEGER,
     &                 BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                 POSITION, COMM, IERR )
        CALL MPI_PACK( ISON, 1, MPI_INTEGER,
     &                 BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                 POSITION, COMM, IERR )
        CALL MPI_PACK( NBROW, 1, MPI_INTEGER,
     &                 BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                 POSITION, COMM, IERR )
        IF (KEEP(50)==0) THEN
        CALL MPI_PACK( LROW, 1, MPI_INTEGER,
     &                 BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                 POSITION, COMM, IERR )
        ELSE
        CALL MPI_PACK( MAX_ROW_LENGTH, 1, MPI_INTEGER,
     &                 BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                 POSITION, COMM, IERR )
        ENDIF
        CALL MPI_PACK( NBROWS_ALREADY_SENT, 1, MPI_INTEGER,
     &                 BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                 POSITION, COMM, IERR )
        CALL MPI_PACK( NBROWS_PACKET, 1, MPI_INTEGER,
     &                 BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                 POSITION, COMM, IERR )
        IF ( PDEST .NE. PDEST_MASTER ) THEN
          IF (KEEP(50)==0) THEN
          CALL MPI_PACK( IW_CBSON( HS + NROW +  NPIV + 1 ), LROW,
     &                 MPI_INTEGER,
     &                 BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                 POSITION, COMM, IERR )
          ELSE
           IF (MAX_ROW_LENGTH > 0) THEN
           CALL MPI_PACK( IW_CBSON( HS + NROW +  NPIV + 1 ),
     &                 MAX_ROW_LENGTH,
     &                 MPI_INTEGER,
     &                 BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                 POSITION, COMM, IERR )
           ENDIF
          ENDIF
        END IF
        DO J=NBROWS_ALREADY_SENT+1,NBROWS_ALREADY_SENT+NBROWS_PACKET
           I = PERM(J)
           INDICE_PERE=MAPROW(I)
           CALL MUMPS_BLOC2_GET_ISLAVE(
     &          KEEP,KEEP8, IPERE, STEP, N, SLAVEF,
     &          ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &
     &          NASS_PERE,
     &          NFRONT_PERE - NASS_PERE,
     &          NSLAVES_PERE,
     &          INDICE_PERE,
     &          NOSLA,
     &          IPOS_IN_SLAVE )
           INDICE_PERE = IPOS_IN_SLAVE
           CALL MPI_PACK( INDICE_PERE, 1, MPI_INTEGER,
     &          BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &          POSITION, COMM, IERR )
        ENDDO
        DO J=NBROWS_ALREADY_SENT+1,NBROWS_ALREADY_SENT+NBROWS_PACKET
           I = PERM(J)
           INDICE_PERE=MAPROW(I)
           CALL MUMPS_BLOC2_GET_ISLAVE(
     &          KEEP,KEEP8, IPERE, STEP, N, SLAVEF,
     &          ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &          
     &          NASS_PERE,
     &          NFRONT_PERE - NASS_PERE,
     &          NSLAVES_PERE,
     &          INDICE_PERE,
     &          NOSLA,
     &          IPOS_IN_SLAVE )
          IF (KEEP(50).ne.0) THEN
            THIS_ROW_LENGTH = LROW + I - LMAP
            CALL MPI_PACK( THIS_ROW_LENGTH, 1, MPI_INTEGER,
     &                      BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                      POSITION, COMM, IERR )
         ELSE
            THIS_ROW_LENGTH = LROW
         ENDIF
         IF (DESC_IN_LU) THEN 
            IF ( COMPRESSCB ) THEN
             IF (NELIM.EQ.0) THEN
               ITMP8 = int(I,8)
             ELSE
               ITMP8 = int(NELIM+I,8)
             ENDIF
             APOS = ITMP8 * (ITMP8-1_8) / 2_8 + 1_8
            ELSE
             APOS = int(I+NELIM-1, 8) * int(LROW,8) + 1_8
            ENDIF
         ELSE
            IF ( COMPRESSCB ) THEN
             IF ( LROW .EQ. NROW )  THEN
               ITMP8 = int(I,8)
               APOS  = ITMP8 * (ITMP8-1_8)/2_8 + 1_8
             ELSE
               ITMP8 = int(I + LROW - NROW,8)
               APOS  = ITMP8 * (ITMP8-1_8)/2_8 + 1_8 -
     &                 int(LROW - NROW, 8) * int(LROW-NROW+1,8) / 2_8
             ENDIF
            ELSE
             APOS = int( I - 1, 8 ) * LDA_SON8 + SHIFTCB_SON + 1_8
            ENDIF
         ENDIF
         CALL MPI_PACK( A_CBSON( APOS ), THIS_ROW_LENGTH,
     &        MPI_DOUBLE_PRECISION,
     &        BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &        POSITION, COMM, IERR )
        ENDDO
      IF (NBROWS_ALREADY_SENT == 0) THEN
        IF (COMPUTE_MAX) THEN
           CALL MPI_PACK(NFS4FATHER,1,
     &          MPI_INTEGER,
     &          BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &          POSITION, COMM, IERR )
           IF(NFS4FATHER .GT. 0) THEN
              BUF_MAX_ARRAY(1:NFS4FATHER) = ZERO
              IF(MAPROW(NROW) .GT. NASS_PERE) THEN
                 DO PS1=1,NROW
                    IF(MAPROW(PS1).GT.NASS_PERE) EXIT
                 ENDDO
                 IF (DESC_IN_LU) THEN
                   IF (COMPRESSCB) THEN
                    APOS = int(NELIM+PS1,8) * int(NELIM+PS1-1,8) /
     &                     2_8 + 1_8
                    NCA  = -44444
                    ASIZE  = int(NROW,8) * int(NROW+1,8)/2_8 -
     &                       int(NELIM+PS1,8) * int(NELIM+PS1-1,8)/2_8
                    LROW1  = PS1 + NELIM
                   ELSE
                    APOS = int(PS1+NELIM-1,8) * int(LROW,8) + 1_8
                    NCA = LROW
                    ASIZE = int(NCA,8) * int(NROW-PS1+1,8)
                    LROW1 = LROW
                   ENDIF
                 ELSE
                    IF (COMPRESSCB) THEN
                      IF (NPIV.NE.0) THEN
         WRITE(*,*) "Error in PARPIV/DMUMPS_BUF_SEND_CONTRIB_TYPE2"
                        CALL MUMPS_ABORT()
                      ENDIF
                      LROW1=LROW-NROW+PS1
                      ITMP8 = int(PS1 + LROW - NROW,8)
                      APOS = ITMP8 * (ITMP8 - 1_8) / 2_8 + 1_8 -
     &                       int(LROW-NROW,8)*int(LROW-NROW+1,8)/2_8
                      ASIZE = int(LROW,8)*int(LROW+1,8)/2_8 -
     &                       ITMP8*(ITMP8-1_8)/2_8
                      NCA   = -555555
                    ELSE
                      APOS = int(PS1-1,8) * LDA_SON8 + 1_8 + SHIFTCB_SON
                      NCA = int(LDA_SON8)
                      ASIZE = SIZFR - (SHIFTCB_SON -
     &                                 int(PS1-1,8) * LDA_SON8)
                      LROW1=-666666
                    ENDIF
                 ENDIF
                 IF ( NROW-PS1+1-KEEP253_LOC .NE. 0 ) THEN
                   CALL DMUMPS_COMPUTE_MAXPERCOL(
     &                A_CBSON(APOS),ASIZE,NCA,
     &                NROW-PS1+1-KEEP253_LOC,
     &                BUF_MAX_ARRAY,NFS4FATHER,COMPRESSCB,LROW1)
                 ENDIF
              ENDIF
              CALL MPI_PACK(BUF_MAX_ARRAY, NFS4FATHER,
     &             MPI_DOUBLE_PRECISION,
     &             BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &             POSITION, COMM, IERR )
           ENDIF
        ENDIF 
      ENDIF  
        CALL MPI_ISEND( BUF_CB%CONTENT( IPOS ), POSITION, MPI_PACKED,
     &                  PDEST, CONTRIB_TYPE2, COMM,
     &                  BUF_CB%CONTENT( IREQ ), IERR )
        IF ( SIZE_PACK.LT. POSITION ) THEN
          WRITE(*,*) ' contniv2: SIZE, POSITION =',SIZE_PACK, POSITION
          WRITE(*,*) ' NBROW, LROW = ', NBROW, LROW
          CALL MUMPS_ABORT()
        END IF
        IF ( SIZE_PACK .NE. POSITION )
     &  CALL BUF_ADJUST( BUF_CB, POSITION )
        NBROWS_ALREADY_SENT=NBROWS_ALREADY_SENT + NBROWS_PACKET
        IF (NBROWS_ALREADY_SENT .NE. NBROW ) THEN
           IERR = -1
        ENDIF
 100    CONTINUE
        RETURN
        END SUBROUTINE DMUMPS_BUF_SEND_CONTRIB_TYPE2
        SUBROUTINE DMUMPS_BUF_SEND_MAPLIG( 
     &                INODE, NFRONT, NASS1, NFS4FATHER,
     &                ISON, MYID, NSLAVES, SLAVES_PERE,
     &                TROW, NCBSON,
     &                COMM, IERR,
     &                DEST, NDEST, SLAVEF, 
     & 
     &                KEEP,KEEP8, STEP, N, 
     &                ISTEP_TO_INIV2, TAB_POS_IN_PERE
     &
     &                                  )
        IMPLICIT NONE
      INTEGER INODE, NFRONT, NASS1, NCBSON, NSLAVES, 
     &          NDEST
      INTEGER SLAVEF, MYID, ISON
      INTEGER TROW( NCBSON )
      INTEGER DEST( NDEST )
      INTEGER SLAVES_PERE( NSLAVES )
      INTEGER COMM, IERR
      INTEGER KEEP(500), N
      INTEGER(8) KEEP8(150)
      INTEGER STEP(N), 
     &        ISTEP_TO_INIV2(KEEP(71)), 
     &        TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
        INTEGER SIZE_AV, IDEST, NSEND, SIZE, NFS4FATHER
        INTEGER TROW_SIZE, POSITION, INDX, INIV2
        INTEGER IPOS, IREQ
        INTEGER IONE
        PARAMETER ( IONE=1 )
        INTEGER NASS_SON
        NASS_SON = -99998
        IERR = 0
        IF ( NDEST .eq. 1 ) THEN
          IF ( DEST(1).EQ.MYID )  GOTO 500
          SIZE = SIZEofINT * ( 7 + NSLAVES + NCBSON )
          IF ( NSLAVES.GT.0 ) THEN
             SIZE = SIZE + SIZEofINT * ( NSLAVES + 1 )
          ENDIF
          CALL BUF_LOOK( BUF_CB, IPOS, IREQ, SIZE, IERR, 
     &                 IONE, DEST
     &                 )
          IF (IERR .LT. 0 ) THEN
             RETURN
          ENDIF
          IF (SIZE.GT.SIZE_RBUF_BYTES ) THEN
             IERR = -3
             RETURN
          END IF
              POSITION = IPOS
              BUF_CB%CONTENT( POSITION ) = INODE
              POSITION = POSITION + 1
              BUF_CB%CONTENT( POSITION ) = ISON
              POSITION = POSITION + 1
              BUF_CB%CONTENT( POSITION ) = NSLAVES
              POSITION = POSITION + 1
              BUF_CB%CONTENT( POSITION ) = NFRONT
              POSITION = POSITION + 1
              BUF_CB%CONTENT( POSITION ) = NASS1
              POSITION = POSITION + 1
              BUF_CB%CONTENT( POSITION ) = NCBSON
              POSITION = POSITION + 1
              BUF_CB%CONTENT( POSITION ) = NFS4FATHER
              POSITION = POSITION + 1
              IF ( NSLAVES.GT.0 ) THEN
                INIV2 = ISTEP_TO_INIV2 ( STEP(INODE) )
                BUF_CB%CONTENT( POSITION: POSITION + NSLAVES )
     &          =  TAB_POS_IN_PERE(1:NSLAVES+1,INIV2)
                POSITION = POSITION + NSLAVES + 1
              ENDIF
              IF ( NSLAVES .NE. 0 ) THEN
                BUF_CB%CONTENT( POSITION: POSITION + NSLAVES - 1 )
     &          = SLAVES_PERE( 1: NSLAVES )
                POSITION = POSITION + NSLAVES
              END IF
              BUF_CB%CONTENT( POSITION:POSITION+NCBSON-1 ) =
     &        TROW( 1: NCBSON )
              POSITION = POSITION + NCBSON
              POSITION = POSITION - IPOS
              IF ( POSITION * SIZEofINT .NE. SIZE ) THEN
                WRITE(*,*) 'Error in DMUMPS_BUF_SEND_MAPLIG :',
     &                     ' wrong estimated size'
                CALL MUMPS_ABORT()
              END IF
              CALL MPI_ISEND( BUF_CB%CONTENT( IPOS ), SIZE,
     &                        MPI_PACKED,
     &                        DEST( NDEST ), MAPLIG, COMM,
     &                        BUF_CB%CONTENT( IREQ ),
     &                        IERR )
        ELSE
          NSEND = 0
          DO IDEST = 1, NDEST
            IF ( DEST( IDEST ) .ne. MYID ) NSEND = NSEND + 1
          END DO
          SIZE = SIZEofINT * 
     &         ( ( OVHSIZE + 7 + NSLAVES )* NSEND + NCBSON )
          IF ( NSLAVES.GT.0 ) THEN
           SIZE = SIZE + SIZEofINT * NSEND*( NSLAVES + 1 )
          ENDIF
          CALL DMUMPS_BUF_SIZE_AVAILABLE( BUF_CB, SIZE_AV )
          IF ( SIZE_AV .LT. SIZE ) THEN
            IERR = -1
            RETURN
          END IF
          DO IDEST= 1, NDEST
            CALL MUMPS_BLOC2_GET_SLAVE_INFO( 
     &                KEEP,KEEP8, ISON, STEP, N, SLAVEF,
     &                ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &                IDEST, NCBSON, 
     &                NDEST, 
     &                TROW_SIZE, INDX  )
            SIZE = SIZEofINT * ( NSLAVES + TROW_SIZE + 7 )
            IF ( NSLAVES.GT.0 ) THEN
             SIZE = SIZE + SIZEofINT * ( NSLAVES + 1 )
            ENDIF
            IF ( MYID .NE. DEST( IDEST ) ) THEN
              CALL BUF_LOOK( BUF_CB, IPOS, IREQ, SIZE, IERR,
     &                       IONE, DEST(IDEST)
     &                     )
              IF ( IERR .LT. 0 )  THEN
                 WRITE(*,*) 'Problem in BUF_LOOK: IERR<0'
                 CALL MUMPS_ABORT()
              END IF
              IF (SIZE.GT.SIZE_RBUF_BYTES) THEN
                 IERR = -3
                 RETURN
              ENDIF
              POSITION = IPOS
              BUF_CB%CONTENT( POSITION ) = INODE
              POSITION = POSITION + 1
              BUF_CB%CONTENT( POSITION ) = ISON
              POSITION = POSITION + 1
              BUF_CB%CONTENT( POSITION ) = NSLAVES
              POSITION = POSITION + 1
              BUF_CB%CONTENT( POSITION ) = NFRONT
              POSITION = POSITION + 1
              BUF_CB%CONTENT( POSITION ) = NASS1
              POSITION = POSITION + 1
              BUF_CB%CONTENT( POSITION ) = TROW_SIZE
              POSITION = POSITION + 1
              BUF_CB%CONTENT( POSITION ) = NFS4FATHER
              POSITION = POSITION + 1
              IF ( NSLAVES.GT.0 ) THEN
                INIV2 = ISTEP_TO_INIV2 ( STEP(INODE) )
                BUF_CB%CONTENT( POSITION: POSITION + NSLAVES )
     &          =  TAB_POS_IN_PERE(1:NSLAVES+1,INIV2)
                POSITION = POSITION + NSLAVES + 1
              ENDIF
              IF ( NSLAVES .NE. 0 ) THEN
                BUF_CB%CONTENT( POSITION: POSITION + NSLAVES - 1 )
     &          = SLAVES_PERE( 1: NSLAVES )
                POSITION = POSITION + NSLAVES
              END IF
              BUF_CB%CONTENT( POSITION:POSITION+TROW_SIZE-1 ) =
     &        TROW( INDX: INDX + TROW_SIZE - 1 )
              POSITION = POSITION + TROW_SIZE
              POSITION = POSITION - IPOS
              IF ( POSITION * SIZEofINT .NE. SIZE ) THEN
               WRITE(*,*) ' ERROR 1 in TRY_SEND_MAPLIG:',
     &          'Wrong estimated size'
               CALL MUMPS_ABORT()
              END IF
              CALL MPI_ISEND( BUF_CB%CONTENT( IPOS ), SIZE,
     &                        MPI_PACKED,
     &                        DEST( IDEST ), MAPLIG, COMM,
     &                        BUF_CB%CONTENT( IREQ ),
     &                        IERR )
            END IF
          END DO
        END IF
 500    CONTINUE
        RETURN
        END SUBROUTINE DMUMPS_BUF_SEND_MAPLIG
        SUBROUTINE DMUMPS_BUF_SEND_BLOCFACTO( INODE, NFRONT,
     &             NCOL, NPIV, FPERE, LASTBL, IPIV, VAL,
     &             PDEST, NDEST, KEEP50, K34, NB_BLOC_FAC,
     &             NSLAVES_TOT,
     &             WIDTH, COMM,
     &
     &             IERR )
      IMPLICIT NONE
        INTEGER, intent(in) :: INODE, NCOL, NPIV, 
     &                         FPERE, NFRONT, NDEST
        INTEGER, intent(in) :: IPIV( NPIV )
        DOUBLE PRECISION, intent(in) :: VAL( NFRONT, * )
        INTEGER, intent(in) :: PDEST( NDEST ) 
        INTEGER, intent(in) :: KEEP50, K34, NB_BLOC_FAC,
     &                         NSLAVES_TOT, COMM, WIDTH
        LOGICAL, intent(in) :: LASTBL
        INTEGER, intent(inout) :: IERR
        INCLUDE 'mpif.h'
        INCLUDE 'mumps_tags.h'
        INTEGER POSITION, IREQ, IPOS, SIZE1, SIZE2, SIZE3, SIZET,
     &          IDEST, IPOSMSG, I
        INTEGER NPIVSENT
        INTEGER SSS
        INTEGER  :: NBMSGS
        INTEGER, ALLOCATABLE, DIMENSION(:) ::  RELAY_INFO
        INTEGER :: LRELAY_INFO, DEST_BLOCFACTO, TAG_BLOCFACTO
        IERR = 0
        LRELAY_INFO = 0
        NBMSGS = NDEST
        IF ( LASTBL ) THEN
          IF ( KEEP50 .eq. 0 ) THEN
            CALL MPI_PACK_SIZE( 4 + NPIV + ( NBMSGS - 1 ) * OVHSIZE +
     &                          1+LRELAY_INFO,
     &                          MPI_INTEGER, COMM, SIZE1, IERR )
          ELSE
            CALL MPI_PACK_SIZE( 6 + NPIV + ( NBMSGS - 1 ) * OVHSIZE + 
     &                          1+LRELAY_INFO,
     &                          MPI_INTEGER, COMM, SIZE1, IERR )
          END IF
        ELSE
          IF ( KEEP50 .eq. 0 ) THEN
          CALL MPI_PACK_SIZE( 3 + NPIV + ( NBMSGS - 1 ) * OVHSIZE + 
     &                        1+LRELAY_INFO,
     &                        MPI_INTEGER, COMM, SIZE1, IERR )
          ELSE
            CALL MPI_PACK_SIZE( 4 + NPIV + ( NBMSGS - 1 ) * OVHSIZE + 
     &                          1+LRELAY_INFO,
     &                          MPI_INTEGER, COMM, SIZE1, IERR )
          END IF
        END IF
        SIZE2 = 0
        IF (NPIV.GT.0) THEN
            CALL MPI_PACK_SIZE( NPIV*NCOL, MPI_DOUBLE_PRECISION,
     &                      COMM, SIZE3, IERR )
            SIZE2 = SIZE2+SIZE3
        ENDIF
        SIZET = SIZE1 + SIZE2 
        IF (LRELAY_INFO.GT.0) THEN
         CALL BUF_LOOK( BUF_CB, IPOS, IREQ, SIZET, IERR, 
     &                 NBMSGS , RELAY_INFO(2)
     &               )
        ELSE
         CALL BUF_LOOK( BUF_CB, IPOS, IREQ, SIZET, IERR, 
     &                 NBMSGS , PDEST
     &               )
        ENDIF
        IF ( IERR .LT. 0 ) THEN
           RETURN
        ENDIF
        IF (SIZET.GT.SIZE_RBUF_BYTES) THEN
          SSS = 0 
          IF ( LASTBL ) THEN
           IF ( KEEP50 .eq. 0 ) THEN
            CALL MPI_PACK_SIZE( 4 + NPIV + 1+LRELAY_INFO,
     &                        MPI_INTEGER, COMM, SSS, IERR )
           ELSE
            CALL MPI_PACK_SIZE( 6 + NPIV + 1+LRELAY_INFO, 
     &                           MPI_INTEGER, COMM, SSS, IERR )
           END IF
          ELSE
           IF ( KEEP50 .eq. 0 ) THEN
            CALL MPI_PACK_SIZE( 3 + NPIV + 1+LRELAY_INFO,
     &                        MPI_INTEGER, COMM, SSS, IERR )
           ELSE
            CALL MPI_PACK_SIZE( 4 + NPIV + 1+LRELAY_INFO,
     &                        MPI_INTEGER, COMM, SSS, IERR )
           END IF
          END IF
          SSS = SSS + SIZE2  
          IF (SSS.GT.SIZE_RBUF_BYTES) THEN
           IERR = -2 
           RETURN
          ENDIF
        ENDIF
        BUF_CB%ILASTMSG = BUF_CB%ILASTMSG + ( NBMSGS - 1 ) * OVHSIZE
        IPOS = IPOS - OVHSIZE
        DO IDEST = 1, NBMSGS - 1
          BUF_CB%CONTENT( IPOS + ( IDEST - 1 ) * OVHSIZE ) =
     &    IPOS + IDEST * OVHSIZE
        END DO
        BUF_CB%CONTENT( IPOS + ( NBMSGS - 1 ) * OVHSIZE ) = 0
        IPOSMSG = IPOS + OVHSIZE * NBMSGS
        POSITION = 0
        CALL MPI_PACK( INODE, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOSMSG ), SIZET,
     &                        POSITION, COMM, IERR )
        NPIVSENT = NPIV
        IF (LASTBL) NPIVSENT = -NPIV
        CALL MPI_PACK( NPIVSENT, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOSMSG ), SIZET,
     &                        POSITION, COMM, IERR )
        IF ( LASTBL .or. KEEP50.ne.0 ) THEN
          CALL MPI_PACK( FPERE, 1, MPI_INTEGER,
     &                   BUF_CB%CONTENT( IPOSMSG ), SIZET,
     &                   POSITION, COMM, IERR )
        END IF
        IF ( LASTBL .AND. KEEP50 .NE. 0 ) THEN
            CALL MPI_PACK( NSLAVES_TOT, 1, MPI_INTEGER,
     &                   BUF_CB%CONTENT( IPOSMSG ), SIZET,
     &                   POSITION, COMM, IERR )
            CALL MPI_PACK( NB_BLOC_FAC, 1, MPI_INTEGER,
     &                   BUF_CB%CONTENT( IPOSMSG ), SIZET,
     &                   POSITION, COMM, IERR )
        END IF
        CALL MPI_PACK( NCOL, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOSMSG ), SIZET,
     &                        POSITION, COMM, IERR )
        IF ( NPIV.GT.0) THEN
          CALL MPI_PACK( IPIV, NPIV, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOSMSG ), SIZET,
     &                        POSITION, COMM, IERR )
            DO I = 1, NPIV
              CALL MPI_PACK( VAL(1,I), NCOL,
     &                        MPI_DOUBLE_PRECISION,
     &                        BUF_CB%CONTENT( IPOSMSG ), SIZET,
     &                        POSITION, COMM, IERR )
            END DO
        ENDIF
        CALL MPI_PACK( LRELAY_INFO, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOSMSG ), SIZET,
     &                        POSITION, COMM, IERR )
        IF ( LRELAY_INFO.GT.0) 
     &    CALL MPI_PACK( RELAY_INFO, LRELAY_INFO, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOSMSG ), SIZET,
     &                        POSITION, COMM, IERR )
        DO IDEST = 1, NBMSGS
          IF (LRELAY_INFO .GT. 0) THEN
            DEST_BLOCFACTO = RELAY_INFO(IDEST+1)
          ELSE
            DEST_BLOCFACTO = PDEST(IDEST)
          ENDIF
          IF ( KEEP50 .EQ. 0) THEN
            TAG_BLOCFACTO = BLOC_FACTO
            CALL MPI_ISEND( BUF_CB%CONTENT( IPOSMSG ), POSITION, 
     &                MPI_PACKED,
     &                DEST_BLOCFACTO, TAG_BLOCFACTO, COMM,
     &                BUF_CB%CONTENT( IREQ + ( IDEST-1 ) * OVHSIZE ),
     &                IERR )
          ELSE
            CALL MPI_ISEND( BUF_CB%CONTENT( IPOSMSG ), POSITION, 
     &                MPI_PACKED,
     &                DEST_BLOCFACTO, BLOC_FACTO_SYM, COMM,
     &                BUF_CB%CONTENT( IREQ + ( IDEST-1 ) * OVHSIZE ),
     &                IERR )
          END IF
        END DO
        SIZET = SIZET - ( NBMSGS - 1 ) * OVHSIZE * SIZEofINT
        IF ( SIZET .LT. POSITION ) THEN
          WRITE(*,*) ' Error sending blocfacto : size < position'
          WRITE(*,*) ' Size,position=',SIZET,POSITION
          CALL MUMPS_ABORT()
        END IF
        IF ( SIZET .NE. POSITION ) CALL BUF_ADJUST( BUF_CB, POSITION )
        RETURN
        END SUBROUTINE DMUMPS_BUF_SEND_BLOCFACTO
        SUBROUTINE DMUMPS_BUF_SEND_BLFAC_SLAVE( INODE,
     &             NPIV, FPERE, IPOSK, JPOSK, UIP21K, NCOLU,
     &             NDEST, PDEST, COMM, 
     &             IERR )
        IMPLICIT NONE
        INTEGER INODE, NCOLU, IPOSK, JPOSK, NPIV, NDEST, FPERE
        DOUBLE PRECISION UIP21K( NPIV, * )
        INTEGER PDEST( NDEST ) 
        INTEGER   COMM, IERR
        INCLUDE 'mpif.h'
        INCLUDE 'mumps_tags.h'
        INTEGER POSITION, IREQ, IPOS, SIZE1, SIZE2, SIZET,
     &          IDEST, IPOSMSG, SSS, SSLR
        IERR = 0
        CALL MPI_PACK_SIZE( 6 + ( NDEST - 1 ) * OVHSIZE,
     &                      MPI_INTEGER, COMM, SIZE1, IERR )
        SIZE2  = 0
        CALL MPI_PACK_SIZE( abs(NPIV)*NCOLU, MPI_DOUBLE_PRECISION,
     &                      COMM, SSLR, IERR )
         SIZE2=SIZE2+SSLR
        SIZET = SIZE1 + SIZE2
        IF (SIZET.GT.SIZE_RBUF_BYTES) THEN
         CALL MPI_PACK_SIZE( 6 ,
     &                      MPI_INTEGER, COMM, SSS, IERR )
         SSS = SSS+SIZE2
         IF (SSS.GT.SIZE_RBUF_BYTES) THEN
           IERR = -2
           RETURN
         ENDIF
        END IF
        CALL BUF_LOOK( BUF_CB, IPOS, IREQ, SIZET, IERR, 
     &                 NDEST, PDEST
     &               )
        IF ( IERR .LT. 0 ) THEN
           RETURN
        ENDIF
        BUF_CB%ILASTMSG = BUF_CB%ILASTMSG + ( NDEST - 1 ) * OVHSIZE
        IPOS = IPOS - OVHSIZE
        DO IDEST = 1, NDEST - 1
          BUF_CB%CONTENT( IPOS + ( IDEST - 1 ) * OVHSIZE ) =
     &    IPOS + IDEST * OVHSIZE
        END DO
        BUF_CB%CONTENT( IPOS + ( NDEST - 1 ) * OVHSIZE ) = 0
        IPOSMSG = IPOS + OVHSIZE * NDEST
        POSITION = 0
        CALL MPI_PACK( INODE, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOSMSG ), SIZET,
     &                        POSITION, COMM, IERR )
        CALL MPI_PACK( IPOSK, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOSMSG ), SIZET,
     &                        POSITION, COMM, IERR )
        CALL MPI_PACK( JPOSK, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOSMSG ), SIZET,
     &                        POSITION, COMM, IERR )
        CALL MPI_PACK( NPIV, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOSMSG ), SIZET,
     &                        POSITION, COMM, IERR )
        CALL MPI_PACK( FPERE, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOSMSG ), SIZET,
     &                        POSITION, COMM, IERR )
        CALL MPI_PACK( NCOLU, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOSMSG ), SIZET,
     &                        POSITION, COMM, IERR )
        CALL MPI_PACK( UIP21K, abs(NPIV) * NCOLU,
     &                        MPI_DOUBLE_PRECISION,
     &                        BUF_CB%CONTENT( IPOSMSG ), SIZET,
     &                        POSITION, COMM, IERR )
        DO IDEST = 1, NDEST
        CALL MPI_ISEND( BUF_CB%CONTENT( IPOSMSG ), POSITION, MPI_PACKED,
     &                  PDEST(IDEST), BLOC_FACTO_SYM_SLAVE, COMM,
     &                  BUF_CB%CONTENT( IREQ + ( IDEST-1 ) * OVHSIZE ),
     &                  IERR )
        END DO
        SIZET = SIZET - ( NDEST - 1 ) * OVHSIZE * SIZEofINT
        IF ( SIZET .LT. POSITION ) THEN
          WRITE(*,*) ' Error sending blfac slave : size < position'
          WRITE(*,*) ' Size,position=',SIZET,POSITION
          CALL MUMPS_ABORT()
        END IF
        IF ( SIZET .NE. POSITION ) CALL BUF_ADJUST( BUF_CB, POSITION )
        RETURN
        END SUBROUTINE DMUMPS_BUF_SEND_BLFAC_SLAVE
        SUBROUTINE DMUMPS_BUF_SEND_CONTRIB_TYPE3( N, ISON,
     &             NBCOL_SON, NBROW_SON, INDCOL_SON, INDROW_SON,
     &             LD_SON, VAL_SON, TAG, SUBSET_ROW, SUBSET_COL,
     &             NSUBSET_ROW, NSUBSET_COL,
     &             NSUPROW, NSUPCOL,
     &             NPROW, NPCOL, MBLOCK, RG2L_ROW, RG2L_COL,
     &             NBLOCK, PDEST, COMM, IERR , 
     &             TAB, TABSIZE, TRANSP, SIZE_PACK,
     &             N_ALREADY_SENT, KEEP, BBPCBP ) 
        IMPLICIT NONE
        INTEGER N, ISON, NBCOL_SON, NBROW_SON, NSUBSET_ROW, NSUBSET_COL
        INTEGER NPROW, NPCOL, MBLOCK, NBLOCK, LD_SON
        INTEGER BBPCBP
        INTEGER PDEST, TAG, COMM, IERR
        INTEGER INDCOL_SON( NBCOL_SON ), INDROW_SON( NBROW_SON )
        INTEGER SUBSET_ROW( NSUBSET_ROW ), SUBSET_COL( NSUBSET_COL )
        INTEGER, DIMENSION(:) :: RG2L_ROW  
        INTEGER, DIMENSION(:) :: RG2L_COL  
        INTEGER NSUPROW, NSUPCOL
        INTEGER(8), INTENT(IN) :: TABSIZE
        INTEGER SIZE_PACK
        INTEGER KEEP(500)
        DOUBLE PRECISION VAL_SON( LD_SON, * ), TAB(*)
        LOGICAL TRANSP
        INTEGER N_ALREADY_SENT
        INCLUDE 'mpif.h'
        INTEGER SIZE1, SIZE2, SIZE_AV, POSITION
        INTEGER SIZE_CBP, SIZE_TMP
        INTEGER IREQ, IPOS, ITAB
        INTEGER ISUB, JSUB, I, J 
        INTEGER ILOC_ROOT, JLOC_ROOT
        INTEGER IPOS_ROOT, JPOS_ROOT
        INTEGER IONE
        LOGICAL RECV_BUF_SMALLER_THAN_SEND
        INTEGER PDEST2(1)
        PARAMETER ( IONE=1 )
        INTEGER N_PACKET
        INTEGER NSUBSET_ROW_EFF, NSUBSET_COL_EFF, NSUPCOL_EFF
        PDEST2(1) = PDEST
        IERR = 0
        IF ( NSUBSET_ROW * NSUBSET_COL .NE. 0 ) THEN
          CALL DMUMPS_BUF_SIZE_AVAILABLE( BUF_CB, SIZE_AV )
          IF (SIZE_AV .LT. SIZE_RBUF_BYTES) THEN
            RECV_BUF_SMALLER_THAN_SEND = .FALSE.
          ELSE
            RECV_BUF_SMALLER_THAN_SEND = .TRUE.
            SIZE_AV = SIZE_RBUF_BYTES
          ENDIF
          SIZE_AV = min(SIZE_AV, SIZE_RBUF_BYTES)
          CALL MPI_PACK_SIZE(8 + NSUBSET_COL,
     &                      MPI_INTEGER, COMM, SIZE1, IERR )
          SIZE_CBP = 0
          IF (N_ALREADY_SENT .EQ. 0 .AND.
     &        min(NSUPROW,NSUPCOL) .GT.0) THEN
            CALL MPI_PACK_SIZE(NSUPROW, MPI_INTEGER, COMM,
     &           SIZE_CBP, IERR)
            CALL MPI_PACK_SIZE(NSUPCOL, MPI_INTEGER, COMM,
     &           SIZE_TMP, IERR)
            SIZE_CBP = SIZE_CBP + SIZE_TMP
            CALL MPI_PACK_SIZE(NSUPROW*NSUPCOL,
     &           MPI_DOUBLE_PRECISION, COMM,
     &           SIZE_TMP, IERR)
            SIZE_CBP = SIZE_CBP + SIZE_TMP
            SIZE1 = SIZE1 + SIZE_CBP
          ENDIF
          IF (BBPCBP.EQ.1) THEN
            NSUBSET_COL_EFF = NSUBSET_COL - NSUPCOL
            NSUPCOL_EFF = 0
          ELSE
            NSUBSET_COL_EFF = NSUBSET_COL
            NSUPCOL_EFF = NSUPCOL
          ENDIF
          NSUBSET_ROW_EFF = NSUBSET_ROW - NSUPROW
          N_PACKET =
     &    (SIZE_AV - SIZE1) / (SIZEofINT + NSUBSET_COL_EFF * SIZEofREAL)
 10       CONTINUE
          N_PACKET = min( N_PACKET,
     &                    NSUBSET_ROW_EFF-N_ALREADY_SENT )
          IF (N_PACKET .LE. 0 .AND.
     &        NSUBSET_ROW_EFF-N_ALREADY_SENT.GT.0) THEN
            IF (RECV_BUF_SMALLER_THAN_SEND) THEN
              IERR=-3
              GOTO 100
            ELSE
              IERR = -1
              GOTO 100
            ENDIF
          ENDIF
          CALL MPI_PACK_SIZE( 8 + NSUBSET_COL_EFF + N_PACKET,
     &                      MPI_INTEGER, COMM, SIZE1, IERR )
          SIZE1 = SIZE1 + SIZE_CBP
          CALL MPI_PACK_SIZE( N_PACKET * NSUBSET_COL_EFF,
     &                      MPI_DOUBLE_PRECISION,
     &                      COMM, SIZE2, IERR )
          SIZE_PACK = SIZE1 + SIZE2
          IF (SIZE_PACK .GT. SIZE_AV) THEN
            N_PACKET = N_PACKET - 1
            IF ( N_PACKET > 0 ) THEN
              GOTO 10
            ELSE
              IF (RECV_BUF_SMALLER_THAN_SEND) THEN
                IERR = -3
                GOTO 100
              ELSE
                IERR = -1
                GOTO 100
              ENDIF
            ENDIF
          ENDIF
#if ! defined(DBG_SMB3)
          IF (N_PACKET + N_ALREADY_SENT .NE. NSUBSET_ROW - NSUPROW
     &      .AND.
     &      SIZE_PACK .LT. SIZE_RBUF_BYTES / 4
     &      .AND. .NOT. RECV_BUF_SMALLER_THAN_SEND)
     &      THEN
            IERR = -1
            GOTO 100
          ENDIF
#endif
        ELSE 
          N_PACKET = 0
          CALL MPI_PACK_SIZE(8,MPI_INTEGER, COMM, SIZE_PACK, IERR )
        END IF
        CALL BUF_LOOK( BUF_CB, IPOS, IREQ, SIZE_PACK, IERR, 
     &                 IONE, PDEST2
     &               )
        IF ( IERR .LT. 0 ) GOTO 100
        IF ( SIZE_PACK.GT.SIZE_RBUF_BYTES ) THEN
             IERR = -3
             GOTO 100
        ENDIF
        POSITION = 0
        CALL MPI_PACK( ISON, 1, MPI_INTEGER,
     &                 BUF_CB%CONTENT( IPOS ),
     &                 SIZE_PACK, POSITION, COMM, IERR )
        CALL MPI_PACK( NSUBSET_ROW, 1, MPI_INTEGER,
     &                 BUF_CB%CONTENT( IPOS ),
     &                 SIZE_PACK, POSITION, COMM, IERR )
        CALL MPI_PACK( NSUPROW, 1, MPI_INTEGER,
     &                 BUF_CB%CONTENT( IPOS ),
     &                 SIZE_PACK, POSITION, COMM, IERR )
        CALL MPI_PACK( NSUBSET_COL, 1, MPI_INTEGER,
     &                 BUF_CB%CONTENT( IPOS ),
     &                 SIZE_PACK, POSITION, COMM, IERR )
        CALL MPI_PACK( NSUPCOL, 1, MPI_INTEGER,
     &                 BUF_CB%CONTENT( IPOS ),
     &                 SIZE_PACK, POSITION, COMM, IERR )
        CALL MPI_PACK( N_ALREADY_SENT, 1, MPI_INTEGER,
     &                 BUF_CB%CONTENT( IPOS ),
     &                 SIZE_PACK, POSITION, COMM, IERR )
        CALL MPI_PACK( N_PACKET, 1, MPI_INTEGER,
     &                 BUF_CB%CONTENT( IPOS ),
     &                 SIZE_PACK, POSITION, COMM, IERR )
        CALL MPI_PACK( BBPCBP, 1, MPI_INTEGER,
     &                 BUF_CB%CONTENT( IPOS ),
     &                 SIZE_PACK, POSITION, COMM, IERR )
        IF ( NSUBSET_ROW * NSUBSET_COL .NE. 0 ) THEN
          IF (N_ALREADY_SENT .EQ. 0 .AND.
     &          min(NSUPROW, NSUPCOL) .GT. 0) THEN
            DO ISUB = NSUBSET_ROW-NSUPROW+1, NSUBSET_ROW
              I =  SUBSET_ROW( ISUB )
              IPOS_ROOT = RG2L_ROW(INDCOL_SON( I ))
              ILOC_ROOT = MBLOCK
     &                 * ( ( IPOS_ROOT - 1 ) / ( MBLOCK * NPROW ) )
     &                 + mod( IPOS_ROOT - 1, MBLOCK ) + 1
              CALL MPI_PACK( ILOC_ROOT, 1, MPI_INTEGER,
     &                      BUF_CB%CONTENT( IPOS ),
     &                      SIZE_PACK, POSITION, COMM, IERR )
            ENDDO
            DO ISUB = NSUBSET_COL-NSUPCOL+1, NSUBSET_COL
               J = SUBSET_COL( ISUB )
               JPOS_ROOT = INDROW_SON( J ) - N
               JLOC_ROOT = NBLOCK
     &                  * ( ( JPOS_ROOT - 1 ) / ( NBLOCK * NPCOL ) )
     &                  + mod( JPOS_ROOT - 1, NBLOCK ) + 1
              CALL MPI_PACK( JLOC_ROOT, 1, MPI_INTEGER,
     &                       BUF_CB%CONTENT( IPOS ),
     &                       SIZE_PACK, POSITION, COMM, IERR )
            ENDDO
            IF ( TABSIZE.GE.int(NSUPROW,8)*int(NSUPCOL,8) ) THEN
              ITAB = 1
              DO JSUB = NSUBSET_ROW - NSUPROW+1, NSUBSET_ROW
                J = SUBSET_ROW(JSUB)
                DO ISUB = NSUBSET_COL - NSUPCOL+1, NSUBSET_COL
                  I = SUBSET_COL(ISUB)
                  TAB(ITAB) = VAL_SON(J, I)
                  ITAB = ITAB + 1
                ENDDO
              ENDDO
              CALL MPI_PACK(TAB(1), NSUPROW*NSUPCOL,
     &         MPI_DOUBLE_PRECISION, 
     &         BUF_CB%CONTENT( IPOS ),
     &         SIZE_PACK, POSITION, COMM, IERR )
            ELSE
              DO JSUB = NSUBSET_ROW - NSUPROW+1, NSUBSET_ROW
                J = SUBSET_ROW(JSUB)
                DO ISUB = NSUBSET_COL - NSUPCOL+1, NSUBSET_COL
                  I = SUBSET_COL(ISUB)
                  CALL MPI_PACK(VAL_SON(J,I), 1,
     &            MPI_DOUBLE_PRECISION, 
     &            BUF_CB%CONTENT( IPOS ),
     &            SIZE_PACK, POSITION, COMM, IERR )
                ENDDO
              ENDDO
            ENDIF
          ENDIF
          IF ( .NOT. TRANSP ) THEN
            DO ISUB = N_ALREADY_SENT+1, N_ALREADY_SENT+N_PACKET
              I         = SUBSET_ROW( ISUB )
              IPOS_ROOT = RG2L_ROW( INDROW_SON( I ) )
              ILOC_ROOT = MBLOCK
     &                 * ( ( IPOS_ROOT - 1 ) / ( MBLOCK * NPROW ) )
     &                 + mod( IPOS_ROOT - 1, MBLOCK ) + 1
              CALL MPI_PACK( ILOC_ROOT, 1, MPI_INTEGER,
     &                      BUF_CB%CONTENT( IPOS ),
     &                      SIZE_PACK, POSITION, COMM, IERR )
            END DO
            DO JSUB = 1, NSUBSET_COL_EFF - NSUPCOL_EFF
              J         = SUBSET_COL( JSUB )
              JPOS_ROOT = RG2L_COL( INDCOL_SON( J ) )
              JLOC_ROOT = NBLOCK
     &                  * ( ( JPOS_ROOT - 1 ) / ( NBLOCK * NPCOL ) )
     &                  + mod( JPOS_ROOT - 1, NBLOCK ) + 1
              CALL MPI_PACK( JLOC_ROOT, 1, MPI_INTEGER,
     &                       BUF_CB%CONTENT( IPOS ),
     &                       SIZE_PACK, POSITION, COMM, IERR )
            END DO
            DO JSUB = NSUBSET_COL_EFF-NSUPCOL_EFF+1, NSUBSET_COL_EFF
               J = SUBSET_COL( JSUB )
               JPOS_ROOT = INDCOL_SON( J ) - N
               JLOC_ROOT = NBLOCK
     &                  * ( ( JPOS_ROOT - 1 ) / ( NBLOCK * NPCOL ) )
     &                  + mod( JPOS_ROOT - 1, NBLOCK ) + 1
              CALL MPI_PACK( JLOC_ROOT, 1, MPI_INTEGER,
     &                       BUF_CB%CONTENT( IPOS ),
     &                       SIZE_PACK, POSITION, COMM, IERR )
            ENDDO
          ELSE
            DO JSUB = N_ALREADY_SENT+1, N_ALREADY_SENT+N_PACKET
              J         = SUBSET_ROW( JSUB )
              IPOS_ROOT = RG2L_ROW( INDCOL_SON( J ) )
              ILOC_ROOT = MBLOCK
     &                 * ( ( IPOS_ROOT - 1 ) / ( MBLOCK * NPROW ) )
     &                 + mod( IPOS_ROOT - 1, MBLOCK ) + 1
              CALL MPI_PACK( ILOC_ROOT, 1, MPI_INTEGER,
     &                       BUF_CB%CONTENT( IPOS ),
     &                       SIZE_PACK, POSITION, COMM, IERR )
            END DO
            DO ISUB = 1, NSUBSET_COL_EFF - NSUPCOL_EFF
              I         = SUBSET_COL( ISUB )
              JPOS_ROOT = RG2L_COL( INDROW_SON( I ) )
              JLOC_ROOT = NBLOCK
     &                  * ( ( JPOS_ROOT - 1 ) / ( NBLOCK * NPCOL ) )
     &                  + mod( JPOS_ROOT - 1, NBLOCK ) + 1
              CALL MPI_PACK( JLOC_ROOT, 1, MPI_INTEGER,
     &                      BUF_CB%CONTENT( IPOS ),
     &                      SIZE_PACK, POSITION, COMM, IERR )
            END DO
            DO ISUB = NSUBSET_COL_EFF - NSUPCOL_EFF + 1, NSUBSET_COL_EFF
              I         = SUBSET_COL( ISUB )
              JPOS_ROOT = INDROW_SON(I) - N
              JLOC_ROOT = NBLOCK
     &                  * ( ( JPOS_ROOT - 1 ) / ( NBLOCK * NPCOL ) )
     &                  + mod( JPOS_ROOT - 1, NBLOCK ) + 1
              CALL MPI_PACK( JLOC_ROOT, 1, MPI_INTEGER,
     &                      BUF_CB%CONTENT( IPOS ),
     &                      SIZE_PACK, POSITION, COMM, IERR )
            ENDDO
          END IF
          IF ( TABSIZE.GE.int(N_PACKET,8)*int(NSUBSET_COL_EFF,8) ) THEN
            IF ( .NOT. TRANSP ) THEN
              ITAB = 1
              DO ISUB = N_ALREADY_SENT+1,
     &                  N_ALREADY_SENT+N_PACKET
                I         = SUBSET_ROW( ISUB )
                DO JSUB = 1, NSUBSET_COL_EFF
                  J              = SUBSET_COL( JSUB )
                  TAB( ITAB )    = VAL_SON(J,I)
                  ITAB           = ITAB + 1
                END DO
              END DO
              CALL MPI_PACK(TAB(1), NSUBSET_COL_EFF*N_PACKET,
     &         MPI_DOUBLE_PRECISION, 
     &         BUF_CB%CONTENT( IPOS ),
     &         SIZE_PACK, POSITION, COMM, IERR )
            ELSE
              ITAB = 1
              DO JSUB = N_ALREADY_SENT+1, N_ALREADY_SENT+N_PACKET
                J = SUBSET_ROW( JSUB )
                DO ISUB = 1, NSUBSET_COL_EFF
                  I         = SUBSET_COL( ISUB )
                  TAB( ITAB ) = VAL_SON( J, I )
                  ITAB = ITAB + 1
                END DO
              END DO
              CALL MPI_PACK(TAB(1), NSUBSET_COL_EFF*N_PACKET,
     &         MPI_DOUBLE_PRECISION, 
     &         BUF_CB%CONTENT( IPOS ),
     &         SIZE_PACK, POSITION, COMM, IERR )
            END IF
          ELSE
            IF ( .NOT. TRANSP ) THEN
              DO ISUB = N_ALREADY_SENT+1, N_ALREADY_SENT+N_PACKET
                I         = SUBSET_ROW( ISUB )
                DO JSUB = 1, NSUBSET_COL_EFF
                  J         = SUBSET_COL( JSUB )
                  CALL MPI_PACK( VAL_SON( J, I ), 1,
     &            MPI_DOUBLE_PRECISION,
     &            BUF_CB%CONTENT( IPOS ),
     &            SIZE_PACK, POSITION, COMM, IERR )
                END DO
              END DO
            ELSE
              DO JSUB = N_ALREADY_SENT+1, N_ALREADY_SENT+N_PACKET
                J = SUBSET_ROW( JSUB )
                DO ISUB = 1, NSUBSET_COL_EFF
                  I         = SUBSET_COL( ISUB )
                  CALL MPI_PACK( VAL_SON( J, I ), 1,
     &            MPI_DOUBLE_PRECISION,
     &            BUF_CB%CONTENT( IPOS ),
     &            SIZE_PACK, POSITION, COMM, IERR )
                END DO
              END DO
            END IF
          ENDIF
        END IF
        CALL MPI_ISEND( BUF_CB%CONTENT( IPOS ), POSITION, MPI_PACKED,
     &                PDEST, TAG, COMM, BUF_CB%CONTENT( IREQ ), IERR )
        IF ( SIZE_PACK .LT. POSITION ) THEN
          WRITE(*,*) ' Error sending contribution to root:Size<positn'
          WRITE(*,*) ' Size,position=',SIZE_PACK,POSITION
          CALL MUMPS_ABORT()
        END IF
        IF ( SIZE_PACK .NE. POSITION )
     &  CALL BUF_ADJUST( BUF_CB, POSITION )
        N_ALREADY_SENT = N_ALREADY_SENT + N_PACKET
        IF (NSUBSET_ROW * NSUBSET_COL .NE. 0) THEN
          IF ( N_ALREADY_SENT.NE.NSUBSET_ROW_EFF ) IERR = -1
        ENDIF
  100   CONTINUE
        RETURN
        END SUBROUTINE DMUMPS_BUF_SEND_CONTRIB_TYPE3
        SUBROUTINE DMUMPS_BUF_SEND_RTNELIND( ISON, NELIM,
     &             NELIM_ROW, NELIM_COL, NSLAVES, SLAVES,
     &             DEST, COMM, IERR )
        INTEGER ISON, NELIM
        INTEGER NSLAVES, DEST, COMM, IERR
        INTEGER NELIM_ROW( NELIM ), NELIM_COL( NELIM )
        INTEGER SLAVES( NSLAVES )
        INCLUDE 'mpif.h'
        INCLUDE 'mumps_tags.h'
        INTEGER SIZE, POSITION, IPOS, IREQ
        INTEGER IONE
        INTEGER DEST2(1)
        PARAMETER ( IONE=1 )
        DEST2(1) = DEST
        IERR = 0
        SIZE = ( 3 + NSLAVES + 2 * NELIM ) * SIZEofINT
        CALL BUF_LOOK( BUF_CB, IPOS, IREQ, SIZE, IERR, 
     &                 IONE, DEST2
     &               )
        IF ( IERR .LT. 0 ) THEN
           RETURN
        ENDIF
          IF (SIZE.GT.SIZE_RBUF_BYTES) THEN
             IERR = -3
             RETURN
          ENDIF
        POSITION = IPOS
        BUF_CB%CONTENT( POSITION ) = ISON
        POSITION = POSITION + 1
        BUF_CB%CONTENT( POSITION ) = NELIM
        POSITION = POSITION + 1
        BUF_CB%CONTENT( POSITION ) = NSLAVES
        POSITION = POSITION + 1
        BUF_CB%CONTENT( POSITION: POSITION + NELIM - 1 ) = NELIM_ROW
        POSITION = POSITION + NELIM
        BUF_CB%CONTENT( POSITION: POSITION + NELIM - 1 ) = NELIM_COL
        POSITION = POSITION + NELIM
        BUF_CB%CONTENT( POSITION: POSITION + NSLAVES - 1 ) = SLAVES
        POSITION = POSITION + NSLAVES
        POSITION = POSITION - IPOS
        IF ( POSITION * SIZEofINT .NE. SIZE ) THEN
          WRITE(*,*) 'Error in DMUMPS_BUF_SEND_ROOT_NELIM_INDICES:',
     &               'wrong estimated size'
           CALL MUMPS_ABORT()
        END IF
        CALL MPI_ISEND( BUF_CB%CONTENT( IPOS ), SIZE, 
     &                  MPI_PACKED,
     &                  DEST, ROOT_NELIM_INDICES, COMM,
     &                  BUF_CB%CONTENT( IREQ ), IERR )
        RETURN
        END SUBROUTINE DMUMPS_BUF_SEND_RTNELIND
        SUBROUTINE DMUMPS_BUF_SEND_ROOT2SON( ISON, NELIM_ROOT,
     &             DEST, COMM, IERR )
        IMPLICIT NONE
        INTEGER ISON, NELIM_ROOT, DEST, COMM, IERR
        INCLUDE 'mpif.h'
        INCLUDE 'mumps_tags.h'
        INTEGER IPOS, IREQ, SIZE
        INTEGER IONE
        INTEGER DEST2(1)
        PARAMETER ( IONE=1 )
        DEST2(1)=DEST
        IERR = 0
        SIZE = 2 * SIZEofINT
        CALL BUF_LOOK( BUF_SMALL, IPOS, IREQ, SIZE, IERR,
     &                 IONE, DEST2
     &               )
        IF ( IERR .LT. 0 ) THEN
          WRITE(*,*) 'Internal error 1 with small buffers '
          CALL MUMPS_ABORT()
        END IF
        IF ( IERR .LT. 0 ) THEN
           RETURN
        ENDIF
        BUF_SMALL%CONTENT( IPOS )     = ISON
        BUF_SMALL%CONTENT( IPOS + 1 ) = NELIM_ROOT
        CALL MPI_ISEND( BUF_SMALL%CONTENT( IPOS ), SIZE, 
     &                  MPI_PACKED,
     &                  DEST, ROOT_2SON, COMM,
     &                  BUF_SMALL%CONTENT( IREQ ), IERR )
        RETURN
        END SUBROUTINE DMUMPS_BUF_SEND_ROOT2SON
        SUBROUTINE DMUMPS_BUF_SEND_ROOT2SLAVE
     &  ( TOT_ROOT_SIZE, TOT_CONT2RECV, DEST, COMM, IERR )
        IMPLICIT NONE
        INTEGER TOT_ROOT_SIZE, TOT_CONT2RECV, DEST, COMM, IERR
        INCLUDE 'mpif.h'
        INCLUDE 'mumps_tags.h'
        INTEGER SIZE, IPOS, IREQ
        INTEGER IONE
        INTEGER DEST2(1)
        PARAMETER ( IONE=1 )
        IERR = 0
        DEST2(1) = DEST
        SIZE = 2 * SIZEofINT
        CALL BUF_LOOK( BUF_SMALL, IPOS, IREQ, SIZE, IERR,
     &                 IONE, DEST2
     &               )
        IF ( IERR .LT. 0 ) THEN
          WRITE(*,*) 'Internal error 2 with small buffers '
           CALL MUMPS_ABORT()
        END IF
        IF ( IERR .LT. 0 ) THEN
           RETURN
        ENDIF
        BUF_SMALL%CONTENT( IPOS     ) = TOT_ROOT_SIZE
        BUF_SMALL%CONTENT( IPOS + 1 ) = TOT_CONT2RECV
        CALL MPI_ISEND( BUF_SMALL%CONTENT( IPOS ), SIZE, 
     &                  MPI_PACKED,
     &                  DEST, ROOT_2SLAVE, COMM,
     &                  BUF_SMALL%CONTENT( IREQ ), IERR )
        RETURN
        END SUBROUTINE DMUMPS_BUF_SEND_ROOT2SLAVE
        SUBROUTINE DMUMPS_BUF_SEND_BACKVEC
     &             ( NRHS, INODE, W, LW, LD_W, DEST,MSGTAG,
     &               JBDEB, JBFIN, COMM,IERR )
        IMPLICIT NONE
        INTEGER NRHS, INODE,LW,COMM,IERR,DEST,MSGTAG, LD_W
        INTEGER, intent(in) :: JBDEB, JBFIN
        DOUBLE PRECISION W(LD_W, *)
        INCLUDE 'mpif.h'
        INTEGER SIZE, SIZE1, SIZE2
        INTEGER POSITION, IREQ, IPOS
        INTEGER IONE, K
        INTEGER DEST2(1)
        PARAMETER ( IONE=1 )
        IERR = 0
        DEST2(1) = DEST
        CALL MPI_PACK_SIZE( 4 , MPI_INTEGER, COMM, SIZE1, IERR )
        CALL MPI_PACK_SIZE( LW*NRHS, MPI_DOUBLE_PRECISION, COMM,
     &                      SIZE2, IERR )
        SIZE = SIZE1 + SIZE2
        CALL BUF_LOOK( BUF_CB, IPOS, IREQ, SIZE, IERR, 
     &                 IONE, DEST2
     &               )
        IF ( IERR .LT. 0 ) THEN
           RETURN
        ENDIF
        POSITION = 0
        CALL MPI_PACK( INODE, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE,
     &                        POSITION, COMM, IERR )
        CALL MPI_PACK( LW   , 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE,
     &                        POSITION, COMM, IERR )
        CALL MPI_PACK( JBDEB   , 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE,
     &                        POSITION, COMM, IERR )
        CALL MPI_PACK( JBFIN   , 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE,
     &                        POSITION, COMM, IERR )
        DO K=1, NRHS
        CALL MPI_PACK( W(1,K), LW, MPI_DOUBLE_PRECISION,
     &                        BUF_CB%CONTENT( IPOS ), SIZE,
     &                        POSITION, COMM, IERR )
        END DO
        CALL MPI_ISEND( BUF_CB%CONTENT( IPOS ), POSITION, MPI_PACKED,
     &                  DEST, MSGTAG, COMM,
     &                  BUF_CB%CONTENT( IREQ ), IERR )
        IF ( SIZE .LT. POSITION ) THEN
          WRITE(*,*) 'Try_update: SIZE, POSITION = ',
     &               SIZE, POSITION
          CALL MUMPS_ABORT()
        END IF
        IF ( SIZE .NE. POSITION ) CALL BUF_ADJUST( BUF_CB, POSITION )
        RETURN
        END SUBROUTINE DMUMPS_BUF_SEND_BACKVEC
        SUBROUTINE DMUMPS_BUF_SEND_UPDATE_LOAD
     &             ( BDC_SBTR,BDC_MEM,BDC_MD, COMM, NPROCS, LOAD,
     &               MEM,SBTR_CUR,
     &               LU_USAGE,
#if ! defined(OLD_LOAD_MECHANISM)
     &               FUTURE_NIV2,
#endif
     &               MYID, IERR)
        IMPLICIT NONE
        INTEGER COMM, NPROCS, MYID, IERR
#if ! defined(OLD_LOAD_MECHANISM)
        INTEGER FUTURE_NIV2(NPROCS)
#endif
        DOUBLE PRECISION LU_USAGE
        DOUBLE PRECISION LOAD
        DOUBLE PRECISION MEM,SBTR_CUR
        LOGICAL BDC_MEM,BDC_SBTR,BDC_MD
        INCLUDE 'mpif.h'
        INCLUDE 'mumps_tags.h'
        INTEGER POSITION, IREQ, IPOS, SIZE1, SIZE2, SIZE
        INTEGER I, NDEST, IDEST, IPOSMSG, WHAT, NREALS
        INTEGER IZERO
        INTEGER MYID2(1)
        PARAMETER ( IZERO=0 )
        IERR = 0
        MYID2(1) = MYID
        NDEST = NPROCS - 1
#if ! defined(OLD_LOAD_MECHANISM)
        NDEST = 0
        DO I = 1, NPROCS
          IF ( I .NE. MYID + 1 .AND. FUTURE_NIV2(I).NE.0) THEN
            NDEST = NDEST + 1
          ENDIF
        ENDDO
#endif
        IF ( NDEST .eq. 0 ) THEN
           RETURN
        ENDIF
        CALL MPI_PACK_SIZE( 1 + (NDEST-1) * OVHSIZE, 
     &                       MPI_INTEGER, COMM,
     &                       SIZE1, IERR )
        NREALS = 1
        IF (BDC_MEM) THEN
          NREALS = 2
        ENDIf
        IF (BDC_SBTR)THEN
          NREALS = 3
        ENDIF
        IF(BDC_MD)THEN
           NREALS=NREALS+1
        ENDIF
        CALL MPI_PACK_SIZE( NREALS, MPI_DOUBLE_PRECISION,
     &                      COMM, SIZE2, IERR )
        SIZE = SIZE1 + SIZE2
        CALL BUF_LOOK( BUF_LOAD, IPOS, IREQ, SIZE, IERR, 
     &                  IZERO, MYID2 
     &               )
        IF ( IERR .LT. 0 ) THEN
           RETURN
        ENDIF
        BUF_LOAD%ILASTMSG = BUF_LOAD%ILASTMSG + ( NDEST - 1 ) * OVHSIZE
        IPOS = IPOS - OVHSIZE
        DO IDEST = 1, NDEST - 1
          BUF_LOAD%CONTENT( IPOS + ( IDEST - 1 ) * OVHSIZE ) =
     &    IPOS + IDEST * OVHSIZE
        END DO
        BUF_LOAD%CONTENT( IPOS + ( NDEST - 1 ) * OVHSIZE ) = 0
        IPOSMSG = IPOS + OVHSIZE * NDEST
        WHAT = 0  
        POSITION = 0
        CALL MPI_PACK( WHAT, 1, MPI_INTEGER,
     &                 BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     &                 POSITION, COMM, IERR )
        CALL MPI_PACK( LOAD, 1, MPI_DOUBLE_PRECISION,
     &                 BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     &                 POSITION, COMM, IERR )
        IF (BDC_MEM) THEN
          CALL MPI_PACK( MEM, 1, MPI_DOUBLE_PRECISION,
     &                   BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     &                   POSITION, COMM, IERR )
        END IF
        IF (BDC_SBTR) THEN
          CALL MPI_PACK( SBTR_CUR, 1, MPI_DOUBLE_PRECISION,
     &                   BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     &                   POSITION, COMM, IERR )
        END IF
        IF(BDC_MD)THEN
           CALL MPI_PACK( LU_USAGE, 1, MPI_DOUBLE_PRECISION,
     &          BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     &          POSITION, COMM, IERR )
        ENDIF
        IDEST = 0
        DO I = 0, NPROCS - 1
#if ! defined(OLD_LOAD_MECHANISM)
        IF ( I .NE. MYID .AND. FUTURE_NIV2(I+1) .NE. 0) THEN
#else
        IF ( I .ne. MYID ) THEN
#endif
            IDEST = IDEST + 1
            CALL MPI_ISEND( BUF_LOAD%CONTENT( IPOSMSG ),
     &                      POSITION, MPI_PACKED, I,
     &                      UPDATE_LOAD, COMM,
     &                      BUF_LOAD%CONTENT( IREQ+(IDEST-1)*OVHSIZE ),
     &                      IERR )
          END IF
        END DO
        SIZE = SIZE - ( NDEST - 1 ) * OVHSIZE * SIZEofINT
        IF ( SIZE .LT. POSITION ) THEN
          WRITE(*,*) ' Error in DMUMPS_BUF_SEND_UPDATE_LOAD'
          WRITE(*,*) ' Size,position=',SIZE,POSITION
          CALL MUMPS_ABORT()
        END IF
        IF ( SIZE .NE. POSITION )
     &  CALL BUF_ADJUST( BUF_LOAD, POSITION )
        RETURN
        END SUBROUTINE DMUMPS_BUF_SEND_UPDATE_LOAD
        SUBROUTINE DMUMPS_BUF_BROADCAST
     &             ( WHAT, COMM, NPROCS, 
#if ! defined(OLD_LOAD_MECHANISM)
     &               FUTURE_NIV2,
#endif
     &               LOAD,UPD_LOAD,
     &               MYID, IERR)
        IMPLICIT NONE
        INTEGER COMM, NPROCS, MYID, IERR, WHAT
        DOUBLE PRECISION LOAD,UPD_LOAD
        INCLUDE 'mpif.h'
        INCLUDE 'mumps_tags.h'
        INTEGER POSITION, IREQ, IPOS, SIZE1, SIZE2, SIZE
        INTEGER I, NDEST, IDEST, IPOSMSG, NREALS
        INTEGER IZERO
        INTEGER MYID2(1)
#if ! defined(OLD_LOAD_MECHANISM)
        INTEGER FUTURE_NIV2(NPROCS)
#endif
        PARAMETER ( IZERO=0 )
        IERR = 0
        IF (WHAT .NE. 2 .AND. WHAT .NE. 3 .AND.
     &       WHAT.NE.6.AND. WHAT.NE.8 .AND.WHAT.NE.9.AND.
     &       WHAT.NE.17) THEN
          WRITE(*,*)
     &  "Internal error 1 in DMUMPS_BUF_BROADCAST",WHAT
        END IF
        MYID2(1) = MYID
        NDEST = NPROCS - 1
#if ! defined(OLD_LOAD_MECHANISM)
        NDEST = 0
        DO I = 1, NPROCS
          IF ( I .NE. MYID + 1 .AND. FUTURE_NIV2(I).NE.0) THEN
            NDEST = NDEST + 1
          ENDIF
        ENDDO
#endif
        IF ( NDEST .eq. 0 ) THEN
           RETURN
        ENDIF
        CALL MPI_PACK_SIZE( 1 + (NDEST-1) * OVHSIZE, 
     &                       MPI_INTEGER, COMM,
     &                       SIZE1, IERR )
        IF((WHAT.NE.17).AND.(WHAT.NE.10))THEN
           NREALS = 1
        ELSE
           NREALS = 2
        ENDIF
        CALL MPI_PACK_SIZE( NREALS, MPI_DOUBLE_PRECISION,
     &                      COMM, SIZE2, IERR )
        SIZE = SIZE1 + SIZE2
        CALL BUF_LOOK( BUF_LOAD, IPOS, IREQ, SIZE, IERR, 
     &                  IZERO, MYID2 
     &               )
        IF ( IERR .LT. 0 ) THEN
           RETURN
        ENDIF
        BUF_LOAD%ILASTMSG = BUF_LOAD%ILASTMSG + ( NDEST - 1 ) * OVHSIZE
        IPOS = IPOS - OVHSIZE
        DO IDEST = 1, NDEST - 1
          BUF_LOAD%CONTENT( IPOS + ( IDEST - 1 ) * OVHSIZE ) =
     &    IPOS + IDEST * OVHSIZE
        END DO
        BUF_LOAD%CONTENT( IPOS + ( NDEST - 1 ) * OVHSIZE ) = 0
        IPOSMSG = IPOS + OVHSIZE * NDEST
        POSITION = 0
        CALL MPI_PACK( WHAT, 1, MPI_INTEGER,
     &                 BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     &                 POSITION, COMM, IERR )
        CALL MPI_PACK( LOAD, 1, MPI_DOUBLE_PRECISION,
     &                 BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     &                 POSITION, COMM, IERR )
        IF((WHAT.EQ.17).OR.(WHAT.EQ.10))THEN
           CALL MPI_PACK( UPD_LOAD, 1, MPI_DOUBLE_PRECISION,
     &          BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     &          POSITION, COMM, IERR )
        ENDIF
        IDEST = 0
        DO I = 0, NPROCS - 1
#if ! defined(OLD_LOAD_MECHANISM)
          IF ( I .NE. MYID .AND. FUTURE_NIV2(I+1) .NE. 0) THEN
#else
          IF ( I .ne. MYID ) THEN
#endif
            IDEST = IDEST + 1
            CALL MPI_ISEND( BUF_LOAD%CONTENT( IPOSMSG ),
     &                      POSITION, MPI_PACKED, I,
     &                      UPDATE_LOAD, COMM,
     &                      BUF_LOAD%CONTENT( IREQ+(IDEST-1)*OVHSIZE ),
     &                      IERR )
          END IF
        END DO
        SIZE = SIZE - ( NDEST - 1 ) * OVHSIZE * SIZEofINT
        IF ( SIZE .LT. POSITION ) THEN
          WRITE(*,*) ' Error in DMUMPS_BUF_BROADCAST'
          WRITE(*,*) ' Size,position=',SIZE,POSITION
          CALL MUMPS_ABORT()
        END IF
        IF ( SIZE .NE. POSITION )
     &  CALL BUF_ADJUST( BUF_LOAD, POSITION )
        RETURN
        END SUBROUTINE DMUMPS_BUF_BROADCAST
        SUBROUTINE DMUMPS_BUF_SEND_FILS
     &             ( WHAT, COMM, NPROCS,
     &               FATHER_NODE,INODE,NCB,K81,
     &               MYID,REMOTE, IERR)
        IMPLICIT NONE
        INTEGER COMM, NPROCS, MYID, IERR, WHAT,REMOTE
        INTEGER FATHER_NODE,INODE
        INCLUDE 'mpif.h'
        INCLUDE 'mumps_tags.h'
        INTEGER POSITION, IREQ, IPOS, SIZE
        INTEGER NDEST, IDEST, IPOSMSG
        INTEGER IZERO,NCB,K81
        INTEGER MYID2(1)
        PARAMETER ( IZERO=0 )
        MYID2(1) = MYID
        NDEST = 1
        IF ( NDEST .eq. 0 ) THEN
           RETURN
        ENDIF
        IF((K81.EQ.2).OR.(K81.EQ.3))THEN
           CALL MPI_PACK_SIZE( 4 + OVHSIZE, 
     &          MPI_INTEGER, COMM,
     &          SIZE, IERR )
        ELSE
           CALL MPI_PACK_SIZE( 2 + OVHSIZE, 
     &          MPI_INTEGER, COMM,
     &          SIZE, IERR )
        ENDIF
        CALL BUF_LOOK( BUF_LOAD, IPOS, IREQ, SIZE, IERR, 
     &                  IZERO, MYID2 
     &               )
        IF ( IERR .LT. 0 ) THEN
           RETURN
        ENDIF
        BUF_LOAD%ILASTMSG = BUF_LOAD%ILASTMSG + ( NDEST - 1 ) * OVHSIZE
        IPOS = IPOS - OVHSIZE
        DO IDEST = 1, NDEST - 1
          BUF_LOAD%CONTENT( IPOS + ( IDEST - 1 ) * OVHSIZE ) =
     &    IPOS + IDEST * OVHSIZE
        END DO
        BUF_LOAD%CONTENT( IPOS + ( NDEST - 1 ) * OVHSIZE ) = 0
        IPOSMSG = IPOS + OVHSIZE * NDEST
        POSITION = 0
        CALL MPI_PACK( WHAT, 1, MPI_INTEGER,
     &                 BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     &                 POSITION, COMM, IERR )
        CALL MPI_PACK( FATHER_NODE, 1, MPI_INTEGER,
     &                 BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     &                 POSITION, COMM, IERR )
        IF((K81.EQ.2).OR.(K81.EQ.3))THEN
           CALL MPI_PACK( INODE, 1, MPI_INTEGER,
     &          BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     &          POSITION, COMM, IERR )
           CALL MPI_PACK( NCB, 1, MPI_INTEGER,
     &          BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     &          POSITION, COMM, IERR )
        ENDIF
        IDEST = 1
        CALL MPI_ISEND( BUF_LOAD%CONTENT( IPOSMSG ),
     &                 POSITION, MPI_PACKED, REMOTE,
     &                 UPDATE_LOAD, COMM,
     &                 BUF_LOAD%CONTENT( IREQ+(IDEST-1)*OVHSIZE ),
     &                 IERR )
        SIZE = SIZE - ( NDEST - 1 ) * OVHSIZE * SIZEofINT
        IF ( SIZE .LT. POSITION ) THEN
          WRITE(*,*) ' Error in DMUMPS_BUF_SEND_FILS'
          WRITE(*,*) ' Size,position=',SIZE,POSITION
          CALL MUMPS_ABORT()
        END IF
        IF ( SIZE .NE. POSITION )
     &  CALL BUF_ADJUST( BUF_LOAD, POSITION )
        RETURN
        END SUBROUTINE DMUMPS_BUF_SEND_FILS
        SUBROUTINE DMUMPS_BUF_SEND_NOT_MSTR( COMM, MYID, NPROCS,
     &  MAX_SURF_MASTER,IERR)
        IMPLICIT NONE
        INCLUDE 'mpif.h'
        INCLUDE 'mumps_tags.h'
        INTEGER IPOS, IREQ, IDEST, IPOSMSG, POSITION, I
        INTEGER COMM, MYID, IERR, NPROCS
        DOUBLE PRECISION MAX_SURF_MASTER
        INTEGER IZERO
        INTEGER MYID2(1)
        PARAMETER ( IZERO=0 )
        INTEGER NDEST, NINTS, NREALS, SIZE, SIZE1, SIZE2
        INTEGER WHAT
        IERR = 0
        MYID2(1) = MYID
        NDEST = NPROCS - 1
        NINTS = 1 + ( NDEST-1 ) * OVHSIZE
        NREALS = 1
        CALL MPI_PACK_SIZE( NINTS, 
     &                       MPI_INTEGER, COMM,
     &                       SIZE1, IERR )
        CALL MPI_PACK_SIZE( NREALS,
     &                       MPI_DOUBLE_PRECISION, COMM,
     &                       SIZE2, IERR )
        SIZE=SIZE1+SIZE2
        CALL BUF_LOOK( BUF_LOAD, IPOS, IREQ, SIZE, IERR,
     &       IZERO, MYID2 )
        IF ( IERR .LT. 0 ) THEN
           RETURN
        ENDIF
        BUF_LOAD%ILASTMSG = BUF_LOAD%ILASTMSG + ( NDEST - 1 ) * OVHSIZE
        IPOS = IPOS - OVHSIZE
        DO IDEST = 1, NDEST - 1
          BUF_LOAD%CONTENT( IPOS + ( IDEST - 1 ) * OVHSIZE ) =
     &    IPOS + IDEST * OVHSIZE
        END DO
        BUF_LOAD%CONTENT( IPOS + ( NDEST - 1 ) * OVHSIZE ) = 0
        IPOSMSG = IPOS + OVHSIZE * NDEST
        POSITION = 0
        WHAT = 4
        CALL MPI_PACK( WHAT, 1, MPI_INTEGER,
     &      BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     &      POSITION, COMM, IERR )
        CALL MPI_PACK( MAX_SURF_MASTER, 1, MPI_DOUBLE_PRECISION,
     &      BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     &      POSITION, COMM, IERR )
        IDEST = 0
        DO I = 0, NPROCS - 1
           IF ( I .ne. MYID ) THEN
              IDEST = IDEST + 1
              CALL MPI_ISEND( BUF_LOAD%CONTENT( IPOSMSG ),
     &             POSITION, MPI_PACKED, I,
     &             UPDATE_LOAD, COMM,
     &             BUF_LOAD%CONTENT( IREQ+(IDEST-1)*OVHSIZE ),
     &             IERR )
           END IF
        END DO
        SIZE = SIZE - ( NDEST - 1 ) * OVHSIZE * SIZEofINT
        IF ( SIZE .LT. POSITION ) THEN
          WRITE(*,*) ' Error in DMUMPS_BUF_BCAST_ARRAY'
          WRITE(*,*) ' Size,position=',SIZE,POSITION
          CALL MUMPS_ABORT()
        END IF
        IF ( SIZE .NE. POSITION )
     &  CALL BUF_ADJUST( BUF_LOAD, POSITION )
        RETURN
        END SUBROUTINE DMUMPS_BUF_SEND_NOT_MSTR
        SUBROUTINE DMUMPS_BUF_BCAST_ARRAY( BDC_MEM,
     &      COMM, MYID, NPROCS,
#if ! defined(OLD_LOAD_MECHANISM)
     &      FUTURE_NIV2,
#endif
     &      NSLAVES,
     &      LIST_SLAVES,INODE,
     &      MEM_INCREMENT, FLOPS_INCREMENT,CB_BAND, WHAT,
     &      IERR )
        IMPLICIT NONE
        INCLUDE 'mpif.h'
        INCLUDE 'mumps_tags.h'
        LOGICAL BDC_MEM
        INTEGER COMM, MYID, NPROCS, NSLAVES, IERR
#if ! defined(OLD_LOAD_MECHANISM)
        INTEGER FUTURE_NIV2(NPROCS)
#endif
        INTEGER LIST_SLAVES(NSLAVES),INODE
        DOUBLE PRECISION MEM_INCREMENT(NSLAVES)
        DOUBLE PRECISION FLOPS_INCREMENT(NSLAVES)
        DOUBLE PRECISION CB_BAND(NSLAVES)
        INTEGER NDEST, NINTS, NREALS, SIZE1, SIZE2, SIZE
        INTEGER IPOS, IPOSMSG, IREQ, POSITION
        INTEGER I, IDEST, WHAT
        INTEGER IZERO
        INTEGER MYID2(1)
        PARAMETER ( IZERO=0 )
        MYID2(1)=MYID
        IERR = 0
#if ! defined(OLD_LOAD_MECHANISM)
        NDEST = 0
        DO I = 1, NPROCS
          IF ( I .NE. MYID + 1 .AND. FUTURE_NIV2(I).NE.0) THEN
            NDEST = NDEST + 1
          ENDIF
        ENDDO
#else
        NDEST = NPROCS - 1
#endif
        IF ( NDEST == 0 ) THEN
           RETURN
        ENDIF
        NINTS = 2 +  NSLAVES + ( NDEST - 1 ) * OVHSIZE + 1
        NREALS = NSLAVES
        IF (BDC_MEM) NREALS = NREALS + NSLAVES
        IF(WHAT.EQ.19) THEN 
           NREALS = NREALS + NSLAVES
        ENDIF
        CALL MPI_PACK_SIZE( NINTS, 
     &                       MPI_INTEGER, COMM,
     &                       SIZE1, IERR )
        CALL MPI_PACK_SIZE( NREALS, MPI_DOUBLE_PRECISION,
     &       COMM, SIZE2, IERR )
        SIZE = SIZE1+SIZE2
        CALL BUF_LOOK( BUF_LOAD, IPOS, IREQ, SIZE, IERR,
     &       IZERO, MYID2 )
        IF ( IERR .LT. 0 ) THEN
           RETURN
        ENDIF
        BUF_LOAD%ILASTMSG = BUF_LOAD%ILASTMSG + ( NDEST - 1 ) * OVHSIZE
        IPOS = IPOS - OVHSIZE
        DO IDEST = 1, NDEST - 1
          BUF_LOAD%CONTENT( IPOS + ( IDEST - 1 ) * OVHSIZE ) =
     &    IPOS + IDEST * OVHSIZE
        END DO
        BUF_LOAD%CONTENT( IPOS + ( NDEST - 1 ) * OVHSIZE ) = 0
        IPOSMSG = IPOS + OVHSIZE * NDEST
        POSITION = 0
        CALL MPI_PACK( WHAT, 1, MPI_INTEGER,
     &      BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     &      POSITION, COMM, IERR )
        CALL MPI_PACK( NSLAVES, 1, MPI_INTEGER,
     &      BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     &      POSITION, COMM, IERR )
        CALL MPI_PACK( INODE, 1, MPI_INTEGER,
     &      BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     &      POSITION, COMM, IERR )
        CALL MPI_PACK( LIST_SLAVES, NSLAVES, MPI_INTEGER,
     &      BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     &      POSITION, COMM, IERR )
        CALL MPI_PACK( FLOPS_INCREMENT, NSLAVES,
     &      MPI_DOUBLE_PRECISION,
     &      BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     &      POSITION, COMM, IERR )
        IF (BDC_MEM) THEN
          CALL MPI_PACK( MEM_INCREMENT, NSLAVES,
     &      MPI_DOUBLE_PRECISION,
     &      BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     &      POSITION, COMM, IERR )
        END IF
        IF(WHAT.EQ.19)THEN
           CALL MPI_PACK( CB_BAND, NSLAVES,
     &          MPI_DOUBLE_PRECISION,
     &          BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     &          POSITION, COMM, IERR )
        ENDIF
        IDEST = 0
        DO I = 0, NPROCS - 1
#if ! defined(OLD_LOAD_MECHANISM)
        IF ( I .NE. MYID .AND. FUTURE_NIV2(I+1) .NE. 0) THEN
#else
        IF ( I .NE. MYID ) THEN
#endif
            IDEST = IDEST + 1
            CALL MPI_ISEND( BUF_LOAD%CONTENT( IPOSMSG ),
     &                      POSITION, MPI_PACKED, I,
     &                      UPDATE_LOAD, COMM,
     &                      BUF_LOAD%CONTENT( IREQ+(IDEST-1)*OVHSIZE ),
     &                      IERR )
          END IF
        END DO
        SIZE = SIZE - ( NDEST - 1 ) * OVHSIZE * SIZEofINT
        IF ( SIZE .LT. POSITION ) THEN
          WRITE(*,*) ' Error in DMUMPS_BUF_BCAST_ARRAY'
          WRITE(*,*) ' Size,position=',SIZE,POSITION
          CALL MUMPS_ABORT()
        END IF
        IF ( SIZE .NE. POSITION )
     &  CALL BUF_ADJUST( BUF_LOAD, POSITION )
        RETURN
        END SUBROUTINE DMUMPS_BUF_BCAST_ARRAY
        SUBROUTINE DMUMPS_BUF_DIST_IRECV_SIZE
     &             ( DMUMPS_LBUFR_BYTES)
        IMPLICIT NONE
        INTEGER DMUMPS_LBUFR_BYTES 
        SIZE_RBUF_BYTES = DMUMPS_LBUFR_BYTES
        RETURN
      END SUBROUTINE DMUMPS_BUF_DIST_IRECV_SIZE
      END MODULE DMUMPS_COMM_BUFFER