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
C
      SUBROUTINE ZMUMPS_ANA_DRIVER(id)
      USE ZMUMPS_LOAD
      USE MUMPS_STATIC_MAPPING
      USE ZMUMPS_STRUC_DEF
      USE TOOLS_COMMON
      USE ZMUMPS_PARALLEL_ANALYSIS
      IMPLICIT NONE
C     
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER IERR, MASTER
      PARAMETER( MASTER = 0 )
C
C     Purpose
C     =======
C
C     Performs analysis and (if required) Max-trans on the master, then
C     broadcasts information to the slaves. Also includes mapping.
C     
C     
C     Parameters
C     ==========
C     
      TYPE(ZMUMPS_STRUC), TARGET :: id
C     
C     Local variables
C     ===============
C     
C     
C     Pointers inside integer array, various data
      INTEGER LIW, IKEEP, FILS, FRERE, PTRAR, NFSIZ
      INTEGER NE, NA
      INTEGER I, allocok
      INTEGER MAXIS1_CHECK
C     Other locals
      INTEGER NB_NIV2, IDEST
      INTEGER :: STATUS(MPI_STATUS_SIZE)
      INTEGER LOCAL_M, LOCAL_N
      INTEGER numroc
      EXTERNAL numroc
      INTEGER IRANK
      INTEGER MP, LP, MPG
      LOGICAL PROK, PROKG, LISTVAR_SCHUR_2BE_FREED, LPOK
      INTEGER SIZE_SCHUR_PASSED
      INTEGER SBUF_SEND, SBUF_REC, TOTAL_MBYTES
      INTEGER(8) SBUF_RECOLD8, MIN_BUF_SIZE8
      INTEGER MIN_BUF_SIZE
      INTEGER(8) MAX_SIZE_FACTOR_TMP
      INTEGER LEAF, INODE, ISTEP, INN, LPTRAR
      INTEGER NBLEAF, NBROOT, MYROW_CHECK, INIV2
C     to store the size of the sequencial peak of stack
C     (or an estimation for not calling REORDER_TREE_N )
      DOUBLE PRECISION PEAK
C     
C     INTEGER WORKSPACE 
C     
      INTEGER, ALLOCATABLE, DIMENSION(:) :: PAR2_NODES
      INTEGER, DIMENSION(:), ALLOCATABLE :: IWtemp
      INTEGER, DIMENSION(:), ALLOCATABLE :: XNODEL, NODEL
      INTEGER, DIMENSION(:), POINTER :: SSARBR
C     Element matrix entry
      INTEGER, POINTER ::  NELT, LELTVAR
      INTEGER, DIMENSION(:), POINTER :: KEEP,INFO, INFOG
      INTEGER(8), DIMENSION(:), POINTER :: KEEP8
      INTEGER(8)                   :: ENTRIES_IN_FACTORS_LOC_MASTERS
      DOUBLE PRECISION, DIMENSION(:), POINTER :: RINFO
      DOUBLE PRECISION, DIMENSION(:), POINTER :: RINFOG
      INTEGER, DIMENSION(:), POINTER :: ICNTL
      LOGICAL I_AM_SLAVE, PERLU_ON, COND
      INTEGER :: OOC_STAT
      INTEGER MUMPS_TYPENODE, MUMPS_PROCNODE
      EXTERNAL MUMPS_TYPENODE, MUMPS_PROCNODE
      INTEGER K,J, IFS
      INTEGER SIZE_TEMP_MEM,SIZE_DEPTH_FIRST,SIZE_COST_TRAV
      LOGICAL IS_BUILD_LOAD_MEM_CALLED
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: TEMP_MEM
      INTEGER, DIMENSION (:,:), ALLOCATABLE :: TEMP_ROOT
      INTEGER, DIMENSION (:,:), ALLOCATABLE :: TEMP_LEAF
      INTEGER, DIMENSION (:,:), ALLOCATABLE :: TEMP_SIZE
      INTEGER, DIMENSION (:), ALLOCATABLE :: DEPTH_FIRST
      INTEGER, DIMENSION (:), ALLOCATABLE :: DEPTH_FIRST_SEQ
      INTEGER, DIMENSION (:), ALLOCATABLE :: SBTR_ID
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: COST_TRAV_TMP
      INTEGER(8) :: TOTAL_BYTES
      INTEGER, POINTER, DIMENSION(:) ::  WORK1PTR, WORK2PTR,
     &     NFSIZPTR,
     &     FILSPTR,
     &     FREREPTR
      ! Used because of multithreaded SIM_NP_
      INTEGER :: locMYID, locMYID_NODES
      LOGICAL, POINTER :: locI_AM_CAND(:)
      INTEGER(kind=8) :: N8, NZ8, LIW8
C
C  Beginning of executable statements
C
      IS_BUILD_LOAD_MEM_CALLED=.FALSE.
      KEEP   => id%KEEP
      KEEP8  => id%KEEP8
      INFO   => id%INFO
      RINFO  => id%RINFO
      INFOG  => id%INFOG
      RINFOG => id%RINFOG
      ICNTL  => id%ICNTL
      NELT    => id%NELT
      LELTVAR => id%LELTVAR
      KEEP8(24) = 0_8  ! reinitialize last used size of WK_USER
C     -------------------------------------
C     Depending on the type of parallelism,
C     the master can now (soon) potentially
C     have the role of a slave
C     -------------------------------------
      I_AM_SLAVE = ( id%MYID .ne. MASTER  .OR.
     &     ( id%MYID .eq. MASTER .AND.
     &     id%KEEP(46) .eq. 1 ) )
      LP  = ICNTL( 1 )
      MP  = ICNTL( 2 )
      MPG = ICNTL( 3 )
C     LP     : errors
C     MP     : INFO
      LPOK  = ((LP.GT.0).AND.(id%ICNTL(4).GE.1))
      PROK  = (( MP  .GT. 0 ).AND.(ICNTL(4).GE.2))
      PROKG = ( MPG .GT. 0 .and. id%MYID .eq. MASTER )
      PROKG = (PROKG.AND.(ICNTL(4).GE.2))
      IF ( PROK ) THEN
         IF ( KEEP(50) .eq. 0 ) THEN
            WRITE(MP, '(A)') 'L U Solver for unsymmetric matrices'
         ELSE IF ( KEEP(50) .eq. 1 ) THEN
            WRITE(MP, '(A)') 
     & 'L D L^T Solver for symmetric positive definite matrices'
         ELSE
            WRITE(MP, '(A)') 
     &           'L D L^T Solver for general symmetric matrices'
         END IF
         IF ( KEEP(46) .eq. 1 ) THEN
            WRITE(MP, '(A)') 'Type of parallelism: Working host'
         ELSE
            WRITE(MP, '(A)') 'Type of parallelism: Host not working'
         END IF
      END IF
      IF ( PROKG .AND. (MP.NE.MPG)) THEN
         IF ( KEEP(50) .eq. 0 ) THEN
            WRITE(MPG, '(A)') 'L U Solver for unsymmetric matrices'
         ELSE IF ( KEEP(50) .eq. 1 ) THEN
            WRITE(MPG, '(A)') 
     & 'L D L^T Solver for symmetric positive definite matrices'
         ELSE
            WRITE(MPG, '(A)') 
     &           'L D L^T Solver for general symmetric matrices'
         END IF
         IF ( KEEP(46) .eq. 1 ) THEN
            WRITE(MPG, '(A)') 'Type of parallelism: Working host'
         ELSE
            WRITE(MPG, '(A)') 'Type of parallelism: Host not working'
         END IF
      END IF
      IF (PROK) WRITE( MP, 110 )
      IF (PROKG .AND. (MPG.NE.MP)) WRITE( MPG, 110 )
C
C     Decode API (ICNTL parameters, mainly)
C     and check consistency of the KEEP array.
C     Note: ZMUMPS_ANA_CHECK_KEEP also sets
C     some INFOG parameters
      CALL ZMUMPS_ANA_CHECK_KEEP(id)
      CALL MUMPS_PROPINFO( ICNTL(1), INFO(1),
     &     id%COMM, id%MYID )
      IF ( INFO(1) .LT. 0 ) RETURN
C     -------------------------------------------
C     Broadcast KEEP(60) since we need to broadcast
C     related information
C     ------------------------------------------
      CALL MPI_BCAST( KEEP(60), 1, MPI_INTEGER, MASTER, id%COMM, IERR )
      IF (id%KEEP(60) .EQ. 2 .or. id%KEEP(60). EQ. 3) THEN
         CALL MPI_BCAST( id%NPROW, 1,
     &        MPI_INTEGER, MASTER, id%COMM, IERR )
         CALL MPI_BCAST( id%NPCOL, 1,
     &        MPI_INTEGER, MASTER, id%COMM, IERR )
         CALL MPI_BCAST( id%MBLOCK, 1,
     &        MPI_INTEGER, MASTER, id%COMM, IERR )
         CALL MPI_BCAST( id%NBLOCK, 1,
     &        MPI_INTEGER, MASTER, id%COMM, IERR )
C     Note that ZMUMPS_INIT_ROOT_ANA will
C     then use that information.
      ENDIF
      CALL MUMPS_PROPINFO( ICNTL(1), INFO(1),
     &     id%COMM, id%MYID )
      IF ( INFO(1) .LT. 0 ) RETURN
C     ----------------------------------------------
C     Broadcast KEEP(54) now to know if the
C     structure of the graph is intially distributed
C     and should be assembled on the master
C     Broadcast KEEP(55) now to know if the
C     matrix is in assembled or elemental format
C     ----------------------------------------------
      CALL MPI_BCAST( KEEP(54), 2, MPI_INTEGER, MASTER, id%COMM, IERR )
C     ----------------------------------------------
C     Broadcast KEEP(69) now to know if
C     we will need to communicate during analysis
C     ----------------------------------------------
      CALL MPI_BCAST( KEEP(69), 1, MPI_INTEGER, MASTER, id%COMM, IERR )
C     ----------------------------------------------
C     Broadcast Out of core strategy (used only on master so far)
C     ----------------------------------------------
      CALL MPI_BCAST( KEEP(201), 1, MPI_INTEGER, MASTER, id%COMM, IERR )
C     ----------------------------------------------
C     Broadcast analysis strategy (used only on master so far)
C     ----------------------------------------------
      CALL MPI_BCAST( KEEP(244), 1, MPI_INTEGER, MASTER, id%COMM, IERR )
C     ---------------------------
C     Fwd in facto
C     Broadcast KEEP(251,252,253) defined on master so far
      CALL MPI_BCAST( KEEP(251), 3, MPI_INTEGER,MASTER,id%COMM,IERR)
C     ----------------------------------------------
C     Broadcast N 
C     ----------------------------------------------
      CALL MPI_BCAST( id%N, 1, MPI_INTEGER, MASTER, id%COMM, IERR )
C     ----------------------------------------------
C     Broadcast NZ for assembled entry
C     ----------------------------------------------
      IF ( KEEP(55) .EQ. 0) THEN
         IF ( KEEP(54) .eq. 3 ) THEN
C     Compute total number of non-zeros
          CALL MPI_ALLREDUCE( id%NZ_loc, id%NZ, 1, MPI_INTEGER, 
     &       MPI_SUM, id%COMM, IERR )
         ELSE
C     Broadcast NZ from the master node
            CALL MPI_BCAST( id%NZ, 1, MPI_INTEGER, MASTER,
     &           id%COMM, IERR )
         END IF
      ELSE
C     Broadcast NA_ELT for elemental entry
         CALL MPI_BCAST( id%NA_ELT, 1, MPI_INTEGER, MASTER,
     &        id%COMM, IERR )
      ENDIF
      IF ( associated(id%MEM_DIST) ) deallocate( id%MEM_DIST )
      allocate( id%MEM_DIST( 0:id%NSLAVES-1 ), STAT=IERR )
      IF ( IERR .GT. 0 ) THEN
         INFO(1) = -7
         INFO(2) = id%NSLAVES
         IF ( LPOK ) THEN
            WRITE(LP, 150) 'MEM_DIST'
         END IF
      END IF
      CALL MUMPS_PROPINFO( ICNTL(1), INFO(1),
     &     id%COMM, id%MYID )
      IF ( INFO(1) .LT. 0 ) RETURN
      id%MEM_DIST(0:id%NSLAVES-1) = 0
      CALL MUMPS_INIT_ARCH_PARAMETERS(
     &     id%COMM,id%COMM_NODES,KEEP(69),KEEP(46),
     &     id%NSLAVES,id%MEM_DIST,INFO)
C     ========================
C     Write problem to a file,
C     if requested by the user
C     ========================
      CALL ZMUMPS_DUMP_PROBLEM(id)
C
C     ====================================================
C     TEST FOR SEQUENTIAL OR PARALLEL ANALYSIS (KEEP(244))
C     ====================================================
      IF (KEEP(244) .EQ. 1) THEN
C     Sequential analysis
         IF ( KEEP(54) .eq. 3 ) THEN
C        -----------------------------------------------
C        Collect on the host -- if matrix is distributed
C        at analysis -- all integer information.
C        -----------------------------------------------
            CALL ZMUMPS_GATHER_MATRIX(id)
         END IF
C        ************************************************
C        BEGINNING OF MASTER CODE FOR SEQUENTIAL ANALYSIS
C        ************************************************
         IF ( id%MYID .eq. MASTER ) THEN
 1234       CONTINUE
C     new case for :
C     - MTRANS scaling
C     - LDLT facto
C     enter if K23 = 7 and scaling required and maxtrans allowed
            IF ( ( (KEEP(23) .NE. 0) .AND.
     &           ( (KEEP(23).NE.7) .OR. KEEP(50).EQ. 2 ) )
     &           .OR.
     &           ( associated(id%A) .AND. KEEP(52) .EQ. 77 .AND.
     &           (KEEP(50).EQ.2))
     &        .OR.
     &           KEEP(52) .EQ. -2 ) THEN
C     MAXIMUM TRANSVERSAL ALGORITHM called on original matrix.
C     KEEP(23) = 7 means that automatic choice
C     of max trans value will be done during Analysis.
C     We compute a permutation of  the original matrix to have zero free diagonal
C     the col. Permutation is held in IS1(1, ...,N).  
C     Max-trans (ZMUMPS_ANA_O) is not used for element entry.
               IF (.not.associated(id%A)) THEN
C     -- If maxtrans is required and A not allocated then reset 
C     -- it to structural based maxtrans.
                  IF (KEEP(23).GT.2) KEEP(23) = 1
               ENDIF
               CALL ZMUMPS_ANA_O(id%N, id%NZ, KEEP(23), id%IS1(1), id,
     &              ICNTL(1), INFO(1))
               IF (INFO(1) .LT. 0) THEN
C     -----------
C     Fatal error
C     -----------
C     Permutation was not computed; reset keep(23)
                  KEEP(23) = 0
                  GOTO 10
               END IF
            END IF
C     END OF MAX-TRANS ON THE MASTER
C
C     **********************************************************
C
C     BEGINNING OF ANALYSIS, STILL ON THE MASTER
C
C     Set up subdivisions of arrays for analysis
C
C     ------------------------------------------------------
C     Define the size of a working array 
C     that will be used as workspace ZMUMPS_ANA_F.
C     For element entry (KEEP(55).ne.0), we do not know NZ, 
C     and so the whole allocation of IW cannot be done at this 
C     point and more workspace is declared/allocated/used
C     inside ZMUMPS_ANA_F_ELT.
C     ------------------------------------------------------
C
            N8=int(id%N,8)
            IF (KEEP(55) .EQ. 0) THEN
C              ----------------
C              Assembled format
C              ----------------
               NZ8=int(id%NZ,8)
               IF ( KEEP(256) .EQ. 1 ) THEN ! KEEP(256) <-- ICNTL(7)
C                 LIW = 2 * id%NZ + 3 * id%N + 2
                  LIW8 = 2_8 * NZ8 + 3_8 * N8 + 2_8
               ELSE
C                 LIW = 2 * id%NZ + 3 * id%N + 2
                  LIW8 = 2_8 * NZ8 + 3_8 * N8 + 2_8
               ENDIF
            ELSE
C              ----------------
C              Elemental format
C              ----------------
C              Only available for AMD, METIS, and given ordering
#if defined(metis) || defined(parmetis) || defined(metis4) || defined(parmetis3)
               COND = (KEEP(60) .NE. 0) .OR. (KEEP(256) .EQ. 5)
#else
               COND = (KEEP(60) .NE. 0)
#endif
               IF( COND ) THEN
C
C
C                 we suppress supervariable detection when Schur
C                 is active or when METIS is applied
C                 Workspaces for FLAG(N), and either LEN(N) or some pointers(N+1)
                  LIW8 = N8 + N8 + 1_8
               ELSE
C                 Spaces FLAG(N), LEN(N), N+3, SVAR(0:N),
                  LIW8 =  N8 + N8 + N8 + 3_8 + N8 + 1_8
               ENDIF
C     
            ENDIF
C           We want to be sure that an array of order 3*N is 
C           available for ZMUMPS_ANA_LNEW
            IF (LIW8.LT.3_8*N8) LIW8 = 3_8*N8
            LIW=int(LIW8)
            IF (LIW8 > huge(LIW)) THEN
               INFO( 1 ) = -51
               INFO( 2 ) = LIW
               IF ( LPOK ) THEN
                  WRITE(LP, 160) 'LIW'
               END IF
               GOTO 10
            ENDIF
C           -----------------------------------------
C           work space of size LIW will be allocated
C           in ZMUMPS_ANA_F_ELT and ZMUMPS_ANA_F
C           -----------------------------------------
            IF (KEEP(23) .NE. 0) THEN
               IKEEP = id%N + 1
            ELSE
               IKEEP = 1
            END IF
            NA      = IKEEP +     id%N
            NE      = IKEEP + 2 * id%N
            FILS    = IKEEP + 3 * id%N
            FRERE   = FILS  +     id%N
            PTRAR   = FRERE +     id%N
            IF (KEEP(55) .EQ. 0) THEN
               NFSIZ  = PTRAR + 4 * id%N
               MAXIS1_CHECK = NFSIZ + id%N - 1
            ELSE
               NFSIZ  = PTRAR + 2 * (NELT + 1)
               MAXIS1_CHECK = NFSIZ + id%N -1
            ENDIF
C     
C     MAXIS1 and MAXIS1_CHECK can be different in the case where
C     max-trans has been called and resulted in the identity permutation.
C     We ensure that MAXIS1 is at least equal to MAXIS1_CHECK.
            IF ( id%MAXIS1 .LT. MAXIS1_CHECK ) THEN
                WRITE(*,*) ' WARNING: MAXIS1 <  MAXIS1_CHECK ',
     &                      'MAXIS1, MAXIS1_CHECK=',id%MAXIS1,
     &                 MAXIS1_CHECK
               CALL MUMPS_ABORT()
            END IF
C     
C     ANALYSIS PHASE
C     Some workspace by ZMUMPS_ANA_F can be reused in subsequent phases.
C       IS(IKEEP) OF LENGTH 3*N
C       IS(PTRAR) OF LENGTH 4*N (for assembled matrices
C         otherwise of size 2*NELT+1 )
C       IS(NFSIZ) OF LENGTH N holds the frontal matrix sizes
C     
            IF ( KEEP(256) .EQ. 1 ) THEN
C     Note that id%PERM_IN has been checked before.
               DO I = 1, id%N
                  id%IS1( IKEEP + I - 1 ) = id%PERM_IN( I )
               END DO
            END IF
            INFOG(1) = 0
            INFOG(2) = 0
C           Initialize structural symmetry value to not yet computed.
            INFOG(8) = -1
C            -------------------------------------------
C            We allocate something only to avoid passing
C            null pointers or arrays of size 0.
C            -------------------------------------------
            IF ( .NOT. associated( id%LISTVAR_SCHUR ) ) THEN
               SIZE_SCHUR_PASSED = 1
               LISTVAR_SCHUR_2BE_FREED=.TRUE.
               allocate( id%LISTVAR_SCHUR( 1 ), STAT=allocok )
               IF ( allocok .GT. 0 ) THEN
                  WRITE(*,*)
     &                 'PB allocating an array of size 1 in Schur '
                  CALL MUMPS_ABORT()
               END IF
            ELSE
               SIZE_SCHUR_PASSED=id%SIZE_SCHUR
               LISTVAR_SCHUR_2BE_FREED = .FALSE.
            END IF
            IF (KEEP(55) .EQ. 0) THEN
               CALL ZMUMPS_ANA_F(id%N, id%NZ, id%IRN(1), id%JCN(1),
     &              LIW, id%IS1(IKEEP),
     &              id%IS1(PTRAR), KEEP(256), id%IS1(NFSIZ),
     &              id%IS1(FILS), id%IS1(FRERE),
     &              id%LISTVAR_SCHUR(1), SIZE_SCHUR_PASSED,
     &              ICNTL(1), INFOG(1), KEEP(1),KEEP8(1),id%NSLAVES, 
     &              id%IS1(1),id)
               IF ( (KEEP(23).LE.-1).AND.(KEEP(23).GE.-6) ) THEN
C                 -- Perform max trans
                  KEEP(23) = -KEEP(23)
                  IF (.NOT. associated(id%A)) KEEP(23) = 1
                  GOTO 1234
               ENDIF
               INFOG(7)     = KEEP(256)
            ELSE
C     -------------------------------------------
C     total space required is LIW
C     -------------------------------------------
               allocate( IWtemp ( 3*id%N ), stat = IERR )
               IF ( IERR .GT. 0 ) THEN
                  INFO( 1 ) = -7
                  INFO( 2 ) = 3*id%N
                  IF ( LPOK ) THEN
                     WRITE(LP, 150) 'IWtemp'
                  END IF
                  GOTO 10
               ENDIF
               allocate( XNODEL ( id%N+1 ), stat = IERR )
               IF ( IERR .GT. 0 ) THEN
                  INFO( 1 ) = -7
                  INFO( 2 ) = id%N + 1
                  IF ( LPOK ) THEN
                     WRITE(LP, 150) 'XNODEL'
                  END IF
                  GOTO 10
               ENDIF
               IF (LELTVAR.ne.id%ELTPTR(NELT+1)-1)  THEN
C                 -- internal error 
                  INFO(1) = -2002
                  INFO(2) = id%ELTPTR(NELT+1)-1
                  GOTO 10
               ENDIF
               allocate( NODEL ( LELTVAR ), stat = IERR )
               IF ( IERR .GT. 0 ) THEN
                  INFO( 1 ) = -7
                  INFO( 2 ) = LELTVAR
                  IF ( LPOK ) THEN
                     WRITE(LP, 150) 'NODEL'
                  END IF
                  GOTO 10
               ENDIF
               CALL ZMUMPS_ANA_F_ELT(id%N, NELT,
     &              id%ELTPTR(1), id%ELTVAR(1), LIW,
     &              id%IS1(IKEEP),
     &              IWtemp(1), KEEP(256), id%IS1(NFSIZ), id%IS1(FILS),
     &              id%IS1(FRERE), id%LISTVAR_SCHUR(1),
     &              SIZE_SCHUR_PASSED, 
     &              ICNTL(1), INFOG(1), KEEP(1),KEEP8(1),
     &              id%ELTPROC(1), id%NSLAVES, 
     &              XNODEL(1), NODEL(1))
               DEALLOCATE(IWtemp)
               INFOG(7)=KEEP(256)
C     
C              XNODEL and NODEL will be deallocated later
C              since they are used in ZMUMPS_FRTELT
C     
            ENDIF
            IF ( LISTVAR_SCHUR_2BE_FREED ) THEN
C              We do not want to have LISTVAR_SCHUR
C              allocated of size 1 if Schur is off. 
               deallocate( id%LISTVAR_SCHUR )
               NULLIFY   ( id%LISTVAR_SCHUR )
            ENDIF
C           ------------------------------
C           Significant error codes should
C           always be in INFO(1/2)
C           ------------------------------
            INFO(1)=INFOG(1)
            INFO(2)=INFOG(2)
C           save statistics in KEEP array.
            KEEP(28) = INFOG(6)
C           Check error during ZMUMPS_ANA_F OR ZMUMPS_ANA_F_ELT
            IF ( INFO(1) .LT. 0 ) THEN
               GO TO 10
            ENDIF
         ENDIF
      ELSE
C     Parallel analysis
         IKEEP   = 1
         NA      = IKEEP +     id%N
         NE      = IKEEP + 2 * id%N
         FILS    = IKEEP + 3 * id%N
         FRERE   = FILS  +     id%N
         PTRAR   = FRERE +     id%N
         NFSIZ   = PTRAR + 4 * id%N
         IF(id%MYID .EQ. MASTER) THEN
            WORK1PTR => id%IS1(IKEEP : IKEEP + 3*id%N-1)
            WORK2PTR => id%IS1(PTRAR : PTRAR + 4*id%N-1)
            NFSIZPTR => id%IS1(NFSIZ : NFSIZ + id%N-1)
            FILSPTR  => id%IS1(FILS  : FILS  + id%N-1)
            FREREPTR => id%IS1(FRERE : FRERE + id%N-1)
         ELSE
C           Because our purpose is to minimize the peak memory consumption,
C           we can afford to allocate IS also on processes other than host
            ALLOCATE(WORK1PTR(3*id%N))
            ALLOCATE(WORK2PTR(4*id%N))
         END IF
         CALL ZMUMPS_ANA_F_PAR(id,
     &        WORK1PTR,
     &        WORK2PTR,
     &        NFSIZPTR,
     &        FILSPTR,
     &        FREREPTR)
         IF(id%MYID .EQ. 0) THEN
            NULLIFY(WORK1PTR, WORK2PTR, NFSIZPTR)
            NULLIFY(FILSPTR, FREREPTR)
         ELSE
            DEALLOCATE(WORK1PTR, WORK2PTR)
         END IF
         KEEP(28) = INFOG(6)
      END IF
 10   CONTINUE
      CALL MUMPS_PROPINFO( ICNTL(1), INFO(1), id%COMM, id%MYID )
      IF ( INFO(1) < 0 ) RETURN
      IF(id%MYID .EQ. MASTER) THEN
C        Save ICNTL(14) value into KEEP(12)
         CALL MUMPS_GET_PERLU(KEEP(12),ICNTL(14),
     &        KEEP(50),KEEP(54),ICNTL(6),KEEP(52))
         CALL ZMUMPS_ANA_R(id%N, id%IS1(FILS), id%IS1(FRERE),
     &        id%IS1(IKEEP+2*id%N), id%IS1(IKEEP+id%N))
C      **********************************************************
C      Continue with CALL to MAPPING routine
C        *********************
C        BEGIN SEQUENTIAL CODE
C        No mapping computed
C        *********************
C
C        In sequential, if no special root
C        reset KEEP(20) and KEEP(38) to 0
C
         IF (id%NSLAVES .EQ. 1) THEN
            id%NBSA = 0
            IF ( (id%KEEP(60).EQ.0).
     &           AND.(id%KEEP(53).EQ.0))  THEN 
C     If Schur is on (keep(60).ne.0)
C     or if RR is on (keep (53) > 0 
C     then we keep root numbers
               id%KEEP(20)=0
               id%KEEP(38)=0
            ENDIF
C     No type 2 nodes:
            id%KEEP(56)=0
C     
            id%PROCNODE = 0
C     It may also happen that KEEP(38) has already been set,
C     in the case of a distributed Schur complement (KEEP(60)=2 or 3).
C     In that case, PROCNODE should be set accordingly and KEEP(38) is
C     not modified.
            IF (id%KEEP(60) .EQ. 2 .OR. id%KEEP(60).EQ.3) THEN
               CALL ZMUMPS_SET_PROCNODE(id%KEEP(38), id%PROCNODE(1),
     &              1+2*id%NSLAVES, id%IS1(FILS),id%N)
            ENDIF
C        *******************
C        END SEQUENTIAL CODE
C        *******************
         ELSE
C        *****************************
C        BEGIN MAPPING WITH CANDIDATES
C        (NSLAVES > 1)
C        *****************************
C     
C     
C      peak is set by default to 1 largest front + One largest CB
       PEAK = dble(id%INFOG(5))*dble(id%INFOG(5)) + ! front matrix
     &        dble(id%KEEP(2))*dble(id%KEEP(2))     ! cb bloc
C     IKEEP(1:N,1) can be used as a work space since it is set
C     to its final state by the SORT_PERM subroutine below.
            SSARBR => id%IS1(IKEEP:IKEEP+id%N-1)
C     Map nodes and assign candidates for dynamic scheduling
            CALL ZMUMPS_DIST_AVOID_COPIES(id%N,id%NSLAVES,ICNTL(1),
     &           INFOG(1),
     &           id%IS1(NE),
     &           id%IS1(NFSIZ),
     &           id%IS1(FRERE),
     &           id%IS1(FILS),
     &           KEEP(1),KEEP8(1),id%PROCNODE(1),
     &           SSARBR(1),id%NBSA,PEAK,IERR
     &           )
            NULLIFY(SSARBR)
            if(IERR.eq.-999) then 
               write(6,*) ' Internal error during static mapping '
               INFO(1) = IERR
               GOTO 11
            ENDIF
            IF(IERR.NE.0) THEN 
               INFO(1) = -135
               INFO(2) = IERR
               GOTO 11
            ENDIF
C     
            CALL ZMUMPS_ANA_R(id%N, id%IS1(FILS),
     &           id%IS1(FRERE), id%IS1(IKEEP+2*id%N),
     &           id%IS1(IKEEP+id%N))
         ENDIF
 11      CONTINUE
      ENDIF
      CALL MUMPS_PROPINFO( ICNTL(1), INFO(1), id%COMM, id%MYID )
      IF ( INFO(1) < 0 ) RETURN
C     The following part is done in parallel
      CALL MPI_BCAST( id%NELT, 1, MPI_INTEGER, MASTER,
     &     id%COMM, IERR )
      IF (KEEP(55) .EQ. 0) THEN
C     Assembled matrix format. Fill up the PTRAR array
C     Broadcast id%SYM_PERM needed to fill up PTRAR
C     postpone to after computation  of id%SYM_PERM 
C     computed after id%DAD_STEPS
         if (associated(id%FRTPTR)) DEALLOCATE(id%FRTPTR)
         if (associated(id%FRTELT)) DEALLOCATE(id%FRTELT)
         allocate( id%FRTPTR(1), id%FRTELT(1) )
      ELSE
C     Element Entry: 
C     -------------------------------
C     COMPUTE THE LIST OF ELEMENTS THAT WILL BE ASSEMBLED
C     AT EACH NODE OF THE ELIMINATION TREE. ALSO COMPUTE
C     FOR EACH ELEMENT THE TREE NODE TO WHICH IT IS ASSIGNED.
C     
C     FRTPTR is an INTEGER array of length N+1 which need not be set by
C     the user. On output, FRTPTR(I) points in FRTELT to first element 
C     in the list of elements assigned to node I in the elimination tree.
C     
C     FRTELT is an INTEGER array of length NELT which need not be set by
C     the user. On output, positions FRTELT(FRTPTR(I)) to
C     FRTELT(FRTPTR(I+1)-1) contain the list of elements assigned to 
C     node I in the elimination tree.
C     
         LPTRAR = id%NELT+id%NELT+2
         CALL MUMPS_REALLOC(id%PTRAR, LPTRAR, id%INFO, LP,
     &        FORCE=.TRUE., STRING='id%PTRAR (Analysis)', ERRCODE=-7)
         CALL MUMPS_REALLOC(id%FRTPTR, id%N+1, id%INFO, LP,
     &        FORCE=.TRUE., STRING='id%FRTPTR (Analysis)', ERRCODE=-7)
         CALL MUMPS_REALLOC(id%FRTELT, id%NELT, id%INFO, LP,
     &        FORCE=.TRUE., STRING='id%FRTELT (Analysis)', ERRCODE=-7)
         CALL MUMPS_PROPINFO( ICNTL(1), INFO(1), id%COMM, id%MYID )
         IF ( INFO(1) < 0 ) RETURN
         IF(id%MYID .EQ. MASTER) THEN
C     In the elemental format case, PTRAR&friends are still
C     computed sequentially and then broadcasted
            CALL ZMUMPS_FRTELT(
     &           id%N, NELT, id%ELTPTR(NELT+1)-1, id%IS1(FRERE),
     &           id%IS1(FILS),
     &           id%IS1(IKEEP+id%N), id%IS1(IKEEP+2*id%N), XNODEL, 
     &           NODEL, id%FRTPTR(1), id%FRTELT(1), id%ELTPROC(1))
            DO I=1, id%NELT+1
               id%PTRAR(id%NELT+I+1)=id%ELTPTR(I)
            ENDDO
            deallocate(XNODEL)
            deallocate(NODEL)
         END IF
         CALL MPI_BCAST( id%PTRAR(id%NELT+2), id%NELT+1, MPI_INTEGER,
     &        MASTER, id%COMM, IERR )
         CALL MPI_BCAST( id%FRTPTR(1), id%N+1, MPI_INTEGER,
     &        MASTER, id%COMM, IERR )
         CALL MPI_BCAST( id%FRTELT(1), id%NELT,  MPI_INTEGER,
     &        MASTER, id%COMM, IERR )
      ENDIF
C     We switch again to sequential computations on the master node
      IF(id%MYID .EQ. MASTER) THEN
         IF ( INFO( 1 ) .LT. 0 ) GOTO 12
         IF ( KEEP(55) .ne. 0 ) THEN
C     ---------------------------------------
C     Build ELTPROC: correspondance between elements and slave numbers. 
C     This is used later in the initial elemental
C     matrix distribution at the beginning of the factorisation phase
C     ---------------------------------------
            CALL ZMUMPS_ELTPROC(id%N, NELT, id%ELTPROC(1),id%NSLAVES,
     &           id%PROCNODE(1))
         END IF
         NB_NIV2 = KEEP(56)
         IF ( NB_NIV2.GT.0 ) THEN
C     
            allocate(PAR2_NODES(NB_NIV2),
     &           STAT=allocok)
            IF (allocok .GT.0) then
               INFO(1)= -7
               INFO(2)= NB_NIV2
               IF ( LPOK ) THEN
                  WRITE(LP, 150) 'PAR2_NODES'
               END IF
               GOTO 12
            END IF
         ENDIF
         IF ((NB_NIV2.GT.0) .AND. (KEEP(24).EQ.0)) THEN
            INIV2 = 0
            DO 777 INODE = 1, id%N
               IF ( ( id%IS1(FRERE+INODE-1) .NE. id%N+1 ) .AND.
     &              ( MUMPS_TYPENODE(id%PROCNODE(INODE),id%NSLAVES)
     &              .eq. 2) ) THEN
                  INIV2 = INIV2 + 1
                  PAR2_NODES(INIV2) = INODE
               END IF
 777        CONTINUE
            IF ( INIV2 .NE. NB_NIV2 ) THEN
               WRITE(*,*) "Internal Error 2 in ZMUMPS_ANA_DRIVER",
     &              INIV2, NB_NIV2
               CALL MUMPS_ABORT()
            ENDIF
         ENDIF
         IF ( (KEEP(24) .NE. 0) .AND. (NB_NIV2.GT.0) ) THEN
C           allocate array to store cadidates stategy
C           for each level two nodes
            IF ( associated(id%CANDIDATES)) deallocate(id%CANDIDATES)
            allocate( id%CANDIDATES(id%NSLAVES+1,NB_NIV2),
     &           stat=allocok)
            if (allocok .gt.0) then
               INFO(1)= -7
               INFO(2)= NB_NIV2*(id%NSLAVES+1)
               IF ( LPOK ) THEN
                  WRITE(LP, 150) 'CANDIDATES'
               END IF
               GOTO 12
            END IF
            CALL MUMPS_RETURN_CANDIDATES
     &           (PAR2_NODES,id%CANDIDATES,IERR)
            IF(IERR.NE.0)  THEN
               INFO(1) = -2002
               GOTO 12
            ENDIF
C     deallocation of variables of module mumps_static_mapping
            CALL MUMPS_END_ARCH_CV()
            IF(IERR.NE.0)  THEN
               INFO(1) = -2002
               GOTO 12
            ENDIF
         ELSE
            IF (associated(id%CANDIDATES)) DEALLOCATE(id%CANDIDATES)
            allocate(id%CANDIDATES(1,1), stat=allocok)
            IF (allocok .NE. 0) THEN
               INFO(1)= -7
               INFO(2)= 1
               IF ( LPOK ) THEN
                  WRITE(LP, 150) 'CANDIDATES'
               END IF
               GOTO 12
            ENDIF
         ENDIF
C*******************************************************************
C     ---------------
 12      CONTINUE
C     ---------------
*     
*     ===============================
*     End of analysis phase on master
*     ===============================
*     
!     blocking factor for multiple RHS for ana_distm
         KEEP(84) = ICNTL(27)
      END IF
      CALL MUMPS_PROPINFO( ICNTL(1), INFO(1), id%COMM, id%MYID )
      IF ( INFO(1) < 0 ) RETURN
C     
C     We now allocate and compute arrays in NSTEPS
C     on the master, as this makes more sense.
C     
C     ==============================
C     PREPARE DATA FOR FACTORIZATION
C     ==============================
C     ------------------
      CALL MPI_BCAST( id%KEEP(1), 110, MPI_INTEGER, MASTER,
     &     id%COMM, IERR )
C     We also need to broadcast KEEP8(21) 
      CALL MUMPS_BCAST_I8( id%KEEP8(21), MASTER,
     &                     id%MYID, id%COMM, IERR)
C     --------------------------------------------------
C     Broadcast KEEP(205) which is outside the first 110
C     KEEP entries but is needed for factorization.
C     --------------------------------------------------
      CALL MPI_BCAST( id%KEEP(205), 1, MPI_INTEGER, MASTER,
     &     id%COMM, IERR )
C     --------------
C     Broadcast NBSA 
      CALL MPI_BCAST( id%NBSA, 1, MPI_INTEGER, MASTER,
     &     id%COMM, IERR )
C     -----------------
C     Global MAXFRT (computed in ZMUMPS_ANA_M)
C     is needed on all the procs during ZMUMPS_ANA_DISTM
C     to evaluate workspace for solve. 
C     We could also recompute it in ZMUMPS_ANA_DISTM
      IF (id%MYID==MASTER) KEEP(127)=INFOG(5)
      CALL MPI_BCAST( id%KEEP(127), 1, MPI_INTEGER, MASTER,
     &     id%COMM, IERR )
C     -----------------
C     Global max panel size KEEP(226)
      CALL MPI_BCAST( id%KEEP(226), 1, MPI_INTEGER, MASTER,
     &     id%COMM, IERR )
C     -----------------
C     Number of leaves not belonging to L0 KEEP(262)
C              and KEEP(263) : inner or outer sends for blocked facto
      CALL MPI_BCAST( id%KEEP(262), 2, MPI_INTEGER, MASTER,
     &     id%COMM, IERR )
C
C
C     ----------------------------------------
C     Allocate new workspace on all processors
C     ----------------------------------------
      CALL MUMPS_REALLOC(id%STEP, id%N, id%INFO, LP, FORCE=.TRUE.,
     &     STRING='id%STEP (Analysis)', ERRCODE=-7)
      IF(INFO(1).LT.0) GOTO 94
      CALL MUMPS_REALLOC(id%PROCNODE_STEPS, id%KEEP(28), id%INFO, LP,
     &     FORCE=.TRUE.,
     &     STRING='id%PROCNODE_STEPS (Analysis)', ERRCODE=-7)
      IF(INFO(1).LT.0) GOTO 94
      CALL MUMPS_REALLOC(id%NE_STEPS, id%KEEP(28), id%INFO, LP, 
     &     FORCE=.TRUE., 
     &     STRING='id%NE_STEPS (Analysis)', ERRCODE=-7)
      IF(INFO(1).LT.0) GOTO 94
      CALL MUMPS_REALLOC(id%ND_STEPS, id%KEEP(28), id%INFO, LP,
     &     FORCE=.TRUE., 
     &     STRING='id%ND_STEPS (Analysis)', ERRCODE=-7)
      IF(INFO(1).LT.0) GOTO 94
      CALL MUMPS_REALLOC(id%FRERE_STEPS, id%KEEP(28), id%INFO, LP,
     &     FORCE=.TRUE., 
     &     STRING='id%FRERE_STEPS (Analysis)', ERRCODE=-7)
      IF(INFO(1).LT.0) GOTO 94
      CALL MUMPS_REALLOC(id%DAD_STEPS, id%KEEP(28), id%INFO, LP, 
     &     FORCE=.TRUE., 
     &     STRING='id%DAD_STEPS (Analysis)', ERRCODE=-7)
      IF(INFO(1).LT.0) GOTO 94
      CALL MUMPS_REALLOC(id%FILS, id%N, id%INFO, LP, FORCE=.TRUE.,
     &     STRING='id%FILS (Analysis)', ERRCODE=-7)
      IF(INFO(1).LT.0) GOTO 94
      CALL MUMPS_REALLOC(id%SYM_PERM, id%N, id%INFO, LP, FORCE=.TRUE.,
     &     STRING='id%SYM_PERM (Analysis)', ERRCODE=-7)
      IF(INFO(1).LT.0) GOTO 94
      IF (KEEP(55) .EQ. 0) THEN
        LPTRAR = id%N+id%N
        CALL MUMPS_REALLOC(id%PTRAR, LPTRAR, id%INFO, LP, FORCE=.TRUE.,
     &       STRING='id%PTRAR (Analysis)', ERRCODE=-7)
        IF(INFO(1).LT.0) GOTO 94
      ENDIF
C     Copy data for factorization and/or solve.
      IF ( associated( id%UNS_PERM ) ) deallocate(id%UNS_PERM)
C     ================================
C     COMPUTE ON THE MASTER, BROADCAST
C     TO OTHER PROCESSES
C     ================================
      IF ( id%MYID == MASTER .AND. id%KEEP(23) .NE. 0 ) THEN
C     This one is only on the master
         allocate(id%UNS_PERM(id%N),stat=allocok)
         IF ( allocok .ne. 0) THEN
            INFO(1) = -7
            INFO(2) = id%N
            IF ( LPOK ) THEN
               WRITE(LP, 150) 'id%UNS_PERM'
            END IF
            GOTO 94
         ENDIF
C     
         DO I=1,id%N
            id%UNS_PERM(I) = id%IS1(I)
         END DO
      ENDIF
 94   CONTINUE
      CALL MUMPS_PROPINFO( ICNTL(1), INFO(1),
     &     id%COMM, id%MYID )
      IF ( id%MYID .EQ. MASTER ) THEN
         DO I=1,id%N
            id%FILS(I) = id%IS1(FILS+I-1)
         ENDDO
      END IF
      IF (id%MYID .EQ. MASTER ) THEN
C     NA -> compressed NA containing only list
C     of leaves of the elimination tree and list of roots 
C     (the two useful informations for factorization/solve).
         IF (id%N.eq.1) THEN
            NBROOT = 1
            NBLEAF = 1
         ELSE IF (id%IS1(NA+id%N-1) .LT.0) THEN
            NBLEAF = id%N
            NBROOT = id%N
         ELSE IF (id%IS1(NA+id%N-2) .LT.0) THEN
            NBLEAF = id%N-1
            NBROOT = id%IS1(NA+id%N-1)
         ELSE
            NBLEAF = id%IS1(NA+id%N-2)
            NBROOT = id%IS1(NA+id%N-1)
         ENDIF
         id%LNA = 2+NBLEAF+NBROOT
      ENDIF
      CALL MPI_BCAST( id%LNA, 1, MPI_INTEGER,
     &     MASTER, id%COMM, IERR )
      CALL MUMPS_REALLOC(id%NA, id%LNA, id%INFO, LP, FORCE=.TRUE., 
     &     STRING='id%NA (Analysis)', ERRCODE=-7)
      IF(INFO(1).LT.0) GOTO 96
      IF (id%MYID .EQ.MASTER ) THEN
C     The structure of NA is the following:
C       NA(1) is the number of leaves.
C       NA(2) is the number of roots.
C       NA(3:2+NA(1)) are the leaves.
C       NA(3+NA(1):2+NA(1)+NA(2)) are the roots.
         id%NA(1) = NBLEAF
         id%NA(2) = NBROOT
C     
C        Initialize NA with the leaves and roots
         LEAF = 3
         IF ( id%N == 1 ) THEN
            id%NA(LEAF) = 1
            LEAF = LEAF + 1
         ELSE IF (id%IS1(NA+id%N-1) < 0) THEN
            id%NA(LEAF) = - id%IS1(NA+id%N-1)-1
            LEAF = LEAF + 1
            DO I = 1, NBLEAF - 1
               id%NA(LEAF) = id%IS1(NA+I-1)
               LEAF = LEAF + 1
            ENDDO
         ELSE IF (id%IS1(NA+id%N-2) < 0 ) THEN
            INODE = - id%IS1(NA+id%N-2) - 1
            id%NA(LEAF) = INODE
            LEAF =LEAF + 1
            IF ( NBLEAF > 1 ) THEN
               DO I = 1, NBLEAF - 1
                  id%NA(LEAF) = id%IS1(NA+I-1)
                  LEAF = LEAF + 1
               ENDDO
            ENDIF
         ELSE
            DO I = 1, NBLEAF
               id%NA(LEAF) = id%IS1(NA+I-1)
               LEAF = LEAF + 1
            ENDDO
         END IF
      END IF
 96   CONTINUE
      CALL MUMPS_PROPINFO( ICNTL(1), INFO(1),
     &     id%COMM, id%MYID )
      IF ( INFO(1).LT.0 ) RETURN
      IF (associated(id%Step2node))      THEN
        DEALLOCATE(id%Step2node)
        NULLIFY(id%Step2node)
      ENDIF
      IF ( id%MYID .EQ. MASTER ) THEN
C     Build array STEP(1:id%N) to hold step numbers in
C     range 1..id%KEEP(28), allowing compression of
C     other arrays from id%N to id%KEEP(28)
C     (the number of nodes/steps in the assembly tree)
         ISTEP = 0
         DO I = 1, id%N
            IF ( id%IS1(FRERE+I-1) .ne. id%N + 1 ) THEN
C     New node in the tree.
c     (Set step( inode_n ) = inode_nsteps for principal
C     variables and -inode_nsteps for internal variables
C     of the node)
               ISTEP = ISTEP + 1
               id%STEP(I)=ISTEP
               INN = id%IS1(FILS+I-1)
               DO WHILE ( INN .GT. 0 )
                  id%STEP(INN) = - ISTEP
                  INN = id%IS1(FILS + INN -1)
               END DO
               IF (id%IS1(FRERE+I-1) .eq. 0) THEN
C     Keep root nodes list in NA
                  id%NA(LEAF) = I
                  LEAF = LEAF + 1
               ENDIF
            ENDIF
         END DO
         IF ( LEAF - 1 .NE. 2+NBROOT + NBLEAF ) THEN
            WRITE(*,*) 'Internal error 2 in ZMUMPS_ANA_DRIVER'
            CALL MUMPS_ABORT()
         ENDIF
         IF ( ISTEP .NE. id%KEEP(28) ) THEN
            write(*,*) 'Internal error 3 in ZMUMPS_ANA_DRIVER'
            CALL MUMPS_ABORT()
         ENDIF
C     ============
C     SET PROCNODE, FRERE, NE
C     ============
         DO I = 1, id%N
            IF (id%IS1(FRERE+I-1) .NE. id%N+1) THEN
               id%PROCNODE_STEPS(id%STEP(I)) = id%PROCNODE( I )
               id%FRERE_STEPS(id%STEP(I))    = id%IS1(FRERE+I-1)
               id%NE_STEPS(id%STEP(I))    = id%IS1(NE+I-1)
               id%ND_STEPS(id%STEP(I))    = id%IS1(NFSIZ+I-1)
            ENDIF
         ENDDO
C     ===============================
C     Algoritme to compute array DAD_STEPS:
C     ----
C       For each node set dad for all of its sons
C       plus, for root nodes set dad to zero.
C     
C     ===============================
         DO I = 1, id%N
C     -- skip non principal nodes
            IF ( id%STEP(I) .LE. 0) CYCLE
C     -- (I) is a principal node
            IF (id%IS1(FRERE+I-1) .eq. 0) THEN
C     -- I is a root node and has no father
               id%DAD_STEPS(id%STEP(I)) = 0
            ENDIF
C     -- Find first son node (IFS)
            IFS = id%IS1(FILS+I-1)
            DO WHILE ( IFS .GT. 0 )
               IFS= id%IS1(FILS + IFS -1)
            END DO
C     -- IFS > 0 if I is not a leave node
C     -- Go through list of brothers of IFS if any
            IFS = -IFS
            DO WHILE (IFS.GT.0) 
C     -- I is not a leave node and has a son node IFS
               id%DAD_STEPS(id%STEP(IFS)) = I
               IFS   = id%IS1(FRERE+IFS-1)
            ENDDO
         END DO
C
C     
C        Following arrays (PROCNODE and IS1) not used anymore 
C        during analysis
         deallocate(id%PROCNODE)
         NULLIFY(id%PROCNODE)
         deallocate(id%IS1)
         NULLIFY(id%IS1)
C     Reorder the tree using a variant of Liu's algorithm. Note that
C     REORDER_TREE MUST always be called since it sorts NA (the list of
C     leaves) in a valid order in the sense of a depth-first traversal.
               CALL ZMUMPS_REORDER_TREE(id%N, id%FRERE_STEPS(1),
     &              id%STEP(1),id%FILS(1), id%NA(1), id%LNA,
     &              id%NE_STEPS(1), id%ND_STEPS(1), id%DAD_STEPS(1), 
     &              id%KEEP(28), .TRUE., id%KEEP(28), id%KEEP(70),
     &              id%KEEP(50), id%INFO(1), id%ICNTL(1),id%KEEP(215),
     &              id%KEEP(234), id%KEEP(55),
     &              id%PROCNODE_STEPS(1),id%NSLAVES,PEAK,id%KEEP(90)
     &              )
            IF(id%KEEP(261).EQ.1)THEN
               CALL MUMPS_SORT_STEP(id%N, id%FRERE_STEPS(1),
     &              id%STEP(1),id%FILS(1), id%NA(1), id%LNA,
     &              id%NE_STEPS(1), id%ND_STEPS(1), id%DAD_STEPS(1), 
     &              id%KEEP(28), .TRUE., id%KEEP(28), id%INFO(1),
     &              id%ICNTL(1),id%PROCNODE_STEPS(1),id%NSLAVES
     &              )
            ENDIF
C     Compute and export some global information on the tree needed by
C     dynamic schedulers during the factorization. The type of
C     information depends on the selected strategy.
         IF ((id%KEEP(76).GE.4).OR.(id%KEEP(76).GE.6).OR.
     &              (id%KEEP(47).EQ.4).OR.((id%KEEP(81).GT.0)
     &              .AND.(id%KEEP(47).GE.2)))THEN
            IS_BUILD_LOAD_MEM_CALLED=.TRUE.
            IF ((id%KEEP(47) .EQ. 4).OR.
     &           (( id%KEEP(81) .GT. 0).AND.(id%KEEP(47).GE.2))) THEN
               IF(id%NSLAVES.GT.1) THEN
C                 NBSA is the total number of subtrees  and
C                 is an upperbound of the local number of
C                 subtrees
                  SIZE_TEMP_MEM = id%NBSA
               ELSE
C                 Only one processor, NA(2) is the number of leaves
                  SIZE_TEMP_MEM = id%NA(2)
               ENDIF
            ELSE
               SIZE_TEMP_MEM = 1
            ENDIF
            IF((id%KEEP(76).EQ.4).OR.(id%KEEP(76).EQ.6))THEN
               SIZE_DEPTH_FIRST=id%KEEP(28)
            ELSE
               SIZE_DEPTH_FIRST=1
            ENDIF
            allocate(TEMP_MEM(SIZE_TEMP_MEM,id%NSLAVES),STAT=allocok) 
            IF (allocok .NE.0) THEN
               INFO(1)= -7
               INFO(2)= SIZE_TEMP_MEM*id%NSLAVES
               IF ( LPOK ) THEN
                  WRITE(LP, 150) 'TEMP_MEM'
               END IF
               GOTO 80
            END IF
            allocate(TEMP_LEAF(SIZE_TEMP_MEM,id%NSLAVES),
     &           stat=allocok) 
            IF (allocok .ne.0) then
               IF ( LPOK ) THEN
                  WRITE(LP, 150) 'TEMP_LEAF'
               END IF
               INFO(1)= -7
               INFO(2)= SIZE_TEMP_MEM*id%NSLAVES
               GOTO 80
            end if
            allocate(TEMP_SIZE(SIZE_TEMP_MEM,id%NSLAVES),
     &           stat=allocok) 
            IF (allocok .ne.0) then
               IF ( LPOK ) THEN
                  WRITE(LP, 150) 'TEMP_SIZE'
               END IF
               INFO(1)= -7
               INFO(2)= SIZE_TEMP_MEM*id%NSLAVES
               GOTO 80
            end if
            allocate(TEMP_ROOT(SIZE_TEMP_MEM,id%NSLAVES),
     &           stat=allocok) 
            IF (allocok .ne.0) then
               IF ( LPOK ) THEN
                  WRITE(LP, 150) 'TEMP_ROOT'
               END IF
               INFO(1)= -7
               INFO(2)= SIZE_TEMP_MEM*id%NSLAVES
               GOTO 80
            end if
            allocate(DEPTH_FIRST(SIZE_DEPTH_FIRST),stat=allocok) 
            IF (allocok .ne.0) then
               IF ( LPOK ) THEN
                  WRITE(LP, 150) 'DEPTH_FIRST'
               END IF
               INFO(1)= -7
               INFO(2)= SIZE_DEPTH_FIRST
               GOTO 80
            end if
            ALLOCATE(DEPTH_FIRST_SEQ(SIZE_DEPTH_FIRST),stat=allocok) 
            IF (allocok .ne.0) then
               IF ( LPOK ) THEN
                  WRITE(LP, 150) 'DEPTH_FIRST_SEQ'
               END IF
               INFO(1)= -7
               INFO(2)= SIZE_DEPTH_FIRST
               GOTO 80
            end if
            ALLOCATE(SBTR_ID(SIZE_DEPTH_FIRST),stat=allocok) 
            IF (allocok .ne.0) then
               IF ( LPOK ) THEN
                  WRITE(LP, 150) 'SBTR_ID'
               END IF
               INFO(1)= -7
               INFO(2)= SIZE_DEPTH_FIRST
               GOTO 80
            end if
            IF(id%KEEP(76).EQ.5)THEN
C     We reuse the same variable as before
               SIZE_COST_TRAV=id%KEEP(28)
            ELSE
               SIZE_COST_TRAV=1
            ENDIF
            allocate(COST_TRAV_TMP(SIZE_COST_TRAV),stat=allocok) 
            IF (allocok .ne.0) then
               IF ( LPOK ) THEN
                  WRITE(LP, 150) 'COST_TRAV_TMP'
               END IF
               INFO(1)= -7
               INFO(2)= SIZE_COST_TRAV
               GOTO 80
            END IF
            IF(id%KEEP(76).EQ.5)THEN
               IF(id%KEEP(70).EQ.0)THEN
                  id%KEEP(70)=5
               ENDIF
               IF(id%KEEP(70).EQ.1)THEN
                  id%KEEP(70)=6
               ENDIF
            ENDIF
            IF(id%KEEP(76).EQ.4)THEN
               IF(id%KEEP(70).EQ.0)THEN
                  id%KEEP(70)=3
               ENDIF
               IF(id%KEEP(70).EQ.1)THEN
                  id%KEEP(70)=4
               ENDIF
            ENDIF
            CALL ZMUMPS_BUILD_LOAD_MEM_INFO(id%N, id%FRERE_STEPS(1),
     &           id%STEP(1),id%FILS(1), id%NA(1), id%LNA,
     &           id%NE_STEPS(1), id%ND_STEPS(1), id%DAD_STEPS(1), 
     &           id%KEEP(28), .TRUE., id%KEEP(28), id%KEEP(70),
     &           id%KEEP(50), id%INFO(1), id%ICNTL(1),id%KEEP(47),
     &           id%KEEP(81),id%KEEP(76),id%KEEP(215),
     &           id%KEEP(234), id%KEEP(55),
     &           id%PROCNODE_STEPS(1),TEMP_MEM,id%NSLAVES, 
     &           SIZE_TEMP_MEM, PEAK,id%KEEP(90),SIZE_DEPTH_FIRST,
     &           SIZE_COST_TRAV,DEPTH_FIRST(1),DEPTH_FIRST_SEQ(1),
     &           COST_TRAV_TMP(1),
     &           TEMP_LEAF,TEMP_SIZE,TEMP_ROOT,SBTR_ID(1)
     &              )
         END IF
         CALL ZMUMPS_SORT_PERM(id%N, id%NA(1), id%LNA,
     &        id%NE_STEPS(1), id%SYM_PERM(1),
     &        id%FILS(1), id%DAD_STEPS(1),
     &        id%STEP(1), id%KEEP(28), id%INFO(1) )
      ENDIF
 80   CONTINUE
C     Broadcast errors
      CALL MUMPS_PROPINFO( ICNTL(1), INFO(1),
     &     id%COMM, id%MYID )
      IF ( INFO(1).LT.0 ) RETURN
C     ---------------------------------------------------
C     Broadcast information computed on the master to
C     the slaves.
C     The matrix itself with numerical values and
C     integer data for the arrowhead/element description 
C     will be received at the beginning of FACTO.
C     ---------------------------------------------------
      CALL MPI_BCAST( id%FILS(1), id%N, MPI_INTEGER,
     &     MASTER, id%COMM, IERR )
      CALL MPI_BCAST( id%NA(1), id%LNA, MPI_INTEGER,
     &     MASTER, id%COMM, IERR )
      CALL MPI_BCAST( id%STEP(1), id%N, MPI_INTEGER,
     &     MASTER, id%COMM, IERR )
      CALL MPI_BCAST( id%PROCNODE_STEPS(1), id%KEEP(28), MPI_INTEGER,
     &     MASTER, id%COMM, IERR )
      CALL MPI_BCAST( id%DAD_STEPS(1), id%KEEP(28), MPI_INTEGER,
     &     MASTER, id%COMM, IERR )
      CALL MPI_BCAST( id%FRERE_STEPS(1), id%KEEP(28), MPI_INTEGER,
     &     MASTER, id%COMM, IERR)
      CALL MPI_BCAST( id%NE_STEPS(1), id%KEEP(28), MPI_INTEGER,
     &     MASTER, id%COMM, IERR )
      CALL MPI_BCAST( id%ND_STEPS(1), id%KEEP(28), MPI_INTEGER,
     &     MASTER, id%COMM, IERR )
      CALL MPI_BCAST( id%SYM_PERM(1), id%N, MPI_INTEGER,
     &     MASTER, id%COMM, IERR )
      IF (KEEP(55) .EQ. 0) THEN
C     Assembled matrix format. Fill up the PTRAR array
C     Broadcast id%SYM_PERM needed to fill up PTRAR
C     At the end of ANA_N_PAR, PTRAR is already on every processor
C     because it is computed in a distributed way.
C     No need to broadcast it again
         CALL ZMUMPS_ANA_N_PAR(id, id%PTRAR(1))
         IF(id%MYID .EQ. MASTER) THEN
C           -----------------------------------
C           For distributed structure on entry,
C           we can now deallocate the complete
C           structure IRN / JCN.
C           -----------------------------------
            IF ( (KEEP(244) .EQ. 1) .AND. (KEEP(54) .EQ. 3) ) THEN
               DEALLOCATE( id%IRN )
               DEALLOCATE( id%JCN )
            END IF
         END IF
      ENDIF
C     
C     Store size of the stack memory for each
C     of the sequential subtree.
      IF((id%KEEP(76).EQ.4).OR.(id%KEEP(76).EQ.6))THEN
         IF(associated(id%DEPTH_FIRST))
     &        deallocate(id%DEPTH_FIRST)
         allocate(id%DEPTH_FIRST(id%KEEP(28)),stat=allocok)
         IF (allocok .ne.0) then
            INFO(1)= -7
            INFO(2)= id%KEEP(28)
            IF ( LPOK ) THEN
               WRITE(LP, 150) 'id%DEPTH_FIRST'
            END IF
            GOTO 87
         END IF
         IF(associated(id%DEPTH_FIRST_SEQ))
     *        DEALLOCATE(id%DEPTH_FIRST_SEQ)
         ALLOCATE(id%DEPTH_FIRST_SEQ(id%KEEP(28)),stat=allocok)
         IF (allocok .ne.0) then
            INFO(1)= -7
            INFO(2)= id%KEEP(28)
            IF ( LPOK ) THEN
               WRITE(LP, 150) 'id%DEPTH_FIRST_SEQ'
            END IF
            GOTO 87
         END IF
         IF(associated(id%SBTR_ID))
     *        DEALLOCATE(id%SBTR_ID)
         ALLOCATE(id%SBTR_ID(id%KEEP(28)),stat=allocok)
         IF (allocok .ne.0) then
            INFO(1)= -7
            INFO(2)= id%KEEP(28)
            IF ( LPOK ) THEN
               WRITE(LP, 150) 'id%DEPTH_FIRST_SEQ'
            END IF
            GOTO 87
         END IF
         IF(id%MYID.EQ.MASTER)THEN
            id%DEPTH_FIRST(1:id%KEEP(28))=DEPTH_FIRST(1:id%KEEP(28))
            id%DEPTH_FIRST_SEQ(1:id%KEEP(28))=
     &           DEPTH_FIRST_SEQ(1:id%KEEP(28))
            id%SBTR_ID(1:KEEP(28))=SBTR_ID(1:KEEP(28))
         ENDIF
         CALL MPI_BCAST( id%DEPTH_FIRST(1), id%KEEP(28), MPI_INTEGER,
     &           MASTER, id%COMM, IERR )         
         CALL MPI_BCAST( id%DEPTH_FIRST_SEQ(1), id%KEEP(28),
     &           MPI_INTEGER,MASTER, id%COMM, IERR )  
         CALL MPI_BCAST( id%SBTR_ID(1), id%KEEP(28),
     &           MPI_INTEGER,MASTER, id%COMM, IERR )  
      ELSE
         IF(associated(id%DEPTH_FIRST))
     &        deallocate(id%DEPTH_FIRST)
         allocate(id%DEPTH_FIRST(1),stat=allocok)
         IF (allocok .ne.0) then
            INFO(1)= -7
            INFO(2)= 1
            IF ( LPOK ) THEN
               WRITE(LP, 150) 'id%DEPTH_FIRST'
            END IF
            GOTO 87
         END IF
         IF(associated(id%DEPTH_FIRST_SEQ))
     *        DEALLOCATE(id%DEPTH_FIRST_SEQ)
         ALLOCATE(id%DEPTH_FIRST_SEQ(1),stat=allocok)
         IF (allocok .ne.0) then
            INFO(1)= -7
            INFO(2)= 1
            IF ( LPOK ) THEN
               WRITE(LP, 150) 'id%DEPTH_FIRST_SEQ'
            END IF
            GOTO 87
         END IF
         IF(associated(id%SBTR_ID))
     *        DEALLOCATE(id%SBTR_ID)
         ALLOCATE(id%SBTR_ID(1),stat=allocok)
         IF (allocok .ne.0) then
            INFO(1)= -7
            INFO(2)= 1
            IF ( LPOK ) THEN
               WRITE(LP, 150) 'id%DEPTH_FIRST_SEQ'
            END IF
            GOTO 87
         END IF
         id%SBTR_ID(1)=0
         id%DEPTH_FIRST(1)=0
         id%DEPTH_FIRST_SEQ(1)=0
      ENDIF
      IF(id%KEEP(76).EQ.5)THEN
         IF(associated(id%COST_TRAV))
     &        deallocate(id%COST_TRAV)
         allocate(id%COST_TRAV(id%KEEP(28)),stat=allocok)
         IF (allocok .ne.0) then
            IF ( LPOK ) THEN
               WRITE(LP, 150) 'id%COST_TRAV'
            END IF
            INFO(1)= -7
            INFO(2)= id%KEEP(28)
            GOTO 87
         END IF
         IF(id%MYID.EQ.MASTER)THEN
            id%COST_TRAV(1:id%KEEP(28))=
     &      dble(COST_TRAV_TMP(1:id%KEEP(28)))
         ENDIF
         CALL MPI_BCAST( id%COST_TRAV(1), id%KEEP(28),
     &        MPI_DOUBLE_PRECISION,MASTER, id%COMM, IERR )         
      ELSE
         IF(associated(id%COST_TRAV))
     &        deallocate(id%COST_TRAV)
         allocate(id%COST_TRAV(1),stat=allocok)
         IF (allocok .ne.0) then
            IF ( LPOK ) THEN
               WRITE(LP, 150) 'id%COST_TRAV(1)'
            END IF
            INFO(1)= -7
            INFO(2)= 1
            GOTO 87
         END IF
         id%COST_TRAV(1)=0.0d0
      ENDIF
      IF (id%KEEP(47) .EQ. 4 .OR.
     &     ((id%KEEP(81) .GT. 0).AND.(id%KEEP(47).GE.2))) THEN
         IF(id%MYID .EQ. MASTER)THEN
            DO K=1,id%NSLAVES
               DO J=1,SIZE_TEMP_MEM
                  IF(TEMP_MEM(J,K) < 0.0D0) GOTO 666 
               ENDDO
 666           CONTINUE
               J=J-1
               IF (id%KEEP(46) == 1) THEN
                  IDEST = K - 1
               ELSE
                  IDEST = K
               ENDIF
               IF (IDEST .NE. MASTER) THEN
                  CALL MPI_SEND(J,1,MPI_INTEGER,IDEST,0,
     &                 id%COMM,IERR)
                  CALL MPI_SEND(TEMP_MEM(1,K),J,MPI_DOUBLE_PRECISION,
     &                 IDEST, 0, id%COMM,IERR)
                  CALL MPI_SEND(TEMP_LEAF(1,K),J,MPI_INTEGER,
     &                 IDEST, 0, id%COMM,IERR)
                  CALL MPI_SEND(TEMP_SIZE(1,K),J,MPI_INTEGER,
     &                 IDEST, 0, id%COMM,IERR)
                  CALL MPI_SEND(TEMP_ROOT(1,K),J,MPI_INTEGER,
     &                 IDEST, 0, id%COMM,IERR)             
               ELSE
                  IF(associated(id%MEM_SUBTREE))
     &                 deallocate(id%MEM_SUBTREE)
                  allocate(id%MEM_SUBTREE(J),stat=allocok)
                  IF (allocok .ne.0) then
                     IF ( LPOK ) THEN
                        WRITE(LP, 150) 'id%MEM_SUBTREE'
                     END IF
                     INFO(1)= -7
                     INFO(2)= J
                     GOTO 87
                  END IF
                  id%NBSA_LOCAL = J
                  id%MEM_SUBTREE(1:J)=TEMP_MEM(1:J,1)
                  IF(associated(id%MY_ROOT_SBTR))
     &                 deallocate(id%MY_ROOT_SBTR)
                  allocate(id%MY_ROOT_SBTR(J),stat=allocok)
                  IF (allocok .ne.0) then
                     IF ( LPOK ) THEN
                        WRITE(LP, 150) 'id%MY_ROOT_SBTR'
                     END IF
                     INFO(1)= -7
                     INFO(2)= J
                     GOTO 87
                  END IF
                  id%MY_ROOT_SBTR(1:J)=TEMP_ROOT(1:J,1)
                  IF(associated(id%MY_FIRST_LEAF))
     &                 deallocate(id%MY_FIRST_LEAF)
                  allocate(id%MY_FIRST_LEAF(J),stat=allocok)
                  IF (allocok .ne.0) then
                     IF ( LPOK ) THEN
                        WRITE(LP, 150) 'id%MY_FIRST_LEAF'
                     END IF
                     INFO(1)= -7
                     INFO(2)= J
                     GOTO 87
                  END IF
                  id%MY_FIRST_LEAF(1:J)=TEMP_LEAF(1:J,1)
                  IF(associated(id%MY_NB_LEAF))
     &                 deallocate(id%MY_NB_LEAF)
                  allocate(id%MY_NB_LEAF(J),stat=allocok)
                  IF (allocok .ne.0) then
                     IF ( LPOK ) THEN
                        WRITE(LP, 150) 'id%MY_NB_LEAF'
                     END IF
                     INFO(1)= -7
                     INFO(2)= J
                     GOTO 87
                  END IF
                  id%MY_NB_LEAF(1:J)=TEMP_SIZE(1:J,1)
               ENDIF
            ENDDO
         ELSE
            CALL MPI_RECV(id%NBSA_LOCAL,1,MPI_INTEGER,
     &           MASTER,0,id%COMM,STATUS, IERR)
            IF(associated(id%MEM_SUBTREE))
     &           deallocate(id%MEM_SUBTREE)
            allocate(id%MEM_SUBTREE(id%NBSA_LOCAL),stat=allocok)
            IF (allocok .ne.0) then
               IF ( LPOK ) THEN
                  WRITE(LP, 150) 'id%MEM_SUBTREE'
               END IF
               INFO(1)= -7
               INFO(2)= id%NBSA_LOCAL
               GOTO 87
            END IF
            IF(associated(id%MY_ROOT_SBTR))
     &           deallocate(id%MY_ROOT_SBTR)
            allocate(id%MY_ROOT_SBTR(id%NBSA_LOCAL),stat=allocok)
            IF (allocok .ne.0) then
               IF ( LPOK ) THEN
                  WRITE(LP, 150) 'id%MY_ROOT_SBTR'
               END IF
               INFO(1)= -7
               INFO(2)= id%NBSA_LOCAL
               GOTO 87
            END IF
            IF(associated(id%MY_FIRST_LEAF))
     &           deallocate(id%MY_FIRST_LEAF)
            allocate(id%MY_FIRST_LEAF(id%NBSA_LOCAL),stat=allocok)
            IF (allocok .ne.0) then
               IF ( LPOK ) THEN
                  WRITE(LP, 150) 'MY_FIRST_LEAF'
               END IF
               INFO(1)= -7
               INFO(2)= id%NBSA_LOCAL
               GOTO 87
            END IF
            IF(associated(id%MY_NB_LEAF))
     &           deallocate(id%MY_NB_LEAF)
            allocate(id%MY_NB_LEAF(id%NBSA_LOCAL),stat=allocok)
            IF (allocok .ne.0) then
               IF ( LPOK ) THEN
                  WRITE(LP, 150) 'MY_NB_LEAF'
               END IF
               INFO(1)= -7
               INFO(2)= id%NBSA_LOCAL
               GOTO 87
            END IF
            CALL MPI_RECV(id%MEM_SUBTREE(1),id%NBSA_LOCAL,
     &           MPI_DOUBLE_PRECISION,MASTER,0,
     &           id%COMM,STATUS,IERR)
            CALL MPI_RECV(id%MY_FIRST_LEAF(1),id%NBSA_LOCAL,
     &           MPI_INTEGER,MASTER,0,
     &           id%COMM,STATUS,IERR)
            CALL MPI_RECV(id%MY_NB_LEAF(1),id%NBSA_LOCAL,
     &           MPI_INTEGER,MASTER,0,
     &           id%COMM,STATUS,IERR)
            CALL MPI_RECV(id%MY_ROOT_SBTR(1),id%NBSA_LOCAL,
     &           MPI_INTEGER,MASTER,0,
     &           id%COMM,STATUS,IERR)
         ENDIF
      ELSE
         id%NBSA_LOCAL = -999999
         IF(associated(id%MEM_SUBTREE))
     &        deallocate(id%MEM_SUBTREE)
         allocate(id%MEM_SUBTREE(1),stat=allocok)
         IF (allocok .ne.0) then
            IF ( LPOK ) THEN
               WRITE(LP, 150) 'id%MEM_SUBTREE(1)'
            END IF
            INFO(1)= -7
            INFO(2)= 1
            GOTO 87
         END IF
         IF(associated(id%MY_ROOT_SBTR))
     &        deallocate(id%MY_ROOT_SBTR)
         allocate(id%MY_ROOT_SBTR(1),stat=allocok)
         IF (allocok .ne.0) then
            IF ( LPOK ) THEN
               WRITE(LP, 150) 'id%MY_ROOT_SBTR(1)'
            END IF
            INFO(1)= -7
            INFO(2)= 1
            GOTO 87
         END IF
         IF(associated(id%MY_FIRST_LEAF))
     &        deallocate(id%MY_FIRST_LEAF)
         allocate(id%MY_FIRST_LEAF(1),stat=allocok)
         IF (allocok .ne.0) then
            IF ( LPOK ) THEN
               WRITE(LP, 150) 'id%MY_FIRST_LEAF(1)'
            END IF
            INFO(1)= -7
            INFO(2)= 1
            GOTO 87
         END IF
         IF(associated(id%MY_NB_LEAF))
     &        deallocate(id%MY_NB_LEAF)
         allocate(id%MY_NB_LEAF(1),stat=allocok)
         IF (allocok .ne.0) then
            IF ( LPOK ) THEN
               WRITE(LP, 150) 'id%MY_NB_LEAF(1)'
            END IF
            INFO(1)= -7
            INFO(2)= 1
            GOTO 87
         END IF
      ENDIF
      IF(id%MYID.EQ.MASTER)THEN
         IF(IS_BUILD_LOAD_MEM_CALLED)THEN 
            deallocate(TEMP_MEM)
            deallocate(TEMP_SIZE)
            deallocate(TEMP_ROOT)
            deallocate(TEMP_LEAF)
            deallocate(COST_TRAV_TMP)
            deallocate(DEPTH_FIRST)
            deallocate(DEPTH_FIRST_SEQ)
            deallocate(SBTR_ID)
         ENDIF
      ENDIF
 87   CONTINUE
      CALL MUMPS_PROPINFO( ICNTL(1), INFO(1),
     &     id%COMM, id%MYID )
      IF ( INFO(1).LT.0 ) RETURN
C     
      NB_NIV2 = KEEP(56)        ! KEEP(1:110) was broadcast earlier
C     NB_NIV2 is now available on all processors.
      IF (  NB_NIV2.GT.0  ) THEN
C        Allocate arrays on slaves
         if (id%MYID.ne.MASTER) then
            IF (associated(id%CANDIDATES)) deallocate(id%CANDIDATES)
            allocate(PAR2_NODES(NB_NIV2),
     &           id%CANDIDATES(id%NSLAVES+1,NB_NIV2),
     &           STAT=allocok)
            IF (allocok .ne.0) then
               INFO(1)= -7
               INFO(2)= NB_NIV2*(id%NSLAVES+1)
               IF ( LPOK ) THEN
                  WRITE(LP, 150) 'PAR2_NODES/id%CANDIDATES'
               END IF
            end if
         end if
         CALL MUMPS_PROPINFO( ICNTL(1), INFO(1),
     &        id%COMM, id%MYID )
         IF ( INFO(1).LT.0 ) RETURN
         CALL MPI_BCAST(PAR2_NODES(1),NB_NIV2,
     &        MPI_INTEGER, MASTER, id%COMM, IERR )
         IF (KEEP(24) .NE.0 ) THEN
            CALL MPI_BCAST(id%CANDIDATES(1,1),
     &           (NB_NIV2*(id%NSLAVES+1)),
     &           MPI_INTEGER, MASTER, id%COMM, IERR )
         ENDIF
      ENDIF
      IF ( associated(id%ISTEP_TO_INIV2)) THEN
         deallocate(id%ISTEP_TO_INIV2)
         NULLIFY(id%ISTEP_TO_INIV2)
      ENDIF
      IF ( associated(id%I_AM_CAND)) THEN
         deallocate(id%I_AM_CAND)
         NULLIFY(id%I_AM_CAND)
      ENDIF
      IF (NB_NIV2.EQ.0) THEN 
C     allocate dummy arrays
C     ISTEP_TO_INIV2 will never be used 
C     Add a parameter SIZE_ISTEP_TO_INIV2 and make
C     it always available in a keep(71)
         id%KEEP(71) = 1
      ELSE
         id%KEEP(71) = id%KEEP(28)
      ENDIF
      allocate(id%ISTEP_TO_INIV2(id%KEEP(71)),
     &     id%I_AM_CAND(max(NB_NIV2,1)),
     &     stat=allocok)
      IF (allocok .gt.0) THEN
         IF ( LPOK ) THEN
            WRITE(LP, 150) 'id%ISTEP_TO_INIV2'
            WRITE(LP, 150) 'id%TAB_POS_IN_PERE'
         END IF
         INFO(1)= -7
         IF (NB_NIV2.EQ.0) THEN
            INFO(2)= 2
         ELSE
            INFO(2)= id%KEEP(28)+NB_NIV2*(id%NSLAVES+2)
         END IF
         GOTO 321
      ENDIF
      IF ( NB_NIV2 .GT.0 ) THEN
         DO INIV2 = 1, NB_NIV2
            INN = PAR2_NODES(INIV2)
            id%ISTEP_TO_INIV2(abs(id%STEP(INN))) = INIV2
         END DO 
         CALL ZMUMPS_BUILD_I_AM_CAND( id%NSLAVES, KEEP(79),
     &        NB_NIV2, id%MYID_NODES,
     &        id%CANDIDATES(1,1), id%I_AM_CAND(1) )
      ENDIF
#if ! defined(OLD_LOAD_MECHANISM)
      IF (associated(id%FUTURE_NIV2)) THEN
         deallocate(id%FUTURE_NIV2)
         NULLIFY(id%FUTURE_NIV2)
      ENDIF
      allocate(id%FUTURE_NIV2(id%NSLAVES), stat=allocok)
      IF (allocok .gt.0) THEN
         IF ( LPOK ) THEN
            WRITE(LP, 150) 'FUTURE_NIV2'
         END IF
         INFO(1)= -7
         INFO(2)= id%NSLAVES
         GOTO 321
      ENDIF
      id%FUTURE_NIV2=0
      DO INIV2 = 1, NB_NIV2
         IDEST = MUMPS_PROCNODE(
     &        id%PROCNODE_STEPS(abs(id%STEP(PAR2_NODES(INIV2)))),
     &        id%NSLAVES)
         id%FUTURE_NIV2(IDEST+1)=id%FUTURE_NIV2(IDEST+1)+1
      ENDDO
#endif
      IF ( I_AM_SLAVE ) THEN
C     Allocate id%TAB_POS_IN_PERE, 
C     TAB_POS_IN_PERE is an array of size (id%NSLAVES+2,NB_NIV2)
C     where NB_NIV2 is the number of type 2 nodes in the tree.
         IF ( associated(id%TAB_POS_IN_PERE)) THEN
            deallocate(id%TAB_POS_IN_PERE)
            NULLIFY(id%TAB_POS_IN_PERE)
         ENDIF
         allocate(id%TAB_POS_IN_PERE(id%NSLAVES+2,max(NB_NIV2,1)),
     &        stat=allocok)
         IF (allocok .gt.0) THEN
            IF ( LPOK ) THEN
               WRITE(LP, 150) 'id%ISTEP_TO_INIV2'
               WRITE(LP, 150) 'id%TAB_POS_IN_PERE'
            END IF
            INFO(1)= -7
            IF (NB_NIV2.EQ.0) THEN
               INFO(2)= 2
            ELSE
               INFO(2)= id%KEEP(28)+NB_NIV2*(id%NSLAVES+2)
            END IF
            GOTO 321
         ENDIF
      END IF
C     deallocate PAR2_NODES  that was computed
C     on master and broadcasted on all slaves
      IF (NB_NIV2.GT.0) deallocate (PAR2_NODES)
 321  CONTINUE
C     ----------------
C     Check for errors
C     ----------------
      CALL MUMPS_PROPINFO( ICNTL(1), INFO(1),
     &     id%COMM, id%MYID )
      IF ( INFO(1).LT.0 ) RETURN
C     ------------------------------
C     Perform again the subdivision of array
C     IS1, both on the master and on
C     the slaves. This is done so to
C     ease the passage to the model
C     where master will work.
C     ------------------------------
C     
      IF ( KEEP(23).NE.0 .and. id%MYID .EQ. MASTER ) THEN
         IKEEP = id%N + 1
      ELSE
         IKEEP = 1
      END IF
      FILS   = IKEEP + 3 * id%N
      NE     = IKEEP + 2 * id%N
      NA     = IKEEP +     id%N
      FRERE  = FILS  + id%N
      PTRAR  = FRERE + id%N
      IF (KEEP(55) .EQ. 0) THEN
         IF ( id%MYID.EQ.MASTER ) THEN
C     The master has PTRAR with length 4 * N
            NFSIZ   = PTRAR  + 4 * id%N
         ELSE
C     The slave only has PTRAR with length 2 * N
            NFSIZ   = PTRAR  + 2 * id%N
         ENDIF
      ELSE
C     For element entry, PTRAR has length 2 * (NELT + 1)
         NFSIZ   = PTRAR  + 2 * (NELT + 1)
      END IF
      IF ( KEEP(38) .NE. 0 ) THEN
C     -------------------------
C     Initialize root structure
C     -------------------------
         CALL ZMUMPS_INIT_ROOT_ANA( id%MYID,
     &        id%NSLAVES, id%N, id%root,
     &        id%COMM_NODES, KEEP( 38 ), id%FILS(1),
     &        id%KEEP(50), id%KEEP(46),
     &        id%KEEP(51)
     &        , id%KEEP(60), id%NPROW, id%NPCOL, id%MBLOCK, id%NBLOCK
     &        )
      ELSE
         id%root%yes = .FALSE.
      END IF
      IF ( KEEP(38) .NE. 0 .and. I_AM_SLAVE ) THEN
C     -----------------------------------------------
C     Check if at least one processor belongs to the
C     root. In the case where all of them have MYROW
C     equal to -1, this could be a problem due to the
C     BLACS. (mpxlf90_r and IBM BLACS).
C     -----------------------------------------------
         CALL MPI_ALLREDUCE(id%root%MYROW, MYROW_CHECK, 1,
     &        MPI_INTEGER, MPI_MAX, id%COMM_NODES, IERR)
         IF ( MYROW_CHECK .eq. -1) THEN
            INFO(1) = -25
            INFO(2) = 0
         END IF
         IF ( id%root%MYROW .LT. -1 .OR.
     &        id%root%MYCOL .LT. -1 ) THEN
            INFO(1) = -25
            INFO(2) = 0
         END IF
         IF ( LPOK .AND. INFO(1) == -25 ) THEN
            WRITE(LP, '(A)')
     &           'Problem with your version of the BLACS.'
            WRITE(LP, '(A)') 'Try using a BLACS version from netlib.'
         ENDIF
      END IF
C     ----------------
C     Check for errors
C     ----------------
      CALL MUMPS_PROPINFO( ICNTL(1), INFO(1),
     &     id%COMM, id%MYID )
      IF ( INFO(1).LT.0 ) RETURN
      IF ( I_AM_SLAVE ) THEN
C     
C     
         IF (KEEP(55) .EQ. 0) THEN
            CALL ZMUMPS_ANA_DIST_ARROWHEADS( id%MYID,
     &           id%NSLAVES, id%N, id%PROCNODE_STEPS(1),
     &           id%STEP(1), id%PTRAR(1),
     &           id%PTRAR(id%N +1),
     &           id%ISTEP_TO_INIV2(1), id%I_AM_CAND(1),
     &           KEEP(1),KEEP8(1), ICNTL(1), id )
         ELSE
            CALL ZMUMPS_ANA_DIST_ELEMENTS( id%MYID,
     &           id%NSLAVES, id%N, id%PROCNODE_STEPS(1),
     &           id%STEP(1),
     &           id%PTRAR(1),
     &           id%PTRAR(id%NELT+2 ),
     &           id%NELT, 
     &           id%FRTPTR(1), id%FRTELT(1),
     &           KEEP(1), KEEP8(1), ICNTL(1), id%KEEP(50) )
         ENDIF
      ENDIF
C     -----------------------------------------
C     Perform some local analysis on the slaves
C     to estimate the size of the working space
C     for factorization
C     -----------------------------------------
      IF ( I_AM_SLAVE ) THEN
         locI_AM_CAND => id%I_AM_CAND
         locMYID_NODES = id%MYID_NODES
         locMYID       = id%MYID
C     
C     Precompute estimates of local_m,local_n
C     (number of rows/columns mapped on each processor)
C     in case of parallel root node.
C     
            IF ( id%root%yes ) THEN
               LOCAL_M = numroc( id%ND_STEPS(id%STEP(KEEP(38))),
     &              id%root%MBLOCK, id%root%MYROW, 0,
     &              id%root%NPROW )
               LOCAL_M = max(1, LOCAL_M)
               LOCAL_N = numroc( id%ND_STEPS(id%STEP(KEEP(38))),
     &              id%root%NBLOCK, id%root%MYCOL, 0,
     &              id%root%NPCOL )
            ELSE
               LOCAL_M = 0
               LOCAL_N = 0
            END IF
            IF  ( KEEP(60) .EQ. 2 .OR. KEEP(60) .EQ. 3 ) THEN
C     Return minimum nb rows/cols to user
               id%SCHUR_MLOC=LOCAL_M
               id%SCHUR_NLOC=LOCAL_N
C     Also store them in root structure for convenience
               id%root%SCHUR_MLOC=LOCAL_M
               id%root%SCHUR_NLOC=LOCAL_N
            ENDIF
               IF ( .NOT. associated(id%CANDIDATES)) THEN
                  ALLOCATE(id%CANDIDATES(id%NSLAVES+1,1))
               ENDIF
               CALL ZMUMPS_ANA_DISTM( locMYID_NODES, id%N,
     &              id%STEP(1), id%FRERE_STEPS(1), id%FILS(1),
     &              id%NA(1), id%LNA, id%NE_STEPS(1), id%DAD_STEPS(1),
     &              id%ND_STEPS(1), id%PROCNODE_STEPS(1),
     &              id%NSLAVES,
     &              KEEP8(11), KEEP(26), KEEP(15),
     &              KEEP8(12),  ! formerly KEEP(16),
     &              KEEP8(14),  ! formerly KEEP(200),
     &              KEEP(224), KEEP(225),
     &              KEEP(27), RINFO(1),
     &              KEEP(1), KEEP8(1), LOCAL_M, LOCAL_N, SBUF_RECOLD8,
     &              SBUF_SEND, SBUF_REC, id%COST_SUBTREES, KEEP(28),
     &              locI_AM_CAND(1), max(KEEP(56),1),
     &              id%ISTEP_TO_INIV2(1), id%CANDIDATES(1,1), 
     &              INFO(1), INFO(2)
     &              ,KEEP8(15)
     &              ,MAX_SIZE_FACTOR_TMP, KEEP8(9) 
     &              ,ENTRIES_IN_FACTORS_LOC_MASTERS
     &              ,id%root%yes, id%root%NPROW, id%root%NPCOL
     &           )
            IF(ASSOCIATED(locI_AM_CAND)) NULLIFY(locI_AM_CAND)
            id%MAX_SURF_MASTER = KEEP8(15)
C
            KEEP8(19)=MAX_SIZE_FACTOR_TMP
            KEEP( 29 ) = KEEP(15) + 2* max(KEEP(12),10)
     &           * ( KEEP(15) / 100 + 1)
C     Relaxed value of size of IS is not needed internally;
C     we save it directly in INFO(19)
            INFO( 19 ) = KEEP(225) + 2* max(KEEP(12),10)
     &           * ( KEEP(225) / 100 + 1)
C     size of S
            KEEP8(13)  = KEEP8(12) + int(KEEP(12),8) *
     &           ( KEEP8(12) / 100_8 + 1_8 )
C     size of S
            KEEP8(17)  = KEEP8(14) + int(KEEP(12),8) *
     &           ( KEEP8(14) /100_8 +1_8)
C     KEEP8( 22 ) is the OLD maximum size of receive buffer 
C     that includes CB related communications.
C     KEEP( 43 ) : min size for send buffer
C     KEEP( 44 ) : min size for receive buffer
C     KEEP(43-44) kept for allocating buffers during
C                 factorization phase
         CALL MUMPS_ALLREDUCEI8 ( SBUF_RECOLD8, KEEP8(22), MPI_MAX,
     &                            id%COMM_NODES )
C     We do a max with KEEP(27)=maxfront because for small
C     buffers, we need at least one row of cb to be sent/
C     received.
         SBUF_SEND = max(SBUF_SEND,KEEP(27))
         SBUF_REC  = max(SBUF_REC ,KEEP(27))
         CALL MPI_ALLREDUCE (SBUF_REC, KEEP(44), 1, 
     &        MPI_INTEGER, MPI_MAX,
     &        id%COMM_NODES, IERR)
         IF (KEEP(48)==5) THEN
            KEEP(43)=KEEP(44)
         ELSE
            KEEP(43)=SBUF_SEND
         ENDIF
C     
         MIN_BUF_SIZE8 = KEEP8(22) / int(KEEP(238),8)
         MIN_BUF_SIZE8 = min( MIN_BUF_SIZE8, int(huge (KEEP(43)),8))
         MIN_BUF_SIZE  = int( MIN_BUF_SIZE8 )
         KEEP(44) = max(KEEP(44), MIN_BUF_SIZE)
         KEEP(43) = max(KEEP(43), MIN_BUF_SIZE)
            IF ( PROK ) THEN
               WRITE(MP,'(A,I10) ') 
     &              ' Estimated INTEGER space for factors         :',
     &              KEEP(26)
               WRITE(MP,'(A,I10) ') 
     &              ' INFO(3), est. complex space to store factors:',
     &              KEEP8(11)
               WRITE(MP,'(A,I10) ') 
     &              ' Estimated number of entries in factors      :',
     &              KEEP8(9)
               WRITE(MP,'(A,I10) ') 
     &              ' Current value of space relaxation parameter :',
     &              KEEP(12)
               WRITE(MP,'(A,I10) ') 
     &              ' Estimated size of IS (In Core factorization):',
     &              KEEP(29)
               WRITE(MP,'(A,I10) ') 
     &              ' Estimated size of S  (In Core factorization):',
     &              KEEP8(13)
               WRITE(MP,'(A,I10) ') 
     &              ' Estimated size of S  (OOC factorization)    :',
     &              KEEP8(17)
            END IF
      ELSE
C     ---------------------
C     Master is not working
C     ---------------------
         ENTRIES_IN_FACTORS_LOC_MASTERS = 0_8
         KEEP8(13) = 0_8
         KEEP(29) = 0
         KEEP8(17)= 0_8
         INFO(19) = 0
         KEEP8(11) = 0_8
         KEEP(26) = 0
         KEEP(27) = 0
         RINFO(1) = 0.0D0
      END IF
C     --------------------------------------
C     KEEP( 13 ) : Real arrowhead size
C     KEEP( 14 ) : Integer arrowhead size
C     INFO(3)/KEEP8( 11 ) : Estimated real space needed for factors
C     INFO(4)/KEEP( 26 )  : Estimated integer space needed for factors
C     INFO(5)/KEEP( 27 )  : Estimated max front size
C     KEEP8(109)          : Estimated number of entries in factor
C                         (based on ENTRIES_IN_FACTORS_LOC_MASTERS computed 
C                          during ZMUMPS_ANA_DISTM, where we assume 
C                          that each master of a node computes
C                          the complete factor size.
C     --------------------------------------
      CALL MUMPS_ALLREDUCEI8( ENTRIES_IN_FACTORS_LOC_MASTERS, 
     &     KEEP8(109), MPI_SUM, id%COMM)
      CALL MUMPS_ALLREDUCEI8( KEEP8(19), KEEP8(119),
     &     MPI_MAX, id%COMM)
      CALL MPI_ALLREDUCE( KEEP(27), KEEP(127), 1,
     &     MPI_INTEGER, MPI_MAX,
     &     id%COMM, IERR)
      CALL MPI_ALLREDUCE( KEEP(26), KEEP(126), 1,
     &     MPI_INTEGER, MPI_SUM,
     &     id%COMM, IERR)
      CALL MUMPS_REDUCEI8( KEEP8(11), KEEP8(111), MPI_SUM,
     &     MASTER, id%COMM )
      CALL MUMPS_SETI8TOI4( KEEP8(111), INFOG(3) )
C     --------------
C     Flops estimate
C     --------------
      CALL MPI_ALLREDUCE( RINFO(1), RINFOG(1), 1,
     &     MPI_DOUBLE_PRECISION, MPI_SUM,
     &     id%COMM, IERR)
      CALL MUMPS_SETI8TOI4( KEEP8(11), INFO(3) )
      INFO ( 4 ) = KEEP(  26 )
      INFO ( 5 ) = KEEP(  27 )
      INFO ( 7 ) = KEEP(  29 )
      CALL MUMPS_SETI8TOI4( KEEP8(13), INFO(8) )
      CALL MUMPS_SETI8TOI4( KEEP8(17), INFO(20) )
      CALL MUMPS_SETI8TOI4( KEEP8(9), INFO(24) )
      INFOG( 4 ) = KEEP( 126 )
      INFOG( 5 ) = KEEP( 127 )
      CALL MUMPS_SETI8TOI4( KEEP8(109), INFOG(20) )
      CALL ZMUMPS_DIAG_ANA(id%MYID, id%COMM, KEEP(1), KEEP8(1),
     &     INFO(1), INFOG(1), RINFO(1), RINFOG(1), ICNTL(1))
C     =========================
C     IN-CORE MEMORY STATISTICS
C     =========================
         OOC_STAT = KEEP(201)
         IF (KEEP(201) .NE. -1) OOC_STAT=0 ! We want in-core statistics
         PERLU_ON = .FALSE.     ! switch off PERLU to compute KEEP8(2)
         CALL ZMUMPS_MAX_MEM( KEEP(1), KEEP8(1),
     &        id%MYID, id%N, id%NELT, id%NA(1), id%LNA, id%NZ,
     &        id%NA_ELT,
     &        id%NSLAVES, TOTAL_MBYTES, .FALSE.,
     &        OOC_STAT, PERLU_ON, TOTAL_BYTES)
         KEEP8(2) = TOTAL_BYTES    
         PERLU_ON  = .TRUE.
         CALL ZMUMPS_MAX_MEM( KEEP(1), KEEP8(1),
     &        id%MYID, id%N, id%NELT, id%NA(1), id%LNA, id%NZ,
     &        id%NA_ELT,
     &        id%NSLAVES, TOTAL_MBYTES, .FALSE.,
     &        OOC_STAT, PERLU_ON, TOTAL_BYTES)
         IF ( PROK ) THEN
            WRITE(MP,'(A,I10) ')
     & ' Estimated space in MBYTES for IC factorization            :',
     &           TOTAL_MBYTES
         END IF
         id%INFO(15) = TOTAL_MBYTES
C     
C     Centralize memory statistics on the host
C     
C     INFOG(16) = after analysis, est. mem size in Mbytes for facto,
C     for the processor using largest memory
C     INFOG(17) = after analysis, est. mem size in Mbytes for facto,
C     sum over all processors
C     INFOG(18/19) = idem at facto.
C     
      CALL MUMPS_MEM_CENTRALIZE( id%MYID, id%COMM,
     &     id%INFO(15), id%INFOG(16), IRANK )
      IF ( PROKG ) THEN
         WRITE( MPG,'(A,I10) ')
     & ' ** Rank of proc needing largest memory in IC facto        :',
     &        IRANK
         WRITE( MPG,'(A,I10) ')
     & ' ** Estimated corresponding MBYTES for IC facto            :',
     &        id%INFOG(16)
         IF ( KEEP(46) .eq. 0 ) THEN
C     Host not working
            WRITE( MPG,'(A,I10) ')
     & ' ** Estimated avg. MBYTES per work. proc at facto (IC)     :'
     &           ,(id%INFOG(17)-id%INFO(15))/id%NSLAVES
         ELSE
            WRITE( MPG,'(A,I10) ')
     & ' ** Estimated avg. MBYTES per work. proc at facto (IC)     :'
     &           ,id%INFOG(17)/id%NSLAVES
         END IF
         WRITE(MPG,'(A,I10) ')
     & ' ** TOTAL     space in MBYTES for IC factorization         :'
     &        ,id%INFOG(17)
      END IF
C        =========================================
C        NOW COMPUTE OUT-OF-CORE MEMORY STATISTICS
C        (except when OOC_STAT is equal to -1 in
C        which case IC and OOC statistics are
C        identical)
C        =========================================
         OOC_STAT = KEEP(201)
#if defined(OLD_OOC_NOPANEL)
         IF (OOC_STAT .NE. -1) OOC_STAT=2
#else
         IF (OOC_STAT .NE. -1) OOC_STAT=1
#endif
         PERLU_ON = .FALSE.     ! PERLU NOT taken into account
C     Used to compute KEEP8(3) (minimum number of bytes for OOC)
         CALL ZMUMPS_MAX_MEM( KEEP(1), KEEP8(1),
     &        id%MYID, id%N, id%NELT, id%NA(1), id%LNA, id%NZ,
     &        id%NA_ELT,
     &        id%NSLAVES, TOTAL_MBYTES, .FALSE.,
     &        OOC_STAT, PERLU_ON, TOTAL_BYTES)
         KEEP8(3) = TOTAL_BYTES
         PERLU_ON  = .TRUE.     ! PERLU taken into account
         CALL ZMUMPS_MAX_MEM( KEEP(1), KEEP8(1),
     &        id%MYID, id%N, id%NELT, id%NA(1), id%LNA, id%NZ,
     &        id%NA_ELT,
     &        id%NSLAVES, TOTAL_MBYTES, .FALSE.,
     &        OOC_STAT, PERLU_ON, TOTAL_BYTES)
         id%INFO(17) = TOTAL_MBYTES
      CALL MUMPS_MEM_CENTRALIZE( id%MYID, id%COMM,
     &     id%INFO(17), id%INFOG(26), IRANK )
      IF ( PROKG  ) THEN
         WRITE( MPG,'(A,I10) ')
     & ' ** Rank of proc needing largest memory for OOC facto      :',
     &        IRANK
         WRITE( MPG,'(A,I10) ')
     & ' ** Estimated corresponding MBYTES for OOC facto           :',
     &        id%INFOG(26)
         IF ( KEEP(46) .eq. 0 ) THEN
C     Host not working
            WRITE( MPG,'(A,I10) ')
     & ' ** Estimated avg. MBYTES per work. proc at facto (OOC)    :'
     &           ,(id%INFOG(27)-id%INFO(15))/id%NSLAVES
         ELSE
            WRITE( MPG,'(A,I10) ')
     & ' ** Estimated avg. MBYTES per work. proc at facto (OOC)    :'
     &           ,id%INFOG(27)/id%NSLAVES
         END IF
         WRITE(MPG,'(A,I10) ')
     & ' ** TOTAL     space in MBYTES for OOC factorization        :'
     &        ,id%INFOG(27)
      END IF
c     #endif
C     -------------------------
C     Define a specific mapping
C     for the user
C     -------------------------
      IF ( id%MYID. eq. MASTER .AND. KEEP(54) .eq. 1 ) THEN
         IF (associated( id%MAPPING))
     &        deallocate( id%MAPPING)
         allocate( id%MAPPING(id%NZ), stat=allocok)
         IF ( allocok .GT. 0 ) THEN
            INFO(1) = -7
            INFO(2) = id%NZ
            IF ( LPOK ) THEN
               WRITE(LP, 150) 'id%MAPPING'
            END IF
            GOTO 92
         END IF
         allocate(IWtemp( id%N ), stat=allocok)
         IF ( allocok .GT. 0 ) THEN
            INFO(1)=-7
            INFO(2)=id%N
            IF ( LPOK ) THEN
               WRITE(LP, 150) 'IWtemp(N)'
            END IF
            GOTO 92
         END IF
         CALL ZMUMPS_BUILD_MAPPING(
     &        id%N, id%MAPPING(1),
     &        id%NZ, id%IRN(1),id%JCN(1), id%PROCNODE_STEPS(1),
     &        id%STEP(1),
     &        id%NSLAVES, id%SYM_PERM(1),
     &        id%FILS(1), IWtemp, id%KEEP(1),id%KEEP8(1),
     &        id%root%MBLOCK, id%root%NBLOCK,
     &        id%root%NPROW, id%root%NPCOL )
         deallocate( IWtemp )
 92      CONTINUE
      END IF
      CALL MUMPS_PROPINFO( ICNTL(1), INFO(1),
     &     id%COMM, id%MYID )
      IF ( INFO(1) .LT. 0 ) RETURN
      RETURN
 110  FORMAT(/' ****** ANALYSIS STEP ********'/)
 150  FORMAT(
     & /' ** FAILURE DURING ZMUMPS_ANA_DRIVER, DYNAMIC ALLOCATION OF',
     &     A30)
 160  FORMAT(
     & /' ** FAILURE DURING ZMUMPS_ANA_DRIVER, INTEGER OVERFLOW OF',
     &     A30)   
      END SUBROUTINE ZMUMPS_ANA_DRIVER
      SUBROUTINE ZMUMPS_DIST_AVOID_COPIES(N,NSLAVES,
     &     ICNTL,INFOG, NE, NFSIZ,
     &     FRERE, FILS,
     &     KEEP,KEEP8,PROCNODE,
     &     SSARBR,NBSA,PEAK,IERR
     &     )
      USE MUMPS_STATIC_MAPPING
      IMPLICIT NONE
      INTEGER N, NSLAVES, NBSA, IERR
      INTEGER ICNTL(40),INFOG(40),KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER NE(N),NFSIZ(N),FRERE(N),FILS(N),PROCNODE(N)
      INTEGER SSARBR(N)
      DOUBLE PRECISION PEAK
      CALL MUMPS_DISTRIBUTE(N,NSLAVES,
     &     ICNTL,INFOG, NE, NFSIZ,
     &     FRERE, FILS,
     &     KEEP,KEEP8,PROCNODE,
     &     SSARBR,NBSA,dble(PEAK),IERR
     &     )
      RETURN
      END SUBROUTINE ZMUMPS_DIST_AVOID_COPIES
      SUBROUTINE ZMUMPS_SET_PROCNODE(INODE, PROCNODE, VALUE, FILS, N)
      INTEGER, intent(in) :: INODE, N, VALUE
      INTEGER, intent(in) :: FILS(N)
      INTEGER, intent(inout) :: PROCNODE(N)
C     
      INTEGER IN
      IN=INODE
      DO WHILE ( IN > 0 )
         PROCNODE( IN ) = VALUE
         IN=FILS( IN )
      ENDDO
      RETURN
      END SUBROUTINE ZMUMPS_SET_PROCNODE
      SUBROUTINE ZMUMPS_ANA_CHECK_KEEP(id)
C     This subroutine decodes the control parameters,
C     stores them in the KEEP array, and performs a
C     consistency check on the KEEP array.
      USE ZMUMPS_STRUC_DEF
      IMPLICIT NONE
      TYPE(ZMUMPS_STRUC)  :: id
C     internal variables
      INTEGER   :: LP, MP, MPG, I
      INTEGER   :: MASTER
      LOGICAL   :: PROK, PROKG
      PARAMETER( MASTER = 0 )
      LP  = id%ICNTL( 1 )
      MP  = id%ICNTL( 2 )
      MPG = id%ICNTL( 3 )
C     LP     : errors
C     MP     : INFO
      PROK  = (( MP  .GT. 0 ).AND.(id%ICNTL(4).GE.2))
      PROKG = ( MPG .GT. 0 .and. id%MYID .eq. MASTER )
      PROKG = (PROKG.AND.(id%ICNTL(4).GE.2))
C     Fwd in facto
      IF (id%MYID.eq.MASTER) THEN
        id%KEEP(256) = id%ICNTL(7) ! copy ordering option
        id%KEEP(252) = id%ICNTL(32)
        IF (id%KEEP(252) < 0 .OR. id%KEEP(252) > 1 ) THEN
          id%KEEP(252) = 0
        ENDIF
C       Which factors to store
        id%KEEP(251) = id%ICNTL(31)
        IF (id%KEEP(251) < 0 .OR. id%KEEP(251) > 2 ) THEN
          id%KEEP(251)=0
        ENDIF
C       For unsymmetric matrices, if forward solve 
C       performed during facto,
C       no reason to store L factors at all. Reset
C       KEEP(251) accordingly... except if the user
C       tells that no solve is needed.
        IF (id%KEEP(50) .EQ. 0 .AND. id%KEEP(252).EQ.1) THEN
          IF (id%KEEP(251) .NE. 1) id%KEEP(251) = 2
        ENDIF
C       Symmetric case, even if no backward needed,
C       store all factors
        IF (id%KEEP(50) .NE.0 .AND. id%KEEP(251) .EQ. 2) THEN
          id%KEEP(251) = 0
        ENDIF
C       Case of solve not needed:
        IF (id%KEEP(251) .EQ. 1) THEN
          id%KEEP(201) = -1
C         In that case, id%ICNTL(22) will
C         be ignored in future phases
        ENDIF
        IF (id%KEEP(252).EQ.1) THEN
          id%KEEP(253) = id%NRHS
          IF (id%KEEP(253) .LE. 0) THEN
            id%INFO(1)=-42
            id%INFO(2)=id%NRHS
            RETURN
          ENDIF
        ELSE
          id%KEEP(253) = 0
        ENDIF
      ENDIF
      IF ( (id%KEEP(24).NE.0) .AND.
     &     id%NSLAVES.eq.1 ) THEN
         id%KEEP(24) = 0
         IF ( PROKG ) THEN
            WRITE(MPG, '(A)')
     &           ' Resetting candidate strategy to 0 because NSLAVES=1'
            WRITE(MPG, '(A)') ' '
         END IF
      END IF
      IF ( (id%KEEP(24).EQ.0) .AND.
     &     id%NSLAVES.GT.1 ) THEN
         id%KEEP(24) = 8
      ENDIF
      IF ( (id%KEEP(24).NE.0)  .AND. (id%KEEP(24).NE.1)  .AND.
     &     (id%KEEP(24).NE.8)  .AND. (id%KEEP(24).NE.10) .AND.
     &     (id%KEEP(24).NE.12) .AND. (id%KEEP(24).NE.14) .AND.
     &     (id%KEEP(24).NE.16) .AND. (id%KEEP(24).NE.18)) THEN
         id%KEEP(24) = 8
         IF ( PROKG ) THEN
            WRITE(MPG, '(A)')
     &           ' Resetting candidate strategy to 8 '
            WRITE(MPG, '(A)') ' '
         END IF
      END IF
C****************************************************
C     
C     The master is doing most of the work
C     
C     NOTE:  Treatment of the errors on the master=
C     Go to the next SPMD part of the code in which
C     the first statement must be a call to PROPINFO
C     
C****************************************************
C     =========================================
C     Check (raise error or modify) some input
C     parameters or KEEP values on the master.
C     =========================================
      id%KEEP8(21) = int(id%KEEP(85),8)
      IF ( id%MYID .EQ. MASTER ) THEN
C     -- OOC/Incore strategy 
        IF (id%KEEP(201).NE.-1) THEN
          id%KEEP(201)=id%ICNTL(22)
          IF (id%KEEP(201) .GT. 0) THEN
#if defined(OLD_OOC_NOPANEL)
            id%KEEP(201)=2
#else
            id%KEEP(201)=1
#endif
          ENDIF
        ENDIF
C     
C     ----------------------------
C     Save id%ICNTL(18) (distributed
C     matrix on entry) in id%KEEP(54)
C     ----------------------------
         id%KEEP(54) = id%ICNTL(18)
         IF ( id%KEEP(54) .LT. 0 .or. id%KEEP(54).GT.3 ) THEN
            IF ( PROKG ) THEN
               WRITE(MPG, *) ' Out-of-range value for id%ICNTL(18).'
               WRITE(MPG, *) ' Used 0 ie matrix not distributed'
            END IF
            id%KEEP(54) = 0
         END IF
C     -----------------------------------------
C     Save id%ICNTL(5) (matrix format) in id%KEEP(55)
C     -----------------------------------------
         id%KEEP(55) = id%ICNTL(5)
         IF ( id%KEEP(55) .LT. 0 .OR. id%KEEP(55) .GT. 1 ) THEN
            IF ( PROKG ) THEN
               WRITE(MPG, *) ' Out-of-range value for id%ICNTL(5).'
               WRITE(MPG, *) ' Used 0 ie matrix is assembled'
            END IF
            id%KEEP(55) = 0
         END IF
         id%KEEP(60) = id%ICNTL(19)
         IF ( id%KEEP( 60 ) .LE. 0 ) id%KEEP( 60 ) = 0
         IF ( id%KEEP( 60 ) .GT. 3 ) id%KEEP( 60 ) = 0
         IF (id%KEEP(60) .NE. 0 .AND. id%SIZE_SCHUR == 0 ) THEN
            WRITE(MPG,'(A)')
     &           ' ** Schur option ignored because SIZE_SCHUR=0'
            id%KEEP(60)=0
         END IF
C        ---------------------------------------
C        Save SIZE_SCHUR in a KEEP, for possible
C        check at factorization and solve phases
C        ---------------------------------------
         IF ( id%KEEP(60) .NE.0 ) THEN
            id%KEEP(116) = id%SIZE_SCHUR
            IF (id%SIZE_SCHUR .LT. 0 .OR. id%SIZE_SCHUR .GE. id%N) THEN
              id%INFO(1)=-49
              id%INFO(2)=id%SIZE_SCHUR
              RETURN
            ENDIF
C           List of Schur variables provided by user.
            IF ( .NOT. associated( id%LISTVAR_SCHUR ) ) THEN
               id%INFO(1) = -22
               id%INFO(2) = 8
               RETURN
            ELSE IF (size(id%LISTVAR_SCHUR)<id%SIZE_SCHUR) THEN
               id%INFO(1) = -22
               id%INFO(2) = 8
               RETURN
            END IF
         ENDIF
         IF (id%KEEP(60) .EQ. 3 .AND. id%KEEP(50).NE.0) THEN
            IF (id%MBLOCK > 0 .AND. id%NBLOCK > 0 .AND.
     &           id%NPROW > 0 .AND. id%NPCOL > 0 ) THEN
               IF (id%NPROW *id%NPCOL .LE. id%NSLAVES) THEN
C     We will eventually have to "symmetrize the
C     Schur complement. For that NBLOCK and MBLOCK
C     must be equal.
                  IF (id%MBLOCK .NE. id%NBLOCK ) THEN
                     id%INFO(1)=-31
                     id%INFO(2)=id%MBLOCK - id%NBLOCK
                     RETURN
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
C     Check the ordering strategy and compatibility with
C     other control parameters
         id%KEEP(244) = id%ICNTL(28)
         id%KEEP(245) = id%ICNTL(29)
#if ! defined(parmetis) && ! defined(parmetis3)        
         IF ((id%KEEP(244) .EQ. 2) .AND. (id%KEEP(245) .EQ. 2)) THEN
            id%INFO(1)  = -38
            IF(id%MYID .EQ.0 ) THEN
               WRITE(LP,'("ParMETIS not available.")')
               WRITE(LP,'("Aborting.")')
               RETURN
            END IF
         END IF
#endif
#if ! defined(ptscotch)         
         IF ((id%KEEP(244) .EQ. 2) .AND. (id%KEEP(245) .EQ. 1)) THEN
            id%INFO(1)  = -38
            IF(id%MYID .EQ.0 ) THEN
               WRITE(LP,'("PT-SCOTCH not available.")')
               WRITE(LP,'("Aborting.")')
               RETURN
            END IF
         END IF
#endif
C     Analysis strategy is set to automatic.
C     Check for availability of parallel ordering tools
C     otherwise set it to sequential.
         IF((id%KEEP(244) .GT. 2) .OR.
     &        (id%KEEP(244) .LT. 0)) id%KEEP(244)=0
         IF(id%KEEP(244) .EQ. 0) THEN
            id%KEEP(244) = 1
         ELSE IF (id%KEEP(244) .EQ. 2) THEN
            IF(id%KEEP(55) .NE. 0) THEN
               id%INFO(1)  = -39
               WRITE(LP,
     &              '("Incompatible values for ICNTL(5), ICNTL(28)")')
               WRITE(LP,
     &              '("Parallel analysis is not possible if the")')
               WRITE(LP,
     &              '("matrix is not assembled")')
               RETURN
            ELSE IF(id%KEEP(60) .NE. 0) THEN
               id%INFO(1)  = -39
               WRITE(LP,
     &              '("Incompatible values for ICNTL(19), ICNTL(28)")')
               WRITE(LP,
     &              '("Parallel analysis is not possible if SCHUR")')
               WRITE(LP,
     &              '("complement must be returned")')
               RETURN
            END IF
C     In the case where there are too few processes to do
C     the parallel analysis we simply revert to sequential version
            IF(id%NSLAVES .LT. 2) THEN
               id%KEEP(244) = 1
               IF(PROKG) WRITE(MPG,
     &              '("Too few processes.
     & Reverting to sequential analysis")',advance='no')
               IF(id%KEEP(245) .EQ. 1) THEN
                  IF(PROKG) WRITE(MPG, '(" with SCOTCH")')
                  id%KEEP(256) = 3
               ELSE IF(id%KEEP(245) .EQ. 2) THEN
                  IF(PROKG) WRITE(MPG, '(" with Metis")')
                  id%KEEP(256) = 5
               ELSE
                  IF(PROKG) WRITE(MPG, '(".")')
                  id%KEEP(256) = 7
               END IF
            END IF
         END IF
         id%INFOG(32) = id%KEEP(244)
         IF ( (id%KEEP(244) .EQ. 1) .AND.
     &        (id%KEEP(256) .EQ. 1) ) THEN
C     ordering given, PERM_IN must be of size N
            IF ( .NOT. associated( id%PERM_IN ) ) THEN
               id%INFO(1) = -22
               id%INFO(2) = 3
               RETURN
            ELSE IF ( size( id%PERM_IN ) < id%N ) THEN
               id%INFO(1) = -22
               id%INFO(2) = 3
               RETURN
            END IF
         ENDIF
C     Check KEEP(9-10) for level 2
         IF (id%KEEP(9) .LE. 1 ) id%KEEP(9) = 500
         IF ( id%KEEP8(21) .GT. 0_8 ) THEN 
            IF ((id%KEEP8(21).LE.1_8) .OR.
     &          (id%KEEP8(21).GT.int(id%KEEP(9),8)))
     &         id%KEEP8(21) = int(min(id%KEEP(9),100),8)
         ENDIF
C     
         IF (id%KEEP(48). EQ. 1 ) id%KEEP(48) = -12345
C     
         IF ( (id%KEEP(48).LT.0) .OR. (id%KEEP(48).GT.5) ) THEN
            id%KEEP(48)=5
         ENDIF
C     Schur 
C     Given ordering must be compatible with Schur variables.
         IF ( (id%KEEP(60) .NE. 0) .AND. (id%KEEP(256) .EQ. 1) ) THEN
            DO I = 1, id%SIZE_SCHUR
               IF (id%PERM_IN(id%LISTVAR_SCHUR(I))
     &              .EQ. id%N-id%SIZE_SCHUR+I)
     &              CYCLE
C              -------------------------------
C              Problem with PERM_IN: -22/3
C              Above constrained explained in
C              doc of PERM_IN in user guide.
C              -------------------------------
               id%INFO(1) = -4
               id%INFO(2) = id%LISTVAR_SCHUR(I)
               RETURN
               IF (MPG.GT.0) THEN
                  WRITE(MPG,'(A)')
     & ' ** Ignoring user-ordering, because incompatible with Schur.'
                  WRITE(MPG,'(A)') ' ** id%ICNTL(7) treated as 0.'
               END IF
               EXIT
            ENDDO
         END IF
C     
C     Note that schur is not compatible with
C     
C     1/Max-trans DONE
C     2/Null space
C     3/Ordering given DONE
C     4/Scaling
C     5/Iterative Refinement
C     6/Error analysis
C     7/Parallel Analysis
*     
*     Graph modification prior to ordering (id%ICNTL(12) option)
*     id%KEEP (95) will hold the eventually modified value of id%ICNTL(12)
*     
         id%KEEP(95) = id%ICNTL(12)
         IF (id%KEEP(50).NE.2) id%KEEP(95) = 1
         IF ((id%KEEP(95).GT.3).OR.(id%KEEP(95).LT.0)) id%KEEP(95) = 0
C     MAX-TRANS
C     
C     id%KEEP (23) will hold the eventually modified value of id%ICNTL(6)
C     (maximum transversal if >= 1)
C     
         id%KEEP(23) = id%ICNTL(6)
C     
C     
C     --------------------------------------------
C     Avoid max-trans unsymmetric permutation in case of
C     ordering is given,
C     or matrix is in element form, or Schur is asked
C     or initial matrix is distributed
C     --------------------------------------------
         IF (id%KEEP(23).LT.0.OR.id%KEEP(23).GT.7) id%KEEP(23) = 7
C        still forbid max trans for LLT
         IF ( id%KEEP(50) .EQ. 1 ) THEN
            IF (id%KEEP(23) .NE. 0) THEN
               IF (MPG.GT.0) THEN
                  WRITE(MPG,'(A)')
     & ' ** Max-trans not compatible with LLT factorization'
               END IF
               id%KEEP(23) = 0
            ENDIF
            IF (id%KEEP(95) .GT. 1) THEN
               IF (MPG.GT.0) THEN
                  WRITE(MPG,'(A)')
     & ' ** ICNTL(12) ignored: not compatible with LLT factorization'
               END IF
            ENDIF
            id%KEEP(95) = 1
         END IF
C     
         IF  (id%KEEP(60) .GT. 0) THEN
            IF (id%KEEP(23) .NE. 0) THEN
               IF (MPG.GT.0) THEN
                  WRITE(MPG,'(A)')
     &                 ' ** Max-trans not allowed because of Schur'
               END IF
               id%KEEP(23) = 0
            ENDIF
            IF (id%KEEP(52).NE.0) THEN
               IF (MPG.GT.0) THEN
                  WRITE(MPG,'(A)')
     & ' ** Scaling during analysis not allowed because of Schur'
               ENDIF
               id%KEEP(52) = 0
            ENDIF
C     also forbid compressed/constrained ordering...
            IF (id%KEEP(95) .GT. 1) THEN
               IF (MPG.GT.0) THEN
                  WRITE(MPG,'(A)')
     & ' ** ICNTL(12) option not allowed because of Schur'
               END IF
            ENDIF
            id%KEEP(95) = 1
         END IF
         IF ( (id%KEEP(23) .NE. 0) .AND. (id%KEEP(256).EQ.1)) THEN
            id%KEEP(23) = 0
            id%KEEP(95) = 1
            IF (MPG.GT.0) THEN
               WRITE(MPG,'(A)')
     & ' ** Max-trans not allowed because ordering is given'
            END IF
         END IF
         IF ( id%KEEP(256) .EQ. 1 ) THEN
            IF (id%KEEP(95) > 1 .AND. MPG.GT.0) THEN
               WRITE(MPG,'(A)')
     & ' ** ICNTL(12) option incompatible with given ordering'
            END IF
            id%KEEP(95) = 1
         END IF
         IF (id%KEEP(54) .NE. 0) THEN
            IF( id%KEEP(23) .NE. 0 ) THEN
               IF (MPG.GT.0) THEN
                  WRITE(MPG,'(A)')
     & ' ** Max-trans not allowed because matrix is distributed'
               END IF
               id%KEEP(23) = 0
            ENDIF
            IF (id%KEEP(52).EQ.-2) THEN
C           Only Ruiz & Bora scaling available for dist format 
C           (Work supported by ANR-SOLSTICE (ANR-06-CIS6-010))
               IF (MPG.GT.0) THEN
                  WRITE(MPG,'(A)')
     & ' ** Scaling during analysis not allowed (matrix is distributed)'
               ENDIF
            ENDIF
            id%KEEP(52) = 0
            IF (id%KEEP(95) .GT. 1 .AND. MPG.GT.0) THEN
               WRITE(MPG,'(A)')
     & ' ** ICNTL(12) option not allowed because matrix is
     &distributed'
            ENDIF
            id%KEEP(95) = 1
         END IF
         IF ( id%KEEP(55) .NE. 0 ) THEN
            IF( id%KEEP(23) .NE. 0 ) THEN
               IF (MPG.GT.0) THEN
                  WRITE(MPG,'(A)')
     & ' ** Max-trans not allowed for element matrix'
               END IF
               id%KEEP(23) = 0
            ENDIF
            IF (MPG.GT.0 .AND. id%KEEP(52).EQ.-2) THEN
               WRITE(MPG,'(A)')
     & ' ** Scaling not allowed at analysis for element matrix'
            ENDIF
            id%KEEP(52) = 0
            id%KEEP(95) = 1
         ENDIF
C     In the case where parallel analysis is done, column permutation
C     is not allowed
         IF(id%KEEP(244) .EQ. 2) THEN
            IF(id%KEEP(23) .EQ. 7) THEN
C     Automatic hoice: set it to 0
               id%KEEP(23) = 0
            ELSE IF (id%KEEP(23) .GT. 0) THEN
               id%INFO(1)  = -39
               id%KEEP(23) = 0
               WRITE(LP,
     &              '("Incompatible values for ICNTL(6), ICNTL(28)")')
               WRITE(LP,
     &              '("Maximum transversal not allowed
     &                 in parallel analysis")')
               RETURN
            END IF
         END IF
C     --------------------------------------------
C     Avoid distributed entry for element matrix.
C     --------------------------------------------
         IF ( id%KEEP(54) .NE. 0 .AND. id%KEEP(55) .NE. 0 ) THEN
            id%KEEP(54) = 0
            IF (MPG.GT.0) THEN
               WRITE(MPG,'(A)')
     & ' ** Distributed entry not available for element matrix'
            END IF
         ENDIF
C     ----------------------------------
C     Choice of symbolic analysis option
C     ----------------------------------
         IF (id%ICNTL(39).NE.1 .and. id%ICNTL(39).NE.2) THEN
            id%KEEP(106)=1
C     Automatic choice leads to new symbolic
C     factorization except(see below) if KEEP(256)==1.
         ELSE
            id%KEEP(106)=id%ICNTL(39)
         ENDIF
C     modify input parameters to avoid incompatible
C     input data between ordering, scaling and maxtrans
C     note that if id%ICNTL(12)/id%KEEP(95) = 0 then
C     the automatic choice will be done in ANA_O
         IF(id%KEEP(50) .EQ. 2) THEN
C     LDLT case
            IF( .NOT. associated(id%A) ) THEN
C     constraint ordering can be computed only if values are
C     given to analysis
               IF(id%KEEP(95) .EQ. 3) THEN
                  id%KEEP(95) = 2
               ENDIF
            ENDIF
            IF(id%KEEP(95) .EQ. 3 .AND. id%KEEP(256) .NE. 2) THEN
C     if constraint and ordering is not AMF then use compress
               IF (PROK) WRITE(MP,*)
     &              'WARNING: ZMUMPS_ANA_O constrained ordering not ', 
     &              'available with selected ordering'
               id%KEEP(95) = 2
            ENDIF 
            IF(id%KEEP(95) .EQ. 3) THEN
C     if constraint ordering required then we need to compute scaling
C     and max trans
C     NOTE that if we enter this condition then
C     id%A is associated because of the test above:
C     (IF( .NOT. associated(id%A) ) THEN)
               id%KEEP(23) = 5
               id%KEEP(52) = -2
            ELSE IF(id%KEEP(95) .EQ. 2 .AND. 
     &              (id%KEEP(23) .EQ. 0 .OR. id%KEEP(23) .EQ. 7) ) THEN
C     compressed ordering requires max trans but not necessary scaling
               IF( associated(id%A) ) THEN
                  id%KEEP(23) = 5
               ELSE
C     we can do compressed ordering without
C     information on the numerical values:
C     a maximum transversal already provides
C     information on the location of off-diagonal
C     nonzeros which can be candidates for 2x2
C     pivots
                  id%KEEP(23) = 1
               ENDIF
            ELSE IF(id%KEEP(95) .EQ. 1) THEN
               id%KEEP(23) = 0
            ELSE IF(id%KEEP(95) .EQ. 0 .AND. id%KEEP(23) .EQ. 0) THEN
C     if max trans desactivated then the automatic choice for type of ord
C     is set to 1, which means that we will use usual ordering
C     (no constraints or compression)
               id%KEEP(95) = 1
            ENDIF
         ELSE
            id%KEEP(95) = 1
         ENDIF
C     --------------------------------
C     Save ICNTL(16) (QR) in KEEP(53)
C     Will be broadcasted to all other
C     nodes in routine ZMUMPS_BDCAST
C     --------------------------------
         id%KEEP(53)=0
         IF(id%KEEP(86).EQ.1)THEN
C     Force the exchange of both the memory and flops information during 
C     the factorization
            IF(id%KEEP(47).LT.2) id%KEEP(47)=2
         ENDIF
         IF(id%KEEP(48).EQ.5)THEN
            IF(id%KEEP(50).EQ.0)THEN
               id%KEEP(87)=50
               id%KEEP(88)=50
            ELSE
               id%KEEP(87)=70
               id%KEEP(88)=70
            ENDIF
         ENDIF
         IF((id%NSLAVES.EQ.1).AND.(id%KEEP(76).GT.3))THEN
            id%KEEP(76)=2
         ENDIF
         IF(id%KEEP(81).GT.0)THEN
            IF(id%KEEP(47).LT.2) id%KEEP(47)=2
         ENDIF
      END IF
      RETURN
      END SUBROUTINE ZMUMPS_ANA_CHECK_KEEP
      SUBROUTINE ZMUMPS_GATHER_MATRIX(id)
C     This subroutine gathers a distributed matrix
C     on the host node 
      USE ZMUMPS_STRUC_DEF
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      TYPE(ZMUMPS_STRUC)  :: id
C     local variables
      INTEGER, ALLOCATABLE :: REQPTR(:,:)
      INTEGER              :: MASTER, IERR, INDX, NRECV
      INTEGER              :: STATUS(MPI_STATUS_SIZE)
      INTEGER              :: LP, MP, MPG, I
      LOGICAL              :: PROK, PROKG
      PARAMETER( MASTER = 0 )
      LP  = id%ICNTL( 1 )
      MP  = id%ICNTL( 2 )
      MPG = id%ICNTL( 3 )
C     LP     : errors
C     MP     : INFO
      PROK  = (( MP  .GT. 0 ).AND.(id%ICNTL(4).GE.2))
      PROKG = ( MPG .GT. 0 .and. id%MYID .eq. MASTER )
      PROKG = (PROKG.AND.(id%ICNTL(4).GE.2))
      IF ( id%KEEP(46) .EQ. 0 .AND. id%MYID .EQ. MASTER ) THEN
C     host-node mode: master has no entries.
         id%NZ_loc = 0
      END IF
      IF ( id%MYID .eq. MASTER ) THEN
C     -----------------------------------
C     Allocate a small array for requests
C     and pointers into IRN/JCN
C     -----------------------------------
         ALLOCATE( REQPTR( id%NPROCS, 3 ), STAT = IERR )
         IF ( IERR .GT. 0 ) THEN
            id%INFO(1) = -7
            id%INFO(2) = 3 * id%NPROCS
            IF ( LP .GT. 0 ) THEN
               WRITE(LP, 150) 'REQPTR'
            END IF
            GOTO 13
         END IF
C     --------------------
C     Allocate now IRN/JCN
C     --------------------
         ALLOCATE( id%IRN( id%NZ ), STAT = IERR )
         IF ( IERR .GT. 0 ) THEN
            id%INFO(1) = -7
            id%INFO(2) = id%NZ
            IF ( LP .GT. 0 ) THEN
               WRITE(LP, 150) 'IRN'
            END IF
            GOTO 13
         END IF
         ALLOCATE( id%JCN( id%NZ ), STAT = IERR )
         IF ( IERR .GT. 0 ) THEN
            id%INFO(1) = -7
            id%INFO(2) = id%NZ
            IF ( LP .GT. 0 ) THEN
               WRITE(LP, 150) 'JCN'
            END IF
            GOTO 13
         END IF
      END IF
 13   CONTINUE
C     Propagate errors
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     id%COMM, id%MYID )
      IF ( id%INFO(1) < 0 ) RETURN
C     -------------------------------------
C     Get numbers of non-zeros for everyone
C     -------------------------------------
      IF ( id%MYID .EQ. MASTER ) THEN
         DO I = 1, id%NPROCS - 1
            CALL MPI_RECV( REQPTR( I+1, 1 ), 1, 
     &           MPI_INTEGER, I,
     &           COLLECT_NZ, id%COMM, STATUS, IERR )
         END DO
         IF ( id%KEEP(46) .eq. 0 ) THEN
            REQPTR( 1, 1 ) = 1
         ELSE
            REQPTR( 1, 1 ) = id%NZ_loc + 1
         END IF
C     --------------
C     Build pointers
C     --------------
         DO I = 2, id%NPROCS
            REQPTR( I, 1 ) = REQPTR( I, 1 ) + REQPTR( I-1, 1 )
         END DO
      ELSE
         CALL MPI_SEND( id%NZ_loc, 1, MPI_INTEGER, MASTER,
     &        COLLECT_NZ, id%COMM, IERR )
      END IF
C     -----------------------------------------------
C     Bottleneck is here master; use synchronous send
C     for slaves, but asynchronous receives on master
C     -----------------------------------------------
      IF ( id%MYID .eq. MASTER ) THEN
         NRECV = 0
         DO I = 1, id%NPROCS - 1
            IF ( REQPTR( I + 1, 1 ) - REQPTR( I, 1 ) .NE. 0 ) THEN
               NRECV = NRECV + 2
               CALL MPI_IRECV( id%IRN( REQPTR( I, 1 ) ),
     &              REQPTR( I + 1, 1 ) - REQPTR( I, 1 ), 
     &              MPI_INTEGER,
     &              I, COLLECT_IRN, id%COMM, REQPTR(I, 2), IERR )
               CALL MPI_IRECV( id%JCN( REQPTR( I, 1 ) ),
     &              REQPTR( I + 1, 1 ) - REQPTR( I, 1 ),   
     &              MPI_INTEGER,
     &              I, COLLECT_JCN, id%COMM, REQPTR(I, 3), IERR )
            ELSE
C     ------------------
C     Nothing to receive
C     ------------------
               REQPTR(I, 2) = MPI_REQUEST_NULL
               REQPTR(I, 3) = MPI_REQUEST_NULL
            END IF
         END DO
      ELSE
C     -----------------------------
C     Send only if size is not zero
C     -----------------------------
         IF ( id%NZ_loc .NE. 0 ) THEN
            CALL MPI_SEND( id%IRN_loc(1), id%NZ_loc, 
     &           MPI_INTEGER, MASTER,
     &           COLLECT_IRN, id%COMM, IERR )
            CALL MPI_SEND( id%JCN_loc(1), id%NZ_loc, 
     &           MPI_INTEGER, MASTER,
     &           COLLECT_JCN, id%COMM, IERR )
         END IF
      END IF
      IF ( id%MYID .eq. MASTER ) THEN
C     ------------------------------------
C     While master is receiving indices
C     from everybody, try to do the local
C     copies for better overlap. If master
C     has other things to do, he could try
C     to do them here.
C     ------------------------------------
         IF ( id%NZ_loc .NE. 0 ) THEN
            DO I=1,id%NZ_loc
               id%IRN(I) = id%IRN_loc(I)
               id%JCN(I) = id%JCN_loc(I)
            ENDDO
         END IF
         REQPTR( id%NPROCS, 2 ) = MPI_REQUEST_NULL
         REQPTR( id%NPROCS, 3 ) = MPI_REQUEST_NULL
         DO I = 1, NRECV
            CALL MPI_WAITANY
     &           ( 2 * id%NPROCS, REQPTR( 1, 2 ), INDX, STATUS, IERR )
         END DO
         DEALLOCATE( REQPTR )
      END IF
      RETURN
 150  FORMAT(
     &/' ** FAILURE DURING ZMUMPS_GATHER_MATRIX, DYNAMIC ALLOCATION OF',
     &     A30)
      END SUBROUTINE ZMUMPS_GATHER_MATRIX
      SUBROUTINE ZMUMPS_DUMP_PROBLEM(id)
      USE ZMUMPS_STRUC_DEF
      IMPLICIT NONE
C
C     Purpose:
C     =======
C
C     If id%WRITE_PROBLEM has been set by the user,
C     possibly on all processors in case of distributed
C     matrix, opens a file and dumps the matrix and/or
C     the right hand side. This subroutine calls
C     ZMUMPS_DUMP_MATRIX and ZMUMPS_DUMP_RHS.
C     The routine should be called on all processors.
C
      INCLUDE 'mpif.h'
C     Arguments
C     =========
      TYPE(ZMUMPS_STRUC)  :: id
C
C     Local variables
C     ===============
C
      INTEGER              :: MASTER, IERR
      INTEGER              :: IUNIT
      LOGICAL              :: IS_ELEMENTAL
      LOGICAL              :: IS_DISTRIBUTED
      INTEGER              :: MM_WRITE
      INTEGER              :: MM_WRITE_CHECK
      CHARACTER(LEN=20)    :: MM_IDSTR
      LOGICAL              :: I_AM_SLAVE, I_AM_MASTER
      PARAMETER( MASTER = 0 )
      IUNIT = 69
      I_AM_SLAVE = ( id%MYID .NE. MASTER  .OR.
     &     ( id%MYID .EQ. MASTER .AND.
     &     id%KEEP(46) .EQ. 1 ) )
      I_AM_MASTER = (id%MYID.EQ.MASTER)
C     Remark: if id%KEEP(54) = 1 or 2, the structure
C     is centralized at analysis. Since ZMUMPS_DUMP_PROBLEM
C     is called at analysis phase, we define IS_DISTRIBUTED
C     as below, which implies that the structure of the problem
C     is distributed in IRN_loc/JCN_loc at analysis.
C     equal to 
      IS_DISTRIBUTED = (id%KEEP(54) .EQ. 3)
      IS_ELEMENTAL   = (id%KEEP(55) .NE. 0)
      IF (id%MYID.EQ.MASTER .AND. .NOT. IS_DISTRIBUTED) THEN
C        ====================
C        Matrix is assembled
C        and centralized
C        ====================
        IF (id%WRITE_PROBLEM(1:20) .NE. "NAME_NOT_INITIALIZED")THEN
          OPEN(IUNIT,FILE=trim(id%WRITE_PROBLEM))
          CALL ZMUMPS_DUMP_MATRIX( id, IUNIT, I_AM_SLAVE, I_AM_MASTER,
     &           IS_DISTRIBUTED,        ! = .FALSE., centralized
     &           IS_ELEMENTAL )         ! Elemental or not
          CLOSE(IUNIT)
        ENDIF
      ELSE IF (id%KEEP(54).EQ.3) THEN
C        =====================
C        Matrix is distributed
C        =====================
         IF (id%WRITE_PROBLEM(1:20) .EQ. "NAME_NOT_INITIALIZED"
     &        .OR. .NOT. I_AM_SLAVE )THEN
            MM_WRITE = 0
         ELSE
            MM_WRITE = 1
         ENDIF
         CALL MPI_ALLREDUCE(MM_WRITE, MM_WRITE_CHECK, 1,
     &        MPI_INTEGER, MPI_SUM, id%COMM, IERR)
C        -----------------------------------------
C        If yes, each processor writes its share
C        of the matrix in a file in matrix market
C        format (otherwise nothing written). We
C        append the process id to the filename.
C        Safer in case all filenames are the
C        same if all processors share the same
C        file system.
C        -----------------------------------------
         IF (MM_WRITE_CHECK.EQ.id%NSLAVES .AND. I_AM_SLAVE) THEN
            WRITE(MM_IDSTR,'(I9)') id%MYID_NODES
            OPEN(IUNIT,
     &           FILE=trim(id%WRITE_PROBLEM)//trim(adjustl(MM_IDSTR)))
            CALL ZMUMPS_DUMP_MATRIX(id, IUNIT, I_AM_SLAVE, I_AM_MASTER,
     &           IS_DISTRIBUTED,           ! =.TRUE., distributed
     &           IS_ELEMENTAL )            ! Elemental or not
            CLOSE(IUNIT)
         ENDIF
C     ELSE ...
C     Nothing written in other cases.
      ENDIF
C     ===============
C     Right-hand side
C     ===============
      IF ( id%MYID.EQ.MASTER .AND.
     &     associated(id%RHS) .AND.
     &     id%WRITE_PROBLEM(1:20)
     &     .NE. "NAME_NOT_INITIALIZED")THEN
        OPEN(IUNIT,FILE=trim(id%WRITE_PROBLEM) //".rhs")
        CALL ZMUMPS_DUMP_RHS(IUNIT, id)
        CLOSE(IUNIT)
      ENDIF
      RETURN
      END SUBROUTINE ZMUMPS_DUMP_PROBLEM
      SUBROUTINE ZMUMPS_DUMP_MATRIX
     & (id, IUNIT, I_AM_SLAVE, I_AM_MASTER,
     &  IS_DISTRIBUTED, IS_ELEMENTAL )
      USE ZMUMPS_STRUC_DEF
      IMPLICIT NONE
C
C  Purpose:
C  =======
C     This subroutine dumps a routine in matrix-market format
C     if the matrix is assembled, and in "MUMPS" format (see
C     example in the MUMPS users'guide, if the matrix is
C     centralized and elemental).
C     The routine can be called on all processors. In case of
C     distributed assembled matrix, each processor writes its
C     share as a matrix market file on IUNIT (IUNIT may have
C     different values on different processors).
C
C  Arguments (input parameters)
C  ============================
C
C     IUNIT: should be set to the Fortran unit where
C            data should be written.
C     I_AM_SLAVE: .TRUE. except on a non working master
C     IS_DISTRIBUTED: .TRUE. if matrix is distributed,
C                     i.e., if IRN_loc/JCN_loc are provided.
C     IS_ELEMENTAL  : .TRUE. if matrix is elemental
C     id            : main MUMPS structure
C
      LOGICAL, intent(in) :: I_AM_SLAVE,
     &                       I_AM_MASTER,
     &                       IS_DISTRIBUTED,
     &                       IS_ELEMENTAL
      INTEGER, intent(in) :: IUNIT
      TYPE(ZMUMPS_STRUC), intent(in)  :: id
C
C  Local variables:
C  ===============
C
      CHARACTER (LEN=10)   :: SYMM
      CHARACTER (LEN=8)    :: ARITH
      INTEGER              :: I
C
C  Executable statements:
C  =====================
      IF (IS_ELEMENTAL) THEN
        RETURN
      ENDIF
      IF (I_AM_MASTER .AND. .NOT. IS_DISTRIBUTED) THEN
C        ==================
C        CENTRALIZED MATRIX
C        ==================
         IF (associated(id%A)) THEN
C     Write header line:
               ARITH='complex'
C
         ELSE
            ARITH='pattern '
         ENDIF
         IF (id%KEEP(50) .eq. 0) THEN
            SYMM="general"
         ELSE
            SYMM="symmetric"
         END IF
         WRITE(IUNIT,FMT=*)'%%MatrixMarket matrix coordinate ',
     &           trim(ARITH)," ",trim(SYMM)
         WRITE(IUNIT,*) id%N, id%N, id%NZ
         IF (associated(id%A)) THEN
            DO I=1,id%NZ
               IF (id%KEEP(50).NE.0 .AND. id%IRN(I).LT.id%JCN(I)) THEN
C              permute upper diag entry
                     WRITE(IUNIT,*) id%JCN(I), id%IRN(I), 
     &                    dble(id%A(I)), aimag(id%A(I))
               ELSE
                     WRITE(IUNIT,*) id%IRN(I), id%JCN(I), 
     &                    dble(id%A(I)), aimag(id%A(I))
               ENDIF
            ENDDO
         ELSE
C           pattern only
            DO I=1,id%NZ
               IF (id%KEEP(50).NE.0 .AND. id%IRN(I).LT.id%JCN(I)) THEN
C                 permute upper diag entry
                  WRITE(IUNIT,*) id%JCN(I), id%IRN(I)
               ELSE
                     WRITE(IUNIT,*) id%IRN(I), id%JCN(I)
               ENDIF
            ENDDO
         ENDIF
      ELSE IF ( IS_DISTRIBUTED .AND. I_AM_SLAVE ) THEN
C        ==================
C        DISTRIBUTED MATRIX
C        ==================
         IF (associated(id%A_loc)) THEN
               ARITH='complex'
         ELSE
               ARITH='pattern '
         ENDIF
         IF (id%KEEP(50) .eq. 0) THEN
            SYMM="general"
         ELSE
            SYMM="symmetric"
         END IF
         WRITE(IUNIT,FMT=*)'%%MatrixMarket matrix coordinate ',
     &           trim(ARITH)," ",trim(SYMM)
         WRITE(IUNIT,*) id%N, id%N, id%NZ_loc
         IF (associated(id%A_loc)) THEN
            DO I=1,id%NZ_loc
               IF (id%KEEP(50).NE.0 .AND.
     &             id%IRN_loc(I).LT.id%JCN_loc(I)) THEN
                     WRITE(IUNIT,*) id%JCN_loc(I), id%IRN_loc(I),
     &                    dble(id%A_loc(I)), aimag(id%A_loc(I))
               ELSE
                     WRITE(IUNIT,*) id%IRN_loc(I), id%JCN_loc(I),
     &                    dble(id%A_loc(I)), aimag(id%A_loc(I))
               ENDIF
            ENDDO
         ELSE
            DO I=1,id%NZ_loc
               IF (id%KEEP(50).NE.0 .AND. 
     &            id%IRN_loc(I).LT.id%JCN_loc(I)) THEN
C                 permute upper diag entry
                  WRITE(IUNIT,*) id%JCN_loc(I), id%IRN_loc(I)
               ELSE
                  WRITE(IUNIT,*) id%IRN_loc(I), id%JCN_loc(I)
               ENDIF
            ENDDO
         ENDIF
      ENDIF
      RETURN
      END SUBROUTINE ZMUMPS_DUMP_MATRIX
      SUBROUTINE ZMUMPS_DUMP_RHS(IUNIT, id)
C
C  Purpose:
C  =======
C     Dumps a dense, centralized,
C     right-hand side in matrix market format on unit
C     IUNIT. Should be called on the host only.
C
      USE ZMUMPS_STRUC_DEF
      IMPLICIT NONE
C  Arguments
C  =========
      TYPE(ZMUMPS_STRUC), intent(in)  :: id
      INTEGER, intent(in)             :: IUNIT
C
C  Local variables
C  ===============
C
      CHARACTER (LEN=8)    :: ARITH
      INTEGER              :: I, J, K, LD_RHS
C
C  Executable statements
C  =====================
C
      IF (associated(id%RHS)) THEN
               ARITH='complex'
        WRITE(IUNIT,FMT=*)'%%MatrixMarket matrix array ',
     &           trim(ARITH),
     &           ' general'
        WRITE(IUNIT,*) id%N, id%NRHS
        IF ( id%NRHS .EQ. 1 ) THEN
           LD_RHS = id%N
        ELSE
           LD_RHS = id%LRHS
        ENDIF
        DO J = 1, id%NRHS
           DO I = 1, id%N
              K=(J-1)*LD_RHS+I
                 WRITE(IUNIT,*) dble(id%RHS(K)), aimag(id%RHS(K))
        ENDDO
        ENDDO
      ENDIF
      RETURN
      END SUBROUTINE ZMUMPS_DUMP_RHS
      SUBROUTINE ZMUMPS_BUILD_I_AM_CAND( NSLAVES, K79, 
     &     NB_NIV2, MYID_NODES,
     &     CANDIDATES, I_AM_CAND )
      IMPLICIT NONE
C
C     Purpose:
C     =======
C     Given  a list of candidate processors per node,
C     returns an array of booleans telling whether the
C     processor is candidate or not for a given node.
C
C     K79 holds splitting strategy (KEEP(79)). If K79>1 then
C     TPYE4,5,6 nodes might have been introduced and 
C     in this case "hidden" slaves should be taken 
C     into account to enable dynamic redistribution 
C     of the hidden slaves while climbing the chain of 
C     split nodes. The master of the first node in the 
C     chain requires a special treatment and is thus here
C     not considered as a slave. 
C     
      INTEGER, intent(in) :: NSLAVES, NB_NIV2, MYID_NODES, K79
      INTEGER, intent(in) :: CANDIDATES( NSLAVES+1, NB_NIV2 )
      LOGICAL, intent(out):: I_AM_CAND( NB_NIV2 )
      INTEGER I, INIV2, NCAND
      IF (K79.GT.0) THEN
C      Because of potential restarting the number of
C      candidates that will be used to distribute 
C      arrowheads have to include all possible candidates.
       DO INIV2=1, NB_NIV2
         I_AM_CAND(INIV2)=.FALSE.
         NCAND = CANDIDATES(NSLAVES+1,INIV2)
C        check if some hidden slaves are there
C        Note that if hidden candidates exists (type 5 or 6 nodes) then
C        in position CANDIDATES (NCAND+1,INIV2) must be the master 
C        of the first node in the chain (type 4) that we skip here because
C        a special treatment (it has to be "considered as a master" for all 
C        nodes in the list) is needed.
         DO I=1, NSLAVES
            IF (CANDIDATES(I,INIV2).LT.0) EXIT ! end of extra slaves
            IF (I.EQ.NCAND+1) CYCLE 
!     skip master of associated TYPE 4 node 
            IF (CANDIDATES(I,INIV2).EQ.MYID_NODES) THEN
               I_AM_CAND(INIV2)=.TRUE.
               EXIT
            ENDIF
         ENDDO
       END DO
      ELSE
       DO INIV2=1, NB_NIV2
         I_AM_CAND(INIV2)=.FALSE.
         NCAND = CANDIDATES(NSLAVES+1,INIV2)
         DO I=1, NCAND
            IF (CANDIDATES(I,INIV2).EQ.MYID_NODES) THEN
               I_AM_CAND(INIV2)=.TRUE.
               EXIT
            ENDIF
         ENDDO
       END DO
      ENDIF
      RETURN
      END SUBROUTINE ZMUMPS_BUILD_I_AM_CAND