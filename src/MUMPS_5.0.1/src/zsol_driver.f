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
      SUBROUTINE ZMUMPS_SOLVE_DRIVER( id)
      USE ZMUMPS_STRUC_DEF
      USE MUMPS_SOL_ES
C
C  Purpose
C  =======
C
C  Performs solution phase (solve), Iterative Refinements
C  and Error analysis.
C
C
C
C
      USE ZMUMPS_COMM_BUFFER
      USE ZMUMPS_OOC
      USE TOOLS_COMMON
      IMPLICIT NONE
C     -------------------
C     Explicit interfaces
C     -------------------
      INTERFACE
      SUBROUTINE ZMUMPS_SIZE_IN_STRUCT( id, NB_INT,NB_CMPLX )
      USE ZMUMPS_STRUC_DEF
      TYPE (ZMUMPS_STRUC) :: id
      INTEGER(8)        :: NB_INT,NB_CMPLX
      END SUBROUTINE ZMUMPS_SIZE_IN_STRUCT
      SUBROUTINE ZMUMPS_CHECK_DENSE_RHS
     &(idRHS, idINFO, idN, idNRHS, idLRHS)
      COMPLEX(kind=8), DIMENSION(:), POINTER :: idRHS
      INTEGER, intent(in)    :: idN, idNRHS, idLRHS
      INTEGER, intent(inout) :: idINFO(:)
      END SUBROUTINE ZMUMPS_CHECK_DENSE_RHS
      END INTERFACE
C
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_headers.h'
      INCLUDE 'mumps_tags.h'
#if defined(V_T)
      INCLUDE 'VT.inc'
#endif
      INTEGER :: STATUS(MPI_STATUS_SIZE)
      INTEGER :: IERR
      INTEGER, PARAMETER :: MASTER = 0
C
C  Parameters
C  ==========
C
      TYPE (ZMUMPS_STRUC), TARGET :: id
C
C  Local variables
C  ===============
C
      INTEGER MP,LP, MPG
      LOGICAL PROK, PROKG, LPOK
      INTEGER MTYPE, ICNTL21
      LOGICAL LSCAL, ERANA, GIVSOL
      INTEGER ICNTL10, ICNTL11
      INTEGER I,K,JPERM, J, II, IZ2
      INTEGER IZ, NZ_THIS_BLOCK
C     pointers in IS
      INTEGER LIW
C     pointers in id%S
      INTEGER(8) :: LA, LA_PASSED
      INTEGER LIW_PASSED
      INTEGER LWCB_MIN, LWCB, LWCB_SOL_C
      INTEGER(8) :: TMP_LWCB8
C     buffer sizes
      INTEGER ZMUMPS_LBUF, ZMUMPS_LBUF_INT
      INTEGER MSG_MAX_BYTES_SOLVE, MSG_MAX_BYTES_GTHRSOL
C     null space
      INTEGER IBEG_ROOT_DEF, IEND_ROOT_DEF,
     &        IBEG_GLOB_DEF, IEND_GLOB_DEF,
     &        IROOT_DEF_RHS_COL1
C
      INTEGER NITREF, NOITER, SOLVET, KASE
      COMPLEX(kind=8) RSOL(1)
C     Meaningful only with tree pruning and sparse RHS
      LOGICAL INTERLEAVE_PAR, DO_PERMUTE_RHS
C     true if ZMUMPS_SOL_C called during postprocessing
      LOGICAL FROM_PP
C
C     TIMINGS
      DOUBLE PRECISION TIMEIT, TIMEEA, TIMEEA1, TIMELCOND
      DOUBLE PRECISION TIME3
      DOUBLE PRECISION TIMEC1,TIMEC2
      DOUBLE PRECISION TIMEGATHER1,TIMEGATHER2
      DOUBLE PRECISION TIMESCATTER1,TIMESCATTER2
C     ------------------------------------------
C     Declarations related to exploit sparsity
C     ------------------------------------------
      INTEGER     :: NRHS_NONEMPTY
      INTEGER     :: STRAT_PERMAM1
      LOGICAL     :: DO_NULL_PIV
      INTEGER, DIMENSION(:), POINTER :: IRHS_PTR_COPY
      INTEGER, DIMENSION(:), POINTER :: IRHS_SPARSE_COPY
      COMPLEX(kind=8), DIMENSION(:), POINTER :: RHS_SPARSE_COPY
      LOGICAL IRHS_SPARSE_COPY_ALLOCATED, IRHS_PTR_COPY_ALLOCATED,
     &        RHS_SPARSE_COPY_ALLOCATED
C
      INTEGER         :: NBCOL, COLSIZE, JBEG_RHS, JEND_RHS, JBEG_NEW,
     &                   NBCOL_INBLOC, IPOS, NBT, IPOSRHSCOMP
      INTEGER     :: PERM_PIV_LIST(max(id%KEEP(112),1))
      INTEGER     :: MAP_PIVNUL_LIST(max(id%KEEP(112),1))
      INTEGER, DIMENSION(:), ALLOCATABLE :: PERM_RHS
      INTEGER, DIMENSION(:), POINTER :: PTR_POSINRHSCOMP_FWD,
     &                                  PTR_POSINRHSCOMP_BWD
      COMPLEX(kind=8), DIMENSION(:), POINTER :: PTR_RHS
      INTEGER :: SIZE_IPTR_WORKING, SIZE_WORKING
C     NRHS_NONEMPTY: holds
C         either the original number of RHS (id%NRHS defined on host)
C         or, when the RHS is sparse, it holds the
C         number of non empty columns.
C         it is computed on master and is
C         then broadcasted on all processes.
C     IRHS_PTR_COPY holds a compress local copy of IRHS_PTR (or points
C              on the master to id%IRHS_PTR if no permutation requested)
C     IRHS_SPARSE_COPY might be allocated or might also point to
C         id%IRHS_SPARSE. To test if we can deallocate it we trace
C         with IRHS_SPARSE_COPY_ALLOCATED when it was effectively
C         allocated.
C     NBCOL_INBLOC  total nb columns to process in this block
c     JBEG_RHS global ptr for starting column requested for this block
c     JEND_RHS global ptr for end column_number requested for this block
C     PERM_PIV_LIST permuted array of pivots
C     MAP_PIVNUL_LIST: mapping of permuted list
C     PERM_RHS -- Permutation of RHS computed on master and broadcasted
C         on all procs (of size id%NRHS orginal)
C         PERM_RHS(k) = i means that i is the kth column to be processed
C         Note that PERM_RHS will be used also in case of interleaving
C     ------------------------------------
      COMPLEX(kind=8) ONE
      COMPLEX(kind=8) ZERO
      PARAMETER( ONE = (1.0D0,0.0D0) )
      PARAMETER( ZERO = (0.0D0,0.0D0) )
      DOUBLE PRECISION RZERO, RONE
      PARAMETER( RZERO = 0.0D0, RONE = 1.0D0 )
C
C     RHS_IR is internal to ZMUMPS and used for iterative refinement
C     or the error analysis section. It either points to the user's
C     RHS (on the host when the solution is centralized or the RHS
C     is dense), or is a workarray allocated inside this routine
C     of size N.
      COMPLEX(kind=8), DIMENSION(:), POINTER :: RHS_IR
      COMPLEX(kind=8), DIMENSION(:), POINTER :: WORK_WCB
      COMPLEX(kind=8), DIMENSION(:), POINTER :: PTR_RHS_ROOT
      INTEGER :: LPTR_RHS_ROOT
C
C  Local workarrays that will be dynamically allocated
C
      COMPLEX(kind=8), ALLOCATABLE :: SAVERHS(:), C_RW1(:),
     &                                 C_RW2(:),
     &                                 SRW3(:), C_Y(:),
     &                                 C_W(:)
      COMPLEX(kind=8), ALLOCATABLE :: CWORK(:)
      DOUBLE PRECISION, ALLOCATABLE :: R_Y(:), D(:)
      DOUBLE PRECISION, ALLOCATABLE :: R_W(:)
C     The 2 following workarrays are temporary local
C     arrays only used for distributed matrix input
C     (KEEP(54) .NE. 0).
      DOUBLE PRECISION,    ALLOCATABLE, DIMENSION(:) :: R_LOCWK54
      COMPLEX(kind=8), ALLOCATABLE, DIMENSION(:) :: C_LOCWK54
      INTEGER   :: NBENT_RHSCOMP, NB_FS_IN_RHSCOMP_F,
     &             NB_FS_IN_RHSCOMP_TOT
      INTEGER, DIMENSION(:), ALLOCATABLE :: UNS_PERM_INV
      INTEGER LIWK_SOLVE, LIWCB
      INTEGER, ALLOCATABLE :: IW1(:), IWK_SOLVE(:), IWCB(:)
C
C  Parameters arising from the structure
C
      INTEGER(8)       :: MAXS
      DOUBLE PRECISION, DIMENSION(:), POINTER :: CNTL
      INTEGER, DIMENSION (:), POINTER :: KEEP,ICNTL,INFO
      INTEGER(8), DIMENSION (:), POINTER :: KEEP8
      INTEGER, DIMENSION (:), POINTER :: IS
      DOUBLE PRECISION, DIMENSION(:),POINTER::   RINFOG
C     ===============================================================
C     SCALING issues:
C       When scaling was performed
C       RHS holds the solution of the scaled system
C       The unscaled second member (b0) was given
C       then we have to scale both rhs adn solution:
C        A(sca) = LU  = D1*A*D2 , with D2 = COLSCA
C                                      D1 = ROWSCA
C        --------------
C        CASE OF A X =B
C        --------------
C         (ICNTL(9)=1 or MTYPE=1)
C           A*x0 = b0
C           b(sca) = D1 * b0 = ROWSCA*S(ISTW3)
C           A(sca) [(D2) **(-1)] x0 = b(sca)
C           so the computed solution by Check y0 of LU *y0 = b(sca)
C           is : y0 =[(D2) **(-1)] x0 and so x0= D2*y0 is modified
C        --------------
C        CASE OF AT X =B
C        --------------
C         (ICNTL(9).NE.0 or MTYPE=0)
C           A(sca) = LU  = D1*A*D2
C           AT*x0 = b0 => D2ATD1 D1-1 x0 = D2b0
C           b(sca) = D2 * b0 = COLSCA*S(ISTW3)
C           A(sca)T [(D1) **(-1)] x0 = b(sca)
C           so the computed solution by Check y0 of LU *y0 = b(sca)
C           is : y0 =[(D1) **(-1)] x0 and so x0= D1*y0 is modified
C
C      In case of distributed RHS we need
C      scaling information on each processor
C
          type scaling_data_t
            SEQUENCE
            DOUBLE PRECISION, dimension(:), pointer :: SCALING
            DOUBLE PRECISION, dimension(:), pointer :: SCALING_LOC
          end type scaling_data_t
          type (scaling_data_t) :: scaling_data
C      To scale on the fly during GATHER SOLUTION
          DOUBLE PRECISION, DIMENSION(:), POINTER :: PT_SCALING
          DOUBLE PRECISION, TARGET                :: Dummy_SCAL(1)
C
C  ==================== END OF SCALING related data ================
C
C  Local variables
C
C     Interval associated to the subblocks of RHS a node has to process
      INTEGER, DIMENSION(:), ALLOCATABLE, TARGET :: RHS_BOUNDS
      INTEGER                            :: LPTR_RHS_BOUNDS
      INTEGER, DIMENSION(:), POINTER     :: PTR_RHS_BOUNDS
      LOGICAL  :: DO_NBSPARSE, NBSPARSE_LOC
      DOUBLE PRECISION ARRET
      COMPLEX(kind=8) C_DUMMY(1)
      DOUBLE PRECISION R_DUMMY(1)
      INTEGER IDUMMY(1), JDUMMY(1), KDUMMY(1), LDUMMY(1), MDUMMY(1)
      INTEGER, TARGET :: IDUMMY_TARGET(1)
      COMPLEX(kind=8), TARGET :: CDUMMY_TARGET(1)
      INTEGER JJ
      INTEGER allocok
      INTEGER NBRHS, NBRHS_EFF, BEG_RHS, NB_RHSSKIPPED,
     &        IBEG, LD_RHS, KDEC,
     &        MASTER_ROOT, MASTER_ROOT_IN_COMM
      INTEGER IPT_RHS_ROOT, SIZE_ROOT, LD_REDRHS
      INTEGER IBEG_REDRHS, IBEG_RHSCOMP, LD_RHSCOMP,
     &        LENRHSCOMP, NCOL_RHS_loc
      INTEGER LD_RHS_loc, IBEG_loc
      INTEGER NB_K133, IRANK, TSIZE
      INTEGER KMAX_246_247
      INTEGER IFLAG_IR, IRStep
      LOGICAL TESTConv
      LOGICAL WORKSPACE_MINIMAL_PREFERRED, WK_USER_PROVIDED
      INTEGER(8) NB_BYTES     !size of data allocated during solve
      INTEGER(8) NB_BYTES_MAX !MAX size of data allocated during solve
      INTEGER(8) NB_BYTES_EXTRA !For Step2Node, which may be freed later
      INTEGER(8) NB_INT, NB_CMPLX, K34_8, K35_8, NB_BYTES_ON_ENTRY
      INTEGER(8) K16_8, ITMP8
#if defined(V_T)
C  Vampir
      INTEGER soln_drive_class, glob_comm_ini, perm_scal_ini, soln_dist,
     &        soln_assem, perm_scal_post
#endif
      LOGICAL I_AM_SLAVE, BUILD_POSINRHSCOMP
      LOGICAL WORK_WCB_ALLOCATED, IS_INIT_OOC_DONE
      LOGICAL STOP_AT_NEXT_EMPTY_COL
      INTEGER  MTYPE_LOC
      INTEGER MUMPS_PROCNODE
      EXTERNAL MUMPS_PROCNODE
C
C  First executable statement
C
#if defined(V_T)
      CALL VTCLASSDEF( 'Soln driver',soln_drive_class,IERR)
      CALL VTFUNCDEF( 'glob_comm_ini',soln_drive_class,
     &     glob_comm_ini,IERR)
      CALL VTFUNCDEF( 'perm_scal_ini',soln_drive_class,
     &     perm_scal_ini,IERR)
      CALL VTFUNCDEF( 'soln_dist',soln_drive_class,soln_dist,IERR)
      CALL VTFUNCDEF( 'soln_assem',soln_drive_class,soln_assem,IERR)
      CALL VTFUNCDEF( 'perm_scal_post',soln_drive_class,
     &     perm_scal_post,IERR)
#endif
C     -- The following pointers xxCOPY might be allocated but then
C     -- the associated xxCOPY_ALLOCATED will be set to
C     -- enable deallocation
      IRHS_PTR_COPY => IDUMMY_TARGET
      IRHS_PTR_COPY_ALLOCATED = .FALSE.
      IRHS_SPARSE_COPY => IDUMMY_TARGET
      IRHS_SPARSE_COPY_ALLOCATED=.FALSE.
      RHS_SPARSE_COPY => CDUMMY_TARGET
      RHS_SPARSE_COPY_ALLOCATED=.FALSE.
      NULLIFY(RHS_IR)
      NULLIFY(WORK_WCB)
      IS_INIT_OOC_DONE   = .FALSE.
      WK_USER_PROVIDED   = .FALSE.
      WORK_WCB_ALLOCATED = .FALSE.
      CNTL =>id%CNTL
      KEEP =>id%KEEP
      KEEP8=>id%KEEP8
      IS   =>id%IS
      ICNTL=>id%ICNTL
      INFO =>id%INFO
C      ASPK =>id%A
C      COLSCA =>id%COLSCA
C      ROWSCA =>id%ROWSCA
      RINFOG =>id%RINFOG
      LP  = ICNTL( 1 )
      MP  = ICNTL( 2 )
      MPG = ICNTL( 3 )
      LPOK  = ((LP.GT.0).AND.(id%ICNTL(4).GE.1))
      PROK  = ((MP.GT.0).AND.(id%ICNTL(4).GE.2))
      PROKG = ( MPG .GT. 0 .and. id%MYID .eq. MASTER )
      PROKG = (PROKG.AND.(id%ICNTL(4).GE.2))
      IF (.not.PROK)  MP =0
      IF (.not.PROKG) MPG=0
      IF ( PROK  ) WRITE(MP,100)
      IF ( PROKG ) WRITE(MPG,100)
      NB_BYTES       = 0_8
      NB_BYTES_MAX   = 0_8
      NB_BYTES_EXTRA = 0_8
      K34_8    = int(KEEP(34), 8)
      K35_8    = int(KEEP(35), 8)
      K16_8    = int(KEEP(16), 8)
      NBENT_RHSCOMP = 0
C     Used by DISTRIBUTED_SOLUTION to skip empty columns
C     that are skipped (case of sparse RHS)
      NB_RHSSKIPPED = 0
C     next 4 initialisations needed in case of error
C     to free space allocated
      LSCAL              = .FALSE.
      WORK_WCB_ALLOCATED = .FALSE.
      ICNTL21  = -99998  ! will be bcasted later to slaves
      IBEG_RHSCOMP =-152525  ! Should not be used
      BUILD_POSINRHSCOMP = .TRUE.
C     Not needed anymore (since new version of gather)
C     LD_RHSCOMP = max(KEEP(89),1)  ! at the nb of pivots eliminated on
                                    ! that proc
      LD_RHSCOMP = 1
      NB_FS_IN_RHSCOMP_TOT = KEEP(89)
!     number of FS var of the pruned tree
!     mapped on this proc
      NB_FS_IN_RHSCOMP_F = NB_FS_IN_RHSCOMP_TOT
C
C     Depending on the type of parallelism,
C     the master can have the role of a slave
      I_AM_SLAVE = ( id%MYID .ne. MASTER  .OR.
     &             ( id%MYID .eq. MASTER .AND.
     &               KEEP(46) .eq. 1 ) )
C
C     Compute the number of integers and nb of reals in the structure
      CALL ZMUMPS_SIZE_IN_STRUCT (id, NB_INT,NB_CMPLX  )
      NB_BYTES = NB_BYTES + NB_INT * K34_8 + NB_CMPLX * K35_8
      NB_BYTES_ON_ENTRY = NB_BYTES !used to check alloc/dealloc count ok
      NB_BYTES_MAX = max(NB_BYTES_MAX,NB_BYTES)
C     ======================================
C     BEGIN CHECK KEEP ENTRIES AND INTERFACE
C     ======================================
C     The checks below used to be in ZMUMPS_DRIVER. It is more better
C     to have them here in ZMUMPS_SOL_DRIVER because this enables
C     more flexibility in the management of priorities between various
C     checks.
      IF (id%MYID .EQ. MASTER) THEN
c     subroutine only because called at facto and solve
          CALL ZMUMPS_SET_K221(id)
          id%KEEP(111) = id%ICNTL(25)
C         For the case of ICNTL(20)=1 one could
C         switch off exploit sparsity when RHS is too dense.
          IF (id%ICNTL(20) .EQ. 1) id%KEEP(235) = -1 !automatic
          IF (id%ICNTL(20) .EQ. 2) id%KEEP(235) = 0  !off
          IF (id%ICNTL(20) .EQ. 3) id%KEEP(235) = 1  !on
          IF (id%ICNTL(20).EQ.1.or.id%ICNTL(20).EQ.2.or.
     &       id%ICNTL(20).EQ.3) THEN
           id%KEEP(248) = 1 !sparse RHS
          ELSE
           id%KEEP(248) = 0 !dense RHS
          ENDIF
          ICNTL21      = id%ICNTL(21)
          IF (ICNTL21 .ne.0.and.ICNTL21.ne.1) ICNTL21=0
          IF ( id%ICNTL(30) .NE.0 ) THEN
C           A-1 is on
            id%KEEP(237) = 1
          ELSE
C           A-1 is off
            id%KEEP(237) = 0
          ENDIF
          IF (id%KEEP(248) .eq.0.and. id%KEEP(237).ne.0) THEN
C            For A-1 we have a sparse RHS in the API.
C            Force KEEP(248) accordingly.
             id%KEEP(248)=1
          ENDIF
          IF ((id%KEEP(221).EQ.2 ).AND.(id%KEEP(248).NE.0) ) THEN
C          -- input RHS is in fact effectively
C          -- stored in REDRHS and RHSCOMP
           id%KEEP(248) = 0
          ENDIF
          IF ((id%KEEP(221).EQ.2 ).AND.(id%KEEP(235).NE.0) ) THEN
C          -- input RHS is in fact effectively
C          -- stored in REDRHS and RHSCOMP
           id%KEEP(235) = 0
          ENDIF
          IF ( (id%KEEP(248).EQ.0).AND.(id%KEEP(111).EQ.0) ) THEN
C           RHS is not sparse and thus exploit sparsity is reset to 0
            id%KEEP(235) = 0
          ENDIF
C         Case of Automatic setting of exploit sparsity (KEEP(235)=-1)
C         (in MUMPS_DRIVER original value of KEEP(235) is reset)
          IF(id%KEEP(111).NE.0) id%KEEP(235)=0
C
          IF (id%KEEP(235).EQ.-1) THEN
            IF (id%KEEP(237).NE.0) THEN
C            for A-1
              id%KEEP(235)=1
            ELSE
              id%KEEP(235)=1
            ENDIF
          ELSE IF (id%KEEP(235).NE.0) THEN
            id%KEEP(235)=1
          ENDIF
C         Setting of KEEP(242) (permute RHS)
          IF ((KEEP(111).NE.0)) THEN
C          In the context of null space, the null pivots
C          are by default permuted to post-order
C          However for null space there is in this case no need to
C          permute null pivots since they are already in correct order.
C          Setting KEEP(242)=1 would just force to go through
C          part of the code permuting to identity.
C          Apart for validation purposes this is not interesting
C          costly (and more risky).
              KEEP(242)   = 0
          ENDIF
          IF ( (KEEP(237) .EQ.0) .and. (KEEP(111).EQ.0) ) THEN
C           Permutation possible so far only with A-1 (KEEP(237).NE.0)
              KEEP(242)   = 0
          ENDIF
C          Interleave (243) possible only
C          when permute RHS (242) is on
           IF (KEEP(242).EQ.0) KEEP(243)=0
          IF (id%KEEP(237).EQ.1) THEN
C           Case of automatic setting of KEEP(243), KEEP(493-498)
C           (exploit sparsity parameters)
           IF (id%NSLAVES.EQ.1) THEN
            IF (id%KEEP(243).EQ.-1) id%KEEP(243)=0  
            IF (id%KEEP(495).EQ.-1) id%KEEP(495)=1
            IF (id%KEEP(497).EQ.-1) id%KEEP(497)=1
           ELSE
            IF (id%KEEP(243).EQ.-1) id%KEEP(243)=1
            IF (id%KEEP(495).EQ.-1) id%KEEP(495)=1
            IF (id%KEEP(497).EQ.-1) id%KEEP(497)=1
           ENDIF
          ELSE
            id%KEEP(243)=0
            id%KEEP(495)=0
            id%KEEP(497)=0
          ENDIF
          MTYPE = id%ICNTL(  9 )
          IF (MTYPE.NE.1) MTYPE=0  ! see interface
          IF ((MTYPE.EQ.0).AND.KEEP(50).NE.0) MTYPE =1
!         suppress option Atx=b for A-1
          IF (id%KEEP(237).NE.0) MTYPE = 1
      ENDIF
      CALL MPI_BCAST(MTYPE,1,MPI_INTEGER,MASTER,
     &               id%COMM,IERR)
      CALL MPI_BCAST( id%KEEP(111), 1, MPI_INTEGER, MASTER, id%COMM,
     &                  IERR )
      CALL MPI_BCAST( id%KEEP(221), 1, MPI_INTEGER, MASTER, id%COMM,
     &                  IERR )
      CALL MPI_BCAST( id%KEEP(235), 1, MPI_INTEGER, MASTER, id%COMM,
     &                  IERR )
      CALL MPI_BCAST( id%KEEP(237), 1, MPI_INTEGER, MASTER, id%COMM,
     &                  IERR )
      CALL MPI_BCAST( id%KEEP(242), 2, MPI_INTEGER, MASTER, id%COMM,
     &                  IERR )
      CALL MPI_BCAST( id%KEEP(248), 1, MPI_INTEGER, MASTER, id%COMM,
     &                  IERR )
      CALL MPI_BCAST( id%KEEP(495), 3, MPI_INTEGER, MASTER, id%COMM,
     &                  IERR )
      CALL MPI_BCAST( ICNTL21, 1, MPI_INTEGER, MASTER, id%COMM, IERR )
C     Broadcast original id%NRHS (used at least for checks on ISOL_loc
C                       and to allocate PERM_RHS in case of exploit sparsity)
      CALL MPI_BCAST( id%NRHS,1, MPI_INTEGER, MASTER, id%COMM,IERR)
C
C     TIMINGS: reset to 0
      TIMEC2=0.0D0
      TIMEGATHER2=0.0D0
      TIMESCATTER2=0.0D0
      id%DKEEP(112)=0.0D0
      id%DKEEP(113)=0.0D0
C     id%DKEEP(114) time for the iterative refinement
C     id%DKEEP(120) time for the error analysis
C     id%DKEEP(121) time for condition number
      id%DKEEP(114)=0.0D0
      id%DKEEP(120)=0.0D0
      id%DKEEP(121)=0.0D0
      id%DKEEP(115)=0.0D0
      id%DKEEP(116)=0.0D0
C     Time for fwd, bwd and sclapack is
C     accumulated in DKEEP(117-119) within SOL_C
C     If requested time for each call to FWD/BWD
C     might be print but on output to solve
C     phase DKEEP will hold on each proc the accumulated time
      id%DKEEP(117)=0.0D0
      id%DKEEP(118)=0.0D0
      id%DKEEP(119)=0.0D0
      CALL MUMPS_SECDEB(TIME3)
C     ------------------------------
C     Check parameters on the master
C     ------------------------------
      IF ( id%MYID .EQ. MASTER ) THEN
          IF ((KEEP(23).NE.0).AND.KEEP(50).NE.0) THEN
C          Maxiumum transversal permutation
C          has not been saved (KEEP(23)>0 and UNS_PERM allocated)
C          when matrix is symmetric.
           IF (PROKG) WRITE(MPG,'(A)')
     &       ' Internal Error 1 in solution driver '
           id%INFO(1)=-444
           id%INFO(2)=KEEP(23)
          ENDIF
C         ------------------------------------
C         Check that factors are available
C         either in-core or on disk, case
C         where factors were discarded during
C         factorization (e.g. useful to simulate
C         an OOC factorization or just get nb of
C         negative pivots or determinant)
C         ------------------------------------
          IF (KEEP(201) .EQ. -1) THEN
             IF (PROKG) WRITE(MPG,'(A)')
     &       ' ERROR: Solve impossible because factors not kept'
             id%INFO(1)=-44
             id%INFO(2)=KEEP(251)
             GOTO 333
          ELSE IF (KEEP(221).EQ.0 .AND. KEEP(251) .EQ. 2
     &            .AND. KEEP(252).EQ.0) THEN
             IF (PROKG) WRITE(MPG,'(A)')
     &       ' ERROR: Solve impossible because factors not kept'
             id%INFO(1)=-44
             id%INFO(2)=KEEP(251)
             GOTO 333
          ENDIF
C         ------------------
          IF (KEEP(252).NE.0 .AND. id%NRHS .NE. id%KEEP(253)) THEN
C            Fwd in facto
C            KEEP(252-253) available on all procs since analysis phase
C            Error: id%NRHS is not allowed to change since analysis
C            because fwd has been performed during facto with
C            KEEP(253) RHS
             IF (PROKG) WRITE(MPG,'(A)')
     &       ' ERROR: id%NRHS not allowed to change when ICNTL(32)=1'
             id%INFO(1)=-42
             id%INFO(2)=id%KEEP(253)
             GOTO 333
          ENDIF
C         Testing MTYPE instead of ICNTL(9)
          IF (KEEP(252).NE.0 .AND. MTYPE.NE.1) THEN
C            Fwd in facto is not compatible with transpose system
             INFO(1) = -43
             INFO(2) = 9
             IF (PROKG) WRITE(MPG,'(A)')
     &       ' ERROR: Transpose system (ICNTL(9).NE.0) not ',
     &       ' compatible with forward performed during',
     &       ' factorization (ICNTL(32)=1)'
             GOTO 333
          ENDIF
          IF (KEEP(248) .NE. 0.AND.KEEP(252).NE.0) THEN
C            Fwd during facto incompatible with sparse RHS
C            Forbid sparse RHS when Fwd performed during facto
C            Sparse RHS may be due to A-1 (ICNTL(30)
             INFO(1) = -43
             IF (KEEP(237).NE.0) THEN
              INFO(2) = 30 ! icntl(30)
              IF (PROKG) WRITE(MPG,'(A)')
     &        ' ERROR: A-1 functionality incompatible with forward',
     &        ' performed during factorization (ICNTL(32)=1)'
             ELSE
              INFO(2) = 20 ! ICNTL(20)
              IF (PROKG) WRITE(MPG,'(A)')
     &        ' ERROR: sparse RHS incompatible with forward',
     &        ' performed during factorization (ICNTL(32)=1)'
             ENDIF
             GOTO 333
          ENDIF
          IF (KEEP(237) .NE. 0 .AND. ICNTL21.NE.0) THEN
             IF (PROKG) WRITE(MPG,'(A)')
     &       ' ERROR: A-1 functionality is incompatible',
     &       ' with distributed solution.'
             INFO(1)=-48
             INFO(2)=21
             GOTO 333
          ENDIF
          IF (KEEP(237) .NE. 0 .AND. KEEP(60) .NE.0) THEN
             IF (PROKG) WRITE(MPG,'(A)')
     &       ' ERROR: A-1 functionality is incompatible',
     &       ' with Schur.'
             INFO(1)=-48
             INFO(2)=19
             GOTO 333
          ENDIF
          IF (KEEP(237) .NE. 0 .AND. KEEP(111) .NE.0) THEN
             IF (PROKG) WRITE(MPG,'(A)')
     &       ' ERROR: A-1 functionality is incompatible',
     &       ' with null space.'
             INFO(1)=-48
             INFO(2)=25
             GOTO 333
          ENDIF
          IF (id%NRHS .LE. 0) THEN
             id%INFO(1)=-45
             id%INFO(2)=id%NRHS
             GOTO 333
          ENDIF
C         Entries of A-1 are stored in place of the input sparse RHS
C         thus no need for RHS to be allocated.
          IF ( (id%KEEP(237).EQ.0) ) THEN
             IF ((id%KEEP(248) == 0 .AND.KEEP(221).NE.2)
     &            .OR. ICNTL21==0) THEN
C              RHS must be of size N on the master either to
C              store the dense centralized RHS, either to store
C              the dense centralized solution.
               CALL ZMUMPS_CHECK_DENSE_RHS
     &         (id%RHS,id%INFO,id%N,id%NRHS,id%LRHS)
               IF (id%INFO(1) .LT. 0) GOTO 333
             ENDIF
          ELSE
C           Check that the constraint NRHS=N is respected
C           Check for valid sparse RHS structure done
            IF (id%NRHS .NE. id%N) THEN
              id%INFO(1)=-47
              id%INFO(2)=id%NRHS
              GOTO 333
            ENDIF
          ENDIF
          IF (id%KEEP(248) == 1) THEN
C           ------------------------------------
C           RHS_SPARSE, IRHS_SPARSE and IRHS_PTR
C           must be allocated of adequate size
C           ------------------------------------
            IF (( id%NZ_RHS .LE.0 ).AND.(KEEP(237).NE.0)) THEN
C             At least one entry of A-1 must be requested
              id%INFO(1)=-46
              id%INFO(2)=id%NZ_RHS
              GOTO 333
            ENDIF
            IF (( id%NZ_RHS .LE.0 ).AND.(KEEP(221).EQ.1)) THEN
C             At least one entry of RHS must be nonzero with
c             Schur reduced RHS option
              id%INFO(1)=-46
              id%INFO(2)=id%NZ_RHS
              GOTO 333
            ENDIF
            IF ( .not. associated(id%RHS_SPARSE) )THEN
              id%INFO(1)=-22
              id%INFO(2)=10
              GOTO 333
            ENDIF
            IF ( .not. associated(id%IRHS_SPARSE) )THEN
              id%INFO(1)=-22
              id%INFO(2)=11
              GOTO 333
            ENDIF
            IF ( .not. associated(id%IRHS_PTR) )THEN
              id%INFO(1)=-22
              id%INFO(2)=12
              GOTO 333
            ENDIF
C
            IF (size(id%IRHS_PTR) < id%NRHS + 1) THEN
              id%INFO(1)=-22
              id%INFO(2)=12
              GOTO 333
            END IF
            IF (id%IRHS_PTR(id%NRHS + 1).ne.id%NZ_RHS+1) THEN
              id%INFO(1)=-27
              id%INFO(2)=id%IRHS_PTR(id%NRHS+1)
              GOTO 333
            END IF
C           compare with dble to prevent overflow
            IF (dble(id%N)*dble(id%NRHS).LT.dble(id%NZ_RHS)) THEN
C             Possible in case of dupplicate entries in Sparse RHS
              IF (PROKG) THEN
                 write(MPG,*)
     &              " WARNING: many dupplicate entries in ", 
     &              " sparse RHS  provided by the user ",
     &              " id%NZ_RHS,id%N,id%NRHS =",
     &              id%NZ_RHS,id%N,id%NRHS
              ENDIF
            END IF
            IF (id%IRHS_PTR(1).ne.1) THEN
              id%INFO(1)=-28
              id%INFO(2)=id%IRHS_PTR(1)
              GOTO 333
            END IF
            IF (size(id%IRHS_SPARSE) < id%NZ_RHS) THEN
              id%INFO(1)=-22
              id%INFO(2)=11
              GOTO 333
            END IF
            IF (size(id%RHS_SPARSE) < id%NZ_RHS) THEN
              id%INFO(1)=-22
              id%INFO(2)=10
              GOTO 333
            END IF
          ENDIF
C         --------------------------------
C         Set null space options for solve
C         --------------------------------
          CALL ZMUMPS_GET_NS_OPTIONS_SOLVE(ICNTL(1),KEEP(1),MPG,INFO(1))
          IF (INFO(1) .LT. 0) GOTO 333
          IF (KEEP(111).eq.-1.AND.id%NRHS.NE.KEEP(112)+KEEP(17))THEN
                INFO(1)=-32
                INFO(2)=id%NRHS
                GOTO 333
          ENDIF
          IF (KEEP(111).gt.0 .AND. id%NRHS .NE. 1) THEN
                INFO(1)=-32
                INFO(2)=id%NRHS
                GOTO 333
          ENDIF
          IF (KEEP(248) .NE.0.AND.KEEP(111).NE.0) THEN
C             Ignore sparse RHS in case we compute
C             vectors of the null space (KEEP(111)).NE.0.)
              IF (PROKG) WRITE(MPG,'(A)')
     &        ' ERROR: ICNTL(20) and ICNTL(30) functionalities ',
     &        ' incompatible with null space'
              INFO(1) = -37
              IF (KEEP(237).NE.0) THEN
                 INFO(2) = 30 ! icntl(30)
                 IF (PROKG) WRITE(MPG,'(A)')
     &           ' ERROR: ICNTL(30) functionality ',
     &           ' incompatible with null space'
              ELSE
                 IF (PROKG) WRITE(MPG,'(A)')
     &           ' ERROR: ICNTL(20) functionality ',
     &           ' incompatible with null space'
                 INFO(2) = 20  ! inclt(20)
              ENDIF
              GOTO 333
          ENDIF
          IF (( KEEP(111) .LT. -1 ) .OR.
     &         (KEEP(111).GT.KEEP(112)+KEEP(17)) .OR.
     &         (KEEP(111) .EQ.-1 .AND. KEEP(112)+KEEP(17).EQ.0))
     &         THEN
             INFO(1)=-36
             INFO(2)=KEEP(111)
             GOTO 333
          ENDIF
          IF (KEEP(221).NE.0.AND.KEEP(111).NE.0) THEN
             INFO(1)=-37
             INFO(2)=26
             GOTO 333
          ENDIF
       END IF                   ! MASTER
C     --------------------------------------
C     Check distributed solution vectors
C     --------------------------------------
      IF (ICNTL21==1) THEN
          IF ( id%MYID .ne. MASTER  .OR.
     &       ( id%MYID .eq. MASTER .AND.
     &               id%KEEP(46) .eq. 1 ) ) THEN
C           (I)SOL_loc should be allocated to hold the
C           distributed solution on exit
            IF ( id%LSOL_loc < id%KEEP(89) ) THEN
              id%INFO(1)= -29
              id%INFO(2)= id%LSOL_loc
              GOTO 333
            ENDIF
            IF (id%KEEP(89) .NE. 0) THEN
              IF ( .not. associated(id%ISOL_loc) )THEN
                id%INFO(1)=-22
                id%INFO(2)=13
                GOTO 333
              ENDIF
              IF ( .not. associated(id%SOL_loc) )THEN
                id%INFO(1)=-22
                id%INFO(2)=14
                GOTO 333
              ENDIF
              IF (size(id%ISOL_loc) < id%KEEP(89) ) THEN
                id%INFO(1)=-22
                id%INFO(2)=13
                GOTO 333
              END IF
              IF (size(id%SOL_loc) <
     &              (id%NRHS-1)*id%LSOL_loc+id%KEEP(89)) THEN
                id%INFO(1)=-22
                id%INFO(2)=14
                GOTO 333
              END IF
            ENDIF
          ENDIF
      ENDIF
      IF (id%MYID .NE. MASTER) THEN
          IF (id%KEEP(248) == 1) THEN
C          RHS should NOT be associated
C          if I am not master since it is
C          not even used to store the solution
           IF ( associated( id%RHS ) ) THEN
             id%INFO( 1 ) = -22
             id%INFO( 2 ) = 7
             GOTO 333
           END IF
           IF ( associated( id%RHS_SPARSE ) ) THEN
             id%INFO( 1 ) = -22
             id%INFO( 2 ) = 10
             GOTO 333
           END IF
           IF ( associated( id%IRHS_SPARSE ) ) THEN
             id%INFO( 1 ) = -22
             id%INFO( 2 ) = 11
             GOTO 333
           END IF
           IF ( associated( id%IRHS_PTR ) ) THEN
             id%INFO( 1 ) = -22
             id%INFO( 2 ) = 12
             GOTO 333
           END IF
          END IF
      ENDIF
      IF (id%MYID.EQ.MASTER) THEN
C         Do some checks (REDRHS), depending on KEEP(221)
          CALL ZMUMPS_CHECK_REDRHS(id)
      END IF ! MYID.EQ.MASTER
      IF (id%INFO(1) .LT. 0) GOTO 333
C     -------------------------
C     Propagate possible errors
C     -------------------------
 333  CONTINUE
      CALL MUMPS_PROPINFO( id%ICNTL(1),
     &                      id%INFO(1),
     &                      id%COMM, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) GO TO 90
C     ====================================
C     END CHECK INTERFACE AND KEEP ENTRIES
C     ====================================
C     ====================================
C     Process case of NZ_RHS = 0  with
C     sparse RHS and General Sparse (NOT A-1)
C     -----------------------------------
      IF ((id%KEEP(248).EQ.1).AND.(id%KEEP(237).EQ.0)) THEN
C
       CALL MPI_BCAST(id%NZ_RHS,1,MPI_INTEGER,MASTER,
     &               id%COMM,IERR)
C
       IF (id%NZ_RHS.EQ.0) THEN
C       We reset solution to zero and we return
C       (first freeing working space at label 90)
        IF ((ICNTL21.EQ.1).AND.(I_AM_SLAVE)) THEN
C       ----------------------
C       SOL_loc reset to zero
C       ----------------------
C         ----------------------
C         Prepare ISOL_loc array
C         ----------------------
          LIW_PASSED=max(1,KEEP(32))
C         Only called if more than 1 pivot
C         was eliminated by the processor.
C         Note that LSOL_loc >= KEEP(89)
          IF (KEEP(89) .GT. 0) THEN
            CALL ZMUMPS_DISTSOL_INDICES( MTYPE, id%ISOL_loc(1),
     &               id%PTLUST_S(1),
     &               id%KEEP(1),id%KEEP8(1),
     &               id%IS(1), LIW_PASSED,id%MYID_NODES,
     &               id%N, id%STEP(1), id%PROCNODE_STEPS(1),
     &               id%NSLAVES, scaling_data, LSCAL )
            DO J=1, id%NRHS
              DO I=1, KEEP(89)
                id%SOL_loc((J-1)*id%LSOL_loc + I) =ZERO
              ENDDO
            ENDDO
          ENDIF
        ENDIF
        IF (ICNTL21.NE.1) THEN  ! centralized solution
C       ----------------------------
C       RHS reset to zero on master
C       ----------------------------
         IF (id%MYID.EQ.MASTER) THEN
            DO J=1, id%NRHS
              DO I=1, id%N
                id%RHS((J-1)*id%LRHS + I) =ZERO
              ENDDO
            ENDDO
         ENDIF
        ENDIF
C
C       print solve phase stats if requested
        IF ( PROKG )  THEN
C          write(6,*) " NZ_LOC is zero "
           WRITE( MPG, 150 )
     &        id%NRHS, ICNTL(27), ICNTL(9), ICNTL(10), ICNTL(11),
     &        ICNTL(20), ICNTL(21), ICNTL(30)
           IF (KEEP(221).NE.0) THEN
            WRITE (MPG, 152) KEEP(221)
           ENDIF
           IF (KEEP(252).GT.0) THEN   ! Fwd during facto
            WRITE (MPG, 153) KEEP(252)
           ENDIF
        ENDIF
C
C       --------
        GOTO 90 ! end of solve deallocate what is needed
C       ====================================
C       END CHECK INTERFACE AND KEEP ENTRIES
C       ====================================
       ENDIF  ! test NZ_RHS.EQ.0
C      --------
      ENDIF ! (id%KEEP(248).EQ.1).AND.(id%KEEP(237).EQ.0)
      INTERLEAVE_PAR   =.FALSE.
      DO_PERMUTE_RHS   =.FALSE.
C
      IF ((id%KEEP(235).NE.0).or.(id%KEEP(237).NE.0)) THEN
C        Case of pruned elimination tree or selected entries in A-1
         IF (id%KEEP(237).NE.0.AND.
     &       id%KEEP(248).EQ.0) THEN
C         When A-1 is requested (keep(237).ne.0)
C         sparse RHS has been forced to be on.
          IF (LPOK) THEN
           WRITE(LP,'(A,I4,I4)')
     &     ' Internal Error 2 in solution driver (A-1) ',
     &       id%KEEP(237), id%KEEP(248)
          ENDIF
          CALL MUMPS_ABORT()
         ENDIF
C        NBT is inout in MUMPS_REALLOC and should be initialized.
         NBT = 0
C        -- Allocate Step2node on each proc
         CALL MUMPS_REALLOC(id%Step2node, id%KEEP(28), id%INFO, LP,
     &        FORCE=.TRUE.,
     &        STRING='id%Step2node (Solve)', MEMCNT=NBT, ERRCODE=-13)
         CALL MUMPS_PROPINFO( ICNTL(1), INFO(1),
     &        id%COMM, id%MYID )
         IF ( INFO(1).LT.0 ) RETURN
C        -- build Step2node on each proc;
C        -- this is usefull to have at each step a unique
C        -- representative node (associated with principal variable of
C        -- that node.
         IF (NBT.NE.0) THEN
          ! Step2node was reallocated and needs be recomputed
          DO I=1, id%N
           IF (id%STEP(I).LE.0) CYCLE  ! nonprincipal variables
           id%Step2node(id%STEP(I)) = I
          ENDDO
C        ELSE
C          we reuse Step2node computed in a previous solve phase
C          Step2node is deallocated each time a new analysis is
C          performed or when job=-2 is called
         ENDIF
         NB_BYTES = NB_BYTES + int(NBT,8)*K34_8
         NB_BYTES_MAX = max(NB_BYTES_MAX,NB_BYTES)
         NB_BYTES_EXTRA = NB_BYTES_EXTRA + int(NBT,8) * K34_8
C        Mapping information used during solve. In case of several
C        facto+solve it has to be recomputed.
C        In case of several solves with the same
C        facto, it is not recomputed.
C        It used to compute the interleaving
C        for A-1, and, in dev_version, passed to sol_c to compute
C        some stats
         IF((KEEP(235).NE.0).OR.(KEEP(237).NE.0)) THEN
           IF(.NOT.associated(id%IPTR_WORKING)) THEN
             CALL ZMUMPS_BUILD_MAPPING_INFO(id)
           END IF
         END IF
      ENDIF
C
C     Initialize SIZE_OF_BLOCK from MUMPS_SOL_ES module
      IF ( I_AM_SLAVE )
     &  CALL MUMPS_SOL_ES_INIT(id%OOC_SIZE_OF_BLOCK, id%KEEP(201))
      DO_NULL_PIV = .TRUE.
      NBCOL_INBLOC = -9998
      NZ_THIS_BLOCK= -9998
      JBEG_RHS  = -9998
c
      IF (id%MYID.EQ.MASTER) THEN ! Compute NRHS_NONEMPTY
C
C       -- Sparse RHS does
        IF( KEEP(111)==0 .AND. KEEP(248)==1
#if defined(RHSCOMP_BYROWS)
C        In case of row storage with reduced right hand side, we
C        do not take into account empty columns during forward.
C        Therefore NRHS_NONEMPTY will simply be set to id%NRHS
     &   .AND. KEEP(221) .NE. 1
#endif
     &   ) THEN
C       -- Note that KEEP(111).NE.0 (null space on)
C       -- and KEEP(248).NE.0 will be made incompatible
C       -- When computing entries of A-1 (or SparseRHS only)
           NRHS_NONEMPTY = 0
           DO I=1, id%NRHS
              IF (id%IRHS_PTR(I).LT.id%IRHS_PTR(I+1))
     &             NRHS_NONEMPTY = NRHS_NONEMPTY+1 !ith col in non empty
           ENDDO
           IF (NRHS_NONEMPTY.LE.0) THEN
C           Internal error: tested before in mumps_driver
            IF (LPOK)
     &        WRITE(LP,*) " Internal Error 3 in solution driver ",
     &                    " NRHS_NONEMPTY= ",
     &        NRHS_NONEMPTY
              CALL MUMPS_ABORT()
           ENDIF
        ELSE
           NRHS_NONEMPTY = id%NRHS
        ENDIF
      ENDIF
C     ------------------------------------
C     If there is a special root node,
C     precompute mapping of root's master
C     ------------------------------------
      SIZE_ROOT   = -33333
      IF ( KEEP( 38 ) .ne. 0 ) THEN
            MASTER_ROOT = MUMPS_PROCNODE(
     &                    id%PROCNODE_STEPS(id%STEP( KEEP(38))),
     &                    id%NSLAVES )
            IF (id%MYID_NODES .eq. MASTER_ROOT) THEN
              SIZE_ROOT = id%root%TOT_ROOT_SIZE
            ELSE IF ((id%MYID.EQ.MASTER).AND.KEEP(60).NE.0) THEN
C             SIZE_ROOT also used for KEEP(221).NE.0
              SIZE_ROOT=id%KEEP(116)
            ENDIF
      ELSE IF (KEEP( 20 ) .ne. 0 ) THEN
            MASTER_ROOT = MUMPS_PROCNODE(
     &                    id%PROCNODE_STEPS(id%STEP(KEEP(20))),
     &                    id%NSLAVES )
            IF (id%MYID_NODES .eq. MASTER_ROOT) THEN
              SIZE_ROOT = id%IS(
     &               id%PTLUST_S(id%STEP(KEEP(20)))+KEEP(IXSZ) + 3)
            ELSE IF ((id%MYID.EQ.MASTER).AND.KEEP(60).NE.0) THEN
C             SIZE_ROOT also used for KEEP(221).NE.0
              SIZE_ROOT=id%KEEP(116)
            ENDIF
      ELSE
            MASTER_ROOT = -44444
      END IF
C     --------------
C     Get block size
C     --------------
C     We work on a maximum of NBRHS at a time.
C     The leading dimension of RHS is id%LRHS on the host process
C     and it is set to N on slave processes.
      IF (id%MYID .eq. MASTER) THEN
        KEEP(84) = ICNTL(27)
        IF (KEEP(252).NE.0) THEN
!       Fwd in facto: all rhs (KEEP(253) need be processed in one pass
          NBRHS = KEEP(253)
        ELSE
          IF (KEEP(201) .EQ. 0 .OR. KEEP(84) .GT. 0) THEN
            NBRHS = abs(KEEP(84))
          ELSE
            NBRHS = -2*KEEP(84)
          END IF
          IF (NBRHS .GT. NRHS_NONEMPTY ) NBRHS = NRHS_NONEMPTY
C
        ENDIF
C  Avoid to have LENRHSCOMP overflowing.
        IF(huge(LENRHSCOMP)/id%N.LT.NBRHS) THEN
          WRITE(*,'(A,I6,A)')'Warning: NBRHS = ',NBRHS,
     &   ' might be too large.'
          NBRHS = huge(LENRHSCOMP)/id%N-1 ! -1 to avoid rounding pbs
          WRITE(*,'(A,I6)')'NBRHS reset to ',NBRHS
        END IF
      ENDIF
#if defined(V_T)
      CALL VTBEGIN(glob_comm_ini,IERR)
#endif
C     NRHS_NONEMPTY needed on all procs to allocate RHSCOMP on slaves
      CALL MPI_BCAST(NRHS_NONEMPTY,1,MPI_INTEGER,MASTER,
     &               id%COMM,IERR)
      CALL MPI_BCAST(NBRHS,1,MPI_INTEGER,MASTER,
     &               id%COMM,IERR)
C
      IF (KEEP(201).GT.0) THEN
C         --- id%KEEP(201) indicates if OOC is on (=1) of not (=0)
C         -- 107: number of buffers
C         Define number of types of files (L, possibly U)
          WORKSPACE_MINIMAL_PREFERRED = .FALSE.
          IF (id%MYID .eq. MASTER) THEN
             KEEP(107) = max(0,KEEP(107))
             IF ((KEEP(107).EQ.0).AND.
     &            (KEEP(204).EQ.0).AND.(KEEP(211).NE.1) ) THEN
C             -- default setting for release 4.8
              ! Case of
              !  -Emmergency buffer only and
              !  -Synchronous mode
              !  -NO_O_DIRECT (because of synchronous choice)
              ! THEN
              !   "Basic system-based version"
              !   We can force to allocate S to a minimal
              !   value.
              WORKSPACE_MINIMAL_PREFERRED=.TRUE.
             ENDIF
          ENDIF
          CALL MPI_BCAST( KEEP(107), 1, MPI_INTEGER,
     &                  MASTER, id%COMM, IERR )
          CALL MPI_BCAST( KEEP(204), 1, MPI_INTEGER,
     &                  MASTER, id%COMM, IERR )
          CALL MPI_BCAST( KEEP(208), 2, MPI_INTEGER,
     &                  MASTER, id%COMM, IERR )
          CALL MPI_BCAST( WORKSPACE_MINIMAL_PREFERRED, 1,
     &                  MPI_LOGICAL,
     &                  MASTER, id%COMM, IERR )
C       --- end of OOC case
      ENDIF
      IF ( I_AM_SLAVE ) THEN
C
C       NB_K133:  Max number of simultaneously processed
C           active fronts.
C         Why more than one active node ?
C         1/ In parallel when we start a level 2 node
C           then we do not know exactly when we will
C           have received all contributions from the
C           slaves.
C           This is very critical in OOC since the
C           size provided to the solve phase is
C           much smaller and since we need
C           to determine the size fo the buffers for IO.
C            We pospone the allocation of the block NFRONT*NB_NRHS
C            and solve the problem.
C
C
C         2/ While processing a node and sending information
C            if we have not enough memory in send buffer
C            then we must receive.
C            We feel that this is not so critical.
C
        NB_K133     = 3
C
C       To this we must add one time KEEP(133) to store
C       the RHS of the root node if the root is local.
C       Furthermore this quantity has to be multiplied by the
C       blocking size in case of multiple RHS.
C
        IF ( KEEP( 38 ) .NE. 0 .OR. KEEP( 20 ) .NE. 0 ) THEN
          IF ( MASTER_ROOT .eq. id%MYID_NODES ) THEN
            IF (
     &          .NOT. associated(id%root%RHS_CNTR_MASTER_ROOT)
     &         ) THEN
                NB_K133 = NB_K133 + 1
            ENDIF
          END IF
        ENDIF
        LWCB_MIN = NB_K133*KEEP(133)*NBRHS
C
C       ---------------------------------------------------------------
C       Set WK_USER_PROVIDED to true when workspace WK_USER is provided
C       by user
C       We can accept WK_USER to be provided on only one proc and
C       different values of WK_USER per processor. Note that we are
C       inside a block "IF (I_AM_SLAVE)"
        WK_USER_PROVIDED = (id%LWK_USER.NE.0)
        IF (id%LWK_USER.EQ.0) THEN
          ITMP8 = 0_8
        ELSE IF (id%LWK_USER.GT.0) THEN
          ITMP8= int(id%LWK_USER,8)
        ELSE
          ITMP8 = -int(id%LWK_USER,8)* 1000000_8
        ENDIF
C       Incore: Check if the provided size is equal to that used during
C               facto (case of ITMP8/=0 and KEEP8(24)/=ITMP8)
C               But also check case of space not provided during solve
C               but was provided during facto
C                (case of ITMP8=0 and KEEP8(24)/=0)
        IF (KEEP(201).EQ.0) THEN  ! incore
C         Compare provided size with previous size
          IF (ITMP8.NE.KEEP8(24)) THEN
C           -- error when reusing space allocated
            INFO(1) = -41
            INFO(2) = id%LWK_USER
            GOTO 99    ! jump to propinfo
                       ! (S is used in between and not allocated)
                       ! NO COMM must occur then before next propinfo
                       ! it happens in Mila's code but only with
                       ! KEEP(209) > 0
          ENDIF
        ELSE
          KEEP8(24)=ITMP8
        ENDIF
C       KEEP8(24) holds the size of WK_USER provided by user.
C
        MAXS = 0_8
        IF (WK_USER_PROVIDED) THEN
           MAXS = KEEP8(24)
           IF (MAXS.LT. KEEP8(20)) THEN
                  INFO(1)= -11
                  ! MAXS should be increased by at least ITMP8
                  ITMP8  = KEEP8(20)+1_8-MAXS
                  CALL  MUMPS_SET_IERROR(ITMP8, INFO(2))
           ENDIF
           IF (INFO(1) .GE. 0 ) id%S => id%WK_USER(1:KEEP8(24))
        ELSE IF (associated(id%S)) THEN
C          Avoid the use of "size(id%S)" because it returns
C          a default integer that may overflow. Also "size(id%S,kind=8)"
C          will only be available with Fortran 2003 compilers.
           MAXS = KEEP8(23)
        ELSE
          ! S not allocated and WK_USER not provided ==> must be in OOC
          IF (KEEP(201).EQ.0) THEN  ! incore
            WRITE(*,*) ' Working array S not allocated ',
     &                ' on entry to solve phase (in core) '
            CALL MUMPS_ABORT()
          ELSE
C         -- OOC and WK_USER not provided:
C            define size (S) and allocate it
C           ---- modify size of MAXS: in a simple
C           ---- system-based version, we want to
C           ---- use a small size for MAXS, to
C           ---- avoid the system pagecache to be
C           ---- polluted by 'our memory'
C
            IF ( KEEP(209).EQ.-1 .AND. WORKSPACE_MINIMAL_PREFERRED)
     &        THEN
C             We need space to load at least the largest factor
              MAXS = KEEP8(20) + 1_8
            ELSE IF ( KEEP(209) .GE.0 ) THEN
C             Use suggested value of MAXS provided in KEEP(209)
              MAXS = max(int(KEEP(209),8), KEEP8(20) + 1_8)
            ELSE
              MAXS  = id%KEEP8(14) ! initial value: do not use more than
                               ! minimum (non relaxed) size of OOC facto
            ENDIF
            ALLOCATE (id%S(MAXS), stat = allocok)
            KEEP8(23)=MAXS
            IF ( allocok .GT. 0 ) THEN
              IF (LPOK) THEN
              WRITE(LP,*) id%MYID,': problem allocation of S at solve'
              ENDIF
              INFO(1) = -13
              CALL MUMPS_SET_IERROR(MAXS, INFO(2))
              NULLIFY(id%S)
              KEEP8(23)=0_8
            ENDIF
            NB_BYTES = NB_BYTES + KEEP8(23) * K35_8
            NB_BYTES_MAX = max(NB_BYTES_MAX,NB_BYTES)
C           --- end of OOC case
          ENDIF
C         -- end of id%S already associated
        ENDIF
C
C       On the slaves, S is divided as follows:
C       S(1..LA) holds the factors,
C       S(LA+1..MAXS) is free workspace
        IF(KEEP(201).EQ.0)THEN
           LA  = KEEP8(31)
        ELSE
C          MAXS has normally be dimensionned to store only factors.
           LA = MAXS
           IF(MAXS.GT.KEEP8(31)+KEEP8(20)*int(KEEP(107)+1,8))THEN
C            If we have a very large MAXS, the size reserved for
C            loading the factors into memory does not need to exceed the
C            total size of factors. The (KEEP8(20)*(KEEP(107)+1)) term
C            is here in order to ensure that even with round-off
C            problems (linked to the number of solve zones) factors can
C            all be stored in-core
             LA=KEEP8(31)+KEEP8(20)*int(KEEP(107)+1,8)
           ENDIF
        ENDIF
C
C       We need to allocate a workspace of size WCB for the solve phase.
C       Either it is available at the end of MAXS, or we perform a
C       dynamic allocation.
        IF ( MAXS-LA .GT. int(LWCB_MIN,8) ) THEN
C          No need for LWCB to overflow.
C          If too much space is available at the
C          end of MAXS, we limit it to huge(LWCB)
           TMP_LWCB8 = min( MAXS - LA, int(huge(LWCB),8) )
           LWCB      = int( TMP_LWCB8, kind(LWCB) )
           WORK_WCB => id%S(LA+1_8:LA+TMP_LWCB8)
           WORK_WCB_ALLOCATED=.FALSE.
        ELSE
           LWCB = LWCB_MIN
           ALLOCATE(WORK_WCB(LWCB_MIN), stat = allocok)
           IF (allocok < 0 ) THEN
                   INFO(1)=-13
                   INFO(2)=LWCB_MIN
           ENDIF
           WORK_WCB_ALLOCATED=.TRUE.
           NB_BYTES = NB_BYTES + int(size(WORK_WCB),8)*K35_8
           NB_BYTES_MAX = max(NB_BYTES_MAX,NB_BYTES)
        ENDIF
      ENDIF ! I_AM_SLAVE
C -----------------------------------
  99  CONTINUE
      CALL MUMPS_PROPINFO( ICNTL(1), INFO(1),
     &                   id%COMM,id%MYID)
      IF (INFO(1) < 0) GOTO 90
C -----------------------------------
      IF ( I_AM_SLAVE ) THEN
        IF (KEEP(201).GT.0) THEN
          CALL ZMUMPS_INIT_FACT_AREA_SIZE_S(LA)
C         -- This includes thread creation
C         -- for asynchronous strategies
          CALL ZMUMPS_OOC_INIT_SOLVE(id)
          IS_INIT_OOC_DONE = .TRUE.
        ENDIF ! KEEP(201).GT.0
      ENDIF
C
      CALL MUMPS_PROPINFO( ICNTL(1), INFO(1),
     &                     id%COMM,id%MYID)
      IF (INFO(1) < 0) GOTO 90
C
      IF (id%MYID.EQ.MASTER) THEN
        IF ( PROKG )  THEN
          WRITE( MPG, 150 )
     &        id%NRHS, NBRHS, ICNTL(9), ICNTL(10), ICNTL(11),
     &        ICNTL(20), ICNTL(21), ICNTL(30)
          IF (KEEP(111).NE.0) THEN
            WRITE (MPG, 151) KEEP(111)
          ENDIF
          IF (KEEP(221).NE.0) THEN
            WRITE (MPG, 152) KEEP(221)
          ENDIF
          IF (KEEP(252).GT.0) THEN   ! Fwd during facto
            WRITE (MPG, 153) KEEP(252)
          ENDIF
        ENDIF
C
C     ====================================
C       Define LSCAL, ICNTL10 and ICNTL11
C     ====================================
C
        LSCAL = (((KEEP(52) .GT. 0) .AND. (KEEP(52) .LE. 8)) .OR. (
     &    KEEP(52) .EQ. -1) .OR. KEEP(52) .EQ. -2)
        IF ((ICNTL(11) .EQ. 1).OR.(ICNTL(11) .EQ. 2)) THEN
          ICNTL11 = ICNTL(11)
        ELSE
          ICNTL11 = 0
          IF (ICNTL(11) .NE. 0 .AND. PROKG) THEN
           WRITE(MPG, *)
           WRITE(MPG,'(A)')' WARNING: ICNTL(11) treated as if set to 0 '
          ENDIF
        ENDIF
        ICNTL10 = ICNTL(10)
        IF ( (KEEP(55).eq.0) .AND. KEEP(54).eq.0 .AND.
     &      .NOT.associated(id%A) ) THEN
C         matrix is assembled centralized but space was released
C         by the master. Error analysis is suppressed.
          ICNTL10 = 0
          ICNTL11 = 0
        ENDIF
C       FORBID ERROR ANALYSIS AND ITERATIVE REFINEMENT
C       IF WE RETURN A NULL SPACE BASIS or compute entries in A-1
C        of Fwd in facto
C       -When only one columns of A-1 is requested then
C        we could try to reactivate IR even if
C          -code need be updated
C          -accuracy could be # when one or more columns are requested
C
        IF ((KEEP(111).NE.0).OR.(KEEP(237).NE.0).OR.
     &       (KEEP(252).NE.0) ) THEN
C         Forbid error analysis and iterative refinement
          IF (ICNTL10 .GT. 0) THEN
            IF (PROKG) WRITE(MPG,'(A)')
     &    ' WARNING: ICNTL(10) treated as if set to 0 '
          ENDIF
          IF (ICNTL11 .GT. 0) THEN
            IF (PROKG) WRITE(MPG,'(A)')
     &    ' WARNING: ICNTL(11) treated as if set to 0 '
          ENDIF
          ICNTL10 = 0
          ICNTL11 = 0
        END IF
C
C       Forbid error analysis and iterative refinement
C       in case of reduced rhs/solution
C
        IF (KEEP(221).NE.0) THEN
C         Forbid error analysis and error refinement
          IF (ICNTL10 .GT. 0) THEN
            IF (PROKG) WRITE(MPG,'(A)')
     &    ' WARNING: ICNTL(10) treated as if set to 0 (reduced RHS))'
          ENDIF
          IF (ICNTL11 .GT. 0) THEN
            IF (PROKG) WRITE(MPG,'(A)')
     &    ' WARNING: ICNTL(11) treated as if set to 0 (reduced RHS)'
          ENDIF
          ICNTL10 = 0
          ICNTL11 = 0
        END IF
C     =====================
C       Define ERANA
C     =====================
C
C       Define ERANA;
C       ICNTL10: iterative refinements,
C       ICNTL11: error analysis
        ERANA = ((ICNTL11 .EQ. 1) .OR.(ICNTL11 .EQ. 2)
     &                             .OR. (ICNTL10 .NE. 0))
C
C       Forbid error analysis and iterative refinement if
C       the solution is distributed.
C
C       Forbid error analysis in the case where nrhs > 1
C       first version: also forbid iterative refinement.
C
        IF ((ERANA .AND. NBRHS > 1) .OR. ICNTL(21) > 0) THEN
          IF (ICNTL11 > 0) THEN
            IF (PROKG) WRITE(MPG,'(A)')
     &     ' WARNING: ICNTL(11) treated as if set to zero'
            ICNTL11=0
          ENDIF
          IF (ICNTL10 > 0) THEN
            IF (PROKG) WRITE(MPG,'(A)')
     &     ' WARNING: ICNTL(10) treated as if set to zero'
            ICNTL10=0
          ENDIF
          ERANA = .FALSE.
        ENDIF
        IF (ERANA) THEN
            ALLOCATE(SAVERHS(id%N*NBRHS),stat = allocok)
            IF ( allocok .GT. 0 ) THEN
              IF (LPOK)
     &        WRITE(LP,*) id%MYID,
     &        ':Problem in solve: error allocating SAVERHS'
              INFO(1) = -13
              INFO(2) = id%N*NBRHS
              GOTO 111
            END IF
            NB_BYTES = NB_BYTES + int(size(SAVERHS),8)*K35_8
            NB_BYTES_MAX = max(NB_BYTES_MAX,NB_BYTES)
        ENDIF
C
C     Forbid entries in a-1, in case of null space computations
c
        IF (KEEP(237).NE.0 .AND.KEEP(111).NE.0) THEN
C         Ignore ENTRIES IN A-1 in case we compute
C         vectors of the null space (KEEP(111)).NE.0.)
C         We should still allocate IRHS_SPARSE
          IF (PROKG) WRITE(MPG,'(A)')
     &    ' WARNING: KEEP(237) treated as if set to 0 (null space)'
          KEEP(237)=0
        ENDIF
C     -- end of test master
      END IF
C     --------------------------------------------------
C     Broadcast information to have all processes do the
C     same thing (error analysis/iterative refinements/
C     scaling/distribution of solution)
C     --------------------------------------------------
      CALL MPI_BCAST(ICNTL10,1,MPI_INTEGER,MASTER,
     &               id%COMM,IERR)
      CALL MPI_BCAST(ICNTL11,1,MPI_INTEGER,MASTER,
     &               id%COMM,IERR)
      CALL MPI_BCAST(ICNTL21,1,MPI_INTEGER,MASTER,
     &               id%COMM,IERR)
      CALL MPI_BCAST(ERANA,1,MPI_LOGICAL,MASTER,
     &               id%COMM,IERR)
      CALL MPI_BCAST(LSCAL,1,MPI_LOGICAL,MASTER,
     &               id%COMM,IERR)
      CALL MPI_BCAST(KEEP(111),1,MPI_INTEGER,MASTER,
     &               id%COMM,IERR)
C     KEEP(248)==1 if not_NullSpace (KEEP(111)=0)
C           and sparse RHS on input (id%ICNTL(20)/KEEP(248)==1)
C     (KEEP(248)==1 implies KEEP(111) = 0, otherwise error was raised)
C     We cant thus isolate the case of
C     sparse RHS associated to Null space computation because
C     in this case preparation is different since
C        -we skip the forward step and
C        -the pattern of the RHS
C         of the bwd is related to null pivot indices found and not
C         to information contained in the sparse rhs input format.
      DO_PERMUTE_RHS = (KEEP(242).NE.0)
C      apply interleaving in parallel (FOR A-1 or Null space only)
      IF ( (id%NSLAVES.GT.1) .AND. (KEEP(243).NE.0)
     &        )  THEN
C       -- Option to interleave RHS only makes sense when
C       -- A-1 option is on or Null space compution are on
C          (note also that KEEP(243).NE.0 only when PERMUTE_RHS is on)
        IF ((KEEP(237).NE.0).or.(KEEP(111).GT.0)) THEN
           INTERLEAVE_PAR= .TRUE.
        ELSE
          IF (PROKG) THEN
            write(MPG,*) ' Warning incompatible options ',
     &      ' interleave RHS reset to false '
          ENDIF
        ENDIF
      ENDIF
C       --------------------------------------
C       Compute an upperbound of message size
C       for forward and backward solutions:
C       --------------------------------------
        MSG_MAX_BYTES_SOLVE =  ( 4 + KEEP(133) ) * KEEP(34) +
     &                           KEEP(133) * NBRHS * KEEP(35)
     &  + 16 * KEEP(34) !   for request id, pointer to next + safety
C       --------------------------------------
C       Compute an upperbound of message size
C       for ZMUMPS_GATHER_SOLUTION
C       --------------------------------------
C
        IF (KEEP(237).EQ.0) THEN
C         Note that for ZMUMPS_GATHER_SOLUTION LBUFR buffer should
C         be larger that MAX_inode(NPIV))*NBRHS + NPIV
C         which is covered by next formula since KMAX_246_247  is larger
C         than  MAX_inode(NPIV))
C                   2 integers packed (npiv and terminaison)
          KMAX_246_247 = max(KEEP(246),KEEP(247))
          MSG_MAX_BYTES_GTHRSOL =  ( ( 2 + KMAX_246_247 ) * KEEP(34) +
     &                               KMAX_246_247 * NBRHS * KEEP(35) )
        ELSE
C         Each message from a slave is of size max 4:
C            2 integers  : I,J
C            1 complex   : (Aij)⁻1
C            1 terminaison
          MSG_MAX_BYTES_GTHRSOL =  (  3  * KEEP(34) + KEEP(35) )
        ENDIF
C       The buffer is used both for solve and for ZMUMPS_GATHER_SOLUTION
        id%LBUFR_BYTES = max(MSG_MAX_BYTES_SOLVE, MSG_MAX_BYTES_GTHRSOL)
        TSIZE = int(min(100_8*int(MSG_MAX_BYTES_GTHRSOL,8),
     &              10000000_8))
        id%LBUFR_BYTES = max(id%LBUFR_BYTES,TSIZE)
        id%LBUFR = ( id%LBUFR_BYTES + KEEP(34) - 1 ) / KEEP(34)
        IF ( associated (id%BUFR) ) THEN
          NB_BYTES = NB_BYTES - int(size(id%BUFR),8)*K34_8
          DEALLOCATE(id%BUFR)
          NULLIFY(id%BUFR)
        ENDIF
        ALLOCATE (id%BUFR(id%LBUFR),stat=allocok)
        IF ( allocok .GT. 0 ) THEN
            IF (LPOK)
     &      WRITE(LP,*) id%MYID,
     &      ' Problem in solve: error allocating BUFR'
            INFO(1) = -13
            INFO(2) = id%LBUFR
            GOTO 111
        ENDIF
        NB_BYTES = NB_BYTES + int(size(id%BUFR),8)*K34_8
        NB_BYTES_MAX = max(NB_BYTES_MAX,NB_BYTES)
      IF ( I_AM_SLAVE ) THEN
C       ------------------------------------------------------
C       Dimension send buffer for small integers, e.g. TRACINE
C       ------------------------------------------------------
        ZMUMPS_LBUF_INT = ( 20 + id%NSLAVES * id%NSLAVES  * 4 )
     &                 * KEEP(34)
        CALL ZMUMPS_BUF_ALLOC_SMALL_BUF( ZMUMPS_LBUF_INT, IERR )
        IF ( IERR .NE. 0 ) THEN
          INFO(1) = -13
          INFO(2) = ZMUMPS_LBUF_INT
          IF ( LPOK) THEN
            WRITE(LP,*) id%MYID,
     &      ':Error allocating small Send buffer:IERR=',IERR
          END IF
          GOTO 111
        END IF
C
C       ---------------------------------------
C       Dimension cyclic send buffer for normal
C       messages, based on largest message
C       size during forward and backward solves
C       ---------------------------------------
        ZMUMPS_LBUF = (MSG_MAX_BYTES_SOLVE + 2*KEEP(34) )*id%NSLAVES
C       Avoid buffers larger than 100 Mbytes ...
        ZMUMPS_LBUF = min(ZMUMPS_LBUF, 100 000 000)
C       ... as long as we can send messages to at least 3
C       destinations simultaneously
        ZMUMPS_LBUF = max(ZMUMPS_LBUF,
     &      (MSG_MAX_BYTES_SOLVE+2*KEEP(34)) * min(id%NSLAVES,3))
        ZMUMPS_LBUF = ZMUMPS_LBUF + KEEP(34)
        CALL ZMUMPS_BUF_ALLOC_CB( ZMUMPS_LBUF, IERR )
        IF ( IERR .NE. 0 ) THEN
          INFO(1) = -13
          INFO(2) = ZMUMPS_LBUF/KEEP(34) + 1
          IF ( LPOK) THEN
            WRITE(LP,*) id%MYID,
     &      ':Error allocating Send buffer:IERR=', IERR
          END IF
          GOTO 111
        END IF
C
C
C     -- end of I am slave
      ENDIF
C
      IF ( ERANA )  THEN
C       When Iterative refinement of error analysis requested
C       Allocate RHS_IR on slave processors
C       (note that on MASTER RHS_IR points to RHS)
        IF ( id%MYID .NE. MASTER ) THEN
C
          ALLOCATE(RHS_IR(id%N),stat=IERR)
          NB_BYTES = NB_BYTES + int(size(RHS_IR),8)*K35_8
          NB_BYTES_MAX = max(NB_BYTES_MAX,NB_BYTES)
          IF ( IERR .GT. 0 ) THEN
            INFO(1)=-13
            INFO(2)=id%N
            IF (LPOK)
     &        WRITE(LP,*) 'ERROR while allocating RHS on a slave'
            GOTO 111
          END IF
        ELSE
          RHS_IR=>id%RHS
        ENDIF
      ENDIF
C
      CALL MPI_BCAST(KEEP(497),1,MPI_INTEGER,MASTER,
     &               id%COMM,IERR)
C     Parallel A-1 and exploit sparsity between columns
      DO_NBSPARSE = (KEEP(237) .NE. 0).AND.(KEEP(497).NE.0)
      IF ( I_AM_SLAVE ) THEN
        IF(DO_NBSPARSE) THEN
c         --- ALLOCATE outside loop RHS_BOUNDS is needed
          LPTR_RHS_BOUNDS = 2*KEEP(28)
          ALLOCATE(RHS_BOUNDS(LPTR_RHS_BOUNDS), STAT=IERR)
          IF (IERR.GT.0) THEN
            INFO(1)=-13
            INFO(2)=LPTR_RHS_BOUNDS
            IF (LPOK)
     &        WRITE(LP,*) 'ERROR while allocating RHS_BOUNDS on a slave'
              GOTO 111
          END IF
          NB_BYTES = NB_BYTES +
     &        int(size(RHS_BOUNDS),8)*K34_8
          NB_BYTES_MAX = max(NB_BYTES_MAX,NB_BYTES)
          PTR_RHS_BOUNDS => RHS_BOUNDS
        ELSE
          LPTR_RHS_BOUNDS = 1
          PTR_RHS_BOUNDS => IDUMMY_TARGET
        ENDIF
      ENDIF
C     --------------------------------------------------
      IF ( I_AM_SLAVE ) THEN
        IF ((KEEP(221).EQ.2 .AND. KEEP(252).EQ.0)) THEN
C          -- RHSCOMP must have been allocated in
C          -- previous solve step (with option KEEP(221)=1)
           IF (.NOT.associated(id%RHSCOMP)) THEN
             INFO(1) = -35
             INFO(2) = 1
             GOTO 111
           ENDIF
C          IF ((KEEP(248).EQ.0) .OR. (id%NRHS.EQ.1)) THEN
C          POSINRHSCOMP_ROW/COL are meaningful and could even be reused
           IF (.NOT.associated(id%POSINRHSCOMP_ROW) ) ! .OR.
!    &        .NOT.(id%POSINRHSCOMP_COL_ALLOC))
     &     THEN
             INFO(1) = -35
             INFO(2) = 2
             GOTO 111
           ENDIF
           LENRHSCOMP = size(id%RHSCOMP)
           IF (.not.id%POSINRHSCOMP_COL_ALLOC) THEN
C             POSINRHSCOMP_COL that is kept from
C             previous call to solve must then (already)
C             point to id%POSINRHSCOMP_ROW
              id%POSINRHSCOMP_COL => id%POSINRHSCOMP_ROW
           ENDIF
        ELSE
C         ----------------------
C         Allocate POSINRHSCOMP_ROW/COL
C         ----------------------
C         The size of POSINRHSCOMP arrays
C         does not depend on the block of RHS
C         POSINRHSCOMP_ROW/COL are initialized in the loop of RHS
          IF (associated(id%POSINRHSCOMP_ROW)) THEN
            NB_BYTES = NB_BYTES -
     &          int(size(id%POSINRHSCOMP_ROW),8)*K34_8
            DEALLOCATE(id%POSINRHSCOMP_ROW)
          ENDIF
          ALLOCATE (id%POSINRHSCOMP_ROW(id%N), stat = allocok)
          IF ( allocok .GT. 0 ) THEN
             INFO(1)=-13
             INFO(2)=id%N
             GOTO 111
          END IF
          NB_BYTES = NB_BYTES +
     &          int(size(id%POSINRHSCOMP_ROW),8)*K34_8
          NB_BYTES_MAX = max(NB_BYTES_MAX,NB_BYTES)
          IF (id%POSINRHSCOMP_COL_ALLOC) THEN
            NB_BYTES = NB_BYTES -
     &          int(size(id%POSINRHSCOMP_COL),8)*K34_8
            DEALLOCATE(id%POSINRHSCOMP_COL)
          ENDIF
C
          IF ((KEEP(50).EQ.0).OR.KEEP(237).NE.0) THEN
           ALLOCATE (id%POSINRHSCOMP_COL(id%N), stat = allocok)
           IF ( allocok .GT. 0 ) THEN
             INFO(1)=-13
             INFO(2)=id%N
             GOTO 111
           END IF
           id%POSINRHSCOMP_COL_ALLOC = .TRUE.
           NB_BYTES = NB_BYTES +
     &          int(size(id%POSINRHSCOMP_COL),8)*K34_8
           NB_BYTES_MAX = max(NB_BYTES_MAX,NB_BYTES)
          ELSE
C          Do no allocate POSINRHSCOMP_COL
           id%POSINRHSCOMP_COL => id%POSINRHSCOMP_ROW
           id%POSINRHSCOMP_COL_ALLOC = .FALSE.
          ENDIF
          IF (KEEP(221).NE.2) THEN
C         -- only in the case of bwd after reduced RHS
C         -- we have to keep "old" RHSCOMP
            IF (associated(id%RHSCOMP)) THEN
             NB_BYTES = NB_BYTES - int(size(id%RHSCOMP),8)*K35_8
             DEALLOCATE(id%RHSCOMP)
             NULLIFY(id%RHSCOMP)
            ENDIF
          ENDIF
        ENDIF
C       ---------------------------
C       Allocate local workspace
C       for the solve (ZMUMPS_SOL_C)
C       ---------------------------
        LIWK_SOLVE = 4 * KEEP(28) + 1
C       KEEP(228)+1 temporary integer positions
C       will be needed in ZMUMPS_SOL_S
        IF (KEEP(201).EQ.1) THEN
          LIWK_SOLVE = LIWK_SOLVE + KEEP(228) + 1
        ELSE
C         Reserve 1 position to pass array of size 1 in routines
          LIWK_SOLVE = LIWK_SOLVE + 1
        ENDIF
        ALLOCATE ( IWK_SOLVE( LIWK_SOLVE), stat = allocok )
        IF (allocok .GT. 0 ) THEN
         INFO(1)=-13
         INFO(2)=LIWK_SOLVE
         GOTO 111
        END IF
        NB_BYTES = NB_BYTES + int(LIWK_SOLVE,8)*K34_8
        NB_BYTES_MAX = max(NB_BYTES_MAX,NB_BYTES)
C       array IWCB used temporarily to hold
C       indices of a front unpacked from a message
C       and to stack (potentially in a recursive call)
C       headers of size 2 positions of CB blocks.
        LIWCB =  20*NB_K133*2 + KEEP(133)
        ALLOCATE ( IWCB( LIWCB), stat = allocok )
        IF (allocok .GT. 0 ) THEN
         INFO(1)=-13
         INFO(2)=LIWCB
         GOTO 111
        END IF
        NB_BYTES = NB_BYTES + int(LIWCB,8)*K34_8
        NB_BYTES_MAX = max(NB_BYTES_MAX,NB_BYTES)
C
C       -- Code for a slave
C       -----------
C       Subdivision
C       of array IS
C       -----------
        LIW = KEEP(32)
C       Define a work array of size maximum global frontal
C       size (KEEP(133)) for the call to ZMUMPS_SOL_C
C       This used to be of size id%N.
        ALLOCATE(SRW3(KEEP(133)), stat = allocok )
        IF ( allocok .GT. 0 ) THEN
          INFO(1)=-13
          INFO(2)=KEEP(133)
          GOTO 111
        END IF
        NB_BYTES = NB_BYTES + int(size(SRW3),8)*K35_8
        NB_BYTES_MAX = max(NB_BYTES_MAX,NB_BYTES)
C     -----------------
C     End of slave code
C     -----------------
      ELSE
C       I am the master with host not working
C
C       LIW is used on master when calling
C       the routine ZMUMPS_GATHER_SOLUTION.
        LIW=0
      END IF
C
C     Precompute inverse of UNS_PERM outside loop
      IF (allocated(UNS_PERM_INV)) DEALLOCATE(UNS_PERM_INV)
      IF ( ( id%MYID .eq. MASTER.AND.(KEEP(23).GT.0) .AND.
     &         (MTYPE .NE. 1).AND.(KEEP(248).NE.0)
     &       )
C          Permute UNS_PERM on master only with
C          sparse RHS (KEEP(248).NE.0 ) when AT x = b is solved
     &      .OR. ( ( KEEP(237).NE.0 ) .AND. (KEEP(23).NE.0)  )
C          When A-1 is active and when the matrix is unsymmetric
C          and a column permutation has been applied (Max transversal)
C          then  we have performed a
C          factorization of a column permuted matrix AQ = LU.
C          In this case,
C          the permuted entry need be used to select the target
C          entries for the BWD (note that a diagonal entry of A-1
C          is not anymore a diagonal of AQ. Thus a diagonal
C          of A-1 does not correspond to the same path
C          in the tree during FWD and BWD steps when MAXTRANS is on
C          and permutation is not identity.)
C          Note that the inverse permutation
C          UNS_PERM_INV need be allocated on each proc
C          since it is used in ZMUMPS_SOL_C routine for pruning.
C          It is allocated only once and its allocation has been
C          migrated outside the blocking on the right hand sides.
     &      )  THEN
           ALLOCATE(UNS_PERM_INV(id%N),stat=allocok)
           if (allocok .GT.0 ) THEN
             INFO(1)=-13
             INFO(2)=id%N
             GOTO 111
           endif
           NB_BYTES = NB_BYTES + int(id%N,8)*K34_8
           NB_BYTES_MAX = max(NB_BYTES_MAX,NB_BYTES)
           IF (id%MYID.EQ.MASTER) THEN
C           Build inverse permutation
            DO I = 1, id%N
              UNS_PERM_INV(id%UNS_PERM(I))=I
            ENDDO
           ENDIF
C
      ELSE
           ALLOCATE(UNS_PERM_INV(1), stat=allocok)
           if (allocok .GT.0 ) THEN
             INFO(1)=-13
             INFO(2)=1
             GOTO 111
           endif
           NB_BYTES = NB_BYTES + 1_8*K34_8
           NB_BYTES_MAX = max(NB_BYTES_MAX,NB_BYTES)
      ENDIF
C
 111  CONTINUE
#if defined(V_T)
      CALL VTEND(glob_comm_ini,IERR)
#endif
C
C     Synchro point + Broadcast of errors
C
C     Propagate errors
      CALL MUMPS_PROPINFO( ICNTL(1), INFO(1),
     &                   id%COMM,id%MYID)
      IF (INFO(1) .LT.0 ) GOTO 90
      IF ( ( KEEP(237).NE.0 ) .AND. (KEEP(23).NE.0)  ) THEN
C       Broadcast UNS_PERM_INV
        CALL MPI_BCAST(UNS_PERM_INV,id%N,MPI_INTEGER,MASTER,
     &                 id%COMM,IERR)
      ENDIF
C     -------------------------------------
C     BEGIN
C     Preparation for distributed solution
C     -------------------------------------
      IF ( ICNTL21==1 ) THEN
        IF (LSCAL) THEN
C         In case of scaling we will need to scale
C         back the RHS. Put the values of the scaling
C         arrays needed to do that on each processor.
          IF (id%MYID.NE.MASTER) THEN
            IF (MTYPE == 1) THEN
              ALLOCATE(id%COLSCA(id%N),stat=allocok)
            ELSE
              ALLOCATE(id%ROWSCA(id%N),stat=allocok)
            ENDIF
            IF (allocok > 0) THEN
              IF (LPOK) THEN
               WRITE(LP,*) 'Error allocating temporary scaling array'
              ENDIF
              INFO(1)=-13
              INFO(2)=id%N
              GOTO 40
            ENDIF
            NB_BYTES = NB_BYTES + int(id%N,8)*K16_8
            NB_BYTES_MAX = max(NB_BYTES_MAX,NB_BYTES)
          ENDIF
          IF (MTYPE == 1) THEN
              CALL MPI_BCAST(id%COLSCA(1),id%N,
     &                       MPI_DOUBLE_PRECISION,MASTER,
     &                       id%COMM,IERR)
              scaling_data%SCALING=>id%COLSCA
          ELSE
              CALL MPI_BCAST(id%ROWSCA(1),id%N,
     &                       MPI_DOUBLE_PRECISION,MASTER,
     &                       id%COMM,IERR)
              scaling_data%SCALING=>id%ROWSCA
          ENDIF
          IF (I_AM_SLAVE) THEN
            ALLOCATE(scaling_data%SCALING_LOC(id%KEEP(89)),
     &               stat=allocok)
            IF (allocok > 0) THEN
              IF (LPOK) THEN
                WRITE(LP,*) 'Error allocating local scaling array'
              ENDIF
              INFO(1)=-13
              INFO(2)=id%KEEP(89)
              GOTO 40
            ENDIF
            NB_BYTES = NB_BYTES + int(id%KEEP(89),8)*K16_8
            NB_BYTES_MAX = max(NB_BYTES_MAX,NB_BYTES)
          ENDIF
        ENDIF
        IF ( I_AM_SLAVE ) THEN
C         ----------------------
C         Prepare ISOL_loc array
C         ----------------------
          LIW_PASSED=max(1,LIW)
C         Only called if more than 1 pivot
C         was eliminated by the processor.
C         Note that LSOL_loc >= KEEP(89)
          IF (KEEP(89) .GT. 0) THEN
            CALL ZMUMPS_DISTSOL_INDICES( MTYPE, id%ISOL_loc(1),
     &               id%PTLUST_S(1),
     &               id%KEEP(1),id%KEEP8(1),
     &               id%IS(1), LIW_PASSED,id%MYID_NODES,
     &               id%N, id%STEP(1), id%PROCNODE_STEPS(1),
     &               id%NSLAVES, scaling_data, LSCAL )
          ENDIF
          IF (id%MYID.NE.MASTER .AND. LSCAL) THEN
C           ---------------------------------
C           Local (small) scaling arrays have
C           been built, free temporary copies
C           ---------------------------------
            IF (MTYPE == 1) THEN
              DEALLOCATE(id%COLSCA)
              NULLIFY(id%COLSCA)
            ELSE
              DEALLOCATE(id%ROWSCA)
              NULLIFY(id%ROWSCA)
            ENDIF
            NB_BYTES = NB_BYTES - int(id%N,8)*K16_8
          ENDIF
        ENDIF
        IF (KEEP(23) .NE. 0 .AND. MTYPE==1) THEN
C         Broadcast the unsymmetric permutation and
C         permute the indices in ISOL_loc
          IF (id%MYID.NE.MASTER) THEN
            ALLOCATE(id%UNS_PERM(id%N),stat=allocok)
            IF (allocok > 0) THEN
              INFO(1)=-13
              INFO(2)=id%N
              GOTO 40
            ENDIF
          ENDIF
        ENDIF
C
C =====================  ERROR handling and propagation ================
 40     CONTINUE
        CALL MUMPS_PROPINFO( ICNTL(1), INFO(1),
     &                   id%COMM,id%MYID)
        IF (INFO(1) .LT.0 ) GOTO 90
C ======================================================================
C
        IF (KEEP(23) .NE. 0 .AND. MTYPE==1) THEN
          CALL MPI_BCAST(id%UNS_PERM(1),id%N,MPI_INTEGER,MASTER,
     &               id%COMM,IERR)
          IF (I_AM_SLAVE) THEN
            DO I=1, KEEP(89)
              id%ISOL_loc(I) = id%UNS_PERM(id%ISOL_loc(I))
            ENDDO
          ENDIF
          IF (id%MYID.NE.MASTER) THEN
            DEALLOCATE(id%UNS_PERM)
            NULLIFY(id%UNS_PERM)
          ENDIF
        ENDIF
      ENDIF
C     --------------------------------------
C     Preparation for distributed solution
C     END
C     --------------------------------------
C     ----------------------------
C     Preparation for reduced RHS
C     ----------------------------
      IF ( ( KEEP(221) .EQ. 1 ) .OR.
     &     ( KEEP(221) .EQ. 2 )
     &   ) THEN
C       -- First compute MASTER_ROOT_IN_COMM proc number in
C          COMM_NODES on which is mapped the master of the root.
         IF (KEEP(46).EQ.1) THEN
             MASTER_ROOT_IN_COMM=MASTER_ROOT
         ELSE
             MASTER_ROOT_IN_COMM =MASTER_ROOT+1
         ENDIF
         IF ( id%MYID .EQ. MASTER ) THEN
C            --------------------------------
C            Avoid using LREDRHS when id%NRHS is
C            equal to 1, as was done for RHS
C            --------------------------------
             IF (id%NRHS.EQ.1) THEN
               LD_REDRHS = id%KEEP(116)
             ELSE
               LD_REDRHS = id%LREDRHS
             ENDIF
         ENDIF
         IF (MASTER.NE.MASTER_ROOT_IN_COMM) THEN
C        -- Make available LD_REDRHS on MASTER_ROOT_IN_COMM
C           This will then be used to test if a single
C           message can be sent
C           (this is possible if LD_REDRHS=SIZE_SCHUR)
            IF ( id%MYID .EQ. MASTER ) THEN
C            -- send LD_REDRHS to MASTER_ROOT_IN_COMM
C               using COMM communicator
             CALL MPI_SEND(LD_REDRHS,1,MPI_INTEGER,
     &       MASTER_ROOT_IN_COMM, 0, id%COMM,IERR)
            ELSEIF ( id%MYID.EQ.MASTER_ROOT_IN_COMM) THEN
C            -- recv LD_REDRHS
             CALL MPI_RECV(LD_REDRHS,1,MPI_INTEGER,
     &       MASTER, 0, id%COMM,STATUS,IERR)
            ENDIF
C          -- other procs not concerned
         ENDIF
      ENDIF
C
      IF ( KEEP(248)==1 ) THEN  ! Sparse RHS (A-1 or general sparse)
!        JBEG_RHS - current starting column within A-1 or sparse rhs
!                      set in the loop below and used to obtain the
!                      global index of the column of the sparse RHS
!                      Also used to get index in global permutation.
!                      It also allows to skip empty columns;
        JEND_RHS = 0 ! last column in current blockin A-1
C
C       Compute and apply permutations
        IF (DO_PERMUTE_RHS) THEN
C          Allocate PERM_RHS
           ALLOCATE(PERM_RHS(id%NRHS),stat=allocok)
           IF (allocok > 0) THEN
                INFO(1) = -13
                INFO(2) = id%NRHS
                GOTO 109
           ENDIF
           NB_BYTES = NB_BYTES +  int(id%NRHS,8)*K34_8
           NB_BYTES_MAX = max(NB_BYTES_MAX,NB_BYTES)
           IF (id%MYID.EQ.MASTER) THEN
C      PERM_RHS is computed on MASTER, it might be modified in case of
C      interleaving and will thus be distributed (BCAST) to all slaves
C      only later.
C      We compute the permutation of the RHS (sparse matrix)
*     - for A-1 is to apply permutation to IRHS_SPARSE ONLY
*     -               (to compute all inverse entries )
C      Note NRHS_NONEMPTY holds the nb of non empty columns in A-1
C       Compute PERM_RHS
C       on output: PERM_RHS(k) = i means that i is the kth column
              STRAT_PERMAM1 = KEEP(242)
                 CALL MUMPS_PERMUTE_RHS_AM1
     &             (STRAT_PERMAM1, id%SYM_PERM(1),
     &             id%IRHS_PTR(1), id%NRHS+1,
     &             PERM_RHS, id%NRHS,
     &             IERR
     &           )
           ENDIF
        ENDIF
      ENDIF
      IF ( .NOT.((id%KEEP(235).EQ.0).and.(id%KEEP(237).EQ.0)) ) THEN
C        Allocate PERM_RHS of size 1 if not allocated
         IF (.NOT. allocated(PERM_RHS)) THEN
C             Note that within ZMUMPS_SOL_C, PERM_RHS is only
C             used for A-1 case (with DO_PERMUTE_RHS OR INTERLEAVE_RHS
C             being tested) to get the column index for the
C             original matrix of RHS (column index in A-1)
C             of the permuted columns that have been selected.
C             We allocate an array PERM_RHS of size 1
              ALLOCATE(PERM_RHS(1),stat=allocok)
              IF (allocok > 0) THEN
                INFO(1) = -13
                INFO(2) = 1
                GOTO 109
              ENDIF
              NB_BYTES = NB_BYTES +  int(size(PERM_RHS),8)*K34_8
              NB_BYTES_MAX = max(NB_BYTES_MAX,NB_BYTES)
         ENDIF
      ENDIF
C     Propagate errors
109   CALL MUMPS_PROPINFO( ICNTL(1), INFO(1),
     &                   id%COMM,id%MYID)
      IF (INFO(1) .LT.0 ) GOTO 90
c --------------------------
c --------------------------
      IF (id%NSLAVES .EQ. 1) THEN
c     - In case of NS/A-1 we may want to permute RHS
c     - for NS thus is to apply permutation to PIVNUL_LIST
*     - before starting loop of NBRHS
       IF (DO_PERMUTE_RHS .AND. KEEP(111).NE.0 ) THEN
C     NOTE:
C        when host not working both master and slaves have
C        in this case the complete list
            WRITE(*,*) id%MYID, ':INTERNAL ERROR 1 : ',
     &                 ' PERMUTE RHS during null space computation ',
     &                 ' not available yet '
            CALL MUMPS_ABORT()
       ENDIF                     ! End Permute_RHS
      ELSE
C   Phase 1 : ZMUMPS_PERMUTE_RHS_NS
C        local permutations to minimize sequential disk access
C        with  chunck of size KEEP(84)/NSLAVES
C   Phase 2 : ZMUMPS_SOL_APPLY_PARPERM
C        parallel redistribution to exploit // disk access feature
         IF (DO_PERMUTE_RHS .AND. KEEP(111).NE.0 ) THEN
C           Phase 1 to be called on each proc
            WRITE(*,*) id%MYID, ':INTERNAL ERROR 2 : ',
     &                 ' PERMUTE RHS during null space computation ',
     &                 ' not available yet '
            CALL MUMPS_ABORT()
C
C
         ENDIF                  !  End DO_PERMUTE_RHS
         IF (INTERLEAVE_PAR) THEN
          IF ( KEEP(111).NE.0 ) THEN
            WRITE(*,*) id%MYID, ':INTERNAL ERROR 3 : ',
     &        ' INTERLEAVE RHS during null space computation ',
     &        ' not available yet '
            CALL MUMPS_ABORT()
          ELSE
C     -   A-1 + Interleave:
C     permute RHS on master
           IF (id%MYID.EQ.MASTER) THEN
C            -- PERM_RHS must have been already set or initialized
C            -- it is then modified in next routine
             SIZE_WORKING = id%IPTR_WORKING(id%NPROCS+1)-1
             SIZE_IPTR_WORKING = id%NPROCS+1
             CALL MUMPS_INTERLEAVE_RHS_AM1(
     &        PERM_RHS, id%NRHS,
     &        id%IPTR_WORKING(1), SIZE_IPTR_WORKING,
     &        id%WORKING(1), SIZE_WORKING,
     &        id%IRHS_PTR(1),
     &        id%STEP(1), id%SYM_PERM(1), id%N, NBRHS,
     &        id%PROCNODE_STEPS(1), KEEP(28), id%NSLAVES,
     &        KEEP(493).NE.0, 
     &        KEEP(495).NE.0, KEEP(496), PROKG, MPG 
     &        )
           ENDIF               !      End Master
          ENDIF                 !    End A-1 / NS
         ENDIF                  !  End INTERLEAVE_PAR
C -------------
      ENDIF                  ! End Parallel Case
c --------------------------
c
      IF (DO_PERMUTE_RHS.AND.(KEEP(111).EQ.0)) THEN
C     --- Distribute PERM_RHS before loop of RHS
C     --- (with null space option PERM_RHS is not allocated / needed
C          to permute the null column pivot list)
        CALL MPI_BCAST(PERM_RHS(1),
     &            id%NRHS,
     &            MPI_INTEGER,
     &            MASTER, id%COMM,IERR)
      ENDIF
C     ==============================
C     BLOCKING ON the number of RHS
C      We work on  a maximum of NBRHS at a time.
C      the leading dimension of RHS is id%LRHS on master
C      and is set to N on slaves
C     ==============================
C  We may want to allow to have NBRHS that varies
C  this is typically the case when a partitionning of
C  the right hand side is performed and leads to
C  irregular partitions.
C  We only have to be sure that the size of each partition
C  is smaller than NBRHS.
      BEG_RHS=1
      DO WHILE (BEG_RHS.LE.NRHS_NONEMPTY)
C       ==========================
C       -- NBRHS     : Original block size
C       -- BEG_RHS   : Column index of the first RHS in the list of
C                      non empty RHS (RHS_LOC) to
C                      be processed during this iteration
C       -- NBRHS_EFF : Effective block size at current iteration
C          In case of sparse RHS (KEEP(248)==1) NBRHS_EFF only refers to
C                  non-empty columns and is used to compute NBCOL_INBLOC
C          -- NBCOL_INBLOC  : the number of columns of sparse RHS needed
C                        to get NBRHS_EFF non empty columns columns of
C                        sparse RHS processed at each step
C
        NBRHS_EFF    = min(NRHS_NONEMPTY-BEG_RHS+1, NBRHS)
C
C       Sparse RHS
C       Free space and reset pointers if needed
        IF (IRHS_SPARSE_COPY_ALLOCATED) THEN
            NB_BYTES =  NB_BYTES -
     &       int(size(IRHS_SPARSE_COPY),8)*K34_8
            DEALLOCATE(IRHS_SPARSE_COPY)
            IRHS_SPARSE_COPY_ALLOCATED=.FALSE.
            NULLIFY(IRHS_SPARSE_COPY)
        ENDIF
        IF (IRHS_PTR_COPY_ALLOCATED) THEN
            NB_BYTES =  NB_BYTES -
     &       int(size(IRHS_PTR_COPY),8)*K34_8
            DEALLOCATE(IRHS_PTR_COPY)
            IRHS_PTR_COPY_ALLOCATED=.FALSE.
            NULLIFY(IRHS_PTR_COPY)
        ENDIF
        IF (RHS_SPARSE_COPY_ALLOCATED) THEN
            NB_BYTES =  NB_BYTES -
     &       int(size(RHS_SPARSE_COPY),8)*K35_8
            DEALLOCATE(RHS_SPARSE_COPY)
            RHS_SPARSE_COPY_ALLOCATED=.FALSE.
            NULLIFY(RHS_SPARSE_COPY)
        ENDIF
C
C       ===============================================================
C       Set LD_RHS and IBEG
C       IBEG might still be updated in case of general sparse
C        and centralized solution to skip empty columns (see next block)
C       ===============================================================
        IF (
C           slave procs
     &      ( id%MYID .NE. MASTER )
C       even on master when RHS not allocated
     &     .or.
C         Case of Master working but with distributed sol and
C            ( sparse RHS or null space )
C         -- Allocate not needed on host not working
     &    ( I_AM_SLAVE .AND. id%MYID .EQ. MASTER .AND.
     &      ICNTL21 .NE.0 .AND.
     &      ( KEEP(248).ne.0 .OR. KEEP(221).EQ.2
     &          .OR. KEEP(111).NE.0 )
     &    )
     &     .or.
C         Case of Master and
C          (compute entries of INV(A))
C         Even when I am a master with host not working I
C         am in charge of gathering solution to scale it
C         and to copy it back in the sparse RHS format
     &    ( id%MYID .EQ. MASTER .AND. (KEEP(237).NE.0) )
C
     &    ) THEN
          LD_RHS = id%N
          IBEG   = 1
        ELSE
          ! (id%MYID .eq. MASTER)
          IF ( associated(id%RHS) ) THEN
C             Leading dimension of RHS on master is id%LRHS
              LD_RHS    = max(id%LRHS, id%N)
          ELSE
C             --- LRHS might not be defined (dont use it)
              LD_RHS    = id%N
          ENDIF
          IBEG      = (BEG_RHS-1) * LD_RHS + 1
        ENDIF
C       JBEG_RHS might also be used in DISTRIBUTED_SOLUTION
C       even when RHS is not sparse on input.
        JBEG_RHS = BEG_RHS
C       ==========================================
C       Shift empty columns in case of sparse RHS
C       ==========================================
        IF ( (id%MYID.EQ.MASTER) .AND.
     &        KEEP(248)==1  ) THEN
C         update position of JBEG_RHS on first non-empty
C         column of this block
          JBEG_RHS = JEND_RHS + 1
#if defined(RHSCOMP_BYROWS)
C         In case RHSCOMP is stored by rows, we need to ensure
C         that the blocks during forward and backward are the
C         same. For that, a simple and safe solution consists in
C         avoiding skipping empty columns during the forward step.
          IF (KEEP(221).NE.1) THEN
#endif
          IF (DO_PERMUTE_RHS.OR.INTERLEAVE_PAR) THEN
            DO WHILE ( id%IRHS_PTR(PERM_RHS(JBEG_RHS)) .EQ.
     &          id%IRHS_PTR(PERM_RHS(JBEG_RHS)+1) )
              IF ((KEEP(237).EQ.0).AND.(ICNTL21.EQ.0).AND.
     &              (KEEP(221).NE.1) ) THEN
C                General sparse RHS (NOT A-1) and centralized solution
C                Set to zero part of the
C                solution corresponding to empty columns
                 DO I=1, id%N
                   id%RHS((PERM_RHS(JBEG_RHS) -1)*LD_RHS+I)
     &                     = ZERO
                 ENDDO
              ENDIF
              JBEG_RHS = JBEG_RHS +1
              CYCLE
            ENDDO
          ELSE
            DO WHILE( id%IRHS_PTR(JBEG_RHS) .EQ.
     &        id%IRHS_PTR(JBEG_RHS+1) )
              IF ((KEEP(237).EQ.0).AND.(ICNTL21.EQ.0).AND.
     &          (KEEP(221).NE.1) ) THEN
C               Case of general sparse RHS (NOT A-1) and
C               centralized solution: set to zero part of
C               the solution corresponding to empty columns
                DO I=1, id%N
                  id%RHS((JBEG_RHS -1)*LD_RHS + I) = ZERO
                ENDDO
              ENDIF
              IF (KEEP(221).EQ.1) THEN
C               Reduced RHS set to ZERO
                DO I = 1, id%SIZE_SCHUR
                  id%REDRHS((JBEG_RHS-1)*LD_REDRHS + I) =  ZERO
                ENDDO
              ENDIF
              JBEG_RHS = JBEG_RHS +1
            ENDDO
          ENDIF                 ! End DO_PERMUTE_RHS.OR.INTERLEAVE_PAR
#if defined(RHSCOMP_BYROWS)
          ENDIF
C         So we will have NBRHS_SKIPPED = 0
C         And we have JBEG_RHS = JEND_RHS+1
          IF (KEEP(221).EQ.1) THEN
            IF ( id%IRHS_PTR(JBEG_RHS) .EQ.
     &          id%IRHS_PTR(JBEG_RHS+1) )  THEN
              DO J=JBEG_RHS, JBEG_RHS + NBRHS_EFF -1
                DO I=1, id%SIZE_SCHUR
                  id%REDRHS((J-1)*LD_REDRHS + I) =  ZERO
                ENDDO
              ENDDO
            ENDIF
          ENDIF
#endif
C         Count nb of RHS columns skipped: useful for
C         * ZMUMPS_DISTRIBUTED_SOLUTION to reset those
C           columns to zero.
C         * in case of reduced right-hand side, to set
C           corresponding entries of RHSCOMP to 0 after
C           forward phase.
          NB_RHSSKIPPED = JBEG_RHS - (JEND_RHS + 1)
          IF ((KEEP(248).EQ.1).AND.(KEEP(237).EQ.0)
     &         .AND. (ICNTL21.EQ.0))
     &         THEN
             ! case of general sparse rhs with centralized solution,
             !set IBEG to shifted columns
             ! (after empty columns have been skipped)
             IBEG      = (JBEG_RHS-1) * LD_RHS + 1
          ENDIF
        ENDIF ! of if (id%MYID.EQ.MASTER) .AND.  KEEP(248)==1
        CALL MPI_BCAST( JBEG_RHS,1, MPI_INTEGER,
     &            MASTER, id%COMM,IERR)
C
C       Shift on REDRHS in reduced RHS functionality
C
        IF (id%MYID.EQ.MASTER .AND. KEEP(221).NE.0) THEN
C         Initialize IBEG_REDRHS
C         Note that REDRHS always has id%NRHS Colmuns
          IBEG_REDRHS= (JBEG_RHS-1)*LD_REDRHS + 1
        ELSE
          IBEG_REDRHS=-142424  ! Should not be used
        ENDIF
C
C       =====================
C       BEGIN
C       Prepare RHS on master
C
#if defined(V_T)
        CALL VTBEGIN(perm_scal_ini,IERR)
#endif
        IF (id%MYID .eq. MASTER) THEN
C         ======================
          IF (KEEP(248)==1) THEN
C         ======================
C
C         Sparse RHS format ( A-1 or sparse input format)
C         is provided as input by the user (IRHS_SPARSE ...)
C         --------------------------------------------------
C         Compute NZ_THIS_BLOCK and NBCOL_INBLOC
C         where
C         NZ_THIS_BLOCK is defined
C         as the number of entries in the next NBRHS_EFF
C         non empty columns (note that since they might be permuted
C         then the following formula is not always valid:
C            NZ_THIS_BLOCK=id%IRHS_PTR(BEG_RHS+NBRHS_EFF)-
C     &                    id%IRHS_PTR(BEG_RHS)
C         anyway NBCOL_INBLOC also need be computed so going through
C         columns one at a time is needed.
C
          NBCOL        = 0
          NBCOL_INBLOC = 0
          NZ_THIS_BLOCK = 0
C         With exploit sparsity we skip empty rows up to reaching
C         the first non empty column; then we process a block of
C         maximum size NBRHS_EFF except if we reach another empty
C         column. (We are not sure to have a copy allocated
C         and thus cannot compress on the fly, as done naturally
C         for A-1).
          STOP_AT_NEXT_EMPTY_COL = .FALSE.
          DO I=JBEG_RHS, id%NRHS
            NBCOL_INBLOC = NBCOL_INBLOC +1
            IF (DO_PERMUTE_RHS.OR.INTERLEAVE_PAR) THEN
C              PERM_RHS(k) = i means that i is the kth
C                            column to be processed
C              PERM_RHS should also be defined for
C                       empty columns i in A-1 (PERM_RHS(K) = i)
               COLSIZE = id%IRHS_PTR(PERM_RHS(I)+1)
     &                    - id%IRHS_PTR(PERM_RHS(I))
            ELSE
               COLSIZE = id%IRHS_PTR(I+1) - id%IRHS_PTR(I)
            ENDIF
            IF ((.NOT.STOP_AT_NEXT_EMPTY_COL).AND.(COLSIZE.GT.0).AND.
     &          (KEEP(237).EQ.0)) THEN
C              -- set STOP_NEXT_EMPTY_COL only for general
C              --  sparse case (not AM-1)
               STOP_AT_NEXT_EMPTY_COL =.TRUE.
            ENDIF
            IF (COLSIZE.GT.0
#if defined(RHSCOMP_BYROWS)
C            In case of forward-only, we do not skip empty RHS.
C            This would cause problems during the backward phase: since
C            each block of RHSCOMP has a row-major storage and inside
C            each block, data is congiguous, blocks must be the same
C            during forward and during backward. Hence NB_RHSSKIPPED
C            will be 0.
C
     &       .OR. KEEP(221) .EQ. 1
#endif
     &      ) THEN
              NBCOL = NBCOL+1
              NZ_THIS_BLOCK = NZ_THIS_BLOCK + COLSIZE
            ELSE IF (STOP_AT_NEXT_EMPTY_COL) THEN
C We have reached an empty column with already selected non empty 
C columns: reduce block size to non empty columns reached so far.
              NBCOL_INBLOC = NBCOL_INBLOC -1
              NBRHS_EFF = NBCOL
              EXIT
            ENDIF
            IF (NBCOL.EQ.NBRHS_EFF) EXIT
          ENDDO
#if defined(RHSCOMP_BYROWS)
          IF (NZ_THIS_BLOCK .eq. 0) THEN
C           Skip block,
C           set REDRHS, RHSCOMP will be set later
            IF (KEEP(221).EQ.1) THEN
              DO J=JBEG_RHS, JBEG_RHS+ NBRHS_EFF  -1
                DO I = 1, id%SIZE_SCHUR
                  id%REDRHS((JBEG_RHS-1)*LD_REDRHS + I) =  ZERO
                ENDDO
              ENDDO
            ELSE
              WRITE(*,*) "Internal error 15 is sol_driver"
              CALL MUMPS_ABORT()
            ENDIF
          ENDIF
#else
          IF (NZ_THIS_BLOCK.EQ.0) THEN
           WRITE(*,*) " Internal Error 16 in sol driver NZ_THIS_BLOCK=",
     &               NZ_THIS_BLOCK
           CALL MUMPS_ABORT()
          ENDIF
#endif
C
          IF (NBCOL.NE.NBRHS_EFF.AND. (KEEP(237).NE.0)
     &        .AND.KEEP(221).NE.1) THEN
C           With exploit sparsity for general sparse RHS (Not A-1)
C           we skip empty rows up to reaching
C           the first non empty column; then we process a block of
C           maximum size NBRHS_EFF except if we reach another empty
C           column. (We are not sure to have a copy allocated
C           and thus cannot compress on the fly, as done naturally
C           for A-1). Thus NBCOL might be smaller than NBRHS_EFF
            WRITE(6,*) ' Internal Error 8 in solution driver ',
     &            NBCOL, NBRHS_EFF
              call MUMPS_ABORT()
          ENDIF
C         -------------------------------------------------------------
C
          IF (NZ_THIS_BLOCK .NE. 0) THEN
C           -----------------------------------------------------------
C           We recall that
C           NBCOL_INBLOC is the number of columns of sparse RHS needed 
C           to get NBRHS_EFF non empty columns:
            ALLOCATE(IRHS_PTR_COPY(NBCOL_INBLOC+1),stat=allocok)
            if (allocok .GT.0 ) then
              INFO(1)=-13
              INFO(2)=NBCOL_INBLOC+1
              GOTO 30
            endif
            IRHS_PTR_COPY_ALLOCATED = .TRUE.
            NB_BYTES =  NB_BYTES +  int(NBCOL_INBLOC+1,8)*K34_8
            NB_BYTES_MAX = max(NB_BYTES_MAX,NB_BYTES)
C
            JEND_RHS =JBEG_RHS + NBCOL_INBLOC - 1
C           -----------------------------------------------------------
C           Initialize IRHS_PTR_COPY
C           compute local copy (compressed) of id%IRHS_PTR on Master
            IF (DO_PERMUTE_RHS.OR.INTERLEAVE_PAR) THEN
              IPOS = 1
              J = 0
              DO I=JBEG_RHS, JBEG_RHS + NBCOL_INBLOC -1
                J = J+1
                IRHS_PTR_COPY(J) = IPOS
                COLSIZE = id%IRHS_PTR(PERM_RHS(I)+1)
     &                    - id%IRHS_PTR(PERM_RHS(I))
                IPOS = IPOS + COLSIZE
              ENDDO
            ELSE
              IPOS = 1
              J = 0
              DO I=JBEG_RHS, JBEG_RHS + NBCOL_INBLOC -1
                J = J+1
                IRHS_PTR_COPY(J) = IPOS
                COLSIZE = id%IRHS_PTR(I+1)
     &                     - id%IRHS_PTR(I)
                IPOS = IPOS + COLSIZE
              ENDDO
            ENDIF                 ! End DO_PERMUTE_RHS.OR.INTERLEAVE_PAR
            IRHS_PTR_COPY(NBCOL_INBLOC+1)= IPOS
            IF ( IPOS-1 .NE. NZ_THIS_BLOCK ) THEN
                WRITE(*,*) "Error in compressed copy of IRHS_PTR"
                IERR = 99
                call MUMPS_ABORT()
            ENDIF
C           -----------------------------------------------------------
C           IRHS_SPARSE : do a copy or point to the original indices
C
C           Check whether IRHS_SPARSE_COPY need be allocated
            IF (KEEP(23) .NE. 0 .and. MTYPE .NE. 1) THEN
C             AP = LU and At x = b ==> b need be permuted
              ALLOCATE(IRHS_SPARSE_COPY(NZ_THIS_BLOCK)
     &               ,stat=allocok)
              if (allocok .GT.0 ) then
                INFO(1)=-13
                INFO(2)=NZ_THIS_BLOCK
                GOTO 30
              endif
              IRHS_SPARSE_COPY_ALLOCATED=.TRUE.
              NB_BYTES = NB_BYTES + int(NZ_THIS_BLOCK,8)*K34_8
              NB_BYTES_MAX = max(NB_BYTES_MAX,NB_BYTES)
            ELSE IF (DO_PERMUTE_RHS.OR.INTERLEAVE_PAR.OR.
     &           (KEEP(237).NE.0)) THEN
C               Columns are not contiguous and need be copied one by one
C               IRHS_SPARSE_COPY will hold a copy of contiguous permuted
C               columns so an explicit copy is needed.
C               IRHS_SPARSE_COPY is also allways allocated with A-1,
C               to enable receiving during mumps_gather_solution
C     .         on the master in any order.
                ALLOCATE(IRHS_SPARSE_COPY(NZ_THIS_BLOCK),
     &                   stat=allocok)
                IF (allocok .GT.0 ) THEN
                    IERR = 99
                    GOTO 30
                ENDIF
                IRHS_SPARSE_COPY_ALLOCATED=.TRUE.
                NB_BYTES = NB_BYTES + int(NZ_THIS_BLOCK,8)*K34_8
                NB_BYTES_MAX = max(NB_BYTES_MAX,NB_BYTES)
C
            ENDIF
C          
C           Initialize IRHS_SPARSE_COPY 
            IF (IRHS_SPARSE_COPY_ALLOCATED) THEN
                IF ( DO_PERMUTE_RHS.OR.INTERLEAVE_PAR ) THEN
                  IPOS = 1
                  DO I=JBEG_RHS, JBEG_RHS + NBCOL_INBLOC -1
                    COLSIZE = id%IRHS_PTR(PERM_RHS(I)+1)
     &                       - id%IRHS_PTR(PERM_RHS(I))
                    IRHS_SPARSE_COPY(IPOS:IPOS+COLSIZE-1) =
     &              id%IRHS_SPARSE(id%IRHS_PTR(PERM_RHS(I)):
     &              id%IRHS_PTR(PERM_RHS(I)+1) -1)
                    IPOS = IPOS + COLSIZE
                  ENDDO
                ELSE
                  IRHS_SPARSE_COPY = id%IRHS_SPARSE(
     &              id%IRHS_PTR(JBEG_RHS):
     &              id%IRHS_PTR(JBEG_RHS)+NZ_THIS_BLOCK-1)
                ENDIF
            ELSE
                IRHS_SPARSE_COPY
c    *                   (1:NZ_THIS_BLOCK)
     &           =>
     &          id%IRHS_SPARSE(id%IRHS_PTR(JBEG_RHS):
     &             id%IRHS_PTR(JBEG_RHS)+NZ_THIS_BLOCK-1)
            ENDIF
            IF (LSCAL.OR.DO_PERMUTE_RHS.OR.INTERLEAVE_PAR.OR.
     &          (KEEP(237).NE.0)) THEN
C             if scaling is on or if columns of the RHS are
C             permuted then a copy of RHS_SPARSE is needed.
C             Also allways allocated with A-1,
c             to enable receiving during mumps_gather_solution
C             on the master in any order.
C
              ALLOCATE(RHS_SPARSE_COPY(NZ_THIS_BLOCK),
     &               stat=allocok)
              IF (allocok .GT.0 ) THEN
                INFO(1)=-13
                INFO(2)=NZ_THIS_BLOCK
                GOTO 30
              ENDIF
              RHS_SPARSE_COPY_ALLOCATED = .TRUE.
              NB_BYTES = NB_BYTES + int(NZ_THIS_BLOCK,8)*K35_8
              NB_BYTES_MAX = max(NB_BYTES_MAX,NB_BYTES)
            ELSE
              IF  ( KEEP(248)==1 ) THEN
                RHS_SPARSE_COPY
c    *            (1:NZ_THIS_BLOCK)
     &          => id%RHS_SPARSE(id%IRHS_PTR(JBEG_RHS):
     &             id%IRHS_PTR(JBEG_RHS)+NZ_THIS_BLOCK-1)
              ELSE
                RHS_SPARSE_COPY
c    *            (1:NZ_THIS_BLOCK)
     &          => id%RHS_SPARSE(id%IRHS_PTR(BEG_RHS):
     &             id%IRHS_PTR(BEG_RHS)+NZ_THIS_BLOCK-1)
              ENDIF
            ENDIF
            IF (DO_PERMUTE_RHS.OR.INTERLEAVE_PAR.OR.
     &          (id%KEEP(237).NE.0)) THEN
              IF (id%KEEP(237).NE.0) THEN
C               --initialized to one; it might be
C               modified if scaling is on (one first entry in each col is scaled)
                RHS_SPARSE_COPY = ONE
              ELSE IF (.NOT. LSCAL) THEN
C               -- Columns are not contiguous and need be copied one by one
C               -- This need not be done if scaling is on because it
C               -- will done and scaled later.
                IPOS = 1
                DO I=JBEG_RHS, JBEG_RHS + NBCOL_INBLOC -1
                  COLSIZE = id%IRHS_PTR(PERM_RHS(I)+1)
     &                    - id%IRHS_PTR(PERM_RHS(I))
                  IF (COLSIZE .EQ. 0) CYCLE
                  RHS_SPARSE_COPY(IPOS:IPOS+COLSIZE-1) =
     &            id%RHS_SPARSE(id%IRHS_PTR(PERM_RHS(I)):
     &            id%IRHS_PTR(PERM_RHS(I)+1) -1)
                  IPOS = IPOS + COLSIZE
                ENDDO
              ENDIF
            ENDIF
C           =========================
            IF (KEEP(23) .NE. 0) THEN
C           =========================
*             maximum transversal was performed
              IF (MTYPE .NE. 1) THEN
*               At x = b is asked while
*               we have AP = LU   where P is the column permutation
*               due to max trans.
*               Therefore we need to modify rhs:
*                 b' = P-1 b   (P-1=Pt)
*               Apply column permutation to the right hand side RHS
*               Column J of the permuted matrix corresponds to
*               column PERMW(J) of the original matrix.
*
C               ==========
C               SPARSE RHS : permute indices rather than values
C               ==========
C               Solve with At X = B should never occur for A-1 
                IPOS = 1
                DO I=1, NBCOL_INBLOC
C                 -- note that IRHS_PTR is compressed
                  COLSIZE = IRHS_PTR_COPY(I+1) - IRHS_PTR_COPY(I)
                  DO K = 1, COLSIZE
                   JPERM = UNS_PERM_INV(IRHS_SPARSE_COPY(IPOS+K-1))
                   IRHS_SPARSE_COPY(IPOS+K-1) = JPERM
                  ENDDO
                  IPOS = IPOS + COLSIZE
                ENDDO
              ENDIF ! MTYPE.NE.1
            ENDIF ! KEEP(23).NE.0
          ENDIF ! NZ_THIS_BLOCK .NE. 0
C         -----
          ENDIF  !  ============ KEEP(248)==1
C         -----
        ENDIF   !  (id%MYID .eq. MASTER)
C
C =====================  ERROR handling and propagation ================
 30     CONTINUE
        CALL MUMPS_PROPINFO( ICNTL(1), INFO(1),
     &                   id%COMM,id%MYID)
        IF (INFO(1) .LT.0 ) GOTO 90
C ======================================================================
C
C       NBCOL_INBLOC depends on loop
        IF (KEEP(248)==1) THEN
            CALL MPI_BCAST( NBCOL_INBLOC,1, MPI_INTEGER,
     &           MASTER, id%COMM,IERR)
        ELSE
          NBCOL_INBLOC = NBRHS_EFF
        ENDIF
        JEND_RHS =JBEG_RHS + NBCOL_INBLOC - 1
        IF ((KEEP(111).eq.0).AND.(KEEP(252).EQ.0)
     &      .AND.(KEEP(221).NE.2 ).AND.(KEEP(248).NE.0) ) THEN
C         ----------------------------
C         -- SPARSE RIGHT-HAND-SIDE
C         ----------------------------
          CALL MPI_BCAST( NZ_THIS_BLOCK,1, MPI_INTEGER,
     &                    MASTER, id%COMM,IERR)
          IF (id%MYID.NE.MASTER .and. NZ_THIS_BLOCK.NE.0) THEN
            ALLOCATE(IRHS_SPARSE_COPY(NZ_THIS_BLOCK),
     &               stat=allocok)
            if (allocok .GT.0 ) then
               INFO(1)=-13
               INFO(2)=NZ_THIS_BLOCK
               GOTO 45
            endif
            IRHS_SPARSE_COPY_ALLOCATED=.TRUE.
C           RHS_SPARSE_COPY is broadcasted
C           for A-1 even if on the slaves the initialisation of the RHS
C           could be only based on the pattern. Doing so we
C           broadcast the scaled version of the RHS (scaling arrays
C           that are not available on slaves).
            ALLOCATE(RHS_SPARSE_COPY(NZ_THIS_BLOCK),
     &               stat=allocok)
            if (allocok .GT.0 ) then
               INFO(1)=-13
               INFO(2)=NZ_THIS_BLOCK
               GOTO 45
            endif
            RHS_SPARSE_COPY_ALLOCATED=.TRUE.
            NB_BYTES = NB_BYTES + int(NZ_THIS_BLOCK,8)*(K34_8+K35_8)
            NB_BYTES_MAX = max(NB_BYTES_MAX,NB_BYTES)
C
            ALLOCATE(IRHS_PTR_COPY(NBCOL_INBLOC+1),stat=allocok)
            if (allocok .GT.0 ) then
               INFO(1)=-13
               INFO(2)=NBCOL_INBLOC+1
               GOTO 45
            endif
            IRHS_PTR_COPY_ALLOCATED = .TRUE.
            NB_BYTES = NB_BYTES + int(NBCOL_INBLOC+1,8)*K34_8
            NB_BYTES_MAX = max(NB_BYTES_MAX,NB_BYTES)
          ENDIF
C
C =====================  ERROR handling and propagation ================
 45       CONTINUE
          CALL MUMPS_PROPINFO( ICNTL(1), INFO(1),
     &                   id%COMM,id%MYID)
          IF (INFO(1) .LT.0 ) GOTO 90
C ======================================================================
          IF (NZ_THIS_BLOCK > 0) THEN
            CALL MPI_BCAST(IRHS_SPARSE_COPY(1),
     &                 NZ_THIS_BLOCK,
     &                 MPI_INTEGER,
     &                 MASTER, id%COMM,IERR)
            CALL MPI_BCAST(IRHS_PTR_COPY(1),
     &             NBCOL_INBLOC+1,
     &             MPI_INTEGER,
     &             MASTER, id%COMM,IERR)
            IF (IERR.GT.0) THEN
               WRITE (*,*)'NOT OK FOR ALLOC PTR ON SLAVES'
               call MUMPS_ABORT()
            ENDIF
          ENDIF
        ENDIF
C
C       ==========================================================
C       INITIALIZE POSINRHSCOMP_ROW/COL, RHSCOMP and related data
C       ==========================================================
        IF ( I_AM_SLAVE ) THEN
C         --------------------------------------------------
C         If I am involved in the solve and if
C         either
C           no null space comput (keep(111)=0) and sparse rhs
C         or
C             null space computation
C         then
C           compute POSINRHSCOMP
C         endif
C
C         Fwd in facto: in this case only POSINRHSCOMP need be computed
C
C         (POSINRHSCOMP_ROW/COL indirection arrays should
C         have been allocated once outside loop)
C         Compute size of RHSCOMP since it might depend
C            on the process index and of the sparsity of the RHS
C            if it is exploited.
C         Initialize POSINRHSCOMP_ROW/COL
C
C         Note that LD_RHSCOMP and LENRHSCOMP
C         are not initialized on the host in
C         the case of a non-working host.
C         Note that POSINRHSCOMP is now always computed in SOL_DRIVER
C         at least during the first block of RHS when sparsity of RHS
C         is not exploited.
C         -------------------------------
C         INITTIALZE POSINRHSCOMP_ROW/COL
C         -------------------------------
C
          IF ( KEEP(221).EQ.2 .AND. KEEP(252).EQ.0
     &      .AND.  (KEEP(248).EQ.0 .OR. (id%NRHS.EQ.1))
     &      ) THEN
C           Reduced RHS was already computed during
C           a previous forward step AND is valid.
C           By valid we mean:
C            -no forward in facto  (KEEP(252)==0) during which
C             POSINRHSCOMP was not computed
C           AND
C            -no exploit sparsity with multiple RHS
C            because in this case POSINRHSCOMP would
C            be valid only for the last block processed during fwd.
C           In those cases since we only perform backward step, we do not
C           need to compute POSINRHSCOMP
            BUILD_POSINRHSCOMP = .FALSE.
          ENDIF
C         ------------------------
C         INITIALIZE POSINRHSCOMP
C         ------------------------
          IF (BUILD_POSINRHSCOMP) THEN
C           -- we first set MTYPE_LOC and
C           -- reset BUILD_POSINRHSCOMP for next iteration in loop
C
C           general case only POSINRHSCOMP is computed
            BUILD_POSINRHSCOMP = .FALSE.
!           POSINRHSCOMP does not change between blocks
            MTYPE_LOC = MTYPE
C
            IF ( (KEEP(111).NE.0) .OR. (KEEP(237).NE.0) .OR.
     &         (KEEP(252).NE.0) ) THEN
C
              IF (KEEP(111).NE.0) THEN
C               -- in the context of null space, we need to
C               -- build RHSCOMP to skip SOL_R. Therefore
C               -- we need to know for each concerned
C               -- row index its position in
C               -- RHSCOMP
C               We use row indices, as these are the ones that
C               were used to detect zero pivots during factorization.
C               POSINRHSCOMP_ROW will allow to find the (row) index of a
C               zero in RHSCOMP before calling ZMUMPS_SOL_S. Then
C               ZMUMPS_SOL_S uses column indices to build the solution
C               (corresponding to null space vectors)
                MTYPE_LOC = 1
              ELSE IF  (KEEP(252).NE.0) THEN
C               -- Fwd in facto: since fwd is skipped we need to build POSINRHSCOMP
                MTYPE_LOC = 1  ! (no transpose)
C     BUILD_POSINRHSCOMP = .FALSE.  ! POSINRHSCOMP does not change between blocks
              ELSE
C               -- A-1 only
                MTYPE_LOC = MTYPE
                BUILD_POSINRHSCOMP = .TRUE.
              ENDIF
            ENDIF
C           -- compute POSINRHSCOMP
            LIW_PASSED=max(1,LIW)
            IF (KEEP(237).EQ.0) THEN
              CALL ZMUMPS_BUILD_POSINRHSCOMP(
     &           id%NSLAVES,id%N,
     &           id%MYID_NODES, id%PTLUST_S(1),
     &           id%KEEP(1),id%KEEP8(1),
     &           id%PROCNODE_STEPS(1), id%IS(1), LIW_PASSED,
     &           id%STEP(1),
     &           id%POSINRHSCOMP_ROW(1), id%POSINRHSCOMP_COL(1),
     &           id%POSINRHSCOMP_COL_ALLOC,
     &           MTYPE_LOC,
     &           NBENT_RHSCOMP, NB_FS_IN_RHSCOMP_TOT )
                NB_FS_IN_RHSCOMP_F = NB_FS_IN_RHSCOMP_TOT
            ELSE
              CALL ZMUMPS_BUILD_POSINRHSCOMP_AM1(
     &           id%NSLAVES,id%N,
     &           id%MYID_NODES, id%PTLUST_S(1), id%DAD_STEPS(1),
     &           id%KEEP(1),id%KEEP8(1),
     &           id%PROCNODE_STEPS(1), id%IS(1), LIW,
     &           id%STEP(1),
     &           id%POSINRHSCOMP_ROW(1), id%POSINRHSCOMP_COL(1),
     &           id%POSINRHSCOMP_COL_ALLOC,
     &           MTYPE_LOC,
     &           IRHS_PTR_COPY(1), NBCOL_INBLOC, IRHS_SPARSE_COPY(1),
     &           NZ_THIS_BLOCK,PERM_RHS, size(PERM_RHS) , JBEG_RHS,
     &           NBENT_RHSCOMP,
     &           NB_FS_IN_RHSCOMP_F, NB_FS_IN_RHSCOMP_TOT,
     &            UNS_PERM_INV, size(UNS_PERM_INV) ! size 1 if not used
     &            )
            ENDIF
          ENDIF   ! BUILD_POSINRHSCOMP=.TRUE.
          IF (KEEP(221).EQ.1) THEN
C           we need save the reduced RHS for all RHS to perform
C           later the backward phase with an updated reduced RHS
C           thus we allocate NRHS_NONEMPTY columns in one shot.
C           Note that RHSCOMP might have been allocated in previous block
C           and RHSCOMP has been deallocated previous to entering loop on RHS
            IF (.not. associated(id%RHSCOMP)) THEN
C             So far we cannot combine this to exploit sparsity
C             so that NBENT_RHSCOMP will not change in the loop
C             and can be used to dimension RHSCOMP
C             Furthermore, during bwd phase the REDRHS provided
C             by the user might also have a different non empty
C             column patter than the sparse RHS provided on input to
C             this phase: thus we need to allocate id%NRHS columns too.
              LD_RHSCOMP = max(NBENT_RHSCOMP,1)
              LENRHSCOMP = LD_RHSCOMP*id%NRHS
              ALLOCATE (id%RHSCOMP(LENRHSCOMP),  stat = allocok )
              IF ( allocok .GT. 0 ) THEN
                 INFO(1)=-13
                 INFO(2)=LENRHSCOMP
                 GOTO 41
              END IF
              NB_BYTES = NB_BYTES + int(size(id%RHSCOMP),8)*K35_8
              NB_BYTES_MAX = max(NB_BYTES_MAX,NB_BYTES)
            ENDIF
          ENDIF
          IF ((KEEP(221).NE.1).AND.
     &        ((KEEP(221).NE.2).OR.(KEEP(252).NE.0))
     &       ) THEN
C           ------------------
C           Allocate RHSCOMP (case of RHSCOMP allocated at each block of RHS)
C           ------------------
C           RHSCOMP allocated per block of maximum size NBRHS
            LD_RHSCOMP = max(NBENT_RHSCOMP, LD_RHSCOMP)
C           NBRHS_EFF could be used instead on NBRHS
            LENRHSCOMP = LD_RHSCOMP*NBRHS
            IF (associated(id%RHSCOMP)) THEN
              IF  ((size(id%RHSCOMP).LT.LENRHSCOMP).OR.
     &             (KEEP(235).NE.0).OR.(KEEP(237).NE.0) ) THEN
                ! deallocate and reallocate if:
                ! _larger array needed
                ! OR 
                ! _exploit sparsity/A-1: since LENRHSCOMP is expected to
                ! vary much in these cases, this should improve locality
                 NB_BYTES = NB_BYTES - int(size(id%RHSCOMP),8)*K35_8
                 DEALLOCATE(id%RHSCOMP)
                 NULLIFY(id%RHSCOMP)
              ENDIF
            ENDIF
            IF (.not. associated(id%RHSCOMP)) THEN
              LD_RHSCOMP = max(NBENT_RHSCOMP, 1)
              LENRHSCOMP = LD_RHSCOMP*NBRHS
              ALLOCATE (id%RHSCOMP(LENRHSCOMP),  stat = allocok )
              IF ( allocok .GT. 0 ) THEN
               INFO(1)=-13
               INFO(2)=LENRHSCOMP
               GOTO 41
              END IF
              NB_BYTES = NB_BYTES + int(size(id%RHSCOMP),8)*K35_8
              NB_BYTES_MAX = max(NB_BYTES_MAX,NB_BYTES)
            ENDIF
          ENDIF
          IF (KEEP(221).EQ.2) THEN
C           RHSCOMP has been allocated (call with KEEP(221).EQ.1)
C           even in the case fwd in facto
            LENRHSCOMP = size(id%RHSCOMP)
            ! Not correct: LD_RHSCOMP = LENRHSCOMP/id%NRHS_NONEMPTY
            LD_RHSCOMP = LENRHSCOMP/id%NRHS
          ENDIF
C
C         Shift on RHSCOMP
C
          IF ( KEEP(221).EQ.0 ) THEN
C            -- RHSCOMP reused in the loop
             IBEG_RHSCOMP= 1
          ELSE
C            Initialize IBEG_RHSCOMP
#if defined(RHSCOMP_BYROWS)
C            Stored by rows but only inside each
C            block. We keep IBEG_RHSCOMP unchanged
C            for locality since both SCATTER_RHS and
C            GATHER_SOLUTION will be done block-by-block?
             IBEG_RHSCOMP= (JBEG_RHS-1)*LD_RHSCOMP + 1
#else
             IBEG_RHSCOMP= (JBEG_RHS-1)*LD_RHSCOMP + 1
#endif
          ENDIF
        ENDIF   ! I_AM_SLAVE
C =====================  ERROR handling and propagation ================
 41     CONTINUE
        CALL MUMPS_PROPINFO( ICNTL(1), INFO(1),
     &                   id%COMM,id%MYID)
        IF (INFO(1) .LT.0 ) GOTO 90
C ======================================================================
C
        IF (id%MYID .eq. MASTER) THEN
C         =========================
          IF (KEEP(23) .NE. 0) THEN
C         =========================
*         maximum transversal was performed
            IF (MTYPE .NE. 1) THEN
*             At x = b is asked while
*             we have AP = LU   where P is the column permutation
*             due to max trans.
*             Therefore we need to modify rhs:
*               b' = P-1 b   (P-1=Pt)
*             Apply column permutation to the right hand side RHS
*             Column J of the permuted matrix corresponds to
*             column PERMW(J) of the original matrix.
*
              IF (KEEP(248)==0) THEN
C               =========
C               DENSE RHS : permute values in RHS
C               =========
                ALLOCATE( C_RW2( id%N ),stat =allocok )
                IF ( allocok .GT. 0 ) THEN
                  INFO(1)=-13
                  INFO(2)=id%N
                  IF ( LPOK) THEN
                    WRITE(LP,*) id%MYID,
     &              ':Error allocating C_RW2 in ZMUMPS_SOLVE_DRIVE'
                  END IF
                  GOTO 30
                END IF
C               We directly permute in id%RHS.
                DO K = 1, NBRHS_EFF
                  KDEC = IBEG+(K-1)*LD_RHS
                  DO I = 1, id%N
                    C_RW2(I)=id%RHS(I-1+KDEC)
                  END DO
                  DO I = 1, id%N
                    JPERM = id%UNS_PERM(I)
                    id%RHS(I-1+KDEC) = C_RW2(JPERM)
                  END DO
                END DO
                DEALLOCATE(C_RW2)
              ENDIF
            ENDIF
          ENDIF
C
          IF (ERANA) THEN
            IF ( KEEP(248) == 0 ) THEN
              DO K = 1, NBRHS_EFF
                KDEC = IBEG+(K-1)*LD_RHS
                DO I = 1, id%N
                  SAVERHS(I+(K-1)*id%N) = id%RHS(KDEC+I-1)
                END DO
              ENDDO
            ELSE
              SAVERHS(:) = ZERO
              DO K = 1, NBRHS
                DO J = id%IRHS_PTR(K), id%IRHS_PTR(K+1)-1
                  I = id%IRHS_SPARSE(J)
                  SAVERHS(I+(K-1)*id%N) = id%RHS_SPARSE(J)
                ENDDO
              ENDDO
            ENDIF
          ENDIF
C
C         RHS is set to scaled right hand side
C
          IF (LSCAL) THEN
C           scaling was performed
            IF (KEEP(248)==0) THEN
C             dense RHS
              IF (MTYPE .EQ. 1) THEN
C               we solve Ax=b, use ROWSCA to scale the RHS
                DO K =1, NBRHS_EFF
                  KDEC = (K-1) * LD_RHS + IBEG - 1
                  DO I = 1, id%N
                    id%RHS(KDEC+I) = id%RHS(KDEC+I) *
     &                               id%ROWSCA(I)
                  ENDDO
                ENDDO
              ELSE
C               we solve Atx=b, use COLSCA to scale the RHS
                DO K =1, NBRHS_EFF
                  KDEC = (K-1) * LD_RHS + IBEG - 1
                  DO I = 1, id%N
                   id%RHS(KDEC+I) = id%RHS(KDEC+I) *
     &                              id%COLSCA(I)
                  ENDDO
                ENDDO
              ENDIF
            ELSE
C             -------------------------
C             KEEP(248)==1 (and MASTER)
C             -------------------------
              KDEC=id%IRHS_PTR(JBEG_RHS)
C             Compute
              IF ((KEEP(248)==1) .AND.
     &          (DO_PERMUTE_RHS.OR.INTERLEAVE_PAR.OR.
     &           (id%KEEP(237).NE.0))
     &        ) THEN
C               -- copy from RHS_SPARSE need be done per
C                  column following PERM_RHS
C               Columns are not contiguous and need be copied one by one
                IPOS = 1
                J    = 0
                DO I=JBEG_RHS, JBEG_RHS + NBCOL_INBLOC -1
                  J = J+1
C                 Note that we work here on compressed IRHS_PTR_COPY
                  COLSIZE = IRHS_PTR_COPY(J+1) - IRHS_PTR_COPY(J)
C                 -- skip empty column
                  IF (COLSIZE .EQ. 0) CYCLE
                  IF (id%KEEP(237).NE.0) THEN
                    IF (DO_PERMUTE_RHS.OR.INTERLEAVE_PAR) THEN
C                     if A-1 only, then, for each non empty target
C                     column PERM_RHS(I), scale in first position
C                     in column the diagonal entry
C                     build the scaled rhs ej on each slave.
                      RHS_SPARSE_COPY(IPOS) =   id%ROWSCA(PERM_RHS(I)) *
     &                       ONE
                    ELSE
                      RHS_SPARSE_COPY(IPOS) =   id%ROWSCA(I) * ONE
                    ENDIF
                  ELSE
                    DO K = 1, COLSIZE
                      II = id%IRHS_SPARSE(id%IRHS_PTR(PERM_RHS(I))+K-1)
                        !! ORIGINAL IRHS_PTR !!
                      IF (MTYPE.EQ.1) THEN
                        RHS_SPARSE_COPY(IPOS+K-1) =
     &                  id%RHS_SPARSE(id%IRHS_PTR(PERM_RHS(I))+K-1)*
     &                  id%ROWSCA(II)
                      ELSE
                        RHS_SPARSE_COPY(IPOS+K-1) =
     &                  id%RHS_SPARSE(id%IRHS_PTR(PERM_RHS(I))+K-1)*
     &                  id%COLSCA(II)
                      ENDIF
                    ENDDO
                  ENDIF
                  IPOS = IPOS + COLSIZE
                ENDDO
              ELSE
              ! general sparse RHS (KEEP(237) .EQ.0)
                IF (MTYPE .eq. 1) THEN
                  DO IZ=1,NZ_THIS_BLOCK
                    I=IRHS_SPARSE_COPY(IZ)
                    RHS_SPARSE_COPY(IZ)=id%RHS_SPARSE(KDEC+IZ-1)*
     &                            id%ROWSCA(I)
                  ENDDO
                ELSE
                  DO IZ=1,NZ_THIS_BLOCK
                    I=IRHS_SPARSE_COPY(IZ)
                    RHS_SPARSE_COPY(IZ)=id%RHS_SPARSE(KDEC+IZ-1)*
     &                                  id%COLSCA(I)
                  ENDDO
                ENDIF
              ENDIF
            ENDIF  ! KEEP(248)==1
          ENDIF  ! LSCAL
        ENDIF  ! id%MYID.EQ.MASTER
#if defined(V_T)
        CALL VTEND(perm_scal_ini,IERR)
#endif
C
C       Prepare RHS on master
C       END
C       =====================
        IF ((KEEP(248).EQ.1).AND.(KEEP(237).EQ.0)) THEN
          ! case of general sparse: in case of empty columns
          ! modifed version of
          ! NBRHS_EFF need be broadcasted since it is used
          ! to update BEG_RHS at the end of the DO WHILE
            CALL MPI_BCAST( NBRHS_EFF,1, MPI_INTEGER,
     &           MASTER, id%COMM,IERR)
            CALL MPI_BCAST(NB_RHSSKIPPED,1,MPI_INTEGER,MASTER,
     &               id%COMM,IERR)
        ENDIF
C       -----------------------------------
C       Two main cases depending on option
C       for null space computation:
C
C       KEEP(111)=0 : use RHS from user
C                     (sparse or dense)
C       KEEP(111)!=0: build an RHS on each
C                     proc for null space
C                     computations
C       -----------------------------------
#if defined(V_T)
        CALL VTBEGIN(soln_dist,IERR)
#endif
        IF ((KEEP(111).eq.0).AND.(KEEP(252).EQ.0)
     &        .AND.(KEEP(221).NE.2 )) THEN
C           ------------------------
C           Use RHS provided by user
C           when not null space and not Fwd in facto
C           ------------------------
            IF (KEEP(248) == 0) THEN
C             ----------------------------
C             -- DENSE RIGHT-HAND-SIDE
C             ----------------------------
              IF(id%MYID.EQ.MASTER) TIMESCATTER1=MPI_WTIME()
              IF ( .NOT.I_AM_SLAVE ) THEN
C               -- Master not working
                CALL ZMUMPS_SCATTER_RHS(id%NSLAVES,id%N, id%MYID,
     &          id%COMM,
     &          MTYPE, id%RHS(IBEG), LD_RHS, NBRHS_EFF,
     &          NBRHS_EFF,
     &          C_DUMMY, 1, 1,
     &          IDUMMY, 0,
     &          JDUMMY, id%KEEP(1), id%KEEP8(1), id%PROCNODE_STEPS(1),
     &          IDUMMY, 1,
     &          id%STEP(1),
     &          id%ICNTL(1),id%INFO(1))
              ELSE
                IF (id%MYID .eq. MASTER) THEN
                  PTR_RHS => id%RHS
                  LD_RHS_loc   = LD_RHS
                  NCOL_RHS_loc = NBRHS_EFF
                  IBEG_loc     = IBEG
                ELSE
                  PTR_RHS => CDUMMY_TARGET
                  LD_RHS_loc     = 1
                  NCOL_RHS_loc   = 1
                  IBEG_loc       = 1
                ENDIF
                LIW_PASSED = max( LIW, 1 )
                CALL ZMUMPS_SCATTER_RHS(id%NSLAVES,id%N, id%MYID,
     &          id%COMM,
     &          MTYPE, PTR_RHS(IBEG_loc),LD_RHS_loc,NCOL_RHS_loc,
     &          NBRHS_EFF,
     &          id%RHSCOMP(IBEG_RHSCOMP), LD_RHSCOMP, NBRHS_EFF,
     &          id%POSINRHSCOMP_ROW(1), NB_FS_IN_RHSCOMP_F,
C
     &          id%PTLUST_S(1), id%KEEP(1), id%KEEP8(1),
     &          id%PROCNODE_STEPS(1),
     &          IS(1), LIW_PASSED,
     &          id%STEP(1),
     &          id%ICNTL(1),id%INFO(1))
              ENDIF
              IF(id%MYID.EQ.MASTER) THEN
                 TIMESCATTER2=MPI_WTIME()-TIMESCATTER1+TIMESCATTER2
              ENDIF
              IF (INFO(1).LT.0) GOTO 90
            ELSE
C             ===  KEEP(248)==1 =========
C             -- SPARSE RIGHT-HAND-SIDE
C             ----------------------------
              IF (NZ_THIS_BLOCK > 0) THEN
                CALL MPI_BCAST(RHS_SPARSE_COPY(1),
     &                     NZ_THIS_BLOCK,
     &                     MPI_DOUBLE_COMPLEX,
     &                     MASTER, id%COMM, IERR)
              ENDIF
C             -- At this point each process has a copy of the
C             -- sparse RHS. We need to store it into RHSCOMP.
C
              IF (KEEP(237).NE.0)  THEN
                IF ( I_AM_SLAVE ) THEN
c                 -----
c                 case of A-1
c                 -----
c                 - Take columns with non-zero entry, say j,
*                 - to build Ej and store it in RHSCOMP
                  K=1              ! Column index in RHSCOMP
                  id%RHSCOMP(1:NBRHS_EFF*LD_RHSCOMP) = ZERO
                  IPOS = 1
                  DO I = 1, NBCOL_INBLOC
                    COLSIZE = IRHS_PTR_COPY(I+1) - IRHS_PTR_COPY(I)
                    IF (COLSIZE.GT.0) THEN
                      ! Find global column index J and set
                      ! column K of RHSCOMP to ej (here IBEG is one)
                      J = I - 1  + JBEG_RHS
                      IF (DO_PERMUTE_RHS.OR.INTERLEAVE_PAR) THEN
                        J = PERM_RHS(J)
                      ENDIF
                      IPOSRHSCOMP = id%POSINRHSCOMP_ROW(J)
C                     IF ( (IPOSRHSCOMP.LE.NB_FS_IN_RHSCOMP_F)
C     &               .AND.(IPOSRHSCOMP.GT.0) ) THEN
                      IF (IPOSRHSCOMP.GT.0)  THEN
C                       Columns J corresponds to ej and thus to variable j
C                       that is on my proc
C                       Note that :
C                       In first entry in column
C                       we have and MUST have already scaled value of diagonal.
C                       This need have been done on master because we do not
C                       have scaling arrays available on slaves.
C                       Furthermore we know that only one entry is
C                       needed the diagonal entry (for the forward with A-1).
C
#if defined(RHSCOMP_BYROWS)
                        id%RHSCOMP((IPOSRHSCOMP-1)*NBRHS_EFF+K) =
     &                      RHS_SPARSE_COPY(IPOS)
#else
                        id%RHSCOMP((K-1)*LD_RHSCOMP+IPOSRHSCOMP) =
     &                       RHS_SPARSE_COPY(IPOS)
#endif
                      ENDIF               ! End of J on my proc
                      K = K + 1
                      IPOS = IPOS + COLSIZE  ! go to next column
                    ENDIF
                  ENDDO
                  IF (K.NE.NBRHS_EFF+1) THEN
                    WRITE(6,*) 'Internal Error 9 in solution driver ',
     &              K,NBRHS_EFF
                    call MUMPS_ABORT()
                  ENDIF
                ENDIF ! I_AM_SLAVE
C               -------
c               END A-1
C               -------
              ELSE
C               --------------
C               General sparse
C               --------------
C               -- reset to zero RHSCOMP for skipped columns (if any)
                IF ((KEEP(221).EQ.1).AND.(NB_RHSSKIPPED.GT.0)
     &            .AND.I_AM_SLAVE) THEN
#if defined(RHSCOMP_BYROWS)
                  WRITE(*,*) "Internal error 17 is sol driver"
                  CALL MUMPS_ABORT()
#else
                  DO K = JBEG_RHS-NB_RHSSKIPPED, JBEG_RHS-1
                    DO I = 1,  LD_RHSCOMP
                      id%RHSCOMP((K-1)*LD_RHSCOMP + I) =  ZERO
                    ENDDO
                  ENDDO
#endif
                ENDIF
#if defined(RHSCOMP_BYROWS)
                IF (I_AM_SLAVE) THEN
                  DO I=1, NBENT_RHSCOMP
                    DO K = 1, NBCOL_INBLOC
C                     NBCOL_INBLOC is equal to NBRHS_EFF in this case
                      id%RHSCOMP(IBEG_RHSCOMP+(I-1)*NBRHS_EFF+K-1)=ZERO
                    ENDDO
                  ENDDO
                ENDIF
C               Test below must be done also on non-working host !!
                IF (NZ_THIS_BLOCK .EQ. 0 .AND. KEEP(221).EQ.1) THEN
C                 Skip the rest, go to next block.
                  GOTO 1000
                ENDIF
                IF (I_AM_SLAVE) THEN
                  DO K = 1, NBCOL_INBLOC
!                   it is equal to NBRHS_EFF in this case
                    KDEC = IBEG_RHSCOMP + K - 1
#else
                IF (I_AM_SLAVE) THEN
                  DO K = 1, NBCOL_INBLOC
!                   it is equal to NBRHS_EFF in this case
                    KDEC = (K-1) * LD_RHSCOMP + IBEG_RHSCOMP - 1
                    id%RHSCOMP(KDEC+1:KDEC+NBENT_RHSCOMP) = ZERO
#endif
                    DO IZ=IRHS_PTR_COPY(K), IRHS_PTR_COPY(K+1)-1
                      I=IRHS_SPARSE_COPY(IZ)
                      IPOSRHSCOMP = id%POSINRHSCOMP_ROW(I)
C                     Since all fully summed variables mapped
C                     on each proc are stored at the beginning
C                     of RHSCOMP, we can compare to KEEP(89)
C                     to know if RHSCOMP should be initialized
C                     So far the tree has not been pruned to exploit
C                     sparsity to compress RHSCOMP so we compare to
C                     NB_FS_IN_RHSCOMP_TOT
                      IF ( (IPOSRHSCOMP.LE.NB_FS_IN_RHSCOMP_TOT)
     &                .AND.(IPOSRHSCOMP.GT.0) ) THEN
C                        ! I is fully summed var mapped on my proc
#if defined(RHSCOMP_BYROWS)
                        id%RHSCOMP(KDEC+(IPOSRHSCOMP-1)*NBRHS_EFF)=
     &                  id%RHSCOMP(KDEC+(IPOSRHSCOMP-1)*NBRHS_EFF) +
     &                  RHS_SPARSE_COPY(IZ)
#else
                        id%RHSCOMP(KDEC+IPOSRHSCOMP)=
     &                  id%RHSCOMP(KDEC+IPOSRHSCOMP) +
     &                  RHS_SPARSE_COPY(IZ)
#endif
                      ENDIF
                    ENDDO
                  ENDDO
                END IF ! I_AM_SLAVE
              ENDIF ! KEEP(237)
            ENDIF  ! ==== KEEP(248)==1 =====
C
        ELSE IF (I_AM_SLAVE) THEN
            ! I_AM_SLAVE AND (null space or Fwd in facto)
            IF (KEEP(111).NE.0) THEN
C             -----------------------
C             Null space computations
C             -----------------------
C
C             We are working on columns BEG_RHS:BEG_RHS+NBRHS_EFF-1
C             of RHS.
C             Columns in 1..KEEP(112):
C                    Put a one in corresponding
C                    position of the right-hand-side,
C                    and zeros in other places.
C             Columns in KEEP(112)+1: KEEP(112)+KEEP(17):
C                    root node => set
C                    0 everywhere and compute the local range
C                    corresponding to IBEG/IEND in root
C                    that will be passed to ZMUMPS_SEQ_SOLVE_ROOT_RR
C                    Also keep track of which part of
C                    ZMUMPS_RHS must be passed to
C                    ZMUMPS_SEQ_SOLVE_ROOT_RR.
C
              IF (KEEP(111).GT.0) THEN
                IBEG_GLOB_DEF = KEEP(111)
                IEND_GLOB_DEF = KEEP(111)
              ELSE
                IBEG_GLOB_DEF = BEG_RHS
                IEND_GLOB_DEF = BEG_RHS+NBRHS_EFF-1
              ENDIF
              IF ( id%KEEP(112) .GT. 0 .AND. DO_NULL_PIV) THEN
                IF (IBEG_GLOB_DEF .GT.id%KEEP(112)) THEN
                  id%KEEP(235) = 0
                  DO_NULL_PIV = .FALSE.
                ENDIF
                IF (IBEG_GLOB_DEF .LT.id%KEEP(112)
     &          .AND. IEND_GLOB_DEF .GT.id%KEEP(112)
     &          .AND. DO_NULL_PIV ) THEN
C                 IEND_GLOB_DEF = id%KEEP(112)
C                 forcing exploit sparsity
C                 - cannot be done at this point
C                 - and is not what the user would have expected the
C                   code to to do anyway !!!!
C                 suppress:  id%KEEP(235) = 1  ! End Block of sparsity ON
                  DO_NULL_PIV = .FALSE.
                ENDIF
              ENDIF
              IF (id%KEEP(235).NE.0) THEN
C               Exploit Sparsity in null space computations
C               We build /allocate the sparse RHS on MASTER
C               based on pivnul_list. Then we broadcast it
C               on the slaves
C               In this case we have ONLY ONE ENTRY per RHS
C
                NZ_THIS_BLOCK=IEND_GLOB_DEF-IBEG_GLOB_DEF+1
                ALLOCATE(IRHS_PTR_COPY(NZ_THIS_BLOCK+1),stat=allocok)
                IF (allocok .GT.0 ) THEN
                  INFO(1)=-13
                  INFO(2)=NZ_THIS_BLOCK
                  GOTO 50
                ENDIF
                IRHS_PTR_COPY_ALLOCATED = .TRUE.
                ALLOCATE(IRHS_SPARSE_COPY(NZ_THIS_BLOCK),stat=allocok)
                IF (allocok .GT.0 ) THEN
                  INFO(1)=-13
                  INFO(2)=NZ_THIS_BLOCK
                  GOTO 50
                ENDIF
                IRHS_SPARSE_COPY_ALLOCATED=.TRUE.
                NB_BYTES = NB_BYTES +
     &                         int(NZ_THIS_BLOCK,8)*(K34_8+K34_8)
     &                         + K34_8
                NB_BYTES_MAX = max(NB_BYTES_MAX,NB_BYTES)
                IF (id%MYID.eq.MASTER) THEN
                  !  compute IRHS_PTR and IRHS_SPARSE_COPY
                  II = 1
                  DO I = IBEG_GLOB_DEF, IEND_GLOB_DEF
                    IRHS_PTR_COPY(I-IBEG_GLOB_DEF+1)      = I
                    IF (INTERLEAVE_PAR .AND.(id%NSLAVES .GT. 1) ) THEN
                      IRHS_SPARSE_COPY(II) = PERM_PIV_LIST(I)
                    ELSE
                      IRHS_SPARSE_COPY(II) = id%PIVNUL_LIST(I)
                    ENDIF
                    II = II +1
                  ENDDO
                  IRHS_PTR_COPY(NZ_THIS_BLOCK+1) = NZ_THIS_BLOCK+1
                ENDIF
C
C =====================  ERROR handling and propagation ================
 50             CONTINUE
                CALL MUMPS_PROPINFO( ICNTL(1), INFO(1),
     &                   id%COMM,id%MYID)
                IF (INFO(1) .LT.0 ) GOTO 90
C ======================================================================
                CALL MPI_BCAST(IRHS_SPARSE_COPY(1),
     &                NZ_THIS_BLOCK,
     &                MPI_INTEGER,
     &                MASTER, id%COMM,IERR)
                CALL MPI_BCAST(IRHS_PTR_COPY(1),
     &                NZ_THIS_BLOCK+1,
     &                MPI_INTEGER,
     &                MASTER, id%COMM,IERR)
C               End IF Exploit Sparsity
              ENDIF
c
C             Initialize RHSCOMP to 0  ! to be suppressed
#if defined(RHSCOMP_BYROWS)
              id%RHSCOMP(1:NBRHS_EFF*LD_RHSCOMP)=ZERO
#else
              DO K=1, NBRHS_EFF
                KDEC = (K-1) *LD_RHSCOMP
                id%RHSCOMP(KDEC+1:KDEC+LD_RHSCOMP)=ZERO
              END DO
#endif
C             Loop over the columns.
C             Note that if ( KEEP(220)+KEEP(109)-1 < IBEG_GLOB_DEF
C             .OR. KEEP(220) > IEND_GLOB_DEF ) then we do not enter
C             the loop.
C             Note that local processor has indices
C             KEEP(220):KEEP(220)+KEEP(109)-1
              IF ((KEEP(235).NE.0).AND.  INTERLEAVE_PAR) THEN
C               When the PIVNUL_LIST has been permuted (in PERM_PIV_LIST)
C               then to exploit sparsity RHSCOMP need be initialized with
c               some care; taking into acount the processor localisation
C               of the indices of the null pivots.
                DO I=IBEG_GLOB_DEF,IEND_GLOB_DEF
C                 Local processor is concerned by I-th column of
C                 global right-hand side.
                  IF (id%MYID_NODES .EQ. MAP_PIVNUL_LIST(I) ) THEN
                    JJ= id%POSINRHSCOMP_ROW(PERM_PIV_LIST(I))
                    IF (JJ.GT.LD_RHSCOMP) THEN
                      WRITE(6,*) ' Internal Error 10 JJ, LD_RHSCOMP=',
     &                JJ, LD_RHSCOMP
                    ENDIF
                    IF (JJ.GT.0) THEN
                      IF (KEEP(50).EQ.0) THEN
C     unsymmetric : always set to fixation used during facto
C     because during factorization we aimed at preserving the
C     sign of the diagonal element, sign here may be different
C     from sign of corresponding diagonal element (not critical)
#if defined(RHSCOMP_BYROWS)
                        id%RHSCOMP(IBEG_RHSCOMP +
     &                  I-IBEG_GLOB_DEF +(JJ-1)* NBRHS_EFF) =
#else
                        id%RHSCOMP(IBEG_RHSCOMP -1+
     &                  (I-IBEG_GLOB_DEF)*LD_RHSCOMP + JJ) =
#endif
     &                  cmplx(abs(id%DKEEP(2)),kind=kind(id%RHSCOMP))
                      ELSE
#if defined(RHSCOMP_BYROWS)
                        id%RHSCOMP(IBEG_RHSCOMP +
     &                  I-IBEG_GLOB_DEF + (JJ-1) *NBRHS_EFF) = ONE
#else
                        id%RHSCOMP(IBEG_RHSCOMP -1+
     &                  (I-IBEG_GLOB_DEF)*LD_RHSCOMP + JJ) = ONE
#endif
                      ENDIF
                    ENDIF
                  ENDIF
                ENDDO
              ELSE
C
C               Computation of null space and computation of backward
C               step incompatible, do one or the other.
                DO I=max(IBEG_GLOB_DEF,KEEP(220)),
     &               min(IEND_GLOB_DEF,KEEP(220)+KEEP(109)-1)
C                 Local processor is concerned by I-th column of
C                 global right-hand side.
                  JJ= id%POSINRHSCOMP_ROW(id%PIVNUL_LIST(I-KEEP(220)+1))
                  IF (JJ.GT.0) THEN
                    IF (KEEP(50).EQ.0) THEN
                      ! unsymmetric : always set to fixation
#if defined(RHSCOMP_BYROWS)
                      id%RHSCOMP( IBEG_RHSCOMP+
     &                   I-IBEG_GLOB_DEF +
     &                   (JJ-1)*NBRHS_EFF )=
#else
                      id%RHSCOMP(IBEG_RHSCOMP-1+
     &                    (I-IBEG_GLOB_DEF)*LD_RHSCOMP +
     &                    JJ)=
#endif
     &                    cmplx(id%DKEEP(2),kind=kind(id%RHSCOMP))
                    ELSE
                      ! Symmetric: always set to one
#if defined(RHSCOMP_BYROWS)
                      id%RHSCOMP(IBEG_RHSCOMP+I-IBEG_GLOB_DEF +
     &                (JJ-1)*NBRHS_EFF)=
#else
                      id%RHSCOMP(IBEG_RHSCOMP-1+
     &                (I-IBEG_GLOB_DEF)*LD_RHSCOMP+JJ)=
#endif
     &                ONE
                    ENDIF
                  ENDIF
                ENDDO
              ENDIF   ! exploit sparsity
              IF ( KEEP(17).NE.0 .AND.
     &             id%MYID_NODES.EQ.MASTER_ROOT) THEN
C               ---------------------------
C               Deficiency of the root node
C               Find range relative to root
C               ---------------------------
C               Among IBEG_GLOB_DEF:IEND_GLOB_DEF, find
C               intersection with KEEP(112)+1:KEEP(112)+KEEP(17)
                IBEG_ROOT_DEF  = max(IBEG_GLOB_DEF,KEEP(112)+1)
                IEND_ROOT_DEF  = min(IEND_GLOB_DEF,KEEP(112)+KEEP(17))
C               First column of right-hand side that must
C               be passed to ZMUMPS_SEQ_SOLVE_ROOT_RR is:
                IROOT_DEF_RHS_COL1 = IBEG_ROOT_DEF-IBEG_GLOB_DEF + 1
C               We look for indices relatively to the root node,
C               substract number of null pivots outside root node
                IBEG_ROOT_DEF = IBEG_ROOT_DEF-KEEP(112)
                IEND_ROOT_DEF = IEND_ROOT_DEF-KEEP(112)
C               Note that if IBEG_ROOT_DEF > IEND_ROOT_DEF, then this
C               means that nothing must be done on the root node
C               for this set of right-hand sides.
              ELSE
                IBEG_ROOT_DEF = -90999
                IEND_ROOT_DEF = -95999
                IROOT_DEF_RHS_COL1= 1
              ENDIF
            ELSE  ! End of null space (test on KEEP(111))
C             case of Fwd in facto
C             id%RHSCOMP need not be initialized. It will be set on the fly
C             to zero for normal fully summed variables of the fronts and
C             to -1 on the roots for the id%N+KEEP(253) variables added
C             to the roots.
            ENDIF ! End of null space (test on KEEP(111))
        ENDIF  ! I am slave
C       -------------------------------------------
C       Reserve space at the end of WORK_WCB on the
C       master of the root node. It will be used to
C       store the reduced RHS.
C       -------------------------------------------
        IF ( I_AM_SLAVE ) THEN
          LWCB_SOL_C = LWCB
          IF ( id%MYID_NODES .EQ. MASTER_ROOT ) THEN
C           This is a special root (otherwise MASTER_ROOT < 0)
            IF ( associated(id%root%RHS_CNTR_MASTER_ROOT) ) THEN
C             RHS_CNTR_MASTER_ROOT may have been allocated
C             during the factorization phase.
              PTR_RHS_ROOT => id%root%RHS_CNTR_MASTER_ROOT
              LPTR_RHS_ROOT = size(id%root%RHS_CNTR_MASTER_ROOT)
            ELSE
C             Otherwise, we use workspace in WCB
              LPTR_RHS_ROOT = NBRHS_EFF * SIZE_ROOT
              IPT_RHS_ROOT  = LWCB - LPTR_RHS_ROOT + 1
              PTR_RHS_ROOT => WORK_WCB(IPT_RHS_ROOT:LWCB)
              LWCB_SOL_C = LWCB_SOL_C - LPTR_RHS_ROOT
            ENDIF
          ELSE
            LPTR_RHS_ROOT = 1
            IPT_RHS_ROOT = LWCB ! Will be passed, but not accessed
            PTR_RHS_ROOT => WORK_WCB(IPT_RHS_ROOT:LWCB)
            LWCB_SOL_C = LWCB_SOL_C - LPTR_RHS_ROOT
          ENDIF
        ENDIF
        IF (KEEP(221) .EQ. 2 ) THEN
C         Copy/send REDRHS in PTR_RHS_ROOT
C         (column by column if leading dimension LD_REDRHS
C         of REDRHS is not equal to SIZE_ROOT).
C         REDRHS was provided on the host
          IF ( ( id%MYID .EQ. MASTER_ROOT_IN_COMM ) .AND.
     &       ( id%MYID .EQ. MASTER ) ) THEN
C         -- Same proc : copy is possible:
            II = 0
            DO K=1, NBRHS_EFF
              KDEC = IBEG_REDRHS+(K-1)*LD_REDRHS -1
              DO I = 1, SIZE_ROOT
                PTR_RHS_ROOT(II+I) = id%REDRHS(KDEC+I)
              ENDDO
              II = II+SIZE_ROOT
            ENDDO
          ELSE
C         -- send REDRHS
            IF ( id%MYID .EQ. MASTER) THEN
C           -- send to MASTER_ROOT_IN_COMM using COMM communicator
C              assert: id%KEEP(116).EQ.SIZE_ROOT
              IF (LD_REDRHS.EQ.SIZE_ROOT) THEN
C              --  One send
                 KDEC = IBEG_REDRHS
                 CALL MPI_SEND(id%REDRHS(KDEC),
     &              SIZE_ROOT*NBRHS_EFF,
     &              MPI_DOUBLE_COMPLEX,
     &              MASTER_ROOT_IN_COMM, 0, id%COMM,IERR)
              ELSE
C             --  NBRHS_EFF sends
                DO K=1, NBRHS_EFF
                  KDEC = IBEG_REDRHS+(K-1)*LD_REDRHS
                  CALL MPI_SEND(id%REDRHS(KDEC),SIZE_ROOT,
     &               MPI_DOUBLE_COMPLEX,
     &               MASTER_ROOT_IN_COMM, 0, id%COMM,IERR)
                ENDDO
              ENDIF
            ELSE IF ( id%MYID .EQ. MASTER_ROOT_IN_COMM ) THEN
C            -- receive from MASTER
              II = 1
              IF (LD_REDRHS.EQ.SIZE_ROOT) THEN
C                -- receive all in on shot
                 CALL MPI_RECV(PTR_RHS_ROOT(II),
     &              SIZE_ROOT*NBRHS_EFF,
     &              MPI_DOUBLE_COMPLEX,
     &              MASTER, 0, id%COMM,STATUS,IERR)
              ELSE
                DO K=1, NBRHS_EFF
                  CALL MPI_RECV(PTR_RHS_ROOT(II),SIZE_ROOT,
     &            MPI_DOUBLE_COMPLEX,
     &            MASTER, 0, id%COMM,STATUS,IERR)
                  II = II + SIZE_ROOT
                ENDDO
              ENDIF
            ENDIF
C         -- other procs are not concerned
          ENDIF
        ENDIF
        TIMEC1=MPI_WTIME()
        IF ( I_AM_SLAVE ) THEN
          LIW_PASSED = max( LIW, 1 )
          LA_PASSED  = max( LA, 1_8 )
C
          IF ((id%KEEP(235).EQ.0).and.(id%KEEP(237).EQ.0) ) THEN
C
C         --- Normal case : we do not exploit sparsity of the RHS
C
            FROM_PP = .FALSE.
            NBSPARSE_LOC = (DO_NBSPARSE.AND.NBRHS_EFF.GT.1)
            PRUNED_SIZE_LOADED = 0_8  ! From ZMUMPS_SOL_ES module
            CALL ZMUMPS_SOL_C(id%root, id%N, id%S(1), LA_PASSED,
     &    IS(1), LIW_PASSED, WORK_WCB(1), LWCB_SOL_C, IWCB, LIWCB,
     &    NBRHS_EFF, id%NA(1),id%LNA,id%NE_STEPS(1), SRW3, MTYPE,
     &    ICNTL(1), FROM_PP, id%STEP(1), id%FRERE_STEPS(1),
     &    id%DAD_STEPS(1), id%FILS(1), id%PTLUST_S(1), id%PTRFAC(1),
     &    IWK_SOLVE, LIWK_SOLVE, id%PROCNODE_STEPS(1),
     &    id%NSLAVES, INFO(1), KEEP(1), KEEP8(1), id%DKEEP(1),
     &    id%COMM_NODES, id%MYID, id%MYID_NODES,
     &    id%BUFR(1), id%LBUFR, id%LBUFR_BYTES,
     &    id%ISTEP_TO_INIV2(1), id%TAB_POS_IN_PERE(1,1),
     &    IBEG_ROOT_DEF, IEND_ROOT_DEF, IROOT_DEF_RHS_COL1,
     &    PTR_RHS_ROOT(1), LPTR_RHS_ROOT, SIZE_ROOT, MASTER_ROOT,
     &    id%RHSCOMP(IBEG_RHSCOMP), LD_RHSCOMP,
     &    id%POSINRHSCOMP_ROW(1), id%POSINRHSCOMP_COL(1)
     &    , 1            , 1       ,  1       , 1
     &    , IDUMMY, 1, JDUMMY, KDUMMY, 1, LDUMMY, 1, MDUMMY
     &    , 1 , 1 , NBSPARSE_LOC, PTR_RHS_BOUNDS(1), LPTR_RHS_BOUNDS
     &    )
          ELSE
C           Exploit sparsity of the RHS (all cases)
C           Remark that JBEG_RHS is already initialized
C
            FROM_PP = .FALSE.
            NBSPARSE_LOC = (DO_NBSPARSE.AND.NBRHS_EFF.GT.1)
            CALL ZMUMPS_SOL_C(id%root, id%N, id%S(1), LA_PASSED,IS(1),
     &    LIW_PASSED, WORK_WCB(1), LWCB_SOL_C, IWCB, LIWCB, NBRHS_EFF,
     &    id%NA(1),id%LNA, id%NE_STEPS(1), SRW3, MTYPE, ICNTL(1),
     &    FROM_PP, id%STEP(1), id%FRERE_STEPS(1), id%DAD_STEPS(1),
     &    id%FILS(1), id%PTLUST_S(1), id%PTRFAC(1), IWK_SOLVE,
     &    LIWK_SOLVE,id%PROCNODE_STEPS(1),id%NSLAVES,INFO(1),KEEP(1),
     &    KEEP8(1), id%DKEEP(1), id%COMM_NODES, id%MYID, id%MYID_NODES,
     &    id%BUFR(1), id%LBUFR, id%LBUFR_BYTES, id%ISTEP_TO_INIV2(1),
     &    id%TAB_POS_IN_PERE(1,1), IBEG_ROOT_DEF, IEND_ROOT_DEF,
     &    IROOT_DEF_RHS_COL1, PTR_RHS_ROOT(1), LPTR_RHS_ROOT, SIZE_ROOT,
     &    MASTER_ROOT, id%RHSCOMP(IBEG_RHSCOMP), LD_RHSCOMP,
     &    id%POSINRHSCOMP_ROW(1), id%POSINRHSCOMP_COL(1),
     &    NZ_THIS_BLOCK, NBCOL_INBLOC, id%NRHS, JBEG_RHS ,
     &    id%Step2node(1), id%KEEP(28), IRHS_SPARSE_COPY(1),
     &    IRHS_PTR_COPY(1), size(PERM_RHS), PERM_RHS, size(UNS_PERM_INV)
C size 1 if not used
     &    , UNS_PERM_INV, NB_FS_IN_RHSCOMP_F, NB_FS_IN_RHSCOMP_TOT
     &    , NBSPARSE_LOC, PTR_RHS_BOUNDS(1), LPTR_RHS_BOUNDS
     &       )
          ENDIF   ! end of exploit sparsity (pruning nodes of the tree)
        END IF
C       -----------------
C       End of slave code
C       -----------------
C
        TIMEC2=MPI_WTIME()-TIMEC1+TIMEC2
C
C       Propagate errors
        CALL MUMPS_PROPINFO( ICNTL(1), INFO(1),
     &                   id%COMM,id%MYID)
C
C       Change error code.
        IF (INFO(1).eq.-2) then
          INFO(1)=-11
          IF (LPOK)
     &    write(LP,*)
     &   ' WARNING : -11 error code obtained in solve'
        END IF
        IF (INFO(1).eq.-3) then
          INFO(1)=-14
          IF (LPOK)
     &    write(LP,*)
     &    ' WARNING : -14 error code obtained in solve'
        END IF
C
C       Return in case of error.
        IF (INFO(1).LT.0) GO TO 90
C
C       ======================================================
C       ONLY FORWARD was performed (case of reduced RHS with Schur
C                               option during factorisation)
C       ======================================================
        IF ( KEEP(221) .EQ. 1 ) THEN ! === Begin OF REDUCED RHS ======
C         --------------------------------------
C         Send (or copy) reduced RHS from PTR_RHS_ROOT located on
C         MASTER_ROOT_IN_COMM to REDRHS located on MASTER (host node).
C         (column by column if leading dimension LD_REDRHS
C         of REDRHS is not equal to SIZE_ROOT)
C         --------------------------------------
          IF ( ( id%MYID .EQ. MASTER_ROOT_IN_COMM ) .AND.
     &        ( id%MYID .EQ. MASTER ) ) THEN
C           -- same proc  --> copy
            II = 0
            DO K=1, NBRHS_EFF
              KDEC = IBEG_REDRHS+(K-1)*LD_REDRHS -1
              DO I = 1, SIZE_ROOT
                id%REDRHS(KDEC+I) = PTR_RHS_ROOT(II+I)
              ENDDO
              II = II+SIZE_ROOT
            ENDDO
          ELSE
C           -- recv in REDRHS
            IF ( id%MYID .EQ. MASTER ) THEN
C             -- recv from MASTER_ROOT_IN_COMM
              IF (LD_REDRHS.EQ.SIZE_ROOT) THEN
C             --  One message to receive
                KDEC = IBEG_REDRHS
                CALL MPI_RECV(id%REDRHS(KDEC),
     &              SIZE_ROOT*NBRHS_EFF,
     &              MPI_DOUBLE_COMPLEX,
     &              MASTER_ROOT_IN_COMM, 0, id%COMM,
     &              STATUS,IERR)
              ELSE
C             --  NBRHS_EFF receives
                DO K=1, NBRHS_EFF
                  KDEC = IBEG_REDRHS+(K-1)*LD_REDRHS
                  CALL MPI_RECV(id%REDRHS(KDEC),SIZE_ROOT,
     &              MPI_DOUBLE_COMPLEX,
     &              MASTER_ROOT_IN_COMM, 0, id%COMM,
     &              STATUS,IERR)
                ENDDO
              ENDIF
            ELSE IF ( id%MYID .EQ. MASTER_ROOT_IN_COMM ) THEN
C           -- send to MASTER
              II = 1
              IF (LD_REDRHS.EQ.SIZE_ROOT) THEN
C               -- send all in on shot
                CALL MPI_SEND(PTR_RHS_ROOT(II),
     &              SIZE_ROOT*NBRHS_EFF,
     &              MPI_DOUBLE_COMPLEX,
     &              MASTER, 0, id%COMM,IERR)
              ELSE
                DO K=1, NBRHS_EFF
                  CALL MPI_SEND(PTR_RHS_ROOT(II),SIZE_ROOT,
     &            MPI_DOUBLE_COMPLEX,
     &            MASTER, 0, id%COMM,IERR)
                  II = II + SIZE_ROOT
                ENDDO
              ENDIF
            ENDIF
C         -- other procs are not concerned
          ENDIF
        ENDIF ! ====== END OF REDUCED RHS (Fwd only performed) ======
C       =======================================================
C       BACKWARD was PERFORMED
C       Postprocess solution that is distributed
        IF ( KEEP(221) .NE. 1 ) THEN  ! BACKWARD was PERFORMED
C       -- KEEP(221).NE.1 => we are sure that backward has been performed
          IF (ICNTL21 == 0) THEN ! CENTRALIZED SOLUTION
C           ========================================================
C           GATHER SOLUTION computed during bwd
C           Each proc holds the pieces of solution corresponding
C           to all fully summed variables mapped on that processor
C           (i.e. corresponding to master nodes mapped on that proc)
C           In case of A-1 we gather directly in RHS_SPARSE
C           the distributed solution.
C           Scaling is done in all case on the fly of the reception
C           Note that when only FORWARD has been performed
C           RSH_MUMPS holds the solution computed during forward step
C           (ZMUMPS_SOL_R)
C           there is no need to copy back in RSH_MUMPS the solution
C           ========================================================
C           centralized solution
            IF (KEEP(237).EQ.0) THEN
C             CWORK not needed for AM1
              ALLOCATE( CWORK(
     &              max(max(KEEP(247),KEEP(246))*NBRHS_EFF,1)),
     &              stat=allocok)
              IF (allocok > 0) THEN
C               Try smaller (it will be less efficient)
                ALLOCATE( CWORK(max(max(KEEP(247),KEEP(246)),1)),
     &                stat=allocok)
                IF (allocok > 0) THEN
                  INFO(1)=-13
                  INFO(2)=max(max(KEEP(247),KEEP(246)),1)
                ENDIF
              ENDIF
            ENDIF
C           Propagate errors
            CALL MUMPS_PROPINFO( ICNTL(1), INFO(1),
     &                       id%COMM,id%MYID)
C           Return in case of error.
            IF (INFO(1).LT.0) GO TO 90
            IF ((id%MYID.NE.MASTER).OR. .NOT.LSCAL) THEN
             PT_SCALING => Dummy_SCAL
            ELSE
              IF (MTYPE.EQ.1) THEN
                PT_SCALING => id%COLSCA
              ELSE
                PT_SCALING => id%ROWSCA
              ENDIF
            ENDIF
            LIW_PASSED = max( LIW, 1 )
            IF(id%MYID.EQ.MASTER) TIMEGATHER1=MPI_WTIME()
            IF ( .NOT.I_AM_SLAVE ) THEN
C             I did not participate to computing part of the solution
C             (id%RHSCOMP not set/allocate) : receive solution, store
C             it and scale it.
              IF (KEEP(237).EQ.0) THEN
C               We need a workspace of minimal size KEEP(247)
C               in order to unpack pieces of the solution.
                CALL ZMUMPS_GATHER_SOLUTION(id%NSLAVES,id%N,
     &              id%MYID, id%COMM, NBRHS_EFF,
     &              MTYPE, id%RHS(IBEG), LD_RHS, NBRHS_EFF,
     &              JDUMMY, id%KEEP(1), id%KEEP8(1),
     &              id%PROCNODE_STEPS(1), IDUMMY, 1,
     &              id%STEP(1), id%BUFR(1), id%LBUFR, id%LBUFR_BYTES,
     &              CWORK(1), size(CWORK),
     &              LSCAL, PT_SCALING(1), size(PT_SCALING),
     &              C_DUMMY, 1 , 1, IDUMMY, 1
     &              )
              ELSE
C               only gather target entries of A-1
                CALL ZMUMPS_GATHER_SOLUTION_AM1(id%NSLAVES,id%N,
     &              id%MYID, id%COMM, NBRHS_EFF,
     &              C_DUMMY, 1, 1,
     &              id%KEEP(1), id%BUFR(1), id%LBUFR, id%LBUFR_BYTES,
     &              LSCAL, PT_SCALING(1), size(PT_SCALING)
C                   --- A-1 related entries
     &             ,IRHS_PTR_COPY(1), size(IRHS_PTR_COPY),
     &              IRHS_SPARSE_COPY(1), size(IRHS_SPARSE_COPY),
     &              RHS_SPARSE_COPY(1), size(RHS_SPARSE_COPY),
     &              UNS_PERM_INV, size(UNS_PERM_INV),
     &              IDUMMY, 1, 0
     &              )
              ENDIF
            ELSE
C             Avoid temporary copy (IS(1)) that some old
C             compilers would do otherwise
              IF (KEEP(237).EQ.0) THEN
                IF (id%MYID.EQ.MASTER) THEN
                  PTR_RHS => id%RHS
                  NCOL_RHS_loc = NBRHS_EFF
                  LD_RHS_loc   = LD_RHS
                  IBEG_loc     = IBEG
                ELSE
                  PTR_RHS => CDUMMY_TARGET
                  NCOL_RHS_loc = 1
                  LD_RHS_loc   = 1
                  IBEG_loc     = 1
                ENDIF
                CALL ZMUMPS_GATHER_SOLUTION(id%NSLAVES,id%N,
     &          id%MYID, id%COMM, NBRHS_EFF,
     &          MTYPE, PTR_RHS(IBEG_loc), LD_RHS_loc, NCOL_RHS_loc,
     &          id%PTLUST_S(1), id%KEEP(1), id%KEEP8(1),
     &          id%PROCNODE_STEPS(1), IS(1), LIW_PASSED,
     &          id%STEP(1), id%BUFR(1), id%LBUFR, id%LBUFR_BYTES,
     &          CWORK(1), size(CWORK),
     &          LSCAL, PT_SCALING(1), size(PT_SCALING),
     &          id%RHSCOMP(IBEG_RHSCOMP), LD_RHSCOMP, NBRHS_EFF,
     &          id%POSINRHSCOMP_COL(1), id%N
     &          )
              ELSE ! only gather target entries of A-1
                CALL ZMUMPS_GATHER_SOLUTION_AM1(id%NSLAVES,id%N,
     &          id%MYID, id%COMM, NBRHS_EFF,
     &          id%RHSCOMP(IBEG_RHSCOMP), LD_RHSCOMP, NBRHS_EFF,
     &          id%KEEP(1), id%BUFR(1), id%LBUFR, id%LBUFR_BYTES,
     &          LSCAL, PT_SCALING(1), size(PT_SCALING)
C               --- A-1 related entries
     &          , IRHS_PTR_COPY(1), size(IRHS_PTR_COPY),
     &          IRHS_SPARSE_COPY(1), size(IRHS_SPARSE_COPY),
     &          RHS_SPARSE_COPY(1), size(RHS_SPARSE_COPY),
     &          UNS_PERM_INV, size(UNS_PERM_INV),
     &          id%POSINRHSCOMP_COL(1), id%N, NB_FS_IN_RHSCOMP_TOT
     &          )
              ENDIF
            ENDIF
            IF (id%MYID.EQ.MASTER) THEN
              TIMEGATHER2=MPI_WTIME()-TIMEGATHER1+TIMEGATHER2
            ENDIF
            IF (KEEP(237).EQ.0)  DEALLOCATE( CWORK )
            IF ( (id%MYID.EQ.MASTER).AND. (KEEP(237).NE.0)
     &        ) THEN
C             Copy back solution from RHS_SPARSE_COPY TO RHS_SPARSE
              IF (DO_PERMUTE_RHS.OR.INTERLEAVE_PAR) THEN
                DO J = JBEG_RHS, JBEG_RHS+NBCOL_INBLOC-1
                  COLSIZE = id%IRHS_PTR(PERM_RHS(J)+1) -
     &                  id%IRHS_PTR(PERM_RHS(J))
                  IF (COLSIZE.EQ.0) CYCLE
                  JJ = J-JBEG_RHS+1
c                 IZ - Column index in Sparse RHS
                  DO IZ= id%IRHS_PTR(PERM_RHS(J)),
     &              id%IRHS_PTR(PERM_RHS(J)+1)-1
                    I = id%IRHS_SPARSE (IZ)
                    DO IZ2 = IRHS_PTR_COPY(JJ),IRHS_PTR_COPY(JJ+1)-1
                      IF (IRHS_SPARSE_COPY(IZ2).EQ.I) EXIT
                      IF (IZ2.EQ.IRHS_PTR_COPY(JJ+1)-1) THEN
                        WRITE(6,*) " Internal Error 13 in solution ",
     &                  " driver, gather "
                        CALL MUMPS_ABORT()
                      ENDIF
                    ENDDO
                    id%RHS_SPARSE(IZ) = RHS_SPARSE_COPY(IZ2)
                  ENDDO
                ENDDO
              ELSE            ! Not (DO_PERMUTE_RHS.OR.INTERLEAVE_PAR)
                DO J = JBEG_RHS, JBEG_RHS+NBCOL_INBLOC-1
                  COLSIZE = id%IRHS_PTR(J+1) - id%IRHS_PTR(J)
                  IF (COLSIZE.EQ.0) CYCLE
                  JJ = J-JBEG_RHS+1
c                 IZ - Column index in Sparse RHS
                  DO IZ= id%IRHS_PTR(J), id%IRHS_PTR(J+1) - 1
                    I = id%IRHS_SPARSE (IZ)
                    DO IZ2 = IRHS_PTR_COPY(JJ),IRHS_PTR_COPY(JJ+1)-1
                      IF (IRHS_SPARSE_COPY(IZ2).EQ.I) EXIT
                      IF (IZ2.EQ.IRHS_PTR_COPY(JJ+1)-1) THEN
                        WRITE(6,*) " Internal Error 14 in solution",
     &                  " driver, gather "
                        CALL MUMPS_ABORT()
                      ENDIF
                    ENDDO
                    id%RHS_SPARSE(IZ) = RHS_SPARSE_COPY(IZ2)
                  ENDDO
                ENDDO
              ENDIF           ! End DO_PERMUTE_RHS.OR.INTERLEAVE_PAR
            ENDIF             ! end A-1 on master
C
C         -- END of backward was performed with centralized solution
          ELSE   ! (KEEP(221).NE.1) .AND.(ICNTL21.NE.0))
C
C           BEGIN of backward performed with distributed solution
C           The non working host should not do this:
            IF ( I_AM_SLAVE ) THEN
              LIW_PASSED = max( LIW, 1 )
C             Only called if more than 1 pivot
C             was eliminated by the processor.
C             Note that LSOL_loc >= KEEP(89)
              IF ( KEEP(89) .GT. 0 ) THEN
                CALL ZMUMPS_DISTRIBUTED_SOLUTION(id%NSLAVES,
     &          id%N, id%MYID_NODES,
     &          MTYPE, id%RHSCOMP(IBEG_RHSCOMP), LD_RHSCOMP,
     &          NBRHS_EFF, id%POSINRHSCOMP_COL(1),
     &          id%ISOL_loc(1),
     &          id%SOL_loc(1), JBEG_RHS-NB_RHSSKIPPED, id%LSOL_loc,
     &          id%PTLUST_S(1), id%PROCNODE_STEPS(1),
     &          id%KEEP(1),id%KEEP8(1),
     &          IS(1), LIW_PASSED,
     &          id%STEP(1), scaling_data, LSCAL, NB_RHSSKIPPED )
              ENDIF
            ENDIF
          ENDIF
C         === BACKWARD was PERFORMED WITH DISTRIBUTED SOLUTION ===
C         ========================================================
        ENDIF ! ==== END of BACKWARD was PERFORMED (KEEP(221).NE.1)
C       note that the main DO-loop on blocks is not ended yet
C
C       ============================================
C       BEGIN
C
C       ITERATIVE REFINEMENT AND/OR ERROR ANALYSIS
C
C       ============================================
        IF ( ICNTL10 > 0 .AND. NBRHS_EFF > 1 ) THEN
C
C         ----------------------------------
C         Multiple RHS: apply a fixed number
C         of iterative refinement steps
C         ----------------------------------
C         DO I = 1, ICNTL10
            write(6,*) ' Internal ERROR 15 in sol_driver '
C           Compute residual:  Y <- SAVERHS - A * RHS
C           Solve RHS <- A^-1 Y, Y modified
C           Assemble in RHS(REDUCE)
C           RHS <- RHS + Y
C         END DO
        END IF
        IF (ERANA) THEN
C
C         SAVERHS holds the original right hand side
C         Sparse rhs are saved in SAVERHS as dense rhs
C
C         * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C         Start iterative refinements. The master is managing the
C         organisation of work, but slaves are used to solve systems of
C         equations and, in case of distributed matrix, perform
C         matrix-vector products. It is more complicated to do this with
C         the SPMD version than it was with the master/slave approach.
C
C         * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c         IF ( PROK  .AND. ICNTL10 .NE. 0 ) WRITE( MP, 270 )
          IF ( PROKG .AND. ICNTL10 .NE. 0 ) WRITE( MPG, 270 )
C         Initializations and allocations
          NITREF = abs(ICNTL10)
          ALLOCATE(R_Y(id%N), stat = allocok)
          IF ( allocok .GT. 0 ) THEN
            INFO(1)=-13
            INFO(2)=id%N
            GOTO 777
          ENDIF
          NB_BYTES = NB_BYTES + int(id%N,8)*K16_8
          ALLOCATE(C_Y(id%N), stat = allocok)
          IF ( allocok .GT. 0 ) THEN
            INFO(1)=-13
            INFO(2)=id%N
            GOTO 777
          ENDIF
          NB_BYTES = NB_BYTES + int(id%N,8)*K35_8
          IF ( id%MYID .EQ. MASTER ) THEN
            ALLOCATE( IW1( 2 * id%N ),stat = allocok )
            IF ( allocok .GT. 0 ) THEN
              INFO(1)=-13
              INFO(2)=2 * id%N
              GOTO 777
            ENDIF
            NB_BYTES = NB_BYTES + int(2*id%N,8)*K34_8
            ALLOCATE( C_W(id%N), stat = allocok )
            IF ( allocok .GT. 0 ) THEN
              INFO(1)=-13
              INFO(2)=id%N
              GOTO 777
            ENDIF
            NB_BYTES = NB_BYTES + int(id%N,8)*K35_8
            ALLOCATE( R_W(2*id%N), stat = allocok )
            IF ( allocok .GT. 0 ) THEN
              INFO(1)=-13
              INFO(2)=id%N
              GOTO 777
            ENDIF
            NB_BYTES = NB_BYTES + int(2*id%N,8)*K16_8
            IF ( PROKG .AND. ICNTL10 .GT. 0 )
     &      WRITE( MPG, 240) 'MAXIMUM NUMBER OF STEPS =', NITREF
C           end allocations on Master
          END IF
          ALLOCATE(C_LOCWK54(id%N),stat = allocok)
          IF ( allocok .GT. 0 ) THEN
            INFO(1)=-13
            INFO(2)=id%N
            GOTO 777
          ENDIF
          NB_BYTES = NB_BYTES + int(id%N,8)*K35_8
          ALLOCATE(R_LOCWK54(id%N),stat = allocok)
          IF ( allocok .GT. 0 ) THEN
            INFO(1)=-13
            INFO(2)=id%N
            GOTO 777
          ENDIF
          NB_BYTES = NB_BYTES + int(id%N,8)*K16_8
          KASE = 0
C         Synchro point with broadcast of errors
 777      CONTINUE
          NB_BYTES_MAX = max(NB_BYTES_MAX,NB_BYTES)
          CALL MUMPS_PROPINFO( ICNTL(1), INFO(1),
     &                       id%COMM,id%MYID)
          IF ( INFO(1) .LT. 0 ) GOTO 90
C         TIMEEA needed if EA and IR with stopping criterium
C         and IR with fixed n.of steps. 
          TIMEEA = 0.0D0
C         TIMEEA1 needed if EA and IR with fixed n.of steps  
          TIMEEA1 = 0.0D0
          CALL MUMPS_SECDEB(TIMEIT)
C         -------------------------
C
C         RHSOL holds the initial guess for the solution
C         We start the loop on the Iterative refinement procedure
C
C
C
C           |-   IRefin.   L O O P   -|
C           V                         V
C
C  =========================================================
C    Computation of the infinity norm of A
C  =========================================================
          IF ((ICNTL11.GT.0).OR.(ICNTL10.GT.0)) THEN
C           We don't get through these lines if ICNTL10<=0 AND ICNTL11<=0
            IF ( KEEP(54) .eq. 0 ) THEN
C             ------------------
C             Centralized matrix
C             ------------------
              IF ( id%MYID .eq. MASTER ) THEN
C               -----------------------------------------
C               Call ZMUMPS_SOL_X outside, if needed,
C               in order to compute w(i,2)=sum|Aij|,j=1:n
C               in vector R_W(id%N+i)
C               -----------------------------------------
                IF (KEEP(55).NE.0) THEN
C                 unassembled matrix and norm of row required
                  CALL ZMUMPS_SOL_X_ELT(MTYPE, id%N,
     &            id%NELT, id%ELTPTR(1),
     &            id%LELTVAR, id%ELTVAR(1),
     &            id%NA_ELT, id%A_ELT(1),
     &            R_W(id%N+1), KEEP(1),KEEP8(1) )
                ELSE
C                 assembled matrix
                  IF ( MTYPE .eq. 1 ) THEN
                    CALL ZMUMPS_SOL_X
     &              ( id%A(1), id%NZ, id%N, id%IRN(1), id%JCN(1),
     &              R_W(id%N+1), KEEP(1),KEEP8(1))
                  ELSE
                    CALL ZMUMPS_SOL_X
     &              ( id%A(1), id%NZ, id%N, id%JCN(1), id%IRN(1),
     &              R_W(id%N+1), KEEP(1),KEEP8(1))
                  END IF
                ENDIF
              ENDIF
            ELSE
C             ---------------------
C             Matrix is distributed
C             ---------------------
              IF ( I_AM_SLAVE .and.
     &             id%NZ_loc .NE. 0 ) THEN
                IF ( MTYPE .eq. 1 ) THEN
                  CALL ZMUMPS_SOL_X(id%A_loc(1),
     &            id%NZ_loc, id%N,
     &            id%IRN_loc(1), id%JCN_loc(1),
     &            R_LOCWK54, id%KEEP(1),id%KEEP8(1) )
                ELSE
                  CALL ZMUMPS_SOL_X(id%A_loc(1),
     &            id%NZ_loc, id%N,
     &            id%JCN_loc(1), id%IRN_loc(1),
     &            R_LOCWK54, id%KEEP(1),id%KEEP8(1) )
                END IF
              ELSE
                R_LOCWK54 = RZERO
              END IF
C             -------------------------
C             Assemble result on master
C             -------------------------
              IF ( id%MYID .eq. MASTER ) THEN
                CALL MPI_REDUCE( R_LOCWK54, R_W( id%N + 1 ),
     &            id%N, MPI_DOUBLE_PRECISION,
     &            MPI_SUM,MASTER,id%COMM, IERR)
              ELSE
                CALL MPI_REDUCE( R_LOCWK54, R_DUMMY,
     &            id%N, MPI_DOUBLE_PRECISION,
     &            MPI_SUM,MASTER,id%COMM, IERR)
              END IF
C           End if KEEP(54)
            END IF
C
            IF ( id%MYID .eq. MASTER ) THEN
C             R_W is available on the master process only
              RINFOG(4) = dble(ZERO)
              DO I = 1, id%N
                RINFOG(4) = max(R_W( id%N +I), RINFOG(4))
              ENDDO
            ENDIF
C           end ICNTL11 =/0  v ICNTL10>0
          ENDIF
C  =========================================================
C    END norm of A
C  =========================================================
C         Initializations for the IR
          NOITER = 0
          IFLAG_IR = 0
          TESTConv = .FALSE.
C         Test of convergence should be made
          IF (( id%MYID .eq. MASTER ).AND.(ICNTL10.GT.0)) THEN
            TESTConv = .TRUE.
            ARRET = CNTL(2)
            IF (ARRET .LT. 0.0D0) THEN
              ARRET = sqrt(epsilon(0.0D0))
            END IF
          ENDIF
C  =========================================================
C         Starting IR
          DO  22 IRStep = 1, NITREF +1
C  =========================================================
C
C  =========================================================
C   Refine the solution starting from the second step of do loop
C  =========================================================
            IF (( id%MYID .eq. MASTER ).AND.(IRStep.GT.1)) THEN
              NOITER = NOITER + 1
              DO I = 1, id%N
                 id%RHS(IBEG+I-1) = id%RHS(IBEG+I-1) + C_Y(I)
              ENDDO
            ENDIF
C  ===========================================
C   Computation of the RESIDUAL and of |A||x|
C  ===========================================
            IF ( KEEP(54) .eq. 0 ) THEN
              IF ( id%MYID .eq. MASTER ) THEN
                IF (KEEP(55).NE.0) THEN
C                 input matrix by element
                  CALL ZMUMPS_ELTYD( MTYPE, id%N,
     &            id%NELT, id%ELTPTR(1), id%LELTVAR,
     &            id%ELTVAR(1), id%NA_ELT, id%A_ELT(1),
     &            SAVERHS, id%RHS(IBEG),
     &            C_Y, R_W, KEEP(50))
                ELSE
                  IF ( MTYPE .eq. 1 ) THEN
                   CALL ZMUMPS_SOL_Y(id%A(1), id%NZ, id%N, id%IRN(1),
     &                   id%JCN(1), SAVERHS,
     &                   id%RHS(IBEG), C_Y, R_W, KEEP(1),KEEP8(1))
                  ELSE
                   CALL ZMUMPS_SOL_Y(id%A(1), id%NZ, id%N, id%JCN(1),
     &                        id%IRN(1), SAVERHS,
     &                        id%RHS(IBEG), C_Y, R_W, KEEP(1),KEEP8(1))
                  ENDIF
                ENDIF
              ENDIF
            ELSE
C             ---------------------
C             Matrix is distributed
C             ---------------------
              CALL MPI_BCAST( RHS_IR(IBEG), id%N,
     &              MPI_DOUBLE_COMPLEX, MASTER,
     &              id%COMM, IERR )
C             --------------------------------------
C             Compute Y = SAVERHS - A * RHS
C             Y, SAVERHS defined only on master
C             --------------------------------------
              IF ( I_AM_SLAVE .and.
     &           id%NZ_loc .NE. 0 ) THEN
                CALL ZMUMPS_LOC_MV( id%N, id%NZ_loc,
     &          id%IRN_loc(1), id%JCN_loc(1), id%A_loc(1),
     &          RHS_IR(IBEG), C_LOCWK54, KEEP(50), MTYPE )
              ELSE
                C_LOCWK54 = ZERO
              END IF
              IF ( id%MYID .eq. MASTER ) THEN
                CALL MPI_REDUCE( C_LOCWK54, C_Y,
     &          id%N, MPI_DOUBLE_COMPLEX,
     &          MPI_SUM,MASTER,id%COMM, IERR)
C              ===========================
                C_Y = SAVERHS - C_Y
C              ===========================
              ELSE
                CALL MPI_REDUCE( C_LOCWK54, C_DUMMY,
     &          id%N, MPI_DOUBLE_COMPLEX,
     &          MPI_SUM,MASTER,id%COMM, IERR)
              END IF
C             --------------------------------------
C             Compute
C             * If MTYPE = 1
C                   W(i) = Sum | Aij | | RHSj |
C                           j
C             * If MTYPE = 0
C                   W(j) = Sum | Aij | | RHSi |
C                           i
C             R_LOCWK54 used as local array for W
C             RHS has been broadcasted
C             --------------------------------------
              IF ( I_AM_SLAVE .and. id%NZ_loc .NE. 0 ) THEN
                CALL ZMUMPS_LOC_OMEGA1( id%N, id%NZ_loc,
     &          id%IRN_loc(1), id%JCN_loc(1), id%A_loc(1),
     &          RHS_IR(IBEG), R_LOCWK54, KEEP(50), MTYPE )
              ELSE
                R_LOCWK54 = RZERO
              END IF
              IF ( id%MYID .eq. MASTER ) THEN
                CALL MPI_REDUCE( R_LOCWK54, R_W,
     &          id%N, MPI_DOUBLE_PRECISION,
     &          MPI_SUM,MASTER,id%COMM, IERR)
              ELSE
                CALL MPI_REDUCE( R_LOCWK54, R_DUMMY,
     &          id%N, MPI_DOUBLE_PRECISION,
     &          MPI_SUM, MASTER, id%COMM, IERR)
              ENDIF
            ENDIF
C  =====================================
C   END computation RESIDUAL and |A||x|
C  =====================================
            IF ( id%MYID .eq. MASTER ) THEN
C
              IF ((ICNTL11.GT.0).OR.(ICNTL10.GT.0)) THEN
C             --------------
C             Error analysis and test of convergence,
C             Compute the sparse componentwise backward error:
C               - at each step if test of convergence of IR is
C                 requested (ICNTL(10)>0)
C               - at step 1 and NITREF+1 if error analysis
C                 to be computed (ICNTL(11)>0) and if ICNTL(10)< 0
                IF (((ICNTL11.GT.0).OR.((ICNTL10.LT.0).AND.
     &               ((IRStep.EQ.1).OR.(IRStep.EQ.NITREF+1)))
     &               .OR.((ICNTL10.EQ.0).AND.(IRStep.EQ.1)))
     &                      .OR.(ICNTL10.GT.0)) THEN
C                 Compute w1 and w2
C                 always if ICNTL10>0 in the other case if ICNTL11>0
C                 -----------------
                  IF (ICNTL10.LT.0) CALL MUMPS_SECDEB(TIMEEA1)
                  CALL ZMUMPS_SOL_OMEGA(id%N,SAVERHS,
     &                id%RHS(IBEG), C_Y, R_W, C_W, IW1, IFLAG_IR,
     &                RINFOG(7), NOITER, TESTConv,
     &                MP, ARRET )
                  IF (ICNTL10.LT.0) THEN
                    CALL MUMPS_SECFIN(TIMEEA1)
                    id%DKEEP(120)=id%DKEEP(120)+TIMEEA1
                  ENDIF
                ENDIF
                IF ((ICNTL11.GT.0).AND.(
     &          (ICNTL10.LT.0.AND.(IRStep.EQ.1.OR.IRStep.EQ.NITREF+1))
     &          .OR.((ICNTL10.GE.0).AND.(IRStep.EQ.1))
     &                           )) THEN
C                 Error analysis before iterative refinement
C                 or for last if icntl10<0
C                 ------------------------------------------
                  CALL MUMPS_SECDEB(TIMEEA)
                  IF (ICNTL10.EQ.0) THEN
C                   No IR : there will be only the EA of the 1st sol.
                    IF ( MPG .GT. 0 ) WRITE( MPG, 170 )
                  ELSEIF (IRStep.EQ.1) THEN
C                   IR :  we print the EA of the 1st sol.
                    IF ( MPG .GT. 0 ) WRITE( MPG, 55 )
                  ELSEIF ((ICNTL10.LT.0).AND.(IRStep.EQ.NITREF+1)) THEN
C                   IR with fixed n. of steps:  we print the EA
C                                               of the last sol.
                    IF ( MPG .GT. 0 ) THEN
                      WRITE( MPG, 81 )
                      WRITE( MPG, * )
                      WRITE( MPG, 141 ) 
     &          'NUMBER OF STEPS OF ITERATIVE REFINEMENT REQUESTED =', 
     &          NOITER
                    ENDIF
                  ENDIF
                  GIVSOL = .TRUE.
                  CALL ZMUMPS_SOL_Q(MTYPE,INFO(1),id%N,id%NZ,
     &            id%RHS(IBEG),
     &            SAVERHS,R_W(id%N+1),C_Y,GIVSOL,
     &            RSOL,RINFOG(4),RINFOG(5),RINFOG(6),MPG,ICNTL(1),
     &            KEEP(1),KEEP8(1))
                  IF ( MPG .GT. 0 ) THEN
C                   Error analysis before iterative refinement
                    WRITE( MPG, 115 )
     &              'RINFOG(7):COMPONENTWISE SCALED RESIDUAL(W1)=',
     &              RINFOG(7)
                    WRITE( MPG, 115 )
     &              '------(8):---------------------------- (W2)=',
     &              RINFOG(8)
                  END IF
                  CALL MUMPS_SECFIN(TIMEEA)
                  id%DKEEP(120)=id%DKEEP(120)+TIMEEA
C                 end EA of the first solution
                END IF
              END IF
C             --------------
              IF (IRStep.EQ.NITREF +1) THEN
C               If we are at the NITREF+1 step , we have refined the
C               solution NITREF times so we have to stop.
                KASE = 0
C               If we test the convergence (ICNTL10.GT.0) and 
C               IFLAG_IR = 0 we set a warning : more than NITREF steps
C               needed
                IF ((ICNTL10.GT.0).AND.(IFLAG_IR.EQ.0))
     &             id%INFO(1) = id%INFO(1) + 8
              ELSE
                IF (ICNTL10.GT.0) THEN
C                 -------------------
C                 Results of the test of convergence.
C                 IFLAG_IR =  0 we should try to improve the solution
C                          =  1 the stopping criterium is satisfied
C                          =  2 the method is diverging, we go back
C                               to the previous iterate
C                          =  3 the convergence is too slow
                  IF (IFLAG_IR.GT.0) THEN
C                   If the convergence criterion is satisfied
C                   or the convergence too slow
C                   we set KASE=0 (end of the Iterative refinement)
                    KASE = 0
C                   If the convergence is not improved,
C                   we go back to the previous iterate.
C                   IFLAG_IR can be equal to 2 only if IRStep >= 2
                    IF (IFLAG_IR.EQ.2)  NOITER = NOITER - 1
                  ELSE
C                   IFLAG_IR=0, try to improve the solution
                    KASE = 2
                  ENDIF
                ELSEIF (ICNTL10.LT.0) THEN
C                 -------------------
                  KASE = 2
                ELSE
C                 ICNTL10 = 0, we want to perform only EA and not IR.
C                 -----------------
                  KASE = 0
                END IF
              ENDIF
C           End Master
            ENDIF
C           --------------
C           Broadcast KASE
C           --------------
            CALL MPI_BCAST( KASE, 1, MPI_INTEGER, MASTER,
     &                      id%COMM, IERR )
C           If Kase= 0 we quit the IR process
            IF (KASE.LE.0) GOTO 666
            IF (KASE.LT.0) THEN
              WRITE(*,*) "Internal error 17 in ZMUMPS_SOL_DRIVER"
            ENDIF
C  =========================================================
C   COMPUTE the solution of Ay = r
C  =========================================================
C           Call internal routine to avoid code duplication
            CALL ZMUMPS_PP_SOLVE()
            IF (INFO(1) .LT. 0) GOTO 90
C           -----------------------
C           Go back to beginning of
C           loop to apply next step
C           of iterative refinement
C           -----------------------
  22      CONTINUE
 666      CONTINUE
C         ************************************************
C
C         End of the iterative refinement procedure
C
C         ************************************************
          CALL MUMPS_SECFIN(TIMEIT)
          IF ( id%MYID .EQ. MASTER ) THEN
            IF ( NITREF .GT. 0 ) THEN
              id%INFOG(15) = NOITER
            END IF
C           id%DKEEP(114) time for the iterative refinement
C           id%DKEEP(120) time for the error analysis
C           id%DKEEP(121) time for condition number
C           these values are meaningful only on the host.
            IF (ICNTL10.EQ.0) THEN
C           No IR has been requested. All the time is needed
C           for computing EA
               id%DKEEP(120)=TIMEIT
            ELSE
C            IR has been requested
             id%DKEEP(114)=TIMEIT - id%DKEEP(120)
            ENDIF
          END IF
          IF ( PROKG ) THEN
              IF (ICNTL10.GT.0) THEN
                WRITE( MPG, 81 )
                WRITE( MPG, * )
                WRITE( MPG, 141 )
     &          'NUMBER OF STEPS OF ITERATIVE REFINEMENTS PERFORMED  =',
     &          NOITER
              ENDIF
          ENDIF
C
C         ==================================================
C         BEGIN
C         Perform error analysis after iterative refinement
C         ==================================================
          IF ((ICNTL11 .GT. 0).AND.(ICNTL10.GT.0)) THEN
C           If IR is requested with test of convergence,
C           the EA of the last step of IR is done here, 
C           otherwise EA of the last step is done at the 
C           end of IR
            CALL MUMPS_SECDEB(TIMEEA)
            KASE = 0
            IF (id%MYID .eq. MASTER ) THEN
C             Test if IFLAG_IR = 2, that is if the the IR was diverging,
C             we went back to the previous iterate
C             We have to do EA on the last computed solution.
              IF (IFLAG_IR.EQ.2) KASE = 2
            ENDIF
C           --------------
C           Broadcast KASE
C           --------------
            CALL MPI_BCAST( KASE, 1, MPI_INTEGER, MASTER,
     &      id%COMM, IERR )
            IF (KASE.EQ.2) THEN
C             We went back to the previous iterate
C             We have to do EA on the last computed solution.
C             Compute the residual in C_Y using IRN, JCN, ASPK
C             and the solution  RHS(IBEG)
C             The norm of the ith row in R_Y(I).
              IF ( KEEP(54) .eq. 0 ) THEN
C               ---------------------
C               Matrix is centralized
C               ---------------------
                IF (id%MYID .EQ. MASTER) THEN
                  IF (KEEP(55).EQ.0) THEN
                    CALL ZMUMPS_QD2( MTYPE, id%N, id%NZ, id%A(1),
     &              id%IRN(1), id%JCN(1),
     &              id%RHS(IBEG), SAVERHS, R_Y, C_Y, KEEP(1),KEEP8(1))
                  ELSE
                    CALL ZMUMPS_ELTQD2( MTYPE, id%N,
     &              id%NELT, id%ELTPTR(1),
     &              id%LELTVAR, id%ELTVAR(1),
     &              id%NA_ELT, id%A_ELT(1),
     &              id%RHS(IBEG), SAVERHS, R_Y, C_Y, KEEP(1),KEEP8(1))
                  ENDIF
                ENDIF
              ELSE
C               ---------------------
C               Matrix is distributed
C               ---------------------
                CALL MPI_BCAST( RHS_IR(IBEG), id%N,
     &              MPI_DOUBLE_COMPLEX, MASTER,
     &              id%COMM, IERR )
C               ----------------
C               Compute residual
C               ----------------
                IF ( I_AM_SLAVE .and.
     &            id%NZ_loc .NE. 0 ) THEN
                  CALL ZMUMPS_LOC_MV( id%N, id%NZ_loc,
     &            id%IRN_loc(1), id%JCN_loc(1), id%A_loc(1),
     &            RHS_IR(IBEG), C_LOCWK54, KEEP(50), MTYPE )
                ELSE
                  C_LOCWK54 = ZERO
                END IF
                IF ( id%MYID .eq. MASTER ) THEN
                  CALL MPI_REDUCE( C_LOCWK54, C_Y,
     &            id%N, MPI_DOUBLE_COMPLEX,
     &            MPI_SUM,MASTER,id%COMM, IERR)
                  C_Y = SAVERHS - C_Y
                ELSE
                  CALL MPI_REDUCE( C_LOCWK54, C_DUMMY,
     &            id%N, MPI_DOUBLE_COMPLEX,
     &            MPI_SUM,MASTER,id%COMM, IERR)
                END IF
              ENDIF
            ENDIF  ! KASE.EQ.2
            IF (id%MYID .EQ. MASTER) THEN
C             Compute which equations are associated to w1 and which
C             ones are associated to w2 in case of IFLAG_IR=2.
C             If IFLAG_IR = 0 or 1 IW1 should be correct
              IF (IFLAG_IR.EQ.2) THEN
                TESTConv = .FALSE.
                CALL ZMUMPS_SOL_OMEGA(id%N,SAVERHS,
     &              id%RHS(IBEG), C_Y, R_W, C_W, IW1, IFLAG_IR,
     &              RINFOG(7), 0, TESTConv,
     &              MP, ARRET )
              ENDIF ! (IFLAG_IR.EQ.2)
c             Compute some statistics for
              GIVSOL = .TRUE.
              CALL ZMUMPS_SOL_Q(MTYPE,INFO(1),id%N,id%NZ,
     &        id%RHS(IBEG),
     &        SAVERHS,R_W(id%N+1),C_Y,GIVSOL,
     &        RSOL,RINFOG(4),RINFOG(5),RINFOG(6),MPG,ICNTL(1),
     &        KEEP(1),KEEP8(1))
            ENDIF ! Master
            CALL MUMPS_SECFIN(TIMEEA)
            id%DKEEP(120)=id%DKEEP(120)+TIMEEA
          ENDIF ! ICNTL11>0 and ICNTL10>0
C  =========================================================
C   Compute the Condition number associated if requested.
C  =========================================================
          CALL MUMPS_SECDEB(TIMELCOND)
          IF (ICNTL11 .EQ. 1) THEN
            IF ( id%MYID .eq. MASTER ) THEN
C             Notice that D is always the identity
              ALLOCATE( D(id%N),stat =allocok )
              IF ( allocok .GT. 0 ) THEN
                INFO(1)=-13
                INFO(2)=id%N
                GOTO 777
              ENDIF
              NB_BYTES = NB_BYTES + int(id%N,8)*K16_8
              DO I = 1, id%N
                D( I ) = RONE
              END DO
            ENDIF
            KASE = 0
 222        CONTINUE
            IF ( id%MYID .EQ. MASTER ) THEN
              CALL ZMUMPS_SOL_LCOND(id%N, SAVERHS,
     &        id%RHS(IBEG), C_Y, D, R_W, C_W, IW1, KASE,
     &        RINFOG(7), RINFOG(9), RINFOG(10),
     &        MP, KEEP(1),KEEP8(1))
            ENDIF
C           --------------
C           Broadcast KASE
C           --------------
            CALL MPI_BCAST( KASE, 1, MPI_INTEGER, MASTER,
     &                      id%COMM, IERR )
C           KASE <= 0
C           We reach the end of iterative method to compute
C           LCOND1 and LCOND2
            IF (KASE.LE.0) GOTO 224
            CALL ZMUMPS_PP_SOLVE()
            IF (INFO(1) .LT. 0) GOTO 90
C           ---------------------------
C           Go back to beginning of
C           loop to apply next step
C           of iterative method
C           -----------------------
            GO TO 222
C    End ICNTL11 = 1
          ENDIF
 224      CONTINUE
          CALL MUMPS_SECFIN(TIMELCOND)
          id%DKEEP(121)=id%DKEEP(121)+TIMELCOND
          IF ((id%MYID .EQ. MASTER).AND.(ICNTL11.GT.0)) THEN
            IF (ICNTL10.GT.0) THEN
C             If ICNTL10<0 these stats have been printed before IR
              IF ( MPG .GT. 0 ) THEN
                WRITE( MPG, 115 )
     &          'RINFOG(7):COMPONENTWISE SCALED RESIDUAL(W1)=',
     &          RINFOG(7)
                WRITE( MPG, 115 )
     &          '------(8):---------------------------- (W2)=',
     &          RINFOG(8)
              ENDIF
            END IF
            IF (ICNTL11.EQ.1) THEN
C           If ICNTL11/=1 these stats haven't been computed
              IF (MPG.GT.0) THEN
               WRITE( MPG, 115 )
     &         '------(9):Upper bound ERROR ...............=',
     &         RINFOG(9)
               WRITE( MPG, 115 )
     &         '-----(10):CONDITION NUMBER (1) ............=',
     &         RINFOG(10)
               WRITE( MPG, 115 )
     &         '-----(11):CONDITION NUMBER (2) ............=',
     &         RINFOG(11)
              END IF
            END IF
          END IF ! MASTER && ICNTL11.GT.0
          IF ( PROKG .AND. abs(ICNTL10) .GT.0 ) WRITE( MPG, 131 )
C===================================================
C Perform error analysis after iterative refinements
C END
C===================================================
C
          IF (id%MYID == MASTER) THEN
            NB_BYTES = NB_BYTES - int(size(C_W),8)*K35_8
            DEALLOCATE(C_W)
            NB_BYTES = NB_BYTES - int(size(R_W),8)*K16_8
     &                        - int(size(IW1),8)*K34_8
            DEALLOCATE(R_W)
            DEALLOCATE(IW1)
            IF (ICNTL11 .EQ. 1) THEN
C             We have used D only for LCOND1,2
              DEALLOCATE(D)
              NB_BYTES = NB_BYTES - int(size(D  ),8)*K16_8
            ENDIF
          ENDIF
          NB_BYTES = NB_BYTES -
     &     (int(size(R_Y),8)+int(size(R_LOCWK54),8))*K16_8
          NB_BYTES = NB_BYTES -
     &     (int(size(C_Y),8)+int(size(C_LOCWK54),8))*K35_8
          DEALLOCATE(R_Y)
          DEALLOCATE(C_Y)
          DEALLOCATE(R_LOCWK54)
          DEALLOCATE(C_LOCWK54)
C End ERANA
        END IF
C============================================
C
C  ITERATIVE REFINEMENT AND/OR ERROR ANALYSIS
C
C  END
C
C============================================
C       ==========================
C       Begin reordering on master
C       corresponding to maximum transversal permutation
C       in case of centralized solution
C       (ICNTL21==0)
C
        IF ( id%MYID .EQ. MASTER .AND. ICNTL21==0
     &     .AND. KEEP(23) .NE. 0.AND.KEEP(237).EQ.0) THEN
C         ((No transpose and backward performed and NO A-1)
C         or null space computation): permutation
C         must be done on solution.
          IF ((KEEP(221).NE.1 .AND. MTYPE .EQ. 1)
     &       .OR. KEEP(111) .NE.0 .OR. KEEP(252).NE.0 ) THEN
C           Permute the solution RHS according to the column
C           permutation held in UNS_PERM
C           Column J of the permuted matrix corresponds to
C           column UNS_PERM(J) of the original matrix.
C           RHS holds the permuted solution
C           Note that id%N>1 since KEEP(23)=0 when id%N=1
C
C
            ALLOCATE( C_RW1( id%N ),stat =allocok )
!           temporary not in NB_BYTES
            IF ( allocok .GT. 0 ) THEN
              INFO(1)=-13
              INFO(2)=id%N
              WRITE(*,*) 'could not allocate ', id%N, 'integers.'
              CALL MUMPS_ABORT()
            END IF
            DO K = 1, NBRHS_EFF
              KDEC = (K-1)*LD_RHS+IBEG-1
              DO I = 1, id%N
                C_RW1(I) = id%RHS(KDEC+I)
              ENDDO
              DO I = 1, id%N
               JPERM = id%UNS_PERM(I)
               id%RHS( KDEC+JPERM ) = C_RW1( I )
              ENDDO
            ENDDO
            DEALLOCATE( C_RW1 ) !temporary not in NB_BYTES
          END IF
        END IF
C
C  End reordering on master
C  ========================
        IF (id%MYID.EQ.MASTER .and.ICNTL21==0.and.KEEP(221).NE.1.AND.
     &     (KEEP(237).EQ.0) ) THEN
*         print out the solution
          IF ( INFO(1) .GE. 0 .AND. ICNTL(4).GE.3 .AND. ICNTL(3).GT.0)
     &    THEN
            K = min0(10, id%N)
            IF (ICNTL(4) .eq. 4 ) K = id%N
            J = min0(10,NBRHS_EFF)
            IF (ICNTL(4) .eq. 4 ) J = NBRHS_EFF
            DO II=1, J
              WRITE(ICNTL(3),110) BEG_RHS+II-1
              WRITE(ICNTL(3),160)
     &      (id%RHS(IBEG+(II-1)*LD_RHS+I-1),I=1,K)
            ENDDO
          END IF
        END IF
#if defined(RHSCOMP_BYROWS)
 1000   CONTINUE
#endif
C     ==========================
C     blocking for multiple RHS (END OF DO WHILE (BEG_RHS.LE.NBRHS)
        IF ((KEEP(248).EQ.1).AND.(KEEP(237).EQ.0)) THEN
          ! case of general sparse: in case of empty columns
          ! NBRHS_EFF might has been updated and broadcasted
          ! and holds the effective size of a contiguous block of
          ! non empty columns
          BEG_RHS = BEG_RHS + NBRHS_EFF ! nb of nonempty columns
        ELSE
          BEG_RHS = BEG_RHS + NBRHS
        ENDIF
      ENDDO
C     DO WHILE (BEG_RHS.LE.id%NRHS)
C     ==========================
C
C     ========================================================
C     Reset RHS to zero for all remaining columns that
C     have not been processed because they were emtpy
C     ========================================================
      IF (   (id%MYID.EQ.MASTER)
     &       .AND. ( KEEP(248).NE.0 )  ! sparse RHS on input
     &       .AND. ( KEEP(237).EQ.0 )  ! No A-1
     &       .AND. ( ICNTL21.EQ.0 )    ! Centralized solution
     &       .AND. ( KEEP(221) .NE.1 ) ! Not Reduced RHS step of Schur
     &       .AND. ( JEND_RHS .LT. id%NRHS )
     &   )
     &         THEN
        JBEG_NEW = JEND_RHS + 1
        IF (DO_PERMUTE_RHS.OR.INTERLEAVE_PAR) THEN
           DO WHILE ( JBEG_NEW.LE. id%NRHS)
              DO I=1, id%N
                 id%RHS((PERM_RHS(JBEG_NEW) -1)*LD_RHS+I)
     &                     = ZERO
              ENDDO
              JBEG_NEW = JBEG_NEW +1
              CYCLE
           ENDDO
        ELSE
          DO WHILE ( JBEG_NEW.LE. id%NRHS)
            DO I=1, id%N
                id%RHS((JBEG_NEW -1)*LD_RHS + I) = ZERO
            ENDDO
            JBEG_NEW = JBEG_NEW +1
          ENDDO
        ENDIF                   ! End DO_PERMUTE_RHS.OR.INTERLEAVE_PAR
      ENDIF
C     ========================================================
C     Reset id%SOL_loc to zero for all remaining columns that
C     have not been processed because they were emtpy
C     ========================================================
      IF ( I_AM_SLAVE .AND. (ICNTL21.NE.0) .AND.
     &      ( JEND_RHS .LT. id%NRHS ) .AND. KEEP(221).NE.1 ) THEN
         JBEG_NEW = JEND_RHS + 1
C
           DO WHILE ( JBEG_NEW.LE. id%NRHS)
            DO I=1, KEEP(89)
                id%SOL_loc((JBEG_NEW -1)*id%LSOL_loc + I) = ZERO
            ENDDO
            JBEG_NEW = JBEG_NEW +1
           ENDDO
      ENDIF
C
C     =====================================================================
C     Reset id%RHSCOMP and id%REDRHS to zero for all remaining columns that
C     have not been processed because they were emtpy
C     ========================================================
      IF ((KEEP(221).EQ.1) .AND.
     &        ( JEND_RHS .LT. id%NRHS ) ) THEN
       IF (id%MYID .EQ. MASTER) THEN
           JBEG_NEW = JEND_RHS + 1
           DO WHILE ( JBEG_NEW.LE. id%NRHS)
            DO I=1,  id%SIZE_SCHUR
              id%REDRHS((JBEG_NEW -1)*LD_REDRHS + I) =  ZERO
            ENDDO
            JBEG_NEW = JBEG_NEW +1
           ENDDO
       ENDIF
       IF (I_AM_SLAVE) THEN
#if defined(RHSCOMP_BYROWS)
           DO I=1,NBENT_RHSCOMP
            JBEG_NEW = JEND_RHS + 1
            DO WHILE ( JBEG_NEW.LE. id%NRHS)
              id%RHSCOMP(JBEG_NEW + (I-1)*NBRHS_EFF) =  ZERO
              JBEG_NEW = JBEG_NEW +1
            ENDDO
           ENDDO
#else
           JBEG_NEW = JEND_RHS + 1
           DO WHILE ( JBEG_NEW.LE. id%NRHS)
            DO I=1,NBENT_RHSCOMP
              id%RHSCOMP((JBEG_NEW -1)*LD_RHSCOMP + I) =  ZERO
            ENDDO
            JBEG_NEW = JBEG_NEW +1
           ENDDO
#endif
       ENDIF
      ENDIF
C
C
C     ! maximum size used on that proc
      id%INFO(26) = int(NB_BYTES_MAX / 1000000_8, 4)
C     Centralize memory statistics on the host
C
C       INFOG(30) = size of mem in bytes for solve
C                   for the processor using largest memory
C       INFOG(31) = size of mem in bytes for solve
C                   sum over all processors
C     ----------------------------------------------------
      CALL MUMPS_MEM_CENTRALIZE( id%MYID, id%COMM,
     &                           id%INFO(26), id%INFOG(30), IRANK )
      IF ( PROKG ) THEN
        WRITE( MPG,'(A,I10) ')
     &  ' ** Rank of processor needing largest memory in solve     :',
     &  IRANK
        WRITE( MPG,'(A,I10) ')
     &  ' ** Space in MBYTES used by this processor for solve      :',
     &  id%INFOG(30)
        IF ( KEEP(46) .eq. 0 ) THEN
        WRITE( MPG,'(A,I10) ')
     &  ' ** Avg. Space in MBYTES per working proc during solve    :',
     &  ( id%INFOG(31)-id%INFO(26) ) / id%NSLAVES
        ELSE
        WRITE( MPG,'(A,I10) ')
     &  ' ** Avg. Space in MBYTES per working proc during solve    :',
     &  id%INFOG(31) / id%NSLAVES
        END IF
      END IF
*===============================
*End of Solve Phase
*===============================
C  Store and print timings
      CALL MUMPS_SECFIN(TIME3)
      id%DKEEP(112)=TIME3
      id%DKEEP(113)=TIMEC2
      id%DKEEP(115)=TIMESCATTER2
      id%DKEEP(116)=TIMEGATHER2
      IF (PROKG) THEN
        WRITE ( MPG, *)
        WRITE( MPG, 434 ) id%DKEEP(115)
        WRITE( MPG, 432 ) id%DKEEP(113)
        WRITE( MPG, 435 ) id%DKEEP(117)
        IF ((KEEP(38).NE.0).OR.(KEEP(20).NE.0))
     &     WRITE( MPG, 437 ) id%DKEEP(119)
        WRITE( MPG, 436 ) id%DKEEP(118)
        WRITE( MPG, 433 ) id%DKEEP(116)
      ELSE IF ( id%MYID .NE. MASTER .and. PROK ) THEN
        WRITE ( MP, *)
        WRITE( MP, 434 ) id%DKEEP(115)
        WRITE( MP, 432 ) id%DKEEP(113)
        WRITE( MP, 435 ) id%DKEEP(117)
        IF ((KEEP(38).NE.0).OR.(KEEP(20).NE.0))
     &     WRITE( MP, 437 ) id%DKEEP(119)
        WRITE( MP, 436 ) id%DKEEP(118)
        WRITE( MP, 433 ) id%DKEEP(116)
      END IF
 90   CONTINUE
      IF (INFO(1) .LT.0 ) THEN
      ENDIF
      IF (KEEP(201).GT.0)THEN
        IF (IS_INIT_OOC_DONE) THEN
          CALL ZMUMPS_OOC_END_SOLVE(IERR)
          IF (IERR.LT.0 .AND. INFO(1) .GE. 0) INFO(1) = IERR
        ENDIF
        CALL MUMPS_PROPINFO( ICNTL(1), INFO(1),
     &         id%COMM,id%MYID)
      ENDIF
C     ------------------------
C     Check allocation before
C     to deallocate (cases of
C     errors that could happen
C     before or after allocate
C     statement)
C
C       Sparse RHS
C       Free space and reset pointers if needed
        IF (IRHS_SPARSE_COPY_ALLOCATED) THEN
             NB_BYTES =  NB_BYTES -
     &        int(size(IRHS_SPARSE_COPY),8)*K34_8
             DEALLOCATE(IRHS_SPARSE_COPY)
             IRHS_SPARSE_COPY_ALLOCATED=.FALSE.
             NULLIFY(IRHS_SPARSE_COPY)
       ENDIF
       IF (IRHS_PTR_COPY_ALLOCATED) THEN
             NB_BYTES =  NB_BYTES -
     &        int(size(IRHS_PTR_COPY),8)*K34_8
             DEALLOCATE(IRHS_PTR_COPY)
             IRHS_PTR_COPY_ALLOCATED=.FALSE.
             NULLIFY(IRHS_PTR_COPY)
       ENDIF
       IF (RHS_SPARSE_COPY_ALLOCATED) THEN
             NB_BYTES =  NB_BYTES -
     &        int(size(RHS_SPARSE_COPY),8)*K35_8
             DEALLOCATE(RHS_SPARSE_COPY)
             RHS_SPARSE_COPY_ALLOCATED=.FALSE.
             NULLIFY(RHS_SPARSE_COPY)
       ENDIF
      IF (allocated(PERM_RHS)) THEN
        NB_BYTES = NB_BYTES - int(size(PERM_RHS),8)*K34_8
        DEALLOCATE(PERM_RHS)
      ENDIF
C     END A-1
      IF (allocated(UNS_PERM_INV)) THEN
        NB_BYTES = NB_BYTES - int(size(UNS_PERM_INV),8)*K34_8
        DEALLOCATE(UNS_PERM_INV)
      ENDIF
      IF (associated(id%BUFR)) THEN
          NB_BYTES = NB_BYTES - int(size(id%BUFR),8)*K34_8
          DEALLOCATE(id%BUFR)
          NULLIFY(id%BUFR)
      ENDIF
      IF ( I_AM_SLAVE ) THEN
        IF (allocated(RHS_BOUNDS)) THEN
          NB_BYTES = NB_BYTES -
     &          int(size(RHS_BOUNDS),8)*K34_8
          DEALLOCATE(RHS_BOUNDS)
        ENDIF
        IF (allocated(IWK_SOLVE)) THEN
          NB_BYTES = NB_BYTES - int(size(IWK_SOLVE),8)*K34_8
          DEALLOCATE( IWK_SOLVE )
        ENDIF
        IF (allocated(IWCB)) THEN
          NB_BYTES = NB_BYTES - int(size(IWCB),8)*K34_8
          DEALLOCATE( IWCB )
        ENDIF
C       ------------------------
C       SLAVE CODE
C       -----------------------
C       Deallocate send buffers
C       -----------------------
        CALL ZMUMPS_BUF_DEALL_CB( IERR )
        CALL ZMUMPS_BUF_DEALL_SMALL_BUF( IERR )
      END IF
C
      IF ( id%MYID .eq. MASTER ) THEN
C       ------------------------
C       SAVERHS may  have been
C       allocated only on master
C       ------------------------
        IF (allocated(SAVERHS)) THEN
         NB_BYTES = NB_BYTES - int(size(SAVERHS),8)*K35_8
         DEALLOCATE( SAVERHS)
        ENDIF
C       Nullify RHS_IR might have been pointing to id%RHS
        NULLIFY(RHS_IR)
      ELSE
C       --------------------
C       Free right-hand-side
C       on slave processors
C       --------------------
        IF (associated(RHS_IR)) THEN
          NB_BYTES = NB_BYTES - int(size(RHS_IR),8)*K35_8
          DEALLOCATE(RHS_IR)
          NULLIFY(RHS_IR)
        END IF
      END IF
      IF (I_AM_SLAVE) THEN
C       Deallocate temporary workspace SRW3
        IF (allocated(SRW3)) THEN
          NB_BYTES = NB_BYTES - int(size(SRW3),8)*K35_8
          DEALLOCATE(SRW3)
        ENDIF
        IF (LSCAL .AND. ICNTL21==1) THEN
C         Free local scaling arrays
          NB_BYTES = NB_BYTES -
     &              int(size(scaling_data%SCALING_LOC),8)*K16_8
          DEALLOCATE(scaling_data%SCALING_LOC)
          NULLIFY(scaling_data%SCALING_LOC)
        ENDIF
C       Free memory until next call to ZMUMPS
        IF (WK_USER_PROVIDED) THEN
C         S points to WK_USER provided by user
C         KEEP8(24) holds size of WK_USER
C         it should be saved and is used
C         in incore to check that size provided is consistent
C         (see error -41)
          NULLIFY(id%S)
        ELSE IF (associated(id%S).AND.KEEP(201).GT.0) THEN
C         OOC: free space for S that was allocated
          NB_BYTES = NB_BYTES - KEEP8(23)*K35_8
          id%KEEP8(23)=0_8
          DEALLOCATE(id%S)
          NULLIFY(id%S)
        ENDIF
        IF (KEEP(221).NE.1) THEN
C       -- After reduction of RHS to Schur variables
C       -- keep compressed RHS generated during FWD step
C       -- to be used for future expansion
         IF (associated(id%RHSCOMP)) THEN
            NB_BYTES = NB_BYTES - int(size(id%RHSCOMP),8)*K35_8
            DEALLOCATE(id%RHSCOMP)
            NULLIFY(id%RHSCOMP)
         ENDIF
         IF (associated(id%POSINRHSCOMP_ROW)) THEN
            NB_BYTES = NB_BYTES -
     &                 int(size(id%POSINRHSCOMP_ROW),8)*K34_8
            DEALLOCATE(id%POSINRHSCOMP_ROW)
            NULLIFY(id%POSINRHSCOMP_ROW)
         ENDIF
         IF (id%POSINRHSCOMP_COL_ALLOC) THEN
           ! IF (associated(id%POSINRHSCOMP_COL)) THEN
            NB_BYTES = NB_BYTES -
     &                 int(size(id%POSINRHSCOMP_COL),8)*K34_8
            DEALLOCATE(id%POSINRHSCOMP_COL)
            NULLIFY(id%POSINRHSCOMP_COL)
           ! ENDIF
          id%POSINRHSCOMP_COL_ALLOC = .FALSE.
         ENDIF
        ENDIF
        IF ( WORK_WCB_ALLOCATED ) THEN
          NB_BYTES = NB_BYTES - int(size(WORK_WCB),8)*K35_8
          DEALLOCATE( WORK_WCB )
        ENDIF
C       Otherwise, WORK_WCB may point to some
C       position inside id%S, nullify it
        NULLIFY( WORK_WCB )
      ENDIF
      RETURN
 55   FORMAT (//' ERROR ANALYSIS BEFORE ITERATIVE REFINEMENT')
 100  FORMAT(//' ****** SOLVE & CHECK STEP ********'/)
 110  FORMAT (//' VECTOR SOLUTION FOR COLUMN ',I12)
 115  FORMAT(1X, A44,1P,D9.2)
 434  FORMAT(' TIME to scatter right-hand-sides =',F15.6)
 432  FORMAT(' TIME in solution step (fwd/bwd)  =',F15.6)
 435  FORMAT('  .. TIME in forward (fwd) step   =   ',F15.6)
 437  FORMAT('  .. TIME in ScaLAPACK root       =   ',F15.6)
 436  FORMAT('  .. TIME in backward (bwd) step  =   ',F15.6)
 433  FORMAT(' TIME to gather solution          =',F15.6)
 150  FORMAT (/' STATISTICS PRIOR SOLVE PHASE     ...........'/
     &        ' NUMBER OF RIGHT-HAND-SIDES                    =',I12/
     &        ' BLOCKING FACTOR FOR MULTIPLE RHS              =',I12/
     &        ' ICNTL (9)                                     =',I12/
     &        '  --- (10)                                     =',I12/
     &        '  --- (11)                                     =',I12/
     &        '  --- (20)                                     =',I12/
     &        '  --- (21)                                     =',I12/
     &        '  --- (30)                                     =',I12)
 151  FORMAT ('  --- (25)                                     =',I12)
 152  FORMAT ('  --- (26)                                     =',I12)
 153  FORMAT ('  --- (32)                                     =',I12)
 160  FORMAT (' RHS'/(1X,1P,5D14.6))
 170  FORMAT (//' ERROR ANALYSIS' )
 240  FORMAT (1X, A42,I4)
 270  FORMAT (//' BEGIN ITERATIVE REFINEMENT' )
  81  FORMAT (/' STATISTICS AFTER ITERATIVE REFINEMENT ')
 131  FORMAT (/' END   ITERATIVE REFINEMENT ')
 141  FORMAT(1X, A52,I4)
      CONTAINS
        SUBROUTINE ZMUMPS_PP_SOLVE()
        IMPLICIT NONE
C
C       Purpose:
C       =======
C       Scatter right-hand side, solve the system,
C       and gather the solution on the host during
C       post-processing.
C       We use an internal subroutine to avoid code
C       duplication without the complication of adding
C       new parameters or local variables. All variables
C       in this routine have the scope of ZMUMPS_SOL_DRIVER.
C
C
        IF (KASE .NE. 1 .AND. KASE .NE. 2) THEN
          WRITE(*,*) "Internal error 1 in ZMUMPS_PP_SOLVE"
          CALL MUMPS_ABORT()
        ENDIF
        IF ( id%MYID .eq. MASTER ) THEN
C         Define matrix B as follows:
C            MTYPE=1 => B=A other values B=At
C         The user asked to solve the system Bx=b
C
C         THEN
C           KASE = 1........ RW1 = INV(TRANSPOSE(B)) * RW1
C           KASE = 2........ RW1 = INV(B) * RW1
          IF ( MTYPE .EQ. 1 ) THEN
            SOLVET = KASE - 1
          ELSE
            SOLVET = KASE
          END IF
C         SOLVET= 1 -> solve A x = B, other values solve Atx=b
C         We force SOLVET to have value either 0 or 1, in order
C         to be able to test both values, and also, be able to
C         test whether SOLVET = MTYPE or not.
          IF ( SOLVET.EQ.2 ) SOLVET = 0
          IF ( LSCAL ) THEN
            IF ( SOLVET .EQ. 1 ) THEN
C             Apply rowscaling
              DO K = 1, id%N
                C_Y( K ) = C_Y( K ) * id%ROWSCA( K )
              END DO
            ELSE
C             Apply column scaling
              DO K = 1, id%N
                C_Y( K ) = C_Y( K ) * id%COLSCA( K )
              END DO
            END IF
          END IF
        END IF ! MYID.EQ.MASTER
C       ------------------------------
C       Broadcast SOLVET to the slaves
C       ------------------------------
        CALL MPI_BCAST( SOLVET, 1, MPI_INTEGER, MASTER,
     &                  id%COMM, IERR)
C       --------------------------------------------
C       Scatter the right hand side C_Y on all procs
C       --------------------------------------------
        IF ( .NOT.I_AM_SLAVE ) THEN
C         -- Master not working
          CALL ZMUMPS_SCATTER_RHS(id%NSLAVES,id%N, id%MYID,
     &      id%COMM,
     &      SOLVET, C_Y(1), id%N, 1,
     &      1,
     &      C_DUMMY, 1, 1,
     &      IDUMMY, 0,
     &      JDUMMY, id%KEEP(1), id%KEEP8(1), id%PROCNODE_STEPS(1),
     &      IDUMMY, 1,
     &      id%STEP(1),
     &      id%ICNTL(1),id%INFO(1))
        ELSE
          IF (SOLVET.EQ.MTYPE) THEN
C           POSINRHSCOMP_ROW is with respect to the
C           original linear system (transposed or not)
            PTR_POSINRHSCOMP_FWD => id%POSINRHSCOMP_ROW
          ELSE
C           Transposed, use column indices of original
C           system (ie, col indices of A or A^T)
            PTR_POSINRHSCOMP_FWD => id%POSINRHSCOMP_COL
          ENDIF
          LIW_PASSED = max( LIW, 1 )
          CALL ZMUMPS_SCATTER_RHS(id%NSLAVES,id%N, id%MYID,
     &      id%COMM,
     &      SOLVET,  C_Y(1), id%N, 1,
     &      1,
     &      id%RHSCOMP(IBEG_RHSCOMP), LD_RHSCOMP, 1,
     &      PTR_POSINRHSCOMP_FWD(1), NB_FS_IN_RHSCOMP_F,
C
     &      id%PTLUST_S(1), id%KEEP(1), id%KEEP8(1),
     &      id%PROCNODE_STEPS(1),
     &      IS(1), LIW_PASSED,
     &      id%STEP(1),
     &      id%ICNTL(1),id%INFO(1))
        ENDIF
        IF (INFO(1).LT.0) GOTO 89
C
C       Solve the system
C
        IF ( I_AM_SLAVE ) THEN
          LIW_PASSED = max( LIW, 1 )
          LA_PASSED = max( LA, 1_8 )
          IF (SOLVET.EQ.MTYPE) THEN
            PTR_POSINRHSCOMP_FWD => id%POSINRHSCOMP_ROW
            PTR_POSINRHSCOMP_BWD => id%POSINRHSCOMP_COL
          ELSE
            PTR_POSINRHSCOMP_FWD => id%POSINRHSCOMP_COL
            PTR_POSINRHSCOMP_BWD => id%POSINRHSCOMP_ROW
          ENDIF
          FROM_PP=.TRUE.
          NBSPARSE_LOC = .FALSE.
          CALL ZMUMPS_SOL_C( id%root, id%N, id%S(1), LA_PASSED,
     & id%IS(1), LIW_PASSED, WORK_WCB(1), LWCB_SOL_C, IWCB, LIWCB,
     & NBRHS_EFF, id%NA(1), id%LNA, id%NE_STEPS(1), SRW3, SOLVET,
     & ICNTL(1), FROM_PP, id%STEP(1), id%FRERE_STEPS(1),
     & id%DAD_STEPS(1), id%FILS(1),
     & id%PTLUST_S(1), id%PTRFAC(1), IWK_SOLVE(1), LIWK_SOLVE,
     & id%PROCNODE_STEPS(1), id%NSLAVES, INFO(1), KEEP(1), KEEP8(1),
     & id%DKEEP(1),
     & id%COMM_NODES, id%MYID, id%MYID_NODES, id%BUFR(1), id%LBUFR,
     & id%LBUFR_BYTES, id%ISTEP_TO_INIV2(1), id%TAB_POS_IN_PERE(1,1),
C         Next 3 arguments are not used in this call
     & IBEG_ROOT_DEF, IEND_ROOT_DEF, IROOT_DEF_RHS_COL1,
C         Case of special root node
     &    PTR_RHS_ROOT(1), LPTR_RHS_ROOT, SIZE_ROOT, MASTER_ROOT,
     &    id%RHSCOMP(IBEG_RHSCOMP), LD_RHSCOMP,
     &    PTR_POSINRHSCOMP_FWD(1), PTR_POSINRHSCOMP_BWD(1)
     &    , 1            , 1       ,  1
     &    , 1
     &    , IDUMMY, 1, JDUMMY, KDUMMY, 1, LDUMMY, 1, MDUMMY
     &    , 1, 1
     &    , NBSPARSE_LOC, PTR_RHS_BOUNDS(1), LPTR_RHS_BOUNDS
     &    )
        END IF
C       ------------------
C       Change error codes
C       ------------------
        IF (INFO(1).eq.-2) INFO(1)=-12
        IF (INFO(1).eq.-3) INFO(1)=-15
C
        IF ( INFO(1) .GE. 0) THEN
C         We need a workspace of minimal size KEEP(247)
C         in order to unpack pieces of the solution during
C         ZMUMPS_GATHER_SOLUTION below
C         - Avoid allocation if error already occurred.
C         - DEALLOCATE called after GATHER_SOLUTION
C         - in case of error with jump to 90 below, CWORK
C           should be automatically deallocated (F95 standard)
C         CWORK not needed for AM1
          ALLOCATE( CWORK(
     &          max(max(KEEP(247),KEEP(246))*NBRHS_EFF,1)),
     &          stat=allocok)
C
          IF (allocok > 0) THEN
C           Try smaller (it will be less efficient)
            ALLOCATE( CWORK(max(max(KEEP(247),KEEP(246)),1)),
     &             stat=allocok)
            IF (allocok > 0) THEN
              INFO(1)=-13
              INFO(2)=max(max(KEEP(247),KEEP(246)),1)
            ENDIF
          ENDIF
        ENDIF
C       -------------------------
C       Propagate possible errors
C       -------------------------
 89     CALL MUMPS_PROPINFO( ICNTL(1), INFO(1),
     &                   id%COMM,id%MYID)
C
C       Return in case of error.
        IF (INFO(1).LT.0) RETURN
C       -------------------------------
C       Assemble the solution on master
C       -------------------------------
C       (Note: currently, if this part of code is executed,
C       then necessarily NBRHS_EFF = 1)
C
C       === GATHER and SCALE solution ==============
C
        IF ((id%MYID.NE.MASTER).OR. .NOT.LSCAL) THEN
          PT_SCALING => Dummy_SCAL
        ELSE
          IF (SOLVET.EQ.1) THEN
            PT_SCALING => id%COLSCA
          ELSE
            PT_SCALING => id%ROWSCA
          ENDIF
        ENDIF
        LIW_PASSED = max( LIW, 1 )
C       Solution computed during ZMUMPS_SOL_C has been stored
C       in id%RHSCOMP and is gathered on the master in C_Y
        IF ( .NOT. I_AM_SLAVE ) THEN
C         I did not participate to computing part of the solution
C         (id%RHSCOMP not set/allocate) : receive solution, store
C         it and scale it.
          CALL ZMUMPS_GATHER_SOLUTION(id%NSLAVES,id%N,
     &      id%MYID, id%COMM, NBRHS_EFF,
     &      SOLVET, C_Y, id%N, NBRHS_EFF,
     &      JDUMMY, id%KEEP(1),id%KEEP8(1), id%PROCNODE_STEPS(1),
     &      IDUMMY, 1,
     &      id%STEP(1), id%BUFR(1), id%LBUFR, id%LBUFR_BYTES,
     &      CWORK(1), size(CWORK),
     &      LSCAL, PT_SCALING(1), size(PT_SCALING),
     &      C_DUMMY, 1 , 1, IDUMMY, 1
     &      )
        ELSE
          CALL ZMUMPS_GATHER_SOLUTION(id%NSLAVES,id%N,
     &      id%MYID, id%COMM, NBRHS_EFF,
     &      SOLVET, C_Y, id%N, NBRHS_EFF,
     &      id%PTLUST_S(1), id%KEEP(1),id%KEEP8(1),
     &      id%PROCNODE_STEPS(1),
     &      IS(1), LIW_PASSED,
     &      id%STEP(1), id%BUFR(1), id%LBUFR, id%LBUFR_BYTES,
     &      CWORK(1), size(CWORK),
     &      LSCAL, PT_SCALING(1), size(PT_SCALING),
     &      id%RHSCOMP(IBEG_RHSCOMP), LD_RHSCOMP, NBRHS_EFF,
     &      PTR_POSINRHSCOMP_BWD(1), id%N)
        ENDIF
        DEALLOCATE( CWORK )
        END SUBROUTINE ZMUMPS_PP_SOLVE
      END SUBROUTINE ZMUMPS_SOLVE_DRIVER