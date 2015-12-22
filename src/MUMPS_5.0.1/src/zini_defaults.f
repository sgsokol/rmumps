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
C**********************************************************************
C
      SUBROUTINE ZMUMPS_SET_TYPE_SIZES( K34, K35, K16, K10 )
      IMPLICIT NONE
C
C     Purpose:
C     =======
C
C     Set the size in bytes of an "INTEGER" in K34
C     Set the size of the default arithmetic (DOUBLE PRECISION, DOUBLE PRECISION,
C     COMPLEX(kind=8) or DOUBLE COMPLEX(kind=8)) in K35
C     Set the size of floating-point types that are real or double
C     precision even for complex versions of MUMPS (DOUBLE PRECISION for S and
C     C versions, DOUBLE PRECISION for D and Z versions)
C     Assuming that the size of an INTEGER(8) is 8, store the ratio
C     nb_bytes(INTEGER(8)) / nb_bytes(INTEGER) = 8 / K34 into K10.
C
C     In practice, we have:
C
C     K35:  Arithmetic   Value    Value for T3E
C              S           4           8
C              D           8          16
C              C           8          16
C              Z          16          32
C
C     K16 = K35 for S and D arithmetics
C     K16 = K35 / 2 for C and Z arithmetics
C
C     K34= 4 and K10 = 2, except on CRAY machines or when compilation
C     flag -i8 is used, in which case, K34 = 8 and K10 = 2
C
C
      INTEGER, INTENT(OUT) :: K34, K35, K10, K16
      INTEGER SIZE_INT, SIZE_REAL_OR_DOUBLE ! Type must match MUMPS_INT 
      INTEGER I(2)
      DOUBLE PRECISION R(2) ! Will be DOUBLE PRECISION if 1
      CALL MUMPS_SIZE_C(I(1),I(2),SIZE_INT)
      CALL MUMPS_SIZE_C(R(1),R(2),SIZE_REAL_OR_DOUBLE)
      K34 = int(SIZE_INT)
      K10 = 8 / K34
      K16 = int(SIZE_REAL_OR_DOUBLE)
      K35 = K16
      K35 = K35 * 2
      RETURN
      END SUBROUTINE ZMUMPS_SET_TYPE_SIZES
C
C**********************************************************************
C
      SUBROUTINE ZMUMPSID( NSLAVES, LWK_USER, CNTL, ICNTL,
     &                    KEEP,KEEP8,
     &                    INFO, INFOG, RINFO, RINFOG, SYM, PAR,
     &                    DKEEP)
      IMPLICIT NONE
C
C  Purpose
C  =======
C
C  The elements of the arrays CNTL and ICNTL control the action of
C  ZMUMPS, ZMUMPS_ANA_DRIVER, ZMUMPS_FAC_DRIVER, ZMUMPS_SOLVE_DRIVER 
C  Default values for the elements are set in this routine.
C
      DOUBLE PRECISION    DKEEP(130)
      DOUBLE PRECISION    CNTL(15), RINFO(40), RINFOG(40)
      INTEGER ICNTL(40), KEEP(500), SYM, PAR, NSLAVES
      INTEGER INFO(40), INFOG(40)
      INTEGER(8) KEEP8(150)
      INTEGER LWK_USER
C
C  Parameters
C  ==========
C===========================================
C       Arrays for control and information
C===========================================
C
C  N  Matrix order
C
C  NELT Number of elements for matrix in ELt format
C
C
C  SYM = 0 ... initializes the defaults for unsymmetric code
C      = 1,2 ... initializes the defaults for symmetric code
C
C
C
C  PAR = 0 ... instance where host is not working
C      = 1 ... instance where host is working as a normal node.
C              (host uses more memory than other processors in
C               the latter case)
C
C  CNTL and the elements of the array ICNTL control the action of
C     ZMUMPS Default values
C     are set by ZMUMPSID.  The elements of the arrays RINFO 
C     and INFO provide information on the action of ZMUMPS.
C  
C  CNTL(1) has default value 0.01 and is used for
C     threshold pivoting. Values greater than 1.0
C     are treated as 1.0, and less than zero as zero.
C     In general, a larger value of CNTL(1) leads to 
C     greater fill-in but a more accurate factorization.
C     If CNTL(1) is nonzero, numerical pivoting will be performed. 
C     If CNTL(1) is zero, no pivoting will be performed and 
C     the subroutine will fail if a zero pivot is encountered. 
C     If the matrix A is diagonally dominant, then
C     setting CNTL(1) to zero will decrease the factorization 
C     time while still providing a stable decomposition.
C
C  CNTL(2) must be set to the tolerance for convegence of iterative
C     refinement. 
C     Default value is sqrt(macheps).
C     Values less than zero are treated as sqrt(macheps).
C
C  CNTL(3) is only used combined with null pivot row
C     detection (ICNTL(24) .eq. 1) and to Rank-Revealing (RR) option.
C     It must be set to the absolute threshold for numerical pivoting.
C     Default value is 0.0. 
C     Let A_{preproc} be the preprocessed matrix to be factored (see
C       equation in the user's guide).
C     A pivot is considered to be null if the infinite norm of its row/column
C     is smaller than a threshold. Let MACHEPS be the machine precision and 
C     ||.|| be the infinite norm.
C    The computed threshold value for postponing pivots in case of RR on root
C     is stored in "SEUIL" and then "SEUIL_LDLT_NIV2"
C     which are identical in current version.
C     This absolute threshold value is stored in DKEEP(9).
C
C     The absolute value to detect a null pivot (when ICNTL(24) .NE.0) 
C     is stored in DKEEP(1) and must be smaller than
C     SEUIL when combined with RR on root.
C
C     IF (ICNTL(16).NE.0) THEN 
C      RR on root is active 
C      IF  (CNTL3 .LT. ZERO) THEN
C          SEUIL = abs(CNTL(3))
C      ELSE IF  (CNTL3 .GT. ZERO) THEN
C          SEUIL = CNTL3*ANORMINF
C      ELSE  !  (CNTL(3) .EQ. ZERO) THEN
C          SEUIL = N*EPS*ANORMINF  ! standard articles
C      ENDIF
C      IF (ICNTL(24).NE.0) THEN
C       null pivot detection
C       IF (CNTL(6).GT.0.AND.CNTL(6).LT.1) THEN
C         we want DKEEP(1) < SEUIL
C         DKEEP(1) = SEUIL*CNTL(6) ! ideally it could be SEUIL*CNTL(6)
C       ELSE 
C         DKEEP(1) = SEUIL* 0.01D0 
C       ENDIF
C      ENDIF
C
C     ELSE (ONLY NULL PIVOT detection is active) 
C         we keep stratgy used in MUMPS_4.10
C      IF CNTL(3)  > 0 THEN
C          DKEEP(1) = CNTL(3)  ||A_{preproc}|| 
C      ELSE IF CNTL(3)  = 0.0 THEN 
C          DKEEP(1) = MACHEPS 10^{-5} ||A_{preproc}||
C      ELSE IF CNTL(3)  < 0 THEN 
C          DKEEP(1) = abs(CNTL(3))   ! this was added for EDF 
C                                    ! in the context of SOLSTICE project
C      ENDIF
C
C  CNTL(4) must be set to value for static pivoting.
C     Default value is -1.0
C     Note that static pivoting is enabled only when 
C     Rank-Revealing and null pivot detection
C     are off (KEEP(19).EQ.0).AND.(KEEP(110).EQ.0).
C     If negative, static pivoting will be set OFF (KEEP(97)=0)
C     If positive, static pivoting is ON (KEEP(97=1) with threshold CNTL(4)
C     If = 0, static pivoting is ON with threshold MACHEPS^1/2 || A ||
C
C  CNTL(5) fixation for null pivots
C     Default value is 0.0
C     Only active if ICNTL(24) = 1
C     If > 0 after finding a null pivot, it is set to CNTL(5) x ||A||
C     (This value is stored in DKEEP(2))
C     If <=  0 then the row/column (except the pivot) is set to zero
C     and the pivot is set to 1
C     Default is 0.
C     Note that in the symmetric parallel case, some elements of the column
C     are not available on the local processor and cannot be set to 0 easily.
C     In such cases, in the current version, 
C            -the corresponding pivot is first set
C                to a large value instead of 1, even when CNTL(5) < 0.
C            -Updating of the off diag block is done with this large 
C                value 
C            -diagonal value is then reset to zero
C
C  CNTL(6) expresses the ratio between 
C        absolute criterion for null pivots and absolute criterion 
C        for posponing pivots before partial pivoting analysis of pivots.
C        Typically 
C          let SEUIL = F(CNTL(3)), and  0 < CNTL(6) < 1
C          SEUIL is stored in DKEEP(9)
C          if ||Pivot row|| < SEUIL*CNTL(6) then
C              null pivot row detected (correct only if LDLT 
C                                 for LU pivot_col must be checked too)
C          else if || Pivot_Row || < SEUIL then
C              pospone pivot
C          else
C              partial threshold pivoting
C          endif
C
C  -----------------------------------------
C
C  ICNTL(1) has default value 6.
C     It is the output stream for error messages.  
C     If it is set to zero, these
C     messages will be suppressed.
C  
C  ICNTL(2) has default value 0.
C     It is the output stream for diagnostic printing and
C     for warning messages that are local to each MPI process.
C     If it is set to zero, these messages are suppressed.
C  
C  ICNTL(3) -- Host only
C            It is the output stream for diagnostic printing
C            and  for warning messages. Default value is 6.
C            If it is set to zero, these messages are suppressed.
C
C  ICNTL(4) is used by ZMUMPS to control printing of error, 
C     warning, and diagnostic messages. It has default value 2.
C     Possible values are:
C    
C    <1       __No messages output.
C     1       __Only error messages printed.
C     2       __Errors and warnings printed.
C     3       __Errors and warnings and terse diagnostics 
C                (only first ten entries
C               of arrays printed).
C     4       __Errors and warnings and all information 
C               on input and output parameters printed.
C
C
C  ICNTL(5) is the format of the input matrix and rhs
C     0: assembled matrix, assembled rhs
C     1: elemental matrix, assembled rhs
C     Default value is 0.
C  
C  ICNTL(6) has default value 7 for unsymmetric and
C      general symmetric matrices, and 0 for SPD matrices.
C      It is only accessed and operational
C      on a call that includes an analysis phase 
C      (JOB = 1, 4, or 6).
C      In these cases, if ICNTL(6)=1, 2, 3, 4, 5, 6 or 7,
C      a column permutation based on algorithms described in 
C      Duff and Koster, 1997, *SIMAX <20>, 4, 889-901, 
C      is applied to the original matrix. Column permutations are
C      then applied to the original matrix to get a zero-free diagonal.
C      Except for ICNTL(6)=1, the numerical values of the 
C      original matrix, id%A(NE), need be provided by the user 
C      during the analysis phase.
C      If ICNTL(6)=7, based on the structural symmetry of the 
C      input matrix the value of ICNTL(6) is automatically chosen.
C     If the ordering is provided by the user 
C     (ICNTL(7)=1) then the value of ICNTL(6) is ignored.
C  
C  ICNTL(7) has default value 7 and must be set by the user to
C     1 if the pivot order in IS is to be used. 
C     Effective value of ordering stored in KEEP(256).
C     Possible values are (depending on the softwares installed)
C       0 AMD: Approximate minimum degree (included in ZMUMPS package)
C       1 Ordering provided by the user
C       2 Approximate minimum fill (included in ZMUMPS package)
C       3 SCOTCH (see http://gforge.inria.fr/projects/scotch/)
C         should be downloaded/installed separately.
C       4 PORD from Juergen Schulze (js@juergenschulze.de) 
C         PORD package is extracted from the SPACE-1.0 package developed at the
C         University of Paderborn by Juergen Schulze 
C         and is provided as a separate package.
C       5 Metis ordering should be downloaded/installed separately.
C       6 Approximate minimum degree with automatic quasi 
C           dense row detection (included in ZMUMPS package). 
C           (to be used when ordering time with AMD is abnormally large)
C       7 Automatic choice done during analysis phase
C     For any other
C     value of ICNTL(7), a suitable pivot order will be 
C     chosen automatically.
C
C  ICNTL(8)  is used to describe the scaling strategy.
C     Default value is 77.
C     Note that scaling is performed only when the numerical
C     factorization step is performed (JOB = 2, 4>, 5>, or 6>).
C     If ICNTL(8) is not equal to
C     any of the values listed below then ICNTL(8) is treated 
C     as if it had its default value of 0 (no scaling).  
C     If the matrix is known to be very badly scaled, 
C     our experience has been that option 6 is the most robust but
C     the best scaling is very problem dependent.
C     If ICNTL(8)=0, COLSCA and ROWSCA are dummy arguments
C     of the subroutine that are not accessed.
C     Possible values of ICNTL(8) are:
C 
C     -2 scaling computed during analysis (and applied during the
C       factorization)
C
C     -1 the user must provide the scaling in arrays 
C        COLSCA and ROWSCA
C  
C     0 no scaling
C  
C     1 Diagonal scaling
C  
C     2 not defined
C
C     3 Column scaling
C
C     4 Row and column scaling
C
C     5,6 not defined
C     7, 8 Scaling based on Daniel Ruiz and Bora Ucar's work done 
C     during the ANR-SOLSTICE project.
C     Reference for this work are:
C     The scaling algorithms are based on those discussed in
C     [1] D. Ruiz, "A scaling algorithm to equilibrate both rows and 
C         columns norms in matrices", Tech. Rep. Rutherford 
C         Appleton Laboratory, Oxon, UK and ENSEEIHT-IRIT, 
C         Toulouse, France, RAL-TR-2001-034 and RT/APO/01/4, 2001.
C     [2] D. Ruiz and B. Ucar, "A symmetry preserving algorithm for
C         matrix scaling", in preparation as of Jan'08.
C     This scaling can work on both centralized and distributed
C     assembled input matrix format. (it works for both symmetric
C     and unsymmetric matrices)
C     Option 8 is similar to 7 but more rigourous and expensive to compute.
C     77 Automatic choice of scaling value done. Proposed algo:
C          if (sym=1) then
C           default = 0
C          else
C           if distributed matrix entry then 
C              default = 7
C           else
C               if (mc64 called or mc77 based matching) then
C                 default=-2 and ordering is computed during analysis
C               else
C                 default = 7
C               endif
C           endif
C          endif
C  
C  ICNTL(9) has default value 1. If ICNTL(9)=1
C     the system of equations A * x = b  is solved. For other
C     values the system A^T *  x = b  is solved.
C     When ICNTL(30) (compute selected entries in A-1) is activated 
C     ICNTL(9) is ignored.
C  
C  ICNTL(10) has default value 0.
C     If ICNTL(10)=0 : iterative refinement is not performed.
C     Values of  ICNTL(10) < 0 : a fix number of steps equal 
C     to ICNTL(10) of IR is done.
C     Values of ICNTL(10) > 0 : mean a maximum of ICNTL(10) number 
C                           of steps of IR is done, and a test of 
C                           convergence is used
C  
C  ICNTL(11) has default value 0.
C     A value equal to 1 will return a backward error estimate in
C     RINFO(4-11).
C     A value equal to 2 will return a backward error estimate in
C     RINFO(4-8). No LCOND 1, 2 and forward error are computed.
C     If ICNTL(11) is negative, zero or greater than 2 no estimate 
C     is returned.
C
C
C  ICNTL(12) has default value 0 and define the strategy for 
C     LDLT orderings
C     0 : automatic choice
C     1 : usual ordering (nothing done)
C     2 : ordering on the compressed graph, available with all orderings
C         except with AMD
C     3 : constraint ordering, only available with AMF, 
C         -> reset to 2 with other orderings
C     Other values are treated as 1 (nothing done).
C     On output KEEP(95) holds the internal value used and INFOG(24) gives 
C     access to KEEP(95) to the user.
C     in LU facto it is always reset to 1
C
C     - ICNTL(12) = 3 has a lower priority than ICNTL(7)
C     thus if ICNTL(12) = 3 and the ordering required is not AMF
C     then ICNTL(12) is set to 2
C
C     - ICNTL(12) = 2 has a higher priority than ICNTL(7)
C     thus if ICNTL(12) = 2 and the ordering required is AMD
C     then the ordering used is QAMD
C
C     - ICNTL(12) has a higher priority than ICNTL(6) and ICNTL(8)
C     thus if ICNTL(12) = 2 then ICNTL(6) is automatically
C     set to a value between 1-6
C     if ICNTL(12) = 3 then ICNTL(6) is automatically
C     set to 5 and ICNTL(6) is set to -2 (we need the scaling  factors
C     to define free and constrained variables)
C
C  ICNTL(13) has default value 0 and allows for selecting Type 3 node.
C            IF ICNTL(13).GT. 0 scalapack is forbidden. Otherwise,
C            scalapack will be activated if the root is large enough.
C            Furthermore
C             IF ((ICNTL(13).GT.0) .AND. (NSLAVES.GT.ICNTL(13),
C               or ICNTL(13)=-1 THEN
C               extra splitting of the root will be activated
C               and is controlled by abs(KEEP(82)).
C               The order of the root node is divided by KEEP(82)
C             ENDIF
C            If ICNTL(13) .EQ. -1 then splitting of the root
C            is done whatever the nb of procs is.
C
C            Authorizing extra root spliting 
C            during analysis might be interesting
C            to further split the root node
C            (combined for example with 
C             null pivot detection option ICNTL(24)=1 OR ICNTL(16))
C             
C            To summarize:
C               -1         : root splitting and scalapack on
C               0  or < -1 : root splitting off and sclalapack on
C               > 0        : scalapack off
C
C  ICNTL(14) has default value 20 (or 30, or 5 depending on NSLAVES,
C            SYM,...) and is the value for memory relaxation 
C            so called "PERLU" in the following.
C
C  ICNTL(18) has default value 0 and is only accessed by the host during
C      the analysis phase if the matrix is assembled (ICNTL(5))= 0).
C      ICNTL(18) defines the strategy for the distributed input matrix.
C      Possible values are: 
C      0: input matrix is centralized on the host. This is the default
C      1: user provides the structure of the matrix on the host at analysis, 
C        ZMUMPS returns
C        a mapping and user should provide the matrix distributed according 
C        to the mapping
C      2:  user provides the structure of the matrix on the host at analysis,  
C          and the
C          distributed matrix on all slave processors at factorization. 
C          Any distribution is allowed
C      3: user directly provides the distributed matrix input both 
C         for analysis and factorization
C     
C      For flexibility and performance issues, option 3 is recommended.
C
C   ICNTL(19) has default value 0 and is only accessed by the host
C    during the analysis phase. If ICNTL(19) \neq 0 then Schur matrix will
C    be returned to the user.
C    The user must set on entry on the host node (before analysis):
C    the integer variable SIZE\_SCHUR to the size fo the Schur matrix,
C    the integer array pointer LISTVAR\_SCHUR to the list of indices
C    of the schur matrix.
C     if = 0      : Schur is off and the root node gets factorized 
C     if = 1      : Schur is on and the Schur complement is returned entirely
C                   on a memory area provided by the user ONLY on the host node
C     if = 2 or 3 : Schur is on and the Schur complement is returned in a 
C                   distributed fashion according to a 2D block-cyclic
C                   distribution. In the case where the matrix is symmetric
C                   the lower part is returned if =2 or the complete 
C                   matrix if =3.
C
C   ICNTL(20) has default value 0 and is only accessed  by the host
C    during the solve phase. If ICNTL(20)=0, the right-hand side must given
C    in dense form in the structure component RHS.
C    If ICNTL(20)=1,2,3, then the right-hand side must be given in sparse form
C    using the structure components IRHS\_SPARSE, RHS\_SPARSE, IRHS\_PTR and
C    NZ\_RHS. 
C    When the right-hand side is provided in sparse form then duplicate entries
C    are summed.
C
C       0 : dense RHS
C       1,2,3 : Sparse RHS
C          1 The decision of exploiting sparsity of the right-hand side to 
C            accelerate the solution phase is done automatically.
C          2 Sparsity of the right-hand sides is NOT exploited 
C            to improve solution phase.
C          3 Sparsity of the right-hand sides is exploited 
C            to improve solution phase.
C    Values different from 0,1, 2,3 are treated as 0.
C    For sparse RHS recommended value is 1.
C
C   ICNTL(21) has default value 0 and is only accessed by the host
C    during the solve phase. If ICNTL(21)=0, the solution vector will be assembled
C    and stored in the structure component RHS, that must have been allocated by
C    the user. If ICNTL(21)=1, the solution vector is kept distributed at the
C    end of the solve phase, and will be available on each slave processor
C    in the structure components ISOL_loc and SOL_loc. ISOL_loc and SOL_loc
C    must then have been allocated by the user and must be of size at least
C    INFO(23), where INFO(23) has been returned by ZMUMPS at the end of the
C    factorization phase.
C    Values of ICNTL(21) different from 0 and 1 are currently treated as 0.
C
C   ICNTL(22) (saved in KEEP(201) controls the OOC setting (0=incore, >0 =OOC,
C    and later -1=discard factors completely)
C    It has default value 0 (incore). 
C    If set before analysis then special setting and massage of the tree 
C    might be done (so far only extra splitting  CUTNODES) is performed.
C    It is then accessed by the host
C    during the factorization phase. If ICNTL(22)=0, then no attempt
C    to use the disks is made. If ICNTL(22)=1, then ZMUMPS will store
C    the computed factors on disk for later use during the solution
C    phase. 
C
C   ICNTL(23) has default value 0 and is accessed by ALL processors
C    at the beginning of the factorization phase.  If positive 
C    it corresponds to the maximum size of the working memory 
C    in MegaBytes that MUMPS can allocate per working processor.
C    If only the host
C    value is non zero, then other processors also use the value on
C    the host. Otherwise, each processor uses the local value
C    provided.
C
C   ICNTL(24) default value is 0
C     if = 0 no null pivot detection (CNTL(5) and CNTL(3) are inactive),
C        = 1 null pivot row detection; CNTL(3) and CNTL(5) are 
C            then used to describe the action taken.
C
C
C   ICNTL(25) has default value 0 and is only accessed by the
C   host during the solution stage. It is only significant if
C   a null space basis was requested during the factorization
C   phase (INFOG(28) .GT. 0); otherwise a normal solution step
C   is performed.
C   If ICNTL(25)=0, then a normal solution step is performed,
C   on the internal problem (excluding the null space). 
C   No special property on the solution (discussion with Serge) 
C   If ICNTL(25)=i, 1 <= i <= INFOG(28), then the i-th vector
C   of the null space basis is computed. In that case, note
C   that NRHS should be set to 1.
C   If ICNTL(25)=-1, then all null space is computed. The
C   user should set NRHS=INFOG(28) in that case.
C   Note that centralized or distributed solutions are
C   applicable in that case, but that iterative refinement,
C   error analysis, etc... are excluded. Note also that the
C   option to solve the transpose system (ICNTL(9)) is ignored.
C
C
C   ICNTL(26) has default value 0 and is accessed on the host only
C     at the beginning of the solution step.
C     It is only effective if the Schur option is ON.
C     (copy in KEEP(221))
C
C
C     During the solution step, a value of 0 will perform a normal
C     solution step on the reduced problem not involving the Schur
C     variables.
C     During the solution step, if ICNTL(26)=1 or 2, then REDRHS
C     should be allocated of size at least LREDRHS*(NRHS-1)+
C     SIZE_SCHUR, where LREDRHS is the leading dimension of
C     LREDRHS (LREDRHS >= SIZE_SCHUR).
C
C     If ICNTL(26)=1, then only a forward substitution is performed,
C     and a reduced RHS will be computed and made available in
C     REDRHS(i+(k-1)*LREDRHS), i=1, ..., SIZE_SCHUR, k=1, ..., NRHS.
C     If ICNTL(26)=2, then REDRHS(i+(k-1)*LREDRHS),i=1, SIZE_SCHUR, k=1,NRHS is
C     considered to be the solution corresponding to the Schur
C     variables. It is injected in ZMUMPS, that computes the solution
C     on the "internal" problem during the backward substitution.
C
C ICNTL(27)  controls the blocking factor for multiple right-hand-sides
C during the solution phase. 
C It influences both the memory used (see INFOG(30-31)) and 
C the solution time
C (Larger values of ICNTL(27) leads to larger memory requirements). 
C Its tuning can be critical when
C the factors are written on disk (out-of core, ICNTL(22)=1). 
C A negative value indicates that  automatic setting is performed by the solver. 
C Default value is -24.  
C
C
C   ICNTL(28) decides whether parallel or sequential analysis should be used. Three
C       values are possible at the moment:
C       0: automatic. This defaults to sequential analysis
C       1: sequential. In this case the ordering strategy is defined by ICNTL(7)
C       2: parallel. In this case the ordering strategy is defined by ICNTL(29)
C
C   ICNTL(29) defines the ordering too to be used during the parallel analysis. Three
C       values are possible at the moment:
C       0: automatic. This defaults to PT-SCOTCH
C       1: PT-SCOTCH.
C       2: ParMetis.
C
C
C   ICNTL(30) controls the activation of functionality A-1. 
C             It has default value 0 and is only accessed by the master 
C             during the solution phase. It enables the solver to 
C             compute entries in the inverse of the original matrix.
C             Possible values are:
C              0  normal solution
C              other values:   compute entries in A-1
C             When ICNTL(30).NE.0 then the user 
C             must describe on entry to the solution phase, 
C             in the sparse right-hand-side 
C             (NZ_RHS, NRHS, RHS_SPARSE, IRHS_SPARSE, IRHS_PTR) 
C             the target entries of A-1 that need be computed.
C             Note that RHS_SPARSE must be allocated but need not be
C             initialized.
C             On output RHS_SPARSE then holds the requested 
C             computed values of A-1.
C             Note that when ICNTL(30).NE.0 then 
C              - sparse right hand side interface is implicitly used 
C                functionality (ICNTL(20)= 1) but RHS need not be 
C                allocated since computed A-1 entries will be stored 
C                in place.
C              - ICNTL(9) option (solve Ax=b or Atx=b) is ignored
C             In case of duplicate entries in the sparse rhs then 
C             on output duplicate entries in the solution are provided 
C             in the same place.
C             This need not be mentioned in the spec since it is a 
C             "natural" extension.
C
C   -----------
C   Fwd in facto
C   -----------
C   ICNTL(31) Must be set before analysis to control storage 
C             of LU factors. Default value is 0. Out of range 
C             values considered as 0. 
C             (copied in KEEP(251) and broadcast, 
C              when setting of ICNTL(31) 
C              results in not factors to be stored then 
C              KEEP(201) = -1, OOC is "suppressed")
C            0 Keep factors needed for solution phase 
C              (when option forward during facto is used then
C               on unsymmetric matrices L factors are not stored)
C            1 Solve not needed (solve phase will never be called).
C              When the user is only interested in the inertia or the 
C              determinant then
C              all factor matrices need not be stored.
C              This can also be useful for testing : 
C              to experiment facto OOC without
C              effective storage of factors on disk.
C            2 L factors not stored: meaningful when 
C              - matrix is unsymmetric and fwd performed during facto
C              - the user is only interested in the null-space basis
C              and thus only need the U factors to be stored.
C              
C   -----------
C   Fwd in facto
C   -----------
C   ICNTL(32) Must be set before analysis to indicate whether 
C             forward is performed during factorization. 
C             Default value is 0 (normal factorization without fwd)
C             (copied in KEEP(252) and broadcast)
C            0 Normal factorization (default value)
C            1 Forward performed during factorization
C
C
C   ICNTL(33) Must be set before the factorization phase to compute
C             the determinant. See also KEEP(258), KEEP(259),
C             DKEEP(6), DKEEP(7), INFOG(34), RINFOG(12), INFOG(34)
C
C            If ICNTL(33)=0 the determinant is not computed
C            For all other values, the determinant is computed. Note that
C            null pivots and static pivots are excluded from the
C            computation of the determinant.
C
C   ICNTL(34-38) are not used in this version
C
C
C
C=========================
C  ARRAYS FOR INFORMATION
C========================
C 
C-----
C INFO is an INTEGER array of length 40 that need not be 
C     set by the user.
C-----
C 
C  INFO(1) is zero if the routine is successful, is negative if an
C     error occurred, and is positive for a warning (see ZMUMPS for
C     a partial documentation and the userguide for a full documentation
C     of INFO(1)).
C  
C  INFO(2)   holds additional information concerning the
C     error (see ZMUMPS).
C  
C ------------------------------------------
C Statistics produced after analysis phase
C ------------------------------------------
C  
C  INFO(3)   Estimated real space needed for factors.
C  
C  INFO(4)   Estimated integer space needed for factors.
C  
C  INFO(5)  Estimated maximum frontal size.
C  
C  INFO(6)  Number of nodes in the tree.
C  
C  INFO(7)  Minimum value of integer working array IS (old MAXIS) 
C     estimated by the analysis phase
C     to run the numerical factorization.
C  
C  INFO(8)  Minimum value of real/complex arry S (old MAXS) 
C     estimated by the analysis phase
C     to run the numerical factorization.
C  
C  INFO(15) Estimated size in MBytes of all ZMUMPS internal data
C      structures to run factorization
C
C  INFO(17) provides an estimation (minimum in Megabytes) 
C  of the total memory required to run 
C  the numerical phases  out-of-core.
C  This memory estimation corresponds to 
C  the least memory consuming out-of-core strategy and it can be
C  used as a lower bound if the user wishes to provide ICNTL(23).
C ---------------------------------------
C Statistics produced after factorization
C ---------------------------------------
C  INFO(9)  Size of the real space used to store the LU factors.
C  
C  INFO(10)  Size of the integer space used to store the LU factors.
C  
C  INFO(11)  Order of largest frontal matrix.
C  
C  INFO(12)  Number of off-diagonal pivots.
C  
C  INFO(13)  Number of uneliminated variables sent to the father.
C  
C  INFO(14)  Number of memory compresses.
C
C  INFO(18)  On exit to factorization: 
C              Local number of null pivots (ICNTL(24)=1) 
C              on the local processor even on master. 
C              (local size of array PIVNUL_LIST).
C              Note that it does not include null pivots 
C              that might have been 
C              further detected on the root (ICNTL(16).NE.0).
C
C  INFO(19) - after analysis: 
C           Estimated size of the main internal integer workarray IS
C     (old MAXIS) to run the numerical factorization out-of-core.
C
C  INFO(21) - after factorization: Effective space used in the main
C           real/complex workarray S -- or in the workarray WK_USER, 
C           in the case where WK_USER is provided.
C
C  INFO(22) - after factorization:
C      Size in millions of bytes of memory effectively used during
C      factorization.
C      This includes the memory effectively used in the workarray
C      WK_USER, in the case where WK_user is provided.
C
C  INFO(23) - after factorization: total number of pivots eliminated
C      on the processor. In the case of a distributed solution (see
C      ICNTL(21)), this should be used by the user to allocate solution
C      vectors ISOL_loc and SOL_loc of appropriate dimensions
C      (ISOL_LOC of size INFO(23), SOL_LOC of size LSOL_LOC * NRHS
C      where LSOL_LOC >= INFO(23)) on that processor, between the
C      factorization and solve steps.
C 
C  INFO(24) - after analysis: estimated number of entries in factors on
C      the processor. If negative, then
C      the absolute value corresponds to {\it millions} of entries
C      in the factors.
C      Note that in the unsymmetric case, INFO(24)=INFO(3).
C      In the symmetric case, however, INFO(24) < INFO(3).
C  INFO(25) - after factorization: number of tiny pivots (number of
C       pivots modified by static pivoting) detected on the processor.
C  INFO(26) - after solution: 
C                 effective size in Megabytes of all working space
C                 to run  the solution phase.
C     (The maximum and sum over all processors are returned 
C    respectively in INFOG(30) and INFOG(31)).
C  INFO(27) - after factorization: effective number of entries in factors
C      on the processor. If negative, then
C      the absolute value corresponds to {\it millions} of entries
C      in the factors.
C      Note that in the unsymmetric case, INFO(27)=INFO(9).
C      In the symmetric case, however, INFO(27) < INFO(9).
C      The total number of entries over all processors is
C      available in INFOG(29).
C
C -------------------------------------------------------------
C -------------------------------------------------------------
C RINFO is a DOUBLE PRECISION/DOUBLE PRECISION array of length 40 that 
C       need not be set by the user. This array supplies 
C       local information on the execution of ZMUMPS.
C
C
C RINFOG is a DOUBLE PRECISION/DOUBLE PRECISION array of length 40 that
C       need not be set by the user. This array supplies
C       global information on the execution of ZMUMPS.
C       RINFOG is only significant on processor 0
C
C
C  RINFO(1) hold the estimated number of floating-point operations
C           for the elimination process on the local processor
C
C  RINFOG(1) hold the estimated number of floating-point operations
C           for the elimination process on all processors
C
C  RINFO(2)  Number of floating-point operations
C     for the assembly process on local processor.
C
C  RINFOG(2) Number of floating-point operations
C     for the assembly process.
C
C  RINFO(3)  Number of floating-point operations
C     for the elimination process on the local processor.
C
C  RINFOG(3)  Number of floating-point operations
C     for the elimination process on all processors.
C
C----------------------------------------------------
C Statistics produced after solve with error analysis
C----------------------------------------------------
C
C  RINFOG(4) Infinite norm of the input matrix.
C
C  RINFOG(5) Infinite norm of the computed solution, where
C
C  RINFOG(6) Norm of scaled residuals
C
C  RINFOG(7), `RINFOG(8) and `RINFOG(9) are used to hold information 
C     on the backward error.
C     We calculate an estimate of the sparse backward error using the
C     theory and measure developed 
C     by Arioli, Demmel, and Duff (1989). The scaled residual w1 
C     is calculated for all equations except those 
C     for which numerator is nonzero and the denominator is small. 
C     For the exceptional equations, w2, is used instead. 
C     The largest scaled residual (w1) is returned in 
C     RINFOG(7) and the largest scaled
C     residual (w2) is returned in `RINFOG(8)>. If all equations are 
C     non exceptional then zero is returned in `RINFOG(8).
C     The upper bound error is returned in `RINFOG(9).
C
C  RINFOG(14) Infinite norm of the residual
C
C  RINFOG(15) 2-norm of the residual
C
C===========================
C DESCRIPTION OF KEEP8 ARRAY
C===========================
C
C   KEEP8 is a 64-bit integer array of length 150 that need not
C   be set by the user
C
C===========================
C DESCRIPTION OF KEEP ARRAY
C===========================
C
C  KEEP is an INTEGER array of length 500 that need not 
C     be set by the user.
C
C
C=============================
C Description of DKEEP array
C=============================
C
C DKEEP internal control array for DOUBLE PRECISION parameters
C     of size 30
C===================================
C Default values for control arrays
C==================================
C     uninitialized values should be 0
      LWK_USER = 0
      KEEP(1:500) = 0
      KEEP8(1:150)= 0_8
      INFO(1:40)  = 0
      INFOG(1:40) = 0
      ICNTL(1:40) = 0
      RINFO(1:40) = 0.0D0
      RINFOG(1:40)= 0.0D0
      CNTL(1:15)  = 0.0D0
      DKEEP(1:130) = 0.0D0
C     ----------------
C     Symmetric code ?
C     ----------------
      KEEP( 50 ) = SYM
C     Check value of SYM
      IF (SYM.EQ.1) THEN
C
C     this option is not available with the complex
C     code on symmetric matrices.
C     We set KEEP(50) to 2 and will exploit symmetry
C     up to the root.
          KEEP(50) = 2
      ENDIF
C     -------------------------------------
C     Only options 0, 1, or 2 are available
C     -------------------------------------
      IF ( KEEP(50).NE.1 .and. KEEP(50).NE.2 ) KEEP( 50 ) = 0
C     threshold value for pivoting
      IF ( KEEP(50) .NE. 1 ) THEN
        CNTL(1)   = 0.01D0
      ELSE
        CNTL(1)   = 0.0D0
      END IF
      CNTL(2) = sqrt(epsilon(0.0D0))
      CNTL(3) = 0.0D0
      CNTL(4) = -1.0D0
      CNTL(5) = 0.0D0
      CNTL(6) = -1.0D0
C     Working host ?
      KEEP(46) = PAR
      IF ( KEEP(46) .NE. 0 .AND.
     &     KEEP(46) .NE. 1 ) THEN
C          ----------------------
C          If out-of-range value,
C          use a working host
C          ----------------------
           KEEP(46) = 1
      END IF
C     control printing
      ICNTL(1)  = 6
      ICNTL(2)  = 0
      ICNTL(3)  = 6
      ICNTL(4)  = 2
C     format of input matrix
      ICNTL(5)  = 0
C     maximum transversal (0=NO, 7=automatic)
      IF (SYM.NE.1) THEN
       ICNTL(6)  = 7
      ELSE
       ICNTL(6)  = 0
      ENDIF
C     Ordering option (icntl(7))
C     Default is automatic choice done during analysis
      ICNTL(7) = 7
C     ask for scaling (0=NO, 4=Row and Column)
C     Default value is 77: automatic choice for analysis
      ICNTL(8)  = 77
C     solve Ax=b (1) or Atx=b (other values)
      ICNTL(9)  = 1
C     Naximum number of IR (0=NO)
      ICNTL(10)  = 0
C     Error analysis (0=NO)
      ICNTL(11)  = 0
C     Control ordering strategy
C     automatic choice
      IF(SYM .EQ. 2) THEN
         ICNTL(12)  = 0
      ELSE
         ICNTL(12)  = 1
      ENDIF
C     Control of the use of ScaLAPACK for root node
C     If null space options asked, ScaLAPACK always
C     used and ICNTL(13) is not significant
C     ICNTL(13) = 0 : Root parallelism on (if size large enough)
C     ICNTL(13) = 1 : Root parallelism off
      ICNTL(13) = 0
C     Default value for the memory relaxation
      IF (SYM.eq.1.AND.NSLAVES.EQ.1) THEN 
        ICNTL(14) = 5  ! it should work with 0
      ELSE IF (NSLAVES .GT. 4) THEN
        ICNTL(14) = 30
      ELSE
        ICNTL(14) = 20
      END IF
C     Minimum size of the null space 
      ICNTL(15) = 0
C     Do not look for rank/null space basis
      ICNTL(16) = 0
C     Max size of null space
      ICNTL(17) = 0
C     Distributed matrix entry
      ICNTL(18) = 0
C     Schur  (default is not active)
      ICNTL(19) = 0
C     dense RHS by default
      ICNTL(20) = 0
C     solution vector centralized on host
      ICNTL(21) = 0
C     out-of-core flag
      ICNTL(22) = 0
C     MEM_ALLOWED (0: not provided)
      ICNTL(23) = 0
C     null pivots
      ICNTL(24) = 0
C     blocking factor for multiple RHS during solution phase
      ICNTL(27) = -24
C     analysis strategy: 0=auto, 1=sequential, 2=parallel
      ICNTL(28) = 1
C     tool used for parallel ordering computation :
C     0 = auto, 1 = PT-SCOTCH, 2 = ParMETIS
      ICNTL(29) = 0
C     --------- Non documented ICNTL options
C     Old or new symbolic factorization
      ICNTL(39) = 1
      ICNTL(40)  = 0 
C===================================
C Default values for some components 
C of KEEP array
C===================================
      KEEP(12) = 0
      KEEP(11) = 2147483646
      KEEP(24) = 18
      KEEP(68) = 0
      KEEP(36) = 1
      KEEP(1) = 5
      KEEP(7)  = 150
      KEEP(8)  = 120
      KEEP(57) = 500
      KEEP(58) = 250
      IF ( SYM .eq. 0 ) THEN
        KEEP(4)  = 32
        KEEP(3)  = 96
        KEEP(5)  = 16
        KEEP(6)  = 32
        KEEP(9)  = 700
        KEEP(85) =  300
        KEEP(62) =  50
      ELSE
        KEEP(4)  = 24
        KEEP(3)  = 96
        KEEP(5)  = 16
        KEEP(6)  = 32
        KEEP(9)  = 400
        KEEP(85) = 100
        KEEP(62) = 50
      END IF
      KEEP(63) = 60
      KEEP(48) = 5
      KEEP(17) = 0
      CALL ZMUMPS_SET_TYPE_SIZES( KEEP(34), KEEP(35),
     &                            KEEP(16), KEEP(10) )
#if defined(SP_)
      KEEP( 51 )  = 70
#else
      KEEP( 51 )  = 48
#endif
      KEEP(37) = max(800, int(sqrt(dble(NSLAVES+1))*dble(KEEP(51))))
      IF ( NSLAVES > 256 ) THEN
        KEEP(39) = 10000
      ELSEIF ( NSLAVES > 128 ) THEN
        KEEP(39) = 20000
      ELSEIF ( NSLAVES > 64 ) THEN
        KEEP(39) = 40000
      ELSEIF ( NSLAVES > 16 ) THEN
        KEEP(39) = 80000
      ELSE
        KEEP(39) = 160000
      END IF
      KEEP(40) = -1 - 456789
      KEEP(45) = 0
      KEEP(47) = 2
      KEEP(64) = 20
      KEEP(69) = 4
C     Forcing proportional is ok with strategy 5
      KEEP(75) = 1
      KEEP(76) = 2
      KEEP(77) = 30
      KEEP(79) = 0  ! old splitting
      !write(6,*) ' TEMPORARY new splitting active, K79=', KEEP(79)
      IF (NSLAVES.GT.4) THEN
          KEEP(78)=max(
     &       int(log(dble(NSLAVES))/log(dble(2))) - 2 
     &       , 0         )
      ENDIF
      KEEP(210) = 2 
      KEEP8(79) = -10_8
      KEEP(80) = 1
      KEEP(81) = 0
      KEEP(82) = 10
      KEEP(83) = min(8,NSLAVES/4)
      KEEP(83) = max(min(4,NSLAVES),max(KEEP(83),1))
      KEEP(86)=1
      KEEP(87)=0
      KEEP(88)=0
      KEEP(90)=1
      KEEP(91)=min(8, NSLAVES)
      KEEP(91) = max(min(4,NSLAVES),min(KEEP(83),KEEP(91)))
      IF(NSLAVES.LT.48)THEN
         KEEP(102)=150
      ELSEIF(NSLAVES.LT.128)THEN
         KEEP(102)=150
      ELSEIF(NSLAVES.LT.256)THEN
         KEEP(102)=200
      ELSEIF(NSLAVES.LT.512)THEN
         KEEP(102)=300
      ELSEIF(NSLAVES.GE.512)THEN
         KEEP(102)=400
      ENDIF
#if defined(OLD_OOC_NOPANEL)
      KEEP(99)=0  ! no panel -> synchronous / no buffer
#else
      KEEP(99)=4  ! new OOC -> asynchronous + buffer
#endif
      KEEP(100)=0
      KEEP(204)=0
      KEEP(205)=0
      KEEP(209)=-1
      KEEP(104) = 16
      KEEP(107)=0
#if ! defined(NO_XXNBPR)
      KEEP(121)=-999999
#endif
      KEEP(211)=2
      IF (NSLAVES .EQ. 2) THEN
        KEEP(213) = 101
      ELSE
        KEEP(213) = 201
      ENDIF
      KEEP(217)=0
      KEEP(215)=0
      KEEP(216)=1
      KEEP(218)=50
      KEEP(219)=1
      IF (KEEP(50).EQ.2) THEN
        KEEP(227)= max(2,32)
      ELSE
        KEEP(227)= max(1,32)
      ENDIF
      KEEP(231) = 1
      KEEP(232) = 3
      KEEP(233) = 0
      KEEP(239) = 1
      KEEP(240) = 10
      DKEEP(4) = -1.0D0
      DKEEP(5) = -1.0D0
      DKEEP(10) = 1000.0D0  ! > 0 : GAP 
      IF(NSLAVES.LE.8)THEN
         KEEP(238)=12
      ELSE
         KEEP(238)=7
      ENDIF
      KEEP(234)= 1
      KEEP(235)=-1
      DKEEP(3)=-5.0D0
      KEEP(242) = 1
      KEEP(243) = -1
      KEEP(250) = 1
      KEEP(261) = 1
      KEEP(262) = 0
      KEEP(263) = 0
      KEEP(488) =  4*KEEP(6)   ! if KEEP(6)=32 then 128
      KEEP(496) = 1
      KEEP(495) = -1
      KEEP(497) = -1
      RETURN
      END SUBROUTINE ZMUMPSID
      SUBROUTINE ZMUMPS_SET_KEEP72(id, LP)
      USE ZMUMPS_STRUC_DEF
      IMPLICIT NONE
      TYPE (ZMUMPS_STRUC) :: id
      INTEGER LP
      IF (id%KEEP(72)==1) THEN
         id%KEEP(37) = 2*id%NSLAVES
         id%KEEP(3)=3
         id%KEEP(4)=2
         id%KEEP(5)=1
         id%KEEP(6)=2
         id%KEEP(9)=3
         id%KEEP(39)=300
         id%CNTL(1)=0.1D0  
         id%KEEP(213) = 101  
         id%KEEP(85)=2
         id%KEEP(85)=-4
         id%KEEP(62) = 2  
         id%KEEP(1)  = 1   
         id%KEEP(51) = 2
         id%KEEP(488) = 4
         id%KEEP(490) = 5
         id%KEEP(491) = 5
      ELSE IF (id%KEEP(72)==2) THEN
         id%KEEP(85)=2            ! default is 
         id%KEEP(85)=-10000         ! default is 160 
         id%KEEP(62) = 10       ! default is 50
         id%KEEP(210) = 1       ! defaults is 0 (automatic)
         id%KEEP8(79) = 160000_8   
         id%KEEP(1) = 2  ! default is 8
         id%KEEP(102) = 110     ! defaults is 150 up to 48 procs
         id%KEEP(213) = 121   ! default is 201
      END IF
      RETURN
      END SUBROUTINE ZMUMPS_SET_KEEP72