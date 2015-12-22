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
C  ===========================
C  FORTRAN 90 Driver for ZMUMPS
C       (MPI based code) 
C  ===========================
C
      SUBROUTINE ZMUMPS( id )
      USE ZMUMPS_OOC
      USE ZMUMPS_STRUC_DEF
      USE ZMUMPS_STATIC_PTR_M ! For Schur pointer
C
      IMPLICIT NONE
C
C  =======
C  Purpose
C  =======
C
C  TO SOLVE a SPARSE SYSTEM OF LINEAR EQUATIONS.
C  GIVEN AN UNSYMMETRIC, SYMMETRIC, OR SYMMETRIC POSITIVE DEFINITE 
C  SPARSE MATRIX A AND AN N-VECTOR B, THIS SUBROUTINE SOLVES THE 
C  SYSTEM A x = b or ATRANSPOSE x = b. 
C
C  List of main functionalities provided by the package:
C  ----------------------------------------------------
C        -Unsymmetric solver with partial pivoting (LU factorization)
C        -Symmetric positive definite solver (LDLT factorization)
C        -General symmetric solver with pivoting
C        -Either elemental or assembled matrix input
C        -Analysis/Factorization/Solve callable separately
C        -Deficient matrices (symmetric or unsymmetric) 
C          -Rank revealing 
C          -Null space basis computation
C          -Solution 
C        -Return the Schur complement matrix while 
C          also providing solution of interior problem
C        -Distributed input matrix and analysis phase
C        -Sequential or parallel MPI version (any number of processors)
C        -Error analysis and iterative refinement
C        -Out-of-Core factorization and solution
C        -Solution phase:
C          -Multiple Right-Hand-sides (RHS)
C          -Sparse RHS
C          -Computation of selected entries of the inverse of original matrix.
C
C Method
C ------
C  The method used is a parallel direct method
C  based on a sparse multifrontal variant
C  of Gaussian elimination with partial numerical pivoting. 
C  An initial ordering for the pivotal sequence
C  is chosen using the pattern of the matrix A + A^T and is
C  later modified for reasons of numerical stability.  Thus this code
C  performs best on matrices whose pattern is symmetric, or nearly so.
C  For symmetric sparse matrices or for very unsymmetric and
C  very sparse matrices, other software might be more appropriate.
C
C
C References :
C -----------
C
C  W. M. Sid-Lakhdar, PhD Thesis from Université de Lyon prepared at ENS Lyon,
C   Scaling the solution of large sparse linear systems using multifrontal
C   methods on hybrid shared-distributed memory architectures (2014).
C
C  P. Amestoy, J.-Y. L'Excellent, W. Sid-Lakhdar,
C   Characterizing asynchronous broadcast trees for multifrontal factorizations,
C   Workshop on Combinatorial Scientific Computing,
C   Lyon, France, July 21-23 (2014).
C
C  P. Amestoy, J.-Y. L'Excellent, F.-H. Rouet, W. Sid-Lakhdar,
C   Modeling 1D distributed-memory dense kernels for an asynchronous
C   multifrontal sparse solver, High-Performance Computing for Computational
C   Science, VECPAR 2014, Eugene, Oregon, USA, June 30 - July 3 (2014).
C
C  J.-Y. L'Excellent and W. M. Sid-Lakhdar,
C   Introduction of shared-memory parallelism in a distributed-memroy
C   multifrontal solver, Parallel Computing (40):3-4, pages 34-46 (2014).
C  
C  C. Weisbecker, PhD Thesis supported by EDF, INPT-IRIT,
C   Improving multifrontal solvers by means of algebraic block low-rank
C   representations (2013).
C
C  E. Agullo, P. Amestoy, A. Buttari, A. Guermouche,  G. Joslin, J.-Y.
C  L'Excellent, X. S. Li, A. Napov, F.-H. Rouet, M. Sid-Lakhdar, S. Wang, C.
C  Weisbecker, I. Yamazaki,
C   Recent Advances in Sparse Direct Solvers, 22nd Conference on Structural
C   Mechanics in Reactor Technology, San Francisco (2013).
C
C  P. Amestoy, A. Buttari, G. Joslin, J.-Y. L'Excellent, W. Sid-Lakhdar, C.
C   Weisbecker, M. Forzan, C. Pozza, R. Perrin, V. Pellissier,
C   Shared memory parallelism and low-rank approximation techniques applied
C   applied to direct solvers in FEM simulation in <i>IEEE Transactions on
C   Magnetics</i>, IEEE, Special issue, Compumag 2013 (2013).
C
C  L. Boucher, P. Amestoy, A, Buttari, F.-H. Rouet and M. Chauvin,
C   INTEGRAL/SPI data segmentation to retrieve sources intensity variations,
C   Astronomy & Astrophysics, Article 52, 20 pages,
C   http://dx.doi.org/10.1051/0004-6361/201219605 (2013).
C
C  F.-H. Rouet, PhD thesis from INPT, Toulouse, France,
C   Memory and Performance issues in parallel multifrontal factorization and
C   triangular solutions with sparse right-hand sides (2014).
C
C  J.-Y. L'Excellent, Habilitation thesis from ENS Lyon,
C   Multifrontal methods: Parallelism, Memory Usage and Numerical
C   Aspects (2012).
C
C  P. Amestoy, I.S. Duff, J.-Y. L'Excellent, Y. Robert, F.H. Rouet
C   and B. Ucar, On computing inverse entries of a sparse matrix in 
C   an out-of-core environment,
C   SIAM J. on Scientific Computing Vol. 34 N. 4, p. 1975-1999 (2012).
C
C  Amestoy, Buttari, Duff, Guermouche, L'Excellent, and Uçar
C   The Multifrontal Method, Encyclopedia of Parallel Computing,
C   editor David Padua, Springer (2011).
C
C  Amestoy, Buttari, Duff, Guermouche, L'Excellent, and Uçar
C   MUMPS, Encyclopedia of Parallel Computing,
C   editor David Padua, Springer (2011).
C
C  Agullo, Guermouche and L'Excellent, Reducing the {I/O} Volume in 
C   Sparse Out-of-core Multifrontal Methods}, SIAM SISC, Vol 31, Nb. 6, 
C   4774-4794 (2010).
C
C  Amestoy, Duff, Guermouche, Slavova, Analysis of the Solution Phase of a
C   Parallel Multifrontal Approach, Parallel Computing, Vol. 36, 3--15 (2010). 
C
C  Tzvetomila Slavova, PhD from INPT prepared at CERFACS,
C   Parallel triangular solution in the out-of-core multifrontal approach 
C   for solving large sparse linear systems, available as CERFACS 
C   Report TH/PA/09/59 (2009).
C
C  Agullo, Guermouche and L'Excellent, A Parallel Out-of-core Multifrontal 
C   Method: Storage of Factors on Disk and Analysis of Models for an 
C   Out-of-core Active Memory, Parallel Computing, Special Issue on Parallel 
C   Matrix Algorithms, Vol. 34, Nb 6-8, 296--317 (2008).
C
C  Emmanuel Agullo, PhD Thesis from LIP-Ecole Normale Superieure de Lyon,
C   On the Out-of-core Factorization of Large Sparse Matrices (Nov 2008).
C
C  Amestoy, Duff, Ruiz, and Ucar, "A parallel
C   matrix scaling algorithm".
C   In proceedings of VECPAR'08-International Meeting-High 
C   Performance Computing for Computational Science, (Jan 2008).
C
C  Guermouche and L'Excellent, Constructing Memory-minimizing Schedules 
C   for Multifrontal Methods, ACM TOMS, Vol. 32, Nb. 1, 17--32 (2006).
C
C  Amestoy, Guermouche, L'Excellent, and Pralet, 
C   Hybrid scheduling for the parallel solution
C   of linear systems.  Vol 32 (2), pp 136-156 (2006).
C
C  Stéphane Pralet, PhD from INPT prepared at CERFACS,
C   Constrained orderings and scheduling for parallel sparse linear algebra,
C   available as CERFACS technical report, TH/PA/04/105, (Sept 2004).
C
C  Abdou Guermouche, PhD Thesis from LIP-Ecole Normale Superieure de Lyon,
C   Etude et optimisation du comportement mémoire dans les méthodes paralleles 
C   de factorisation de matrices creuses (2004).
C
C  Guermouche, L'Excellent and Utard, Impact of Reordering on the Memory of a
C   Multifrontal Solver, Parallel Computing, Vol. 29, Nb. 9, 1191--1218 (2003).
C
C  Amestoy, Duff, L'Excellent and Xiaoye S. Li, Impact of the Implementation 
C   of MPI Point-to-Point Communications on the Performance of Two General 
C   Sparse Solvers, Parallel Computing, Vol. 29, Nb 7, 833--847 (2003).
C
C  Amestoy, Duff, L'Excellent and Xiaoye S. Li, Analysis and Comparison of 
C   Two General Sparse Solvers for Distributed Memory Computers, ACM TOMS,
C   Vol. 27, Nb 4, 388--421 (2001).
C
C  Amestoy, Duff, Koster and  L'Excellent (2001), 
C   A fully asynchronous multifrontal solver using distributed dynamic
C   scheduling, SIAM Journal of Matrix Analysis and Applications,
C   Vol 23, No 1, pp 15-41 (2001).
C
C  Amestoy, Duff and  L'Excellent (2000),
C   Multifrontal parallel distributed symmetric and unsymmetric solvers,
C   Comput. Methods in Appl. Mech. Eng., 184,  501-520 (2000)
C
C  Amestoy, Duff and L'Excellent (1998),
C   Parallelisation de la factorisation LU de matrices
C   creuses non-symmetriques pour des architectures a memoire distribuee,
C   Calculateurs Paralleles Reseaux et systemes repartis, 
C   Vol 10(5), 509-520 (1998).
C
C  PARASOL Deliverable D2.1d (final report), 
C   ZMUMPS Version 3.1, A MUltifrontal Massively Parallel Solver,
C   PARASOL project, EU ESPRIT IV LTR project 20160, (June 1999).
C
C  Jacko Koster, PhD from INPT prepared at CERFACS, On the parallel solution 
C   and the reordering of unsymmetric sparse linear systems (1997).
C
C  Vincent Espirat, Master's thesis from INPT(ENSEEIHT)-IRIT, Développement 
C   d'une approche multifrontale pour machines à mémoire distribuée et 
C   réseau hétérogene de stations de travail (1996).
C
C  Patrick Amestoy, PhD from INPT prepared at CERFACS, Factorization of large 
C  sparse matrices based on a multifrontal approach in a multiprocessor 
C  environment, Available as CERFACS report TH/PA/91/2 (1991).
C
C============================================
C Argument lists and calling sequences
C============================================
C
C There is only one entry:
*
*  A Fortran 90 driver subroutine ZMUMPS has been designed as a user
*   friendly interface to the multifrontal code. 
*   This driver, in addition to providing the 
*   normal functionality of a sparse solver, incorporates some
*   pre- and post-processing.
*   This driver enables the user to preprocess the matrix to obtain a 
*   maximum
*   transversal so that the permuted matrix has a zero-free diagonal,
*   to perform prescaling
*   of the original matrix (a choice of scaling strategies is provided),
*   to use iterative refinement to improve the solution,
*   and finally to perform error analysis.
* 
* The driver routine ZMUMPS offers similar functionalities to other 
* sparse direct solvers, depending on the value of one of 
* its parameters (JOB).  These are:
*
* (i)  JOB = -1 
C    initializes an instance of the package. This must be
C    called before any other call to the package concerning that instance.
C    It sets default values for other
C    components of ZMUMPS_STRUC, which may then be altered before
C    subsequent calls to ZMUMPS.
C    Note that three components of the structure must always be set by the
C    user (on all processors) before a call with JOB=-1. These are
C        id%COMM,
C        id%SYM, and
C        id%PAR.
C    CNTL, ICNTL can then be modified (see documentation) by the user.²
C
* A value of JOB = -1 cannot be combined with other values for JOB
*
* (ii) JOB = 1 accepts the pattern of matrix A and chooses pivots
* from the diagonal using a selection criterion to
* preserve sparsity.  It uses the pattern of A + A^T 
* but ignores numerical values. It subsequently constructs subsidiary
* information for the actual factorization by a call with JOB_=_2.  
* An option exists for the user to
* input the pivot sequence, in which case only the necessary
* information for a JOB = 2 entry will be generated.  We call the JOB=1
* entry, the analysis phase.
C The following components of the structure define the centralized matrix 
C pattern and must be set by the user (on the host only) 
C before a call with JOB=1:
C   --- id%N, id%NZ,
C       id%IRN, and id%JCN
C       if the user wishes to input the structure of the
C       matrix in assembled format (ICNTL(5)=0, and ICNTL(18) $\neq$ 3),
C   --- id%ELTPTR, and id%ELTVAR
C       if the user wishes to input the matrix in elemental
C       format (ICNTL(5)=1).
C A distributed matrix format is also available (see documentation)
C
* (iii) JOB = 2 factorizes a matrix A using the information
* from a previous call with JOB = 1. The actual pivot sequence
* used may differ slightly from that of this earlier call if A is not
* diagonally dominant.
*
* (iv) JOB = 3 uses the factors generated by a JOB = 2 call to solve
* a system of equations A X = B or A^T X =B, where X and B are matrices
* that can be either dense or sparse.
* The sparsity of B is exploited to limit the number of operations 
* performed during solution. When only part of the solution is
* also needed (such as when computing selected entries of A^1) then
* further reduction of the number of operations is performed.
* This is particularly beneficial in the context of an 
* out-of-core factorization.
*
* A call with JOB=3 must be preceded by a call with JOB=2,
* which in turn must be preceded by a call with JOB=1.  Since
* the information passed from one call to the next is not
* corrupted by the second, several calls with JOB=2 for matrices
* with the same sparsity pattern but different values may follow
* a single call with JOB=1, and similarly several calls with JOB=3 
* can be used for different right-hand sides.
* Other values for the parameter JOB can invoke combinations of these
* three basic operations.
C
*********
C     --------------------------------------
C     Explicit interface needed for routines
C     using a target argument if they appear
C     in the same compilation unit.
C     --------------------------------------
      INTERFACE
      SUBROUTINE ZMUMPS_CHECK_DENSE_RHS
     &(idRHS, idINFO, idN, idNRHS, idLRHS)
      COMPLEX(kind=8), DIMENSION(:), POINTER :: idRHS
      INTEGER, intent(in)    :: idN, idNRHS, idLRHS
      INTEGER, intent(inout) :: idINFO(:)
      END SUBROUTINE ZMUMPS_CHECK_DENSE_RHS
      SUBROUTINE ZMUMPS_ANA_DRIVER( id )
      USE ZMUMPS_STRUC_DEF
      TYPE (ZMUMPS_STRUC), TARGET :: id
      END SUBROUTINE ZMUMPS_ANA_DRIVER
      SUBROUTINE ZMUMPS_FAC_DRIVER( id )
      USE ZMUMPS_STRUC_DEF
      TYPE (ZMUMPS_STRUC), TARGET :: id
      END SUBROUTINE ZMUMPS_FAC_DRIVER
      SUBROUTINE ZMUMPS_SOLVE_DRIVER( id )
      USE ZMUMPS_STRUC_DEF
      TYPE (ZMUMPS_STRUC), TARGET :: id
      END SUBROUTINE ZMUMPS_SOLVE_DRIVER
      SUBROUTINE ZMUMPS_PRINT_ICNTL(id, LP)
      USE ZMUMPS_STRUC_DEF
      TYPE (ZMUMPS_STRUC), TARGET, INTENT(IN) :: id
      INTEGER  :: LP
      END SUBROUTINE ZMUMPS_PRINT_ICNTL
      END INTERFACE
*  MPI
*  ===
      INCLUDE 'mpif.h'
      INTEGER MASTER, IERR
      PARAMETER( MASTER = 0 )
*
*  ==========
*  Parameters
*  ==========
      TYPE (ZMUMPS_STRUC) :: id
C
C  Main components of the structure are:
C  ------------------------------------
C
C   (see documentation for a complete description)
C
C  JOB is an INTEGER variable which must be set by the user to
C    characterize the factorization step.  Possible values of JOB
C    are given below
C
C     1   Analysis: Ordering and symbolic factorization steps.
C     2   Scaling and Numerical Factorization
C     3   Solve and Error analysis
C     4   Analysis followed by numerical factorization
C     5   Numerical factorization followed by Solving step
C     6   Analysis, Numerical factorization and Solve
C
C  N is an INTEGER variable which must be set by the user to the
C    order n of the matrix A.  It is not altered by the
C     subroutine.  
C
C  NZ is an INTEGER variable which must be set by the user to the
C    number of entries being input.  It is not altered by the
C    subroutine. Only used if ICNTL(5).eq.0 (assembled matrix entry).
C    Restriction: NZ > 0.
C
C  NELT is an INTEGER variable which must be set by the user to the
C    number of elements being input.  It is not altered by the
C    subroutine. Only used if ICNTL(5).ne.0 (elemental matrix entry).
C    Restriction: NELT > 0.
C
C  IRN and JCN  are INTEGER  arrays of length NZ. 
C    IRN(k) and JCN(k), k=1..NZ must be set on entry to hold 
C    the row and column indices respectively.
C    They are not altered by the subroutine except when ICNTL(6) = 1.
C    (in which case only the column indices are modified).
C    The arrays are only used if ICNTL(5).eq.0 (assembled entry).
C
C  ELTPTR is an INTEGER array of length NELT+1. 
C  ELTVAR is an INTEGER array of length ELTPTR(NELT+1)-1.
C    ELTPTR(I) points in ELTVAR to the first variable in the list of
C    variables that correspond to element I. ELTPTR(NELT+1) points
C    to the first unused location in ELTVAR.
C    The positions ELTVAR(I) .. ELTPTR(I+1)-1 contain the variables
C    for element I. No free space is allowed between variable lists.
C    ELTPTR/ELTVAR are not altered by the subroutine.
C    The arrays are only used if ICNTL(5).ne.0 (element entry).
C
C  A is a COMPLEX(kind=8) array of length NZ. 
C     The user must set A(k) to the value 
C     of the entry in row IRN(k) and column JCN(k) of the matrix.
C     It is not altered by the subroutine.
C     (Note that the matrix can also be provided in a distributed 
C      assembled input format)
C
C  RHS is a COMPLEX(kind=8) array of length N that is only accessed when
C    JOB = 3, 5, or 6. On entry, RHS(i)
C     must hold the i th component of the right-hand side of the
C     equations being solved.
C     On exit, RHS(i) will hold the i th component of the
C     solution vector.  For other values of JOB, RHS is not accessed and
C     can be declared to have size one.
C     RHS should only be available on the host processor. If
C     it is associated on other processors, an error is raised.
C     (Note that the right-hand sides can also be provided in a 
C      sparse format).
C
C COLSCA, ROWSCA are DOUBLE PRECISION
C     arrays of length N that are used to hold
C     the values used to scale the columns and the rows
C     of the original matrix, respectively. 
C     These arrays need to be set by the user
C     only if ICNTL(8) is set to -1. If ICNTL(8)=0,
C     COLSCA and ROWSCA are not accessed and 
C     so can be declared to have size one.
C     For any other values of ICNTL(8),
C     the scaling arrays are computed before
C     numerical factorization.  The factors of the scaled matrix
C     diag(ROWSCA(i)) <A diag(COLSCA(i)) are computed.
C 
C  The workspace is automatically allocated by the package.
C  At the beginning of the numerical phase. If the user wants to increase
C   the allocated workspace (typically, numerical pivoting that leads to extra
C   storage, or previous call to MUMPS that failed because of 
C   a lack of allocated memory), 
C   we describe in the following how the user can modify the size 
C   of the workspace:
C    1/ The memory relaxation parameter
C       ICNTL(14) is designed to control the increase, with respect to the 
C       estimations performed during analysis, in the size of the workspace 
C       allocated during the numerical phase.
C    2/ The user can also provide 
C       a unique parameter,  ICNTL(23),  holding the maximum size of the total 
C       workspace (in Megabytes) that the package is allowed to use internally.
C       In this case we try as much as possible to follow the indication given
C       by the relaxation parameter (ICNTL(14)).
C
C   If ICNTL(23) is greater than 0 
C   then MUMPS automatically computes the size of the internal working arrays
C   such that the storage for all MUMPS internal data is equal to ICNTL(23).
C   The relaxation ICNTL(14) is first applied to
C   the internal integer working array and communication buffer sizes;
C   the remaining available space is given to the real/complex 
C   internal working arrays.
C   A lower bound of ICNTL(23) (if ICNTL(14) has not
C   been modified since the analysis) is given by INFOG(26).
C   
C   If ICNTL(23) is left to its default value 0 
C   then each processor will allocate workspace based on
C   the estimates computed during the analysis (INFO(17)
C   if ICNTL(14) has not been modified since analysis,
C   or larger if ICNTL(14) was increased). 
C   Note that these estimates are accurate in the sequential
C   version of {\tt MUMPS}, but that they can be inaccurate
C   in the parallel case. Therefore, in parallel, we recommend
C   to use ICNTL(23) and provide a value significantly larger
C   than INFOG(26).
C --------------------------------------------------------------------------------
C    
C CNTL is a DOUBLE PRECISION array of length 15
C  that contains control parameters and must be set by the user. Default
C  values for the components may be set by a call to ZMUMPS(JOB=-1)
C  Details of the control parameters are given in ZMUMPSID.
C
C ICNTL is an INTEGER array of length 40
C  that contains control parameters and must be set by the user. Default
C  values for the components may be set by a call to ZMUMPS(JOB=-1)
C  Details of the control parameters are given in ZMUMPSID.
C
C INFO is an INTEGER array of length 40 that need not be set by the
C  user.  On return from ZMUMPS, a value of zero for INFO(1)
C  indicates that the subroutine has performed successfully. 
C  Details of the control parameters are given in ZMUMPSID.
C 
C RINFO is a DOUBLE PRECISION  array of length 40 that need not be set by the
C  user.  This array supplies information on the execution of ZMUMPS.
C  Details of the control parameters are given in ZMUMPSID.
C
C
*
*
*   ====================
*    .. Error Return ..
*   ====================
*
C MUMPS uses the following mechanism to process errors that
C may occur during the parallel execution of the code. 
C If, during a call to MUMPS, an error occurs on a processor, 
C this processor informs all the other processors before they
C return from the call.
C In parts of the code where messages are sent asynchronously 
C (for example the factorization and solve phases), 
C the processor on which the error occurs sends a message 
C to the other processors with a specific error tag. 
C On the other hand, if the error occurs in a subroutine that
C does not use asynchronous communication, the processor propagates 
C the error to the other processors.
C On successful completion, a call to MUMPS will exit with the 
C parameter id%INFOG(1) set to zero.
C A negative value for id%INFOG(1) indicates that an 
C error has been detected on one of the processors.
C For example, if processor s returns with
C INFO(1)= -8 and INFO(2)=1000, then processor s ran out of integer 
C workspace during the factorization and the size of the workspace 
C should be increased by 1000 at least. 
C The other processors are informed about this error and return with
C INFO(1)=-1 (i.e., an error occurred on another processor) and 
C INFO(2)=s (i.e., the error occurred on processor s).
C If several processors raised an error, those processors do not overwrite 
C INFO(1), i.e., only processors that did not produce an error will set 
C INFO(1) to -1 and INFO(2) to the rank of the processor having the most 
C negative error code.
C
C The behaviour is slightly different for the global information
C parameters INFOG(1) and INFOG(2):
C in the previous example, all processors would return with
C INFOG(1)=-8 and INFOG(2)=1000.
C
C The possible error codes returned in INFO(1) (and INFOG(1))
C are fully described in the documentation.
C
C A positive value of INFO(1) is associated with a warning message 
C which  will be output on unit ICNTL(2) (see documentation).
C
C
C      .. Local variables ..
C
      INTEGER JOBMIN, JOBMAX, OLDJOB
      INTEGER I, J, MP, LP, MPG, KEEP235SAVE, KEEP242SAVE,
     &        KEEP243SAVE, KEEP495SAVE, KEEP497SAVE
      LOGICAL LANA, LFACTO, LSOLVE, PROK, LPOK, FLAG, PROKG
      LOGICAL NOERRORBEFOREPERM
      LOGICAL UNS_PERM_DONE
C     Saved communicator (pb of interference)
      INTEGER COMM_SAVE
C     Local copies of fields JOB, N, NZ, NELT in the structure
      INTEGER JOB, N, NZ, NELT
      INTEGER, PARAMETER :: ICNTL18DIST_MIN = 1
      INTEGER, PARAMETER :: ICNTL18DIST_MAX = 3
      INTEGER, DIMENSION(:), ALLOCATABLE :: UNS_PERM_INV
C     TIMINGS
      DOUBLE PRECISION TIMEG, TIMETOTAL
      NOERRORBEFOREPERM = .FALSE.
      UNS_PERM_DONE = .FALSE.
      JOB  = id%JOB
      N    = id%N
      NZ   = id%NZ
      NELT = id%NELT
C
C     Initialize error return codes to 0.
      id%INFO(1) = 0
      id%INFO(2) = 0
C     -----------------------------------
C     Check that MPI has been initialized
C     -----------------------------------
      CALL MPI_INITIALIZED( FLAG, IERR )
      IF ( .NOT. FLAG ) THEN
        id%INFO(1) = -23
        id%INFO(2) =   0
        WRITE(6,990)
 990  FORMAT(' Unrecoverable Error in ZMUMPS initialization: ',
     &       ' MPI is not running.')
        RETURN               
      END IF
C     ---------------------------
C     Duplicate user communicator
C     to avoid communications not
C     related to ZMUMPS
C     ---------------------------
       COMM_SAVE = id%COMM
       CALL MPI_COMM_DUP( COMM_SAVE, id%COMM, IERR )
C     -------------------------
C     Check if value of JOB is
C     the same on all processes
C     -------------------------
      CALL MPI_ALLREDUCE(JOB,JOBMIN,1,MPI_INTEGER,MPI_MAX,
     &                   id%COMM,IERR)
      CALL MPI_ALLREDUCE(JOB,JOBMAX,1,MPI_INTEGER,MPI_MIN,
     &                   id%COMM,IERR)
      IF ( JOBMIN .NE. JOBMAX ) THEN
        id%INFO(1) = -3 
        id%INFO(2) = JOB
        GOTO 499
      END IF
C
C     Default setting for printing
      LP = 6
      MP = 0
      MPG = 0
      LPOK  = .TRUE.
      PROK  = .FALSE.
      PROKG = .FALSE.
C   
C     Check value of JOB and previous value of JOB
C
      IF (JOB.LT.-2.OR.JOB.EQ.0.OR.JOB.GT.6) THEN
C       Out of range value
        id%INFO(1) = -3 
        id%INFO(2) = JOB
        GOTO 499
      END IF
      IF (JOB.NE.-1) THEN
C      Check the previous value of JOB
C      One should be able to test for old job value
       OLDJOB = id%KEEP( 40 ) + 456789
       IF (OLDJOB.NE.-1.AND.OLDJOB.NE.-2.AND.
     &    OLDJOB.NE.1.AND.OLDJOB.NE.2.AND.
     &    OLDJOB.NE.3) THEN
        id%INFO(1) = -3 
        id%INFO(2) = JOB
        GOTO 499
       END IF
      END IF
C
C
      IF (JOB.EQ.-2.OR.JOB.EQ.1.OR.JOB.EQ.2.OR.JOB.EQ.3.OR.
     &    JOB.EQ.4.OR.JOB.EQ.5.OR.JOB.EQ.6) THEN
C       Correct value of JOB
C       ICNTL should have been initialized and can be used
        LP      = id%ICNTL(1)
        MP      = id%ICNTL(2)
        MPG     = id%ICNTL(3)
        LPOK    = ((LP.GT.0).AND.(id%ICNTL(4).GE.1))
        PROK    = ((MP.GT.0).AND.(id%ICNTL(4).GE.2))
        PROKG   = ( MPG .GT. 0 .and. id%MYID .eq. MASTER )
        PROKG   = (PROKG.AND.(id%ICNTL(4).GE.2))
        IF (PROKG) THEN
C         Print basic information on matrix
           IF (id%ICNTL(5) .NE. 1) THEN
              IF (id%ICNTL(18) .EQ. 0) THEN
                WRITE(MPG,'(A,A,A,I4,I12,I15)') 
     &                'Entering ZMUMPS ',
     &                trim(adjustl(id%VERSION_NUMBER)),
     &                ' driver with JOB, N, NZ =', JOB,N,NZ
              ELSE
                WRITE(MPG,'(A,A,A,I4,I12,I15)') 
     &                'Entering ZMUMPS ',
     &                trim(adjustl(id%VERSION_NUMBER)),
     &                ' driver with JOB, N =', JOB,N
              ENDIF
           ELSE
              WRITE(MPG,'(A,A,A,I4,I12,I15)') 
     &             'Entering ZMUMPS ',
     &             trim(adjustl(id%VERSION_NUMBER)),
     &             ' driver with JOB, N, NELT =', JOB,N
     &             ,NELT
           ENDIF
        ENDIF
      END IF
C
C----------------------------------------------------------------
C
C     JOB = -1 : START INITIALIZATION PHASE
C                (NEW INSTANCE)
C
C     JOB = -2 : TERMINATE AN INSTANCE
C----------------------------------------------------------------
C
      IF ( JOB .EQ. -1 ) THEN
C
C       ------------------------------------------
C       Check that we have called (JOB=-2), ie
C       that the previous JOB is not 1 2 or 3,
C       before calling the initialization routine.
C       --------------------------------------------
        id%INFO(1)=0
        id%INFO(2)=0
        OLDJOB = id%KEEP( 40 ) + 456789
        IF ( OLDJOB .EQ. 1 .OR.
     &       OLDJOB .EQ. 2 .OR.
     &       OLDJOB .EQ. 3  ) THEN
          IF ( id%N > 0 ) THEN
           id%INFO(1)=-3
           id%INFO(2)=JOB
          ENDIF
        ENDIF
C       Initialize id%MYID now because it is
C       required by MUMPS_PROPINFO. id%MYID
C       used to be initialized inside ZMUMPS_INI_DRIVER,
C       leading to an uninitialized access here.
        CALL MPI_COMM_RANK(id%COMM, id%MYID, IERR)
        CALL MUMPS_PROPINFO( id%ICNTL(1),
     &                       id%INFO(1),
     &                       id%COMM, id%MYID )
        IF ( id%INFO(1) .LT. 0 ) THEN
C
C         If there was an error, then initialization
C         was already called and we can rely on the null
C         or non null value of the pointers related to OOC
C         stuff.
C         We use ZMUMPS_CLEAN_OOC_DATA that should work even
C         on the master. Note that KEEP(201) was also
C         initialized in a previous call to Mumps.
C
C         If ZMUMPS_END_DRIVER or ZMUMPS_FAC_DRIVER is called after
C         this error, then ZMUMPS_CLEAN_OOC_DATA will be called
C         a second time, though.
C
           IF (id%KEEP(201).GT.0) THEN
             CALL ZMUMPS_CLEAN_OOC_DATA(id, IERR)
           ENDIF
           GOTO 499
        ENDIF
C       ----------------------------------------
C       Initialization ZMUMPS_INI_DRIVER 
C       ----------------------------------------
C       - Default values for ICNTL, KEEP,KEEP8, CNTL
C       - Attach emission buffer for buffered Send
C       - Nullify pointers in the structure
C       - Get rank and size of the communicator
C       ----------------------------------------
        CALL ZMUMPS_INI_DRIVER( id )
        GOTO 500
      END IF
      IF ( JOB .EQ. -2 ) THEN
C       -------------------------------------
C       Deallocation of the instance id
C       -------------------------------------
        id%KEEP(40)= -2 - 456789
        CALL ZMUMPS_END_DRIVER( id )
        GOTO 500
      END IF
C
C----------------------------------------------------------------
C
C     MAIN DRIVER
C     OTHER VALUES OF JOB : 1 to 6
C
C----------------------------------------------------------------
C
C     Check some input parameters first
C     ---------------------------------
C
C     Check N on the master only
      IF (id%MYID.EQ.MASTER) THEN
        IF ((N.LE.0).OR.((N+N+N)/3.NE.N)) THEN
          id%INFO(1) = -16
          id%INFO(2) = N
        END IF
C       Checking only for centralized matrix
C       This is done by testing non-equality to
C       distributed options (1..3)
        IF ( id%ICNTL(18) .LT. ICNTL18DIST_MIN
     &       .OR. id%ICNTL(18) .GT. ICNTL18DIST_MAX ) THEN
          IF (id%ICNTL(5).NE.1) THEN
            IF (NZ.LE.0) THEN
              id%INFO(1) = -2
              id%INFO(2) = NZ
            END IF
          ELSE
C           Element entry: check NELT on the master only
            IF (NELT.LE.0) THEN
              id%INFO(1) = -24
              id%INFO(2) = NELT
            END IF
          ENDIF
        END IF
C       -----------------------------
C       Check incompatibility between
C       par (=0) and nprocs (=1)
C       -----------------------------
        IF ( (id%KEEP(46).EQ.0).AND.(id%NPROCS.LE.1) ) 
     &     THEN
          id%INFO(1) = -21
          id%INFO(2) = id%NPROCS
        ENDIF
      END IF
C
C     Propagate possible error to all nodes
      CALL MUMPS_PROPINFO( id%ICNTL(1),
     &                    id%INFO(1),
     &                    id%COMM, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) GOTO 499
C
C     ----------------------------------
C     Initialize, LANA, LFACTO, LSOLVE
C     LANA indicates if analysis must be performed
C     LFACTO indicates if factorization must be performed
C     LSOLVE indicates if solution must be performed
C     ----------------------------------
      LANA  = .FALSE.
      LFACTO = .FALSE.
      LSOLVE = .FALSE.
      IF ((JOB.EQ.1).OR.(JOB.EQ.4).OR.
     &    (JOB.EQ.6))               LANA  = .TRUE.
      IF ((JOB.EQ.2).OR.(JOB.EQ.4).OR.
     &    (JOB.EQ.5).OR.(JOB.EQ.6)) LFACTO = .TRUE.
      IF ((JOB.EQ.3).OR.(JOB.EQ.5).OR.
     &    (JOB.EQ.6))               LSOLVE = .TRUE.
C
C
C     Print ICNTL and KEEP
C
      IF (PROK) CALL ZMUMPS_PRINT_ICNTL(id, MP)
C-----------------------------------------------------------------------
C
C           CHECK SEQUENCE
C
C-----------------------------------------------------------------------
C     TIMINGS
      IF (id%MYID .eq. MASTER) THEN
         id%DKEEP(70)=0.0D0
         CALL MUMPS_SECDEB(TIMETOTAL)
      END IF 
      OLDJOB = id%KEEP( 40 ) + 456789
      IF ( LANA ) THEN
        IF ( PROKG .AND. OLDJOB .EQ. -1 ) THEN
C         Print compilation options at first call to analysis
          CALL  MUMPS_PRINT_IF_DEFINED(MPG)
        ENDIF
C
C       User wants to perform analysis. Previous value of
C       JOB must be -1, 1, 2 or 3.
C
        IF ( OLDJOB .EQ. 0 .OR. OLDJOB .GT. 3 .OR. OLDJOB .LT. -1 ) THEN
          id%INFO(1) = -3
          id%INFO(2) = JOB
          GOTO 499
        END IF
        IF ( OLDJOB .GE. 2 ) THEN
C         -----------------------------------------
C         Previous step was factorization or solve.
C         As analysis is now performed, deallocate
C         at least some big arrays from facto.
C         -----------------------------------------
          IF (associated(id%IS)) THEN
            DEALLOCATE  (id%IS)
            NULLIFY     (id%IS)
          END IF
          IF (associated(id%S)) THEN
            DEALLOCATE  (id%S)
            NULLIFY     (id%S)
          END IF
        END IF   
      END IF
      IF ( LFACTO ) THEN
C        ------------------------------------
C        User wants to perform factorization.
C        Analysis must have been performed.
C        ------------------------------------
         IF ( OLDJOB .LT. 1 .and. .NOT. LANA ) THEN
            id%INFO(1) = -3
            id%INFO(2) = JOB
            GOTO 499
         END IF
      END IF
      IF ( LSOLVE ) THEN
C        -------------------------------
C        User wants to perform solve.
C        Facto must have been performed.
C        -------------------------------
         IF ( OLDJOB .LT. 2 .AND. .NOT. LFACTO ) THEN
            id%INFO(1) = -3
            id%INFO(2) = JOB
            GOTO 499
         END IF
      END IF
C     ------------------------------------------
C     Permute JCN on entry to JOB if no analysis
C     to be performed and IRN/JCN are needed.
C     (facto: arrowheads + solve: iterative
C      refinement and error analysis)
C     ------------------------------------------
#if ! defined (LARGEMATRICES)
      NOERRORBEFOREPERM =.TRUE.
      UNS_PERM_DONE=.FALSE.
      IF (id%MYID .eq. MASTER .AND. id%KEEP(23) .NE. 0) THEN
        IF ( id%JOB .EQ. 2 .OR. id%JOB .EQ. 5 .OR.
     &       (id%JOB .EQ. 3 .AND. (id%ICNTL(10) .NE.0 .OR.
     &        id%ICNTL(11).NE. 0))) THEN
          UNS_PERM_DONE = .TRUE.
          ALLOCATE(UNS_PERM_INV(id%N),stat=IERR)
          IF (IERR .GT. 0) THEN
C             --------------------------------
C             Exit with an error.
C             We are not able to permute
C             JCN correctly after a MAX-TRANS
C             permutation resulting from a
C             previous call to ZMUMPS.
C             --------------------------------
              id%INFO(1)=-13
              id%INFO(2)=id%N
              IF (LPOK) WRITE(LP,99993)
              GOTO 510
          ENDIF
          DO I = 1, id%N
            UNS_PERM_INV(id%UNS_PERM(I))=I
          END DO
          DO I = 1, id%NZ
            J = id%JCN(I)
C           -- skip out-of range (that are ignored in ANA_O)
            IF (J.LE.0.OR.J.GT.id%N) CYCLE
            id%JCN(I)=UNS_PERM_INV(J)
          END DO
          DEALLOCATE(UNS_PERM_INV)
        END IF
      END IF
#endif
C
C       Propagate possible error
        CALL MUMPS_PROPINFO( id%ICNTL(1),
     &                    id%INFO(1),
     &                    id%COMM, id%MYID )
        IF ( id%INFO( 1 ) .LT. 0 ) GO TO 499
*
*********
* MaxTrans-Analysis-Distri, Scale-Arrowhead-factorize, and
* Solve-IR-Error_Analysis (depending on the value of JOB)
*********
*
C
      IF ( LANA ) THEN
C-----------------------------------------------------
C-
C-       ANALYSIS : Max-Trans, Analysis, Distribution
C-
C-----------------------------------------------------
C
C        Allocations
C
C        IS1 :allocated on the master now, will be allocated on
C             the slaves later
C        IS : will be allocated on the slaves later
C        PROCNODE : on the master only,
C             because slave does not know N yet.
C             Will be allocated in analysis for the slave.
C
C        For assembled entry: 
C        IRN, JCN : check that they have been allocated by the
C             user on the master, and if their size is adequate
C
C        For element entry:
C        ELTPTR, ELTVAR : check that they have been allocated by the
C             user on the master, and if their size is adequate
C       ----------------------------
C       Reset KEEP(40) to -1 for the
C       case where an error occurs
C       ----------------------------
        id%KEEP(40)=-1 -456789
C
        IF (id%MYID.EQ.MASTER) THEN
C     -- initialize values of respectively
C     icntl(6), (7) and (12) to not done/chosen
          id%INFOG(7) = -9999
          id%INFOG(23) = 0
          id%INFOG(24) = 1
          IF (associated(id%IS1)) DEALLOCATE(id%IS1)
C         -------------------------------------------
C         Allocate array IS1 for analysis of size:
C          - assembled entry: 10 * N or 11 * N
C                             depending on max-trans
C          - element entry: 7 * N + 3 * NELT + 3
C                           max-trans not allowed
C         -------------------------------------------
          IF ( id%ICNTL(5) .NE. 1 ) THEN ! assembled matrix
            IF ( id%KEEP(50) .NE. 1 
     &           .AND. (
     &           (id%ICNTL(6) .NE. 0 .AND. id%ICNTL(7) .NE.1)
     &           .OR.
     &           id%ICNTL(12) .NE. 1) ) THEN
              id%MAXIS1 = 11 * N
            ELSE
              id%MAXIS1 = 10 * N
            END IF
          ELSE
            id%MAXIS1 = 6 * N + 2 * NELT + 2
          ENDIF
          ALLOCATE( id%IS1(id%MAXIS1), stat=IERR )
          IF (IERR.gt.0) THEN
            id%INFO(1) = -7
            id%INFO(2) = id%MAXIS1
            IF ( LPOK ) WRITE(LP,'(A)')
     &      ' Problem in allocating work array for analysis'
            GO TO 100
          END IF
C
C         ----------------------
C         Allocate PROCNODE(1:N)
C         ----------------------
          IF ( associated( id%PROCNODE ) )
     &         DEALLOCATE( id%PROCNODE )
          ALLOCATE( id%PROCNODE(id%N), stat=IERR )
          IF (IERR.gt.0) THEN
            id%INFO(1) = -7
            id%INFO(2) = id%N
            IF ( LPOK ) WRITE(LP,'(A)')
     &        'Problem in allocating work array PROCNODE'
            GOTO 100
          END IF
          id%PROCNODE(1:id%N) = 0
C         ---------------------------------------
C         Element entry: allocate ELTPROC(1:NELT)
C         ---------------------------------------
          IF ( id%ICNTL(5) .EQ. 1 ) THEN ! Elemental matrix
            IF ( associated( id%ELTPROC ) )
     &           DEALLOCATE( id%ELTPROC )
            ALLOCATE( id%ELTPROC(id%NELT), stat=IERR )
            IF (IERR.gt.0) THEN
              id%INFO(1) = -7
              id%INFO(2) = id%NELT
              IF ( LPOK ) WRITE(LP,'(A)')
     &          'Problem in allocating work array ELTPROC'
              GOTO 100
            END IF
          END IF
C         ---------------------------------------------------
C         Assembled centralized entry: check input parameters
C         IRN/JCN
C         Element entry: check input parameters ELTPTR/ELTVAR
C         ---------------------------------------------------
          IF ( id%ICNTL(5) .NE. 1 ) THEN ! Assembled matrix
            id%NA_ELT=0
            IF ( id%ICNTL(18) .LT. ICNTL18DIST_MIN
     &           .OR. id%ICNTL(18) .GT. ICNTL18DIST_MAX ) THEN
              IF ( .not. associated( id%IRN ) ) THEN
                id%INFO(1) = -22
                id%INFO(2) = 1
              ELSE IF ( size( id%IRN ) < id%NZ ) THEN
                id%INFO(1) = -22
                id%INFO(2) = 1
              ELSE IF ( .not. associated( id%JCN ) ) THEN
                id%INFO(1) = -22
                id%INFO(2) = 2
              ELSE IF ( size( id%JCN ) < id%NZ ) THEN
                id%INFO(1) = -22
                id%INFO(2) = 2
              END IF
            END IF
            IF ( id%INFO( 1 ) .eq. -22 ) THEN
              IF ( LPOK ) WRITE(LP,'(A)')
     &           'Error in analysis: IRN/JCN badly allocated.'
            END IF
          ELSE
            IF ( .not. associated( id%ELTPTR ) ) THEN
              id%INFO(1) = -22
              id%INFO(2) = 1
            ELSE IF ( size( id%ELTPTR ) < id%NELT+1 ) THEN
              id%INFO(1) = -22
              id%INFO(2) = 1
            ELSE IF ( .not. associated( id%ELTVAR ) ) THEN
              id%INFO(1) = -22
              id%INFO(2) = 2
            ELSE 
              id%LELTVAR = id%ELTPTR( id%NELT+1 ) - 1
              IF ( size( id%ELTVAR ) < id%LELTVAR ) THEN 
                id%INFO(1) = -22
                id%INFO(2) = 2
              ELSE
C               If no error, we compute NA_ELT, required
C               for ZMUMPS_MAX_MEM already in analysis, and
C               then later during facto to check the size of A_ELT
                id%NA_ELT = 0
                IF ( id%KEEP(50) .EQ. 0 ) THEN
C                 Unsymmetric elements (but symmetric structure)
                  DO I = 1,NELT
                    J = id%ELTPTR(I+1) - id%ELTPTR(I)
                    J = (J * J)
                    id%NA_ELT = id%NA_ELT + J
                  ENDDO
                ELSE
C                 Symmetric elements
                  DO I = 1,NELT
                    J = id%ELTPTR(I+1) - id%ELTPTR(I)
                    J = (J * (J+1))/2
                    id%NA_ELT = id%NA_ELT + J
                  ENDDO
                ENDIF
              ENDIF
            END IF
            IF ( id%INFO( 1 ) .eq. -22 ) THEN
              IF ( LPOK ) WRITE(LP,'(A)')
     &           'Error in analysis: ELTPTR/ELTVAR badly allocated.'
            END IF
          ENDIF
 100      CONTINUE
        END IF
C
C       Propagate possible error
        CALL MUMPS_PROPINFO( id%ICNTL(1),
     &                    id%INFO(1),
     &                    id%COMM, id%MYID )
        IF ( id%INFO( 1 ) .LT. 0 ) GO TO 499
C       -----------------------------------------
C       Call analysis procedure ZMUMPS_ANA_DRIVER
C       -----------------------------------------
        IF (id%MYID .eq. MASTER) THEN
          id%DKEEP(71)=0.0D0
          CALL MUMPS_SECDEB(TIMEG)
        END IF 
C        -------------------------------------------
C        Set scaling option for analysis in KEEP(52)
C        -------------------------------------------
        id%KEEP(52) = id%ICNTL(8)
C        Out-of-range values => automatic choice
        IF ( id%KEEP(52) .GT. 8 .OR. id%KEEP(52).LT.-2)
     &       id%KEEP(52) = 77
        IF ( id%KEEP(52) .EQ. 2 .OR. id%KEEP(52).EQ.5 
     &       .OR. id%KEEP(52) .EQ. 6 )
     &       id%KEEP(52) = 77
        IF ((id%KEEP(52).EQ.77).AND.(id%KEEP(50).EQ.1)) THEN
         ! for SPD matrices default is no scaling
          id%KEEP(52) = 0
        ENDIF
        IF ( id%KEEP(52).EQ.77 .OR. id%KEEP(52).LE.-2) THEN 
C          -- suppress scaling computed during analysis 
C          -- if centralized matrix is not associated
          IF (.not.associated(id%A)) id%KEEP(52) = 0
        ENDIF
C deactivate analysis scaling if scaling given
        IF(id%KEEP(52) .EQ. -1) id%KEEP(52) = 0
        CALL ZMUMPS_ANA_DRIVER( id )
C Save scaling option in INFOG(33)
        IF (id%MYID .eq. MASTER) THEN
          IF (id%KEEP(52) .NE. 0) THEN
            id%INFOG(33)=id%KEEP(52)
          ELSE
            id%INFOG(33)=id%ICNTL(8)
          ENDIF
        ENDIF
C       return value of ICNTL(12) effectively used
C       that was saved on the master in KEEP(95)
        IF (id%MYID .eq. MASTER) id%INFOG(24)=id%KEEP(95)
C       TIMINGS:
        IF (id%MYID .eq. MASTER) THEN
          CALL MUMPS_SECFIN(TIMEG)
          id%DKEEP(71) = TIMEG
        ENDIF
        IF (PROKG) THEN
          WRITE( MPG,'(A,F12.4)')
     &         ' ELAPSED TIME IN ANALYSIS DRIVER= ', TIMEG
        END IF 
C       -----------------------
C     Return in case of error
C     -----------------------
        IF ( id%INFO( 1 ) .LT. 0 ) GO TO 499
        id%KEEP(40) = 1 -456789
      END IF
C
      IF ( LFACTO ) THEN
         IF (id%MYID .eq. MASTER) THEN
            id%DKEEP(91)=0.0D0
            CALL MUMPS_SECDEB(TIMEG)
         END IF 
C        ----------------------
C        Reset KEEP(40) to 1 in
C        case of error in facto
C        ----------------------
         id%KEEP(40) = 1 - 456789
C
C-------------------------------------------------------
C-
C-      CHECKS, SCALING, ARROWHEAD + FACTORIZATION PHASE
C-
C-------------------------------------------------------
C
        IF ( id%MYID .EQ. MASTER ) THEN
C         -------------------------
C         Check if Schur complement
C         is allocated.
C         -------------------------
          IF (id%KEEP(60).EQ.1) THEN
             IF ( associated( id%SCHUR_CINTERFACE)) THEN
C              Called from C interface...
C              Since id%SCHUR_CINTERFACE is of size 1,
C              instruction below which causes bound check
C              errors should be avoided. We cheat by first
C              setting a static pointer with a routine with
C              implicit interface, and then copying this pointer
C              into id%SCHUR.
               CALL ZMUMPS_SET_TMP_PTR(id%SCHUR_CINTERFACE(1),
     &         id%SIZE_SCHUR*id%SIZE_SCHUR)
               CALL ZMUMPS_GET_TMP_PTR(id%SCHUR)
             ENDIF
             IF ( .NOT. associated (id%SCHUR)) THEN
              IF (LP.GT.0) 
     &        write(LP,'(A)') 
     &                      ' SCHUR not associated'
              id%INFO(1)=-22
              id%INFO(2)=9
             ELSE IF ( size(id%SCHUR) .LT.
     &                id%SIZE_SCHUR * id%SIZE_SCHUR ) THEN
                IF (LP.GT.0) 
     &          write(LP,'(A)') 
     &                ' SCHUR allocated but too small' 
                id%INFO(1)=-22
                id%INFO(2)=9
             END IF
          END IF
C         ------------------------------------------
C         Assembled entry: check input parameter A
C         Element entry: check input parameter A_ELT
C         ------------------------------------------
          IF ( id%KEEP(55) .EQ. 0 ) THEN
           IF ( id%KEEP(54).eq.0 ) THEN
            IF ( .not. associated( id%A ) ) THEN
              id%INFO( 1 ) = -22
              id%INFO( 2 ) = 4
            ELSE IF ( size( id%A ) < id%NZ ) THEN
              id%INFO( 1 ) = -22
              id%INFO( 2 ) = 4
            END IF
           END IF
          ELSE
            IF ( .not. associated( id%A_ELT ) ) THEN
              id%INFO( 1 ) = -22
              id%INFO( 2 ) = 4
            ELSE 
              IF ( size( id%A_ELT ) < id%NA_ELT ) THEN
                id%INFO( 1 ) = -22
                id%INFO( 2 ) = 4
              ENDIF
            END IF
          ENDIF
C         ----------------------
C         Get the value of PERLU
C         ----------------------
          CALL MUMPS_GET_PERLU(id%KEEP(12),id%ICNTL(14),
     &         id%KEEP(50),id%KEEP(54),id%ICNTL(6),id%ICNTL(8))
C
C         ----------------------
C         Get null space options
C         Note that nullspace is forbidden in case of Schur complement
C         ----------------------
          CALL ZMUMPS_GET_NS_OPTIONS_FACTO(N,id%KEEP(1),id%ICNTL(1),MPG)
C         ========================================
C         Decode and set scaling options for facto
C         ========================================
          IF( id%KEEP(52) .EQ. -2 .AND. id%ICNTL(8) .NE. -2 .AND.
     &        id%ICNTL(8).NE. 77 ) THEN
             IF ( MPG .GT. 0 ) THEN
                WRITE(MPG,'(A)') ' ** WARNING : SCALING'
                WRITE(MPG,'(A)') 
     &               ' ** scaling already computed during analysis'
                WRITE(MPG,'(A)') 
     &               ' ** keeping the scaling from the analysis'
             ENDIF
          ENDIF
          IF (id%KEEP(52) .NE. -2) THEN
            id%KEEP(52)=id%ICNTL(8)
          ENDIF
          IF ( id%KEEP(52) .GT. 8 .OR. id%KEEP(52).LT.-2)
     &    id%KEEP(52) = 77
          IF ( id%KEEP(52) .EQ. 2 .OR. id%KEEP(52).EQ.5 
     &        .OR. id%KEEP(52) .EQ. 6 )
     &        id%KEEP(52) = 77
          IF (id%KEEP(52).EQ.77) THEN
            IF (id%KEEP(50).EQ.1) THEN
              ! for SPD matrices the default is "no scaling"
              id%KEEP(52) = 0
            ELSE
              ! SYM .ne. 1  the default is cheap SIMSCA
              id%KEEP(52) = 7 
            ENDIF
          ENDIF
          IF (id%KEEP(23) .NE. 0 .AND. id%ICNTL(8) .EQ. -1) THEN
             IF ( MPG .GT. 0 ) THEN
                WRITE(MPG,'(A)') ' ** WARNING : SCALING'
                WRITE(MPG,'(A)') 
     &               ' ** column permutation applied:'
                WRITE(MPG,'(A)') 
     &               ' ** column scaling has to be permuted'
             ENDIF 
          ENDIF
C
          IF ( id%KEEP( 19 ) .ne. 0 .and. id%KEEP( 52 ).ne. 0 ) THEN
            IF ( MPG .GT. 0 ) THEN
              WRITE(MPG,'(A)') ' ** Warning: Scaling not applied.'
              WRITE(MPG,'(A)') ' ** (incompatibility with null space)'
            END IF
            id%KEEP(52) = 0
          END IF
C         ------------------------
C         If Schur has been asked
C         for, scaling is disabled
C         ------------------------
          IF ( id%KEEP(60) .ne. 0 .and. id%KEEP(52) .ne. 0 ) THEN
            id%KEEP(52) = 0
            IF ( MPG .GT. 0 .AND. id%ICNTL(8) .NE. 0 ) THEN
              WRITE(MPG,'(A)') ' ** Warning: Scaling not applied.'
              WRITE(MPG,'(A)') ' ** (incompatibility with Schur)'
            END IF
          END IF
C         -------------------------------
C         If matrix is distributed on
C         entry, only options 7 and 8
C         of scaling are allowed.
C         -------------------------------
          IF (id%KEEP(54) .NE. 0 .AND. 
     &        id%KEEP(52).NE.7 .AND. id%KEEP(52).NE.8 .AND.
     &        id%KEEP(52) .NE. 0 ) THEN
             id%KEEP(52) = 0
             IF ( MPG .GT. 0 .and. id%ICNTL(8) .ne. 0 ) THEN
               WRITE(MPG,'(A)')
     &         ' ** Warning: This scaling option not available'
               WRITE(MPG,'(A)') ' ** for distributed matrix entry'
             END IF
          END IF
C         ------------------------------------
C         If matrix is symmetric, only scaling
C         options -1 (given scaling), 1
C         (diagonal scaling), and 7 (SIMSCALING)
C         ------------------------------------
          IF ( id%KEEP(50) .NE. 0 ) THEN
             IF ( id%KEEP(52).ne.  1 .and.
     &            id%KEEP(52).ne. -1 .and.
     &            id%KEEP(52).ne.  0 .and.
     &            id%KEEP(52).ne.  7 .and.
     &            id%KEEP(52).ne.  8 .and.
     &            id%KEEP(52).ne. -2 .and.
     &            id%KEEP(52).ne. 77) THEN
              IF ( MPG .GT. 0 ) THEN
                WRITE(MPG,'(A)')
     &  ' ** Warning: Scaling option n.a. for symmetric matrix'
              END IF
              id%KEEP(52) = 0
            END IF
          END IF
C         ----------------------------------
C         If matrix is elemental on entry, 
C         automatic scaling is now forbidden
C         ----------------------------------
          IF (id%KEEP(55) .NE. 0 .AND. 
     &        ( id%KEEP(52) .gt. 0 ) ) THEN
            id%KEEP(52) = 0
            IF ( MPG .GT. 0 ) THEN
              WRITE(MPG,'(A)') ' ** Warning: Scaling not applied.'
              WRITE(MPG,'(A)')
     &        ' ** (only user scaling av. for elt. entry)'
            END IF
          END IF
C         --------------------------------------
C         Check input parameters ROWSCA / COLSCA
C         --------------------------------------
          IF ( id%KEEP(52) .eq. -1 ) THEN
            IF ( .not. associated( id%ROWSCA ) ) THEN
              id%INFO(1) = -22
              id%INFO(2) = 5
            ELSE IF ( size( id%ROWSCA ) < id%N ) THEN
              id%INFO(1) = -22
              id%INFO(2) = 5
            ELSE IF ( .not. associated( id%COLSCA ) ) THEN
              id%INFO(1) = -22
              id%INFO(2) = 6
            ELSE IF ( size( id%COLSCA ) < id%N ) THEN
              id%INFO(1) = -22
              id%INFO(2) = 6
            END IF
          END IF
C
C  Allocate -- if required,
C  ROWSCA and COLSCA on the master
C
C  Allocation of scaling arrays.
C  If ICNTL(8) == -1, ROWSCA and COLSCA must have been associated and
C  filled by the user. If ICNTL(8) is >0 and <= 8, the scaling is
C  computed at the beginning of ZMUMPS_FAC_DRIVER and is allocated now.
C
          IF (id%KEEP(52).GT.0 .AND.
     &        id%KEEP(52) .LE.8) THEN
            IF ( associated(id%COLSCA))
     &             DEALLOCATE( id%COLSCA )
            IF ( associated(id%ROWSCA))
     &             DEALLOCATE( id%ROWSCA )
            ALLOCATE( id%COLSCA(N), stat=IERR)
            IF (IERR .GT.0) id%INFO(1)=-13
            ALLOCATE( id%ROWSCA(N), stat=IERR)
            IF (IERR .GT.0) id%INFO(1)=-13
          END IF
C
C         Allocate scaling arrays of size 1 if
C         they are not used to avoid problems
C         when passing them in arguments
C
          IF (.NOT. associated(id%COLSCA)) THEN
            ALLOCATE( id%COLSCA(1), stat=IERR)
          END IF
          IF (IERR .GT.0) id%INFO(1)=-13
          IF (.NOT. associated(id%ROWSCA))
     &    ALLOCATE( id%ROWSCA(1), stat=IERR)
          IF (IERR .GT.0) id%INFO(1)=-13
          IF ( id%INFO(1) .eq. -13 ) THEN
            IF ( LPOK ) WRITE(LP,'(A)')
     &         'Problems in allocations before facto'
            GOTO 200
          END IF
          IF (id%KEEP(252) .EQ. 1) THEN
             CALL ZMUMPS_CHECK_DENSE_RHS
     &       (id%RHS,id%INFO,id%N,id%NRHS,id%LRHS)
C            Sets KEEP(221) and do some checks
             CALL ZMUMPS_SET_K221(id)
             CALL ZMUMPS_CHECK_REDRHS(id)
          ENDIF
 200      CONTINUE
        END IF        ! End of IF (MYID .eq. MASTER)
C       KEEP(221) was set in ZMUMPS_SET_K221 but not broadcast
        CALL MPI_BCAST( id%KEEP(221), 1, MPI_INTEGER, MASTER, id%COMM,
     &                  IERR )
C
C  Check Schur complement on all processors.
C  ZMUMPS_PROPINFO will be called right after those checks.
C
        IF (id%KEEP(60).EQ.2.OR.id%KEEP(60).EQ.3) THEN
          IF ( id%root%yes ) THEN
            IF ( associated( id%SCHUR_CINTERFACE )) THEN
C             Called from C interface...
C             
              id%SCHUR=>id%SCHUR_CINTERFACE
     &          (1:id%SCHUR_LLD*(id%root%SCHUR_NLOC-1)+
     &          id%root%SCHUR_MLOC)
            ENDIF
C           Check that SCHUR_LLD is large enough
            IF (id%SCHUR_LLD < id%root%SCHUR_MLOC) THEN
              IF (LP.GT.0) write(LP,*) 
     &          ' SCHUR leading dimension SCHUR_LLD ', 
     &          id%SCHUR_LLD, 'too small with respect to', 
     &          id%root%SCHUR_MLOC
              id%INFO(1)=-30
              id%INFO(2)=id%SCHUR_LLD
            ELSE IF ( .NOT. associated (id%SCHUR)) THEN
              IF (LP.GT.0) write(LP,'(A)') 
     &                      ' SCHUR not associated'
              id%INFO(1)=-22
              id%INFO(2)=9
            ELSE IF (size(id%SCHUR) <
     &          id%SCHUR_LLD*(id%root%SCHUR_NLOC-1)+
     &          id%root%SCHUR_MLOC) THEN
              IF (LP.GT.0) THEN 
                write(LP,'(A)') 
     &                      ' SCHUR allocated but too small'
                write(LP,*) id%MYID, ' : Size Schur=', 
     &          size(id%SCHUR), 
     &          ' SCHUR_LLD= ', id%SCHUR_LLD, 
     &          ' SCHUR_MLOC=', id%root%SCHUR_NLOC, 
     &          ' SCHUR_NLOC=', id%root%SCHUR_NLOC
              ENDIF
              id%INFO(1)=-22
              id%INFO(2)= 9
            ELSE
C              We initialize the pointer that
C              we will use within ZMUMPS here.
               id%root%SCHUR_LLD=id%SCHUR_LLD
               IF (id%root%SCHUR_NLOC==0) THEN
                 ALLOCATE(id%root%SCHUR_POINTER(1))
               ELSE
                id%root%SCHUR_POINTER=>id%SCHUR
               ENDIF
            ENDIF
          ENDIF
        ENDIF
C       -------------------------
C       Propagate possible errors
C       -------------------------
        CALL MUMPS_PROPINFO( id%ICNTL(1),
     &                      id%INFO(1),
     &                      id%COMM, id%MYID )
        IF ( id%INFO(1) .LT. 0 ) GO TO 499
C       -----------------------------------------------
C       Call factorization procedure ZMUMPS_FAC_DRIVER
C       -----------------------------------------------
        CALL ZMUMPS_FAC_DRIVER(id)
C       Save scaling in INFOG(33)
        IF (id%MYID .eq. MASTER) id%INFOG(33)=id%KEEP(52)
C
C       In the case of Schur, free or not associated
C       id%root%SCHUR_POINTER now rather than in end_driver.F
C       (Case of repeated factorizations).
        IF (id%KEEP(60).EQ.2.OR.id%KEEP(60).EQ.3) THEN
          IF (id%root%yes) THEN
            IF (id%root%SCHUR_NLOC==0) THEN
               DEALLOCATE(id%root%SCHUR_POINTER)
               NULLIFY(id%root%SCHUR_POINTER)
            ELSE
               NULLIFY(id%root%SCHUR_POINTER)
            ENDIF
          ENDIF
        ENDIF
        IF (id%MYID .eq. MASTER) THEN
           CALL MUMPS_SECFIN(TIMEG)
           id%DKEEP(91) = TIMEG
        ENDIF
        IF (PROKG) THEN
            WRITE( MPG,'(A,F12.4)')
     &         ' ELAPSED TIME IN FACTORIZATION DRIVER= ', TIMEG
        END IF 
C
C       Check for errors after FACTO
C       (it was propagated inside)
        IF ( id%INFO(1) .LT. 0 ) GO TO 499
C
C       Update last successful step
C
        id%KEEP(40) = 2 - 456789
      END IF
      IF (LSOLVE) THEN
        IF (id%MYID .eq. MASTER) THEN
           id%DKEEP(111)=0.0D0
           CALL MUMPS_SECDEB(TIMEG)
        END IF 
C       ---------------------
C       Reset KEEP(40) to 2.
C       (last successful step
C       was facto)
C       ---------------------
        id%KEEP(40) = 2 -456789
C       ------------------------------------------
C       Call solution procedure ZMUMPS_SOLVE_DRIVER
C       ------------------------------------------
        IF (id%MYID .eq. MASTER) THEN
           KEEP235SAVE = id%KEEP(235)
           KEEP242SAVE = id%KEEP(242)
           KEEP243SAVE = id%KEEP(243)
           KEEP495SAVE = id%KEEP(495)
           KEEP497SAVE = id%KEEP(497)
           ! if no permutation of RHS asked then suppress request
           ! to interleave the RHS
           ! to interleave the RHS on ordering given then 
           ! using option to set permutation to identity should be 
           ! used (note though that 
           ! they # with A-1/sparseRHS and Null Space)
           IF (id%KEEP(242).EQ.0) id%KEEP(243)=0
        ENDIF
        CALL ZMUMPS_SOLVE_DRIVER(id)
        IF (id%MYID .eq. MASTER) THEN
            CALL MUMPS_SECFIN(TIMEG)
            id%DKEEP(111) = TIMEG
        ENDIF
        IF (PROKG) THEN
            WRITE( MPG,'(A,F12.4)')
     &         ' ELAPSED TIME IN SOLVE DRIVER= ', TIMEG
        END IF 
        IF (id%MYID .eq. MASTER) THEN
           id%KEEP(235) = KEEP235SAVE
           id%KEEP(242) = KEEP242SAVE
           id%KEEP(243) = KEEP243SAVE
           id%KEEP(495) = KEEP495SAVE
           id%KEEP(497) = KEEP497SAVE
        ENDIF
        IF (id%INFO(1).LT.0) GOTO 499
C       ---------------------------
C       Update last successful step
C       ---------------------------
        id%KEEP(40) = 3 -456789
      ENDIF
C
C  What was actually done is saved in KEEP(40)
C
      IF (PROK) CALL ZMUMPS_PRINT_ICNTL(id, MP)
      GOTO 500
*
*=================
* ERROR section
*=================
  499 CONTINUE
*     Print error message if PROK
      IF (LPOK) WRITE (LP,99995) id%INFO(1)
      IF (LPOK) WRITE (LP,99994) id%INFO(2)
*
500   CONTINUE
#if ! defined(LARGEMATRICES)
C     ---------------------------------
C     Permute JCN on output to ZMUMPS if
C     KEEP(23) is different from 0.
C     ---------------------------------
      IF (id%MYID .eq. MASTER .AND. id%KEEP(23) .NE. 0
     &    .AND. NOERRORBEFOREPERM) THEN
C       -------------------------------
C       IF JOB=3 and PERM was not
C       done (no iterative refinement/
C       error analysis), then we do not
C       permute JCN back.
C       -------------------------------
        IF (id%JOB .NE. 3 .OR. UNS_PERM_DONE) THEN
          DO I = 1, id%NZ
            J=id%JCN(I)
C           -- skip out-of range (that are ignored in ANA_O)
            IF (J.LE.0.OR.J.GT.id%N) CYCLE
            id%JCN(I)=id%UNS_PERM(J)
          END DO
        END IF
      END IF
#endif
 510  CONTINUE
C     ---------------------------------
C     Modify INFOG(1:2) for concordance
C     with specification sheets 3.2, ie
C     return the same significant error
C     code on all processors.
C     And broadcast other INFOG entries
C     ---------------------------------
      CALL ZMUMPS_SET_INFOG(id%INFO(1), id%INFOG(1), id%COMM, id%MYID)
C
C     --------------------------------
C     Broadcast RINFOG entries to make
C     them available on all procs.
C     --------------------------------
      CALL MPI_BCAST( id%RINFOG(1), 40, MPI_DOUBLE_PRECISION, MASTER,
     &                    id%COMM, IERR )
      IF (id%INFOG(1).GE.0 .AND. JOB.NE.-1  
     &     .AND. JOB.NE.-2 ) THEN
         IF (id%MYID .eq. MASTER) THEN
            CALL MUMPS_SECFIN(TIMETOTAL)
            id%DKEEP(70) = TIMEG
         ENDIF
      ENDIF
*===============
* ERRORG section
*===============
      IF (id%MYID.EQ.MASTER.and.MPG.GT.0.and.
     & id%INFOG(1).lt.0) THEN
        WRITE(MPG,'(A,I12)') ' On return from ZMUMPS, INFOG(1)=',
     &      id%INFOG(1)
        WRITE(MPG,'(A,I12)') ' On return from ZMUMPS, INFOG(2)=',
     &      id%INFOG(2)
      END IF
C     -------------------------
C     Restore user communicator
C     -------------------------
       CALL MPI_COMM_FREE( id%COMM, IERR )
       id%COMM = COMM_SAVE
      RETURN
*
99995 FORMAT (' ** ERROR RETURN ** FROM ZMUMPS INFO(1)=', I3)
99994 FORMAT (' ** INFO(2)=', I10)
99993 FORMAT (' ** Allocation error: could not permute JCN.')
      END SUBROUTINE ZMUMPS
*
      SUBROUTINE ZMUMPS_SET_INFOG( INFO, INFOG, COMM, MYID )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
C
C  Purpose:
C  =======
C
C  If one proc has INFO(1).lt.0 and INFO(1) .ne. -1,
C  puts INFO(1:2) of this proc on all procs in INFOG
C
C  Arguments:
C  =========
C
      INTEGER INFO(40), INFOG(40), COMM, MYID
C
C  Local variables
C  ===============
C
      INTEGER TMP1(2),TMP(2)
      INTEGER ROOT, IERR
      INTEGER MASTER
      PARAMETER (MASTER=0)
C
C
      IF ( INFO(1) .ge. 0  .and. INFO(2) .ge. 0 ) THEN
C
C       This can only happen if the phase was successful
C       on all procs. If one proc failed, then all other
C       procs would have INFO(1)=-1.
C
        INFOG(1) = INFO(1)
        INFOG(2) = INFO(2)
      ELSE
C       ---------------------
C       Find who has smallest
C       error code INFO(1)
C       ---------------------
        INFOG(1) = INFO(1)
C        INFOG(2) = MYID
        TMP1(1) = INFO(1)
        TMP1(2) = MYID
        CALL MPI_ALLREDUCE(TMP1,TMP,1,MPI_2INTEGER,
     &                     MPI_MINLOC,COMM,IERR )
        INFOG(2) = INFO(2)
        ROOT = TMP(2)
        CALL MPI_BCAST( INFOG(1), 1, MPI_INTEGER, ROOT, COMM, IERR )
        CALL MPI_BCAST( INFOG(2), 1, MPI_INTEGER, ROOT, COMM, IERR )
      END IF
C
C    Make INFOG available on all procs:
C
      CALL MPI_BCAST(INFOG(3), 38, MPI_INTEGER, MASTER, COMM, IERR )
      RETURN
      END SUBROUTINE ZMUMPS_SET_INFOG
      SUBROUTINE ZMUMPS_PRINT_ICNTL(id, LP)
      USE ZMUMPS_STRUC_DEF
*
*  ==========
*  Parameters
*  ==========
      TYPE (ZMUMPS_STRUC), TARGET, INTENT(IN) :: id
      INTEGER  :: LP
** Local Variables
      INTEGER, POINTER :: JOB 
      INTEGER,DIMENSION(:),POINTER::ICNTL
      INTEGER MASTER
      PARAMETER( MASTER = 0 )
      IF (LP.LE.0) RETURN
      JOB=>id%JOB
      ICNTL=>id%ICNTL
      IF (id%MYID.EQ.MASTER) THEN
         SELECT CASE (JOB)
         CASE(1);
           WRITE (LP,980) 
           WRITE (LP,990) ICNTL(1),ICNTL(2),ICNTL(3),ICNTL(4)
           WRITE (LP,991) ICNTL(5),ICNTL(6),ICNTL(7),ICNTL(12),
     &          ICNTL(13),ICNTL(18),ICNTL(19),ICNTL(22)
           IF ((ICNTL(6).EQ.5).OR.(ICNTL(6).EQ.6).OR.
     &          (ICNTL(12).NE.1) )  THEN
              WRITE (LP,992) ICNTL(8)
           ENDIF   
           IF (id%ICNTL(19).NE.0)
     &      WRITE(LP,998) id%SIZE_SCHUR
           WRITE (LP,993) ICNTL(14)
         CASE(2);
           WRITE (LP,980) 
           WRITE (LP,990) ICNTL(1),ICNTL(2),ICNTL(3),ICNTL(4)
           WRITE (LP,992) ICNTL(8)
           WRITE (LP,993) ICNTL(14)
         CASE(3);
           WRITE (LP,980)
           WRITE (LP,990) ICNTL(1),ICNTL(2),ICNTL(3),ICNTL(4)
           WRITE (LP,995)
     &     ICNTL(9),ICNTL(10),ICNTL(11),ICNTL(20),ICNTL(21)
         CASE(4);
           WRITE (LP,980) 
           WRITE (LP,990) ICNTL(1),ICNTL(2),ICNTL(3),ICNTL(4)
           WRITE (LP,992) ICNTL(8)
           IF (id%ICNTL(19).NE.0)
     &      WRITE(LP,998) id%SIZE_SCHUR
           WRITE (LP,993) ICNTL(14)
         CASE(5);
           WRITE (LP,980) 
           WRITE (LP,990) ICNTL(1),ICNTL(2),ICNTL(3),ICNTL(4)
           WRITE (LP,991) ICNTL(5),ICNTL(6),ICNTL(7),ICNTL(12),
     &          ICNTL(13),ICNTL(18),ICNTL(19),ICNTL(22)
           WRITE (LP,992) ICNTL(8)
           WRITE (LP,993) ICNTL(14)
           WRITE (LP,995)
     &     ICNTL(9),ICNTL(10),ICNTL(11),ICNTL(20),ICNTL(21)
         CASE(6);
           WRITE (LP,980)
           WRITE (LP,990) ICNTL(1),ICNTL(2),ICNTL(3),ICNTL(4)
           WRITE (LP,991) ICNTL(5),ICNTL(6),ICNTL(7),ICNTL(12),
     &          ICNTL(13),ICNTL(18),ICNTL(19),ICNTL(22)
           IF (id%ICNTL(19).NE.0)
     &      WRITE(LP,998) id%SIZE_SCHUR
           WRITE (LP,992) ICNTL(8)
           WRITE (LP,995)
     &     ICNTL(9),ICNTL(10),ICNTL(11),ICNTL(20),ICNTL(21)
           WRITE (LP,993) ICNTL(14)
        END SELECT
      ENDIF
 980  FORMAT (/'***********CONTROL PARAMETERS (ICNTL)**************'/)
 990  FORMAT (
     &     'ICNTL(1)   Output stream for error messages        =',I10/
     &     'ICNTL(2)   Output stream for diagnostic messages   =',I10/
     &     'ICNTL(3)   Output stream for global information    =',I10/
     &     'ICNTL(4)   Level of printing                       =',I10)
 991  FORMAT (
     &     'ICNTL(5)   Matrix format  ( keep(55) )             =',I10/
     &     'ICNTL(6)   Maximum transversal  ( keep(23) )       =',I10/
     &     'ICNTL(7)   Ordering                                =',I10/
     &     'ICNTL(12)  LDLT ordering strat ( keep(95) )        =',I10/
     &     'ICNTL(13)  Parallel root (0=on, 1=off)             =',I10/
     &     'ICNTL(18)  Distributed matrix  ( keep(54) )        =',I10/
     &     'ICNTL(19)  Schur option ( keep(60) 0=off,else=on ) =',I10/
     &     'ICNTL(22)  Out-off-core option (0=Off, >0=ON)      =',I10)
 992  FORMAT (
     &     'ICNTL(8)   Scaling strategy                        =',I10)
 993  FORMAT (
     &     'ICNTL(14)  Percent of memory increase              =',I10)
 995  FORMAT (
     &     'ICNTL(9)   Solve A x=b (1) or A''x = b (else)       =',I10/
     &     'ICNTL(10)  Max steps iterative refinement          =',I10/
     &     'ICNTL(11)  Error analysis ( 0= off, else=on)       =',I10/
     &     'ICNTL(20)  Dense (0) or sparse (1) RHS             =',I10/
     &     'ICNTL(21)  Gathered (0) or distributed(1) solution =',I10)
 998  FORMAT (
     &     '      Size of SCHUR matrix (SIZE_SHUR)             =',I10)
      END SUBROUTINE ZMUMPS_PRINT_ICNTL
C--------------------------------------------------------------------
      SUBROUTINE ZMUMPS_PRINT_KEEP(id, LP)
      USE ZMUMPS_STRUC_DEF
*
*  ==========
*  Parameters
*  ==========
      TYPE (ZMUMPS_STRUC), TARGET, INTENT(IN) :: id
      INTEGER ::LP
** Local Variables
      INTEGER, POINTER :: JOB 
      INTEGER,DIMENSION(:),POINTER::ICNTL, KEEP
      INTEGER MASTER
      PARAMETER( MASTER = 0 )
      IF (LP.LE.0) RETURN
      JOB=>id%JOB
      ICNTL=>id%ICNTL
      KEEP=>id%KEEP
      IF (id%MYID.EQ.MASTER) THEN
         SELECT CASE (JOB)
         CASE(1);
           WRITE (LP,980) 
           WRITE (LP,990) ICNTL(1),ICNTL(2),ICNTL(3),ICNTL(4)
           WRITE (LP,991) KEEP(55),KEEP(23),ICNTL(7),KEEP(95),
     &          ICNTL(13),KEEP(54),KEEP(60),ICNTL(22)
           IF ((KEEP(23).EQ.5).OR.(KEEP(23).EQ.6))THEN
              WRITE (LP,992) KEEP(52)
           ENDIF   
           WRITE (LP,993) KEEP(12)
         CASE(2);
           WRITE (LP,980)
           WRITE (LP,990) ICNTL(1),ICNTL(2),ICNTL(3),ICNTL(4)
           IF (KEEP(23).EQ.0)THEN
              WRITE (LP,992) KEEP(52)
           ENDIF   
           WRITE (LP,993) KEEP(12)
         CASE(3);
           WRITE (LP,980)
           WRITE (LP,990) ICNTL(1),ICNTL(2),ICNTL(3),ICNTL(4) 
           WRITE (LP,995)
     &     ICNTL(9),ICNTL(10),ICNTL(11),ICNTL(20),ICNTL(21)
         CASE(4);
           WRITE (LP,980) 
           WRITE (LP,990) ICNTL(1),ICNTL(2),ICNTL(3),ICNTL(4)
           IF (KEEP(23).NE.0)THEN
              WRITE (LP,992) KEEP(52)
           ENDIF  
           WRITE (LP,991) KEEP(55),KEEP(23),ICNTL(7),KEEP(95),
     &          ICNTL(13),KEEP(54),KEEP(60),ICNTL(22)
           WRITE (LP,995)
     &     ICNTL(9),ICNTL(10),ICNTL(11),ICNTL(20),ICNTL(21)
           WRITE (LP,993) KEEP(12)
         CASE(5);
           WRITE (LP,980) 
           WRITE (LP,990) ICNTL(1),ICNTL(2),ICNTL(3),ICNTL(4)
           WRITE (LP,991) KEEP(55),KEEP(23),ICNTL(7),KEEP(95),
     &          ICNTL(13),KEEP(54),KEEP(60),ICNTL(22)
           IF ((KEEP(23).EQ.5).OR.(KEEP(23).EQ.6)
     &       .OR. (KEEP(23).EQ.7)) THEN
              WRITE (LP,992) KEEP(52)
           ENDIF              
           IF (KEEP(23).EQ.0)THEN
              WRITE (LP,992) KEEP(52)
           ENDIF   
           WRITE (LP,993) KEEP(12)
         CASE(6);
           WRITE (LP,980)
           WRITE (LP,990) ICNTL(1),ICNTL(2),ICNTL(3),ICNTL(4)
           WRITE (LP,991) KEEP(55),KEEP(23),ICNTL(7),KEEP(95),
     &          ICNTL(13),KEEP(54),KEEP(60),ICNTL(22)
           IF ((KEEP(23).EQ.5).OR.(KEEP(23).EQ.6)
     &       .OR. (KEEP(23).EQ.7)) THEN
              WRITE (LP,992) KEEP(52)
           ENDIF   
           IF (KEEP(23).EQ.0)THEN
              WRITE (LP,992) KEEP(52)
           ENDIF   
           WRITE (LP,995)
     &     ICNTL(9),ICNTL(10),ICNTL(11),KEEP(248),ICNTL(21)
           WRITE (LP,993) KEEP(12)
        END SELECT
      ENDIF
 980  FORMAT (/'******INTERNAL VALUE OF PARAMETERS (ICNTL/KEEP)****'/)
 990  FORMAT (
     &     'ICNTL(1)   Output stream for error messages        =',I10/
     &     'ICNTL(2)   Output stream for diagnostic messages   =',I10/
     &     'ICNTL(3)   Output stream for global information    =',I10/
     &     'ICNTL(4)   Level of printing                       =',I10)
 991  FORMAT (
     &     'ICNTL(5)   Matrix format  ( keep(55) )             =',I10/
     &     'ICNTL(6)   Maximum transversal  ( keep(23) )       =',I10/
     &     'ICNTL(7)   Ordering                                =',I10/
     &     'ICNTL(12)  LDLT ordering strat ( keep(95) )        =',I10/
     &     'ICNTL(13)  Parallel root (0=on, 1=off)             =',I10/
     &     'ICNTL(18)  Distributed matrix  ( keep(54) )        =',I10/
     &     'ICNTL(19)  Schur option ( keep(60) 0=off,else=on ) =',I10/
     &     'ICNTL(22)  Out-off-core option (0=Off, >0=ON)      =',I10)
 992  FORMAT (
     &     'ICNTL(8)   Scaling strategy ( keep(52) )           =',I10)
 993  FORMAT (
     &     'ICNTL(14)  Percent of memory increase ( keep(12) ) =',I10)
 995  FORMAT (
     &     'ICNTL(9)   Solve A x=b (1) or A''x = b (else)       =',I10/
     &     'ICNTL(10)  Max steps iterative refinement          =',I10/
     &     'ICNTL(11)  Error analysis ( 0= off, else=on)       =',I10/
     &     'ICNTL(20)  Dense (0) or sparse (1) RHS             =',I10/
     &     'ICNTL(21)  Gathered (0) or distributed(1) solution =',I10)
      END SUBROUTINE ZMUMPS_PRINT_KEEP
      SUBROUTINE ZMUMPS_CHECK_DENSE_RHS
     &       (idRHS, idINFO, idN, idNRHS, idLRHS)
      IMPLICIT NONE
C
C  Purpose:
C  =======
C
C     Check that the dense RHS is associated and of
C     correct size. Called on master only, when dense
C     RHS is supposed to be allocated. This can be used
C     either at the beginning of the solve phase or
C     at the beginning of the factorization phase
C     if forward solve is done during factorization
C     (see ICNTL(32)) ; idINFO(1), idINFO(2) may be
C     modified.
C
C
C  Arguments:
C  =========
C
C     id* : see corresponding components of the main
C     MUMPS structure.
C
      COMPLEX(kind=8), DIMENSION(:), POINTER :: idRHS
      INTEGER, intent(in)    :: idN, idNRHS, idLRHS
      INTEGER, intent(inout) :: idINFO(:)
      IF ( .not. associated( idRHS ) ) THEN
              idINFO( 1 ) = -22
              idINFO( 2 ) = 7
      ELSE IF (idNRHS.EQ.1) THEN
               IF ( size( idRHS ) < idN ) THEN
                  idINFO( 1 ) = -22
                  idINFO( 2 ) = 7
               ENDIF
      ELSE IF (idLRHS < idN) 
     &            THEN
                  idINFO( 1 ) = -26
                  idINFO( 2 ) = idLRHS
      ELSE IF 
     &      (size(idRHS)<(idNRHS*idLRHS-idLRHS+idN)) 
     &            THEN
                  idINFO( 1 ) = -22
                  idINFO( 2 ) = 7
      END IF
      RETURN
      END SUBROUTINE ZMUMPS_CHECK_DENSE_RHS
C
      SUBROUTINE ZMUMPS_SET_K221(id)
      USE ZMUMPS_STRUC_DEF
      IMPLICIT NONE
C
C     Purpose:
C     =======
C
C     Sets KEEP(221) on master.
C     Constraint: must be called before ZMUMPS_CHECK_REDRHS.
C     Can be called at factorization or solve phase
C
      TYPE (ZMUMPS_STRUC) :: id
      INTEGER MASTER
      PARAMETER( MASTER = 0 )
      IF (id%MYID.EQ.MASTER) THEN
        id%KEEP(221)=id%ICNTL(26)
        IF (id%KEEP(221).ne.0 .and. id%KEEP(221) .NE.1
     &      .AND.id%KEEP(221).ne.2) id%KEEP(221)=0
      ENDIF
      RETURN
      END SUBROUTINE ZMUMPS_SET_K221
C
      SUBROUTINE ZMUMPS_CHECK_REDRHS(id)
      USE ZMUMPS_STRUC_DEF
      IMPLICIT NONE
C
C  Purpose:
C  =======
C
C  * Decode API related to REDRHS and check REDRHS
C  * Can be called at factorization or solve phase
C  * Constraints:
C    - Must be called after solve phase.
C    - KEEP(60) must have been set (ok to check
C    since KEEP(60) was set during analysis phase)
C  * Remark that during solve phase, ICNTL(26)=1 is
C    forbidden in case of fwd in facto.
C
      TYPE (ZMUMPS_STRUC) :: id
      INTEGER MASTER
      PARAMETER( MASTER = 0 )
      IF (id%MYID .EQ. MASTER) THEN
          IF ( id%KEEP(221) == 1 .or. id%KEEP(221) == 2 ) THEN
            IF (id%KEEP(221) == 2 .and. id%JOB == 2) THEN
              id%INFO(1)=-35
              id%INFO(2)=id%KEEP(221)
              GOTO 333
            ENDIF
            IF (id%KEEP(221) == 1 .and. id%KEEP(252) == 1
     &          .and. id%JOB == 3) THEN
              id%INFO(1)=-35
              id%INFO(2)=id%KEEP(221)
            ENDIF
            IF ( id%KEEP(60).eq. 0 .or. id%SIZE_SCHUR.EQ.0 ) THEN
              id%INFO(1)=-33
              id%INFO(2)=id%KEEP(221)
              GOTO 333
            ENDIF
            IF ( .NOT. associated( id%REDRHS)) THEN
              id%INFO(1)=-22
              id%INFO(2)=15
              GOTO 333
            ELSE IF (id%NRHS.EQ.1) THEN
              IF (size(id%REDRHS) < id%SIZE_SCHUR ) THEN
                id%INFO(1)=-22
                id%INFO(2)=15
                GOTO 333
              ENDIF
            ELSE IF (id%LREDRHS < id%SIZE_SCHUR) THEN
              id%INFO(1)=-34
              id%INFO(2)=id%LREDRHS
              GOTO 333
            ELSE IF
     &      (size(id%REDRHS)<
     &         id%NRHS*id%LREDRHS-id%LREDRHS+id%SIZE_SCHUR)
     &      THEN
              id%INFO(1)=-22
              id%INFO(2)=15
              GOTO 333
            ENDIF
          ENDIF
      ENDIF
 333  CONTINUE
C     Error is not propagated. It should be propagated outside.
C     The reason to propagate it outside is that there can be
C     one call to PROPINFO instead of several ones.
      RETURN
      END SUBROUTINE ZMUMPS_CHECK_REDRHS