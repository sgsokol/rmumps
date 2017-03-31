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
      MODULE DMUMPS_ANA_LR
      USE DMUMPS_LR_CORE
      USE DMUMPS_LR_STATS
      USE MUMPS_LR_COMMON
      USE MUMPS_ANA_ORD_WRAPPERS
      IMPLICIT NONE
      CONTAINS
      SUBROUTINE GET_CUT(IWR, NASS, NCB, LRGROUPS, NPARTSCB,
     &   NPARTSASS, CUT)
      INTEGER, INTENT(IN) :: NASS, NCB
      INTEGER, INTENT(IN) :: IWR(*)
      INTEGER, INTENT(IN), DIMENSION(:) :: LRGROUPS
      INTEGER, INTENT(OUT) :: NPARTSCB, NPARTSASS
      INTEGER, POINTER, DIMENSION(:) :: CUT
      INTEGER :: I, CURRENT_PART, CUTBUILDER
      INTEGER, ALLOCATABLE, DIMENSION(:) :: BIG_CUT
      ALLOCATE(BIG_CUT(max(NASS,1)+NCB+1))
      CURRENT_PART = LRGROUPS(IWR(1))
      BIG_CUT(1) = 1
      BIG_CUT(2) = 2
      CUTBUILDER = 2
      NPARTSASS = 0
      NPARTSCB  = 0 
      DO I = 2,NASS + NCB
        IF (LRGROUPS(IWR(I)) == CURRENT_PART) THEN
          BIG_CUT(CUTBUILDER) = BIG_CUT(CUTBUILDER) + 1
        ELSE
          CUTBUILDER = CUTBUILDER + 1
          BIG_CUT(CUTBUILDER) = BIG_CUT(CUTBUILDER-1) + 1
          CURRENT_PART = LRGROUPS(IWR(I))
        END IF
        IF (I == NASS) NPARTSASS = CUTBUILDER - 1
      END DO
      IF (NASS.EQ.1) NPARTSASS= 1
      NPARTSCB = CUTBUILDER - 1 - NPARTSASS
      ALLOCATE(CUT(max(NPARTSASS,1)+NPARTSCB+1))
      IF (NPARTSASS.EQ.0) THEN
        CUT(1) = 1
        CUT(2:2+NPARTSCB) = BIG_CUT(1:1+NPARTSCB)
      ELSE
        CUT = BIG_CUT(1:NPARTSASS+NPARTSCB+1)
      ENDIF 
      if(allocated(BIG_CUT)) DEALLOCATE(BIG_CUT)
      END SUBROUTINE
      SUBROUTINE SEP_GROUPING(NV, VLIST, N, NZ, LRGROUPS, NBGROUPS, IW, 
     &   LW, IPE, LEN, GROUP_SIZE, HALO_DEPTH, TRACE, WORKH, NODE,
     &   GEN2HALO, K482, K472, K469, SEP_SIZE, 
     &   KEEP10, LP, LPOK, IFLAG, IERROR)
        INTEGER(8), INTENT(IN) :: NZ, LW
        INTEGER, INTENT(IN)    :: NV, N, GROUP_SIZE, HALO_DEPTH
        INTEGER, INTENT(IN)    :: IW(LW), LEN(N), NODE, K482
        INTEGER(8), INTENT(IN) :: IPE(N+1)
        INTEGER, INTENT(IN)    :: K472, K469, SEP_SIZE, KEEP10, LP
        LOGICAL                :: LPOK
        INTEGER, INTENT(INOUT) :: NBGROUPS, WORKH(N)
        INTEGER, INTENT(INOUT) :: LRGROUPS(N), VLIST(NV), TRACE(N)
        INTEGER, INTENT(INOUT) :: GEN2HALO(N)
        INTEGER, INTENT(INOUT) :: IFLAG, IERROR
        INTEGER(8), ALLOCATABLE, DIMENSION(:) :: IPTRHALO
        INTEGER,    ALLOCATABLE, DIMENSION(:) :: PARTS, JCNHALO
        INTEGER(8) :: HALOEDGENBR
        INTEGER    :: NHALO,
     &       NBGROUPS_KWAY, I, GROUP_SIZE2, LRGROUPS_SIGN, IERR
#if defined (metis) || defined (parmetis) || defined (metis4) || defined (parmetis3)
        INTEGER :: METIS_IDX_SIZE
#endif
#if defined (scotch) || defined (ptscotch)
        INTEGER :: SCOTCH_IDX_SIZE
#endif
        CALL COMPUTE_BLR_VCS(K472, GROUP_SIZE2, GROUP_SIZE, NV)
        NBGROUPS_KWAY = MAX(NINT(dble(NV)/dble(GROUP_SIZE2)),1)
        IF (NV .GE. SEP_SIZE) THEN
          LRGROUPS_SIGN = 1
        ELSE
          LRGROUPS_SIGN = -1
        ENDIF
        IF (NBGROUPS_KWAY > 1) THEN 
          IF (K469.EQ.3) THEN
!$OMP CRITICAL(gethalo_cri)
            CALL GETHALONODES(N, IW, LW, IPE, VLIST, NV, HALO_DEPTH,
     &                 NHALO, TRACE, WORKH, NODE, LEN, HALOEDGENBR,
     &                 GEN2HALO)
            ALLOCATE(PARTS(NHALO), IPTRHALO(NHALO+1),
     &      JCNHALO(HALOEDGENBR), STAT=IERR)
            IF (IERR.GT.0) THEN
              IF (LPOK) WRITE(LP,*) 
     &        " Error allocate integer array of size: ", 
     &         int(NHALO+(KEEP10*(NHALO+1)),8) + HALOEDGENBR
              IFLAG  = -7
              CALL MUMPS_SET_IERROR 
     &         (int(NHALO+(KEEP10*(NHALO+1)),8) + HALOEDGENBR,
     &          IERROR)
            ENDIF
            CALL GETHALOGRAPH(WORKH, NHALO, N, IW, LW, IPE, IPTRHALO,
     &     JCNHALO, HALOEDGENBR,TRACE,NODE, GEN2HALO)
!$OMP END CRITICAL(gethalo_cri)
            IF (IFLAG.LT.0) RETURN
          ELSE
            CALL GETHALONODES(N, IW, LW, IPE, VLIST, NV, HALO_DEPTH,
     &                 NHALO, TRACE, WORKH, NODE, LEN, HALOEDGENBR,
     &                 GEN2HALO)
            ALLOCATE(PARTS(NHALO), IPTRHALO(NHALO+1),
     &      JCNHALO(HALOEDGENBR), STAT=IERR)
            IF (IERR.GT.0) THEN
              IF (LPOK) WRITE(LP,*) 
     &        " Error allocate integer array of size: ", 
     &         int(NHALO+(KEEP10*(NHALO+1)),8) + HALOEDGENBR
              IFLAG  = -7
              CALL MUMPS_SET_IERROR 
     &         (int(NHALO+(KEEP10*(NHALO+1)),8) + HALOEDGENBR,
     &          IERROR)
              RETURN
            ENDIF
            CALL GETHALOGRAPH(WORKH, NHALO, N, IW, LW, IPE, IPTRHALO, 
     &     JCNHALO, HALOEDGENBR,TRACE,NODE, GEN2HALO)
          ENDIF
          IF (K482.EQ.1) THEN
#if defined (metis) || defined (parmetis) || defined (metis4) || defined (parmetis3)
             CALL MUMPS_METIS_IDXSIZE(METIS_IDX_SIZE)
             IF (METIS_IDX_SIZE .EQ. 64) THEN
                CALL MUMPS_METIS_KWAY_MIXEDto64(NHALO, HALOEDGENBR,
     &               IPTRHALO,
     &               JCNHALO,
     &               NBGROUPS_KWAY, PARTS, LP, LPOK, KEEP10,
     &               IFLAG, IERROR)
             ELSE
               IF (KEEP10.EQ.1) THEN
                IFLAG  = -52
                IERROR = 1
               ELSE
                CALL MUMPS_METIS_KWAY_MIXEDto32(NHALO, HALOEDGENBR,
     &               IPTRHALO,
     &               JCNHALO,
     &               NBGROUPS_KWAY, PARTS, LP, LPOK, KEEP10,
     &               IFLAG, IERROR)
               ENDIF
             ENDIF
#endif
          ELSE IF (K482.EQ.2) THEN
#if defined (scotch) || defined (ptscotch)
             CALL MUMPS_SCOTCH_INTSIZE(SCOTCH_IDX_SIZE)
             IF (SCOTCH_IDX_SIZE .EQ. 32) THEN
               IF (KEEP10.EQ.1) THEN
                IFLAG  = -52
                IERROR = 2
               ELSE
                CALL MUMPS_SCOTCH_KWAY_MIXEDto32(
     &               NHALO, HALOEDGENBR, IPTRHALO, JCNHALO,
     &               NBGROUPS_KWAY, PARTS, LP, LPOK, KEEP10,
     &               IFLAG, IERROR)
               ENDIF
             ELSE
                CALL MUMPS_SCOTCH_KWAY_MIXEDto64(
     &               NHALO, HALOEDGENBR, IPTRHALO, JCNHALO,
     &               NBGROUPS_KWAY, PARTS, LP, LPOK, KEEP10,
     &               IFLAG, IERROR)
             END IF
#endif
          ELSE
            WRITE(6,*) " Internal ERROR K482=", K482
            CALL MUMPS_ABORT()
          END IF
          IF (IFLAG.LT.0) GOTO 500
          CALL GET_GLOBAL_GROUPS(PARTS,VLIST,NV,
     &     NBGROUPS_KWAY, LRGROUPS, N, NBGROUPS, LRGROUPS_SIGN)
        ELSE 
!$OMP CRITICAL(lrgrouping_cri)
          DO I=1,NV
            LRGROUPS(VLIST(I)) = LRGROUPS_SIGN*(NBGROUPS + 1)
          END DO
            NBGROUPS = NBGROUPS + 1
!$OMP END CRITICAL(lrgrouping_cri)
        END IF
  500   IF (allocated(IPTRHALO)) then
           DEALLOCATE(IPTRHALO)
        ENDIF
        IF (allocated(PARTS)) then
           DEALLOCATE(PARTS) 
        ENDIF
        IF (allocated(JCNHALO)) then
           DEALLOCATE(JCNHALO   )
        ENDIF
        RETURN
      END SUBROUTINE SEP_GROUPING
      SUBROUTINE GET_GLOBAL_GROUPS(PARTS, SEP, NSEP, NPARTS, 
     &                    LRGROUPS, N, NBGROUPS, LRGROUPS_SIGN)
      INTEGER,INTENT(IN) :: NSEP, N, LRGROUPS_SIGN
      INTEGER :: PARTS(:)
      INTEGER,DIMENSION(:),INTENT(INOUT) :: SEP 
      INTEGER, INTENT(INOUT) :: NPARTS 
      INTEGER, INTENT(INOUT) :: NBGROUPS
      INTEGER, INTENT(INOUT) :: LRGROUPS(N)
      INTEGER::I,CNT,NB_PARTS_WITHOUT_SEP_NODE
      INTEGER,DIMENSION(:),ALLOCATABLE::SIZES, RIGHTPART
      INTEGER,DIMENSION(:),ALLOCATABLE::PARTPTR
      INTEGER,DIMENSION(:),ALLOCATABLE :: NEWSEP
      ALLOCATE( NEWSEP(NSEP),
     &          SIZES(NPARTS),
     &          RIGHTPART(NPARTS),
     &          PARTPTR(NPARTS+1))
      NB_PARTS_WITHOUT_SEP_NODE = 0
      RIGHTPART = 0
      SIZES = 0
      DO I=1,NSEP
        SIZES(PARTS(I)) = SIZES(PARTS(I)) + 1
      END DO
      PARTPTR(1)=1
      CNT = 0
      DO I=2,NPARTS+1
        PARTPTR(I) = PARTPTR(I-1) + SIZES(I-1)
        IF (SIZES(I-1)==0) THEN
          NB_PARTS_WITHOUT_SEP_NODE = NB_PARTS_WITHOUT_SEP_NODE + 1
        ELSE
          CNT = CNT + 1
          RIGHTPART(I-1) = CNT
        END IF
      END DO
      NPARTS = NPARTS - NB_PARTS_WITHOUT_SEP_NODE 
!$OMP CRITICAL(lrgrouping_cri)
      DO I=1,NSEP
        NEWSEP(PARTPTR(PARTS(I))) = SEP(I)  
        LRGROUPS(SEP(I)) = LRGROUPS_SIGN*(RIGHTPART(PARTS(I)) 
     &    + NBGROUPS)
        PARTPTR(PARTS(I)) = 
     &    PARTPTR(PARTS(I)) + 1
      END DO
      NBGROUPS = NBGROUPS + NPARTS
!$OMP END CRITICAL(lrgrouping_cri)
      SEP = NEWSEP
      DEALLOCATE(NEWSEP,SIZES,RIGHTPART,PARTPTR)
      END SUBROUTINE GET_GLOBAL_GROUPS
      SUBROUTINE GETHALONODES(N, IW, LW, IPE, IND, NIND, PMAX, 
     &                        NHALO, TRACE, WORKH, NODE, LEN, CNT,
     &                        GEN2HALO)
        INTEGER,DIMENSION(:),INTENT(IN) :: IND
        INTEGER(8), INTENT(IN)          :: LW
        INTEGER, INTENT(IN)             :: N, NODE
        INTEGER, INTENT(IN)             :: IW(LW), LEN(N)
        INTEGER(8), INTENT(IN)          :: IPE(N+1)
        INTEGER, INTENT(IN)             :: PMAX,NIND
        INTEGER, INTENT(OUT)            :: NHALO
        INTEGER, INTENT(INOUT)          :: TRACE(N), WORKH(N)
        INTEGER                         :: GEN2HALO(N)
        INTEGER(8), INTENT(OUT)         :: CNT
        INTEGER                         :: DEPTH, I, LAST_LVL_START
        INTEGER                         :: HALOI
        INTEGER(8)                      :: J
        WORKH(1:NIND) = IND
        LAST_LVL_START = 1
        NHALO = NIND
        CNT = 0
        DO I=1,NIND
          HALOI = WORKH(I)
          GEN2HALO(HALOI) = I
          IF (TRACE(HALOI) .NE. NODE) THEN
             TRACE(HALOI) = NODE
          END IF
          DO J=IPE(HALOI),IPE(HALOI+1)-1
            IF (TRACE(IW(J)).EQ.NODE) THEN
              CNT = CNT + 2
            END IF
          END DO
        END DO
        DO DEPTH=1,PMAX
          CALL NEIGHBORHOOD(WORKH, NHALO, N, IW, LW, IPE,
     &                      TRACE, NODE, LEN, CNT, LAST_LVL_START,
     &                      DEPTH, PMAX, GEN2HALO)
        END DO
      END SUBROUTINE
      SUBROUTINE NEIGHBORHOOD(HALO, NHALO, N, IW, LW, IPE,
     &                        TRACE, NODE, LEN, CNT, LAST_LVL_START,
     &                        DEPTH, PMAX, GEN2HALO)
        INTEGER, INTENT(IN)                 :: N, NODE, DEPTH, PMAX
        INTEGER,INTENT(INOUT)               :: NHALO, GEN2HALO(N)
        INTEGER, INTENT(INOUT)              :: LAST_LVL_START
        INTEGER(8), INTENT(INOUT)           :: CNT
        INTEGER,DIMENSION(:),INTENT(INOUT)  :: HALO
        INTEGER(8), INTENT(IN)              :: LW
        INTEGER(8), INTENT(IN)              :: IPE(N+1)
        INTEGER, TARGET, INTENT(IN)         :: IW(LW)
        INTEGER, INTENT(IN)                 :: LEN(N)
        INTEGER,DIMENSION(:)                :: TRACE
        INTEGER :: I,INEI,NADJI,NEWNHALO, NEIGH
        INTEGER, DIMENSION(:), POINTER :: ADJI
        INTEGER(8) :: J
        NEWNHALO = 0
        DO I=LAST_LVL_START,NHALO
          NADJI = LEN(HALO(I))
          ADJI => IW(IPE(HALO(I)):IPE(HALO(I)+1)-1)
          DO INEI=1,NADJI
            IF (TRACE(ADJI(INEI)) .NE. NODE) THEN
              NEIGH = ADJI(INEI)
              TRACE(NEIGH) = NODE
              NEWNHALO = NEWNHALO + 1
              HALO(NHALO+NEWNHALO) = NEIGH
              GEN2HALO(NEIGH) = NHALO + NEWNHALO
              DO J=IPE(NEIGH),IPE(NEIGH+1)-1
                IF (TRACE(IW(J)).EQ.NODE) THEN
                  CNT = CNT + 2
                END IF
              END DO
            END IF
          END DO
        END DO
        LAST_LVL_START = NHALO + 1
        NHALO = NHALO + NEWNHALO
      END SUBROUTINE
      SUBROUTINE GETHALOGRAPH(HALO,NHALO,N,IW,LW,IPE,IPTRHALO,JCNHALO,
     &     HALOEDGENBR,TRACE,NODE, GEN2HALO)
      INTEGER, INTENT(IN) :: N
      INTEGER,INTENT(IN):: NHALO, NODE
      INTEGER,INTENT(IN):: GEN2HALO(N)
      INTEGER,DIMENSION(NHALO),INTENT(IN) :: HALO
      INTEGER(8), INTENT(IN)              :: LW
      INTEGER(8), INTENT(IN)              :: IPE(N+1)
      INTEGER, INTENT(IN)     :: IW(LW), TRACE(N)
      INTEGER(8),INTENT(IN)   :: HALOEDGENBR
      INTEGER(8), INTENT(OUT) :: IPTRHALO(NHALO+1)
      INTEGER, INTENT(OUT)    :: JCNHALO(HALOEDGENBR)
      INTEGER::I,IPTR_CNT,JCN_CNT,HALOI
      INTEGER(8) :: J, CNT
      CNT = 0
      IPTR_CNT = 2
      JCN_CNT = 1
      IPTRHALO(1) = 1
      DO I=1,NHALO
         HALOI = HALO(I)
         DO J=IPE(HALOI),IPE(HALOI+1)-1
            IF (TRACE(IW(J))==NODE) THEN
               CNT = CNT + 1
               JCNHALO(JCN_CNT) = GEN2HALO(IW(J))
               JCN_CNT = JCN_CNT + 1
            END IF
         END DO
         IPTRHALO(IPTR_CNT) = CNT + 1
         IPTR_CNT = IPTR_CNT + 1
      END DO
      END SUBROUTINE
      SUBROUTINE GET_GROUPS(NHALO,PARTS,SEP,NSEP,NPARTS,
     &     CUT,NEWSEP,PERM,IPERM)
      INTEGER,INTENT(IN) :: NHALO,NSEP
      INTEGER,DIMENSION(:),INTENT(IN) :: SEP
      INTEGER,POINTER,DIMENSION(:)::PARTS
      INTEGER,POINTER,DIMENSION(:)::CUT,NEWSEP,PERM,
     &   IPERM
      INTEGER,INTENT(INOUT) :: NPARTS 
      INTEGER::I,CNT,NB_PARTS_WITHOUT_SEP_NODE
      INTEGER,DIMENSION(:),ALLOCATABLE::SIZES
      INTEGER,DIMENSION(:),ALLOCATABLE::PARTPTR
      ALLOCATE(NEWSEP(NSEP))
      ALLOCATE(PERM(NSEP))
      ALLOCATE(IPERM(NSEP))
      ALLOCATE(SIZES(NPARTS))
      ALLOCATE(PARTPTR(NPARTS+1))
      NB_PARTS_WITHOUT_SEP_NODE = 0
      SIZES = 0
      DO I=1,NSEP
        SIZES(PARTS(I)) = 
     &     SIZES(PARTS(I))+1
      END DO
      PARTPTR(1)=1
      DO I=2,NPARTS+1
        PARTPTR(I) = PARTPTR(I-1) + SIZES(I-1)
        IF (SIZES(I-1)==0) THEN
          NB_PARTS_WITHOUT_SEP_NODE = NB_PARTS_WITHOUT_SEP_NODE + 1
        END IF
      END DO
      ALLOCATE(CUT(NPARTS-NB_PARTS_WITHOUT_SEP_NODE+1))
      CUT(1) = 1
      CNT = 2
      DO I=2,NPARTS+1
        IF (SIZES(I-1).NE.0) THEN
          CUT(CNT) = PARTPTR(I) 
          CNT = CNT + 1
        END IF
      END DO
      NPARTS = NPARTS - NB_PARTS_WITHOUT_SEP_NODE 
      CUT(NPARTS+1) = NSEP+1
      DO I=1,NSEP
        NEWSEP(PARTPTR(PARTS(I))) = SEP(I)  
        PERM(PARTPTR(PARTS(I))) = I
        IPERM(I) = PARTPTR(PARTS(I))
        PARTPTR(PARTS(I)) = 
     &     PARTPTR(PARTS(I)) + 1
      END DO
      DEALLOCATE(SIZES,PARTPTR)
      END SUBROUTINE
      END MODULE DMUMPS_ANA_LR
