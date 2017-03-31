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
      MODULE DMUMPS_FAC_LR
      USE DMUMPS_LR_CORE
      USE DMUMPS_LR_TYPE
      USE DMUMPS_LR_STATS
      USE DMUMPS_ANA_LR
      IMPLICIT NONE
      CONTAINS
      SUBROUTINE DMUMPS_BLR_UPDATE_TRAILING_LDLT(
     &        A, LA, POSELT, IFLAG, IERROR, NFRONT,
     &        BEGS_BLR, NB_BLR, CURRENT_BLR, BLR_L,
     &        NELIM, IW2, BLOCK,
     &        MAXI_CLUSTER, NPIV, NIV, 
     &        COMPRESS_MID_PRODUCT, TOLEPS, KPERCENT)
!$    USE OMP_LIB
      INTEGER(8), intent(in)       :: LA
      INTEGER(8), intent(in)       :: POSELT 
      INTEGER, intent(in)          :: NFRONT, NB_BLR, CURRENT_BLR,
     &   NELIM, MAXI_CLUSTER, NPIV, NIV
      INTEGER, intent(out)         :: IFLAG, IERROR
      DOUBLE PRECISION, intent(inout)    :: A(LA)
      TYPE(LRB_TYPE),intent(in) :: BLR_L(NB_BLR-CURRENT_BLR)
      DOUBLE PRECISION, INTENT(INOUT), TARGET :: BLOCK(:,:)
      INTEGER, intent(in) :: IW2(*)
      INTEGER, POINTER, DIMENSION(:) :: BEGS_BLR
      INTEGER,intent(in) :: COMPRESS_MID_PRODUCT, KPERCENT
      DOUBLE PRECISION,intent(in) :: TOLEPS
      INTEGER :: I, NB_BLOCKS_PANEL, J, MID_RANK
      DOUBLE PRECISION, POINTER, DIMENSION(:) :: BLOCK_PTR
      LOGICAL :: BUILDQ
      INTEGER :: OMP_NUM
      INTEGER :: IBIS
#if defined(BLR_MT)
      INTEGER :: CHUNK
#endif
      INTEGER(8) :: POSELTT, POSELTD
      DOUBLE PRECISION :: ONE, MONE, ZERO
      PARAMETER (ONE = 1.0D0, MONE=-1.0D0)
      PARAMETER (ZERO=0.0D0)
      NB_BLOCKS_PANEL = NB_BLR-CURRENT_BLR
      POSELTD = POSELT + int(NFRONT,8) * int(BEGS_BLR(CURRENT_BLR)-1,8)
     &          + int(BEGS_BLR(CURRENT_BLR) - 1,8)
      OMP_NUM = 0
      BLOCK_PTR => BLOCK(1:MAXI_CLUSTER,1)
#if defined(BLR_MT) 
      CHUNK = 1
!$OMP DO SCHEDULE(DYNAMIC, CHUNK)
!$OMP& PRIVATE(I, J, POSELTT, OMP_NUM, BLOCK_PTR,
!$OMP&         MID_RANK, BUILDQ)
#endif
      DO IBIS = 1, (NB_BLOCKS_PANEL*(NB_BLOCKS_PANEL+1)/2) 
        IF (IFLAG.LT.0) CYCLE
        I = CEILING((1.0D0+SQRT(1.0D0+8.0D0*dble(IBIS)))/2.0D0)-1
        J = IBIS - I*(I-1)/2
#if defined(BLR_MT)         
!$      OMP_NUM = OMP_GET_THREAD_NUM() 
        BLOCK_PTR => BLOCK(1:MAXI_CLUSTER,OMP_NUM*MAXI_CLUSTER+1)
#endif
            POSELTT = POSELT + int(NFRONT,8) *
     &                int(BEGS_BLR(CURRENT_BLR+I)-1,8)
     &           + int(BEGS_BLR(CURRENT_BLR+J) - 1, 8)
            CALL DMUMPS_LRGEMM3('N', 'T', MONE,
     &            BLR_L(J),BLR_L(I), ONE, A, LA, 
     &            POSELTT, NFRONT, 1, NIV, IFLAG, IERROR, 
     &            COMPRESS_MID_PRODUCT, TOLEPS, KPERCENT,
     &            MID_RANK, BUILDQ,
     &            POSELTD, NFRONT, 
     &            IW2,
     &            BLOCK_PTR, 
     &            MAXI_CLUSTER)
            IF (IFLAG.LT.0) CYCLE
            CALL UPDATE_FLOP_STATS_LRB_PRODUCT(BLR_L(J), BLR_L(I), 'N',
     &          'T', NIV, COMPRESS_MID_PRODUCT, MID_RANK, BUILDQ
     &          , (I.EQ.J)
     &           )
         ENDDO
#if defined(BLR_MT) 
!$OMP END DO
#endif
      END SUBROUTINE DMUMPS_BLR_UPDATE_TRAILING_LDLT
      SUBROUTINE DMUMPS_SLAVE_BLR_UPD_TRAIL_LDLT(A, LA, POSELT, 
     &        IFLAG, IERROR, NCOL, NROW, POSBLOCFACTO, LD_BLOCFACTO, 
     &        BEGS_BLR_LM, NB_BLR_LM, BLR_LM, ISHIFT_LM,
     &        BEGS_BLR_LS, NB_BLR_LS, BLR_LS, ISHIFT_LS,
     &        CURRENT_BLR_LM, CURRENT_BLR_LS,
     &        IW2, BLOCK,
     &        MAXI_CLUSTER,
     &        COMPRESS_MID_PRODUCT, TOLEPS, KPERCENT
     &        )
!$    USE OMP_LIB      
      INTEGER(8), intent(in)  :: LA, POSBLOCFACTO
      DOUBLE PRECISION, intent(inout)  :: A(LA)
      INTEGER(8), intent(in)  :: POSELT 
      INTEGER, intent(out)    :: IFLAG, IERROR
      INTEGER, intent(in)     :: NCOL, NROW, IW2(*), 
     &                           MAXI_CLUSTER, LD_BLOCFACTO
      INTEGER, intent(in)     :: NB_BLR_LM, NB_BLR_LS, 
     &                           ISHIFT_LM, ISHIFT_LS, 
     &                           CURRENT_BLR_LM, CURRENT_BLR_LS
      DOUBLE PRECISION, INTENT(INOUT) :: BLOCK(MAXI_CLUSTER)
      INTEGER, POINTER, DIMENSION(:) :: BEGS_BLR_LM, BEGS_BLR_LS
      TYPE(LRB_TYPE),intent(in) :: BLR_LM(NB_BLR_LM-CURRENT_BLR_LM),
     &                             BLR_LS(NB_BLR_LS-CURRENT_BLR_LS)
      INTEGER,intent(in) :: COMPRESS_MID_PRODUCT, KPERCENT
      DOUBLE PRECISION,intent(in) :: TOLEPS
      INTEGER :: I, NB_BLOCKS_PANEL_LM, NB_BLOCKS_PANEL_LS, J, MID_RANK
      LOGICAL :: BUILDQ
      INTEGER :: IBIS
#if defined(BLR_MT)
      INTEGER :: CHUNK
#endif
      INTEGER(8) :: POSELTT, POSELTD
      DOUBLE PRECISION :: ONE, MONE, ZERO
      PARAMETER (ONE = 1.0D0, MONE=-1.0D0)
      PARAMETER (ZERO=0.0D0)
      NB_BLOCKS_PANEL_LM = NB_BLR_LM-CURRENT_BLR_LM
      NB_BLOCKS_PANEL_LS = NB_BLR_LS-CURRENT_BLR_LS
      POSELTD = POSBLOCFACTO 
#if defined(BLR_MT) 
      CHUNK = 1
!$OMP DO SCHEDULE(DYNAMIC,CHUNK)
!$OMP& PRIVATE(I, J, POSELTT, MID_RANK, BUILDQ) 
#endif
      DO IBIS = 1, (NB_BLOCKS_PANEL_LS*NB_BLOCKS_PANEL_LM) 
        IF (IFLAG.LT.0) CYCLE
        I = (IBIS-1)/NB_BLOCKS_PANEL_LM+1
        J = IBIS - (I-1)*NB_BLOCKS_PANEL_LM
            POSELTT = POSELT 
     &           + int(NCOL,8) * 
     &             int((BEGS_BLR_LS(CURRENT_BLR_LS+I)+ISHIFT_LS-1),8)
     &           + int((BEGS_BLR_LM(CURRENT_BLR_LM+J)+ISHIFT_LM-1),8)
            CALL DMUMPS_LRGEMM3('N', 'T', MONE,
     &            BLR_LM(J),BLR_LS(I), ONE, A, LA, 
     &            POSELTT, NCOL, 
     &            1, 2, IFLAG, IERROR, 
     &            COMPRESS_MID_PRODUCT, TOLEPS, KPERCENT,
     &            MID_RANK, BUILDQ,
     &            POSELTD, LD_BLOCFACTO, 
     &            IW2,
     &            BLOCK, 
     &            MAXI_CLUSTER)
            IF (IFLAG.LT.0) CYCLE
            CALL UPDATE_FLOP_STATS_LRB_PRODUCT(BLR_LM(J), BLR_LS(I),
     &           'N','T', 2, COMPRESS_MID_PRODUCT, MID_RANK, BUILDQ,
     &            .FALSE.) 
         ENDDO
#if defined(BLR_MT) 
!$OMP END DO
         IF (IFLAG.LT.0) RETURN
!$OMP DO SCHEDULE(DYNAMIC,CHUNK)
!$OMP& PRIVATE(I, J, POSELTT, MID_RANK, BUILDQ) 
#endif
         DO IBIS = 1, (NB_BLOCKS_PANEL_LS*(NB_BLOCKS_PANEL_LS+1)/2)
          IF (IFLAG.LT.0) CYCLE
          I = CEILING((1.0D0+SQRT(1.0D0+8.0D0*dble(IBIS)))/2.0D0)-1
          J = IBIS - I*(I-1)/2
          POSELTT = POSELT 
     &        + int(NCOL,8) * 
     &          int((BEGS_BLR_LS(CURRENT_BLR_LS+I)+ISHIFT_LS-1),8)
     &        + int((NCOL-NROW+(BEGS_BLR_LS(CURRENT_BLR_LS+J)-1)),8)
          CALL DMUMPS_LRGEMM3('N', 'T', MONE,
     &            BLR_LS(J),BLR_LS(I), ONE, A, LA, 
     &            POSELTT, NCOL, 
     &            1, 2, IFLAG, IERROR, 
     &            COMPRESS_MID_PRODUCT, TOLEPS, KPERCENT,
     &            MID_RANK, BUILDQ,
     &            POSELTD, LD_BLOCFACTO, 
     &            IW2,
     &            BLOCK, 
     &            MAXI_CLUSTER)
            IF (IFLAG.LT.0) CYCLE
            CALL UPDATE_FLOP_STATS_LRB_PRODUCT(BLR_LS(J), BLR_LS(I),
     &           'N','T', 2, COMPRESS_MID_PRODUCT, MID_RANK, BUILDQ,
     &            (I.EQ.J)) 
      ENDDO
#if defined(BLR_MT) 
!$OMP END DO
#endif
      END SUBROUTINE DMUMPS_SLAVE_BLR_UPD_TRAIL_LDLT
      SUBROUTINE DMUMPS_BLR_UPDATE_NELIM_VAR(
     &        A, LA, POSELT, IFLAG, IERROR, NFRONT,
     &        BEGS_BLR_L, BEGS_BLR_U, CURRENT_BLR, BLR_L, NB_BLR_L, 
     &        FIRST_BLOCK, NELIM, LBANDSLAVE, ISHIFT, NIV, SYM)
!$    USE OMP_LIB
      INTEGER(8), intent(in)       :: LA
      INTEGER(8), intent(in)       :: POSELT 
      INTEGER, intent(in)          :: NFRONT, NB_BLR_L, CURRENT_BLR,
     &                                NELIM,  SYM, NIV, FIRST_BLOCK
      LOGICAL, intent(in)          :: LBANDSLAVE
      INTEGER, intent(out)         :: IFLAG, IERROR
      INTEGER, intent(in)          :: ISHIFT
      DOUBLE PRECISION, TARGET, intent(inout) :: A(LA)
      TYPE(LRB_TYPE),TARGET,intent(in) :: BLR_L(NB_BLR_L-CURRENT_BLR)
      INTEGER, POINTER, DIMENSION(:) :: BEGS_BLR_L, BEGS_BLR_U
      INTEGER :: I, NB_BLOCKS_PANEL_L, KL, ML, NL, IS
      INTEGER :: allocok
      INTEGER(8) :: POSELT_INCB, POSELT_TOP
      DOUBLE PRECISION, ALLOCATABLE,DIMENSION(:,:) :: TEMP_BLOCK
      DOUBLE PRECISION :: ONE, MONE, ZERO
      PARAMETER (ONE = 1.0D0, MONE=-1.0D0)
      PARAMETER (ZERO=0.0D0)
      NB_BLOCKS_PANEL_L = NB_BLR_L-CURRENT_BLR
      IF (LBANDSLAVE) THEN
       IS = ISHIFT
      ELSE
       IS = 0
      ENDIF
#if defined(BLR_MT) 
!$OMP SINGLE
#endif
      IF (NELIM.NE.0) THEN 
        DO I = FIRST_BLOCK-CURRENT_BLR, NB_BLOCKS_PANEL_L
          KL = BLR_L(I)%K 
          ML = BLR_L(I)%M 
          NL = BLR_L(I)%N 
          IF (BLR_L(I)%ISLR) THEN
             IF (KL.GT.0) THEN
               allocate(TEMP_BLOCK( NELIM, KL ), stat=allocok )
               IF (allocok .GT. 0) THEN
                 IFLAG  = -13
                 IERROR = NELIM * KL
                 write(*,*) 'Allocation problem in BLR routine 
     &         DMUMPS_BLR_UPDATE_NELIM_VAR: ',
     &         'not enough memory? memory requested = ', IERROR
                 GOTO 100
               ENDIF
               POSELT_TOP  = POSELT 
     &         + int(NFRONT,8) * int((BEGS_BLR_U(CURRENT_BLR)-1),8)
     &         + int(BEGS_BLR_U(CURRENT_BLR+1) + IS - NELIM - 1,8)
               POSELT_INCB = POSELT 
     &         + int(NFRONT,8) * int((BEGS_BLR_L(CURRENT_BLR+I)-1),8)
     &         + int(BEGS_BLR_U(CURRENT_BLR+1)+IS-NELIM-1,8)
               CALL dgemm('N' , 'T' , NELIM, KL, NL , ONE ,
     &              A(POSELT_TOP) , NFRONT , BLR_L(I)%R(1,1) , KL ,
     &              ZERO , TEMP_BLOCK , NELIM) 
               CALL dgemm('N' , 'T' , NELIM , ML , KL , MONE ,
     &              TEMP_BLOCK , NELIM , BLR_L(I)%Q(1,1) , ML ,
     &              ONE , A(POSELT_INCB) , NFRONT) 
               deallocate(TEMP_BLOCK)
             ENDIF
          ELSE
            IF (SYM.EQ.0) THEN
              POSELT_TOP  = POSELT 
     &         + int(NFRONT,8) * int((BEGS_BLR_L(CURRENT_BLR)-1),8)
     &         + int(BEGS_BLR_U(CURRENT_BLR+1)+IS-NELIM-1,8)
              POSELT_INCB = POSELT 
     &         + int(NFRONT,8) * int((BEGS_BLR_L(CURRENT_BLR+I)-1),8)
     &         + int(BEGS_BLR_U(CURRENT_BLR+1) + IS - NELIM - 1, 8)
              CALL dgemm('N' , 'T' , NELIM, ML, NL , MONE ,
     &              A(POSELT_TOP) , NFRONT , BLR_L(I)%Q(1,1) , ML ,
     &              ONE , A(POSELT_INCB) , NFRONT) 
            ELSE
              POSELT_TOP  = POSELT + int(NFRONT,8) 
     &         * int(BEGS_BLR_U(CURRENT_BLR+1)+IS-NELIM-1,8)
     &         + int((BEGS_BLR_L(CURRENT_BLR)-1),8)
              POSELT_INCB = POSELT 
     &         + int(NFRONT,8) * int((BEGS_BLR_L(CURRENT_BLR+I)-1),8)
     &         + int(BEGS_BLR_U(CURRENT_BLR+1) + IS - NELIM - 1, 8)
              CALL dgemm('T' , 'T' , NELIM, ML, NL , MONE ,
     &              A(POSELT_TOP) , NFRONT , BLR_L(I)%Q(1,1) , ML ,
     &              ONE , A(POSELT_INCB) , NFRONT) 
            ENDIF
          ENDIF
        ENDDO
      ENDIF
 100  CONTINUE
#if defined(BLR_MT) 
!$OMP END SINGLE
#endif
      END SUBROUTINE DMUMPS_BLR_UPDATE_NELIM_VAR
      SUBROUTINE DMUMPS_BLR_UPDATE_TRAILING(
     &        A, LA, POSELT, IFLAG, IERROR, NFRONT,
     &        BEGS_BLR_L, BEGS_BLR_U, CURRENT_BLR, BLR_L, NB_BLR_L, 
     &        BLR_U,
     &        NB_BLR_U, NELIM, LBANDSLAVE, ISHIFT, NIV, SYM, K470,
     &        COMPRESS_MID_PRODUCT, TOLEPS, KPERCENT)
!$    USE OMP_LIB
      INTEGER(8), intent(in)       :: LA
      INTEGER(8), intent(in)       :: POSELT 
      INTEGER, intent(in)          :: NFRONT, NB_BLR_L, NB_BLR_U, 
     &                                CURRENT_BLR, K470,
     &                                NELIM, NIV, SYM
      INTEGER, intent(out)         :: IFLAG, IERROR
      LOGICAL, intent(in)          :: LBANDSLAVE
      INTEGER, intent(in)          :: ISHIFT
      DOUBLE PRECISION, TARGET, intent(inout) :: A(LA)
      TYPE(LRB_TYPE),TARGET,intent(in) :: BLR_U(NB_BLR_U-CURRENT_BLR)
      TYPE(LRB_TYPE),TARGET,intent(in) :: BLR_L(NB_BLR_L-CURRENT_BLR)
      INTEGER, POINTER, DIMENSION(:) :: BEGS_BLR_L, BEGS_BLR_U
      INTEGER,intent(in) :: COMPRESS_MID_PRODUCT, KPERCENT
      DOUBLE PRECISION,intent(in) :: TOLEPS
      INTEGER :: I, NB_BLOCKS_PANEL_L, NB_BLOCKS_PANEL_U, 
     &           KL, ML, NL, J, IS, MID_RANK
      INTEGER :: allocok
      LOGICAL :: BUILDQ
      INTEGER :: OMP_NUM
      CHARACTER(len=1) :: TRANSB1
      INTEGER :: IBIS
#if defined(BLR_MT)
      INTEGER :: CHUNK
#endif
      INTEGER(8) :: POSELT_INCB, POSELT_TOP
      DOUBLE PRECISION, ALLOCATABLE,DIMENSION(:,:) :: TEMP_BLOCK
      DOUBLE PRECISION :: ONE, MONE, ZERO
      PARAMETER (ONE = 1.0D0, MONE=-1.0D0)
      PARAMETER (ZERO=0.0D0)
      NB_BLOCKS_PANEL_L = NB_BLR_L-CURRENT_BLR
      NB_BLOCKS_PANEL_U = NB_BLR_U-CURRENT_BLR
      IF (LBANDSLAVE) THEN
       IS = ISHIFT
      ELSE
       IS = 0
      ENDIF
#if defined(BLR_MT)
!$OMP SINGLE
#endif
      IF (NELIM.NE.0) THEN 
         DO I = 1, NB_BLOCKS_PANEL_L
            KL = BLR_L(I)%K 
            ML = BLR_L(I)%M 
            NL = BLR_L(I)%N 
            IF (BLR_L(I)%ISLR) THEN
               IF (KL.GT.0) THEN
               allocate(TEMP_BLOCK( NELIM, KL ), stat=allocok )
               IF (allocok .GT. 0) THEN
                 IFLAG  = -13
                 IERROR = NELIM * KL
                 write(*,*) 'Allocation problem in BLR routine 
     &           DMUMPS_BLR_UPDATE_TRAILING: ',
     &           'not enough memory? memory requested = ', IERROR
                 GOTO 100
               ENDIF
               POSELT_TOP  = POSELT 
     &           + int(NFRONT,8) * int((BEGS_BLR_U(CURRENT_BLR)-1),8)
     &           + int(BEGS_BLR_U(CURRENT_BLR+1) + IS - NELIM - 1,8)
               POSELT_INCB = POSELT 
     &           + int(NFRONT,8) * int((BEGS_BLR_L(CURRENT_BLR+I)-1),8)
     &           + int(BEGS_BLR_U(CURRENT_BLR+1)+IS-NELIM-1,8)
               CALL dgemm('N' , 'T' , NELIM, KL, NL , ONE ,
     &                A(POSELT_TOP) , NFRONT , BLR_L(I)%R(1,1) , KL ,
     &                ZERO , TEMP_BLOCK , NELIM) 
               CALL dgemm('N' , 'T' , NELIM , ML , KL , MONE ,
     &                TEMP_BLOCK , NELIM , BLR_L(I)%Q(1,1) , ML ,
     &                ONE , A(POSELT_INCB) , NFRONT) 
               deallocate(TEMP_BLOCK)
               ENDIF
            ELSE
              POSELT_TOP  = POSELT 
     &           + int(NFRONT,8) * int((BEGS_BLR_L(CURRENT_BLR)-1),8)
     &           + int(BEGS_BLR_U(CURRENT_BLR+1)+IS-NELIM-1,8)
              POSELT_INCB = POSELT 
     &           + int(NFRONT,8) * int((BEGS_BLR_L(CURRENT_BLR+I)-1),8)
     &           + int(BEGS_BLR_U(CURRENT_BLR+1) + IS - NELIM - 1, 8)
               CALL dgemm('N' , 'T' , NELIM, ML, NL , MONE ,
     &                A(POSELT_TOP) , NFRONT , BLR_L(I)%Q(1,1) , ML ,
     &                ONE , A(POSELT_INCB) , NFRONT) 
            ENDIF
         ENDDO
      ENDIF
 100  CONTINUE
#if defined(BLR_MT) 
!$OMP END SINGLE
#endif
      IF (IFLAG.LT.0) GOTO 200
      OMP_NUM = 0
#if defined(BLR_MT)
      CHUNK = 1
!$OMP DO SCHEDULE(DYNAMIC,CHUNK) 
!$OMP& PRIVATE(I, J, POSELT_INCB, MID_RANK, BUILDQ)
#endif
      DO IBIS = 1, (NB_BLOCKS_PANEL_L*NB_BLOCKS_PANEL_U) 
        IF (IFLAG.LT.0) CYCLE
        I = (IBIS-1)/NB_BLOCKS_PANEL_U+1
        J = IBIS - (I-1)*NB_BLOCKS_PANEL_U
            POSELT_INCB = POSELT 
     &           + int(NFRONT,8) * int((BEGS_BLR_L(CURRENT_BLR+I)-1),8)
     &           + int(BEGS_BLR_U(CURRENT_BLR+J) +IS - 1,8)
         IF (SYM.EQ.0) THEN
           IF (K470.EQ.1) THEN
             TRANSB1 = 'N'
           ELSE
             TRANSB1 = 'T'
           ENDIF
            CALL DMUMPS_LRGEMM3(TRANSB1, 'T', MONE, BLR_U(J),
     &            BLR_L(I), ONE, A, LA, POSELT_INCB,
     &            NFRONT, 0, NIV, IFLAG, IERROR, 
     &            COMPRESS_MID_PRODUCT, TOLEPS,
     &            KPERCENT, MID_RANK, BUILDQ)
            IF (IFLAG.LT.0) CYCLE
            CALL UPDATE_FLOP_STATS_LRB_PRODUCT(BLR_U(J), BLR_L(I),
     &          TRANSB1,
     &          'T', NIV, COMPRESS_MID_PRODUCT, MID_RANK, BUILDQ)
         ELSE
            CALL DMUMPS_LRGEMM3('N', 'T', MONE, BLR_U(J),
     &            BLR_L(I), ONE, A, LA, POSELT_INCB,
     &            NFRONT, 0, NIV, IFLAG, IERROR, 
     &            COMPRESS_MID_PRODUCT, TOLEPS,
     &            KPERCENT, MID_RANK, BUILDQ)
            IF (IFLAG.LT.0) CYCLE
            CALL UPDATE_FLOP_STATS_LRB_PRODUCT(BLR_U(J), BLR_L(I), 'N',
     &          'T', NIV, COMPRESS_MID_PRODUCT, MID_RANK, BUILDQ)
         ENDIF
         ENDDO
#if defined(BLR_MT) 
!$OMP END DO
#endif
 200  CONTINUE
      END SUBROUTINE DMUMPS_BLR_UPDATE_TRAILING
      SUBROUTINE DMUMPS_DECOMPRESS_PANEL(A, LA, POSELT, NFRONT,
     &        COPY_DENSE_BLOCKS,
     &        BEGS_BLR_DIAG, BEGS_BLR_FIRST_OFFDIAG,
     &        NB_BLR, BLR_PANEL, CURRENT_BLR, DIR,
     &        LD_OR_NPIV, K470,
     &        BEG_I_IN, END_I_IN)
!$    USE OMP_LIB 
      INTEGER(8), intent(in)       :: LA
      DOUBLE PRECISION, intent(inout)       :: A(LA)
      INTEGER(8), intent(in)       :: POSELT 
      LOGICAL, intent(in)          :: COPY_DENSE_BLOCKS  
      INTEGER, intent(in)          :: NFRONT, NB_BLR, CURRENT_BLR
      INTEGER, intent(in)          :: BEGS_BLR_DIAG, 
     &                                BEGS_BLR_FIRST_OFFDIAG
      TYPE(LRB_TYPE),intent(inout) :: BLR_PANEL(NB_BLR-CURRENT_BLR)
      CHARACTER(len=1) :: DIR
      INTEGER, intent(in) :: LD_OR_NPIV, K470
      INTEGER,OPTIONAL,intent(in) :: BEG_I_IN, END_I_IN
      INTEGER :: IP, M, N, BIP, BEG_I, END_I
#if defined(BLR_MT)
      INTEGER :: LAST_IP, CHUNK
#endif
      INTEGER :: K, I
      INTEGER(8) :: POSELT_BLOCK, NFRONT8, LD_BLK_IN_FRONT
      DOUBLE PRECISION :: ONE, ALPHA, ZERO
      PARAMETER (ONE = 1.0D0, ALPHA=-1.0D0)
      PARAMETER (ZERO = 0.0D0)
      IF(present(BEG_I_IN)) THEN
        BEG_I = BEG_I_IN
      ELSE
        BEG_I = CURRENT_BLR+1
      ENDIF
      IF(present(END_I_IN)) THEN
        END_I = END_I_IN
      ELSE
        END_I = NB_BLR
      ENDIF
      NFRONT8         = int(NFRONT,8)
      LD_BLK_IN_FRONT = NFRONT8
      BIP             = BEGS_BLR_FIRST_OFFDIAG
#if defined(BLR_MT)
      LAST_IP = BEG_I
      CHUNK = 1
!$OMP PARALLEL DO PRIVATE(POSELT_BLOCK, M, N, K, I)
!$OMP&         FIRSTPRIVATE(BIP, LAST_IP) SCHEDULE(DYNAMIC, CHUNK)
#endif
      DO IP = BEG_I, END_I
#if defined(BLR_MT)
        DO I = 1, IP - LAST_IP
           IF (DIR .eq. 'V') THEN
              BIP  = BIP +  BLR_PANEL(LAST_IP-CURRENT_BLR+I-1)%M
           ELSE 
             IF (K470.EQ.1) THEN
               BIP  = BIP +  BLR_PANEL(LAST_IP-CURRENT_BLR+I-1)%M
             ELSE
               BIP  = BIP +  BLR_PANEL(LAST_IP-CURRENT_BLR+I-1)%N
             ENDIF
           ENDIF  
        ENDDO
        LAST_IP = IP
#endif
        IF (DIR .eq. 'V') THEN
           IF (BIP .LE. LD_OR_NPIV) THEN
             POSELT_BLOCK = POSELT + NFRONT8*int(BIP-1,8) + 
     &              int(BEGS_BLR_DIAG - 1,8)
           ELSE
             POSELT_BLOCK = POSELT +NFRONT8*int(LD_OR_NPIV,8)+
     &              int(BEGS_BLR_DIAG - 1,8)
             POSELT_BLOCK = POSELT_BLOCK +
     &                      int(LD_OR_NPIV,8)*int(BIP-1-LD_OR_NPIV,8)
             LD_BLK_IN_FRONT=int(LD_OR_NPIV,8)
           ENDIF
        ELSE 
         POSELT_BLOCK = POSELT +
     &              NFRONT8*int(BEGS_BLR_DIAG-1,8) +
     &              int(BIP - 1,8)
        ENDIF
           M = BLR_PANEL(IP-CURRENT_BLR)%M
           N = BLR_PANEL(IP-CURRENT_BLR)%N
           K = BLR_PANEL(IP-CURRENT_BLR)%K
        IF ((BLR_PANEL(IP-CURRENT_BLR)%ISLR).AND.
     &      (BLR_PANEL(IP-CURRENT_BLR)%LRFORM.EQ.1)) THEN
        IF (K.EQ.0) THEN
          IF (K470.NE.1.OR.DIR .eq. 'V') THEN
            DO I = 1, M
            A(POSELT_BLOCK+int(I-1,8)*LD_BLK_IN_FRONT :
     &                 POSELT_BLOCK+int(I-1,8)*LD_BLK_IN_FRONT
     &                                                + int(N-1,8)) 
     &               = ZERO
            ENDDO
          ELSE 
            DO I = 1, N
            A(POSELT_BLOCK+int(I-1,8)*NFRONT8:
     &                  POSELT_BLOCK+int(I-1,8)*NFRONT8 + int(M-1,8))
     &               = ZERO
            ENDDO
          ENDIF
          GOTO 1800
        ENDIF
        IF (K470.NE.1.OR.DIR .eq. 'V') THEN
           CALL dgemm('T', 'T', N, M, K, ONE ,
     &          BLR_PANEL(IP-CURRENT_BLR)%R(1,1) , K, 
     &          BLR_PANEL(IP-CURRENT_BLR)%Q(1,1) , M, 
     &          ZERO, A(POSELT_BLOCK), int(LD_BLK_IN_FRONT))
        ELSE 
           CALL dgemm('N', 'N', M, N, K, ONE ,
     &          BLR_PANEL(IP-CURRENT_BLR)%Q(1,1) , M, 
     &          BLR_PANEL(IP-CURRENT_BLR)%R(1,1) , K, 
     &          ZERO, A(POSELT_BLOCK), NFRONT)
        ENDIF
         ELSE  IF (COPY_DENSE_BLOCKS) THEN
             IF (K470.NE.1.OR.DIR .eq. 'V') THEN
                 DO I = 1, M
                     A(POSELT_BLOCK+int(I-1,8)*LD_BLK_IN_FRONT :
     &                 POSELT_BLOCK+int(I-1,8)*LD_BLK_IN_FRONT
     &                                                + int(N-1,8)) 
     &               = BLR_PANEL(IP-CURRENT_BLR)%Q(I,1:N)
                 ENDDO
             ELSE 
                 DO I = 1, N
                     A(POSELT_BLOCK+int(I-1,8)*NFRONT8:
     &                  POSELT_BLOCK+int(I-1,8)*NFRONT8 + int(M-1,8))
     &               = BLR_PANEL(IP-CURRENT_BLR)%Q(1:M,I)
                 ENDDO
             ENDIF
        ENDIF
 1800   CONTINUE
#if !defined(BLR_MT)
        IF (DIR .eq. 'V') THEN
           BIP  = BIP +  BLR_PANEL(IP-CURRENT_BLR)%M
        ELSE 
          IF (K470.EQ.1) THEN
            BIP  = BIP +  BLR_PANEL(IP-CURRENT_BLR)%M
          ELSE
            BIP  = BIP +  BLR_PANEL(IP-CURRENT_BLR)%N
          ENDIF
        ENDIF
#endif
      END DO 
#if defined(BLR_MT)
!$OMP END PARALLEL DO
#endif
      END SUBROUTINE DMUMPS_DECOMPRESS_PANEL
      SUBROUTINE DMUMPS_FAKE_COMPRESS_CB(A, LA, POSELT, NFRONT,
     &        BEGS_BLR_L, NB_BLR_L, 
     &        BEGS_BLR_U,  NB_BLR_U, NPARTSASS_U,
     &        TOLEPS, NASS, NROW,
     &        SYM, WORK, TAU, JPVT, LWORK, RWORK,
     &        BLOCK, MAXI_CLUSTER, INODE, NIV, 
     &        LBANDSLAVE, ISHIFT,KPERCENT) 
      INTEGER(8), intent(in)       :: LA
      DOUBLE PRECISION, intent(inout)       :: A(LA)
      INTEGER(8), intent(in)       :: POSELT 
      INTEGER, intent(in)          :: NFRONT, INODE
      INTEGER, INTENT(IN)          :: NIV, NROW, KPERCENT
      INTEGER                      :: MAXI_CLUSTER, LWORK, SYM, NASS, 
     &                                NB_BLR_L, NB_BLR_U, NPARTSASS_U
      DOUBLE PRECISION,intent(in)              :: TOLEPS
      LOGICAL, intent(in)          :: LBANDSLAVE
      INTEGER, intent(in)          :: ISHIFT
      INTEGER, POINTER, DIMENSION(:) :: BEGS_BLR_L,  BEGS_BLR_U
      DOUBLE PRECISION :: BLOCK(MAXI_CLUSTER,MAXI_CLUSTER)
      DOUBLE PRECISION,DIMENSION(:) :: RWORK
      DOUBLE PRECISION, DIMENSION(:) :: WORK, TAU
      INTEGER, DIMENSION(:) :: JPVT
      INTEGER :: M, N, NCB, BEGLOOP, RANK, MAXRANK, FRONT_CB_BLR_SAVINGS
      INTEGER :: INFO, I, J, JJ, IB, JDEB, IS
      INTEGER :: allocok, MREQ
      INTEGER(8) :: POSELT_BLOCK
      DOUBLE PRECISION :: HR_COST, BUILDQ_COST, CB_DEMOTE_COST, 
     &                                          CB_PROMOTE_COST
      INTEGER T1, T2, COUNT_RATE
      DOUBLE PRECISION :: LOC_PROMOTING_TIME
      DOUBLE PRECISION :: LOC_CB_DEMOTING_TIME
      DOUBLE PRECISION, ALLOCATABLE :: R(:,:)
      DOUBLE PRECISION :: ONE, ZERO
      PARAMETER (ONE = 1.0D0)
      PARAMETER (ZERO = 0.0D0)
      LOC_PROMOTING_TIME = 0.0D0 
      LOC_CB_DEMOTING_TIME = 0.0D0 
      CB_DEMOTE_COST = 0.0D0
      CB_PROMOTE_COST = 0.0D0
      allocate(R(MAXI_CLUSTER,MAXI_CLUSTER),stat=allocok)
      IF (allocok .GT. 0) THEN
          MREQ=MAXI_CLUSTER*MAXI_CLUSTER
          write(*,*) 'Allocation problem in BLR routine 
     &    DMUMPS_FAKE_COMPRESS_CB: ',
     &    'not enough memory? memory requested = ', MREQ
          CALL MUMPS_ABORT()
      ENDIF
      FRONT_CB_BLR_SAVINGS = 0
      NCB = NFRONT - NASS 
      IF (NCB.LE.0) RETURN
      IF (LBANDSLAVE) THEN
       IS = ISHIFT
      ELSE
       IS = 0
      ENDIF
      DO J =  NPARTSASS_U+1, NB_BLR_U
         IF (NIV.EQ.1) THEN
           IF (SYM.GT.0) THEN
             BEGLOOP = J
           ELSE
             BEGLOOP = NPARTSASS_U + 1
           ENDIF
         ELSE
           BEGLOOP = 2
         ENDIF
         IF ((BEGS_BLR_U(J+1)+IS).LE.NASS+1) CYCLE
         JDEB = max(BEGS_BLR_U(J)+IS,NASS+1)
          N = BEGS_BLR_U(J+1)+IS-JDEB
         DO I = BEGLOOP, NB_BLR_L
           CALL SYSTEM_CLOCK(T1)
           JPVT = 0
           M = BEGS_BLR_L(I+1)-BEGS_BLR_L(I)
           POSELT_BLOCK = POSELT 
     &           + int(NFRONT,8) * int((BEGS_BLR_L(I)-1),8)
     &           + int(JDEB - 1,8)
           DO IB=1,M
               IF((I.EQ.J).AND.(SYM.GT.0).AND.(NIV.EQ.1)) THEN
                   BLOCK(IB,1:IB) =
     &             A( POSELT_BLOCK+int((IB-1),8)*int(NFRONT,8) :
     &                POSELT_BLOCK+
     &                int((IB-1),8)*int(NFRONT,8)+int(IB-1,8)  )
                   BLOCK(1:IB-1,IB) = BLOCK(IB,1:IB-1)
               ELSE
                   BLOCK(IB,1:N) =
     &              A( POSELT_BLOCK+int((IB-1),8)*int(NFRONT,8) :
     &             POSELT_BLOCK+int((IB-1),8)*int(NFRONT,8)+int(N-1,8) )
               ENDIF
           END DO  
           MAXRANK = floor(dble(M*N)/dble(M+N))
           MAXRANK = max (1, int((MAXRANK*KPERCENT/100)))
           CALL DMUMPS_TRUNCATED_RRQR( M, N, BLOCK(1,1),
     &          MAXI_CLUSTER, JPVT(1), TAU(1), WORK(1), N, 
     &          RWORK(1), TOLEPS, RANK, MAXRANK, INFO )
           CALL SYSTEM_CLOCK(T2,COUNT_RATE)
           LOC_CB_DEMOTING_TIME = LOC_CB_DEMOTING_TIME 
     &          + DBLE(T2-T1)/DBLE(COUNT_RATE)
            IF (INFO < 0) THEN
               WRITE(*,*) " PROBLEM IN ARGUMENT NUMBER ",INFO,
     &                    " OF TRUNCATED_RRQR WHILE COMPRESSING A BLOCK
     &                      IN CB (FAKE COMPRESSION anyway) "
                CALL MUMPS_ABORT()
            END IF
            HR_COST =  4.0D0*dble(RANK)*dble(RANK)*dble(RANK)/3.0D0 
     &               + 4.0D0*dble(RANK)*dble(M)*dble(N) 
     &               - 2.0D0*dble((M+N))*dble(RANK)*dble(RANK) 
            IF (RANK.LE.MAXRANK) THEN 
               CALL SYSTEM_CLOCK(T1)
               DO JJ=1, N
                  R(1:MIN(RANK,JJ),JPVT(JJ)) =
     &                 BLOCK(1:MIN(RANK,JJ),JJ)
                  IF(JJ.LT.RANK) R(MIN(RANK,JJ)+1:
     &                 RANK,JPVT(JJ))= ZERO
               END DO
               CALL dorgqr(M, RANK, RANK,
     &              BLOCK(1,1), MAXI_CLUSTER,
     &              TAU(1), WORK(1), LWORK, INFO)
               CALL dgemm('T', 'T', N, M, RANK, ONE ,
     &              R , MAXI_CLUSTER, 
     &              BLOCK(1,1) , MAXI_CLUSTER, 
     &              ZERO, A(POSELT_BLOCK), NFRONT)
               CALL SYSTEM_CLOCK(T2,COUNT_RATE)
               LOC_PROMOTING_TIME = LOC_PROMOTING_TIME + 
     &              DBLE(T2-T1)/DBLE(COUNT_RATE)
               BUILDQ_COST = 4.0D0*dble(RANK)*dble(RANK)*dble(M)
     &                      - dble(RANK)*dble(RANK)*dble(RANK)
     &                       
               CB_DEMOTE_COST = CB_DEMOTE_COST +
     &                  (HR_COST+BUILDQ_COST) 
               CB_PROMOTE_COST = CB_PROMOTE_COST +
     &                      2.0D0*dble(RANK)*dble(M)*dble(N)
               FRONT_CB_BLR_SAVINGS = FRONT_CB_BLR_SAVINGS + 
     &                      (M-RANK)*(N-RANK)-RANK*RANK
            ELSE 
               CB_DEMOTE_COST = CB_DEMOTE_COST + HR_COST 
            END IF
         END DO 
      END DO
      deallocate(R)
      CALL STATS_COMPUTE_MRY_FRONT_CB(NCB, NROW, SYM, NIV, INODE,
     &                        FRONT_CB_BLR_SAVINGS)
      CALL UPDATE_FLOP_STATS_CB_DEMOTE(CB_DEMOTE_COST, NIV)
      CALL UPDATE_FLOP_STATS_CB_PROMOTE(CB_PROMOTE_COST, NIV)
      CALL UPDATE_CB_DEMOTING_TIME(INODE, LOC_CB_DEMOTING_TIME)
      CALL UPDATE_PROMOTING_TIME(INODE, LOC_PROMOTING_TIME)
      END SUBROUTINE DMUMPS_FAKE_COMPRESS_CB
      SUBROUTINE DMUMPS_COMPRESS_PANEL(
     &        A, LA, POSELT, IFLAG, IERROR, NFRONT,
     &        BEGS_BLR, NB_BLR, TOLEPS, K473, BLR_PANEL, CURRENT_BLR,
     &        DIR, WORK, TAU, JPVT, 
     &        LWORK, RWORK, BLOCK,
     &        MAXI_CLUSTER, NELIM, 
     &        LBANDSLAVE, NPIV, ISHIFT, NIV, KPERCENT, 
     &        K470, KEEP8, K480,
     &        BEG_I_IN, END_I_IN
     &        )
!$    USE OMP_LIB
      INTEGER(8), intent(in)       :: LA
      INTEGER(8), intent(in)       :: POSELT 
      INTEGER, intent(in)          :: NFRONT, NB_BLR, CURRENT_BLR, NIV
      INTEGER, intent(out)          :: IFLAG, IERROR
      TYPE(LRB_TYPE), intent(inout) :: BLR_PANEL(NB_BLR-CURRENT_BLR)
      DOUBLE PRECISION, intent(inout)       :: A(LA)
      DOUBLE PRECISION, TARGET, DIMENSION(:) :: RWORK
      DOUBLE PRECISION, TARGET, DIMENSION(:,:) :: BLOCK
      DOUBLE PRECISION, TARGET, DIMENSION(:) :: WORK, TAU
      INTEGER, TARGET, DIMENSION(:) :: JPVT
      INTEGER, POINTER :: BEGS_BLR(:)
      INTEGER(8) :: KEEP8(150)
      INTEGER, OPTIONAL, intent(in) :: K480
      INTEGER,OPTIONAL,intent(in) :: BEG_I_IN, END_I_IN
      INTEGER, intent(in)          :: NPIV, ISHIFT, KPERCENT, K473, K470
      LOGICAL, intent(in)          :: LBANDSLAVE
      INTEGER                      :: MAXI_CLUSTER, LWORK, NELIM
      DOUBLE PRECISION,intent(in)              :: TOLEPS
      CHARACTER(len=1) :: DIR
      INTRINSIC maxval
      INTEGER :: IP, NB_BLOCKS_PANEL, M, N, RANK, MAXRANK
      INTEGER :: INFO, I, J, IS, BEG_I, END_I
      INTEGER(8) :: POSELT_BLOCK
      LOGICAL :: ISLR
      DOUBLE PRECISION :: ONE, ALPHA, ZERO
      PARAMETER (ONE = 1.0D0, ALPHA=-1.0D0)
      PARAMETER (ZERO = 0.0D0)
      INTEGER :: OMP_NUM
      DOUBLE PRECISION, POINTER, DIMENSION(:) :: RWORK_THR
      DOUBLE PRECISION, POINTER, DIMENSION(:,:) :: BLOCK_THR
      DOUBLE PRECISION, POINTER, DIMENSION(:) :: WORK_THR, TAU_THR
      INTEGER, POINTER, DIMENSION(:) :: JPVT_THR
#if defined(BLR_MT) 
      INTEGER :: CHUNK
#endif
      IF(present(BEG_I_IN)) THEN
        BEG_I = BEG_I_IN
      ELSE
        BEG_I = CURRENT_BLR+1
      ENDIF
      IF(present(END_I_IN)) THEN
        END_I = END_I_IN
      ELSE
        END_I = NB_BLR
      ENDIF
      IF (LBANDSLAVE) THEN
       IS = ISHIFT
      ELSE
       IS=0
      ENDIF
      IF (DIR .eq. 'V') THEN
         IF (LBANDSLAVE) THEN
          N = NPIV
         ELSE
          N = BEGS_BLR(CURRENT_BLR+1)-BEGS_BLR(CURRENT_BLR)-NELIM
         ENDIF
      ELSE IF (DIR .eq. 'H') THEN
        IF (K470.EQ.1) THEN
          N = BEGS_BLR(CURRENT_BLR+1)-BEGS_BLR(CURRENT_BLR)-NELIM
        ELSE
          M = BEGS_BLR(CURRENT_BLR+1)-BEGS_BLR(CURRENT_BLR)-NELIM
        ENDIF
      ELSE
         WRITE(*,*) " WRONG ARGUMENT IN DMUMPS_COMPRESS_PANEL "
         CALL MUMPS_ABORT()
      END IF
      NB_BLOCKS_PANEL = NB_BLR-CURRENT_BLR
      OMP_NUM = 0
#if defined(BLR_MT) 
      CHUNK = 1
!$OMP DO PRIVATE(INFO, POSELT_BLOCK, RANK, MAXRANK, I, J, OMP_NUM)
!$OMP&   SCHEDULE(DYNAMIC,CHUNK)
#endif
      DO IP = BEG_I, END_I
        IF (IFLAG.LT.0) CYCLE
#if defined(BLR_MT)         
!$      OMP_NUM = OMP_GET_THREAD_NUM()
#endif
        BLOCK_THR => BLOCK(1:MAXI_CLUSTER,OMP_NUM*MAXI_CLUSTER+1:
     &                  (OMP_NUM+1)*MAXI_CLUSTER) 
        JPVT_THR  => JPVT(OMP_NUM*MAXI_CLUSTER+1:
     &                 (OMP_NUM+1)*MAXI_CLUSTER) 
        TAU_THR   => TAU(OMP_NUM*MAXI_CLUSTER+1:
     &                 (OMP_NUM+1)*MAXI_CLUSTER) 
        WORK_THR  => WORK(OMP_NUM*LWORK+1:
     &                 (OMP_NUM+1)*LWORK) 
        RWORK_THR => RWORK(OMP_NUM*2*MAXI_CLUSTER+1:
     &                  (OMP_NUM+1)*2*MAXI_CLUSTER) 
        IF (DIR .eq. 'V') THEN
            M = BEGS_BLR(IP+1)-BEGS_BLR(IP)
            POSELT_BLOCK = POSELT + 
     &              int(NFRONT,8) * int(BEGS_BLR(IP)-1,8) + 
     &              int(BEGS_BLR(CURRENT_BLR) + IS - 1,8)
        ELSE 
          IF (K470.EQ.1) THEN
            M = BEGS_BLR(IP+1)-BEGS_BLR(IP)
          ELSE
            N = BEGS_BLR(IP+1)-BEGS_BLR(IP)
          ENDIF
          POSELT_BLOCK = POSELT +
     &              int(NFRONT,8)*int(BEGS_BLR(CURRENT_BLR)-1,8) + 
     &              int( BEGS_BLR(IP) - 1,8)
        END IF
        JPVT_THR(1:MAXI_CLUSTER) = 0
        IF (K473.EQ.1) THEN
            MAXRANK = 1
            RANK = MAXRANK+1
            INFO = 0
            GOTO 3800
        ENDIF
        IF (K470.NE.1.OR.DIR .eq. 'V') THEN
            DO I=1,M
                BLOCK_THR(I,1:N)=
     &          A( POSELT_BLOCK+int(I-1,8)*int(NFRONT,8) :
     &          POSELT_BLOCK+int(I-1,8)*int(NFRONT,8)+int(N-1,8) )
            END DO  
        ELSE 
            DO I=1,N
                BLOCK_THR(1:M,I)=
     &          A( POSELT_BLOCK+int(I-1,8)*int(NFRONT,8) :
     &          POSELT_BLOCK+int(I-1,8)*int(NFRONT,8)+int(M-1,8) )
            END DO  
        END IF
        MAXRANK = floor(dble(M*N)/dble(M+N))
        MAXRANK = max (1, int((MAXRANK*KPERCENT/100)))
        CALL DMUMPS_TRUNCATED_RRQR( M, N,
     &       BLOCK_THR(1,1),
     &       MAXI_CLUSTER, JPVT_THR(1), 
     &       TAU_THR(1), 
     &       WORK_THR(1), N, 
     &       RWORK_THR(1), 
     &       TOLEPS, RANK, MAXRANK, INFO)
 3800 CONTINUE
        IF (INFO < 0) THEN
           WRITE(*,*) " PROBLEM IN ARGUMENT NUMBER ",INFO,
     &                 " OF TRUNCATED_RRQR WHILE COMPRESSING A BLOCK "
            CALL MUMPS_ABORT()
        END IF
        ISLR = ((RANK.LE.MAXRANK).AND.(M.NE.0).AND.(N.NE.0))     
        CALL ALLOC_LRB(BLR_PANEL(IP-CURRENT_BLR), RANK, RANK, 
     &                 M, N, ISLR, IFLAG, IERROR, KEEP8)
        IF (IFLAG.LT.0) CYCLE
        IF (ISLR) THEN 
           IF (RANK .EQ. 0) THEN 
           ELSE 
               BLR_PANEL(IP-CURRENT_BLR)%Q = ZERO
               DO I=1,RANK
                 BLR_PANEL(IP-CURRENT_BLR)%Q(I,I) = ONE
               END DO
               CALL dormqr 
     &           ('L', 'N', M, RANK, RANK,
     &           BLOCK_THR(1,1), 
     &           MAXI_CLUSTER, TAU_THR(1), 
     &           BLR_PANEL(IP-CURRENT_BLR)%Q(1,1),
     &           M, WORK_THR(1), LWORK, INFO )
               IF (INFO < 0) THEN
                 WRITE(*,*) " PROBLEM IN ARGUMENT NUMBER ",INFO,
     &                     " OF CUNMQR WHILE COMPRESSING A BLOCK "
                 CALL MUMPS_ABORT()
               END IF
               DO J=1,N
                 BLR_PANEL(IP-CURRENT_BLR)%R(1:MIN(RANK,J),
     &               JPVT_THR(J)) =
     &               BLOCK_THR(1:MIN(RANK,J),J)
                 IF(J.LT.RANK) BLR_PANEL(IP-CURRENT_BLR)%
     &               R(MIN(RANK,J)+1:RANK,JPVT_THR(J))= ZERO
               ENDDO
             CALL UPDATE_FLOP_STATS_DEMOTE(
     &               BLR_PANEL(IP-CURRENT_BLR), NIV)
          END IF
        ELSE 
           IF (K470.NE.1.OR.DIR .eq. 'V') THEN
               DO I=1,M
                   BLR_PANEL(IP-CURRENT_BLR)%Q(I,1:N) =
     &             A( POSELT_BLOCK+int((I-1),8)*int(NFRONT,8) :
     &               POSELT_BLOCK+int((I-1),8)*int(NFRONT,8)
     &                           +int(N-1,8) )
               END DO  
           ELSE 
               DO I=1,N
                   BLR_PANEL(IP-CURRENT_BLR)%Q(1:M,I) =
     &             A( POSELT_BLOCK+int((I-1),8)*int(NFRONT,8) :
     &               POSELT_BLOCK+int((I-1),8)*int(NFRONT,8)
     &                           +int(M-1,8) )
               END DO  
           END IF
           IF (K473.EQ.0) THEN
               CALL UPDATE_FLOP_STATS_DEMOTE(BLR_PANEL(IP-CURRENT_BLR), 
     &                                   NIV)
           ENDIF
           BLR_PANEL(IP-CURRENT_BLR)%K = -1
        END IF
      END DO 
#if defined(BLR_MT) 
!$OMP END DO NOWAIT
#endif
      END SUBROUTINE DMUMPS_COMPRESS_PANEL
      END MODULE DMUMPS_FAC_LR
