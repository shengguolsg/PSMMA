!
      SUBROUTINE PSMMA_CAUCHY( N,K,ALPHA,A,IA,JA,DESCA,B,LDB, &
                               BETA,Q,IQ,JQ,DESCQ,WORK,REDIST )
!
       implicit none
    
      include 'mpif.h'
!
!  -- HSSPACK computational routine (version 1.0.0) --
!  -- HSSPACK is a software package provided by Nat. Univ. of Defense Tech. --
!     September 2019
!
!     .. Scalar Arguments ..
      INTEGER            N, K, IA, JA, IQ, JQ, LDB
      DOUBLE PRECISION   ALPHA, BETA
      LOGICAL            REDIST
!     ..
!     .. Array Arguments ..
      INTEGER          ::   DESCA(*), DESCQ( * )
      DOUBLE PRECISION ::   A(*), Q(*), WORK(*)
      DOUBLE PRECISION, INTENT(IN) :: B(LDB,4) 
!     ..
!
!  Purpose
!  =======
!
!  PSMMA_CAUCHY computes the multiplication of a matrix A with an HSS matrix, 
!  Q = alpha* A * H + beta * Q, where A is N-by-K, and H is K-by-K,
!  H(i,j) = u(i)*z(j) / ( d(i)-lambda(j) ).
!  It uses structured Cauchy-like parallel matrix-matrix multiplication  
!  algorithms, which implements the SPMMA algorithm.  
!  
!  The matrix A is already distributed in the block-cyclic form, and the 
!  generators of H are stored globally, all the entries can be
!  constructed locally. In this routine, we can redistribute A to use a larger
!  NB or use the original NB. We compare these two methods and test which one
!  is more efficient.     
!
!  PSMMA works for block-cyclic from. It nearly does not relate with NB except
!  that the off-diagonal ranks may be different. It prefers to using a larger NB
!  too. But the performance of PDGEMM depends more on NB.  
! 
!  We assume A is N-by-K, and H is K-by-K. 
! 
! ===========================
! Written by S.G. Li, NUDT, July. 2020
! =======================================================
!
!     .. Parameters ..
!
      INTEGER            BLOCK_CYCLIC_2D, DLEN_, DTYPE_, CTXT_, M_, N_, &
                         MB_, NB_, RSRC_, CSRC_, LLD_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,  &
                         CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,   &
                         RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            ICTXT, info, IAT, IQT, IWEND, &
                         LDAT,LDQT,MYCOL,MYROW,NB,NPCOL, NPROW,&
                         MBT, NBT, NQT, NPT, MB, LDA, LDQ, ierr
      REAL(8)            time1, time2
!     ..
!     .. Local Arrays ..
      INTEGER       ::   DESCAT(DLEN_),DESCQT(DLEN_)
!
!      DOUBLE PRECISION, ALLOCATABLE :: AT(:), QT(:)
!     ..
!     .. External Functions ..
      INTEGER            NUMROC
      EXTERNAL           NUMROC
!     ..
!     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, DESCINIT
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     This is just to keep ftnchek and toolpack/1 happy
      IF( BLOCK_CYCLIC_2D*CSRC_*CTXT_*DLEN_*DTYPE_*LLD_*MB_*M_*NB_*N_* &
          RSRC_.LT.0 ) RETURN
!
!     Test the input parameters.
!
      CALL BLACS_GRIDINFO( DESCQ( CTXT_ ), NPROW, NPCOL, MYROW, MYCOL )
!
      ICTXT = DESCQ( CTXT_ )
      LDA = DESCA( LLD_ )
      LDQ = DESCQ( LLD_ )
      MB    = DESCQ( MB_ )
      NB    = DESCQ( NB_ )

      ! Use the normal PSMMA algorithm without redistribution
      IF( .NOT. REDIST ) THEN
            time1 = MPI_Wtime()
            CALL PSCAUCHY_COMPUTE( N,K,ALPHA,A,LDA,DESCA,B, LDB,&
                  BETA, Q, LDQ, 0, 0, WORK )
            time2 = MPI_Wtime()
            IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
                write(*,*) 'The PSMMA without redistribution costs', time2-time1
            END IF     
            GOTO 90
      END IF
!
!     Construct the descriptor for the reshaped matrix A and Q
      MBT = N / NPROW 
      NBT = K / NPCOL
      IF( MBT*NPROW .NE. N )  MBT = MBT+1
      IF( NBT*NPCOL .NE. K )  NBT = NBT+1 
      NPT = NUMROC( N, MBT, MYROW, 0, NPROW )
      NQT = NUMROC( K, NBT, MYCOL, 0, NPCOL )

      !LDAT = MAX( NPT,MBT )   ! The locate leading dimension of AT
      LDAT = NPT
      LDQT = LDAT 
!
!     Here we assume N == K
      CALL DESCINIT( DESCAT, N, K, MBT, NBT, 0, 0, ICTXT, LDAT, &
                     INFO )
      CALL DESCINIT( DESCQT, N, K, MBT, NBT, 0, 0, ICTXT, LDQT, &
                     INFO )
!
!      ALLOCATE( AT(MBT*NBT), QT(MBT*NBT), stat=ierr )
!      IF( ierr .ne. 0 ) THEN
!         WRITE(*,*) 'Allocation failed in pdhssevc1, and return'
!         RETURN
!      END IF

      IAT = 1
      IQT = 1+ NPT*NBT
      IWEND= IQT+NPT*NBT ! IWEND must left NPT*NQT
       
      ! Redistribute matrix A to blocked form. 
      time1 = MPI_Wtime()
      CALL PDGEMR2D( N,K,A,1,1,DESCA,WORK(IAT),1,1,DESCAT,ICTXT ) 
      time2 = MPI_Wtime()
      IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
           write(*,*) 'The first redistribution costs', time2-time1, MBT, NBT
      END IF 

      ! Call PSMMA Algorithm for Cauchy-like matrices
      time1 = MPI_Wtime()
      CALL PSCAUCHY_COMPUTE( N,K, ALPHA,WORK(IAT),LDAT,DESCAT,B, LDB,&
                             BETA, WORK(IQT), LDQT, 0, 0, WORK(IWEND) ) 
      time2 = MPI_Wtime()
      IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
           write(*,*) 'The PSMMA with redistribution costs', time2-time1
      END IF 

!      CALL BLACS_BARRIER(ictxt,'All')
      ! Redistribute the matrix Q back
      time1 = MPI_Wtime()
      CALL PDGEMR2D( N,K,WORK(IQT),1,1,DESCQT,Q,IQ,JQ,DESCQ,ICTXT ) 
      time2 = MPI_Wtime()
      IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
           write(*,*) 'The second redistribution costs', time2-time1
      END IF
       
      ! Clean-up
!      DEALLOCATE( AT,QT )

  90    RETURN

      END SUBROUTINE PSMMA_CAUCHY
