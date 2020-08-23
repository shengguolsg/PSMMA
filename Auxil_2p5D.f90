
    SUBROUTINE ConstIndex( N, NB, MYCOL, NPCOL, CINdex, Length )
!
    IMPLICIT NONE
!
!  ---- ******* ---
    INTEGER  ::   N, NB, MYCOL, NPCOL, Length
    INTEGER  ::   CIndex(*)
!    INTEGER, INTENT(out), OPTIONAL :: Length
!
!   ========
! N     INTEGER (input)
!       The column dimension of matrix
!
! NB    INTEGER (input)
!       The block size for block cyclic distribution
!
! MYCOL INTEGER (input)
!        The column label of current process, it starts from 0
!
! NPCOL INTEGER (input)
!        The column size of process grid
!
! CIndex INTEGER (output)
!        The column index of current process contains.
!
! ==========
!  Written by Shengguo Li
!  It is written for Cauchy matrix multiplications.
! =====================================
!
    INTEGER :: I, J, KNB, TMP, NFSeg, Cstart, Gap, Len, Left
!
    INTEGER, EXTERNAL :: ICEIL

    !KNB = ICEIL(N,NB)    ! The total number of blocks or minus 1
    KNB = N/NB
    TMP = 1
    Length = 0

    NFSeg = KNB/NPCOL
    Cstart = MYCOL*NB+1
    Gap = NPCOL*NB

    DO I=1, NFSeg        ! Each segment has NB continous index
        CIndex(TMP:TMP+NB-1) = (/ (J, J=Cstart, Cstart+NB-1) /)
        Cstart = Cstart+Gap
        TMP = TMP+NB
        Length = TMP-1
    END DO

    ! Deal with the last segment, check whether is one more block or part
    Left = N-NFSeg*Gap-MYCOL*NB
    IF( Left .GT. 0 ) THEN
        Len = MIN(NB, Left)
        CIndex(TMP:TMP+Len-1) = (/ (J, J=Cstart, Cstart+Len-1) /)
        Length = TMP+Len-1
    END IF

!    IF( present(Length) )   THEN
!        Length= TMP+Len-1
!    END IF 

    RETURN

END SUBROUTINE ConstIndex

SUBROUTINE ConstCauchy( LM,LN,A,LDA,B,RIndex,CIndex )
!
        use cauchylowrank
        IMPLICIT NONE
!
!     HSSPACK
!     S.G. Li,  National University of Defense Technology
!     May 11th, 2020.  .
!
!     .. Scalar Arguments ..
        INTEGER            LM, LN, LDA
!     ..
!     .. Array Arguments ..
        DOUBLE PRECISION   A( LDA,* ), B( * )
        INTEGER ::         RIndex(*), CIndex(*)
!  ..
!  Purpose
!  =======
!
!  This routine constructs a Cauchy-like matrix locally with generators stored
!  in matrix A, which has four or five columns. B is used as a workspace to
!  store the compute Cauchy-like submatrix, and the indexes of row generators are 
!  stored in RIndex and the columns indexes are stored in CIndex. 
!      
!  Arguments
!  =========
!
!  M      - (input) INTEGER
!           The row dimension of locally matrix B, to be constructed.
!
!  N      - (input) INTEGER
!           The column dimension of locally matrix B, to be constructed.
!
!  A      - (input) DOUBLE PRECISION Array, with dimension (LDA, 4 or 5)
!           Stores the generators of a Cauchy-like matrix. A =[U V D W]
!
!  LDA    - (input) INTEGER
!           Leading dimension of A
!
!  B      - (output) DOUBLE PRECISION,
!           Stores the constructed Cauchy-like matrix with dimension (M, N).
!           B is defined as B_{ij} = U_i*V_j / (D_i -W_j), for i,j=1,...,M.
!
!           B is approximated by two low-rank matrices, B=X*Y, X is M-by-Rk,
!           and Y is Rk-by-N. X is stored in the first part of B, and Y is
!           stored in the second part. B is used as 1D array.
!
! RInd_start  - (input) INTEGER
!               The starting position for row generators
!
! CInd_Start  - (input) INTEGER
!               The starting position for column generators
!
!  Rk     - (output) INTEGER
!           It records the rank of this off-diagonal block.
!
! ======================================================================

!     ..
!     .. Local Scalars ..
!        DOUBLE PRECISION :: time, time1
!     ..
!     .. Local Arrays ..
        DOUBLE PRECISION, ALLOCATABLE :: D(:), F(:), U(:), V(:)

!     ..
!     .. Execution Parts ..

        ALLOCATE( D(LM),F(LN),U(LM),V(LN) )

        D(1:LM) = A(RIndex(1:LM),1)
        F(1:LN) = A(CIndex(1:LN),2)
        U(1:LM) = A(RIndex(1:LM),3)
        V(1:LN) = A(CIndex(1:LN),4)

        !call cpu_time(time)
        call Cauchylike( B,D,F,U,V,LM,LN )
        !call cpu_time(time1)
        !time = time1 - time
        !write(*,*) 'Construct Cauchylike costs', time, M, N

        DEALLOCATE( D,F,U,V )

END SUBROUTINE ConstCauchy


    SUBROUTINE ConstIndex1( N, NB, MYCOL, NPCOL, CINdex,PIndex )
!
    IMPLICIT NONE
!
!  ---- ******* ---
    INTEGER  ::   N, NB, MYCOL, NPCOL
    INTEGER  ::   CIndex(*),PIndex(*)
!
!   ========
! N     INTEGER (input)
!       The column dimension of matrix
!
! NB    INTEGER (input)
!       The block size for block cyclic distribution
!
! MYCOL INTEGER (input)
!        The column label of current process, it starts from 0
!
! NPCOL INTEGER (input)
!        The column size of process grid
!
! CIndex INTEGER (output)
!        The column index of current process contains.
!
! ==========
!  Written by Shengguo Li
!  It is written for Cauchy matrix multiplications.
! =====================================
!
    INTEGER :: I, KNB, TMP, NFSeg, Cstart, Gap
!
    INTEGER, EXTERNAL :: ICEIL

    KNB = ICEIL(N,NB)    ! The total number of blocks or minus 1
    TMP = 1

    NFSeg = KNB/NPCOL
    Cstart = MYCOL*NB+1
    Gap = NPCOL*NB

    DO I=1, NFSeg        ! Each segment has NB continous index
        CIndex(TMP:TMP+NB-1) = PIndex( Cstart:Cstart+NB-1 ) 
        Cstart = Cstart+Gap
        TMP = TMP+NB
    END DO

            
END SUBROUTINE ConstIndex1

! =====================================
!
SUBROUTINE is_intersect( M, SetA, N, SetB, is_inters )
!
! This function tests whether SetA and SetB have common elements.
! If they are intersected, return .TRUE., ELSE, return .FALSE. 
!
! SetA and SetB are 1D arrays, and their lengths are M and N, respectively.
! 
! ========================
 
    INTEGER   :: M, N
    LOGICAL   :: is_inters
!
    INTEGER   :: SetA(*), SetB(*)
!
! ..
! .. Local parameters ..
    INTEGER  :: MN, I

    MN = min(M,N)

    is_inters = .TRUE.

    IF( M > N ) THEN
        DO I=1, MN
            if( any(SetA(1:M) == SetB(I)) ) then
                is_inters = .FALSE.
                return 
            end if
        END DO  
    ELSE
        DO I=1, MN
            if( any(SetB(1:N) == SetA(I)) ) then
                is_inters = .FALSE.
                return 
            end if
        END DO
    END IF  

    RETURN 

END SUBROUTINE is_intersect

! ===================================
! 
    SUBROUTINE ConstCauchylowrank( M,N,A,LDA,H,RIndex,CIndex,Rk ) 
!
        use cauchylowrank
        IMPLICIT NONE
!
!     HSSPACK 
!     S.G. Li,  National University of Defense Technology
!     August 31th, 2019.  .
!
!     .. Scalar Arguments ..
        INTEGER            M, N, LDA, Rk
!     ..
!     .. Array Arguments ..
        DOUBLE PRECISION   A( LDA,* ), H( * )
        INTEGER            RIndex(*), CIndex(*)
!
!     ..
!
!  Purpose
!  =======
!
!  This routine constructs a Cauchy-like matrix locally with generators stored
!  in matrix A, which has four or five columns. B is used as a workspace to 
!  store the compute low-rank approximation of locally Cauchy-like submatrix.  
!
!  Arguments
!  =========
!
!
! ====================================================================== 

!     ..
!     .. Local Scalars ..
        DOUBLE PRECISION :: tol
        DOUBLE PRECISION :: time, time1
!     ..
!     .. Local Arrays ..
        DOUBLE PRECISION, ALLOCATABLE :: D(:), F(:), U(:), V(:)

!     ..
!     .. Execution Parts ..

        tol = 1.0E-17
! 
        ALLOCATE( D(M),F(N),U(M),V(N) )

        D(1:M) = A(RIndex(1:M),1) 
        F(1:N) = A(CIndex(1:N),2) 
        U(1:M) = A(RIndex(1:M),3) 
        V(1:N) = A(CIndex(1:N),4) 

        !call cpu_time(time) 
        call compress_cauchy( 'R',D,F,U,V,tol,M,N,H,Rk ) 
        !call cpu_time(time1) 
        !time = time1 - time
        !write(*,*) 'Construct low-rank approx. costs', time, M, N, Rk

        DEALLOCATE( D,F,U,V ) 

        END SUBROUTINE ConstCauchyLowrank
