module Cauchylowrank
  implicit none

contains

!!!!!!
    subroutine rrluCauchy(D,F,A,V,tol,U,W,PL,PU,M,N,Rk)
!
! .. Scalar parameter ..
    DOUBLE PRECISION  TOL
    INTEGER           M,N,Rk
! .. Array parameters ..
    DOUBLE PRECISION D(*),F(*),A(*),V(*),U(*),W(*)
    INTEGER          PL(*),PU(*)
!
! Purpose
! =========
! This routine computes a Rank Revealing Schur Complement, RRSC, for a Cauchy-like matrix
! via structured matrices. This matrix C has dimensions M-by-N, with generators D,F,A and V.
! C = (A(i)*V(j) / (D(i)-F(j)) )_{i,j}. 
!
! ============
! Written by Sheng-Guo Li, On Sept. 18th, 2012
! ============

! .. Local Scalars
    double precision rho, Amax, Ukk, zero, one, Nmax
    integer          nswap,j,k,mn,prd,flgl,ii,jj
    parameter        (zero = 0.0D0, one=1.0D0)

!  .. Intrinsic Functions ..
    intrinsic    max, abs, maxval,mod

    U(1:M)  = A(1:M)
    Rk = min(M,N)
    mn = Rk
    PL(1:M) = (/ (j, j=1,M) /)
    PU(1:N) = (/ (j, j=1,N) /)
    nswap = 0
    Amax = zero
    rho = 1.12D0
    prd = 10
    
    do k = 1, mn
       call CauchyPivtG(D,F,U,V,PL,PU,k,A,M,N)

       Ukk = u(k)*V(k)/( D(k)-F(k))
       Amax = max(Amax,abs(Ukk) )
       if (abs(Ukk) .lt. tol*Amax ) then  ! first step converged
          call CauchyPivtG_CP(D,F,U,V,PL,PU,k,A,M,N)
          Ukk = u(k)*V(k)/( D(k)-F(k))
          Amax = max(Amax,abs(Ukk) )
          if (abs(Ukk) .lt. tol*Amax ) then  ! final converged
             Rk = k -1
             exit
          end if
       end if
       
       U(k+1:M) = U(k+1:M)* (D(k+1:M)-D(k)) / (D(k+1:M)-F(k))
       V(k+1:N) = V(k+1:N)* (F(k+1:N)-F(k)) / (F(k+1:N)-D(k))
       do j = 1, k-1
          W(j) = W(j) * ( (F(k)-D(j))/ (D(k)-D(j)) )
       end do
       W(k) = (D(k)-F(k))/A(k)
       do j = 1, k-1
          W(k) = W(k) * (F(j)-D(k)) / ( D(j)-D(k) )
       end do

       ! swap
       flgl = mod(k,prd)
       do while(flgl .lt. 1)
          flgl = 1
          call searchMax2(U(k+1),W,D(k+1),D,M-k,k,ii,jj,Nmax)
          
          if(Nmax .gt. rho) then
             nswap = nswap + 1
             flgl = 0
             jj = jj + k 
             V(k+1:N)    = V(k+1:N) * ( (F(k+1:N)-D(ii)) / (F(k+1:N)-D(jj)) )
             U(k+1:jj-1) = U(k+1:jj-1) * ( (D(k+1:jj-1)-D(jj)) / (D(k+1:jj-1)-D(ii)) )
             U(jj+1:M)   = U(jj+1:M) * ( (D(jj+1:M)-D(jj)) / (D(jj+1:M)-D(ii)) )
             W(1:ii-1)   = W(1:ii-1) * ( (D(1:ii-1)-D(ii)) / (D(1:ii-1)-D(jj)) )
             W(ii+1:k)   = W(ii+1:k) * ( (D(ii+1:k)-D(ii)) / (D(ii+1:k)-D(jj)) )
             U(jj)       = A(ii) * ( (D(ii)-D(jj)) / (D(ii)-F(ii)) )  
             W(ii)       = (D(jj)-F(ii)) / A(jj)
             do j = 1, ii-1
                U(jj) = U(jj) * ( (D(ii)-D(j)) / (D(ii)-F(j)) )
                W(ii) = W(ii) * ( (F(j)-D(jj)) / (D(j)-D(jj)) )
             end do
             do j = ii+1, k
                U(jj) = U(jj) * ( (D(ii)-D(j)) / (D(ii)-F(j)) )
                W(ii) = W(ii) * ( (F(j)-D(jj)) / (D(j)-D(jj)) )
             end do
             call iswap(PL,ii,jj)
             call dswap(D,ii,jj)
             call dswap(A,ii,jj)
!             write(*,*) 'swap once '
          end if ! Nmax
          
       end do ! while

    end do ! main loop

  end subroutine rrluCauchy

!!!!!!
    subroutine searchMax2(U,W,D2,D1,LDU,LDW,ii,jj,Nmax)

    double precision U(*),W(*),D2(*),D1(*)
    integer LDU,LDW,ii,jj

    double precision zero,junk,Nmax
    integer j,jjL
    parameter(zero = 0.0D0)
    
    call CauchyMax(D2,D1(1),u,LDU,junk,jjL)
    junk = junk*abs(w(1))       
    Nmax = junk
    ii = 1
    jj = jjL
    
    do j = 2, LDW   ! ii: col, jj: row
       call CauchyMax(D2,D1(j),u,LDU,junk,jjL)
       junk = junk*abs(w(j))       
       if(junk .gt. Nmax) then
          Nmax = junk
          ii = j
          jj = jjL
       end if
    end do

  end subroutine searchMax2

!!!!!!
  subroutine CauchyMax(D,F,u,N,junk,jj)

    double precision D(*),u(*),junk,F
    integer N,jj
    integer temp(1)

!  .. Intrinsic Functions ..
    intrinsic    Maxloc, ABS,MAXVAL

    double precision, allocatable :: LK(:) ! change it to LK(N)

    allocate( LK(N) )

    Lk = u(1:N)/ ( D(1:N)-F )
    Lk = abs(Lk)
    junk = maxval(Lk(1:N))
    temp = maxloc(Lk(1:N))
    jj = temp(1)
    
    deallocate(LK)
    
  end subroutine CauchyMax

!!!!!!
  subroutine CauchyPivtG( D,F,u,v,PL,PU,k,A,M,N )
!
    integer k,M,N
    double precision D(*),F(*),u(*),v(*),A(*)
    integer PL(*), PU(*)

!  .. Intrinsic Functions ..
    intrinsic    ABS

!   local parameters 
    integer jjL,jjU,flg,jj
    double precision junkL,junkU, Piv, pivot,zero,one  
    parameter    (zero = 0.0D0, one=1.0D0)

    call CauchyMax(D(k),F(k),u(k),M-k+1,junkL,jjL)
    call CauchyMax(F(k),D(k),v(k),N-k+1,junkU,jjU)
    junkL = junkL * abs(v(k))
    junkU = junkU * abs(u(k))
    Piv   = abs(u(k)*v(k)/ (D(k)-F(k)) )
    
    if (junkL .le. Piv .and. junkU .le. Piv ) then
       return
    end if
    
    pivot = zero
    flg = 0
    if(junkL > junkU) flg = 1
    
    do while (1 < 2)
       pivot = pivot +1
       if (flg == 1) then
          jj = jjL
          call dswap(D,k,jj+k-1)
          call dswap(u,k,jj+k-1)
          call dswap(A,k,jj+k-1)
          call iswap(PL,k,jj+k-1)
          call CauchyMax( F(k),D(k),v(k),N-k+1,junkU,jjU ) ! N-k+1
          if(jjU == 1) return
          
          flg = 0
          continue
       end if
       jj = jjU
       call dswap(F,k,jj+k-1)
       call dswap(v,k,jj+k-1)
       call iswap(PU,k,jj+k-1)
       call CauchyMax( D(k),F(k),u(k),M-k+1,junkL,jjL )
       if(jjL == 1) return 

       flg = 1
    end do

  end subroutine CauchyPivtG

!!!!!!!!
    subroutine CauchyPivtG_CP( D,F,u,v,PL,PU,k,A,M,N )
!
!   This routine chooses the largest entry for the k-th Schur Complement.
!   It uses complete pivoting. If this routine could be faster than Matrix-Matrix
!   multiplication, we can try to modify the codes for bidiagonal SVD too. 
! 
    integer k,M,N
    double precision D(*),F(*),u(*),v(*),A(*)
    integer PL(*), PU(*)

!  .. Intrinsic Functions ..
    intrinsic    ABS
! ..
! .. local parameters 
    integer ii,jj,j
    double precision junk,zero,one
    parameter    (zero = 0.0D0, one=1.0D0)
!   ii : the row index; jj : the column index;
! ..
! .. local array ..
    integer iiU(N-k+1)    ! the row index for each column
    double precision junkU(N-k+1)  ! the largest entry for each column
!
!$OMP PARALLEL DO    
    DO j =1, N-k+1
       call CauchyMax(D(k),F(k+j-1),u(k),M-k+1,junkU(j),iiU(j) )
       junkU(j) = abs( junkU(j)*V(k+j-1) )
    END DO
!$OMP END PARALLEL DO

    junk = junkU(1)
    ii = iiU(1)
    jj = 1
    DO j = 2, N-k+1
       IF(junkU(j) > junk ) THEN
          junk = junkU(j)
          ii = iiU(j)
          jj = j
       END IF
    END DO
        
    call dswap(D,k,ii+k-1)
    call dswap(u,k,ii+k-1)
    call dswap(A,k,ii+k-1)
    call iswap(PL,k,ii+k-1)
    
    call dswap(F,k,jj+k-1)
    call dswap(v,k,jj+k-1)
    call iswap(PU,k,jj+k-1)
    
  end subroutine CauchyPivtG_CP

!!!!!!
  subroutine CauchyPivtG_VP( D,F,u,v,PL,PU,k,A,M,N )
!
    integer k,M,N
    double precision D(*),F(*),u(*),v(*),A(*)
    integer PL(*), PU(*)

!  .. Intrinsic Functions ..
    intrinsic    ABS, MIN

!   local parameters 
    integer flg,jj,MN,ii,j,ierr
    double precision junk, zero,one  
    parameter    (zero = 0.0D0, one=1.0D0)
!
!   local array
    integer, allocatable :: iiU(:)
    double precision, allocatable :: junkU(:)

    flg = 0
    MN = MIN(10, N-K+1)
    allocate( iiU(MN), junkU(MN), stat=ierr )
    if(ierr .ne. 0 ) then
       write(*,*) 'Allocate failed in PivtG_VP'
    end if

    do while ( flg .ne. 1 )
       call CauchyPivtG(D,F,u,v,PL,PU,k,A,M,N)

       ! maximun of first 10 columns
       DO j =1, MN
          call CauchyMax( D(k),F(k+j-1),u(k),M-k+1,junkU(j),iiU(j) )
          junkU(j) = abs( junkU(j)*V(k+j-1) ) 
       END DO

       junk = junkU(1)
       ii = iiU(1)
       jj = 1
       DO j = 2, MN
          IF(junkU(j) > junk ) THEN
             junk = junkU(j)
             ii = iiU(j)
             jj = j
          END IF
       END DO

       if( jj .eq. 1 ) then
          flg = 1
          return
       else
          call dswap(D,k,ii+k-1)
          call dswap(u,k,ii+k-1)
          call dswap(A,k,ii+k-1)
          call iswap(PL,k,ii+k-1)

          call dswap(F,k,jj+k-1)
          call dswap(v,k,jj+k-1)
          call iswap(PU,k,jj+k-1)
       end if

    end do

    deallocate(iiU, junkU)
  end subroutine CauchyPivtG_VP


!!!!!!
  subroutine dswap(D,mk,nk)
! double 1D array

    double precision D(*), tt
    integer mk,nk
    tt = D(mk)
    D(mk) = D(nk)
    D(nk) = tt
    
  end subroutine dswap

  subroutine iswap(D,mk,nk)
! integer 1D array

    integer D(*), tt, mk, nk
    tt = D(mk)
    D(mk) = D(nk)
    D(nk) = tt
    
  end subroutine iswap

!!!!!!
  subroutine compress_cauchy(rowcol,D,F,U,V,tol,M,N,H,Rk)
!
! Scalar parameters
    integer M,N,Rk
    double precision tol
    character(len=1) rowcol
! Array parameters
    double precision D(*),F(*),U(*),V(*),H(*)
!
! Purpose
! ========
! This routine computes a low rank approximation of some off-diagonal blocks of a
! Cauchy-like matrix. 
! 
! =========
! Written by S.G. Li, on Sept. 1st, 2019
!
! Note that this version only works for 'Row' compression case. 
! =========

! local scalars 
    integer mn,i,pnh
    logical  CR
    
! local arrays
    double precision zero, one
    parameter (zero = 0.0D0, one = 1.0D0)
    integer, allocatable :: PL(:), PU(:)
    double precision, allocatable :: Q(:,:),Z(:),W(:),LD(:),LF(:),LU(:),LV(:)

!  .. Intrinsic Functions ..
    intrinsic    Min
    
!  .. External Functions
      LOGICAL     LSAME
      EXTERNAL    LSAME
      
      pnh = 1
      mn = min(M,N)
      allocate(Z(M),W(mn),LD(M),LF(N),LU(M),LV(N),PL(M),PU(N) )
      LD(1:M) = D(1:M)
      LU(1:M) = U(1:M)
      LF(1:N) = F(1:N)
      LV(1:N) = V(1:N)

      Z = zero
      W = zero
      CR = lsame(rowcol,'r')
      if( CR )  then
         call rrluCauchy(LD,LF,LU,LV,tol,Z,W,PL,PU,M,N,Rk)
         call invp(PL,M)       ! PL --> InvPL
         allocate( Q(M,Rk) )
         if (Rk .lt. M) then            
            call Cauchylike2(Q,LD(Rk+1),LD,Z(Rk+1),W,M,Rk)    
            Q(1:M,1:Rk) = Q(PL,1:Rk)
            call dlacpy('A',M,Rk,Q,M,H(pnh),M)        ! copy tall matrix WH to H
            pnh = pnh+M*Rk
            call Cauchylike( H(pnh),LD,F,LU,V,Rk,N )    ! copy fat matrix QT to H            
         else
            ! copy identity matrix to generators
            WRITE(*,*) "Warning: This off-diagonal block is full rank."
         end if

      else ! block col
         
         call rrluCauchy(F,D,V,U,tol,Z,W,PL,PU,M,N,Rk)
         call invp(PL,M)  ! PL --> InvPL
         allocate( Q(Rk,M) )
         if (Rk .lt. M) then            
            call Cauchylike2(Q,F(Rk+1),F,Z(Rk+1),W,Rk,M)             
            Q(1:Rk,1:M) = Q(1:Rk,PL)
            call dlacpy('A',Rk,M,Q,Rk,H(pnh),Rk)      ! copy V to H
            !call hssexpmm2('V',PH,pnh,Rk,M,nodi,info) ! copy Q to generators            
         else
            ! copy identity matrix to generators
            Q(1:M,1:M) = zero
            do i = 1,M
               Q(i,i) = one
            end do
            Q(1:M,1:M) = Q(1:M,PL)            
            call dlacpy('A',M,M,Q,M,H(pnh),M)         ! copy V to H
            !call hssexpmm2('V',PH,pnh,Rk,M,nodi,info) ! copy Q to generators
         end if
         
      end if ! compr type

      deallocate(Z,W,Q,LD,LF,LU,LV,PL,PU)

    end subroutine compress_cauchy

!!!!!!
    subroutine invp(P,N)
      integer N,P(*)
      integer i
      integer :: IP(N)

      IP = 0

      do i = 1, N
         IP( P(i) ) = i
      end do
      P(1:N) = IP(1:N)

    end subroutine invp


subroutine Cauchylike(A,D,F,U,V,M,N)
!   A(i,j) = U(i)V(j) / (D(i)-F(j))

    integer M,N,i,j
    double precision A(*), D(*),F(*),u(*),v(*)

    do j = 1,N
       do i = 1,M
          A((j-1)*M+i) = U(i)*V(j) / (D(i)-F(j))
       end do
    end do

  end subroutine Cauchylike

!!!!!!
  subroutine Cauchylike2(A,D,F,u,v,M,N)
! A(i,j) = u(i)v(j) / (D(i)-F(j))

    integer M,N,i,j,Rk
    double precision A(M,*), D(*),F(*),u(*),v(*)
    double precision zero, one
    parameter(zero = 0.0D0, one = 1.0D0)
    
    Rk = M-N
    A(1:M,1:N) = zero
    
    if(Rk .eq. 0) then
       do j = 1, M
          A(j,j) = one
       end do
    else
       if(Rk .gt. 0) then
          do j = 1,N
             A(j,j) = one
             do i = N+1,M
                A(i,j) = u(i-N)*v(j) / (D(i-N)-F(j))
             end do
          end do
       else
          do j = 1, M
             A(j,j) = one
          end do
          do j = M+1, N
             do i = 1, M
                A(i,j) = v(i)*u(j-M) / (D(j-M)-F(i))
             end do
          end do
       end if
    end if
    
  end subroutine Cauchylike2
!
!!!!!!
       function comp_time()

         real(8) :: comp_time
         integer :: time0(8)
         call date_and_time(values=time0)
         comp_time = time0(5) * 3600 + time0(6)*60 + time0(7) +0.001 * time0(8)

       end function comp_time
!

end module Cauchylowrank
