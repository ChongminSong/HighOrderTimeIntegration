      MODULE MKLTimeHighOrderWrapper

      USE MKLpardiso_solver

      IMPLICIT NONE

      !      INCLUDE 'mkl.fi'
!      include 'mkl_spblas.fi'
      INTERFACE
      subroutine mkl_dcsrsymv( uplo, m,                                 &
     &      a, ia, ja,  x, y)
            character          uplo
            integer            m
            integer 		 ia(*), ja(*)
            double precision   a(*)
            double precision   y(*), x(*)
      END
      subroutine mkl_zcsrsymv( uplo, m,                                 &
     &a, ia, ja,  x, y)
      character          uplo
      integer            m
      integer 		 ia(*), ja(*)
      double complex   a(*)
      double complex   y(*), x(*)
      END

      SUBROUTINE DGEEV(JOBVL,JOBVR,N,A,LDA,WR,WI,VL,LDVL,VR,LDVR,WORK,  &
     &                 LWORK,INFO)
      CHARACTER          JOBVL,JOBVR
      INTEGER            INFO,LDA,LDVL,LDVR,LWORK,N
      DOUBLE PRECISION   A(LDA,*),VL(LDVL,*),VR(LDVR,*),WI(*),WORK(*),  &
     &                   WR(*)
      END

      END INTERFACE

      INTEGER, PRIVATE, PARAMETER :: dp = KIND(1.0D0)

      CONTAINS
!      
!     ***********************************************************
!
      SUBROUTINE TimeSolver01(aInfty, n, ia, ja, K, M, C, F,          &
     &                          dt, ft, pDOF, t, dsp, vel, acc, u0, v0)
 
      IMPLICIT NONE

      !..   Internal solver memory pointer 
      TYPE (MKLPardisoType) :: pM, pKd

      INTEGER, PARAMETER :: p = 1
      
      REAL(KIND=DP) :: aInfty
      INTEGER :: n, nnz, n2 
      INTEGER :: ia(:)
      INTEGER :: ja(:)
      REAL(KIND=DP) :: K(:),  M(:),  C(:)
      REAL(KIND=DP) :: F(:)
      
      REAL(KIND=DP) :: u0(:), v0(:)     

      INTEGER        :: npt, ns
      REAL(KIND=DP)  :: dt
      REAL(KIND=DP)  :: ft( :, : )
      INTEGER        :: npDOF
      INTEGER        :: pDOF( : )
      REAL(KIND=DP)  :: t( : ), dsp( :, : ), vel( :, : ), acc( :, : )
      
      REAL(KIND=DP), ALLOCATABLE :: Fb( : )

      REAL(KIND=DP), ALLOCATABLE :: f0( : ), Pf( :,: ), pt( : )
      REAL(KIND=DP), ALLOCATABLE :: b( : )
      REAL(KIND=DP), ALLOCATABLE :: x( :,: )
      
      REAL(KIND=DP), ALLOCATABLE :: pcoe( : )
      REAL(KIND=DP), ALLOCATABLE :: qcoe( : )
      REAL(KIND=DP), ALLOCATABLE :: cf( :,: )
      COMPLEX(KIND=DP) rs(1)
      REAL(KIND=DP) r
      
      REAL(KIND=DP), ALLOCATABLE ::  Kd( : )
      REAL(KIND=DP), ALLOCATABLE :: ftmp( : )

      INTEGER       :: i, it, is
      REAL(KIND=DP) :: fn
      REAL(4) WallClock, solverWallClock, solverTotalWallClock 
      
      n2 = n + n
      nnz = ia(n)
      ns = size(ft,1) + 1
      npt = min( size(ft,2), p+1 )
      npDOF = size(pDOF)
      
      ALLOCATE( pcoe(p+1) )
      ALLOCATE( qcoe(p+1) )
      ALLOCATE( cf(p+1,p) )
      call TimeIntgCoeff(p, aInfty, pcoe, qcoe, rs, cf)
      r = dreal(rs(1))
      
      call MKLpardiso_SetParameter()
     
      pM%mtype = 2
      pKd%mtype = 2

      WallClock = SECNDS(0.0)
      call MKLpardiso_factor_real(pM, n, ia, ja, M)
      WallClock = SECNDS(wallClock)
      write(*,*) '   Wall Clock time for factorization of M (sec) : ', WallClock

      allocate( Fb(n) ) 
      call MKLpardiso_solve_real(pM, n, ia, ja, M, F, Fb)
      
      allocate( f0(n2), Pf(n2, npt) , pt(npt) )
      f0 = 0.0_DP
      f0(1:n) = Fb
      Pf = 0.0_DP
      DO i = 1, npt
         Pf(:,i) = Pf(:,i) + cf(i,1)*f0
         pt(i) = 0.5_DP**(i-1)
      ENDDO
      
      ALLOCATE( b(n2), ftmp(n2) )
      ALLOCATE( x(n2,p+1) )
      x   = 0.0_DP
      
      ! initial conditions
      it = 1
      t(it) = 0.0_DP
      x(1:n,1)    = v0 
      x(n+1:n2,1) = u0 
      call MKLpardiso_Ax(pM, n, ia, ja, K, C, M, x(:,1), x(:,2))
      ! ft at t=0
      fn = ft(it,1)
      do i = 2, npt
         fn = fn + ft(it,i) * ( (-0.5_DP)**(i-1) )
      end do
      DO i = 1, npDOF
         dsp( it, i ) = u0( pDOF(i) )
         vel( it, i ) = v0( pDOF(i) )
         acc( it, i ) = x( pDOF(i),2 ) + fn*Fb( pDOF(i) ) 
      END DO
!      dsp(1,:) = 0.0_DP
!      vel(1,:) = 0.0_DP
!      acc(1,:) = 0.0_DP

           
      allocate( Kd(nnz) )
      Kd = (r*r)*M + r*C + K
      WallClock = SECNDS(0.0)
      call MKLpardiso_factor_real(pKd, n, ia, ja, Kd)
      WallClock = SECNDS(wallClock)
      write(*,*) '   Wall Clock time for factorization of Kd (sec) : ', WallClock
      
      WallClock = SECNDS(0.0)
      solverTotalWallClock = 0.0

      DO it = 2, ns
         is = it - 1
         
         b = x(:,1)*pcoe(1)
         do i = 2, p+1
            b = b + x(:,i)*pcoe(i)
         end do
         do i = 1, npt
            b = b + Pf(:,i)*ft(is,i)
         end do
         call mkl_dcsrsymv('U', n, M, ia, ja, b(1:n),    ftmp(1:n))
         call mkl_dcsrsymv('U', n, K, ia, ja, b(n+1:n2), ftmp(n+1:n2))
         ftmp(1:n) = r*ftmp(1:n) - ftmp(n+1:n2)
         
         solverWallClock = SECNDS(0.0)
         call MKLpardiso_solve_real(pKd, n, ia, ja, Kd, ftmp(1:n), ftmp(n+1:n2))
         solverWallClock = SECNDS(solverWallClock)
         solverTotalWallClock = solverTotalWallClock + solverWallClock

         x(1:n,1)    = ftmp(n+1:n2)
         x(n+1:n2,1) = ( ftmp(n+1:n2) +  b(n+1:n2) )/r
         
         x(:,2) = (b - qcoe(1)*x(:,1) )/qcoe(2)
         
         fn = 0.0_DP
         do i = 1, npt
            fn = fn + ft(is,i)*pt(i)
         end do
         
         t(it) = t(it-1) + dt
         DO i = 1, npDOF
            dsp( it, i ) = x( n+pDOF(i),1 )
            vel( it, i ) = x( pDOF(i),1 )
            acc( it, i ) = x( pDOF(i),2 ) + fn*Fb( pDOF(i) ) 
         END DO
         
      END DO
      
      write(*,*) '   Wall Clock time for back-substitution (sec) : ', solverTotalWallClock 
      
      WallClock = SECNDS(wallClock)
      write(*,*) '   Wall Clock time for stepping (sec) : ', WallClock

      vel = vel/dt
      acc = acc/(dt*dt)
      
      CALL MKLpardiso_delete(pM, n)
      CALL MKLpardiso_delete(pKd, n)
      
      DEALLOCATE( Fb )
      
      DEALLOCATE( f0, Pf, pt )
      DEALLOCATE( b )
      DEALLOCATE( x )
      
      DEALLOCATE( pcoe, qcoe, cf )
      DEALLOCATE( Kd, ftmp)

      
      RETURN
      END SUBROUTINE TimeSolver01
      
!      
!     ***********************************************************
!
      SUBROUTINE TimeSolver12(aInfty, n, ia, ja, K, M, C, F,          &
     &                          dt, ft, pDOF, t, dsp, vel, acc, u0, v0)
 
      IMPLICIT NONE

      !..   Internal solver memory pointer
      TYPE (MKLPardisoType) :: pM, pKd

      INTEGER, PARAMETER :: p = 2
      
      REAL(KIND=DP)  aInfty
      INTEGER :: n, nnz, n2 
      INTEGER :: ia(:)
      INTEGER :: ja(:)
      REAL(KIND=DP) :: K(:),  M(:),  C(:)
      REAL(KIND=DP) :: F(:)

      REAL(KIND=DP) :: u0(:), v0(:)     

      INTEGER        :: npt, ns
      REAL(KIND=DP)  :: dt
      REAL(KIND=DP)  :: ft( :, : )
      INTEGER        :: npDOF
      INTEGER        :: pDOF( : )
      REAL(KIND=DP)  :: t( : ), dsp( :, : ), vel( :, : ), acc( :, : )

      REAL(KIND=DP), ALLOCATABLE :: Fb( : )
      REAL(KIND=DP), ALLOCATABLE :: f0( :,: ), Pf( :,: ), pt( : )
      REAL(KIND=DP), ALLOCATABLE :: b( : )
      REAL(KIND=DP), ALLOCATABLE :: x( :,: )
      
      
      REAL(KIND=DP), ALLOCATABLE :: pcoe( : )
      REAL(KIND=DP), ALLOCATABLE :: qcoe( : )
      REAL(KIND=DP), ALLOCATABLE :: cf( :,: )
      COMPLEX(KIND=DP) r, rs(1)
      
      COMPLEX(KIND=DP), ALLOCATABLE ::  Kd( : )
      COMPLEX(KIND=DP), ALLOCATABLE ::  y( : )
      REAL(KIND=DP),    ALLOCATABLE :: ftmp( : )
      COMPLEX(KIND=DP), ALLOCATABLE :: ctmp( : )

      INTEGER       :: i, it, is
      REAL(KIND=DP) :: fn
      REAL(4) WallClock, solverWallClock, solverTotalWallClock 
      
      n2 = n + n
      nnz = ia(n)
      ns = size(ft,1) + 1
      npt = min( size(ft,2), p+1 )
      npDOF = size(pDOF)
      
      ALLOCATE( pcoe(p+1) )
      ALLOCATE( qcoe(p+1) )
      ALLOCATE( cf(p+1,p) )
      call TimeIntgCoeff(p, aInfty, pcoe, qcoe, rs, cf)
      r = rs(1)

      call MKLpardiso_SetParameter()
     
!.. Initialize the internal solver memory pointer. This is only
! necessary for the FIRST call of the PARDISO solver.

      pM%mtype = 2
      pKd%mtype = 6

      WallClock = SECNDS(0.0)
      call MKLpardiso_factor_real(pM, n, ia, ja, M)
      WallClock = SECNDS(wallClock)
      write(*,*) '   Wall Clock time for factorization of M (sec) : ', WallClock

      allocate( Fb(n) ) 
      call MKLpardiso_solve_real(pM, n, ia, ja, M, F, Fb)
      
      allocate( f0(n2,p), Pf(n2, npt) , pt(npt) )
      f0 = 0.0_DP
      f0(1:n,1) = Fb
      call  MKLpardiso_Ax(pM, n, ia, ja, K, C, M, f0(:,1), f0(:,2))
      Pf = 0.0_DP
      DO i = 1, npt
         Pf(:,i) = Pf(:,i) + (cf(i,1)*f0(:,1) + cf(i,2)*f0(:,2))
         pt(i) = 0.5_DP**(i-1)
      ENDDO
      
      ALLOCATE( x(n2,p+1) )
      x   = 0.0_DP
      
      ! initial conditions
      it = 1
      t(it) = 0.0_DP
      x(1:n,1)    = v0 
      x(n+1:n2,1) = u0 
      call MKLpardiso_Ax(pM, n, ia, ja, K, C, M, x(:,1), x(:,2))
      call MKLpardiso_Ax(pM, n, ia, ja, K, C, M, x(:,2), x(:,3))
      ! ft at t=0
      fn = ft(it,1)
      do i = 2, npt
         fn = fn + ft(it,i) * ( (-0.5_DP)**(i-1) )
      end do
      DO i = 1, npDOF
         dsp( it, i ) = u0( pDOF(i) )
         vel( it, i ) = v0( pDOF(i) )
         acc( it, i ) = x( pDOF(i),2 ) + fn*Fb( pDOF(i) ) 
      END DO
      
      ALLOCATE( b(n2), y(n2), ftmp(n2), ctmp(n2) )
      
      allocate( Kd(nnz) )
      Kd = (r*r)*M + r*C + K
      WallClock = SECNDS(0.0)
      call MKLpardiso_factor_complex(pKd, n, ia, ja, Kd)
      WallClock = SECNDS(wallClock)
      write(*,*) '   Wall Clock time for factorization of Kd (sec) : ', WallClock
      
      WallClock = SECNDS(0.0)
      solverTotalWallClock = 0.0

      DO it = 2, ns
         is = it - 1
         
         b = x(:,1)*pcoe(1)
         do i = 2, p+1
            b = b + x(:,i)*pcoe(i)
         end do
         do i = 1, npt
            b = b + Pf(:,i)*ft(is,i)
         end do
         call mkl_dcsrsymv('U', n, M, ia, ja, b(1:n),    ftmp(1:n))
         call mkl_dcsrsymv('U', n, K, ia, ja, b(n+1:n2), ftmp(n+1:n2))
         ctmp(1:n) = r*ftmp(1:n) - ftmp(n+1:n2)
         
         solverWallClock = SECNDS(0.0)
         call MKLpardiso_solve_complex(pKd, n, ia, ja, Kd, ctmp(1:n), ctmp(n+1:n2))

!         if (it.eq.ns/2) then
!         call mkl_zcsrsymv('U', n, Kd, ia, ja, ctmp(n+1:n2), y(1:n))
!         ftmp(1:n) = abs(ctmp(1:n)-y(1:n))
!         write(*,*) 'error: ', maxval(ftmp(1:n))/maxval(abs(ctmp(1:n))) 
!         endif

         
         solverWallClock = SECNDS(solverWallClock)
         solverTotalWallClock = solverTotalWallClock + solverWallClock
         y(1:n) = ctmp(n+1:n2)
         y(n+1:n2) = ( ctmp(n+1:n2) +  b(n+1:n2) )/r
         
         x(:,1) = -dimag(y)/dimag(r)
         x(:,2) = -dimag(r*y)/dimag(r)
         x(:,3) = (b - qcoe(1)*x(:,1) - qcoe(2)*x(:,2))/qcoe(3)
         
         fn = 0.0_DP
         do i = 1, npt
            fn = fn + ft(is,i)*pt(i)
         end do
         
         t(it) = t(it-1) + dt
         DO i = 1, npDOF
            dsp( it, i ) = x( n+pDOF(i),1 )
            vel( it, i ) = x( pDOF(i),1 )
            acc( it, i ) = x( pDOF(i),2 ) + fn*Fb( pDOF(i) ) 
         END DO
         
      END DO
      vel = vel/dt
      acc = acc/(dt*dt)

      write(*,*) '   Wall Clock time for back-substitution (sec) : ', solverTotalWallClock 
      
      WallClock = SECNDS(wallClock)
      write(*,*) '   Wall Clock time for stepping (sec) : ', WallClock
      
      CALL MKLpardiso_delete(pM, n)
      CALL MKLpardiso_delete(pKd, n)
      
      DEALLOCATE( Fb )

      DEALLOCATE( f0, Pf, pt )
      DEALLOCATE( b )
      DEALLOCATE( x )
      
      DEALLOCATE( pcoe, qcoe, cf )
      DEALLOCATE( Kd, y, ftmp, ctmp )

      
      RETURN
      END SUBROUTINE TimeSolver12
      
!      
!     ***********************************************************
!
      SUBROUTINE TimeSolver23(aInfty, n, ia, ja, K, M, C, F,          &
     &                          dt, ft, pDOF, t, dsp, vel, acc, u0, v0)
 
      IMPLICIT NONE

      !..   Internal solver memory pointer 
      TYPE (MKLPardisoType) :: pM, pKd1, pKd2

      INTEGER, PARAMETER :: p = 3
      
      REAL(KIND=DP) aInfty
      INTEGER :: n, nnz, n2 
      INTEGER :: ia(:)
      INTEGER :: ja(:)
      REAL(KIND=DP) :: K(:),  M(:),  C(:)
      REAL(KIND=DP) :: F(:)
      
      REAL(KIND=DP) :: u0(:), v0(:)     

      INTEGER        :: npt, ns
      REAL(KIND=DP)  :: dt
      REAL(KIND=DP)  :: ft( :, : )
      INTEGER        :: npDOF
      INTEGER        :: pDOF( : )
      REAL(KIND=DP)  :: t( : ), dsp( :, : ), vel( :, : ), acc( :, : )

      REAL(KIND=DP), ALLOCATABLE :: Fb( : )
      REAL(KIND=DP), ALLOCATABLE :: f0( :,: ), Pf( :,: ), pt( : )
      REAL(KIND=DP), ALLOCATABLE :: b( : )
      REAL(KIND=DP), ALLOCATABLE :: x( :,: )
      
      REAL(KIND=DP), ALLOCATABLE :: pcoe( : )
      REAL(KIND=DP), ALLOCATABLE :: qcoe( : )
      REAL(KIND=DP), ALLOCATABLE :: cf( :,: )
      COMPLEX(KIND=DP) rs(2)
      REAL(KIND=DP) r1
      COMPLEX(KIND=DP) r2
      
      REAL(KIND=DP),    ALLOCATABLE ::  Kd1( : )
      REAL(KIND=DP),    ALLOCATABLE ::  x1( : )
      COMPLEX(KIND=DP), ALLOCATABLE ::  Kd2( : )
      COMPLEX(KIND=DP), ALLOCATABLE ::  y( : )
      REAL(KIND=DP),    ALLOCATABLE :: ftmp( : )
      COMPLEX(KIND=DP), ALLOCATABLE :: ctmp( : )

      INTEGER       :: i, it, is
      REAL(KIND=DP) :: fn
      REAL(4) WallClock, solverWallClock, solverTotalWallClock 
      
      n2 = n + n
      nnz = ia(n)
      ns = size(ft,1) + 1
      npt = min( size(ft,2), p+1 )
      npDOF = size(pDOF)
      
      ALLOCATE( pcoe(p+1) )
      ALLOCATE( qcoe(p+1) )
      ALLOCATE( cf(p+1,p) )
      call TimeIntgCoeff(p, aInfty, pcoe, qcoe, rs, cf)
      r1 = dreal(rs(1))
      r2 = rs(2)
      
      call MKLpardiso_SetParameter()
     
!.. Initialize the internal solver memory pointer. This is only
! necessary for the FIRST call of the PARDISO solver.

      pM%mtype = 2
      pKd1%mtype = 2
      pKd2%mtype = 6

      WallClock = SECNDS(0.0)
      call MKLpardiso_factor_real(pM, n, ia, ja, M)
      WallClock = SECNDS(wallClock)
      write(*,*) '   Wall Clock time for factorization of M (sec) : ', WallClock

      allocate( Fb(n) ) 
      call MKLpardiso_solve_real(pM, n, ia, ja, M, F, Fb)
      
      allocate( f0(n2,p), Pf(n2, npt) , pt(npt) )
      f0 = 0.0_DP
      f0(1:n,1) = Fb
      call  MKLpardiso_Ax(pM, n, ia, ja, K, C, M, f0(:,1), f0(:,2))
      call  MKLpardiso_Ax(pM, n, ia, ja, K, C, M, f0(:,2), f0(:,3))
      Pf = 0.0_DP
      DO i = 1, npt
         Pf(:,i) = Pf(:,i) + (cf(i,1)*f0(:,1) + cf(i,2)*f0(:,2) + cf(i,3)*f0(:,3))
         pt(i) = 0.5_DP**(i-1)
      ENDDO
      
      ALLOCATE( x(n2,p+1) )
      x   = 0.0_DP
      
      ! initial conditions
      it = 1
      t(it) = 0.0_DP
      x(1:n,1)    = v0 
      x(n+1:n2,1) = u0 
      do i = 1, p
         call MKLpardiso_Ax(pM, n, ia, ja, K, C, M, x(:,i), x(:,i+1))
      enddo
      ! ft at t=0
      fn = ft(it,1)
      do i = 2, npt
         fn = fn + ft(it,i) * ( (-0.5_DP)**(i-1) )
      end do
      DO i = 1, npDOF
         dsp( it, i ) = u0( pDOF(i) )
         vel( it, i ) = v0( pDOF(i) )
         acc( it, i ) = x( pDOF(i),2 ) + fn*Fb( pDOF(i) ) 
      END DO

      ALLOCATE( b(n2), x1(n2), y(n2), ftmp(n2), ctmp(n2) )

      allocate( Kd1(nnz), Kd2(nnz) )
      Kd1 = (r1*r1)*M + r1*C + K
      WallClock = SECNDS(0.0)
      call MKLpardiso_factor_real(pKd1, n, ia, ja, Kd1)
      WallClock = SECNDS(wallClock)
      write(*,*) '   Wall Clock time for factorization of Kd1 (sec) : ', WallClock
      Kd2 = (r2*r2)*M + r2*C + K
      WallClock = SECNDS(0.0)
      call MKLpardiso_factor_complex(pKd2, n, ia, ja, Kd2)
      WallClock = SECNDS(wallClock)
      write(*,*) '   Wall Clock time for factorization of Kd2 (sec) : ', WallClock
      
      WallClock = SECNDS(0.0)
      solverTotalWallClock = 0.0
      
      DO it = 2, ns
         is = it - 1
         
         b = x(:,1)*pcoe(1)
         do i = 2, p+1
            b = b + x(:,i)*pcoe(i)
         end do
         do i = 1, npt
            b = b + Pf(:,i)*ft(is,i)
         end do
         call mkl_dcsrsymv('U', n, M, ia, ja, b(1:n),    ftmp(1:n))
         call mkl_dcsrsymv('U', n, K, ia, ja, b(n+1:n2), ftmp(n+1:n2))
         ftmp(1:n) = r1*ftmp(1:n) - ftmp(n+1:n2)
         solverWallClock = SECNDS(0.0)
         call MKLpardiso_solve_real(pKd1, n, ia, ja, Kd1, ftmp(1:n), ftmp(n+1:n2))
         solverWallClock = SECNDS(solverWallClock)
         solverTotalWallClock = solverTotalWallClock + solverWallClock
        x1(1:n)    = ftmp(n+1:n2)
         x1(n+1:n2) = ( ftmp(n+1:n2) +  b(n+1:n2) )/r1
         
         call mkl_dcsrsymv('U', n, M, ia, ja, x1(1:n),    ftmp(1:n))
         call mkl_dcsrsymv('U', n, K, ia, ja, x1(n+1:n2), ftmp(n+1:n2))
         ctmp(1:n) = r2*ftmp(1:n) - ftmp(n+1:n2)
         
         solverWallClock = SECNDS(0.0)
         call MKLpardiso_solve_complex(pKd2, n, ia, ja, Kd2, ctmp(1:n), ctmp(n+1:n2))
         solverWallClock = SECNDS(solverWallClock)
         solverTotalWallClock = solverTotalWallClock + solverWallClock
         y(1:n) = ctmp(n+1:n2)
         y(n+1:n2) = ( ctmp(n+1:n2) +  x1(n+1:n2) )/r2
         
         x(:,1) = -dimag(y)/dimag(r2)
         x(:,2) = -dimag(r2*y)/dimag(r2)
         x(:,3) = x1 - dreal(r2*dconjg(r2))*x(:,1) + 2*dreal(r2)*x(:,2)
         x(:,4) = ( b - qcoe(1)*x(:,1) - qcoe(2)*x(:,2) - qcoe(3)*x(:,3) )/qcoe(4)
         
         fn = 0.0_DP
         do i = 1, npt
            fn = fn + ft(is,i)*pt(i)
         end do
         
         t(it) = t(it-1) + dt
         DO i = 1, npDOF
            dsp( it, i ) = x( n+pDOF(i),1 )
            vel( it, i ) = x( pDOF(i),1 )
            acc( it, i ) = x( pDOF(i),2 ) + fn*Fb( pDOF(i) ) 
         END DO
         
      END DO
      vel = vel/dt
      acc = acc/(dt*dt)

      write(*,*) '   Wall Clock time for back-substitution (sec) : ', solverTotalWallClock 
      
      WallClock = SECNDS(wallClock)
      write(*,*) '   Wall Clock time for stepping (sec) : ', WallClock
      
      CALL MKLpardiso_delete(pM, n)
      CALL MKLpardiso_delete(pKd1, n)
      CALL MKLpardiso_delete(pKd2, n)
      
      DEALLOCATE( Fb )

      DEALLOCATE( f0, Pf, pt )
      DEALLOCATE( b )
      DEALLOCATE( x )
      
      DEALLOCATE( pcoe, qcoe, cf )
      DEALLOCATE( Kd1, Kd2, x1, y, ftmp, ctmp )
      
      RETURN
      END SUBROUTINE TimeSolver23
!      
!     ***********************************************************
!
      SUBROUTINE TimeSolver34(aInfty, n, ia, ja, K, M, C, F,          &
     &                          dt, ft, pDOF, t, dsp, vel, acc, u0, v0)
 
      IMPLICIT NONE

      !..   Internal solver memory pointer 
      TYPE (MKLPardisoType) :: pM, pKd1, pKd2

      INTEGER, PARAMETER :: p = 4
      
      REAL(KIND=DP) :: aInfty
      INTEGER :: n, nnz, n2 
      INTEGER :: ia(:)
      INTEGER :: ja(:)
      REAL(KIND=DP) :: K(:),  M(:),  C(:)
      REAL(KIND=DP) :: F(:)
      
      REAL(KIND=DP) :: u0(:), v0(:)     

      INTEGER        :: npt, ns
      REAL(KIND=DP)  :: dt
      REAL(KIND=DP)  :: ft( :, : )
      INTEGER        :: npDOF
      INTEGER        :: pDOF( : )
      REAL(KIND=DP)  :: t( : ), dsp( :, : ), vel( :, : ), acc( :, : )

      REAL(KIND=DP), ALLOCATABLE :: Fb( : )
      REAL(KIND=DP), ALLOCATABLE :: f0( :,: ), Pf( :,: ), pt( : )
      REAL(KIND=DP), ALLOCATABLE :: b( : )
      REAL(KIND=DP), ALLOCATABLE :: x( :,: )
      
      REAL(KIND=DP), ALLOCATABLE :: pcoe( : )
      REAL(KIND=DP), ALLOCATABLE :: qcoe( : )
      REAL(KIND=DP), ALLOCATABLE :: cf( :,: )
      COMPLEX(KIND=DP) rs(2)
      COMPLEX(KIND=DP) r1, r2
      
      COMPLEX(KIND=DP), ALLOCATABLE ::  Kd1( : )
      REAL(KIND=DP), ALLOCATABLE    ::  x1( : )
      COMPLEX(KIND=DP), ALLOCATABLE ::  y1( : )
      COMPLEX(KIND=DP), ALLOCATABLE ::  Kd2( : )
      COMPLEX(KIND=DP), ALLOCATABLE ::  y2( : )
      REAL(KIND=DP),    ALLOCATABLE :: ftmp( : )
      COMPLEX(KIND=DP), ALLOCATABLE :: ctmp( : )

      INTEGER       :: i, j, it, is
      REAL(KIND=DP) :: fn
      REAL(4) WallClock, solverWallClock, solverTotalWallClock 
      
      n2 = n + n
      nnz = ia(n)
      ns = size(ft,1) + 1
      npt = min( size(ft,2), p+1 )
      npDOF = size(pDOF)
      
      ALLOCATE( pcoe(p+1) )
      ALLOCATE( qcoe(p+1) )
      ALLOCATE( cf(p+1,p) )

      call TimeIntgCoeff(p, aInfty, pcoe, qcoe, rs, cf)
      r1 = rs(1)
      r2 = rs(2)

      call MKLpardiso_SetParameter()
     
!.. Initialize the internal solver memory pointer. This is only
! necessary for the FIRST call of the PARDISO solver.

      pM%mtype = 2
      pKd1%mtype = 6
      pKd2%mtype = 6

      WallClock = SECNDS(0.0)
      call MKLpardiso_factor_real(pM, n, ia, ja, M)
      WallClock = SECNDS(wallClock)
      write(*,*) '   Wall Clock time for factorization of M (sec) : ', WallClock

      allocate( Fb(n) ) 
      call MKLpardiso_solve_real(pM, n, ia, ja, M, F, Fb)
      
      allocate( f0(n2,p), Pf(n2, npt) , pt(npt) )
      f0 = 0.0_DP
      f0(1:n,1) = Fb
      DO i = 2, p
         call  MKLpardiso_Ax(pM, n, ia, ja, K, C, M, f0(:,i-1), f0(:,i))
      ENDDO
      Pf = 0.0_DP
      DO i = 1, npt
         DO j = 1, p
            Pf(:,i) = Pf(:,i) + cf(i,j)*f0(:,j)
         ENDDO
         pt(i) = 0.5_DP**(i-1)
      ENDDO
      
      ALLOCATE( x(n2,p+1) )
      x   = 0.0_DP
      
      ! initial conditions
      it = 1
      t(it) = 0.0_DP
      x(1:n,1)    = v0 
      x(n+1:n2,1) = u0 
      do i = 1, p
         call MKLpardiso_Ax(pM, n, ia, ja, K, C, M, x(:,i), x(:,i+1))
      enddo
      ! ft at t=0
      fn = ft(it,1)
      do i = 2, npt
         fn = fn + ft(it,i) * ( (-0.5_DP)**(i-1) )
      end do
      DO i = 1, npDOF
         dsp( it, i ) = u0( pDOF(i) )
         vel( it, i ) = v0( pDOF(i) )
         acc( it, i ) = x( pDOF(i),2 ) + fn*Fb( pDOF(i) ) 
      END DO

      ALLOCATE( b(n2), x1(n2), y1(n2), y2(n2), ftmp(n2), ctmp(n2) )
      
      allocate( Kd1(nnz), Kd2(nnz) )
!      Kd1 = (r1*r1)*M(1:nnz) + r1*C(1:nnz) + K(1:nnz)
      Kd1 = (r1*r1)*M + r1*C + K
      WallClock = SECNDS(0.0)
      call MKLpardiso_factor_complex(pKd1, n, ia, ja, Kd1)
      WallClock = SECNDS(wallClock)
      write(*,*) '   Wall Clock time for factorization of Kd1 (sec) : ', WallClock
      Kd2 = (r2*r2)*M + r2*C + K
      WallClock = SECNDS(0.0)
      call MKLpardiso_factor_complex(pKd2, n, ia, ja, Kd2)
      WallClock = SECNDS(wallClock)
      write(*,*) '   Wall Clock time for factorization of Kd2 (sec) : ', WallClock
      
      WallClock = SECNDS(0.0)
      solverTotalWallClock = 0.0
      
      DO it = 2, ns
         is = it - 1
         
         b = x(:,1)*pcoe(1)
         do i = 2, p+1
            b = b + x(:,i)*pcoe(i)
         end do
         do i = 1, npt
            b = b + Pf(:,i)*ft(is,i)
         end do
         call mkl_dcsrsymv('U', n, M, ia, ja, b(1:n),    ftmp(1:n))
         call mkl_dcsrsymv('U', n, K, ia, ja, b(n+1:n2), ftmp(n+1:n2))
         ctmp(1:n) = r1*ftmp(1:n) - ftmp(n+1:n2)
         solverWallClock = SECNDS(0.0)
         call MKLpardiso_solve_complex(pKd1, n, ia, ja, Kd1, ctmp(1:n), ctmp(n+1:n2))
         solverWallClock = SECNDS(solverWallClock)
         solverTotalWallClock = solverTotalWallClock + solverWallClock
         y1(1:n)    = ctmp(n+1:n2)
         y1(n+1:n2) = ( ctmp(n+1:n2) +  b(n+1:n2) )/r1
         x1 =  -dimag(y1)/dimag(r1);

         
         call mkl_dcsrsymv('U', n, M, ia, ja, x1(1:n),    ftmp(1:n))
         call mkl_dcsrsymv('U', n, K, ia, ja, x1(n+1:n2), ftmp(n+1:n2))
         ctmp(1:n) = r2*ftmp(1:n) - ftmp(n+1:n2)
         solverWallClock = SECNDS(0.0)
         call MKLpardiso_solve_complex(pKd2, n, ia, ja, Kd2, ctmp(1:n), ctmp(n+1:n2))
         solverWallClock = SECNDS(solverWallClock)
         solverTotalWallClock = solverTotalWallClock + solverWallClock
         y2(1:n) = ctmp(n+1:n2)
         y2(n+1:n2) = ( ctmp(n+1:n2) +  x1(n+1:n2) )/r2
         
         x(:,1) = -dimag(y2)/dimag(r2)
         x(:,2) = -dimag(r2*y2)/dimag(r2)
         x(:,3) = x1 - dreal(r2*dconjg(r2))*x(:,1) + 2.0_DP*dreal(r2)*x(:,2)
         x(:,4) = dreal(dconjg(r1)*x1 - y1 - (r2*dconjg(r2)*x(:,2)     &
      &                 - 2.0_DP*dreal(r2)*x(:,3)));
         x(:,5) = b
         DO j = 1, p
            x(:,5) = x(:,5) -  qcoe(j)*x(:,j)
         END DO
         x(:,5) = x(:,5)/qcoe(p+1)
         
         fn = 0.0_DP
         do i = 1, npt
            fn = fn + ft(is,i)*pt(i)
         end do
         
         t(it) = t(it-1) + dt
         DO i = 1, npDOF
            dsp( it, i ) = x( n+pDOF(i),1 )
            vel( it, i ) = x( pDOF(i),1 )
            acc( it, i ) = x( pDOF(i),2 ) + fn*Fb( pDOF(i) ) 
         END DO
         
      END DO
      vel = vel/dt
      acc = acc/(dt*dt)

      write(*,*) '   Wall Clock time for back-substitution (sec) : ', solverTotalWallClock 
      
      WallClock = SECNDS(wallClock)
      write(*,*) '   Wall Clock time for stepping (sec) : ', WallClock
      
      CALL MKLpardiso_delete(pM, n)
      CALL MKLpardiso_delete(pKd1, n)
      CALL MKLpardiso_delete(pKd2, n)
      
      DEALLOCATE( Fb )

      DEALLOCATE( f0, Pf, pt )
      DEALLOCATE( b )
      DEALLOCATE( x )
      
      DEALLOCATE( pcoe, qcoe, cf )
      DEALLOCATE( Kd1, Kd2, x1, y1, y2, ftmp, ctmp )

      
      RETURN
      END SUBROUTINE TimeSolver34
!      
!     ***********************************************************
!
      SUBROUTINE TimeSolver45(aInfty, n, ia, ja, K, M, C, F,          &
     &                          dt, ft, pDOF, t, dsp, vel, acc, u0, v0)
 
      IMPLICIT NONE

      !..   Internal solver memory pointer 
      TYPE (MKLPardisoType) :: pM, pKd1, pKd2, pKd3

      INTEGER, PARAMETER :: p = 5
      

      REAL(KIND=DP) aInfty
      INTEGER :: n, nnz, n2 
      INTEGER :: ia(:)
      INTEGER :: ja(:)
      REAL(KIND=DP) :: K(:),  M(:),  C(:)
      REAL(KIND=DP) :: F(:)
      
      REAL(KIND=DP) :: u0(:), v0(:)     

      INTEGER        :: npt, ns
      REAL(KIND=DP)  :: dt
      REAL(KIND=DP)  :: ft( :, : )
      INTEGER        :: npDOF
      INTEGER        :: pDOF( : )
      REAL(KIND=DP)  :: t( : ), dsp( :, : ), vel( :, : ), acc( :, : )

      REAL(KIND=DP), ALLOCATABLE :: Fb( : )
      REAL(KIND=DP), ALLOCATABLE :: f0( :,: ), Pf( :,: ), pt( : )
      REAL(KIND=DP), ALLOCATABLE :: b( : )
      REAL(KIND=DP), ALLOCATABLE :: x( :,: )
      
      REAL(KIND=DP), ALLOCATABLE :: pcoe( : )
      REAL(KIND=DP), ALLOCATABLE :: qcoe( : )
      REAL(KIND=DP), ALLOCATABLE :: cf( :,: )
      COMPLEX(KIND=DP) rs(3)
      REAL(KIND=DP)    r1
      COMPLEX(KIND=DP) r2, r3
      
      REAL(KIND=DP), ALLOCATABLE    ::  Kd1( : )
      REAL(KIND=DP), ALLOCATABLE    ::  x1( : )
      COMPLEX(KIND=DP), ALLOCATABLE ::  y1( : )
      COMPLEX(KIND=DP), ALLOCATABLE ::  Kd2( : )
      REAL(KIND=DP), ALLOCATABLE    ::  x2( : )
      COMPLEX(KIND=DP), ALLOCATABLE ::  y2( : )
      COMPLEX(KIND=DP), ALLOCATABLE ::  Kd3( : )
      COMPLEX(KIND=DP), ALLOCATABLE ::  y3( : )
      REAL(KIND=DP),    ALLOCATABLE :: ftmp( : )
      COMPLEX(KIND=DP), ALLOCATABLE :: ctmp( : )

      INTEGER       :: i, j, it, is
      REAL(KIND=DP) :: fn
      REAL(4) WallClock, solverWallClock, solverTotalWallClock 
      
      n2 = n + n
      nnz = ia(n)
      ns = size(ft,1) + 1
      npt = min( size(ft,2), p+1 )
      npDOF = size(pDOF)
      
      ALLOCATE( pcoe(p+1) )
      ALLOCATE( qcoe(p+1) )
      ALLOCATE( cf(p+1,p) )
      call TimeIntgCoeff(p, aInfty, pcoe, qcoe, rs, cf)
      r1 = dreal( rs(1) )
      r2 = rs(2)
      r3 = rs(3)
      
      call MKLpardiso_SetParameter()
     
!.. Initialize the internal solver memory pointer. This is only
! necessary for the FIRST call of the PARDISO solver.

      pM%mtype = 2
      pKd1%mtype = 2
      pKd2%mtype = 6
      pKd3%mtype = 6

      WallClock = SECNDS(0.0)
      call MKLpardiso_factor_real(pM, n, ia, ja, M)
      WallClock = SECNDS(wallClock)
      write(*,*) '   Wall Clock time for factorization of M (sec) : ', WallClock

      allocate( Fb(n) ) 
      call MKLpardiso_solve_real(pM, n, ia, ja, M, F, Fb)
      
      allocate( f0(n2,p), Pf(n2, npt) , pt(npt) )
      f0 = 0.0_DP
      f0(1:n,1) = Fb
      DO i = 2, p
         call  MKLpardiso_Ax(pM, n, ia, ja, K, C, M, f0(:,i-1), f0(:,i))
      ENDDO
      Pf = 0.0_DP
      DO i = 1, npt
         DO j = 1, p
            Pf(:,i) = Pf(:,i) + cf(i,j)*f0(:,j)
         ENDDO
         pt(i) = 0.5_DP**(i-1)
      ENDDO
      
      ALLOCATE( x(n2,p+1) )
      x   = 0.0_DP
      
      ! initial conditions
      it = 1
      t(it) = 0.0_DP
      x(1:n,1)    = v0 
      x(n+1:n2,1) = u0 
      do i = 1, p
         call MKLpardiso_Ax(pM, n, ia, ja, K, C, M, x(:,i), x(:,i+1))
      enddo
      ! ft at t=0
      fn = ft(it,1)
      do i = 2, npt
         fn = fn + ft(it,i) * ( (-0.5_DP)**(i-1) )
      end do
      DO i = 1, npDOF
         dsp( it, i ) = u0( pDOF(i) )
         vel( it, i ) = v0( pDOF(i) )
         acc( it, i ) = x( pDOF(i),2 ) + fn*Fb( pDOF(i) ) 
      END DO

      ALLOCATE( b(n2), x1(n2), y1(n2), x2(n2), y2(n2), y3(n2), ftmp(n2), ctmp(n2) )
      
      allocate( Kd1(nnz), Kd2(nnz), Kd3(nnz) )
      Kd1 = (r1*r1)*M + r1*C + K
      WallClock = SECNDS(0.0)
      call MKLpardiso_factor_real(pKd1, n, ia, ja, Kd1)
      WallClock = SECNDS(wallClock)
      write(*,*) '   Wall Clock time for factorization of Kd1 (sec) : ', WallClock
      Kd2 = (r2*r2)*M + r2*C + K
      WallClock = SECNDS(0.0)
      call MKLpardiso_factor_complex(pKd2, n, ia, ja, Kd2)
      WallClock = SECNDS(wallClock)
      write(*,*) '   Wall Clock time for factorization of Kd2 (sec) : ', WallClock
      WallClock = SECNDS(0.0)
      Kd3 = (r3*r3)*M + r3*C + K
      call MKLpardiso_factor_complex(pKd3, n, ia, ja, Kd3)
      WallClock = SECNDS(wallClock)
      write(*,*) '   Wall Clock time for factorization of Kd3 (sec) : ', WallClock
      
      WallClock = SECNDS(0.0)
      solverTotalWallClock = 0.0
      
      DO it = 2, ns
         is = it - 1
         
         b = x(:,1)*pcoe(1)
         do i = 2, p+1
            b = b + x(:,i)*pcoe(i)
         end do
         do i = 1, npt
            b = b + Pf(:,i)*ft(is,i)
         end do
         call mkl_dcsrsymv('U', n, M, ia, ja, b(1:n),    ftmp(1:n))
         call mkl_dcsrsymv('U', n, K, ia, ja, b(n+1:n2), ftmp(n+1:n2))
         ftmp(1:n) = r1*ftmp(1:n) - ftmp(n+1:n2)
         solverWallClock = SECNDS(0.0)
         call MKLpardiso_solve_real(pKd1, n, ia, ja, Kd1, ftmp(1:n), ftmp(n+1:n2))
         solverWallClock = SECNDS(solverWallClock)
         solverTotalWallClock = solverTotalWallClock + solverWallClock
         x1(1:n)    = ftmp(n+1:n2)
         x1(n+1:n2) = ( ftmp(n+1:n2) +  b(n+1:n2) )/r1
         
         call mkl_dcsrsymv('U', n, M, ia, ja, x1(1:n),    ftmp(1:n))
         call mkl_dcsrsymv('U', n, K, ia, ja, x1(n+1:n2), ftmp(n+1:n2))
         ctmp(1:n) = r2*ftmp(1:n) - ftmp(n+1:n2)
         solverWallClock = SECNDS(0.0)
         call MKLpardiso_solve_complex(pKd2, n, ia, ja, Kd2, ctmp(1:n), ctmp(n+1:n2))
         solverWallClock = SECNDS(solverWallClock)
         solverTotalWallClock = solverTotalWallClock + solverWallClock
         y2(1:n) = ctmp(n+1:n2)
         y2(n+1:n2) = ( ctmp(n+1:n2) +  x1(n+1:n2) )/r2
         x2 =  -dimag(y2)/dimag(r2)
         
         call mkl_dcsrsymv('U', n, M, ia, ja, x2(1:n),    ftmp(1:n))
         call mkl_dcsrsymv('U', n, K, ia, ja, x2(n+1:n2), ftmp(n+1:n2))
         ctmp(1:n) = r3*ftmp(1:n) - ftmp(n+1:n2)
         solverWallClock = SECNDS(0.0)
         call MKLpardiso_solve_complex(pKd3, n, ia, ja, Kd3, ctmp(1:n), ctmp(n+1:n2))
         solverWallClock = SECNDS(solverWallClock)
         solverTotalWallClock = solverTotalWallClock + solverWallClock
         y3(1:n)    = ctmp(n+1:n2)
         y3(n+1:n2) = ( ctmp(n+1:n2) +  x2(n+1:n2) )/r3
         
         x(:,1) = -dimag(y3)/dimag(r3)
         x(:,2) = -dimag(r3*y3)/dimag(r3)
         x(:,3) = x2 - dreal(r3*dconjg(r3))*x(:,1) + 2.0_DP*dreal(r3)*x(:,2)
         x(:,4) = dreal(dconjg(r2)*x2 - y2 - (r3*dconjg(r3)*x(:,2)      &
     &                 - 2.0_DP*dreal(r3)*x(:,3)));
         x(:,5) = dreal(x1 - r2*dconjg(r2)*x2 +                         &
     &                  (  2.0_DP*dreal(r2)*r3*dconjg(r3)*x(:,2) -      &
     &        ( r3*dconjg(r3) + 4.0_DP*dreal(r2)*dreal(r3)   )*x(:,3) + & 
     &           ( 2.0_DP*dreal(r2)+2.0_DP*dreal(r3) )*x(:,4) )   )
         x(:,6) = b
         DO j = 1, p
            x(:,6) = x(:,6) -  qcoe(j)*x(:,j)
         END DO
         x(:,6) = x(:,6)/qcoe(p+1)
         
         fn = 0.0_DP
         do i = 1, npt
            fn = fn + ft(is,i)*pt(i)
         end do
         
         t(it) = t(it-1) + dt
         DO i = 1, npDOF
            dsp( it, i ) = x( n+pDOF(i),1 )
            vel( it, i ) = x( pDOF(i),1 )
            acc( it, i ) = x( pDOF(i),2 ) + fn*Fb( pDOF(i) ) 
         END DO
         
      END DO

      write(*,*) '   Wall Clock time for back-substitution (sec) : ', solverTotalWallClock 
      
      WallClock = SECNDS(wallClock)
      write(*,*) '   Wall Clock time for stepping (sec) : ', WallClock
      
      vel = vel/dt
      acc = acc/(dt*dt)
      
      CALL MKLpardiso_delete(pM, n)
      CALL MKLpardiso_delete(pKd1, n)
      CALL MKLpardiso_delete(pKd2, n)
      
      DEALLOCATE( Fb )

      DEALLOCATE( f0, Pf, pt )
      DEALLOCATE( b )
      DEALLOCATE( x )
      
      DEALLOCATE( pcoe, qcoe, cf )
      DEALLOCATE( Kd1, Kd2, x1, y1, y2, ftmp, ctmp )

      
      RETURN
      END SUBROUTINE TimeSolver45
!      
!     ***********************************************************
!
      SUBROUTINE OutputResponse(ich, t, dsp, vel, acc)
 
      IMPLICIT NONE
    
      INTEGER        :: ich
      REAL(KIND=DP)  :: t( : )
      REAL(KIND=DP)  :: dsp( :, : )
      REAL(KIND=DP)  :: vel( :, : )
      REAL(KIND=DP)  :: acc( :, : )

      INTEGER        :: i(2)

      i = shape(dsp)

      WRITE(ich, *) 'ns'
      WRITE(ich, *) i(1)
      WRITE(ich, *) 'npDOF'
      WRITE(ich, *) i(2)
      WRITE(ich, *) 't'
      WRITE(ich, *) t
      WRITE(ich, *) 'dsp'
      WRITE(ich, *) dsp
      WRITE(ich, *) 'vel'
      WRITE(ich, *) vel
      WRITE(ich, *) 'acc'
      WRITE(ich, *) acc
      
      RETURN
      END SUBROUTINE OutputResponse
!      
!     ***********************************************************
!
      SUBROUTINE  MKLpardiso_Ax(p, n, ia, ja, K, C, M, x, y)
      
      IMPLICIT NONE
      
      INTEGER, PARAMETER :: dp = KIND(1.0D0)

!..   Internal solver memory pointer 
      TYPE (MKLPardisoType) :: p
  
      INTEGER, INTENT(IN) :: n
      INTEGER, INTENT(IN) :: ia( : ), ja( : )
      REAL(KIND=DP), INTENT(IN) :: K( : )
      REAL(KIND=DP), INTENT(IN) :: C( : )
      REAL(KIND=DP), INTENT(IN) :: M( : )
      REAL(KIND=DP), INTENT(IN) :: x( : )
      REAL(KIND=DP), INTENT(OUT) :: y( : )

      INTEGER :: n2
      REAL(KIND=DP), ALLOCATABLE :: ftmp( : )

      n2 = n + n

      allocate( ftmp(n2) )

      call mkl_dcsrsymv('U', n, C, ia, ja, x(1:n),    ftmp(1:n))
      call mkl_dcsrsymv('U', n, K, ia, ja, x(n+1:n2), ftmp(n+1:n2))
      ftmp(1:n) = - ftmp(1:n) -  ftmp(n+1:n2) 
      
      call MKLpardiso_solve_real(p, n, ia, ja, M, ftmp(1:n), y(1:n))
      y(n+1:n2) = x(1:n)

      
      END SUBROUTINE  MKLpardiso_Ax
      
!      
!     ***********************************************************
!
      SUBROUTINE TimeIntgCoeff(M, aInfty, p, q, r, C)

      IMPLICIT NONE

      integer, intent(in) :: M
      double precision, intent(in) :: aInfty

      double precision, intent(out) :: p(0:M), q(0:M)
      double complex,   intent(out) :: r(M) 
      double precision, intent(out) :: C(M+1,M)


      integer M1, L, i, j, info
      double precision, allocatable :: p1(:), q1(:)
      double precision, allocatable :: p2(:), q2(:)
      double precision, allocatable :: fML(:), fM(:) 

      double precision f1, f2
      double precision, allocatable :: a(:,:), wr(:), wi(:), fv1(:) 


      M1 =  M + 1
      L = M - 1
      allocate (p1(0:M), q1(0:M), p2(0:L), q2(0:M))
      allocate (fML(-1:M), fM(0:M))

      fM(0) = 1.d0
      fML(0) = 1.d0
      do i = 1, M
         fM(i) = fM(i-1)*dble(i)
         fML(i) = fML(i-1)*(M+i)
      enddo
      fML(-1) = fM(L)
      fML(0:M) = fML(0:M)*fM(M)
      p1 = fML(M:0:-1)/(fM(0:M)*fM(M:0:-1))
      q1 = p1
      q1(1:M:2) = - q1(1:M:2)

      p2 = fML(L:0:-1)/(fM(0:L)*fM(L:0:-1))
      q2 = (fM(M)/fM(L))*(fML(L:-1:-1)/(fM(0:M)*fM(M:0:-1)))
      q2(1:M:2) = - q2(1:M:2)

      p = aInfty*p1
      p(0:L) = p(0:L) + (1.d0-aInfty)*p2
      q = aInfty*q1 + (1.d0-aInfty)*q2

      deallocate(fM, fML, p1, q1, p2, q2)

      if (M .eq. 1) then
            r = -q(0)/q(1)
      else if (M .eq. 2) then
            r = dcmplx(-q(1),dsqrt(-q(1)*q(1)+4.d0*q(0)*q(2)))/(q(2)+q(2))
      else
!ltx  calculate eignvalues and eigenvectors
      allocate(wr(M), wi(M), a(M,M), fv1(4*min(4,M)) )
      a = 0.d0
      a(M,1:M) = -q(0:M-1)/q(M)
      forall (i = 1:M-1) a(i,i+1) = 1.d0

      call DGEEV('N','N',M,a,M,wr,wi,fv1,1,fv1,1,fv1,4*min(4,M),info)
      if(info.ne.0) stop ' calculation of eigenvalues failed'
      j = 1
      do while (j.le.M1/2)    
          i = minloc(wi, 1, MASK = wi.gt.-1.d-6)
          r(j) = dcmplx( wr(i), wi(i))
          wi(i) = -1.d0     
         j = j + 1
      enddo
      deallocate(a, wr, wi, fv1) 
      end if

      allocate( fv1(M1) )
      fv1 = p - q
      C = 0.d0
      C(1,:) = fv1(2:M1)
      f1 = 1.d0
      f2 = 1.d0
      do i = 1, M
         f1 = -f1/2.d0
         f2 = -f2
         fv1 = f1*(p-f2*q)
         fv1(1:M) = fv1(1:M) + i*C(i,:)
         C(i+1,:) = fv1(2:M1)
      end do 
      deallocate(fv1)

      return
      end  SUBROUTINE TimeIntgCoeff

      END MODULE MKLTimeHighOrderWrapper