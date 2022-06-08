   INCLUDE 'mkl_pardiso.f90'
   MODULE MKLpardiso_solver
      
      USE mkl_pardiso

      IMPLICIT NONE

      INTEGER, PRIVATE, PARAMETER :: dp = KIND(1.0D0)

      TYPE MKLPardisoType
         TYPE(MKL_PARDISO_HANDLE) :: pt( 64 )
!         INTEGER*8 pt(64)
         INTEGER mtype
                        ! mtype  = -2 ! symmetric, indefinite
                        ! mtype  = 6 ! complex and symmetric
                        ! mtype  = 2 ! real and symmetric positive definite
         INTEGER,  PRIVATE:: iparm(64)
!         REAL(KIND=DP),  PRIVATE::  dparm(64) 
      END TYPE MKLPardisoType

      INTEGER nproc

      INTEGER,  PRIVATE:: mtype
      ! mtype  = -2 ! symmetric, indefinite
      ! mtype  = 6 ! complex and symmetric
      ! mtype  = 2 ! real and symmetric positive definite

      INTEGER,  PRIVATE:: maxfct, mnum, phase, nrhs, error, msglvl
!      INTEGER,  PRIVATE:: iparm(64)

      INTEGER,  PRIVATE:: idum(1)
      REAL(KIND=DP),  PRIVATE:: ddum(1)
      COMPLEX(KIND=DP),  PRIVATE:: zdum(1)
      
      CONTAINS

!      
!     ***********************************************************
!
      SUBROUTINE MKLpardiso_SetParameter()

      IMPLICIT NONE
      
      character val*40
      integer leng, status
      
      nrhs = 1 
      maxfct = 1 
      mnum = 1

      error  = 0 ! initialize error flag
      msglvl = 0 ! = 1 print statistical information

      call get_environment_variable ('OMP_NUM_THREADS', val, leng, status)
      if ( status.eq.0 .and. leng.gt.0 ) then
            read(val,*) nproc
            write(*,*) 'OMP_NUM_THREADS = ', nproc
      !else
      !      call get_environment_variable ('NUMBER_OF_PROCESSORS', val, leng, status)
      !      if ( status.eq.0 .and. leng.gt.0 ) then
      !            read(val,*) nproc
      !            write(*,*) 'NUMBER_OF_PROCESSORS = ', nproc
      !            nproc = nproc/2
      !            write(*,*) 'OMP_NUM_THREADS chosen = ', nproc
      !      endif
      endif
      nproc = 0
      write (*,*) 'Using MKL solver, revision 1 11-July'
      
      END SUBROUTINE MKLpardiso_SetParameter
!      
!     ***********************************************************
!
      SUBROUTINE MKLpardiso_factor_real(p, n, ia, ja, a)
    
      IMPLICIT NONE

      !.. Internal solver memory pointer
      TYPE (MKLPardisoType) :: p

      INTEGER, INTENT(IN) :: n
      INTEGER, INTENT(IN) :: ia( : )
      INTEGER, INTENT(IN) :: ja( : )
      REAL(KIND=DP), INTENT(IN) :: a( : )

      INTEGER i

!.. Initialize the internal solver memory pointer. This is only
! necessary for the FIRST call of the PARDISO solver.
      DO i = 1, 64
         p%pt(i)%DUMMY =  0
      END DO

!..
!.. Set up PARDISO control parameter
!..
      p.iparm = 0
      p.iparm(1) = 1 ! no solver default
      p.iparm(2) = 2 ! fill-in reordering from METIS
      p.iparm(3) = nproc ! numbers of processors, value of OMP_NUM_THREADS
!      p.iparm(8) = 10 ! numbers of iterative refinement steps
      p.iparm(10) = 13 ! perturb the pivot elements with 1E-13
      p.iparm(11) = 1 ! use nonsymmetric permutation and scaling MPS
      p.iparm(18) = -1 ! Output: number of nonzeros in the factor LU
      p.iparm(19) = -1 ! Output: Mflops for LU factorization
      p.iparm(24) = 1 ! two-level factorization algorithm

!.. Reordering and Symbolic Factorization, This step also allocates
! all memory that is necessary for the factorization

      phase = 11 ! only reordering and symbolic factorization
      CALL pardiso (p%pt, maxfct, mnum, p%mtype, phase, n, a, ia, ja, &
                    idum, nrhs, p%iparm, msglvl, ddum, ddum, error)
    
      WRITE(*,*) 'Reordering completed ... '
      IF (error /= 0) THEN
         WRITE(*,*) 'The following ERROR was detected: ', error
         GOTO 1000
      END IF
      WRITE(*,*) 'Number of nonzeros in factors = ',p%iparm(18)
      WRITE(*,*) 'Number of factorization MFLOPS = ',p%iparm(19)

!.. Factorization.
      phase = 22 ! only factorization
      CALL pardiso (p%pt, maxfct, mnum, p%mtype, phase, n, a, ia, ja, &
                    idum, nrhs, p%iparm, msglvl, ddum, ddum, error)
      WRITE(*,*) 'Factorization  of real matrix completed ... '
      IF (error /= 0) THEN
         WRITE(*,*) 'The following ERROR was detected: ', error
         GOTO 1000
      ENDIF

1000  CONTINUE

      RETURN

      END  SUBROUTINE MKLpardiso_factor_real

!      
!     ***********************************************************
!
      SUBROUTINE MKLpardiso_solve_real(p, n, ia, ja, a, b, x)
    
      IMPLICIT NONE
!..   Internal solver memory pointer 
      TYPE (MKLPardisoType) :: p
  
      INTEGER, INTENT(IN) :: n
      INTEGER, INTENT(IN) :: ia( : )
      INTEGER, INTENT(IN) :: ja( : )
      REAL(KIND=DP), INTENT(IN) :: a( : )
      REAL(KIND=DP) b( * )
      REAL(KIND=DP), INTENT(OUT) :: x( : )


!.. Back substitution and iterative refinement
      phase = 33 ! only solving
      CALL pardiso (p%pt, maxfct, mnum, p%mtype, phase, n, a, ia, ja, &
                    idum, nrhs, p%iparm, msglvl, b, x, error)

      IF (error /= 0) THEN
         WRITE(*,*) 'The following ERROR was detected: ', error
      ENDIF

!       WRITE(*,*) 'Number of iterative refinements: ', p%iparm(7)
      
      RETURN

      END SUBROUTINE MKLpardiso_solve_real

!      
!     ***********************************************************
!
      SUBROUTINE MKLpardiso_factor_complex(p, n, ia, ja, a)
    
      IMPLICIT NONE

      !.. Internal solver memory pointer 
      TYPE (MKLPardisoType) :: p

      INTEGER, INTENT(IN) :: n
      INTEGER, INTENT(IN) :: ia( : )
      INTEGER, INTENT(IN) :: ja( : )
      COMPLEX(KIND=DP), INTENT(IN) :: a( : )

      INTEGER i

!.. Initialize the internal solver memory pointer. This is only
! necessary for the FIRST call of the PARDISO solver.
      DO i = 1, 64
         p%pt(i)%DUMMY =  0
      END DO
!..
!.. Set up PARDISO control parameter
!..
      p.iparm = 0
      p.iparm(1) = 1 ! no solver default
      p.iparm(2) = 2 ! fill-in reordering from METIS
      p.iparm(3) = nproc ! numbers of processors, value of OMP_NUM_THREADS
!      p.iparm(8) = 10 ! numbers of iterative refinement steps
      p.iparm(10) = 13 ! perturb the pivot elements with 1E-13
      p.iparm(11) = 1 ! use nonsymmetric permutation and scaling MPS
      p.iparm(18) = -1 ! Output: number of nonzeros in the factor LU
      p.iparm(19) = -1 ! Output: Mflops for LU factorization


!..Reordering and Symbolic Factorization, This step also allocates
! all memory that is necessary for the factorization

      phase = 11 ! only reordering and symbolic factorization
      CALL pardiso (p%pt, maxfct, mnum, p%mtype, phase, n, a, ia, ja, &
                    idum, nrhs, p%iparm, msglvl, zdum, zdum, error)
    
      WRITE(*,*) 'Reordering completed ... '
      IF (error /= 0) THEN
         WRITE(*,*) 'The following ERROR was detected: ', error
         GOTO 1000
      END IF
      WRITE(*,*) 'Number of nonzeros in factors = ',p%iparm(18)
      WRITE(*,*) 'Number of factorization MFLOPS = ',p%iparm(19)

!.. Factorization.
      phase = 22 ! only factorization
      CALL pardiso (p%pt, maxfct, mnum, p%mtype, phase, n, a, ia, ja, &
                    idum, nrhs, p%iparm, msglvl, zdum, zdum, error)
      WRITE(*,*) 'Factorization of complex matrix completed ... '
      IF (error /= 0) THEN
         WRITE(*,*) 'The following ERROR was detected: ', error
         GOTO 1000
      ENDIF

1000  CONTINUE

      RETURN

      END  SUBROUTINE MKLpardiso_factor_complex
      

!      
!     ***********************************************************
!
      SUBROUTINE MKLpardiso_solve_complex(p, n, ia, ja, a, b, x)
    
      IMPLICIT NONE

!.. Internal solver memory pointer 
      TYPE (MKLPardisoType) :: p
!.. All other variables

      INTEGER, INTENT(IN) :: n
      INTEGER, INTENT(IN) :: ia( : )
      INTEGER, INTENT(IN) :: ja( : )
      COMPLEX(KIND=DP), INTENT(IN) :: a( : )
      COMPLEX(KIND=DP) :: b( : )
      COMPLEX(KIND=DP), INTENT(OUT) :: x( : )

!.. Back substitution and iterative refinement
      phase = 33 ! only solving
      CALL pardiso (p%pt, maxfct, mnum, p%mtype, phase, n, a, ia, ja, &
                    idum, nrhs, p%iparm, msglvl, b, x, error)

      IF (error /= 0) THEN
         WRITE(*,*) 'The following ERROR was detected: ', error
      ENDIF

!     WRITE(*,*) 'Number of iterative refinements: ', p%iparm(7)

      
 1000 CONTINUE

      RETURN
      END SUBROUTINE MKLpardiso_solve_complex
      
!      
!     ***********************************************************
!
      SUBROUTINE MKLpardiso_delete(p, n)
    
      IMPLICIT NONE
!..   Internal solver memory pointer 
      TYPE (MKLPardisoType) :: p
  
      INTEGER, INTENT(IN) :: n

      !.. Termination and release of memory
      phase = -1 ! release internal memory
      CALL pardiso (p%pt, maxfct, mnum, p%mtype, phase, n, ddum, idum, idum, &
                    idum, nrhs, p%iparm, msglvl, ddum, ddum, error)


      IF (error /= 0) THEN
         WRITE(*,*) 'The following ERROR on release stage was detected: ', error
         STOP 1
      ENDIF

      RETURN

      END SUBROUTINE MKLpardiso_delete

      
      END MODULE MKLpardiso_solver
