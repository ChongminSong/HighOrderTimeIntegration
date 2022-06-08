      PROGRAM HighOrderTimeIntg

      USE TextInput
      USE MKLTimeHighOrderWrapper
      !USE pardiso6TimeHighOrderWrapper

      IMPLICIT NONE
      INTEGER, PARAMETER :: DP = KIND(1.0D0)

      INTEGER nm

      INTEGER  n, nnz
      INTEGER  porder
      INTEGER, ALLOCATABLE :: ia( : )
      INTEGER, ALLOCATABLE :: ja( : )
      REAL(KIND=DP) aInfty
      REAL(KIND=DP), ALLOCATABLE :: stff(:), mass(:), damping(:)

      REAL(KIND=DP), ALLOCATABLE :: F( : )
      REAL(KIND=DP), ALLOCATABLE :: u0( : ),  v0( : )

      INTEGER                    :: npt, ns
      REAL(KIND=DP)              :: dt
      REAL(KIND=DP), ALLOCATABLE :: ft( :,: )

      INTEGER              :: npDOF
      INTEGER, ALLOCATABLE :: pDOF( : )
      
      REAL(KIND=DP), ALLOCATABLE :: t( : )
      REAL(KIND=DP), ALLOCATABLE :: dsp( :,: )
      REAL(KIND=DP), ALLOCATABLE :: vel( :,: )
      REAL(KIND=DP), ALLOCATABLE :: acc( :,: )

      REAL(4) WallClock

      integer ich, ichb
      INTEGER  ichOut
      character*80 scratchFolder
      character*80 fNameResponse
      

      integer narg
      integer(2) iarg
      character*80 argBuffer
      integer i

      integer, parameter :: ncLinebuffer = 255
      character*ncLinebuffer linebuffer ! buffer for reading
      character*2   kword
      logical endoffile
      integer ierror

      write(*,*) 'revision 2, 11/09/2021'
     
      narg = NARGS( ) 
      if(narg.eq.2) then      
         iarg = 1
         call GETARG(iarg, argBuffer)
         scratchFolder = adjustl(trim(argBuffer))
         i = len_trim(scratchFolder)
         if (scratchFolder(i:i) .ne. '\') then
            scratchFolder(i+1:i+1) = '\'
         endif
      else
         scratchFolder = 'c:\tmp\'
       endif 
      
!      call ReadInputFiles(scratchFolder, fNameResponse)
      ich = 7
      open(ich,file=trim(scratchFolder)//'cntrl.txt',status='OLD')


      READ(ich,*)
      READ(ich,*) porder, aInfty
      porder = porder - porder/10*10

      ichb = 8
      do
          call ReadLine(ich,linebuffer,endoffile,ierror)
          	if(ierror.eq.1) then
                write(*,*) 'a line in the input file is longer than line buffer'
                stop
	        endif

          kword = adjustl(linebuffer(1:2))
          call uppercase(kword)
          select case (kword)
          case ('KM')
              open(ichb,file=trim(scratchFolder)//adjustl(trim(linebuffer(4:60))),status='OLD')
              READ(ichb,*)
              READ(ichb,*) n, nnz, nm
              ALLOCATE(ia(n + 1))
              ALLOCATE(ja(nnz))
              ALLOCATE(stff(nnz))
    
              READ(ichb,*)
              READ(ichb,*) ia

              READ(ichb,*)
              READ(ichb,*) ja

              READ(ichb,*)
              READ(ichb,*) stff

              if (nm .ge. 2) then
                  ALLOCATE( mass(nnz), damping(nnz) ) 
                  READ(ichb,*)
                  READ(ichb,*) mass
                  if (nm .ge. 3) then
                     READ(ichb,*)
                     READ(ichb,*) damping
                   else
                     damping = 0.0_DP
                   endif
              end if
              close(ichb)
              
              ALLOCATE( F(n) , u0(n) , v0(n) )
              F =  0.0_DP
              u0 = 0.0_DP
              v0 = 0.0_DP
      
          case ('FR')
              open(ichb,file=trim(scratchFolder)//adjustl(trim(linebuffer(4:60))),status='OLD')
              READ(ichb,*)
              READ(ichb,*) F
              close(ichb)
          case ('IC')
              open(ichb,file=trim(scratchFolder)//adjustl(trim(linebuffer(4:60))),status='OLD')
              READ(ichb,*)
              READ(ichb,*) u0
              READ(ichb,*)
              READ(ichb,*) v0
              close(ichb)
          case ('TH')
              open(ichb,file=trim(scratchFolder)//adjustl(trim(linebuffer(4:60))),status='OLD')
              READ(ichb,*)
              READ(ichb,*) dt
              READ(ichb,*)
              READ(ichb,*) ns, npt !ns: number of integration steps
              ALLOCATE( ft(ns, npt) ) ! polynomial expansion of force history 
                                      ! at middle points of time steps 
              READ(ichb,*)
              READ(ichb,*) ft
              close(ichb)
          case ('OC')
              open(ichb,file=trim(scratchFolder)//adjustl(trim(linebuffer(4:60))),status='OLD')
              READ(ichb,*)
              READ(ichb,*) npDOF
              ALLOCATE( pDOF(npDOF) )
              READ(ichb,*) pDOF
              close(ichb)
          case('RS')
              fNameResponse = adjustl(trim(linebuffer(4:60)))
          case('EN')
              exit
            case default
          end select
      end do

      ns = ns + 1
      ALLOCATE( t(ns), dsp(ns,npDOF), vel(ns,npDOF), acc(ns,npDOF) )

      WallClock = SECNDS(0.0)
      if (porder.eq.1) then
            call TimeSolver01(aInfty, n, ia, ja, stff, mass, damping, F,     &
           &                      dt, ft, pDOF, t, dsp, vel, acc, u0, v0)
      else if (porder.eq.2) then
            call TimeSolver12(aInfty, n, ia, ja, stff, mass, damping, F,     &
           &                      dt, ft, pDOF, t, dsp, vel, acc, u0, v0)
      else if (porder.eq.3) then
            call TimeSolver23(aInfty, n, ia, ja, stff, mass, damping, F,     &
           &                      dt, ft, pDOF, t, dsp, vel, acc, u0, v0)
      else if (porder.eq.4) then
            call TimeSolver34(aInfty, n, ia, ja, stff, mass, damping, F,     &
           &                      dt, ft, pDOF, t, dsp, vel, acc, u0, v0)
      else if (porder.eq.5) then
            call TimeSolver45(aInfty, n, ia, ja, stff, mass, damping, F,     &
           &                      dt, ft, pDOF, t, dsp, vel, acc, u0, v0)
      else
            write(*,*) 'not supported, porder = ',porder
      endif
      WallClock = SECNDS(wallClock)
      write(*,*) ' Wall Clock time (sec) : ', WallClock
     
      ichOut = 9
      open(ichOut,file=trim(scratchFolder)//adjustl(trim(fNameResponse)),action='WRITE')
      call OutputResponse(ichOut, t, dsp, vel, acc)
      close(ichOut)

      open(10,file=trim(scratchFolder)//'timing.txt',action='WRITE')
      write(10,*) WallClock, '%Time spent on time integration'
      close(10)

END PROGRAM HighOrderTimeIntg
