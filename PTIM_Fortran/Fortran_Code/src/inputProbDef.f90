      MODULE TextInput


      IMPLICIT NONE
      
      CONTAINS
!      
!     ***********************************************************
!
       subroutine ReadLine(ch_inp,linebuffer,endoffile,ierror)
!

!      read one line of an input file
      implicit none

      integer ch_inp
      integer, parameter :: ncLinebuffer = 255
      character*ncLinebuffer linebuffer ! buffer for reading

      logical endoffile
      integer i,ierror

      endoffile=.FALSE.
      linebuffer=''
      ierror=0
      i=1
      do while(i.le.ncLinebuffer) !loop over characters
         read(ch_inp,fmt='(a)', advance='NO', eor=100, end=101, err=1001) linebuffer(i:i)
         if(iachar(linebuffer(i:i)).eq.9) then
!	        expand tab to 6 blank spaces
            linebuffer(i:i+5)=" "
            i = i + 6
         else
            i = i + 1
	     end if
      end do !loop over characters
         ierror=1
 100     continue ! finished reading one line

      return

 101  endoffile=.TRUE.
      return

 1001 write(*,*) 'input error'
      stop

      END subroutine ReadLine

!      
!     ***********************************************************
!
      subroutine uppercase(string)

      implicit none

      character (*) string
      integer i

      do i=1, len_trim(string)
         if(ichar(string(i:i)).gt.97.and.ichar(string(i:i)).le.122) &
     &                string(i:i) = char( ichar(string(i:i))-32 )
      end do
      
      end subroutine uppercase
  
      END MODULE TextInput