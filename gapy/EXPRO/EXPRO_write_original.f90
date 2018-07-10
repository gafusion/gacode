!--------------------------------------------------------
! EXPRO_write_original.f90
!
! PURPOSE:
!  Rewrite original quantities to input.profiles
!
!  NOTE:
!   q is kept positive definite.  All information
!   about sign of q is retained in IPCCW and BTCCW.  
!--------------------------------------------------------

subroutine EXPRO_write_original(io1,datafile1,io2,datafile2,tag)

  use EXPRO_globals
  use EXPRO_interface

  implicit none

  integer, intent(in) :: io1,io2
  character (len=*), intent(in) :: datafile1,datafile2
  character (len=*), intent(in) :: tag

  integer :: ierr
  character (len=80) :: line
  
  open(unit=io1,file=trim(path)//trim(datafile1),status='old')
  open(unit=io2,file=trim(path)//trim(datafile2),status='replace')

  write(io2,'(a)') '# '//tag
  write(io2,'(a)') '# '
  do 
     read(io1,'(a)',iostat=ierr) line
     if (ierr < 0) exit
     if (line(2:4) /= 'rho') then
        write(io2,'(a)') line
     else
        exit
     endif
  enddo
  call EXPRO_write(io2)
  
  close(io1)
  close(io2)

end subroutine EXPRO_write_original
