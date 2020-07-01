!-------------------------------------------------------------------
! timer_lib.f90
!
! PURPOSE:
!  A library to manage a series of up to 64 code timers in Fortran.
!
! NOTES:
!  ! Initialize timer
!  call timer_lib_init('mytimer')
!
!  ! Start 
!  call timer_lib_in('mytimer')
!   ...
!  ! Stop 
!  call timer_lib_out('mytimer')
!
!  ! Get value
!  print *,timer_lib_time('mytimer')
!-------------------------------------------------------------------

module timer_lib

  implicit none 

  public :: timer_lib_init
  public :: timer_lib_in
  public :: timer_lib_out

  real, dimension(64) :: timer_cpu=0.0
  real, dimension(64) :: timer_cpu_in
  character(len=19), dimension(64) :: timer_cpu_tag
  integer :: timer_cpu_maxindx=0

contains 

  subroutine timer_lib_init(tag)

    implicit none
    character(len=*), intent(in) :: tag

    timer_cpu_maxindx = timer_cpu_maxindx+1
    timer_cpu_tag(timer_cpu_maxindx) = trim(tag)

  end subroutine timer_lib_init

  subroutine timer_lib_in(tag)

    use mpi

    implicit none
    character(len=*), intent(in) :: tag
    integer :: indx

    do indx=1,timer_cpu_maxindx
       if (trim(tag) == trim(timer_cpu_tag(indx))) then
          timer_cpu_in(indx) = MPI_Wtime()
       endif
    enddo

  end subroutine timer_lib_in

  subroutine timer_lib_out(tag)

    use mpi

    implicit none
    character(len=*), intent(in) :: tag
    integer :: indx

    do indx=1,timer_cpu_maxindx
       if (trim(tag) == trim(timer_cpu_tag(indx))) then
          timer_cpu(indx) = timer_cpu(indx)+MPI_Wtime()-timer_cpu_in(indx)
       endif
    enddo

  end subroutine timer_lib_out

  real function timer_lib_time(tag)

    use mpi

    implicit none
    character(len=*), intent(in) :: tag
    integer :: indx

    do indx=1,timer_cpu_maxindx
       if (trim(tag) == trim(timer_cpu_tag(indx))) then
          timer_lib_time = timer_cpu(indx)
       endif
    enddo

  end function timer_lib_time

end module timer_lib
