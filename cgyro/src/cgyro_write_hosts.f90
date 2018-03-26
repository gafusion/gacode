!------------------------------------------------
! cgyro_write_hosts.f90
!
! PURPOSE:
!  Write the hostnames used by the application
!------------------------------------------------

subroutine cgyro_write_hosts

  use iso_c_binding
  use mpi
  use cgyro_globals
  use cgyro_io

  !----------------------------------------------
  implicit none
  !
  character(kind=C_CHAR) :: hostnameraw(128)
  character(LEN=128) :: hostname
  character(LEN=192) :: line
  integer :: i,j
  !
  ! Required for MPI-IO: 
  !
  integer :: filemode
  integer :: finfo
  integer :: fhv
  integer :: fstatus(MPI_STATUS_SIZE)
  integer(kind=MPI_OFFSET_KIND) :: disp
  integer(kind=MPI_OFFSET_KIND) :: offset1
  !----------------------------------------------

  interface
    subroutine gethostname(hname,len) bind(C, name='gethostname')
        use iso_c_binding
        implicit none
        character(kind = C_CHAR) :: hname(*)
        integer(kind = C_INT), VALUE :: len
    end subroutine gethostname
  end interface

  call gethostname(hostnameraw,128)
  ! must convert from null terminated to fortran string
  do i = 1, 128
    if (hostnameraw(i) == char(0) ) exit
    hostname(i:i) = hostnameraw(i)
  end do
  do j = i, 128
    hostname(j:j) = ' '
  end do

  do i = 1, 191
    line(i:i) = ' '
  end do
  write (line, "(A,I8,A,I6,A,I6,A,A)") "RANK G ",i_proc," C1 ",i_proc_1," C2 ",i_proc_2," Host ",hostname
  line(192:192) = NEW_LINE('A')


  !-----------------------------------------------
  ! Dump h and blending coefficients:
  !
  filemode = IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE)
  disp     = 0

  offset1 = 192*i_proc

  call MPI_INFO_CREATE(finfo,i_err)

  call MPI_INFO_SET(finfo,"striping_factor","4",i_err)

  call MPI_FILE_OPEN(CGYRO_COMM_WORLD,&
          trim(path)//runfile_hosts,&
          filemode,&
          finfo,&
          fhv,&
          i_err)

  call MPI_FILE_SET_VIEW(fhv,&
          disp,&
          MPI_CHAR,&
          MPI_CHAR,&
          'native',&
          finfo,&
          i_err)

  call MPI_FILE_WRITE_AT(fhv,&
          offset1,&
          line,&
          192,&
          MPI_CHAR,&
          fstatus,&
          i_err)
  if (i_err /= 0) then
     call cgyro_error('ERROR: (CGYRO) MPI_FILE_WRITE_AT in cgyro_write_hosts failed')
     return
  endif

  call MPI_FILE_SYNC(fhv,i_err)
  call MPI_FILE_CLOSE(fhv,i_err)
  call MPI_INFO_FREE(finfo,i_err)

 !---------------------------------------------------------

end subroutine cgyro_write_hosts

