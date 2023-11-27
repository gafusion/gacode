!--------------------------------------------------------------
! prgen_read_genf.f90
!
! PURPOSE:
!  Read profiles from General Fusion Statefile
!--------------------------------------------------------------

subroutine prgen_read_genf

  use prgen_globals

  implicit none

  integer :: ios,i,ncol
  character(len=1) :: a

  integer, parameter :: cmax=13
  character(len=1) :: b(cmax)
  real, dimension(cmax) :: vec

  open(unit=1,file=file_state,status='old')
  nx = -1
  do
     read(1,*,iostat=ios) a
     if (ios < 0) exit
     nx = nx+1
  enddo
  close(1)

  call prgen_allocate('')

  ! Determine number of columns
  open(unit=1,file=file_state,status='old')
  read(1,*) b(:)
  do ncol=1,cmax
     if (iachar(b(ncol)) < 49) exit
  enddo
  ncol = ncol-1

  rewind(1)
  read(1,*) 
  do i=1,nx

     read(1,*) vec(1:ncol)
     
     dpsi(i)     = vec(2)
     p_tot(i)    = vec(4)
     q(i)        = vec(5)
     te_kev(i)   = vec(6)
     ti_kev(i)   = vec(7)
     ne_e19m3(i) = vec(8)
     ni_e19m3(i) = vec(9)
     qohm(i)     = vec(12)

     if (ncol > 12) then
        omega0(i) = vec(13)
     endif
     
  enddo

  dpsi = (-dpsi+dpsi(1))/(2*pi)

  ! Units conversions
  te_kev = te_kev*1e-3
  ti_kev = ti_kev*1e-3

  ne_e19m3 = ne_e19m3*1e-19
  ni_e19m3 = ni_e19m3*1e-19

  zeff = 1.0

end subroutine prgen_read_genf
