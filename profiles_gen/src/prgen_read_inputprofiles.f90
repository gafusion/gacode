!--------------------------------------------------------------
! prgen_read_inputprofiles.f90
!
! PURPOSE:
!  Read input.gacode
!--------------------------------------------------------------

subroutine prgen_read_inputprofiles

  use prgen_globals
  use expro

  implicit none
  integer :: i
  integer :: nexp,nion
  character(len=99) :: line
  real :: x(5)
  real :: b_ref,arho,rvbv

  expro_ctrl_quasineutral_flag = 0
  expro_ctrl_numeq_flag = 0 

  !----------------------------------------------

  open(unit=1,file='input.profiles',status='old')
  do while (line(1:2) /= '#r')
     read(1,'(a)') line
     if (line(1:5) == 'N_EXP') then
        read(line(7:),*) expro_n_exp
     endif
     if (line(1:5) == 'N_ION') then
        read(line(7:),*) expro_n_ion
     endif
     if (line(1:6) == 'BT_EXP') then
        read(line(8:),*) b_ref
     endif
     if (line(1:6) == 'IP_EXP') then
        read(line(8:),*) current
     endif
     if (line(1:4) == 'RVBV') then
        read(line(6:),*) rvbv
     endif
     if (line(1:8) == 'ARHO_EXP') then
        read(line(10:),*) arho
     endif
  enddo
  torfluxa = 0.5*b_ref*arho**2

  call expro_init(1)

  nexp = expro_n_exp
  nion = expro_n_ion

  nx = expro_n_exp

  call prgen_allocate

  ! 1
  do i=1,nexp
     read(1,*) x
     rho(i)        = x(1)
     expro_rmin(i) = x(2)
     dpsi(i)       = x(3)
     q(i)          = x(4)
     expro_w0(i)   = x(5)
  enddo

  read(1,'(a)') line
  read(1,'(a)') line

  ! 2
  do i=1,nexp
     read(1,*) x
     expro_rmaj(i)  = x(1)
     zmag(i)  = x(2)
     kappa(i) = x(3)
     delta(i) = x(4)
     zeta(i)  = x(5)
  enddo

  read(1,'(a)') line
  read(1,'(a)') line

  ! 3
  do i=1,nexp
     read(1,*) x
     expro_ne(i)    = x(1)
     expro_te(i)    = x(2)
     expro_ptot(i)  = x(3)
     expro_z_eff(i) = x(4)
  enddo

  read(1,'(a)') line
  read(1,'(a)') line

  ! 4
  do i=1,nexp
     read(1,*) x
     expro_ni(1:nion,i) = x(1:nion)
  enddo

  read(1,'(a)') line
  read(1,'(a)') line

  ! 5 (assume < 6 ions, so skip)
  do i=1,nexp
     read(1,*) x
  enddo

  read(1,'(a)') line
  read(1,'(a)') line

  ! 6
  do i=1,nexp
     read(1,*) x
     expro_ti(1:nion,i) = x(1:nion)
  enddo

  read(1,'(a)') line
  read(1,'(a)') line

  ! 7 (assume < 6 ions, so skip)
  do i=1,nexp
     read(1,*) x
  enddo

  read(1,'(a)') line
  read(1,'(a)') line

  ! 8
  do i=1,nexp
     read(1,*) x
     expro_vtor(1:nion,i) = x(1:nion)
  enddo

  read(1,'(a)') line
  read(1,'(a)') line

  ! 9 (assume < 6 ions, so skip)
  do i=1,nexp
     read(1,*) x
  enddo

  read(1,'(a)') line
  read(1,'(a)') line

  ! 10
  do i=1,nexp
     read(1,*) x
     expro_vpol(1:nion,i) = x(1:nion)
  enddo

  read(1,'(a)') line
  read(1,'(a)') line

  ! 11 (assume < 6 ions, so skip)
  do i=1,nexp
     read(1,*) x
  enddo
  close(1)

  ! Needed for diagnostic printing
  rmin(:) = expro_rmin(:)
  rmaj(:) = expro_rmaj(:)

  ! Missing stuff
  expro_name = 'unknown'
  expro_type = 'unknown'

  bcentr = b_ref
  rcentr = rvbv/b_ref
  expro_rho(:) = rho

  open(unit=1,file='profile_shot',status='old')
  read(1,*) expro_shot
  read(1,*) expro_time
  close(1)

  open(unit=1,file='profile_header',status='old')
  do i=1,nion
     read(1,*) expro_z(i),expro_mass(i),expro_type(i)
     call prgen_ion_name(nint(expro_mass(i)),nint(expro_z(i)),expro_name(i))
  enddo
  close(1)

end subroutine prgen_read_inputprofiles
