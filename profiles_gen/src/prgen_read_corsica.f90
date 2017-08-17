!--------------------------------------------------------------
! prgen_read_corsica.f90
!
! PURPOSE:
!  Extract corsica data from native TEXT file.
!--------------------------------------------------------------

subroutine prgen_read_corsica

  use prgen_globals

  implicit none

  integer :: i
  real :: index
  character(len=16) :: buf

  ! profile lengths are always the same (I think)
  ! so arrays can be allocated ahead of time
  nx = corsica_nvals
  call allocate_corsica_vars

  open(unit=1, file=raw_data_file, action='read', status='old')

  ! skip header
  do 
     read(1,'(A)',err=10) buf
     if(buf.eq.'#New time slice ') goto 99
  end do
10 continue
  print *, 'Error: could not find time slice data'
  close(1)
  return

99 continue
  ! read scalar data
100 format (1e12.5)
  read(1, 100) corsica_time
  read(1, 100) corsica_current
  read(1, 100) corsica_loop_voltage
  read(1, 100) corsica_fusion_power
  read(1, 100) corsica_bootstrap_fraction
  read(1, 100) corsica_betat
  read(1, 100) corsica_betan
  read(1, 100) corsica_taue
  read(1, 100) corsica_q95
  read(1, 100) corsica_r
  read(1, 100) corsica_a
  read(1, 100) corsica_kappa
  read(1, 100) corsica_delta
  read(1, 100) corsica_fusion_gain
  read(1, 100) corsica_betap
  read(1, 100) corsica_li3

  ! skip comments
  do i=1, 4
     read(1,*)
  end do

  ! read profiles
101 format (15e13.5)
  do i=1, corsica_nvals
     read(1, 101) index, corsica_rho(i), corsica_r_a(i), corsica_psin(i), &
          corsica_vl(i), corsica_te(i), corsica_ti(i), corsica_ne(i), &
          corsica_ndt(i), corsica_nz(i), corsica_nalpha(i), &
          corsica_zeff(i), corsica_q(i), corsica_j(i), corsica_jbs(i)
  end do

  close(1)


  call allocate_internals

  dpsi(:)   = corsica_psin(:)
  rmin(:)   = 0.0
  rmaj(:)   = 0.0
  rho(:)    = 0.0
  kappa(:)  = 0.0
  delta(:)  = 0.0
  zmag(:)   = 0.0
  zeta(:)   = 0.0
  omega0(:) = 0.0

end subroutine prgen_read_corsica
