subroutine expromake_read_input

  use expromake_globals

  implicit none

  integer :: i

  open(unit=1,file='input.expromake.gen',status='old')
  read(1,*) nx
  read(1,*) exm_b_ref
  read(1,*) exm_arho
  read(1,*) exm_a_meters
  read(1,*) exm_rmaj
  read(1,*) exm_q
  read(1,*) exm_kappa
  read(1,*) exm_delta
  read(1,*) exm_te_axis
  read(1,*) exm_ne_axis
  read(1,*) exm_mu_e
  read(1,*) exm_z_eff
  read(1,*) exm_w0
  read(1,*) exm_zeta
  read(1,*) exm_zmag
  do i=1,nions_max
     read(1,*) exm_ni_axis(i)
     read(1,*) exm_ti_axis(i)
     read(1,*) exm_z(i)
     read(1,*) exm_mu(i)
  enddo

  read(1,*) set_exm_b_ref
  read(1,*) set_exm_arho
  read(1,*) set_exm_rho
  read(1,*) set_exm_rmin
  read(1,*) set_exm_rmaj
  read(1,*) set_exm_q
  read(1,*) set_exm_kappa
  read(1,*) set_exm_delta
  read(1,*) set_exm_te
  read(1,*) set_exm_ne
  read(1,*) set_exm_z_eff
  read(1,*) set_exm_w0
  read(1,*) set_exm_zeta
  read(1,*) set_exm_zmag
  do i=1,nions_max
     read(1,*) set_exm_ni(i)
     read(1,*) set_exm_ti(i)
  enddo

  read(1,*) exm_te_model
  read(1,*) exm_alte
  read(1,*) exm_ne_model
  read(1,*) exm_alne
  do i=1,nions_max
     read(1,*) exm_ti_model(i)
     read(1,*) exm_alti(i)
     read(1,*) exm_ni_model(i)
     read(1,*) exm_alni(i)
  enddo

  close(1)

end subroutine expromake_read_input
