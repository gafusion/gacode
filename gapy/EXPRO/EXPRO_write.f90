!--------------------------------------------------------
! EXPRO_write.f90
!
! PURPOSE:
!  Rewrite tag/variables in input.profiles format.  Unit
!  must already be open.
!--------------------------------------------------------

subroutine EXPRO_write(io)

  use EXPRO_globals
  use EXPRO_interface

  implicit none

  integer, intent(in) :: io
  integer :: i,nx

  nx = EXPRO_n_exp
  
  write(io,20) '# '
  write(io,20) '#',EXPRO_rho_tag,EXPRO_rmin_tag,EXPRO_polflux_tag,EXPRO_q_tag,EXPRO_w0_tag
  do i=1,nx
     write(io,10) EXPRO_rho(i),EXPRO_rmin(i),EXPRO_polflux(i),EXPRO_q(i),EXPRO_w0(i)
  enddo
  write(io,20) '# '
  write(io,20) '#',EXPRO_rmaj_tag,EXPRO_zmag_tag,EXPRO_kappa_tag,EXPRO_delta_tag,EXPRO_zeta_tag
  do i=1,nx
     write(io,10) EXPRO_rmaj(i),EXPRO_zmag(i),EXPRO_kappa(i),EXPRO_delta(i),EXPRO_zeta(i)
  enddo
  write(io,20) '# '
  write(io,20) '#',EXPRO_ne_tag,EXPRO_te_tag,EXPRO_ptot_tag,EXPRO_z_eff_tag,EXPRO_null_tag
  do i=1,nx
     write(io,10) EXPRO_ne(i),EXPRO_te(i),EXPRO_ptot(i),EXPRO_z_eff(i),0.0
  enddo
  write(io,20) '# '
  write(io,20) '#',EXPRO_ni_tag(1:5)
  do i=1,nx
     write(io,10) EXPRO_ni(1:5,i)
  enddo
  write(io,20) '# '
  write(io,20) '#',EXPRO_ni_tag(6:10)
  do i=1,nx
     write(io,10) EXPRO_ni(6:10,i)
  enddo
  write(io,20) '# '
  write(io,20) '#',EXPRO_ti_tag(1:5)
  do i=1,nx
     write(io,10) EXPRO_ti(1:5,i)
  enddo
  write(io,20) '# '
  write(io,20) '#',EXPRO_ti_tag(6:10)
  do i=1,nx
     write(io,10) EXPRO_ti(6:10,i)
  enddo
  write(io,20) '# '
  write(io,20) '#',EXPRO_vtor_tag(1:5)
  do i=1,nx
     write(io,10) EXPRO_vtor(1:5,i)
  enddo
  write(io,20) '# '
  write(io,20) '#',EXPRO_vtor_tag(6:10)
  do i=1,nx
     write(io,10) EXPRO_vtor(6:10,i)
  enddo
  write(io,20) '# '
  write(io,20) '#',EXPRO_vpol_tag(1:5)
  do i=1,nx
     write(io,10) EXPRO_vpol(1:5,i)
  enddo
  write(io,20) '# '
  write(io,20) '#',EXPRO_vpol_tag(6:10)
  do i=1,nx
     write(io,10) EXPRO_vpol(6:10,i)
  enddo
  write(io,20) '# '
  write(io,20) '#',EXPRO_flow_beam_tag,EXPRO_flow_wall_tag,EXPRO_flow_mom_tag,EXPRO_null_tag,EXPRO_null_tag
  do i=1,nx
     write(io,10) EXPRO_flow_beam(i),EXPRO_flow_wall(i),EXPRO_flow_mom(i),0.0,0.0
  enddo
  write(io,20) '# '
  write(io,20) '#',EXPRO_pow_e_tag,EXPRO_pow_i_tag,EXPRO_pow_ei_tag,EXPRO_pow_e_aux_tag,EXPRO_pow_i_aux_tag
  do i=1,nx
     write(io,10) EXPRO_pow_e(i),EXPRO_pow_i(i),EXPRO_pow_ei(i),EXPRO_pow_e_aux(i),EXPRO_pow_i_aux(i)
  enddo
  write(io,20) '# '
  write(io,20) '#',EXPRO_pow_e_fus_tag,EXPRO_pow_i_fus_tag,EXPRO_pow_e_sync_tag,EXPRO_pow_e_brem_tag,EXPRO_pow_e_line_tag
  do i=1,nx
     write(io,10) EXPRO_pow_e_fus(i),EXPRO_pow_i_fus(i),EXPRO_pow_e_sync(i),EXPRO_pow_e_brem(i),EXPRO_pow_e_line(i)
  enddo
  write(io,20) '# '
  write(io,20) '#',EXPRO_sbeame_tag,EXPRO_sbcx_tag,EXPRO_sscxl_tag,EXPRO_null_tag,EXPRO_null_tag
  do i=1,nx
     write(io,10) EXPRO_sbeame(i),EXPRO_sbcx(i),EXPRO_sscxl(i),0.0,0.0
  enddo

10 format(5(1pe14.7,2x))
20 format(6(a))

end subroutine EXPRO_write
