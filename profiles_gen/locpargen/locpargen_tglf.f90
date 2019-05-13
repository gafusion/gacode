subroutine locpargen_tglf

  use locpargen_globals
  use EXPRO_interface
  use EXPRO_locsim_interface

  implicit none

  integer :: isr
  
  write(1,'(a)') '# Geometry (Miller)'
  write(1,10) 'RMIN_LOC=',r0
  write(1,10) 'RMAJ_LOC=',rmaj_loc
  write(1,*)
  write(1,10) 'DRMAJDX_LOC=',shift_loc
  write(1,10) 'ZMAJ_LOC=',zmag_loc
  write(1,10) 'DZMAJDX_LOC=',dzmag_loc
  write(1,10) 'Q_LOC=',abs(q_loc)
  write(1,10) 'Q_PRIME_LOC=',(q_loc/r0)**2*s_loc
  write(1,10) 'KAPPA_LOC=',kappa_loc
  write(1,10) 'S_KAPPA_LOC=',s_kappa_loc
  write(1,10) 'DELTA_LOC=',delta_loc
  write(1,10) 'S_DELTA_LOC=',s_delta_loc
  write(1,10) 'ZETA_LOC=',zeta_loc
  write(1,10) 'S_ZETA_LOC=',s_zeta_loc
  write(1,10) 'SIGN_BT=',-real(EXPRO_signb)
  write(1,10) 'SIGN_IT=',-real(EXPRO_signq*EXPRO_signb)
  write(1,*)
  write(1,'(a)') '# Collisions and pressure'
  write(1,10) 'ZEFF=',z_eff_loc
  write(1,10) 'XNUE=',nu_ee*a/cs_loc
  write(1,10) 'DEBYE=',7.43*sqrt(1e3*temp_loc(ise)/(1e13*dens_loc(ise)))/abs(rhos_loc)
  write(1,10) 'BETAE=',betae_unit
  write(1,10) 'P_PRIME_LOC=',(abs(q_loc)/r0)*(-beta_star_loc/(8*pi))
  write(1,*)
  write(1,'(a)')  '# Rotation'
  write(1,10) 'VEXB_SHEAR=',gamma_e_loc*a/cs_loc

  !---------------------------------------------------------
  ! Species data
  
  write(1,*)
  write(1,'(a)') '# Species data'
  write(1,11) 'NS=',ise

  do is=1,ise
     if (is == ise) then
        isr = 1
     else
        isr = is+1
     endif
     write(1,*)
     write(1,11) 'ZS_'//tag(isr)//'=',int(z_loc(is))
     ! Deuteron mass normalization
     write(1,10) 'MASS_'//tag(isr)//'=',mass_loc(is)/2.0
     write(1,10) 'AS_'//tag(isr)//'=',dens_loc(is)/dens_loc(ise)
     write(1,10) 'TAUS_'//tag(isr)//'=',temp_loc(is)/temp_loc(ise)
     write(1,10) 'RLNS_'//tag(isr)//'=',dlnndr_loc(is)
     write(1,10) 'RLTS_'//tag(isr)//'=',dlntdr_loc(is)
     write(1,10) 'VPAR_SHEAR_'//tag(isr)//'=',gamma_p_loc*a/cs_loc
     write(1,10) 'VPAR_'//tag(isr)//'=',mach_loc/cs_loc
  enddo
  close(1)

10 format(a,sp,1pe12.5)
11 format(a,i0)

end subroutine locpargen_tglf
