subroutine locpargen_tglf

  use locpargen_globals
  use expro_locsim_interface

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
  write(1,*)
  write(1,'(a)')  '# Geometry (Advanced)'
  write(1,10) 'SHAPE_SIN3=',shape_sin3_loc
  write(1,10) 'SHAPE_S_SIN3=',shape_s_sin3_loc
  write(1,10) 'SHAPE_SIN4=',shape_sin4_loc
  write(1,10) 'SHAPE_S_SIN4=',shape_s_sin4_loc
  write(1,10) 'SHAPE_SIN5=',shape_sin5_loc
  write(1,10) 'SHAPE_S_SIN5=',shape_s_sin5_loc
  write(1,10) 'SHAPE_SIN6=',shape_sin6_loc
  write(1,10) 'SHAPE_S_SIN6=',shape_s_sin6_loc
  write(1,10) 'SHAPE_COS0=',shape_cos0_loc
  write(1,10) 'SHAPE_S_COS0=',shape_s_cos0_loc
  write(1,10) 'SHAPE_COS1=',shape_cos1_loc
  write(1,10) 'SHAPE_S_COS1=',shape_s_cos1_loc
  write(1,10) 'SHAPE_COS2=',shape_cos2_loc
  write(1,10) 'SHAPE_S_COS2=',shape_s_cos2_loc
  write(1,10) 'SHAPE_COS3=',shape_cos3_loc
  write(1,10) 'SHAPE_S_COS3=',shape_s_cos3_loc 
  write(1,10) 'SHAPE_COS4=',shape_cos4_loc
  write(1,10) 'SHAPE_S_COS4=',shape_s_cos4_loc 
  write(1,10) 'SHAPE_COS5=',shape_cos5_loc
  write(1,10) 'SHAPE_S_COS5=',shape_s_cos5_loc 
  write(1,10) 'SHAPE_COS6=',shape_cos6_loc
  write(1,10) 'SHAPE_S_COS6=',shape_s_cos6_loc 
  write(1,*)
  write(1,10) 'SIGN_BT=',btccw
  write(1,10) 'SIGN_IT=',ipccw
  write(1,*)
  write(1,'(a)') '# Collisions and pressure'
  write(1,10) 'ZEFF=',z_eff_loc
  write(1,10) 'XNUE=',nu_ee*a/cs_loc
  write(1,10) 'DEBYE=',7.43*sqrt(1e3*temp_loc(ise)/(1e13*dens_loc(ise)))/abs(rhos_loc)
  write(1,10) 'BETAE=',betae_unit
  write(1,10) 'P_PRIME_LOC=',(abs(q_loc)/r0)*(-beta_star_loc/(8*pi))
  write(1,*)
  write(1,'(a)')  '# Rotation'
  if(q_loc < 0.0) then
     write(1,10) 'VEXB_SHEAR=',ipccw*gamma_e_loc*a/cs_loc
  else
     write(1,10) 'VEXB_SHEAR=',-ipccw*gamma_e_loc*a/cs_loc
  endif
     
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
     write(1,10) 'VPAR_SHEAR_'//tag(isr)//'=',-ipccw*gamma_p_loc*a/cs_loc
     write(1,10) 'VPAR_'//tag(isr)//'=',-ipccw*mach_loc/cs_loc
  enddo
  close(1)

10 format(a,sp,1pe12.5)
11 format(a,i0)

end subroutine locpargen_tglf
