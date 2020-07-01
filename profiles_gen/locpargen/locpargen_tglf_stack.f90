subroutine locpargen_tglf_stack

  use locpargen_globals
  use expro_locsim_interface

  implicit none

  integer :: isr
  
  write(1,10) 'RMIN_LOC=',r0, 'RMAJ_LOC=',rmaj_loc, 'DRMAJDX_LOC=',shift_loc, &
       'ZMAJ_LOC=',zmag_loc, 'DZMAJDX_LOC=',dzmag_loc, 'Q_LOC=',abs(q_loc), &
       'Q_PRIME_LOC=',(q_loc/r0)**2*s_loc, 'KAPPA_LOC=',kappa_loc, 'S_KAPPA_LOC=',s_kappa_loc, &
       'DELTA_LOC=',delta_loc, 'S_DELTA_LOC=',s_delta_loc, 'ZETA_LOC=',zeta_loc, &
       'S_ZETA_LOC=',s_zeta_loc, 'SIGN_BT=',btccw, 'SIGN_IT=',ipccw, &
       'ZEFF=',z_eff_loc, 'XNUE=',nu_ee*a/cs_loc, &
       'DEBYE=',7.43*sqrt(1e3*temp_loc(ise)/(1e13*dens_loc(ise)))/abs(rhos_loc), &
       'BETAE=',betae_unit, 'P_PRIME_LOC=',(abs(q_loc)/r0)*(-beta_star_loc/(8*pi)), &
       'VEXB_SHEAR=',gamma_e_loc*a/cs_loc, &
       'NS=',ise, &
       'ZS_2=',int(z_loc(1)), &
       'MASS_2=',mass_loc(1)/2.0,  &
       'AS_2=',dens_loc(1)/dens_loc(ise), &
       'TAUS_2=',temp_loc(1)/temp_loc(ise), &
       'RLNS_2=',dlnndr_loc(1), &
       'RLTS_2=',dlntdr_loc(1), &
       'VPAR_SHEAR_2=',gamma_p_loc*a/cs_loc, &
       'VPAR_2=',mach_loc/cs_loc, &
       'ZS_3=',int(z_loc(2)), &
       'MASS_3=',mass_loc(2)/2.0,  &
       'AS_3=',dens_loc(2)/dens_loc(ise), &
       'TAUS_3=',temp_loc(2)/temp_loc(ise), &
       'RLNS_3=',dlnndr_loc(2), &
       'RLTS_3=',dlntdr_loc(2), &
       'VPAR_SHEAR_3=',gamma_p_loc*a/cs_loc, &
       'VPAR_3=',mach_loc/cs_loc, &
       'ZS_4=',int(z_loc(3)), &
       'MASS_4=',mass_loc(3)/2.0,  &
       'AS_4=',dens_loc(3)/dens_loc(ise), &
       'TAUS_4=',temp_loc(3)/temp_loc(ise), &
       'RLNS_4=',dlnndr_loc(3), &
       'RLTS_4=',dlntdr_loc(3), &
       'VPAR_SHEAR_4=',gamma_p_loc*a/cs_loc, &
       'VPAR_4=',mach_loc/cs_loc, &
       'ZS_1=',int(z_loc(4)), &
       'MASS_1=',mass_loc(4)/2.0,  &
       'AS_1=',dens_loc(4)/dens_loc(ise), &
       'TAUS_1=',temp_loc(4)/temp_loc(ise), &
       'RLNS_1=',dlnndr_loc(4), &
       'RLTS_1=',dlntdr_loc(4), &
       'VPAR_SHEAR_1=',gamma_p_loc*a/cs_loc, &
       'VPAR_1=',mach_loc/cs_loc
  close(1)

10 format(sp,21(a,1pe12.5,2x),a,s,i0,2x,4(a,s,i0,2x,sp,7(a,1pe12.5,2x)))

end subroutine locpargen_tglf_stack
