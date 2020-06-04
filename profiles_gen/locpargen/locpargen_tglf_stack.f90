subroutine locpargen_tglf_stack

  use locpargen_globals
  use expro_locsim_interface

  implicit none
  
  if (appendflag == 0) then
     write(1,'(54a16)') 'RMIN_LOC', 'RMAJ_LOC', 'DRMAJDX_LOC', 'ZMAJ_LOC', &
          'DZMAJDX_LOC', 'Q_LOC', 'Q_PRIME_LOC', 'KAPPA_LOC', 'S_KAPPA_LOC', &
          'DELTA_LOC', 'S_DELTA_LOC', 'ZETA_LOC', 'S_ZETA_LOC', 'SIGN_BT', &
          'SIGN_IT', 'ZEFF', 'XNUE', 'DEBYE', 'BETAE', 'P_PRIME_LOC', &
          'VEXB_SHEAR', 'NS', &
          'ZS_2', 'MASS_2', 'AS_2', 'TAUS_2', 'RLNS_2', 'RLTS_2', 'VPAR_SHEAR_2', 'VPAR_2', &
          'ZS_3', 'MASS_3', 'AS_3', 'TAUS_3', 'RLNS_3', 'RLTS_3', 'VPAR_SHEAR_3', 'VPAR_3', &
          'ZS_4', 'MASS_4', 'AS_4', 'TAUS_4', 'RLNS_4', 'RLTS_4', 'VPAR_SHEAR_4', 'VPAR_4', &
          'ZS_1', 'MASS_1', 'AS_1', 'TAUS_1', 'RLNS_1', 'RLTS_1', 'VPAR_SHEAR_1', 'VPAR_1'
 endif

 write(1,'(54g16.8)') r0,rmaj_loc,shift_loc, zmag_loc, &
      dzmag_loc, abs(q_loc), (q_loc/r0)**2*s_loc,kappa_loc,s_kappa_loc, &
      delta_loc,s_delta_loc,zeta_loc, s_zeta_loc,btccw, &
      ipccw, z_eff_loc,nu_ee*a/cs_loc, &
      7.43*sqrt(1e3*temp_loc(ise)/(1e13*dens_loc(ise)))/abs(rhos_loc), &
      betae_unit,(abs(q_loc)/r0)*(-beta_star_loc/(8*pi)), &
      gamma_e_loc*a/cs_loc, &
      ise, &
      int(z_loc(1)), mass_loc(1)/2.0, dens_loc(1)/dens_loc(ise), temp_loc(1)/temp_loc(ise), &
      dlnndr_loc(1), dlntdr_loc(1), gamma_p_loc*a/cs_loc, mach_loc/cs_loc, &
      int(z_loc(2)), mass_loc(2)/2.0, dens_loc(2)/dens_loc(ise), temp_loc(2)/temp_loc(ise), &
      dlnndr_loc(2), dlntdr_loc(2), gamma_p_loc*a/cs_loc, mach_loc/cs_loc, &
      int(z_loc(3)), mass_loc(3)/2.0, dens_loc(3)/dens_loc(ise),  temp_loc(3)/temp_loc(ise), &
      dlnndr_loc(3), dlntdr_loc(3), gamma_p_loc*a/cs_loc, mach_loc/cs_loc, &
      int(z_loc(4)), mass_loc(4)/2.0, dens_loc(4)/dens_loc(ise), temp_loc(4)/temp_loc(ise), &
      dlnndr_loc(4), dlntdr_loc(4), gamma_p_loc*a/cs_loc, mach_loc/cs_loc
  close(1)

end subroutine locpargen_tglf_stack
