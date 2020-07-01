subroutine locpargen_neo

  use locpargen_globals
  use expro_locsim_interface

  implicit none

  !---------------------------------------------------------
  ! Resolution info (default)

  write(1,'(a)')  '# Resolution'
  write(1,11) 'N_ENERGY=', 6
  write(1,11) 'N_XI=', 17
  write(1,11) 'N_THETA=', 17
  write(1,11) 'N_RADIAL=', 1
  write(1,*)
  
  !---------------------------------------------------------
  ! Geometry data 

  write(1,'(a)')  '# Geometry (Miller)'
  write(1,10) 'RMIN_OVER_A=',r0
  write(1,10) 'RMAJ_OVER_A=',rmaj_loc
  write(1,10) 'RHO_STAR=', rhos_loc/a
  write(1,*)
  write(1,11) 'EQUILIBRIUM_MODEL=', 2  
  write(1,10) 'SHIFT=',shift_loc
  write(1,10) 'ZMAG_OVER_A=',zmag_loc
  write(1,10) 'S_ZMAG=',dzmag_loc
  write(1,10) 'Q=',abs(q_loc)
  write(1,10) 'SHEAR=',s_loc
  write(1,10) 'KAPPA=',kappa_loc
  write(1,10) 'S_KAPPA=',s_kappa_loc
  write(1,10) 'DELTA=',delta_loc
  write(1,10) 'S_DELTA=',s_delta_loc
  write(1,10) 'ZETA=',zeta_loc
  write(1,10) 'S_ZETA=',s_zeta_loc
  write(1,10) 'SHAPE_SIN3=',shape_sin3_loc
  write(1,10) 'SHAPE_S_SIN3=',shape_s_sin3_loc
  write(1,10) 'SHAPE_COS0=',shape_cos0_loc
  write(1,10) 'SHAPE_S_COS0=',shape_s_cos0_loc
  write(1,10) 'SHAPE_COS1=',shape_cos1_loc
  write(1,10) 'SHAPE_S_COS1=',shape_s_cos1_loc
  write(1,10) 'SHAPE_COS2=',shape_cos2_loc
  write(1,10) 'SHAPE_S_COS2=',shape_s_cos2_loc
  write(1,10) 'SHAPE_COS3=',shape_cos3_loc
  write(1,10) 'SHAPE_S_COS3=',shape_s_cos3_loc
  write(1,11) 'IPCCW=',int(ipccw)
  write(1,11) 'BTCCW=',int(btccw)
  write(1,*)

  !---------------------------------------------------------
  ! Rotation data

  write(1,'(a)')  '# Rotation (Sonic)'
  write(1,11) 'ROTATION_MODEL=',2   
  write(1,10) 'OMEGA_ROT=',mach_loc/rmaj_loc/cs_loc
  write(1,10) 'OMEGA_ROT_DERIV=',-gamma_p_loc*a/cs_loc/rmaj_loc
  write(1,*)  

  !---------------------------------------------------------
  ! Collision data
  write(1,'(a)')  '# Collisions (FP)'
  write(1,11) 'COLLISION_MODEL=', 4     
  write(1,10) 'NU_1=',nu_ee*a/cs_loc &
       * (z_loc(1)**4) * (dens_loc(1)/dens_loc(ise)) &
       * (mass_loc(ise)/mass_loc(1))**0.5 * (temp_loc(ise)/temp_loc(1))**1.5
  write(1,*)
  
  !---------------------------------------------------------
  ! Species data

  write(1,'(a)')  '# Species'
  write(1,11) 'N_SPECIES=',ise

  do is=1,ise
     write(1,*)
     write(1,11) 'Z_'//tag(is)//'=',int(z_loc(is))
     ! Deuteron mass normalization
     write(1,10) 'MASS_'//tag(is)//'=',mass_loc(is)/2.0
     write(1,10) 'DENS_'//tag(is)//'=',dens_loc(is)/dens_loc(ise)
     write(1,10) 'TEMP_'//tag(is)//'=',temp_loc(is)/temp_loc(ise)
     write(1,10) 'DLNNDR_'//tag(is)//'=',dlnndr_loc(is)
     write(1,10) 'DLNTDR_'//tag(is)//'=',dlntdr_loc(is)
  enddo
  close(1)

10 format(a,sp,1pe12.5)
11 format(a,i0)
  
end subroutine locpargen_neo
