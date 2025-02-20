subroutine locpargen_cgyro

  use locpargen_globals
  use expro_locsim_interface

  implicit none

  !---------------------------------------------------------
  ! Resolution info (default)

  write(1,'(a)')  '# Resolution'
  write(1,11) 'N_ENERGY=', 8
  write(1,11) 'N_XI=', 16
  write(1,11) 'N_THETA=', 24
  write(1,11) 'N_RADIAL=', 6
  write(1,11) 'N_TOROIDAL=', 1
  write(1,11) 'NONLINEAR_FLAG=', 0
  write(1,*)
  write(1,11) 'BOX_SIZE=', 1
  write(1,10) 'KY=', 0.3
  write(1,*)

  !---------------------------------------------------------
  ! Time Stepping (default)
  write(1,10) 'DELTA_T=', 0.008
  write(1,10) 'MAX_TIME=', 100.0
  write(1,11) 'PRINT_STEP=', 100
  write(1,*)
  
  !---------------------------------------------------------
  ! Geometry data
  
  write(1,'(a)')  '# Geometry (Miller)'
  write(1,10) 'RMIN=',r0
  write(1,10) 'RMAJ=',rmaj_loc
  write(1,*)
  write(1,11) 'EQUILIBRIUM_MODEL=', 2  
  write(1,10) 'IPCCW=',ipccw
  write(1,10) 'BTCCW=',btccw
  write(1,*)
  write(1,10) 'SHIFT=',shift_loc
  write(1,10) 'ZMAG=',zmag_loc
  write(1,10) 'DZMAG=',dzmag_loc
  write(1,10) 'Q=',abs(q_loc)
  write(1,10) 'S=',s_loc
  write(1,10) 'KAPPA=',kappa_loc
  write(1,10) 'S_KAPPA=',s_kappa_loc
  write(1,10) 'DELTA=',delta_loc
  write(1,10) 'S_DELTA=',s_delta_loc
  write(1,10) 'ZETA=',zeta_loc
  write(1,10) 'S_ZETA=',s_zeta_loc
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
  
  !---------------------------------------------------------
  ! Rotation 
  
  write(1,'(a)')  '# Rotation (Sonic)'
  write(1,11) 'ROTATION_MODEL=', 2   
  write(1,10) 'GAMMA_E=',gamma_e_loc*a/cs_loc
  write(1,10) 'GAMMA_P=',gamma_p_loc*a/cs_loc
  write(1,10) 'MACH=',mach_loc/cs_loc
  write(1,*)

  !---------------------------------------------------------
  ! Collisions and data
  write(1,'(a)')  '# Collisions (Sugama)'
  write(1,11) 'COLLISION_MODEL=', 4    
  write(1,10) 'NU_EE=',nu_ee*a/cs_loc
  write(1,*)
  write(1,11) 'N_FIELD=', 2
  write(1,10) 'BETAE_UNIT=',betae_unit
  write(1,10) 'LAMBDA_STAR=',lambda_star
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
     write(1,10) 'SDLNNDR_'//tag(is)//'=',sdlnndr_loc(is)
     write(1,10) 'SDLNTDR_'//tag(is)//'=',sdlntdr_loc(is)
  enddo
  write(1,10) 'SBETA=',sbeta_loc
  close(1)

10 format(a,sp,1pe12.5)
11 format(a,i0)

end subroutine locpargen_cgyro
