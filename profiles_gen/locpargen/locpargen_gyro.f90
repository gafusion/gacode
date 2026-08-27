subroutine locpargen_gyro

  use locpargen_globals
  use expro_locsim_interface

  implicit none

  !---------------------------------------------------------
  ! Resolution info (default)

  write(1,'(a)')  '# Resolution'
  write(1,11) 'RADIAL_GRID=', 6
  write(1,11) 'PASS_GRID=', 4
  write(1,11) 'TRAP_GRID=', 4
  write(1,11) 'BLEND_GRID=', 6
  write(1,11) 'ENERGY_GRID=', 6
  write(1,11) 'N_TOROIDAL=', 1
  write(1,11) 'NONLINEAR_FLAG=', 0
  write(1,*)
  write(1,10) 'BOX_MULTIPLIER=',1.0
  write(1,*)

  !---------------------------------------------------------
  ! Time Stepping (default)
  write(1,10) 'TIME_STEP=',0.2
  write(1,10) 'TIME_MAX=', 15.0
  write(1,11) 'TIME_SKIP=', 5
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
  write(1,10) 'SAFETY_FACTOR=',abs(q_loc)
  write(1,10) 'SHEAR=',s_loc
  write(1,10) 'KAPPA=',kappa_loc
  write(1,10) 'S_KAPPA=',s_kappa_loc
  write(1,10) 'DELTA=',delta_loc
  write(1,10) 'S_DELTA=',s_delta_loc
  write(1,10) 'ZETA=',zeta_loc
  write(1,10) 'S_ZETA=',s_zeta_loc
  write(1,*)
  
  !---------------------------------------------------------
  ! Rotation 
  
  write(1,'(a)')  '# Rotation'
  write(1,10) 'GAMMA_E=',gamma_e_loc*a/cs_loc
  write(1,10) 'GAMMA_P=',gamma_p_loc*a/cs_loc
  write(1,10) 'MACH=',mach_loc/cs_loc
  write(1,*)

  !---------------------------------------------------------
  ! Collisions and beta
  write(1,'(a)')  '# Collisions'
  write(1,10) 'NU_EE=',nu_ee*a/cs_loc
  write(1,10) 'BETAE_UNIT=',betae_unit
  write(1,*)
  
  !---------------------------------------------------------
  ! Species data
  
  write(1,'(a)')  '# Species'

  do is=1,ise
     if (is == 1) then
        mytag = ''
     else
        mytag = '_'//tag(is)
     endif
     write(1,*)
     write(1,11) 'Z'//mytag//'=',int(z_loc(is))
     ! Deuteron mass normalization
     write(1,10) 'MU'//mytag//'=',mass_loc(is)/2.0
     write(1,10) 'DENS_'//tag(is)//'=',dens_loc(is)/dens_loc(ise)
     write(1,10) 'TEMP_'//tag(is)//'=',temp_loc(is)/temp_loc(ise)
     write(1,10) 'DLNNDR_'//tag(is)//'=',dlnndr_loc(is)
     write(1,10) 'DLNTDR_'//tag(is)//'=',dlntdr_loc(is)
  enddo
  close(1)

10 format(a,sp,1pe12.5)
11 format(a,i0)

end subroutine locpargen_gyro
