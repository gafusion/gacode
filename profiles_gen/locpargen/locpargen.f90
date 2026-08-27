!------------------------------------------------------------
! locpargen.f90
!
! PURPOSE:
!  Generate list of local parameters given input radius or 
!  input rho.
!------------------------------------------------------------

program locpargen

  use locpargen_globals
  use expro
  use expro_locsim_interface
  use geo

  implicit none

  integer :: i
  character(len=80) :: path

  path = './'
  
  open(unit=1,file='input.locpargen',status='old')
  read(1,*) r0
  read(1,*) rho0
  read(1,*) psi0
  read(1,*) qnflag
  read(1,*) appendflag
  read(1,*) ntheta
  read(1,*) nion
  close(1)
  
  expro_ctrl_quasineutral_flag = qnflag

  call expro_read('input.gacode',0)
  
  ! Minor radius
  a = expro_rmin(expro_n_exp)

  ! Number of ions to retain
  if (nion <= 0) then
     nion = expro_n_ion
  endif

  ! Electron index
  ise = nion+1
  
  if (rho0 > 0.0) then

     ! Use local rho

     x(1) = rho0
     call cub_spline(expro_rho,expro_rmin/a,expro_n_exp,x,y,1)
     r0 = y(1)

  else if (psi0 > 0.0) then

     ! Use local psi_N

     x(1) = psi0*abs(expro_polflux(expro_n_exp))
     call cub_spline(abs(expro_polflux),expro_rmin/a,expro_n_exp,x,y,1)
     r0 = y(1)

  endif
  
  call expro_locsim_profiles(&
       qnflag,&
       nion+1,&
       r0,&
       btccw,&
       ipccw,&
       a,&
       path,0)
       
  if (ntheta < 0) then
     goto 100
  endif
     
  print '(a,i2,a,i2,a)','INFO: (locpargen) Keeping',nion,' of',expro_n_ion,' ions'
  print 10,'INFO: (locpargen) rmin/a   =',rmin_loc
  print 10,'INFO: (locpargen) rhos/a   =',rhos_loc/a
  !print 10,'rhoi/a   =',rhos_loc/a*sqrt(temp_loc(ise)/temp_loc(1))
  print 10,'INFO: (locpargen) Te [keV] =',temp_loc(ise)
  print 10,'INFO: (locpargen) Ti [keV] =',temp_loc(1)
  print 10,'INFO: (locpargen) Bunit    =',b_unit_loc
  print 10,'INFO: (locpargen) beta_*   =',beta_star_loc
  print 10,'INFO: (locpargen) ky*rhos (n=1) =',abs(q_loc/rmin_loc*rhos_loc/a)

  betae_unit = betae_loc

  ! Compute collision frequency
  !
  !          4 pi ne e^4 ne ln(Lambda)
  !  nu_ei = ------------------------- 
  !            me^(1/2) (2 Te)^(3/2)

  pi = 4.0*atan(1.0)
  cc = sqrt(2.0)*pi*charge_norm_fac**4/(4.0*pi*8.8542)**2 &
       *1e9/(sqrt(mass_deuterium)*temp_norm_fac**1.5)

  loglam = 24.0-log(sqrt(dens_loc(ise)*1e13)/(temp_loc(ise)*1e3))
  nu_ee  = cc*loglam*dens_loc(ise)/(sqrt(mass_loc(ise)/2.0)*temp_loc(ise)**1.5)

  lambda_star = 7.43 * sqrt((1e3*temp_loc(ise))/(1e13*dens_loc(ise)))/rhos_loc

  tag(:) = (/'1','2','3','4','5','6','7','8','9'/)

  ! Write data for banana width calculation
  open(unit=1,file='out.locpargen',status='replace')
  write(1,*) q_loc
  write(1,*) r0
  write(1,*) rmaj_loc
  write(1,*) shift_loc
  write(1,*) kappa_loc
  write(1,*) q_loc*rhos_loc*sqrt(temp_loc(ise)/temp_loc(1))/a
  close(1)

  call fileopen('input.cgyro.locpargen') ; call locpargen_cgyro
  call fileopen('input.gyro.locpargen') ; call locpargen_gyro
  call fileopen('input.tglf.locpargen')  ; call locpargen_tglf
  call fileopen('input.tglf.locpargen_stack') ; call locpargen_tglf_stack
  call fileopen('input.neo.locpargen')   ; call locpargen_neo
  if (qnflag == 0) then 
     print 10,'INFO: (locpargen) Quasineutrality NOT enforced.'
  else
     print 10,'INFO: (locpargen) Quasineutrality enforced.'
  endif
  print 10,'INFO: (locpargen) Wrote input.*.locpargen'

  !---------------------------------------------------------------------------
  ! GEO output
  if (ntheta > 0) then
     allocate(theta(ntheta))
     do i=1,ntheta
        theta(i) = -pi+(i-1)*2*pi/ntheta
     enddo
     geo_model_in = 0
     geo_rmin_in = rmin_loc
     geo_rmaj_in = rmaj_loc
     geo_drmaj_in = shift_loc
     geo_zmag_in = zmag_loc
     geo_dzmag_in = dzmag_loc
     geo_q_in = q_loc
     geo_s_in = s_loc
     geo_kappa_in = kappa_loc
     geo_s_kappa_in = s_kappa_loc
     geo_delta_in = delta_loc
     geo_s_delta_in = s_delta_loc
     geo_zeta_in = zeta_loc
     geo_s_zeta_in = s_zeta_loc

     geo_shape_cos0_in = shape_cos0_loc
     geo_shape_s_cos0_in = shape_s_cos0_loc
     geo_shape_cos1_in = shape_cos1_loc
     geo_shape_s_cos1_in =  shape_s_cos1_loc
     geo_shape_cos2_in = shape_cos2_loc
     geo_shape_s_cos2_in = shape_s_cos2_loc
     geo_shape_cos3_in = shape_cos3_loc
     geo_shape_s_cos3_in = shape_s_cos3_loc
     geo_shape_sin3_in = shape_sin3_loc
     geo_shape_s_sin3_in = shape_s_sin3_loc
     geo_shape_cos4_in = shape_cos4_loc
     geo_shape_s_cos4_in = shape_s_cos4_loc
     geo_shape_sin4_in = shape_sin4_loc
     geo_shape_s_sin4_in = shape_s_sin4_loc
     geo_shape_cos5_in = shape_cos5_loc
     geo_shape_s_cos5_in = shape_s_cos5_loc
     geo_shape_sin5_in = shape_sin5_loc
     geo_shape_s_sin5_in = shape_s_sin5_loc
     geo_shape_cos6_in = shape_cos6_loc
     geo_shape_s_cos6_in = shape_s_cos6_loc
     geo_shape_sin6_in = shape_sin6_loc
     geo_shape_s_sin6_in = shape_s_sin6_loc

     geo_beta_star_in = beta_star_loc
     geo_beta_star_1_in = 0.0 
     geo_beta_star_2_in = 0.0

     call geo_interp(ntheta,theta,.true.)

     open(unit=1,file='bin.locpargen.theta',status='replace',access='stream')
     write(1) real(theta,kind=4)
     write(1) real(geo_grad_r,kind=4)
     write(1) real(geo_bigr,kind=4)
     write(1) real(geo_bigz,kind=4)
     write(1) real(geo_gsin,kind=4)
     write(1) real(geo_gcos1,kind=4)
     write(1) real(geo_gcos2,kind=4)
     write(1) real(geo_bt,kind=4)
     write(1) real(geo_bp,kind=4)
     write(1) real(geo_g_theta,kind=4)
     write(1) real(geo_gq,kind=4)
     write(1) real(geo_captheta,kind=4)
     close(1)
     deallocate(theta)
  endif

10 format(a,1x,f7.5)

100 continue
  
end program locpargen

subroutine fileopen(fname)

  use locpargen_globals

  implicit none

  character(len=*) :: fname

  if (appendflag == 0) then
     open(unit=1,file=trim(fname),status='replace')
  else
     open(unit=1,file=trim(fname),position='append')
     write(1,*) '#---------------------------------------------------------'
  endif

end subroutine fileopen
