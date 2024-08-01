!-----------------------------------------------------------
! tgyro_init_profiles.f90
!
! PURPOSE:
!  Manage generation of profiles on transport grid.
!----------------------------------------------------------

subroutine tgyro_init_profiles

  use mpi
  use tgyro_globals
  use tgyro_ped
  use expro

  implicit none

  integer :: i_ion
  integer :: i
  real :: tmp_ped
  real :: p_ave
  real :: x0(1),y0(1)
  real :: dx,dx_min,dx_max

  !------------------------------------------------------
  ! PHYSICAL CONSTANTS
  !
  pi = 4.0*atan(1.0)
  !
  e_alpha = 3.5e6*1.6022e-12 ! eV*(erg/eV)
  e       = 4.8032e-10 ! statcoul
  k       = 1.6022e-12 ! erg/eV
  me      = 9.1094e-28 ! g
  md      = expro_mass_deuterium ! g
  malpha  = 2*md       ! g
  c       = 2.9979e10  ! cm/s
  !
  mu_0    = 4*pi*1e-7  ! N/A^2
  !------------------------------------------------------

  !------------------------------------------------------
  ! Construct vector of thermal ions
  !
  ! For example, if ions 1 and 3 are thermal, we have
  !
  !  therm_vec = (1,3)  
  !
  ! Size of therm_vec is number of thermal ions.
  !
  i=0
  do i_ion=1,loc_n_ion
     if (therm_flag(i_ion) == 1) then
        i = i+1
        therm_vec(i) = i_ion
     endif
  enddo
  !------------------------------------------------------

  !-----------------------------------------------------------
  ! Initialize integrated power vectors (may not be necessary)
  !
  p_i_fus(:) = 0.0
  p_e_fus(:) = 0.0
  p_exch(:)  = 0.0
  p_brem(:)  = 0.0
  p_sync(:)  = 0.0
  p_brem(:)  = 0.0
  p_expwd(:) = 0.0
  f_he_fus(:) = 0.0
  !------------------------------------------------------------

  !----------------------------------------------
  ! Generate radial vector:
  !
  if (tgyro_rmin > 0.0 .and. n_r > 2) then

     r(:) = 0.0

     do i=2,n_r
        ! Normalized r (dimensionless)
        r(i) = tgyro_rmin+(i-2)/(n_r-2.0)*(tgyro_rmax-tgyro_rmin)
     enddo

  else

     do i=1,n_r
        ! Normalized r (dimensionless)
        r(i) = (i-1)/(n_r-1.0)*tgyro_rmax
     enddo
     
  endif

  ! Overwrite radii with special values
  do i=2,n_r
     if (inputrads(i-1) > 0.0) r(i) = inputrads(i-1)
  enddo

  if (tgyro_use_rho == 1) then
     ! Using equally-spaced rho grid, not default r grid.  This is 
     ! useful for benchmarking with other codes.
     rho = r
  endif
  !----------------------------------------------

  !----------------------------------------------
  ! Check for highly-nonuniform grid
  !
  dx_min = 1.0
  dx_max = 0.0
  do i=2,n_r
     dx = r(i)-r(i-1) 
     if (dx < dx_min) dx_min = dx
     if (dx > dx_max) dx_max = dx
  enddo
  if (dx_max/dx_min > 10.0) then
     use_trap = 1
  else
     use_trap = 0
  endif
  !----------------------------------------------
  
  expro_ctrl_n_ion = loc_n_ion
  expro_ctrl_quasineutral_flag = 0
  expro_ctrl_numeq_flag = 0

  call expro_read('input.gacode',MPI_COMM_WORLD) 

  ! Check for acceptable number of ions
  if (expro_n_ion < loc_n_ion) then
     call tgyro_catch_error('ERROR: (tgyro_init_profiles) LOC_N_ION > expro_n_ion')
  endif

  ! Mass and charge taken from input.gacode
  mi_vec(:) = expro_mass(1:loc_n_ion)
  zi_vec(:) = expro_z(1:loc_n_ion)

  ! Get ion names (truncated to 3 characters)
  ion_name = expro_name(1:loc_n_ion)
  
  shot = 0
  
  n_exp = expro_n_exp

  ! r_min in m:
  r_min = expro_rmin(n_exp)

  ! Aspect ratio
  aspect_rat = expro_rmaj(n_exp)/expro_rmin(n_exp)

  ! Is this a DT plasma
  dt_flag = 0
  if (loc_n_ion > 1) then
     if (ion_name(1) == 'D' .and. ion_name(2) == 'T') dt_flag = 1
     if (ion_name(1) == 'T' .and. ion_name(2) == 'D') dt_flag = 1
  else if (loc_n_ion .eq. 1) then
     if (ion_name(1) == 'D' .and. expro_name(2) == 'T') dt_flag = 2
     if (ion_name(1) == 'T' .and. expro_name(2) == 'D') dt_flag = 2
  endif
  
  !------------------------------------------------------
  ! Convert dimensionless mass to grams.
  mi(:) = mi_vec(:)*(md*0.5)
  !------------------------------------------------------
  
  !------------------------------------------------------------------------------------------
  ! Direct input of simple profiles:
  !
  if (tgyro_use_rho == 1) then
     ! Equally-spaced in rho
     call cub_spline(expro_rho(:),expro_rmin(:)/r_min,n_exp,rho,r,n_r)
  else  
     ! Equally-spaced in r (default)
     call cub_spline(expro_rmin(:)/r_min,expro_rho(:),n_exp,r,rho,n_r)
  endif
  call cub_spline(expro_rmin(:)/r_min,expro_q(:),n_exp,r,q,n_r)
  call cub_spline(expro_rmin(:)/r_min,expro_s(:),n_exp,r,s,n_r)
  call cub_spline(expro_rmin(:)/r_min,expro_kappa(:),n_exp,r,kappa,n_r)
  call cub_spline(expro_rmin(:)/r_min,expro_delta(:),n_exp,r,delta,n_r)
  call cub_spline(expro_rmin(:)/r_min,expro_skappa(:),n_exp,r,s_kappa,n_r)
  call cub_spline(expro_rmin(:)/r_min,expro_sdelta(:),n_exp,r,s_delta,n_r)
  ! Convert r_maj to cm (from m):
  call cub_spline(expro_rmin(:)/r_min,100*expro_rmaj(:),n_exp,r,r_maj,n_r)
  call cub_spline(expro_rmin(:)/r_min,expro_drmaj(:),n_exp,r,shift,n_r)
  ! Convert zmag to cm (from m):
  call cub_spline(expro_rmin(:)/r_min,100*expro_zmag(:),n_exp,r,zmag,n_r)
  call cub_spline(expro_rmin(:)/r_min,expro_dzmag(:),n_exp,r,dzmag,n_r)
  call cub_spline(expro_rmin(:)/r_min,expro_zeta(:),n_exp,r,zeta,n_r)
  call cub_spline(expro_rmin(:)/r_min,expro_szeta(:),n_exp,r,s_zeta,n_r)

  ! Miller Extended Harmonic (MXH) geometry
  call cub_spline(expro_rmin(:)/r_min,expro_shape_sin3(:),n_exp,r,shape_sin3,n_r)
  call cub_spline(expro_rmin(:)/r_min,expro_shape_ssin3(:),n_exp,r,shape_ssin3,n_r)
  call cub_spline(expro_rmin(:)/r_min,expro_shape_cos0(:),n_exp,r,shape_cos0,n_r)
  call cub_spline(expro_rmin(:)/r_min,expro_shape_scos0(:),n_exp,r,shape_scos0,n_r)
  call cub_spline(expro_rmin(:)/r_min,expro_shape_cos1(:),n_exp,r,shape_cos1,n_r)
  call cub_spline(expro_rmin(:)/r_min,expro_shape_scos1(:),n_exp,r,shape_scos1,n_r)
  call cub_spline(expro_rmin(:)/r_min,expro_shape_cos2(:),n_exp,r,shape_cos2,n_r)
  call cub_spline(expro_rmin(:)/r_min,expro_shape_scos2(:),n_exp,r,shape_scos2,n_r)
  call cub_spline(expro_rmin(:)/r_min,expro_shape_cos3(:),n_exp,r,shape_cos3,n_r)
  call cub_spline(expro_rmin(:)/r_min,expro_shape_scos3(:),n_exp,r,shape_scos3,n_r)

  ! Convert psi in Weber to Maxwell (1 Weber = 10^8 Maxwell)  
  call cub_spline(expro_rmin(:)/r_min,expro_polflux(:)*1e8,n_exp,r,polflux,n_r)

  ! b_ref in Gauss (used for wce in Synchroton rad)
  call cub_spline(expro_rmin(:)/r_min,1e4*expro_bt0(:),n_exp,r,b_ref,n_r)
  
  ! Convert ptot to Ba from Pascals (1 Pa = 10 Ba)
  call cub_spline(expro_rmin(:)/r_min,expro_ptot(:)*10.0,n_exp,r,ptot,n_r)

  ! Convert fpol to Gauss-cm from T-m (1 T-m = 10^6 G-cm)
  call cub_spline(expro_rmin(:)/r_min,expro_fpol(:)*1e6,n_exp,r,fpol,n_r)

  ! Convert V and dV/dr from m^3 to cm^3
  call cub_spline(expro_rmin(:)/r_min,expro_vol(:)*1e6,n_exp,r,vol,n_r)
  call cub_spline(expro_rmin(:)/r_min,expro_volp(:)*1e4,n_exp,r,volp,n_r)

  call cub_spline(expro_rmin(:)/r_min,expro_ave_grad_r(:),n_exp,r,ave_grad_r,n_r)

  ! Convert B to Gauss (from T):
  call cub_spline(expro_rmin(:)/r_min,1e4*expro_bunit(:),n_exp,r,b_unit,n_r)

  ! Convert T to eV (from keV) and length to cm (from m):
  call cub_spline(expro_rmin(:)/r_min,1e3*expro_te(:),n_exp,r,te,n_r)
  call cub_spline(expro_rmin(:)/r_min,1e13*expro_ne(:),n_exp,r,ne,n_r)
  call cub_spline(expro_rmin(:)/r_min,expro_dlnptotdr(:)/100.0,n_exp,r,dlnptotdr,n_r)
  call cub_spline(expro_rmin(:)/r_min,expro_dlnnedr(:)/100.0,n_exp,r,dlnnedr,n_r)
  call cub_spline(expro_rmin(:)/r_min,expro_dlntedr(:)/100.0,n_exp,r,dlntedr,n_r)
  do i_ion=1,loc_n_ion
     call cub_spline(expro_rmin(:)/r_min,1e3*expro_ti(i_ion,:),n_exp,r,ti(i_ion,:),n_r)
     call cub_spline(expro_rmin(:)/r_min,1e13*expro_ni(i_ion,:),n_exp,r,ni(i_ion,:),n_r)
     call cub_spline(expro_rmin(:)/r_min,expro_dlntidr(i_ion,:)/100.0,n_exp,r,dlntidr(i_ion,:),n_r)
     call cub_spline(expro_rmin(:)/r_min,expro_dlnnidr(i_ion,:)/100.0,n_exp,r,dlnnidr(i_ion,:),n_r)
     ! Define default ratios (these will change if tgyro_ped_ratio < 0.0)
     n_ratio(i_ion) = ni(i_ion,n_r)/ne(n_r)
     t_ratio(i_ion) = ti(i_ion,n_r)/te(n_r)
     if (ni(i_ion,1) < 1e-10) then
        call tgyro_catch_error('ERROR (tgyro_init_profiles) Zero density in ion '//ion_name(i_ion))
     endif
  enddo

  if (tgyro_consistent_flag == 1) then
     ! Exact recovery of TGYRO-grid gradients
     call math_zfind(n_r,te,r*(100*r_min),dlntedr)
     call math_zfind(n_r,ne,r*(100*r_min),dlnnedr)
     do i_ion=1,loc_n_ion
        call math_zfind(n_r,ti(i_ion,:),r*(100*r_min),dlntidr(i_ion,:))
        call math_zfind(n_r,ni(i_ion,:),r*(100*r_min),dlnnidr(i_ion,:))
     enddo
  endif

  if (tgyro_ptot_flag == 1) then

     ! Total pressure correction from included species

     pr(:) = ne(:)*k*te(:)
     do i_ion=1,loc_n_ion
        pr(:) = pr(:)+ni(i_ion,:)*k*ti(i_ion,:)
     enddo
     pext(:) = ptot(:)-pr(:)
     pr(:)   = ptot(:) 

     dlnpdr(:) = ne(:)*k*te(:)*(dlnnedr(:)+dlntedr(:))/pr(:)
     do i_ion=1,loc_n_ion
        dlnpdr(:) = dlnpdr(:)+&
             ni(i_ion,:)*k*ti(i_ion,:)*(dlnnidr(i_ion,:)+dlntidr(i_ion,:))/pr(:)
     enddo
     dpext(:) = pr(:)*(dlnptotdr(:)-dlnpdr(:))
     dlnpdr(:) = dlnptotdr(:)

  else

     pext(:)  = 0.0
     dpext(:) = 0.0

  endif
  !------------------------------------------------------------------------------------------

  !------------------------------------------------------------------------------------------
  ! Quasineutrality:
  !
  if (loc_lock_profile_flag == 0) then

     ! Update some species (evo_e=1) based on quasineutrality
     call tgyro_quasigrad

     ! Reintegrate density profiles
     do i_ion=1,loc_n_ion
        ! ni in 1/cm^3
        call math_scaleintv(dlnnidr(i_ion,:),r,n_r,ni(i_ion,:),'log')
     enddo

  endif
  !------------------------------------------------------------------------------------------

  !------------------------------------------------------------------------------------------
  ! Helium ash detection (diagnostic only -- set alpha source with evo_e=2)
  !
  i_ash   = 0
  i_alpha = 0
  do i=1,loc_n_ion
     if (trim(ion_name(i)) == 'He' .and. therm_flag(i) == 1) then
        i_ash = i
     endif
     if (trim(ion_name(i)) == 'He' .and. therm_flag(i) == 0) then
        i_alpha = i
     endif
  enddo
  !------------------------------------------------------------------------------------------

  !------------------------------------------------------------------------------------------
  ! Rotation and rotation shear:
  !
  ! w0 (1/s)
  call cub_spline(expro_rmin(:)/r_min,expro_w0(:),n_exp,r,w0,n_r)
  ! w0p = d(w0)/dr (1/s/cm)
  call cub_spline(expro_rmin(:)/r_min,expro_w0p(:)/100.0,n_exp,r,w0p,n_r)
  ! v_pol
  do i_ion=1,loc_n_ion
    call cub_spline(expro_rmin(:)/r_min,1e2*expro_vpol(i_ion,:),n_exp,r,v_pol(i_ion,:),n_r)
  enddo  
  ! vtor_p = d(vtor)/dr 
  call bound_deriv(vtor_p, w0*r_maj, r, n_r)
  ! vpol_p = d(vpol)/dr 
  call bound_deriv(vpol_p, v_pol, r, n_r)

  !------------------------------------------------------------------------------------------

  !------------------------------------------------------------------------------------------
  ! Field orientation
  !
  ! signb = -btccw          OR     btccw = -signb  
  ! signq = ipccw*btccw            ipccw = -signb*signq
  !
  signb = expro_signb
  signq = expro_signq
  !------------------------------------------------------------------------------------------

  !------------------------------------------------------------------------------------------
  ! Apply rescaling factors 
  !
  ne(:) = tgyro_input_den_scale*ne(:)
  te(:) = tgyro_input_te_scale*te(:)
  do i_ion=1,loc_n_ion
     ni(i_ion,:) = tgyro_input_den_scale*ni(i_ion,:)
     ti(i_ion,:) = tgyro_input_ti_scale*ti(i_ion,:)
  enddo
  w0(:) = tgyro_input_w0_scale*w0(:)
  w0p(:) = tgyro_input_w0_scale*w0p(:)

  dlntedr(:)   = dlntedr(:)  *tgyro_input_dlntdr_scale
  dlntidr(:,:) = dlntidr(:,:)*tgyro_input_dlntdr_scale
  !------------------------------------------------------------------------------------------

  !------------------------------------------------------------------------------------------
  ! Z_eff:
  !
  if (loc_zeff_flag == 1) then
     ! Set based on data in input.gacode
     call cub_spline(expro_rmin(:)/r_min,expro_z_eff(:),n_exp,r,z_eff,n_r)
  else
     ! Set to unity
     z_eff = 1.0
  endif
  !------------------------------------------------------------------------------------------

  !------------------------------------------------------------------------------------------
  ! SOURCES:
  !
  ! (1) Power -- convert powers to erg/s from MW:
  !
  ! Integrated TOTAL electron and ion powers
  call cub_spline(expro_rmin(:)/r_min,expro_pow_e(:)*1e13,n_exp,r,p_e_in,n_r)
  call cub_spline(expro_rmin(:)/r_min,expro_pow_i(:)*1e13,n_exp,r,p_i_in,n_r)
  !
  ! Collisional exchange power
  call cub_spline(expro_rmin(:)/r_min,expro_pow_ei(:)*1e13,n_exp,r,p_exch_in,n_r)
  !
  ! (1a) Detailed powers for reactor simulation
  !
  ! Integrated fusion powers
  call cub_spline(expro_rmin(:)/r_min,expro_pow_e_fus(:)*1e13,n_exp,r,p_e_fus_in,n_r)
  call cub_spline(expro_rmin(:)/r_min,expro_pow_i_fus(:)*1e13,n_exp,r,p_i_fus_in,n_r)
  !
  ! Integrated radiated powers 
  call cub_spline(expro_rmin(:)/r_min,expro_pow_e_sync(:)*1e13,n_exp,r,p_sync_in,n_r)
  call cub_spline(expro_rmin(:)/r_min,expro_pow_e_brem(:)*1e13,n_exp,r,p_brem_in,n_r)
  call cub_spline(expro_rmin(:)/r_min,expro_pow_e_line(:)*1e13,n_exp,r,p_line_in,n_r)
  !  
  ! Integrated auxiliary heating powers (NB + RF + Ohmic)
  call cub_spline(expro_rmin(:)/r_min,expro_pow_e_aux(:)*1e13,n_exp,r,p_e_aux_in,n_r)
  call cub_spline(expro_rmin(:)/r_min,expro_pow_i_aux(:)*1e13,n_exp,r,p_i_aux_in,n_r)
  call cub_spline(expro_rmin(:)/r_min,expro_pow_e_ohmic(:)*1e13,n_exp,r,p_e_ohmic_in,n_r)
  
  ! Apply auxiliary power rescale
  ! 1. subtract off
  p_e_in = p_e_in-p_e_aux_in
  p_i_in = p_i_in-p_i_aux_in
  ! 2. Rescale
  p_e_aux_in = tgyro_input_paux_scale*p_e_aux_in
  p_i_aux_in = tgyro_input_paux_scale*p_i_aux_in 
  ! 3. Add back in
  p_e_in = p_e_in + p_e_aux_in 
  p_i_in = p_i_in + p_i_aux_in
  !
  ! (2) Particle flow -- convert to 1/s from MW/keV
  !
  call cub_spline(expro_rmin(:)/r_min,expro_flow_beam(:),n_exp,r,f_b_in,n_r)
  call cub_spline(expro_rmin(:)/r_min,expro_flow_wall(:),n_exp,r,f_w_in,n_r)
  !
  ! (3) Angular momentum flow -- convert to erg (dyne-cm) from N-m.
  !
  call cub_spline(expro_rmin(:)/r_min,expro_flow_mom(:)*1e7,n_exp,r,mf_in,n_r)
  !------------------------------------------------------------------------------------------

  !-----------------------------------------------------------------
  ! Parameters for EPED pedestal model 
  ! ** BEWARE: these are NOT all in CGS units.
  !
  allocate(exp_te(n_exp))
  allocate(exp_ne(n_exp))
  allocate(exp_ti(loc_n_ion,n_exp))
  allocate(exp_ni(loc_n_ion,n_exp))
  allocate(exp_w0(n_exp))
  allocate(exp_nu_exch(n_exp))
  ! exp_ne, exp_ni: [1/cm^3]
  exp_ne = expro_ne*1e13
  exp_ni(1:loc_n_ion,:) = expro_ni(1:loc_n_ion,:)*1e13
  ! exp_te, exp_ti: [eV]
  exp_te = expro_te*1e3
  exp_ti(1:loc_n_ion,:) = expro_ti(1:loc_n_ion,:)*1e3
  ! exp_w0 [1/s]
  exp_w0 = expro_w0
  exp_nu_exch = 0.0

  allocate(volp_exp(n_exp))
  volp_exp = expro_volp
  allocate(ptot_exp(n_exp))

  ! Compute pressure: ptot_exp
  call tgyro_pressure

  ! Volume average (p_ave)
  call tgyro_volume_ave(ptot_exp,expro_rmin,volp_exp,p_ave,n_exp)
  !
  allocate(rmin_exp(n_exp))
  rmin_exp = expro_rmin*100.0
  allocate(psi_exp(n_exp))
  ! Psi_norm
  psi_exp = expro_polflux/expro_polflux(n_exp)
  allocate(dpsidr_exp(n_exp))
  ! d (Psi_norm)/dr in units of 1/cm
  dpsidr_exp = expro_bunit*expro_rmin/expro_q/expro_polflux(n_exp)/100.0

  ! Check for sanity of psi_exp profile
  do i=2,n_exp
     if (abs(psi_exp(i)-psi_exp(i-1)) < 1e-8) then
        call tgyro_catch_error('ERROR: (tgyro_init_profiles) Poloidal flux profile (polflux) has two equal values.')
     endif
  enddo
  
  if (tgyro_ped_model > 1) then
     ! a [m]
     a_in = r_min
     ! Bt on axis (EFIT BCENTR) [T]
     if (abs(expro_bcentr) < 1e-10) then
        call tgyro_catch_error('expro_bcentr = 0.  Check input.gacode')
     endif
     bt_in = expro_bcentr
     ! Plasma current (EFIT CURRENT) [MA]
     if (abs(expro_current) < 1e-10) then
        call tgyro_catch_error('expro_current = 0. Check input.gacode')
     endif
     ip_in = expro_current
     ! betan [%] = betat/In*100 where In = Ip/(a Bt) 
     betan_in = abs(( p_ave/(0.5*bt_in**2/mu_0) ) / ( ip_in/(a_in*bt_in) ) * 100.0)
     ! Triangularity [-]
     delta_in = expro_delta(n_exp-3)  
     ! Elongation [-]
     kappa_in = expro_kappa(n_exp-3) 
     ! Main ion mass [0.5*md]
     m_in = mi_vec(1)
     ! R0(a) [m]
     r_in = expro_rmaj(n_exp)
     !
     ! Pedestal density
     if (tgyro_neped < 0.0) then
        ! Set neped to ne(psi_0), where psi_0=-tgyro_neped
        x0(1) = -tgyro_neped
        call cub_spline(psi_exp,expro_ne(:),n_exp,x0,y0,1)
        tgyro_neped = y0(1)
     endif
     !
     ! Pedestal zeff
     if (tgyro_zeffped < 0.0) then
        ! Set zeffped to zeff(psi_0), where psi_0=-tgyro_zeffped
        x0(1) = -tgyro_zeffped
        call cub_spline(psi_exp,expro_z_eff(:),n_exp,x0,y0,1)
        tgyro_zeffped = y0(1)
     endif

     ! Pedestal density/temperature ratios 
     if (tgyro_ped_ratio < 0.0) then
        do i_ion=1,loc_n_ion
           x0(1) = -tgyro_ped_ratio

           call cub_spline(psi_exp,expro_ne(:),n_exp,x0,y0,1)
           tmp_ped = y0(1)
           call cub_spline(psi_exp,expro_ni(i_ion,:),n_exp,x0,y0,1)
           n_ratio(i_ion) = y0(1)/tmp_ped

           call cub_spline(psi_exp,expro_te(:),n_exp,x0,y0,1)
           tmp_ped = y0(1)
           call cub_spline(psi_exp,expro_ti(i_ion,:),n_exp,x0,y0,1)
           t_ratio(i_ion) = y0(1)/tmp_ped
        enddo
     endif

     ! Quantities used to compute (ne,ni,Te,Ti) from <n>, <T>.
     n_frac = 2.0/(1.0+sum(n_ratio(1:loc_n_ion)))
     t_frac = (1.0+sum(n_ratio(1:loc_n_ion)))/(1.0+sum(n_ratio(1:loc_n_ion)*t_ratio(1:loc_n_ion)))

  endif
  !-----------------------------------------------------------------

  ! Convert r_min to cm (from m):
  r_min = r_min*100.0

  !----------------------------------------------------------
  ! Restore dimension to r (cm)
  r = r*r_min
  !----------------------------------------------------------

  !----------------------------------------------------------
  ! Primitive quantity to evolve for rotation/Er evolution is 
  !
  ! f_rot [1/cm] = w0p/w0_norm
  !
  ! w0_norm = c_s/R_maj at r=0.
  !
  w0_norm = sqrt(k*te(1)/md)/r_maj(1)
  !
  f_rot(:) = w0p(:)/w0_norm
  !----------------------------------------------------------
  
  !----------------------------------------------------------
  if (loc_restart_flag == 1) then
     quasifix = 1
     call tgyro_restart
     w0p(:) = w0_norm*f_rot(:)
  else
     i_tran = 0
  endif
  !----------------------------------------------------------

  ! Axis boundary conditions
  call tgyro_init_profiles_axis

end subroutine tgyro_init_profiles


!------------------------------------------------------------
! Subroutine to compute z that corresponds to profile
!------------------------------------------------------------
subroutine math_zfind(n,p,r,z)

  implicit none

  integer, intent(in) :: n
  real, intent(inout) :: p(n)
  real, intent(in) :: r(n)
  real, intent(inout) :: z(n)
  
  real, dimension(n) :: rat

  integer :: i

  ! Assume zero gradient at r=0
  z(1) = 0.0

  ! Modify axis temperature for smooth profile
  p(1) = p(2)*exp(0.5*z(2)*(r(2)-r(1)))

  rat = log(p/p(n))
  do i=2,n
     z(i) = 2*(rat(i)-rat(i-1))/(r(i)-r(i-1))-z(i-1)
  enddo
  z = -z
  
end subroutine math_zfind
 
