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
  use EXPRO_interface

  implicit none

  integer :: i_ion
  integer :: i
  integer :: n
  real :: arho
  real :: tmp_ped
  real :: p_ave
  real :: x0(1),y0(1)

  !------------------------------------------------------
  ! PHYSICAL CONSTANTS
  !
  pi = 4.0*atan(1.0)
  !
  e_alpha = 3.5e6*1.6022e-12 ! eV*(erg/eV)
  e       = 4.8032e-10 ! statcoul
  k       = 1.6022e-12 ! erg/eV
  me      = 9.1094e-28 ! g
  mp      = 1.6726e-24 ! g
  malpha  = 4*mp       ! g
  c       = 2.9979e10  ! cm/s
  !
  mu_0    = 4*pi*1e-7  ! N/A^2
  !------------------------------------------------------

  !------------------------------------------------------
  ! Convert dimensionless mass to grams.
  mi(:) = mi_vec(:)*mp
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
  if (tgyro_rmin > 0.0) then

     r(1) = 0.0
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
  ! Radius where profiles will be matched.
  !
  i_bc = n_r-loc_bc_offset
  !----------------------------------------------

  EXPRO_ctrl_n_ion = loc_n_ion
  EXPRO_ctrl_quasineutral_flag = 0
  EXPRO_ctrl_z = 0.0
  EXPRO_ctrl_z(1:loc_n_ion) = zi_vec(1:loc_n_ion)
  EXPRO_ctrl_numeq_flag = loc_num_equil_flag

  call EXPRO_palloc(MPI_COMM_WORLD,'./',1) 
  call EXPRO_pread

  shot = EXPRO_shot

  n_exp = EXPRO_n_exp

  ! r_min in m:
  r_min = EXPRO_rmin(n_exp)

  ! arho in cm
  arho  = 100*EXPRO_arho

  ! b_ref in Gauss
  b_ref = 1e4*EXPRO_b_ref

  ! Aspect ratio
  aspect_rat = EXPRO_rmaj(n_exp)/EXPRO_rmin(n_exp)

  !------------------------------------------------------------------------------------------
  ! Direct input of simple profiles:
  !
  if (tgyro_use_rho == 1) then
     ! Equally-spaced in rho
     call cub_spline(EXPRO_rho(:),EXPRO_rmin(:)/r_min,n_exp,rho,r,n_r)
  else  
     ! Equally-spaced in r (default)
     call cub_spline(EXPRO_rmin(:)/r_min,EXPRO_rho(:),n_exp,r,rho,n_r)
  endif
  call cub_spline(EXPRO_rmin(:)/r_min,EXPRO_q(:),n_exp,r,q,n_r)
  call cub_spline(EXPRO_rmin(:)/r_min,EXPRO_s(:),n_exp,r,s,n_r)
  call cub_spline(EXPRO_rmin(:)/r_min,EXPRO_kappa(:),n_exp,r,kappa,n_r)
  call cub_spline(EXPRO_rmin(:)/r_min,EXPRO_delta(:),n_exp,r,delta,n_r)
  call cub_spline(EXPRO_rmin(:)/r_min,EXPRO_skappa(:),n_exp,r,s_kappa,n_r)
  call cub_spline(EXPRO_rmin(:)/r_min,EXPRO_sdelta(:),n_exp,r,s_delta,n_r)
  call cub_spline(EXPRO_rmin(:)/r_min,EXPRO_drmaj(:),n_exp,r,shift,n_r)
  call cub_spline(EXPRO_rmin(:)/r_min,EXPRO_zmag(:),n_exp,r,zmag,n_r)
  call cub_spline(EXPRO_rmin(:)/r_min,EXPRO_dzmag(:),n_exp,r,dzmag,n_r)
  call cub_spline(EXPRO_rmin(:)/r_min,EXPRO_zeta(:),n_exp,r,zeta,n_r)
  call cub_spline(EXPRO_rmin(:)/r_min,EXPRO_szeta(:),n_exp,r,s_zeta,n_r)

  ! Convert ptot to Ba from Pascals (1 Pa = 10 Ba)
  call cub_spline(EXPRO_rmin(:)/r_min,EXPRO_ptot(:)*10.0,n_exp,r,ptot,n_r)

  ! Convert V and dV/dr from m^3 to cm^3
  call cub_spline(EXPRO_rmin(:)/r_min,EXPRO_vol(:)*1e6,n_exp,r,vol,n_r)
  call cub_spline(EXPRO_rmin(:)/r_min,EXPRO_volp(:)*1e4,n_exp,r,volp,n_r)

  call cub_spline(EXPRO_rmin(:)/r_min,EXPRO_ave_grad_r(:),n_exp,r,ave_grad_r,n_r)

  ! Convert B to Gauss (from T):
  call cub_spline(EXPRO_rmin(:)/r_min,1e4*EXPRO_bunit(:),n_exp,r,b_unit,n_r)

  ! Convert r_maj to cm (from m):
  call cub_spline(EXPRO_rmin(:)/r_min,100*EXPRO_rmaj(:),n_exp,r,r_maj,n_r)

  ! Convert T to eV (from keV) and length to cm (from m):
  call cub_spline(EXPRO_rmin(:)/r_min,1e3*EXPRO_te(:),n_exp,r,te,n_r)
  call cub_spline(EXPRO_rmin(:)/r_min,1e13*EXPRO_ne(:),n_exp,r,ne,n_r)
  call cub_spline(EXPRO_rmin(:)/r_min,EXPRO_dlnptotdr(:)/100.0,n_exp,r,dlnptotdr,n_r)
  call cub_spline(EXPRO_rmin(:)/r_min,EXPRO_dlnnedr(:)/100.0,n_exp,r,dlnnedr,n_r)
  call cub_spline(EXPRO_rmin(:)/r_min,EXPRO_dlntedr(:)/100.0,n_exp,r,dlntedr,n_r)
  do i_ion=1,loc_n_ion
     call cub_spline(EXPRO_rmin(:)/r_min,1e3*EXPRO_ti(i_ion,:),n_exp,r,ti(i_ion,:),n_r)
     call cub_spline(EXPRO_rmin(:)/r_min,1e13*EXPRO_ni(i_ion,:),n_exp,r,ni(i_ion,:),n_r)
     call cub_spline(EXPRO_rmin(:)/r_min,EXPRO_dlntidr(i_ion,:)/100.0,n_exp,r,dlntidr(i_ion,:),n_r)
     call cub_spline(EXPRO_rmin(:)/r_min,EXPRO_dlnnidr(i_ion,:)/100.0,n_exp,r,dlnnidr(i_ion,:),n_r)
     ! Define default ratios (these will change if tgyro_ped_ratio < 0.0)
     n_ratio(i_ion) = ni(i_ion,n_r)/ne(n_r)
     t_ratio(i_ion) = ti(i_ion,n_r)/te(n_r)
  enddo

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
  ! Helium ash option
  !
  i_ash = 0
  do i=1,loc_n_ion
     if (nint(zi_vec(i)) == 2 .and. nint(mi_vec(i)) == 4 .and. therm_flag(i) == 1) then
        i_ash = i
     endif
  enddo
  !------------------------------------------------------------------------------------------

  !------------------------------------------------------------------------------------------
  ! Rotation and rotation shear:
  !
  ! w0 (1/s)
  call cub_spline(EXPRO_rmin(:)/r_min,EXPRO_w0(:),n_exp,r,w0,n_r)
  ! w0p = d(w0)/dr (1/s/cm)
  call cub_spline(EXPRO_rmin(:)/r_min,EXPRO_w0p(:)/100.0,n_exp,r,w0p,n_r)
  !------------------------------------------------------------------------------------------

  !------------------------------------------------------------------------------------------
  ! Field orientation
  !
  ! signb = -btccw          OR     btccw = -signb  
  ! signq = ipccw*btccw            ipccw = -signb*signq
  !
  signb = EXPRO_signb
  signq = EXPRO_signq
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
     ! Set based on data in input.profiles
     call cub_spline(EXPRO_rmin(:)/r_min,EXPRO_z_eff(:),n_exp,r,z_eff,n_r)
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
  call cub_spline(EXPRO_rmin(:)/r_min,EXPRO_pow_e(:)*1e13,n_exp,r,p_e_in,n_r)
  call cub_spline(EXPRO_rmin(:)/r_min,EXPRO_pow_i(:)*1e13,n_exp,r,p_i_in,n_r)
  !
  ! Collisional exchange power
  call cub_spline(EXPRO_rmin(:)/r_min,EXPRO_pow_ei(:)*1e13,n_exp,r,p_exch_in,n_r)
  !
  ! (1a) Detailed powers for reactor simulation
  !
  ! Integrated fusion powers
  call cub_spline(EXPRO_rmin(:)/r_min,EXPRO_pow_e_fus(:)*1e13,n_exp,r,p_e_fus_in,n_r)
  call cub_spline(EXPRO_rmin(:)/r_min,EXPRO_pow_i_fus(:)*1e13,n_exp,r,p_i_fus_in,n_r)
  !
  ! Integrated radiated powers 
  call cub_spline(EXPRO_rmin(:)/r_min,EXPRO_pow_e_sync(:)*1e13,n_exp,r,p_sync_in,n_r)
  call cub_spline(EXPRO_rmin(:)/r_min,EXPRO_pow_e_brem(:)*1e13,n_exp,r,p_brem_in,n_r)
  call cub_spline(EXPRO_rmin(:)/r_min,EXPRO_pow_e_line(:)*1e13,n_exp,r,p_line_in,n_r)
  !  
  ! Integrated auxiliary heating powers (NB + RF + Ohmic)
  call cub_spline(EXPRO_rmin(:)/r_min,EXPRO_pow_e_aux(:)*1e13,n_exp,r,p_e_aux_in,n_r)
  call cub_spline(EXPRO_rmin(:)/r_min,EXPRO_pow_i_aux(:)*1e13,n_exp,r,p_i_aux_in,n_r)

  ! Apply auxiliary power rescale
  p_e_in = tgyro_input_paux_scale*p_e_aux_in + (p_e_in-p_e_aux_in)
  p_i_in = tgyro_input_paux_scale*p_i_aux_in + (p_i_in-p_i_aux_in)
  !
  ! (2) Particle flow -- convert to 1/s from MW/keV
  !
  call cub_spline(EXPRO_rmin(:)/r_min,EXPRO_flow_beam(:)*1e22/1.6022,n_exp,r,f_b_in,n_r)
  call cub_spline(EXPRO_rmin(:)/r_min,EXPRO_flow_wall(:)*1e22/1.6022,n_exp,r,f_w_in,n_r)
  !
  ! (3) Angular momentum flow -- convert to erg (dyne-cm) from N-m.
  !
  call cub_spline(EXPRO_rmin(:)/r_min,EXPRO_flow_mom(:)*1e7,n_exp,r,mf_in,n_r)
  !------------------------------------------------------------------------------------------

  !------------------------------------------------------------------------------------------
  ! Fourier coefficients for plasma shape
  if (EXPRO_nfourier > 0) then

     if (i_proc_global == 0) then
        open(unit=1,file=trim(runfile),position='append')
        write(1,*) 'INFO: (TGYRO) Passing input.profiles.geo information to components'
        write(1,*)
        close(1)
     endif

     n_fourier_geo = EXPRO_nfourier

     do n=0,n_fourier_geo
        do i=1,4  

           ! aR_n = EXPRO_geo(1,n,:)
           ! bR_n = EXPRO_geo(2,n,:)
           ! aZ_n = EXPRO_geo(3,n,:)
           ! bZ_n = EXPRO_geo(4,n,:)
           ! d(aR_n)/dr
           ! d(bR_n)/dr
           ! d(aZ_n)/dr
           ! d(bZ_n)/dr

           call cub_spline(EXPRO_rmin(:)/r_min,EXPRO_geo(i,n,:)/r_min,&
                n_exp,r,a_fourier_geo(i,n,:),n_r)
           call cub_spline(EXPRO_rmin(:)/r_min,EXPRO_dgeo(i,n,:),&
                n_exp,r,a_fourier_geo(i+4,n,:),n_r)

        enddo
     enddo

  else

     ! Numerical equilibrium not available

     n_fourier_geo      = 0
     loc_num_equil_flag = 0

  endif
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
  ! exp_ne, exp_ni: [1/cm^3]
  exp_ne = EXPRO_ne*1e13
  exp_ni(1:loc_n_ion,:) = EXPRO_ni(1:loc_n_ion,:)*1e13
  ! exp_te, exp_ti: [eV]
  exp_te = EXPRO_te*1e3
  exp_ti(1:loc_n_ion,:) = EXPRO_ti(1:loc_n_ion,:)*1e3
  ! exp_w0 [1/s]
  exp_w0 = EXPRO_w0

  allocate(volp_exp(n_exp))
  volp_exp = EXPRO_volp
  allocate(ptot_exp(n_exp))

  ! Compute pressure: ptot_exp
  call tgyro_pressure

  ! Volume average (p_ave)
  call tgyro_volume_ave(ptot_exp,EXPRO_rmin,volp_exp,p_ave,n_exp)
  !
  allocate(rmin_exp(n_exp))
  rmin_exp = EXPRO_rmin*100.0
  allocate(psi_exp(n_exp))
  ! Psi_norm
  psi_exp = EXPRO_polflux/EXPRO_polflux(n_exp)
  allocate(dpsidr_exp(n_exp))
  ! d (Psi_norm)/dr in units of 1/cm
  dpsidr_exp = EXPRO_bunit*EXPRO_rmin/EXPRO_q/EXPRO_polflux(n_exp)/100.0

  if (tgyro_ped_model > 1) then
     ! a [m]
     a_in = r_min
     ! Bt on axis [T]
     if (EXPRO_rvbv == 0.0) then
        call tgyro_catch_error('EXPRO_rvbv = 0.  Check input.profiles')
     endif
     bt_in = EXPRO_rvbv/EXPRO_rmaj(n_exp)
     ! Plasma current Ip [Ma]
     if (EXPRO_ip_exp == 0.0) then
        call tgyro_catch_error('EXPRO_ip_exp = 0. Check input.profiles')
     endif
     ip_in = 1e-6*EXPRO_ip_exp
     ! betan [%] = betat/In*100 where In = Ip/(a Bt) 
     betan_in = ( p_ave/(0.5*bt_in**2/mu_0) ) / ( ip_in/(a_in*bt_in) ) * 100.0
     ! Triangularity [-]
     delta_in = EXPRO_delta(n_exp-3)  
     ! Elongation [-]
     kappa_in = EXPRO_kappa(n_exp-3) 
     ! Main ion mass [mp]
     m_in = mi_vec(1)
     ! R0(a) [m]
     r_in = EXPRO_rmaj(n_exp)
     !
     ! Pedestal density
     if (tgyro_neped < 0.0) then
        ! Set neped to ne(psi_0), where psi_0=-tgyro_neped
        x0(1) = -tgyro_neped
        call cub_spline(psi_exp,EXPRO_ne(:),n_exp,x0,y0,1)
        tgyro_neped = y0(1)
     endif
     !
     ! Pedestal zeff
     if (tgyro_zeffped < 0.0) then
        ! Set zeffped to zeff(psi_0), where psi_0=-tgyro_zeffped
        x0(1) = -tgyro_zeffped
        call cub_spline(psi_exp,EXPRO_z_eff(:),n_exp,x0,y0,1)
        tgyro_zeffped = y0(1)
     endif

     ! Pedestal density/temperature ratios 
     if (tgyro_ped_ratio < 0.0) then
        do i_ion=1,loc_n_ion
           x0(1) = -tgyro_ped_ratio

           call cub_spline(psi_exp,EXPRO_ne(:),n_exp,x0,y0,1)
           tmp_ped = y0(1)
           call cub_spline(psi_exp,EXPRO_ni(i_ion,:),n_exp,x0,y0,1)
           n_ratio(i_ion) = y0(1)/tmp_ped

           call cub_spline(psi_exp,EXPRO_te(:),n_exp,x0,y0,1)
           tmp_ped = y0(1)
           call cub_spline(psi_exp,EXPRO_ti(i_ion,:),n_exp,x0,y0,1)
           t_ratio(i_ion) = y0(1)/tmp_ped
        enddo
     endif

     ! Quantities used to compute (ne,ni,Te,Ti) from <n>, <T>.
     n_frac = 2.0/(1.0+sum(n_ratio(1:loc_n_ion)))
     t_frac = (1.0+sum(n_ratio(1:loc_n_ion)))/(1.0+sum(n_ratio(1:loc_n_ion)*t_ratio(1:loc_n_ion)))

  endif
  !-----------------------------------------------------------------

  call EXPRO_palloc(MPI_COMM_WORLD,'./',0)

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
  w0_norm = sqrt(k*te(1)/mi(1))/r_maj(1)
  !
  f_rot(:) = w0p(:)/w0_norm
  !----------------------------------------------------------

  !----------------------------------------------------------
  if (loc_restart_flag == 1) then
     call tgyro_restart
     w0p(:) = w0_norm*f_rot(:)
  else
     i_tran = 0
  endif
  !----------------------------------------------------------

  ! Axis boundary conditions
  call tgyro_init_profiles_axis


end subroutine tgyro_init_profiles
