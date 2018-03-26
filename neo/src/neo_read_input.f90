subroutine neo_read_input

  use neo_globals

  implicit none

  integer :: is
  integer :: stat
  character (len=1) :: cdummy

  open(unit=1,file=trim(path)//'input.neo.gen',status='old')
  read(1,*) n_energy
  read(1,*) n_xi
  read(1,*) n_theta
  read(1,*) n_radial
  read(1,*) matsz_scalefac
  read(1,*) rmin_1_in
  read(1,*) rmin_2_in
  read(1,*) rmaj_in
  read(1,*) silent_flag
  read(1,*) sim_model
  read(1,*) equilibrium_model
  read(1,*) collision_model
  read(1,*) profile_model
  read(1,*) profile_erad0_model
  read(1,*) profile_equilibrium_model
  read(1,*) ipccw_in
  read(1,*) btccw_in
  read(1,*) te_ade_in
  read(1,*) ne_ade_in
  read(1,*) dlntdre_ade_in
  read(1,*) dlnndre_ade_in

  read(1,*) rotation_model
  read(1,*) omega_rot_in
  read(1,*) omega_rot_deriv_in
  read(1,*) rbf_dir

  read(1,*) spitzer_model
  read(1,*) epar0_spitzer

  read(1,*) coll_uncoupledei_model
  read(1,*) coll_uncoupledaniso_model

  read(1,*) n_species

  read(1,*) nu_1_in

  do is=1,11
     read(1,*) z_in(is)
     read(1,*) mass_in(is)
     read(1,*) dens_in(is)
     read(1,*) temp_in(is)
     read(1,*) dlnndr_in(is)
     read(1,*) dlntdr_in(is)
     read(1,*) aniso_model_in(is)
     read(1,*) temp_para_in(is)
     read(1,*) dlntdr_para_in(is)
     read(1,*) temp_perp_in(is)
     read(1,*) dlntdr_perp_in(is)
     read(1,*) profile_dlnndr_scale(is)
     read(1,*) profile_dlntdr_scale(is)
  enddo

  read(1,*) dphi0dr_in
  read(1,*) epar0_in
  read(1,*) q_in
  read(1,*) rho_in
  read(1,*) shear_in
  read(1,*) shift_in     
  read(1,*) kappa_in    
  read(1,*) s_kappa_in   
  read(1,*) delta_in        
  read(1,*) s_delta_in 
  read(1,*) zeta_in        
  read(1,*) s_zeta_in
  read(1,*) zmag_in        
  read(1,*) s_zmag_in
  read(1,*) beta_star_in
  read(1,*) profile_delta_scale
  read(1,*) profile_zeta_scale
  read(1,*) profile_zmag_scale

  read(1,*) subroutine_flag

  read(1,*) threed_model
  read(1,*) threed_exb_model
  read(1,*) threed_exb_dphi0dr
  read(1,*) threed_drift_model
  read(1,*) laguerre_method
  read(1,*) write_cmoments_flag
  read(1,*) threed_hyperxi

  close(1)

  ! GEO fourier coefficients are not yet available to read-in
  geo_ny_in = 0
  geo_yin_in(:,:) = 0.0
  if (subroutine_flag == 0 .and. equilibrium_model == 3) then
     open(unit=1,file=trim(path)//'input.geo',status='old')
     ! header skip
     do
        read(1,'(a)') cdummy
        if (cdummy /= '#') exit
     enddo
     backspace 1
     ! n_fourier
     read(1,*) geo_ny_in
     ! fourier coefficients
     read(1,*) geo_yin_in(:,0:geo_ny_in)
     close(1)
  endif

  ! 3D equilibrium from LE3:
  if (threed_model == 1) then

     open(unit=1,file=trim(path)//'out.le3.geoscalar',status='old',iostat=stat)
     if (stat /= 0) then
        call neo_error('ERROR: (NEO) le3 files not available')
     else
        read(1,*) n_tptheta
        read(1,*) n_tpvarphi
        read(1,*) tpmatsize
        read(1,*) indx_c00
        close(1)
     endif

  endif

end subroutine neo_read_input
