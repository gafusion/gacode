subroutine neo_check

  use neo_globals

  implicit none

  integer :: ir, is

  if(silent_flag == 0 .and. i_proc == 0) then
     open(unit=io_neoout,file=trim(path)//runfile_neoout,&
          status='old',position='append')
  endif

  !-----------------------------------------------------------
  ! Grid parameter checks
  !
  if (modulo(n_theta,2) == 0) then 
     call neo_error('ERROR: (NEO) n_theta must be odd')
     return
  endif
  !
  if(n_species > 6) then
     call neo_error('ERROR: (NEO) max n_species is 6')
     return
  endif

  if(rho_in < 0) then
     call neo_error('ERROR: (NEO) rho_unit must be positive')
     return
  endif

  !-----------------------------------------------------------
  if(silent_flag == 0 .and. i_proc == 0) then
     write(io_neoout,*) 'CONTROL PARAMETERS'
     write(io_neoout,*) '------------------'
  end if
  !------------------------------------------------------------

  ! Simulation model
  select case (sim_model)
  case(0)
     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,30) 'sim_model','THEORY'
     endif
  case(1)
     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,30) 'sim_model','NUMERICAL (with theory and nclass)'
     endif
  case(2)
     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,30) 'sim_model','NUMERICAL (with theory)'
     endif
  case default   
     call neo_error('ERROR: (NEO) invalid sim_model')
     return
  end select

  ! Collision model
  !
  select case (collision_model)  

  case (1) 

     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,30) 'collision_model','CONNOR'
     endif

  case (2) 

     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,30) 'collision_model','REDUCED HIRSHMAN-SIGMAR'
     end if

  case (3) 

     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,30) 'collision_model','FULL HIRSHMAN-SIGMAR'
     endif

  case (4) 

     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,30) 'collision_model','FULL LINEARIZED FOKKER-PLANCK'
     end if

  case (5) 

     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,30) 'collision_model','FULL LINEARIZED FOKKER-PLANCK WITH AD-HOC FIELD PARTICLE TERMS' 
     end if

  case default

     call neo_error('ERROR: (NEO) invalid collision_model')
     return

  end select

  select case (coll_uncoupledei_model) 

  case (0) 
     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,30) 'collision_model','Full e-i collisional coupling'
     endif

  case (1)
     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,30) 'collision_model','Reduced e-i collisional coupling'
     endif

  case (2)
     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,30) 'collision_model','Reduced e-i collisional coupling'
     endif

  case default

     call neo_error('ERROR: (NEO) invalid coll_uncoupledei_model')
     return

  end select

  !------------------------------------------------------------

  !------------------------------------------------------------
  ! Equilibrium model
  !
  select case (equilibrium_model)  


  case (0) 

     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,30) 'equilibrium_model','S-ALPHA'
     end if

  case (1) 

     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,30) 'equilibrium_model','LARGE-ASPECT-RATIO'
     end if

  case (2) 

     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,30) 'equilibrium_model','MILLER'
     end if

  case (3) 

     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,30) 'equilibrium_model','GENERAL'
     end if

     if(geo_ny <= 0) then
        call neo_error('ERROR: (NEO) geometry coefficients missing')
        return
     endif

  case default

     call neo_error('ERROR: (NEO) equilibrium_model invalid')
     return

  end select
  !------------------------------------------------------------

  !------------------------------------------------------------
  ! Profile model
  !
  select case (profile_model)  
  case (1) 
     if(n_radial > 1) then
        call neo_error('ERROR: (NEO) profile_model=1 must be run with n_radial = 1')
        return
     endif
     ir=1
     do is=1,n_species
        if(dens(is,ir) <= 0.0) then
           call neo_error('ERROR: (NEO) density must be positive')
           return
        end if
        if(temp(is,ir) <= 0.0) then
           call neo_error('ERROR: (NEO) temperature must be positive')
           return
        end if
        if(nu(is,ir) <= 0.0) then
           call neo_error('ERROR: (NEO) collision frequency must be positive')
           return
        end if
        if(z(is) == 0.0) then
           call neo_error('ERROR: (NEO) charge must be non-zero')
           return
        end if
     enddo

     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,30) 'profile_model','LOCAL'
     end if

  case (2) 

     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,30) 'profile_model','GLOBAL PROFILE'
     end if

     select case (profile_erad0_model)
     case(0)
        if(silent_flag == 0 .and. i_proc == 0) then
           write(io_neoout,*) 'GLOBAL PROFILE profile_erad0_model: ERAD0 NOT INCLUDED'
        end if
     case (1)
        if(silent_flag == 0 .and. i_proc == 0) then
           write(io_neoout,*) 'GLOBAL PROFILE profile_erad0_model: ERAD0 INCLUDED'
        end if
     case default
        call neo_error('ERROR: (NEO) invalid profile_erad0_model')
        return
     end select

     select case (profile_equilibrium_model)
     case(0)
        if(silent_flag == 0 .and. i_proc == 0) then
           write(io_neoout,*) 'GLOBAL PROFILE profile_equilibrium_model: WITH S-ALPHA GEOMETRY'
        end if
     case (1)
        if(silent_flag == 0 .and. i_proc == 0) then
           write(io_neoout,*) 'GLOBAL PROFILE profile_equilibrium_model: WITH MILLER GEOMETRY'
           if(abs(profile_delta_scale-1.0) > epsilon(0.) ) then
              write(io_neoout,*) 'GLOBAL PROFILE profile_equilibrium_model: DELTA AND S_DELTA ARE RE-SCALED'
           endif
           if(abs(profile_zeta_scale-1.0) > epsilon(0.) ) then
              write(io_neoout,*) 'GLOBAL PROFILE profile_equilibrium_model: ZETA AND S_ZETA ARE RE-SCALED'
           endif
           if(abs(profile_zmag_scale-1.0) > epsilon(0.) ) then
              write(io_neoout,*) 'GLOBAL PROFILE profile_equilibrium_model: ZMAG AND S_MAG ARE RE-SCALED'
           endif
        end if
     case (2)
        if(silent_flag == 0 .and. i_proc == 0) then
           write(io_neoout,*) 'GLOBAL PROFILE profile_equilibrium_model: WITH GENERAL GEOMETRY'
        endif
     case default
        call neo_error('ERROR: (NEO) invalid profile_equilibrium_model')
        return
     end select

     if(n_species >= 1) then
        if(abs(profile_dlnndr_scale(1)-1.0) > epsilon(0.) ) then
           write(io_neoout,*) 'GLOBAL PROFILE profile_equilibrium_model: DLNNDR_1 IS RE-SCALED'
        endif
        if(abs(profile_dlntdr_scale(1)-1.0) > epsilon(0.) ) then
           write(io_neoout,*) 'GLOBAL PROFILE profile_equilibrium_model: DLNTDR_1 IS RE-SCALED'
        endif
     endif
     if(n_species >= 2) then
        if(abs(profile_dlnndr_scale(2)-1.0) > epsilon(0.) ) then
           write(io_neoout,*) 'GLOBAL PROFILE profile_equilibrium_model: DLNNDR_2 IS RE-SCALED'
        endif
        if(abs(profile_dlntdr_scale(2)-1.0) > epsilon(0.) ) then
           write(io_neoout,*) 'GLOBAL PROFILE profile_equilibrium_model: DLNTDR_2 IS RE-SCALED'
        endif
     endif
     if(n_species >= 3) then
        if(abs(profile_dlnndr_scale(3)-1.0) > epsilon(0.) ) then
           write(io_neoout,*) 'GLOBAL PROFILE profile_equilibrium_model: DLNNDR_3 IS RE-SCALED'
        endif
        if(abs(profile_dlntdr_scale(3)-1.0) > epsilon(0.) ) then
           write(io_neoout,*) 'GLOBAL PROFILE profile_equilibrium_model: DLNTDR_3 IS RE-SCALED'
        endif
     endif
     if(n_species >= 4) then
        if(abs(profile_dlnndr_scale(4)-1.0) > epsilon(0.) ) then
           write(io_neoout,*) 'GLOBAL PROFILE profile_equilibrium_model: DLNNDR_4 IS RE-SCALED'
        endif
        if(abs(profile_dlntdr_scale(4)-1.0) > epsilon(0.) ) then
           write(io_neoout,*) 'GLOBAL PROFILE profile_equilibrium_model: DLNTDR_4 IS RE-SCALED'
        endif
     endif
     if(n_species >= 5) then
        if(abs(profile_dlnndr_scale(5)-1.0) > epsilon(0.) ) then
           write(io_neoout,*) 'GLOBAL PROFILE profile_equilibrium_model: DLNNDR_5 IS RE-SCALED'
        endif
        if(abs(profile_dlntdr_scale(5)-1.0) > epsilon(0.) ) then
           write(io_neoout,*) 'GLOBAL PROFILE profile_equilibrium_model: DLNTDR_5 IS RE-SCALED'
        endif
     endif
     if(n_species >= 6) then
        if(abs(profile_dlnndr_scale(6)-1.0) > epsilon(0.) ) then
           write(io_neoout,*) 'GLOBAL PROFILE profile_equilibrium_model: DLNNDR_6 IS RE-SCALED'
        endif
        if(abs(profile_dlntdr_scale(6)-1.0) > epsilon(0.) ) then
           write(io_neoout,*) 'GLOBAL PROFILE profile_equilibrium_model: DLNTDR_6 IS RE-SCALED'
        endif
     endif

  case default

     call neo_error('ERROR: (NEO) invalid profile_model')
     return

  end select

  !-----------------------------------------------------------
  ! Sign of B checks
  if(sign_q > 0.0) then
     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,30) 'sign_q','POSITIVE'
     end if
  else
     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,30) 'sign_q','NEGATIVE'
     end if
  end if

  if(sign_bunit > 0.0) then
     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,30) 'sign_bunit','POSITIVE (BT CW)'
     end if
  else
     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,30) 'sign_bunit','NEGATIVE (BT CCW)'
     end if
  end if


  !-----------------------------------------------------------
  ! Rotation parameter checks
  !
  select case (rotation_model)     
  case (1)
     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,30) 'rotation model','ROTATION EFFECTS NOT INCLUDED'
     end if
  case (2)
     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,30) 'rotation model','ROTATION EFFECTS INCLUDED'
     end if
  case default
     call neo_error('ERROR: (NEO) invalid rotation_model')
     return
  end select

  select case (spitzer_model)
  case(0)
     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,30) 'spitzer_model:','OFF (NEOCLASSICAL TRANSPORT)'
     end if
  case (1)
     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,30) 'spitzer_model','ON (SOLVE SPITZER PROBLEM)'
     end if
  case default
     call neo_error('ERROR: (NEO) invalid spitzer_model')
     return
  end select

  select case(threed_model)
  case(0)
     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,30) 'threed_model','AXISYMMETRIC EQUILIBRIUM'
     end if
  case (1)
     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,30) 'threed_model','NON-AXISYMMETRIC EQUILIBRIUM (LE3)'
     end if

     select case(threed_exb_model)
     case(0)
        if(silent_flag == 0 .and. i_proc == 0) then
           write(io_neoout,30) 'threed_exb_model','NO HIGHER-ORDER V_EXB DRIFT'
        end if
     case (1)
        if(silent_flag == 0 .and. i_proc == 0) then
           write(io_neoout,30) 'threed_exb_model','HIGHER-ORDER V_EXB DRIFT'
        end if
     case default
        call neo_error('ERROR: (NEO) invalid threed_exb_model')
        return
     end select

  case default
     call neo_error('ERROR: (NEO) invalid threed_model')
     return
  end select

  !-----------------------------------------------------------
  ! Anisotropic species checks
  !
  select case (aniso_model)     
  case (1)
  case (2)
     if(profile_model == 2) then
        call neo_error('ERROR: (NEO) aniso_model not available with global profiles')
     endif
     if(rotation_model == 1) then
        call neo_error('ERROR: (NEO) aniso_model requires rotation_method=2')
     endif
     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,30) 'aniso model','ANISOTROPIC SPECIES INCLUDED'
     end if
  case default
     call neo_error('ERROR: (NEO) invalid aniso_model')
     return
  end select

  !------------------------------------------------------------

  if (silent_flag == 0 .and. i_proc == 0) then

     write(io_neoout,*)
     write(io_neoout,*) 'GRID DIMENSIONS'
     write(io_neoout,*) '---------------'
     write(io_neoout,10) 'n_radial',n_radial
     write(io_neoout,10) 'n_energy',n_energy
     write(io_neoout,10) 'n_xi',n_xi
     if (threed_model == 0) then
        write(io_neoout,10) 'n_theta',n_theta
     else
        write(io_neoout,10) 'n_theta (modes)',n_tptheta
        write(io_neoout,10) 'n_phi   (modes)',n_tpvarphi
        write(io_neoout,10) 'n_tp',tpmatsize
        write(io_neoout,"(t2,a,t19,': ',i5)") 'n_row',n_row
     endif

     do ir=1,n_radial
        write(io_neoout,*) 
        write(io_neoout,*) 'PHYSICS PARAMETERS'
        write(io_neoout,*) '------------------'
        write(io_neoout,20) 'r/R',r(ir)/rmaj(ir)
        write(io_neoout,20) 'dphi0/dr',dphi0dr(ir)
        write(io_neoout,20) 'omega_rot',omega_rot(ir)
        write(io_neoout,20) 'omega_rot_deriv',omega_rot_deriv(ir)
        write(io_neoout,20) 'q',q(ir)
        write(io_neoout,20) 's',shat(ir)
        write(io_neoout,20) 'shift',shift(ir)
        write(io_neoout,20) 'kappa',kappa(ir)
        write(io_neoout,20) 's_kappa',s_kappa(ir)
        write(io_neoout,20) 'delta',delta(ir)
        write(io_neoout,20) 's_delta',s_delta(ir)
        write(io_neoout,20) 'zeta',zeta(ir)
        write(io_neoout,20) 's_zeta',s_zeta(ir)
        write(io_neoout,20) 'zmag',zmag(ir)
        write(io_neoout,20) 's_zmag',s_zmag(ir)

        do is=1,n_species
           write(io_neoout,*) 
           write(io_neoout,'(t2,a,i1)') 'Species ',is
           write(io_neoout,*) '----------'
           write(io_neoout,10) 'Z',z(is)
           write(io_neoout,20) 'dens',dens(is,ir)
           write(io_neoout,20) 'temp',temp(is,ir)
           write(io_neoout,20) 'mass',mass(is)
           write(io_neoout,20) 'a/Ln',dlnndr(is,ir)
           write(io_neoout,20) 'a/LT',dlntdr(is,ir)
           write(io_neoout,20) 'nu',nu(is,ir)
        enddo
     enddo

     if (aniso_model == 2) then
        write(io_neoout,*) 
        write(io_neoout,10) 'Z_aniso',z_aniso
        write(io_neoout,20) 'dens',dens_aniso
        write(io_neoout,20) 'temp_parallel',temp_para_aniso
        write(io_neoout,20) 'temp_perp',temp_perp_aniso
        write(io_neoout,20) 'mass',mass_aniso
     endif

     write(io_neoout,*)

  endif

  if(silent_flag == 0 .and. i_proc == 0) then
     close(io_neoout)
  endif

10 format(t2,a,t19,': ',i3)
20 format(t2,a,t19,':',1pe12.4) 
30 format(t2,a,t22,': ',a)

end subroutine neo_check
